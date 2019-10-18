# coding=utf-8
# the above line is for python 2 compatibility
# write suite suite [1.0] in progress: inserted suite load, asym as a matrix, fft on matrix CHECK!!!
######################################################
# Gui tabs correspond to distinct gui methods with independes scopes and additional local methods
# gui attributes: 
#     entities that must be shared between tabs 
#     including variables passed to functions  outside mugui
###############################################################
#             Implementation of multi fit
#                 logical         logical        list of lists of musr2py instances
# type of fit   | self._global_ |self._single_ | len(self._the_runs_)
#-------------------------------------------------------
# single        | False         | True          | 1 (e.g. [[run234]] or [[run234,run235]] (the latter adds data of two runs)
# non global    | False         | False         | >1
# global (TODO) | True          | False         | >1
#-------------------------------------------------------
# asymmetry loads 
# single runs, both for single and for non global fit
# run suites, both for global fit and for multiplot
###############################################################
            # another dialog to explore
            #import tkinter
            #from tkinter import simpledialog
            #move = simpledialog.askstring("Pause","hit return when ready")
            #simpledialog.mainloop(0)


class mugui(object):
    '''
    The main method is gui, that designs a tabbed gui. Planned to be run in a Jupyter notebook. 

    Public methods correspond to the main tabs (fft, fit, setup, suite, ...). 
    Other public methods produce the fundamental data structure (asymmetry) 
    or common functions (plot). Their nested methods are not documented by sphynx.
    One tab handles (badly) the text output, mixing error messages and fit results.

    function __init__(self)

      * initiates an instance and a few attributes, 
      * launches the gui.
      * Use as follows::

             from mugui import mugui as MG
             MuJPy = MG() # instance is MuJPy

    '''

##########################
# INIT
##########################
    def __init__(self):
        '''
        * Initiates an instance and a few attributes, 
        * launches the gui.
        * Use a template notebook, basically as follows::

             from mugui import mugui as MG
             MuJPy = MG() # instance is MuJPy

        '''
        import numpy as np
        import os
        from scipy.constants import physical_constants as C
        from mujpy import __file__ as MuJPyName
        from IPython.display import display

# check where constants are needed and how to define them
        self.TauMu_mus = 2.1969811 # numbers are from Particle Data Group 2017
        self.TauPi_ns = 2.6033 # numbers are from Particle Data Group 2017
        self.gamma_Mu_MHzperT = 3.183345142*C['proton gyromag. ratio over 2 pi'][0]  # numbers are from Particle Data Group 2017
        self.gamma_e_MHzperT = C['electron gyromag. ratio over 2 pi'][0]
# end check constants

# initializations
        self.offset0 = 7 # initial value
        self.offset = [] # this way the first run load calculates get_totals with self.offset0
        self.firstbin = 0
        self.second_plateau = 100
        self.peakheight = 100000.
        self.peakwidth = 1.   # rough guesses for default
        self.histoLength = 7900 # initialize
        self.bin_range0 = '0,500' # initialize (for plots, counter inspection)
        self.nt0_run = []
        self._global_ = False
        self.thermo = 1 # sample thermometer is 1 on gps (check or adapt to other instruments)

        self.binwidth_ns = [] # this way the first call to asymmetry(_the_runs_) initializes self.time
        self.grouping = {'forward':np.array([1]),'backward':np.array([0])} # normal dict
        self._the_runs_ = []  # if self._the_runs_: is False, to check whether the handle is created
        self.first_t0plot = True
        self.fitargs= [] # initialize

# mujpy paths
        self.__path__ = os.path.dirname(MuJPyName)
        self.__logopath__ = os.path.join(self.__path__,"logo")
        self.__startuppath__ = os.getcwd() # working directory, in which, in case, to find mujpy_setup.pkl 
# mujpy layout
        self.button_color = 'lightgreen'
        self.button_color_off = 'gray'
        self.newtlogdir = True

        #####################################
        # actually produces the gui interface
        #####################################

        self.gui()
        self.output() # this calls output(self) that defines self._output_
        self.setup()
        self.suite()
        self.fit()
        self.fft()
        self.plots()

        #####################################
        # static figures 
        ##########################        self.fig_fit = [] # initialize to false, it will become a pyplot.subplots instance 
###########
        self.fig_fit = [] # initialize to false, it will become a pyplot.subplots instance 
        self.fig_fft = [] # initialize to false, it will become a pyplot.subplots instance
        self.fig_multiplot = []
        self.fig_counters = []
        self.fig_tlog = []
#        self.graph_fit = Toplevel()
#        canvas_fit = FigureCanvas(self.fig_fit, master=self.graph_fit)
#        canvas_fit.get_tk_widget().pack()
#        self.toolbar_fit = Navig(canvas_fit,self.graph_fit)
#        self.toolbar_fit.update()
#        self.graph_fit.withdraw()  # self.graph_fit.deiconify() brings it back up


        # self fft
        # self.tlog
        self.about()
        try:
            whereamI = get_ipython().__class__
            if not str(whereamI).find('erminal')+1:
                display(self.gui) # you are in a Jupyter notebook
            else:
                print(str(wheremI)) # you are in an ipython terminal
        except:
                print('Python test script') # other option?

        def _eval(string):
            '''
            yes, I know eval is evil, but
            mujpy users with jupyter have already full control
            of the machine, hence I do not care!
            **** (UNUSED! should go) *****
            BTW exec is used even by python for good: see input!
            '''
            try: 
                return eval(string)
            except Exception as e: 
                print(e)

##########################
# ABOUT
##########################
    def about(self):
        '''
        about tab:

        - a few infos (version and authors)

        '''
        from ipywidgets import Textarea, Layout

        _version = 'MuJPy          version '+'1.0' # increment while progressing
        _authors = '\n\n  Authors: Roberto De Renzi, Pietro Bonfà (*)'
        _blahblah = ('\n\n  A Python MuSR data analysis graphical interface.'+
                     '\n  Based on classes, designed for jupyter.'+
                     '\n  Released under the MIT licence')
        _pronounce = ('\n  See docs in ReadTheDocs'+
                      '\n  Pronounce it as mug + pie')
        _additional_credits_ = ('\n ---------------------\n (*) dynamic Kubo-Toyabe algorithm by G. Allodi\n     MuSR_td_PSI by A. Amato and A.-R. Raselli \n     acme algorithm code from NMRglue, by Jonathan J. Helmus')
        _about_text = _version+_blahblah+_pronounce+_authors+_additional_credits_
        _about_area = Textarea(value=_about_text,
                                   placeholder='Info on MuJPy',
                                   layout=Layout(width='100%',height='170px'),
                                   disabled=True)
        # now collect the handles of the three horizontal frames to the main fit window (see tabs_contents for index)
        self.mainwindow.children[6].children = [_about_area] # add the list of widget handles as the third tab, fit
       
##########################
# ASYMMETRY
##########################
    def asymmetry(self):   
        """
        * defines self.time, generates asymmetry end error without rebinning, 

        * all are 2D numpy arrays, shape of self.time is (1,nbins), of self.asymm, self.asyme are (nruns, nbins)

        * self._the_runs_ is a list of lists of musr2py (psi loading routine) instances 

        * inner list is for adding runs

        * outer list is suites of runs

        may be treated equally by a double for loop (single, no addition implies that k,j=0,0)
        
        - returns 0 if ok and -1 if an error is encountered 

        """ 
        import numpy as np

        # no checks, consistency in binWidth and numberHisto etc are done with run loading
        self.numberHisto = self._the_runs_[0][0].get_numberHisto_int()
        self.histoLength = self._the_runs_[0][0].get_histoLength_bin() - self.nt0.max() - self.offset.value # max available bins on all histos
        self.firstrun  = True
        self.binwidth_ns = self._the_runs_[0][0].get_binWidth_ns() 
        time = (np.arange(self.histoLength) + self.offset.value +
                     np.mean(self.dt0 [np.append(self.grouping['forward'],self.grouping['backward'])] ) 
                     )*self.binwidth_ns/1000. # in microseconds, 1D np.array
        self.time = np.array([time]) # in microseconds, 2D np.array
    ##################################################################################################
    # Time definition: 
    # 1) Assume the prompt is entirely in bin self.nt0. (python convention, the bin index is 0,...,n,... 
    # The content of bin self.nt0 will be the t=0 value for this case and self.dt0 = 0.
    # The center of bin self.nt0 will correspond to time t = 0, time = (n-self.nt0 + self.offset.value + self.dt0)*mufit.binWidth_ns/1000.
    # 2) Assume the prompt is equally distributed between n and n+1. Then self.nt0 = n and self.dt0 = 0.5, the same formula applies
    # 3) Assume the prompt is 0.45 in n and 0.55 in n+1. Then self.nt0 = n+1 and self.dt0 = -0.45, the same formula applies.
    ##################################################################################################
        # calculate asymmetry in y and error in ey
        for k, runs in enumerate(self._the_runs_):
            yforw = np.zeros(time.shape[0]) # counts with background substraction
            cforw = np.zeros(time.shape[0]) # pure counts for Poisson errors
            ybackw = np.zeros(time.shape[0]) # counts with background substraction
            cbackw = np.zeros(time.shape[0]) # pure counts for Poisson errors
            for j, run in enumerate(runs): 

                for counter in self.grouping['forward']:
                    n1, n2 = self.nt0[counter]+self.offset.value, self.nt0[counter]+self.offset.value+self.histoLength
                    histo = run.get_histo_array_int(counter)
                    background = np.mean(histo[self.firstbin:self.lastbin])
                    yforw += histo[n1:n2]-background
                    cforw += histo[n1:n2]

                for counter in self.grouping['backward']:
                    n1, n2 = self.nt0[counter]+self.offset.value, self.nt0[counter]+self.offset.value+self.histoLength
                    histo = run.get_histo_array_int(counter)
                    background = np.mean(histo[self.firstbin:self.lastbin])
                    ybackw += histo[n1:n2]-background
                    cbackw += histo[n1:n2]

            yplus = yforw + self.alpha.value*ybackw
            x = np.exp(-time/self.TauMu_mus)
            enn0 = np.polyfit(x,yplus,1)
            enn0 = enn0[0] # initial rate per ns
            y = (yforw-self.alpha.value*ybackw)/enn0*np.exp(time/self.TauMu_mus)  # since self.time is an np.arange, this is a numpy array
            ey = np.sqrt(cforw + self.alpha.value**2*cbackw)*np.exp(time/self.TauMu_mus)/enn0 # idem
            ey[np.where(ey==0)] = 1 # substitute zero with one in ey
            if self._single_:  # len(self._the_runs_)==1 and k=0 
                self.asymm = np.array([y]) # 2D np.array
                self.asyme = np.array([ey]) # 2D np.array
                self.nrun = [runs[0].get_runNumber_int()]
            else: # the first call of the suite the master, resets binwidth_ns, hence self.firstrun=True
                if self.firstrun:
                    self.asymm = y # 1D np.array
                    self.asyme = ey # idem
                    self.firstrun = False
                    self.nrun = [runs[0].get_runNumber_int()]
                else:
                    self.asymm = np.row_stack((self.asymm, y)) # columns are times, rows are successive runs (for multiplot and global)
                    self.asyme = np.row_stack((self.asyme, ey))
                    self.nrun.append(runs[0].get_runNumber_int())  # this is a list
                ######################################################
                # self.nrun contains only the first run in case of run addition
                # used by save_fit (in file name), 
                #         write_csv (first item is run number)
                #         animate (multiplot) 
                ######################################################
        # self.console('1st, last = {},{} bin for bckg subtraction'.format(self.firstbin,self.lastbin)) #dbg

#########
# CONSOLE
#########

    def console(self,string,end = '\n'):
        '''
        printed output method::

            if self._output_ == self._outputtab_ prints on tab 
                   and sets self.mainwindow.selected_index = 3
            otherwise writes on terminal 

        ''' 
        if self._output_ == self._outputtab_:
            with self._output_:
                print(string)
                self.mainwindow.selected_index = 3
        else:
            p = open(self._output_,'w') 
            p.write(''.join([string,end]))
            p.close()

###################
# CONSOLE TRACEBACK
###################

    def console_trcbk(self):
        '''
        method to print traceback in console

        ''' 
        import sys, traceback
        if self._output_ == self._outputtab_:
            with self._output_:
                raise 
                self.mainwindow.selected_index = 3
        else:
            p = open(self._output_,'w') 
            formatted_lines = traceback.format_exc().splitlines()
            for string in formatted_lines:
                p.write(''.join([string,'\n']))
            p.close()
        
##########################
# CREATE_RUNDICT
##########################
    def create_rundict(self,k=0):
        '''
        creates a dictionary rundict to identify and compare runs;
        refactored for adding runs

        '''
        rundict={}
        instrument = self.filespecs[0].value.split('_')[2] # valid for psi with standard names 'deltat_tdc_gpd_xxxx.bin'
        for j,run in enumerate(self._the_runs_[k]): # more than one: add sequence 
            rundict0 = {}
            rundict0.update({'nhist':run.get_numberHisto_int()})
            rundict0.update({'histolen':run.get_histoLength_bin()})
            rundict0.update({'binwidth':run.get_binWidth_ns()})
            rundict0.update({'instrument':instrument})
            if not rundict: # rundict contains only the first run of an add sequence (ok also for no add)
                rundict = rundict0
            elif rundict0!=rundict: # trying to add runs with different nhist, histolen, binwidth, instrument?
                rundict.update({'error':run.get_runNumber_int()})
                break     
        rundict.update({'nrun':self._the_runs_[k][0].get_runNumber_int()}) 
        rundict.update({'date':self._the_runs_[k][0].get_timeStart_vector()})
        return rundict

##########################
# GUI
##########################
    def gui(self):
        '''
        Main gui layout. Based on ipywidgets.
        Executed only once,
        it designs 
        
        * an external frame
        
        * the logo and title header
        
        * the tab structure. Tabs mostly correspond to public mugui methods. Their nested methods cannot be documented by spynx so they are replicated under each public method

        At the end (Araba.Phoenix) the gui method redefines self.gui 
        as a Vbox named 'whole', containing the entire gui structure

        '''
        from ipywidgets import Image, Text, Layout, HBox, Output, VBox, Tab
        import os                      

        file = open(os.path.join(self.__logopath__,"logo.png"), "rb")
        image = file.read()
        logo = Image(value=image,format='png',width=132,height=132)
        self.title = Text(description='run title', value='none yet',layout=Layout(width='55%'),disabled=True)
        self._the_runs_display = Text(description='run number',value='no run',layout=Layout(width='45%'),disabled=True)
        title_content = [self._the_runs_display, self.title]
        titlerow = HBox(description='Title')
        titlerow.children = title_content
        comment_box = HBox(description='comment',layout=Layout(width='100%'))
        self.comment_handles = [Text(description='Comment',layout=Layout(width='46%'),disabled=True),
                                Text(description='Start date',layout=Layout(width='27%'),disabled=True),
                                Text(description='Stop date',layout=Layout(width='27%'),disabled=True)]
        comment_box.children = self.comment_handles
        counts = ['Total counts', 'Group counts','ns/bin'] # needs an HBox with three Text blocks
        self.totalcounts = Text(value='0',description='Total counts',layout=Layout(width='30%'),disabled=True)
        self.groupcounts = Text(value='0',description='Group counts',layout=Layout(width='30%'),disabled=True)
        self.nsbin = Text(description='ns/bin',layout=Layout(width='20%'),disabled=True)
        self.maxbin = Text(description='Last bin',layout=Layout(width='20%'),disabled=True)
        secondrow = HBox(description='counts',layout=Layout(width='100%'))
        secondrow.children = [self.totalcounts, self.groupcounts, self.nsbin, self.maxbin]
        titlewindow = VBox()
        titlewindow_content = [titlerow, comment_box, secondrow] # ,thirdrow (moved to 4th tab)
        titlewindow.children = titlewindow_content
        titlelogowindow = HBox()
        titlelogowindow_content = [logo, titlewindow]
        titlelogowindow.children = titlelogowindow_content

        # main layout: tabs
        tabs_contents = ['setup', 'suite', 'fit', 'output', 'fft', 'plots', 'about']
        tabs = [VBox(description=name,layout=Layout(border='solid')) for name in tabs_contents]
        self.mainwindow = Tab(children = tabs,layout=Layout(width='99.8%')) # '99.6%' works

        self.mainwindow.selected_index = 0 # to stipulate that the first display is on tab 0, setup
        for i in range(len(tabs_contents)):
            self.mainwindow.set_title(i, tabs_contents[i])
        # Araba.Phoenix:
        self.gui = VBox(description='whole',layout=Layout(width='100%'))
        self.gui.children = [titlelogowindow, self.mainwindow]

##########################
# FFT
##########################
    def fft(self): 
        '''
        fft tab of mugui.

        on_fft_request(b) performs fft and plot (WARNING)
        * two options: (partial) residues or full asymmetry
        * two modes: real amplitude or power
        * vectorized: range(len(self.fitargs)) is (0,1) or (0,n>1) for single or suite  

        '''
        def on_fft_request(b):
            '''
            perform fft and plot 
            * two options: (partial) residues or full asymmetry
            * two modes: real amplitude or power
            * vectorized: range(len(self.fitargs)) is (0,1) or (0,n>1) for single or suite  

            WARNING: relies on self._the_model_._add_ or self._the_model_._fft_add_
            to produce the right function for each run (never checked yet)
            insert expected noise level (see bottom comment)

            '''
            import numpy as np
            from mujpy.aux.aux import derange, autops, ps, _ps_acme_score, _ps_peak_minima_score, plotile, get_title
            from copy import deepcopy
            import matplotlib.pyplot as P
            from matplotlib.path import Path
            import matplotlib.patches as patches
            import matplotlib.animation as animation

            ###################
            # PYPLOT ANIMATIONS
            ###################
            def animate(i):
                '''
                anim function
                update fft data, fit fft and their color 

                '''
                # color = next(ax_fft._get_lines.prop_cycler)['color']
                self.ax_fft.set_title(str(self._the_runs_[i][0].get_runNumber_int())+': '+get_title(self._the_runs_[i][0]))
                marks.set_ydata(ap[i])
                marks.set_color(color[i])
                line.set_ydata(apf[i])
                line.set_color(color[i])
                top = fft_e[i] 
                errs.set_facecolor(color[i])
                return line, marks, errs,


            def init():
                '''
                anim init function
                blitting (see wikipedia)
                to give a clean slate 

                '''
                self.ax_fft.set_title(str(self._the_runs_[0][0].get_runNumber_int())+': '+get_title(self._the_runs_[0][0]))
                marks.set_ydata(ap[0])
                marks.set_color(color[0])
                line.set_ydata(apf[0])
                line.set_color(color[0])
                top = fft_e[0] 
                errs.set_facecolor(color[0])
                return line, marks, errs, 

            def fft_std():
                '''
                Returns fft_e, array, one fft std per bin value per run index k  
                using time std ey[k] and filter filter_apo.
                The data slice is equivalent (not equal!) to 
                *       y[k] = yf[k] + ey[k]*np.random.randn(ey.shape[1])
                
                It is composed of l data plus l zero padding (n=2*l).

                Here we deal only with the first l data bins (no padding). 
                Assuming that the frequency noise is uniform, 
                the f=0 value of the filtered fft(y) is

                *      ap[k] = (y[k]*filter_apo).sum()

                and the j-th sample of the corresponding noise is

                *      eapj[k] = ey[k]*np.random.randn(ey.shape[1])*filter_apo).sum()

                Repeat n times to average the variance, 

                *  eapvar[k] = [(eapj[k]**2 for j in range(n)]

                *  fft_e = np.sqrt(eapvar.sum()/n)

                '''
                n = 10
                fft_e = np.empty(ey.shape[0])
                for k in range(ey.shape[0]):
                    eapvariance = [((ey[k]*np.random.randn(ey.shape[1])*filter_apo).sum())**2 for j in range(n)]
                    fft_e[k] = np.sqrt(sum(eapvariance)/n)
                return fft_e

            # ON_FFT_REQUEST STARTS HERE
            #################################
            # retrieve self._the_model_, pars, 
            #          fit_start,fit_stop=rangetup[0], rangetup[1]
            #          with rangetup = derange(self.fit_range.value), 

            if not self._the_model_._alpha_:
                # this is used as a check that an appropriate fit was performed
                # so we can assume that fit_range.value has already been checked
                self.console('No fit yet. Please first produce a fit attempt.')
                return
            if self._global_:
                print('not yet!')
            else:           
                ####################
                # setup fft
                ####################
                dt = self.time[0,1]-self.time[0,0]
                fmax = 0.5/dt  # max frequancy available
                rangetup = derange(self.fit_range.value,self.histoLength)
                # no checks, it has already been used in fit
                fit_start, fit_stop = int(rangetup[0]), int(rangetup[1]) # = self.time[fit_start]/dt, self.time[fit_stop]/dt
                # print('fit_start, fit_stop = {}, {}'.format(fit_start, fit_stop))
                l = fit_stop-fit_start # dimension of data
                df = 1/(dt*l)
                n = 2*l # not a power of 2, but surely even
                filter_apo = np.exp(-(dt*np.linspace(0,l-1,l)*float(fft_filter.value))**3) # hypergaussian filter mask
                                                 # is applied as if first good bin were t=0
                filter_apo = filter_apo/sum(filter_apo)/dt # approx normalization
                # try hypergauss n=3, varying exponent
                dfa = 1/n/dt         # digital frequency resolution

                #####################################################################################
                # asymm, asyme and the model are a row arrays if _single_ and matrices if not _single_
                #####################################################################################

                ##########################################
                # zero padding, apodization [and residues]
                ##########################################
                y = np.zeros((self.asymm.shape[0],n)) # for data zero padded to n
                ey = np.zeros((self.asyme.shape[0],l)) #  for errors, l bins, non zero padded
                yf = np.zeros((self.asymm.shape[0],n)) # for fit function zero padded to n
                for k in range(len(self.fitargs)):
                    pars = [self.fitargs[k][name] for name in self.minuit_parameter_names]
                    yf[k,0:l] = self._the_model_._add_(self.time[0,fit_start:fit_stop],*pars) # full fit zero padded, 
                if residues_or_asymmetry.value == 'Residues': 
                    fft_include_components = []
                    fft_include_da = False
                    for j,dic in enumerate(self.model_components):
                        if dic['name']=='da' and self.fftcheck[j].value:
                            fft_include_da = True # flag for "da is a component" and "include it"
                        elif dic['name']!='da': # fft_include_components, besides da, True=include, False=do not  
                            fft_include_components.append(self.fftcheck[j].value) # from the gui FFT checkboxes
                        self._the_model_._fft_init(fft_include_components,fft_include_da) # sets _the_model_ in fft 
                # t = deepcopy(self.time[fit_start:fit_stop])
                # print('self.time.shape = {}, t.shape = {}, range = {}'.format(self.time.shape,t.shape,fit_stop-fit_start))
                for k in range(len(self.fitargs)):
                    y[k,0:l] = self.asymm[k,fit_start:fit_stop] # zero padded data
                    ey[k] = self.asyme[k,fit_start:fit_stop] #  slice of time stds
                    
                    # print('yf.shape = {}, the_model.shape = {}'.format(yf[k,0:l].shape,t.shape))
                    ############################################
                    # if Residues
                    # subtract selected fit components from data
                    ############################################ 
                    if residues_or_asymmetry.value == 'Residues': 
                                     # fft partial subtraction mode: only selected components are subtracted
                        pars = [self.fitargs[k][name] for name in self.minuit_parameter_names]
                        y[k,0:l] -= self._the_model_._add_(self.time[0,fit_start:fit_stop],*pars)
                    y[k,0:l] *= filter_apo # zero padded, filtered data or residues
                    yf[k,0:l] *= filter_apo # zero padded, filtered full fit function
                #################################################
                # noise in the FFT: with scale=1 noise in n data bins, one gets sqrt(n/2) noise per fft bin, real and imag
                # generalising to scale=sigma noise in n bins -> sqrt(0.5*sum_i=1^n filter_i)
                #################################################
                fft_e = fft_std() # array of fft standard deviations per bin for each run

                fft_amplitude = np.fft.fft(y)  # amplitudes (complex), matrix with rows fft of each run
                fftf_amplitude = np.fft.fft(yf)  # amplitudes (complex), same for fit function
                #################
                # frequency array
                #################
                nf = np.hstack((np.linspace(0,l,l+1,dtype=int), np.linspace(-l+1,-1,l-2,dtype=int)))
                f = nf*dfa  # all frequencies, l+1 >=0 followed by l-1 <0
                rangetup = derange(fft_range.value,fmax,int_or_float='float') # translate freq range into bins
                fstart, fstop = float(rangetup[0]), float(rangetup[1]) 
                start, stop = int(round(fstart/dfa)), int(round(fstop/dfa))
                f = deepcopy(f[start:stop]) # selected slice
                
                ########################
                # build or recall Figure
                ########################
                if self.fig_fft: # has been set to a handle once
                    self.fig_fft.clf()
                    self.fig_fft,self.ax_fft = P.subplots(num=self.fig_fft.number)
                else: # handle does not exist, make one
                    self.fig_fft,self.ax_fft = P.subplots(figsize=(6,4))
                    self.fig_fft.canvas.set_window_title('FFT')
                self.ax_fft.set_xlabel('Frequency [MHz]')
                self.ax_fft.set_title(get_title(self._the_runs_[0][0]))
                xm, xM = f.min(),f.max()                    
                self.ax_fft.set_xlim(xm,xM)
                if real_or_power.value=='Real part':
                    ########################
                    # REAL PART
                    # APPLY PHASE CORRECTION
                    # try acme
                    ########################
                    with self._outputtab_:
                        fftf_amplitude[0], p0, p1 = autops(fftf_amplitude[0],'acme') # fix phase on theory 
                    fft_amplitude[0] = ps(fft_amplitude[0], p0=p0 , p1=p1).real # apply it to data
                    for k in range(1,fft_amplitude.shape[0]):
                        fft_amplitude[k] = ps(fft_amplitude[k], p0=p0 , p1=p1)
                        fftf_amplitude[k] = ps(fftf_amplitude[k], p0=p0 , p1=p1)
                    ap = deepcopy(fft_amplitude[:,start:stop].real)
                    apf = deepcopy(fftf_amplitude[:,start:stop].real)
                    label = 'Real part'
                else:
                    ##################
                    # POWER
                    ##################
                    ap = fft_amplitude.real[:,start:stop]**2+fft_amplitude.imag[:,start:stop]**2
                    apf = fftf_amplitude.real[:,start:stop]**2+fftf_amplitude.imag[:,start:stop]**2
                    label = 'Power'

                ########
                # tile 
                ########
                if not anim_check.value or self._single_:  # TILES: creates matrices for offset multiple plots      
                    foffset = 0 # frequency offset
                    yoffset = 0.1*apf.max()  # add offset to each row, a fraction of the function maximum
                    f, ap, apf = plotile(f,xdim=ap.shape[0],offset=foffset),\
                                 plotile(ap,offset=yoffset),\
                                 plotile(apf,offset=yoffset) 
                    # f, ap, apf are (nrun,nbins) arrays

                #############
                # animation
                #############
                if anim_check.value and not self._single_: # a single cannot be animated
                    ##############
                    # initial plot
                    ##############
                    color = []
                    for k in range(ap.shape[0]):
                        color.append(next(self.ax_fft._get_lines.prop_cycler)['color'])
                    yM = 1.02*max(ap.max(),apf.max())
                    ym = min(0,1.02*ap.min(),1.02*apf.min())
                    line, = self.ax_fft.plot(f,apf[0],'-',lw=1,color=color[0],alpha=0.8)
                    marks, = self.ax_fft.plot(f,ap[0],'o',ms=2,color=color[0],alpha=0.8)
                    self.ax_fft.set_ylim(ym,yM)
                    left, bottom, right, top = f[0],0.,f[-1],fft_e[0]
                    verts = [
                            (left, bottom), # left, bottom
                            (left, top), # left, top
                            (right, top), # right, top
                            (right, bottom), # right, bottom
                            (0., 0.), # ignored
                            ]

                    codes = [Path.MOVETO,
                             Path.LINETO,
                             Path.LINETO,
                             Path.LINETO,
                             Path.CLOSEPOLY,
                             ]
                    path = Path(verts, codes)
                    errs = patches.PathPatch(path, facecolor=color[0], lw=0, alpha=0.3)
                    self.ax_fft.add_patch(errs)
                    #######
                    # anim
                    #######
                    self.anim_fft = animation.FuncAnimation(self.fig_fft, animate, 
                                                            np.arange(0,len(self.fitargs)),
                                                            init_func=init,
                                                            interval=anim_delay.value,
                                                            blit=False)

                ###############################
                # single and tiles with offset
                ###############################
                else:
                    # print('f.shape = {}, ap.shape = {}'.format(f.shape,ap.shape))  
                    color = []                  
                    for k in range(ap.shape[0]): 
                        color.append(next(self.ax_fft._get_lines.prop_cycler)['color'])
                        self.ax_fft.plot(f[k],ap[k],'o',ms=2,alpha=0.5,color=color[k]) # f, ap, apf are plotiled!
                        self.ax_fft.plot(f[k],apf[k],'-',lw=1,alpha=0.5,color=color[k])
                        self.ax_fft.fill_between([f[0,0],f[0,-1]],[k*yoffset,k*yoffset],[k*yoffset+fft_e[k],k*yoffset+fft_e[k]],facecolor=color[k],alpha=0.2)
                    ###################
                    # errors, alpha_version for single
                    ################### 
#                    if self._single_: 
                    self.ax_fft.relim(), self.ax_fft.autoscale_view()
                    ym,yM = self.ax_fft.get_ylim()
                    xm,xM = self.ax_fft.get_xlim()
                    ytext = yM-(ap.shape[0]+1)*yoffset
                    xtext = xM*0.90
                    for k in range(ap.shape[0]):   
                        ytext = ytext+yoffset                  
                        self.ax_fft.text(xtext,ytext,str(self._the_runs_[k][0].get_runNumber_int()),color=color[k])              

                if residues_or_asymmetry.value == 'Residues':
                    self.ax_fft.set_ylabel('FFT '+label+' (Residues/Fit)')
                    self._the_model_._include_all_() # usual _the_model_ mode: all components included
                else:
                    self.ax_fft.set_ylabel('FFT '+label+' (Asymmetry/Fit)')

            self.fig_fft.canvas.manager.window.tkraise()
            P.draw()

        def on_filter_changed(change):
            '''
            observe response of fit tab widgets:
            validate float   
     
            '''
            string = change['owner'].value # description is three chars ('val','fun','flg') followed by an integer nint
                                               # iterable in range(ntot), total number of internal parameters
            try: 
                float(string)
            except:
                change['owner'].value = '{:.4f}'.format(filter0)

        def on_range(change):
            '''
            observe response of FFT range widgets:
            check for validity of function syntax

            '''
            from mujpy.aux.aux import derange

            fmax = 0.5/(self.time[0,1]-self.time[0,0])
            returnedtup = derange(fft_range.value,fmax,int_or_float='float') # errors return (-1,-1),(-1,0),(0,-1), good values are all positive
            if sum(returnedtup)<0:
                fft_range.background_color = "mistyrose"
                fft_range.value = fft_range0
            else:
                fft_range.background_color = "white"


        def on_start_stop(change):
            if anim_check.value: 
                if change['new']:
                    self.anim_fft.event_source.start()
                else:
                    self.anim_fft.event_source.stop()

        # begins fft gui
        import numpy as np
        from ipywidgets import HBox, VBox, Layout, Button, FloatText, Text, IntText, Dropdown, Checkbox, ToggleButton

        # must inherit/retrieve self._the_model_, pars, fit_range = range(fit_start,fit_stop)
        # layout a gui to further retrieve  
        # fft_range (MHz), lb (us-1), real_amplitude (True/False) if False then power, autophase (True/False)

        # Layout gui
        fft_button = Button(description='Do FFT',layout=Layout(width='12%'))
        fft_button.style.button_color = self.button_color
        fft_button.on_click(on_fft_request)

        filter0 = 0.3
        fft_filter = Text(description='Filter ($\mu s^{-1}$)',
                               value='{:.4f}'.format(filter0),
                               layout=Layout(width='20%'),
                               continuous_update=False) # self.filter.value
        fft_filter.observe(on_filter_changed,'value')

        fft_range0 = '0,50'
        fft_range = Text(description='fit range\nstart,stop\n  (MHz)',
                         value=fft_range0,
                         layout=Layout(width='28%'),
                         continuous_update=False)
        fft_range.style.description_width='60%'
        fft_range.observe(on_range,'value')

        real_or_power = Dropdown(options=['Real part','Power'], 
                                 value='Real part',
                                 layout=Layout(width='12%'))

        residues_or_asymmetry = Dropdown(options=['Residues','Asymmetry'], 
                                         value='Residues',
                                         layout=Layout(width='13%'))

        autophase = Checkbox(description='Autophase',
                             value=True,
                             layout=Layout(width='15%'))
        autophase.style.description_width='10%'

        anim_check = Checkbox(description='Animate',value=False, layout=Layout(width='12%'))
        anim_check.style.description_width = '1%'

        anim_delay = IntText(description='Delay (ms)',value=1000, layout=Layout(width='20%'))

        anim_stop_start = ToggleButton(description='start/stop',value=True)
        anim_stop_start.observe(on_start_stop,'value')
        anim_stop_start.style.button_color = self.button_color
        fft_frame_handle = VBox(description='FFT_bar',children=[HBox(description='first_row',children=[fft_button,
                                                                fft_filter,
                                                                fft_range,
                                                                real_or_power,
                                                                residues_or_asymmetry,
                                                                autophase]),
                                                                HBox(description='second_row',children=[anim_check,
                                                                anim_delay, anim_stop_start])])
        # now collect the handles of the three horizontal frames to the main fit window 
        self.mainwindow.children[4].children = [fft_frame_handle] 
                             # add the list of widget handles as the third tab, fit

##########################
# FIT
##########################
    def fit(self, model_in = 'daml'): # self.fit(model_in = 'mgmgbl') produces a different layout
        '''
        fit tab of mugui, used 

        - to set: self.alpha.value, self.offset.value, forw and backw groups fit and plot ranges, model version           
        - to display: model name 
        - to activate: fit, plot and update buttons
        - to select and load model (load from folder is missing/deprecated in python)
        - to select parameters value, fix, function, fft subtract

        ::

        # the calculation is performed in independent class mucomponents
        # the methods are "inherited" by mugui 
        # via the reference instance self._the_model_, initialized in steps: 
        #     __init__ share initial attributes (constants) 
        #     _available_components_ automagical list of mucomponents 
        #     clear_asymmetry: includes reset check when suite is implemented
        #     create_model: lay out self._the_model_
        #     delete_model: for a clean start
        #     functions use eval, evil but needed, checked by muvalid, 
        #     iminuit requires them to be formatted as fitarg by int2min
        #     help  
        #     load 
        #     save_fit/load_ft save results in mujpy format (dill)
        #     write_csv produces a qtiplot/origin loadable summary
        # 
        # Three fit types: single, suite no global, suite global.
        # Suite non global iterates a single fit over several runs
        # Suite global performs a single fit over many runs,
        # with common (global) and run dependent (local) parameters

        '''


        from mujpy.mucomponents.mucomponents import mumodel
        import numpy as np

        def _available_components_():
            from iminuit import describe
            '''
            Method, returns a template tuple of dictionaries (one per fit component):

            Each dictionary contains 'name' and 'pars', 
            the latter in turns is a list of dictionaries, one per parameter, 'name','error,'limits'

            ::  ({'name':'bl','pars':[{'name':'asymmetry','error':0.01,'limits'[0,0]},

                                      {'name':'Lor_rate','error':0.01,'limits'[0,0]}}, 
                                      ...)

            retreived magically from the mucomponents class.

            '''
            _available_components = [] # is a list, mutable
            # generates a template of available components.
            for name in [module for module in dir(mumodel()) if module[0]!='_']: # magical extraction of component names
                pars = describe(mumodel.__dict__[name])[2:]            # the [2:] is because the first two arguments are self and x
                _pars = []
                # print('pars are {}'.format(pars))
                for parname in pars:
                # The style is like iminuit fitargs, but not exactly,
                # since the latter is a minuit instance:
                # it will contain parameter name: parname+str(k)[+'_'+str(nrun)]
                # error_parname, fix_parname (False/True), limits_parname, e.g.
                #   {'amplitude1_354':0.154,'error_amplitude1_354':0.01,'fix_amplitude1_354':False,'limits_amplitude1_354':[0, 0]
                # 
                # In this template only
                #   {'name':'amplitude','error':0.01,'limits':[0, 0]}
                    error, limits = 0.01, [0, 0] # defaults
                    if parname == 'field' or parname == 'phase' or parname == 'dipfield': error = 1.0
                    if parname == 'beta': error,limits = 0.05, [1.e-2, 1.e2]
                    # add here special cases for errors and limits, e.g. positive defined parameters
                    _pars.append({'name':parname,'error':error,'limits':limits})
                _available_components.append({'name':name,'pars':_pars})
            self.available_components = (_available_components) # these are the mucomponents method directories
                                                                # transformed in tuple, immutable
            self.component_names = [self.available_components[i]['name'] 
                                    for i in range(len(self.available_components))] # list of just mucomponents method names

        def addcomponent(name,label):
            '''
            myfit = MuFit()
            addcomponent('ml') # adds e.g. a mu precessing, lorentzian decay, component
            this method adds a component selected from self.available_component, tuple of directories
            with zeroed values, stepbounds from available_components, flags set to '~' and empty functions

            '''
            from copy import deepcopy
            if name in self.component_names:
                k = self.component_names.index(name)
                npar = len(self.available_components[k]['pars']) # number of pars
                pars = deepcopy(self.available_components[k]['pars']) # list of dicts for 
                # parameters, {'name':'asymmetry','error',0.01,'limits',[0, 0]}

                # now remove parameter name degeneracy                   
                for j, par in enumerate(pars):
                    pars[j]['name'] = par['name']+label
                    pars[j].update({'value':0.0})
                    pars[j].update({'flag':'~'})
                    pars[j].update({'function':''}) # adds these three keys to each pars dict
                    # they serve to collect values in mugui
                self.model_components.append({'name':name,'pars':pars})
                return True # OK code
            else:
                self.console('\nWarning: '+name+' is not a known component. Not added.\n'+
                           'With myfit = mufit(), type myfit.help to see the available components')
                return False # error code

        def create_model(name):
            '''
            myfit = MuFit()       
            myfit.create_model('daml') # adds e.g. the two component 'da' 'ml' model
            this method adds a model of components selected from the available_component tuple of  
            directories

            with zeroed values, stepbounds from available_components, flags set to '~' and empty functions

            '''
            import string
    # name 2_0_mlml_blbl for 2 global parameters (A0 R), 0 kocal parameters (B end T) and two models 
    # e.g. alpha fit with a WTF and a ZF run, with two muon fractions of amplitude A0*R and A0*(1-R) respectively
    # find the three underscores in name by
    # [i for i in range(len(name)) if name.startswith('_', i)]

            components = checkvalidmodel(name)
            if components: # exploits the fact that [] is False and ['da'] is true
                self.model = name
                self.model_components = [] # start from empty model
                for k,component in enumerate(components):
                    label = string.ascii_uppercase[k]
                    if not addcomponent(component,label):
                        return False
                return True
            else:
                return False

        def checkvalidmodel(name):
            '''
            checkvalidmodel(name) checks that name is a 
            ::      2*component string of valid component names, e.g.
                                  'daml' or 'mgmgbl'
            
            '''
            components = [name[i:i+2] for i in range(0, len(name), 2)]
            for component in components:            
                if component not in self.component_names:
                    self.console('Warning: '+component+' is not a known component. Not added.\n'+
                               'With myfit = mufit(), type myfit.help to see the available components')
                    return [] # error code
            return components

        def chi(t,y,ey,pars):
            '''
            stats for the right side of the plot 

            '''
            nu = len(t) - self.freepars # degrees of freedom in plot
            # self.freepars is calculated in int2min
            self._the_model_._load_data_(t,y,int2_int(),self.alpha.value,e=ey)
            f = self._the_model_._add_(t,*pars) # f for histogram
            chi2 = self._the_model_._chisquare_(*pars)/nu # chi2 in plot
            return nu,f,chi2
            
        def fitplot(guess=False,plot=False):
            '''
            Plots fit results in external Fit window
            guess=True plot dash guess values
            guess=False plot best fit results
            plot=False best fit, invoke write_csv
            plot=True do not
            
            This is a complex routine that allows for
            ::  - single, multiple or global fits
                - fit range different form plot range
                - either
                    one plot range, the figure is a subplots((2,2))
                        plot ax_fit[(0,0), chi2_prints ax_fit[(0,-1)]
                        residues ax_fit[(1,0)], chi2_histograms ax_fit[(1,-1)]
                    two plot ranges, early and late, the figure is a subplots((3,2))
                        plot_early ax_fit[(0,0)], plot_late ax_fit[(0,1)], chi2_prints ax_fit[(0,-1)]
                        residues_early ax_fit[(1,0)], residues_late ax_fit[(1,1)], chi2_histograms ax_fit[(1,-1)]

            If multi/globalfit, it also allows for either
            ::  - anim display 
                - offset display 

            '''  
            import matplotlib.pyplot as P
            from mujpy.aux.aux import derange, rebin, get_title, plotile, set_bar
            from scipy.stats import norm
            from scipy.special import gammainc
            import matplotlib.path as path
            import matplotlib.patches as patches
            import matplotlib.animation as animation

            font = {'family':'Ubuntu','size':10}
            P.rc('font', **font)

            ###################
            # PYPLOT ANIMATIONS
            ###################
            def animate(i):
                '''
                anim function
                update errorbar data, fit, residues and their color,
                       chisquares, their histograms 

                '''
                # from mujpy.aux.aux import get_title
                # print('animate')
                # nufit,ffit,chi2fit = chi(tfit[0],yfit[i],eyfit[i],pars[i])
                # nu,dum,chi2plot = chi(t[0],y[i],ey[i],pars[i])
                # color = next(self.ax_fit[(0,0)]._get_lines.prop_cycler)['color']
                line.set_ydata(y[i]) # begin errorbar
                line.set_color(color[i])
                line.set_markerfacecolor(color[i])
                line.set_markeredgecolor(color[i])
                segs = [np.array([[q,w-a],[q,w+a]]) for q,w,a in zip(t[0],y[i],ey[i])]
                ye[0].set_segments(segs) 
                ye[0].set_color(color[i]) # end errorbar
                fline.set_ydata(f[i]) # fit
                fline.set_color(color[i])
                res.set_ydata(y[i]-fres[i]) # residues
                res.set_color(color[i])
                # self.ax_fit[(0,0)].relim(), self.ax_fit[(0,0)].autoscale_view()
       
                if len(returntup)==5:
                    linel.set_ydata(ylate[i]) # begin errorbar
                    linel.set_color(color[i])
                    linel.set_markerfacecolor(color[i])
                    linel.set_markeredgecolor(color[i])
                    segs = [np.array([[q,w-a],[q,w+a]]) for q,w,a in zip(tlate[0],ylate[i],eylate[i])]
                    yel[0].set_segments(segs) 
                    yel[0].set_color(color[i]) # end errorbar
                    flinel.set_ydata(fl[i]) # fit
                    flinel.set_color(color[i])
                    resl.set_ydata(ylate[i]-flres[i]) # residues
                    resl.set_color(color[i])
                    # self.ax[(0,1)].relim(), self.ax[(0,1)].autoscale_view()

                self.ax_fit[(0,0)].set_title(get_title(self._the_runs_[i][0]))
                nhist,dum = np.histogram((yfit[i]-ffit[i])/eyfit[i],xbin)
                top = bottomf + nhist
                vertf[1::5, 1] = top
                vertf[2::5, 1] = top

                nhist,dum = np.histogram((y[i]-fres[i])/ey[i],xbin,weights=nufit[i]/nu[i]*np.ones(t.shape[1]))
                top = bottomp + nhist
                vertp[1::5, 1] = top
                vertp[2::5, 1] = top
                patchplot.set_facecolor(color[i])
                patchplot.set_edgecolor(color[i])
                nufitplot.set_ydata(nufit[i]*yh)

                string = '$\chi^2_f=$ {:.4f}\n ({:.2f}-{:.2f})\n$\chi^2_c=$ {:.4f}\n{} dof\n'.format(chi2fit[i],
                                                                            lc[i],hc[i],gammainc(chi2fit[i],nufit[i]),nufit[i])
                if  len(returntup)==5: 
                    nulate,dum,chi2late = chi(tlate[0],ylate[i],eylate[i],pars[i])
                    string += '$\chi^2_e=$ {:.4f}\n$\chi^2_l=$ {:.4f}'.format(chi2plot[i],chi2late)
                else:
                    string += '$\chi^2_p=$ {:.4f}'.format(chi2plot[i])
                text.set_text('{}'.format(string))

                if  len(returntup)==5: 
                    return line, ye[0], fline, res, linel, yel[0], flinel, resl, patchfit, patchplot, nufitplot, text 
                else:
                    return line, ye[0], fline, res, patchfit, patchplot, nufitplot, text

            def init():
                '''
                anim init function
                blitting (see wikipedia)
                to give a clean slate 

                '''
                from mujpy.aux.aux import get_title
                # nufit,ffit,chi2fit = chi(tfit[0],yfit[0],eyfit[0],pars[0])
                # nu,dum,chi2plot = chi(t[0],y[0],ey[0],pars[0])
                # color = next(self.ax_fit[(0,0)]._get_lines.prop_cycler)['color']
                line.set_ydata(y[0]) # begin errorbar
                line.set_color(color[0])
                segs = [np.array([[q,w-a],[q,w+a]]) for q,w,a in zip(t[0],y[0],ey[0])]
                ye[0].set_segments(segs)
                ye[0].set_color(color[0]) # end errorbar
                fline.set_ydata(f[0]) # fit
                fline.set_color(color[0])
                res.set_ydata(y[0]-fres[0]) # residues
                res.set_color(color[0])

                if len(returntup)==5:
                    linel.set_ydata(ylate[0]) # begin errorbar
                    linel.set_color(color[0])
                    segs = [np.array([[q,w-a],[q,w+a]]) for q,w,a in zip(tlate[0],ylate[0],eylate[0])]
                    yel[0].set_segments(segs)
                    yel[0].set_color(color[0]) # end errorbar
                    flinel.set_ydata(fl[0]) # fit
                    flinel.set_color(color[0])
                    resl.set_ydata(ylate[0]-flres[0]) # residues
                    resl.set_color(color[0])

                self.ax_fit[(0,0)].set_title(get_title(self._the_runs_[0][0]))
                nhist,dum = np.histogram((yfit[0]-ffit[0])/eyfit[0],xbin)
                top = bottomf + nhist
                vertf[1::5, 1] = top
                vertf[2::5, 1] = top

                nhist,dum = np.histogram((y[0]-fres[0])/ey[0],xbin,weights=nufit[0]/nu[0]*np.ones(t.shape[1]))
                top = bottomp + nhist
                vertp[1::5, 1] = top
                vertp[2::5, 1] = top
                patchplot.set_facecolor(color[0])
                patchplot.set_edgecolor(color[0])
                nufitplot.set_ydata(nufit[0]*yh)
                string = '$\chi^2_f=$ {:.4f}\n ({:.2f}-{:.2f})\n$\chi^2_c=$ {:.4f}\n{} dof\n'.format(chi2fit[0],
                                                                            lc[0],hc[0],gammainc(chi2fit[0],nufit[0]),nufit[0])
                if  len(returntup)==5: 
                    nulate,dum,chi2late = chi(tlate[0],ylate[0],eylate[0],pars[0])
                    string += '$\chi^2_e=$ {:.4f}\n$\chi^2_l=$ {:.4f}'.format(chi2plot[0],chi2late)
                else:
                    string += '$\chi^2_p=$ {:.4f}'.format(chi2plot[0])
                text.set_text('{}'.format(string))
                # print('init')
                if  len(returntup)==5: 
                    return line, ye[0], fline, res,  linel, yel[0], flinel, resl, patchfit, patchplot, nufitplot, text 
                else:
                    return line, ye[0], fline, res, patchfit, patchplot, nufitplot, text


            #     FITPLOT BEGINS HERE
            ######################################################
            # pars is a list of lists of best fit parameter values
            # self.time is a 1D array
            # self.asymm, self.asyme are 2D arrays 
            # y, ey, f, fres, ylate, eylate, fl, flres, yfit, eyfit, ffit are 2D arrays
            # tf, tfl, tlate, tfit are 1D array 

            ##############################
            # plot according to plot_range
            ##############################
            
            # must provide maximum available bins = self.histoLength 
            returntup = derange(self.plot_range.value,self.histoLength) 
            if sum(n<0 for n in returntup)>0:  # illegal values self.plot_range.value
                                               # signalled by returntup = (-1,-1);
                                               # restore defaults
                tmp = self.plot_range.value
                self.plot_range.value = self.plot_range0
                self.plot_range.background_color = "mistyrose"
                self.console('Wrong plot range: {}'.format(tmp))
                return
            self.asymmetry() # prepare asymmetry, 
            ############################################
            #  choose pars for first/single fit function
            ############################################ 
            fitarg  = int2min(return_names=True) # from dash, fitarg is a list of dictionaries
            # print('fitarg = {}\nself.minuit_parameter_names = {}'.format(fitarg,self.minuit_parameter_names))
            if guess: # from dash, for plot guess
                pars = [[fitarg[k][name] for name in self.minuit_parameter_names] for k in range(len(fitarg))]
                ###############################################################
                # mock data loading to set alpha and global in self._the_model_
                # in case no fit was done yet
                ###############################################################
                if not self._the_model_._alpha_:  # False if no _load_data_ yet 
                    if self._global_:
                        self.console('global, mumodel load_data') 
                        self._the_model_._load_data_(self.time[0],self.asymm,int2_int(),
                                                     self.alpha.value,e=self.asyme,
                                                     _nglobal_=self.nglobals,
                                                     _locals_=self.locals)                
                    else:
                        # self.console('no global, mumodel load_data') # debug
                        self._the_model_._load_data_(self.time[0],self.asymm[0],int2_int(),self.alpha.value,e=self.asyme[0]) # pok = is debug
                        # string = 'Load data iok = {} [0=fine]. Time last = {:.2e}, A[0] = {:.2f}, alpha = {:.3f}'.format(pok,self._the_model_._x_[-1],self._themodel_._y_[0],self._the_model_._alpha_) # debug
                        # self.console(string)
            else:  # from lastfit, for best fit and plot best fit
                pars = [[self.fitargs[k][name] for name in self.minuit_parameter_names] for k in range(len(self.fitargs))]
            ##########################################
            # now self.time is a 1D array
            # self.asymm, self.asyme are 1D or 2D arrays
            # containing asymmetry and its std, 
            # for either single run or suite of runs
            # pars[k] is the k-th par list for the fit curve of the k-th data row
            ##########################################

            ###############################################
            # rebinnig for plot (different packing from fit)
            ###############################################
            # early and late plots
            ######################
            if len(returntup)==5: # start stop pack=packearly last packlate
                start, stop, pack, last, packlate = returntup
                tlate,ylate,eylate = rebin(self.time,self.asymm,[stop,last],packlate,e=self.asyme)
                tfl,dum = rebin(self.time,self.asymm,[stop,last],1)
                ncols, width_ratios = 3,[2,2,1]

            ###################
            # single range plot
            ###################
            else:
                pack = 1
                ncols, width_ratios = 2,[4,1]
                if len(returntup)==3: # plot start stop pack
                    start, stop, pack = returntup
                elif len(returntup)==2: # plot start stop
                    start, stop = returntup

            t,y,ey = rebin(self.time,self.asymm,[start,stop],pack,e=self.asyme)
            tf,dum = rebin(self.time,self.asymm,[start,stop],1)
            yzero = y[0]-y[0]

            #############################
            # rebinning of data as in fit
            #############################
            fittup = derange(self.fit_range.value,self.histoLength) # range as tuple
            fit_pack =1
            if len(fittup)==3: # plot start stop pack
                fit_start, fit_stop, fit_pack = fittup[0], fittup[1], fittup[2]
            elif len(fittup)==2: # plot start stop
                fit_start, fit_stop = fittup[0], fittup[1]
            # if not self._single_ each run is a row in 2d ndarrays yfit, eyfit

            tfit,yfit,eyfit = rebin(self.time,self.asymm,[fit_start,fit_stop],fit_pack,e=self.asyme)
            # print('pars = {}'.format(pars))
            # print('t = {}'.format(t))

            f = np.array([self._the_model_._add_(tf[0],*pars[k]) for k in range(len(pars))]) # tf,f for plot curve
            fres = np.array([self._the_model_._add_(t[0],*pars[k]) for k in range(len(pars))]) # t,fres for residues
            ffit = np.array([self._the_model_._add_(tfit[0],*pars[k]) for k in range(len(pars))]) # t,fres for residues

            if len(returntup)==5:
            ##############################################
            # prepare fit curves for second window, if any
            ##############################################         
                fl = np.array([self._the_model_._add_(tfl[0],*pars[k]) for k in range(len(pars))]) # tfl,fl for plot curve 
                flres = np.array([self._the_model_._add_(tlate[0],*pars[k]) for k in range(len(pars))]) # tlate,flate for residues 
            ###############################
            #  set or recover figure, axes 
            ###############################
            if self.fig_fit: # has been set to a handle once
                self.fig_fit.clf()
                self.fig_fit,self.ax_fit = P.subplots(2,ncols,sharex = 'col', 
                             gridspec_kw = {'height_ratios':[3,1],'width_ratios':width_ratios},
                                            num=self.fig_fit.number)
                self.fig_fit.subplots_adjust(hspace=0.05,top=0.90,bottom=0.12,right=0.97,wspace=0.03)
            else: # handle does not exist, make one
                self.fig_fit,self.ax_fit = P.subplots(2,ncols,figsize=(6,4),sharex = 'col',
                             gridspec_kw = {'height_ratios':[3, 1],'width_ratios':width_ratios})
                self.fig_fit.canvas.set_window_title('Fit')
                self.fig_fit.subplots_adjust(hspace=0.05,top=0.90,bottom=0.12,right=0.97,wspace=0.03)

            ##########################
            #  plot data and fit curve
            ##########################
            #############
            # animation
            #############
            if anim_check.value and not self._single_: # a single cannot be animated
                # THIS BLOCK TAKE CARE OF THE FIRST ROW OF DATA (errobars, fit curve, histograms and all)
                # pars[k] are the parameters to the run of the FIRST row, both for global and multi fits 
                #      in anim therefore FIT CURVES  (f, fres, fl, flres) ARE ALWAYS 1D ARRAYS
                # animate must take care of updating parameters and producing correct fit curves
                ##############
                # initial plot
                ##############
                nufit,dum,chi2fit = chi(tfit[0],yfit[0],eyfit[0],pars[0])
                color = []
                for k in range(len(self.fitargs)):
                    color.append(next(self.	ax_fit[(0,0)]._get_lines.prop_cycler)['color'])
                line, xe, ye, = self.ax_fit[(0,0)].errorbar(t[0],y[0],yerr=ey[0],
                                            fmt='o',elinewidth=1.0,ms=2.0,
                                            mec=color[0],mfc=color[0],ecolor=color[0],alpha=0.5) # data
                fline, = self.ax_fit[(0,0)].plot(tf[0],f[0],'-',lw=1.0,color=color[0],alpha=0.5) # fit
                res, = self.ax_fit[(1,0)].plot(t[0],y[0]-fres[0],'-',lw=1.0,color=color[0],alpha=0.5) # residues
                self.ax_fit[(1,0)].plot(t[0],yzero,'k-',lw=0.5,alpha=0.3) # zero line
                ym,yM =  y.min()*1.02,y.max()*1.02
                rm,rM =  (y-fres).min()*1.02,(y-fres).max()*1.02
                ym,rm = min(ym,0), min(rm,0)

                ############################
                # plot second window, if any
                ############################
                if len(returntup)==5:
                    linel, xel, yel, = self.ax_fit[(0,1)].errorbar(tlate[0],ylate[0],yerr=eylate[0],
                                                fmt='o',elinewidth=1.0,ms=2.0,alpha=0.5,
                                                mec=color[0],mfc=color[0],ecolor=color[0]) # data
                    flinel, = self.ax_fit[(0,1)].plot(tfl[0],fl[0],'-',lw=1.0,alpha=0.5,color=color[0]) # fit
                    self.ax_fit[(0,1)].set_xlim(tlate[0,0], tlate[0,-1])
                    # plot residues
                    resl, = self.ax_fit[(1,1)].plot(tlate[0],ylate[0]-flres[0],'-',lw=1.0,alpha=0.5,color=color[0]) # residues
                    self.ax_fit[(1,1)].plot(tlate[0],ylate[0]-ylate[0],'k-',lw=0.5,alpha=0.3) # zero line
                    self.ax_fit[(0,1)].set_xlim(tlate.min(),tlate.max()) # these are the global minima
                    self.ax_fit	[(1,1)].set_xlim(tlate.min(),tlate.max())
                    self.ax_fit	[(1,1)].set_xlim(tlate.min(),tlate.max())
                    self.ax_fit	[(1,1)].set_xlim(tlate.min(),tlate.max())
                    yml,yMl =  ylate.min()*1.02,ylate.max()*1.02
                    rml,rMl =  (ylate-flres).min()*1.02,(ylate-flres).max()*1.02
                    ym,yM,rm,rM = min(ym,yml),max(yM,yMl),min(rm,rml),max(rM,rMl)
                    self.ax_fit[(0,1)].set_ylim(ym,yM)
                    self.ax_fit[(1,1)].set_ylim(rm,rM)
                    self.ax_fit[(0,1)].set_yticklabels([])
                    self.ax_fit[(1,1)].set_yticklabels([])
                ###############################
                # set title, labels
                ###############################
                # print('title = {}'.format(get_title(self._the_runs_[0][0])))
                self.ax_fit[(0,0)].set_title(get_title(self._the_runs_[0][0]))
                self.ax_fit[(0,0)].set_xlim(0,t.max())
                self.ax_fit[(0,0)].set_ylim(ym,yM)
                self.ax_fit[(1,0)].set_ylim(rm,rM)
                self.ax_fit[(1,0)].set_xlim(0,t.max())
                self.ax_fit[(0,0)].set_ylabel('Asymmetry')
                self.ax_fit[(1,0)].set_ylabel('Residues')
                self.ax_fit[(1,0)].set_xlabel(r'Time [$\mu$s]')
                self.ax_fit[(1,-1)].set_xlabel("$\sigma$")
                self.ax_fit[(1,-1)].set_yticklabels(['']*len(self.ax_fit[(1,-1)].get_yticks()))    
                self.ax_fit[(1,-1)].set_xlim([-5., 5.])
                self.ax_fit[(0,-1)].axis('off')

                ########################
                # chi2 distribution: fit
                ########################
                xbin = np.linspace(-5.5,5.5,12)
                nhist,dum = np.histogram((yfit[0]-ffit[0])/eyfit[0],xbin) # fc, lw, alpha set in patches
                vertf, codef, bottomf, xlimf = set_bar(nhist,xbin) 

                barpathf = path.Path(vertf, codef)
                patchfit = patches.PathPatch(
                    barpathf, facecolor='w', edgecolor= 'k', alpha=0.5,lw=0.7)
                self.ax_fit[(1,-1)].add_patch(patchfit)  #hist((yfit-ffit)/eyfit,xbin,rwidth=0.9,fc='w',ec='k',lw=0.7)

                self.ax_fit[(1,-1)].set_xlim(xlimf[0],xlimf[1])
                # self.ax_fit[(1,-1)].set_ylim(0, 1.15*nhist.max())

                #########################################
                # chi2 distribution: plots, scaled to fit
                #########################################
                nu,dum,chi2plot = chi(t[0],y[0],ey[0],pars[0])
                nhist,dum = np.histogram((y[0]-fres[0])/ey[0],xbin,weights=nufit/nu*np.ones(t.shape[1]))
                vertp, codep, bottomp, xlimp = set_bar(nhist,xbin) # fc, lw, alpha set in patches

                barpathp = path.Path(vertp, codep)
                patchplot = patches.PathPatch(
                    barpathp, facecolor=color[0], edgecolor= color[0], alpha=0.5,lw=0.7)
                self.ax_fit[(1,-1)].add_patch(patchplot)  # hist((y[0]-f/ey[0],xbin,weights=nufit/nu*np.ones(t.shape[0]),rwidth=0.9,fc=color,alpha=0.2)

                ###############################
                # chi2 dist theo curve & labels 
                ###############################
                xh = np.linspace(-5.5,5.5,23)        # static
                yh = norm.cdf(xh+1)-norm.cdf(xh)     # static

                nufitplot, = self.ax_fit[(1,-1)].plot(xh+0.5,nufit*yh,'r-') # nufit depends on k
                mm = round(nufit/4)              # nu, mm, hb, cc, lc, hc depend on k
                hb = np.linspace(-mm,mm,2*mm+1)
                cc = gammainc((hb+nufit)/2,nufit/2) # muchi2cdf(x,nu) = gammainc(x/2, nu/2);
                lc = 1+hb[min(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nufit
                hc = 1+hb[max(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nufit
                string = 'F-B: {} - {}\n'.format(self.group[0].value,self.group[1].value)
                string += r'$\alpha=$ {}'.format(self.alpha.value)
                string += '\n$\chi^2_f=$ {:.4f}\n ({:.2f}-{:.2f})\n$\chi^2_c=$ {:.4f}\n{} dof\n'.format(chi2fit,
                                                                            lc,hc,gammainc(chi2fit,nufit),nufit)
                if  len(returntup)==5: 
                    nulate,dum,chi2late = chi(tlate[0],ylate[0],eylate[0],pars[0])
                    string += '$\chi^2_e=$ {:.4f}\n$\chi^2_l=$ {:.4f}'.format(chi2plot,chi2late)
                else:
                    string += '$\chi^2_p=$ {:.4f}'.format(chi2plot)
                text = self.ax_fit[(0,-1)].text(-4,0.3,string)

                self.fig_fit.canvas.manager.window.tkraise()
                # save all chi2 values now
                nufit,chi2fit,nu,chi2plot,lc,hc = [nufit],[chi2fit],[nu],[chi2plot],[lc],[hc] # initialize lists with k=0 value
                for k in range(1,len(self.fitargs)):
                    nufitk,dum,chi2fitk = chi(tfit[0],yfit[k],eyfit[k],pars[k])
                    nufit.append(nufitk)
                    chi2fit.append(chi2fitk)
                    nuk,dum,chi2plotk = chi(t[0],y[k],ey[k],pars[k])
                    nu.append(nuk)
                    chi2plot.append(chi2plotk)
                    mm = round(nufitk/4)              # nu, mm, hb, cc, lc, hc depend on k
                    hb = np.linspace(-mm,mm,2*mm+1)
                    cc = gammainc((hb+nufitk)/2,nufitk/2) # muchi2cdf(x,nu) = gammainc(x/2, nu/2);
                    lc.append(1+hb[min(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nufitk)
                    hc.append(1+hb[max(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nufitk)
                if not plot:
                    for k in range(len(self.fitargs)):
                        write_csv(chi2fit[k],lc[k],hc[k],k)  # writes csv file
                        path = save_fit(k)  # saves .fit file
                        if path.__len__ != 2:
                            # assume path is a path string, whose length will be definitely > 2
                            self.console('chi2r = {:.4f} ({:.4f} - {:.4f}), saved in  {}'.format(chi2fit[k],lc[k],hc[k],path))
                        else:
                            # assume path is a tuple containing (path, exception)
                            self.console('Could not save results in  {}, error: {}'.format(path[0],path[1]))

                # print('len(self.fitargs)) = {}'.format(len(self.fitargs)))
                #########################################################
                # animate (see): TAKES CARE OF i>0 PLOTS IN 2D ARRAYS
                #                DOES ALSO UPDATE pars AND FIT CURVES
                #########################################################
                self.anim_fit = animation.FuncAnimation(self.fig_fit, 
                                                        animate, 
                                                        range(0,len(self.fitargs)),init_func=init,
                                                        interval=anim_delay.value,repeat=True,
                                                        blit=False) # 


            ###############################
            # single and tiles with offset
            ###############################
            ########
            # tile 
            ########
            else:  # TILES: creates matrices for offset multiple plots (does nothing on single)
                ##############################
                # THIS BLOCK TAKES CARE OF ALL ROWS OF DATA AT ONCE (errobars, fit curve, histograms and all)
                # pars must refer to the run of the FIRST row, both for global and multi fits 
                #      in anim therefore FIT CURVES  (f, fres, fl, flres) ARE ALWAYS 1D ARRAYS
                # animate must take care of updating parameters and producing correct fit curves               
                ##############
                # initial plot
                ##############
                yoffset = 0.05
                ymax = yoffset*fres.max()
                rmax = 0.3*(y-fres).max()
                xoffset = 0.
                # print ('fres = {}'.format(fres.shape))
                ttile, ytile, yres = plotile(t,y.shape[0],offset=xoffset), plotile(y,offset=ymax), plotile(y-fres,offset=rmax)  # plot arrays, full suite
                tftile, ftile  = plotile(tf,y.shape[0],offset=xoffset), plotile(f,offset=ymax)  
                # print('ttile.shape = {}, ytile.shape= {}, yres.shape = {}, tftile.shape = {}, ftile.shape = {}'.format(ttile.shape,ytile.shape,yres.shape,tftile.shape,ftile.shape))
                # print('f_tile = {}'.format(f_tile[0,0:50]))
                #############################
                # plot first (or only) window
                #############################
                # print(color)
                # errorbar does not plot multidim
                t1 = t.max()/20.
                t0 = np.array([0,t1])
                y0 = np.array([0,0])
                for k in range(y.shape[0]):
                    color = next(self.ax_fit[0,0]._get_lines.prop_cycler)['color']
                    self.ax_fit[(0,0)].errorbar(ttile[k],
                                                ytile[k],
                                                yerr=ey[k],
                                                fmt='o',
                                                elinewidth=1.0,ecolor=color,mec=color,mfc=color,
                                                ms=2.0,alpha=0.5) # data
                    self.ax_fit[(0,0)].plot(t0,y0,'-',lw=0.5,alpha=0.3,color=color)
                    if not self._single_:
                        self.ax_fit[(0,0)].text(t1,y0.max(),str(self.nrun[k]))
                    self.ax_fit[(1,0)].plot(ttile[k],yres[k],'-',lw=1.0,alpha=0.3,zorder=2,color=color) # residues 

                    self.ax_fit[(0,0)].plot(tftile[k],ftile[k],'-',lw=1.5,alpha=0.5,zorder=2,color=color) # fit
                    y0 = y0 + ymax
                self.ax_fit[(1,0)].plot(t[0],yzero,'k-',lw=0.5,alpha=0.3,zorder=0) # zero line

                ############################
                # plot second window, if any
                ############################
                if len(returntup)==5:
                    tltile, yltile, ylres  = plotile(tlate,xdim=ylate.shape[0],offset=xoffset), plotile(ylate,offset=xoffset), plotile(ylate-flres,offset=rmax)  # plot arrays, full suite
                    tfltile, fltile  = plotile(tfl,offset=xoffset),plotile(fl,
                                                          offset=foffset)  # res offset is 0.03 = 0.1-0.07
                    for k in range(y.shape[0]):
                        color = next(self.ax_fit[0,1]._get_lines.prop_cycler)['color']
                        self.ax_fit[(0,1)].errorbar(tltile[k],yltile[k],yerr=eylate[k],
                                                    fmt='o',elinewidth=1.0,
                                                    mec=color,mfc=color,ecolor=color,ms=2.0,alpha=0.5) # data
                        self.ax_fit[(1,1)].plot(tltile[k],ylres[k],'-',lw=1.0,alpha=0.3,zorder=2,color=color) # residues 
                        self.ax_fit[(0,1)].plot(tfl[0],fl_tile,'-',lw=1.5,alpha=0.5,zorder=2,color=color) # fit
                        self.ax_fit[(0,1)].set_xlim(tlate[0,0], tlate_tile[-1,-1])
                    self.ax_fit[(1,1)].plot(tlate,tlate-tlate,'k-',lw=0.5,alpha=0.3,zorder=0) # zero line

                ###############################
                # set title, labels
                ###############################
                self.ax_fit[(0,0)].set_ylabel('Asymmetry')
                self.ax_fit[(1,0)].set_ylabel('Residues')
                self.ax_fit[(1,0)].set_xlabel(r'Time [$\mu$s]')

                if self._single_:
                    self.ax_fit[(0,0)].set_title(str(self.nrun[0])+': '+self.title.value)
                    ########################
                    # chi2 distribution: fit
                    ########################
                    nufit,dum,chi2fit = chi(tfit[0],yfit[0],eyfit[0],pars[0])
                    nu,f,chi2plot = chi(t[0],y[0],ey[0],pars[0])
                    self.ax_fit[(0,0)].plot(t[0],f,'g--',lw=1.5	,alpha=1,zorder=2)#,color=color) # fit
                    xbin = np.linspace(-5.5,5.5,12)
                    self.ax_fit[(1,-1)].hist((yfit[0]-ffit[0])/eyfit[0],xbin,rwidth=0.9,fc='w',ec='k',lw=0.7)
                    # self.ax_fit[(1,-1)].set_ylim(0, 1.15*nhist.max())

                    #########################################
                    # chi2 distribution: plots, scaled to fit
                    #########################################
                    self.ax_fit[(1,-1)].hist((y[0]-fres[0])/ey[0],xbin,weights=nufit/nu*np.ones(t.shape[1]),rwidth=0.9,fc=color,alpha=0.2)

                    ###############################
                    # chi2 dist theo curve & labels 
                    ###############################
                    xh = np.linspace(-5.5,5.5,23)
                    yh = norm.cdf(xh+1)-norm.cdf(xh)
                    self.ax_fit[(1,-1)].plot(xh+0.5,nufit*yh,'r-')
                    self.ax_fit[(1,-1)].set_xlabel("$\sigma$")
                    self.ax_fit[(1,-1)].set_yticklabels(['']*len(self.ax_fit[(1,-1)].get_yticks()))    
                    self.ax_fit[(1,-1)].set_xlim([-5.5, 5.5])
                    mm = round(nu/4)
                    hb = np.linspace(-mm,mm,2*mm+1)
                    cc = gammainc((hb+nu)/2,nu/2) # muchi2cdf(x,nu) = gammainc(x/2, nu/2);
                    lc = 1+hb[min(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nufit
                    hc = 1+hb[max(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nufit
                    if not plot:
                        write_csv(chi2fit,lc,hc,0)  # writes csv file
                        path = save_fit(0)  # saves .fit file
                        if path.__len__ != 2:
                            # assume path is a path string, whose length will be definitely > 2
                            self.console('chi2r = {:.4f} ({:.4f} - {:.4f}), saved in  {}'.format(chi2fit,lc,hc,path))
                        else:
                            # assume path is a tuple containing (path, exception)
                            self.console('Could not save results in  {}, error: {}'.format(path[0],path[1]))
                    string = 'F-B: {} - {}\n'.format(self.group[0].value,self.group[1].value)
                    string += r'$\alpha=$ {}'.format(self.alpha.value)
                    string += '\n$\chi^2_f=$ {:.4f}\n ({:.2f}-{:.2f})\n$\chi^2_c=$ {:.4f}\n{} dof\n'.format(chi2fit,lc,hc,gammainc(chi2fit,nufit),nufit)
                    if  len(returntup)==5: 
                        string += '$\chi^2_e=$ {:.4f}\n$\chi^2_l=$ {:.4f}'.format(chi2plot,chi2late)
                    else:
                        string += '$\chi^2_p=$ {:.4f}'.format(chi2plot)
                    self.ax_fit[(0,-1)].text(-4.,0.3,string)
                else:
                    self.ax_fit[(0,0)].set_title(self.title.value)
                    ########################
                    # chi2 distribution: fit
                    ########################
                    fittup = derange(self.fit_range.value,self.histoLength) # range as tuple
                    fit_pack =1
                    if len(fittup)==3: # plot start stop pack
                        fit_start, fit_stop, fit_pack = fittup[0], fittup[1], fittup[2]
                    elif len(fittup)==2: # plot start stop
                        fit_start, fit_stop = fittup[0], fittup[1]
                    # if not self._single_ each run is a row in 2d ndarrays yfit, eyfit
                    # tfit,yfit,eyfit = rebin(self.time,self.asymm,[fit_start,fit_stop],fit_pack,e=self.asyme)
                    ychi = 0.
                    for k in range(len(pars)):
                        #########################################
                        # chi2 distribution: plots, scaled to fit
                        #########################################
                        nufit,ffit,chi2fit = chi(tfit[0],yfit[k],eyfit[k],pars[k])
                        nu,f,chi2plot = chi(t[0],y[k],ey[k],pars[k])
                        mm = round(nufit/4)
                        hb = np.linspace(-mm,mm,2*mm+1)
                        cc = gammainc((hb+nu)/2,nu/2) # muchi2cdf(x,nu) = gammainc(x/2, nu/2);
                        lc = 1+hb[min(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nufit
                        hc = 1+hb[max(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nufit
                        if not plot:
                            write_csv(chi2fit,lc,hc,k)  # writes csv file
                            path = save_fit(k)  # saves .fit file
                            if path.__len__ != 2:
                                # assume path is a path string, whose length will be definitely > 2
                                self.console('chi2r = {:.4f} ({:.4f} - {:.4f}), saved in  {}'.format(chi2fit,lc,hc,path))
                            else:
                                # assume path is a tuple containing (path, exception)
                                self.console('Could not save results in  {}, error: {}'.format(path[0],path[1]))
                        pedice = '_{'+str(self.nrun[k])+'}'
                        string = 'F-B: {} - {}\n'.format(self.group[0].value,self.group[1].value)
                        string += r'$\alpha=$ {}'.format(self.alpha.value)
                        string += '\n$\chi^2'+pedice+'=$ {:.3f}'.format(chi2fit)
                        self.ax_fit[(0,-1)].text(0.03,ychi,string)
                        ychi += ymax
                    self.ax_fit[(1,-1)].axis('off')
                    self.ax_fit[(0,-1)].set_ylim(self.ax_fit[(0,0)].get_ylim())
                self.ax_fit[(0,-1)].axis('off')
                # self.mainwindow.selected_index = 3  # focus on output tab

            self.fig_fit.canvas.manager.window.tkraise()
            P.draw()
                
        def int2_int():
            '''
            From internal parameters to the minimal representation 
            for the use of mucomponents._add_.
            Invoked just before submitting minuit 

            '''
            from mujpy.aux.aux import translate

#_components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value            
# refactor : this routine has much in common with min2int
            ntot = sum([len(self.model_components[k]['pars']) for k in range(len(self.model_components))])
            lmin = [-1]*ntot
            nint = -1 # initializec\c   
            nmin = -1 # initialize
            _int = []
            for k in range(len(self.model_components)):  # scan the model
                name = self.model_components[k]['name']
                # print('name = {}, model = {}'.format(name,self._the_model_))
                bndmthd = [] if name=='da' else self._the_model_.__getattribute__(name)
                # to set dalpha apart
                keys = []
                isminuit = []
                for j in range(len(self.model_components[k]['pars'])): # 
                    nint += 1  # internal parameter incremente always   
                    if self.flag[nint].value == '=': #  function is written in terms of nint
                        # nint must be translated into nmin 
                        string = translate(nint,lmin,self.function)
                        keys.append(string) # the function will be eval-uated, eval(key) inside mucomponents
                        isminuit.append(False)
                    else:
                        nmin += 1
                        keys.append('p['+str(nmin)+']')                        
                        lmin[nmin] = nint # lmin contains the int number of the minuit parameter
                        isminuit.append(True)
                _int.append([bndmthd,keys]) #,isminuit])  # ([component_dict,keys])
           # for k in range(len(_int)):
           #     print(_int[k])      

            return _int

        def int2min(return_names=False):
            '''
            From internal parameters to minuit parameters.
            Invoked just before submitting minuit 
            Internal are numbered progressively according to the display:
               first global parameters not belonging to components - e.g. A0, R, 
                        such as for asymmetry1 = A0*R and asymmetry2= A0*(1.-R)
               then local parameters not belonging to components - e.g. B and T
                        from the data file headers
               then the list of components' parameters
            Minuit parameters are the same, including fixed ones, but 
               the ones defined by functions or sharing
            Each parameter requires name=value, error_name=value, fix_name=value, limits_name=value,value
            [plus
               the local replica of the non global component parameters
               to be implemented]

            New version for suite of runs
                fitarg becomes a list of dictionaries

            '''
            ntot = sum([len(self.model_components[k]['pars']) for k in range(len(self.model_components))])
            ntot -= sum([1 for k in range(ntot) if self.flag[k]=='=']) # ntot minus number of functions 
            fitarg = [] # list of dictionaries
            parameter_names = []
##########################################
# single produces a list of one dictionary 
#            with keys 'par_name':guess_value,'error_par_name':step,...
# suite no global produces one dictionary per run
#       furthermore, if flag == 'l', each run may have a different guess value
# suite global produces a single dictionary but has
#       local, run dependent parameters, that may have flag=='l'
##########################################
            if not self._global_:

                for lrun in range(len(self._the_runs_)):
                    lmin = [-1]*ntot
                    nint = -1 # initialize
                    nmin = -1 # initialize
                    free = -1
                    fitargs= {}
                    for k in range(len(self.model_components)):  # scan the model
                        component_name = self.model_components[k]['name'] # name of component
                        keys = []
                        for j, par in enumerate(self.model_components[k]['pars']): # list of dictionaries, par is a dictionary
                            nint += 1  # internal parameter incremented always   
                            if self.flag[nint].value == '~': #  skip functions, they are not new minuit parameter
                                keys.append('~')
                                nmin += 1
                                free += 1
                                lmin[nmin] = nint # correspondence between nmin and nint, is it useful?
                                fitargs.update({par['name']:float(self.parvalue[nint].value)})
                                parameter_names.append(par['name'])
                                fitargs.update({'error_'+par['name']:float(par['error'])})
                                if not (par['limits'][0] == 0 and par['limits'][1] == 0):
                                    fitargs.update({'limit_'+par['name']:par['limits']})
                            elif self.flag[nint].value == 'l':
                                keys.append('~')
                                nmin += 1
                                free += 1
                                lmin[nmin] = nint # correspondence between nmin and nint, is it useful?
                                fitargs.update({par['name']:muvalue(lrun,self.function[nint].value)})
                                parameter_names.append(par['name'])
                                fitargs.update({'error_'+par['name']:float(par['error'])})
                                if not (par['limits'][0] == 0 and par['limits'][1] == 0):
                                    fitargs.update({'limit_'+par['name']:par['limits']})                            
                            elif self.flag[nint].value == '!':
                                nmin += 1
                                lmin[nmin] = nint # correspondence between nmin and nint, is it useful?
                                fitargs.update({par['name']:float(self.parvalue[nint].value)})
                                parameter_names.append(par['name'])
                                fitargs.update({'fix_'+par['name']:True})
                    fitarg.append(fitargs)
                    self.freepars = free

            else: # global
                # to be done
                fitarg.append(fitargs)


            # print('fitargs= {}'.format(fitargs))
            if return_names:
                self.minuit_parameter_names = tuple(parameter_names)
            return fitarg

        def load_fit(b):
            '''
            loads fit values such that the same fit can be reproduced on the same data

            '''
            import dill as pickle
            import os

            path_and_filename = path_file_dialog(self.paths[2].value) # returns the full path and filename
            if path_and_filename == '':
                return
            try:            
                with open(path_and_filename,'rb') as f:
                    fit_dict = pickle.load(f)
                try: 
                    del self._the_model_
                    self.fitargs = []
                except:
                    pass
                model.value = fit_dict['model.value']
                self.fit(model.value) # re-initialize the tab with a new model
                self.version.value = fit_dict['version']
                self.offset.value = fit_dict['self.offset.value']
                self.model_components = fit_dict['self.model_components']
                self.grouping = fit_dict['self.grouping']
                set_group()
                self.alpha.value = fit_dict['self.alpha.value']
                self.offset.value = fit_dict['self.offset.value']
                nint = fit_dict['nint']
                self.fit_range.value = fit_dict['self.fit_range.value']
                self.plot_range.value = fit_dict['self.plot_range.value'] # keys
                for k in range(nint+1):
                    self.parvalue[k].value = fit_dict['_parvalue['+str(k)+']']  
                    self.flag[k].value = fit_dict['_flag['+str(k)+    ']']  
                    self.function[k].value = fit_dict['_function['+str(k)+']']
                self.fitargs = fit_dict['self.fitargs'] 
                self.load_handle[0].value = fit_dict['self.load_handle[0].value']

            except Exception as e:
                self.console('Problems with reading {} file\n\nException: {}'.format(path_and_filename,e))

        def min2int(fitargs):
            '''
            From minuit parameters to internal parameters,
            see int2min for a description   
            Invoked just after minuit convergence for save_fit, [on_update]
            '''
            # refactor : this routine has much in common with int2_int
            # initialize
            from mujpy.aux.aux import translate

            ntot = sum([len(self.model_components[k]['pars']) for k in range(len(self.model_components))])
            _parvalue =  []
            lmin = [-1]*ntot 
            p = [0.0]*ntot 
            nint = -1
            nmin = -1
            for k in range(len(self.model_components)):  # scan the model
                keys = []
                for j, par in enumerate(self.model_components[k]['pars']): # list of dictionaries, par is a dictionary
                    nint += 1  # internal parameter incremented always   
                    if self.flag[nint].value != '=': #  skip functions, they are not new minuit parameter
                        nmin += 1
                        p[nmin] = fitargs[par['name']] # needed also by functions
                        _parvalue.append('{:4f}'.format(p[nmin]))   # _parvalue item is a string                     
                        lmin[nint] = nmin # number of minuit parameter
                    else: # functions, calculate as such
                        # nint must be translated into nmin 
                        string = translate(nint,lmin,self.function) # 
                        _parvalue.append('{:4f}'.format(eval(string))) # _parvalue item is a string
            return _parvalue
            
        def min2print(fitargs):
            '''
            From minuit parameters to plain print in xterminal ,
            see int2min for a description   
            Invoked just after minuit convergence for output 
            '''
            # refactor : this routine has much in common with int2_int
            # initialize
            from mujpy.aux.aux import translate, value_error

            ntot = sum([len(self.model_components[k]['pars']) for k in range(len(self.model_components))])
            _parvalue =  []
            lmin = [-1]*ntot 
            p = [0.0]*ntot 
            e = [0.0]*ntot 
            nint = -1
            nmin = -1
            for k in range(len(self.model_components)):  # scan the model
                keys = []
                for j, par in enumerate(self.model_components[k]['pars']): # list of dictionaries, par is a dictionary
                    nint += 1  # internal parameter incremented always   
                    if self.flag[nint].value != '=': #  skip functions, they are not new minuit parameter
                        nmin += 1
                        string = ' {}-{} = {}'.format(nmin,par['name'],value_error(fitargs[par['name']],fitargs['error_'+par['name']]))
                        self.console(string,end='') # needed also by functions
                        p[nmin] = fitargs[par['name']] # needed also by functions
                        lmin[nint] = nmin # number of minuit parameter
                        e[nmin] = fitargs['error_'+par['name']]
                    else: # functions, calculate as such
                        # nint must be translated into nmin 
                        string = translate(nint,lmin,self.function) # 
                        stringe =string.replace('p','e')
                        stringa = '  -{} = {}'.format(par['name'],value_error(eval(string),eval(stringe)))
                        self.console(stringa,end='') # needed also by functions
                self.console('')
            return

        def on_alpha_changed(change):
            '''
            observe response of fit tab widgets:
            validate float        
            '''
            string = change['owner'].value # description is three chars ('val','fun','flg') followed by an integer nint
                                               # iterable in range(ntot), total number of internal parameters
            try: 
                float(string)
            except:
                change['owner'].value = '{:.4f}'.format(alpha0)

        def on_fit_request(b):
            '''
            retrieve data from the gui dashboard:
            parameters values (parvalue[nint].value), flags (flag[nint].value), 
            errors, limits, functions (function[nint].value), self.alpha.value, range and pack
            pass _int, generated by int2_int. to mumodel._add_ (distribute minuit parameters)
            obtain fitargs dictionary, needed by migrad, either from self.fitargs or from min2int
            pass them to minuit
            call fit_plot
            save fit file in save_fit
            write summary in write_csv 

            '''
            from iminuit import Minuit as M
            from mujpy.aux.aux import derange, rebin, get_title
            # from time import localtime, strftime

    ###################
    # error: no run yet
    ###################
            if not self._the_runs_:
                self.console('No run loaded yet! Load one first (select suite tab).')
            else:
    ###################
    # run loaded
    ###################
                self.asymmetry() # prepare asymmetry
                # self.time is 1D asymm, asyme can 
                pack = 1 # initialize default
                returntup = derange(eval('self.fit_range.value'),self.histoLength) 
                # self.console('values  = {}'.format(returntup)) # debug
                if len(returntup)==3: # 
                    start, stop, pack = returntup
                elif len(returntup)==0:
                    self.console('Empty ranges. Choose fit/plot range')
                else:
                    start, stop = returntup
                time,asymm,asyme = rebin(self.time,self.asymm,[start,stop],pack,e=self.asyme)
                level = 0
                self.fitargs = []
                fitarg = int2min(return_names=True) # from dash
                if self._global_:
                    self._the_model_._load_data_(time[0],asymm,int2_int(),
                                                 self.alpha.value,e=asyme,
                                                 _nglobal_=self.nglobals,_locals_=self.locals) # pass all data to model                    
                    ##############################
                    # actual global migrad call
                    lastfit = M(self._the_model_._chisquare_,
                                     pedantic=False,
                                     forced_parameters=self.minuit_parameter_names,
                                     print_level=level,**fitarg[0])
                    # self.console('{} *****'.format([self.nrun[k] for k in range(len(self.nrun))]))
                    lastfit.migrad()
                    self.fitargs.append(lastfit.fitarg)
                    min2print(lastfit.fitarg) 
                    #    lastfit[0].hesse()
                    ##############################
                else:
                    # self.console('_single is {}'.format(self._single_)) # debug
                    if self._single_:
                        self._the_model_._load_data_(time[0],asymm[0],int2_int(),self.alpha.value,e=asyme[0]) # pass data to model, one at a time
                        ##############################
                        # actual single migrad calls
                        lastfit = M(self._the_model_._chisquare_,
                                         pedantic=False,
                                         forced_parameters=self.minuit_parameter_names,
                                         print_level=level,**fitarg[0])
                        lastfit.migrad()
                        # if lastfit.migrad_ok():
                        self.fitargs.append(lastfit.fitarg)                            
                        nu = len(time[0]) - self.freepars # degrees of freedom in plot
                        # self.freepars is calculated in int2min
                        chi2 = lastfit.fval/nu  
                        # self.console('{}: {} ***** chi2 = {:.3f} ***** {}'.format(self.nrun[0],get_title(self._the_runs_[0][0]),chi2,strftime("%d %b %Y %H:%M:%S", localtime())))
                        min2print(lastfit.fitarg)
                    else:
                        for k in range(len(self._the_runs_)): 
                            self._the_model_._load_data_(time[0],asymm[k],int2_int(),self.alpha.value,e=asyme[k]) # pass data to model, one at a time
                            ##############################
                            # actual single migrad calls
                            lastfit = M(self._the_model_._chisquare_,
                                             pedantic=False,
                                             forced_parameters=self.minuit_parameter_names,
                                             print_level=level,**fitarg[k])
                            lastfit.migrad()
                            # if lastfit.migrad_ok():
                            self.fitargs.append(lastfit.fitarg)
                            nu = len(time[0]) - self.freepars # degrees of freedom in plot
                            # self.freepars is calculated in int2min
                            chi2 = lastfit.fval/nu  
                            # self.console('{}: {} ***** chi2 = {:.3f} ***** {}'.format(self.nrun[k],get_title(self._the_runs_[k][0]),chi2,strftime("%d %b %Y %H:%M:%S", localtime())))
                            min2print(lastfit.fitarg)
                            #    lastfit.hesse()
                            ##############################
                    if lastfit.migrad_ok():
                        fitplot() # plot the best fit results 
                    else:
                        self.console('Fit not converged')
#                        string = 'Time last = {:.2e}, A[0] = {:.2f}, alpha = {:.3f}'.format(self._the_model_._x_[-1],self._themodel_._y_[0],self._the_model_._alpha_) # debug
#                        self.console(string)
         
        def on_flag_changed(change):
            '''
            observe response of fit tab widgets:
            set disabled on corresponding function (True if flag=='!' or '~', False if flag=='=') 
            '''
            dscr = change['owner'].description # description is three chars ('val','fun','flg') followed by an integer nint
                                               # iterable in range(ntot), total number of internal parameters
            n = int(dscr[4:]) # description='flag'+str(nint), skip 'flag'
            self.function[n].disabled=False if change['new']=='=' else True

        def on_function_changed(change):
            '''
            observe response of fit tab widgets:
            check for validity of function syntax

            '''
            from mujpy.aux.aux import muvalid

            dscr = change['owner'].description # description is three chars ('val','fun','flg') followed by an integer nint
                                               # iterable in range(ntot), total number of internal parameters
            n = int(dscr[4:]) # description='func'+str(nint), skip 'func'
            error_message = muvalid(change['new'])
            if error_message!='':
                self.function[n].value = ''
                self.console(error_message) 
  
        def on_group_changed(change):
            '''
            observe response of setup tab widgets:

            '''
            from mujpy.aux.aux import get_grouping
            name = change['owner'].description
            groups = ['forward','backward']
    #       now parse groupcsv shorthand
            self.grouping[name] = get_grouping(self.group[groups.index(name)].value) # stores self.group shorthand in self.grouping dict, grouping is the python based address of the counters
            # self.console('Group {} -> grouping {}'.format(self.group[groups.index(name)].value, get_grouping(self.group[groups.index(name)].value))) # debug
            if self.grouping[name][0]==-1:
                self.console('Wrong group syntax: {}'.format(self.group[groups.index(name)].value))
                self.group[groups.index(name)].value = ''
                self.grouping[name] = np.array([])
            else:
                self.get_totals()
        def on_integer(change):
            name = change['owner'].description
            if name == 'offset':
                if self.offset.value<0: # must be positive
                   self.offset.value = self.offset0 # standard value

        def on_load_model(change):
            '''
            observe response of fit tab widgets:
            check that change['new'] is a valid model
            relaunch MuJPy.fit(change['new'])
            '''
            if checkvalidmodel(change['new']): # empty list is False, non empty list is True
                try:
                    del self._the_model_
                    self.fitargs=[] # so that plot understands that ther is no previous minimization
                except:
                    pass
                self.fit(change['new']) # restart the gui with a new model
                self.mainwindow.selected_index = 2
            else:
                loadmodel.value=''

        def on_parvalue_changed(change):
            '''
            observe response of fit tab widgets:
            check for validity of function syntax
            '''
            dscr = change['owner'].description # description is three chars ('val','fun','flg') followed by an integer nint
                                               # iterable in range(ntot), total number of internal parameters
            n = int(dscr[5:]) # description='value'+str(nint), skip 'func'
            try:
                float(self.parvalue[n].value)
                self.parvalue[n].background_color = "white"
            except:
                self.parvalue[n].value = '0.0' 
                self.parvalue[n].background_color = "mistyrose"
 
        def on_plot_request(b):
            '''
            plot wrapper
            '''
            if not guesscheck.value and not self._the_model_._alpha_:
                self.console('No best fit yet, to plot the guess function tick the checkbox')

            else:
                fitplot(guess=guesscheck.value,plot=True) # 

        def on_range(change):
            '''
            observe response of FIT, PLOT range widgets:
            check for validity of function syntax
            '''
            from mujpy.aux.aux import derange

            fit_or_plot = change['owner'].description[0] # description is a long sentence starting with 'fit range' or 'plot range'
            if fit_or_plot=='f':
                name = 'fit'
            else:
                name = 'plot'
            returnedtup = derange(change['owner'].value,self.histoLength) 
            # print('sum = {}'.format(sum(returnedtup)))
            if sum(returnedtup)<0: # errors return (-1,-1), good values are all positive
                                   # go back to defaults
                if name == 'fit':
                    self.fit_range.value = '0,'+str(self.histoLength)
                    self.fit_range.background_color = "mistyrose"
                else:
                    self.plot_range.value = self.plot_range0
                    self.plot_range.background_color = "mistyrose"
            else:
                if name == 'fit':
                    self.fit_range.background_color = "white"
                    if len(returnedtup)==5:
                        if returnedtup[4]>self.histoLength:
                            change['owner'].value=str(returnedtup[:-1],self.histoLength)        
                    if returnedtup[1]>self.histoLength:
                        change['owner'].value=str(returnedtup[0],self.histoLength) if len(returnedtup)==2 else str(returnedtup[0],self.histoLength,returnedtup[2:])         
                else:
                    self.plot_range.background_color = "white"
                    if returnedtup[1]>self.histoLength:
                        change['owner'].value=str(returnedtup[0],self.histoLength) if len(returnedtup)==2 else str(returnedtup[0],self.histoLength,returnedtup[2])         

        def on_start_stop(change):
            if anim_check.value: 
                if change['new']:
                    self.anim_fit.event_source.start()
                else:
                    self.anim_fit.event_source.stop()

        def on_update(b):
            '''
            update parvalue[k].value with last best fit results
            '''
            if self.fitargs:
                _parvalue = min2int(self.fitargs[0]) # best fit parameters (strings)
                for k in range(len(_parvalue)):
                    self.parvalue[k].value = _parvalue[k]                  

        def path_file_dialog(path):
            import tkinter
            from tkinter import filedialog
            import os
            here = os.getcwd()
            os.chdir(path)
            tkinter.Tk().withdraw() # Close the root window
            in_path = filedialog.askopenfilename(filetypes=(('.fit','*.fit'),('all','*.*')))
            os.chdir(here)
            return in_path

        def save_fit(k):
            '''
            saves fit values such that load_fit can reproduce the same fit
            includes fit of suite of runs and global fits
            '''
            import dill as pickle
            import os

            version = str(self.version.value)
            fittype = '' # single run fit
            if self._global_: # global fit of run suite
                fittype = '.G.'
            elif not self._single_: # sequential fit of run suite
                fyttype = '.S.'
            strgrp = self.group[0].value.replace(',','_')+'-'+self.group[1].value.replace(',','_')
            path_fit = os.path.join(self.paths[2].value, model.value+'.'+version+fittype+'.'+str(self.nrun[k])+'.'+strgrp+'.fit')
            # create dictionary setup_dict to be pickled 
            # the inclusion of self.load_handle[0] will reload the data upon load_fit (?)
            names = ['self.alpha.value','self.offset.value',
                     'self.grouping','model.value',
                     'self.model_components','self.load_handle[0].value',
                     'version','nint',
                     'self.fit_range.value','self.plot_range.value',
                     'self.fitargs','self._global_','self._single_'] # keys
            fit_dict = {}
            for k,key in enumerate(names):
               fit_dict[names[k]] = eval(key) # key:value
            _parvalue = min2int(self.fitargs[0]) # starting values from first bestfit
            for k in range(nint+1):
                fit_dict['_parvalue['+str(k)+']'] = _parvalue[k] # either fit or dashboard
                fit_dict['_flag['+str(k)+    ']'] = self.flag[k].value # from fit tab
                fit_dict['_function['+str(k)+']'] = self.function[k].value # from fit tab
                
            with open(path_fit,'wb') as f:
                try:
                    # print ('dictionary to be saved: fit_dict = {}'.format(fit_dict))
                    pickle.dump(fit_dict, f) 
                except Exception as e:
                    return path_fit, e
            return path_fit

        def set_group():
            """
            return shorthand csv out of grouping
            name = 'forward' or 'backward'
            grouping[name] is an np.array wth counter indices
            group.value[k] for k=0,1 is a shorthand csv like '1:3,5' or '1,3,5' etc.
            """
            import numpy as np

            # two shorthands: either a list, comma separated, such as 1,3,5,6 
            # or a pair of integers, separated by a colon, such as 1:3 = 1,2,3 
            # only one column is allowed, but 1, 3, 5 , 7:9 = 1, 3, 5, 7, 8, 9 
            # or 1:3,5,7 = 1,2,3,5,7  are also valid
            #       get the shorthand from the gui Text 
            groups = ['forward','backward']
            for k, name in enumerate(groups):
                s = ''
                aux = np.split(self.grouping[name],np.where(np.diff(self.grouping[name]) != 1)[0]+1)
                # this finds the gaps in the array, 
                # e.g. 1,2,0,3 yields [array([1,2]),array([0]),array([3])]  
                for j in aux:                    
                    s += str(j[0]+1) # e.g '2', convention is from 1, python is from 0
                    if len(j)>1:
                        s += ':'+str(j[-1]+1) # e.g. '2:3'
                    s += ','
                s = s[:-1] # remove last ','
                self.group[k].value = s

        def write_csv(chi2,lowchi2,hichi2,k):
            '''
            writes a csv-like file of best fit parameters 
            that can be imported by qtiplot
            or read by python to produce figures::

                refactored for adding runs 
                and for writing one line per run
                in run suite, both local and global
                do not use csv, in order to control format (precision)

            '''
            import os
            #import csv
            from mujpy.aux.aux import get_title, spec_prec
            from time import localtime, strftime

            # print('k = {}, self.nrun = {}'.format(k,[j for j in self.nrun]))
            version = str(self.version.value)
            strgrp = self.group[0].value.replace(',','_')+'-'+self.group[1].value.replace(',','_')
            path_csv = os.path.join(self.paths[2].value,model.value+'.'+version+'.'+strgrp+'.csv')

            TsTc, eTsTc = self._the_runs_[k][0].get_temperatures_vector(), self._the_runs_[k][0].get_devTemperatures_vector()
            Bstr = self._the_runs_[k][0].get_field()
            n1,n2 = spec_prec(eTsTc[0]),spec_prec(eTsTc[1]) # calculates format specifier precision
            form = '{} {:.'
            form += '{}'.format(n1)
            form += 'f}  {:.'
            form += '{}'.format(n1)
            form += 'f}  {:.'
            form += '{}'.format(n2)
            form += 'f}  {:.'
            form += '{}'.format(n2)
            form += 'f} {}' #".format(value,most_significant)'
            # string = "{} {:.1f} {:.1f} {:.1f} {:.1f} {}".format(self.nrun[k], TsTc[0],eTsTc[0],TsTc[1],eTsTc[1], Bstr[:Bstr.find('G')]) #debug
            # self.console(string) # debug
            # self.console(form) #debug
            row = form.format(self.nrun[k], TsTc[0],eTsTc[0],TsTc[1],eTsTc[1], Bstr[:Bstr.find('G')])

            for name in self.minuit_parameter_names:
                value, error = self.fitargs[k][name], self.fitargs[k]['error_'+name]
                n1 = spec_prec(error) # calculates format specifier precision
                form = ' {:.'
                form += '{}'.format(n1)
                form += 'f}  {:.'
                form += '{}'.format(n1)
                form += 'f}'
                row += form.format(value,error)
            echi = max(chi2-lowchi2,hichi2-chi2)
            n1 = spec_prec(echi) # calculates format specifier precision
            form = ' {:.'
            form += '{}'.format(n1)
            form += 'f}  {:.'
            form += '{}'.format(n1)
            form += 'f}  {:.'
            form += '{}'.format(n1)
            form += 'f} {:.3f} {}'
            row += form.format(chi2,chi2-lowchi2,hichi2-chi2,self.alpha.value,self.offset.value)
            row += ' {}'.format(strftime("%d %b %Y %H:%M:%S", localtime()))
            form =' {} {:.2f}'
            for j in range(len(self.nt0)):
                row += form.format(self.nt0[j],self.dt0[j])
            row += '\n'
            # row is fromatted with appropriate rounding, write directly
            # self.console(row)
            header = ['Run','T_cryo[K]','e_T_cryo[K]','T_sample[K]','e_T_sample[K]','B[G]']
            for j,name in enumerate(self.minuit_parameter_names):
                header.append(name)
                header.append('e_'+name)
            header.append('chi2_r')
            header.append('e_chi2_low')
            header.append('e_chi2_hi')
            header.append('alpha')
            header.append('offset time')
            for j in range(len(self.nt0)):
                header.append('nt0{}'.format(j))
                header.append('dt0{}'.format(j))
            header = ' '.join(header)+'\n'

            try: # the file exists
                lineout = [] # is equivalent to False
                with open(path_csv,'r') as f_in:
                    notexistent = True
                    for nline,line in enumerate(f_in.readlines()):
                        if nline==0:
                            if header!=line: # different headers
                                break
                            lineout = [header]
                        elif int(line.split(" ")[0]) < self.nrun[k]:
                            lineout.append(line)
                        elif int(line.split(" ")[0]) == self.nrun[k]:
                            lineout.append(row) # substitute an existing fit
                            notexistent = False
                        else:
                            if notexistent:
                                lineout.append(row) # insert before last existing fit
                                notexistent = False
                            lineout.append(line)
                    if notexistent:
                        lineout.append(row) # append at the end
                        notexistent = False

                if not lineout:
                    raise # if headers were different this is an exception
                with open(path_csv,'w') as f_out:                 
                    for line in lineout:
                        f_out.write(line)
                self.console('Run {}: {} *** {}\nbest fit logged in {}'.format(self.nrun[k],get_title(self._the_runs_[k][0]),strftime("%d %b %Y %H:%M:%S", localtime()),path_csv))
            except: # write a new file
                with open(path_csv,'w') as f:
                    f.write(header)
                    f.write(row)
                self.console('Run {}: {} *** {}\nbest fit logged in NEW {}'.format(self.nrun[k],get_title(self._the_runs_[k][0]),strftime("%d %b %Y %H:%M:%S", localtime()),path_csv))

        ######### here starts the fit method of MuGui
        # no need to observe parvalue, since their value is a perfect storage point for the latest value
        # validity check before calling fit
        from ipywidgets import FloatText, Text, IntText, Layout, Button, HBox, \
                               Checkbox, VBox, Dropdown, ToggleButton, Label
        _available_components_() # creates tuple self.available_components automagically from mucomponents
        self._the_model_ = mumodel() # local instance, need a new one each time a fit tab is reloaded (on_load_model)
        try:
            alpha0 = self.alpha.value
        except:
            alpha0 = 1.01 # generic initial value
        try:
            self.offset0 = self.offset.value
        except:
            self.offset0 = 7 # generic initial value 
        loadbutton = Button(description='Load fit',layout=Layout(width='8%'))
        loadbutton.style.button_color = self.button_color
        loadbutton.on_click(load_fit)

        self.alpha = FloatText(description='alpha',value='{:.4f}'.format(alpha0),
                                layout=Layout(width='12%'),continuous_update=False) # self.alpha.value
        self.alpha.observe(on_alpha_changed,'value')
        self.offset = IntText(description='offset',value=self.offset0,
                                layout=Layout(width='11%'),continuous_update=False) # offset, is an integer
        # initialized to 7, only input is from an IntText, integer value, or saved and reloaded from mujpy_setup.pkl
        self.alpha.style.description_width='32%' 
        self.offset.style.description_width='38%' 
        # group and grouping: csv shorthand 
        self.group = [Text(description='forward',layout=Layout(width='16%'),
                                    continuous_update=False),
                      Text(description='backward',layout=Layout(width='16%'),
                                    continuous_update=False)]
        set_group() # inserts shorthand from self.grouping into seld.group[k].value, k=0,1
        self.group[0].observe(on_group_changed,'value')
        self.group[0].style.description_width='40%'
        self.group[1].style.description_width='40%'
        self.group[1].observe(on_group_changed,'value')
        guesscheck = Checkbox(description='guess',value=False, layout=Layout(width='8%'))
        guesscheck.style.description_width='2%'
        # end moved

        model = Text(description = '', layout=Layout(width='10%'), disabled = True) # this is static, empty description, next to loadmodel
        model.value = model_in
        version0 = 1
        loadmodel = Text(description='loadmodel',layout=Layout(width='20%'),continuous_update=False) # this is where one can input a new model name
        loadmodel.observe(on_load_model,'value')
        loadmodel.style.description_width='37%'

        try:
            version0 = self.version.value 
        except:
            version0 = 1

        try:
            self.plot_range0 = self.plot_range.value 
        except:
            if not self._the_runs_:
                self.plot_range0 = ''

        try:
            fit_range0 = self.fit_range.value 
        except:
            fit_range0 = self.plot_range0

        self.version = IntText(description='version',value=version0,layout=Layout(width='11%',indent=False)) # version.value is an int
        self.version.style.description_width='48%'

        fit_button = Button (description='Fit',layout=Layout(width='6%'))
        fit_button.style.button_color = self.button_color
        fit_button.on_click(on_fit_request)

        self.fit_range = Text(description='fit range\nstart,stop[,pack]',value=fit_range0,layout=Layout(width='22%'),continuous_update=False)
        self.fit_range.style.description_width='36%'
        self.fit_range.observe(on_range,'value')

        plot_button = Button (description='Plot',layout=Layout(width='6%'))
        plot_button.style.button_color = self.button_color
        plot_button.on_click(on_plot_request)

        self.plot_range = Text(description='plot range\nstart,stop\n[,pack]\n[last,pack]',value=self.plot_range0,
                               layout=Layout(width='22%'),continuous_update=False)
        self.plot_range.style.description_width='36%'
        self.plot_range.observe(on_range,'value')

        update_button = Button (description='Update',layout=Layout(width='8%'))
        update_button.style.button_color = self.button_color
        update_button.on_click(on_update)

        anim_stop_start = ToggleButton(description='stop/start',value=True,layout=Layout(width='10%'))   
        anim_stop_start.observe(on_start_stop,'value')

        label_delay = Label(value='Delay (ms)', layout=Layout(width='8%'))
        anim_delay = IntText(value=1000, layout=Layout(width='8%'))

        anim_check = Checkbox(description='Animate',value=False, layout=Layout(width='10%'))
        anim_check.style.description_width = '1%'
        
        topframe_handle = HBox(description = 'Model', children=[update_button,loadmodel,self.offset,
                                                                fit_button,self.group[1],
                                                                self.fit_range,anim_check,anim_delay])  #
        alphaframe_handle = HBox(description = 'Alpha', children=[loadbutton,guesscheck,self.alpha,
                                                                self.version,plot_button,self.group[0],
                                                                self.plot_range,anim_stop_start,
                                                                label_delay]) # 
        bottomframe_handle = HBox(description = 'Components', layout=Layout(width='100%',border='solid')) #

        try: 
            create_model(model.value) # this may be not a valid model, e.g. after fit('da#ò')
        except:
            self.fit()  # this starts over, producing model = 'daml', which is valid.

        leftframe_list, rightframe_list = [],[]
  
        words  = ['#','name','value','~!=','function']
        nint = -1 # internal parameter count, each widget its unique name
        ntot = np.array([len(self.model_components[k]['pars']) 
                         for k in range(len(self.model_components))]).sum()
        self.parvalue, self.flag, self.function = [], [], [] # lists, index runs according to internal parameter count nint
        self.compar = {} # dictionary: key nint corresponds to a list of two values, c (int index of component) and p (int index of parameter)
                # use: self.compar[nint] is a list of two integers, the component index k and its parameter index j
        self.fftcheck = []            

        for k in range(len(self.model_components)):  # scan the model
            self.fftcheck.append(Checkbox(description='FFT',value=True))
            header = HBox([ Text(value=self.model_components[k]['name'],disabled=True,layout=Layout(width='8%')),
                            self.fftcheck[k]]) # list of HBoxes, the first is the header for the component
                                                                     # composed of the name (e.g. 'da') and the FFT flag
                                                                     # fft will be applied to a 'residue' where only checked components
                                                                     # are subtracted
            componentframe_list = [header] # list of HBoxes, header and pars
            componentframe_handle = VBox()
            for j in range(len(self.model_components[k]['pars'])): # make a new par for each parameter 
                                                                                # and append it to component_frame_content
                nint += 1      # all parameters are internal parameters, first is pythonically zero 
                self.compar.update({nint:[k,j]}) # stores the correspondence between nint and component,parameter
                nintlabel_handle = Text(value=str(nint),layout=Layout(width='10%'),disabled=True)
                parname_handle = Text(value=self.model_components[k]['pars'][j]['name'],layout=Layout(width='22%'),disabled=True)
                # parname can be overwritten, not important to store

                self.parvalue.append(Text(value='{:.4}'.format(self.model_components[k]['pars'][j]['value']),
                                     layout=Layout(width='18%'),description='value'+str(nint),continuous_update=False))
                self.parvalue[nint].style.description_width='0%'
                try:
                    self.parvalue[nint].value = _parvalue[nint]
                except:
                    pass
                # parvalue handle must be unique and stored at position nint, it will provide the initial guess for the fit

                self.function.append(Text(value=self.model_components[k]['pars'][j]['function'],
                                     layout=Layout(width='33%'),description='func'+str(nint),continuous_update=False))
                self.function[nint].style.description_width='0%'
                try:
                    self.function[nint].value = _function[nint]
                except:
                    pass
                # function handle must be unique and stored at position nint, it will provide (eventually) the nonlinear relation 

                fdis = False if self.model_components[k]['pars'][j]['flag']=='=' else True 
                self.function[nint].disabled = fdis # enabled only if flag='='
                self.flag.append(Dropdown(options=['~','!','='], 
                                     value=self.model_components[k]['pars'][j]['flag'],
                                     layout=Layout(width='10%'),description='flag'+str(nint)))
                self.flag[nint].style.description_width='0%'
                try:
                    self.flag[nint].value = _flag[nint]
                except:
                    pass
                 # flag handle must be unique and stored at position nint, it will provide (eventually) the nonlinear relation to be evaluated
 
                # now put this set of parameter widgets for the new parameter inside an HBox
                par_handle = HBox([nintlabel_handle, parname_handle, self.parvalue[nint], self.flag[nint], self.function[nint]])
                           # handle to an HBox of a list of handles; notice that parvalue, flag and function are lists of handles
                
                # now make value flag and function active 
                self.parvalue[nint].observe(on_parvalue_changed,'value')
                self.flag[nint].observe(on_flag_changed,'value') # when flag[nint] is modified, function[nint] is z(de)activated
                self.function[nint].observe(on_function_changed,'value') # when function[nint] is modified, it is validated

                componentframe_list.append(par_handle) # add par widget to the frame list

            componentframe_handle.children = componentframe_list # add full component to the frame
            if k%2==0:                                         # and ...
                leftframe_list.append(componentframe_handle) # append it to the left if k even
            else:
                rightframe_list.append(componentframe_handle) # or to the right if k odd   

        # end of model scan, ad two vertical component boxes to the bottom frame
        bottomframe_handle.children = [VBox(leftframe_list),VBox(rightframe_list)]  # list of handles 

        # backdoors
        self._load_fit = load_fit
        self._fit = fitplot
        self._int2_int = int2_int

        # now collect the handles of the three horizontal frames to the main fit window 
        self.mainwindow.children[2].children = [alphaframe_handle, topframe_handle ,bottomframe_handle] 
                             # add the list of widget handles as the third tab, fit

##########################
# OUTPUT
##########################
    def output(self):
        '''
        create an Output widget in fourth tab       
        select by 
        self.mainwindow.selected_index = 3
        '''
        import os
        import platform
        from shutil import which
        from subprocess import Popen, PIPE
        from ipywidgets import Output, HBox, Layout 
                     # Output(layout={'height': '100px', 'overflow_y': 'auto', 'overflow_x': 'auto'})


        self._outputtab_ = Output(layout={'height': '300px','width':'100%','overflow_y':'auto','overflow_x':'auto'})
        _output_box = HBox([self._outputtab_],layout=Layout(width='100%')) # x works y does scroll
        self.mainwindow.children[3].children = [_output_box] 
                             # add the list of widget handles as the fourth tab, output
        if platform.system()=='Linux':
            if which('gnome-terminal') is not None:
                self._output_ = "/tmp/mujpy_pipe"
                if not os.path.exists(self._output_):
                    os.mkfifo(self._output_)
                terminal = 'gnome-terminal'
                settitle = "PROMPT_COMMAND='echo -ne \"\033]0;mujpy console\007\"'"
                # this is to open the xterm mujpy console on linux
                self._xterm_ = Popen([terminal,'shell=False',settitle,'-e', 'tail -f %s' % self._output_])
            else:
               self._output_ = self._outputtab_ # all in the tab
        elif platform.system()=='Windows': # put here also a cmd output option
            # this Windows and the Mac part bust be still checked
            self._xterm_ = Popen(["cmd.exe","/k"],stdout=PIPE,stderr=PIPE)
#            self._output_ = self._outputtab_ # all in the tab
        elif platform.system()=='Darwin':  # put here also a xterm option
            self._xterm_ = Popen(['open', '-a', 'Terminal', '-n'],stdout=PIPE,stderr=PIPE)
#            self._output_ = self._outputtab_ # all in the tab
        self.console(''.join(['***************************************************\n',
                              '* This is where the output and the error messages *\n'
                              '* will be displayed. Keep an eye. DO NOT KILL!    *\n',
                              '***************************************************\n']))
            

################
# PLOTS
################
    def plots(self):
        '''
        tlog plot
        multi plot (if not _single_)
        '''

        def on_counter(b):
            '''
            check syntax of counter_range
            '''
            from mujpy.aux.aux import get_grouping
            from numpy import array
            # abuse of get_grouping: same syntax here
            if counter_range.value == '':
                return
            counters = get_grouping(counter_range.value)
            ok = 0
            for k in range(counters.shape[0]):
                if counters[k]<0 or counters[k]>=self._the_runs_[0][0].get_numberHisto_int():
                    # print('k = {}, counters[k] = {}, numberHisto = {}'.format(k,counters[k],self._the_runs_[0][0].get_numberHisto_int()))
                    ok = -1
            if counters[0] == -1 or ok == -1:
                    #print('Wrong counter syntax or counters out of range: {}'.format(counter_range.value))
                self.console('Wrong counter syntax or counters out of range: {}'.format(counter_range.value))
                counter_range.value = ''
                counters = array([])

        def on_counterplot(b):
            '''
            COUNTERPLOT:
            produce plot
            '''
            from numpy import zeros, arange
            from mujpy.aux.aux import get_grouping, derange
            import matplotlib.pyplot as P

            font = {'family':'Ubuntu','size':8}
            P.rc('font', **font)
            dpi = 100.

            if not self._the_runs_:
                self.console('No run loaded yet! Load one first (select suite tab).')

            ############
            #  bin range
            ############
            returntup = derange(self.counterplot_range.value,self.histoLength) # 
            start, stop = returntup
            # abuse of get_grouping: same syntax here
            counters = get_grouping(counter_range.value) # already tested
            # now counters is an np.array of counter indices

            #############
            # load histos
            #############
            histo = zeros((self._the_runs_[0][0].get_numberHisto_int(),stop-start),dtype=int)
            bins = arange(start,stop,dtype=int)

            # 4x4, 3x3 or 2x3 counters
            ncounters = counters.shape[0] # self.numberHisto
            screen_x, screen_y = P.get_current_fig_manager().window.wm_maxsize() # screen size in pixels
            y_maxinch = float(screen_y)/dpi -0.5 # maximum y size in inches, 1 inch for window decorations 
            fx, f, f1 = 1., 4./5., 16./25. # fraction of screen for 
            if ncounters > 9:
                nrows,ncols = 4,4
                x,y = fx*y_maxinch, y_maxinch
            elif ncounters > 6:
                nrows,ncols = 3,3
                x,y = fx*y_maxinch*f, y_maxinch*f
            elif ncounters > 4:
                nrows,ncols = 2,3
                x,y = fx*y_maxinch*f, y_maxinch*f1
            elif ncounters > 1:
                nrows,ncols = 2,2
                x,y = fx*y_maxinch*f1, y_maxinch*f1
            else:
                nrows,ncols = 1,1

            ##############################
            #  set or recover figure, axes 
            ##############################
            if self.fig_counters:
                self.fig_counters.clf()
                self.fig_counters,self.ax_counters = P.subplots(nrows,ncols,figsize=(x,y),num=self.fig_counters.number) 
            else:
                # residual problem: when this is the first pyplot window a Figure 1 is opened that nobody ordered
                self.fig_counters,self.ax_counters = P.subplots(num=10,nrows=nrows,ncols=ncols,figsize=(x,y),dpi=dpi,squeeze=False)
                self.fig_counters.subplots_adjust(hspace=0.1,top=0.95,bottom=0.11,right=0.98,wspace=0.28)
                self.fig_counters.canvas.set_window_title('Counters')

            nplots = nrows*ncols    
            for kk,nrun in enumerate(self.nrun):
                if nrun==int(self.choose_nrun.value):
                    this_run = kk # if not single selects which run to plot counters for 
            for k in range(nplots):

                if k < counters.shape[0]: 
                    counter = counters[k] # already an index 0:n-1
                    for run in self._the_runs_[this_run]: # allow for add runs 
                        histo[counter] += run.get_histo_array_int(counter)[start:stop]
                    ymax = histo[counter].max()
                    if stop-start<100:
                        self.ax_counters[divmod(counter,ncols)].bar(bins,histo[counter,:],edgecolor='k',color='silver',alpha=0.7,lw=0.7)
                    else:
                        self.ax_counters[divmod(counter,ncols)].plot(bins,histo[counter,:],'k-',lw=0.7)
                    if divmod(counter,ncols)[0]==counters.shape[0]/ncols-1:
                        self.ax_counters[divmod(counter,ncols)].set_xlabel('bins')
                        # does not print
                    if divmod(counter,ncols)[1]==0:
                        self.ax_counters[divmod(counter,ncols)].set_ylabel('counts')
                    self.ax_counters[divmod(counter,ncols)].text(start+(stop-start)*0.9, ymax*0.9,'# '+str(counter+1)) # from index to label
                else:
                    self.ax_counters[divmod(k,ncols)].cla()
                    self.ax_counters[divmod(k,ncols)].axis('off')
            P.show()          
            
        def on_multiplot(b):
            '''
            MULTIPLOT:
            produce plot
            '''
            import matplotlib.pyplot as P
            from numpy import array
            from mujpy.aux.aux import derange, rebin, get_title
            import matplotlib.animation as animation

            ###################
            # PYPLOT ANIMATIONS
            ###################
            def animate(i):
                '''
                anim function
                update multiplot data and its color 
                '''
                line.set_ydata(asymm[i])
                line.set_color(color[i])
                self.ax_multiplot.set_title(str(self.nrun[i])+': '+get_title(self._the_runs_[0][0]))
                return line, 


            def init():
                '''
                anim init function
                to give a clean slate 
                '''
                line.set_ydata(asymm[0])
                line.set_color(color[0])
                self.ax_multiplot.set_title(str(self.nrun[0])+': '+get_title(self._the_runs_[0][0]))
                return line, 

            dpi = 100.
            ############
            #  bin range
            ############
            returntup = derange(self.multiplot_range.value,self.histoLength) # 
            pack = 1
            if len(returntup)==3: # plot start stop packearly last packlate
                start, stop, pack = returntup
            else:
                start, stop = returntup

            ####################
            # load and rebin 
            #     time,asymm are 2D arrays, 
            #     e.g. time.shape = (1,25000), 
            #     asymm.shape = (nruns,25000) 
            ###################
            self.asymmetry() # prepare asymmetry
            time,asymm = rebin(self.time,self.asymm,[start,stop],pack)
            nruns,nbins = asymm.shape

            #print('start, stop, pack = {},{},{}'.format(start,stop,pack))
            #print('shape time {}, asymm {}'.format(time.shape,asymm.shape))
            y = 4. # normal y size in inches
            x = 6. # normal x size in inches
            my = 12. # try not to go beyond 12 run plots

            ##############################
            #  set or recover figure, axes 
            ##############################
            if self.fig_multiplot:
                self.fig_multiplot.clf()
                self.fig_multiplot,self.ax_multiplot = P.subplots(figsize=(x,y),num=self.fig_multiplot.number) 
            else:
                self.fig_multiplot,self.ax_multiplot = P.subplots(figsize=(x,y),dpi=dpi)
                self.fig_multiplot.canvas.set_window_title('Multiplot')
            screen_x, screen_y = P.get_current_fig_manager().window.wm_maxsize() # screen size in pixels
            y_maxinch = float(screen_y)/float(self.fig_multiplot.dpi) # maximum y size in inches

            ########## note that "inches" are conventional, dince they depend on the display pitch  
            # print('your display is y_maxinch = {:.2f} inches'.format(y_maxinch))
            ########## XPS 13 is 10.5 "inches" high @160 ppi (cfr. conventional self.fig_multiplot.dpi = 100)
            bars = 1. # overhead y size(inches) for three bars (tools, window and icons)
            dy = 0. if anim_check.value else (y_maxinch-y-1)/my   # extra y size per run plot
            y = y + nruns*dy if nruns < 12 else y + 12*dy # size, does not dilate for anim 
            # self.fig_multiplot.set_size_inches(x,y, forward=True)

            ##########################
            #  plot data and fit curve
            ##########################
            color = []
            for run in range(nruns):
                color.append(next(self.ax_multiplot._get_lines.prop_cycler)['color'])

            if anim_check.value and not self._single_:
            #############
            # animation
            #############
                ##############
                # initial plot
                ##############
                ylow, yhigh = asymm.min()*1.02, asymm.max()*1.02
                line, = self.ax_multiplot.plot(time[0],asymm[0],'o-',ms=2,lw=0.5,color=color[0],alpha=0.5,zorder=1)
                self.ax_multiplot.set_title(str(self.nrun[0])+': '+get_title(self._the_runs_[0][0]))
                self.ax_multiplot.plot([time[0,0],time[0,-1]],[0,0],'k-',lw=0.5,alpha=0.3)
                self.ax_multiplot.set_xlim(time[0,0],time[0,-1])
                self.ax_multiplot.set_ylim(ylow,yhigh)
                self.ax_multiplot.set_ylabel('Asymmetry')
                self.ax_multiplot.set_xlabel(r'time [$\mu$s]')
                #######
                # anim
                #######
                self.anim_multiplot = animation.FuncAnimation(self.fig_multiplot, animate, nruns, init_func=init,
                      interval=anim_delay.value, blit=False)

            ###############################
            # tiles with offset
            ###############################
            else: 
                aoffset = asymm.max()*float(multiplot_offset.value)*array([[run] for run in range(nruns)])
                asymm = asymm + aoffset # exploits numpy broadcasting
                ylow,yhigh = min([0,asymm.min()+0.01]),asymm.max()+0.01
                for run in range(nruns):
                    self.ax_multiplot.plot(time[0],asymm[run],'o-',lw=0.5,ms=2,alpha=0.5,color=color[run],zorder=1)
                    self.ax_multiplot.plot([time[0,0],time[0,-1]],
                                           [aoffset[run],aoffset[run]],'k-',lw=0.5,alpha=0.3,zorder=0)
                    self.ax_multiplot.text(time[0,-1]*1.025,aoffset[run],self._the_runs_[run][0].get_runNumber_int())
                self.ax_multiplot.set_title(get_title(self._the_runs_[0][0]))
                self.ax_multiplot.set_xlim(time[0,0],time[0,-1]*9./8.)
                self.ax_multiplot.set_ylim(ylow,yhigh)
                # print('axis = [{},{},{},{}]'.format(time[0,0],time[0,-1]*9./8.,ylow,yhigh))
                self.ax_multiplot.set_ylabel('Asymmetry')
                self.ax_multiplot.set_xlabel(r'time [$\mu$s]')
                # self.fig_multiplot.tight_layout()
            self.fig_multiplot.canvas.manager.window.tkraise()
            P.show()

        def on_range(change):
            '''
            observe response of MULTIPLOT range widgets:
            check for validity of function syntax
            on_range (PLOTS, FIT, FFT) perhaps made universal and moved to aux
            '''
            from mujpy.aux.aux import derange

            # change['owner'].description
            if change['owner'].description[0:2] == 'fft':
                fmax = 0.5/(self.time[0,1]-self.time[0,0])
                returnedtup = derange(change['owner'].value,fmax,int_or_float='float')
            else:
                returnedtup = derange(change['owner'].value,self.histoLength) # errors return (-1,-1),(-1,0),(0,-1), good values are all positive
            if sum(returnedtup)<0:
                change['owner'].background_color = "mistyrose"
                if change['owner'].description[0:3] == 'plot':
                    change['owner'].value = self.plot_range0
                else:
                    change['owner'].value = self.bin_range0
            else:
                change['owner'].background_color = "white"
                if returnedtup[1]>self.histoLength:
                    change['owner'].value=str(returnedtup[0],self.histoLength) if len(returnedtup)==2 else str(returnedtup[0],self.histoLength,returnedtup[2])         

        def on_nexttlog(b):
            '''
            select next run tlog
            '''
            runs = self.choose_tlogrun.options.keys()
            runs = sorted([int(run) for run in runs])
            runindex = runs.index(self.choose_tlogrun.options.index(self.choose_tlogrun.value))
            if runindex > 0:
                self.choose_tlogrun.value = str(runs.index(runindex-1))
                on_tlogdisplay([])

        def on_prevtlog(b):
            '''
            select prev run tlog
            '''
            runs = self.choose_tlogrun.options.keys()
            runs = sorted([int(run) for run in runs])
            runindex = runs.index(self.choose_tlogrun.options.index(self.choose_tlogrun.value))
            if runindex < len(runs):
                self.choose_tlogrun.value = str(runs.index(runindex+1))
                on_tlogdisplay([])

        def on_start_stop(change):
            if anim_check.value: 
                if change['new']:
                    self.anim_multiplot.event_source.start()
                    anim_step.style.button_color = self.button_color_off
                    anim_step.disabled=True
                else:
                    self.anim_multiplot.event_source.stop()
                    anim_step.style.button_color = self.button_color
                    anim_step.disabled=False

        def on_step(b):
            '''
            step when stop animate
            ''' 
            if not anim_step.disabled:
                self.anim_multiplot.event_source.start()

        def on_tlogdisplay(b):
            '''
            display a PSI tlog if the files exist
            '''
            import os
            import matplotlib.pyplot as P
            from mujpy.aux.aux import muzeropad
            from numpy import array, mean, std
            import csv
            from matplotlib.dates import HourLocator, MinuteLocator, DateFormatter
            import datetime
            from datetime import datetime as d

            ################
            # load tlog file
            ################
            pathfile = os.path.join(self.paths[1].value,'run_'+muzeropad(self.choose_tlogrun.value)+'.mon')
            with open(pathfile,'r') as f:
                reader=csv.reader(f)
                header, t, T1,T2,pause,go = [],[],[],[],[],[]
                for k in range(9):
                    header.append(next(reader))
                # print(header[7][0][2:22])
                starttime = header[7][0][2:22]
                start = d.strptime(starttime, "%d-%b-%Y %H:%M:%S")
                # print(start)
                for row in reader:
                    # print(row)
                    if row[0][0]!='!':
                        row = row[0].split('\\')
                        stop = d.strptime(row[0], "%H:%M:%S")
                        row = row[2].split()
                        T1.append(float(row[0]))
                        T2.append(float(row[1]))
                        # print('record = {}, stop = {}'.format(row[0][0:8],stop.time()))
                        time = stop.time()
                        t.append(time.hour*60.+time.minute+time.second/60.)
                        # print('{}'.format(t))
                    else:
                        # print(row)
                        if row[0][24:29]=='Paused':
                            pause.append(d.strptime(row[0][2:22], "%d-%b-%Y %H:%M:%S"))
                        else:
                            go.append(d.strptime(row[0][2:22], "%d-%b-%Y %H:%M:%S"))
            ##############################
            #  set or recover figure, axes 
            ##############################
            if self.fig_tlog:
                self.fig_tlog.clf()
                self.fig_tlog,self.ax_tlog = P.subplots(num=self.fig_tlog.number) 
            else:
                self.fig_tlog,self.ax_tlog = P.subplots()
                self.fig_tlog.canvas.set_window_title('Tlogger')
            T1,T2 = array(T1), array(T2)
            self.ax_tlog.plot(t,T1,'r-',label='T1')
            self.ax_tlog.plot(t,T2,'b-',label='T2')
            tlim,Tlim = self.ax_tlog.get_xlim(), self.ax_tlog.get_ylim()
            T1ave, T1std, T2ave, T2std = mean(T1),std(T1),mean(T2),std(T2)
            self.ax_tlog.plot(tlim,[T1ave, T1ave],'r-',lw=0.5,alpha=0.8)
            self.ax_tlog.fill_between(tlim, [T1ave-T1std,T1ave-T1std ],[T1ave, T1ave],facecolor='r',alpha=0.2)
            self.ax_tlog.fill_between(tlim, [T1ave+T1std,T1ave+T1std ],[T1ave, T1ave],facecolor='r',alpha=0.2)
            self.ax_tlog.plot(tlim,[T2ave, T2ave],'b-',lw=0.5,alpha=0.8)
            self.ax_tlog.fill_between(tlim, [T2ave-T2std,T2ave-T2std ],[T2ave, T2ave],facecolor='b',alpha=0.2)
            self.ax_tlog.fill_between(tlim, [T2ave+T2std,T2ave+T2std ],[T2ave, T2ave],facecolor='b',alpha=0.2)
            self.ax_tlog.set_title('Run '+self.choose_tlogrun.value+' start at '+starttime)
            self.ax_tlog.set_xlabel('time [min]')
            self.ax_tlog.set_ylabel('T [k]')
            self.ax_tlog.legend()
            if Tlim[1]-Tlim[0] < 1.:
                T0 = (Tlim[0]+Tlim[1])/2
                self.ax_tlog.set_ylim(T0-1.,T0+0.5)
            y1,y2 = self.ax_tlog.get_ylim()
            for x1,x2 in zip(pause,go):
                self.ax_tlog.fill_between([x1,x2], [y1,y1 ],[y2,y2],facecolor='k',alpha=0.5)
                self.ax_tlog.text(x1*0.9+x2*0.1,y1*0.9+y2*0.1,'PAUSE',color='w')
            P.show()

        from ipywidgets import HBox, VBox, Button, Text, Textarea, Accordion, Layout, Checkbox, IntText, ToggleButton, Label, Dropdown
        ###########
        # multiplot
        ###########
        multiplot_button = Button(description='Multiplot',layout=Layout(width='10%'))
        multiplot_button.on_click(on_multiplot)
        multiplot_button.style.button_color = self.button_color
        anim_check = Checkbox(description='Animate',value=False, layout=Layout(width='10%'))
        anim_check.style.description_width = '1%'
        anim_delay = IntText(description='Delay (ms)',value=1000, layout=Layout(width='20%'))
        anim_delay.style.description_width = '45%'
        anim_stop_start = ToggleButton(description='start/stop',value=True,layout={'width':'12%'})
        anim_stop_start.observe(on_start_stop,'value')
        # anim_stop_start.style.button_color = self.button_color
        anim_step = Button(description='step',layout={'width':'10%'})
        anim_step.on_click(on_step)
        anim_step.style.button_color = self.button_color_off

        self.multiplot_range = Text(description='plot range\nstart,stop[,pack]',
                               value=self.plot_range0,layout=Layout(width='26%'),
                               continuous_update=False)
        self.multiplot_range.style.description_width='43%'
        self.multiplot_range.observe(on_range,'value')
        multiplot_offset0 = '0.1'
        multiplot_offset = Text(description='offset',
                               value=multiplot_offset0,layout=Layout(width='12%'),
                               continuous_update=False)
        multiplot_offset.style.description_width='35%'
        # self.tlog_accordion.layout.height='10'

        multibox = HBox([multiplot_button,anim_check,anim_delay,anim_stop_start,self.multiplot_range,multiplot_offset,
                         Label(layout=Layout(width='3%'))],layout=Layout(width='100%',border='solid'))
        ###################
        # counters inspect
        ###################
        counterlabel = Label(value='Inspect',layout=Layout(width='7%'))# count
        counterplot_button = Button(description='counters',layout=Layout(width='10%'))
        counterplot_button.on_click(on_counterplot)
        counterplot_button.style.button_color = self.button_color
        self.counternumber = Label(value='{} counters per run'.format(' '),layout=Layout(width='15%'))
        counter_range = Text(description='counters',
                               value='', continuous_update=False,layout=Layout(width='20%'))
        counter_range.style.description_width='33%'
        counter_range.observe(on_counter,'value')
        self.counterplot_range = Text(description='bins: start,stop',
                               value=self.bin_range0,continuous_update=False,layout=Layout(width='25%'))
        self.counterplot_range.style.description_width='40%'
        self.counterplot_range.observe(on_range,'value')
        self.choose_nrun = Dropdown(options=[], description='run', layout=Layout(width='15%'))
        self.choose_nrun.style.description_width='25%'
        counterbox = HBox([Label(layout=Layout(width='3%')),counterlabel,counterplot_button,Label(layout=Layout(width='3%')),self.counternumber,
                          counter_range,self.counterplot_range,self.choose_nrun],layout=Layout(width='100%',border='solid'))

        ##########
        # TLOG PSI
        ##########
        tloglabel = Label(value='Tlog',layout=Layout(width='7%'))
        tlog_button = Button(description='display',layout=Layout(width='10%'))
        tlog_button.on_click(on_tlogdisplay)
        tlog_button.style.button_color = self.button_color
        options = {} # empty slot to start with
        self.choose_tlogrun = Dropdown(options=options,description='Tlog run', layout=Layout(width='15%')) #
        self.choose_tlogrun.style.description_width='25%'
        nexttlog_button = Button(description='Next',layout=Layout(width='10%'))
        nexttlog_button.on_click(on_nexttlog)
        nexttlog_button.style.button_color = self.button_color
        prevtlog_button = Button(description='Prev',layout=Layout(width='10%'))
        prevtlog_button.on_click(on_prevtlog)
        prevtlog_button.style.button_color = self.button_color      

        self.tlog_accordion = Accordion(children=[Textarea(layout={'width':'100%','height':'200px',
                                                 'overflow_y':'auto','overflow_x':'auto'})])
        self.tlog_accordion.set_title(0,'run: T(eT)')
        self.tlog_accordion.selected_index = None

        tlogbox = HBox([Label(layout=Layout(width='3%')),tloglabel,tlog_button,Label(layout=Layout(width='3%')),self.choose_tlogrun,
                         nexttlog_button,prevtlog_button, self.tlog_accordion],layout=Layout(width='100%',border='solid'))

        vbox = VBox()
        vbox.children = [multibox, counterbox, tlogbox]
        self.mainwindow.children[5].children = [vbox]

##########################i
# SETUP
##########################
    def setup(self):
        '''
        setup tab of mugui
        used to set: 

            paths, fileprefix and extension
            prepeak, postpeak (for prompt peak fit)
            prompt plot check, 

        to activate: 

            fit, save and load setup buttons

        '''

        def load_setup(b):
            """
            when user presses this setup tab widget: mugui
            loads mujpy_setup.pkl with saved attributes
            and replaces them in setup tab Text widgets
            """
            import dill as pickle
            import os
                            
            path = os.path.join(self.__startuppath__,'mujpy_setup.pkl')
            # print('loading {}, presently in {}'.format(path,os.getcwd()))
            try:
                with open(path,'rb') as f:
                    mujpy_setup = pickle.load(f) 
            except:
                self.console('File {} not found'.format(path))
            # _paths_content = [ self.paths[k].value for k in range(3) ] # should be 3 ('data','tlag','analysis')
            # _filespecs_content = [ self.filespecs[k].value for k in range(2) ] # should be 2 ('fileprefix','extension')
            # _prepostpk = [self.prepostpk[k].value for k in range(2)] # 'pre-prompt bin','post-prompt bin' len(bkg_content)
            # _nt0 = self.nt0 # numpy array
            # _dt0 = self.dt0 # numpy array
            try:
                for k in range(3):  # len(paths_contents)
                    self.paths[k].value =  mujpy_setup['_paths_content'][k] # should be 3 ('data','tlag','analysis')
                for k in range(2):  # len(filespecs.content)
                    self.filespecs[k].value = mujpy_setup['_filespecs_content'][k] # should be 2 ('fileprefix','extension')
                for k in range(2):  # len(bkg_content)
                    self.prepostpk[k].value = mujpy_setup['_prepostpk'][k] # 'pre-prompt bin','post-prompt bin' 
                self.nt0 = mujpy_setup['self.nt0'] # bin of peak, nd.array of shape run.get_numberHisto_int()
                self.dt0 = mujpy_setup['self.dt0'] # fraction of bin, nd.array of shape run.get_numberHisto_int()
                self.lastbin = mujpy_setup['self.lastbin'] # fraction of bin, nd.array of shape run.get_numberHisto_int()
                self.nt0_run = mujpy_setup['self.nt0_run'] # dictionary to identify runs belonging to the same setup
                return 0
            except Exception as e:
                self.console('Error in load_setup: {}'.format(e)) 
                return -1   
 

        def on_introspect(b):
            '''
            print the unclean list of class attributes
            '''
            self.introspect() 

        def on_paths_changed(change):
            '''
            when user changes this setup tab widget mugui
            checks that paths exist, in case it creates log path
            '''
            import os

            path = change['owner'].description # description is paths[k] for k in range(len(paths)) () 
            k = paths_content.index(path) # paths_content.index(path) is 0,1,2 for paths_content = 'data','tlog','analysis'
            if k==1: # tlog
                self.newtlogdir = True
            directory = self.paths[k].value # self.paths[k] = handles of the corresponding Text
            if not os.path.isdir(directory):
                if k==2: # analysis, if it does not exist mkdir
                    # eventualmente togli ultimo os.path.sep = '/' in directory
                    dire=directory
                    if dire.rindex(os.path.sep)==len(dire):
                        dire=dire[:-1]
                    # splitta all'ultimo os.path.sep = '/'
                    prepath=dire[:dire.rindex(os.path.sep)+1]
                    # controlla che prepath esista
                    # print('prepath for try = {}'.format(prepath))
                    try:
                        os.stat(prepath)
                        os.mkdir(dire+os.path.sep)
                        self.console('Analysis path {} created\n'.format(directory)) 
                    except:
                        self.paths[k].value = os.path.curdir
                        self.console('Analysis path {} does not exist and cannot be created\n'.format(directory))
                else:
                    self.paths[k].value = os.path.curdir
                    self.console('Path {} does not exist, reset to .\n'.format(directory)) 
 
        def on_prompt_fit_click(b):
            '''
            when user presses this setup tab widget mugui
            executes prompt fits
            '''
            promptfit(mplot=self.plot_check.value) # mprint we leave always False

        def promptfit(mplot = False, mprint = False):
            '''
            launches t0 prompts fit::

                fits peak positions 
                prints migrad results
                plots prompts and their fit (if plot checked)
                stores bins for background and t0

            refactored for run addition and
            suite of runs

            WARNING: this module is for PSI only        
            '''
            import numpy as np
            from iminuit import Minuit, describe
            import matplotlib.pyplot as P
            from mujpy.mucomponents.muprompt import muprompt

            font = {'family' : 'Ubuntu','size'   : 8}
            P.rc('font', **font)
            dpi = 100.

            if not self._the_runs_:
                self.console('No run loaded yet! Load one first (select suite tab).')
            else:
                ###################################################
                # fit a peak with different left and right plateaus
                ###################################################

                #############################
                # guess prompt peak positions
                ############################# 
                npeaks = []
                for counter in range(self._the_runs_[0][0].get_numberHisto_int()):
                    histo = np.empty(self._the_runs_[0][0].get_histo_array_int(counter).shape)
                    for k in range(len(self._the_runs_[0])): # may add runs
                        histo += self._the_runs_[0][k].get_histo_array_int(counter)
                    npeaks.append(np.where(histo==histo.max())[0][0])
                npeaks = np.array(npeaks)

                ###############
                # right plateau
                ###############
                nbin =  max(npeaks) + self.second_plateau # this sets a counter dependent second plateau bin interval
                x = np.arange(0,nbin,dtype=int) # nbin bins from 0 to nbin-1
                self.lastbin, np3s = npeaks.min() - self.prepostpk[0].value, npeaks.max() + self.prepostpk[1].value # final bin of first and 

                if mplot:

                    ##############################
                    #  set or recover figure, axes 
                    ##############################
                    if self.fig_counters:
                        self.fig_counters.clf()
                        self.fig_counters,self.ax_counters = P.subplots(2,3,figsize=(7.5,5),num=self.fig_counters.number) 
                    else:
                        self.fig_counters,self.ax_counters = P.subplots(2,3,figsize=(7.5,5),dpi=dpi)
                        self.fig_counters.canvas.set_window_title('Prompts fit')
                    screen_x, screen_y = P.get_current_fig_manager().window.wm_maxsize() # screen size in pixels
                    y_maxinch = float(screen_y)/dpi #  maximum y size in inches

                    prompt_fit_text = [None]*self._the_runs_[0][0].get_numberHisto_int()    

                    for counter in range(self._the_runs_[0][0].get_numberHisto_int(),sum(self.ax_counters.shape)):
                        self.ax_counters[divmod(counter,3)].cla()
                        self.ax_counters[divmod(counter,3)].axis('off')

                x0 = np.zeros(self._the_runs_[0][0].get_numberHisto_int()) # for center of peaks
                for counter in range(self._the_runs_[0][0].get_numberHisto_int()):
                    # prepare for muprompt fit
                    histo = np.empty(self._the_runs_[0][0].get_histo_array_int(counter).shape)
                    for k in range(len(self._the_runs_[0])): # may add runs
                        histo += self._the_runs_[0][k].get_histo_array_int(counter)
                    p = [ self.peakheight, float(npeaks[counter]), self.peakwidth, 
                          np.mean(histo[self.firstbin:self.lastbin]), 
                          np.mean(histo[np3s:nbin])]
                    y = histo[:nbin]
                    ##############
                    # guess values
                    ##############
                    pars = dict(a=p[0],error_a=p[0]/100,x0=p[1]+0.1,error_x0=p[1]/100,
                         dx=1.1,error_dx=0.01,    
                         ak1=p[3],error_ak1=p[3]/100,ak2=p[4],error_ak2=p[4]/100)
                    level = 1 if mprint else 0
                    mm = muprompt()
                    mm._init_(x,y)
                    m = Minuit(mm,pedantic=False,print_level=level,**pars)
                    m.migrad()
                    A,X0,Dx,Ak1,Ak2 = m.args
                    x0[counter] = X0 # store float peak bin position (fractional)     
                    if mplot:    

                        n1 = npeaks[counter]-50
                        n2 = npeaks[counter]+50
                        x3 = np.arange(n1,n2,1./10.)
                        # with self.t0plot_container:
#                       if self.first_t0plot:
                        self.ax_counters[divmod(counter,3)].cla()
                        self.ax_counters[divmod(counter,3)].plot(x[n1:n2],y[n1:n2],'.')
                        self.ax_counters[divmod(counter,3)].plot(x3,mm.f(x3,A,X0,Dx,Ak1,Ak2))
                        x_text,y_text = npeaks[counter]+10,0.8*max(y)
                        prompt_fit_text[counter] = self.ax_counters[divmod(counter,3)].text(x_text,y_text,'Det #{}\nt0={}bin\n$\delta$t0={:.2f}'.format
                                                              (counter+1,x0.round().astype(int)[counter],x0[counter]-x0.round().astype(int)[counter]))
            if mplot:
                P.draw()
 
               ##################################################################################################
                # Simple cases: 
                # 1) Assume the prompt is entirely in bin nt0. (python convention, the bin index is 0,...,n,... 
                # The content of bin nt0 will be the t=0 value for this case and dt0 = 0.
                # The center of bin nt0 will correspond to time t = 0, time = (n-nt0 + mufit.offset + mufit.dt0)*mufit.binWidth_ns/1000.
                # 2) Assume the prompt is equally distributed between n and n+1. Then nt0 = n and dt0 = 0.5, the same formula applies
                # 3) Assume the prompt is 0.45 in n and 0.55 in n+1. Then nt0 = n+1 and dt0 = -0.45, the same formula applies.
                ##################################################################################################

            # these three are the sets of parameters used by other methods
            self.nt0 = x0.round().astype(int) # bin of peak, nd.array of shape run.get_numberHisto_int() 
            self.dt0 = x0-self.nt0 # fraction of bin, nd.array of shape run.get_numberHisto_int() 
            self.lastbin = self.nt0.min() - self.prepostpk[0].value # nd.array of shape run.get_numberHisto_int() 
            
            self.nt0_run = self.create_rundict()
            nt0.children[0].children[0].value = ' '.join(map(str,self.nt0.astype(int)))
            dt0.children[0].children[0].value = ' '.join(map('{:.2f}'.format,self.dt0))
                                                   # refresh, they may be slightly adjusted by the fit
            # self.t0plot_results.clear_output()

           #with self.t0plot_results:
            #    print('\n\n\n\nRun: {}'.format(self._the_runs_[0][0].get_runNumber_int()))
            #    print(' Bin nt0')
            #    for counter in range(self._the_runs_[0][0].get_numberHisto_int()):
            #        print('#{}: {}'.format(counter,self.nt0[counter]))
            #    print('\n\n dt0 (bins)')
            #    for counter in range(self._the_runs_[0][0].get_numberHisto_int()):
            #        print('#{}: {:.2f}'.format(counter,self.dt0[counter]))
            ##################################################################################################

        def save_log(b):
            """
            when user presses this setup tab button mugui
            saves ascii file .log with run list in data directory
            """
            import os
            from glob import glob
            from mujpy.musr2py.musr2py import musr2py as muload
            from mujpy.aux.aux import value_error

            run_files = sorted(glob(os.path.join(self.paths[0].value, '*.bin')))
            run = muload()
            run.read(run_files[0])
            filename=run.get_sample()+'.log'
            nastychar=list(' #%&{}\<>*?/$!'+"'"+'"'+'`'+':@')
            for char in nastychar:
                filename = "_".join(filename.split(char))
            path_file = os.path.join(self.paths[0].value, filename) # set to [2] for analysis
            with open (path_file,'w') as f:
                        #7082  250.0  250.0(1)   3   4.8  23:40:52 17-DEC-12  PSI8KMnFeF  Powder   PSI 8 K2.5Mn2.5Fe2.5F15, TF cal 30G, Veto ON, SR ON
                f.write("Run\tT_nom/T_meas(K)\t\tB(mT)\tMev.\tStart Time & Date\tSample\t\tOrient.\tComments\n\n")
                for run_file in run_files:
                    run.read(run_file)
                    TdT = value_error(run.get_temperatures_vector()[self.thermo],
                                      run.get_devTemperatures_vector()[self.thermo])
                    tsum = 0
                    for counter in range(run.get_numberHisto_int()):
                        histo = run.get_histo_array_int(counter).sum()
                        tsum += histo
                    BmT = float(run.get_field().strip()[:-1])/10. # convert to mT, avoid last chars 'G '
                    Mev = float(tsum)/1.e6
                            #Run T  TdT  BmT     Mev    Date   sam or com
                    f.write('{}\t{}/{}\t{:.1f}\t{:.1f}\t{}\t{}\t{}\t{}\n'.format(run.get_runNumber_int(),
                                      run.get_temp(), TdT, BmT, Mev, run.get_timeStart_vector(),
                                      run.get_sample().strip(), run.get_orient().strip(), run.get_comment().strip() ))
            self.console('Saved logbook {}'.format(path_file)) 

        def save_setup(b):
            """
            when user presses this setup tab button mugui
            saves mujpy_setup.pkl with setup tab values
            """
            import dill as pickle
            import os

            path = os.path.join(self.__startuppath__, 'mujpy_setup.pkl')
            # create dictionary setup_dict to be pickled 
            _paths_content = [ self.paths[k].value for k in range(3) ] # should be 3 ('data','tlag','analysis')
            _filespecs_content = [ self.filespecs[k].value for k in range(2) ] # should be 2 ('fileprefix','extension')
            _prepostpk = [self.prepostpk[k].value for k in range(2)] # 'pre-prompt bin','post-prompt bin' len(bkg_content)
            names = ['_paths_content','_filespecs_content',
                      '_prepostpk','self.nt0','self.dt0','self.lastbin','self.nt0_run'] # keys
            setup_dict = {}
            for k,key in enumerate(names):
               setup_dict[names[k]] = eval(key) # key:value
            with open(path,'wb') as f:
                pickle.dump(setup_dict, f) # according to __getstate__()
            self.console('Saved {}'.format(os.path.join(self.__startuppath__,'mujpy_setup.pkl'))) 

        from ipywidgets import HBox, Layout, VBox, Text, Textarea, IntText, Checkbox, Button, Output, Accordion   
        from numpy import array 

        # first tab: setup for things that have to be set initially (paths, t0, etc.)
        # the tab is self.mainwindow.children[0], a VBox 
        # containing a setup_box of three HBoxes: path, and t0plot 
        # path is made of a firstcolumn, paths, and a secondcolumns, filespecs, children of setup_box[0]
        # agt0 is made of three 
        setup_contents = ['path','promptfit','nt0_dt0','t0plot'] # needs two VBoxes
        setup_hbox = [HBox(description=name,layout=Layout(border='solid',)) for name in setup_contents]
        self.mainwindow.children[0].children = setup_hbox # first tab (setup)
        # first path
        paths_content = ['data','tlog','analysis'] # needs a VBox with three Text blocks
        paths_box = VBox(description='paths',layout=Layout(width='60%'))
        self.paths = [Text(description=paths_content[k],layout=Layout(width='90%'),continuous_update=False) for k in range(len(paths_content))]
        # self.paths[k].value='.'+os.path.sep # initial value, 
        paths_box.children = self.paths
        for k in range(len(paths_content)):
            self.paths[k].observe(on_paths_changed,'value')
        filespecs_content = ['fileprefix','extension']
        filespecs_box = VBox(description='filespecs',layout=Layout(width='40%'))
        self.filespecs = [Text(description=filespecs_content[k],layout=Layout(width='90%'),continuous_update=False) 
                          for k in range(len(filespecs_content))]
        filespecs_box.children = self.filespecs
        #        for k in range(len(filespecs)):  # not needed, only check that data and tlog exixt
        #            self.filespecs_list[k].observe(on_filespecs_changed,'value')
        # paths finished
        # now agt0
        self.prepostpk = [IntText(description='prepeak',value = 7, layout=Layout(width='20%'),
                              continuous_update=False), 
                          IntText(description='postpeak',value = 7, layout=Layout(width='20%'),    
                              continuous_update=False)]
        self.prepostpk[0].style.description_width='60%'
        self.prepostpk[1].style.description_width='60%'
        self.plot_check = Checkbox(description='prompt plot',value=True,layout=Layout(width='15%'))
        self.plot_check.style.description_width='10%'
        fit_button = Button(description='prompt fit',layout=Layout(width='15%'))
        fit_button.on_click(on_prompt_fit_click)
        fit_button.style.button_color = self.button_color
        save_button = Button(description='save setup',layout=Layout(width='15%'))
        save_button.style.button_color = self.button_color
        load_button = Button(description='load setup',layout=Layout(width='15%'))
        load_button.style.button_color = self.button_color
        prompt_fit = [self.prepostpk[0], self.prepostpk[1], self.plot_check, fit_button ,save_button, load_button] 
        # fit bin range is [self.binrange[0].value:self.binrange[1].value]
        save_button.on_click(save_setup)
        load_button.on_click(load_setup)
        nt0 = Accordion(font_size=10,children=[VBox(children=[Text(description='t0 [bins]',layout={'width':'99%'})])],layout={'width':'35%'})#,layout={'height':'22px'})
        dt0 = Accordion(font_size=10,children=[VBox(children=[Text(description='dt0 [bins]',layout={'width':'99%'})])],layout={'width':'35%'})#,layout={'height':'22px'})
        nt0.children[0].children[0].style.description_width='20%'
        nt0.set_title(0,'t0 bins')
        nt0.selected_index =  None
        dt0.children[0].children[0].style.description_width='20%'
        dt0.set_title(0,'t0 remainders')
        dt0.selected_index =  None
        # self.tots_all = Textarea(description='All',layout={'width':'100%','height':'200px',
        #                                         'overflow_y':'auto','overflow_x':'auto'},disabled=True)
        # self.tots_all.style.description_width='15%'
        # self.tots_group = Textarea(description='Group',layout={'width':'100%','height':'200px',
        #                                         'overflow_y':'auto','overflow_x':'auto'},disabled=True)
        # self.tots_group.style.description_width='30%'
        # tots = Accordion(font_size=10,children=[HBox(children=[self.tots_all, self.tots_group])],layout={'width':'35%'})
        # tots.set_title(0,'Total counts')
        # tots.selected_index =  None

        introspect_button = Button(description='Introspect',layout=Layout(width='15%'))
        introspect_button.on_click(on_introspect)
        introspect_button.style.button_color = self.button_color
        log_button = Button(description='Data log',layout=Layout(width='15%'))
        log_button.on_click(save_log)
        log_button.style.button_color = self.button_color 

        #self.t0plot_container = Output(layout=Layout(width='85%'))     
        self.t0plot_results = Output(layout=Layout(width='15%')) 
        setup_hbox[0].children = [paths_box, filespecs_box]
        setup_hbox[1].children = prompt_fit
        setup_hbox[2].children = [nt0, dt0,introspect_button,log_button]
        #setup_hbox[3].children = [self.t0plot_container,self.t0plot_results]
        self.nt0,self.dt0 = array([0.]),array([0.])
        load_setup([])
        nt0.children[0].children[0].value = ' '.join(map(str,self.nt0.astype(int)))
        dt0.children[0].children[0].value = ' '.join(map('{:.2f}'.format,self.dt0))
        if not self.nt0_run:
            self.console('WARNING: you must fix t0 = 0, please do a prompt fit from the setup tab')

######################
# GET_TOTALS
######################
    def get_totals(self):
        '''
        calculates the grand totals and group totals 
        after a single run 
        or a run suite are read

        '''
        import numpy as np
        # called only by self.suite after having loaded a run or a run suite

        ###################
        # grouping set 
        # initialize totals
        ###################
        gr = set(np.concatenate((self.grouping['forward'],self.grouping['backward'])))
        ts,gs =  [],[]

        if self.offset:  # True if self.offset is already created by self.fit()
            offset_bin = self.offset.value # self.offset.value is int
        else:     # should be False if self.offset = [], as set in self.__init__()
            offset_bin = self.offset0  # temporary parking

        for k,runs in enumerate(self._the_runs_):
            tsum, gsum = 0, 0
            for j,run in enumerate(runs): # add values for runs to add
                for counter in range(run.get_numberHisto_int()):
                    n1 = offset_bin+self.nt0[counter] 
                    # if self.nt0 not yet read it is False and totals include prompt
                    histo = run.get_histo_array_int(counter)[n1:].sum()
                    tsum += histo
                    if counter in gr:
                        gsum += histo
            ts.append(tsum)
            gs.append(gsum)
            # print('In get totals inside loop,k {}, runs {}'.format(k,runs))

        #######################
        # strings containing 
        # individual run totals
        #######################
        # self.tots_all.value = '\n'.join(map(str,np.array(ts)))
        # self.tots_group.value = '       '.join(map(str,np.array(gs)))

        # print('In get totals outside loop, ts {},gs {}'.format(ts,gs))
        #####################
        # display values for self._the_runs_[0][0] 
        self.totalcounts.value = str(ts[0])
        self.groupcounts.value = str(gs[0])
        # self.console('Updated Group Total for group including counters {}'.format(gr)) # debug 
        self.nsbin.value = '{:.3}'.format(self._the_runs_[0][0].get_binWidth_ns())
        self.maxbin.value = str(self.histoLength)


##########################
# SUITE
##########################
    def suite(self):
        '''
        suite tab of mugui
        used to select: 

            run (single/suite)
            load next previous, add next previous
             
        and to print: 

            run number, title, 
            total counts, group counts, ns/bin
            comment, start stop date, next run, last add                           

        '''

        def run_headers(k):
            '''
            Stores and displays
            title, comments and histoLength only for master run
            Saves T, dT and returns 0
            '''
            import numpy as np
            from mujpy.aux.aux import get_title, value_error
            if k==0:
                try:
                    dummy = self.nt0.sum() # fails if self.nt0 does not exist yet
                except: # if self.nt0 does not exist, guess from the first in self._the_runs_
                    self.nt0 = np.zeros(self._the_runs_[0][0].get_numberHisto_int(),dtype=int)
                    self.dt0 = np.zeros(self._the_runs_[0][0].get_numberHisto_int(),dtype=float)
                    for j in range(self._the_runs_[0][0].get_numberHisto_int()):
                        self.nt0[j] = np.where(self._the_runs_[0][0].get_histo_array_int(j)==
                                               self._the_runs_[0][0].get_histo_array_int(j).max())[0][0]

                # self.nt0 exists
                self.title.value = get_title(self._the_runs_[0][0])                
                self.comment_handles[0].value = self._the_runs_[0][0].get_comment() 
                self.comment_handles[1].value = self._the_runs_[0][0].get_timeStart_vector() 
                self.comment_handles[2].value = self._the_runs_[0][0].get_timeStop_vector()
                self._the_runs_display.value = str(self.load_handle[0].value)
                # but if it is not compatible with present first run issue warning 
                if len(self.nt0)!=self._the_runs_[0][0].get_numberHisto_int(): # reset nt0,dt0
                    self.nt0 = np.zeros(self._the_runs_[0][0].get_numberHisto_int(),dtype=int)
                    self.dt0 = np.zeros(self._the_runs_[0][0].get_numberHisto_int(),dtype=float)
                    for j in range(self._the_runs_[0][0].get_numberHisto_int()):
                        self.nt0[j] = np.where(self._the_runs_[0][0].get_histo_array_int(j)==
                                               self._the_runs_[0][0].get_histo_array_int(j).max())[0][0]
                    self.console('WARNING! Run {} mismatch in number of counters, rerun prompt fit'.format(self._the_runs_[0][0].get_runNumber_int())) 

                # store max available bins on all histos
                self.histoLength = self._the_runs_[0][0].get_histoLength_bin() - self.nt0.max() - self.offset.value 
                self.counternumber.value = '    {} counters per run'.format(self._the_runs_[0][0].get_numberHisto_int())
                self.plot_range0 = '0,{},100'.format(self.histoLength)
                self.multiplot_range.value = self.plot_range0
                if self.plot_range.value == '':
                    self.plot_range.value = self.plot_range0
                if self.fit_range.value == '':
                    self.fit_range.value = self.plot_range0                
                npk = float(self.nt0.sum())/float(self.nt0.shape[0])
                self.bin_range0 = '{},{}'.format(int(0.9*npk),int(1.1*npk))
                self.counterplot_range.value = self.bin_range0
            else:  # k > 0 
                self._single_ = False
                ok = [self._the_runs_[k][0].get_numberHisto_int() == self._the_runs_[0][0].get_numberHisto_int(),
                      self._the_runs_[k][0].get_binWidth_ns() == self._the_runs_[0][0].get_binWidth_ns()]
                if not all(ok): 
                    self._the_runs_=[self._the_runs_[0]] # leave just the first one
                    # self.load_handle[1].value=''  # just loaded a single run, incompatible with suite
                    self.console('\nFile {} has wrong histoNumber or binWidth'.format(path_and_filename))
                    return -1  # this leaves the first run of the suite
            TdT = value_error(*t_value_error(k))
            self.tlog_accordion.children[0].value += '{}: '.format(self._the_runs_[k][0].get_runNumber_int())+TdT+' K\n'
            # print('3-run_headers')

            return 0

        def check_next():
            '''
            Checks if next run file exists
            '''
            import os
            from mujpy.aux.aux import muzeropad

            runstr = str(self.nrun[0] +1)
            if len(runstr)>4:
                self.console('Too long run number!')
                next_label.value = ''
            else:
                filename = ''
                filename = filename.join([self.filespecs[0].value,
                                  muzeropad(runstr),
                                  '.',self.filespecs[1].value]) 
                            # data path + filespec + padded run rumber + extension)
                next_label.value = runstr if os.path.exists(os.path.join(self.paths[0].value,filename)) else ''
                        
        def check_runs(k):
            '''
            Checks nt0, etc.
            Returns -1 with warnings  
            for severe incompatibility
            Otherwise calls run_headers to store and display 
            title, comments, T, dT, histoLength [,self._single]
            '''
            from copy import deepcopy
            from dateutil.parser import parse as dparse
            import datetime
 
            if self.nt0_run: # either freshly produced or loaded from load_setup
                nt0_experiment = deepcopy(self.nt0_run) # needed to preserve the original from the pops
                nt0_experiment.pop('nrun')
                nt0_days = dparse(nt0_experiment.pop('date'))
                try:
                    this_experiment = self.create_rundict(k) # disposable, no deepcopy, for len(runadd)>1 check they are all compatible
                    # print('check - {}'.format(self._the_runs_[k][0].get_runNumber_int()))
                    rn = this_experiment.pop('nrun') # if there was an error with files to add in create_rundict this will fail
                except:
                    self.console('\nRun {} not added. Non existent or incompatible'.format(this_experiment('errmsg')))
                    return -1 # this leaves the previous loaded runs n the suite 
                
                this_date = this_experiment.pop('date') # no errors with add, pop date then
                dday = abs((dparse(this_date)-nt0_days).total_seconds())
                if nt0_experiment != this_experiment or abs(dday) > datetime.timedelta(7,0).total_seconds(): # runs must have same binwidth etc. and must be within a week
                    self.console('Warning: mismatch in histo length/time bin/instrument/date\nConsider refitting prompt peaks (in setup)')
            # print('2-check_runs, {} loaded '.format(rn))            
            return run_headers(k)

        def add_runs(k,runs):
            '''
            Tries to load one or more runs to be added together
            by means of murs2py. 
            runs is a list of strings containing integer run numbers provided by aux.derun
            Returns -1 and quits if musr2py complains
            If not invokes check_runs an returns its code   
            '''
            import os
            from mujpy.musr2py.musr2py import musr2py as muload
            from mujpy.aux.aux import muzeropad
            read_ok = 0
            runadd = []
            options = self.choose_tlogrun.options.copy()
            for j,run in enumerate(runs): # run is a string containing a single run number
                if len(run)>4:
                    self.console('Too long run number!')
                    read_ok = 99
                    break
                else:
                    filename = ''
                    filename = filename.join([self.filespecs[0].value,
                                          muzeropad(str(run)),
                                          '.',self.filespecs[1].value]) 
                    path_and_filename = os.path.join(self.paths[0].value,filename)
                                    # data path + filespec + padded run rumber + extension)
                    if os.path.exists(os.path.join(self.paths[0].value,filename)):
                        runadd.append(muload())  # this adds to the list in j-th position a new instance of muload()
                        read_ok += runadd[j].read(path_and_filename) # THE RUN DATA FILE IS LOADED HERE
                        if read_ok==0:
                        # print('tlog dropdown position {} run {}'.format(str(j),str(run)))
                            options.update({str(run):str(run)}) # adds this run to the tlog display dropdown, on_loads_changed checks that tlog exists before value selection
            if read_ok==0: # no error condition, set by musr2py.cpp or
                self.choose_tlogrun.options = options
                # ('self.choose_tlogrun.options = {}'.format(options))
                self._the_runs_.append(runadd) # 
                self.nrun.append(runadd[0].get_runNumber_int())
            else:
                self.console('\nFile {} not read. Check paths, filespecs and run rumber on setup tab'.format(path_and_filename))
                return -1 # this leaves the previous loaded runs n the suite 
            return check_runs(k)

        def on_load_nxt(b):
            '''
            load next run (if it exists)
            '''
            if self._single_:
                # print('self.nrun[0] = {}'.format(self.nrun[0]))
                self.load_handle[0].value=str(self.nrun[0]+1)
                # print('self.nrun[0] = {}'.format(self.nrun[0]))
            else:
                self.console('Cannot load next run (multiple runs loaded)')
                return -1 # this leaves the multiple runs of the suite 
                
        def on_load_prv(b):
            '''
            load previous run (if it exists)
            '''
            if self._single_:
                self.load_handle[0].value=str(self.nrun[0]-1)
            else:
                self.console('Cannot load previous run (multiple runs loaded)')
                return -1 # this leaves the multple runs of the suite 
                
        def on_add_nxt(b):
            '''
            add next run (if it exists)
            '''
            if self._single_:
                load_single(self.nrun[0]+1)
                self.get_totals()
            else:
                self.console('Cannot add next run (multiple runs loaded)')
                return -1 # this leaves the multiple loaded runs of the suite 
                
        def on_add_prv(b):
            '''
            add previous run (if it exists)
            '''
            if self._single_:
                load_single(self.nrun[0]-1)
                self.get_totals()
            else:
                self.console('Cannot add previous run (multiple runs loaded)')
                return -1 # this leaves the multiple loaded runs of the suite 
                

        def on_loads_changed(change):
            '''
            observe response of suite tab widgets:
            load a run via musrpy 
            single run and run suite unified in a list
               clears run suite
               loads run using derun parsing of a string
                  csv, n:m for range of runs 
                  [implement n+m+... for run addition]
               sets _single_ to True if single
            plan: derun must recognize '+', e.g.
                  '2277:2280,2281+2282,2283:2284'
                  and produce 
                  run = [['2277'],['2278'],['2279'],['2280'],['2281','2282'],['2283'],['2284']]
                  Then the loop must subloop on len(run) to recreate the same list structure in self._the_runs_
                  and all occurencies of self._the_runs_ must test to add data from len(self._the_runs_[k])>1
                    check also asymmetry, create_rundict, write_csv, get_totals, promptfit, on_multiplot
            '''
            from mujpy.aux.aux import derun, tlog_exists
 
            # rm: run_or_runs = change['owner'].description # description is either 'Single run' or 'Run  suite'
            if self.load_handle[0].value=='': # either an accitental empty text return, or reset due to derun error
                return
            self._single_ = True
            self._the_runs_ = []  # it will be a list of muload() runs
            self.nrun = [] # it will contain run numbers (the first in case of run add)
            self.tlog_accordion.children[0].value=''
 
            #######################
            # decode the run string
            #######################           
            runs, errormessage = derun(self.load_handle[0].value) # runs is a list of lists of run numbers (string)
            if errormessage is not None: # derun error
                self.console('Run syntax error: {}. You typed: {}'.format(errormessage,self.load_handle[0].value))
                self.load_handle[0].value=''
                return

            ##################################
            # load a single run or a run suite
            ##################################
            read_ok = 0
            for k,runs_add in enumerate(runs):#  rs can be a list of run numbers (string) to add
                read_ok += add_runs(k,runs_add)                
                # print('on_loads_change, inside loop, runs {}'.format(self._the_runs_))
            if read_ok == 0:
                self.choose_nrun.options  = [str(n) for n in self.nrun]
                self.choose_nrun.value = str(self.nrun[0])
                options = self.choose_tlogrun.options.copy()
                runs = self.choose_tlogrun.options.keys() 
                # scheme for popping items from list without altering index count
                kk = 0
                for run in runs: # use original to iterate
                    if not tlog_exists(self.paths[1].value,run): # loading tlogs is optional
                        options.pop(run) # pop runs that do not have a tlog file
                self.choose_tlogrun.options = options
                try: 
                    self.choose_tlogrun.value = str((sorted(list(options.keys())))[0])
                except:
                    if self.newtlogdir:
                        self.console("No tlogs")
                        self.newtlogdir = False
                self.get_totals() # sets totalcounts, groupcounts and nsbin
                
            if self._single_:
                check_next()
            if not self.nt0_run:
                self.console('WARNING: you must fix t0 = 0, please do a prompt fit from the setup tab')

        def t_value_error(k):
            '''
            calculates T and eT values also for runs to be added
            silliy, but it works also for single run 
            '''
            from numpy import sqrt

            m = len(self._the_runs_[k])
            weight = [float(sum(self._the_runs_[k][j].get_histo_array_int(2))) for j in range(m)]
            weight = [w/sum(weight) for k,w in enumerate(weight)]
            t_value = sum([self._the_runs_[k][j].get_temperatures_vector()[self.thermo]*weight[j] for j in range(m)])
            t_error = sqrt(sum([(self._the_runs_[k][j].get_devTemperatures_vector()[self.thermo]*weight[j])**2 for j in range(m)]))
            return t_value, t_error
                                       
        from ipywidgets import HBox, Layout, Text, Button

        # second tab: select run or suite of runs (for sequential or global fits)
        # the tab is self.mainwindow.children[1], a VBox 
        # containing three HBoxes, loads_box, comment_box, speedloads_box
        # path is made of a firstcolumn, paths, and a secondcolumns, filespecs, children of setup_box[0]
        # rm: loads = ['Single run','Run suite'] 
        speedloads = ['Next run' 'Load next', 'Load previous', 'Add next', 'Add previous', 'Last added']
        load_box = HBox(description='loads',layout=Layout(width='100%'))
        speedloads_box = HBox(description='speedloads',layout=Layout(width='100%'))
        self.load_handle = [Text(description='Run[run suite]:       \nsingle run\ne.g 431\nor run suites\ne.g. 431, 435:439, 443+444',
                              layout=Layout(width='100%'),continuous_update=False)]
        self.load_handle[0].style.description_width='11%'
        self.load_handle[0].observe(on_loads_changed,'value')
        # the following doesn't work yet
        Ln_button = Button(description='Load nxt')
        Ln_button.on_click(on_load_nxt)
        Ln_button.style.button_color = self.button_color
        Lp_button = Button(description='Load prv')
        Lp_button.on_click(on_load_prv)
        Lp_button.style.button_color = self.button_color
        An_button = Button(description='Add nxt')
        An_button.on_click(on_add_nxt)
        An_button.style.button_color = self.button_color
        Ap_button = Button(description='Add prv')
        Ap_button.on_click(on_add_prv)
        Ap_button.style.button_color = self.button_color
        next_label = Text(description='Next run',disabled=True)
        self.speedloads_handles = [next_label,
                                   Ln_button, Lp_button, An_button, Ap_button, 
                                   Text(description='Last add',disabled=True)]
        load_box.children = self.load_handle
        speedloads_box.children = self.speedloads_handles

        self.mainwindow.children[1].children = [load_box, speedloads_box] # second tab (suite)

    def introspect(self):
        '''
        print updated attributes of the class mugui 
        after each fit in file "mugui.attributes.txt" 
        in self.__startuppath__
        '''
        import os
        import pprint
        from ipywidgets import VBox, HBox, Image, Text, Textarea, Layout, Button, IntText, Checkbox, Output, Accordion, Dropdown, FloatText, Tab
        # trick to avoid printing the large mujpy log image binary file 
        image = self.__dict__['gui'].children[0].children[0].value
        self.__dict__['gui'].children[0].children[0].value=b''
        with open(os.path.join(self.__startuppath__,"mugui.attributes.txt"),'w') as f:
            pprint.pprint('**************************************************',f)
            pprint.pprint('*               Mugui attribute list:            *',f)
            pprint.pprint('**************************************************',f)
            pprint.pprint(self.__dict__,f)
        self.__dict__['gui'].children[0].children[0].value=image

