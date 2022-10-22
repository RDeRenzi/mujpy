'''A musr analysis class
Notes

    #   NB mujpy is python3 only
    Starts from init, the rest is operated by gui ipywidgets
    output for linux; Popen _xterm_ and tail -f /tmp/mujpy_pipe
    then self.console(string) takes care of writing strings to pipe
    Jupyterlab provides a way of writing stderr to console
    description_tooltip in ipywidgets >7.3 should for all items with description
    as of 7.5 buttons still require "tooltip" instead of "description_tooltip"
    
    Gui tabs correspond to distinct gui methods with independes scopes and additional local methods
    gui attributes: 
         entities that must be shared between tabs 
         including variables passed to functions  outside mugui
    

    Centered around a minuit fit tab. The fit is performed only in asymmetry mode, 
    and its details are derived from the syntax of the conventional model name.
    Hitting the FIT button after the selection on the fit tab dashboard loads 
    asymmetry data and errors as asymm, asyme, by a method switchyard, in either 
                            single runs [tbin] 
                            suites of runs [run,tbin]
                            multi groups [group, tbin]
                            multi groups suites [run, group, tbin]           
    This corresponds to equivalent switchyards in the mucomponents class (mumodel), 
    in the int2min parameter value translations from dash to minuit 
    and min2int viceversa get rid of
    in the int2_int component method selection, 
    in the fitargs dictionary used in save_fit, load_fit methods 
    in the write_csv method 
                

'''

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
        # only constant used here, all others in mujpy.mucomponents
        self.TauMu_mus = 2.1969811 # number is from Particle Data Group 2017 (not in scipy.constants.physical_constants) 
        # print('__init__ now ...')
        self.initialize_mujpy()
        
    def initialize_mujpy(self):
        from mujpy import __file__ as MuJPyName
        import numpy as np
        import os
        from IPython.display import display
       
# some initializations

        self.offset0 = 7 # initial value
        self.offset = [] # this way the first run load calculates get_totals with self.offset0
        self.firstbin = 0
        # These three RELOCATED TO promptfit local variables
        # self.second_plateau = 100 
        # self.peakheight = 100000.
        # self.peakwidth = 1.   # rough guesses for default 
        self.histoLength = 7900 # initialize
        self.bin_range0 = '0,500' # initialize (for plots, counter inspection)
        self.nt0_run = []
        # self._global_ = False
        self.thermo = 1 # sample thermometer is 1 on gps (check or adapt to other instruments)

        self.binwidth_ns = [] # this way the first call to asymmetry(_the_runs_) initializes self.time
        self.grouping = {'forward':np.array([1]),'backward':np.array([0])} # normal dict
        self._the_runs_ = []  # if self._the_runs_: is False, to check whether the handle is created
        self.first_t0plot = True
        self.fitargs = [] # initialize, 
        # fitargs was iminuit 1.0 shuttle, dictionary containing names, values, errors
        # reprocuce it! still useful in fft of residues (print subtracted components)
        # write_csv to print names 
        # save_fit and load_fit

# mujpy paths
        self.__path__ = os.path.dirname(MuJPyName)
        self.__logopath__ = os.path.join(self.__path__,"logo")
        self.__startuppath__ = os.getcwd() # working directory, in which, in case, to find mujpy_setup.pkl 
# mujpy layout
        self.button_color = 'lightblue'
        self.button_color_off = 'lightgray'
        self.border_color = 'dodgerblue'
        self.newtlogdir = True
        # print('...finshed initializing')

        ##########################################
        # This actually produces the gui interface
        # with tabs: 
        # suite fit fft plots setup [output] about
        ##########################################
        
        # print('Got here: create gui')
        self.output() # this calls output(self) that defines self._output_
        # must be called first! 
        # either it launches a detached os terminal or it adds a tab
        
        self.gui()
        self.setup()

        #####################################
        # static figures 
        #####################################
        self.fig_fft = [] # initialize to false, it will become a pyplot.subplots instance
        self.fig_multiplot = []
        self.fig_counters = []

        self.fit()
        self.fft()
        self.plots()
        self.about()


        #self.fig_fit = [] # (initialize to false),  now not ititialized
        try:
            whereamI = get_ipython().__class__
            if not str(whereamI).find('erminal')+1:
                display(self.gui) # you are in a Jupyter notebook
            else:
                print(str(wheremI)) # you are in an ipython terminal
        except:
                print('Python test script') # other option?



##########################
# ABOUT
##########################
    def about(self):
        '''
        about tab:

        - a few infos (version and authors)

        '''
        from ipywidgets import Textarea, Layout, HTML, VBox

        _version = 'MuJPy          version '+'1.0' # increment while progressing
        _authors = '\n\n  Authors: Roberto De Renzi, Pietro Bonfà (*)'
        _blahblah = ('\n\n  A Python MuSR data analysis graphical interface.'+
                     '\n  Based on classes, designed for jupyter.'+
                     '\n  Released under the MIT licence')
        _html = ('\n  See docs in the <a href="https://mujpy.readthedocs.io/en/latest/Tutorial.html"target="_blank">Tutorial</a>')
        _pronounce = ('\n  Pronounce it as mug + pie')
        _additional_credits_ = ('\n ---------------------\n (*) dynamic Kubo-Toyabe algorithm by G. Allodi\n     MuSR_td_PSI by A. Amato and A.-R. Raselli \n     acme algorithm code from NMRglue, by Jonathan J. Helmus')
        _about_text = _version+_blahblah+_pronounce+_authors+_additional_credits_
        _about_area = Textarea(value=_about_text,
                                   description='Info on MuJPy',
                                   layout=Layout(width='98%',height='250px'),
                                   disabled=True)
        _about_html = HTML(value=_html,
                                   description='ReadTheDocs',
                                   layout=Layout(width='98%',height='250px'),
                                   disabled=True)
        # now collect the handles of the three horizontal frames to the main fit window (see tabs_contents for index)
        if self._output_==self._outputtab_:
            self.mainwindow.children[4].children = self._park_ # was parked here in self.output
            self.mainwindow.children[4].layout = Layout(border = '2px solid dodgerblue',width='100%')
            self.mainwindow.children[5].children = [VBox([_about_area,_about_html])] # add the list of widget handles as the third tab, fit
            self.mainwindow.children[5].layout = Layout(border = '2px solid dodgerblue',width='100%')
        else:
            self.mainwindow.children[4].children = [VBox([_about_area,_about_html])]  # add the list of widget handles as the third tab, fit
            self.mainwindow.children[4].layout = Layout(border = '2px solid dodgerblue',width='100%')
            
            
##########################
# ASYMMETRY
##########################
    def timebase(self):
        """
        * initializes self numberHisto histoLength binwidth_ns
        * returns time array
        * 1D numpy array
        """ 
        import numpy as np

        # no checks, consistency in binWidth and numberHisto etc are done with run loading
        self.numberHisto = self._the_runs_[0][0].get_numberHisto_int()
        self.histoLength = self._the_runs_[0][0].get_histoLength_bin() - self.nt0.max() - self.offset.value # max available bins on all histos
        self.binwidth_ns = self._the_runs_[0][0].get_binWidth_ns() 
   
    ##################################################################################################
    # Time definition: 
    # 1) Assume the prompt is entirely in bin self.nt0. (python convention, the bin index is 0,...,n,... 
    # The content of bin self.nt0 will be the t=0 value for this case and self.dt0 = 0.
    # The center of bin self.nt0 will correspond to time t = 0, 
    #          time = (n-self.nt0 + self.offset.value + self.dt0)*mufit.binWidth_ns/1000.
    # 2) Assume the prompt is equally distributed between n and n+1. 
    #    Then self.nt0 = n and self.dt0 = 0.5, the same formula applies
    # 3) Assume the prompt is 0.45 in n and 0.55 in n+1. 
    #    Then self.nt0 = n+1 and self.dt0 = -0.45, the same formula applies.
    ##################################################################################################
        
        if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu': # only for PSI, not for ISIS
            time = (np.arange(self.histoLength) + self.offset.value +
                         np.mean(self.dt0 [np.append(self.grouping['forward'],self.grouping['backward'])] ) 
                         )*self.binwidth_ns/1000. # in microseconds, 1D np.array
        else:  # ISIS
            time = (np.arange(self.histoLength) + self.offset.value +
                         self.dt0)*self.binwidth_ns/1000. # in microseconds, 1D np.array       
        return time

    def single_asymmetry(self,runs,grouping,alpha):
        """
        * returns asymmetry and error, 
        * all are 1D numpy arrays
        """
        import numpy as np

#       initialize to zero 
        yforw = np.zeros(self.time.shape[0]) # counts with background substraction
        ybackw = np.zeros(self.time.shape[0]) # counts with background substraction
        if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu': # only for PSI, not for ISIS
            cforw = np.zeros(self.time.shape[0]) # corrected counts for Poisson errors
            cbackw = np.zeros(self.time.shape[0]) # corrected counts for Poisson errors
            M = self.lastbin-self.firstbin+1 # number of pre-prompt bins for background average

        n1, n2 = self.nt0[0] + \
                 self.offset.value, self.nt0[0] + \
                                    self.offset.value + \
                                    self.histoLength 
                                 # ISIS: set in suite.run_headers

        for j, run in enumerate(runs): # This is for adding runs

            for counter in grouping['forward']: # first good bin, last good bin from data attay start
                histo = run.get_histo_array_int(counter) # counter data array in forw group
                if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu': # only for PSI, not for ISIS
                    n1, n2 = self.nt0[counter] + \
                             self.offset.value, self.nt0[counter] + \
                                                self.offset.value + \
                                                self.histoLength 
                                         # set in suite.run_headers as max common good data length
                    background = np.mean(histo[self.firstbin:self.lastbin])  # from prepromt estimate
                    yforw += histo[n1:n2]-background           #  background subtracted add to others
                    cforw += abs(histo[n1:n2]-background*(1-1/M))  #  add with correction, see error eval box
                else:
                    yforw += histo[n1:n2]
            for counter in grouping['backward']: # first good bin, last good bin from data attay start
                histo = run.get_histo_array_int(counter) # counter data array in back group
                if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu': # only for PSI, not for ISIS
                    n1, n2 = self.nt0[counter] + \
                             self.offset.value, self.nt0[counter] + \
                                                self.offset.value + \
                                                self.histoLength 
                                         # set in suite.run_headers as max common good data length
                    background = np.mean(histo[self.firstbin:self.lastbin])  # from prepromt estimate
                    ybackw += histo[n1:n2]-background   #  background subtracted add to others
                    cbackw += abs(histo[n1:n2]-background*(1-1/M)) #  add with correction, see error eval box
                else:
                    ybackw += histo[n1:n2]

        yplus = yforw + float(alpha)*ybackw # denominator
        x = np.exp(-self.time/self.TauMu_mus)        # pure muon decay
        enn0 = np.polyfit(x,yplus,1)            # fit of initial count rate plus average background vs. decay
        enn0 = enn0[0]                          # initial rate per ns    
        
            # A(t) = y  with  corrected counts Nfc(t) = Nf(t) - bf = yforw, Nbc(t) = Nb(t) - bb = ybackw
            # until 11-2020 error Poissonian in Nf, Nb, wrong since Njc  and bj  j=f,b are independent   
            # from  11-2020 correction, thanks to Giuseppe Allodi for pointing this out
            # and Andreas Suter, inspiring through musrfit bitbucker code
            # Nf(t) = cforw ,   Nb(t) = cbackw            # whereas  yield 
            
                        #########################################################
                        #                  error evaluation:                    #
                        # A = d [Nfc(t) - alpha Nbc(t)] = d a(t)     with       #
                        #   d = 1/[2N0e(-t/tau)], Nfc = Nf -bf, Nbc = Nb - bb   #
                        # eNfc = sqrt(|Nf-bf|), eNbc = sqrt(|Nb-bb|)            #
                        # ebf = sqrt(bf/M), ebb= sqrt(bb/M)  b = M sum_i^M bi   #
                        #   sum_i^M b_i = Mb  e_sum = sqrt(Mb) eb = 1/M e_sum   #
                        #                  error propagation:                   #
                        # |da/dNfc| = |da/dbf| = 1, |da/dNbc| = |da/dbb|= alpha #
                        # eA/d = sqrt[(|da/dNfc|eNfc)^2 + (|da/dbf|ebf)^2) +    #
                        #        sqrt[(|da/dNb|eNb)^2 + (|da/dbb|ebb)^2)        #
                        # eA = d sqrt[Nf - bf(1-1/M) + alpha^2(Nb - bb(1-1/M))] #
                        #########################################################
                        
        # the correction happens in cforw, cbackw above
        if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu': # PSI
                
            ey = np.sqrt(cforw + float(alpha)**2*cbackw) \
                        *np.exp(self.time/self.TauMu_mus)/enn0 #
        else: # ISIS
            ey = np.sqrt(yforw + float(alpha)**2*ybackw) \
                        *np.exp(self.time/self.TauMu_mus)/enn0 #
#         ey definition, in both cforw, cbackw, contains difference
#         of independent terms: histo - background*(1-1/M), could be <=0
#         Denominator of fcn sum_bin ((y_bin-y_th)/ey)**2 cannot be zero 
#         Set to 1 the minimum error (weights less very few points closer to ~ zero) 
        ey[np.where(ey==0)] = 1 # substitute 0 with 1 in ey 

        return (yforw-float(alpha)*ybackw)/enn0*np.exp(self.time/self.TauMu_mus), ey

    def multi_asymmetry(self,runs):
        """
        * returns asymmetry and error, 
        * all are 2D numpy arrays 
        """

        for k, grouping in enumerate(self.grouping):
            y, ey = single_asymmetry(runs,grouping,self.alpha[k])
            if ~k:
                asymm, asyme = [y], [ey]
            else:
                asymm.append(y)
                asyme.append(ey)
        return np.array(asymm), np.array(asyme)

    def multisuite_asymmetry(self,runs):
        """
        * returns asymmetry and error, 
        * all are 3D numpy arrays 
        """
        for j, runs in enumerate(self._the_runs_):
            y, ey = multi_asymmetry(runs)
            if ~j:
                asymm, asyme = [y], [ey]
            else:
                asymm.append(y)
                asyme.append(ey )
        return np.array(asymm), np.array(asyme)   
    
    def asymmetry_single(self):
        """
        * defines self.time, generates asymmetry end error without rebinning, 
        * all are 1D numpy arrays
        * self._the_runs_ is either one or a list (run add) of musr2py instances (psi/isis load routine)  
        * list is for adding runs
        """ 

        # calculate time
        self.time = self.timebase() # in microseconds, 1D np.array 
        # calculate asymmetry and error; self._the_runs_[0] is a run or a list of runs to be added 
        self.asymm, self.asyme = self.single_asymmetry(self._the_runs_[0],
                                                  self.grouping,
                                                  self.alpha.value)                    
        self.nrun = [self._the_runs_[0][0].get_runNumber_int()]
        ######################################################
        # self.nrun contains only the first run in case of run addition
        # used by save_fit (in file name), 
        #         write_csv (first item is run number)
        #         animate_fit (multiplot) 
        ######################################################
        # self.console('1st, last = {},{} bin for bckg subtraction'.format(self.firstbin,self.lastbin)) #dbg


    def asymmetry_suite(self):
        """
        * defines self.time, generates asymmetry end error without rebinning, 
        * first is 1D numpy, the others 2D arrays
        * self._the_runs_ is a list of lists of musr2py instances (psi/isis load routine)
        * inner list is for adding runs
        * outer list is suites of runs
        """ 
        # calculate time
        self.time = self.timebase() # in microseconds, 1D np.array 
        # calculate asymmetry and error
        for k, runs in enumerate(self._the_runs_): # initialize to zero    
            y, ey = self.single_asymmetry(runs, 
                                     self.grouping,
                                     self.alpha.value)    
            if ~k:
                self.asymm = y # 1D np.array
                self.asyme = ey # idem
                self.nrun = [runs[0].get_runNumber_int()]
            else:
                self.asymm = np.row_stack((self.asymm, y)) # columns are times, rows are successive runs (for multiplot and global)
                self.asyme = np.row_stack((self.asyme, ey))
                self.nrun.append(runs[0].get_runNumber_int())  # this is a list
                ######################################################
                # self.nrun list contains only the first run in case of run addition
                # used by save_fit (in file name), 
                #         write_csv (first item is run number)
                #         animate_fit (multiplot) 
                ######################################################
        # self.console('1st, last = {},{} bin for bckg subtraction'.format(self.firstbin,self.lastbin)) #dbg

    def asymmetry_multi(self):
        """
        * defines self.time, generates asymmetry end error without rebinning, 
        * first is 1D numpy, the others 2D arrays
        * self._the_runs_ is either one or a list of musr2py instances (psi/isis load routine)
        * list is for adding runs
        """ 
        # calculate time
        self.time = self.timebase() # in microseconds, 1D np.array 
        # calculate asymmetry and error
        self.asymm, self.asyme = self.multi_asymmetry(runs)   
        self.nrun = [runs[0].get_runNumber_int()]
                ######################################################
                # self.nrun list contains only the first run in case of run addition
                # used by save_fit (in file name), 
                #         write_csv (first item is run number)
                #         animate_fit (multiplot) 
                ######################################################
        # self.console('1st, last = {},{} bin for bckg subtraction'.format(self.firstbin,self.lastbin)) #dbg
    

    def asymmetry_multisuite(self):
        """
        * defines self.time, generates asymmetry end error without rebinning, 
        * first is 1D numpy, the others 3D arrays
        * self._the_runs_ is a list of lists of musr2py instances (psi/isis load routine)
        * inner list is for adding runs
        * outer list is suites of runs
        """ 
        # calculate time
        self.time = self.timebase() # in microseconds, 1D np.array 
        # calculate asymmetry and error
        for k, runs in enumerate(self._the_runs_): # initialize to zero    
            y, ey = self.multi_asymmetry(runs)   
            if ~k:
                self.asymm = y # 1D np.array
                self.asyme = ey # idem
                self.nrun = [runs[0].get_runNumber_int()]
            else:
                self.asymm = np.row_stack((self.asymm, y)) # columns are times, rows are successive runs (for multiplot and global)
                self.asyme = np.row_stack((self.asyme, ey))
                self.nrun.append(runs[0].get_runNumber_int())  # this is a list
                ######################################################
                # self.nrun list contains only the first run in case of run addition
                # used by save_fit (in file name), 
                #         write_csv (first item is run number)
                #         animate_fit (multiplot) 
                ######################################################
    

#########
# CONSOLE
#########

    def console(self,string,end = '\n'):
        '''
        printed output method::

            if self._output_ == self._outputtab_ prints on tab 
                   and sets self.mainwindow.selected_index = 5
            otherwise writes on terminal 

        ''' 
        import logging
        logger = logging.getLogger()

        if self._output_ == self._outputtab_:
            with self._output_:
                print(string)
                self.mainwindow.selected_index = 5
        else:
#            p = open(self._output_,'w') 
#            p.write(''.join([string,end]))
#            p.close()
            logger.info(''.join([string,end]))

###################
# CONSOLE TRACEBACK (to be deleted)
####################

#    def console_trcbk(self):
#        '''
#        method to print traceback in console

#        ''' 
#        import sys, traceback
#        if self._output_ == self._outputtab_:
#            with self._output_:
#                raise 
#                self.mainwindow.selected_index = 5
#        else:
#            p = open(self._output_,'w') 
#            formatted_lines = traceback.format_exc().splitlines()
#            for string in formatted_lines:
#                p.write(''.join([string,'\n']))
#            p.close()
        
##########################
# CREATE_RUNDICT
##########################
    def create_rundict(self,k=0):
        '''
        creates a dictionary rundict to identify and compare runs;
        refactored for adding runs

        '''
        rundict={}
        instrument = self.filespecs[0].value.split('_')[2] if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu' else self.filespecs[0].value # valid for PSI or ISIS
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
# FIT
##########################
    def fit(self, model_in = ''): # '' is reset to a 'daml' default before create_model(model_in')
                                  #  self.fit(model_in = 'mgmgbl') produces a different layout
        '''
        fit tab of mugui, used 
        - switchyard version: only single fit is active
        - to set: self.alpha.value, self.offset.value, forw and backw groups fit and plot ranges, model version           
        - to display: model name (new switchyard starts from here) 
        - to activate: fit, plot and update buttons
        - to select and load model (load from folder is missing/deprecated in python)
        - to select parameters value, fix, function, fft subtract

        ::

        # the calculation is performed in independent class mumodel, in library mucomponents
        # the methods are "inherited" by mugui 
        # via the reference instance self._the_model_, initialized in steps: 
        #     __init__ share initial attributes (constants) 
        #     _available_components_ automagical list of mucomponents 
        #     clear_asymmetry: includes reset check when suite is implemented
        #     create_model: lay out self._the_model_
        #     delete_model: for a clean start
        #     functions use eval, evil but needed, checked by muvalid, 
        #     int2_int to pass instances of mumodel methods
        #     iminuit2 requires them to be formatted by int2min as 
        #              numpy arrays fitval, fiterr, fitfix, fitlim, stored here as
        #              fitarg={"fitval":fitval,'"fiterr":fiterr, "fitfix":fitfix, "fitlim":fitlim}
        #     help  
        #     load 
        #     save_fit/load_ft save results in mujpy format (dill)
        #     write_csv produces a qtiplot/origin loadable summary
        # 
        # Six fit types: 
        #       single A1, multi single A2, 
        #       suite B1, multi suite B2, 
        #       suite C1, multi suite C2.
        # Multi optimizes a single chi2 over multiple groups of detectors
        # Single is the standard fit of one time-differential asymmetry, A1
        # Multi single optimizes a single chi2 over multiple groups, A2.
        # Suite either iterates a single fit over several runs, B1,
        #       or optimizes a single chi2 over the same runs, C1.
        # Multi suite either iterates a single fit over several runs with multiple groups, B2,
        #             or optimizes a single chi2 over the same runs, C2.
        # with global (one distinct value for all asymmetries, all runs and all groups), 
        #      glocal (one distinct value per group for all runs)
        #      or rlocal (one distinct value for all groups and each run) a parameters
        # another distinction is between 'internal' parameters from mumodel components
        # and 'external' parameters, both those extracted from each run header, such as T,B,theta
        # and those defined by the model, 
        # such as e.g. total asymmetry, fractions, T=0 order parameters, Tc etc
        '''


        from mujpy.mucomponents.mucomponents import mumodel
        import numpy as np


        def addcomponent(name,label):
            '''
            myfit = MuFit()
            addcomponent('ml') # adds e.g. a mu precessing, lorentzian decay, component
            this method adds a component selected from self.available_component, tuple of directories
            with zeroed values, stepbounds from available_components, flags set to '~' and empty functions
            plan also addgroupcomponents and addruncomponents (for A2, B2, C1, C2)
            '''
            from copy import deepcopy
            if name in self.component_names:
                k = self.component_names.index(name)
                npar = len(self.available_components[k]['pardicts']) # number of pars
                pars = deepcopy(self.available_components[k]['pardicts']) # list of dicts for 
                # parameters, {'name':'asymmetry','error':0.01,'limits':[0, 0]}

                # now remove parameter name degeneracy                   
                for j, par in enumerate(pars):
                    pars[j]['name'] = par['name']+label
                    pars[j].update({'value':0.0})
                    pars[j].update({'flag':'~'})
                    pars[j].update({'function':''}) # adds these three keys to each pars dict
                    # they serve to collect values in mugui
                self.model_components.append({'name':name,'pardicts':pars})
                return True # OK code
            else:
                self.console('\nWarning: '+name+' is not a known component. Not added.\n'+
                           'With myfit = mufit(), type myfit.help to see the available components')
                return False # error code

        def create_model(name):
            '''
            myfit = MuFit()       
            myfit.create_model('daml') # adds e.g. the two component 'da' 'ml' model
            this method 
            does not check syntax (prechecked by on_load_model)
            separates nexternafls and ngroups numbers from model name (e.g. '3mgml2' -> 'mgml', 3, 2)
            starts switchyard for A1,A1,B1, B2, C1, C2 fits
            adds a model of components selected from the available_component tuple of  
            directories
            with zeroed values, stepbounds from available_components, flags set to '~' and empty functions
            '''
            import string
            from mujpy.aux.aux import modelstrip
    # name 3mgbl2 for 3 global parameters (A0 f ), 0 kocal parameters (B end T) and two models 
    # e.g. alpha fit with a WTF and a ZF run, with two muon fractions of amplitude A0*R and A0*(1-R) respectively
    # find the three underscores in name by
    # [i for i in range(len(name)) if name.startswith('_', i)]

            model, nexternals, ngroups = modelstrip(name)
            self.console('created {}, with {} externals and {} groups'.format(model,nexternals,ngroups))
            components = [model[i:i+2] for i in range(0, len(model), 2)]
            # self.console('components are {}: {}'.format(len(components),components))
            self.model_components = [] # start from empty model
            for k,component in enumerate(components):
                label = string.ascii_lowercase[k] # was uppercase[k]
                if not addcomponent(component,label):
                    return False
                # self.console('create model added {}'.format(component+label))
            return True

#        def checkvalidmodel(name):
#            '''
#            checkvalidmodel(name) checks that name is either  
#            ::      A1, B1: 2*component string of valid component names, e.g.
#                                'daml' or 'mgmgbl'
#                                                                          
#            ::      or A2, B2: same, ending with 1 digit, number of groups (max 9 groups), 
#                                'daml2' or 'mgmgml2' (2 groups)
#            ::      or C1: same, beginning with 1 digit, number of external minuit parameters (max 9)
#                                '3mgml' (3 external parameters e.g. A, f, phi)
#            ::      or C2: same, both previous options
#                                '3mgml2' (3 external parameters, 2 groups)  
#            '''
#            from mujpy.aux.aux import modelstrip
#            try:
#                name, nexternals, ngroup = modelstrip(name)
#            except:
#                self.console('name error: '+name+' contains too many externals or groups (max 9 each)')
#                return [] 
#            # decode model
#            numberofda = 0
#            components = [name[i:i+2] for i in range(0, len(name), 2)]
#            for component in components: 
#                if component == 'da':
#                    numberofda += 1           
#                if numberofda > 1:
#                    self.console('name error: '+name+' contains too many da. Not added.')
#                    return [] # error code
#                if component not in self.component_names:
#                    self.console('name error: '+component+' is not a known component. Not added.\n'+
#                               'With myfit = mufit(), type myfit.help to see the available components')
#                    return [] # error code
#            return components

        def chi(t,y,ey,pars):
            '''
            stats for the right side of the plot 

            '''
            nu = len(t) - self.freepars # degrees of freedom in plot
            # self.freepars is calculated in int2min
            self._the_model_._load_data_(t,y,int2_int(),
                                         float(self.alpha.value),e=ey) 
                                         # int2_int() returns a list of methods to calculate the components
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
            from mujpy.aux.aux import derange, rebin, get_title
            from mujpy.aux.plot import plotile, set_bar#, animate_fit, init_animate_fit
            from scipy.stats import norm
            from scipy.special import gammainc
            import matplotlib.path as path
            import matplotlib.patches as patches
            import matplotlib.animation as animation
            from matplotlib import ticker

            font = {'family':'Ubuntu','size':10}
            P.rc('font', **font)

            ###################
            # PYPLOT ANIMATIONS
            ###################
            def animate_fit(i):
                '''
                anim function
                update errorbar data, fit, residues and their color,
                       chisquares, their histograms 

                '''
                # from mujpy.aux.aux import get_title
                # print('animate_fit')
                # nufit,ffit,chi2fit = chi(tfit[0],yfit[i],eyfit[i],pars[i])
                # nu,dum,chi2plot = chi(t[0],y[i],ey[i],pars[i])
                # color = next(self.ax_fit[(0,0)]._get_lines.prop_cycler)['color']
                line.set_ydata(y[i]) # begin errorbar
                line.set_color(color[i])
                line.set_markerfacecolor(color[i])
                line.set_markeredgecolor(color[i])
                segs = [np.array([[q,w-a],[q,w+a]]) for q,w,a in zip(t,y[i],ey[i])]
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
                    segs = [np.array([[q,w-a],[q,w+a]]) for q,w,a in zip(tlate,ylate[i],eylate[i])]
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

            def init_animate_fit():
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
                segs = [np.array([[q,w-a],[q,w+a]]) for q,w,a in zip(t,y[0],ey[0])]
                ye[0].set_segments(segs)
                ye[0].set_color(color[0]) # end errorbar
                fline.set_ydata(f[0]) # fit
                fline.set_color(color[0])
                res.set_ydata(y[0]-fres[0]) # residues
                res.set_color(color[0])

                if len(returntup)==5:
                    linel.set_ydata(ylate[0]) # begin errorbar
                    linel.set_color(color[0])
                    segs = [np.array([[q,w-a],[q,w+a]]) for q,w,a in zip(tlate,ylate[0],eylate[0])]
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
            # asymm[0] contains data from bin corresponding to
            # the prompt peak bin (nt0) + self.offset bins (fit_start = 0 to start from there)
            # pars is a list of lists of best fit parameter values
            # self.time, self.asymm, self.asyme are 2D arrays (self.times.shape[0] is always 1)
            # y, ey, f, fres, ylate, eylate, fl, flres, yfit, eyfit, ffit are 2D arrays
            # tf, tfl, tlate, tfit are 2D array like self.time (e.g. tf.shape[0] is always 1)

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
            [fitvalues,fiterrors,fitfixed,fitlimits,parameter_names] = int2min(dashboard)

            # fitarg  = int2min(return_names=True) # from dash, fitarg is a list of lists
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
                        #  self._the_model_ is an instance of mucomponents, and _load_data_ is one of its methods
                        self._the_model_._load_data_(self.time[0],self.asymm,int2_int(),
                                                     float(self.alpha.value),e=self.asyme)     # int2_int() returns a list of methods to calculate the components           
                    else:
                        # self.console('no global, mumodel load_data') # debug
                        self._the_model_._load_data_(self.time[0],
                                                     self.asymm[0],
                                                     int2_int(),
                                                     float(self.alpha.value),
                                                     e=self.asyme[0]) 
                        # int2_int() returns a list of methods to calculate the components
                        # string = 'Load data iok = {} [0=fine]. Time last = {:.2e}, A[0] = {:.2f}, alpha = {:.3f}'.format(pok,self._the_model_._x_[-1],self._themodel_._y_[0],self._the_model_._alpha_) # debug
                        # self.console(string)
            else:  # from self.lastfit, for best fit and plot best fit
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

#                if self.fig_tlog:
#                    self.fig_tlog.clf()
#                    self.fig_tlog,self.ax_tlog = P.subplots(num=self.fig_tlog.number) 
#            
#                self.fig_tlog,self.ax_tlog = P.subplots()
                   self.fig_fit.clf()
                   self.fig_fit,self.ax_fit = P.subplots(2,ncols,sharex = 'col', 
                             gridspec_kw = {'height_ratios':[3,1],'width_ratios':width_ratios},
                                            num=self.fig_fit.number)
                   self.fig_fit.subplots_adjust(hspace=0.05,top=0.90,bottom=0.12,right=0.97,wspace=0.03)
            except: # handle does not exist, make one
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
                # animate_fit must take care of updating parameters and producing correct fit curves
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
                    self.ax_fit[(0,1)].yaxis.set_major_formatter(ticker.NullFormatter())#set_yticklabels([])#
                    self.ax_fit[(1,1)].yaxis.set_major_formatter(ticker.NullFormatter())#set_yticklabels([])#
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
                self.ax_fit[(1,-1)].yaxis.set_major_formatter(ticker.NullFormatter())
                       #set_yticklabels(['']*len(self.ax_fit[(1,-1)].get_yticks()))#set_yticklabels([])#   
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


                # print('len(self.fitargs)) = {}'.format(len(self.fitargs)))
                #########################################################
                # animate_fit (see): TAKES CARE OF i>0 PLOTS IN 2D ARRAYS
                #                DOES ALSO UPDATE pars AND FIT CURVES
                #########################################################
                self.anim_fit = animation.FuncAnimation(self.fig_fit, 
                                                        animate_fit, 
                                                        range(0,len(self.fitargs)),init_func=init_animate_fit,
                                                        interval=anim_delay.value,repeat=True,
                                                        blit=False,
                                                        fargs=(self.ax_fit[(0,0)],
                                                        get_title(self._the_runs_[i][0]),) )# 


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
                # animate_fit must take care of updating parameters and producing correct fit curves               
                ##############
                # initial plot
                ##############f*
                yoffset = 0.4
                ymax = yoffset*fres.max() # this is an offset in the y direction, to be refined
                # self.console('ymax = {:.2f}'.format(ymax))
                rmax = 0.3*(y-fres).max()
                xoffset = 0.
                # print ('fres = {}'.format(fres.shape))
#                ttile, ytile, yres = plotile(t,y.shape[0],offset=xoffset), plotile(y,offset=ymax), plotile(y-fres,offset=rmax)  # plot arrays, ytile is shifted by ymax
#                tftile, ftile  = plotile(tf,y.shape[0],offset=xoffset), plotile(f,offset=ymax)  
                ytile, yres, ftile = plotile(y,offset=ymax), plotile(y-fres,offset=rmax), plotile(f,offset=ymax)   # y plot arrays shifted by ymax
                # self.console('after plotile ymax = {:.2f}'.format(ymax))
                # print('ttile.shape = {}, ytile.shape= {}, yres.shape = {}, tftile.shape = {}, ftile.shape = {}'.format(ttile.shape,ytile.shape,yres.shape,tftile.shape,ftile.shape))
                # print('f_tile = {}'.format(f_tile[0,0:50]))
                #############################
                # plot first (or only) window
                #############################
                # print(color)
                # errorbar does not plot multidim
                t1 = t.max()#/10.
                t0 = np.array([0.8*t1,t1])
                y0 = np.array([0.,0.])
                for k in range(y.shape[0]):
                    color = next(self.ax_fit[0,0]._get_lines.prop_cycler)['color']
                    self.ax_fit[(0,0)].errorbar(t[0],
                                                ytile[k],
                                                yerr=ey[k],
                                                fmt='o',
                                                elinewidth=1.0,ecolor=color,mec=color,mfc=color,
                                                ms=2.0,alpha=0.5) # data
                    self.ax_fit[(0,0)].plot(t0,y0,'-',lw=0.5,alpha=1,color=color) 
                    self.ax_fit[(1,0)].plot(t[0],yres[k],'-',lw=1.0,alpha=0.3,zorder=2,color=color) # residues 

                    self.ax_fit[(0,0)].plot(tf[0],ftile[k],'-',lw=1.5,alpha=1,zorder=2,color=color) # fit
                    y0 = y0 + ymax
                    self.ax_fit[(1,0)].plot(t[0],yzero,'k-',lw=0.5,alpha=0.3,zorder=0) # zero line
                    yzero = yzero + rmax
                ############################
                # plot second window, if any
                ############################
                y0 = np.array([0.,0.])
                if len(returntup)==5:
#                    tltile, yltile, ylres  = plotile(tlate,xdim=ylate.shape[0],
#                                             offset=xoffset),plotile(ylate,
#                                             offset=ymax),plotile(ylate-flres,offset=rmax)  # plot arrays, full suite
#                    tfltile, fltile  = plotile(tfl,offset=xoffset),plotile(fl,offset=ymax)  # res offset is 0.03                     
                    yltile, ylres, fltile = plotile(ylate,offset=ymax), plotile(ylate-flres,
                                                          offset=rmax),plotile(fl,offset=ymax)   
                                                          #  y arrays, late
                    for k in range(y.shape[0]):
                        color = next(self.ax_fit[0,1]._get_lines.prop_cycler)['color']
                        self.ax_fit[(0,1)].errorbar(tlate[0],
                                                    yltile[k],
                                                    yerr=eylate[k],
                                                    fmt='o',
                                                    elinewidth=1.0,ecolor=color,mec=color,mfc=color,
                                                    ms=2.0,alpha=0.5) # data
                        self.ax_fit[(1,1)].plot(tlate[0],ylres[k],'-',lw=1.0,alpha=0.3,zorder=2,color=color) # residues 
                        self.ax_fit[(0,1)].plot(tfl[0],fltile[k],'-',lw=1.5,alpha=1,zorder=2,color=color) # fit
                    # self.ax_fit[(1,1)].plot(tlate,tlate-tlate,'k-',lw=0.5,alpha=0.3,zorder=0) # zero line
                    self.ax_fit[(0,1)].set_xlim(tlate[0,0], tlate[-1,-1])

                ###############################
                # set title, labels
                ###############################
                ym,yM,rm,rM = ytile.min()-0.05,ytile.max()+0.01,yres.min()-0.005,yres.max()+0.005
                self.ax_fit[(0,0)].set_xlim(0,t.max())
                self.ax_fit[(0,0)].set_ylim(ym,yM)
                self.ax_fit[(1,0)].set_ylim(rm,rM)
                self.ax_fit[(1,0)].set_xlim(0,t.max())
                if len(returntup)==5:
                    self.ax_fit[(0,1)].set_xlim(tlate[0,0], tlate[-1,-1])
                    self.ax_fit[(0,1)].set_ylim(ym,yM)
                    self.ax_fit[(1,1)].set_xlim(tlate[0,0], tlate[-1,-1])
                    self.ax_fit[(1,1)].set_ylim(rm,rM)
                    self.ax_fit[(0,1)].yaxis.set_major_formatter(ticker.NullFormatter())
                    self.ax_fit[(1,1)].yaxis.set_major_formatter(ticker.NullFormatter())
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
                    self.ax_fit[(1,-1)].yaxis.set_major_formatter(ticker.NullFormatter())# set_yticklabels(['']*len(self.ax_fit[(1,-1)].get_yticks()))    #
                    self.ax_fit[(1,-1)].set_xlim([-5.5, 5.5])
                    mm = round(nu/4)
                    hb = np.linspace(-mm,mm,2*mm+1)
                    cc = gammainc((hb+nu)/2,nu/2) # muchi2cdf(x,nu) = gammainc(x/2, nu/2);
                    lc = 1+hb[min(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nufit
                    hc = 1+hb[max(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nufit
                    string = 'F-B: {} - {}\n'.format(self.group[0].value,self.group[1].value)
                    string += r'$\alpha=$ {}'.format(self.alpha.value)
                    string += '\n$\chi^2_f=$ {:.4f}\n ({:.2f}-{:.2f})\n$\chi^2_c=$ {:.4f}\n{} dof\n'.format(chi2fit,lc,hc,gammainc(chi2fit,nufit),nufit)
                    if  len(returntup)==5: 
                        string += '$\chi^2_e=$ {:.4f}\n$\chi^2_l=$ {:.4f}'.format(chi2plot,chi2late)
                    else:
                        string += '$\chi^2_p=$ {:.4f}'.format(chi2plot)
                    self.ax_fit[(0,-1)].text(-4.,0.3,string)
                else:
                    self.ax_fit[(0,0)].set_title(str(self.nrun[0])+': '+self.title.value)
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
                    ychi = yM
                    string = 'F-B: {} - {}\n'.format(self.group[0].value,self.group[1].value)
                    string += r'$\alpha=$ {}'.format(self.alpha.value)
                    self.ax_fit[(0,-1)].text(0.03,ychi,string)
                    dychi = (yM-ym)/(len(pars)+2) # trying to separate chi2
                    ychi -= 2*dychi
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
                        pedice = '_{'+str(self.nrun[k])+'}'
                        string = '$\chi^2'+pedice+'=$ {:.3f}'.format(chi2fit)
                        self.ax_fit[(0,-1)].text(0.03,ychi,string)
                        ychi -= dychi
                    self.ax_fit[(1,-1)].axis('off')
                    self.ax_fit[(0,-1)].set_ylim(self.ax_fit[(0,0)].get_ylim())
                self.ax_fit[(0,-1)].axis('off')

            self.fig_fit.canvas.manager.window.tkraise()# fig.canvas.manager.window.raise_()
            P.draw()

        def int2csv(k):
            '''
            Never used!
            translates nint into order numbers in the csv list, plus component number
            accepts single int, np.array of int, list of int
            e.g. model daml, da is nint=0, ncomp=0
                 in daml.1.3-4.csv the columns are
                 Run T1 eT1 T2 eT2 B da eda ...
                 0   1  2   3  4   5 6  7
                 int2csv(0) returns [[6,7,0]]
            '''
            from numpy import array
            ncomp = len(self.model_components)
            npar = [len(self.model_components[k]['pardicts']) for k in range(ncomp)]
            ntot = sum(npar)
            lmin = [-1]*ntot
            lcomp = [-1]*ntot
            nint = -1 # initialize   
            nmin = -1 # initialize
            
            for k in range(ncomp):  # scan the model
                for j in range(npar[k]): 
                    nint += 1  # internal parameter incremente always  
                    if self.flag[nint].value != '=': #  free or fixed par
                        nmin += 1
                        lmin[nint] = nmin # correspondence nint -> nmin 
                        lcomp[nint] = k # correspondence nint -> ncomp
            out = []
            lint = array(k) # to grant its ndarray
            if lint.shape: # to avoid single numbers
                for k in lint:
                    out.append([5+2*lmin[k]+1,5+2*lmin[k]+2,lcomp[k]])
            return out # lists of lists the inner list contains 
                       # the csv indices of par and its error, 
                       # and the internal index of their component
                    
                
        def int2_int(dashboard,the_model):
            '''
            input: the dashboard dict structure and the fit model instance
            output: a list of methods, in the order of the model components 
                    for the use of mumodel._add_.
            Invoked by the iMinuit initializing call
                     self._the_model_._load_data_, 
            just before submitting migrad, 
            self._the_model_ is an instance of mumodel 
             
            This function applies aux.translate to the parameter numbers in formulas
            since on the dash the numbers are internal numebers (each parameter of each component gets one)
            and in minuit some are missed because they are shared or determined by formulas 
            '''
            from mujpy.aux.aux import translate

            model_components = dashboard['model_components']

            ntot = sum([len(model_components[k]['pardicts']) for k in range(len(model_components))])
            lmin = [-1]*ntot # initialize the minuit parameter index of dashboard function indices 
            nint = -1 # initialize the number of internal parameters
            nmin = -1 # initialize the number of minuit parameters
            _int = []
            for k in range(len(model_components)):  # scan the model
                name = model_components[k]['name']
                # print('name = {}, model = {}'.format(name,self._the_model_))
                bndmthd = [] if name=='da' else the_model.__getattribute__(name) # this is the method to calculate a component, to set dalpha apart
                keys = []
                isminuit = []
                flag = [item.['flag'] for item in model_components[k]['pardicts']]
                function = [item.['function'] for item in model_components[k]['pardicts']]                
                for j in range(len(model_components[k]['pardicts'])): 
                    nint += 1  # internal parameter incremente always   
                    if flag[j] == '=': #  function is written in terms of nint
                        # nint must be translated into nmin 
                        string = translate(nint,lmin,function) # here is where lmin is used
                        keys.append(string) # the function will be eval-uated, eval(key) inside mucomponents
                        isminuit.append(False)
                    else:
                        nmin += 1
                        keys.append('p['+str(nmin)+']')                        
                        lmin[nmin] = nmin # 
                        isminuit.append(True)
                _int.append([bndmthd,keys]) 
            return _int

        def int2min(return_names=False):
            '''
            From internal parameters (dashboard) to minuit parameters.
            New version for iminuit >= 2.0
            Returns fitarg, a list lists:  
                [fitvalues,fiterrors,fitfixed,fitlimits] =           
                                   int2min(return_names=True) 
                    self.lastfit = Minuit(self._the_model_._chisquare_,
                                     name=self.minuit_parameter_names,
                                     *fitvalues)                                        
                    self.lastfit.errors = fiterrors
                    self.lastfit.limits = fitlimits
                    self.lastfit.fixed = fitfixed
                    self.lastfit.migrad()
            If returm_names = True, 
            self.minuit_parameter_names and self.minuit_component_names are updates.
            '''
            from mujpy.aux.aux import muvalue
            
            ntot = sum([len(self.model_components[k]['pardicts']) for k in range(len(self.model_components))])
            ntot -= sum([1 for k in range(ntot) if self.flag[k]=='=']) # ntot minus number of functions 
            # this is the number of minuit parameters
            fitval, fiterr, fitfix, fitlim = [], [], [], [] 
            
            fitarg = [] # list of dictionaries 
            parameter_name = []
            component_name = []
##########################################################
# this produces only the values, the names are obtained by 
##########################################################

#            for lrun in range(len(self._the_runs_)):
            nint = -1 # initialize
            # nmin = -1 # initialize
            free = 0
            # fitargs= {} # may not be used, check
            for k in range(len(self.model_components)):  # scan the model
                for j, par in enumerate(self.model_components[k]['pardicts']): # list of dictionaries, par is a dictionary
                    nint += 1  # internal parameter incremented always   
                    component_name.append(self.model_components[k]['name']) # name of component
                    if self.flag[nint].value == '~': #  skip functions, they are not new minuit parameter
                        # nmin += 1
                        free += 1
                        # lmin[nmin] = nint # correspondence between nmin and nint, is it useful?
                        # fitargs.update({par['name']:float(self.parvalue[nint].value)}) 
                        fitval.append(float(self.parvalue[nint].value))
                        parameter_name.append(par['name']) 
                        # fitargs.update({'error_'+par['name']:float(par['error'])})  
                        fiterr.append(float(par['error']))
                        if not (par['limits'][0] == 0 and par['limits'][1] == 0):
                            # fitargs.update({'limit_'+par['name']:par['limits']})
                            fitlim.append(par['limits'])
                        else:
                            fitlim.append((None,None))
                        #    fitargs.update({'limit_'+par['name']:None}) # not needed, defaut is None
                        fitfix.append(False)
                    elif self.flag[nint].value == '!':
                        # nmin += 1
                        # lmin[nmin] = nint # correspondence between nmin and nint, is it useful?
                        # fitargs.update({par['name']:float(self.parvalue[nint].value)})
                        fitval.append(float(self.parvalue[nint].value))
                        parameter_name.append(par['name'])
                        # fitargs.update({'fix_'+par['name']:True})
                        fitfix.append(True)
            fitarg.append(fitval)
            fitarg.append(fiterr)
            fitarg.append(fitfix)
            fitarg.append(fitlim)
            fitarg.append(component_name)
            fitarg.append(parameter_name)
            # self.console('fitval = {}\nfiterr = {}\nfitfix = {}\nfitlim = {}\ncomp name = {},\npar name = {} '.format(fitval,fiterr,fitfix,fitlim,component_name,parameter_name)) 
            self.freepars = free
            if free != ntot: self.console('check ntot = {} and free = {} in int2min, not equal!'.format(ntot,free)) 

 

            # print('fitargs= {}'.format(fitargs))
            if return_names:
                self.minuit_parameter_names = parameter_name
                self.minuit_component_names = component_name
            return fitarg
            # return [fitval, fiterr, fitfix, fitlim]
            # this routine returns only fitarg, a list of dict 
            # (must return lists of arrays fitvalues, fiterrors, fitfixed, fitlimits)
            # sets self.minuit_parameter_names, self.freepars (probably equal to Minuit.nfit)
            # a number of internal variables seem useless: lmin ntot (used to dimension lmin), nint, nmin
            # furthermore the name list is the same for all runs, can be done once, and then arrays filled
              
            
        def load_fit(b):
            '''
            loads fit values such that the same fit can be reproduced on the same data
            '''
            import dill as pickle
            from mujpy.aux.aux import path_file_dialog
            
            pre = self.paths[1].value
            path_and_filename = path_file_dialog(pre,'fit') # returns the full path and filename
            if path_and_filename == '':
                return
#            else:
#                self.console('Trying to read {}'.format(path_and_filename)) 
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
                # put self.lastfit after this one and it should be consistent:
                # initialized to [] = None at startup, equal to saved fit on load fit,
                # equal to last fit after a minimization
                self.lastfit = fit_dict['self.lastfit']
                self.version.value = fit_dict['version']
                self.offset.value = fit_dict['self.offset.value']
                self.model_components = fit_dict['self.model_components']
                self.grouping = fit_dict['self.grouping']
                set_group()
                self.alpha.value = fit_dict['self.alpha.value']
                # self.alphavalue = float(self.alpha.value)
                self.offset.value = fit_dict['self.offset.value']
                nint = fit_dict['nint']
                self.fit_range.value = fit_dict['self.fit_range.value']
                self.plot_range.value = fit_dict['self.plot_range.value'] # keys
                for k in range(nint+1):
                    self.parvalue[k].value = fit_dict['_parvalue['+str(k)+']']  
                    self.flag[k].value = fit_dict['_flag['+str(k)+    ']']  
                    self.function[k].value = fit_dict['_function['+str(k)+']']
                # self.fitargs = fit_dict['self.fitargs'] 
                # iminuit 2.0  get rid of
                self.load_handle.value = fit_dict['self.load_handle.value']

            except Exception as e:
                self.console('Problems with reading {} file\n\nException: {}'.format(path_and_filename,e))

#       DELETE  def min2fitarg(lastfit):
#            '''
#            reproduce fitarg (single fit) as from iminuit 1.0
#            as an ordered directory 
#            containing           true parameter name and its value
#                                 'error_'+true parameter name and its error
#                                 (not fixed and limits)
#            '''
#            from collections import OrderedDict
#            fitarg = OrderedDict()
#            value = lastfit.values
#            error = lastfit.errors
#            if len(value) != len(self.minuit_parameter_names):
#                self.console('len(value) = {} != len(name) in  min2fitarg'.format(
#                                 len(value),len(self.minuit_parameter_names)))
#                return fitarg
#            fullname = list(self.minuit_parameter_names) # a tuple made by int2min
#            for k in range(len(value)):
#                fitarg[fullname[k]] = value[k]
#                fitarg['error_'+fullname[k]] = error[k]
#                # name = getname(fullname[k]) detects the initial common part of the name
#            return fitarg
 
#  used only after fit, equivalent to lastfit.values, get rid
#           
        def min2int(fitargs):
            '''
            From minuit parameters to internal parameters,
            see int2min for a description   
            Invoked just after minuit convergence for save_fit, [on_update]
            '''
            # 
            # initialize
            #
            from mujpy.aux.aux import translate

            ntot = sum([len(self.model_components[k]['pardicts']) for k in range(len(self.model_components))])
            # ntot is number of 
            parvalue =  []
            lmin = [-1]*ntot 
            p = [0.0]*ntot 
            nint = -1
            nmin = -1
            for k in range(len(self.model_components)):  # scan the model
                keys = []
                for j, par in enumerate(self.model_components[k]['pardicts']): # list of dictionaries, par is a dictionary
                    nint += 1  # internal parameter incremented always   
                    if self.flag[nint].value != '=': #  skip functions, they are not new minuit parameter
                        nmin += 1
                        p[nmin] = fitargs[par['name']] # needed also by functions
                        parvalue.append('{:4f}'.format(p[nmin]))   # parvalue item is a string
                        lmin[nint] = nmin # number of minuit parameter
                    else: # functions, calculate as such
                        # nint must be translated into nmin 
                        string = translate(nint,lmin,self.function) # 
                        parvalue.append('{:4f}'.format(eval(string))) # parvalue item is a string
            return parvalue # list of parameter values 
            # iminuit 2.0 does not provide fitargs, use min2int to reproduce it 
            
        def min2print():
            '''
            From minuit parameters to plain print in console,
            see int2min for a description   
            Invoked just after minuit convergence for output 
            Does not need fitargs any more. Uses self.lastfit instead
            '''
            # 
            # initialize
            from mujpy.aux.aux import translate, translate_nint, value_error

            ntot = sum([len(self.model_components[k]['pardicts']) for k in range(len(self.model_components))])
            _parvalue =  []
            lmin = [-1]*ntot 
            p = [0.0]*ntot 
            e = [0.0]*ntot 
            nint = -1
            nmin = -1
            if isinstance(self._the_runs_[0], list):
                run = self._the_runs_[0][0].get_runNumber_int()
            else:
                run = self._the_runs_
            self.console('*****************************************************\n'+
                         '* Run: {} * {}*'.format(self.nrun[0],self.title.value))
            for k in range(len(self.model_components)):  # scan the model
                keys = []
                for j, par in enumerate(self.model_components[k]['pardicts']): # list of dictionaries, 
                                                                           # par is a dictionary
                    nint += 1  # internal parameter incremented always   
                    if self.flag[nint].value != '=': #  skip functions, they are not new minuit parameter
                        nmin += 1
                        p[nmin] = self.lastfit.values[nmin] # needed also by functions 
                        lmin[nint] = nmin              # number of minuit parameter
                        e[nmin] = self.lastfit.errors[nmin] # needed also by functions
                        # par['name'] is fullname (generic name + label, e.g. "λb"), which is unique
                        string = ' {}-{} = {}'.format(nmin,par['name'],value_error(p[nmin],e[nmin]))
                        self.console(string,end='')   # needed also by functions
                    else: # functions, calculate as such
                        # nint must be translated into nmin 
                        string = translate(nint,lmin,self.function) # 
                        stringe = string.replace('p','e')
                        nfix = translate_nint(nint,lmin,self.function)
                        stringa = '{}-{} = {}'.format(nfix,par['name'],value_error(eval(string),eval(stringe)))
                        self.console(stringa,end='') # needed also by functions
                self.console('')
            return
            # iminuit 2.0 must be adapted to fitargs not a dict any more

        def on_alpha_changed(change):
            '''
            observe response of fit tab widgets:
            validate float        
            '''
            string = change['owner'].value # description is three chars ('val','fun','flg') followed by an integer nint
                                               # iterable in range(ntot), total number of internal parameters
            try: 
                float(string) # self.alphavalue = 
            except:
                change['owner'].value = alpha0 # a string
            
#        def on_anim_check(change): 
#            '''
#            toggles multiplot and animations for suite fit
#            no effect on single run fit
#            '''
#            if change['owner'].value:
#                change['owner'].value=False
#            else:
#                change['owner'].value=True

        def on_fit_request(b):
            '''
            this is the entry to iminuit and fitplot (plot and log method) triggered by Fit Button
            retrieve data from the gui dashboard:
            self.alpha.value, range and pack with derange
            int2min obtains parameters values (parvalue[nint].value), flags (flag[nint].value), 
            errors, limits, functions (function[nint].value) in the dictionary format used by iminuit
            wrapped inside a list, to allow suite fits of more runs
            pass _int, generated by int2_int. to mumodel._add_ (distribute minuit parameters)
            obtain fitargs dictionary, needed by migrad, either from self.fitargs or from min2int
            pass them to minuit
            call fit_plot
            save fit file in save_fit
            write summary in write_csv 
            
            must be adapted to major change in iminuit 2.0, see 
            https://iminuit.readthedocs.io/en/stable/changelog.html#december-7-2020

            '''
            from iminuit import Minuit, cost
            from mujpy.aux.aux import derange, rebin, get_title, chi2std
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
                self.asymmetry_single() # prepare asymmetry
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
                # self.console('time.shape = {}, asymm.shape = {}, asyme.shape = {}'.format(time.shape,asymm.shape,asyme.shape))
                level = 0 # 0 (quiet),1,2
                self.fitargs = []  # CORREGGERE
                [fitvalues,fiterrors,fitfixed,fitlimits] = int2min(return_names=True) # from dash, 
                # with iminuit 2.0  they are four lists of arrays of (guess) values, one array per item 
                # relate to a single minimization, single run for the time being
                # self.console('start parameters:\n{}'.format(fitarg))
                
                if self._single_: # works fine
                    self._the_model_._load_data_(time,asymm,int2_int(),
                                                 float(self.alpha.value),
                                                 e=asyme) 
                                                 # pass data to model, one at a time
                    ############################## int2_int() returns a list of methods to calculate the components
                    # actual single migrad calls
                    self.lastfit = Minuit(self._the_model_._chisquare_,
                                     name=self.minuit_parameter_names,
                                     *fitvalues)                                        
                    self.lastfit.errors = fiterrors
                    self.lastfit.limits = fitlimits
                    self.lastfit.fixed = fitfixed
                    self.lastfit.migrad()
                    self.lastfit.hesse()
                    if self.lastfit.valid:
                        # self.console('self.lastfit.values[1] = {}'.format(self.lastfit.values[1]))
                        nu = len(time) - self.freepars # degrees of freedom in plot
                        # self.freepars is calculated in int2min
                        chi2 = self.lastfit.fval/nu  
                        # self.console('{}: {} ***** chi2 = {:.3f} ***** {}'.format(self.nrun[0],get_title(self._the_runs_[0][0]),chi2,strftime("%d %b %Y %H:%M:%S", localtime())))
                        min2print()
                         
                        lc,hc = chi2std(nu)
                        pathfitfile = save_fit(0)  # saves .fit file
                        if pathfitfile.__len__ != 2:
                            # assume pathfitfile is a path string, whose length will be definitely > 2
                            self.console('chi2r = {:.3f} ({:.3f} - {:.3f}) on {:.3f}-{:.3f} mus, saved in  {}'.format(chi2,lc,hc,time[0],time[-1],pathfitfile))
                        else:
                            # assume pathfitfile is a tuple containing (path, exception)
                            self.console('Could not save results in  {}, error: {}'.format(pathfitfile[0],pathfitfile[1]))
                        write_csv(chi2,lc,hc,0)  # writes csv file
                        self.console('*****************************************************')
                    else:
                        self.console('Fit not converged')
                else: # suites
                    valid = True
                    for k in range(len(self._the_runs_)): # loop over runs
                        self._the_model_._load_data_(time[0],asymm[k],int2_int(),
                                                     float(self.alpha.value),
                                                     e=asyme[k]) 
                                    # pass data to model, one run at a time
                        ############################## int2_int() returns a list of methods to calculate the components
                        # actual single migrad calls
                        # self.console('just before calling iminuit')
                        self.lastfit = Minuit(self._the_model_._chisquare_,
                                         name=self.minuit_parameter_names,
                                         *fitvalues[k])
                        self.lastfit.print_level = level
                        self.lastfit.errors = fiterrors[k]
                        self.lastfit.limits = fitlimits[k]
                        self.lastfit.fixed = fitfixed[k]
                        self.lastfit.migrad()
                        self.lastfit.hesse()
                        if not lastfit.valid:
                            valid = False
                        nu = len(time[0]) - self.freepars # degrees of freedom in plot
                        # self.freepars is calculated in int2min
                        chi2 = self.lastfit.fval/nu  
                        # self.console('{}: {} ***** chi2 = {:.3f} ***** {}'.format(self.nrun[k],get_title(self._the_runs_[k][0]),chi2,strftime("%d %b %Y %H:%M:%S", localtime())))
                        if k==0:
                            self.console('Fit on time interval ({:.3f},{:.3f}) mus'.format(time[0][0],time[0][-1]))
                        min2print()
                        lc,hc = chi2std(nu)
                        pathfitfile = save_fit(k)  # saves .fit file
                        if pathfitfile.__len__ != 2:
                            # assume pathfitfile is a path string, whose length will be definitely > 2
                            self.console('chi2r = {:.4f} ({:.4f} - {:.4f}), saved in  {}'.format(chi2,lc,hc,pathfitfile))
                        else:
                            # assume pathfitfile is a tuple containing (path, exception)
                            self.console('Could not save results in  {}, error: {}'.format(pathfitfile[0],pathfitfile[1]))
                        write_csv(chi2,lc,hc,k)  # writes csv file
                        ##############################
                    self.console('*****************************************************')
                    if not valid:
                        self.console('* CHECK: some fits did not converge!\n*****************************************************')
                    self.console('I am about to plot the bestfit')
                fitplot() # plot the best fit results 
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
            from mujpy.aux.aux import checkvalidmodel
            if checkvalidmodel(change['new']): # empty list is False, non empty list is True
                try:
                    del self._the_model_
                    self.fitargs=[] # so that plot understands that there is no previous minimization # CORREGGERE
                    model.value = loadmodel.value
                except:
                    pass
                self.fit(change['new']) # restart the gui with a new model
                self.mainwindow.selected_index = 0
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
 
        def on_plot_par(change):
            '''
            plot par wrapper
            '''
            string = change.description
            plotpar(x=string[9]) # 'B' or 'T'
            
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
            '''
            toggle start stop animation of suite fits
            '''
            if anim_check.value and not self._single_: 
                if change['owner'].value:
                    self.anim_fit.event_source.start()
                else:
                    self.anim_fit.event_source.stop()
        def on_update(b):
            '''
            update parvalue[k].value with last best fit results
            '''
            if self.lastfit: 
                # parvalue = min2int(self.fitargs[0]) # best fit parameters (strings)
                for k in range(len(parvalue)):
                    self.parvalue[k].value = self.lastfit.values[k]
            else:
                self.console('No fit yet! Cannot update')                 


        def print_nint(self):
            '''
            debug, delete me
            '''
            print(nint)

        def plotpar(x='X'): 
            '''
            plots parameters that have plotflag set to an integer 0<n<7
            x is 'T' or 'B'
            '''
            from mujpy.aux.aux import component, get_title
            from mujpy.aux.plot import plot_parameters
            import os
            from numpy import array
            import matplotlib.pyplot as P

            min2csv = lambda k: 2*k+6 # min2csv( kmin) returns the parameter column in the csv file 
                        
            # matplotlib initialization will provide same color for same component
            colors = P.rcParams['axes.prop_cycle'].by_key()['color']
            markers = ['o','d','^','V']

            nsubplots, lmin, lcomp, kminplot = [], [], [], [] # initialize service lists
            ylabels = [1]*6 # max number of labels
            kmin = -1
            # first loop, the second will be on 
            # parameters with plotflag, find them from dash with kint iterator
            # kmin is their minuit index
            # if plotflag <0 or excees 6 warning on console
            for kint in range(nint+1): # start from the dashboard parameters
                # discover which component it is (for color coherence)
                kcomp = component(self.model_components,kint) # aux method returns index of component 
                                                              # that kint belongs to
                if self.flag[kint]!='=': # this is a minuit parameter
                    kmin += 1 # count again minuit parameters (all but '=')
                    lmin.append(kint) # stored replica for finding kint from kmin 
                    if self.plotflag[kint].value>0 and self.plotflag[kint].value<=6:
                    # legal subplot (this is also a syntax checker).
                        # hence this minuit parameter must be plotted
                        kminplot.append(kmin) # store which kmin (for a plotted parameters iteration)
                        nsub = self.plotflag[kint].value  - 1 # which axis pythonized
                        nsubplots.append(nsub) # store in which axis, pythonic
                        lcomp.append(kcomp) # to remember which component kmin belongs to
                        if not isinstance(ylabels[nsub],str): # skip if ylabel already allocated
                            ylabels[nsub] = self.parname[kint].value[:-1].capitalize() # name without A,B,C,...
                          # ylabel: minuit pars, parname: dashboard parameters
                            if ylabels[nsub][-4:]=='rate': ylabels[nsub]+=' [mus-1]'
                            if ylabels[nsub][-4:]=='hase': ylabels[nsub]+=' [deg]'
                            if ylabels[nsub][-4:]=='ield': ylabels[nsub]+=' [mT]'
                    elif self.plotflag[k].value<0 or self.plotflag[k].value>6: # excluded, issue Warning
                        self.console('Warning: Parameter {} plot flag, {}, must be 0<=n<=6'.format
                                    (self.parname[k].value,self.plotflag[k].value))
                        
            if nsubplots:
                nsubs = int(array(nsubplots).max())+1 # effective subplots (nsubplots is 3 for flag 4)
            else:
                self.console('No plots requested: get fit results and select Panels')
                return
            nmin = kmin +1 # number of minuit parameters    
            # create dictionary labels, for plot_parameters (organizing subplots layout)
            #    the sample title
            #    the xlabel
            #    the ylabels from self.parname ù
            if x=='T':
                xlab = x+' [K]'
                title = get_title(self._the_runs_[0][0],notemp=True)
            elif x=='B':
                xlab = x+' [mT]'
                title = get_title(self._the_runs_[0][0],nofield=True)
            else:
                xlab = 'index'  
                title = get_title(self._the_runs_[0][0])
            labels ={'xlabel':xlab,
                    'title':title}
            labels['ylabels'] = ylabels # list of yaxis labels
            try:
                self.fig_pars.clf()
                self.fig_pars, ax = plot_parameters(nsubs,labels,self.fig_pars)
            except:
                self.fig_pars, ax = plot_parameters(nsubs,labels) # default fig=None

            version = str(self.version.value)
            strgrp = self.group[0].value.replace(',','_')+'-'+self.group[1].value.replace(',','_')
            path_csv = os.path.join(self.paths[1].value,model.value+'.'+version+'.'+strgrp+'.csv')
            p = np.genfromtxt(path_csv,comments='R')
            tb = p[:,1] if x=='T' else p[:,5] if x=='B' else np.arange(p.shape[0])
            etb = p[:,2] if x=='T' else np.zeros(p.shape[0])

            for kmp in range(len(kminplot)):   # do plotted parameters scan
                kmin = kminplot[kmp] # which minuit parameter is this
                kcsv = min2csv(kmin) # column of minuit par in p
                kcomp = lcomp[kmp]  # component of kmin parameter
                l = nsubplots[kmp] # plotflags are not pythonic (start from 1)
                ax[l].errorbar(tb,p[:,kcsv],yerr=p[:,kcsv+1],xerr=etb,fmt=markers[kcomp],color=colors[kcomp])
            if nsubs==5: ax[5].axis('off') # the last of 6
            for l in range(nsubs):
                if ax[l].get_ylabel()=='1':
                    ax[l].axis('off')
                ym,yM = ax[l].get_ylim()
                if ym>0: 
                    ax[l].set_ylim(0.,yM)
            if nsubs==2:
                if ax[0].axison: P.setp(ax[0].get_yticklabels()[0], visible=False) 
            elif nsubs==3:
                if ax[0].axison: P.setp(ax[0].get_yticklabels()[0], visible=False)  
                if ax[1].axison: P.setp(ax[1].get_yticklabels()[0], visible=False)  
            elif nsubs==4:
                if ax[0].axison: P.setp(ax[0].get_yticklabels()[0], visible=False)
                if ax[2].axison: P.setp(ax[2].get_yticklabels()[0], visible=False)  
            else: # 6
                if ax[0].axison: P.setp(ax[0].get_yticklabels()[0], visible=False)  
                if ax[1].axison: P.setp(ax[1].get_yticklabels()[0], visible=False)  
                if ax[3].axison: P.setp(ax[3].get_yticklabels()[0], visible=False)  
                if ax[4].axison: P.setp(ax[4].get_yticklabels()[0], visible=False)  
            self.fig_pars.canvas.manager.window.tkraise()
            P.draw()
            
        def save_fit(k):
            '''
            saves fit values such that load_fit can reproduce the same fit
            back to single run
            saves also unfitted models (with self.lastfit = [] equiv to None)
            '''
            import dill as pickle
            import os

            version = str(self.version.value)
            fittype = '' # single run fit
            # if self._global_: # global fit of run suite
            #    fittype = '.G.'
            # elif not self._single_: # sequential fit of run suite
            #    fyttype = '.S.'
            strgrp = self.group[0].value.replace(',','_')+'-'+self.group[1].value.replace(',','_')
            path_fit = os.path.join(self.paths[1].value,model.value+'.'+version+fittype+'.'+str(self.nrun[k])+'.'+strgrp+'.fit')
            # create dictionary setup_dict to be pickled 
            # the inclusion of self.load_handle will reload the data upon load_fit (?)
            names = ['self.alpha.value','self.offset.value',
                     'self.grouping','model.value',
                     'self.model_components','self.load_handle.value',
                     'version','nint', 'self.lastfit',
                     'self.fit_range.value','self.plot_range.value'] # keys
                     # iminuit 2.0 self.fitargs: eventually get rid (also min2fitargs)
            # self.fitargs
            fit_dict = {}
            for k,key in enumerate(names):
               fit_dict[names[k]] = eval(key) # key:value
            # _parvalue = min2int(self.fitargs[0]) # starting values from first bestfit
            for k in range(nint+1):
                if self.lastfit:
                    fit_dict['_parvalue['+str(k)+']'] = self.lastfit.values[k] # either fit 
                else:   
                    fit_dict['_parvalue['+str(k)+']'] = self.parvalue[k].value  # or dashboard
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
            from mujpy.aux.aux import get_title, spec_prec
            from time import localtime, strftime

            # print('k = {}, self.nrun = {}'.format(k,[j for j in self.nrun]))
            version = str(self.version.value)
            strgrp = self.group[0].value.replace(',','_')+'-'+self.group[1].value.replace(',','_')
            path_csv = os.path.join(self.paths[1].value,model.value+'.'+version+'.'+strgrp+'.csv')

            Bstr = self._the_runs_[k][0].get_field()
            if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu':
                TsTc, eTsTc = self._the_runs_[k][0].get_temperatures_vector(), self._the_runs_[k][0].get_devTemperatures_vector()
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
                row = form.format(self.nrun[k], TsTc[0],eTsTc[0],TsTc[1],eTsTc[1], Bstr[:Bstr.find('G')])
            elif self.filespecs[1].value=='nxs':
                Ts = self._the_runs_[k][0].get_temperatures_vector()
                n1 = '1'       
                form = '{} {:.'
                form += '{}'.format(n1)
                form += 'f} {}' #".format(value)
                row = form.format(self.nrun[k], Ts[0], Bstr[:Bstr.find('G')])

            k = -1
            for name in self.minuit_parameter_names:
                k += 1
                value, error = self.fitargs[k][name] , self.fitargs[k]['error_'+name]
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
            form += 'f} {} {}'
            row += form.format(chi2,chi2-lowchi2,hichi2-chi2,self.alpha.value,self.offset.value)
            row += ' {}'.format(strftime("%d %b %Y %H:%M:%S", localtime()))
            form =' {} {:.2f}'
            if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu':           
                for j in range(len(self.nt0)):
                    row += form.format(self.nt0[j],self.dt0[j])
            elif self.filespecs[1].value=='nxs': 
                row += form.format(self.nt0[0],self.dt0[0])
            row += '\n'
            # row is formatted with appropriate rounding, write directly
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


######### here starts the fit method of MuGui ####################################################
# no need to observe parvalue, since their value is a perfect storage point for the latest value #
# validity check before calling fit                                                              #
##################################################################################################
        from ipywidgets import Text, IntText, Layout, Button, HBox, \
                               Checkbox, VBox, Dropdown, ToggleButton, Label
        from mujpy.aux.aux import _available_components_
        
        self.lastfit = [] # initialize to None                               
        self.available_components =_available_components_() # creates list automagically from mucomponents
        self.component_names = [self.available_components[i]['name'] 
                                    for i in range(len(self.available_components))]
        # list of just mucomponents method name, used just in addcomponent
        self._the_model_ = mumodel() # local instance, need a new one each time a fit tab is reloaded (on_load_model)
        self._the_model_._alpha_ = [] # initialize to empty for plot guess sake

#------------------------- oneframe

        fit_button = Button(description='Fit',
                             tooltip='Execute the Minuit fit',
                             layout=Layout(width='15%'))                      # 15%
        fit_button.style.button_color = self.button_color
        fit_button.on_click(on_fit_request)

        loadmodel = Text(description='model',
                         description_tooltip='The name is the game: 2 letter-codes,\ne.g. daml, for da + ml,\n(see About)',
                         value='',
                         layout=Layout(width='12%'),continuous_update=False)  # this is where one can input a new model name
                                                                              # 27%
        loadmodel.observe(on_load_model,'value')
        loadmodel.style.description_width='37%'
        model = Text(value=model_in) # This is shown in the mainwindow, but is used to store present model
        
        # group[0] in oneframe , group[1] really ends up in in twopframe
        self.group = [Text(description='forward',
                           description_tooltip='list here all detectors\nin the forward group\ne.g. 1, 3, 5 , 7:9\nonly one : allowed',
                           layout=Layout(width='16%'),
                           continuous_update=False),                          # 43%
                      Text(description='backward',
                           description_tooltip='list here all detectors\nin the backward group\ne.g. 2, 4, 6 , 10:12\nonly one : allowed',
                           layout=Layout(width='16%'),
                           continuous_update=False)]
        # group and grouping: csv shorthand 
        set_group() # inserts shorthand from self.grouping into seld.group[k].value, k=0,1
        self.group[0].observe(on_group_changed,'value')
        self.group[0].style.description_width='38%'
        self.group[1].observe(on_group_changed,'value')
        self.group[1].style.description_width='38%'

        # self.alpha.value is Text and self.alfavalue is its float
        try:
            alpha0 = self.alpha.value
        except:
            alpha0 = '1.01' # generic initial value
        # self.alphavalue=float(alpha0)            
        self.alpha = Text(description='alpha',
                               description_tooltip='N_forward/N_backward\nfor the chosen grouping',
                               value=alpha0,
                               layout=Layout(width='10%'),
                               continuous_update=False) # self.alpha.value    # 53%
        self.alpha.observe(on_alpha_changed,'value')
        self.alpha.style.description_width='40%' 

        try:
            fit_range0 = self.fit_range.value 
        except:
            if not self._the_runs_:
                fit_range0 = '' # will be fixed at first data load
        self.fit_range = Text(description='fit range',
                              description_tooltip='start,stop[,pack]\ne.g. 0,20000,10\n(start from first good bin (see offset)\nuse 20000 bins\n pack them 10 by  10',
                              value=fit_range0,
                              layout=Layout(width='21.7%'),
                              continuous_update=False)                        # 75%
        self.fit_range.style.description_width='30%'
        self.fit_range.observe(on_range,'value')



        guesscheck = Checkbox(description='guess',
                              description_tooltip='Tick to plot\nstarting guess function',
                              value=False,
                              layout=Layout(width='8%'))                      # 83%
        guesscheck.style.description_width='2%'

        anim_delay = IntText(value=1000,
                             description='Delay (ms)',
                             description_tooltip='between successive runs\n(animation frames)',
                             layout=Layout(width='15%'))                       # 98%

        


#------------------------------- twoframe

        loadbutton = Button(description='Load fit',
                            tooltip='Opens GUI to load one of the existing fit templates',
                            layout=Layout(width='7.1%'))                         #  7.3%
        loadbutton.style.button_color = self.button_color
        loadbutton.on_click(load_fit)

        update_button = Button (description='Update',
                                tooltip='Update parameter starting guess\nfrom latest fit\n(must have fitted this model once).',
                                layout=Layout(width='7.3%'))                     # 15%
        update_button.style.button_color = self.button_color
        update_button.on_click(on_update)

        try:
            version0 = self.version.value 
        except:
            version0 = '1'
        self.version = Text(description='version',
                               description_tooltip='String to distinguish among model output files',
                               value=version0,
                               layout=Layout(width='11.5%'))#,indent=False)) # version.value is an int
                                                                               # 27%
        self.version.style.description_width='40%'

        # group[1] really ends up here   (width=16%)                           # 43%

        try:
            self.offset0 = self.offset.value
        except:
            self.offset0 = 7 # generic initial value 
        self.offset = IntText(description='offset',
                              description_tooltip='First good bin\n(number of bins skipped\nafter prompt peak)',
                              value=self.offset0,
                              layout=Layout(width='10%'),
                              continuous_update=False) # offset, is an integer
                                                                               # 53%
        self.offset.style.description_width='38%' 
        # initialized to 7, only input is from an IntText, integer value, or saved and reloaded from mujpy_setup.pkl

        try:
            self.plot_range0 = self.plot_range.value 
        except:
            if not self._the_runs_:
                self.plot_range0 = ''
        self.plot_range = Text(description='plot range',
                               description_tooltp='start,stop[,pack][,last,pack]\n0,20000,10 see fit range\n0,2000,10,20000,100 pack 10 up to bin 2000\npCK 100 from bin 2001 to bin 20000',
                               value=self.plot_range0,
                               layout=Layout(width='22%'),
                               continuous_update=False)                        # 75%
        self.plot_range.style.description_width='30%'
        self.plot_range.observe(on_range,'value')
         
        plot_button = Button (description='Plot',
                              tooltip='Generate a plot\n(see guess tooltip)',
                              layout=Layout(width='7.6%'))                       # 82.6%
        plot_button.style.button_color = self.button_color
        plot_button.on_click(on_plot_request)

        anim_check = Checkbox(description='Anim',
                         description_tooltip='Activate animation with Play button\ntoggle on/off movie with Loop button\nregulate frame Delay.',
                         value=False,
                         layout=Layout(width='7.6%'))                               # 92.2%
       # anim_check.observe(on_anim_check)
        anim_check.style.description_width='2%'
              
        anim_start = ToggleButton(description='off/on',
                                 tooltip='stop/start movie of fit plots\nfor run suite (see Run:)\nwith frame Delay\nwhen Anim is ticked.',
                                 value = False,
                                 layout=Layout(width='7.6%'))                   #  99.8%
        anim_start.observe(on_start_stop)
#        anim_start_value = False
#        anim_start.style.button_color = self.button_color
        
        
   
#-------------------- create the model template
# arrive here in three ways: model_in='' at start; model typed in loadmodel text box, checked by checkvalidmodel;
# load_fit
        if model_in == '':
            create_model('daml')
            model.value = 'daml'     # this sets the initial default model
        else:
            create_model(model_in)
        # self.console('in Fit main: components are {}'.format(len(self.model_components)))
        # this should always be a valid model, it is harwdired, model_in = 'daml'; from now on loadmodel is used and checked

#-------------------- fill model template into input widgets, two columns

#
#                    12345678901234567890123456789012345678901234567890123456789012345678901234567890
        s_n,s_nam,s_val,s_flag,s_func,s_plot ='Par n','Name','Value','Fit flag','Function','Panel'
        dashhead = HBox([Label(s_n,layout={'width':'8%','height':'16pt'},description_tooltip='Number to be used in Functions.'),
                        Label(s_nam,layout={'width':'18%','height':'16pt'}),
                        Label(s_val,layout={'width':'15%','height':'16pt'},description_tooltip='initial guess.'),
                        Label(s_flag,layout={'width':'11%','height':'16pt'},description_tooltip='~ is free\n! is fixed,= activates Function.'),
                        Label(s_func,layout={'width':'36%','height':'16pt'},description_tooltip='P[0] replaces by par 0 (quote only previous num).\nSimple algebra is allowed, e.g.\n0.5*P[0]+0.5*P[4].'),
                        Label(s_plot,layout={'width':'9%','height':'16pt'},description_tooltip='multipanel plot. \nChoose panel\nUse 1<=n<=6\nPars can share a panel.'),
                        ])
        leftframe_list, rightframe_list = [],[]
  
        words  = ['#','name','value','~!=','function']
        nint = -1 # internal parameter count, each widget its unique name
        ntot = np.array([len(self.model_components[k]['pardicts']) 
                         for k in range(len(self.model_components))]).sum()
        self.parname, self.parvalue, self.flag, self.function, self.plotflag = [] , [], [], [], [] # lists, index runs according to internal parameter count nint
        # self.compar = {} # dictionary: key nint corresponds to a list of two values, c (int index of component) and p (int index of parameter)
                # never use: self.compar[nint] is a list of two integers, the component index k and its parameter index j
        self.fftcheck = []            
        nleft,nright = 0,0
        # self.console('n of components {}'.format(len(self.model_components)))
        for k in range(len(self.model_components)):  # scan the model
            self.fftcheck.append(Checkbox(description='FFT',
                                          description_tooltip='uncheck for showing this component\nin the FFT of residues',
                                          value=True,
                                          layout=Layout(width='16%')))              
            self.fftcheck[k].style.description_width='2%'
            header = HBox([ Text(value=self.model_components[k]['name'],
                                 disabled=True,
                                 layout=Layout(width='8%')),                        #  8 %
                            self.fftcheck[k]]) # component header HBox              # 24 %
                                               # composed of the name (e.g. 'da') and the FFT flag
                                               # fft will be applied to a 'residue' where only checked components
                                               # are subtracted
            if k%2==0:                         # and ...
                leftframe_list.append(header)  # append it to the left if k even
                if k==0:
                      leftframe_list.append(dashhead)
            else:
                rightframe_list.append(header) # or to the right if k odd 
                if k==1:
                      rightframe_list.append(dashhead)
  
                                               # list of HBoxes, headers and pars
            nleftorright = 0                                   
            for j in range(len(self.model_components[k]['pardicts'])): # make a new par for each parameter 
                                                                                # and append it to component_frame_content
                nint += 1      # all parameters are internal parameters, first is pythonically zero 
                nleftorright += 1
                # self.compar.update({nint:[k,j]}) # stores the correspondence between nint and component,parameter
                nintlabel_handle = Text(value=str(nint),
                                        layout=Layout(width='7%'),
                                        disabled=True)                               # 7%
                name = self.model_components[k]['pardicts'][j]['name']
                baloon = ''
                if 'field' in name:
                    baloon = '[mT]'
                elif 'rate' in name:
                    baloon = '[mus-1]'
                elif 'phase' in name:
                    baloon = '[deg]'
                self.parname.append( 
                                    Text(value=name,
                                      description_tooltip=baloon,
                                      layout=Layout(width='16%'),
                                      disabled=True))                                 # 23%
                # self.console('\n{} - comp {} par {} appended'.format(nint,self.model_components[k]['name'],name))

                self.parvalue.append(
                                  Text(value='{:.4}'.format(self.model_components[k]['pardicts'][j]['value']),
                                  layout=Layout(width='15%'),
                                  description='value'+str(nint),
                                  continuous_update=False))                          # 38%
                self.parvalue[nint].style.description_width='0%'
                try:
                    self.parvalue[nint].value = _parvalue[nint]
                except:
                    pass
                # parvalue handle must be unique and stored at position nint, it will provide the initial guess for the fit

                self.flag.append(Dropdown(options=['~','!','='], 
                                 value=self.model_components[k]['pardicts'][j]['flag'],
                                 layout=Layout(width='11%'),
                                 description='flag'+str(nint)))                      # 49%
                self.flag[nint].style.description_width='0%'
                try:
                    self.flag[nint].value = _flag[nint]
                except:
                    pass
                 # flag handle must be unique and stored at position nint, it will provide (eventually) the nonlinear relation to be evaluated

                self.function.append(Text(value=self.model_components[k]['pardicts'][j]['function'],
                                     layout=Layout(width='36%'),
                                     description='func'+str(nint),
                                     description_tooltip='multipanel plot. \nChoose panel\nUse 1<=n<=6\nPars can share a panel.',
                                     continuous_update=False))                       # 85%
                self.function[nint].style.description_width='0%'

                self.plotflag.append(IntText(value=0,
                                             layout=Layout(width='9%')))      # 100%
                                             # self.plotflag[k].value is the subplots axis
                try:
                    self.function[nint].value = _function[nint]
                except:
                    pass
                # function handle must be unique and stored at position nint, it will provide (eventually) the nonlinear relation 

                fdis = False if self.model_components[k]['pardicts'][j]['flag']=='=' else True 
                self.function[nint].disabled = fdis # enabled only if flag='='
 
                # now put this set of parameter widgets for the new parameter inside a parameter HBox
                par_handle = HBox([nintlabel_handle, self.parname[nint], self.parvalue[nint], self.flag[nint], self.function[nint],self.plotflag[nint]],layout=Layout(width='100%'))
                           # handle to an HBox of a list of handles; notice that parvalue, flag and function are lists of handles
                
                # now make value flag and function active 
                self.parvalue[nint].observe(on_parvalue_changed,'value')
                self.flag[nint].observe(on_flag_changed,'value') # when flag[nint] is modified, function[nint] is z(de)activated
                self.function[nint].observe(on_function_changed,'value') # when function[nint] is modified, it is validated
 
                if k%2==0:                                         # and ...
                    leftframe_list.append(par_handle) # append it to the left if k even
                    nleft += nleftorright
                else:
                    rightframe_list.append(par_handle) # or to the right if k odd  
                    nright += nleftorright 

        if model_in == '': # at startup, model daml, 0 da, 1 A ml, 2 B ml, 3 ph ml, 4 lam ml
            self.parvalue[1].value = '0.2'  # asymmetry
            self.parvalue[2].value = '3.0'  # mT  
            self.parvalue[4].value = '0.2'  # mus-1 
        
        PT = Button(description='Plot vs. T',
                    tooltip='First perform suite vs. T fits.\nThis creates a csv file in the log folder\nSelect panel for pars you want plotted.',
                    layout=Layout(width='24%'))
        PT.on_click(on_plot_par)
        PT.style.button_color = self.button_color
        PB = Button(description='Plot vs. B',
                    tooltip='First perform suite vs. B fits\nThis creates a csv file in the log folder\nSelect panel for pars you want plotted.',
                    layout=Layout(width='24%'))
        PB.on_click(on_plot_par)
        PB.style.button_color = self.button_color

        PTB = VBox([HBox([Label('Show parameter plots',layout={'width':'34%'}),PT,PB],layout=Layout(width='100%'))],
                     layout={'border':'2px solid dodgerblue','align_items':'stretch'})
        blank = VBox([Label('',layout={'height':'32pt'})])
        if nleft<=nright:
            leftframe_list.append(blank)
            leftframe_list.append(PTB)
        else:
            rightframe_list.append(blank)
            rightframe_list.append(PTB)
        width = '99.8%'
        widthhalf = '100%'
        leftframe_handle = VBox(layout=Layout(width=widthhalf),
                                children=leftframe_list)#,layout=Layout(width='100%')
        rightframe_handle = VBox(layout=Layout(width=widthhalf),
                                 children=rightframe_list)# ,layout=Layout(width='100%')

#------------------ include frames in boxes        
        oneframe_handle = HBox(layout=Layout(width=width),
                               children=[fit_button,loadmodel,self.fit_range,guesscheck,anim_delay,self.group[0],
                               self.alpha])
        twoframe_handle = HBox(layout=Layout(width=width),
                               children=[loadbutton,update_button,self.version,self.plot_range,
                               plot_button,anim_check,anim_start,self.group[1],self.offset])
        bottomframe_handle = HBox(layout=Layout(width=width,border='1px solid dodgerblue'),
                                  children=[leftframe_handle,rightframe_handle])

        # end of model scan, ad two vertical component boxes to the bottom frame
        # backdoors
        self._load_fit = load_fit
        self._fit = fitplot
        self._int2_int = int2_int
        # now collect the handles of the three horizontal frames to the main fit window 
        self.mainwindow.children[0].children = [VBox([oneframe_handle,
                                                      twoframe_handle,
                                                      bottomframe_handle],layout=Layout(width=width))]#      
                                                # WHY DOES THE BOX EXTEND TO '100%' ?
                             # add the list of widget handles as the third tab, fit
        self.mainwindow.children[0].layout = Layout(border = '2px solid dodgerblue',width='100%')

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
            from mujpy.aux.aux import derange, autops, ps, _ps_acme_score, _ps_peak_minima_score, get_title#, animate_fft, init_animate_fft
            from mujpy.aux.plot import plotile
            from copy import deepcopy
            import matplotlib.pyplot as P
            from matplotlib.path import Path
            import matplotlib.patches as patches
            import matplotlib.animation as animation

            ###################
            # PYPLOT ANIMATIONS
            ###################
            def animate_fft(i):
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


            def init_animate_fft():
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
                    fftf_amplitude[0][start:stop], p0, p1, out = autops(fftf_amplitude[0][start:stop],'acme') # fix phase on theory 
                    if self.prntps.value:
                        self.console(out)
                    fft_amplitude[0][start:stop] = ps(fft_amplitude[0][start:stop], p0=p0 , p1=p1).real # apply it to data
                    for k in range(1,fft_amplitude.shape[0]):
                        fft_amplitude[k][start:stop] = ps(fft_amplitude[k][start:stop], p0=p0 , p1=p1)
                        fftf_amplitude[k][start:stop] = ps(fftf_amplitude[k][start:stop], p0=p0 , p1=p1)
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
                    self.anim_fft = animation.FuncAnimation(self.fig_fft, animate_fft, 
                                                            np.arange(0,len(self.fitargs)),
                                                            init_func=init_animate_fft,
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
        from ipywidgets import HBox, VBox, Layout, Button, Label, Text, IntText, Dropdown, Checkbox, ToggleButton

        # must inherit/retrieve self._the_model_, pars, fit_range = range(fit_start,fit_stop)
        # layout a gui to further retrieve  
        # fft_range (MHz), lb (us-1), real_amplitude (True/False) if False then power, autophase (True/False)

        # Layout gui
        fft_button = Button(description='Do FFT',layout=Layout(width='12%'))  # 12%
        fft_button.style.button_color = self.button_color
        fft_button.on_click(on_fft_request)

        filter0 = 0.3
        fft_filter = Text(description='Filter (mus-1)',
                               value='{:.4f}'.format(filter0),
                               layout=Layout(width='20%'),
                               continuous_update=False) # self.filter.value   # 32%
        fft_filter.observe(on_filter_changed,'value')

        fft_range0 = '0,50'
        fft_range = Text(description='fit range\nstart,stop\n  (MHz)',
                         value=fft_range0,
                         layout=Layout(width='28%'),
                         continuous_update=False)                              # 60%
        fft_range.style.description_width='60%'
        fft_range.observe(on_range,'value')

        real_or_power = Dropdown(options=['Real part','Power'], 
                                 value='Real part',
                                 layout=Layout(width='12%'))                    # 72%

        residues_or_asymmetry = Dropdown(options=['Residues','Asymmetry'], 
                                         value='Residues',
                                         layout=Layout(width='13%'))            # 85%

        autophase = Checkbox(description='Autophase',
                             value=True,
                             layout=Layout(width='15%'))                        #100%
        autophase.style.description_width='10%'

        anim_check = Checkbox(description='Animate',
                              value=False,
                              layout=Layout(width='12%'))                       # 12%
        anim_check.style.description_width = '1%'

        anim_delay = IntText(description='Delay (ms)',
                             value=1000,
                             description_tooltip='between frames',
                             layout=Layout(width='20%'))                        # 32%
                             
        anim_stop_start = ToggleButton(description='start/stop',
                                       description_tooltip='toggle animation loop',
                                       value=True,
                                       layout=Layout(width='12%'))               # 44%
        anim_stop_start.observe(on_start_stop,'value')
        
        empty = Label(layout=Layout(width='30%'))                                # 74%
        
        prntpslabel = Label('print autophase',
                            description_tooltip='verbose fit',
                            layout=Layout(width='10%'))     # 92%
        self.prntps = Checkbox(description=' ',
                      value=False,
                      layout=Layout(width='12%'))     # 100%

        fft_frame_handle = VBox(description='FFT_bar',children=[HBox(description='first_row',children=[fft_button,
                                                                fft_filter,
                                                                fft_range,
                                                                real_or_power,
                                                                residues_or_asymmetry,
                                                                autophase]),
                                     HBox(description='second_row',children=[anim_check, 
                                                                anim_delay,
                                                                anim_stop_start, 
                                                                empty,
                                                                self.prntps, 
                                                                prntpslabel])])
        # now collect the handles of the three horizontal frames to the main fit window 
        self.mainwindow.children[1].children = [fft_frame_handle] 
                             # add the list of widget handles as the third tab, fit
        self.mainwindow.children[1].layout = Layout(border = '2px solid dodgerblue',width='100%')
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
        
        Integrated with suite tab of mugui in v. 1.05
        used to select: 

            run (single/suite)
            load next previous, add next previous
             
        and to print: 

            run number, title, 
            total counts, group counts, ns/bin
            comment, start stop date, next run, last add                           
        '''

##########################
# SUITE
##########################

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
                    if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu':
                        self.nt0 = np.zeros(self._the_runs_[0][0].get_numberHisto_int(),dtype=int)
                        self.dt0 = np.zeros(self._the_runs_[0][0].get_numberHisto_int(),dtype=float)
                    elif self.filespecs[1].value=='nxs':
                        self.nt0, self.dt0 = np.array([0]), np.array([0])
                       
                    for j in range(self._the_runs_[0][0].get_numberHisto_int()):
                        self.nt0[j] = np.where(self._the_runs_[0][0].get_histo_array_int(j)==
                                               self._the_runs_[0][0].get_histo_array_int(j).max())[0][0]

                # self.nt0 exists
                self.title.value = get_title(self._the_runs_[0][0])                
                self.comment.value = self._the_runs_[0][0].get_comment() 
                self.start_date.value = self._the_runs_[0][0].get_timeStart_vector() 
                self.stop_date.value = self._the_runs_[0][0].get_timeStop_vector()
                self._the_runs_display.value = str(self.load_handle.value)
                # but if it is not compatible with present first run issue warning 
                if len(self.nt0)!=self._the_runs_[0][0].get_numberHisto_int(): # reset nt0,dt0
                    if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu':
                        self.nt0 = np.zeros(self._the_runs_[0][0].get_numberHisto_int(),dtype=int)
                        self.dt0 = np.zeros(self._the_runs_[0][0].get_numberHisto_int(),dtype=float)
                        for j in range(self._the_runs_[0][0].get_numberHisto_int()):
                            self.nt0[j] = np.where(self._the_runs_[0][0].get_histo_array_int(j)==
                                                   self._the_runs_[0][0].get_histo_array_int(j).max())[0][0]
                    elif self.filespecs[1].value=='nxs':
                        self.nt0, self.dt0 = np.array([0]), np.array([0])
                        
                     # self.console('WARNING! Run {} mismatch in number of counters, rerun prompt fit'.format(self._the_runs_[0][0].get_runNumber_int())) 

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
                #self.bin_range0 = '{},{}'.format(int(0.9*npk),int(1.1*npk))
                self.counterplot_range.value = self.bin_range0
            else:  # k > 0 
                self._single_ = False
                ok = [self._the_runs_[k][0].get_numberHisto_int() == self._the_runs_[0][0].get_numberHisto_int(),
                      self._the_runs_[k][0].get_binWidth_ns() == self._the_runs_[0][0].get_binWidth_ns()]
                if not all(ok): 
                    self._the_runs_=[self._the_runs_[0]] # leave just the first one
                    self.console('\nFile {} has wrong histoNumber or binWidth'.format(path_and_filename))
                    return -1  # this leaves the first run of the suite
            if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu':
                TdT = value_error(*t_value_error(k))
            else: # no error
                TdT = '{:.1f}'.format(t_value_error(k)[0])
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
            if ((self.filespecs[1].value=='bin' and len(runstr)>4) or (self.filespecs[1].value=='mdu' and len(runstr)>5)) or (self.filespecs[1].value=='nxs' and len(runstr)>8):
                self.console('Too long run number! {} {}'.format(self.filespecs[1].value,runstr))
                next_label.value = ''
                self.load_handle.value =''
            else:
                filename = ''
                nzeros = 4 if self.filespecs[1]=='bin' else (8 if self.filespecs[1]=='nxs' else 6)
                filename = filename.join([self.filespecs[0].value,
                                  muzeropad(runstr,nzeros),
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
            self.counternumber.value
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
            from mujpy.musr2py.musr2py import musr2py as psiload
            from mujpy.muisis2py.muisis2py import muisis2py as isisload
            from mujpy.aux.aux import muzeropad
            read_ok = 0
            runadd = []
            options = self.choose_tlogrun.options.copy() # existing dict (initialized to empty dict on dropdown creation) 
            if self.filespecs[1].value=='bin':            
                ndigits = 4 
            elif self.filespecs[1].value=='mdu':        
                ndigits = 5
                self.counternext_button.disable = True
                self.counternext_button.style.button_color = self.button_color_off
                self.counterprev_button.disable = True
                self.counternext_button.style.button_color = self.button_color_off                
            elif self.filespecs[1].value=='nxs': 
                ndigits = 8
                self.counternext_button.disable = False
                self.counterprev_button.disable = False
            for j,run in enumerate(runs): # run is a string containing a single run number
                if ((self.filespecs[1].value=='bin' and len(run)>4) or (self.filespecs[1].value=='mdu' and len(run)>5)) or (self.filespecs[1].value=='nxs' and len(run)>8):
                    self.console('Too long run number! {} {}'.format(self.filespecs[1].value,runstr))
                    read_ok = 99
                    break
                else:
                    filename = ''
                    filename = filename.join([self.filespecs[0].value,
                                          muzeropad(str(run),ndigits),
                                          '.',self.filespecs[1].value]) 
                    path_and_filename = os.path.join(self.paths[0].value,filename)
                    # data path + filespec + padded run rumber + extension)
                    if os.path.exists(os.path.join(self.paths[0].value,filename)):
                        if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu':            
                            runadd.append(psiload())  # adds a new instance of psiload() 
                            read_ok += runadd[j].read(path_and_filename) # THE RUN DATA FILE IS LOADED HERE
                        elif  self.filespecs[1].value=='nxs': 
                            try: 
                                #self.console('{} ns'.format(isisload(path_and_filename,'r').get_binWidth_ns()))
                                runadd.append(isisload(path_and_filename,'r'))  # adds a new instance of isisload
                                #self.console('path = {}'.format(path_and_filename))
                            except:
                                read_ok += 1
                        if read_ok==0:
                           	self.console('Run {} loaded'.format(path_and_filename))
                            # print('tlog dropdown position {} run {}'.format(str(j),str(run)))
                           	options.update({str(run):str(run)})
                           	# adds this run to the tlog display dropdown, 
                           	# on_loads_changed checks that tlog exists before value selection
                    else:
                   	    self.console('\nRun file: {} does not exist'.format(os.path.join(self.paths[0].value,filename)))
                   	    self.console('            if reading from afs check that klog is not expired!')
                   	    return -1
            if read_ok==0: # no error condition, set by musr2py.cpp or
                self.choose_tlogrun.options = options
                # ('self.choose_tlogrun.options = {}'.format(options))
                self._the_runs_.append(runadd) # 
                self.nrun.append(runadd[0].get_runNumber_int())
            else:
                try:
                    self.console('\nFile {} not read. Check paths, filespecs and run rumber on setup tab'.format(path_and_filename))
                except:
                    pass
                return -1 # this leaves the previous loaded runs n the suite 
            return check_runs(k)

        def on_load_nxt(b):
            '''
            load next run (if it exists)
            '''
            if self._single_:
                # print('self.nrun[0] = {}'.format(self.nrun[0]))
                self.load_handle.value=str(self.nrun[0]+1)
                # print('self.nrun[0] = {}'.format(self.nrun[0]))
            else:
                self.console('Cannot load next run (multiple runs loaded)')
                return -1 # this leaves the multiple runs of the suite 
                
        def on_load_prv(b):
            '''
            load previous run (if it exists)
            '''
            if self._single_:
                self.load_handle.value=str(self.nrun[0]-1)
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
                
        def on_exit(b):
            '''
            ditto
            '''
            from mujpy.aux.aux import exit_safe
            import os
            import signal
            from IPython.display import clear_output
            
            self.console('Exit request, called exit_safe')
            if exit_safe():
                if self.xterm: 
                    os.killpg(os.getpgid(self.xterm.pid),signal.SIGTERM)
                    clear_output()
                    self.gui.close()
                    
        def on_load_file(b):
            '''
            loads data from GUI and calls on loads_changed

            '''
            from mujpy.aux.aux import path_file_dialog, get_run_number_from
            
            pre,ext = self.paths[0].value,self.filespecs[1].value
            filename = path_file_dialog(pre,ext)
            if filename != '': # if empty, no choice made, nothing should happen
                run = get_run_number_from(filename,[self.filespecs[0].value,self.filespecs[1].value]) # path_file_dialog returns the full path and filename
                if run =='-1':
                    self.console('File mismatch: path = {}, filespecs = {}, {}'.format(p,f[0],f[1]))
                else:
                    self.load_handle.value = run
                    # self.console('Read {}'.format(self.load_handle.value))
                    #on_loads_changed('value')
            
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
            if self.load_handle.value==None:
                self.load_handle.value==''
            if self.load_handle.value=='': # either an accitental empty text return, or reset due to derun error
                return
            self._single_ = True
            self._the_runs_ = []  # it will be a list of muload() runs
            self.nrun = [] # it will contain run numbers (the first in case of run add)
            self.tlog_accordion.children[0].value=''
 
            #######################
            # decode the run string
            #######################           
            runs, errormessage = derun(self.load_handle.value) # runs is a list of lists of run numbers (string)
            if errormessage is not None: # derun error
                self.console('Run syntax error: {}. You typed: {}'.format(errormessage,self.load_handle.value))
                self.load_handle.value=''
                return

            ##################################
            # load a single run or a run suite
            ##################################
            read_ok = 0
            for k,runs_add in enumerate(runs):#  rs can be a list of run numbers (string) to add
                read_ok += add_runs(k,runs_add)                
                # print('on_loads_change, inside loop, runs {}'.format(self._the_runs_))
                if read_ok !=0:
                    return
#            if read_ok == 0:
                self.choose_nrun.options  = [str(n) for n in self.nrun]
                self.choose_nrun.value = str(self.nrun[0])
                options = self.choose_tlogrun.options.copy()
                runs = list(self.choose_tlogrun.options.keys())
                # scheme for popping items from list without altering index count
                kk = 0
                if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu':
                    for run in runs: # use original to iterate
                        if not tlog_exists(self.paths[2].value,run,4): # loading tlogs is optional
                            options.pop(run) # pop runs that do not have a tlog file
                    self.choose_tlogrun.options = options
                    try: 
                        self.choose_tlogrun.value = str((sorted(list(options.keys())))[0])
                    except:
                        if self.newtlogdir:
                            self.console("No tlogs")
                            self.newtlogdir = False
                self.get_totals() # sets totalcounts, groupcounts and nsbin
                # self.console('nextnext check next, self_single_ = {}'.format(self._single_))
            if self._single_:
                # self.console('next check next')
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
            if sum(weight)>0:
                weight = [w/sum(weight) for k,w in enumerate(weight)]
                t_value = sum([self._the_runs_[k][j].get_temperatures_vector()[self.thermo]*weight[j] for j in range(m)])
                t_error = sqrt(sum([(self._the_runs_[k][j].get_devTemperatures_vector()[self.thermo]*weight[j])**2 for j in range(m)])) if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu' else 0
            else:
                t_value, t_error = 0, 0
            return t_value, t_error
                                       
        from ipywidgets import Image, Text, Layout, HBox, Output, VBox, Tab
        from ipywidgets import Button, Label
        import os                      

        file = open(os.path.join(self.__logopath__,"logo.png"), "rb")
        image = file.read()
        logo = Image(value=image,format='png',width=132,height=132)
        
        ######################
        mujpy_width = '1024px'# '70%'  #
        ######################
        
        width = '100%'
#----------------- firstrow
        self.maxbin = Text(description='Bins',
                           description_tooltip='per histogram',
                           layout=Layout(width='12%'),
                           disabled=True)                              # 12%
        self.maxbin.style.description_width = '25%'                           

        self.comment = Text(description='Comment',
                                     layout=Layout(width='70%'),
                                     disabled=True)                    # 82%
        self.comment.style.description_width = '10%'


        self.nsbin = Text(description='ns/bin',
                          layout=Layout(width='18%'),
                          disabled=True)                               # 100%
        self.nsbin.style.description_width = '40%'


#----------------- secondtrow
        self._the_runs_display = Text(description='#',
                                 description_tooltip='Run number',
                                 value='none',
                                 layout=Layout(width='12%'),
                                 disabled=True)                        # 12%
        self._the_runs_display.style.description_width = '30%'

        self.title = Text(description='Title', 
                          value='none yet',
                          layout=Layout(width='48%'),
                          disabled=True)                               # 60%
        self.title.style.description_width = '14.5%'

        self.totalcounts = Text(value='0',
                                description='Total counts',
                                layout=Layout(width='20%'),
                                disabled=True)                         # 80%

        self.groupcounts = Text(value='0',
                                description='Group counts',
                                layout=Layout(width='20%'),
                                disabled=True)                         # 100%

#----------------- thirdrow
        self.start_date = Text(description='Start date',
                                     layout=Layout(width='28%'),
                                     disabled=True)                    # 28%
        self.start_date.style.description_width = '27%'
        exit_button = Button(description='Exit',
                           tooltip='cleanly exit mujpy',
                           layout=Layout(width='8%'))                   # 36%
        exit_button.on_click(on_exit)                   
        exit_button.style.button_color = self.button_color
        spacer = Label(description='     ',layout=Layout(width='1%'))   # x2 38%
        self.stop_date = Text(description='Stop date',
                                     layout=Layout(width='28%'),
                                     disabled=True)                    # 66%
        self.stop_date.style.description_width = '24%'
        
        Ap_button = Button(description='<Add',
                           tooltip='Add previous run\n(refers to Last+\nor to Next#)',
                           layout=Layout(width='7%'))                   # 73%
        Ap_button.on_click(on_add_prv)
        Ap_button.style.button_color = 'darkgray'#self.button_color
        
        An_button = Button(description='Add>',
                           tooltip='Add next run\n(refers to Last+\nor to Next#)',
                           layout=Layout(width='7%'))                   # 80%
        An_button.on_click(on_add_nxt)
        An_button.style.button_color = 'darkgray'#self.button_color

        last_add = Text(description='Last +',
                        description_tooltip='Last added run',
                        disabled=True,
                        layout=Layout(width='18%'))                     # 100%
        last_add.style.description_width = '27%'

                       # 100%
#------------------- fourthrow
        small = Label(layout=Layout(width='3.6%'))                      #  4%
        Ld_button = Button(description='Load',
                           tooltip='Opens a run file GUI\nin the data path (see setup tab)',
                           layout=Layout(width='8%'))                   # 12% 
        Ld_button.style.button_color = self.button_color
        Ld_button.on_click(on_load_file)

        self.load_handle = Text(description='Run:',
                      description_tooltip='Single run\ne.g. 431\nor run suites\ne.g. 431, 435:439\n or 431, 443+444',
                      layout=Layout(width='54%'),
                      continuous_update=False)                          # 66%
        self.load_handle.style.description_width='11%'
        self.load_handle.observe(on_loads_changed,'value') 

        Lp_button = Button(description='<Load',
                           tooltip='Loads previous run\n(if it exists)',
                           layout=Layout(width='8%'))                   # 74%
        Lp_button.on_click(on_load_prv)
        Lp_button.style.button_color = self.button_color
        Ln_button = Button(description='Load>',
                           tooltip='Loads Next # run\n(if it exists)',
                           layout=Layout(width='8%'))                   # 82%
        Ln_button.on_click(on_load_nxt)
        Ln_button.style.button_color = self.button_color
        next_label = Text(description='Next #',
                          description_tooltip='Next run to load',
                          disabled=True,
                          layout=Layout(width='18%'))                   # 100%
        next_label.style.description_width = '30%'


        firstrow = HBox(layout=Layout(width=width))
        firstrow.children = [self.maxbin, self.comment, self.nsbin]
        secondrow = HBox(layout=Layout(width=width))
        secondrow.children = [self._the_runs_display, self.title, self.totalcounts, self.groupcounts]
        thirdrow = HBox(layout=Layout(width=width))
        thirdrow.children = [self.start_date, spacer, exit_button, 
                             spacer, self.stop_date, spacer, Ap_button, 
                             spacer, An_button, spacer, last_add]
        fourthrow = HBox(layout=Layout(width=width))
        fourthrow.children = [small, Ld_button, self.load_handle,
                                   Lp_button, Ln_button, next_label]

        titlewindow = VBox(layout=Layout(width='100%')) 
        titlewindow.children = [firstrow, secondrow, thirdrow, fourthrow] 
        
        titlelogowindow = HBox(layout=Layout(width=mujpy_width))
        titlelogowindow.children = [logo, titlewindow]

        # main layout: tabs
        if self._output_==self._outputtab_:
            tabs_contents = ['fit', 'fft','setup',  'plots', 'output', 'about']# 'suite', 
        else:
            tabs_contents = [ 'fit', 'fft','setup',  'plots', 'about']# 'suite',
            
        tabs = [VBox(description=name,layout=Layout(border='2px solid dodgerblue')) for name in tabs_contents]
        self.mainwindow = Tab(children = tabs,layout=Layout(width=mujpy_width,border='2px solid dodgerblue'))# 

        self.mainwindow.selected_index = 2 # to stipulate that the first display is on tab 2, setup
        for i in range(len(tabs_contents)):
            self.mainwindow.set_title(i, tabs_contents[i])
            
            
        # Araba.Phoenix:
        self.gui = VBox(description='whole',layout=Layout(width='auto'))
        self.gui.children = [titlelogowindow, self.mainwindow]

        # This is the ex-suite tab to select run or suite of runs (for sequential or global fits)
        # moved to titlewindow thirdrow, an HBox containing various run load methods ['Single run','Run suite'] 


##########################
# OUTPUT
##########################
    def output(self):
        '''
        create an Output in terminal,
        or in widget in sixth tab, if needed      
        select by 
        self.mainwindow.selected_index = 5
        '''
        import os
        import platform
        import datetime
        import locale
        from shutil import which
        from subprocess import Popen, PIPE
        from ipywidgets import Output, HBox, Layout 
                     # Output(layout={'height': '100px', 'overflow_y': 'auto', 'overflow_x': 'auto'})

        self.xterm = []
        self._output_=''
        term_title = 'mujpy console -  '+datetime.datetime.today().strftime(locale.nl_langinfo(locale.D_T_FMT))
        if platform.system()=='Linux':
            terminal = ''
            if which('gnome-terminal') is not None:
                # terminal = 'gnome-terminal' # works on Ubuntu
                terminal = 'xterm'
#            if which('x-terminal-emulator') is not None:
#                terminal = 'x-terminal-emulator' # works on Debian, Ubuntu
            if terminal:
                self._output_ = "/tmp/mujpy_pipe"
                if not os.path.exists(self._output_):
                    os.mkfifo(self._output_)
                else:
                    os.unlink(self._output_)
                    os.mkfifo(self._output_)
                    # should write: 'gnome-terminal --title "mujpy console" -- bash -c "tail -f "'+self._output_
                # this is to open the xterm mujpy console on linux           
                self.xterm = Popen([terminal,'-T',term_title,'-e','tail -f %s' % self._output_],preexec_fn=os.setpgrp)
                #self.xterm = Popen([terminal,'--title',term_title,'shell=False','-e','tail -f %s' % self._output_],preexec_fn=os.setpgrp)
                # xterm  opens a terminal and executes tail -f on the open pipe self._output_
                self._outputtab_=[] # for Linux with terminal, if self._output_==self._outputtab_: 
                                    # Linux no terminal is redirected to tab 5 as other OSs

        elif platform.system()=='Windows': # put here also a cmd output option
            # this Windows and the Mac part bust be still checked
            # xterm = Popen(["cmd.exe","/k"],stdout=PIPE,stderr=PIPE)
            pass
        elif platform.system()=='Darwin':  # put here also a xterm option
            # xterm = Popen(['open', '-a', 'Terminal', '-n'],stdout=PIPE,stderr=PIPE)
            pass
        if not self._output_:
            self._outputtab_ = Output(layout={'height': '300px','width':'auto','overflow_y':'auto','overflow_x':'auto'})
            _output_box = HBox([self._outputtab_],layout=Layout(width='100%')) # x works y does scroll
            self._park_ = [_output_box] 
            # add the list of widget handles as the last tab, output
            self._output_ = self._outputtab_ # all in the tab
               

        self.console(''.join(['*****************************************************\n',
                              '* The output of mujpy is displayed here. Watch out! *\n',
                              '*                     DO NOT CLOSE!                 *\n',
                              '*****************************************************\n']))
                       

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

        def on_counterplot(change):
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
            self.numberHisto = self._the_runs_[0][0].get_numberHisto_int()
            returntup = derange(self.counterplot_range.value,self.histoLength) # 
            start, stop = returntup
            # abuse of get_grouping: same syntax here
            # counters = arange(self.numberHisto) 
            # now counters is an np.array of counter indices

            #############
            # load histos
            #############
            histo = zeros((self.numberHisto,stop-start),dtype=int)
            bins = arange(start,stop,dtype=int)

            # 4x4, 3x3 or 2x3 counters
            screen_x, screen_y = P.get_current_fig_manager().window.wm_maxsize() # screen size in pixels
            y_maxinch = float(screen_y)/dpi -0.5 # maximum y size in inches, 1 inch for window decorations 
            fx, f, f1 = 1., 4./5., 16./25. # fraction of screen for 
            P.close(fig=1)          
            if self.numberHisto > 9:
                nrows,ncols = 4,4
                x,y = fx*y_maxinch, y_maxinch 
                left,right,bottom,top,hspace,wspace = 0.06,0.97,0.05,0.95,0.02,0.02
            elif self.numberHisto > 6:
                nrows,ncols = 3,3
                x,y = fx*y_maxinch*f, y_maxinch*f
                left,right,bottom,top,hspace,wspace = 0.06,0.97,0.05,0.95,0.02,0.02
            elif self.numberHisto > 4:
                nrows,ncols = 2,3
                x,y = fx*y_maxinch*f, y_maxinch*f1
                left,right,bottom,top,hspace,wspace = 0.06,0.97,0.05,0.95,0.02,0.02
            elif self.numberHisto > 1:
                nrows,ncols = 2,2
                x,y = fx*y_maxinch*f1, y_maxinch*f1
                left,right,bottom,top,hspace,wspace = 0.06,0.97,0.05,0.95,0.02,0.02
            else:
                nrows,ncols = 1,1
                left,right,bottom,top,hspace,wspace = 0.06,0.97,0.05,0.95,0.02,0.02

                
            nplots = nrows*ncols    
            if change.description=='counters':  # this lets
            # buttons <> switch display of 16 counter banks
                self.bank = 0
                self.counterprev_button.disable = True
                self.counterprev_button.style.button_color = self.button_color_off
                self.counternext_button.disable = True                   
                self.counternext_button.style.button_color = self.button_color
            elif change.description=='>' and self.bank-1 < self.numberHisto/nplots:
                self.bank += 1
                self.counterprev_button.disable = False
                self.counterprev_button.style.button_color = self.button_color
                if (self.bank+1)*16>self.numberHisto:
                    self.counternext_button.disable = True
                    self.counternext_button.style.button_color = self.button_color_off
                else:
                    self.counternext_button.disable = False
                    self.counternext_button.style.button_color = self.button_color
            elif change.description=='<' and self.bank>0:
                self.bank -= 1
                self.counternext_button.disable = False                
                if self.bank==0:
                    self.counterprev_button.disable = True
                    self.counterprev_button.style.button_color = self.button_color_off
                else:
                    self.counterprev_button.disable = False                
                    self.counterprev_button.style.button_color = self.button_color 
            else:
                self.console('... requesting non existent detector banks')
                return         
            ##############################
            #  set or recover figure, axes 
            ##############################

#            if self.fig_counters:
#                self.fig_counters.clf()
#                self.fig_counters,self.ax_counters = P.subplots(nrows,ncols,figsize=(x,y),num=self.fig_counters.number) 
#            else:
            self.fig_counters,self.ax_counters = P.subplots(num=10,
                                                            nrows=nrows,
                                                            ncols=ncols,
                                                            figsize=(x,y),
                                                            dpi=dpi,
                                                            squeeze=False,
                                                            sharex=True,
                                                            sharey=False)
            self.fig_counters.canvas.set_window_title('Counters')
                
            for krun,nrun in enumerate(self.nrun):
                if nrun==int(self.choose_nrun.value):
                    this_run = krun # if not single selects which run to plot counters for 
            for k in range(nplots):
                kk = k+self.bank*16
                if kk < self.numberHisto: 
                    counter = k # already an index 0:n-1
                    for run in self._the_runs_[this_run]: # allow for add runs 
                        histo[counter] += run.get_histo_array_int(kk)[start:stop]
                    if k==0: ymax = histo[counter].max()
                    if stop-start<100:
                        self.ax_counters[divmod(counter,ncols)].bar(bins,
                                         histo[counter,:],edgecolor='k',color='silver',alpha=0.7,lw=0.7)
                        # self.console('Now histogramming histo {}'.format(kk))
                    else:
                        self.ax_counters[divmod(counter,ncols)].plot(bins,histo[counter,:],'k-',lw=0.7)
                        # self.console('Now plotting histo {}'.format(kk))
                    if counter>=(nrows-1)*ncols:
                        self.ax_counters[divmod(counter,ncols)].set_xlabel('bins')
                        P.setp(self.ax_counters[divmod(counter,ncols)].get_xticklabels()[-1], visible=False)
                    if divmod(counter,ncols)[1]==0:
                        self.ax_counters[divmod(counter,ncols)].set_ylabel('counts')
                    self.ax_counters[divmod(counter,ncols)].text(start+(stop-start)*0.7, ymax*0.9,'# '+str(kk+1)) # from index to label
                else:
                    self.ax_counters[divmod(k,ncols)].cla()
                    self.ax_counters[divmod(k,ncols)].axis('off')
            self.fig_counters.subplots_adjust(right=right,
                                              left=left,
                                              bottom=bottom,
                                              top=top,
                                              wspace=wspace,
                                              hspace=hspace)
# an unwanted Figure 1
            self.fig_counters.canvas.manager.window.tkraise()
            P.draw()          

            
        def on_multiplot(b):
            '''
            MULTIPLOT:
            produce plot
            '''
            import matplotlib.pyplot as P
            from numpy import array
            from mujpy.aux.aux import derange, rebin, get_title#, animate_multiplot, init_animate_multiplot
            import matplotlib.animation as animation

            ###################
            # PYPLOT ANIMATIONS
            ###################
            def animate_multiplot(i):
                '''
                anim function
                update multiplot data and its color 
                '''
                line.set_ydata(asymm[i])
                line.set_color(color[i])
                self.ax_multiplot.set_title(str(self.nrun[i])+': '+get_title(self._the_runs_[0][0]))
                return line, 


            def init_animate_multiplot():
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
            #  set figure, axes 
            ##############################
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
                self.anim_multiplot = animation.FuncAnimation(self.fig_multiplot, animate_multiplot, nruns, init_func=init_animate_multiplot,
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
            P.draw()

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
                # change['owner'].background_color = "mistyrose"
                if change['owner'].description[0:3] == 'plot':
                    change['owner'].value = self.plot_range0
                else:
                    change['owner'].value = self.bin_range0
            else:
                # change['owner'].background_color = "white"
                if returnedtup[1]>self.histoLength:
                    change['owner'].value=str(returnedtup[0],self.histoLength) if len(returnedtup)==2 else str(returnedtup[0],self.histoLength,returnedtup[2])         

        def on_nexttlog(b):
            '''
            select next run tlog 
            '''
            runs = list(self.choose_tlogrun.options.keys())
            runs = sorted([int(run) for run in runs])
            runindex = self.choose_tlogrun.index
            if runindex < len(runs):
                self.choose_tlogrun.value = str(runs[runindex+1])
                on_tlogdisplay([])

        def on_prevtlog(b):
            '''
            select prev run tlog
            '''
            runs = list(self.choose_tlogrun.options.keys())
            runs = sorted([int(run) for run in runs])
            runindex = self.choose_tlogrun.index
            if runindex > 0:
                self.choose_tlogrun.value = str(runs[runindex-1])
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
            if not self._the_runs_:
                self.console('Cannot plot temperatures: first load the runs!')
                return
            # the plotter requires one or two numpy arrays for as many sensors
            # plus a time array (in seconds)
            if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu':
                pathfile = os.path.join(self.paths[2].value,'run_'+muzeropad(self.choose_tlogrun.value)+'.mon')
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
            elif self.filespecs[1].value=='nxs':
                for k,run in enumerate(self.nrun):
                    if self.choose_tlogrun.value ==str(run):
                        break  # self.nrun is a list of int
                T1 = self._the_runs_[k][0].get_temperatures_vector()
                t = self._the_runs_[k][0].get_timeTemperature_vector() # in seconds
                t = t/60 # in minutes (float)
                starttime = self._the_runs_[k][0].get_timeStart_vector()
            ##############################
            #  set or recover figure, axes 
            ##############################
            try:
                if self.fig_tlog:
                    self.fig_tlog.clf()
                    self.fig_tlog,self.ax_tlog = P.subplots(num=self.fig_tlog.number) 
            except:
                self.fig_tlog,self.ax_tlog = P.subplots()
                self.fig_tlog.canvas.set_window_title('Tlogger')
            if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu':
                T1,T2 = array(T2), array(T1) # notice inversion! T1 is sample
                self.ax_tlog.plot(t,T2,'b-',label=r'$T_{\rm diff}$')
                T2ave, T2std = mean(T2), std(T2)
                self.ax_tlog.plot(tlim,[T2ave, T2ave],'b-',lw=0.5,alpha=0.8,label=r'$\langle T_{\rm sample}\rangle$')
                self.ax_tlog.fill_between(tlim, [T2ave-T2std,T2ave-T2std ],[T2ave, T2ave],facecolor='b',alpha=0.2)
                self.ax_tlog.fill_between(tlim, [T2ave+T2std,T2ave+T2std ],[T2ave, T2ave],facecolor='b',alpha=0.2)
            self.ax_tlog.plot(t,T1,'r-',lw=1.5,label=r'$T_{\rm sample}$')
            tlim,Tlim = self.ax_tlog.get_xlim(), self.ax_tlog.get_ylim()
            T1ave, T1std = mean(T1), std(T1)
            self.ax_tlog.plot(tlim,[T1ave, T1ave],'r-',lw=0.5,alpha=0.8,label=r'$\langle T_{\rm diff}\rangle$')
            self.ax_tlog.fill_between(tlim, [T1ave-T1std,T1ave-T1std ],[T1ave, T1ave],facecolor='r',alpha=0.1)
            self.ax_tlog.fill_between(tlim, [T1ave+T1std,T1ave+T1std ],[T1ave, T1ave],facecolor='r',alpha=0.1)
            self.ax_tlog.set_title('Run '+self.choose_tlogrun.value+' started at '+starttime)
            self.ax_tlog.set_xlabel('time [min]')
            self.ax_tlog.set_ylabel('T [K]')
            self.ax_tlog.set_xlim(min(t),max(t))
            if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu':
                self.ax_tlog.legend()
            if Tlim[1]-Tlim[0] < 1.:
                T0 = (Tlim[0]+Tlim[1])/2
                self.ax_tlog.set_ylim(T0-1.,T0+0.5)
            y1,y2 = self.ax_tlog.get_ylim()
            if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu':
                for x1,x2 in zip(pause,go):
                    self.ax_tlog.fill_between([x1,x2], [y1,y1 ],[y2,y2],facecolor='k',alpha=0.5)
                    self.ax_tlog.text(x1*0.9+x2*0.1,y1*0.9+y2*0.1,'PAUSE',color='w')
            self.fig_tlog.canvas.manager.window.tkraise()
            P.draw()

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
                         Label(layout=Layout(width='3%'))],layout=Layout(width='100%',border='2px solid dodgerblue'))
        ###################
        # counters inspect
        ###################
        counterlabel = Label(value='Inspect',layout=Layout(width='7%'))# count
        counterplot_button = Button(description='counters',layout=Layout(width='10%'))
        counterplot_button.on_click(on_counterplot)
        counterplot_button.style.button_color = self.button_color
        self.counternumber = Label(value='{} counters per run'.format(' '),layout=Layout(width='15%'))
        self.counternext_button = Button(description='>',layout=Layout(width='10%'))
        self.counternext_button.on_click(on_counterplot)
        self.counternext_button.style.button_color = self.button_color_off
        self.counterprev_button = Button(description='<',layout=Layout(width='10%'))
        self.counterprev_button.on_click(on_counterplot)
        self.counternext_button.style.button_color = self.button_color_off
        self.counterprev_button.style.button_color = self.button_color_off
        self.counternext_button.disable = True
        self.counterprev_button.disable = True


#        counter_range = Text(description='counters',
#                               value='', continuous_update=False,layout=Layout(width='20%'))
#        counter_range.style.description_width='33%'
#        counter_range.observe(on_counter,'value')
        self.counterplot_range = Text(description='bins: start,stop',
                               value=self.bin_range0,continuous_update=False,layout=Layout(width='25%'))
        self.counterplot_range.style.description_width='40%'
        self.counterplot_range.observe(on_range,'value')
        self.choose_nrun = Dropdown(options=[], description='run', layout=Layout(width='15%'))
        self.choose_nrun.style.description_width='25%'
        counterbox = HBox([Label(layout=Layout(width='3%')),
                           counterlabel,
                           counterplot_button,
                           Label(layout=Layout(width='3%')),
                           self.counternumber,
                           self.counterprev_button,
                           self.counternext_button,
                           self.counterplot_range,
                           self.choose_nrun],
                           layout=Layout(width='100%',border='2px solid dodgerblue'))

        ##########
        # TLOG PSI
        ##########
        spacer = Label(layout=Layout(width='3%'))
        tloglabel = Label(value='Tlog',layout=Layout(width='7%'))
        tlog_button = Button(description='display',layout=Layout(width='10%'))
        tlog_button.on_click(on_tlogdisplay)
        tlog_button.style.button_color = self.button_color
        options = {} # empty slot to start with
        self.choose_tlogrun = Dropdown(options=options,description='Tlog run', layout=Layout(width='15%')) #
        self.choose_tlogrun.style.description_width='35%'
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

        tlogbox = HBox([spacer,tloglabel,
                      tlog_button,
                      self.choose_tlogrun,
                      nexttlog_button,prevtlog_button,
                      self.tlog_accordion],layout=Layout(width='100%',border='2px solid dodgerblue'))

        vbox = VBox(layout=Layout(border='2px solid dodgerblue'))
        vbox.children = [multibox, counterbox, tlogbox]
        self.mainwindow.children[3].children = [vbox]
        self.mainwindow.children[3].layout = Layout(border = '2px solid dodgerblue',width='100%')
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
            Also invoked at gui layout. 
            If mujpy_setup.pkl does not exist, produces one.
            """
            import dill as pickle
            import os
            from numpy import array
                                        
            path = os.path.join(self.__startuppath__,'mujpy_setup.pkl')
            # self.console('loading {}, presently in {}'.format(path,os.getcwd()))
            
            iok = 0
            try: # if mujpy_setup.pkl exists and is valid setup file
                with open(path,'rb') as f:
                    mujpy_setup = pickle.load(f) 
            except:
                self.console('mujpy_setup file not found or unreadable')
                iok = -1
            try:
                self.nt0 = mujpy_setup['self.nt0'] # peak bin, shape run.get_numberHisto_int()
                self.dt0 = mujpy_setup['self.dt0'] # bin fractions, shape run.get_numberHisto_int()           
            except:
                self.console('mujpy_setup file not found or unreadable')
                iok = -1
               
            try:
                # now find out which facility
                self.filespecs[0].options = mujpy_setup['_filespecs_options']
                for k in range(2):  # len(filespecs.content)
                    self.filespecs[k].value = mujpy_setup['_filespecs_content'][k] # ('fileprefix','extension')
            except:
                self.console('unrecognized file specs: {}xxxx.{}'.format(mujpy_setup['_filespecs_content'][0],mujpy_setup['_filespecs_content'][1]))
                iok = -1
            try:
                self.nt0_run = mujpy_setup['self.nt0_run'] # dictionary, checks run belongs to same setup
            except:
                self.console('nt0_run dictionary missing in mujpy_setup')
                iok = -1
                
            if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu': # PSI bulk musr
                try:
                    for k in range(3):  # len(self.paths)=len(mujpy_setup['_paths_content'])
                        self.paths[k].value =  mujpy_setup['_paths_content'][k] 
                    datapath.description_tooltip += '\n'+self.paths[0].value
                    logpath.description_tooltip += '\n'+self.paths[1].value
                    tlogpath.description_tooltip += '\n'+self.paths[2].value
                except:
                    self.console('unrecognized paths: data {}'.format(mujpy_setup['_paths_content'][0]))
                    iok = -1               
                try:
                    for k in range(2):  # len(bkg_content)
                        self.prepostpk[k].value = mujpy_setup['_prepostpk'][k] # pre-prompt, post-prompt bins 
                    self.lastbin = mujpy_setup['self.lastbin'] # upper edge for background
                    fit_button.tooltip = 't0 is fitted as prompt centre\n(peak on step function)'
                    self.prepostpk[0].disabled=False
                    self.prepostpk[1].disabled=False
                except:
                    self.console('prepost, etc. missing in mujpy_setup')
                    iok = -1
            elif self.filespecs[1].value=='nxs': # ISIS
                try:
                    for k in range(2):  # len(self.paths)=len(mujpy_setup['_paths_content'])
                        self.paths[k].value =  mujpy_setup['_paths_content'][k] 
                    self.paths[2].value = 'None'
                    fit_button.tooltip = 't0 is fitted on sum\nof all detectors\nas centre of sigmoid\n convoluted with pion decay'
                    self.prepostpk[0].value=0
                    self.prepostpk[1].value=0
                    self.prepostpk[0].disabled=True
                    self.prepostpk[1].disabled=True
                except:
                    self.console('ISIS mujpy_setup.pkl is corrupted')
                    self.console('unrecognized paths: data {}'.format(mujpy_setup['_paths_content'][0]))
                    iok = -1

            if iok!=0: # if mujpy_setup.pkl doesn't exist create a mock PSI setup
                self.console('File {} not found'.format(path))
                self.console('Creating a PSI GPS setup. To switch to ISIS choose file extension nxs')
                self.filespecs[0].options=['deltat_tdc_gps_',
                                           'deltat_tdc_ltf_',
                                           'deltat_tdc_dol_',
                                           'deltat_tdc_gpd_']
                for k in range(3):  # len(paths_contents)
                    self.paths[k].value =  './' 
                self.filespecs[0].value = 'deltat_tdc_gps_' # 'fileprefix'
                self.filespecs[1].value = 'bin' #,'extension')
                self.prepostpk[0].value = 7 # 'pre-prompt bin'
                self.prepostpk[1].value = 7 # 'post-prompt bin' 
                self.nt0,self.dt0 = array([0.]),array([0.])
                self.lastbin = [] # upper edge for background
                self.nt0_run = {} # dictionary to identify runs belonging to the same setup

        def on_extension_changed(change):
            '''
            when user changes this setup tab widget mugui
            changes the optiond for file prefixes
            '''
            if self.filespecs[1].value == 'bin':
                self.filespecs[0].options=['deltat_tdc_gps_',
                                           'deltat_tdc_ltf_',
                                           'deltat_tdc_dol_',
                                           'deltat_tdc_gpd_']
                self.prepostpk[0].disabled=False       
                self.prepostpk[1].disabled=False      
            if self.filespecs[1].value == 'mdu':
                self.filespecs[0].options=['tdc_hifi_2021_',
                                           'tdc_hifi_2020_',
                                           'tdc_hifi_2019_',
                                           'tdc_hifi_2018_',
                                           'tdc_hifi_2017_',
                                           'tdc_hifi_2016_',
                                           'tdc_hifi_2015_',
                                           'tdc_hifi_2022_']
                self.prepostpk[0].disabled=False       
                self.prepostpk[1].disabled=False      
            elif self.filespecs[1].value == 'nxs':
                self.filespecs[0].options=['EMU',
                               'MUS',
                               'HIFI',
                               'ARGUS']
                self.prepostpk[0].disabled=True        
                self.prepostpk[1].disabled=True        

        def on_paths_changed(change):
            '''
            when user changes this setup tab widget mugui
            checks that paths exist, in case it creates log path
            '''
            import os

            path = change['owner'].description # description is paths[k] for k in range(len(paths)) () 
            k = paths_content.index(path) # paths_content.index(path) is 0,1,2 for paths_content = 'data','analysis','tlog'
            if k==2: # tlog
                self.newtlogdir = True
            directory = self.paths[k].value # self.paths[k] = handles of the corresponding Text
            if not os.path.isdir(directory):
                if k==1: # analysis, if it does not exist mkdir
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
            # 59 is the end of the standard message
            self.path[k].description = self.path[k].description[:59]+'\n'+self.path[k].value
 
        def on_prompt_fit_click(b):
            '''
            when user presses this setup tab widget mugui
            executes prompt fits
            '''
            if self._the_runs_:
                promptfit() # mplot mprint we leave always True False, respectively
            else:
                self.console('Cannot fit prompts without data set. Please load the data first')

        def promptfit(mplot = True, mprint = False):
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
            from iminuit import Minuit, cost
            
            import matplotlib.pyplot as P
            from mujpy.mucomponents.muprompt import muprompt
            from mujpy.mucomponents.muedge import muedge
            from mujpy.aux.aux import scanms, step 

            font = {'family' : 'Ubuntu','size'   : 8}
            P.rc('font', **font)
            dpi = 100.
            if self.filespecs[1].value=='bin':
                second_plateau = 100
                peakheight = 100000.
                peakwidth = 1.
            elif self.filespecs[1].value=='mdu':
                first_plateau = - 500
                second_plateau = 1500
            
            if not self._the_runs_:
                self.console('No run loaded yet! Load one first (select suite tab).')
            else:
                if self.filespecs[1].value=='bin':  # PSI
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
                    nbin =  max(npeaks) + second_plateau # this sets a counter dependent second plateau bin interval
                    x = np.arange(0,nbin,dtype=int) # nbin bins from 0 to nbin-1
                    self.lastbin, np3s = npeaks.min() - self.prepostpk[0].value, npeaks.max() + self.prepostpk[1].value # final bin for background average, first bin for right plateau estimate (last is nbin)


                    x0 = np.zeros(self._the_runs_[0][0].get_numberHisto_int()) # for center of peaks
                    if mplot:

                        ###################
                        #  set figure, axes 
                        ###################
                        self.fig_counters,self.ax_counters = P.subplots(2,3,figsize=(7.5,5),dpi=dpi)
                        self.fig_counters.canvas.set_window_title('Prompts t0 fit')
                        screen_x, screen_y = P.get_current_fig_manager().window.wm_maxsize() # screen size in pixels
                        y_maxinch = float(screen_y)/dpi #  maximum y size in inches

                        prompt_fit_text = [None]*self._the_runs_[0][0].get_numberHisto_int()   
                             
                    for counter in range(self._the_runs_[0][0].get_numberHisto_int(),sum(self.ax_counters.shape)):
                        self.ax_counters[divmod(counter,3)].cla()
                        self.ax_counters[divmod(counter,3)].axis('off')

                    for counter in range(self._the_runs_[0][0].get_numberHisto_int()):
                        # prepare for muprompt fit
                        histo = np.empty(self._the_runs_[0][0].get_histo_array_int(counter).shape)
                        for k in range(len(self._the_runs_[0])): # may add runs
                            histo += self._the_runs_[0][k].get_histo_array_int(counter)
                        p = [ peakheight, float(npeaks[counter]), peakwidth, 
                              np.mean(histo[self.firstbin:self.lastbin]), 
                              np.mean(histo[np3s:nbin])]
                        y = histo[:nbin]
                        ##############
                        # guess values
                        ##############
                        mm = muprompt()
                        mm._init_(x,y)
                        m = Minuit(mm,a=p[0],x0=p[1],dx=p[2],ak1=p[3],ak2=p[4])
                        # m.values = p
                        m.errors = (p[0]/100,p[1]/100,0.01,p[3]/100,p[4]/100)
                        # m.print_level = 1 if mprint else 0
                        m.migrad()
                        A,X0,Dx,Ak1,Ak2 = m.values
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
                            prompt_fit_text[counter] = self.ax_counters[divmod(counter,3)].text(x_text,y_text,'Det #{}\nt0={}bin\n$\delta$t0={:.2f}'.format(counter+1,x0.round().astype(int)[counter],x0[counter]-x0.round().astype(int)[counter]))
 
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
                    nt0.children[0].value = ' '.join(map(str,self.nt0.astype(int)))
                    dt0.children[0].value = ' '.join(map('{:.2f}'.format,self.dt0))
                                                   # refresh, they may be slightly adjusted by the fit
                                                   
                elif self.filespecs[1].value=='mdu': # PSI HIFI
                #############################
                # very rough guess of histo start bin
                # then 
                # fit a step
                ############################# 
                    ncounters = self._the_runs_[0][0].get_numberHisto_int()
                    npeaks = []
                    a = 0.5*np.ones(ncounters)
                    b = 30*np.ones(ncounters)
                    dn = 5*np.ones(ncounters)
                    for counter in range(ncounters):
                        histo = np.empty(self._the_runs_[0][0].get_histo_array_int(counter).shape)
                        for k in range(len(self._the_runs_[0])): # may add runs
                            histo += self._the_runs_[0][k].get_histo_array_int(counter)
                        npeakguess = scanms(histo,100) # simple search for a step pattern
                        if npeakguess>0:
                            npeaks.append(npeakguess)
                        elif counter != 0:
                            self.console('**** WARNING: step in hifi detector {} not found'.format(counter))
                            self.console('     set to arbitrary bin 20000')
                            npeaks.append(20000)
                        else:
                            npeaks.append(np.where(histo==histo.max())[0][0])
                        # print(npeaks)
                        ###############
                        # now fit it
                        ###############
                        if counter != 0:
                            n2 = npeaks[counter] + second_plateau # counter dependent bin interval
                            n1 = npeaks[counter] + first_plateau
                            x = np.arange(n1,n2+1,dtype=int) # n2-n1+1 bins from n1 to n2 included for plotting
                            y = histo[n1:n2+1]
                            # next will try likelihood
                            c = cost.LeastSquares(x,y,1,step)
                            m = Minuit(c,a=a[counter],n=npeaks[counter],dn=dn[counter],b=b[counter])
                            # m.errors(1.,10.,1.)
                            m.migrad()
                            a[counter],n,dn[counter],b[counter] = m.values
                            if m.valid:                               
                                npeaks.pop()
                                npeaks.append(n)
                            else:
                                self.console('     step fit not converged for detector {}'.format(counter))
                    x0 = np.array(npeaks).astype(int)
                    self.lastbin = x0.min() - self.prepostpk[0].value # final bin for background average 

                 ############################
                 # just show where this is and save parameters
                 ############################
                    if mplot:

                        ###################
                        #  set figure, axes (8  real counters, python number 1 2 3 4 5 6 7 8
                        ###################
                        self.fig_counters,self.ax_counters = P.subplots(3,3,figsize=(9.5,9.5),dpi=dpi)
                        self.fig_counters.canvas.set_window_title('HIFI Start histo guess')
                        screen_x, screen_y = P.get_current_fig_manager().window.wm_maxsize() # screen size in pixels
                        y_maxinch = float(screen_y)/dpi #  maximum y size in inches

                        prompt_fit_text = [None]*ncounters   
                        n2 = x0.max() + second_plateau # counter independent bin interval
                        n1 = x0.min() + first_plateau
                        for counter in range(ncounters):
                            self.ax_counters[divmod(counter,3)].cla()
                            # self.ax_counters[divmod(counter,3)].axis('off')
                            histo = np.empty(self._the_runs_[0][0].get_histo_array_int(counter).shape)
                            for k in range(len(self._the_runs_[0])): # may add runs
                                histo += self._the_runs_[0][k].get_histo_array_int(counter)
                            x = np.arange(n1,n2+1,dtype=int) # n2-n1+1 bins from n1 to n2 included for plotting
                            y = histo[n1:n2+1]
                            x3 = np.arange(n1,n2)
                            # with self.t0plot_container:
            #                       if self.first_t0plot:
                            # self.ax_counters[divmod(counter,3)].cla()
                            self.ax_counters[divmod(counter,3)].plot(x,y,'.')
                            self.ax_counters[divmod(counter,3)].plot(x,step(x,a[counter],npeaks[counter],dn[counter],b[counter]),'r-')
                            x_text,y_text = npeaks[counter]+10,0.8*histo.max()
                            prompt_fit_text[counter] = self.ax_counters[divmod(counter,3)].text(x_text,y_text,'Det #{}\nt0={}bin'.format(counter+1,x0[counter]))
                    self.nt0 = x0 # bin of peak, nd.array of shape run.get_numberHisto_int() 
                    self.dt0 = np.zeros(x0.shape) # fraction of bin, nd.array of shape run.get_numberHisto_int() 
                    
                    self.nt0_run = self.create_rundict()
                    nt0.children[0].value = ' '.join(map(str,self.nt0))
                    dt0.children[0].value = ' '.join(map('{:.2f}'.format,self.dt0))
                                                   # refresh, they may be slightly adjusted by the fit


                elif self.filespecs[1].value=='nxs': # ISIS
                    if mplot:

                        ###################
                        #  set figure, axes 
                        ###################
                        self.fig_counters,self.ax_counters = P.subplots(figsize=(5,4),dpi=dpi)
                        self.fig_counters.canvas.set_window_title('Edge t0 fit')
                        screen_x, screen_y = P.get_current_fig_manager().window.wm_maxsize() # screen size in pixels
                        y_maxinch = float(screen_y)/dpi #  maximum y size in inches

                    histo = np.empty(self._the_runs_[0][0].get_histo_array_int(0).shape[0])
                    for counter in range(self._the_runs_[0][0].get_numberHisto_int()):
                        for k in range(len(self._the_runs_[0])): # may add runs
                            histo += self._the_runs_[0][k].get_histo_array_int(counter)
                    error = np.sqrt(histo)
                    error[np.where(error==0)]=1
                    dh = histo[1:]-histo[:-1]
                    kt0 = np.where(dh==dh.max())[0] # [0]
                    musbin = float(self.nsbin.value)/1e3
                    t0 = kt0*musbin
                    N = histo[int(kt0)+10]*self.TauMu_mus
                    D = 0.080
                    n1 = 0
                    n2 = 101
                    t = musbin*np.linspace(n1,n2-1,n2)
                    mm = muedge()
                    mm._init_(t,histo[n1:n2])
                    m = Minuit(mm,t00=t0,N=N,D=D)
                    m.errors=(t0/100,N/100,0.8)
                    m.print_level = 1 if mprint else 0                   
#                    if mplot:    
#                        # with self.t0plot_container:
#        #                       if self.first_t0plot:
#                        self.ax_counters.plot(t,histo[n1:n2],'.')
#                        self.ax_counters.plot(t,mm.f(t,t0,N,D))
#                        x_text,y_text = t[int(2*n2/3)],0.2*max(histo[n1:n2])
#                        self.ax_counters.text(x_text,y_text,'t0 = {:.1f} mus'.format(t0))
#                        self.fig_counters.canvas.manager.window.tkraise()
#                        P.draw()
                    m.migrad()
                    t0,N,D = m.values
                    if mplot:    
                        # with self.t0plot_container:
        #                       if self.first_t0plot:
                        self.ax_counters.plot(t,histo[n1:n2],'.')
                        self.ax_counters.plot(t,mm.f(t,t0,N,D))
                        x_text,y_text = t[int(2*n2/3)],0.2*max(histo[n1:n2])
                        self.ax_counters.text(x_text,y_text,'t0 = {:.1f} mus'.format(t0))
                    self.nt0 = np.array([t0/float(self.nsbin.value)]).round().astype(int) # bin of peak, nd.array of shape run.get_numberHisto_int() 
                    self.dt0 = np.array(t0-self.nt0) # fraction of bin, nd.array of shape run.get_numberHisto_int()                     
                    self.nt0_run = self.create_rundict()
                    nt0.children[0].value = ' '.join(map(str,self.nt0.astype(int)))
                    dt0.children[0].value = ' '.join(map('{:.2f}'.format,self.dt0))
            if mplot:                        
                self.fig_counters.canvas.manager.window.tkraise()
                P.draw()


        def save_run_list(b):
            """
            when user presses this setup tab button mugui
            saves ascii file .log with run list in data directory
            """
            import os
            from glob import glob
            from mujpy.musr2py.musr2py import musr2py as muload
            # from psibin import MuSR_td_PSI_bin as muload
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
            _paths_content = [ self.paths[k].value for k in range(3) ] # should be 3 ('data','analysis','tlog')
            _filespecs_content = [ self.filespecs[k].value for k in range(2) ] # should be 2 ('fileprefix','extension')
            _filespecs_options = self.filespecs[0].options
            setup_dict = {}
            if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu': # PSI
                facility = 'PSI'
                _prepostpk = [self.prepostpk[k].value for k in range(2)] # pre-prompt, post-prompt bins
                names = ['_paths_content','_filespecs_content','_filespecs_options',
                         '_prepostpk','self.nt0','self.dt0','self.lastbin','self.nt0_run'] # keys
            elif  self.filespecs[1].value=='nxs': # ISIS:
                facility = 'ISIS'
                names = ['_paths_content','_filespecs_content','_filespecs_options',
                         'self.nt0','self.dt0','self.nt0_run'] # keys
            #self.console('ext = {}, names = {}'.format(self.filespecs[1],names))
            for k,key in enumerate(names):
               setup_dict[names[k]] = eval(key) # key:value
            with open(path,'wb') as f:
                pickle.dump(setup_dict, f) # according to __getstate__()
            self.console('Saved {} {}'.format(facility,os.path.join(self.__startuppath__,'mujpy_setup.pkl'))) 

        def select_paths(b):
            '''
            select paths where data, logs, [tlogs] are stored
            '''
            from mujpy.aux.aux import path_dialog, shorten
            import os
            
            tpath = ''
            # shorten provides './subpath' if path is subpath of startuppath
            thispath = self.paths[0].value if os.path.isdir(self.paths[0].value) else self.__startuppath__
            dpath = path_dialog(thispath,'Select path to data files') # do not shorten (maybe no common root)
            thispath = self.paths[1].value if os.path.isdir(self.paths[1].value) else self.__startuppath__
            lpath = path_dialog(thispath,'Select path to log files') # do not shorten (maybe no common root)
            if self.filespecs[1]=='bin':
                tpath = shorten(path_dialog(thispath,'Select path to T mon files'),self.__startuppath__)
            if dpath == '' or lpath == ' ':
                self.console('Paths not selected')
                return
            else:
                self.paths[0].value, self.paths[1].value = dpath, lpath
                if tpath: self.paths[2].value = tpath
            
        from ipywidgets import HBox, Layout, VBox, Text, Textarea, IntText, Checkbox, Button, Output, Accordion, Dropdown   
        from numpy import array 

        # setup for things that have to be set initially (paths, t0, etc.)
        # the tab is self.mainwindow.children[2], a VBox 
        # containing a setup_box of three HBoxes: path, and t0plot 
        # path is made of a firstcolumn, paths, and a secondcolumns, filespecs, children of setup_box[0]
        # agt0 is made of three 


        width = '100%'
        height = 'auto'

        paths_content = ['data','analysis','tlog']   

 # ---------------firstrow

        load_button = Button(description='load setup',
                            tooltip = 'from file mujpy_setup.pkl',
                            layout=Layout(width='11%')) # 11%
        load_button.style.button_color = self.button_color
        load_button.on_click(load_setup)
        #          012345678901234567890123456789012345678901234567890123456789
        balloon = 'data path\n(./ is is the folder\nwhere you launched jupyter)'
        datapath = Text(description='data',
                        description_tooltip=balloon,
                        disabled=True,
                        layout=Layout(width='22%'))     # 33%
        datapath.style.description_width='15%'

        balloon = 'logs path\n(./ is is the folder\nwhere you launched jupyter)'
        logpath = Text(description='logs',
                       description_tooltip=balloon,
                       disabled=True,
                       layout=Layout(width='22%'))      # 55%
        logpath.style.description_width='15%'
        
        balloon = 'tlog path\n(./ is is the folder\nwhere you launched jupyter)'
        tlogpath = Text(description='tlog',
                        description_tooltip=balloon,
                        disabled=True,
                        layout=Layout(width='22%'))     # 77%
        tlogpath.style.description_width='15%'
        
        self.paths = [datapath, logpath, tlogpath]
        
        select_paths_button = Button(description='Select paths',
                                     tooltip='GUI to select where\ndat, logs [and tlogs]\nare stored',
                                     layout={'width':'11%'}) # 88%
        select_paths_button.on_click(select_paths)
        select_paths_button.style.button_color = self.button_color
        
        run_list_button = Button(description='Run list',
                            tooltip = 'Writes titles of runs\nin list file\nin data path',
                            layout=Layout(width='11%'))    # 99%
        run_list_button.on_click(save_run_list)
        run_list_button.style.button_color = self.button_color 
        
             

        firstrow = HBox([load_button, datapath, logpath, tlogpath,
                         select_paths_button, run_list_button],
                                                 layout=Layout(width=width))


#----------------- secondrow

        save_button = Button(description='save setup',
                            tooltip = 'in file mujpy_setup.pkl\nloaded at startup',
                            layout=Layout(width='11%')) # 11%
        save_button.style.button_color = self.button_color
        save_button.on_click(save_setup)

        extension = Dropdown(options=['bin','mdu','nxs'],description='File extension',
                        description_tooltip='bin=PSI\nmdu=PSI HIFI\nnxs=ISIS',
                        layout=Layout(width='22%'))     # 33%
        extension.style.description_width='45%'
        extension.observe(on_extension_changed)

        fileprefix = Dropdown(options=['deltat_tdc_gps_',
                               'deltat_tdc_ltf_',
                               'deltat_tdc_dol_',
                               'deltat_tdc_gpd_',
                               'tdc_hifi_'],
                               description='prefix',
                               value='deltat_tdc_gps_',
                               description_tooltip='e.g. deltat_tdc_gps_',
                               disabled=False,
                               layout=Layout(width='22%'))     # 55%
        fileprefix.style.description_width='30%'
        
        self.filespecs = [fileprefix, extension] 

        prepeak = IntText(description='prepeak',
                          description_tooltip='bins before prompt\nfor left plateau',
                          value = 7,
                          layout=Layout(width='16.5%'), 
                          continuous_update=False)     # 72%
        prepeak.style.description_width='50%'
        
        postpeak = IntText(description='postpeak',
                          description_tooltip='bins after prompt\nfor right plateau',
                          value = 7,
                          layout=Layout(width='16.5%'),
                          continuous_update=False)     # 89%
        postpeak.style.description_width='50%'

        self.prepostpk = [prepeak, postpeak]

        fit_button = Button(description='t0 fit',
                            layout=Layout(width='11%')) # 100%
        if self.filespecs[1].value=='bins': # PSI
            fit_button.tooltip = 't0 is fitted as prompt centre\n(peak on step function)'
        else:
            fit_button.tooltip = 't0 is fitted on sum\nof all detectors\nas centre of sigmoid\n convoluted with pion decay'
        fit_button.on_click(on_prompt_fit_click)
        fit_button.style.button_color = self.button_color

        load_setup([]) # invokes load_setup to load mujpy_setup.pkl 
        # provides for non existent file!

        secondrow = HBox([save_button, self.filespecs[1], self.filespecs[0], 
                          prepeak, postpeak,
                          fit_button],layout=Layout(width=width))

#--------------thirdrow

        nt0 = Accordion(children=[Textarea(description='t0 [bin]',disabled=True,
                                    layout={'width':'99%'})],
                        layout={'width':'35%','height':height})
        nt0.children[0].style.description_width='15%'
        nt0.set_title(0,'t0 bin')
        nt0.selected_index =  None
        
        dt0 = Accordion(children=[Textarea(description='dt0 [ns]',disabled=True,
                                    layout={'width':'99%'})],
                        layout={'width':'35%','height':height})
        dt0.children[0].style.description_width='20%'
        dt0.set_title(0,'fraction from bin center [ns]')
        dt0.selected_index =  None

        nt0.children[0].value = ','.join(str(l)+'\n' *(n % 4 == 3) for n, l in enumerate(self.nt0))
        dt0.children[0].value = ','.join(map('{:.2f}'.format,self.dt0))

        thirdrow = HBox([nt0,dt0],layout=Layout(width=width,height=height))
        if not self.nt0_run:
            self.console('WARNING: you must fix t0 = 0, please do a t0 fit')
        setup_hbox = [VBox([firstrow,secondrow,thirdrow],
        layout=Layout(width=width))]
        self.mainwindow.children[2].children = setup_hbox # first tab (setup)
        self.mainwindow.children[2].layout=Layout(border='2px solid dodgerblue',
                                           width=width,height='auto')

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
                n1 = offset_bin+self.nt0[0]
                for counter in range(run.get_numberHisto_int()):
                    if self.filespecs[1].value=='bin' or self.filespecs[1].value=='mdu':
                        n1 = offset_bin+self.nt0[counter] 
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


#    def introspect(self):# mark for deletion: use MuJPy.__dict__ directly below GUI
#        '''
#        print updated attributes of the class mugui 
#        after each fit in file "mugui.attributes.txt" 
#        in self.__startuppath__
#        '''
#        import os
#        import pprint
#        # from ipywidgets import VBox, HBox, Image, Text, Textarea, Layout, Button, IntText, Checkbox, Output, Accordion, Dropdown, FloatText, Tab
#        # trick to avoid printing the large mujpy log image binary file 
#        image = self.__dict__['gui'].children[0].children[0].value
#        self.__dict__['gui'].children[0].children[0].value=b''
#        with open(os.path.join(self.__startuppath__,"mugui.attributes.txt"),'w') as f:
#            pprint.pprint('**************************************************',f)
#            pprint.pprint('*               Mugui attribute list:            *',f)
#            pprint.pprint('**************************************************',f)
#            pprint.pprint(self.__dict__,f)
#        self.__dict__['gui'].children[0].children[0].value=image

