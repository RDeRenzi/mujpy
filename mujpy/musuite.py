class suite(object):
    """
    A suite class, file read through musr2py (bin, mdu), muroot2py (root, old and MusrRoot), muisis2py (nxs)
    
      muroot2py provides methods equivalent to all MuSR_td_PSI_bin(), same syntax (but data are always numpy arrays) 
      
      Different Instruments/Facilities require different methods for setting t0, specifically:
            HAL not known
            GPS, GPD, DOLLY, FLAME, LTF, determine a prompt position by identifying the count maximum
                                         and fit a prompt peak over a given interval around the peak
            LEM trust the instrument tof calibration encoded in the header
            NXS fit the edge ISIS function, common to all instrument hence to the nxs file spec
      Instrument identification as of 2024 is BY FILE SPEC: mdu ->   HAL (hifi), bin -> GPS, FLAME, GPD, DOLLY, nxs -> ISIS
      However, from 2023 on root can be any instrument 
      Instrument identification must be performed early in the method. Root contains it. Bin can be oly gps, gpd, flame, dolly, ltf so it's ok.
      mdu is only hifi. nxs is also fine (and it may contain  
    self._init_ reads input fixed variables
        console, runlist, datafile 
        grp_calib, offset 
        (console(string) must print string somewhere)
        t0 parameters are fixed depending on bin mdu root spec
    
    self.groups = grp.calib contains a list of dictionaries for forward, backward groups and their alpha values
        where groups may be in shorthand notation
        self.alpha value used for normal fit 
        as opposed to calibration fits where alpha is a fit parameter
    self.grouping is a list of dictionaries for forward, backward groups and their alpha values
        where groups as np array of 0 based indices of detectors (tools get_grouping does the translation)
    imports musr2py or muisis2py and loads it as instance for each data set
    which can be an individual run or the sum of several runs

    self._the_runs_ is a list of lists of musr2py instances
    self._the_runs_[k]  is a list of musr2py instances to be added

    invokes prompfit(self) calls and calculates t = 0 parameters. modify to accommodate root

    self.timebase returns time, always 1d array
    self.single_for_back_counts(runs,grouping) acts on runs
            yforw, ybackw are sums, 
            background_forw, background_back their backgrounds (PSI) or zero (ISIS) 
            yfm, ybm are <(yi -background_i)*exp(t/tau)> i = forw, backw (works also for ISIS)
            allow on the fly asymmetry in Minuit, with alpha fit parameter
    self.asymmetry_single calculates time, asymmetry and asymmetry error, 1d arrays,  with given alpha
            invoking single_for_back_counts
    self.asymmetry_multirun calculates time (1d), asymmetry  and asymmetry error, 2d arrays, with given alpha
            invokes self.single_for_back_counts(runs,grouping) for runs in range(len(self._the_runs_))
    self.etc methods for suite, multi suite,  global suite and global multi suite
            to be done
    Notes   NB mujpy is python3 only
        Output ends up in notebook, below cell, and 
        sys.__stdout__ is <_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>
        sys.__stderr__ is <_io.TextIOWrapper name='<stderr>' mode='w' encoding='utf-8'>
    """ 

    def __init__(self, datafile , runlist , grp_calib , offset , startuppath, console = 'print',mplot=False):
                 
        """
        __init__ for musuite class

        * inputs: 
            the suite_input_file a dict containing
                    datafile,   containing the full path 
                    runlist, run number or list of runs
                    grp_calib, grouping and alpha paramter dictionary, see below
                    offset, first good bin
                    startuppath, path where mudash is lauched                      
        * grp_calib is a list of dictionaries (minimum one)
          {'forward':stringfw,'backward':stringbw,'alpha':alpha}
          strings are translated into np.arrays of histograms by tools get_grouping
        * upon initialization automatically
                checks input, stores paths
                load_runs(): stores data load instance(s) in self._the_runs_
                store_groups() stores list of dicts in self.grouping, each containing
                               'forward' and 'backward' lists of detector indices
                promptfit(): determines t0
                timebase(): stores self.time (1d)  
        if console = 'print' self.console will exec print(string)
        if mudashed calls suite with console = 'self.log', a method of mudashed  
              self.console will exec self.log(string), i.e. write on board output                           
        """

        from mujpy.tools.tools import derun
        import json
        import os
        from mujpy import __file__ as MuJPyName
        from mujpy._version import __version__
        self.__version__ = __version__
        # print('__init__ now ...')
        self.firstbin = 5 # common to all data file specs, to evaluate background for subtraction 
        self.loadfirst = False  # non initialized trivialities in the present instance of self.suite
        self._the_runs_ = [] # initialized to empty
#        with open(suite_input_file,"r") as f:
#            suite_input_json = f.read()
#            suite_input = json.loads(suite_input_json)
        self.console_method = console # either 'print' for scripts or 'self.log' for mudashed calls, eg. from jupyter-lab or ipython
        try:
            if self.console_method == 'print': # not a call from a mudash instance
                self.console('******************* SUITE *********************')
                # self.console('debug')
            elif isinstance(console,object):
                self.console('**** SUITE: using caller log ****')
            else:
                self.console('')
        except ValueError as e:
            print('Except: {}\nno suite.console!'.format(e))
            return  # with no console error message, it means no console
        # determine number of runs, filenames etc.
            #######################
            # decode the run string
            #######################           
        self.runs, errormessage = derun(runlist) # self.runs is a list of lists of run numbers (string)
        if errormessage is not None: # derun error
            self.console('Run syntax error? {}. You typed: {}'.format(errormessage,runlist))
            return  # with console error message

        self.nruns = len(self.runs) # vanished on 23 July 2026 at 17:00 ?!??!  reinserted here
        if os.path.isfile(datafile):
            self.datafile = datafile
            if self.datafile[-4:]=='root':
                self.thermo = 0 # sample thermometer is always a single value list
            elif self.datafile[-3:]=='bin': 
                self.thermo = 1 # sample thermometer is 1 on gps (check or adapt to other instruments)
                self.prepostpk = [50, 50]
            elif self.datafile[-3:]=='mdu':
                self.thermo = 1 # sample thermometer is 1 on gps (check for hifi)
                self.prepostpk = [70, 70]
            self.loadfirst = True
        else:
            self.console('* File {} not found'.format(datafile))
#            self.console('* CHECK YOUR ACCESS (e.g. klog)')
            self.console('**************************************')
            return None  # with console error message
        self.__startuppath__ = os.getcwd()
        # implement new folder policy, see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.Dashboard
        self.__path__ = os.path.dirname(MuJPyName) # mujpy path
        # REMEMBER! ../ does not work if it points outside the directory where jupyterlab is launched from
        pre_mujpy_path = os.path.dirname(MuJPyName)[:os.path.dirname(MuJPyName).rfind(os.sep)]
        self.__grppath__ = pre_mujpy_path+os.sep+'groups'+os.sep
        self.__templatepath__ = pre_mujpy_path+os.sep+'templates'+os.sep
        if not os.path.exists(self.__startuppath__+os.sep+'fit'+os.sep):
            os.mkdir(self.__startuppath__+os.sep+'fit'+os.sep)
        self.__fitpath__ = self.__startuppath__+os.sep+'fit'+os.sep
        if not os.path.exists(self.__startuppath__+os.sep+'csv'+os.sep):
            os.mkdir(self.__startuppath__+os.sep+'csv'+os.sep)
        self.__csvpath__ = self.__startuppath__+os.sep+'csv'+os.sep
        if not os.path.exists(self.__startuppath__+os.sep+'cache'+os.sep):
            os.mkdir(self.__startuppath__+os.sep+'cache'+os.sep)
        self.__cachepath__ = self.__startuppath__+os.sep+'cache'+os.sep
        # reverse of tools get_grouping is tools get_group
        self.offset = int(offset) # offset belongs to suite, that needs it for asymmetries
        if self.load_runs(): #           load data instances in self._the_runs_
            self.groups = grp_calib
        # grp_calib is a list of dictionaries, one per group, with keys 'forward,'backward','alpha'
        # where groups may be in shorthand notation
            self.grouping = [] # reproduce the same with arrays of histogram numbers, done by tools get_grouping
            self.store_groups() #        in self.grouping        self.promptfit(mplot)   #    to be done: make switch for ISIS
            self.promptfit(mplot)   #    to be done: switch for ISIS broken and root LEM wip
            self.timebase()
        # self.console('... end of initialize suite')
        else:
            self.console('*********************** suite exits without loading data ****************************')
    
    def console(self,string):
        """
        writes string to initiated self.console_method

           when suite invoked without console, self.console_method defaults top 'print'
           when invoked from mudashes, self.console_method = mudashed.log
        """

        #print('debug string = '+string)
        #print('debug console_method = {}'.format(self.console_method))
        if self.console_method == 'print': 
            print(string)
        else:
            exec('self.console_method'+'("'+string+'")')
 
    def add_runs(self,k):
        """
        Tries to load one or more runs to be added together

        by means of murs2py. 
        self.runs is a list of strings containing integer run numbers 
        Returns -1 and quits if musr2py or muisis2py complain, 0 otherwise
        """

        from mujpy.muroot2py.muroot2py import muroot2py as rootload 
        from musr2py import MuSR_td_PSI_bin as psiload
        from mujpy.muisis2py.muisis2py import munxs2py as isisload
        # muisis2py has the same methods as musr2py
        from mujpy.tools.tools import get_datafilename, get_title, short_path
        
        def instrument(run):
            instrument_set = set(['gps','gpd','dolly','flame','ltf','hifi'])
            return list(set(run.Filename().split('_')).intersection(instrument_set))[0].upper()

        read_ok = True
        runadd = []
        for j,run in enumerate(self.runs[k]): # run is a single run number
            path_and_filename =  get_datafilename(self.datafile,run) 
            # print('path_and_filename in suite: {}'.format(path_and_filename))
            if self.datafile[-4:]=='root':
                try:
                    runadd.append(rootload(path_and_filename))
                    self._the_instrument_ = runadd[-1].get_instrument()
                    self._the_facility_ = 'PSI'
                except:
                    read_ok = False
            elif self.datafile[-3:]=='bin' or self.datafile[-3:]=='mdu':
                runadd.append(psiload())
                if runadd[j].read(path_and_filename) != 0:
                    read_ok = False # THE RUN DATA FILE IS LOADED HERE
                else:
                    self._the_instrument_ = instrument(runadd[j])
                    self._the_facility_ = 'PSI'

            elif self.datafile[-3:]=='nxs': 
                try: 
                    #self.console('{} ns'.format(isisload(path_and_filename,'r').get_binWidth_ns()))
                    runadd.append(isisload(path_and_filename,'r'))  # adds a new instance of isisload
                    #self.console('path = {}'.format(path_and_filename))
                    self._the_instrument_ = runadd[-1].get_instrument()
                    self._the_facility_ = 'ISIS'
                except:
                    read_ok = False
            if read_ok==True:
                self._the_runs_.append(runadd) # 
                    
                self.console('{} loaded:'.format(short_path(path_and_filename,self.__startuppath__)))
                self.console('Run {}: {}'.format(run,get_title(self._the_runs_[-1][0])))
#                import sys
#                self.console('sys.__stdout__ is {}, sys.__stderr__ is {}'.format(sys.__stdout__,sys.__stderr__))

                if k>0:
                    ok = [self._the_runs_[k][0].get_numberHisto_int() == 
                          self._the_runs_[0][0].get_numberHisto_int(),
                          self._the_runs_[k][0].get_binWidth_ns() == 
                          self._the_runs_[0][0].get_binWidth_ns()]
                    if not all(ok): 
                        self.console(r'\nFile {} has wrong histoNumber or binWidth'.format(
                                  get_datafilename(self.datafile,self._the_runs_[k][0].get_runNumber_int())))
                        return -1  # this leaves the first run of the suite

                return True
            else:
                # print(path_and_filename)
                self.console(r'\nRun file: {} does not exist'.format(path_and_filename))
           	    #self.console('            if reading from afs check that klog is not expired!')
                return False

    def load_runs(self):
        """
        load musr2py. muroot2py or muisis2py instances stored as a list of lists [[runs to add], ...]

        self._the_runs_[0][0] a single run, or the first of a suite 
        self._the_runs_[k][0] the k-th run of a run suite
        
        Invoked after creating a suite instance, typically as
            the_suite = suite('log/input.suite') # the_suite implements the class suite according to input.suite
            if the_suite.load_runs():            # this and the following two statements load data
                if the_suite.store_groups():     #                                       define groups
                    the_suite.promptfit(mplot=False)    #                                fit t0 = 0
        """

        read_ok = True
        for k,runs_add in enumerate(self.runs):#  runs_add can be a list of run numbers (string) to add
            read_ok = read_ok and self.add_runs(k)                
            # print('on_loads_change, inside loop, runs {}'.format(self._the_runs_))
        
        return read_ok # False with console error message in add_runs

    def check_group(self,group):
        """
        rough check that this is a group of existing detectors
        """

        # numberHisto_int is the number of physics detectors
        # RedGreen mode has more than one period (in Isis parlance) for different stimuli ON or OFF 
        # for this purpose nexus detector counts have three indices: period, detector, bin
        # MusrRoot instead has offsets, typically [0,20,40,80] or [0,10,20,30], reflected in the label Histo No
        # whereas the list self._histo.counts() has two indices: histogram, bin, and histogram numbers are contiguous
        # e.g. on lem for numberHisto_int = 8 and offsets [0,10] indices 0,...,7 are the original histos and 8,...,15 are PPC
        numberHisto = self._the_runs_[0][0].get_numberHisto_int()
        periods = 1
        if 'get_RedGreen_offsets' in self._the_runs_[0][0].__dir__(): # RedGreen may be present
            periods = len(self._the_runs_[0][0].get_RedGreen_offsets()) # are there RedGreen copies (periods>1)?
        if 'get_beamline' in self._the_runs_[0][0].__dir__(): # PSI only
            numberHisto = numberHisto*periods           
        #print('group = {} nH =  {} p = {}'.format(group,numberHisto,periods))
        return (group>=0).all()*(group<numberHisto).all()

    def store_groups(self):
        """
        from self.groups dashboard shorthand dict to self.grouping dict alpha, lists of histogram numbers  
        """

        from mujpy.tools.tools import get_grouping
        for k,group in enumerate(self.groups):
            fgroup, bgroup, alpha = get_grouping(group['forward']), get_grouping(group['backward']), group['alpha']
            if alpha>0 and self.check_group(fgroup) and self.check_group(bgroup): # checks legal grpcalib_file
            #a,b = self.check_group(fgroup), self.check_group(bgroup)
            # print('fwd = {} bwd = {}'.format(a,b))
            #if alpha>0 and a and b:
                if k==0: self.grouping=[]
                self.grouping.append({'forward':fgroup, 'backward':bgroup, 'alpha':alpha})
                # fgroup bgroup are two np.arrays of integers
            else:
                self.console('forw {}, backw {}, alpha {:.2f}, Nhisto = {}'.format(fgroup,bgroup,alpha,
                                self._the_runs_[0][0].get_numberHisto_int()))
                self.console('Groups calibration file corrupted')
                return False
        return True

              
    def t_value_error(self,k):
        """
        calculates T and eT values also for runs to be added
        
        sillily, but it works also for single run 
        """

        from numpy import sqrt

        m = len(self._the_runs_[k]) # number of added runs
        weight = [float(sum(self._the_runs_[k][j].get_histo_vector(counter,1))) for counter in range(self._the_runs_[0][0].get_numberHisto_int()) for j in range(m)]
        if sum(weight)>0:
            weight = [w/sum(weight) for k,w in enumerate(weight)]
            t_value = sum([self._the_runs_[k][j].get_temperatures_vector()[self.thermo]
                              *weight[j] for j in range(m)])
            t_error = sqrt(sum([(self._the_runs_[k][j].get_devTemperatures_vector([self.thermo])
                              *weight[j])**2 for j in range(m)])) if self.datafile[-3:] in ['oot','bin','mdu'] else 0
        else:
            t_value, t_error = 0, 0
        return t_value, t_error

    def promptfit(self,mplot, mprint = False):
        """
        indentifies t0 and stores self.nt0 array (all intruments) [ISIS not yet]

        t0 prompt fit method for PSI gps, gpd, dolly, ltf, flame identified by muroot2py.get_instrument() or by musr2py.readbin
        t0 bin value for PSI lem by muroot2py.get_t0_int()
        t0 edge method for ISIS (all), now broken
        t0 guess method for PSI hifi 

        refactored for run addition and
        suite of runs
        WARNING: this module is tenporarily for PSI only, included root        
        """

        from numpy import array, where, arange, zeros, zeros_like, mean, ones, sqrt, linspace
        from iminuit import Minuit, cost
        import matplotlib.pyplot as P
        from mujpy.mucomponents.muprompt import muprompt
        from mujpy.mucomponents.muedge import muedge
        from mujpy.tools.tools import TauMu_mus, scanms, step
        from mujpy.tools.plot import set_fig 
    
        if mplot:  # setup figure window if a plot of the prompt peak or edge fits is required
            font = {'family' : 'Ubuntu','size'   : 8}
            P.rc('font', **font)
            dpi = 100. # conventional screen dpi
            num = 0 # unique window number
            if self.datafile[-3:] == 'bin': 
                nrow, ncol = 2,3
                kwargs = {'figsize':(7.5,5),'dpi':dpi}
                title = 'Prompts t0 fit'
                prompt_fit_text = [None]*self._the_runs_[0][0].get_numberHisto_int()
            elif self.datafile[-3:] =='mdu': # PSI HIFI        
                nrow, ncol = 3,3
                ###################
                #  set figure, axes (8  real counters, python number 1 2 3 4 5 6 7 8
                ###################
                # fig_counters,ax_counters = P.subplots(3,3,figsize=(9.5,9.5),dpi=dpi)
                kwargs = {'figsize':(9.5,9.5),'dpi':dpi}
                title = 'HIFI start histo guess'
            elif self.datafile[-3:]=='nxs': # ISIS
                nrow, ncol = 3,3
                ###################
                #  set figure, axes 
                ###################
                kwargs = {'figsize':(5,4),'dpi':dpi}
                title = 'Edge t0 fit'
            fig_counters,ax_counters = set_fig(num,nrow,ncol,title,**kwargs)
            fig_counters.canvas.set_window_title(title)
        
        if self._the_instrument_  == 'LEM':
            # t0 calibration encoded in file, override only by setting self.nt0 by hand
            self.nt0 = array([int(self._the_runs_[0][0].get_t0_int(k)) for k in range(self._the_runs_[0][0].get_numberHisto_int())])
            self.dt0 = zeros_like(self.nt0) # fraction of bin, ignored with tdc resolution of 195 pds  
            self.lastbin = self.nt0.min() # unused for lem

            
        elif self._the_instrument_ in ['GPS','GPD','DOLLY','FLAME','LTF']: # self.datafile[-3:] == 'bin': # PSI gps, flame, ltf, dolly, gpd
            second_plateau = 100
            peakheight = 100000.
            peakwidth = 1.
            ###################################################
            # fit a peak with different left and right plateaus
            ###################################################

            #############################
            # guess prompt peak positions
            ############################# 
            npeaks = []
            for counter in range(self._the_runs_[0][0].get_numberHisto_int()):
                histo = zeros(self._the_runs_[0][0].get_histoLength_bin())
                for k in range(len(self._the_runs_[0])): # may add runs
                    histo += array(self._the_runs_[0][k].get_histo_vector(counter,1))
                binpeak = where(histo==histo.max())[0][0]
                npeaks.append(binpeak)
            npeaks = array(npeaks)

            ###############
            # right plateau
            ###############
            nbin =  int(max(npeaks) + second_plateau) # this sets a counter dependent second plateau bin interval
            x = arange(0,nbin,dtype=int) # nbin bins from 0 to nbin-1
            self.lastbin, np3s = npeaks.min() - self.prepostpk[0], int(npeaks.max() + self.prepostpk[1])
            # final bin for background average, first bin for right plateau estimate (last is nbin)

            x0 = zeros(self._the_runs_[0][0].get_numberHisto_int()) # for center of peaks

            if mplot:
                for counter in range(self._the_runs_[0][0].get_numberHisto_int(),sum(ax_counters.shape)):
                    ax_counters[divmod(counter,3)].cla()
                    ax_counters[divmod(counter,3)].axis('off')

            for counter in range(self._the_runs_[0][0].get_numberHisto_int()):
                # prepare for muprompt fit
                histo = zeros(self._the_runs_[0][0].get_histoLength_bin())
                for k in range(len(self._the_runs_[0])): # may add runs
                    histo += self._the_runs_[0][k].get_histo_vector(counter,1)
                p = [ peakheight, float(npeaks[counter]), peakwidth, 
                      mean(histo[self.firstbin:self.lastbin]), 
                      mean(histo[np3s:nbin])]
                y = histo[:nbin]
                ##############
                # guess values
                ##############
                mm = muprompt()
                mm._init_(x,y)
                mm.errordef = Minuit.LEAST_SQUARES
                m = Minuit(mm,a=p[0],x0=p[1],dx=p[2],ak1=p[3],ak2=p[4])
                # m.values = p
                m.errors = (p[0]/100,p[1]/100,0.01,p[3]/100,p[4]/100)
                m.migrad()
                A,X0,Dx,Ak1,Ak2 = m.values
                x0[counter] = X0 # store float peak bin position (fractional)  

                if mplot:    # do plot
                    n1 = npeaks[counter]-50
                    n2 = npeaks[counter]+50
                    x3 = arange(n1,n2,1./10.)
                    ax_counters[divmod(counter,3)].cla()
                    ax_counters[divmod(counter,3)].plot(x[n1:n2],y[n1:n2],'.')
                    ax_counters[divmod(counter,3)].plot(x3,mm.f(x3,A,X0,Dx,Ak1,Ak2))
                    x_text,y_text = npeaks[counter]+10,0.8*max(y)
                    prompt_fit_text[counter] = ax_counters[
                                                  divmod(counter,3)].text(x_text,y_text,
                     r'Det #{}\nt0={}bin\n$\delta$t0={:.2f}'.format(counter+1,
                     x0.round().astype(int)[counter],x0[counter]-x0.round().astype(int)[counter]))

                ##############################################################################
                # Simple cases:                                                              #
                # 1) Assume the prompt is entirely in bin nt0.                               #
                #   (python convention, the bin index is 0,...,n,...                         #
                # The content of bin nt0 will be the t=0 value for this case and dt0 = 0.    #
                # The center of bin nt0 will correspond to time t = 0,                       #
                #         time = (n-nt0 + mufit.offset + mufit.dt0)*mufit.binWidth_ns/1000.  #
                # 2) Assume the prompt is equally distributed between n and n+1.             #
                #    Then nt0 = n and dt0 = 0.5, the same formula applies                    #
                # 3) Assume the prompt is 0.45 in n and 0.55 in n+1.                         #
                #    Then nt0 = n+1 and dt0 = -0.45, the same formula applies.               #
                ##############################################################################

                # these three are the sets of parameters used by other methods
            self.nt0 = x0.round().astype(int) # bin of peak, nd.array of shape run.get_numberHisto_int() 
            self.dt0 = x0-self.nt0 # fraction of bin, nd.array of shape run.get_numberHisto_int() 
            self.lastbin = self.nt0.min() - self.prepostpk[0] # nd.array of shape run.get_numberHisto_int() 

        elif self._the_instrument_ == 'HIFI' and self._the_facility_ == 'PSI': # self.datafile[-3:] =='mdu': # PSI HIFI
            first_plateau = - 500
            second_plateau = 1500
            #############################
            # very rough guess of histo start bin
            # then 
            # fit a step
            ############################# 
            ncounters = self._the_runs_[0][0].get_numberHisto_int()
            npeaks = []
            a = 0.5*ones(ncounters)
            b = 30*ones(ncounters)
            dn = 5*ones(ncounters)
            for counter in range(ncounters):
                histo = zeros(self._the_runs_[0][0].get_histoLength_bin())
                for k in range(len(self._the_runs_[0])): # may add runs
                    histo += self._the_runs_[0][k].get_histo_vector(counter,1)
                npeakguess = scanms(histo,100) # simple search for a step pattern
                if npeakguess>0:
                    npeaks.append(npeakguess)
                elif counter != 0:
                    self.console('**** WARNING: step in hifi detector {} not found'.format(counter))
                    self.console('     set to arbitrary bin 20000')
                    npeaks.append(20000)
                else:
                    npeaks.append(where(histo==histo.max())[0][0])
                ###############
                # now fit it
                ###############
                if counter != 0:
                    n2 = npeaks[counter] + second_plateau # counter dependent bin interval
                    n1 = npeaks[counter] + first_plateau
                    x = arange(n1,n2+1,dtype=int) # n2-n1+1 bins from n1 to n2 included for plotting
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
                        self.console('****   step fit not converged for detector {}'.format(counter))
            x0 = array(npeaks).astype(int)
            self.lastbin = x0.min() - self.prepostpk[0].value # final bin for background average 

            ############################
            # just show where this is and save parameters
            ############################
            if mplot:     # do plot
                prompt_fit_text = [None]*ncounters   
                n2 = x0.max() + second_plateau # counter independent bin interval
                n1 = x0.min() + first_plateau
                for counter in range(ncounters):
                    ax_counters[divmod(counter,3)].cla()
                    # ax_counters[divmod(counter,3)].axis('off')
                    histo = zeros(self._the_runs_[0][0].get_histoLength_bin())
                    for k in range(len(self._the_runs_[0])): # may add runs
                        histo += self._the_runs_[0][k].get_histo_vector(counter,1)
                    x = arange(n1,n2+1,dtype=int) # n2-n1+1 bins from n1 to n2 included for plotting
                    y = histo[n1:n2+1]
                    x3 = arange(n1,n2)
                    ax_counters[divmod(counter,3)].plot(x,y,'.')
                    ax_counters[divmod(counter3)].plot(x,
                                step(x,a[counter],npeaks[counter],dn[counter],b[counter]),'r-')
                    x_text,y_text = npeaks[counter]+10,0.8*histo.max()
                    prompt_fit_text[counter] = ax_counters[divmod(counter,3)].text(x_text,
                                  y_text,'Det #{}\nt0={}bin'.format(counter+1,x0[counter]))
            self.nt0 = x0 # bin of peak, nd.array of shape run.get_numberHisto_int() 
            self.dt0 = zeros(x0.shape) # fraction of bin, nd.array of shape run.get_numberHisto_int()

        elif self._the_facility_ == 'ISIS': # self.datafile[-3:]=='nxs': # ISIS
            histo = zeros(self._the_runs_[0][0].get_histoLength_bin())
            for counter in range(self._the_runs_[0][0].get_numberHisto_int()):
                for k in range(len(self._the_runs_[0])): # may add runs
                    histo += self._the_runs_[0][k].get_histo_vector(counter,1)
            error = sqrt(histo)
            error[where(error==0)]=1
            dh = histo[1:]-histo[:-1]
            kt0 = where(dh==dh.max())[0] # [0]
            musbin = float(self.nsbin.value)/1e3
            t0 = kt0*musbin
            N = histo[int(kt0)+10]*TauMu_mus()
            D = 0.080
            n1 = 0
            n2 = 101
            t = musbin*linspace(n1,n2-1,n2)
            mm = muedge()
            mm._init_(t,histo[n1:n2])
            m = Minuit(mm,t00=t0,N=N,D=D)
            m.errors=(t0/100,N/100,0.8)
            m.print_level = 1 if mprint else 0                   
            m.migrad()
            t0,N,D = m.values
            
            
            if mplot:    # do plot
                ax_counters.plot(t,histo[n1:n2],'.')
                ax_counters.plot(t,mm.f(t,t0,N,D))
                x_text,y_text = t[int(2*n2/3)],0.2*max(histo[n1:n2])
                ax_counters.text(x_text,y_text,'t0 = {:.1f} mus'.format(t0))
            self.nt0 = array([t0/float(self.nsbin.value)]).round().astype(int) # bin of peak, 
                                             # nd.array of shape run.get_numberHisto_int() 
            self.dt0 = array(t0-self.nt0) # fraction of bin, in ns


        if mplot:   # show results                  
            fig_counters.canvas.manager.window.tkraise()
            P.draw()
#            self.console('Succesfully completed prompt Minuit fit, check plots')
        else:
            pass
#            self.console('Succesfully completed prompt Minuit fit, check nt0, dt0 ')
#        self.console('****************END OF SUITE*****************')

##########################
# ASYMMETRY
##########################
    def mean_dt0(self):
        """
        PSI only, calculates average of dt0 over histograms in self.grouping       
        """

        from numpy import append, mean
        histos = append(self.grouping[0]['forward'],self.grouping[0]['backward'])
        if len(self.grouping)>1:
        # self.grouping[:]['forward'] or ['backward'] are np.arrays
            for k in range(len(self.grouping)): # find the mean of dt0 over all histos of all groups
                histos = append(histos,self.grouping[k]['forward'])
                histos = append(histos,self.grouping[k]['backward'])
                # a list of np.arrays
        return mean(self.dt0[histos])

    def timebase(self):
        """
        generates self.time

        * initializes self histoLength 
        * fills self.time. 1D numpy array
        * all histogram selects common time
        * PSI has different t0 per histogram
        * and must standardize to a common length 
        
        # Time definition for center of bin n: 
        #          time = (n - self.nt0 + self.offset + self.dt0)*binWidth_ns/1000.
        # 1) Assume the prompt is entirely in bin self.nt0. (python convention, the bin index is 0,...,n,... 
        # The content of bin self.nt0 will be the t=0 value for this case and self.dt0 = 0.
        # The center of bin self.nt0 will correspond to time t = 0
        # 2) Assume the prompt is equally distributed between n and n+1. 
        #    Then self.nt0 = n and self.dt0 = 0.5, the same formula applies
        # 3) Assume the prompt is 0.45 in n and 0.55 in n+1. 
        #    Then self.nt0 = n+1 and self.dt0 = -0.45, the same formula applies.
        """

        import numpy as np
   
        ##################################################################################################
        # self histoLength = self._the_runs_[0][0].get_histoLength_bin() - self.nt0.max() - self.offset
        # needed to set a
        ##################################################################################################

        time_bins = np.arange(self.offset,self._the_runs_[0][0].get_histoLength_bin() - 
                               self.nt0.max())   # 1D np.array
        binwidth_mus = self._the_runs_[0][0].get_binWidth_ns()/1000.
        self.histoLength = self._the_runs_[0][0].get_histoLength_bin() - self.nt0.max() - self.offset

        self.time = time_bins*binwidth_mus 
        
        if self.datafile[-3:]=='bin' or self.datafile[-3:]=='mdu': # PSI
            self.time += self.mean_dt0()*binwidth_mus # mean dt0 correction (fraction of a bin, probaby immaterial)

            
#        # Asymmetry as 
##        denominator = yfm_ + alpha*self._ybm_)*exp(-x/TauMu_mus()) # f+b normalization count
##        A = (yf - alpha*yb - (bf - alpha*bb)) / (yfm+alpha*ybm)*exp(-x/TauMu_mus()) 
##        errexp = sqrt(yf + alpha**2*yb) # equivalent f+b counts
##        errexp[where(errexp==0)] = 1  #   set to 1 the minimum error for zero equivalent f+b counts
##        self._e_ = errexp / denominator 

#        # Error eval box:
#        # Nf(t), Nb(t) are row counts from the two groupings, forward and backward
#        # A(t) = y  with background corrected counts Nfc(t) = Nf(t) - bf = yf, 
#        #                                            Nbc(t) = Nb(t) - bb = yb
#        # errors eA(T) with renormalized counts
#        #                                 Nf(t) = cf,
#        #                                 Nb(t) = cb

#        #############################################################
#        #  ISIS)         Error evaluation, no backgrounds:          #
#        # Brewers trick to avoid double error propagation:          #
#        # the denominator is evaluated as an average                #
#        #       A = [yf(t) - alpha yb(t)]/d          with           #
#        #    d = (<yf e^t/tau> + alpha <yb e^t/tau>)e-t/tau         #
#        # yf = sum_{i in f) Ni        yb = sum_{i in b} Ni          #
#        # ef^2 = yf                   eb^2 = yb                     #
#        #          eA = sqrt(yf + alpha^2 yb)/d                     #
#        #-----------------------------------------------------------#
#        # PSI)          With background                             #
#        # evaluate bf, bb before prompt for yf, yb respectively     #
#        #                                                           #
#        #  A = [yf-bf - alpha(yb-bb)]/[yf-bf + alpha(yb-bb)]        #
#        #           ef^2, eb^2 are the same                         #
#        # d = [<(yf-bf)e^t/tau)> + alpha <(yb-bb)e^t/tau>] e^-t/tau #
#        #   =   [<yfbe>     +    alpha     <ybbe>] e^-t/tau         #
#        #                                                           #
#        #     A = [yf - alpha yb - (bf - alpha bb)]/d               #
#        #                                                           #
#        #    eA = sqrt(yf + alpha^2 yb)/d                           #
#        #-----------------------------------------------------------#
#        # if alpha is a paramter                                    #
#        # compute and return yf, yb, bf, bb, <yfbe>, <ybbe>         #
#        # mumodel must compute   d,  A, eA                          #
#        #############################################################
#        # for ISIS the PSI formula work with bb and bf zero         #
#        #############################################################
#        # rebin eA works for ISIS, must be corrected for PSI        #
#        # yfm, ybm depend on binning   

    def single_for_back_counts(self,runs,grouping):
        """
        basic method for single group single run count arrays

        * input: 
        *         runs, runs to add 
        *         grouping, {'forward':[3],'backward':[4]} for 3-4 
        * output:
        *    for PSI, with b_j = mean(data_j before muon arrival)
        *         yfc, ybc     = sum_{j in for or back} (data_j - b_j)
        *    for ISIS (b_j = 0)
        *                      = sum_{j in for or back} data_j) 
        *         eyfc, eyfc   = sqrt(sum_{j in for or back}(data_j + (b_j + p(0,b)/n)))
        *                 with error e = sqrt(N+(b+p(0,b))/n), 
        *                 where p(0,b) is probability for 0 count, 
        *                              either
        *                 p_0b = normal probability , ~Poisson std=sqrt(b)
        *                      = exp(-b/2)/sqrt(2*pi*b), or
        *                 P_0b = true Poisson = exp(-b) 
        * Method used both in self.asymmetry_single normal fits
        * and in calib (alpha minuit parameter)
        #        NOTE: in calib fits alpha is a fit parameter with an error unknown while minimizing
        #              at the minimum e_alpha is typically 2e-3 on PSI calib runs
        #              Neglected here, could recalculate chi_square at minimum including error from non corrected chi_square
        #        eA = 2*alpha/(self._yfc_ - alpha*self._ybc_)**2 * sqrt((self._ybc_*self._eyfc_)**2 + (self._yfc_*self._eybc_)**2) 
        *
        * all output objects are 1D numpy arrays
        """

        from numpy import zeros, array, mean, exp, where, sqrt
        from mujpy.tools.tools import TauMu_mus

        filespec = self.datafile[-3:] # 'bin', 'mdu', "oot" or 'nsx'
        if self.loadfirst:
            
#            n1 = self.nt0[0] + self.offset # ISIS and lem PSI root
#            n2 = n1 + self.histoLength # ISIS and lem PSI root
#            print('musuite single_for_back_counts debug: n1 {}, n2 {}, self.histoLength {}'.format(n1,n2,self.histoLength))
            yfc, eyfc = zeros(self.histoLength), zeros(self.histoLength) # counts 
            ybc, eybc = zeros(self.histoLength), zeros(self.histoLength) # counts 
            back_f, back_b = 0., 0.             # background initialized
            n = self.lastbin-self.firstbin+1    # n of bins in back estimation 
                           
            for j, run in enumerate(runs): # Add runs
                #print(run)
                for counter in grouping['forward']: 
                
                    histo = array(run.get_histo_vector(counter,1)) # counter data array in forw group
#                    if array(where(histo==0)).size:
#                        print('musuite single_for_back_counts debug: run {} counter fwd {} contains {}  zeros'.format(run.get_runNumber_int(),counter,array(where(histo==0)).size))
                    if filespec =='bin' or filespec=='mdu' or filespec=='oot': # PSI, counter specific range                  
                        n1 = self.nt0[counter] + self.offset # first good bin
                        n2 = n1 + self.histoLength           # last good bin
                        back_f += mean(histo[self.firstbin:self.lastbin]) # before muon arrival
                    yfc += histo[n1:n2] # counts 
                p_0b =  exp(-back_f) if back_f > 0 else 0. # true Poisson probability for zero counts
                eyfc = sqrt(yfc+(back_f+p_0b)/n) # error accounting also for zero counts
                yfc -= back_f # counts corrected for background, valid also for ISIS
                yfc[where(yfc<0)] = 0
                        
                for counter in grouping['backward']: # first good bin, last good bin from data attay start
                
                    histo = array(run.get_histo_vector(counter,1)) # counter data array in back group
#                    if array(where(histo==0)).size:
#                        print('musuite single_for_back_counts debug: run {} counter bkw {} contains {}  zeros'.format(run.get_runNumber_int(),counter,array(where(histo==0)).size))
                    if filespec=='bin' or filespec=='mdu' or filespec=='oot': #  PSI, counter specific range  
                        n1 = self.nt0[counter] + self.offset
                        n2 = n1 + self.histoLength 
                        back_b += mean(histo[self.firstbin:self.lastbin])  # before muon arrival
                    ybc += histo[n1:n2]
                p_0b =  exp(-back_b) if back_b > 0 else 0. # true Poisson probability for zero counts
                eybc = sqrt(ybc+(back_b+p_0b)/n)
                ybc -= back_b
                ybc[where(ybc<0)] = 0
             
            return yfc, ybc, eyfc, eybc
        else:
            return None, None, None, None
        #-----------------------------------------------------------#
        # compute and return yfc, ybc, eyfc, eybc                   #
        # asymmetry computes A, eA if alpha is constant             #
        # mucomponents computes A, eA if alpha is a fit parameter   #
        # eb = sqrt(sum_i y_i + p(0)*n)/n = sqrt((b+p(o))/n)        #
        # errors eyc = sqrt(eN**2+eb**2) = sqrt(N+(b+p(0))/n) <---- #
        #-----------------------------------------------------------#
    def slice_for_back_counts(self,krun,kgroup):
        """
        basic methods for count array slices
        
        input:
            krun = index in range(len(self._the_runs_)) or -1 for [:len(self._the_runs_)]
            kgroup = index in range(len(self.groups)) or -1 for [:len(self.groups)]
        output
            the corresponding yf, yb, eyf, eyb slices
        cases:
            slice==suite: A1, A21, C1, C2 slice is (-1,-1) but it is not really needed 
            slice<suite:  A20 1d (0,kgroup), B1 1d (krun,0), B20 1d (krun,kgroup), B21 2d (krun,-1)  
        """

        from numpy import array, vstack
        if krun==-1 and kgroup==-1: # 3d
            for krun,run in enumerate(self._the_runs_):
                for kgroup,grouping in enumerate(self.grouping):
                    yforw, ybackw, ey_forw, ey_backw = self.single_for_back_counts(run,grouping)
                    if kgroup:
                        yf, yb = vstack((yf,array([yforw]))), vstack((yb,array([ybackw])))
                        eyf, eyb = vstack((eyf,array([ey_forw]))), vstack((eyb,array([ey_backw])))
                    else: # kgroup = 0
                        yf, yb = array([yforw]), array([ybackw])
                        eyf, eyb = array([ey_forw]), array([ey_backw])
                if krun: 
                    yfr, ybr  = vstack((yfr,array([yf]))), vstack((ybr,array([yb])))
                    eyfr, eybr = vstack((eyfr,array([eyf]))),vstack((eybr,array([eyb])))
                else: # krun = 0
                    yfr, ybr = array([yf]),array([yb]) # 3rd dimension runs
                    eyfr, eybr = array([eyf]),array([eyb])

        elif krun==-1: # 2d
            for krun,run in enumerate(self._the_runs_):
                yforw, ybackw, ey_forw, ey_backw = self.single_for_back_counts(run,self.grouping[kgroup])
                if krun:
                    yfr, ybr = vstack((yfr,array([yforw]))), vstack((ybr,array([ybackw])))
                    eyfr, eybr = vstack((eyfr,array([ey_forw]))), vstack((eybr,array([ey_backw])))
                else: # krun = 0
                    yfr, ybr = array([yforw]), array([ybackw])
                    eyfr, eybr = array([ey_forw]), array([ey_backw])
        elif kgroup==-1: # 2d
            for kgroup,grouping in enumerate(self.grouping):
                yforw, ybackw, ey_forw, ey_backw = self.single_for_back_counts(self._the_runs_[krun],grouping)
                if kgroup:
                    yfr, ybr = vstack((yfr,array([yforw]))), vstack((ybr,array([ybackw])))
                    eyfr, eybr = vstack((eyfr,array([ey_forw]))), vstack((eybr,array([ey_backw])))
                else: # kgroup=0
                    yfr, ybr = array([yforw]), array([ybackw])
                    eyfr, eybr = array([ey_forw]), array([ey_backw])
        else: # 1d, only bins
            yfr, ybr, eyfr, eybr = self.single_for_back_counts(self._the_runs_[krun],self.grouping[kgroup])
        return yfr, ybr, eyfr, eybr

    def asymmetry_single(self,the_run,kgroup):
        """
        basic method for plain asymmetry and errors arrays

        input:
            the_run = list containing the instance[s] of the run[s to be added]
            k = index of self.grouping, a list of dicts 
                self.grouping[k]['forward'] and ['backward'] (py-index, i.e "counter 1-Backw" is  0)
                containing the respective lists of detectors
        * run instances from musr2py/muisis2py  (psi/isis load routine) 
        *
        outputs: 
            # can be A1 fit, but is also invoked by all others
            asymmetry and asymmetry error (1d)
        """

        from numpy import exp, sqrt, where

        if self.loadfirst:
            # print(the_run)
            alpha = self.grouping[kgroup]['alpha']            
            yf, yb, eyf, eyb = self.single_for_back_counts(the_run,self.grouping[kgroup])
            
            denominator = (yf + alpha*yb)
            denominator[where(denominator==0)] = 1  
            asymm = (yf - alpha*yb) / denominator  # 1d array
            
            #   ey_i in fcn sum_i ((y_i-y_th)/ey_i)**2 cannot be zero, but errexp can
            asyme = 2*alpha/(denominator)**2*sqrt((yb*eyf)**2+(yf*eyb)**2) 
            asyme[where(asyme==0)] = 1 # could be zero, yielding inf chi_square 
            return asymm, asyme
        else:
            return None, None
        #############################################################
        #  ISIS)         Error evaluation, no backgrounds:          #
        #-----------------------------------------------------------#
        # PSI)          With background                             #
        # evaluate yfc, ybc, eyfc, eybc (corrected for background)  #
        #                                                           #
        #  A = [yfc - alpha ybc)]/[yfc + alpha ybc]                 #
        #  eA = 2 alpha/(yfc + alpha*ybc)**2*sqrt((ybc*eyfc)**2 +   #
        #                                         (yfcn*eybc)**2)   #
        #############################################################
        # for ISIS the PSI formula works with zero background       #
        #############################################################

    def asymmetry_slice(self,krun,kgroup):
        """
        asymmetry for mufit (-1 index lingo), can yield slice <= suite

        input:
            krun = index in range(len(self._the_runs_)) or -1 for [:len(self._the_runs_)]
            kgroup = index in range(len(self.groups)) or -1 for [:len(self.groups)]
        output
            the corresponding asymmetry slice 
        cases:
            slice==suite: A1  1d (0,0), A21 2d (0,-1), C1 2d (-1,0), C2 3d (-1,-1)
                        the same is obtained with asymmetry_multirun_multigroup(self)
            slice<suite:  A20 1d (0,kgroup), B1 1d (krun,0), B20 1d (krun,kgroup), B21 2d (krun,-1)
        """

        from numpy import array, vstack
        if krun==-1 and kgroup==-1: # 3d
            for krun,run in enumerate(self._the_runs_):
                for kgroup in range(len(self.groups)):
                    a,e = self.asymmetry_single(run,kgroup)
                    if kgroup:
                        asy,ase = vstack((asy,a)), vstack((ase,e)) # groups are vstacked in 2nd dimension
                    else: # kgroup = 0
                        asy,ase = a,e # 1d dimension bins
                if krun:
                    asymm, asyme  = vstack((asymm,array([asy]))), vstack((asyme,array([ase])))
                else: # krun=0
                    asymm, asyme = array([asy]),array([ase]) # 3rd dimension runs
        elif krun==-1: # 2d
            for krun,run in enumerate(self._the_runs_):
                a,e = self.asymmetry_single(run,kgroup)
                if krun:
                    asymm, asyme  = vstack((asymm,a)), vstack((asyme,e)) # runs are vstacked in 2nd dimension
                else: # krun=0
                    asymm, asyme = a,e # 1nd dimension bins
        elif kgroup==-1:
            for kgroup in range(len(self.groups)):
                a,e = self.asymmetry_single(self._the_runs_[krun],kgroup)
                if kgroup:
                    asymm,asyme = vstack((asymm,a)), vstack((asyme,e)) # groups are vstacked in 2nd dimension
                else: # kgroup = 0
                    asymm,asyme = a,e# 1d dimension bins
        else: # 1d, only bins
            asymm,asyme = self.asymmetry_single(self._the_runs_[krun],kgroup)
        return asymm, asyme
     
    def single(self):
        """
        True if len(self.runs)==1

        output:
            True if there is a single run (fit type A)
            False if there are many runs (fit types B and C)
        """

        try:
            test = len(self._the_runs_)==1
            
        except:
            self.console('Warning: data are not available: access expired?')
        return test
            
    def multi_groups(self):
        """
        False if A1 or B1 or C1, True otherwise

        output:
            True if more groups (fits A2, B2, C2)
            False if just one group (fits A1, B1, C1)
        """
        # print('multi_group suite debug: self.grouping {} len {}'.format(self.grouping,len(self.grouping)))
        return len(self.grouping)>1     
        
    def multirun(self):
        """
        True if B1, B20, B21, C1, C2, False otherwise
        """

        return len(self._the_runs_)>1

    def scan(self):
        """
        returns suite scan string, "T[K] " , "B[mT]", "[deg]", "#   ", False if self.single()

        output
            False if single
            'B[mT]' if it's a B scan
            'T[K] ' if it's a T scan
            '[deg]' if it's an angle scan 
            '#    ' if it's another scan
        """

        if self.single(): return False
        elif [run[0].get_temp() for run in self._the_runs_]!=[self._the_runs_[0][0].get_temp()]*len(self._the_runs_): return 'T[K] '
        elif [run[0].get_field() for run in self._the_runs_]!=[self._the_runs_[0][0].get_field()]*len(self._the_runs_): return 'B[mT]'
        elif [run[0].get_orient() for run in self._the_runs_]!=[self._the_runs_[0][0].get_orient()]*len(self._the_runs_): return '[deg]'
        else: return '#    '

