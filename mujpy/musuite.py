class suite(object):
    '''
    A suite class, 
    self._init_ reads input from JSON suite_input_file, that can be edited, including
        runlist, datatemplatefile, logpath, 
        groups_calibration file, offset and a console method 
        (console(string) must print string somewhere)
        t0 parameters are fixed depending on bin or mdu spec
    
    loads grp.calib JSON file containing a list of dictionaries for groups and their alpha values
        (future global fits with different asymmetry functions, to start with only single)
        self.alpha value for normal fit (in future can be list for multi fits)
        as opposed to fits producing calibration files

    imports musr2py or muisis2py and loads it as instance for each data set
    which can be an individual run or the sum of several runs

    self._the_runs_ is a list of lists of musr2py instances
    self._the_runs_[k]  is a list of musr2py instances to be added

    invokes prompfit(self) calls and calculates t = 0 parameters

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
    ''' 
    def __init__(self,suite_input_file = 'input.suite',mplot=True):
        # print('__init__ now ...')
                 
        self.firstbin = 0
        self._initialise_suite(suite_input_file,mplot)       

    def _initialise_suite(self,suite_input_file,mplot):
        '''
        * Initiates an instance of suite, 
        * inputs: 
            the suite_input_file a dict containing
                    runlist, datatemplatefile, logpath, offset
                    the group[s]_calibration file, console method  
        * the group[s] calibration file is a JSON file containing
          a list of dictionaries (minimum one)
          {'forward':stringfw,'backward':stringbw,'alpha':alpha}
          strings are translated into np.arrays of histograms by aux get_grouping
        * upon initialization automatically
                checks input
                load_runs(): stores data load instance(s) in self._the_runs_
                store_groups() stores list of dicts in self.grouping, each containing
                               'forward' and 'backward' lists of detector indices
                promptfit(): determines t0
                timebase(): stores self.time (1d)                               
        '''
        from mujpy.aux.aux import derun
        import json
        import os
        # read suite_input
        self.loadfirst = False
        self._the_runs_ = [] # initialized 
        self.grouping = [] # reproduce the same with arrays of histogram numbers
        with open(suite_input_file,"r") as f:
            suite_input_json = f.read()
            suite_input = json.loads(suite_input_json)
        self.console_method = suite_input["console"]
        try:
            self.console('******************* SUITE *********************')
        except:
            return False # with no console error message, it means no console
        # determine number of runs, filenames etc.
            #######################
            # decode the run string
            #######################           
        self.runs, errormessage = derun(suite_input['runlist']) # self.runs is a list of lists of run numbers (string)
        if errormessage is not None: # derun error
            self.console('Run syntax error: {}. You typed: {}'.format(errormessage,suite_input['runlist']))
            return False  # with console error message
        if os.path.isfile(suite_input['datafile']):
            self.datafile = suite_input['datafile']
            if self.datafile[-3:]=='bin': 
                self.thermo = 1 # sample thermometer is 1 on gps (check or adapt to other instruments)
                self.prepostpk = [7, 7]
            elif self.datafile[-3:]=='mdu':
                self.thermo = 1 # sample thermometer is 1 on gps (check for hifi)
                self.prepostpk = [70, 70]
            self.loadfirst = True
        else:
            self.console('File {} not found'.format(suite_input['datafile']))
            return False  # with console error message
        if os.path.isdir(suite_input['logpath']):
            self.logpath = suite_input['logpath']
        else:
            self.console('Log path {} not found'.format(suite_input['logpath']))   
            return False # with console error message
        if os.path.isfile(self.logpath+suite_input['groups calibration']):
            grpcalib_file = self.logpath+suite_input['groups calibration']
        else:
            self.console('Calibration file {} not found'.format(suite_input['groups calibration']))  
            return False  # with console error message             
        # load groups calibration file
        with open(grpcalib_file,"r") as f:
            grpin = f.read()
            try:
                self.groups = json.loads(grpin)
            except:
                self.console('Groups calibration file {} is not proper JSON'.format(suite_input['groups calibration']))  
                return False  #  with console error message                             
            # groups_in is in dashboard jargon
        self.logpath = suite_input_file[:suite_input_file.rfind('/')+1]
        self.offset = suite_input['offset']
        # initialization automatically implies also loading data and groups
        #try:
        self.load_runs() #               load data instances in self._the_runs_
        self.store_groups() #        in self.grouping
        self.promptfit(mplot)   # make switch for ISIS
        self.timebase()
        #except Exception as e:
        #    self.console('*** suite: something went wrong, ask RDR ***')
        #    self.console('* {} *'.format(e))
        # self.console('... end of initialize suite')

    def add_runs(self,k):
        '''
        Tries to load one or more runs to be added together
        by means of murs2py. 
        self.runs is a list of strings containing integer run numbers 
        Returns -1 and quits if musr2py or muisis2py complain, 0 otherwise
        '''
        
        from mujpy.musr2py.musr2py import musr2py as psiload
        from mujpy.muisis2py.muisis2py import muisis2py as isisload
        # muisis2py has the same methods as musr2py
        from mujpy.aux.aux import get_datafilename, get_title

        read_ok = True
        runadd = []
        for j,run in enumerate(self.runs[k]): # run is a single run number
            path_and_filename =  get_datafilename(self.datafile,run)
            if self.datafile[-3:]=='bin' or self.datafile[-3:]=='mdu':
                runadd.append(psiload())
                if runadd[j].read(path_and_filename) != 0:
                    read_ok = False # THE RUN DATA FILE IS LOADED HERE
            elif self.datafile[-3:]=='nxs': 
                try: 
                    #self.console('{} ns'.format(isisload(path_and_filename,'r').get_binWidth_ns()))
                    runadd.append(isisload(path_and_filename,'r'))  # adds a new instance of isisload
                    #self.console('path = {}'.format(path_and_filename))
                except:
                    read_ok = False
            if read_ok==True:
                self._the_runs_.append(runadd) # 
               	self.console('Run {} loaded'.format(path_and_filename))
               	self.console(' {}'.format(get_title(self._the_runs_[-1][0])))
                if k>0:
                    ok = [self._the_runs_[k][0].get_numberHisto_int() == 
                          self._the_runs_[0][0].get_numberHisto_int(),
                          self._the_runs_[k][0].get_binWidth_ns() == 
                          self._the_runs_[0][0].get_binWidth_ns()]
                    if not all(ok): 
                        self._the_runs_=[self._the_runs_[0]] # leave just the first one
                        self.console('\nFile {} has wrong histoNumber or binWidth'.format(
                                  get_datafilename(self.datafile,self._the_runs_[k][0].get_runNumber_int())))
                        return -1  # this leaves the first run of the suite

                return True
            else:
           	    self.console('\nRun file: {} does not exist'.format(path_and_filename))
           	    self.console('            if reading from afs check that klog is not expired!')
           	    return False

    def load_runs(self):
        '''
        load musr2py or muisis2py instances
        stored as a list of lists 
        self._the_runs_[0][0] a single run, or the first of a suite 
        self._the_runs_[k][0] the k-th run of a run suite
        
        Invoked after creating a suite instance, typically as
            the_suite = suite('log/input.suite') # the_suite implements the class suite according to input.suite
            if the_suite.load_runs():            # this and the following two statements load data
                if the_suite.store_groups():     #                                       define groups
                    the_suite.promptfit(mplot=False)    #                                fit t0 = 0
        '''
        read_ok = True
        for k,runs_add in enumerate(self.runs):#  runs_add can be a list of run numbers (string) to add
            read_ok = read_ok and self.add_runs(k)                
            # print('on_loads_change, inside loop, runs {}'.format(self._the_runs_))
        
        return read_ok # False with console error message in add_runs

    def store_groups(self):
        '''
        reads groups dictionary in dashboard shortnote
        and appends lists of histogram numbers, alphas to self.grouping 
        '''
        from mujpy.aux.aux import get_grouping
        for k,group in enumerate(self.groups):
            fgroup, bgroup, alpha = get_grouping(group['forward']), get_grouping(group['backward']),       group['alpha']
            ok = 1
            ok *= alpha>0*(fgroup>=0).all()*(fgroup<self._the_runs_[0][0].get_numberHisto_int()).all()
            ok *= (bgroup>=0).all()*(bgroup<self._the_runs_[0][0].get_numberHisto_int()).all()
            if ok: # checks for errors in grpcalib_file
                self.grouping.append({'forward':fgroup, 'backward':bgroup, 'alpha':alpha})
                # fgroup bgroup are two np.arrays of integers
            else:
                self.console('forw {}, backw {}, alpha {:.2f}, Nhisto = {}'.format(fgroup,bgroup,alpha,
                                self._the_runs_[0][0].get_numberHisto_int()))
                self.console('Groups calibration file corrupted')
                return False
        return True

    def console(self,string):
        exec(self.console_method+'("'+string+'")')
               
    def t_value_error(self,k):
        '''
        calculates T and eT values also for runs to be added
        sillily, but it works also for single run 
        '''
        from numpy import sqrt

        m = len(self._the_runs_[k]) # number of added runs
        weight = [float(sum(self._the_runs_[k][j].get_histo_array_int(2))) for j in range(m)]
        if sum(weight)>0:
            weight = [w/sum(weight) for k,w in enumerate(weight)]
            t_value = sum([self._the_runs_[k][j].get_temperatures_vector()[self.thermo]
                              *weight[j] for j in range(m)])
            t_error = sqrt(sum([(self._the_runs_[k][j].get_devTemperatures_vector([self.thermo])
                              *weight[j])**2 for j in range(m)])) if self.datafile[-3:]=='bin' or self.datafile[-3:]=='mdu' else 0
        else:
            t_value, t_error = 0, 0
        return t_value, t_error


    def promptfit(self,mplot, mprint = False):
        '''
        launches t0 prompts fit::

            fits peak positions 
            prints migrad results
            plots prompts and their fit (if plot checked, mprint not implemented)
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
        from mujpy.aux.aux import TauMu_mus, scanms, step, set_fig 
    
        if mplot:  # setup figure window
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
            elif self.filespecs[1].value=='nxs': # ISIS
                nrow, ncol = 3,3
                ###################
                #  set figure, axes 
                ###################
                kwargs = {'figsize':(5,4),'dpi':dpi}
                title = 'Edge t0 fit'
            fig_counters,ax_counters = set_fig(num,nrow,ncol,title,**kwargs)
            fig_counters.canvas.set_window_title(title)
                
        if self.datafile[-3:] == 'bin':  # PSI gps, flame, dolly, gpd
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
            self.lastbin, np3s = npeaks.min() - self.prepostpk[0], npeaks.max() + self.prepostpk[1] 
            # final bin for background average, first bin for right plateau estimate (last is nbin)

            x0 = np.zeros(self._the_runs_[0][0].get_numberHisto_int()) # for center of peaks
                   
            if mplot:
                for counter in range(self._the_runs_[0][0].get_numberHisto_int(),sum(ax_counters.shape)):
                    ax_counters[divmod(counter,3)].cla()
                    ax_counters[divmod(counter,3)].axis('off')

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
                    x3 = np.arange(n1,n2,1./10.)
                    ax_counters[divmod(counter,3)].cla()
                    ax_counters[divmod(counter,3)].plot(x[n1:n2],y[n1:n2],'.')
                    ax_counters[divmod(counter,3)].plot(x3,mm.f(x3,A,X0,Dx,Ak1,Ak2))
                    x_text,y_text = npeaks[counter]+10,0.8*max(y)
                    prompt_fit_text[counter] = ax_counters[
                                                  divmod(counter,3)].text(x_text,y_text,
                     'Det #{}\nt0={}bin\n$\delta$t0={:.2f}'.format(counter+1,
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

        elif self.datafile[-3:] =='mdu': # PSI HIFI
            first_plateau = - 500
            second_plateau = 1500
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
                        self.console('****   step fit not converged for detector {}'.format(counter))
            x0 = np.array(npeaks).astype(int)
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
                    histo = np.empty(self._the_runs_[0][0].get_histo_array_int(counter).shape)
                    for k in range(len(self._the_runs_[0])): # may add runs
                        histo += self._the_runs_[0][k].get_histo_array_int(counter)
                    x = np.arange(n1,n2+1,dtype=int) # n2-n1+1 bins from n1 to n2 included for plotting
                    y = histo[n1:n2+1]
                    x3 = np.arange(n1,n2)
                    ax_counters[divmod(counter,3)].plot(x,y,'.')
                    ax_counters[divmod(counter3)].plot(x,
                                step(x,a[counter],npeaks[counter],dn[counter],b[counter]),'r-')
                    x_text,y_text = npeaks[counter]+10,0.8*histo.max()
                    prompt_fit_text[counter] = ax_counters[divmod(counter,3)].text(x_text,
                                  y_text,'Det #{}\nt0={}bin'.format(counter+1,x0[counter]))
            self.nt0 = x0 # bin of peak, nd.array of shape run.get_numberHisto_int() 
            self.dt0 = np.zeros(x0.shape) # fraction of bin, nd.array of shape run.get_numberHisto_int()

        elif self.filespecs[1].value=='nxs': # ISIS
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
            N = histo[int(kt0)+10]*TauMu_mus()
            D = 0.080
            n1 = 0
            n2 = 101
            t = musbin*np.linspace(n1,n2-1,n2)
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
            self.nt0 = np.array([t0/float(self.nsbin.value)]).round().astype(int) # bin of peak, 
                                             # nd.array of shape run.get_numberHisto_int() 
            self.dt0 = np.array(t0-self.nt0) # fraction of bin, in ns


        if mplot:   # show results                  
            fig_counters.canvas.manager.window.tkraise()
            P.draw()
            self.console('Succesfully completed prompt Minuit fit, check plots')
        else:
            self.console('Succesfully completed prompt Minuit fit, check nt0, dt0 ')
        self.console('****************END OF SUITE*****************')

##########################
# ASYMMETRY
##########################
    def mean_dt0(self):
        '''
        PSI only
        calculates average of dt0 over histograms in self.grouping       
        '''
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
        * initializes self histoLength 
        * fills self.time. 1D numpy array
        * all histogram selects common time
        * PSI has different t0 per histogram
        * and must standardize to a common length 
        
        # Time definition for center of bin n: 
        #          time = (n - self.nt0 + self.offset.value + self.dt0)*binWidth_ns/1000.
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
                               self.nt0.max(),dtype=float)   # 1D np.array
        binwidth_mus = self._the_runs_[0][0].get_binWidth_ns()/1000.
        self.histoLength = self._the_runs_[0][0].get_histoLength_bin() - self.nt0.max() - self.offset

        self.time = time_bins*binwidth_mus 
        
        if self.datafile[-3:]=='bin' or self.datafile[-3:]=='mdu': # PSI
            self.time += self.mean_dt0()*binwidth_mus # mean dt0 correction (fraction of a bin, probaby immaterial)
             
    def single_for_back_counts(self,runs,grouping):
        """
        * input: 
        *         runs, runs to add
        *         grouping, dict with list of detectors 
                            grouping['forward'] and grouping['backward']
        * output:
        *         yforw, ybackw  
        *                        = sum_{i=for or back}(data_i - background_i), PSI, 
        *                        = sum_{i=for or back}data_i, ISIS
        *         background_forw  =
        *         background_backw = average backgrounds
        *         yfbmean        = mean of (yforw-bf)*exp(t/TauMu)
        *         ybackw         = mean of (ybackw-bb)*exp(t/TauMu)
        * used both by self.asymmetry_single (see) in normal fits (alpha is fixed)
        * and directly in calibration fits (alpha is a parameter)
        * all are 1D numpy arrays
        """
        import numpy as np
        from mujpy.aux.aux import TauMu_mus

        filespec = self.datafile[-3:] # 'bin', 'mdu' or 'nsx'
        if self.loadfirst:
            
    #       initialize to zero self.histolength
            n1 = self.nt0[0] + self.offset # ISIS
            n2 = n1 + self.histoLength # ISIS
            yforw, ybackw = np.zeros(self.histoLength), np.zeros(self.histoLength) # counts 
            background_forw, background_backw = 0., 0. # background estimate
                           
            for j, run in enumerate(runs): # Add runs
                #print(run)
                for counter in grouping['forward']: # first good bin, last good bin from data array start
                
                    histo = run.get_histo_array_int(counter) # counter data array in forw group
                    if filespec =='bin' or filespec=='mdu': # PSI, counter specific range                  
                        n1 = self.nt0[counter] + self.offset
                        n2 = n1 + self.histoLength 
                        background_forw += np.mean(histo[self.firstbin:self.lastbin])  # from prepromt estimate
                    yforw += histo[n1:n2]
                        
                for counter in grouping['backward']: # first good bin, last good bin from data attay start
                
                    histo = run.get_histo_array_int(counter) # counter data array in back group
                    if filespec=='bin' or filespec=='mdu': #  PSI, counter specific range  
                        n1 = self.nt0[counter] + self.offset
                        n2 = n1 + self.histoLength 
                        background_backw += np.mean(histo[self.firstbin:self.lastbin])  # from prepromt estimate
                    ybackw += histo[n1:n2]              

            x = np.exp(self.time/TauMu_mus())
            yfm = np.mean((yforw-background_forw)*x)
            ybm = np.mean((ybackw-background_backw)*x)
            return yforw, ybackw, background_forw, background_backw, yfm, ybm
        else:
            return None, None, None, None, None, None
            
        # Error eval box:
        # Nf(t), Nb(t) are row counts from the two groupings, forward and backward
        # A(t) = y  with background corrected counts Nfc(t) = Nf(t) - bf = yf, 
        #                                            Nbc(t) = Nb(t) - bb = yb
        # errors eA(T) with renormalized counts
        #                                 Nf(t) = cf,
        #                                 Nb(t) = cb

        #############################################################
        #  ISIS)         Error evaluation, no backgrounds:          #
        # Brewers trick to avoid double error propagation:          #
        # the denominator is evaluated as an average                #
        #       A = [yf(t) - alpha yb(t)]/d          with           #
        #    d = (<yf e^t/tau> + alpha <yb e^t/tau>)e-t/tau         #
        # yf = sum_{i in f) Ni        yb = sum_{i in b} Ni          #
        # ef^2 = yf                   eb^2 = yb                     #
        #          eA = sqrt(yf + alpha^2 yb)/d                     #
        #-----------------------------------------------------------#
        # PSI)          With background                             #
        # evaluate bf, bb before prompt for yf, yb respectively     #
        #                                                           #
        #  A = [yf-bf - alpha(yb-bb)]/[yf-bf + alpha(yb-bb)]        #
        #           ef^2, eb^2 are the same                         #
        # d = [<(yf-bf)e^t/tau)> + alpha <(yb-bb)e^t/tau>] e^-t/tau #
        #   =   [<yfbe>     +    alpha     <ybbe>] e^-t/tau         #
        #                                                           #
        #     A = [yf - alpha yb - (bf - alpha bb)]/d               #
        #                                                           #
        #    eA = sqrt(yf + alpha^2 yb)/d                           #
        #-----------------------------------------------------------#
        # if alpha is a paramter                                    #
        # compute and return yf, yb, bf, bb, <yfbe>, <ybbe>         #
        # mumodel must compute   d,  A, eA                          #
        #############################################################
        # for ISIS the PSI formula work with bb and bf zero         #
        #############################################################
        # rebin eA works for ISIS, must be corrected for PSI        #
        # yfm, ybm depend on binning   

    def single_multigroup_for_back_counts(self,runs,groupings):
        """
        * input: 
        *         runs, runs to add
        *         grouping, dict with list of detectors 
                            grouping['forward'] and grouping['backward']
        * output:
        *         yforw, ybackw  
        *                        = sum_{i=for or back}(data_i - background_i), PSI, 
        *                        = sum_{i=for or back}data_i, ISIS
        *         background_forw  =
        *         background_backw = average backgrounds
        *         yfbmean        = mean of (yforw-bf)*exp(t/TauMu)
        *         ybackw         = mean of (ybackw-bb)*exp(t/TauMu)
        * used only by calib multigroups
        * yforw, ybackw are 2D numpy arrays, the last four output items are arrays 
        """
        from numpy import vstack,array
        bf,bb,yfm,ybm = [],[],[],[]
        for k,grouping in enumerate(groupings):
            yforw, ybackw, background_forw, background_backw, yforwm, ybackwm = self.single_for_back_counts(runs,grouping)
            bf.append(background_forw)
            bb.append(background_backw)
            yfm.append(yforwm)
            ybm.append(ybackwm)
            if k:
                yf = vstack((yf,yforw))
                yb = vstack((yb,ybackw))
            else:
                yf = yforw
                yb = ybackw
        bf = array(bf)
        bb = array(bb)
        yfm = array(yfm)
        ybm = array(ybm)
        return yf,yb,bf,bb,yfm,ybm
        

    def asymmetry_single(self,the_run,kgroup):
        """
        input:
            the_run = list containing the instance[s] of the run[s to be added]
            k = index of self.grouping, a list of dicts 
                self.grouping[k]['forward'] and ['backward']
                containing the respective lists of detectors
        * run instances from musr2py/muisis2py  (psi/isis load routine) 
        *
        outputs: 
            asymmetry and asymmetry error (1d)
         """
        from numpy import exp, sqrt, where
        from mujpy.aux.aux import TauMu_mus

        if self.loadfirst:
            # print(the_run)
            alpha = self.grouping[kgroup]['alpha']
            
            yf, yb, bf, bb, yfm, ybm = self.single_for_back_counts(the_run,self.grouping[kgroup])
            
            # calculate asymmetry and error
            denominator = (yfm + alpha*ybm)*exp(-self.time/TauMu_mus())
            asymm = (yf - alpha*yb - (bf - alpha*bb)) / denominator 
            
            #   ey_i in fcn sum_i ((y_i-y_th)/ey_i)**2 cannot be zero, but errexp can
            errexp = sqrt(yf + alpha**2*yb)
            errexp[where(errexp==0)] = 1                                                               #
            #   set to 1 the minimum error (weights less very few points closer to ~ zero) 
            asyme = errexp / denominator 

            return asymm, asyme
        else:
            return None, None
            
    def asymmetry_multirun(self,kgroup):
        """
        input:
                kgroup, index forward - backward pair 
                    self.grouping[kgroup]['forward'] and ['backward']
                    containing the respective lists of detectors
        * uses the suite of run instances from musr2py/muisis2py  (psi/isis load routine) 
        *
        outputs: 
            asymmetry and asymmetry error (2d)
                 also generates self.time (1d)
        """
        from numpy import vstack

        if self.loadfirst:
            for k,run in enumerate(self._the_runs_):
                a,e = self.asymmetry_single(run,kgroup)
                if k==0:
                    asymm, asyme  = a, e
                else:
                    asymm, asyme = vstack((asymm,a)), vstack((asyme,e))
            return asymm, asyme
        else:
            return None, None

    def asymmetry_multigroup(self):
        """
        input: none
            calls self.asymmetry_single which calls self.single_for_back_counts
        outputs: 
            asymmetry and asymmetry error (2d)
         """
        from numpy import vstack

        if self.loadfirst:
            if not self.multi_groups():
                self.console('** ONLY ONE GROUP! Use asymmetry_single instead') 
            run = self._the_runs_[0]   # must be only one run, switch brings here only if self.suite.single     
            if not self.single():
                self.console('** You are programmatically invoking asymmetry_multigroup with a multi-run suite')
                self.console('*  Only the first run in the suite will be analysed') 
            for kgroup in range(len(self.grouping)):
                a,e = self.asymmetry_single(run,kgroup)
                # self.console('Loaded run {}, group {} ({}), alpha = {}'.format(run[0].get_runNumber_int(), kgroup, 
#                                                  self.groups[kgroup]['forward']+'-'+self.groups[kgroup]['backward'],
#                                                  self.groups[kgroup]["alpha"]))  
                if kgroup==0:
                    asymm, asyme  = a, e
                else:
                    asymm, asyme = vstack((asymm,a)), vstack((asyme,e))
            return asymm, asyme
        else:
            self.console('** CHECK ACCESS to database (or load runs first)') 
            return None, None
                
    def single(self):
        '''
        True if there is a single run (fit type A)
        False if there are many runs (fit types B and C)
        Usage:
            if single:
                do something with single run
            else:
                do something with multiple runs
        '''
        try:
            test = len(self._the_runs_)==1
            
        except:
            self.console('Warning: data are not available: access expired?')
        return test
            
    def multi_groups(self):
        '''
        True if more groups
        False if just one group
        Usage:
            if multi_groups:
                do something with single group
            else:
                do something with multi groups
        '''
        # print('multi_group suite debug: self.grouping {} len {}'.format(self.grouping,len(self.grouping)))
        return len(self.grouping)>1        
        
