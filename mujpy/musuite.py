class suite(object):
    '''
    A suite class, 
    self._init_ reads input fixed variables
        console, runlist, datafile 
        grp_calib, offset 
        (console(string) must print string somewhere)
        t0 parameters are fixed depending on bin or mdu spec
    
    Self.groups = grp.calib contains a list of dictionaries for forward, backward groups and their alpha values
        where groups may be in shorthand notation
        self.alpha value used for normal fit 
        as opposed to calibration fits where alpha is a fit parameter
    self.grouping is a list of dictionaries for forward, backward groups and their alpha values
        where groups as np array of 0 based indices of detectors (tools get_grouping does the translation)
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
        Output ends up in notebook, below cell, and 
        sys.__stdout__ is <_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>
        sys.__stderr__ is <_io.TextIOWrapper name='<stderr>' mode='w' encoding='utf-8'>
    ''' 
    def __init__(self, datafile , runlist , grp_calib , offset , startuppath, console = 'print',dash=None, mplot=False):
                 
        '''
        * Initiates an instance of suite, 
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
        if dash = None and cosole = 'print' self.console will exec print(string)
        if dash = self in mudash calls to suite, and console = 'self.dash.log'   
              self.console will exec self.dash.log(string), i.e. write on board output                           
        '''
        from mujpy.tools.tools import derun
        import json
        import os
        from mujpy import __file__ as MuJPyName
        from mujpy._version import __version__
        self.__version__ = __version__
        # print('__init__ now ...')
        self.firstbin = 0
        self.loadfirst = False  
        self._the_runs_ = [] # initialized 
#        with open(suite_input_file,"r") as f:
#            suite_input_json = f.read()
#            suite_input = json.loads(suite_input_json)
        self.dash = dash
        self.console_method = console
        try:
            if not dash: # not a call from a mudash instance
                self.console('******************* SUITE *********************')
            else:
                self.console('')
        except:
            print('no suite.console!')
            return  # with no console error message, it means no console
        # determine number of runs, filenames etc.
            #######################
            # decode the run string
            #######################           
        self.runs, errormessage = derun(runlist) # self.runs is a list of lists of run numbers (string)
        if errormessage is not None: # derun error
            self.console('Run syntax error: {}. You typed: {}'.format(errormessage,runlist))
            return  # with console error message
        if os.path.isfile(datafile):
            self.datafile = datafile
            if self.datafile[-3:]=='bin': 
                self.thermo = 1 # sample thermometer is 1 on gps (check or adapt to other instruments)
                self.prepostpk = [50, 50]
            elif self.datafile[-3:]=='mdu':
                self.thermo = 1 # sample thermometer is 1 on gps (check for hifi)
                self.prepostpk = [70, 70]
            self.loadfirst = True
        else:
            self.console('* File {} not found'.format(datafile))
            self.console('* CHECK YOUR ACCESS (e.g. klog)')
            self.console('**************************************')
            return None  # with console error message
        try: # working directory, in which, in case, to find mujpy_setup.pkl 
            self.__startuppath__ = startuppath 
        except:
            self.__startuppath__ = os.getcwd()
        # implement new folder policy, see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.Dashboard
        self.__path__ = os.path.dirname(MuJPyName) # mujpy path
        self.__templatepath__ = os.path.abspath(os.path.join(os.path.dirname(MuJPyName), '../', 'templates/'))
        if not os.path.exists(self.__startuppath__+'/fit/'):
            os.mkdir(self.__startuppath__+'/fit/')
        self.__fitpath__ = self.__startuppath__+'/fit/'
        if not os.path.exists(self.__startuppath__+'/csv/'):
            os.mkdir(self.__startuppath__+'/csv/')
        self.__csvpath__ = self.__startuppath__+'/csv/'
        if not os.path.exists(self.__startuppath__+'/cache/'):
            os.mkdir(self.__startuppath__+'/cache/')
        self.__cachepath__ = self.__startuppath__+'/cache/'
        # reverse of tools get_grouping is tools get_group
        self.offset = int(offset) # offset belongs to suite, that needs it for asymmetries
        self.load_runs() #           load data instances in self._the_runs_
        self.groups = grp_calib
        # grp_calib is a list of dictionaries, one per group, with keys 'forward,'backward','alpha'
        # where groups may be in shorthand notation
        self.grouping = [] # reproduce the same with arrays of histogram numbers, done by tools get_grouping
        self.store_groups() #        in self.grouping        self.promptfit(mplot)   #    to be done: make switch for ISIS
        self.promptfit(mplot)   #    to be done: make switch for ISIS
        self.timebase()
        # self.console('... end of initialize suite')

    def add_runs(self,k):
        '''
        Tries to load one or more runs to be added together
        by means of murs2py. 
        self.runs is a list of strings containing integer run numbers 
        Returns -1 and quits if musr2py or muisis2py complain, 0 otherwise
        '''
        
        from musr2py import MuSR_td_PSI_bin as psiload
        from mujpy.muisis2py.muisis2py import muisis2py as isisload
        # muisis2py has the same methods as musr2py
        from mujpy.tools.tools import get_datafilename, get_title

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
               	self.console('{} loaded'.format(path_and_filename))
               	self.console('Run {}: {}'.format(run,get_title(self._the_runs_[-1][0])))
#                import sys
#                self.console('sys.__stdout__ is {}, sys.__stderr__ is {}'.format(sys.__stdout__,sys.__stderr__))

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

    def check_group(self,group):
        '''
        check that this is a group of existing detectors
        '''
        return (group>=0).all()*(group<self._the_runs_[0][0].get_numberHisto_int()).all()

    def store_groups(self):
        '''
        reads groups dictionary in dashboard shortnote
        and appends lists of histogram numbers, alphas to self.grouping 
        '''
        from mujpy.tools.tools import get_grouping
        for k,group in enumerate(self.groups):
            fgroup, bgroup, alpha = get_grouping(group['forward']), get_grouping(group['backward']), group['alpha']
            if alpha>0 and self.check_group(fgroup) and self.check_group(bgroup): # checks legal grpcalib_file
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
        weight = [float(sum(self._the_runs_[k][j].get_histo_vector(counter,1))) for counter in range(self._the_runs_[0][0].get_numberHisto_int()) for j in range(m)]
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
        from numpy import array, where, arange, zeros, mean, ones, sqrt, linspace
        from iminuit import Minuit, cost
        
        import matplotlib.pyplot as P
        from mujpy.mucomponents.muprompt import muprompt
        from mujpy.mucomponents.muedge import muedge
        from mujpy.tools.tools import TauMu_mus, scanms, step, set_fig 
    
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
            elif self.datafile[-3:]=='nxs': # ISIS
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

        elif self.datafile[-3:]=='nxs': # ISIS
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
        from numpy import zeros, array, mean, exp, where
        from mujpy.tools.tools import TauMu_mus

        filespec = self.datafile[-3:] # 'bin', 'mdu' or 'nsx'
        if self.loadfirst:
            
    #       initialize to zero self.histoLength, maximum available good bins valid for all histos 
            n1 = self.nt0[0] + self.offset # ISIS
            n2 = n1 + self.histoLength # ISIS
#            print('musuite single_for_back_counts debug: n1 {}, n2 {}, self.histoLength {}'.format(n1,n2,self.histoLength))
            yforw, ybackw = zeros(self.histoLength), zeros(self.histoLength) # counts 
            background_forw, background_backw = 0., 0. # background estimate
                           
            for j, run in enumerate(runs): # Add runs
                #print(run)
                for counter in grouping['forward']: # first good bin, last good bin from data array start
                
                    histo = array(run.get_histo_vector(counter,1)) # counter data array in forw group
#                    if array(where(histo==0)).size:
#                        print('musuite single_for_back_counts debug: run {} counter fwd {} contains {}  zeros'.format(run.get_runNumber_int(),counter,array(where(histo==0)).size))
                    if filespec =='bin' or filespec=='mdu': # PSI, counter specific range                  
                        n1 = self.nt0[counter] + self.offset
                        n2 = n1 + self.histoLength 
                        background_forw += mean(histo[self.firstbin:self.lastbin])  # from prepromt estimate
                    yforw += histo[n1:n2]
                        
                for counter in grouping['backward']: # first good bin, last good bin from data attay start
                
                    histo = array(run.get_histo_vector(counter,1)) # counter data array in back group
#                    if array(where(histo==0)).size:
#                        print('musuite single_for_back_counts debug: run {} counter bkw {} contains {}  zeros'.format(run.get_runNumber_int(),counter,array(where(histo==0)).size))
                    if filespec=='bin' or filespec=='mdu': #  PSI, counter specific range  
                        n1 = self.nt0[counter] + self.offset
                        n2 = n1 + self.histoLength 
                        background_backw += mean(histo[self.firstbin:self.lastbin])  # from prepromt estimate
                    ybackw += histo[n1:n2]              

            x = exp(self.time/TauMu_mus())
            yfm = mean((yforw-background_forw)*x)
            ybm = mean((ybackw-background_backw)*x)
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

#    def single_multigroup_for_back_counts(self,runs,groupings):
#        """
#        # unused?!
#        * input: 
#        *         runs, runs to add
#        *         grouping, dict with list of detectors 
#                            grouping['forward'] and grouping['backward']
#        * output:
#        *         yforw, ybackw  
#        *                        = sum_{i=for or back}(data_i - background_i), PSI, 
#        *                        = sum_{i=for or back}data_i, ISIS
#        *         background_forw  =
#        *         background_backw = average backgrounds
#        *         yfbmean        = mean of (yforw-bf)*exp(t/TauMu)
#        *         ybackw         = mean of (ybackw-bb)*exp(t/TauMu)
#        * used only by calib multigroups
#        * yforw, ybackw are 2D numpy arrays, the last four output items are arrays 
#        """
#        from numpy import vstack,array
#        bf,bb,yfm,ybm = [],[],[],[]
#        for k,grouping in enumerate(groupings):
#            yforw, ybackw, background_forw, background_backw, yforwm, ybackwm = self.single_for_back_counts(runs,grouping)
#            bf.append(background_forw)
#            bb.append(background_backw)
#            yfm.append(yforwm)
#            ybm.append(ybackwm)
#            if k:
#                yf = vstack((yf,yforw))
#                yb = vstack((yb,ybackw))
#            else:
#                yf = yforw
#                yb = ybackw
#        bf = array(bf)
#        bb = array(bb)
#        yfm = array(yfm)
#        ybm = array(ybm)
#        return yf,yb,bf,bb,yfm,ybm
        
    def asymmetry_single(self,the_run,kgroup):
        """
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
        from numpy import exp, sqrt, where, array, intersect1d, finfo
        from mujpy.tools.tools import TauMu_mus

        if self.loadfirst:
            # print(the_run)
            alpha = self.grouping   [kgroup]['alpha']
            
            yf, yb, bf, bb, yfm, ybm = self.single_for_back_counts(the_run,self.grouping[kgroup])
            
            # yf, yb >=0 are flaot(int) counts; bf, bb are float average background counts
            # yfm, ybm are  <(y-b)*exp(t/tau>, could be negative?
            # Now calculate asymmetry and error, adding inner list runs
            # self.grouping py-index, i.e "counter 1-Backw" is  0
            if (yfm + alpha*ybm)<=0:
                self.console('Too low counts on run {} group {}: negative asymmetry denominator: exiting'.format(the_run[0].get_numberRun_int(),kgroup))
                return None, None
            denominator = (yfm + alpha*ybm)*exp(-self.time/TauMu_mus()) # yfm, ybm are mean, not == 0
            asymm = (yf - alpha*yb - (bf - alpha*bb)) / denominator  # 1d array
            
            #   ey_i in fcn sum_i ((y_i-y_th)/ey_i)**2 cannot be zero, but errexp can
            norma = yf + alpha**2*yb # could be zero, yielding zero error 
            # norma represents (corrected) total counts, 
            # its minimum physical value is 1 
            # a) can be zero when both yf and yb are zeros
            # b) can be 1 when yf = 1 and yb = 0
            # c) can be alpha**" when yf = 1 and yb = 0, which are both fine
            # cure by setting to 1 only case a), others take care of themselves
            norma[norma==0] = 1.
            asyme = sqrt(norma) / denominator
            
            # further check (effects of funny denominator? Maybe can be removed
            sqepsi = sqrt(finfo('d').eps)
            if (array(where(asyme<sqepsi))).size: # should never happen now
                self.console('suite asymmetry_single debug: run {}, kgroup {} asyme contains {} value(s) < 1.5e-8 This should NOT happen'.format(the_run[0].get_runNumber_int(),kgroup ,(array(where(asyme<sqepsi))).size))
                asyme[where(asyme<epsi)] = 2*sqepsi

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
        # can be B1, C1 fits 
        outputs: 
            asymmetry and asymmetry error (2d)
                 also generates self.time (1d)
        """
        from numpy import vstack

        if self.loadfirst:
            for k,run in enumerate(self._the_runs_):
                a,e = self.asymmetry_single(run,kgroup)
                if a is None: 
                    return None, None
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
            # can be A20, A21 fits 
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
                if a is None: 
                    return None, None
                if kgroup==0:
                    asymm, asyme  = a, e
                else:
                    asymm, asyme = vstack((asymm,a)), vstack((asyme,e))
            return asymm, asyme
        else:
            self.console('** CHECK ACCESS to database (or load runs first)') 
            return None, None
            
    def asymmetry_multirun_multigroup(self):
        '''
        input: none
            calls self.asymmetry_single which calls self.single_for_back_counts
        outputs: 
            # can be B20, B21 or C2 fit 
            asymmetry and asymmetry error (3d)
            for run in runs:
                for group in groups:
                    np.vstack # axis=1
                np.vstack # axis=0
        '''
        from numpy import array, vstack
        if self.loadfirst:
            for krun,run in enumerate(self._the_runs_):
                for kgroup in range(len(self.grouping)):
                    if kgroup:
                        a,e = self.asymmetry_single(run,kgroup)
                        if a is None: 
                            return None, None
                        asy,ase = vstack((asy,a)), vstack((ase,e)) # groups are vstacked
                    else: # kgroup = 0
                        asy,ase = self.asymmetry_single(run,kgroup)
                if krun:
                    asymm, asyme  = vstack((asymm,array([asy]))), vstack((asyme,array([ase])))
                else: # krun=0
                    asymm, asyme = array([asy]),array([ase]) # into 2nd dimension
            return asymm, asyme
        else:
            self.console('** CHECK ACCESS to database (or load runs first)') 
            return None, None
               
    def single(self):
        '''
        Boole test
        output:
            True if there is a single run (fit type A)
            False if there are many runs (fit types B and C)
        '''
        try:
            test = len(self._the_runs_)==1
            
        except:
            self.console('Warning: data are not available: access expired?')
        return test
            
    def multi_groups(self):
        '''
        Boole test
        output:
            True if more groups (fits A2, B2, C2)
            False if just one group (fits A1, B1, C1)
        '''
        # print('multi_group suite debug: self.grouping {} len {}'.format(self.grouping,len(self.grouping)))
        return len(self.grouping)>1     
        
    def scan(self):
        '''
        output
            False if single
            'B[mT]' if it's a B scan
            'T[K] ' if it's a T scan
            '[deg]' if it's an angle scan 
            '#    ' if it's another scan
        '''  
        if self.single(): return False
        elif [run[0].get_temp() for run in self._the_runs_]!=[self._the_runs_[0][0].get_temp()]*len(self._the_runs_): return 'T[K] '
        elif [run[0].get_field() for run in self._the_runs_]!=[self._the_runs_[0][0].get_field()]*len(self._the_runs_): return 'B[mT]'
        elif [run[0].get_orient() for run in self._the_runs_]!=[self._the_runs_[0][0].get_orient()]*len(self._the_runs_): return '[deg]'
        else: return '#    '
        
        
