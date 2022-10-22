class mufitplot(object):
    '''
    input 
       plot_range: 2 to 5 values
          None (skip time plot)
          start stop
          start stop pack
          start stopearly packearly stop packlate
       the_suite: a suite instance
       the_fit: pass direclty the mufit instance, 
          with all its methods, including .suite 
       guess: True, False (default)
       rotating_frame_frequencyMHZ: 0.0 default
       fit_range: 1 or 3 values
          None default
          start, stop, sigma (LB) in MHz  
    mufitplot class 
    produces the fitplot
    multiple runs sequential produce anim
    can plot in the rotating frame
    optional fft of residues toggles
    fig handle is now handled automatically 
    '''
# A bit of discussion to be removed when done
# plot needs to access a version of mumodel, 
# that already stored in the_fit,
# with different ranges for rebinned data, 
# denser function plot, original data
# mumodel() is needed for generating model function
# and calculating chi2.
# everything can be obtained from the_fit, including dashboard initial guess
# Maybe make another version that does multiplots from the dashboards
#
#
#  the function that performs the plot is mujpy.aux.aux import set_single_fit, self_sequence_fit 
#        invoked by self.plot_run in 
#        plot_singlerun or plot_singlerun_multigroup_sequential or plot_singlerun_multigroup_userpar

    def __init__(self,plot_range,the_fit,guess=False,rotating_frame_frequencyMHz = 0.0,fft_range = None,real=True,fig_fit=None,fig_fft=None):
        '''
        input: 
            either
                suite instance 
                dashboard best fit json file 
            or 
               fit instance
            fig fig_handle
        '''
#        from mujpy.mucomponents.mucomponents import mumodel
        
# reassigned, for backward compatibility  
        # print('mufitplot __init__ debug: inside')      
        self.lastfits = the_fit.lastfits
        self._the_model_ = the_fit._the_model_                           
        self.suite = the_fit.suite
        self.log = the_fit.log
        if not self.suite.loadfirst: 
            self.log('Sorry, no access to data... quitting mufitplot')
            ok = False
        else:
            ok = True 
#        Reassigns dashboard for backward compatibility.
#        Dashboard, from the_fit, is just needed for the initial guess 
#        Results in fft are now passed directly through lastfit 
            self.dashboard = the_fit.dashboard
        #print('mufitplot __init__ debug: the_fit.dashboard =\n{}'.format(self.dashboard))
        self.fig = fig_fit
        self.fig_fft = fig_fft
        if ok:
            self.rotating_frame_frequencyMHz = rotating_frame_frequencyMHz
            self.guess = guess
            self.model = "model_guess" 
            if not self.guess:
                if "model_result" in self.dashboard:
                    self.model = "model_result"
                elif "userpardicts_guess" in self.dashboard.keys():   
                    self.model = "model_result"  
                else:            
                    # self.log('Sorry no fit results yet, plotting guess instead')
                    self.guess = True
            # print('__init__ mufitplot debug: model = {}, guess = {}'.format(self.model,self.guess))    
           
            if plot_range is not None:
                self.chooseplot(plot_range)
                if fft_range is not None:
                    self.choosefftplot(fft_range,real)
            else:
                if fft_range is None:
                    self.log('Exiting mufitplot: nothing to do.')    
                              
  
    def chooseplot(self,plot_range):
        '''
            switch for single (A1), sequential (B1),
            ... 
        #       multi_groups suite.single userpars multigroup_in_components userlocals
        # A1     False          True        False       False               False *
        # A20    True           True        False       False               False
        # A21    True           True        True        True                False *
        # B1     False          False       False       False               False
        # B20    True           False       False       False               False does not exist
        # B21    True           False       True        True                False
        # C1     False          False       True        False               True  *
        # C2     True           False       True        True                True  *
        '''
        from mujpy.aux.aux import multigroup_in_components, userpars, userlocals
        if self.suite.single(): # A1, A20, A21
            if self.suite.multi_groups(): # A20 A21
                self.log('Multigroup fit animation: toggle pause/resume by clicking on the plot')
                if sum(multigroup_in_components(self.dashboard)): # A21
                    ok, msg = self.plot_singlerun_multigroup_userpar(plot_range)
                else:                          # A20
                    ok, msg = self.plot_singlerun_multigroup_sequential(plot_range)
            else:                              # A1
                ok, msg = self.plot_singlerun(plot_range)
        else: # B1, B2, C1, C2, no calib in this lot
            self.log('Multirun fit animation: toggle pause/resume by clicking on the plot')
            if self.suite.multi_groups(): # B2, C2
                if userlocals(self.dashboard): # C2
                    self.log('No C2 yet. Exiting mufitplot without a plot')
                    ok, msg = False, 'C2'
                    pass
                else: # B2 (means B21 group global, B20 group sequential does not exist)
                    ok, msg = self.plot_multirun_multigrup_userpar(plot_range)
                    pass
            else: # B1, C1
                if userpars(self.dashboard): # C1
                    self.log('No B2, C2 yet. Exiting mufitplot without a plot')
                    ok, msg = False, 'C1'
                    pass
                else: # B1
                    ok, msg = self.plot_multirun_singlegroup_sequential(plot_range)       
        if not ok:
            self.log('Exiting mufitplot without a plot: '+msg)
    
    def plot_singlerun(self,plot_range):
        '''
        input plot_range, passed to self.plot_run together with
            pars, list of one list of fit parameter values
            asymm, asyme 1d        
        calib pops first par and passes it as alpha to standard plot
        '''
        from mujpy.aux.aux import int2min, mixer, calib
        from numpy import cos, pi
#        self.log('muplotfit plot_singlerun debug')
        kgroup = 0 # default single group 
        pars,_,_,_,_,_ = int2min(self.dashboard[self.model])
        if calib(self.dashboard):
            self.suite.grouping[kgroup]['alpha'] = pars[0] # from fit parameter to standard asymmetry mode
            self.suite.groups[kgroup]['alpha'] = pars[0]
        asymm, asyme = self.suite.asymmetry_single(self.suite._the_runs_[0],0)
        if self.rotating_frame_frequencyMHz:
            self.rrf_asymm = mixer(self.suite.time,asymm,self.rotating_frame_frequencyMHz)
            self.rrf_asyme = mixer(self.suite.time,asyme,self.rotating_frame_frequencyMHz)
        #self.log('mufitplot: Inside single plot; debug mode')
        return self.plot_run(plot_range,pars,asymm,asyme)        

    def plot_singlerun_multigroup_userpar(self,plot_range):
        '''
        input plot_range, passed to self.plot_run together with
            pars, list of one list of fit parameter values
            # A21
            asymm, asyme 2d        
        '''
        from mujpy.aux.aux import int2min_multigroup, mixer, calib
        from numpy import cos, pi, vstack
        userpardicts = (self.dashboard["userpardicts_guess"] if self.guess else 
                        self.dashboard["userpardicts_result"])
        pardict = self.dashboard["model_guess"][0]["pardicts"][0]
        pars,_,_,_,_,_ = int2min_multigroup(userpardicts)
        p = pars
        if calib(self.dashboard):
            for kgroup,group in enumerate(self.suite.groups):
                group['alpha'] = eval(pardict["function_multi"][kgroup])
                self.suite.groups[kgroup]["alpha"] = eval(pardict["function_multi"][kgroup])
        asymm, asyme = self.suite.asymmetry_multigroup()        
        if self.rotating_frame_frequencyMHz:
            # print('mufitplot single_plot_multi debug: asymm.shape[0] = {}'.format(asymm.shape[0]))
            for k in range(asymm.shape[0]):
                if not k:
                    time = self.suite.time
                else:
                    time = vstack((time,self.suite.time))
            self.rrf_asymm = mixer(time,asymm,self.rotating_frame_frequencyMHz)
            self.rrf_asyme = mixer(time,asyme,self.rotating_frame_frequencyMHz)
        #self.log('mufitplot: Inside single plot; debug mode')
        return self.plot_run(plot_range,pars,asymm,asyme)        

    def plot_multirun_multigrup_userpar(self,plot_range):
        '''
         inputs: 
            plot_range
        # B2 fit
        asymm, asyme d
        '''
        from mujpy.aux.aux import int2min_multigroup, mixer, calib
        from numpy import cos, pi, vstack
        userpardicts = (self.dashboard["userpardicts_guess"] if self.guess else 
                        self.dashboard["userpardicts_result"])
        pardict = self.dashboard["model_guess"][0]["pardicts"][0]
        pars,_,_,_,_,_ = int2min_multigroup(userpardicts)
# no calib mode!
#        p = pars
#        if calib(self.dashboard):
#            for kgroup,group in enumerate(self.suite.groups):
#                group['alpha'] = eval(pardict["function_multi"][kgroup])
#                self.suite.groups[kgroup]["alpha"] = eval(pardict["function_multi"][kgroup])
        asymm, asyme = self.suite.asymmetry_multirun_multigroup()
        if self.rotating_frame_frequencyMHz:
            for krun in range(asymm.shape[0]):
                for kgroup in range(asymm.shape[1]):
                    if not kgroup: # kgroup=0
                        tim = self.suite.time
                    else:
                        tim = vstack((time,self.suite.time))
                if not krun: # krun=0
                    time = array([tim])
                else:
                    time = vstack(time,array([tim]))                      
            self.rrf_asymm = mixer(time,asymm,self.rotating_frame_frequencyMHz)
            self.rrf_asyme = mixer(time,asyme,self.rotating_frame_frequencyMHz)
        #self.log('mufitplot: Inside single plot; debug mode')
        return self.plot_run(plot_range,pars,asymm,asyme)        
   
    def plot_singlerun_multigroup_sequential(self,plot_range):
        '''
        inputs: 
            plot_range, passed to self.plot_run together with
            pars, list of lists  of fit parameter values
            asymm, asyme 2d
            # A20 fit
        if model_result
            self.lastfit is a list of lists 
        if model_guess
            reproduces the same from a single guesses            
        '''
        from mujpy.aux.aux import int2min, mixer
        from numpy import cos, pi, vstack
        # dashboard is a multi_sequential thing: each sequential fit has its own
        asymm, asyme = self.suite.asymmetry_multigroup()
        if self.rotating_frame_frequencyMHz:
            for k in range(asymm.shape[0]):
                if not k:
                    time = self.suite.time
                else:
                    time = vstack((time,self.suite.time))  ### TIME BECOMES MULTIDIM, NEEDED IN RRF?
            self.rrf_asymm = mixer(time,asymm,self.rotating_frame_frequencyMHz)
            self.rrf_asyme = mixer(time,asyme,self.rotating_frame_frequencyMHz)
        pars = []
        # self.log('mufitplot single_plot_multi_sequential debug')
        for kgroup in range(asymm.shape[0]):
            if self.model == "model_result":
                #print('single_plot_multi_sequential mufitplot debug: self model = {}',format(self.dashboard[self.model][kgroup]))
                values = self.lastfits[kgroup].values                
            else:
                #print('single_plot_multi_sequential mufitplot debug: self model = {}',format(self.dashboard[self.model]))
                values,_,_,_,_,_ = int2min(self.dashboard[self.model])                
            pars.append(values) # here pars is list of list!!
        return self.plot_run(plot_range,pars,asymm,asyme)        

    def plot_multirun_singlegroup_sequential(self, plot_range):
        '''
        input plot_range, passed to self.plot_run together with
            pars, list of lists  of fit parameter values
            asymm, asyme 2d
            # B1 fit
        '''
        from mujpy.aux.aux import int2min, mixer
        from numpy import cos, pi, vstack
        
        kgroup = 0 # default single 
        # dashboard must become a suite thing: each sequential fit has its own
        pars = []
        # self.log('mufitplot: Inside sequential plot; debug mode')
        for k,lastfit in enumerate(self.lastfits):
            values = lastfit.values
            
            #print('plot_run muplotfit debug: pars = {}'.format(pars)) run {},  values = {}'.format(run[0].get_runNumber_int(),values))
            pars.append(values)
        asymm, asyme = self.suite.asymmetry_multirun(kgroup) # 
        if self.rotating_frame_frequencyMHz:
            for k in range(asymm.shape[0]):
                if not k:
                    time = self.suite.time
                else:
                    time = vstack((time,self.suite.time))
            self.rrf_asymm = mixer(time,asymm,self.rotating_frame_frequencyMHz)
            self.rrf_asyme = mixer(time,asyme,self.rotating_frame_frequencyMHz)
        return self.plot_run(plot_range,pars,asymm,asyme)        
     
    def plot_run(self,plot_range,pars,asymm,asyme):
        '''
        input :
            plot_range
                (start,stop)
                (start,stop,pack)
                (start_early,stop_early,pack_early,stop_late,pack_late)
            pars list (single, calib) or list of lists (sequence)
            asymm, asyme  1d (single, calib) or 2d (sequence)
        calls either 
            self.chi_1 or self.chi_2 (stats and model function, 
                                      deals also with plotting model function)
            set_single_fit or set_sequence_fit, from aux.plot
                                    (produces actual figures)
        Draws either static or anim figures using aux.plot functions 
        '''
        from mujpy.aux.aux import derange, rebin, multigroup_in_components
        from mujpy.aux.aux import get_run_title, userpars, mixer
        from mujpy.aux.plot import set_single_fit, set_sequence_fit
        from iminuit import Minuit
        from numpy import ones 

        # print('plot_run muplotfit debug: pars = {}'.format(pars))
        # chi_1: single chisquare, different asymm.shape, 
        # chi_2: many chisquare, different asymm.shape, 
        chi = self.chi_1 if self.single_chi() else self.chi_2
        # print('plot_run mufitplot debug: single chi is {}'.format(self.single_chi()))

        run_title = get_run_title(self.suite)    # always a list, even for single 
        string = 'global ' if userpars(self.dashboard) else ''  
                  
        if self.guess:
            run_title = [title + ": "+string+"guess values" for title in run_title]
        else:
            run_title = [title + ": "+string+"fit results" for title in run_title]
        plottup,ermsg = derange(plot_range,self.suite.histoLength)
#        self.log('mufitplot plot_run debug: plot_range {}, plottup {}'.format(plot_range,plottup))
        #############################
        # rebinning of data as in fit 
        # this works for single 
        # and for sequential
        #############################
        fittup,ermsg = derange(self.dashboard["fit_range"],self.suite.histoLength) # range as tuple
        fit_pack = 1
        if len(fittup)==3: # plot start stop pack
            fit_start, fit_stop, fit_pack = fittup[0], fittup[1], fittup[2]
        elif len(fittup)==2: # plot start stop
            fit_start, fit_stop = fittup[0], fittup[1]

        # load modules and reproduce fit
        t_fit,y_fit,ey_fit = rebin(self.suite.time,asymm,[fit_start,fit_stop],fit_pack,e=asyme)
        # single slices as in fit 
        # print('muplotfit plot_run debug: calling chi for f_fit')
        nu_fit, f_fit, chi_fit = chi(t_fit,y_fit,ey_fit,pars)
        if self.rotating_frame_frequencyMHz: # (t_fit,y_fit,f_fit,ey_fit) the last three must be transformed to rrf
            _,y_fit,ey_fit = rebin(self.suite.time,self.rrf_asymm,[fit_start,fit_stop],fit_pack,e=self.rrf_asyme)
            f_fit = mixer(t_fit,f_fit,self.rotating_frame_frequencyMHz)
        if not nu_fit:
            return False, 'plot_run y_fit was None'
        # function as in fit

        if len(plottup)==5: 
        ###################
        # double range plot
        ###################
            early_late = True
            start, stop, pack, last, packlate = plottup
            
            t_late,y_late,ey_late = rebin(self.suite.time,asymm,[stop,last],packlate,e=asyme)
            # rebinned late slices 
            packfit = int(packlate/2)
            tfl,dum = rebin(self.suite.time,asymm,[stop,last],packfit)
            # late time slice for function
            # tfl,fl for late plot curve 
            _, f_late_res,_  = chi(t_late,None,None,pars) # for f_late_res calculating residues (on data points)
            _,fl,_ = chi(tfl,None,None,pars)            # t_late,fls for plotting residues  
            if self.rotating_frame_frequencyMHz: # (t_late,y_late,f_late_res,ey_late) the last three must be transformed to rrf
                                                 # (tfl,fl) non rebinned, the last must be transformed to rrf
                _,y_late,ey_late = rebin(self.suite.time,self.rrf_asymm,[start,stop],packlate,e=self.rrf_asyme)
                f_late_res = mixer(t_late, f_late_res,self.rotating_frame_frequencyMHz) 
                fl = mixer(tfl, fl,self.rotating_frame_frequencyMHz) 
            
            fit_late_start = int(stop/packlate*fit_pack) # divides fit_range in early and late
            
            # redo the same with the original fit binning fit_pack, only for histo and chi2
            t_fit_late,y_fit_late,ey_fit_late = rebin(self.suite.time,asymm,        
                                                [fit_late_start,fit_stop],fit_pack,e=asyme)
            # no rotating frame version for this!
            nu_fit_late, f_fit_late, chi_fit_late = chi(t_fit_late,y_fit_late,ey_fit_late,pars)
            if self.rotating_frame_frequencyMHz: # (t_fit_late,y_fit_late,f_fit_late,ey_fit_late) the last three must be transformed to rrf
                _,y_fit_late,ey_fit_late = rebin(self.suite.time,self.rrf_asymm,[fit_start,fit_late_start],fit_pack,e=self.rrf_asyme)
                f_fit_late = mixer(t_fit_late, f_fit_late,self.rotating_frame_frequencyMHz) 
            # print('mufitplot plot_run debug: y_fit_late ey_fit_late f_fit_late shape = {},  {},  {}'.format(y_fit_late.shape, ey_fit_late.shape, f_fit_late.shape))

            if not nu_fit_late:
                return False, 'plot_run y_fit_late was None'
            # function as in fit

            t_fit_early,y_fit_early,ey_fit_early = rebin(self.suite.time,asymm,        
                                                [fit_start,fit_late_start],fit_pack,e=asyme)
            # single slices as in fit 
            nu_fit_early, f_fit_early, chi_fit_early = chi(t_fit_early,y_fit_early,ey_fit_early,pars)
            if self.rotating_frame_frequencyMHz:
                _,y_fit_early,ey_fit_early = rebin(self.suite.time,self.rrf_asymm,[fit_start,fit_late_start],fit_pack,e=self.rrf_asyme)
                f_fit_early = mixer(t_fit_early, f_fit_early,self.rotating_frame_frequencyMHz) 

            if not nu_fit_early:
                return False,  'plot_run y_fit_early was None'
            # function as in fit

        else:
        ###################
        # single range plot
        ###################
            early_late = False
            pack = 1
            chi_fit_late, nu_fit_late, chi_fit_early, nu_fit_early = None, None, None, None
            if len(plottup)==3: # plot start stop pack
                start, stop, pack = plottup
            elif len(plottup)==2: # plot start stop
                start, stop = plottup
#        self.log('plot_range= {}'.format(plot_range))
#        self.log('start, stop, pack = {},{},{}'.format(start, stop, pack))

        t,y,ey = rebin(self.suite.time,asymm,[start,stop],pack,e=asyme)
        # rebinned single or early slices 
        # print('muplotfit plot_run debug: calling chi for f_res')
        nudum, f_res, _ = chi(t,None,None,pars)
        packfit = int(pack/2)
        tf,dum = rebin(self.suite.time,asymm,[start,stop],packfit)
        # single or early time slice for plot function 
        # print('muplotfit plot_run debug: calling chi for f (None)')
        _,f,_ = chi(tf,None,None,pars)
        # print('mufitplot plot_run debug: min max f_res {}, min max f {}'.format([f_res.min(),f_res.max()],[f.min(),f.max()])) 
        if self.rotating_frame_frequencyMHz: # (t,y,f_res,ey) the last three must be transformed to rrf
                                             # (tf,f) non rebinned, the last must be transformed to rrf
            _,y,ey = rebin(self.suite.time,self.rrf_asymm,[start,stop],pack,e=self.rrf_asyme)
            f_res = mixer(t, f_res,self.rotating_frame_frequencyMHz) 
            f = mixer(tf, f,self.rotating_frame_frequencyMHz) 
        
        # print('nu = {}, chi_res = {}'.format(nu,chi_res))
        # t,fres calculated on plot points for residues


        
        # assume self.suite.single() 
        # prepare figure  
        fgroup = [self.suite.groups[k]['forward'] for k in range(len(self.suite.groups))]
        bgroup = [self.suite.groups[k]['backward'] for k in range(len(self.suite.groups))]
        alpha = [self.suite.groups[k]['alpha'] for k in range(len(self.suite.groups))]
        group = [fgroup, bgroup, alpha]
        dy_fit = (y_fit-f_fit)/ey_fit

        data = [t, y, ey, f_res, tf, f, dy_fit]

        chi_dof =[nu_fit,chi_fit]

        if self.suite.single() and not self.suite.multi_groups():
            set_figure_fit = set_single_fit
        else:
            set_figure_fit = set_sequence_fit
        # single draws a static plot, many runs are dealt by animation
        # invoke ploting of dat, fit and residues

        if early_late:
            dy_fit_early = (y_fit_early-f_fit_early)/ey_fit_early
            dy_fit_late = (y_fit_late-f_fit_late)/ey_fit_late            
            w_early  = nu_fit*ones(t_fit_early.shape[0])/nu_fit_early
            w_late = nu_fit*ones(t_fit_late.shape[0])/nu_fit_late
            data_late = [t_late,y_late,ey_late,f_late_res,tfl,fl,dy_fit_early,dy_fit_late]
            chi_dof_late = [nu_fit_early,nu_fit_late,\
                            chi_fit_early,chi_fit_late,
                            w_early,w_late]
        else:
            data_late, chi_dof_late = None, None
        # self.log('Debug-mufitplot: set_figure  = {}'.format(set_figure_fit))
        self.fig = set_figure_fit(  self.fig,self.model_name(),
                                    early_late,
                                    data,
                                    group,
                                    run_title,
                                    chi_dof,
                                    data_late,
                                    chi_dof_late,rrf=self.rotating_frame_frequencyMHz)
        return True, ''

    def chi_1(self,t,yin,eyin,pars):
        '''
        input:
            t, yin, eyin are time, asymmetry, error all 1d 
            unless multigroup_in_components(self.dashboard) is True (yin eyin 2d)
            pars list of fit parameter values for (single run or calib) 
               (fitvalues as obtained by int2min)
            kgroup index of grouping assumed 0
        output:
            number_dof, fit function (1d), chi2_r  scalar for single and calib
        '''
        from mujpy.aux.aux import multigroup_in_components, int2_multigroup_method_key
        from mujpy.aux.aux import _nparam, calib, int2_method_key
        
        # print('chi_1 mufitplot debug: {}'.format(pars)) 
        if sum(multigroup_in_components(self.dashboard)):
            pardicts = self.dashboard['userpardicts_guess']
            # print('mufitplot chi_1 Minuit pars {}'.format([[k,pardict["name"]] for k,pardict in enumerate(pardicts)]))
            parfixed = sum([1 if pardict['flag'] =='!' else 0 for pardict in pardicts]) 
            freepars = len(pardicts) - parfixed
            methods_keys = int2_multigroup_method_key(self.dashboard,self._the_model_) 
        else:        
            kgroup = 0
            _,_, freepars = _nparam(self.dashboard['model_guess'])
            methods_keys = int2_method_key(self.dashboard,self._the_model_)
        if yin is None:
            nu = None
            # print('mufitplot chi_1 debug: yin and nu are None') 
            f = self._the_model_._add_(t,*pars)
            chi2 = None
        else: # yin is not None: # data are not None
            nu = yin.size - freepars # degrees of freedom in plot
                  
             # int2_method_key() returns a list of methods and keys to calculate the components
             # int2calib__method_key() removes 'al' first parameter and renumbers accordingly the parameters
            if sum(multigroup_in_components(self.dashboard)):
                if calib(self.dashboard):
                # must be handles as a standard global multigroup plot
                    # print('mufitplot chi_1 debug: _load_data_multigroup_calib_')
                    methods_keys = methods_keys[1:]
                ok, msg = self._the_model_._load_data_multigroup_(t,
                                                                      yin,
                                                                      methods_keys,
                                                                      e=eyin)
            elif calib(self.dashboard):
                methods_keys = methods_keys[1:]
                ok, msg = self._the_model_._load_data_calib_(t,
                                                             yin,
                                                             methods_keys,
                                                             e=eyin) 
            else:
                ok, msg = self._the_model_._load_data_(t,
                                                       yin,
                                                       methods_keys,
                                                       e=eyin)             
            if ok:
                # print('mufitplot chi_1 debug: pars = {}'.format(pars))
                # self.debug(pars) # works only for global fits
                f = self._the_model_._add_(t,*pars)
                chi2 = self._the_model_._chisquare_(*pars)/nu # chi2 in plot
                if sum(multigroup_in_components(self.dashboard)):
                    chicchi = []
                    nu = nu*0.5
                    for k in range(yin.shape[0]):
                        chicchi.append(chi2)
                    chi2=chicchi
            else: 
                self.log(msg)
                return None, None, None            
        return nu,f,chi2
        
    def chi_2(self,t,yin,eyin,pars):
        '''
        input:
            t, yin, eyin are time, asymmetry, error 2d
            pars list of lists  of fit parameter values for sequence 
               (fitvalues as obtained by int2min)
            kgroup index of grouping is 0 if run sequence
            otherwise sequential or global groups 
        output:
            number_dof, fit function (2d), chi2_r (list)
        deals also with cases where only f is needed (plot function)
        #  B1 2d (runs) sequential, pars is a list, uses _load_data_
        #  A20 2d (groups) sequential, pars is a list, uses _load_data_
        # to be implemented
        #  B20 3d (group,runs) sequential, pars is a list, uses _load_data_multigroup_
        #  B21 3d multigroup userpar, to be done
        '''
        from mujpy.aux.aux import _nparam, int2_method_key, calib
        
        # print('chi_2 muplotfit debug: pars = {}'.format(pars))
        # really the next if should go 
        #if len(self.suite.grouping)==1 or self.model == "model_guess":
        _,_, freepars = _nparam(self.dashboard[self.model])
            # model_guess and model_result are lists
        #else:
        #    _,_, freepars = _nparam(self.dashboard[self.model])# [0])
            # multigroup model_results is a list of lists, both share the same number of free fit parameters 
        chi2 = []
        nu = len(t) - freepars # degrees of freedom in plot
        # print('chi_2 mufitplot debug: nu {}, freepars {}'.format(nu,freepars))
        if yin is not None: # data are not None
            mthdk = int2_method_key(self.dashboard,self._the_model_)
            if calib(self.dashboard):
                mthdk = mthdk[1:]
            ok, msg = self._the_model_._load_data_(t,yin,mthdk,e=eyin) 
            if ok:
                # print('mufitplot chi_2 debug: pars {}'.format(pars))
                f = self.fstack(t,*pars)
            else:
                self.log(msg)
                return None, None, None
            # print('mufitplot chi_2 debug: len(pars) {}'.format(len(pars)))
            for k in range(len(self.lastfits)): # works only for 2d cases
                chi2.append(self.lastfits[k].fval/nu) # chi2 in plot
                # print('mufitplot chi_2 debug: k {}, chi2 = {}'.format(k,chi2[-1]))
        else: # only f is really needed, this can be called only after a call with yin not None
            f = self.fstack(t,*pars)
            chi2 = [None for k in range(len(pars))]
        return nu,f,chi2
                
    def fstack(self,t,*pars):
        from numpy import vstack
        # print('fstack mufitplot debug: pars = {}'.format(pars))
        for k,par in enumerate(pars):
            # print('fstak mufitplot debug: par = {}'.format(par))
            fin = self._the_model_._add_(t,*par)
            if k==0: # f for histogram
               f = fin
            else:
               f = vstack((f,fin))
        return f

    def choosefftplot(self,fft_range,real):    
        '''
        distinguishes  
          single-multi run
          single-multi group
          sequential-global  
        '''
        from mujpy.aux.aux import multigroup_in_components, userlocals
        from mujpy.aux.aux import get_nruns, int2min, int2min_multigroup
        from numpy import array

        # single - multi run sequential  A1 B1 both produce list of pars 
        # single - multi group sequential A1 A20 both produce list of pars 
        # multi run sequential multi group global B2 produces list of pars
        # multi group global A21 produces par
        # multi run global C1 produces par
        # multi run multi group global C2 produces par
        if self.suite.single():
            if self.suite.multi_groups(): # A2
                # print('mufitplot choosefftplot debug: A2')
                self.log('Multigroup fft animation: toggle pause/resume by clicking on the plot')
                if sum(multigroup_in_components(self.dashboard)): # A21 single chi2
                    userpars = "userpardicts_guess" if self.guess else "userpardicts_results"
                    pardicts = self.dashboard[userpars]
                    pars,_,_,_,_,_ = int2min_multigroup(pardicts)
                    # ok = self.single_fft_plot_multi_global(fft_range,pars,real)
                else: # A20 as many chi2 as groups now results are saved in single group dashboards
                    if self.model=="model_result":
                        pars = [array(lastfit.values) for lastfit in self.lastfits]
                        # print('mufitplot choosefftplot debug: A20 fit pars = {}'.format(pars))
                    else:
                        par,_,_,_,_,_ = int2min(self.dashboard[self.model])                    
                        pars = [array(par) for k in range(len(self.suite.groups))]
                        # print('mufitplot choosefftplot debug: A20 guess pars = {}'.format(pars))
                    # ok = self.single_fft_plot_multi_sequential(fft_range,pars,real)
                asymm, asyme = self.suite.asymmetry_multigroup()
            else:  # A1 simple single plot
                pars,_,_,_,_,_ = int2min(self.dashboard[self.model])
                # ok = self.single_fft_plot(fft_range,pars,real)
                asymm, asyme = self.suite.asymmetry_single(self.suite._the_runs_[0],0)
        else: 
            if userpars(self.dashboard):
                if not self.suite.multi_groups(): # C1, userpardicts no multigroup
                    self.log('Multirun fft animation: toggle pause/resume by clicking on the plot')
                    pass# not yet
                else: 
                    self.log('Multi group and run fft animation: toggle pause/resume by clicking on the plot')
                    if userlocals(self.dashboard): # C2
                        pass# not yet
                    else: # B21 userpardicts, multigroup_in_components no tilde_incomponents
                        pars = [] 
                        # not yet
#                        for run in get_nruns(self.suite):
#                            par = 
#                            pars.append(par)
            else:
                pars = [] 
                if self.suite.multi_groups(): 
                    for run in get_nruns(self.suite): # B1 no multigroup
                        for group in self.suite.groups: # or B20 no userpardicts multigroup
                            par,_,_,_,_,_ = int2min(self.dashboard[self.model])
                            pars.append(par)
                            asymm, asyme = self.suite.asymmetry_multirun_multigroup()
                else:
                    for run in get_nruns(self.suite): # B1 no multigroup
                        par,_,_,_,_,_ = int2min(self.dashboard[self.model])
                        pars.append(par)
                    asymm, asyme = self.suite.asymmetry_multirun(0)
                # ok = self.sequential_fft_plot(fft_rangepars,real)

        # print('mufitplot choosefftplot debug: pars = {}'.format(pars))
        self.plot_fft(fft_range,pars,asymm,asyme,real)
 
    def plot_fft(self,fft_range,pars,asymm,asyme,real):
        '''
        input:
            fft_range
                start,stop,sigma MHz
            pars list (single, calib) or list of lists (sequence)
            asymm, asyme  1d (single, calib) or 2d (sequence)
        uses data as dictated by self.dashboard["fit_range"]
        
        calls either 
            self.chi_fft (only for model function) 
            set_single_fft or set_sequence_fft, from aux.plot
            first version single only
                                    (produces actual figures
                                     using aux.plot functions)
        '''
        from mujpy.aux.aux import derange, rebin, autops, ps
        from mujpy.aux.aux import get_run_title, userpars
        from mujpy.aux.plot import set_figure_fft
        from copy import deepcopy
        from numpy import ones, exp, linspace, sqrt, mean, fft
        from numpy import hstack, linspace, zeros, mgrid

        # print('plot_fft muplotfit debug: pars = {}'.format(pars))

        run_title = get_run_title(self.suite)    # always a list, even for single 
        string = 'global ' if userpars(self.dashboard) else ''  
        if len(self.suite.groups)>1:
            strgrps = [groups['forward']+'-'+groups['backward'] for groups in self.suite.groups]
        else:
            strgrps = self.suite.groups['forward']+'-'+self.suite.groups['backward']
        if self.guess:
            if len(run_title)==len(self.suite.groups):
                run_title = [runtitle + " ("+string+"guess) group" + strgrp for runtitle,strgrp in zip(run_title,strgrps)]
            else:
                run_title = [runtitle + " ("+string+"guess)" for runtitle in run_title]
        else:
            if len(run_title)==len(self.suite.groups):
                run_title = [runtitle + " ("+string+"fit) group" + strgrp for runtitle,strgrp in zip(run_title,strgrps)]
            else:
                run_title = [runtitle + " ("+string+"fit)" for runtitle in run_title]
        #############################
        # rebinning of data as in fit 
        # this works for single 
        # and for sequential
        #############################
        fittup,ermsg = derange(self.dashboard["fit_range"],self.suite.histoLength) # range as tuple
        fit_pack = 1
        if len(fittup)==3: # plot start stop pack
            fit_start, fit_stop, fit_pack = fittup[0], fittup[1], fittup[2]
        elif len(fittup)==2: # plot start stop
            fit_start, fit_stop = fittup[0], fittup[1]

        t_fit,y_fit,ey_fit = rebin(self.suite.time,asymm,[fit_start,fit_stop],fit_pack,e=asyme)

        f_fit_res, f_fit = self.chi_fft(t_fit,y_fit,*pars) # returns both partial and full model
        
        fgroup = [self.suite.groups[k]['forward'] for k in range(len(self.suite.groups))]
        bgroup = [self.suite.groups[k]['backward'] for k in range(len(self.suite.groups))]
        alpha = [self.suite.groups[k]['alpha'] for k in range(len(self.suite.groups))]
        group = [fgroup, bgroup, alpha]
        
        dt = t_fit[1]-t_fit[0]  # time bin
        fmax = 0.5/dt  # max frequancy available
        l = (fit_stop-fit_start)//fit_pack # floor division, dimension of data
        df = 1/(dt*l)
        n = 2*l # not a power of 2, but surely even
        nf = hstack((linspace(0,l,l+1,dtype=int), linspace(-l+1,-1,l-2,dtype=int)))
        dfa = 1/n/dt         # digital frequency resolution
        f = nf*dfa  # all frequencies, l+1 >=0 followed by l-1 <0

        fstart, fstop, fsigma = float(fft_range[0]), float(fft_range[1]), float(fft_range[2])
        start, stop = int(round(fstart/dfa)), int(round(fstop/dfa))
        f = deepcopy(f[start:stop]) # selected slice

# asymmetry has 1,2,3 dimensions to distinguish run from group
# here they are on same par: the data are reshaped into a 2 d array
        if len(asymm.shape)==1:
            y = zeros(n) # for data zero padded to n
            yf = zeros(n) # fit function for rephasing, zero padded to n
            ey = zeros(n)
            y[0:l] = y_fit[fit_start:fit_stop] - f_fit_res # zero padded partial residues
            yf[0:l] = f_fit
            ey[0:l] = ey_fit[fit_start:fit_stop] #  slice of time stds
            t = dt*linspace(0,n-1,n)
        elif len(asymm.shape)==2:
            # print('mufitplot plot_fft debug: nruns/ngroups {}, nbins {}'.format(asymm.shape[0],y_fit.shape[-1]))
            y = zeros((y_fit.shape[0],n)) # for data zero padded to n
            yf = zeros((y_fit.shape[0],n)) # fit function for rephasing, zero padded to n
            ey = zeros((ey_fit.shape[0],n))
            # print('mufitplot plot_fft debug: shapes y {}, asymm {}, f_fit_res {}'.format(y[:,0:l].shape,y_fit[:,fit_start:fit_stop].shape,f_fit_res.shape))
            y[:,0:l] = y_fit[:,fit_start:fit_stop]- f_fit_res # zero padded partial residues
            yf[:,0:l] = f_fit
            ey[:,0:l] = ey_fit[:,fit_start:fit_stop] #  slice of time stds
            _,t = dt*mgrid[0:asymm.shape[0],0:n]
        else: # 3 is reshaped to 2d
            nruns,ngroups = asymm.shape[0], asymm.shape[1]
            y = zeros((nruns*ngroups,n)) # for data zero padded to n
            yf = zeros((nruns*ngroups,n)) # fit function for rephasing, zero padded to n
            ey = zeros((nruns*ngroups,n))
            y[:,:,0:l] = y_fit[:,:,fit_start:fit_stop]- f_fit_res # zero padded partial residues
            yf[:,:,0:l] = f_fit # zero padded full fit function
            ey[:,:,0:l] = ey_fit[:,:,fit_start:fit_stop] #  slice of time stds
            y = y.reshape(nruns*ngroups,n)  # reshaped to 2d
            yf = yf.reshape(nruns*ngroups,n)  # reshaped to 2d
            ey = ey.reshape(nruns*ngroups,n)  # reshaped to 2d
            _,t = dt*mgrid[0:nruns*ngroups,0:n] # time for zero padded data
        # print('mufitplot plot_fft debug: t.min t.max = {},{}'.format(t.min(),t.max()))
        filter_apo = exp(-(t*fsigma)**3) # hypergaussian filter mask
                                           # is applied as if first good bin were t=0
        filter_apo = filter_apo/filter_apo.sum(axis=1)[0]/dt # approx normalization
        dfa = 1/n/dt         # digital frequency resolution
        y *= filter_apo # zero padded, filtered partial residues
        yf *= filter_apo # zero padded, filtered full fit function
        ey *= filter_apo # filteres 
        fft_e = (sqrt(mean(fft.fft(ey)**2,axis=-1))).real#.reshape(ey.shape[0],1)*ones(y.shape) # one per asymmetry
        fft_amplitude = fft.fft(y)  # amplitudes (complex), matrix with rows fft of each run
        fftf_amplitude = fft.fft(yf)  # amplitudes (complex), same for fit function
        if real:
            ########################
            # REAL PART
            # APPLY PHASE CORRECTION
            # try acme
            ########################
            if len(y_fit.shape)==1:
                fftf_amplitude[start:stop], p0, p1, out = autops(fftf_amplitude[start:stop],'acme') # fix phase on theory 
                out = out[:out.index('\n')]
                self.log('Autophase: '+out)
                fft_amplitude[start:stop] = ps(fft_amplitude[start:stop], p0=p0 , p1=p1).real 
                ap = deepcopy(fft_amplitude[start:stop].real)
                apf = deepcopy(fftf_amplitude[start:stop].real)            
            else:
                for k in range(fftf_amplitude.shape[0]):
                    fftf_amplitude[k,start:stop], p0, p1, out = autops(
                                              fftf_amplitude[k,start:stop],'acme') # fix phase on theory 
                    out = out[:out.index('\n')]
                    self.log('Autophase {}: '.format(k)+out)
                    # fft_amplitude[k,start:stop], p0, p1, out = autops(fft_amplitude[k,start:stop],'acme') # fix phase on theory 
                    fft_amplitude[k,start:stop] = ps(fft_amplitude[k,start:stop],p0=p0,p1=p1).real # same correction as theory
                ap = deepcopy(fft_amplitude[:,start:stop].real)
                apf = deepcopy(fftf_amplitude[:,start:stop].real)            
            ylabel = 'FFT Real part [arb. units]'
        else:
            ##################
            # POWER
            ##################
            if len(asymm.shape)==1:
                ap = fft_amplitude.real[start:stop]**2+fft_amplitude.imag[start:stop]**2
                apf = fftf_amplitude.real[start:stop]**2+fftf_amplitude.imag[start:stop]**2
            else:
                ap = fft_amplitude.real[:,start:stop]**2+fft_amplitude.imag[:,start:stop]**2
                apf = fftf_amplitude.real[:,start:stop]**2+fftf_amplitude.imag[:,start:stop]**2
            ylabel = 'FFT Power [arb. units]'
        self.fig_fft = set_figure_fft(self.fig_fft,self.model_name(),ylabel,
                                    f,
                                    ap,
                                    apf,
                                    fft_e,
                                    group,
                                    run_title)
#        self._the_model_._include_all_() # usual _the_model_ mode: all components included

        
    def chi_fft(self,t_fit,y_fit,*pars):
        '''
        returns f for the fft of residues
        for the moment works only for single or sequential fits
        '''
        from mujpy.aux.aux import int2_method_key
        from numpy import vstack
        method_keys = int2_method_key(self.dashboard,self._the_model_)
        self._the_model_._load_data_(t_fit,y_fit,method_keys)
        fft_include_components = []      
        for component in self.dashboard["model_guess"]:
            if "fft" in component.keys():
                fft_include_components.append(component["fft"]) # conditionally subtract component
            else:
                fft_include_components.append('False') # fft of data if fft is not in model dict
        self._the_model_._fft_init(fft_include_components) # sets 
        # distinguish list from list of lists 
        # this applies to A20 single run multigroup sequential fits
        # print('chi_fft mufitplot debug: pars = {}'.format(pars))
        fres = self._the_model_._add_fft_(t_fit,y_fit,*pars)
        if isinstance(pars,tuple):
            # print('mufitplot chi_fft debug: pars is tuple')
            f = self.fstack(t_fit,*pars)
        else: 
            # print('mufitplot chi_fft debug: pars is not tuple')
            self._the_model_._add_(t_fit,*pars)
        # debug
#        from matplotlib.pyplot import draw, subplots 
#        fig,ax = subplots()
#        for k in range(f.shape[0]):
#            ax.plot(t_fit,f[k,:])
#            ax.plot(t_fit,fres[k,:],alpha=0.5)
#        draw()    
        return fres, f
    
    def model_name(self):
        '''
        output:
            model name (e.g. 'mlbg')
        '''
        model = self.dashboard["model_guess"]
        return ''.join([component["name"] for component in model]) 
        
    def single_chi(self):
        '''
        output:
            True if chi_1 is required (single cost function)
            False if chi_2 is required (multi cost finctions)
        '''
        #       multi_groups suite.single userpars multigroup_in_components userlocals
        # A1     False          True        False       False               False *
        # A20    True           True        False       False               False
        # A21    True           True        True        True                False *
        # B1     False          False       False       False               False
        # B20    True           False       False       False               False
        # B21    True           False       True        True                False
        # C1     False          False       True        False               True  *
        # C2     True           False       True        True                True  *
        #  True : A1 1d, A21 2d multigroup userpar, C1 2d userpar, C2 3d multgrup userpar
        #  False: B1 2d, A20, B20 2d multigroup, B21 2d multigroup userpar
        from mujpy.aux.aux import userlocals, userpars
        return (userlocals(self.dashboard) or 
                   (self.suite.single() and 
                        (not self.suite.multi_groups() or userpars(self.dashboard)))) 
                    
    def debug(self,pars):
        # can be deleted after debugging stage of mujpy 2.0
        from mujpy.aux.aux import min2int_multigroup
        names, ps, _ = min2int_multigroup(self.dashboard,pars,pars)
        print('muplotfit debug: internal dashboard parameters {}'.format([[name,p] for name,p in  zip(names,ps)]))
