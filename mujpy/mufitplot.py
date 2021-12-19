class mufitplot(object):
    '''
    plotfit class 
    produces the fitplot
    multiple runs sequential produce anim
    '''
# a bit of discussion to be removed when done
# plot needs to access a version of mumodel, 
# with different ranges for rebinned data, 
# denser function plot, original data
# mumodel() is needed for generating model function
# and calculating chi2.
# everything can be obtained from dashboard 
# essentially by replicating what mufit does
    def __init__(self,plot_range,the_suite,dashboard_file,guess=False,fig=None,rotating_frame_frequencyMHz = 0.0):
        '''
        input: 
            either
                suite instance 
                dashboard best fit json file 
            or 
               fit instance
            fig fig_handle
        '''
        from mujpy.mucomponents.mucomponents import mumodel
        import json
        
        if the_suite:
            if dashboard_file==None: # dashboard is missing
                the_suite.console('Cannot plot: best fit json file missing')
                ok = False
            else: # everything else fine
                self._the_model_ = mumodel()                            
                self.suite = the_suite
                if not self.suite.loadfirst: 
                    self.suite.console('Sorry, no access to data... quitting mufitplot')
                    ok = False
                else:
                    ok = self.get_dashboards(dashboard_file)
        else: # the_suite was missing
                print('Cannot plot: suite instance missing')
                ok = False
        self.fig = fig  
        if ok:
            self.rotating_frame_frequencyMHz = rotating_frame_frequencyMHz
            self.guess = guess
            self.model = "model_guess" 
            if not self.guess:
                if "model_result" in self.dashboard[0]:
                    self.model = "model_result"
                elif "userpardicts_guess" in self.dashboard[0].keys():   
                    self.model = "model_result"  
                else:            
                    self.suite.console('Sorry no fit results yet, plotting guess instead')
                    self.guess = True
            # print('__init__ mufitplot debug: model = {}, guess = {}'.format(self.model,self.guess))    
           
            self.chooseplot(plot_range)

    def get_nruns(self):
        '''
        get nrun strings
        '''
        nruns = []
        print
        for k,run in enumerate(self.suite._the_runs_):
            nruns.append(str(run[0].get_runNumber_int()))
        return nruns
        
    def get_dashboards(self,file):
        '''
        load dashboard from json files 
        into the list self.dashboard 
        '''
        import json
        prefile = ''
        self.dashboard = []
        runs = self.get_nruns()
        if len(runs)==1:
            self.suite.console('Plotting best fit of run {}'.format(runs[0]))
        for run in runs:
            n = file.find(run) 
            if n>0:
                if len(runs)>1:
                    self.suite.console('Plotting best fit of runs {}'.format(runs))
                    self.suite.console('Animation: toggle pause/resume by clicking on the plot')
                m =len(run)
                prefile = file[:n]
                postfile = file[n+m:]
                break  

        if prefile:                
            for nrun in runs:
                dash_file = prefile+nrun+postfile
                try:
                    with open(dash_file,'r') as f:
                        dash = json.load(f)
                        self.dashboard.append(dash) 
                        # self.suite.console('get_dashboard mufitplot debug: appended {} to dashboard list'.format(dash_file))
                except Exception as e:
                    self.suite.console('Exception: {}'.format(e))
                    self.suite.console('Cannot json.load(open({}))'.format(dash_file))
                    return False
            # print('get_dashboard mufitplot debug: first dash result {}'.format(self.dashboard[0]['model_result']))
        else:
            self.suite.console('Warning: run list {} mismatch with dashboard file {}'.format(runs,dash_file))
        return True	

    def get_run_title(self):
        '''
        write run_title strings
        '''
        from mujpy.aux.aux import get_title
        run_title = []
        if self.suite.single():
            run = self.suite._the_runs_[0]
            if self.multigroup:
                for k in range(len(self.suite.grouping)):
                    run_title.append(str(run[0].get_runNumber_int())+'-'+get_title(run[0]))
            else:
                run_title.append(str(run[0].get_runNumber_int())+'-'+get_title(run[0]))
        else:
            for k,run in enumerate(self.suite._the_runs_):
                run_title.append(str(run[0].get_runNumber_int())+'-'+get_title(run[0]))
        return run_title
        
              
    def chooseplot(self,plot_range):
        '''
            switch for single (A1), sequential (B1),
            ... 
        '''
        from mujpy.aux.aux import multigroup_in_components
        if self.suite.single():
            if self.suite.multi_groups():
                self.suite.console('Multigroup animation: toggle pause/resume by clicking on the plot')
                if sum(multigroup_in_components(self.dashboard[0])): # single chi2
                    ok = self.single_plot_multi(plot_range)
                else:                          # as many chi2 as groups
                    ok = self.single_plot_multi_sequential(plot_range)
            else:                              # simple single plot
                ok = self.single_plot(plot_range)
        else:
            ok = self.sequential_plot(plot_range)       
        if not ok:
            self.suite.console('Exiting mufitplot without a plot')
    
    def single_plot(self,plot_range):
        '''
        input plot_range, passed to self.plot_run together with
            pars, list of one list of fit parameter values
            asymm, asyme 1d        
        calib pops first par and passes it as alpha to standard plot
        '''
        from mujpy.aux.aux import int2min, mixer
        from numpy import cos, pi
        kgroup = 0 # default single group 
        pars,_,_,_,_ = int2min(self.dashboard[0][self.model])
        if self.calib():
            self.suite.grouping[kgroup]['alpha'] = pars[0] # from fit parameter to standard asymmetry mode
        asymm, asyme = self.suite.asymmetry_single(self.suite._the_runs_[0],0)
        if self.rotating_frame_frequencyMHz:
            self.rrf_asymm = mixer(self.suite.time,asymm,self.rotating_frame_frequencyMHz)
            self.rrf_asyme = mixer(self.suite.time,asyme,self.rotating_frame_frequencyMHz)
        #self.suite.console('mufitplot: Inside single plot; debug mode')
        self.multigroup = False
        return self.plot_run(plot_range,pars,asymm,asyme)        

    def single_plot_multi(self,plot_range):
        '''
        input plot_range, passed to self.plot_run together with
            pars, list of one list of fit parameter values
            asymm, asyme 1d        
        '''
        from mujpy.aux.aux import int2min_multigroup, mixer
        from numpy import cos, pi, vstack
        userpardicts = (self.dashboard[0]["userpardicts_guess"] if self.guess else 
                        self.dashboard[0]["userpardicts_result"])
        pardict = self.dashboard[0]["model_guess"][0]["pardicts"][0]
        pars,_,_,_,_ = int2min_multigroup(userpardicts)
        p = pars
        if self.calib():
            for kgroup,group in enumerate(self.suite.grouping):
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
        self.multigroup = True
        #self.suite.console('mufitplot: Inside single plot; debug mode')
        return self.plot_run(plot_range,pars,asymm,asyme)        
    
    def single_plot_multi_sequential(self,plot_range):
        '''
        input plot_range, passed to self.plot_run together with
            pars, list of lists  of fit parameter values
            asymm, asyme 2d
        if model_result
            dashboard contain a "model_result" key which is a list of lists
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
                    time = vstack((time,self.suite.time))
            self.rrf_asymm = mixer(time,asymm,self.rotating_frame_frequencyMHz)
            self.rrf_asyme = mixer(time,asyme,self.rotating_frame_frequencyMHz)
        pars = []
        # self.suite.console('mufitplot: Inside sequential plot; debug mode')
        for kgroup in range(asymm.shape[0]):
            if self.model == "model_result":
                #print('single_plot_multi_sequential mufitplot debug: self model = {}',format(self.dashboard[0][self.model][kgroup]))
                values,_,_,_,_ = int2min(self.dashboard[0][self.model][kgroup])                
            else:
                #print('single_plot_multi_sequential mufitplot debug: self model = {}',format(self.dashboard[0][self.model]))
                values,_,_,_,_ = int2min(self.dashboard[0][self.model])                
            pars.append(values) # here pars is list of list!!
        self.multigroup = True   
        return self.plot_run(plot_range,pars,asymm,asyme)        

    def sequential_plot(self, plot_range):
        '''
        input plot_range, passed to self.plot_run together with
            pars, list of lists  of fit parameter values
            asymm, asyme 2d
        '''
        from mujpy.aux.aux import int2min, mixer
        from numpy import cos, pi, vstack
        
        kgroup = 0 # default single 
        # dashboard must become a suite thing: each sequential fit has its own
        pars = []
        # self.suite.console('mufitplot: Inside sequential plot; debug mode')
        for k,run in enumerate(self.suite._the_runs_):
            values,_,_,_,_ = int2min(self.dashboard[k][self.model])
            
            #print('plot_run muplotfit debug: pars = {}'.format(pars)) run {},  values = {}'.format(run[0].get_runNumber_int(),values))
            pars.append(values)
        asymm, asyme = self.suite.asymmetry_multirun(kgroup) # 
        self.multigroup = False
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
        Draws either a 2x2 or 2x3 subplots figure 
        using aux.plot functions 
        '''
        from mujpy.aux.aux import derange, rebin, multigroup_in_components, userpars, mixer
        from mujpy.aux.plot import set_single_fit, set_sequence_fit
        from iminuit import Minuit
        from numpy import ones 

        # print('plot_run muplotfit debug: pars = {}'.format(pars))
        # chi_1 is for a single run single group (kgroup = 0) asymm 1d
        # chi_2 is for 2d asymm either sequential run or multigroup
        chi = (self.chi_1 if (self.suite.single() and not self.multigroup) or 
                    sum(multigroup_in_components(self.dashboard[0])) else self.chi_2)
        # print('plot_run mufitplot debug: single not multigroup {}, multigroup in comp {}'.format(self.suite.single() and not self.multigroup,bool(sum(multigroup_in_components(self.dashboard[0])))))

        run_title = self.get_run_title()    # always a list, even for single 
        string = 'global ' if userpars(self.dashboard[0]) else ''  
                  
        if self.guess:
            run_title = [title + ": "+string+"guess values" for title in run_title]
        else:
            run_title = [title + ": "+string+"fit results" for title in run_title]
        plottup = derange(plot_range,self.suite.histoLength)
        #############################
        # rebinning of data as in fit 
        # this works for single 
        # and for sequential
        #############################
        fittup = derange(self.dashboard[0]["fit_range"],self.suite.histoLength) # range as tuple
        fit_pack = 1
        if len(fittup)==3: # plot start stop pack
            fit_start, fit_stop, fit_pack = fittup[0], fittup[1], fittup[2]
        elif len(fittup)==2: # plot start stop
            fit_start, fit_stop = fittup[0], fittup[1]

        # load modules and reproduce fit
        t_fit,y_fit,ey_fit = rebin(self.suite.time,asymm,[fit_start,fit_stop],fit_pack,e=asyme)
        # single slices as in fit 
        nu_fit, f_fit, chi_fit = chi(t_fit,y_fit,ey_fit,pars)
        if self.rotating_frame_frequencyMHz: # (t_fit,y_fit,f_fit,ey_fit) the last three must be transformed to rrf
            _,y_fit,ey_fit = rebin(self.suite.time,self.rrf_asymm,[fit_start,fit_stop],fit_pack,e=self.rrf_asyme)
            f_fit = mixer(t_fit,f_fit,self.rotating_frame_frequencyMHz)
        if not nu_fit:
            return False
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
            print('mufitplot plot_run debug: y_fit_late ey_fit_late f_fit_late shape = {},  {},  {}'.format(y_fit_late.shape, ey_fit_late.shape, f_fit_late.shape))

            if not nu_fit_late:
                return False
            # function as in fit

            t_fit_early,y_fit_early,ey_fit_early = rebin(self.suite.time,asymm,        
                                                [fit_start,fit_late_start],fit_pack,e=asyme)
            # single slices as in fit 
            nu_fit_early, f_fit_early, chi_fit_early = chi(t_fit_early,y_fit_early,ey_fit_early,pars)
            if self.rotating_frame_frequencyMHz:
                _,y_fit_early,ey_fit_early = rebin(self.suite.time,self.rrf_asymm,[fit_start,fit_late_start],fit_pack,e=self.rrf_asyme)
                f_fit_early = mixer(t_fit_early, f_fit_early,self.rotating_frame_frequencyMHz) 

            if not nu_fit_early:
                return False
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
#        self.suite.console('plot_range= {}'.format(plot_range))
#        self.suite.console('start, stop, pack = {},{},{}'.format(start, stop, pack))

        t,y,ey = rebin(self.suite.time,asymm,[start,stop],pack,e=asyme)
        # rebinned single or early slices 
        nudum, f_res, _ = chi(t,y,ey,pars)
        packfit = int(pack/2)
        tf,dum = rebin(self.suite.time,asymm,[start,stop],packfit)
        # single or early time slice for plot function 
        _,f,_ = chi(tf,None,None,pars)
        if self.rotating_frame_frequencyMHz: # (t,y,f_res,ey) the last three must be transformed to rrf
                                             # (tf,f) non rebinned, the last must be transformed to rrf
            _,y,ey = rebin(self.suite.time,self.rrf_asymm,[start,stop],pack,e=self.rrf_asyme)
            f_res = mixer(t, f_res,self.rotating_frame_frequencyMHz) 
            f = mixer(tf, f,self.rotating_frame_frequencyMHz) 
        
        if not nudum:
            return False
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

        if self.suite.single() and not self.multigroup:
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
        # self.suite.console('Debug-mufitplot: set_figure  = {}'.format(set_figure_fit))
        self.fig = set_figure_fit(  self.fig,self.model_name(),
                                    early_late,
                                    data,
                                    group,
                                    run_title,
                                    chi_dof,
                                    data_late,
                                    chi_dof_late,rrf=self.rotating_frame_frequencyMHz)
        return True

    def chi_1(self,t,yin,eyin,pars):
        '''
        input:
            t, yin, eyin are time, asymmetry, error all 1d 
            unless multigroup_in_components(self.dashboard[0]) is True (yin eyin 2d)
            pars list of fit parameter values for (single run or calib) 
               (fitvalues as obtained by int2min)
            kgroup index of grouping assumed 0
        output:
            number_dof, fit function (1d), chi2_r  scalar for single and calib
        '''
        from mujpy.aux.aux import _nparam, int2_method_key, int2_calib_method_key      
        from mujpy.aux.aux import multigroup_in_components, int2_multigroup_method_key   
        
        # print('chi_1 mufitplot debug: {}'.format(pars)) 
        if sum(multigroup_in_components(self.dashboard[0])):
            pardicts = self.dashboard[0]['userpardicts_guess']
            parfixed = sum([1 if pardict['flag'] =='!' else 0 for pardict in pardicts]) 
            freepars = len(pardicts) - parfixed
            methods_keys = int2_multigroup_method_key(self.dashboard[0],self._the_model_) 
        else:         
            kgroup = 0
            _,_, freepars = _nparam(self.dashboard[0]['model_guess'])
            methods = int2_method_key(self.dashboard[0],self._the_model_)
        if yin is None:
            nu = None
        else:
            nu = yin.size - freepars # degrees of freedom in plot
                  
             # int2_method_key() returns a list of methods and keys to calculate the components
             # int2calib__method_key() removes 'al' first parameter and renumbers accordingly the parameters
        if yin is not None: # data are not None
        # (int2_calib_method_key(self.dashboard[0],self._the_model_) if self.calib() else
            if sum(multigroup_in_components(self.dashboard[0])):
                if self.calib():
                    # print('mufitplot chi_1 debug: y.shape = {}, e.shape = {}'.format(eyin.shape,yin.shape))
                    ok, msg = self._the_model_._load_data_multigroup_calib_(t,yin,methods_keys,e=eyin) 
                else:
                    ok, msg = self._the_model_._load_data_multigroup_(t,yin,methods_keys,e=eyin)
            elif self.calib():
                ok, msg = self._the_model_._load_data_calib_(t,yin,methods,self.suite.grouping[kgroup]['alpha'],e=eyin) 
            else:
                ok, msg = self._the_model_._load_data_(t,yin,methods,self.suite.grouping[kgroup]['alpha'],e=eyin)             
            if ok:
                f = self._the_model_._add_(t,*pars)
                chi2 = self._the_model_._chisquare_(*pars)/nu # chi2 in plot
                if sum(multigroup_in_components(self.dashboard[0])):
                    chicchi = []
                    nu = nu*0.5
                    for k in range(yin.shape[0]):
                        chicchi.append(chi2)
                    chi2=chicchi
            else: 
                self.suite.console(msg)
                return None, None, None            
        else: # only f is really needed
            f = self._the_model_._add_(t,*pars)
            chi2 = None
        return nu,f,chi2
        
    def chi_2(self,t,yin,eyin,pars):
        '''
        input:
            t, y, ey are time, asymmetry, error 2d
            pars list of lists  of fit parameter values for sequence 
               (fitvalues as obtained by int2min)
            kgroup index of grouping is 0 if run sequence
            otherwise sequential or global groups 
        output:
            number_dof, fit function (2d), chi2_r (list)
        deals also with cases where only f is needed (plot function)
        '''
        from mujpy.aux.aux import _nparam, int2_method_key
        
        # print('chi_2 muplotfit debug: pars = {}'.format(pars)) 
        if len(self.suite.grouping)==1 or self.model == "model_guess":
            _,_, freepars = _nparam(self.dashboard[0][self.model])
            # model_guess and model_result are lists
        else:
            _,_, freepars = _nparam(self.dashboard[0][self.model][0])
            # multigroup model_results is a list of lists, both share the same number of free fit parameters 
        chi2 = []
        nu = len(t) - freepars # degrees of freedom in plot
        if yin is not None: # data are not None
            for k,(y,ey) in enumerate(zip(yin,eyin)): # single and calib have k = 0
                kgroup = 0 if len(self.suite.grouping)==1 else k
                # always group 0 multirun if not multigroup, sequential groups otherwise
                # print('chi_2 mufitplot debug: group {}: t:{}-{}mus, alpha {:.4f}'.format(kgroup,t[0],t[-1],self.suite.groups[kgroup]['alpha']))
                if self.calib():
                    ok, msg = self._the_model_._load_data_calib_(t,y,
                                         int2_method_key(self.dashboard[0],self._the_model_),
                                         self.suite.groups[kgroup]['alpha'],e=ey) 
                                         # int2_method_key() returns a list of 
                                         # methods and keys to calculate the components
                else:
                    ok, msg = self._the_model_._load_data_(t,y,
                                         int2_method_key(self.dashboard[0],self._the_model_),
                                         self.suite.groups[kgroup]['alpha'],e=ey) 
                                         # int2_method_key() returns a list of 
                                         # methods and keys to calculate the components
                if ok:
                    f = self.fstack(t,*pars)
                    chi2.append(self._the_model_._chisquare_(*pars[k])/nu) # chi2 in plot
                else:
                    self.suite.console(msg)
                    return None, None, None
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

    def calib(self):
        '''
        True if the first component is 'al'
        '''
        # print(self.dashboard)
        return self.dashboard[0]['model_guess'][0]['name']=='al'
        
    def model_name(self):
        '''
        output:
            model name (e.g. 'mlbg')
        '''
        model = self.dashboard[0]["model_guess"]
        return ''.join([component["name"] for component in model])    
