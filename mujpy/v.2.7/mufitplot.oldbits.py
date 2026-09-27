
##############
# Deprecated
##############
                              
    def reload(self,tup):
        """
        Deprecated
        reloads data through mufit.mumodel and mufit.suite

            utilizing self.pars, self.dashboard, self.suite
        into self.x, self.y, self.e 
        that overloads mucomponents 
            self._x_, self._y_, self._e_, required by the cost function
        """

        from mujpy.tools.tools import set_alpha, rebin, rshp
        if self.calib: 
            #self.fit._alpha = set_alpha(self.pars,len(self.suite.runs),len(self.suite.groups))
            self.model._reload_calib_(tup)
        else:
            self.model._reload_(tup) # 
        self.x, self.y, self.e = rshp(self.model._x_), reshaplot(self.model._y_), reshaplot(self.model._e_)
 
    def single_chi(self):
        """
        unused
        output:
            True if chi_1 is required (single cost function in global mufitplot: A1,A21,C1,C2 adn their calib)
            False if chi_2 is required (multi cost finctions in sequential mufitplots: A20,B1,B20,B21 and their calib)
        """
        #       multi_groups suite.single userpars function_multi_in_components locals  single_cost_function
        # A1     False          True        False       False               False        *
        # A20    True           True        False       False               False
        # A21    True           True        True        True                False        *
        # B1     False          False       False       False               False
        # B20    True           False       False       False               False
        # B21    True           False       True        True                False
        # C1     False          False       True        False               True         *
        # C2     True           False       True        True                True         *
        #  True : A1 1d, A21 2d multigroup userpar, C1 2d userpar, C2 3d multgrup userpar
        #  False: B1 2d, A20, B20 2d multigroup, B21 2d multigroup userpar
#        from mujpy.tools.tools import locals
#        return (locals(self.dashboard) or 
#                   (self.suite.single() and 
#                        (not self.suite.multi_groups() or userpars(self.dashboard)))) 
        return self.fit.A1() or self.fit.A21() or self.fit.A1_calib() or self.fit.A21_calib() or self.fit.C1() or self.fit.C1_calib() or self.fit.C2() or self.fit.C2_calib()

    def chooseplot(self,plot_range):
        """
        deprecated
            switch for single (A1), sequential (B1),
            ... 
        #       multi_groups suite.single userpars function_multi_in_components :w

        # A1     False          True        False       False               False *
        # A20    True           True        False       False               False
        # A21    True           True        True        True                False *
        # B1     False          False       False       False               False
        # B20    True           False       False       False               False does not exist
        # B21    True           False       True        True                False
        # C1     False          False       True        False               True  *
        # C2     True           False       True        True                True  *
        """
        from mujpy.tools.tools import function_multi_in_components

        A1, A20, A21, B1, B20, B21, C1, C2 = self.fit_types()
        if A1:
            ok, msg = self.plot_singlerun(plot_range)
        elif A20 or A21:
            self.log('Multigroup fit animation: toggle pause/resume by clicking on the plot')
            ok, msg = self.plot_singlerun_multigroup(plot_range)
        else:
            self.log('Multirun fit animation: toggle pause/resume by clicking on the plot')
            if B1():
                ok, msg = self.plot_multirun_singlegroup_sequential(plot_range)
            elif B21():
                ok, msg = self.plot_multirun_multigroup_userpar(plot_range)
            elif C1():
                self.log('No C1 yet. Exiting mufitplot without a plot')
                ok, msg = False, 'C1'
            else: #  C2()
                ok, msg = self.plot_multirun_multigroup_userpar(plot_range)
       
        if not ok:
            self.log('Exiting mufitplot without a plot: '+msg)
    
    def plot_singlerun(self,plot_range):
        """
        input plot_range a csv string, start, stop [,pack [start_r, stop_r]]
              obtains asymm, asyme = self.the_model._asymmetry_(alpha)
        A1
        """
        from mujpy.tools.tools import int2min, mixer, calib
        from numpy import cos, pi
#        self.log('muplotfit plot_singlerun debug')
        kgroup = 0 # default single group
        pars = self.lastfit.values
        #pars,_,_,_,_,_ = int2min(self.dashboard,self.suite.runs)
        if calib(self.dashboard): # A1 calib
            self.suite.grouping[kgroup]['alpha'] = pars[0] 
            self.suite.groups[kgroup]['alpha'] = pars[0]
            asymm, asyme = self.the_model._asymmetry_(pars[0])
        else: # A1
            asymm, asyme = self.suite.asymmetry_single(self.suite._the_runs_[0],0)
        if self.rotating_frame_frequencyMHz:
            self.rrf_asymm = mixer(self.suite.time,asymm,self.rotating_frame_frequencyMHz)
            self.rrf_asyme = mixer(self.suite.time,asyme,self.rotating_frame_frequencyMHz)
        #self.log('mufitplot: Inside single plot; debug mode')
        return self.plot_run(plot_range,pars,asymm,asyme)        

    def plot_singlerun_multigroup(self,plot_range):
        """
        remove!
        method chosen to standardize, all done
        input plot_range, passed to self.plot_run together with
            pars, list of one list of fit parameter values
            asymm, asyme 2d # A20, A20_calib, A21, A21 calib
        calib retrieves alpha as [par[0],par[1],...] and
              obtains asymm, asyme = self.the_model._asymmetry_(alpha)
        A20 A21
        """
        from mujpy.tools.tools import mixer, calib, int2min_multigroup, function_multi_in_components, int2_method_key, int2_multigroup_method_key, set_alpha, rebin
        from numpy import vstack, array
        
       

    def plot_multirun_multigroup_userpar(self,plot_range):
        """
         inputs: 
            plot_range
        # C2 fit 
        # asymm, asyme are 3d
        """
        from mujpy.tools.tools import int2min_multirun_multigroup, mixer, calib
        from numpy import cos, pi, vstack
#        userpardicts = (self.dashboard["userpardicts_guess"] if self.guess else 
#                        self.dashboard["userpardicts_result"])
#        pardict = self.dashboard["model_guess"][0]["pardicts"][0]
        pars,_,_,_,_,_ = int2min_multirun_multigroup(self.dashboard,self.suite.runs)
# no calib mode!
#        p = pars
#        if calib(self.dashboard):
#            for kgroup,group in enumerate(self.suite.groups):
#                group['alpha'] = eval(pardict["function_multi"][kgroup])
#                self.suite.groups[kgroup]["alpha"] = eval(pardict["function_multi"][kgroup])
        asymm, asyme = self.suite.asymmetry_multirun_multigroup()
        if self.rotating_frame_frequencyMHz:
            for krun in range(asymm.shapcalibe[0]):
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


    def plot_multirun_singlegroup_sequential(self, plot_range):
        """
        input plot_range, passed to self.plot_run together with
            pars, list of lists  of fit parameter values
            asymm, asyme 2d
            # B1 fit
        """
        from mujpy.tools.tools import mixer
        from numpy import vstack
        
        kgroup = 0 # default single 
        # dashboard must become a suite thing: each sequential fit has its own
        pars = []
        # self.log('mufitplot: Inside sequential plot; debug mode')
        for k,lastfit in enumerate(self.fit.lastfits): # loops over runs
            values = lastfit.values   # list of values for a run
            
            #print('plot_run muplotfit debug: pars = {}'.format(pars)) run {},  values = {}'.format(run[0].get_runNumber_int(),values))
            pars.append(values) # list of lists one per run
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
  
    def plot_multirun_singlegroup_user(self,plot_range):  
        """
        input plot_range, passed to self.plot_run together with
            pars, list of lists of fit parameter values
            asymm, asyme 2d
            # C1 fit (presented to plot_run like A21 singlerun_multigroup_userpar)
            the difference is that A21 all Minuit pars are userpardicts, 
            while here there is a mixed situation
        """
        from mujpy.tools.tools import mixer
        from numpy import vstack
        
        kgroup = 0 # default
        # dashboard must become a suite thing: each sequential fit has its own
        pars = list(self.fit.lastfit.values)  # self.lastfit.values is the damn ValueView
        # self.log('mufitplot: Inside sequential plot; debug mode')
        # special arrangement for C1 fits plotted here as B1 fits
        # self.fit.lastfits[0] is guess, self.fit.lastfits[1] is result
        # both are Minuit instances with data loaded as multirun sequential
        asymm, asyme = self.suite.asymmetry_multirun(kgroup) # 
        if self.rotating_frame_frequencyMHz

