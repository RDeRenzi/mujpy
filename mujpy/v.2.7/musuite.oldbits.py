###############
# Deprecated
###############

    def single_multigroup_for_back_counts(self,runs,groupings):
        """
        Deprecated

        * input: 
        *         runs, runs-to-add
        *         grouping, {'forward':[3],'backward':[4]] for 3-4
        *         uses self.single_for_back_counts
        * output:
        *         yf, yb, eyf, eyf = vstacks over groups
        *         2D numpy arrays, used in calib multigroups
        """

        from numpy import vstack,array
        for k,grouping in enumerate(groupings):
            yforw, ybackw, ey_forw, ey_backw = self.single_for_back_counts(runs,grouping)
            #        all are 1D numpy arrays
            if k:
                yf = vstack((yf,yforw))
                yb = vstack((yb,ybackw))
                eyf = vstack((eyf,ey_forw))
                eyb = vstack((eyb,ey_backw))
            else:
                yf = yforw
                yb = ybackw
                eyf = ey_forw
                eyb = ey_backw
        return yf,yb,eyf,eyb

    def multirun_multigroup_for_back_counts(self,runs,groupings):
        """
        Deprecated

        * input: 
        *         runs, list of list, 
                        [[run0 runs to add],[run1 runs to add], ...] 
        *         grouping, {'forward':[3],'backward':[4]] for 3-4
        *         uses self.single_multigroup_for_back_counts
        * output:
        *         yf, yb, eyf, eyf = vstacks over runs, groups
        *         3D numpy arrays, used in calib multirun multigroups
        """

        from numpy import vstack,array
        for k,run in enumerate(runs):
            yforw, ybackw, ey_forw, ey_backw = self.single_multigroup_for_back_counts(run,groupings)
            #        all are 1D numpy arrays
            if k:
                yf = vstack((yf,array([yforw])))
                yb = vstack((yb,array([ybackw])))
                eyf = vstack((eyf,array([ey_forw])))
                eyb = vstack((eyb,array([ey_backw])))
            else:
                yf = array([yforw])
                yb = array([ybackw])
                eyf = array([ey_forw])
                eyb = array([ey_backw])
        return yf,yb,eyf,eyb
                        
    def asymmetry_multirun(self,kgroup):
        """
        Deprecated

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
        Deprecated

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
 
    def asymmetry_multirun_multigroup(self,multirun,multigroup): # 
        """
        Deprecated

        input: 
            multirun True/False 
            multigroup True/False
            calls self.asymmetry_single which calls self.single_for_back_counts
        outputs: 
            asymmetry and asymmetry error (3d,2d/1d)
            for run in runs:
                for group in groups:
                    np.vstack # axis=1
                np.vstack # axis=0
        A1 -> False False self.asymmetry_single
        A20 -> False True vstack of ngroups, iterate kgroup 
        A21 -> False True vstack of ngroups
        B1 -> True False vstack of nruns, iterate krun
        B20 ->  True True vstack of vstacks, iterate krun. kgroup
        B21 -> True True vstack of vstacks, iterate krun
        C1 -> True False vstack of nruns
        C2 -> True True vstack of vstacks
        """

        from numpy import array, vstack

        if self.loadfirst:
            if multirun and multigroup: # 3d
            
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
            elif multirun: # 2d
                for krun,run in enumerate(self._the_runs_):
                    a,e = self.asymmetry_single(run,0)
                    if krun:
                        asymm, asyme  = vstack((asymm,a)), vstack((asyme,e)) # runs are vstacked in 2nd dimension
                    else: # krun=0
                        asymm, asyme = a,e # 1nd dimension bins
            elif multigroup:
                for kgroup in range(len(self.groups)):
                    a,e = self.asymmetry_single(self._the_runs_[0],kgroup)
                    if kgroup:
                        asymm,asyme = vstack((asymm,a)), vstack((asyme,e)) # groups are vstacked in 2nd dimension
                    else: # kgroup = 0
                        asymm,asyme = a,e# 1d dimension bins
            else: # 1d, only bins
                asymm,asyme = self.asymmetry_single(self._the_runs_[0],0)
            return asymm, asyme
        else:
            self.console('** CHECK ACCESS to database (or load runs first)') 
            return None, None

