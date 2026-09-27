    def _chisquare_single_(self,*argv,k=0,l=None):
        """
        Deprecated

        input:
            argv ar single run single group fit parameters
            k[, l] are indices of _y_ and _e_ multidimensional arrays
            if l==None k = either krun or kgroup 
        Used in mufit.prepare_csv_row and (maybe?) mufitplot
        Provides partial chisquares over individual runs or groups
        """ 

       # print('_chisquare_ mucomponents debug: {} {} {}'.format(self._x_.shape,self._y_.shape,self._e_.shape))
#        self._ndata = self._y_[k,:].size
        if l is None:
            
            return sum(  ( (self._add_(self._x_,*argv) - self._y_[k,:]) /self._e_[k,:])**2 )
        else:
            return sum(  ( (self._add_(self._x_,*argv) - self._y_[k,l,:]) /self._e_[k,l,:])**2 )
            
    from iminuit import Minuit as _M
    _chisquare_.errordef = _M.LEAST_SQUARES

# ----- old stuff -------------------------
 
    def _rebin_calib_(self,*argv):
        """
        Deprecated
        sets _alpha, calls _asymmetry_, rebins(_start_,_stop_,_pack_) into _x_, _y_, _e_

        input: 
            *argv   passed as a variable number of parameter values 
                    alpha,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [alpha,val1,val2,val3,val4,val5, ...]
        """

        from mujpy.tools.tools import rebin, set_alpha
        from numpy import set_printoptions
        set_printoptions(precision=2)
        print('mucomponents _rebin_calib_ argv = {}'.format(argv))
        self._alpha = set_alpha(argv,self._slice_nruns_,self._slice_ngroups_)
        y,e = self._asymmetry_()
        return rebin(self._t_,y,[self._start_,self._stop_],self._pack_,e=e)


    def _slice_(self,krun,kgroup):
        """
        Deprecated, suite.slice_for_back_counts is used directly
        was used to select a mufitplot data slice (self. implied for _)
        assumed alpha are already set
        input:
            krun, kgroup both aither 0 or -1
            0 this dimension is singleton
            -1 this dimension spans the suite 
        output:
            loads slices of _a_, _ea_, by self._asymmetry_ if _calib 
        """
        from mujpy.tools.tools import slice, rebin
        if self._calib:
            a,eq = self._asymmetry_plot_(
        self._a_,self._ea_ = slice(self._asymmetry_(),krun.kgroup) if self._calib else slice(self._aa_,self._eaa_,krun,kgroup)


    def _load_calib_(self,suite,returntup,components,multigroup=False,multirun=False):
        """
        Deprecated
        works for A1, A20, A21, B1, B20, B21, C1, C2 calib
        fit with alpha as free parameter
        input:
            (not homogeneous with _load_, mufit.load_ tests calib)
            suite is the musuite instance that loads the data 
            returntup is a tuple of integer bins start, stop, pack
            (_reload_calib_(returntup) 
                            sets alpha
                            calculates _asymmetry_ 
                            invokes tools rebin)
            components is a list [[method,[key...]],...,[method,[key...]]], 
                produced by the appropriate mujpy.tools.tools int2_method() 
        """
        from mujpy.tools.tools import rebin
        from mujpy.tools.tools import list_depth
        # self._suite_ = suite # remove
        self._t_ = suite.time  # pristine time array, unbinned
        krun, kgroup = 0, 0 # for the single run and/or single group cases
        self._ntruecomponents_ = len(components) # remove? still used in fft
        # this blockdefined self._lambdakeys_ and self._keys for use in self._add_
        key_depth = list_depth(components) # 3,4,5 respectively, for single, multirg, or multir_multig
        # _add_ loops for mthd, keys in zip(self._methods_,self._keys_): 
        #                 pars = keys(p,keys)
        #                 f += methd(t,*pars)
        self._methods_ = [mthd_keys[0] for mthd_keys in components[1:]] # extracts methods, skip 'al' 
        keys_ = [mthd_keys[1] for mthd_keys in components[1:]] # extracts corresponding keys
        self._lambdakeys_ = []
        self._keys_ = []
        for keys in keys_:
            lkeys = lambda p,kys: [key(p) for key in kys] if key_depth==3 else (
                    lambda p,kys: [[key(p) for key in rgkey] for rgkey in kys]) if key_depth==4 else (
                    lambda p,kys: [[[key(p) for key in gkey] for gkey in rkey] for rkey in kys]) # key_depth==4
            self._lambdakeys_.append(lkeys)
            self._keys_.append(keys)

        # wrappers
        self._run, self._grouping = (suite._the_runs_,suite.grouping) if (
                multigroup and multirun) else (suite._the_runs_,suite.grouping[kgroup]) if (
                               multirun) else (suite._the_runs_[krun],suite.grouping) if (
                multigroup             ) else (suite._the_runs_[krun],suite.grouping[kgroup])

        self._for_back_counts = suite.multirun_multigroup_for_back_count if (
                multigroup and multirun) else suite.multirun_for_back_counts if (
                               multirun) else suite.single_multigroup_for_back_counts if (
                multigroup             ) else suite.single_for_back_counts
        self._reload_calib_(returntup)
        self._ndata = self._yf_.size

        self._add_ = self._add_calib_
        self._nruns_, self._ngroups_ = (self._yf_.shape[0],self._yf_.shape[1]) if multirun and multigroup else (
                                        self._yf_.shape[0],                1)  if                multirun else (
                                                        1,self._yf_.shape[0])  if              multigroup else (
                                                        1,                1)  

    def _reload_calib_(self,returntup):
        """
        Deprecated
        mucomponent _load_  defines a suite 
        this method redefines range, returntup = (start, stop, pack)
        for wrapper self._for_back_counts_ in the following calib fits 
            A1 single    
            A20 multigroup sequential
            A21 multigroup global
            B1  multirun sequential
            B20 multirun multigroup sequential
            B21 multirun sequential multigroup global 
            C1 multirun global
            C2 multirun multigroup global 
        """
        from mujpy.tools.tools import rebin
        self._include_all_() # to rcover from possible fft mode
        self._start_, self._stop_, self._pack_ = returntup
        self._yf_,self._yb_,self._eyf_,self._eyb_ = self._for_back_counts(self._run,self._grouping)
        y,e = self._asymmetry_()
        self._x_,self._y_,self._e_ = rebin(self._t_,y,[self._start_,self._stop_],self._pack_,e=e) # rebinned time array
        #print('********* entered mucomponents _reload_calib_ self._x_.shape,self._y_.shape {},{}'.format(self._x_.shape,self._y_.shape))



    def _add_single_(self,x,*argv): 
        """
        Deprecated
            *argv   passed as a variable number of parameter values 
                    val1,val2,val3,val4,val5,val6, ... at this iteration 
                    argv is a list of values [val1,val2,val3,val4,val5,val6, ...]
        _add_single_ DISTRIBUTES THESE PARAMETER VALUES::
              asymmetry fit with fixed alpha
              order driven by model e.g. blml
        """      
        f = zeros_like(x)  # initialize a 1D array
        p = argv
        for j in range(self._n0truecomponents_,self._ntruecomponents_): 
            component = self._components_[j][0]
            keys = self._components_[j][1] 
            pars = [key(p) for key in keys] 
            f += component(x,*pars)  
        return f

    def _add_multigroup_(self,x,*argv):   
        """
        Deprecated
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    val0,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [val0,val1,val2,val3,val4,val5, ...]
        _add_multigroup_ DISTRIBUTES THESE PARAMETER VALUES::
              asymmetry fit with fixed alpha
              order driven by 
              first global parameters 
              then local run parameters, first local and then "~","!" pars in model, e.g. mgbl 
        method is vectorized as a vstack over groups, whose number n = y.shape[0]
        and produce a n-valued np.array function f, f[k,:] for y[k,:],e[k,:] 
        """      
        # from numpy import savez,set_printoptions
        # set_printoptions(precision = 2)
        f = zeros((self._y_.shape[0],x.shape[0]))  # initialize a 2D array shape (groups,bins)   
        p = argv 
        
        for method, keys in self._components_[self._n0truecomponents_:]: # works also for calib [self._n0truecomponents_ skips alphas]
            pars = [[key(p) for key in groups_key] for groups_key in keys]
            f += method(x,*pars)  
        return f
          
    def _add_multirun_(self,x,*argv):   
        """
        Deprecated
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    val0,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [val0,val1,val2,val3,val4,val5, ...]
        _add_multirun_ DISTRIBUTES THESE PARAMETER VALUES::
              asymmetry fit with fixed alpha
              order driven by 
              first global parameters 
              then local run parameters, first local and then "~","!" pars in model, e.g. mgbl 
        method is vectorized as a vstack over runs, whose number n = y.shape[0]
        and produce a n-valued np.array function f, f[k,:] for y[k,:],e[k,:] 
        """    
        f = zeros((self._y_.shape[0],x.shape[0]))  # initialize a 2D array shape (groups,bins)   
        p = argv 
        for method, keys in self._components_: #
            pars = [[key(p) for key in groups_key] for groups_key in keys]
            f += method(x,*pars) 
        return f

    def _add_multirun_multigroup_(self,x,*argv):   
        """
        Deprecated
         input: 
            x       time array
            *argv   minuit pars in multirun multigroup fit
        (mgbl models: called directly
         almgbl models call _add_calib_multirun_multigroup_
            that EXTRACTS alpha VALUES FROM argv)
            and then _add_
        output: asymmetry fit model
        """
        # from numpy import savez,set_printoptions
        # set_printoptions(precision = 2)
        f = zeros((self._y_.shape[0],self._y_.shape[1],x.shape[0]))  # initialize a 3D array with shape (runs,groups,bins)   
        p = argv
        #print('debug mucomponent _add_multirun_multigroup_ n0 {}'.format(self._n0truecomponents_))
        for method, keys in self._components_[self._n0truecomponents_:]: # works also for calib [self._n0truecomponents_ skips alphas]
            if key_depth==1:
                pars = [key(p) for key in keys]
            elif key_depth==2:
                pars = [[key(p) for key in run_or_group_key] for run_or_group_key in keys]
            else:
                pars = [[[key(p) for key in group_key] for group_key in run_key] for run_key in keys]
         #   ff =  method(x,*pars)
         #   print('debug mucomponent _add_multirun_multigroup_ method {}, keys.shape {} f.shape {} ff.shape {}'.format(method, array(keys).shape,f.shape,ff.shape))
            f += method(x,*pars)
        return f

    def _add_calib_multigroup_(self,x,*argv):   
        """
        Deprecated
         input: 
            x       time array
            *argv   list of minuit parameters in al two group fits
                    par[0], par[1], ... are alpha VALUE
            CALCULATES asymmetry fit with fixed alpha
            then rebins it and stores in self._y_ self._e_
            CALLS appropriate _add_(t,*argv) for the model 
        """      
        from mujpy.tools.tools import set_alpha, rebin
        self._alpha = set_alpha(argv,self._nruns_,self._ngroups_)
        y,e = self._asymmetry_()
        t,self._y_,self._e_ = rebin(x,y,[self._start_,self._stop_],self._pack_,e=e)
        return self._add_multigroup_(t,*argv)

    def _load_multirun_user_(self,x,y,components,e=1):
        """
        Deprecated
        input: 
            x, y, e are numpy arrays, y, e are 2d 
            e = 1 or missing yields unitary errors 
            components is a list [[method,[key,...,key]],...,[method,[key,...,key]]], 
                produced by int2_multirun_user_method_key() from mujpy.tools.tools
                where method is an instantiation of a component, e.g. self.ml 
                and value = eval(key) produces the parameter value
            _add_multirun_ must produce a 2d function f.shape(ngroup,nbins)
            therefore _components_ must be a np.vectorize of nrun copies of the method 
            method da not allowed here, no need for alpha
            no fft of residues   
        """
        self._x_ = x
        self._y_ = y        # self._global_ = true if _nglobals_ is not none else false
        self._ndata = y.size
        self._components_ = components
        self._ntruecomponents_ = len(components)
        self._add_ = self._add_multirun_
        self._add_norebin_ = self._add_multirun_
        # print('mucomponents _load_multirun_user_ mucomponents debug: {} components, n of parameters/component: {}'.format(len(components),[len(par) for group in components for par in group[1]]))
        try:
            if isinstance(e,int):
                self._e_ = ones((y.shape))
            else:
                if len(y.shape)>1:
                    # print('_load_multirun_user_ mucomponents debug: x,y,e not e=1')
                    if e.shape!=y.shape or x.shape[0]!=y.shape[-1]:
                        # print('_load_multirun_user_ mucomponents debug: x,y,e different shape[0]>1')
                        raise valueerror('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                       y.shape,
                                                                                       e.shape))          
                elif e.shape!=y.shape or e.shape[0]!=x.shape[0]:
                    # print('_load_multirun_user_ mucomponents debug: x,y,e different shape[0]=1')
                    raise valueerror('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                           y.shape,
                                                                                           e.shape))          
            # print('_load_multigroup_ mucomponents debug: defining self._e_')
            self._e_ = e
            # print('mucomponents _load_multigroup_ debug: self._x_ {}, self._y_ {}, self._e_  {}shape'.format(self._x_,self._y_,self._e_) 
        except ValueError as e:
            return False, e       
        return True, '' # no error

    def _add_calib_multirun_multigroup_(self,x,*argv):   
        """
         Deprecated
         input: 
            x       time array
            *argv   list of minuit parameters in al two group fits
                    par[0], par[1], ... are alpha VALUE
            CALCULATES unbinned _asymmetry_ with fixed alpha
            rebins and stores in self._x_ self._y_ self._e_
            self._x_ is [only] used in mufitplot
            CALLED by _add_(t,*argv)
            _asymmetry_ works for any model defined in _load_calib_
        """      
        from mujpy.tools.tools import set_alpha, rebin
        self._alpha = set_alpha(argv,self._nruns_,self._ngroups_)
        y,e = self._asymmetry_() # works for multirun and/or multigroup when alphoa = set_alpha
        self._x_,self._y_,self._e_ = rebin(x,y,[self._start_,self._stop_],self._pack_,e=e)
        return self._add_multirun_multigroup_(t,*argv)



    def _load_multigroup_(self,asymm,asyme,returntup,components): #self,x,y,,e=1): 

        """
        Deprecated, remove
        input: 
            suite
            components is a list [[method,[[key,...],...,[key,...]]],...,[method,[[key...],...,[key,...]]]], 
                produced by int2_method() from mujpy.tools.tools
                where method is  e.g. self.ml 
                and value = eval(key) produces one parameter value
            _add_multigroup_ must produce a 2d function f.shape(ngroup,nbins)
            therefore _components_ must be a np.vectorize of ngroup copies of the method 
        """
        from mujpy.tools.tools import rebin
        self._start_, self._stop_, self._pack_ = returntup
        self._x_, self._y_,self._e_ = rebin(self._suite_.time,asymm,[self._start_,self._stop_],self._pack_,e=asyme) 
        self._components_ = components
        self._ntruecomponents_ = len(components)
        self._ndata = self._y_.size

        self._add_ = self._add_multigroup_
        self._add_norebin_ = self._add_multigroup_

        self._include_all_() # to rcover from possible fft mode

        return True, '' # no error remove?

    def _load_calib_multigroup_(self,suite,returntup,krun,kgroup,components): 
        """
        deprecated, remove
        fit with alpha as free parameter
        input: 
            x, yf, ybm eyfm eyb  are numpy arrays
            returntup = [start,stop],pack  for rebinning of asymmetry, created in _add_calib_multigroup_
        #         A = (yf-alpha*yb)/(yf+alpha*yb)
        #        eA = 2*alpha/(self._yfc_ - alpha*self._ybc_)**2 * sqrt((self._ybc_*self._eyfc_)**2 + (self._yfc_*self._eybc_)**2) 
        # must be calculated in _add_calib_multigroup_ with the minuit value of fit parameter alpha
            components is a list [[method,[[key,...],...,[key,...]]],...,[method,[[key...],...,[key,...]]]], 
                produced by int2_method() from mujpy.tools.tools
                where method is an instantiation of a component, e.g. self.ml 
                and value = eval(key) produces the parameter value (the inner key list is for different groups)
            e = True, default, for caluclasting eA, provide e=False for simple plot
            self._add_ = self._add_calib_multigroup_ 
            no fft of residues
        """
         # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
        # self._ntruecomponents_ = number of components apart from dalpha        
        from mujpy.tools.tools import rebin
        self._suite_ = suite
        self._start_, self._stop_, self._pack_ = returntup
        self._yf_,self._yb_,self._eyf_,self._eyb_ = self._suite_.single_multigroup_for_back_counts(
                                            self._suite_._the_runs_[krun],self._suite_.grouping) 
        self._t_ = self._suite_.time  # pristine time array, unbinned 

        self._x_,_ = rebin(self._t_,self._yf_,[self._start_,self._stop_],self._pack_)
        self._components_ = components
        self._ntruecomponents_ = len(components)
        self._n0truecomponents_ = 1
        self._ndata = self._yf_.size

        self._add_ = self._add_calib_multigroup_
        self._add_norebin_ = self._add_calib_multigroup_norebin_
        # print('debug mucomponents._load_calib_multigroup_: components = {}'.format(self._components_[self._n0truecomponents_:]))
#        for k, val in enumerate(components):
#            self._ntruecomponents_ += 1
#            # print('mucomponents _load_calib_multigroup_ debug: val = {}'.format(val))
#            self._components_.append(val) # store again [method, [key,...,key]], ismin
        return True, '' # no error remove?

    def _load_multirun_multigroup_():
        """
         deprecated, remove
         see calib below

        """

    def _load_calib_multirun_multigroup_(self,suite,returntup,krun,kgroup,components):
        """
        deprecated, remove
        input
            suite
            loas self._t_, self._x_, self._yf_, self._yb_, self._eyf_, self._eyb_ 
            returntup = [start,stop],pack  stored in their self versions
            components = [method,[[[key,...],...],[[key,...],...],...],method,[[[key...],...],[[key,...],...]],...] 
                produced by int2_multirun_multigroup_method_key() mujpy.tools.tools
            stored in self._components_
        """
        from mujpy.tools.tools import rebin
        self.suite = suite
        self._start_, self._stop_, self._pack_ = returntup
        self._yf_,self._yb_,self._eyf_,self._eyb_ = self.suite.multirun_multigroup_for_back_counts(
                                                        self.suite._the_runs_,self.suite.grouping)
        self._t_ = self.suite.time  # pristine time array, unbinned

        self._x_,_ = rebin(self._t_,self._yf_,[self._start_,self._stop_],self._pack_) # rebinned time array
        self._components_ = components
        self._ntruecomponents_ = len(components)
        self._n0truecomponents_ = 1
        self._ndata = self._yf_.size

        self._add_ = self._add_calib_multirun_multigroup_
        self._add_norebin_ = self._add_calib_multirun_multigroup_norebin_
        self._nruns_, self._, _ = self._yf_.shape 
        # print('debug mucomponents._load_calib_multigroup_: components = {}'.format(self._components_[self._n0truecomponents_:]))
#        for k, val in enumerate(components):
#            self._ntruecomponents_ += 1
#            # print('mucomponents _load_calib_multigroup_ debug: val = {}'.format(val))
#            self._components_.append(val) # store again [method, [key,...,key]], ismin
        return True, '' # no error, remove?
