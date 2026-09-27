    def load_(self,returntup,values,krun,kgroup,mr,mg):
        """
        Deprecated! self.dofit_ calls directly self.the_model._load_

        loads slice of unbinned data into self.the_model

        input:
            returntup = (start, stop, pack)
            values = Minuit start guess, required only if self.calib()
            krun, kgroup = indices in suite run/group slice if >=0
                            -1 is understood as [:]
        """

        from mujpy.tools.tools import set_alpha

# try removing args mr,mg and using instead the last two output values of suite.slice and suite.asymmetry_slice

        if self.calib():
            #yf,yb,eyf,eyb,multirun,multigroup = self.suite.slice_for_back_counts(krun,kgroup)
            #print('mufit load_ mr = {}'.format(mr))
            self.the_model._alpha = set_alpha(values,len(self.suite.runs),len(self.suite.grouping))
            self.the_model._load_calib_(self.suite,returntup,self.methods_keys,multirun=mr,multigroup=mg) # was yf,yb,eyf,eyb
            # self.the_model._asymmetry_() # loads data into self.the_model._y_, self.the_model._e_
        else:
            asymm,asyme,multirun,multigroup = self.suite.asymmetry_slice(krun,kgroup)
            self.the_model._load_(self.suite.time,asymm,asyme,returntup,self.methods_keys,multirun=mr,multigroup=mg)

    def write_multirun_glob_csv(self,file_csv,scan=None):
        """
        Deprecated

        reuse for writing glob csv
        this is a one-shot csv write after a global fit
        input :
            the_runs is suite _the_runs_
            file_csv = full path/filename to csv file 
        """

        from mujpy.tools.tools import get_title, min2int_multirun
        from datetime import datetime

    # prepare_csv_row writes a header: # column-index-name 
    # plos a list of rows, a row per each run, composed of 
    # run number T eT B 
    # local parameters and their errors (columns)
    # global parameters and their errors (columns with repeated values)
    # chi2 end their errors (partial), chi2 global and its error repeated) 
        names, values, errors = min2int_multirun(self.dashboard,
							            self.lastfit.values,self.lastfit.errors,self.suite._the_runs_)
        header, rows = self.prepare_csv_row()    
        with open(file_csv,'w') as f_out:  
            f_out.write(header)               
            for line in rows:
                f_out.write(line)
        file_csv = file_csv[file_csv.rfind('/')+1:]
        nrun0 = self.suite._the_runs_[0][0].get_runNumber_int()
        nrun1 = self.suite._the_runs_[-1][0].get_runNumber_int()
        self.log('Global fit of runs {}-{}:'.format(nrun0,nrun1)+
                ' values and errors saved in {}'.format(file_csv))


    def dofit_singlerun_singlegroup(self,returntup):  
        """
        returntup is a tuple of integers, start stop pack
        performs fit on single run, single group
        (A1) tested
        """
        from mujpy.tools.tools import int2min, int2_method_key, rebin

        krun, kgroup = 0, 0 
        a,e = self.suite.asymmetry_single(self.suite._the_runs_[0],0) 
        start, stop, pack = returntup
        time,asymm,asyme = rebin(self.suite.time,a,[start,stop],pack,e=e)
        #dt = time[1]-time[0]

        values,errors,fixed,limits,names, pospar = int2min(self.dashoard,self.suite.runs)
        self.methods_keys = int2_method_key(self.dashboard,self.the_model)
        self.the_model._load_(self.suite,returntup,methods_keys) # this is multirun = False, multigroup = False -> single fit

        cost = self.the_model._chisquare_
        summary = self.summary
        savefit = self.save_fit
        self.execute_log_save_fit(cost,
                                  names,
                                  values,
                                  errors,
                                  limits,
                                  fixed,
                                  pospar,
                                  summary,
                                  start,
                                  stop,
                                  krun,
                                  kgroup,
                                  savefit)

    def dofit_calib_singlerun_singlegroup(self,kgroup,returntup):
        """
        returntup is a tuple of integers, start stop pack
        performs calib fit on single run, single group 
        (A1-calib) tested
        input 
            kgroup is group index in suitegrouping
        """
        from mujpy.tools.tools import int2min, int2_method_key 

        krun, kgroup = 0, 0
        yf,yb,eyf,eyb = self.suite.single_for_back_counts(self.suite._the_runs_[0],self.suite.grouping[kgroup]) 
        start, stop, _ = returntup
        #dt = self.suite.time[1]-self.suite.time[0]

        values,errors,fixed,limits,names,pospar = int2min(self.dashboard,self.suite.runs)
        self.the_model._load_calib_(self.suite,returntup,krun,kgroup,
                                                  int2_method_key(self.dashboard,self.the_model))
        
        cost = self.the_model._chisquare_calib_
        summary = self.summary
        savefit = self.save_fit
        self.execute_log_save_fit(cost,
                                  names,
                                  values,
                                  errors,
                                  limits,
                                  fixed,
                                  pospar,
                                  summary,
                                  start,
                                  stop,
                                  krun,
                                  kgroup,
                                  savefit)
               
    def dofit_singlerun_multigroup_sequential(self,returntup):
        """
        returntup is a tuple of integers, start stop pack
        performs fit on single run, multi-group data sequentially
        (A20) tested
        """
        from iminuit import Minuit
        from mujpy.tools.tools import int2min, int2_method_key, rebin        
        
        krun = 0  #  single run!!
        a,e = self.suite.asymmetry_multigroup() # the second dimension is group
        start, stop, pack = returntup
        time,asymm,asyme = rebin(self.suite.time,a,[start,stop],pack,e=e)
        
        values,errors,fixed,limits,names,pospar = int2min(self.dashboard,self.suite.runs)

        #for kgroup,(a,e) in enumerate(zip(asymm,asyme)):
        self.methods_keys = int2_method_key(self.dashboard,self.the_model)
        self.the_model._load_(self.suite,returntup,methods_keys) # this is multirun = False, multigroup = False -> single fit 
        # self.methods_keys if for single fit, _add_ is for single fit, need a clever way to iterate groups in mucomponents
        # i.e. know which group and reload asymmetry slice
        cost = self.the_model._chisquare_
        summary = self.summary
        savefit = self.save_fit
        self.execute_log_save_fit(cost,
                                  names,
                                  values,
                                  errors,
                                  limits,
                                  fixed,
                                  pospar,
                                  summary,
                                  start,
                                  stop,
                                  krun,
                                  kgroup,
                                  savefit)

    def dofit_calib_singlerun_multigroup_sequential(self,returntup):
        """
        performs calib fit on single run, multiple groups sequentially
        returntup is a tuple of integers, start stop pack
        (A20-calib) tested
        """
        from mujpy.tools.tools import int2min, int2_method_key

        krun = 0 # single run
        #dt = self.suite.time[1]-self.suite.time[0]
        for kgroup,group in enumerate(self.suite.grouping):
            yf,yb,eyf,eyb = self.suite.single_for_back_counts(self.suite._the_runs_[0],group) 
            start, stop, pack = returntup

            values,errors,fixed,limits,names,pospar = int2min(self.dashboard,self.suite.runs)
            self.the_model._load_calib_(self.suite,returntup,krun,kgroup,
                                                  int2_method_key(self.dashboard,self.the_model))

            cost = self.the_model._chisquare_calib_
            summary = self.summary
            savefit = self.save_fit
            self.execute_log_save_fit(cost,
                                      names,
                                      values,
                                      errors,
                                      limits,
                                      fixed,
                                      pospar,
                                      summary,
                                      start,
                                      stop,
                                      krun,
                                      kgroup,
                                      savefit)

    def dofit_singlerun_multigroup_globpardicts(self,returntup):
        """
        returntup is a tuple of integers, start stop pack
        performs fit on single run, global multi-group data
        (A21) tested
        All minuit parameters predefined as globpardicts
        All component parameters assigned by functions to globpardicts
        (absence of omponents' parameters "flag":"~" identifies this fit)
        """
        from mujpy.tools.tools import rebin, int2min_multigroup, int2_multigroup_method_key
        
        krun, kgroup = 0, 0  #  single run, group ignored for A21 in execute_log_save_fit
        a,e = self.suite.asymmetry_multigroup() # first dim group, last bins
        start, stop, pack = returntup
        time,asymm,asyme = rebin(self.suite.time,a,[start,stop],pack,e=e) 
        #dt = self.suite.time[1]-self.suite.time[0]

        values,errors,fixed,limits,names,pospar = int2min_multigroup(self.dashboard["globpardicts_guess"])
        
        self.methods_keys = int2_multigroup_method_key(self.dashboard,self.the_model)
        self.the_model._load_(self.suite,returntup,multigroup=True)
        cost = self.the_model._chisquare_
        summary = self.summary_global
        savefit = self.save_fit_multigroup
        self.execute_log_save_fit(cost,
                                  names,
                                  values,
                                  errors,
                                  limits,
                                  fixed,
                                  pospar,
                                  summary,
                                  start,
                                  stop,
                                  krun,
                                  kgroup,
                                  savefit)
 
    def dofit_calib_singlerun_multigroup_globpardicts(self,returntup):
        """
        returntup is a tuple of integers, start stop pack
        performs calib fit on single run, multiple groups global
        (A21_calib) tested (non optimized fit model mufit works)
        """
        from mujpy.tools.tools import int2min_multigroup, int2_multigroup_method_key 
        
        krun, kgroup = 0, 0  #  single run, group ignored for A21 in execute_log_save_fit
        start, stop, pack = returntup
        yf,yb,eyf,eyb = self.suite.single_multigroup_for_back_counts(self.suite._the_runs_[0],self.suite.grouping) 
        #dt = self.suite.time[1]-self.suite.time[0]

        values,errors,fixed,limits,names, pospar = int2min_multigroup(self.dashboard["globpardicts_guess"])

        self.methods_keys = int2_multigroup_method_key(self.dashboard,self.the_model) 
#        p = [ 0.13,0.14,0.3,0.2,0.3,10.1,30,0.1,0.7]
#        pars = [k(p) for m,keys in methods_keys for key in keys for k in key]
#        print('dofit cal mg glob debug p = {}'.format(p))
#        print('dofit cal mg glob debug pars = {}'.format(pars))
        ok,errmsg = self.the_model._load_calib_multigroup_(self.suite,returntup,krun,kgroup,methods_keys)  #self.suite.time,yf,yb,eyf,eyb,returntup,) 
        if not ok:
            self.log('Error in _load_multigroup_: '+errmsg)
            self.log('mufit stops here')            
            return

        cost = self.the_model._chisquare_calib_
        summary = self.summary_global
        savefit = self.save_fit_multigroup
        self.execute_log_save_fit(cost,
                                  names,
                                  values,
                                  errors,
                                  limits,
                                  fixed,
                                  pospar,
                                  summary,
                                  start,
                                  stop,
                                  krun,
                                  kgroup,
                                  savefit)

    def dofit_multirun_singlegroup_sequential(self,returntup):
        """
        returntup is a tuple of integers, start stop pack
        performs sequential fit on many-run, single-group data
        (B1) tested
        """
        #from iminuit import Minuit
        from mujpy.tools.tools import int2min, int2_method_key, rebin

        #a, e = self.suite.asymmetry_multirun(0)
        # a, e are 2d: (run,timebin) 
        #start, stop, pack = returntup
        #time,asymms,asymes = rebin(self.suite.time,a,[start,stop],pack,e=e)

        values,errors,fixed,limits,names,pospar = int2min(self.dashboard,self.suite.runs)

        cost = self.the_model._chisquare_
        summary = self.summary_sequential
        savefit = self.save_fit
        kgroup = 0
        krun = -1
        #for asymm, asyme in zip(asymms,asymes): 
        self.methods_keys = int2_method_key(self.dashboard,self.the_model)
        self.the_model._load_(self.suite,returntup,methods_keys) 

        self.execute_log_save_fit(cost,
                                    names,
                                    values,
                                    errors,
                                    limits,
                                    fixed,
                                    pospar,
                                    summary,
                                    start,
                                    stop,
                                    krun,
                                    kgroup,
                                    savefit)

    def dofit_multiruns_sequential_multigroup_globpardicts(self,returntup):
        """
        returntup is a tuple of integers, start stop pack
        performs fit on sequential mani-run, global multi-group data
        (B2) testing, not fully converted
        """
        # from iminuit import Minuit
        from mujpy.tools.tools import rebin, int2min_multigroup, int2_multigroup_method_key 
        # from mujpy.tools.tools import stringify_groups, write_csv, version_flag
        # from numpy import where, array, finfo, sqrt
        
        # from matplotlib.pyplot import subplots, draw 

        # print('dofit_multirun_singlegroup_sequential mufit debug')
        # self.log('In sequential single')   
        a, e = self.suite.asymmetry_multirun_multigroup() # runs to loaded, group index
        # a, e are 2d: (run,timebin) 
        # print('mufit dofit_multiruns_sequential_multigroup_globpardicts debug: shape asymm, asyme = {}, {}'.format(a.shape,e.shape))
        start, stop, pack = returntup
        time,asymmrg,asymerg = rebin(self.suite.time,a,[start,stop],pack,e=e)
        #dt = self.suite.time[1]-self.suite.time[0]

        #zer = array(where(asymerg<2e-162))
        # time (1d): (timebin)    asymms, asymes (2d): (run,timebin) 
        values,errors,fixed,limits,names,pospar = int2min_multigroup(
                                            self.dashboard["globpardicts_guess"])

#        for k in range(len(fitvalues)):
#            self.log('{} = {}, step = {}, fix = {}, limits ({},{})'.format(names[k],values[k],errors[k],fixed[k],limits[k][0],limits[k][1]))
        # print('mufit dofit_multiruns_sequential_multigroup_globpardicts debug: Minuit inputs')
        #j = -1
        #for ns,vs,es,fx,lm in zip(names,values,errors,fixed,limits):
        #    j +=1
        #    print('{} {} = {}({}), {}, {} '.format(j,ns,vs,es,fx,lm))
        
        self.methods_keys = int2_multigroup_method_key(self.dashboard,self.the_model)
        # print('mufit dofit_multiruns_sequential_multigroup_globpardicts debug: methods_keys contains {} methods with{} keys/method'.format(len(methods_keys),[len(c) for g in methods_keys for c in g[1]]))
        #krun = -1
        
        
        #if self.dofit:
        #    fig,ax = subplots()
        #    da, ms, lw = 0.2, 0.1, 0.3
        
        #for asymm, asyme in zip(asymmrg,asymerg): 
        #    krun += 1
            #for kg in range(asyme.shape[0]):
            #    
            #    if array(where(asyme[kg,:]==0)).sum():
            #        print('mufit dofit_multiruns_sequential_multigroup_globpardicts debug: check asyme[{},{}] contains zeros!'.format(krun,kg))
            #        print('mufit dofit_multiruns_sequential_multigroup_globpardicts debug: asymm.shape {} '.format(asymm.shape))
            # asymm is 2d (group, bins)
        self.the_model._load_(self.suite,returntup,multigroup=True) 
                                    # int2_int() returns a list of methods to calculate the components remove?

#            if self.dofit:
#                fs = self.the_model._add_(time,*values)
#                print('mufit dofit_multiruns_sequential_multigroup_globpardicts debug: fs.shape  {} '.format(fs.shape)) 
#                kk, line, fmt,  = -1, ['b-','g-'],['r.','m.']
#                for a,e,f in zip(asymm,asyme,fs):
#                    kk+=1
#                    ax.errorbar(time,a+krun*da,yerr=e,fmt=fmt[kk],ms=ms,alpha=0.3)
#                    ax.plot(time,f+krun*da,line[kk],lw=lw,alpha=0.8)

        self.lastfit = Minuit(self.the_model._chisquare_,
                          name=names,
                          *values)                                        
        self.lastfit.errors = errors
        self.lastfit.limits = limits
        self.lastfit.fixed = fixed
            # self.freepars = self.lastfit.nfit
        self.number_dof = asymm.size - self.lastfit.nfit
            # print('mufit dofit_multiruns_sequential_multigroup_globpardicts debug: name value error limits fixed {}'.format([[name,value,error,limit,fix] for name,value,error,limit,fix in zip(names,values,errors,limits,fixed)]))
        if self.dofit:            
#                print('mufit dofit_multiruns_sequential_multigroup_globpardicts debug: limits {}'.format(self.lastfit.limits))
            self.lastfit.migrad()
            # check if some parameters are positive parity 
            if pospar:
                for k in pospar:
                    self.lastfit.limits[k] = [None,None]                    
                self.lastfit.migrad()
                self.lastfit.hesse()
            self.lastfits.append(self.lastfit) #  muplotfit compatibility with multifits
        # write summary on console and log
            self.summary_global(start,stop,time[1]-time[0],krun)
            print('mufit dofit_multiruns_sequential_multigroup_globpardicts debug: fval {}, ndof {}'.format(self.lastfit.fval,self.number_dof))

            version = self.dashboard["version"]+'_'+version_flag(self)
            strgrp = stringify_groups(self.suite.groups)
            modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
            file_csv = self.suite.__csvpath__+modelname+'.'+strgrp+'.'+version+'.csv'
            the_run = self.suite._the_runs_[krun][0]
            filespec = self.suite.datafile[-3:]
            header, row = self.prepare_csv_row(krun=krun)
            string1, string2 = write_csv(header,row,the_run,file_csv,filespec,scan=self.scan)
            self.log(string1)
                #self.log(string2)
            self.save_fit_multigroup(krun,string2)

            if (chain):
                values = self.lastfit.values
                
#        if self.dofit:                
#            ax.set_xlim(0,4)
#            ax.set_ylim(-0.5,2.7)
#            draw()

    def dofit_multirun_singlegroup_globpardicts(self,returntup):          
        """
        broken, must be standardized
        returntup is a tuple of integers, start stop pack
        performs global fit of many-run single-group data, not tested
        (C1) WIP, strategy:
        globpardicts is a list of glob parameter dictionaries
            each composed of keys 
            "name", "value", "error" (step), "limits" (default [None,None]), "pospar" (positive parity)
                     , "local" (default False) is DEPRECATED, use flag = "#" in model parameter instead 
                     daughter parameters as there are runs in the suite (musrfit-style)
                     'local' is used only in save_fit_multirun, which is not in the latest DONE python3 scripts check
        component parameters can be 
            equal to a global parameter or a function of global parameters
            equal to a previous local parameter or a function of global and local parameters
                both cases do not introduce a new minuit parameter and are dealt with by functions
            active, therefore local 
                i.e. the parent component parameter generates automatically as many 
                daughter parameters as there are runs in the suite
        """
        from iminuit import Minuit
        from mujpy.mucomponents.mucomponents import mumodel
        from mujpy.tools.tools import int2min_multirun, int2_multirun_glob_method_key 
        from mujpy.tools.tools import int2_multirun_grad_method_key
        from mujpy.tools.tools import minglobal2sequential, int2min, int2_method_key, version_flag
        from mujpy.tools.tools import rebin, stringify_groups #, _available_gradients_
        from numpy import array
        from time import time as timeit 

        kgroup = 0 # only one group
        a,e = self.suite.asymmetry_multirun(kgroup) # the second dimension is run 
        start, stop, pack = returntup
        time,asymm,asyme = rebin(self.suite.time,a,[start,stop],pack,e=e)

        values_in,errors,fixed,limits,names,pospar = int2min_multirun(self.dashboard,self.suite.runs)

        string = []
        method_key = int2_multirun_glob_method_key(self.dashboard,self.the_model,self.suite.nruns)
        ok, errmsg = self.the_model._load_multirun_glob_(
                                    time,asymm,method_key,e=asyme) 
        #gradmthd_key = int2_multirun_grad_mthdkey(self.dashboard,self.the_model,self.suite.nruns)
        if not ok:
            self.log(repr(errmsg))
            return
        # print('debug mufit dofit_multirun_singlegroup_globpardicts: names =\n{}\npospar =\n{}\nvalues_in =\n{}'.format(names,pospar,values_in))
        if self.grad:
            self.the_model._load_multirun_grad_(int2_multirun_grad_method_key(self.dashboard,self.the_model,self.suite.nruns))
            self.lastfit = minuit(self.the_model._chisquare_,
                              name=names,
                              grad = self.the_model._add_multirun_grad_,                          
                              *values_in
                                )                                        
        else:
            self.lastfit = minuit(self.the_model._chisquare_,
                              name=names,
                              *values_in)                                        
        self.lastfit.print_level = 0
        self.lastfit.errors = errors
        self.lastfit.limits = limits
        self.lastfit.fixed = fixed
        self.number_dof = asymm.size - self.lastfit.nfit
        
        if self.dofit:  # do the fit
            tic = timeit()
            self.lastfit.migrad()
            toc =timeit()-tic
            self.log('migrad converged in {} s, {} calls, {} grads'.format(toc,self.lastfit.nfcn,self.lastfit.ngrad))
            # check if some parameters are positive parity 
        else:
            if self.grad:
                grad = self.the_model._add_multirun_grad_(*values_in)
                from numpy import set_printoptions as npopt, array
                npopt(threshold=1000)
                print('debug grad components as per minuit internal parameter index:') 
                print(grad)
        if self.dofit and self.lastfit.valid: 
            if pospar:
                self.log('... now redo without limits')
                #print('debug mufit dofit_multirun_singlegroup_globpardicts pospar = {}'.format(pospar))
                for k in pospar:
                    #print('debug mufit dofit_multirun_singlegroup_globpardicts k = {}, par = {}'.format(k,names[k]))
                    self.lastfit.limits[k] = [none,none]                    
                tic = timeit()                
                self.lastfit.migrad()
                toc =timeit()-tic
                self.log('migrad no limits redone in {} s, {} calls, {} grads'.format(toc,self.lastfit.nfcn,self.lastfit.ngrad))
        if self.dofit and self.lastfit.valid: 
            tic = timeit()                
            self.lastfit.hesse()
            tuc =timeit()-tic
            self.log('hesse in {} s, {} calls, {} grads'.format(tuc,self.lastfit.nfcn,self.lastfit.ngrad))
            
            n_runs = self.suite.nruns

            # write summary on console and log
            self.summary_multirun_global(start,stop,time[1]-time[0])

            # record result in csv file
            version = self.dashboard["version"]+'_'+version_flag(self)
            strgrp = stringify_groups(self.suite.groups)
            modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])

            # this is a one-shot csv, not incremental
            file_csv = self.suite.__csvpath__+modelname+'.'+strgrp+'.'+version+'.csv'
            self.write_multirun_glob_csv(file_csv,scan=self.scan)
            self.save_fit_multirun()
            self.lastfits.append(self.lastfit) #  muplotfit compatibility with multifits
        if self.dofit and not self.lastfit.valid:
            self.log('**** minuit did not converge! ****')
            print(self.lastfit)


           
    def dofit_multirun_multigroup_globpardicts(self,returntup):
        """
        returntup is a tuple of integers, start stop pack
        performs global fit of many-run, many-group data
        (c2) not yet
        """
        from mujpy.tools.tools import int2min_multirun, int2_multigroup_method_key

        a, e = self.suite.asymmetry_multirun_multigroup()
#       a and e are [[[run 0 grp 0 ... bins],
#                     [run 0 grp 1 ... bins], 
#                            ...            ],
#                    [[run 1 grp 0 ... bins],i
#                     [run 1 grp 1 ... bins]
#                            ...            ], 
#                      ...                   ] 

        self.log('not yet!')
        return

    def dofit_calib_multirun_multigroup_globpardicts(self,returntup):
        """
        returntup is a tuple of integers, start stop pack
        perform calib global fit of many-run, many-group data
        (c2 calib) in progress
        """
#        from numpy import array
        from mujpy.tools.tools import int2min_multirun_multigroup, int2_global_method_key

        krun, kgroup = 0,0 # krun used to generate file names, kgroup not used, overwritten in execute_...
        start, stop, pack = returntup
        # does not exist yet:
        yf,yb,eyf,eyb = self.suite.multirun_multigroup_for_back_counts(self.suite._the_runs_,self.suite.grouping) 
#       yf, etc are [[[run 0 grp 0 ... bins],
#                     [run 0 grp 1 ... bins], 
#                            ...            ],
#                    [[run 1 grp 0 ... bins],i
#                     [run 1 grp 1 ... bins]
#                            ...            ], 
#                      ...                   ] 
        #dt = self.suite.time[1]-self.suite.time[0]
        # from calib_singlerun_multigroup_globpardicts
        # adapt int2min and int2_ _method_key to mutirun_multigroup_globpardicts
        values,errors,fixed,limits,names, pospar = int2min_multirun_multigroup(self.dashboard,self.suite.runs)
#        print('debug mufit multirun multigroup: names {}, values {} errors {}'.format(names[0],values[0],errors[0]))
#        for k in range(1,len(values)):
#            print('                                       {},       {},       {}'.format(names[k],values[k],errors[k]))
#        print('debug mufit: len values errors = {}  {}'.format(len(values),len(errors)))
        self.methods_keys = int2_global_method_key(self.dashboard,self.the_model,self.suite.runs) 
#        p = [0.13,0.14,0.3,0.2,0.3,34.1,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7,10,0.1,11,0.7]

        load = (self.the_model._load_calib_multirun_multigroup_ if self.C2() else  
               self.the_model._load_calib_multigroup_ if self.A21() or self.B21() else 
               self.the_model._load_calib_) # implement later
        ok,errmsg = self.the_model._load_calib_multirun_multigroup_(self.suite,returntup,krun,kgroup,methods_keys)  #self.suite.time,yf,yb,eyf,eyb,returntup,self.methods_keys) 
        if not ok:
            self.log('Error in _load_multigroup_: '+errmsg)
            self.log('mufit stops here')            
            return
# buid an execute_log_save_fit for this case
        cost = self.the_model._chisquare_calib_
        summary = self.summary_global
        savefit = self.save_fit_multigroup
        self.execute_log_save_fit(cost,
                                  names,
                                  values,
                                  errors,
                                  limits,
                                  fixed,
                                  pospar,
                                  summary,
                                  start,
                                  stop,
                                  krun,
                                  kgroup,
                                  savefit)

 
    def summary_sequential(self, start, stop, kgroup, krun=0):
        """
        input: k is index in _the_runs_, default 0
        initial version: prints single fit single group result
        used by B1 multirun sequential singlegroup fits
        """
        from mujpy.tools.tools import get_title, chi2std, len_print_components, print_components, min2int, version_flag
        from mujpy.tools.tools import print_csv_components, write_csv
        from datetime import datetime

        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        version = self.dashboard["version"]+'_'+version_flag(self)
        the_run = self.suite._the_runs_[krun][0]
        nrun = the_run.get_runNumber_int()
        title = get_title(the_run)
        group = self.suite.groups[kgroup] # assumes only one group
        fgroup, bgroup, alpha = group['forward'],\
    					        group['backward'],\
    					        group['alpha']
        strgrp = fgroup.replace(',','_')+'-'+bgroup.replace(',','_')
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")  
        dt = (self.suite.time[1] - self.suite.time[0])*self.pack
        start, stop = self.suite.time[start]*1000, self.suite.time[stop]
        if krun==0:
            self.log('|'+77*'-'+'|') 
            fit_string = '| Fit [{:.2f}ns, {:.2}µs, {:.2f}ns/bin] on group: {} - {}  α = {:.3f}'
            self.log(fit_string.format(start,stop,dt,fgroup,bgroup,alpha)+8*' '+'|')
            fit_string = '|'+50*' '+'@{} |'
            self.log(fit_string.format(dt_string))
            self.log('|'+77*'-'+'|') 
        chi = self.lastfit.fval/self.number_dof#/self.lastfit.ndof #/self.number_dof 
        lowchi, highchi = chi2std(self.number_dof)
        file_log = self.suite.__cachepath__+modelname+'.'+str(nrun)+'.'+strgrp+'.'+version+'.log'
        names, values, errors = min2int(self.dashboard["model_guess"],
							        self.lastfit.values,self.lastfit.errors)
        with open(file_log,'w') as f:
            f.write(' '+85*'_'+'\n')
            f.write('| Run {}: {}              on group: {} - {}       α = {:.3f}'.format(nrun,
		                                 title,fgroup,bgroup,alpha)+4*' '+'|\n')

            self.log('| Run {}: {}         χᵣ² = {:.3f}({:.3f},{:.3f})'.format(nrun,
		                             title,chi,lowchi,highchi))
            f.write('| χᵣ² = {:.3f}({:.3f},{:.3f}), fit on [{:.2f}ns, {:.2}µs, {:.2f}ns/bin]   @{} \n'.format(chi,
		                                 lowchi,highchi,start,stop,dt*1000,dt_string))
            if not self.lastfit.valid:
                self.log('')
                self.log(27*'*'+' Minuit did not converge! '+27*'*')
                self.log('')
                f.write('')
                f.write(27*'*'+' Minuit did not converge! '+27*'*')
                f.write('')
            f.write('|'+85*'-'+'|\n') 
            self.log('|'+77*'-'+'|') 
            maxlen = 0
            par_err_str = ''
            for name,value,error in zip(names,values,errors): 
                maxlen = max(maxlen, len_print_components(name, value, error))   
            for name,value,error in zip(names,values,errors): 
                f.write('| '+print_components(name, value, error,maxlen)+'\n')
                par_err_str += print_csv_components(value,error)
            zip_forw, zip_backw = zip(names,values,errors), zip(reversed(names),reversed(values),reversed(errors))
            zip_nam_val_err = zip_forw if self.suite.console_method=='print' else zip_backw
            for name,value,error in zip_nam_val_err:
                self.log('| '+print_components(name, value, error,maxlen))
            f.write('|'+85*'_'+'|\n')
            self.log('|'+77*'-'+'|') 
            if self.verbose and not self.lastfit.valid:
                self.log(self.lastfit)

       # record result in csv file
#        version = self.dashboard["version"]+'_'+version_flag(self)
#        group = self.suite.groups[kgroup] # assumes only one group
#        fgroup, bgroup, alpha = group['forward'],\
#                                group['backward'],\
#                                group['alpha']
#        strgrp = fgroup.replace(',','_')+'-'+bgroup.replace(',','_')
#        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        file_csv = self.suite.__csvpath__+modelname+'.'+strgrp+'.'+version+'.csv'
#        the_run = self.suite._the_runs_[k][0]
        filespec = self.suite.datafile[-3:]
        header, row = self.prepare_csv_row(par_err_str,krun=k) 
        string1, string2 = write_csv(header,row,the_run,file_csv,filespec,scan=self.scan) 
        return string2

    def summary_multirun_global(self,start,stop):
        """
        deprecated
        print summary on Output and log file
        multirun glob version
        """
        from mujpy.tools.tools import get_title, chi2std, stringify_groups, value_error, version_flag
        from mujpy.tools.tools import len_print_components_multirun, print_components_multirun, min2int_multirun
        from datetime import datetime

        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        version = self.dashboard["version"]+'_'+version_flag(self)
        nrun0 = self.suite._the_runs_[0][0].get_runNumber_int()
        nrun1 = self.suite._the_runs_[-1][0].get_runNumber_int()
        title = get_title(self.suite._the_runs_[0][0])
        strgrp = stringify_groups(self.suite.groups)
        chi = self.lastfit.fval /self.number_dof 
        lowchi, highchi = chi2std(self.number_dof)
        dt = (self.suite.time[1] - self.suite.time[0])*self.pack
        start, stop = self.suite.time[start]*1000, self.suite.time[stop]
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")  

        nruns = str(nrun0)+'-'+str(nrun1)
        file_log = self.suite.__cachepath__+modelname+'.'+nruns+'.'+strgrp+'.g.'+version+'.log'
        # n_runs = self.suite._the_runs_
        names, values, errors = min2int_multirun(self.dashboard,
							        self.lastfit.values,self.lastfit.errors,self.suite._the_runs_)
        #print('debug mufit summary_multirun_global: names = {}\nvalues= {},errors = {}'.format(names,values,errors))
        fg1,bg1,al1 = self.suite.groups[0]['forward'], self.suite.groups[0]['backward'], self.suite.groups[0]['alpha'] 
        sumlength = 123
        with open(file_log,'w') as f:
            nch = sumlength - 2
            f.write(' '+nch*'_'+'\n')
            self.log(' '+nch*'_')
            string = '| Runs {}-{}: {}  Global fit {} on group: {} - {}   α = {:.3f}   '.format(nrun0,nrun1,title,dt_string,fg1,bg1,al1)
            nch = sumlength - 2 - len(string) if sumlength-len(string) - 2 >=0 else sumlength - 2
            f.write(string+nch*' '+' |\n')
#            print('debug summary_multirun len(string) {}'.format(len(string)))
            self.log(string+nch*' '+' |')
            string = '| χᵣ² = {:.3f}({:.3f},{:.3f}) ,    on [{:.2f}ns, {:.2}µs, {:.2f}ns/bin]'.format(chi,lowchi,highchi,start,stop,dt*1000)
            nch = sumlength - 2 - len(string) if sumlength-len(string) - 2 >=0 else sumlength - 2

            if self.verbose and not self.lastfit.valid:
                self.log(self.lastfiti)
    # check! and place here as well the Minuit did not converge! string

            f.write(string+nch*' '+'  |\n')
            self.log(string+nch*' '+' |')
            nparperrow = 10
            maxlen = 0     
            scan = self.suite.scan()
            for k,(nam,val,err) in enumerate(zip(names,values,errors)):   
                for na,va,er in zip([nam[i:i+nparperrow] for i in range(0, len(nam), nparperrow)],
                [val[i:i+nparperrow] for i in range(0, len(val), nparperrow)],
                [err[i:i+nparperrow] for i in range(0, len(err), nparperrow)]):
                    maxlen = max(maxlen,len_print_components_multirun(na, va, er))
                    if k==0: na0,va0,er0 = na,va,er
            namstring, _ = print_components_multirun(na,va,er,maxlen)
            
            nam0string, val0string = print_components_multirun(na0,va0,er0,maxlen)
            prestring = 'Run     '
            nrunstr = len(prestring)
            prestring += scan+'   ' # len(scan) = 4 + len blanks = 3 is 7 
            nbk = sumlength-len(namstring)-3-len(prestring)
            nbk0 = sumlength-len(nam0string)-3
            for k,(nam,val,err) in enumerate(zip(names,values,errors)):   # k=0 globals, k=1...nruns+1 run parameters, including locals
                for na,va,er in zip([nam[i:i+nparperrow] for i in range(0, len(nam), nparperrow)],
                [val[i:i+nparperrow] for i in range(0, len(val), nparperrow)],
                [err[i:i+nparperrow] for i in range(0, len(err), nparperrow)]): # na va er include k=0 globals (not used) 
                    if k==0:
                        # these are the global glob parameters
                        f.write('| '+nam0string+nbk0*' '+'|\n')
                        self.log('| '+nam0string+nbk0*' '+'|')
                        f.write('| '+val0string+nbk0*' '+'|\n')
                        self.log('| '+val0string+nbk0*' '+'|')
                        nch = sumlength - 2
                        f.write('|'+nch*'.'+'|\n')
                        self.log('|'+nch*'.'+'|')
                        f.write('| '+prestring+namstring+nbk*' '+'|\n')
                        self.log('| '+prestring+namstring+nbk0*' '+'|')
                    else:
                        # these are the run parameters and k=1 is run[0]
                        runscan = str(self.suite._the_runs_[k-1][0].get_runNumber_int())
                        runscan += (nrunstr-len(runscan))*' '
                        if scan[0]=='B':
                            field = self.suite._the_runs_[k-1][0].get_field()
                            fieldstring = '{:.0f}'.format(float(field[:field.index('G')])/10)
                            runscan += fieldstring + (7-len(fieldstring))*' '
                        elif scan[0]=='T':
                            TsTc, eTsTc = self.suite._the_runs_[k-1][0].get_temperatures_vector(), self.suite._the_runs_[k-1][0].get_devTemperatures_vector()
                            Tstring = value_error(TsTc[1],eTsTc[1])
                            runscan += Tstring + (7-len(Tstring))*' '
                        elif scan[0]=='[':
                            orientstring = self.suite._the_runs_[k-1][0].get_orient() 
                            runscan += orientstring + (7-len(orientstring))*' '
                        else:
                            runscan += 7*' '
                        _, valstring = print_components_multirun(na,va,er,maxlen)
                        nbk = sumlength-len(valstring)-4-len(runscan)
                        f.write('| '+runscan+valstring+nbk*' '+'|\n')
                        self.log('| '+runscan+valstring+nbk*' '+'|')                    
            f.write('|'+nch*'_'+'|\n')
            nch = sumlength - 2
            self.log('|'+nch*'_'+'|')

    def summary_multirun_multigroup_global(self,start,stop):
        """
        deprecated
        input: 
            start, stop: initial, final bin
        output: C2 C2_calib fit results on: log console print, cache/ .log file and fit/ .csv saves
        is called by self.execute_log_csv_fit, no return tuple
 
        """
        from mujpy.tools.tools import get_title, chi2std, stringify_groups, version_flag
        from mujpy.tools.tools import len_print_components, print_components, min2int_multirun_multigroup
        from mujpy.tools.tools import print_csv_components, write_csv
        from datetime import datetime

        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        version = self.dashboard["version"]+'_'+version_flag(self)
        the_run = self.suite._the_runs_[krun][0]
        nrun = the_run.get_runNumber_int()
        title = get_title(the_run)
        strgrp = stringify_groups(self.suite.groups)
        chi = self.lastfit.fval /self.number_dof 
        lowchi, highchi = chi2std(self.number_dof)
        dt = (self.suite.time[1] - self.suite.time[0])*self.pack
        start, stop = self.suite.time[start]*1000, self.suite.time[stop] # now times, in [ns] and [us]
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")  

        file_log = self.suite.__cachepath__+modelname+'.'+str(nrun)+'.'+strgrp+'.'+version+'.log'
        names, values, errors = min2int_multirun_multigroup(self.dashboard,
							        self.lastfit.values,self.lastfit.errors)
        file_csv = self.suite.__csvpath__+modelname+'.'+strgrp+'.'+version+'.csv'
        filespec = self.suite.datafile[-3:]
        string2 = []
 
# list (groups) of lists (omponents) of lists (parameters)
        sumlength = 100
        with open(file_log,'w') as f:
            f.write(' '+96*'_'+'\n')
            nch = sumlength - 2
            self.log(' '+nch*'_')
            string = '| Run {}: {}    Global fit of {}'.format(nrun,title,dt_string)
            f.write(string+21*' '+'|\n')
            self.log(string+24*' '+'|')
            string = '| χᵣ² = {:.3f}({:.3f},{:.3f}) ,    on [{:.2f}ns, {:.2}µs, {:.2f}ns/bin]'.format(chi,lowchi,highchi,start,stop,dt*1000)
            f.write(string+33*' '+'|\n')
            nch = sumlength - 1 - len(string) if sumlength-len(string) - 1 >=0 else sumlength - 1
            self.log(string+nch*' '+'|')
            if not self.lastfit.valid:
                self.log('')
                self.log(27*'*'+' Minuit did not converge! '+27*'*')
                self.log('')
                f.write('')
                f.write(27*'*'+' Minuit did not converge! '+27*'*')
                f.write('')
            for g1,n1,v1,e1,g2,n2,v2,e2 in zip(self.suite.groups[::2],names[::2],values[::2],errors[::2],
                                               self.suite.groups[1::2],names[1::2],values[1::2],errors[1::2]):
                fg1,bg1,al1 = g1['forward'], g1['backward'], g1['alpha'] 
                fg2,bg2,al2 = g2['forward'], g2['backward'], g2['alpha'] 

                string = ' on group: {} - {}   α = {:.3f}   |'.format(fg1,bg1,al1)
                nch = sumlength - 1 - len(string) if sumlength-len(string) - 1 >=0 else sumlength - 1
                f.write('|'+(nch-3)*'-'+string+'\n')
                self.log('|'+nch*'-'+string)
                maxlen = 0
                par_err_str = ''
                for nam,val,err in zip(n1,v1,e1):
                    maxlen = max(maxlen, len_print_components(nam, val, err))
                for nam,val,err in zip(n1,v1,e1):
                    f.write('| '+print_components (nam, val, err,maxlen)+'\n')
                    par_err_str += print_csv_components(val,err)
                zip_forw, zip_backw = zip(n1,v1,e1), zip(reversed(n1),reversed(v1),reversed(e1))
                zip_nam_val_err = zip_forw if self.suite.console_method=='print' else zip_backw
                for nam,val,err in zip_nam_val_err:
                    self.log('| '+print_components(nam, val, err,maxlen))
                kgroup = 0
                header, row = self.prepare_csv_row(par_err_str,kgroup=kgroup) 
                string1, string = write_csv(header,row,the_run,file_csv,filespec,scan=self.scan) 
                string2.append(string)
                
                string = ' on group: {} - {}   α = {:.3f}   |'.format(fg2,bg2,al2)
                nch = sumlength - 1 - len(string) if sumlength-len(string) - 1 >=0 else sumlength - 1
                f.write('|'+(nch-3)*'-'+string+'\n') 
                self.log('|'+nch*'-'+string)
                if self.verbose and not self.lastfit.valid:
                    self.log(self.lastfiti)
                maxlen = 0
                par_err_str = ''
                for nam,val,err in zip(n2,v2,e2):
                    maxlen = max(maxlen, len_print_components(nam, val, err))
                for nam,val,err in zip(n2,v2,e2):
                    f.write('| '+print_components (nam, val, err,maxlen)+'\n')
                    par_err_str += print_csv_components(val,err)
                zip_forw, zip_backw = zip(n2,v2,e2), zip(reversed(n2),reversed(v2),reversed(e2))
                zip_nam_val_err = zip_forw if self.suite.console_method=='print' else zip_backw
                for nam,val,err in zip_nam_val_err:
                    self.log('| '+print_components(nam, val, err,maxlen))
                kgroup = 1
                header, row = self.prepare_csv_row(par_err_str,kgroup=kgroup) 
                string1, string = write_csv(header,row,the_run,file_csv,filespec,scan=self.scan) 
                string2.append(string)
            nch = sumlength - 5
            f.write('|'+nch*'_'+'|\n')
        nch = sumlength - 2
        self.log('|'+nch*'_'+'|')
        return string2

    def save_fit_multirun(self):
        """
        is this in use?
        fit is multirun global (C1, C1_calib)
            saves a dashboard json adding the bestfit parameters as "globpardicts_result"
            and "model_result"
        to be consistent a single-run model_result is saved
        with lists of values, one per run, in place of single values as in the model_guess 
        filename is __cachepath__ + modelname + nruns + srtgrp + version.json
        nruns = shorthand for runNumbers, strgrp = shorthand for allgroups
        """
        from mujpy.tools.tools import stringify_groups, min2int_multirun, version_flag
        import json
        import os
        from copy import deepcopy
        
        # file name composition        
        # print('save_fit_multirun mufit debug: dashboard version {}'.format(self.dashboard['version']))
        version = self.dashboard["version"]+'_'+version_flag(self)
        strgrp = stringify_groups(self.suite.groups)
        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        the_runs = self.suite._the_runs_[:][0]
        nruns = str(the_runs[0].get_runNumber_int())+'-'+str(the_runs[-1].get_runNumber_int())
        file_json = self.suite.__fitpath__+modelname+'.'+nruns+'.'+strgrp+'.'+version+'_fit.json'
        model_result = deepcopy(self.dashboard["model_guess"])
        names, values, errors = min2int_multirun(self.dashboard,
							        self.lastfit.values,self.lastfit.errors,self.suite._the_runs_)
        # names, values, errors are list of lists, the first list is for the global parameters
        # the others lists are one for each run in the suite, and refer to the local parameters
        n_locals = 0
        n_globals = 0
        digits = '0123456789'
        for k, pardict in enumerate(self.dashboard['globpardicts_guess']):
            if pardict['local'] or type(pardict['value'])==list:
                n_locals += 1 # number of local glob parameters
        self.n_locals = n_locals
        globpardicts = []
        # model indices and names for local component parameters 
        componentindex = [k for k,component in enumerate(model_result) for pardict in component['pardicts'] if pardict['flag']=='~']
        parname =[pardict["name"] for component in model_result for pardict in component['pardicts']  if pardict['flag']=='~']

        for nam,val,err in zip(names[0],values[0],errors[0]): # global parameters
            globpardicts.append({'name':nam,'value':val,'std':err, 'local':False})
            n_globals += 1
        for j,nam in enumerate(names[1]): # names of minuit parameters for first run
            # print('debug mufit save_fit_multirun j = {}, nam = {}, n_locals = {}'.format(j, nam, n_locals))
            # the first n_locals appended to globpardicts
            if j<n_locals: # first ones are glob locals
                na = nam.rstrip(digits).rstrip('_') # stripped of run number
                va = [vals[j] for vals in values[1:]] # vals is a run list and val[j] is a glob local  
                er = [errs[j] for errs in errors[1:]] # errs is a run list and err[j] is its error
                # va and er ar lists over runs  
                globpardicts.append({'name':na,'value':va,'std':er,'label':'','local':True})
            for component in model_result:
                for pardict in component["pardicts"]:
                    if pardict["flag"] !="=":
                        pardict["name"] = nam.rstrip(digits).rstrip('_') # stripped of run number
                        pardict["value"] = [vals[j] for vals in values[1:]] # vals is a run list and val[j] is a component par  
                        pardict["std"] = [errs[j] for errs in errors[1:]] # errs is a run list and err[j] is its error
                #self.log('debug mufit save_fit_multirun: minuit name = {}, parname = {}'.format(na,model_result[index]["pardicts"]["name"]))
        self.dashboard["globpardicts_result"] = globpardicts
        self.dashboard["model_result"] = model_result
        self.dashboard["chi2"] = self.lastfit.fval /self.number_dof
        if os.path.isfile(file_json): 
            os.rename(file_json,file_json+'~')
        with open(file_json,"w") as f:
            json.dump(self.dashboard,f, indent=2,ensure_ascii=False)
        string_in = 'Best fit saved in {} '.format(file_json)
        self.log(string_in)

    """
       old flow:
            choosefit identifies cases A1, A1_calib, ... (see self.A1, ...)
            dofit_fittype executes each type
            
                dofit_singlerun_singlegroup
                    self.suite.asymmetry_single      A1
                    rebin
                    int2min
                    int2_method_key
                    the_model_._load_
                    self.summary
                    self.save_fit
                dofit_calib_singlerun_singlegroup    A1_calib
                    self.suite.single_for_back_counts
                    int2min
                    int2_method_key
                    the_model__load_calib_
                    self.summary
                    self.save_fit
                dofit_singlerun_multigroup_sequential   A20
                    self.suite.asymmetry_multigroup
                    rebin
                    int2min
                    int2_method_key
                    the_model_._load_
                    self.summary
                    self.save_fit
                dofit_calib_singlerun_multigroup_sequential A20_calib
                    self.suite.single_for_back_counts
                    int2min
                    int2_method_key
                    the_model._load_calib_
                    self.summary
                    self.save_fit
                dofit_singlerun_multigroup_globpardicts    A21
                    self.suite.asymmetry_multigroup
                    rebin
                    int2min_multigroup
                    int2_multigroup_method_key
                    the_model_._load_multigroup_
                    self.summary_global
                    self.save_fit_multigroup
                dofit_calib_singlerun_multigroup_globpardicts    A21_calib
                    self.suite.single_multigroup_for_back_counts
                    int2min_multigroup
                    int2_multigroup_method_key
                    the_model._load_calib_multigroup_
                    self.summary_global
                    self.save_fit_multigroup
                dofit_multirun_singlegroup_sequential    B1
                    self.suite.asymmetry_multirun
                    rebin
                    int2min
                    int2_method_key
                    the_model_._load_
                    self.summary
                    min2int
                    self.save_fit
        modes to do 
                dofit_multiruns_sequential_multigroup_globpardicts     
                # B20, B21, C1, C2, B20_calib, B21_calib, C1_calib
                    self.suite.multirun_multigroup_for_back_counts DONE
                    int2min_multirun_multigroup TODO

                dofit_calib_multirun_multigroup_globpardicts # C2_calib
                    self.suite.multirun_multigroup_for_back_counts DONE
                    int2min_multirun_multigroup TODO
                    int2_global_method_key TODO
                    the_model._load_calib_multirun_multigroup TODO
                    self.summary_global MUST UPGRADE
                    self.save_fit_?
    """

    def show_calib(self,plot_range):
        """
        Deprecated

        input:
            plot_range = '0,2000,40'
        output:
            t time 
            a asymmetry
            e asymmetry error
            f guess fit function for calib mode
        for debugging single run calibs  remove?
        """

        from mujpy.tools.tools import int2_method_key, int2min
        run = self.suite._the_runs_[0]
        yf, yb, eyf, eyb = self.suite.single_for_back_counts(run,self.suite.grouping[0])
        t = self.suite.time
        returntup,_ = derange(plot_range,self.suite.histoLength)
        par,_,_,_,name = int2min(self.dashboard,self.suite.runs)
        self.the_model._load_calib_(self.suite,returntup,0,0,
                                                  int2_method_key(self.dashboard,self.the_model))
        f = self.the_model._add_calib_(t,*par)
        return t,self.the_model._y_,e,f

