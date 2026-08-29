class mufit(object):
    """
    fit class, creates a mumodel from a json file, executes Minuit, produces logs and json files

    reads from a dashboard file
        will be interfaced to a dash gui
    """

    def __init__(self,suite,dashboard_file,chain=True,dash_log=None,no_fit=False, grad=False, scan = None, verbose = False):
        """
        accepts suite instance and dashboard file [options]

        input
            suite is the instance of the runs
            dashboard_file is a JSON file of a dictionary structure
        """

        from mujpy._version import __version__
        from mujpy.mucomponents.mucomponents import mumodel
        from mujpy.tools.tools import derange 

        self.__version__ = __version__
        self.suite = suite
        self.log = dash_log if dash_log else self.suite.console # input log prevails
        self.nodash = True
        self.nodata = True
        self.dofit = not no_fit # initialize_only True just loads guess values, no Minuit
        self.grad = grad 
        self.scan = scan  # is handled automatically if caller sets scan = the_suite.scan()[0]
        self.scan = 'θ' if self.scan == '[' else self.scan
        self.verbose = verbose
        self.chain = chain

        if self._dash_load_(dashboard_file):# model syntax and data/range ok
            if not self.suite.loadfirst: # is suite loaded?
                self.log('******* no data in musuite')
                self.log('* check access to database')
            else:
                self.nodata = False
                #print('mufit __init__ self.nodata {}'.format(self.nodata))
                returntup,errmsg = derange(self.dashboard["fit_range"],self.suite.histoLength) 
                if returntup[1]>0:
                    start, stop, self.pack = returntup # self.pack used in summary
                else:
                    self.log('fit range: {}, histoLength: {}, errmsg {},{}, check syntax!'.format(
                                      self.dashboard["fit_range"],self.suite.histoLength,returntup[0],returntup[1]))
                    self.nodata = True # misuse to stop execution here below           self.the_model = mumodel()
            self.lastfits = [] # lastfit initialization: 
            self.dofit_(returntup) # execute fit
        if self.nodata:# 
            self.log('     mufit stops here')
       
    def _dash_load_(self,dashboard_file):
        """
        try load dashboard file
        
        input:
            json dashboard_file produces a dict structure

        output:
            True/False if dashboard checks
                json syntax
                tools check_function
                multigroup global has function_multi
            are passed/failed 
        """

        from mujpy.tools.tools import check_function, short_path #
        import json


        try: # is dashboard_file readable?  
            self.log('Reading fit guess from: '+short_path(dashboard_file,self.suite.__startuppath__))
            with open(dashboard_file,"r") as f:
                self.dashboard = json.load(f)
                self.nodash = False                
        except Exception as e:
            print('Log file {}'.format(dashboard_file))  
            #self.log(getattr(e, 'message', repr(e)))
            self.log(getattr(e, 'message', str(e)))
             
            self.log('******* log file not found or corrupted')
            self.log('* {}'.format(e))

        kc, kp, kgroup = check_function(self.dashboard,self.suite.groups) # minimal model syntax check
        if kc >= 0 and kp > 0: # failed json dashboard syntax check
            self.log('++++ Check file {} model_guess: wrong syntax in component {}, parameter {}, group {}'.format(dashboard_file, kc, kp, kgroup))
            return False
        elif kc < 0 and kp > 0:
            self.log('++++ Check file {} globpardicts_guess: no function/function_multi allowed in parameter {}'.format(dashboard_file, kp))
            return False
        elif self.global_fit() and len(self.suite.groups)>1:
            # check if multigroup global have at least one 'function_multi' parameter key
            # at present int2_global_method_key and mufitplot relie on this feature for multigroup recognition
            if not any(['function_multi' in pardict.keys() for component in self.dashboard['model_guess'] for pardict in component['pardicts']]):
                self.log('++++ Check file {} globpardicts_guess or grouping: no function_multi while more groups') 
                return False
        return True

    def dofit_(self,returntup):
        """
        main method, auto invoked by _init_

        returntup is a tuple of integers, start stop pack
        switchyard for all fits
        to be tested
        """

        from mujpy.tools.tools import int2min, int2min_global, int2_method_key, int2_global_method_key 
        from mujpy.tools.tools import version_flag, list_depth, get_title#, chunk
        from mujpy.mucomponents.mucomponents import mumodel
        from iminuit import Minuit
        from numpy import printoptions
        from matplotlib.pyplot import subplots, show
        #from matplotlib.pyplot import subplots, show

        self.the_model = mumodel()
        start, stop, pack = returntup

        #######################
        # wrappers definition #
        #######################
        A1, A20, A21, B1, B20, B21, C1, C2 = self.fit_types()
        #mr, mg = (True, True) if C2 or B20 or B21 else (True, False) if C1 or B1 else (False, True) if A20 or A21 else (False, False)
        # whether suite is multirun (mr), multigroup (mg) 
        #sr, sg = (True, True) if B20 else (True, False) if B1 or B21 else (False, True) if A20 else (False, False)
        # whether fits are run sequential (sr), group sequential (sg)

        # set the wrappers, but use them inside the if switchyard
        int2min_ = int2min_global if C1 or C2 or A21 or B21 else int2min # if A1 or B1 or A20, or B20
        int2_method_key_ = int2_global_method_key if A21 or B21 or C1 or C2 else int2_method_key # if A1, A20, B1, B20
        cost = self.the_model._chisquare_calib_ if self.calib() else self.the_model._chisquare_
        switch = 'C2' if C2 else 'A21 or B21' if A21 or B21 else 'C1' if C1 else ''
        summary = self.summary if A1 or A20 or B1 or B20 else self.summary_global # if A21 or B21 or C1 or C2  
        savefit = self.save_fit if A1 or A20 or B1 or B20 else self.save_fit_global # if A21 else        
        runs = self.suite.runs[0] if B1 or B20 or B21 else self.suite.runs # the first value for B fits
        nruns = 1 if A21 or B21 else self.suite.nruns
        self.last_summary = False
 
        # use wrappers, replicate inside if switchyard
        methods_keys = int2_method_key_(self.dashboard,self.the_model,nruns)
        #n,npar = list_depth(methods_keys[1][1]),len(methods_keys[1][1])
        #print('mufit dofit_ list_depth(methods_keys[1][1]) {}, len(methods_keys[1][1]) {}'.format(n,npar))
        if not methods_keys: # if translate, called by int2_method_keys gets an exception
            self.log('     mufit stops here')
            return
        values,errors,fixed,limits,names,pospar = int2min_(self.dashboard,runs) 
        self.lastfit = Minuit(cost,
                              name=names,
                              *values)  # the only instance   
#       print(self.lastfit)  # debug to see input minuit parameters
        self.lastfit.errors = errors
        self.lastfit.fixed = fixed

        ############################
        # scan the sequentials     #
        # B1, B20, B21, A20        #
        # and execute the singles  #
        # A1, A21, C1, C2          #
        # works also for calib     #
        ############################
        if B20:
            for krun in range(self.suite.nruns):
                for kgroup in range(len(self.suite.grouping)):
                    self.the_model._load_(self.suite,methods_keys,krun,kgroup,self.fit_types) # loads slice, 
                    # the_model detects calib: first 'al' component, methods_keys[0][0] = None
                    self.the_model._rebin_(returntup,*values) # updates alphas if _calib and rebins 
                    if not krun and not kgroup: # only the first time
                        self.ndof_(switch)
                        
                    self.execute_fit(limits,pospar) # checks if self.dofit

                    if self.dofit: 
                        if self.calib():
                            self.suite.groups[kgroup]["alpha"] = self.lastfit.values[0]
                            self.suite.grouping[kgroup]["alpha"] = self.lastfit.values[0]
                        summary(start, stop, krun, kgroup) # self.summary
                        savefit(krun,kgroup) # self.save_fit
                        self.lastfit.values = values if not self.chain else self.lastfit.values
            nr, ng, nt = self.suite.nruns,len(self.suite.groups),self.suite.nruns*len(self.suite.groups)
            if self.dofit: self.log('                                       {} runs x {} groups = {} best fits'.format(nr,ng,nt))

        elif B1 or B21:
            for krun in range(self.suite.nruns):
                kgroup = 0 if B1 else -1 # if B21
                self.the_model._load_(self.suite,methods_keys,krun,kgroup,self.fit_types) # loads slice, 
                self.the_model._rebin_(returntup,*values) # updates alphas if _calib and rebins
                if not krun: # only the first time
                    self.ndof_(switch)
                self.execute_fit(limits,pospar)

                if self.dofit:
                    lastgroup = len(self.suite.grouping) - 1
                    for kgroup in range(len(self.suite.grouping)): # loop for summary only, but not for B1
                        if self.calib():
                            self.suite.groups[kgroup]["alpha"] = self.lastfit.values[kgroup]
                            self.suite.grouping[kgroup]["alpha"] = self.lastfit.values[kgroup]
                        self.last_summary = True if kgroup == lastgroup else False
                        #summary(start, stop, krun, kgroup) # self.summary_global for B21, self.summary for B1
                        summary(start, stop, 0, kgroup) # self.summary_global for B21, self.summary for B1
                    savefit(krun,0) # save_fit for B1, save_fit_global for B21
                    self.lastfit.values = values if not self.chain else self.lastfit.values
            nr, ng = self.suite.nruns,len(self.suite.groups)
            if self.dofit: self.log('                                       {} runs, {} group[s], {} best fits'.format(nr,ng,nr))

        elif A20:
            krun = 0
            for kgroup in range(len(self.suite.grouping)):
                self.the_model._load_(self.suite,methods_keys,krun,kgroup,self.fit_types) # loads slice, 
                self.the_model._rebin_(returntup,*values) # updates alphas if _calib and rebins
                if not kgroup: # only the first time
                    self.ndof_(switch)
                self.execute_fit(limits,pospar)
                if self.dofit: 

                    if self.calib():
                        self.suite.groups[kgroup]["alpha"] = self.lastfit.values[0]
                        self.suite.grouping[kgroup]["alpha"] = self.lastfit.values[0]
                    summary(start, stop, krun, kgroup) # self.summary
                    savefit(krun,kgroup) # self.save_fit
                    self.lastfit.values = values if not self.chain else list(self.lastfit.values)

        else: # A1, A21, C1, C2 
            krun,kgroup = (-1,-1) if A1 else (0,-1) if A21 else (-1,0) if C1 else (-1,-1)
            self.the_model._load_(self.suite,methods_keys,krun,kgroup,self.fit_types) # loads slice, 
            self.the_model._rebin_(returntup,*values) # updates alphas if _calib and rebins
            self.ndof_(switch)
            self.execute_fit(limits,pospar)
            if self.dofit: 

                last = (len(self.suite.runs)-1,len(self.suite.grouping)-1)
                for krun in range(len(self.suite.runs)): # loops for summary only, but not for A1
                    for kgroup in range(len(self.suite.grouping)): # no loop for A1, C1
                        if self.calib():
                            self.suite.groups[kgroup]["alpha"] = self.lastfit.values[kgroup]
                            self.suite.grouping[kgroup]["alpha"] = self.lastfit.values[kgroup]
                        self.last_summary = True if (krun,kgroup) == last else False
                        summary(start, stop, krun, kgroup) # self.summary for A1, self.summary_global for A21, C1, C2
                savefit(0,0) # self.save_fit for A1, self.save_fit_global for A21, C1, C2
            nr, ng = self.suite.nruns,len(self.suite.groups)
            if C1 or C2 and self.dofit:
                self.log('                                       {} runs , {} group[s], 1 best fit'.format(nr,ng))


        if not self.dofit:
            self.log('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            self.log('+     Minuit initialized with guess values and 1 fcn call       +')
            self.log('+           mufitplot will show the fit function[s]             +')
            self.log('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            self.log('')

    def execute_fit(self,limits,pospar):
        """
        Minuit execution standard, auto invoked by dofit_

        input: 
            limits, pospar, see int2min
        executes Minuit (twice for pospars) 
        appends to lastfits
        both plain and calib 
        """

        # print('mr is multirun {}. mg is multigroup {}'.format(mr == multirun, mg == multigroup))
        # if true mr,mg can be omitted from execute ..
        self.lastfit.limits = limits
        values = list(self.lastfit.values) # a deep copy of guess values
        if self.dofit: 
            #print('mufit execute_fit self.lastfit.values = {}'.format(self.lastfit.values))
            self.lastfit.migrad()
            if pospar: # positive defined parameters
                for k in pospar:
                    self.lastfit.limits[k] = [None,None]   
                self.lastfit.migrad()
            self.lastfit.hesse()
        else:
            self.lastfit.migrad(ncall=1) # loads and changes Minuit values in 1 iteration
            self.lastfit.values = values # restores guess values for plot
        # Minuit parameters for mufitplot (self.lastfit may change)
        values = list(self.lastfit.values) # saves a result deep copy
        self.lastfits.append(values) # stores result deep copy 

    def ndof_(self,switch):
        """
        calculates model number of degrees of freedom

        input:
            switch = 'C2' or 'A21 or B21' or 'C1'
            requires a completed self.lastfit
        saves self.number_dof
        """

        # required in summary and save_fit
        self.number_dof =  self.the_model._x_.size 
        #print('mufit ndof_ _x_.size {}'.format(self.the_model._x_.size))
        self.number_dof *= len(self.suite.groups) * self.suite.nruns if switch=='C2' else (
                           len(self.suite.groups)) if switch=='A21 or B21' else ( 
                           self.suite.nruns) if switch=='C1'     else 1
        self.number_dof -= self.lastfit.nfit

    def summary(self,start, stop, krun, kgroup):
        """
        single run, single group A1 A20 B1 B20 & calib results logged and saved (fit, csv)

        input: 
            krun index in runs
            kgroup index in grouping
        output:
                  saved on cache/file.log 
                  printed on self.log
                  saved in csv/file.csv 
                  saved in fit/file.json
        invoked by self.dofit_ inside a krun and/or kgroup loops
        """

        from mujpy.tools.tools import get_title, chi2std, stringify_groups, version_flag
        from mujpy.tools.tools import print_components, min2int
        from mujpy.tools.tools import print_csv_components, write_csv
        from datetime import datetime
        import os

        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        version = self.dashboard["version"]+'_'+version_flag(self)
        the_run = self.suite._the_runs_[krun][0]
        nrun = the_run.get_runNumber_int() # = int(the_run)
        title = get_title(the_run)
        strgrp = stringify_groups(self.suite.groups) 
        chi = self.lastfit.fval /self.number_dof
        # print('mufit summary self.lastfit.fval, ndof = {},{}'.format(self.lastfit.fval,self.number_dof))
        lowchi, highchi = chi2std(self.number_dof)
        dt = (self.suite.time[1] - self.suite.time[0])*self.pack
        start, stop = self.suite.time[start]*1000, self.suite.time[stop]
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")  

        # logs  model+grp+run+grp+version for each individual fit (A1, A20, B1, B20)
        #     ignore if bloated, goes in cache!
        file_log = self.suite.__cachepath__+modelname+'.'+str(nrun)+'.'+strgrp+'.'+version+'.log'
        nruns = 1 # dummy
        names, values, errors, shared = min2int(self.dashboard,
							                    self.lastfit.values,self.lastfit.errors,
                                                krun,kgroup,self.the_model._slice_nruns_,self.the_model._slice_ngroups_)
        # names values, errors, shared are internal values referred to mudel_guess out of a single-run-sungle-group fit
        # model csv: adds one line for each individual fit (A1, A20, B1, B20
        strg = '' if self.A20() or self.A20_calib() or self.B20() or self.B20_calib() else strgrp+'.'
        # run is in csv and for 20 fits grp is too 
        file_csv = self.suite.__csvpath__+modelname+'.'+strg+version+'.csv'
        group = self.suite.groups[kgroup] # assumes only one group [at a time]
        fgroup, bgroup, alpha = (group['forward'], 
    					        group['backward'],
    					        group['alpha'])

        wa = 'w' if (krun,kgroup)==(0,0) else 'a' # appends if Krun or kgroup > 0  
        if wa=='w' and os.path.isfile(file_log): # create prev version 
            os.rename(file_log,file_log+'~')
        # A1 A20 B1 B20 always 'w' 
        with open(file_log,wa) as f:
            f.write(' '+85*'_'+'\n')
            self.log('|'+77*'-'+'|') 
            f.write('| Run {}: {}       group: {} - {}      offs = {}  α = {:.3f}'.format(nrun,
		                                 title,fgroup,bgroup,self.suite.offset,alpha)+3*' '+'|\n')

            f.write('| χᵣ² = {:.3f}({:.3f},{:.3f}), fit on [{:.2f}ns, {:.2}µs, {:.2f}ns/bin]   @{} \n'.format(chi,
		                                 lowchi,highchi,start,stop,dt*1000,dt_string))
            self.log('| Run {}: {}     α = {:.3f}  offs = {}  group: {} - {} |'.format(nrun,
                                                                                        title,alpha,self.suite.offset,fgroup,bgroup))
            self.log('|χᵣ² = {:.3f}({:.3f},{:.3f}), on [{:.2f}ns, {:.2}µs] {:.2f}ns/bin @{}|'.format(chi,
		                                 lowchi,highchi,start,stop,dt*1000,dt_string))
#            if not self.lastfit.valid:
#                self.log('')
#                self.log(27*'*'+' Minuit did not converge! '+27*'*')
#                self.log('')
#                f.write('')
#                f.write(27*'*'+' Minuit did not converge! '+27*'*')
#                f.write('')
            f.write('|'+85*'-'+'|\n') 
            self.log('|'+77*'-'+'|')
            maxlen = 0
            par_err_str = ''
            for name,value,error,share in zip(names,values,errors,shared):
                f.write('| '+print_components(name, value, error, share)+'\n')
                par_err_str += print_csv_components(value,error)
            zip_forw, zip_backw = zip(names,values,errors,shared), zip(reversed(names),reversed(values),reversed(errors),reversed(shared))
            zip_nam_val_err = zip_forw if self.suite.console_method=='print' else zip_backw
            for name,value,error,share in zip_nam_val_err:
                self.log('| '+print_components(name, value, error, share))
            f.write('|'+85*'-'+'|')
            self.log('|'+77*'_'+'|')
            self.if_not_converged(f)

        filespec = self.suite.datafile[-3:]
        header, row = self.prepare_csv_row(par_err_str,krun=krun,kgroup=kgroup) 
        string1, string2 = write_csv(header,row,the_run,file_csv,filespec,scan=self.scan)
        self.log(string2)
        return

    def summary_global(self,start,stop,krun,kgroup):
        """
        global A21 B21 C1 C2 fits & calib results logged and saved (fit, csv)

        input:
            start, stop, fit bins
            krun index in run
            kgroup index in grouping
        called
            A21 ngroup times, C1 nrun times, B21, C2 nrun*ngroup times 
            shows single fit, single group results      in self.log
            adds them                to a unique cache/file_log.log
            adds a row                 to a unique csv/file_csv.csv
        if (krun,kgroup)==(0,0) 
           C1,C2     print a unique Minuit cache/U_file_log.log
           A21, B21 adds a Minuit row to a unique csv/U_file_csv.csv           
        invoked by self.dofit_ 
        """

        from mujpy.tools.tools import get_title, chi2std, stringify_groups, version_flag
        from mujpy.tools.tools import print_components
        from mujpy.tools.tools import min2int_global
        from mujpy.tools.tools import value_error, chunk
        from mujpy.tools.tools import print_csv_components, write_csv#, write_glob_csv
        from datetime import datetime
        import os

        #   print('mufit summary_global krun, kgroup {} {}'.format(krun, kgroup))
        _,_,A21,_,_,B21,C1,C2 = self.fit_types()
        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        version = self.dashboard["version"]+'_'+version_flag(self)
        the_run = self.suite._the_runs_[krun][0] # this is an instance 
        title = get_title(the_run)
        nrun = the_run.get_runNumber_int()
        #print('mufit summary_global nrun {}'.format(nrun))a # this nrun is correc
        strgrp = stringify_groups(self.suite.groups) if C1 else stringify_groups(self.suite.groups,joinch='+') 
        chi = self.lastfit.fval/self.number_dof 
        lowchi, highchi = chi2std(self.number_dof)
        dt = (self.suite.time[1] - self.suite.time[0])*self.pack
        start, stop = self.suite.time[start]*1000, self.suite.time[stop]
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")  

        strun = str(nrun)+'.' if A21 or B21 else '' # C1, C2  have run in the name
        #   bloated, but goes to csv
        # strg = strgrp+'.' if C1 else ''
        file_log = self.suite.__cachepath__+modelname+'.'+strun+strgrp+version+'.log'
        # csv write group and run in each row
        file_csv = self.suite.__csvpath__+modelname+'.'+version+'.csv'
        #print('mufit summary_global file_csv {}'.format(file_csv))

        scan = '' if not self.scan else 'θ' if self.scan[0] == '[' else self.scan[0]
        cnames, cvalues, cstds, cnonhash = min2int_global(self.dashboard,self.lastfit.values,self.lastfit.errors,
                                                      krun,kgroup,self.the_model._slice_nruns_,self.the_model._slice_ngroups_) # comp list of par lists
        # print('mufit summary_global cnames {}'.format(cnames))
        # names values, errors, honhash are model_guess internal ones obtained by a global fit
        filespec = self.suite.datafile[-3:]
        string2 = []
 
# list (groups) of lists (omponents) of lists (parameters)
        sumlength = 100
        wa = 'w' if (krun,kgroup)==(0,0) else 'a' # appends if Krun or kgroup > 0  
        if wa=='w' and os.path.isfile(file_log): # create prev version 
            os.rename(file_log,file_log+'~')
        # A21 B21 C1 C2 always 'w' 
        with open(file_log,wa) as f:
            if wa=='w' or (B21 and kgroup==0):  
            # write Run chi2 header
                f.write(' '+96*'_'+'\n')
                nch = sumlength - 2
                self.log(' '+nch*'_')
                if A21 or B21: 
                    string = '| Run {}: {}    Global fit of {}'.format(nrun,title,dt_string)
                    f.write(string+21*' '+'|\n')
                    self.log(string+24*' '+'|')
                    string = '| χᵣ² = {:.3f}({:.3f},{:.3f}) ,    on [{:.2f}ns, {:.2}µs, {:.2f}ns/bin]'.format(chi,lowchi,highchi,start,stop,dt*1000)
                    f.write(string+33*' '+'|\n')
                else:
                    string = '| Global fit χᵣ² = {:.3f} ({:.3f},{:.3f}) ,    on [{:.2f}ns, {:.2}µs, {:.2f}ns/bin]'.format(chi,lowchi,highchi,start,stop,dt*1000)
                    f.write(string+21*' '+'|\n')
                nch = sumlength - 2 - len(string) if sumlength-len(string) - 1 >=0 else sumlength - 1
                self.log(string+nch*' '+'|')
                string= '|'
                if all([cn for component in cnonhash for cn in component]): # A21 and B21
                    n, v, e = [pardict['name'] for pardict in self.dashboard['globpardicts_guess']], self.lastfit.values, self.lastfit.errors
                    for nam, val, err in zip(n,v,e):
                        string += ' '+nam+'='+value_error(val,err)
                    f.write(string+'\n')
                    self.log(string)
                    # add row[s] to a global csv/'U'+file_csv
                    #if wa=='w' and (A21 or B21):
                    # write unique Minuit
                    csv_format = '{:.4f}'
                    csv_large = '{:.4e}'
                    large = 8
                    minstr = [csv_format.format(par) if len(csv_format.format(par))<= large  else (
                              csv_large.format(par)) for par in self.lastfit.values] # minuit values in {:.4f}/{:.4e} format
                    eminstr = [csv_format.format(std) if len(csv_format.format(par))<= large else (
                               csv_large.format(std)) for par,std in zip(self.lastfit.values,self.lastfit.errors)] # minuit stds in {:.4f}/{:.4e} format
                    par_err_str = ','.join([x for xs in zip(minstr,eminstr) for x in xs])+',' 
                    n = [pardict['name'] for pardict in self.dashboard['globpardicts_guess']] # just plain globpardicts names
                    # reconstruct run replicas

                    file = self.suite.__csvpath__+'U_'+modelname+'.'+version+'.csv'
                    header, row = self.prepare_globcsv_row(n,par_err_str,krun)
                    string1, string0 = write_csv(header,row,the_run,file,filespec,scan=scan) 
                    self.log(string0)
                else: # C1 and C2
                    # a sinle log file with une line run
                    start = len(self.dashboard['globpardicts_guess'])
                    nhash = len([k for k,pardir in enumerate(self.dashboard['globpardicts_guess']) if pardir['flag']=='#'])
                    globnames = [self.lastfit.params[i].name for i in range(self.lastfit.npar)]
                    n, v, e  = globnames[:start], self.lastfit.values[:start], self.lastfit.errors[:start]
                    strun = str(nrun)
                    for nam, val, err in zip(n,v,e):
                        if strun in nam:
                            f.write(string+'\n')
                            self.log(string)
                            string= '|'+20*' '
                            strun = 'strun'
                        string += ' '+nam+'='+value_error(val,err)
                    f.write(string+'\n')
                    self.log(string)
                    n, v, e  = globnames[start:], self.lastfit.values[start:], self.lastfit.errors[start:]
                    for nam, val, err in zip(chunk(n,nhash),chunk(v,nhash),chunk(e,nhash)):
                        string = '|'+20*' '
                        for na,va,er in zip(nam,val,err):
                            string += ' '+na+'='+value_error(va,er)
                        f.write(string+'\n')
                        self.log(string)
###################################################################
# this is a single krun, single kgroup global summary             #
#           each:                printed on self.log              #
#                                saved in csv/file.csv rows       #
#           krun,kgroup=(0,0)#   saved in a single cache/file.log #
#                                saved in a single  fit/file.json #
###################################################################
            runs = self.suite.runs[krun] 
            # names is list of lists of lists of lists of parameters
            # for krn,runadd in enumerate(runs): # A21, B21 no loop  
#                for g1,n1,v1,e1 in zip(self.suite.groups[::2],names[::2],values[::2],errors[::2]),
#                                               self.suite.groups[1::2],names[1::2],values[1::2],errors[1::2]):
            g = self.suite.groups[kgroup]
            fg,bg,al = g['forward'], g['backward'], g['alpha'] 
            string = ' on group: {} - {}   α = {:.3f}   |'.format(fg,bg,al)
            nch = sumlength - 1 - len(string) if sumlength-len(string) - 1 >=0 else sumlength - 1
            if not (self.C2_calib()): 
                f.write('|'+(nch-3)*'-'+string+'\n')
                self.log('|'+nch*'-'+string)
            par_err_str = ''
            for kc,(n,v,e,nh) in enumerate(zip(cnames,cvalues,cstds,cnonhash)): # names, ... model_guess [component[parameter]] for krun, kgroup
                #print('mufit summary_global len(n) len(nh) {} {}'.format(len(n),len(nh))) 
                if not (self.C2_calib() and kc==0):
                    f.write('| '+print_components(n, v, e, nh,global_fit=True)+'\n')
                    self.log('| '+print_components(n, v, e, nh, global_fit=True))
                else:
                    f.write('| '+print_components(n, v, e, nh,global_fit=True)+3*' '+(nch-26)*'-'+string+'\n')
                    self.log('| '+print_components(n, v, e, nh, global_fit=True)+3*' '+(nch-22)*'-'+string)
                par_err_str += print_csv_components(v,e)
            nch = sumlength - 2
            header, row = self.prepare_csv_row(par_err_str,krun=krun,kgroup=kgroup) 
            string1, string_out = write_csv(header,row,the_run,file_csv,filespec,scan=scan)
            run='+'.join(self.suite.runs[krun])
            get_value = the_run.get_temp() if scan == 'T' else the_run.get_field() if scan == 'B' else the_run.get_orient() if scan == 'θ' else ''
            string = '|    Run '+run+': '+scan+'='+get_value
            f.write(string+'\n')
            self.log(string+'  '+string_out)
            if self.last_summary:
                f.write('|'+(nch-4)*'_'+'|\n')
                self.log('|'+nch*'_'+'|')
            #string2.append(string)

    #       close f

        if self.last_summary: # only last summary_global call  
            self.if_not_converged(f)

        # now write Minuit parameter csv
  
 #           else: # C1 or C2 # write a global log/'U_'+file_log
 #               file = self.suite.__cachepath__+'U_'+modelname+'.'+strun+strgrp+version+'.log'
 #               if os.path.isfile(file_log): # create prev version 
 #                   os.rename(file_log,file_log+'~')
 #               with open(file,'w') as ff:
 #                   ff.write(' '+96*'_'+'\n')
 #                   nch = sumlength - 2
 #                   string = '| Run {}: {}    Global fit of {}'.format(nrun,title,dt_string)
 #                   ff.write(string+21*' '+'|\n')
 #                   string = '| χᵣ² = {:.3f}({:.3f},{:.3f}) ,    on [{:.2f}ns, {:.2}µs, {:.2f}ns/bin]'.format(chi,lowchi,highchi,start,stop,dt*1000)
 #                   ff.write(string+33*' '+'|\n')
 #                   nch = sumlength - 1 - len(string) if sumlength-len(string) - 1 >=0 else sumlength - 1
 #                   string = 'Run: '+run+'  '+scan + ' = ' + get_value
 #                   string = ' on group: {} - {}   α = {:.3f}   |'.format(fg,bg,al)
 #                   nch = sumlength - 1 - len(string) if sumlength-len(string) - 1 >=0 else sumlength - 1
 #                   ff.write('|'+(nch-3)*'-'+string+'\n'+'|'+nch*'_'+'|\n')
 #                   for nam,val,err in zip(n,v,e): # one line per Minuit parameter
 #                       ff.write('| {}:  {} '.format(n,value_error(v,e)))
 #                       if not self.lastfit.valid:
 #                           ff.write('')
 #                           ff.write(27*'*'+' Minuit did not converge! '+27*'*')
 #                           ff.write('')
 #                   ff.write('|'+nch*'_'+'|\n')
        return 

    def if_not_converged(self,f):
        """
        prints varoious things if Minuit not converged 
        """

        if not self.lastfit.valid:
            self.log('')
            self.log(27*'*'+' Minuit did not converge! '+27*'*')
            self.log('')
            f.write('')
            f.write(27*'*'+' Minuit did not converge! '+27*'*')
            f.write('')
        # how many groups?
        if self.verbose and not self.lastfit.valid:
            strout = str(self.lastfit).split('\n')
            for k in range(len(strout)):
                self.log(strout[k])
    
    def save_fit(self,krun,kgroup):
        """
        saves indiviudual A1, A20, B1, B20 results in new json .fit files, with fit type label

        input:
            krun is index in self.suite._the_runs_
            kgroup is indek in self.suite.groups
        saves a dashboard file adding the bestfit parameters as "model_result"
        filename is __cachepath__ + modelname + nrun  + strgrp + version .json
        nrun = runNumber, strgrp = shorthand for group
        These saves are  for individual runs: A1, A20, B1
            use "version" as additional label to qualify fit
        """

        from mujpy.tools.tools import min2int, version_flag
        import json
        import os
        from copy import deepcopy
                
        version = self.dashboard["version"]+'_'+version_flag(self)
        group = self.suite.groups[kgroup] # assumes only one group
        fgroup, bgroup, alpha = group['forward'],\
				                group['backward'],\
				                group['alpha']
        strgrp = fgroup.replace(',','_')+'-'+bgroup.replace(',','_')
        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        the_run = self.suite._the_runs_[krun][0]
        nrun = str(the_run.get_runNumber_int())
        file_json = self.suite.__fitpath__+modelname+'.'+nrun+'.'+strgrp+'.'+version+'.json'

        # replace (guess) values with Minuit.values, leave error as step, add fit_range, std and chi2
        # do not replace names, they are autogenerated by mufit
        names, values, errors,_ = min2int(self.dashboard,self.lastfit.values,self.lastfit.errors,
                                          krun,kgroup,self.the_model._slice_nruns_,self.the_model._slice_ngroups_)
        # print(self.dashboard["model_guess"])
        self.dashboard["model_result"] = deepcopy(self.dashboard["model_guess"])
        for k,component in enumerate(self.dashboard["model_result"]):
            value, std = values[k], errors[k]
            for j,pardict in enumerate(component['pardicts']):
                pardict["value"] = value[j]
                pardict["std"] = std[j]
        self.dashboard["chi2"] = self.lastfit.fval /self.number_dof
        self.dashboard["grp_calib"] = self.suite.groups
        if os.path.isfile(file_json): 
            os.rename(file_json,file_json+'~')
        with open(file_json,"w") as f:
            json.dump(self.dashboard,f, indent=2,ensure_ascii=False) # ,object_pairs_hook=OrderedDict)
        short_json = file_json.replace(self.suite.__startuppath__,'.')              
        self.log('{}  saved.  '.format(short_json))

    def save_fit_global(self,krun,kgroup):
        """
        saves A21, B21, C1, C2 results in json .fit files, with fit type label

        input:
            krun index in self.suite._the_runs_
            kgroup is dummy for compatibility with save_fit
        fit is global
            saves one dashboard json adding the bestfit parameters as "globpardicts_result"
        Use "version" as additional label to qualify fit (auto 'g_
        filename is __cachepath__ + modelname + nrun  + srtgrp0 + strgrp...  + version .json
        nrun = runNumber, strgrp0,1,... = shorthand for allgroups
        """

        from mujpy.tools.tools import stringify_groups, version_flag, min2int_global
        import json
        import os
        from copy import deepcopy
        
        # file name composition        
        # print('save_fit_multigroup mufit debug: dashboard version {}'.format(self.dashboard['version']))
        version = self.dashboard["version"]+'_'+version_flag(self)
        strgrp = stringify_groups(self.suite.groups)
        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        the_run = self.suite._the_runs_[krun][0]
        nrun = str(the_run.get_runNumber_int())
        # print('debug {}.{}.{}.{}'.format(modelname,nrun,strgrp,version))
        file_json = self.suite.__fitpath__+modelname+'.'+nrun+'.'+strgrp+'.'+version+'.json'
        self.dashboard["model_result"] = deepcopy(self.dashboard["model_guess"])
        self.dashboard["globpardicts_result"] = deepcopy(self.dashboard["globpardicts_guess"])
        _, _, _, nonhash = min2int_global(self.dashboard,self.lastfit.values,self.lastfit.errors,
                                          krun,kgroup,self.the_model._slice_nruns_,self.the_model._slice_ngroups_) # comp list of par lists
        chi2 = []
#        names = [parameter["name"] for parameter in self.dashboard["globpardicts_guess"]]
        if all([nh for nhc in nonhash for nh in nhc]): # no hash parameters A21 or B21, all minuit parameters are globpars
            #print('mufit save_fit_global nonhash {}'.format(nonhash))
            j = -1
            for value,std in zip(self.lastfit.values,self.lastfit.errors):
                j += 1
                self.dashboard["globpardicts_result"][j]['value'] = value
                self.dashboard["globpardicts_result"][j]['std'] = std
        else:  # C1 or C2, save globpars and define rest as hash_results
            for j,pardict in enumerate(self.dashboard['globpardicts_result']):
                pardict['value'] = self.lastfit.values[j]
                pardict['std'] = self.lastfit.errors[j]
            start = len(self.dashboard['globpardicts_result'])
            self.dashboard['hash_result'] = []
            names = [self.lastfit.params[i].name for i in range(start,self.lastfit.npar)]
            for name,value,std in zip(names,self.lastfit.values[start:],self.lastfit.errors[start:]):
                self.dashboard['hash_result'].append({'name':name,'value':value,'std':std})
# must add here also the "model_guess" with values replaced by fit results.
# for multigroup_globpardicts these are the model 
        self.dashboard["chi2"] = self.lastfit.fval /self.number_dof
        if os.path.isfile(file_json): 
            os.rename(file_json,file_json+'~')
        with open(file_json,"w") as f:
            json.dump(self.dashboard,f, indent=2,ensure_ascii=False)
        string = '{} saved'.format(file_json)
#        print('mufit save_fit_multigroup debug, string_in:\n{}'.format(string_in))
# ,object_pairs_hook=OrderedDict)    
        self.log(string)

    def prepare_csv_row(self,par_err_str,krun = 0,kgroup=0):
        """
        writes the same header row, model single row for all fit types

        input: 
            par_err_str text string containing csv of parameter values, std 
            krun, kgroup are indices in self.suite._the_runs_[krun][0], self.suite.groups[kgroup]
        output: 
            header, the model specific csv header 
                    to compare with that of the csv file
            row, the line to be added to the csv file
        prepares a csv-like row of best fit parameters 
        that can be imported to produce figures
        """

        from mujpy.tools.tools import chi2std, init_csv_row, stringify_group #, get_title, spec_prec
        from mujpy.tools.tools import chi2_csv, min2int_names
        #from mujpy.tools.tools import min2_csv, min2_global_csv #, list_depth

        
        # A1, A20, A21, B1, B20, B21, C1, C2 = self.fit_types()
        filespec = self.suite.datafile[-3:]
        lowchi, hichi = chi2std(self.number_dof)

        # if A21, B21, C1, C2: # includes calib
        #    nruns = len(self.suite.runs) # is 1 for A21
        #    ngroups = len(self.suite.groups) # is 1 for C1
        #    vstr, estr, mvstr, mestr = min2_global_csvself.dashboard,self.lastfit.values,self.lastfit.errors,self.suite.nruns)
        #    # A21, B21: single run multi group write row and append in group loop 
        #    # C1: multi run single group write row in and append in group loop is fine, only one append per run
        #    # C2: multi both write row and append in group is fine
        #    # component parameters appended sequentially at the same level per each run and group 
        #    #minparam2_csv(self.dashboard,self.lastfit.values,self.lastfit.errors,multirun=self.suite._the_runs_)
        #    #nest = list_depth(vstr)i# always 3
        #    rows = []
        #    dof_single = self.number_dof/(nruns*ngroups)
        #    for krun in range(nruns): # one line per run, only model parameter csv
        #        for kgroup in range(ngroups): # one line per group 
        #            the_run = self.suite._the_runs_[krun][0]
        #            k, l = (krun, None) if C1 else (kgroup, None) if B21 or A21 else (krun,kgroup) # if C2  
        #            chi_run = self.the_model._chisquare_single_(*self.lastfit.values,k=k,l=l)/= dof_single
        #            group = self.suite.groups[kgroup] 
        #            fgroup, bgroup, alpha = group['forward'],\
        #                                    group['backward'],\
        #                                    group['alpha']
        #            row = init_csv_row(the_run.get_field(), filespec, the_run, group = stringify_group(group) )
        #            for v,e in zip(vstr[krun][kgroup][:],estr[krun][kgroup][:])
        #                row +=  v+','+e,+','
        #        row += chi2_csv(chi_run,lowchi,hichi,alpha,self.suite.offset)   
        #        row += '\n'
        #        rows.append(row)
        #    row = rows # in this case, C1, row is a list of rows
        # else: # A1 or A20 or B1 or B20, includes calib
            # A1, A20, B1, B20 are all single single when it comes to the summary
        the_run = self.suite._the_runs_[krun][0]
        group = self.suite.groups[kgroup] # assumes only one group
        fgroup, bgroup, alpha = group['forward'],\
                                group['backward'],\
                                group['alpha']
        row = init_csv_row(filespec, the_run, group=stringify_group(group)) # run T B valid for all fits
        # group is always included in the csv
        model =  self.dashboard["model_guess"]
        chi = self.lastfit.fval/self.number_dof # true fit reduced chi square, fval is cost (chi2) at the minimum
        row += par_err_str 
        row += chi2_csv(chi,lowchi,hichi,alpha,self.suite.offset)
        row += '\n'
        # row is formatted with appropriate rounding, write directly
        # self.console(row)
        if filespec == 'bin' or filespec == 'mdu' or filespec == 'root':
            header = ['#0.Run','1.T_cryo[K]','2.e_T_cryo[K]','3.T_sample[K]','4.e_T_sample[K]','5.B[G]','6.group']
            k = 6
        else: # ISIS NeXus
            header = ['#0.Run','1.T[K]','2.eT[K]','3.B[G]','4.group']
            k = 4
        # now model parameter names for header
        nruns = 1
        names = min2int_names(self.dashboard)
        # works for all fits to produce model component list of parameter names lists
        for parnames in names: 
            # print('mufit prepare_csv_row debug: parlist: {}'.format(parlist))
            for parname in parnames:
                # print('mufit prepare_csv_row debug: parnames: {}'.format(parname))
                k += 1
                header.append(str(k)+'.'+parname)
                k += 1
                header.append(str(k)+'.'+'e_'+parname)
        k += 1
        header.append(str(k)+'.'+'chi2_r')
        k +=1
        header.append(str(k)+'.'+'chi2_low')
        k += 1
        header.append(str(k)+'.'+'chi2_hi')
#        if self.suite.multi_groups():
#            for jgroup,group in enumerate(self.suite.groups):
#                header.append('alpha{}'.format(jgroup))
#        else:
        k += 1
        header.append(str(k)+'.'+'alpha')
        k += 1
        header.append(str(k)+'.'+'offset_bin')
        header.append('timestring\n')
        return ','.join(header), row

    def prepare_globcsv_row(self,names,par_err_str,krun):
        """
        writes the same header row, Minuit paramenter row for B21,A21 fits

        input: 
            par_err_str text string containing csv of Minuit parameter values, std 
            krun index in run list
        output: 
            header, a globper specific csv header 
                    with a string for all runs and a string for all groups
            row, the line to be added to the csv file
        prepares a csv-like row of best fit parameters 
        that can be imported to produce figures
        """

        from mujpy.tools.tools import chi2std, init_csv_row, stringify_group #, get_title, spec_prec
        from mujpy.tools.tools import chi2_csv, min2int_names
        #from mujpy.tools.tools import min2_csv, min2_global_csv #, list_depth

        
        # A1, A20, A21, B1, B20, B21, C1, C2 = self.fit_types()
        filespec = self.suite.datafile[-3:]
        lowchi, hichi = chi2std(self.number_dof)
        the_run = self.suite._the_runs_[krun][0]
        row = init_csv_row(filespec, the_run) # run T B valid for all fits, group = False default 
        # group is always included in the csv
        #model =  self.dashboard["model_guess"]
        if self.calib():
            alpha = 1.00
        else:
            alpha = [group['alpha'] for group in self.suite.groups]
        chi = self.lastfit.fval/self.number_dof # true fit reduced chi square, fval is cost (chi2) at the minimum
        print('debug prepare_globcsv_row\nrow = {}\npar_err_str = {}'.format(row,par_err_str))
        row += par_err_str 
        row += chi2_csv(chi,lowchi,hichi,alpha,self.suite.offset) # now accepts also lists of alphas
        row += '\n'
        # row is formatted with appropriate rounding, write directly
        # self.console(row)
        if filespec == 'bin' or filespec == 'mdu' or filespec == 'root':
            header = ['#0.Run','1.T_cryo[K]','2.e_T_cryo[K]','3.T_sample[K]','4.e_T_sample[K]','5.B[G]']
            k = 5
        else: # ISIS NeXus
            header = ['#0.Run','1.T[K]','2.eT[K]','3.B[G]']
            k = 3
        # now model parameter names for header
        # works for all fits to produce model component list of parameter names lists
        for name in names:
            # print('mufit prepare_csv_row debug: parnames: {}'.format(parname))
            k += 1
            header.append(str(k)+' '+name)
            k += 1
            header.append(str(k)+' '+'e_'+name)
        k += 1
        header.append(str(k)+' '+'chi2_r')
        k +=1
        header.append(str(k)+' '+'chi2_low')
        k += 1
        header.append(str(k)+' '+'chi2_hi')
#        if self.suite.multi_groups():
#            for jgroup,group in enumerate(self.suite.groups):
#                header.append('alpha{}'.format(jgroup))
#        else:
        if self.calib():
            k += 1
            header.append(str(k)+' '+'alpha')
        else:
            for a in alpha:
                k += 1
                header.append(str(k)+' '+'alpha')
        k += 1
        header.append(str(k)+' '+'offset_bin')
        header.append('timestring\n')
        return ','.join(header), row

    def calib(self):
        """
        True if first model component is 'al'
        """
        
        return self.dashboard['model_guess'][0]['name']=='al'

    def globpar(self): # this is a global fit using globpardicts
        """
        True if global fit A21, B21, C1, C2

        self.globpar() synonym of self.global_fit() 
        """
        return "globpardicts_guess" in self.dashboard.keys()

    def tilde_in_component(self): # this dashboard has minuit parameters in the model components
        """
        Deprecated

        see table at the bottom of https://musr-nmr.unipr.it/dispense/pmwiki.phpdd?n=Mujpy.GlobalSwitch
        """

        # empty list, no ~ flags, is equivalent to False, non empty list is True
        return any([par['flag']=='~' for component in self.dashboard["model_guess"]  for par in component['pardicts']]) 

    def hash_in_component(self): # this dashboard has minuit parameters in the model components
        """
        True for C1 and C2 fits

        global fits with 'flag':'#' in 'globpardicts_guess'
        these parameters generate a labelled Minuit replica for each run in the suite 
        """

        # empty list, no ~ flags, is equivalent to False, non empty list is True
        return any([pardict['flag']=='#' for pardict in self.dashboard["globpardicts_guess"]])

    def A1(self): # single run singlegroup 
        """
        True for A1 fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """
        return self.suite.single() and not (self.calib() or self.suite.multi_groups() or self.globpar())
            
    def A1_calib(self): # single run calib singlegroup 
        """
        True for A1 calib fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return self.suite.single() and self.calib() and not (self.suite.multi_groups() or self.globpar())
            
    def A20(self): # single run multigroup sequential 
        """
        True for A20 fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return self.suite.single() and self.suite.multi_groups() and not self.globpar() and not self.calib()

    def A20_calib(self): # single run calib multigroup sequential 
        """
        True for A20 calib fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return self.suite.single() and self.suite.multi_groups() and not self.globpar() and self.calib()
                
    def A21(self): # single run multigroup global 
        """
        True for A21 fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return self.suite.single() and self.suite.multi_groups() and self.globpar() and not self.calib()
                
    def A21_calib(self): # single run calib multigroup global 
        """
        True for A21 calib fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return self.suite.single() and self.suite.multi_groups() and self.globpar() and self.calib()
                
    def B1(self): # multirun sequential singlegroup 
        """
        True for B2 fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return not (self.suite.single() or self.suite.multi_groups() or self.globpar() or self.calib())

    def B1_calib(self): # multirun sequential singlegroup 
        """
        True for B1 calib fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return not (self.suite.single() or self.suite.multi_groups() or self.globpar()) and self.calib()

    def B20(self): # multirun sequential multigroup global
        """
        True for B20 fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return self.suite.multi_groups() and not (self.suite.single() or self.globpar() or self.calib())

    def B20_calib(self):
        """
        True for B20 calib fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return self.suite.multi_groups() and not (self.suite.single() or self.globpar()) and self.calib()

    def B21(self): # multirun sequential multigroup global
        """
         True for B21 fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return self.suite.multi_groups() and self.globpar() and not self.suite.single() and not self.hash_in_component() and not self.calib()

    def B21_calib(self):
        """
        True for B21 calib fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return self.suite.multi_groups() and self.globpar() and not self.suite.single() and not self.hash_in_component() and self.calib()

    def C1(self): # multirun global singlegroup
        """
        True for C1 fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """
        return not self.suite.multi_groups() and not self.suite.single() and self.globpar() and self.hash_in_component() and not self.calib()

    def C1_calib(self): # multirun global singlegroup
        """
         True for C1 calib fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return not self.suite.multi_groups() and not self.suite.single() and self.globpar() and self.hash_in_component() and self.calib()

    def C2(self): # multirun global multigroup global
        """
        True for C2 fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return self.suite.multi_groups() and not self.suite.single() and self.globpar() and self.hash_in_component() and not self.calib()

    def C2_calib(self): # multirun global multigroup global
        """
        True for C2 calib fit

        see https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy
        """

        return self.suite.multi_groups() and not self.suite.single() and self.globpar() and self.hash_in_component() and self.calib()

    def check_fit(self):
        """
        check that fit verifications are mutually exclusive
        """

        def single_true(iterable):
            i = iter(iterable)
            return any(i) and not any(i)

        return single_true([self.A1(), self.A1_calib(),
                            self.A20(), self.A20_calib(), 
                            self.A21(), self.A21_calib(), 
                            self.B20(), self.B20_calib(),
                            self.B21(), self.B21_calib(),
                            self.C1(), self.C1_calib(),
                            self.C2(), self.C2_calib()])
               
    def fit_types(self):
        """
        list all A1, ... , C2 booleans, in order

        Xyz is True for Xyz or Xyz_calib
        """

        A1 = self.A1() or self.A1_calib()
        A20 = self.A20() or self.A20_calib()
        A21 = self.A21() or self.A21_calib()
        B1 = self.B1() or self.B1_calib()
        B20 = self.B20() or self.B20_calib()
        B21 = self.B21() or self.B21_calib()
        C1 = self.C1() or self.C1_calib()
        C2 = self.C2() or self.C2_calib()
        return A1, A20, A21, B1, B20, B21, C1, C2

    def global_fit(self):
        """
        True for global fits

        True for global fits A21,B21, C1, C2 & calib 
        False for plain fits A1, A20, B1, B20 & calib
        """

        return "globpardicts_guess" in self.dashboard.keys()

#    def which_true(self):
#        """
#        tell which fit it may be
#        """
#        string = ['A1','A1_calib',
#                  'A20','A20_calib',
#                  'A21','A21_calib',
#                  'B1','B1_calib',
#                  'B20','B20_calib',
#                  'B21','B21_calib',
#                  'C1','C1_calib',
#                  'C2','C2_calib']
#        fits = [self.A1(), self.A1_calib(),
#                self.A20(), self.A20_calib(), 
#                self.A21(), self.A21_calib(), 
#                self.B20(), self.B20_calib(),
#                self.B21(), self.B21_calib(),
#                self.C1(), self.C1_calib(),
#                self.C2(), self.C2_calib()]
#        return [string[k] for k,fit in enumerate(fits) if fit]
 
