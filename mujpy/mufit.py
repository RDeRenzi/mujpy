class mufit(object):
    '''
    fit class 
    reads from a dashboard file
        can also be used to generate the gui
    '''
    def __init__(self,suite,dashboard_file,chain=False,dash=None,initialize_only=False, grad=False, scan = None):
        '''
        input
            suite is the instance of the runs
            dashboard_file is a JSON file of a dictionary structure
        '''
        from mujpy._version import __version__
        self.__version__ = __version__
        self.suite = suite
        self.dash = dash
        self.log = self.dash.log if self.dash else self.suite.console
        self.nodash = True
        self.nodata = True
        self.n_locals = [] # used only by C1 and C2
        self.dofit = not initialize_only # initialize_only True just loads guess values, no Minuit
        self.grad = grad 
        self.scan = scan  # set Scan = 'T' to order csv by temperature, 'B' to order csv by field  
        self._initialise_fit(dashboard_file,chain)   
        
    def _initialise_fit(self,dashboard_file,chain):
        '''
        input:
            json dashboard_file produces a dict structure
        '''
        # from collections import OrderedDict
        from mujpy.tools.tools import _available_components_
        from mujpy.mucomponents.mucomponents import mumodel
        import json

        if not self.suite.loadfirst:
            self.log('******* no data in musuite')
            self.log('* check access to database')
        else:
            self.nodata = False
        try:  
            with open(dashboard_file,"r") as f:
                self.dashboard = json.load(f) # ,object_pairs_hook=OrderedDict)
                # print('mufit _initialise_fit debug: dash {}'.format(self.dashboard['model_guess']))
                self.nodash = False                
                # dashboard is a dict structure, not an Ordered Dictionary structure              
        except Exception as e:
            print('Log file {}'.format(dashboard_file))  
            self.log(getattr(e, 'message', repr(e)))
            self.log(getattr(e, 'message', str(e)))
             
            self.log('******* log file not found or corrupted')
            self.log('* {}'.format(e))
            
# DELETE?            
#        self.available_components = _available_components_() 
#        # list of templates dictionaries 'name' 'pardicts' 
#        # now each pardict contains only 'name', 'error', 'limits'
#        # e.g. 'name':'A','error':0.01, 'limits':(0,0)
        
        # self.log("* Check dashboard method *")
        self.component_names = [item['name'] for item in _available_components_()]
        self._the_model_ = mumodel()
        self._the_model_single_ = mumodel()
        self.lastfits = [] # lastfit initialization: 
                           # A1, A21 and calib, C1, C2 will add a single Minuit instance
                           # A20 and calib, B1, B2 will add a sequence of instances 
                           # self.lastfit survives for backward compatibility
                           # and is always the last
                           # it is also simply appended to self.lastfits, a list
#        self.log("**** Fit initialized *****")
            
        if self.choosefit(chain):
            self.log('     mufit stops here')
#        else:
#            print('mufitplot debug: Really finished')

        #self.log('{}'.format(self.lastfit.params))       
    
    def choosefit(self,chain):
        '''
        select type of fit (Ai,Bi,Ci, i=1,2)
            i = 1,2 single or multi groups (single cost function)
            A, B single or multiple sequential runs (single or multiple cost functions)
            C global (single cost function)   
        '''
        from iminuit import Minuit
        from mujpy.tools.tools import derange, add_step_limits_to_model, checkvalidmodel
        from mujpy.tools.tools import model_name, userpars, userlocals, multigroup_in_components
        
        if self.nodata or self.nodash:
            return True        

        ok, msg = checkvalidmodel(model_name(self.dashboard),self.component_names)

        if not ok:
            self.log(msg)
            self.log('*** Check your dashboard:')
            return True
        else:
            # buffer to transfer sequential fits to save_fit_multigroups
            #   used only by A20, left empty otherwise, can be used to detect A20
            self.names = []
            self.values = []
            self.stds = [] 
            self.fvals = []
            # add errors and limits to dashboard            
            # self.dashboard = add_step_limits_to_model(self.dashboard)
            # require six switches: suite.single, suite.multi_groups, calib, userpardicts, sequential_runs 

            #####################
            # begins switchyard #
            #####################
            # in case remove previous global tags
            version = self.dashboard["version"]
            while len(version)>3 and version[0]=='g' and version[2]=='_':
                if version[1] in {'r','g','G'}:
                    version = version[3:]
                else:
                    break 
            # print('choosefit mufit debug:    self.suite.single() = {}'.format(self.suite.single()))  
               
            returntup,errmsg = derange(self.dashboard["fit_range"],self.suite.histoLength) 
                                        # histoLength set after asymmetry_single
            if returntup[1]>0:
                start, stop, pack = returntup
            else:
                self.log('fit range: {}, histoLength: {}, errmsg {},{}, check syntax!'.format(
                                  self.dashboard["fit_range"],self.suite.histoLength,returntup[0],returntup[1]))
                return  True # = stop here|
                            
            if self.A1():                       # A1 singlerun singlegroup DONE
                self.dofit_singlerun_singlegroup(returntup)
                # print('mufit choosefit debug: should be finished!')       
                
            elif self.A1_calib():               # A1 calib singlerun singlegoup DONE 
                self.dofit_calib_singlerun_singlegroup(0,returntup)  

            elif self.A20():                   # A2.1 singlerun multigroup sequential DONE
                self.dofit_singlerun_multigroup_sequential(returntup)   
            
            elif self.A20_calib():             # A2.0 calib singlerun multigroup sequential DONE
                self.dofit_calib_singlerun_multigroup_sequential(returntup)  # DONE (?)
                
            elif self.A21():                   # A2.1 singlerun multigroup global DONE
                self.dashboard["version"]=('gg_'+version if 
                                self.dashboard["version"][0:3]!='gg_' else version)
                self.dofit_singlerun_multigroup_userpardicts(returntup) 

            elif self.A21_calib():             # A2.1 calib singlerun multigroup global DONE
                self.dashboard["version"]=('gg_'+version if 
                                self.dashboard["version"][0:3]!='gg_' else version)        
                self.dofit_calib_singlerun_multigroup_userpardicts(returntup)

            elif self.B1():                     # B1 multirun singlegroup sequential DONE
#                self.log('Doing a multirun singlegroup sequential fit')
#                self.log('singlerun = {}, calib = {}, multi_groups = {}, userpar = {}, tilde = {}'.format(self.suite.single(),self.calib(),self.suite.multi_groups(),self.userpar(),self.tild_in_component()))
                self.dofit_multirun_singlegroup_sequential(returntup,chain)

            elif self.B20():                   # B2.0 multirun multigroup sequential Does not exist
                self.dashboard["version"]=('gg_'+version if 
                                self.dashboard["version"][0:3]!='gg_' else version)        
                # self.dofit_multirun_multigroup_sequential(returntup,chain)
                print("multirun multigroup sequential fit doesn't exist yet")

            elif self.B21():                   # B2.1 multirun_sequential multigroup_global DONE
                self.dashboard["version"]=('gg_'+version if 
                                self.dashboard["version"][0:3]!='gg_' else version)
                self.dofit_sequentialrun_multigroup_userpardicts(returntup,chain) 

            elif self.C1():                     # C1 multirun singlegroup global DOING
#                self.log('Doing a multirun singlegroup global fit')
                self.dofit_multirun_singlegroup_userpardicts(returntup) 

            elif self.C2():                     # C2 multirun singlegroup global NOT YET
                self.dashboard["version"]=('gg_'+version if 
                                self.dashboard["version"][0:3]!='gg_' else version)
                self.dofit_multirun_multigroup_userpardicts(returntup)      
            else:
                self.log('Not clear which fit!')
                return True
        return False

    def dofit_calib_singlerun_singlegroup(self,kgroup,returntup):
        '''
        performs calib fit on single run, single group, 
        (A1-calib) tested
        input 
            kgroup is group index in suitegrouping
        '''
        from iminuit import Minuit
        from mujpy.tools.tools import int2min, int2_method_key, rebin_decay, write_csv
        from mujpy.mucomponents.mucomponents import mumodel
        _the_model_ = mumodel()                
        # self.log('In single calib.')      

        yf,yb,bf,bb,yfm,ybm = self.suite.single_for_back_counts(self.suite._the_runs_[0],self.suite.grouping[kgroup]) 
                              # the second dimension is group
        start, stop, pack = returntup
        t,yf,yb,bf,bb,yfm,ybm = rebin_decay(self.suite.time,yf,yb,bf,bb,[start,stop],pack)

        fitvalues,fiterrors,fitfixed,fitlimits,names,pospar = int2min(self.dashboard["model_guess"])

        # print('dofit_calib_singlerun_singlegroup mufit debug: fitvalues = {}'.format(fitvalues))
#        for k in range(len(fitvalues)):
#            self.log('{} = {}, step = {}, fix = {}, limits ({},{})'.format(names[k], values[k],errors[k],fixed[k],limits[k][0],limits[k][1]))

        _the_model_._load_calib_single_data_(t,yf,yb,bf,bb,yfm,ybm,
                                                  int2_method_key(self.dashboard,_the_model_))
                                             # int2_int() returns a list of methods to calculate the components
        self.lastfit = Minuit(_the_model_._chisquare_,
                              name=names,
                              *fitvalues) 
        # print('dofit_calib_singlerun_singlegroup mufit debug: fitvalues = {}'.format(fitvalues))                                       
        self.lastfit.errors = fiterrors
        self.lastfit.limits = fitlimits
        self.lastfit.fixed = fitfixed
        # self.freepars = self.lastfit.nfit
        self.number_dof = len(t) - self.lastfit.nfit
        if self.dofit:
            self.lastfit.migrad()
            # check if some parameters are positive parity 
            if pospar:
                for k in pospar:
                    self.lastfit.limits[k] = [None,None]                    
                self.lastfit.migrad()
            self.lastfit.hesse()
        self.lastfits.append(self.lastfit) # not really needed here

        if self.dofit:
            kgroup = 0
            if self.lastfit.valid:
                self.suite.groups[0]["alpha"] = self.lastfit.values[0]
                self.suite.grouping[0]["alpha"] = self.lastfit.values[0]
                # write summary on console

                self.summary(start, stop, t[1]-t[0],kgroup) 
                # record result in csv file
                version = self.dashboard["version"]
                group = self.suite.groups[kgroup] # assumes only one group
                fgroup, bgroup, alpha = group['forward'],\
					                    group['backward'],\
					                    group['alpha']
                strgrp = fgroup.replace(',','_')+'-'+bgroup.replace(',','_')
                modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
                file_csv = self.suite.__csvpath__+modelname+'.'+version+'.'+strgrp+'.csv'
                the_run = self.suite._the_runs_[0][0]
                filespec = self.suite.datafile[-3:]
                
                header, row = self.prepare_csv() 
                
                string1, string2 = write_csv(header,row,the_run,file_csv,filespec,scan=self.scan) 

                # self.log(string1)
                #self.log(string2)
                krun = 0
                
                self.save_fit(krun,0,string2) 
               
            else:
                self.log('**** Minuit did not converge! ****')
                print(self.lastfit)
            return True

    def dofit_singlerun_singlegroup(self,returntup):  
        '''
        performs fit on single run, single group
        (A1) tested
        '''
        from iminuit import Minuit
        from mujpy.tools.tools import int2min, int2_method_key, rebin, write_csv
          
        #self.log('In single single: grouping is {}'.format(self.suite.grouping)) 

        a,e = self.suite.asymmetry_single(self.suite._the_runs_[0],0) # runs to be added, group index
        start, stop, pack = returntup
        time,asymm,asyme = rebin(self.suite.time,a,[start,stop],pack,e=e)

#        for cmp in self.dashboard["model_guess"]:
#            for pd in cmp["pardicts"]:
#                print('dofit_singlerun_singlegroup debug: {} = {}, step = {}, fix = {}, limits ({},{})'.format(pd['name'], pd['value'],pd['error'],pd['flag'],pd['limits'][0],pd['limits'][1]))

        values,errors,fixed,limits,names, pospar = int2min(self.dashboard["model_guess"])
#        self.log('in dofit_singlerun_singlegroup:\nnames = {}\nvalues = {}\nerrors = {}\n fixed = {}\nlimits = {}'.format(names,values,errors,fixed,limits))                                  


        self._the_model_._load_data_(time,asymm,int2_method_key(self.dashboard,self._the_model_),
                                     e=asyme) 
                                     # pass data to model, one at a time
        ############################## int2_int() returns a list of methods to calculate the components
        # actual single migrad calls
        # print('dofit_singlerun_singlegroup mufit debug: names {}, values {}'.format(names,values))
        self.lastfit = Minuit(self._the_model_._chisquare_,
                              name=names,
                              *values)     
        self.lastfit.errors = errors
        self.lastfit.limits = limits
        self.lastfit.fixed = fixed
        # self.freepars = self.lastfit.nfit
        self.number_dof = len(asymm) - self.lastfit.nfit
        # print('mufit dofit_singlerun_singlegroup debug: Minuit name value error limits fixed {}'.format([[name,value,error,limit,fix] for name,value,error,limit,fix in zip(names,values,errors,limits,fixed)]))
        if self.dofit:
            self.lastfit.migrad()
            # check if some parameters are positive parity 
            if pospar:
                for k in pospar:
                    self.lastfit.limits[k] = [None,None]   
                    self.log('....fit again wih no lims')                 
                self.lastfit.migrad()
            self.lastfit.hesse()
        self.lastfits.append(self.lastfit) # not really needed by this fit

        if self.dofit:
            # write summary on console
            kgroup = 0
            self.summary(start, stop, time[1]-time[0],kgroup)

            # record result in csv file
            version = self.dashboard["version"]
            group = self.suite.groups[kgroup] # assumes only one group
            fgroup, bgroup, alpha = group['forward'],\
					                group['backward'],\
					                group['alpha']
            strgrp = fgroup.replace(',','_')+'-'+bgroup.replace(',','_')
            modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
            file_csv = self.suite.__csvpath__+modelname+'.'+version+'.'+strgrp+'.csv'
            the_run = self.suite._the_runs_[0][0]
            filespec = self.suite.datafile[-3:]
            header, row = self.prepare_csv()
            string1, string2 = write_csv(header,row,the_run,file_csv,filespec,scan=self.scan)
            # self.log(string1)
            #self.log(string2)
            krun,kgroup = 0,0
            self.save_fit(krun,kgroup,string2)
                  
    def dofit_singlerun_multigroup_sequential(self,returntup):
        '''
        performs fit on single run, multi-group data sequentially
        (A20) tested
        '''
        from iminuit import Minuit
        from mujpy.tools.tools import int2min, int2_method_key, min2int
        from mujpy.tools.tools import rebin, write_csv
        
        a,e = self.suite.asymmetry_multigroup() # the second dimension is group
        start, stop, pack = returntup
        time,asymm,asyme = rebin(self.suite.time,a,[start,stop],pack,e=e)
        # print('dofit_singlerun_multigroup_sequential mufit debug: {} {} {}'.format(time.shape,asymm.shape,asyme.shape))

        values_in,errors,fixed,limits,names,pospar = int2min(self.dashboard["model_guess"])

        # print('dofit_singlerun_multigroup_sequential mufit debug: names {}, values_in {}'.format(names,values_in))#        for name,value,error,fix,limit in zip(names,values,errors,fixed,limits):
#            self.log('dofit_singlerun_multigroup_sequential debug{} = {}, step = {}, fix = {}, limits ({},{})'.format(name, value,error,fix,limit[0],limit[1]))

        krun = 0  #  single run!!
        string = []
        for kgroup,(a,e) in enumerate(zip(asymm,asyme)):
            # print('dofit_singlerun_multigroup_sequential mufit debug in loop: names {}, values_in {}'.format(names,values_in))
            ok, errmsg = self._the_model_._load_data_(
                                        time,a,
                                        int2_method_key(self.dashboard,self._the_model_),
                                        e=e) 
            if not ok:
                self.log(repr(errmsg))
                break
            # actual single migrad calls

            self.lastfit = Minuit(self._the_model_._chisquare_,
                                  name=names,
                                  *values_in)                                        
            self.lastfit.errors = errors
            self.lastfit.limits = limits
            self.lastfit.fixed = fixed
            self.number_dof = len(a) - self.lastfit.nfit
            if self.dofit:
                self.lastfit.migrad()
                # check if some parameters are positive parity 
                if pospar:
                    for k in pospar:
                        self.lastfit.limits[k] = [None,None]                    
                    self.lastfit.migrad()
                self.lastfit.hesse()
            self.lastfits.append(self.lastfit)

            if self.dofit:
               # write summary on console
                self.summary(start, stop, time[1]-time[0],kgroup)

                # record result in csv file
                version = self.dashboard["version"]
                group = self.suite.groups[kgroup] # assumes only one group
                fgroup, bgroup, alpha = group['forward'],\
					                    group['backward'],\
					                    group['alpha']
                strgrp = fgroup.replace(',','_')+'-'+bgroup.replace(',','_')
                modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
                file_csv = self.suite.__csvpath__+modelname+'.'+version+'.'+strgrp+'.csv'
                the_run = self.suite._the_runs_[0][0]
                filespec = self.suite.datafile[-3:]
                header, row = self.prepare_csv()
                string1, string2 = write_csv(header,row,the_run,file_csv,filespec,scan=self.scan)
                # self.log(string1)
                #self.log(string2)
                #self.save_fit(krun,kgroup,string2)
                string.append(string2)
                names, values, stds = min2int(self.dashboard["model_guess"],
						                self.lastfit.values,self.lastfit.errors)
                # in summary values are obtained from self.lastfit.values 
                # within the loop (one at a time) and produce right output 
                # in save_fit_multigroup to fill component "values" like in
                
                self.names.append(names)
                self.values.append(values)
                self.stds.append(stds)
                self.fvals.append(self.lastfit.fval)
                # print('dofit_singlerun_multigroup_sequential mufit debug: self.values {}'.format(self.values))
                self.save_fit_multigroup(krun,string)

    def dofit_singlerun_multigroup_userpardicts(self,returntup):
        '''
        performs fit on single run, global multi-group data
        (A21) tested
        All minuit parameters must be predefined as user pardicts
        All component parameters must be assigned by functions to the previous
        (the absence of "flag":"~" parameters identifies this type of fit)
        It could be also run sequentially (next devel)
        '''
        from iminuit import Minuit
        from mujpy.tools.tools import rebin, derange, write_csv, stringify_groups
        from mujpy.tools.tools import int2min_multigroup, int2_multigroup_method_key
        
        # self.log('Single run, global multigroup')      
        

        a,e = self.suite.asymmetry_multigroup() # the first dimension  is group  the last is time bins
        start, stop, pack = returntup
        time,asymm,asyme = rebin(self.suite.time,a,[start,stop],pack,e=e) #  same slice for all groups

        values,errors,fixed,limits,names,pospar = int2min_multigroup(self.dashboard["userpardicts_guess"])
        # values, errors etc. corrispond to Minuit parameters, only the user defined ones
        #
        # dashboard must contain "userpardicts_guess":lstpars, a list of pardicts, one per user defined parameter
        # The layout consists in defining 
        #     all Minuit parameters, including fixed ones (e.g. phase = 0 for ZF) as user parameters
        #                  each with pardict {"name":_,"value":_,"flag":_,"error":_,"limits":[_,_]} 
        #     all usual components with parameters defined by "flag":"=" and "function":"eval command"
        #     if "function":"" another key must be present
        #        "function_multi":["eval command 0", ... ], a different command for each group
        #        command can reference only the user parameters
         

#        for k in range(len(fitvalues)):
#            self.log('dofit_singlerun_multigroup_userpardicts mufit DEBUG:')
#            self.log(' {} = {}, step = {}, fix = {}, limits ({},{})'.format(names[k],values[k],errors[k],fixed[k],limits[k][0],limits[k][1]))
           
        krun = 0  #  single run!!
        string = []
        methods_keys = int2_multigroup_method_key(self.dashboard,self._the_model_)
        # as many as the total component parameters for one group
        if not methods_keys:
            self.log('Dashboard incompatible with single run multi group single chi2 fit')
            self.log('Check that component parameters are all defined through ''function''')
            self.log('                           with at least one ''function_multi''')
            return
        ok, errmsg = self._the_model_._load_data_multigroup_(time,asymm,methods_keys,e=asyme)
                                              # time 1d, asymm 2d, alpha list
                                              # pass data to model, one at a time
        ############################## int2_multigroup_method_key(dash,model,nruns) returns 
                                     #   list of [methods, keys] to calculate the 2d components
                                     #   in the single migrad call
        if not ok:
            self.log('Error in _load_data_multigroup_: '+errmsg)
            self.log('mufit stops here')            
            return

        self.lastfit = Minuit(self._the_model_._chisquare_,
                              name=names,
                              *values)                                        
        self.lastfit.errors = errors
        self.lastfit.limits = limits
        self.lastfit.fixed = fixed
#        self.lastfit.grad = self._the_model_._add_multirun_grad_
        self.number_dof = asymm.size - self.lastfit.nfit
        # print('mufit dofit_singlerun_multigroup_userpardicts debug: name value error limits fixed {}'.format([[name,value,error,limit,fix] for name,value,error,limit,fix in zip(names,values,errors,limits,fixed)]))
        if self.dofit:
            self.lastfit.migrad()
            # check if some parameters are positive parity 
            if pospar:
                for k in pospar:
                    self.lastfit.limits[k] = [None,None]                    
                self.lastfit.migrad()
            self.lastfit.hesse()
        self.lastfits.append(self.lastfit) # not really needed here

        if self.dofit:
            # write summary on console
            self.summary_global(start,stop)

            # record result in csv file
            version = self.dashboard["version"]
            strgrp = stringify_groups(self.suite.groups)
            modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
            file_csv = self.suite.__csvpath__+modelname+'.'+version+'.'+strgrp+'.csv'
            the_run = self.suite._the_runs_[0][0]
            filespec = self.suite.datafile[-3:]
            header, row = self.prepare_csv()
            string1, string2 = write_csv(header,row,the_run,file_csv,filespec,scan=self.scan)
            # self.log(string1)
            #self.log(string2)
            krun = 0
            self.save_fit_multigroup(krun,string2)
            if not self.lastfit.valid:
                self.log('**** Minuit did not converge! ****')
                print(self.lastfit)

    def dofit_calib_singlerun_multigroup_userpardicts(self,returntup):
        '''
        performs calib fit on single run, multiple groups global
        (A21-calib) tested
        '''
        from iminuit import Minuit
        from mujpy.tools.tools import int2min_multigroup, int2_multigroup_method_key 
        from mujpy.tools.tools import rebin_decay, write_csv, stringify_groups
        
        # self.log('Multigroup calib global: does not work yet')      

        yf,yb,bf,bb,_,_ = self.suite.single_multigroup_for_back_counts(self.suite._the_runs_[0],self.suite.grouping) 
                              # the second dimension is group
        start, stop, pack = returntup
        t,yf,yb,bf,bb,yfm,ybm = rebin_decay(self.suite.time,yf,yb,bf,bb,[start,stop],pack)
        # print('mufit dofit_calib_singlerun_multigroup_userpardicts debug: yfm {}, ybm {}'.format(yfm,ybm))

        values,errors,fixed,limits,names, pospar = int2min_multigroup(self.dashboard["userpardicts_guess"])

        methods_keys = int2_multigroup_method_key(self.dashboard,self._the_model_) 
        # as many as the total component parameters for one group, includes "al"!

        # print('dofit_calib_singlerun_singlegroup mufit debug: fitvalues = {}'.format(fitvalues))
#        for k in range(len(fitvalues)):
#            self.log('{} = {}, step = {}, fix = {}, limits ({},{})'.format(names[k], values[k],errors[k],fixed[k],limits[k][0],limits[k][1]))

        self._the_model_._load_data_calib_multigroup_(t,yf,yb,bf,bb,yfm,ybm,methods_keys)

        self.lastfit = Minuit(self._the_model_._chisquare_,
                              name=names,
                              *values) 
        # print('dofit_calib_singlerun_singlegroup mufit debug: fitvalues = {}'.format(fitvalues))                                       
        self.lastfit.errors = errors
        self.lastfit.fixed = fixed
        self.lastfit.limits = limits
        # self.freepars = self.lastfit.nfit
        self.number_dof = len(t) - self.lastfit.nfit
        if self.dofit:
            self.lastfit.migrad()
            # check if some parameters are positive parity 
            if pospar:
                for k in pospar:
                    self.lastfit.limits[k] = [None,None]                    
                self.lastfit.migrad()
            self.lastfit.hesse()
        self.lastfits.append(self.lastfit) # not really needed

        if self.dofit and self.lastfit.valid:
            pardict = self.dashboard["model_guess"][0]["pardicts"][0]
            p = self.lastfit.values
            for kgroup,group in enumerate(self.suite.grouping):
                group["alpha"] = eval(pardict["function_multi"][kgroup])
                self.suite.groups[kgroup]["alpha"] = eval(pardict["function_multi"][kgroup])            
            # write summary on console
            self.summary_global(start,stop)

            # record result in csv file
            version = self.dashboard["version"]
            strgrp = stringify_groups(self.suite.groups)
            modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
            file_csv = self.suite.__csvpath__+modelname+'.'+version+'.'+strgrp+'.csv'
            the_run = self.suite._the_runs_[0][0]
            filespec = self.suite.datafile[-3:]
            header, row = self.prepare_csv()
            string1, string2 = write_csv(header,row,the_run,file_csv,filespec,scan=self.scan)
            # self.log(string1)
            #self.log(string2)
            krun = 0
            self.save_fit_multigroup(krun,string2)   
         
        elif self.dofit:
            self.log('**** Minuit did not converge! ****')
            print(self.lastfit)
    
    def dofit_calib_singlerun_multigroup_sequential(self,returntup):
        '''
        performs calib fit on single run, multiple groups sequentially
        (A20-calib) tested
        '''
        from iminuit import Minuit
        from mujpy.tools.tools import int2min, int2_method_key, rebin_decay, write_csv, min2int
        
        # self.log('Multigroup calib: does not work yet')
        string = []
        for kgroup,group in enumerate(self.suite.grouping):
            yf,yb,bf,bb,yfm,ybm = self.suite.single_for_back_counts(self.suite._the_runs_[0],group) 
                                  # the second dimension is group
            start, stop, pack = returntup
            t,yf,yb,bf,bb,yfm,ybm = rebin_decay(self.suite.time,yf,yb,bf,bb,[start,stop],pack)

            values,errors,fixed,limits,names,pospar = int2min(self.dashboard["model_guess"])

            # print('dofit_calib_singlerun_singlegroup mufit debug: fitvalues = {}'.format(fitvalues))
    #        for k in range(len(fitvalues)):
    #            self.log('{} = {}, step = {}, fix = {}, limits ({},{})'.format(names[k], values[k],errors[k],fixed[k],limits[k][0],limits[k][1]))

            self._the_model_._load_calib_single_data_(t,yf,yb,bf,bb,yfm,ybm,
                                    int2_method_key(self.dashboard,self._the_model_))
                                                 # int2_int() returns a list of methods to calculate the components

            self.lastfit = Minuit(self._the_model_._chisquare_,
                                  name=names,
                                  *values) 
            # print('dofit_calib_singlerun_singlegroup mufit debug: fitvalues = {}'.format(fitvalues))                                       
            self.lastfit.errors = errors
            self.lastfit.fixed = fixed
            self.lastfit.limits = limits
            self.number_dof = len(t) - self.lastfit.nfit
            if self.dofit:
                self.lastfit.migrad()
            # check if some parameters are positive parity 
                if pospar:
                    for k in pospar:
                        self.lastfit.limits[k] = [None,None]                    
                    self.lastfit.migrad()
                self.lastfit.hesse()
            self.lastfits.append(self.lastfit)

            if self.dofit and self.lastfit.valid:
                self.suite.groups[kgroup]["alpha"] = self.lastfit.values[0]
                self.suite.grouping[kgroup]["alpha"] = self.lastfit.values[0]
                # write summary on console

                self.summary(start, stop, t[1]-t[0],kgroup)  

                # record result in csv file
                version = self.dashboard["version"]
                group = self.suite.groups[kgroup] # assumes only one group
                fgroup, bgroup, alpha = group['forward'],\
				                        group['backward'],\
				                        group['alpha']
                strgrp = fgroup.replace(',','_')+'-'+bgroup.replace(',','_')
                modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
                file_csv = self.suite.__csvpath__+modelname+'.'+version+'.'+strgrp+'.csv'
                the_run = self.suite._the_runs_[0][0]
                filespec = self.suite.datafile[-3:]
                
                header, row = self.prepare_csv() # DEBUG
                
                string1, string2 = write_csv(header,row,the_run,file_csv,filespec,scan=self.scan) # DEBUG

                # self.log(string1)
                #self.log(string2)
                krun = 0
                names, values, stds = min2int(self.dashboard["model_guess"],
						            self.lastfit.values,self.lastfit.errors)
                self.names.append(names)
                self.values.append(values)
                self.stds.append(stds)
                self.fvals.append(self.lastfit.fval)
                string.append(string2)
                # just a check (maybe can be removed):    
            elif self.dofit:
                self.log('**** Minuit did not converge! ****')
                print(self.lastfit)

        if self.dofit:
            self.save_fit_multigroup(krun,string)  # DEBUG

    def dofit_multirun_singlegroup_userpardicts(self,returntup):          
        '''
        performs global fit of many-run single-group data
        (C1) WIP, strategy:
        userpardicts is a list of user parameter dictionaries
            each composed of keys 
            "name", "value", "error" (step), "limits" (default [None,None]), "pospar" (positive parity), "local" (default False)
        the default userpardict "Local":False is global
        if "local":True the user parameter (parent) generates automatically as many 
            daughter parameters as there are runs in the suite (musrfit-style)
        component parameters can be 
            equal to a global parameter or a function of global parameters
            equal to a previous local parameter or a function of global and local parameters
                both cases do not introduce a new minuit parameter and are dealt with by functions
            active, therefore local 
                i.e. the parent component parameter generates automatically as many 
                daughter parameters as there are runs in the suite
        '''
        from iminuit import Minuit
        from mujpy.mucomponents.mucomponents import mumodel
        from mujpy.tools.tools import int2min_multirun, int2_multirun_user_method_key 
        from mujpy.tools.tools import int2_multirun_grad_method_key
        from mujpy.tools.tools import minglobal2sequential, int2min, int2_method_key
        from mujpy.tools.tools import rebin, stringify_groups #, _available_gradients_
        from numpy import array
        from time import time as timeit 

        kgroup = 0 # only one group
        a,e = self.suite.asymmetry_multirun(kgroup) # the second dimension is run 
        start, stop, pack = returntup
        time,asymm,asyme = rebin(self.suite.time,a,[start,stop],pack,e=e)

        values_in,errors,fixed,limits,names,pospar = int2min_multirun(self.dashboard,self.suite.runs)

        string = []
        method_key = int2_multirun_user_method_key(self.dashboard,self._the_model_,len(self.suite._the_runs_))
        ok, errmsg = self._the_model_._load_data_multirun_user_(
                                    time,asymm,method_key,e=asyme) 
        #gradmthd_key = int2_multirun_grad_mthdkey(self.dashboard,self._the_model_,len(self.suite._the_runs_))
        if not ok:
            self.log(repr(errmsg))
            return
        # print('debug mufit dofit_multirun_singlegroup_userpardicts: names =\n{}\npospar =\n{}\nvalues_in =\n{}'.format(names,pospar,values_in))
        if self.grad:
            self._the_model_._load_data_multirun_grad_(int2_multirun_grad_method_key(self.dashboard,self._the_model_,len(self.suite._the_runs_)))
            self.lastfit = Minuit(self._the_model_._chisquare_,
                              name=names,
                              grad = self._the_model_._add_multirun_grad_,                          
                              *values_in
                                )                                        
        else:
            self.lastfit = Minuit(self._the_model_._chisquare_,
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
                grad = self._the_model_._add_multirun_grad_(*values_in)
                from numpy import set_printoptions as npopt, array
                npopt(threshold=1000)
                print('debug grad components as per minuit internal parameter index:') 
                print(grad)
        if self.dofit and self.lastfit.valid: 
            if pospar:
                self.log('... now redo without limits')
                #print('debug mufit dofit_multirun_singlegroup_userpardicts pospar = {}'.format(pospar))
                for k in pospar:
                    #print('debug mufit dofit_multirun_singlegroup_userpardicts k = {}, par = {}'.format(k,names[k]))
                    self.lastfit.limits[k] = [None,None]                    
                tic = timeit()                
                self.lastfit.migrad()
                toc =timeit()-tic
                self.log('migrad no limits redone in {} s, {} calls, {} grads'.format(toc,self.lastfit.nfcn,self.lastfit.ngrad))
        if self.dofit and self.lastfit.valid: 
            tic = timeit()                
            self.lastfit.hesse()
            tuc =timeit()-tic
            self.log('hesse in {} s, {} calls, {} grads'.format(tuc,self.lastfit.nfcn,self.lastfit.ngrad))
            
            n_runs = len(self.suite._the_runs_)

            # write summary on console
            self.summary_multirun_global(start,stop,time[1]-time[0])

            # record result in csv file
            version = self.dashboard["version"]
            strgrp = stringify_groups(self.suite.groups)
            modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])

            # this is a one-shot csv, not incremental
            file_csv = self.suite.__csvpath__+modelname+'.'+version+'.'+strgrp+'.csv'
            self.write_multirun_user_csv(file_csv,scan=self.scan)
            self.save_fit_multirun()
            self.lastfits.append(self.lastfit) # not really needed, for consistency
        if self.dofit and not self.lastfit.valid:
            self.log('**** Minuit did not converge! ****')
            print(self.lastfit)


    def dofit_multirun_singlegroup_sequential(self,returntup,chain):
        '''
        performs sequential fit on many-run, single-group data
        (B1) tested
        '''
        from iminuit import Minuit
        from mujpy.tools.tools import int2min, int2_method_key, rebin, write_csv

        # print('dofit_multirun_singlegroup_sequential mufit debug')
        # self.log('In sequential single')   
        a, e = self.suite.asymmetry_multirun(0) # runs to loaded, group index
        # a, e are 2d: (run,timebin) 
        start, stop, pack = returntup
        time,asymms,asymes = rebin(self.suite.time,a,[start,stop],pack,e=e)
        # time (1d): (timebin)    asymms, asymes (2d): (run,timebin) 

        values,errors,fixed,limits,names,pospar = int2min(self.dashboard["model_guess"])

#        for k in range(len(fitvalues)):
#            self.log('{} = {}, step = {}, fix = {}, limits ({},{})'.format(names[k],values[k],errors[k],fixed[k],limits[k][0],limits[k][1]))

        kgroup = 0
        krun = -1
        for asymm, asyme in zip(asymms,asymes): 
            krun += 1
            self._the_model_._load_data_(time,asymm,int2_method_key(self.dashboard,self._the_model_),
                                     e=asyme) 
                                    # int2_int() returns a list of methods to calculate the components

            self.lastfit = Minuit(self._the_model_._chisquare_,
                              name=names,
                              *values)                                        
            self.lastfit.errors = errors
            self.lastfit.limits = limits
            self.lastfit.fixed = fixed
            # self.freepars = self.lastfit.nfit
            self.number_dof = len(asymm) - self.lastfit.nfit
            if self.dofit:
                self.lastfit.migrad()
                # check if some parameters are positive parity 
                if pospar:
                    for k in pospar:
                        self.lastfit.limits[k] = [None,None]                    
                    self.lastfit.migrad()
                self.lastfit.hesse()
            self.lastfits.append(self.lastfit)

            if self.dofit:
        # write summary on console
                self.summary_sequential(start, stop, time[1]-time[0],k=krun)

            # record result in csv file
                version = self.dashboard["version"]
                group = self.suite.groups[kgroup] # 
                fgroup, bgroup, alpha = group['forward'],\
        					            group['backward'],\
        					            group['alpha']
                strgrp = fgroup.replace(',','_')+'-'+bgroup.replace(',','_')
                modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
                file_csv = self.suite.__csvpath__+modelname+'.'+version+'.'+strgrp+'.csv'
                the_run = self.suite._the_runs_[krun][0]
                # print(the_run.get_runNumber_int())
                filespec = self.suite.datafile[-3:]
                header, row = self.prepare_csv(krun=krun)
                string1, string2 = write_csv(header,row,the_run,file_csv,filespec,scan=self.scan)
                #self.log(string1)
                kgroup = 0
                self.save_fit(krun,kgroup,string2)
                if (chain):
                    values = self.lastfit.values
           
    def dofit_multirun_multigroup_userpardicts(self,returntup):
        '''
        performs single global fit of many-run, many-group data
        (C2) not yet
        '''
            
    def dofit_sequentialrun_multigroup_userpardicts(self,returntup,chain):
        '''
        performs sequential mani-run, multigroup-global fits over a suite of runs
        (B2) testing
        '''
        from iminuit import Minuit
        from mujpy.tools.tools import int2min_multigroup, int2_multigroup_method_key, rebin
        from mujpy.tools.tools import stringify_groups, write_csv
        from numpy import where, array, finfo, sqrt
        
        from matplotlib.pyplot import subplots, draw 

        # print('dofit_multirun_singlegroup_sequential mufit debug')
        # self.log('In sequential single')   
        a, e = self.suite.asymmetry_multirun_multigroup() # runs to loaded, group index
        # a, e are 2d: (run,timebin) 
        print('mufit dofit_sequentialrun_multigroup_userpardicts debug: shape asymm, asyme = {}, {}'.format(a.shape,e.shape))
        start, stop, pack = returntup
        time,asymmrg,asymerg = rebin(self.suite.time,a,[start,stop],pack,e=e)
        zer = array(where(asymerg<2e-162))
        # time (1d): (timebin)    asymms, asymes (2d): (run,timebin) 
        values,errors,fixed,limits,names,pospar = int2min_multigroup(
                                            self.dashboard["userpardicts_guess"])

#        for k in range(len(fitvalues)):
#            self.log('{} = {}, step = {}, fix = {}, limits ({},{})'.format(names[k],values[k],errors[k],fixed[k],limits[k][0],limits[k][1]))
        print('mufit dofit_sequentialrun_multigroup_userpardicts debug: Minuit inputs')
        j = -1
        for ns,vs,es,fx,lm in zip(names,values,errors,fixed,limits):
            j +=1
            print('{} {} = {}({}), {}, {} '.format(j,ns,vs,es,fx,lm))
        
        methods_keys = int2_multigroup_method_key(self.dashboard,self._the_model_)
        # print('mufit dofit_sequentialrun_multigroup_userpardicts debug: methods_keys contains {} methods with{} keys/method'.format(len(methods_keys),[len(c) for g in methods_keys for c in g[1]]))
        krun = -1
        
        
        if self.dofit:
            fig,ax = subplots()
            da, ms, lw = 0.2, 0.1, 0.3
        
        for asymm, asyme in zip(asymmrg,asymerg): 
            krun += 1
            for kg in range(asyme.shape[0]):
                
                if array(where(asyme[kg,:]==0)).sum():
                    print('mufit dofit_sequentialrun_multigroup_userpardicts debug: check asyme[{},{}] contains zeros!'.format(krun,kg))
                    print('mufit dofit_sequentialrun_multigroup_userpardicts debug: asymm.shape {} '.format(asymm.shape))
            # asymm is 2d (group, bins)
            self._the_model_._load_data_multigroup_(time,asymm,methods_keys,e=asyme) 
                                    # int2_int() returns a list of methods to calculate the components

            if self.dofit:
                fs = self._the_model_._add_(time,*values)
                print('mufit dofit_sequentialrun_multigroup_userpardicts debug: fs.shape  {} '.format(fs.shape)) 
                kk, line, fmt,  = -1, ['b-','g-'],['r.','m.']
                for a,e,f in zip(asymm,asyme,fs):
                    kk+=1
                    ax.errorbar(time,a+krun*da,yerr=e,fmt=fmt[kk],ms=ms,alpha=0.3)
                    ax.plot(time,f+krun*da,line[kk],lw=lw,alpha=0.8)

            self.lastfit = Minuit(self._the_model_._chisquare_,
                              name=names,
                              *values)                                        
            self.lastfit.errors = errors
            self.lastfit.limits = limits
            self.lastfit.fixed = fixed
            # self.freepars = self.lastfit.nfit
            self.number_dof = asymm.size - self.lastfit.nfit
            # print('mufit dofit_sequentialrun_multigroup_userpardicts debug: name value error limits fixed {}'.format([[name,value,error,limit,fix] for name,value,error,limit,fix in zip(names,values,errors,limits,fixed)]))
            if self.dofit:            
                print('mufit dofit_sequentialrun_multigroup_userpardicts debug: limits {}'.format(self.lastfit.limits))
                self.lastfit.migrad()
                # check if some parameters are positive parity 
                if pospar:
                    for k in pospar:
                        self.lastfit.limits[k] = [None,None]                    
                    self.lastfit.migrad()
                self.lastfit.hesse()
            self.lastfits.append(self.lastfit)

            if self.dofit:
        # write summary on console
                self.summary_global(start,stop,krun)
                print('mufit dofit_sequentialrun_multigroup_userpardicts debug: fval {}, ndof {}'.format(self.lastfit.fval,self.number_dof))

                version = self.dashboard["version"]
                strgrp = stringify_groups(self.suite.groups)
                modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
                file_csv = self.suite.__csvpath__+modelname+'.'+version+'.'+strgrp+'.csv'
                the_run = self.suite._the_runs_[krun][0]
                filespec = self.suite.datafile[-3:]
                header, row = self.prepare_csv(krun=krun)
                string1, string2 = write_csv(header,row,the_run,file_csv,filespec,scan=self.scan)
                self.log(string1)
                #self.log(string2)
                self.save_fit_multigroup(krun,string2)


                if (chain):
                    values = self.lastfit.values
                
        if self.dofit:                
            ax.set_xlim(0,4)
            ax.set_ylim(-0.5,2.7)
            draw()

    def global_fit(self):
        '''
        True for fit type C
        False for fit types A and B
        '''
        return self.dashboard['model_guess'][0]['pardicts'][0].__contains__('local')
        
    def summary(self,start, stop, dt, kgroup, krun=0):
        '''
        input: k is index in _the_runs_, default 0
        initial version: prints single fit single group result
        '''
        # strategy: two outputs, archive and instant
        # Archive, to replicate the converged fit
        # save a run specific json dashboard with added 
        #          best fit values, 
        #          local chisquare, 
        #          grp_calib dictionary
        # Instant summary   
        # explore: write to a text file cache/log.txt (need self.suite.__cachedir__) 
        #          open a lab terminal
        #          display  it side-by-side (open it from lab, drag its tab to the right)
        #          cd to cache  !!! os dependent
        #          execute tail -f -n +1 log.txt  !!! os dependent
        # The Log Console option would be much better if one could drag it to the right         
        #          it would be os independent
        from mujpy.tools.tools import get_grouping, get_title, chi2std
        from mujpy.tools.tools import len_print_components, print_components, min2int, value_error
        from datetime import datetime

        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        version = self.dashboard["version"]
        the_run = self.suite._the_runs_[krun][0]
        nrun = the_run.get_runNumber_int()
        title = get_title(the_run)
        group = self.suite.groups[kgroup] # assumes only one group
        fgroup, bgroup, alpha = group['forward'],\
						        group['backward'],\
						        group['alpha']

        strgrp = fgroup.replace(',','_')+'-'+bgroup.replace(',','_')
        chi = self.lastfit.fval /self.number_dof 
        # print('summary mufit debug FCN = {}, number of DOF = {}'.format(self.lastfit.fval,self.number_dof))
        lowchi, highchi = chi2std(self.number_dof)
        start, stop = self.suite.time[start]*1000, self.suite.time[stop]
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")  

        file_log = self.suite.__cachepath__+modelname+'.'+str(nrun)+'.'+strgrp+'.'+version+'.log'

        names, values, errors = min2int(self.dashboard["model_guess"],
							        self.lastfit.values,self.lastfit.errors)
        #print('mufit summary debug: values, errors {}'.format([[name,value,error] for name,value,error in zip(names,values,errors)]))
        #print('mufit summary debug: minvalues, minerrors {}'.format([[name,value,error] for name,value,error in zip(self.lastfit.parameters,self.lastfit.values,self.lastfit.errors)]))        
        with open(file_log,'w') as f:
            f.write(' '+86*'_'+'\n')
            f.write('| Run {}: {}      α = {:.3f} on group: {} - {}                      |\n'.format(nrun,
		                                 title,alpha,fgroup,bgroup))
            f.write('| χᵣ² = {:.3f}({:.3f},{:.3f}), fit on [{:.2f}ns, {:.2}µs] {:.2f}ns/bin                           |\n'.format(chi,lowchi,highchi,start,stop,dt*1000))
            f.write('|'+86*'-'+'|\n') 
            # for k,name in enumerate(names): # replaced by:
            maxlen = 0
            for name,value,error in zip(reversed(names),reversed(values),reversed(errors)):
                maxlen = max(maxlen, len_print_components(name, value, error))    
            for name,value,error in zip(names,values,errors):
                f.write('| '+print_components (name, value, error,maxlen)+'\n')
#                self.log('mufit summary debug: {} = {}, error = {}'.format(name, value,error))
            f.write(86*'-'+'\n')
            f.write('             Best fit performed on {}'.format(dt_string))
            self.log(78*'-')
            maxlen = 0
            for name,value,error in zip(reversed(names),reversed(values),reversed(errors)):
                maxlen = max(maxlen, len_print_components(name, value, error))    
            for name,value,error in zip(reversed(names),reversed(values),reversed(errors)):
                self.log('| '+print_components(name, value, error,maxlen))
            self.log('|'+76*'-'+'|') 
            self.log('|χᵣ² = {:.3f}({:.3f},{:.3f}), on [{:.2f}ns, {:.2}µs] {:.2f}ns/bin @{}|'.format(chi,
		                                 lowchi,highchi,start,stop,dt*1000,dt_string))
            self.log('| Run {}: {}         α = {:.3f} on group: {} - {}         |'.format(nrun,
		                                 title,alpha,fgroup,bgroup))
            self.log(' '+77*'_')

    def summary_sequential(self, start, stop, dt, k=0):
        '''
        input: k is index in _the_runs_, default 0
        initial version: prints single fit single group result
        '''
        from mujpy.tools.tools import get_title, chi2std, len_print_components, print_components, min2int
        from datetime import datetime

        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        version = self.dashboard["version"]
        the_run = self.suite._the_runs_[k][0]
        nrun = the_run.get_runNumber_int()
        title = get_title(the_run)
        group = self.suite.groups[0] # assumes only one group
        fgroup, bgroup, alpha = group['forward'],\
    					        group['backward'],\
    					        group['alpha']
        strgrp = fgroup.replace(',','_')+'-'+bgroup.replace(',','_')
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")  
        start, stop = self.suite.time[start]*1000, self.suite.time[stop]
        if k==0:
            self.log(' '+71*'_')
            fit_string = '| Fit on [{:.2f}ns, {:.2}µs, {:.2f}ns/bin] on group: {} - {}  α = {:.3f}, {} |'
            self.log(fit_string.format(start,stop,dt,fgroup,bgroup,alpha,dt_string))
        chi = self.lastfit.fval /self.number_dof 
        lowchi, highchi = chi2std(self.number_dof)
        file_log = self.suite.__cachepath__+modelname+'.'+str(nrun)+'.'+strgrp+'.'+version+'.log'
        names, values, errors = min2int(self.dashboard["model_guess"],
							        self.lastfit.values,self.lastfit.errors)
        with open(file_log,'w') as f:
            f.write(' '+85*'_'+'\n')
            f.write('| Run {}: {}     on group: {} - {}  α = {:.3f}     |\n'.format(nrun,
		                                 title,fgroup,bgroup,alpha))
            self.log('| Run {}: {}         χᵣ² = {:.3f}({:.3f},{:.3f})'.format(nrun,
		                             title,chi,lowchi,highchi,))
            f.write('| χᵣ² = {:.3f}({:.3f},{:.3f}), fit on [{:.2f}ns, {:.2}µs, {:.2f}ns/bin], at {}|\n'.format(chi,
		                                 lowchi,highchi,start,stop,dt*1000,dt_string))
            f.write('|'+85*'-'+'|\n') 
            self.log('|'+71*'-'+'|') 
            maxlen = 0
            for name,value,error in zip(names,values,errors): 
                maxlen = max(maxlen, len_print_components(name, value, error))    
            for name,value,error in zip(names,values,errors): 
                f.write('| '+print_components(name, value, error,maxlen)+'\n')
                self.log('| '+print_components(name, value, error,maxlen))
#                self.log('summary_sequential mufit debug: {} = {}, error = {}'.format(name, value,error))
            f.write('|'+85*'_'+'|\n')
            f.write('             Best fit performed on {}'.format(dt_string))
            self.log('|'+71*'_'+'|')

    def summary_global(self,start,stop,krun=0):
        '''
        input: krun is index in _the_runs_, default 0
        initial version: prints multigroup globa fit result
        '''
        from mujpy.tools.tools import get_title, chi2std, stringify_groups
        from mujpy.tools.tools import len_print_components, print_components, min2int_multigroup
        from datetime import datetime

        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        version = self.dashboard["version"]
        the_run = self.suite._the_runs_[krun][0]
        nrun = the_run.get_runNumber_int()
        title = get_title(the_run)
        strgrp = stringify_groups(self.suite.groups)
        chi = self.lastfit.fval /self.number_dof 
        lowchi, highchi = chi2std(self.number_dof)
        start, stop = self.suite.time[start]*1000, self.suite.time[stop]
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")  

        file_log = self.suite.__cachepath__+modelname+'.'+str(nrun)+'.'+strgrp+'.'+version+'.log'
        names, values, errors = min2int_multigroup(self.dashboard,
							        self.lastfit.values,self.lastfit.errors,nruns)
        # list (groups) of lists (omponents) of lists (parameters)
        sumlength = 122
        with open(file_log,'w') as f:
            f.write(' '+96*'_'+'\n')
            nch = sumlength - 2
            self.log(' '+nch*'_')
            string = '| Run {}: {}    Global fit of {}'.format(nrun,title,dt_string)
            f.write(string+39*' '+'|\n')
            self.log(string+14*' '+'|')
            string = '| χᵣ² = {:.3f}({:.3f},{:.3f}) ,    on [{:.2f}ns, {:.2}µs, {:.2f}ns/bin]'.format(chi,lowchi,highchi,start,stop,dt*1000)
            f.write(string+46*' '+'|\n')
            nch = sumlength - 2 - len(string) if sumlength-len(string) - 2 >=0 else sumlength - 2
            self.log(string+nch*' '+'|')
            for g1,n1,v1,e1,g2,n2,v2,e2 in zip(self.suite.groups[::2],names[::2],values[::2],errors[::2],
                                               self.suite.groups[1::2],names[1::2],values[1::2],errors[1::2]):
                fg1,bg1,al1 = g1['forward'], g1['backward'], g1['alpha'] 
                fg2,bg2,al2 = g2['forward'], g2['backward'], g2['alpha'] 
                string = ' on group: {} - {}   α = {:.3f}   |'.format(fg1,bg1,al1)
                nch = sumlength - 1 - len(string) if sumlength-len(string) - 2 >=0 else sumlength - 2
                f.write('|'+nch*'-'+string+'\n')
                self.log('|'+nch*'-'+string)
            # for k,name in enumerate(names): # replaced by:
                mexlen = 0
                for nam,val,err in zip(n1,v1,e1):
                    maxlen = max(maxlen, len_print_components(nam, val, err))
                for nam,val,err in zip(n1,v1,e1):
                    f.write('| '+print_components (nam, val, err,maxlen)+'\n')
                    self.log('| '+print_components(nam, val, err,maxlen))
#                self.log('mufit summary debug: {} = {}, error = {}'.format(name, value,error))
                string = ' on group: {} - {}   α = {:.3f}   |'.format(fg2,bg2,al2)
                nch = sumlength - 1 - len(string) if sumlength-len(string) - 2 >=0 else sumlength - 2
                f.write('|'+nch*'-'+string+'\n') 
                self.log('|'+nch*'-'+string)
                mexlen = 0
                for nam,val,err in zip(n2,v2,e2):
                    maxlen = max(maxlen, len_print_components(nam, val, err))
                for nam,val,err in zip(n2,v2,e2):
                    f.write('| '+print_components (nam, val, err,maxlen)+'\n')
                    self.log('| '+print_components(nam, val, err,maxlen))
#                self.log('mufit summary debug: {} = {}, error = {}'.format(name, value,error))
            nch = sumlength - 2
            f.write('|'+nch*'_'+'|\n')
        nch = sumlength - 2
        self.log('|'+nch*'_'+'|')

    def summary_multirun_global(self,start,stop,dt):
        '''
        print summary on Output and log file
        multirun user version
        '''
        from mujpy.tools.tools import get_title, chi2std, stringify_groups, value_error
        from mujpy.tools.tools import len_print_components_multirun, print_components_multirun, min2int_multirun
        from datetime import datetime

        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        version = self.dashboard["version"]
        nrun0 = self.suite._the_runs_[0][0].get_runNumber_int()
        nrun1 = self.suite._the_runs_[-1][0].get_runNumber_int()
        title = get_title(self.suite._the_runs_[0][0])
        strgrp = stringify_groups(self.suite.groups)
        chi = self.lastfit.fval /self.number_dof 
        lowchi, highchi = chi2std(self.number_dof)
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
                        # these are the global user parameters
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


    def prepare_csv(self,krun = 0):
        '''
        input: 
            k is index in _the_runs_, default k = 0
        output: 
            header, the model specific csv header 
                    to compare with that of the csv file
            row, the line to be added to the csv file
        prepares a csv-like row of best fit parameters 
        that can be imported to produce figures
        Identifies multigroup as dashboard = False in minparam2_csv::
        '''
        from mujpy.tools.tools import get_title, spec_prec, chi2std, initialize_csv 
        from mujpy.tools.tools import minparam2_csv, chi2_csv, min2int, min2int_multigroup
        from mujpy.tools.tools import min2int_multirun, multigroup_in_components, userpars

        # print('k = {}, self.nrun = {}'.format(k,[j for j in self.nrun]))


        filespec = self.suite.datafile[-3:]
        lowchi, hichi = chi2std(self.number_dof)
        chi = self.lastfit.fval /self.number_dof # fval is cost (chi2) at the minimum
        if self.C1():
            rowpars = minparam2_csv(self.dashboard,self.lastfit.values,self.lastfit.errors,multirun=self.suite._the_runs_)
            rows = []
            nruns = len(self.suite._the_runs_)
            dof_run = self.number_dof/nruns
            for run in range(nruns):
                the_run = self.suite._the_runs_[run][0]
                chi_run = self._the_model_._chisquare_single_(*self.lastfit.values,run)/dof_run
                row = initialize_csv(the_run.get_field(), filespec, the_run)
                row += rowpars[run] 
                row += chi2_csv(chi_run,lowchi,hichi,self.suite.groups,self.suite.offset)   
                row += '\n'
                rows.append(row)
            row = rows # in this case, C1, row is a list of rows
        else:
            the_run = self.suite._the_runs_[krun][0]
            # print('mufit prepare_csv debug, krun ={} run={}'.format(krun,the_run.get_runNumber_int()))
            row = initialize_csv(the_run.get_field(), filespec, the_run) # run T B valid for all fits
            dashcsv = False if self.suite.multi_groups() else self.dashboard["model_guess"]
            min2dash = min2int_multigroup if (self.suite.multi_groups() and userpars(self.dashboard)) else min2int
            dashboard =  self.dashboard if (self.suite.multi_groups() and userpars(self.dashboard)) else self.dashboard["model_guess"]        
            row += minparam2_csv(dashcsv,self.lastfit.values,self.lastfit.errors) # multirun=0
            row += chi2_csv(chi,lowchi,hichi,self.suite.groups,self.suite.offset)
            row += '\n'
        # row is formatted with appropriate rounding, write directly
        # self.console(row)

        if filespec == 'bin' or filespec == 'mdu':
            header = ['#0.Run','1.T_cryo[K]','2.e_T_cryo[K]','3.T_sample[K]','4.e_T_sample[K]','5.B[G]']
            k = 5
        else:
            header = ['#0.Run','1.T[K]','2.eT[K]','2.B[G]']
            k = 2
        # now component names for header
        # print('prepare_csv mufit debug:  sum(multigroup_in_components(self.dashboard)) {}'.format(bool(sum(multigroup_in_components(self.dashboard)))))
        if sum(multigroup_in_components(self.dashboard)) and self.suite.single():
            pardicts = self.dashboard['userpardicts_guess']
            for pardict in pardicts:
                k += 1
                header.append(str(k)+'.'+pardict['name'])
                k += 1
                header.append(str(k)+'.'+'e_'+pardict['name'])
        elif self.C1(): # multirun global 
            parlists, _, _ = min2int_multirun(self.dashboard,
							        self.lastfit.values,self.lastfit.errors,self.suite._the_runs_)
            for name in parlists[1]: # list of local parameter names for a run            
                k += 1
                header.append(str(k)+'.'+name)
                k += 1
                header.append(str(k)+'.'+'e_'+name)
            for name in parlists[0]: # list of global parameter names
                k += 1
                header.append(str(k)+'.'+name)
                k += 1
                header.append(str(k)+'.'+'e_'+name)                
        else:            
            parlists, _, _ = min2dash(dashboard,
							        self.lastfit.values,self.lastfit.errors)
            parlists = parlists[0] if isinstance(parlists[0][0],list) else parlists
            for parlist in parlists:
                # print('mufit prepare_csv debug: parlist: {}'.format(parlist))
                for parname in parlist:
                    # print('mufit prepare_csv debug: parnames: {}'.format(parname))
                    k += 1
                    header.append(str(k)+'.'+parname)
                    k += 1
                    header.append(str(k)+'.'+'e_'+parname)
        k += 1
        header.append(str(k)+'.'+'chi2_r')
        header.append('e_chi2_low')
        header.append('e_chi2_hi')
        if self.suite.multi_groups():
            for jgroup,group in enumerate(self.suite.groups):
                header.append('alpha{}'.format(jgroup))
        else:
            header.append('alpha')
        header.append('offset_bin')
        header.append('timestring\n')
        return ','.join(header), row

    def write_multirun_user_csv(self,file_csv,scan=None):
        '''
        this is a one-shot csv write after a global fit
        input :
            the_runs is suite _the_runs_
            file_csv = full path/filename to csv file 
        '''
        from mujpy.tools.tools import get_title, min2int_multirun
        from datetime import datetime

    # prepare_csv writes a header: # column-index-name 
    # plos a list of rows, a row per each run, composed of 
    # run number T eT B 
    # local parameters and their errors (columns)
    # global parameters and their errors (columns with repeated values)
    # chi2 end their errors (partial), chi2 global and its error repeated) 
        names, values, errors = min2int_multirun(self.dashboard,
							            self.lastfit.values,self.lastfit.errors,self.suite._the_runs_)
        header, rows = self.prepare_csv()    
        with open(file_csv,'w') as f_out:  
            f_out.write(header)               
            for line in rows:
                f_out.write(line)
        file_csv = file_csv[file_csv.rfind('/')+1:]
        nrun0 = self.suite._the_runs_[0][0].get_runNumber_int()
        nrun1 = self.suite._the_runs_[-1][0].get_runNumber_int()
        self.log('Global fit of runs {}-{}:'.format(nrun0,nrun1)+
                ' values and errors saved in {}'.format(file_csv))
    
    def save_fit(self,krun,kgroup,string):
        '''
        input:
            krun is index in self.suite._the_runs_
            kgroup is indek in self.suite.groups
        saves a dashboard file adding the bestfit parameters as "model_result"
        These saves are individual runs
        use "version" as additional label to qualify fit
        filename is __cachepath__ + modelname + nrun  + strgrp + version .json
        nrun = runNumber, strgrp = shorthand for group
        '''
        from mujpy.tools.tools import min2int
        import json
        import os
        from copy import deepcopy
                
        version = self.dashboard["version"]
        group = self.suite.groups[kgroup] # assumes only one group
        fgroup, bgroup, alpha = group['forward'],\
				                group['backward'],\
				                group['alpha']
        strgrp = fgroup.replace(',','_')+'-'+bgroup.replace(',','_')
        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        the_run = self.suite._the_runs_[krun][0]
        nrun = str(the_run.get_runNumber_int())
        file_json = self.suite.__fitpath__+modelname+'.'+nrun+'.'+strgrp+'.'+version+'_fit.json'

        # replace (guess) values with Minuit.values, leave error as step, add fit_range, std and chi2
        # do not replace names, they are autogenerated by mufit
        names, values, errors = min2int(self.dashboard["model_guess"],self.lastfit.values,self.lastfit.errors)
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
        self.log('{}  saved'.format(short_json)+string)

    def save_fit_multigroup(self,krun,string_in):
        '''
        input:
            krun is index in self.suite._the_runs_
            string_in is the result of the write_csv process
        if fit is global
            saves a dashboard json adding the bestfit parameters as "userpardicts_result"
        else if sequential
            saves as many single run single group "model_result" dashboard json 
        Use "version" as additional label to qualify fit (auto 'gg_
        filename is __cachepath__ + modelname + nrun  + srtgrp0 + strgrp...  + version .json
        nrun = runNumber, strgrp0,1,... = shorthand for allgroups
        '''
        from mujpy.tools.tools import stringify_groups
        import json
        import os
        from copy import deepcopy
        
        # file name composition        
        # print('save_fit_multigroup mufit debug: dashboard version {}'.format(self.dashboard['version']))
        version = self.dashboard["version"] 
        strgrp = stringify_groups(self.suite.groups)
        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        the_run = self.suite._the_runs_[krun][0]
        nrun = str(the_run.get_runNumber_int())
        # print('debug {}.{}.{}.{}'.format(modelname,nrun,strgrp,version))
        file_json = self.suite.__fitpath__+modelname+'.'+nrun+'.'+strgrp+'.'+version+'_fit.json'
        # dashboard result-augmented
        chi2 = []
        if isinstance(string_in, list): # sequential fit
        # MUST refactor into saving one _fit json dashboard file per group
        # reflects the fact that as many chi2 as groups were minimized
        # use just a "model_result" list for each instead of: 
        # remove =  [], split stringify_groups in bits, update a file_json per group 
        # and bring with open(file_json,"w") within the if
            # use "model_result" as a list of lists
            
            self.dashboard["model_result"] = []
            for string in string_in:
                self.log(string)
            for kgroup in range(len(string_in)):
                fgroup, bgroup = self.suite.groups[kgroup]['forward'],\
                                 self.suite.groups[kgroup]['backward']
                strgrp = fgroup.replace(',','_') + '-' + bgroup.replace(',','_')
                file_json = self.suite.__fitpath__+modelname+'.'+nrun+'.'+strgrp+'.'+version+'_fit.json'
                group_results = deepcopy(self.dashboard["model_guess"])
                value, std = self.values[kgroup], self.stds[kgroup]
                # print('save_fit_multigroup mufit debug: value {}'.format(value))
                for k,component in enumerate(group_results):
                    # print('save_fit_multigroup mufit debug: n components {}'.format(len(component["pardicts"]))) 
                    for j,pardict in enumerate(component['pardicts']):
                        pardict["value"] = value[k][j]
                        pardict["std"] = std[k][j]               
                chi2 = self.fvals[kgroup] /self.number_dof
                self.dashboard["model_result"] = group_results
                self.dashboard["chi2"]=chi2
                with open(file_json,"w") as f:
                    json.dump(self.dashboard,f, indent=2,ensure_ascii=False)   
                short_json = file_json.replace(self.suite.__startuppath__,'./')              
                self.log('{} saved'.format(short_json))
        else: # global
            # userpardicts fit
            names = [parameter["name"] for parameter in self.dashboard["userpardicts_guess"]]
            userpardicts = []
            for name,value,std in zip(names,self.lastfit.values,self.lastfit.errors):
                 userpardicts.append({'name':name,'value':value,'std':std})
            self.dashboard["userpardicts_result"] = userpardicts
            self.dashboard["chi2"] = self.lastfit.fval /self.number_dof
            if os.path.isfile(file_json): 
                os.rename(file_json,file_json+'~')
            with open(file_json,"w") as f:
                json.dump(self.dashboard,f, indent=2,ensure_ascii=False)
            string_in = 'Best fit saved in {} '.format(file_json)+string_in
            self.log(string_in)
# ,object_pairs_hook=OrderedDict)    

    def save_fit_multirun(self):
        '''
        fit is multirun global
            saves a dashboard json adding the bestfit parameters as "userpardicts_result"
            and "model_result"
        to be consistent a single-run model_result is saved
        with lists of values, one per run, in place of single values as in the model_guess 
        filename is __cachepath__ + modelname + nruns + srtgrp + version.json
        nruns = shorthand for runNumbers, strgrp = shorthand for allgroups
        '''
        from mujpy.tools.tools import stringify_groups, min2int_multirun
        import json
        import os
        from copy import deepcopy
        
        # file name composition        
        # print('save_fit_multigroup mufit debug: dashboard version {}'.format(self.dashboard['version']))
        version = self.dashboard["version"] 
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
        for k, pardict in enumerate(self.dashboard['userpardicts_guess']):
            if pardict['local'] or type(pardict['value'])==list:
                n_locals += 1 # number of local user parameters
        self.n_locals = n_locals
        userpardicts = []
        # model indices and names for local component parameters 
        componentindex = [k for k,component in enumerate(model_result) for pardict in component['pardicts'] if pardict['flag']=='~']
        parname =[pardict["name"] for component in model_result for pardict in component['pardicts']  if pardict['flag']=='~']

        for nam,val,err in zip(names[0],values[0],errors[0]): # global parameters
            userpardicts.append({'name':nam,'value':val,'std':err, 'local':False})
            n_globals += 1
        for j,nam in enumerate(names[1]): # names of minuit parameters for first run
            # print('debug mufit save_fit_multirun j = {}, nam = {}, n_locals = {}'.format(j, nam, n_locals))
            # the first n_locals appended to userpardicts
            if j<n_locals: # first ones are user locals
                na = nam.rstrip(digits).rstrip('_') # stripped of run number
                va = [vals[j] for vals in values[1:]] # vals is a run list and val[j] is a user local  
                er = [errs[j] for errs in errors[1:]] # errs is a run list and err[j] is its error
                # va and er ar lists over runs  
                userpardicts.append({'name':na,'value':va,'std':er,'label':'','local':True})
            for component in model_result:
                for pardict in component["pardicts"]:
                    if pardict["flag"] !="=":
                        pardict["name"] = nam.rstrip(digits).rstrip('_') # stripped of run number
                        pardict["value"] = [vals[j] for vals in values[1:]] # vals is a run list and val[j] is a component par  
                        pardict["std"] = [errs[j] for errs in errors[1:]] # errs is a run list and err[j] is its error
                #self.log('debug mufit save_fit_multirun: minuit name = {}, parname = {}'.format(na,model_result[index]["pardicts"]["name"]))
        self.dashboard["userpardicts_result"] = userpardicts
        self.dashboard["model_result"] = model_result
        self.dashboard["chi2"] = self.lastfit.fval /self.number_dof
        if os.path.isfile(file_json): 
            os.rename(file_json,file_json+'~')
        with open(file_json,"w") as f:
            json.dump(self.dashboard,f, indent=2,ensure_ascii=False)
        string_in = 'Best fit saved in {} '.format(file_json)
        self.log(string_in)

    def show_calib(self):
        '''
        output:
            t time 
            a asymmetry
            e asymmetry error
            f guess fit function for calib mode
        for degugging single run calibs
        '''
        from mujpy.tools.tools import int2_method_key, int2min
        run = self.suite._the_runs_[0]
        yf, yb, bf, bb, yfm, ybm = self.suite.single_for_back_counts(run,self.suite.grouping[0])
        t = self.suite.time
        a, e = self.suite.asymmetry_single(run,0) # returned for errorbar of data, in case
        par,_,_,_,name = int2min(self.dashboard["model_guess"])
        self._the_model_._load_calib_single_data_(t,yf,yb,bf,bb,yfm,ybm,
                                                  int2_method_key(self.dashboard,self._the_model_))
        f = self._the_model_._add_calib_single_(t,*par)
        return t,a,e,f

    def calib(self):
        '''
        True if the first component is 'al'
        '''
        return self.dashboard['model_guess'][0]['name']=='al'

    def userpar(self): # this is a global fit using userpardicts
        '''
        see table at the bottom of https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.GlobalSwitch
        '''
        return "userpardicts_guess" in self.dashboard.keys()

    def tilde_in_component(self): # this dashboard has minuit parameters in the model components
        '''
        see table at the bottom of https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.GlobalSwitch
        '''
        # empty list, no ~ flags, is equivalent to False, non empty list is True
        return any([par['flag']=='~' for component in self.dashboard["model_guess"]  for par in component['pardicts']]) 

    def A1(self): # single run singlegroup 
        '''
        see table at the bottom of https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.GlobalSwitch
        '''
        return self.suite.single() and not (self.calib() or self.suite.multi_groups() or self.userpar())
            
    def A1_calib(self): # single run calib singlegroup 
        '''
        see table at the bottom of https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.GlobalSwitch
        '''
        return self.suite.single() and self.calib() and not (self.suite.multi_groups() or self.userpar())
            
    def A20(self): # single run multigroup sequential 
        '''
        see table at the bottom of https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.GlobalSwitch
        '''
        return self.suite.single() and self.suite.multi_groups and not self.userpar()

    def A20_calib(self): # single run calib multigroup sequential 
        '''
        see table at the bottom of https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.GlobalSwitch
        '''
        return self.suite.single() and self.suite.multi_groups and self.tilde_in_component() and self.calib()
                
    def A21(self): # single run multigroup global 
        '''
        see table at the bottom of https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.GlobalSwitch
        '''
        return self.suite.single() and self.suite.multi_groups and self.userpar() and not self.calib()
                
    def A21_calib(self): # single run calib multigroup global 
        '''
        see table at the bottom of https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.GlobalSwitch
        '''
        return self.suite.single() and self.suite.multi_groups and self.userpar() and self.calib()
                
    def B1(self): # multirun sequential singlegroup 
        '''
        see table at the bottom of https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.GlobalSwitch
        '''
        return not (self.suite.multi_groups() or self.userpar() or self.suite.single()) and self.tilde_in_component()

    def B20(self): # multirun sequential multigroup sequential
        '''
        see table at the bottom of https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.GlobalSwitch
        '''
        return not (self.suite.multi_groups() or self.userpar() or self.suite.single()) and self.tilde_in_component()

    def B21(self): # multirun sequential multigroup global
        '''
        see table at the bottom of https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.GlobalSwitch
        '''
        return not (self.suite.multi_groups() or self.userpar() or self.suite.single()) and self.tilde_in_component()

    def C1(self): # multirun global singlegroup
        '''
        see table at the bottom of https://musr-nmr.unipr.it/dispense/pmwiki.php?n=Mujpy.GlobalSwitch
        '''
        return not self.suite.multi_groups() and self.userpar() and self.tilde_in_component()

    def C2(self): # multirun global multigroup global
        return not self.suite.single() and self.userpar() and self.tilde_in_component()
        
        
        
