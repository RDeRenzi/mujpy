class mufit(object):
    '''
    fit class 
    reads from a dashboard file
        can also be used to generate the gui
    '''
    def __init__(self,suite,dashboard_file,chain=False,dash=None,initialize_only=False):
        '''
        input
            suite is the instance of the runs
            dashboard_file is a JSON file of a dictionary structure
        '''
        self.suite = suite
        self.dash = dash
        self.log = self.dash.log if self.dash else self.suite.console
        self.nodash = True
        self.nodata = True
        self.dofit = not initialize_only # initialize_only True just loads guess values, no Minuit
        self._initialise_fit(dashboard_file,chain)       
        
    def _initialise_fit(self,dashboard_file,chain):
        '''
        input:
            json dashboard_file produces a dict structure
        '''
        # from collections import OrderedDict
        from mujpy.aux.aux import _available_components_
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
        from mujpy.aux.aux import derange, add_step_limits_to_model, checkvalidmodel
        from mujpy.aux.aux import model_name, userpars, userlocals, multigroup_in_components
        
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
                
            if not self.suite.multi_groups(): # A1, B1, C1 single group fits
            
                if self.suite.single():       # A1
                
                    if self.calib():          # A1 calib DONE  (True if the first component is 'al')
                        self.dofit_calib_singlerun_singlegroup(0,returntup)  
                    else:                     # A1 DONE
                        self.dofit_singlerun_singlegroup(returntup)
                        # print('mufit choosefit debug: should be finished!')       
                else:                         # B1, C1  
                    if userpars(self.dashboard): # C1
                        self.dashboard["version"]=('gr_'+version if 
                                        self.dashboard["version"][0:3]!='gr_' else version)
                        self.dofit_multirun_singlegroup_userpardicts(returntup) 
                    else:                     # B1 DONE
                        self.dofit_multirun_singlegroup_sequential(returntup,chain)

            else:                             # A2, B2, C2 multi groups fits 
                if self.suite.single():       # A2 (single run  multi group)
                                              # can be calib 
                                              # can be sequential or global
                    if self.calib():          # A2-calib, True if the first component is 'al'
                        if sum(multigroup_in_components(self.dashboard)): # single run calib 
                        
                            self.dashboard["version"]=('gg_'+version if 
                                            self.dashboard["version"][0:3]!='gg_' else version)                           
                            self.dofit_calib_singlerun_multigroup_userpardicts(returntup)  
                                                                          # multi group global A2_01-calib
                        else: 
                            self.dofit_calib_singlerun_multigroup_sequential(returntup)  # DONE (?)
                                                                          # multi group sequential A2_00-calib 
                    else:                     # A2 DONE multi group single run
                        if sum(multigroup_in_components(self.dashboard)):   # A2_1 DONE 
                                              # A2 multi group single run single chi2 global parameters
                            self.dashboard["version"]=('gg_'+version if 
                                            self.dashboard["version"][0:3]!='gg_' else version)
                            self.dofit_singlerun_multigroup_userpardicts(returntup) # True, True, False, True
                        else:                 # A2_0 DONE multi-group-sequentially single run 
                            self.dofit_singlerun_multigroup_sequential(returntup)   # True, True, False, False, True
                else:                         # B2, C2 multi group multirun 
                    if userlocals(self.dashboard): # C2 multi group multirun global 
                        # this is a fit with global parameters
                        self.dashboard["version"]=('gg_'+version if 
                                        self.dashboard["version"][0:3]!='gg_' else version)
                        self.dofit_multirun_multigroup_userpardicts(returntup)      # False, True, False, True, False
                    else:                     # B2 multigroup, sequential multirun (sequential all does not exist)
                        self.dashboard["version"]=('gg_'+version if 
                                        self.dashboard["version"][0:3]!='gg_' else version)
                        self.dofit_sequentialrun_multigroup_userpardicts(returntup,chain) 
                        # pe.g TF sequential of multigroups with shared parameters
        return False

    def dofit_calib_singlerun_singlegroup(self,kgroup,returntup):
        '''
        performs calib fit on single run, single group, 
        (A1-calib) tested
        input 
            kgroup is group index in suitegrouping
        '''
        from iminuit import Minuit
        from mujpy.aux.aux import int2min, int2_method_key, rebin_decay, write_csv
        
        # self.log('In single calib.')      

        yf,yb,bf,bb,yfm,ybm = self.suite.single_for_back_counts(self.suite._the_runs_[0],self.suite.grouping[kgroup]) 
                              # the second dimension is group
        start, stop, pack = returntup
        t,yf,yb,bf,bb,yfm,ybm = rebin_decay(self.suite.time,yf,yb,bf,bb,[start,stop],pack)

        fitvalues,fiterrors,fitfixed,fitlimits,names,pospar = int2min(self.dashboard["model_guess"])

        # print('dofit_calib_singlerun_singlegroup mufit debug: fitvalues = {}'.format(fitvalues))
#        for k in range(len(fitvalues)):
#            self.log('{} = {}, step = {}, fix = {}, limits ({},{})'.format(names[k], values[k],errors[k],fixed[k],limits[k][0],limits[k][1]))

        self._the_model_._load_calib_single_data_(t,yf,yb,bf,bb,yfm,ybm,
                                                  int2_method_key(self.dashboard,self._the_model_))
                                             # int2_int() returns a list of methods to calculate the components
        self.lastfit = Minuit(self._the_model_._chisquare_,
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
        self.lastfits.append(self.lastfit)

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
                
                string1, string2 = write_csv(header,row,the_run,file_csv,filespec) 

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
        from mujpy.aux.aux import int2min, int2_method_key, rebin, write_csv
          
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
        self.lastfits.append(self.lastfit)

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
            string1, string2 = write_csv(header,row,the_run,file_csv,filespec)
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
        from mujpy.aux.aux import int2min, int2_method_key, min2int
        from mujpy.aux.aux import rebin, write_csv
        
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
            # self.freepars = self.lastfit.nfit
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
                string1, string2 = write_csv(header,row,the_run,file_csv,filespec)
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
        from mujpy.aux.aux import rebin, derange, write_csv, stringify_groups
        from mujpy.aux.aux import int2min_multigroup, int2_multigroup_method_key
        
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
        ############################## int2_multigroup_method_key() returns 
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
        # self.freepars = self.lastfit.nfit
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
        self.lastfits.append(self.lastfit)

        if self.dofit:
            # write summary on console
            self.summary_global(start, stop, time[1]-time[0])

            # record result in csv file
            version = self.dashboard["version"]
            strgrp = stringify_groups(self.suite.groups)
            modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
            file_csv = self.suite.__csvpath__+modelname+'.'+version+'.'+strgrp+'.csv'
            the_run = self.suite._the_runs_[0][0]
            filespec = self.suite.datafile[-3:]
            header, row = self.prepare_csv()
            string1, string2 = write_csv(header,row,the_run,file_csv,filespec)
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
        from mujpy.aux.aux import int2min_multigroup, int2_multigroup_method_key 
        from mujpy.aux.aux import rebin_decay, write_csv, stringify_groups
        
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
        self.lastfits.append(self.lastfit)

        if self.dofit and self.lastfit.valid:
            pardict = self.dashboard["model_guess"][0]["pardicts"][0]
            p = self.lastfit.values
            for kgroup,group in enumerate(self.suite.grouping):
                group["alpha"] = eval(pardict["function_multi"][kgroup])
                self.suite.groups[kgroup]["alpha"] = eval(pardict["function_multi"][kgroup])            
            # write summary on console
            self.summary_global(start, stop, t[1]-t[0])

            # record result in csv file
            version = self.dashboard["version"]
            strgrp = stringify_groups(self.suite.groups)
            modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
            file_csv = self.suite.__csvpath__+modelname+'.'+version+'.'+strgrp+'.csv'
            the_run = self.suite._the_runs_[0][0]
            filespec = self.suite.datafile[-3:]
            header, row = self.prepare_csv()
            string1, string2 = write_csv(header,row,the_run,file_csv,filespec)
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
        from mujpy.aux.aux import int2min, int2_method_key, rebin_decay, write_csv, min2int
        
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
                
                string1, string2 = write_csv(header,row,the_run,file_csv,filespec) # DEBUG

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
        (C1) not yet
        '''
    def dofit_multirun_singlegroup_sequential(self,returntup,chain):
        '''
        performs sequential fit on many-run, single-group data
        (B1) tested
        '''
        from iminuit import Minuit
        from mujpy.aux.aux import int2min, int2_method_key, rebin, write_csv

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
                string1, string2 = write_csv(header,row,the_run,file_csv,filespec)
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
        from mujpy.aux.aux import int2min_multigroup, int2_multigroup_method_key, rebin
        from mujpy.aux.aux import stringify_groups, write_csv
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
            self.lastfit.fixed = fixedthe_fit
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
                self.summary_global(start, stop, time[1]-time[0],krun=krun)
                print('mufit dofit_sequentialrun_multigroup_userpardicts debug: fval {}, ndof {}'.format(self.lastfit.fval,self.number_dof))

                version = self.dashboard["version"]
                strgrp = stringify_groups(self.suite.groups)
                modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
                file_csv = self.suite.__csvpath__+modelname+'.'+version+'.'+strgrp+'.csv'
                the_run = self.suite._the_runs_[krun][0]
                filespec = self.suite.datafile[-3:]
                header, row = self.prepare_csv(krun=krun)
                string1, string2 = write_csv(header,row,the_run,file_csv,filespec)
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
        from mujpy.aux.aux import get_grouping, get_title, chi2std, print_components, min2int, value_error
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
            for name,value,error in zip(names,values,errors):
                f.write('| '+print_components (name, value, error)+'\n')
#                self.log('mufit summary debug: {} = {}, error = {}'.format(name, value,error))
            f.write(86*'-'+'\n')
            f.write('             Best fit performed on {}'.format(dt_string))
            self.log(78*'-')
            for name,value,error in zip(reversed(names),reversed(values),reversed(errors)):
                self.log('| '+print_components(name, value, error))
            self.log('|'+76*'-'+'|') 
            self.log('|χᵣ² = {:.3f}({:.3f},{:.3f}), on [{:.2f}ns, {:.2}µs] {:.2f}ns/bin @{}|'.format(chi,
		                                 lowchi,highchi,start,stop,dt*1000,dt_string))
            self.log('| Run {}: {}         α = {:.3f} on group: {} - {}         |'.format(nrun,
		                                 title,alpha,fgroup,bgroup))
            self.log(' '+77*'_')

    def summary_sequential(self, start, stop, dt, k):
        '''
        input: k is index in _the_runs_, default 0
        initial version: prints single fit single group result
        '''
        from mujpy.aux.aux import get_grouping, get_title, chi2std, print_components, min2int
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
        file_log = self.suite.logpath+modelname+'.'+str(nrun)+'.'+strgrp+'.'+version+'.log'
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
            for name,value,error in zip(names,values,errors): 
                f.write('| '+print_components(name, value, error)+'\n')
                self.log('| '+print_components(name, value, error))
#                self.log('summary_sequential mufit debug: {} = {}, error = {}'.format(name, value,error))
            f.write('|'+85*'_'+'|\n')
            f.write('             Best fit performed on {}'.format(dt_string))
            self.log('|'+71*'_'+'|')

    def summary_global(self,start, stop, dt, krun=0):
        '''
        input: k is index in _the_runs_, default 0
        initial version: prints single fit single group result
        '''
        from mujpy.aux.aux import get_grouping, get_title, chi2std, stringify_groups
        from mujpy.aux.aux import print_components, min2int_multigroup, value_error
        from datetime import datetime

        modelname = ''.join([component["name"] for component in self.dashboard['model_guess']])
        version = self.dashboard["version"]
        the_run = self.suite._the_runs_[krun][0]
        nrun = the_run.get_runNumber_int()
        title = get_title(the_run)
        #group = self.suite.groups[kgroup] # assumes only one group
        #fgroup, bgroup, alpha = group['forward'],\
		#				        group['backward'],\
		#				        group['alpha']
        strgrp = stringify_groups(self.suite.groups)
        chi = self.lastfit.fval /self.number_dof 
        # print('summary_global mufit debug FCN = {}, number of DOF = {}'.format(self.lastfit.fval,self.number_dof))
        lowchi, highchi = chi2std(self.number_dof)
        start, stop = self.suite.time[start]*1000, self.suite.time[stop]
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")  

        file_log = self.suite.logpath+modelname+'.'+str(nrun)+'.'+strgrp+'.'+version+'.log'
        names, values, errors = min2int_multigroup(self.dashboard,
							        self.lastfit.values,self.lastfit.errors)
        # list (groups) of lists (omponents) of lists (parameters)
        with open(file_log,'w') as f:
            f.write(' '+96*'_'+'\n')
            self.log(' '+71*'_')
            string = '| Run {}: {}    Global fit of {}'.format(nrun,title,dt_string)
            f.write(string+39*' '+'|\n')
            self.log(string+14*' '+'|')
            string = '| χᵣ² = {:.3f}({:.3f},{:.3f}) ,    on [{:.2f}ns, {:.2}µs, {:.2f}ns/bin]'.format(chi,lowchi,highchi,start,stop,dt*1000)
            f.write(string+46*' '+'|\n')
            self.log(string+21*' '+'|')
            for g1,n1,v1,e1,g2,n2,v2,e2 in zip(self.suite.groups[::2],names[::2],values[::2],errors[::2],
                                               self.suite.groups[1::2],names[1::2],values[1::2],errors[1::2]):
                fg1,bg1,al1 = g1['forward'], g1['backward'], g1['alpha'] 
                fg2,bg2,al2 = g2['forward'], g2['backward'], g2['alpha'] 
                f.write('|'+65*'-'+' on group: {} - {}   α = {:.3f}   |\n'.format(fg1,bg1,al1)) 
                self.log('|'+40*'-'+' on group: {} - {}   α = {:.3f}   |'.format(fg1,bg1,al1))
            # for k,name in enumerate(names): # replaced by:
                for nam,val,err in zip(n1,v1,e1):
                    f.write('| '+print_components (nam, val, err)+'\n')
                    self.log('| '+print_components(nam, val, err))
#                self.log('mufit summary debug: {} = {}, error = {}'.format(name, value,error))
                string = ' on group: {} - {}   α = {:.3f}   |'.format(fg2,bg2,al2)
                f.write('|'+65*'-'+string+'\n') 
                self.log('|'+40*'-'+string)
                for nam,val,err in zip(n2,v2,e2):
                    f.write('| '+print_components (nam, val, err)+'\n')
                    self.log('| '+print_components(nam, val, err))
#                self.log('mufit summary debug: {} = {}, error = {}'.format(name, value,error))
            f.write('|'+96*'_'+'|\n')
        self.log('|'+71*'_'+'|')



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
        from mujpy.aux.aux import get_title, spec_prec, chi2std, initialize_csv 
        from mujpy.aux.aux import minparam2_csv, chi2_csv, min2int, min2int_multigroup
        from mujpy.aux.aux import multigroup_in_components, userpars

        # print('k = {}, self.nrun = {}'.format(k,[j for j in self.nrun]))


        filespec = self.suite.datafile[-3:]
        the_run = self.suite._the_runs_[krun][0]
        # print('mufit prepare_csv debug, krun ={} run={}'.format(krun,the_run.get_runNumber_int()))
        chi = self.lastfit.fval /self.number_dof 
        lowchi, hichi = chi2std(self.number_dof)

        row = initialize_csv(the_run.get_field(), filespec, the_run) 

        dashcsv = False if self.suite.multi_groups() else self.dashboard["model_guess"]
        min2dash = min2int_multigroup if (self.suite.multi_groups() and userpars(self.dashboard)) else min2int
        dashboard =  self.dashboard if (self.suite.multi_groups() and userpars(self.dashboard)) else self.dashboard["model_guess"]
        
        row += minparam2_csv(dashcsv,self.lastfit.values,self.lastfit.errors)
        row += chi2_csv(chi,lowchi,hichi,self.suite.groups,self.suite.offset)
        row += '\n'
        # row is formatted with appropriate rounding, write directly
        # self.console(row)
        if filespec == 'bin' or filespec == 'mdu':
            header = ['Run','T_cryo[K]','e_T_cryo[K]','T_sample[K]','e_T_sample[K]','B[G]']
        else:
            header = ['Run','T[K]','eT[K]','B[G]']
        # now component names for header
        # print('prepare_csv mufit debug:  sum(multigroup_in_components(self.dashboard)) {}'.format(bool(sum(multigroup_in_components(self.dashboard)))))
        if sum(multigroup_in_components(self.dashboard)) and self.suite.single():
            pardicts = self.dashboard['userpardicts_guess']
            for pardict in pardicts:
                header.append(pardict['name'])
                header.append('e_'+pardict['name'])
        else:
            
            parlists, _, _ = min2dash(dashboard,
							        self.lastfit.values,self.lastfit.errors)
            parlists = parlists[0] if isinstance(parlists[0][0],list) else parlists
            for parlist in parlists:
                # print('mufit prepare_csv debug: parlist: {}'.format(parlist))
                for parname in parlist:
                    # print('mufit prepare_csv debug: parnames: {}'.format(parname))
                    header.append(parname)
                    header.append('e_'+parname)
        header.append('chi2_r')
        header.append('e_chi2_low')
        header.append('e_chi2_hi')
        if self.suite.multi_groups():
            for jgroup,group in enumerate(self.suite.groups):
                header.append('alpha{}'.format(jgroup))
        else:
            header.append('alpha')
        header.append('offset_bin')
        header.append('timestring\n')
        return ' '.join(header), row

    def save_fit(self,krun,kgroup,string):
        '''
        input:
            krun is index in self.suite._the_runs_
            kgroup is indek in self.suite.groups
        saves a dashboard file adding the bestfit parameters as "model_result"
        These saves are individual runs
        use "version" as additional label to qualify fit
        filename is logpath + modelname + nrun  + strgrp + version .json
        nrun = runNumber, strgrp = shorthand for group
        '''
        from mujpy.aux.aux import min2int
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
        filename is logpath + modelname + nrun  + srtgrp0 + strgrp...  + version .json
        nrun = runNumber, strgrp0,1,... = shorthand for allgroups
        '''
        from mujpy.aux.aux import stringify_groups
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


    def show_calib(self):
        '''
        output:
            t time 
            a asymmetry
            e asymmetry error
            f guess fit function for calib mode
        for degugging single run calibs
        '''
        from mujpy.aux.aux import int2_method_key, int2min
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
        
