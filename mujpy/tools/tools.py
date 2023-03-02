########################
# FFT AUTO PHASE METHODS
########################
def autops(data, fn, p0=0.0, p1=0.0):
    """
    Automated phase correction from NMRglue by https://github.com/jjhelmus
    These functions provide support for automatic phasing of NMR data. 


    Automatic linear phase correction

    Parameters

        data : ndarray

             Array of NMR data.

        fn : str or function

             Algorithm to use for phase scoring. Built in functions can be
             specified by one of the following strings: "acme", "peak_minima"

        p0 : float

            Initial zero order phase in degrees.

        p1 : float

            Initial first order phase in degrees.

    Returns

        ndata : ndarray

            Phased NMR data.

    """

    import numpy as np
    import scipy.optimize
    from io import StringIO # Python3 use: from io import StringIO
    from contextlib import redirect_stdout

    
    if not callable(fn):
        fn = {
            'peak_minima': _ps_peak_minima_score,
            'acme': _ps_acme_score,
        }[fn]
    
    opt = [p0, p1]
    with StringIO() as buf, redirect_stdout(buf):   
        opt = scipy.optimize.fmin(fn, x0=opt, args=(data, ))
        mystdout = buf.getvalue()
    return ps(data, p0=opt[0], p1=opt[1]), opt[0], opt[1], mystdout


def _ps_acme_score(ph, data):
    """
    Phase correction using ACME algorithm by Chen Li et al.
    Journal of Magnetic Resonance 158 (2002) 164-168

    Parameters
    * pd : tuple, current p0 and p1 values
    * data : ndarray, array of NMR data.

    Returns
    * score : float, value of the objective function (phase score)

    """
    import numpy as np

    stepsize = 1

    phc0, phc1 = ph

    s0 = ps(data, p0=phc0, p1=phc1)
    data = np.real(s0)

    # Calculation of first derivatives
    ds1 = np.abs((data[1:]-data[:-1]) / (stepsize*2))
    p1 = ds1 / np.sum(ds1)

    # Calculation of entropy
    p1[p1 == 0] = 1

    h1 = -p1 * np.log(p1)
    h1s = np.sum(h1)

    # Calculation of penalty
    pfun = 0.0
    as_ = data - np.abs(data)
    sumas = np.sum(as_)

    if sumas < 0:
        pfun = pfun + np.sum((as_/2) ** 2)

    p = 1000 * pfun

    return h1s + p


def _ps_peak_minima_score(ph, data):
    """
    Phase correction using simple minima-minimisation around highest peak
    This is a naive approach but is quick and often achieves reasonable
    results.  The optimisation is performed by finding the highest peak in the
    spectra (e.g. TMSP) and then attempting to reduce minima surrounding it.
    Parameters
    * pd : tuple, current p0 and p1 values
    * data : ndarray, array of NMR data.

    Returns
    * score : float, value of the objective function (phase score)

    """

    phc0, phc1 = ph

    s0 = ps(data, p0=phc0, p1=phc1)
    data = np.real(s0)

    i = np.argmax(data)
    mina = np.min(data[i-100:i])
    minb = np.min(data[i:i+100])

    return np.abs(mina - minb)

def ps(data, p0=0.0, p1=0.0, inv=False):
    """
    Linear phase correction

    Parameters

        data : ndarray

            Array of NMR data.

        p0 : float

            Zero order phase in degrees.

        p1 : float

            First order phase in degrees.

        inv : bool, optional

            True for inverse phase correction

    Returns

        ndata : ndarray

            Phased NMR data.

    """
    import numpy as np

    p0 = p0 * np.pi / 180.  # convert to radians
    p1 = p1 * np.pi / 180.
    size = data.shape[-1]
    apod = np.exp(1.0j * (p0 + (p1 * np.arange(size) / size))).astype(data.dtype)
    if inv:
        apod = 1 / apod
    return apod * data

##############
# MU FIT AUX
##############

def TauMu_mus():
    '''
    muon mean lifetime in microsecond
    from Particle Data Group 2017 
    (not present in scipy.constants)
    '''
    return 2.1969811 
    
def _errors_(component,available_components):
    '''
    inputs: one legal mucomponent name contained 
    in the _available_components_(), which must be the second input
    output: a list of errors (steps), one for each parameter of this component
    '''
    #print(component,available_components)
    k = [item['name'] for item in available_components].index(component)
    return [pardict["error"] for pardict in available_components[k]['pardicts']] 

def _limits_(component,available_components):
    '''
    inputs: one legal mucomponent name contained 
    in the _available_components_(), which must be the second input
    output: a list of lists of limits (low, high), one for each parameter of this component
    '''
    k = [item['name'] for item in available_components].index(component)
    return [pardict["limits"] for pardict in available_components[k]['pardicts']] 

def add_step_limits_to_model(dash_in):
    '''
    input: original dashboard dash_in, already checked 
    output: dash_out is a deepcopy including 'error' and 'limits'
    '''
    from copy import deepcopy
    from mujpy.tools.tools import _available_components_, _errors_, _limits_
    
    available_components = _available_components_()
    dash_out = deepcopy(dash_in)   
    # these lists contain all parameter values in the dashboard, including their error steps and limits

    for component in dash_out['model_guess']:
        steps = _errors_(component['name'],available_components)             
        limits = _limits_(component['name'],available_components) 
        for j,pardict in enumerate(component['pardicts']):
            pardict['error'] = steps[j]                  
            pardict['limits'] = limits[j]                         
    return dash_out
   
def _available_components_():
    '''
    returns a list of template dictionaries (one per fit component):
    retreived magically from the mucomponents mumodel class.

    Each dictionary contains 'name' and 'pardicts', 
           'pardicts' = list of parameter dictionaries, 
                        keys: 
                          'name',
                          'error,
                          'limits'
           errore are used by minuit as initial steps
           limits are 
               [None,None] for uncostrained parameters A,B,φ,λ
               [0,None] for positive parity parameters Δ,σ
                        and for positive defined parameters 'α','β','Λ','ν'
               [0,0] for fake parameter BL
    ::  ({'name':'bl','pardicts':[{'name':'A','error':0.01,'limits'[None,None]},
                                  {'name':'λ','error':0.01,'limits'[None,None]}}, 
                                  ...)
    '''
    from mujpy.mucomponents.mucomponents import mumodel
    from iminuit import describe
    
    available_components = [] # generates the template of available components.
    for name in [module for module in dir(mumodel()) if module[0]!='_']: # magical extraction of component names
        pars = describe(mumodel.__dict__[name])[2:]            #  [12:] because the first two arguments are self, x
        _pars = [] 
        # print('pars are {}'.format(pars))
        tip = eval('mumodel.'+name+'.__doc__')
        positive_defined = ['α','β','Λ','ν']
        positive_parity = ['Δ','σ']
        for parname in pars:
        # parname, error, limits
        # In this template only
        #   {'name':'amplitude','error':0.01,'limits':[0, 0]}
        # parameter name will get a label later 
            error, limits = 0.002, [None, None] # defaults for 'A', 'λ', 'Γ'
            if parname == 'B' or parname == 'Bd': error = 0.05
            if parname == 'BL': error, limits = 0, [0,0]
            if parname == 'φ': error = 1.0
            if parname in positive_defined+positive_parity: limits = [0., None]
            # add here special cases for errors and limits, e.g. positive defined parameters
            _pars.append({'name':parname,'error':error,'limits':limits})
        available_components.append({'name':name,'pardicts':_pars,'tip':tip})
    # [available_components[i]['name'] for i in range(len(available_components))] 
    # list of just mucomponents method names
    return available_components
    
#def _available_gradients_(component):
#    '''
#    returns True if the component has an analytic gradient
#    i.e. for component name xx in the mucomponents mumodel class, 
#    a method _grad_xx_ in the same class.
#    '''
#    from mujpy.mucomponents.mucomponents import mumodel
#    
#    methods_with_grad = [module[6:8] for module in dir(mumodel()) if module[0:6]=='_grad_']: # magical extraction of component names
#    return component in methods_with_grad
    
def validmodel(model):
    '''
    checks valid simple name "almlmg"
    '''
    from mujpy.tools.tools import _available_components_
    # print('validmodel: {}'.format(model))
    available_components =_available_components_() # creates list automagically from mucomponents
    component_names = [available_components[i]['name'] 
                            for i in range(len(available_components))]
    components = [model[i:i+2] for i in range(0, len(model), 2)]
    # print('valid model, available components: ',*component_names)
    if not components: # empty model
        return False
    for component in components: 
        if component in component_names:
            pass
        else:
            return False
    if 'al' in components: # check that model has only one 'al' at the beginning
        if model.count('al')>1 or model.index('al')>0:
            return False      
    return True

def get_fit_range(string):
    '''
    transform a valid string for fit_range
    into a list of integers
    '''
    fit_range = []
    for chan in string.split(','):
        fit_range.append(int(chan))
    return fit_range

def checkvalidmodel(name,component_names):
    '''
    checkvalidmodel(name) checks that name is either  
    ::      A1, B1: 2*component string of valid component names, e.g.
                        'daml' or 'mgmgbl'
                                                                  
    ::      or A2, B2: same, ending with 1 digit, number of groups (max 9 groups), 
                        'daml2' or 'mgmgml2' (2 groups)
    ::      or C1: same, beginning with 1 digit, number of external minuit parameters (max 9)
                        '3mgml' (3 external parameters e.g. A, f, phi)
    ::      or C2: same, both previous options
                        '3mgml2' (3 external parameters, 2 groups)  
    '''
    from mujpy.tools.tools import modelstrip
    
    try:
        name, nexternals = modelstrip(name)
    except:
        # self.console('name error: '+name+' contains too many externals or groups (max 9 each)')
        error_msg = 'name error: '+name+' contains too many externals or groups (max 9 each)'
        return False, error_msg # err code mess
    # decode model
    numberofda = 0
    components = [name[i:i+2] for i in range(0, len(name), 2)]
    for component in components: 
        if component == 'da':
            numberofda += 1           
        if component == 'al':
            numberofda += 1           
        if numberofda > 1:
            # self.console('name error: '+name+' contains too many da. Not added.')
            error_msg = 'name error: '+name+' contains too many da/al. Not added.'
            return False, error_msg # error code, message
        if component not in component_names:
            # self.console()
            error_msg = 'name error: '+component+' is not a known component. Not added.'
            return False, error_msg # error code, message
    return True, None

######################
# GET_TOTALS
######################
def get_totals(suite):
    '''
    calculates the grand totals and group totals 
    of a single run 
    iput is self.suite of class musuite
    returns strings totalcounts groupcounts nsbin maxbin

    '''
    import numpy as np
    # called only by self.suite after having loaded a run or a run suite

    ###################
    # grouping set 
    # suite.grouping['forward'] and suite.grouping['backward'] are np.arrays of integers
    # initialize totals
    ###################
    
    for k,d in enumerate(suite.grouping):
        if not k:
            gr = np.concatenate((suite.grouping[k]['forward'],suite.grouping[k]['backward']))
        else:
            gr = np.concatenate((gr,np.concatenate((suite.grouping[k]['forward'],suite.grouping[k]['backward']))))
    ts,gs =  [],[]

    for k,runs in enumerate(suite._the_runs_):
        tsum, gsum = 0, 0
        for j,run in enumerate(runs): # add values for runs to add
            n1 = suite.offset+suite.nt0[0]
            for counter in range(run.get_numberHisto_int()):
                if suite.datafile[-3:]=='bin' or suite.datafile[-3:]=='mdu':
                    n1 = suite.offset+suite.nt0[counter] 
                histo = np.array(run.get_histo_vector(counter,1)).sum() 
                tsum += histo
                if counter in gr:
                    gsum += histo
        ts.append(tsum)
        gs.append(gsum)
        # print('In get totals inside loop,k {}, runs {}'.format(k,runs))

    #######################
    # strings containing 
    # individual run totals
    #######################
    # self.tots_all.value = '\n'.join(map(str,np.array(ts)))
    # self.tots_group.value = '       '.join(map(str,np.array(gs)))

    # print('In get totals outside loop, ts {},gs {}'.format(ts,gs))
    #####################
    # display values for self._the_runs_[0][0] 
#        self.totalcounts.value = str(ts[0])
#        self.groupcounts.value = str(gs[0])
        # self.console('Updated Group Total for group including counters {}'.format(gr)) # debug 
#        self.nsbin.value = '{:.3}'.format(self._the_runs_[0][0].get_binWidth_ns())
#        self.maxbin.value = str(self.histoLength)
    return str(int(ts[0])), str(int(gs[0])), '{:.3}'.format(suite._the_runs_[0][0].get_binWidth_ns()), str(suite.histoLength)


def _nparam(model):
    '''
    input: dashboard['model_guess']
    output: ntot, nmintot, nfree
    '''
    number_components = len(model)
    # print('_nparam tools debug: model {}'.format(model))
    ntot = sum([len(model[k]['pardicts']) 
                                 for k in range(number_components)]) # total number of component parameters
    flag = [pardict['flag'] for component in model for pardict in component['pardicts']]
    nmintot = ntot - sum([1 for k in range(ntot) if flag[k]=='=']) # ntot minus number of functions 
    nfree = nmintot - sum([1 for k in range(ntot) if flag[k]=='!']) # ntot minus number of fixed parameters 
    return ntot, nmintot, nfree
    
##################################################################
# int2min methods: generate guess values, errors and limits
#                  of minuit parameters
#  int2min : 
#  int2min_multigroup : assumes all parameters are in userpardicts
#  int2min_multirun : assumes both userpardicts and active parameters
#                    must generate daughters for local parameters
##################################################################

def int2min(model):
    '''
    input: 
        model 
            either dashboard["model_guess"] (after add_step_limits_to_model)
            or  dashboard["model_guess"] both lists of dicts
    output: a list of lists:  
        values: minuit parameter values, either guess of result
        errors: their steps
        fixed: True/False for each
        limits: [low, high] limits for each or [None,None]  
        names: name of parameter 'x_label' for each parameter
        pospar: parameter for which component is positive parity, eg s in e^{-(s*t)^2/2}
    '''
    from mujpy.tools.tools import _nparam

    dum, ntot, dum  = _nparam(model)
    
    #####################################################
    # the following variables contain the same as input #
    # parameters to iMinuit, removing '='s (functions)  #
    #####################################################
    
    positive_parity = ['Δ','σ']                                                    
    val, err, fix, lim = [], [], [], []           
    names = []
    pospar = [] # contains index of positive parity parameters, to rerun with no limits

    for component in model:  # scan the model components
        label = component['label']
        for k,pardict in enumerate(component['pardicts']):  # list of dictionaries
            if pardict['flag'] != '=': #  skip functions, only new minuit parameters
                if pardict["name"] in positive_parity: pospar.append(k)
                if pardict['flag'] == '~':
                    fix.append(False)
                elif pardict['flag'] == '!':
                    fix.append(True)
                val.append(float(pardict['value']))
                names.append(pardict['name']+'_'+label) 
                err.append(float(pardict['error']))
                lim.append(pardict['limits'])
                # print('tools int2min debug: pardict name {} limits {}'.format(names[-1],lim[-1]))
    # self.console('val = {}\nerr = {}\nfix = {}\nlim = {},\npar name = {} '.format(val,err,fix,lim, names)) 

    return val, err, fix, lim, names, pospar

def int2min_multigroup(pardicts):
    '''
    input: 
        pardicts 
            either dashboard["userpardicts_guess"] if guess = True
            or  dashboard["userpardicts_result"] if guess = False
    output: a list of lists:  
        values: minuit parameter values, either guess of result
        errors: their steps
        fixed: True/False for each
        limits: [low, high] limits for each or [None,None]  
        name: name of parameter 'x_label' for each parameter
        pospar: parameter for which component is positive parity, eg s in e^{-(s*t)^2/2}
    this works for A2 single fit, multigroup with userpardicts parameters = Minuit parameters
    '''
    
    #####################################################
    # the following variables contain the same as input #
    # parameters to iMinuit, removing '='s (functions)  #
    #####################################################
                                                        
    val, err, fix, lim = [], [], [], []           
    name = []
    pospar = [] # contains index of positive parity parameters, to rerun with no limits

    for k,pardict in enumerate(pardicts):  # scan the model components
        if 'positive_parity' in pardict.keys(): pospar.append(k)
        errstd = 'error' if 'error' in pardict.keys() else 'std'
        val.append(float(pardict['value']))
        name.append(pardict['name']) 
        err.append(float(pardict[errstd]))
        if 'error' in pardict.keys():
            lim.append(pardict['limits'])
            # print('tools debug: par {} limits {}'.format(pardict['name'],pardict['limits']))
        if 'flag' in pardict.keys():
            if pardict['flag'] == '!':
                fix.append(True)
            elif pardict['flag'] == '~':
                fix.append(False)
            else:
                return False,_,_,_,_,_,_
        # self.console('val = {}\nerr = {}\nfix = {}\nlim = {}\python list with more repeated valuesncomp name = {},\npar name = {} '.format(val,err,fix,lim,name)) 
    return val, err, fix, lim, name, pospar

def int2min_multirun(dashboard,runs):
    '''
    input: 
        dashboar for single group multirun
            containins both
               dashboard["userpardicts_guess"] or dashboard["userpardicts_result"]
               dashboard["model_guess"] (after add_step_limits_to_model) or dashboard["model_result"]
        runs = list of run numbers
    output: a list of lists, each list containing minuit internal parameters for a run 
        values: minuit parameter values, either guess of result
        errors: their steps
        fixed: True/False for each
        limits: [low, high] limits for each or [None,None]  
        names: name of parameter 'x_label' for each parameter
        pospar: parameter for which component is positive parity, eg s in e^{-(s*t)^2/2}
    Note: this method knows the model, hence it can generate only minuit parameters
          mucomponents, i.e. the cost function, needs to feed also non minuit parameters to the components
          e.g. mg with A global B1 local phi global s1 global will have only two minuit parameters per run
          , the second and the fourth, and it must know which value to use for the first and third
          This is accomplished by _add_multirun_singlegroup_ in mucomponents, using self._components_
    '''
    #####################################################
    # the following variables contain the same as input #
    # parameters to iMinuit, removing '='s (functions)  #
    #####################################################
                                                        
    positive_parity = ['Δ','σ']                                                    
    val, err, fix, lim = [], [], [], []           
    name = []
    user_local = []
    pospar = [] # contains index of positive parity parameters, to rerun with no limits
    pospar_loc = []
    nlocals = 0
    # first scan the global and local user parameters 
    
    pardicts = dashboard["userpardicts_result"] if "userpardicts_result" in dashboard.keys() else dashboard["userpardicts_guess"]    
    # REMIND: insert plot guess option for _result dashboard
    model = dashboard["model_result"] if "model_result" in dashboard.keys() else dashboard["model_guess"]
    for k,pardict in enumerate(pardicts):  # scan the model components
        if 'positive_parity' in pardict.keys(): 
            if pardict["positive_parity"]:
                pospar_loc.append(k)
                pardict['limits'][0] = 0. 
            # print('debug tools int2min_multirun: pospar {} lim({}) = {}'.format(pardict["name"],k, pardict['limits']))
        errstd = 'error' if 'error' in pardict.keys() else 'std'
        if pardict["local"] or type(pardict["value"])==list: # the local key is set to False by default
            nlocals += 1
            user_local.append(pardict) # set aside this parameter for the loop over runs
            pardict['local'] = True
        else: # the first n parameters are the global user parameters
            val.append(pardict['value'])
            name.append(pardict['name'])
            err.append(pardict[errstd])
            if 'error' in pardict.keys():
                lim.append(pardict['limits'])
                # print('tools debug: par {} limits {}'.format(pardict['name'],pardict['limits']))
            if 'flag' in pardict.keys():
                if pardict['flag'] == '!':
                    fix.append(True)
                elif pardict['flag'] == '~':
                    fix.append(False)
                else:
                    return False,_,_,_,_,_,_     
                    

    kloc = k-nlocals
    # print('debug tools int2min_multirun: pospar_local {}\nuser_local = {}'.format(pospar_loc,user_local))
#    print('debug tools int2min_multirun: kloc = {}, len(fix) is {}'.format(kloc,len(fix)))
#    print('first the userpars\nval = {}\nerr = {}\nfix = {}\nlim = {}\ncomp name = {},\npar name = {} '.format(val,err,fix,lim,name)) 
        
    # now scan the runs andcreate as many replicas of the local paratmeters
    for krun,run in enumerate(runs): # run[0] is a string with the run number
        # "value" may be a single guess value for all or a list of guess values, one per run, checked at start
        for kusr,usr in enumerate(user_local): # first the local user parameter names 
            kloc += 1
            fix.append(False) # can only be not-fixed 
            if kusr in pospar_loc: 
                pospar.append(kloc) # this parameter is run version of a positive parity user local par
            if type(usr["value"])==list: 
                # print('list = {}, krun = {}'.format(usr["value"],krun))
                val.append(usr["value"][krun])
            else:
                val.append(usr["value"])
            name.append(usr["name"]+'_'+run[0])
            errstd = 'error' if 'error' in usr.keys() else 'std'
            if type(usr[errstd])==list: 
                err.append(usr[errstd][krun]) 
            else: 
                err.append(usr[errstd])
            lim.append(usr['limits'])
            
        for component in model:  # then scan the model components and add only non "="-flag parameters
            label = component['label']
            for k,pardict in enumerate(component['pardicts']):  # list of dictionaries
                if pardict['flag'] != '=': # minuit parameter
                    kloc += 1
                    if pardict["name"][0] in positive_parity: 
                        pospar.append(kloc)
                        pardict['limits'][0] = 0.
                        # print('debug tools int2min_multirun: pospar {} lim({}) = {}'.format(pardict["name"],k, pardict['limits']))
                    if pardict['flag'] == '~':
                        fix.append(False)
                    elif pardict['flag'] == '!':
                        fix.append(True)
#                    else:
#                        print('debug tools int2min_multirun: kloc = {}, pardict["flag"] is {}'.format(kloc,pardict['flag']))
                    if type(pardict["value"])==list: 
                        #print('val = {}, krun = {}'.format(pardict["value"],krun))
                        val.append(pardict["value"][krun]) 
                    else: 
                        val.append(pardict["value"])
                    name.append(pardict['name']+'_'+label+'_'+run[0]) 
                    errstd = 'error' if 'error' in pardict.keys() else 'std'
                    if type(pardict[errstd])==list: 
                        err.append(pardict[errstd][krun])
                    else: 
                        err.append(pardict[errstd])
                    lim.append(pardict['limits'])
                    pre = 0
                    for k in pospar_loc:
                        if k not in pospar: 
                            pospar.insert(pre,k)
                            pre += 1
#    print('debug tools int2min_multirun: kloc = {}, len(fix) is {}'.format(kloc,len(fix)))
    return val, err, fix, lim, name, pospar # all simple lists of sequential parameters, minuit order 

def int2fft(model):
    '''
    input: 
        model 
            dashboard["model_guess"] 
    output: 
        fft_subtract: a list of boolean values, one per model component
            fft flag True, component subtracted in residues 
    '''
    from mujpy.tools.tools import _nparam
    fft_flag = []
    fft_name = []
    for componentdict in model:  # scan the model components
        if "fft" not in componentdict.keys():
            append(False)
        else:
            append(componentdict["fft"])
        fft_name.append(componentdict["name"])
    return fft_name, fft_flag
    
##################################
# method and key methods: provide component methods 
#                           and parameter key for eval(key) in _add_
#   int2_method_key :                single run single group 
#   int2_multigroup_method_key :     single run multi group
#   int2_multirun_user_method_key :  multirun single group user
#   int2_multirun_grad_method_key :  same with grad

def int2_method_key(dashboard,the_model):
    '''
    input: 
       dashboard, the dashboard dict structure 
       the_model,  a fit model instance (not necessarily loaded)
    output: 
       a list of lists, the inner lists contain each
         method,  a mumodel component method, in the order of the model components
                   for the use of mumodel._add_.
         keys,   a list of as many lambda functions as the parameters of teh component
                 hard coding the translated "function" string for fast evaluation.
    This function applies tools.translate to the parameter numbers in formulas:
    dashboard "function" is written in terms of the internal parameter index,
    while Minuit parameter index skips shared or formula-determined ('=') parameters  
    '''
    from mujpy.tools.tools import translate, set_key

    model_guess = dashboard['model_guess']  # guess surely exists

    ntot = sum([len(model_guess[k]['pardicts']) for k in range(len(model_guess))])
    lmin = [] # initialize the minuit parameter index of dashboard function indices 
    nint = -1 # initialize the number of internal parameters
    nmin = -1 # initialize the number of minuit parameters
    method_key = []
    function = [pardict['function'] for component in model_guess for pardict in component['pardicts']]
    for k in range(len(model_guess)):  # scan the model
        name = model_guess[k]['name']
        # print('name = {}, model = {}'.format(name,self._the_model_))
        is_al_da = name=='al' or name=='da'
        bndmthd = [] if is_al_da else the_model.__getattribute__(name) 
                            # this is the method to calculate a component, to set alpha, dalpha apart
        keys = []
        # isminuit = [] not used
        flag = [item['flag'] for item in model_guess[k]['pardicts']]
        for j,pardict in enumerate(model_guess[k]['pardicts']): 
            nint += 1  # internal parameter incremente always   
            if flag[j] == '=': #  function is written in terms of nint
                # nint must be translated into nmin 
                string = translate(nint,lmin,pardict['function']) # here is where lmin is used
                # translate substitutes lmin[n] where n is the index read in the function (e.g. p[3])
                key_as_lambda = set_key(string) # NEW! calculates simple functions and speedup
                keys.append(key_as_lambda) # the function key in keys will be evaluated, key(p), inside mucomponents
                # isminuit.append(False)
                lmin.append(0)
            else:# flag[j] == '~' or flag[j] == '!'
                nmin += 1
                key_as_lambda = set_key('p['+str(nmin)+']') # NEW! calculates simple functions and speedup
                keys.append(key_as_lambda) # the function key in keys will be evaluated, key(p), inside mucomponents
                lmin.append(nmin) # 
                # isminuit.append(True)
        # print('int2_method tools debug: bndmthd = {}, keys = {}'.format(bndmthd,keys))
        method_key.append([bndmthd,keys]) 
    return method_key

def int2_multigroup_method_key(dashboard,the_model):
    '''
    input: 
        dashboard, the dashboard dict structure
        fit._the_model_ is an instance of mumodel 
            (the number of groups is obtained from dashboard)
    output: a list of methods and keys, in the order of the model components 
            for the use of mumodel._add_multigroup_.
            method is a 2d vector function 
            accepting time and a variable number of lists of (component) parameters
                e.g if one component is mumodel.bl(x,A,λ)
                the corresponding component for a two group fit 
                accepts the following argument (t,[A1, A2],[λ1,λ2])               
            keys is a list of lists of strings 
            they are evaluated to produce the method parameters, 
            there are 
                ngroups strings per parameter (inner list)
                npar parametes per component (outer list)
    Invoked by the iMinuit initializing call
             self._the_model_._load_data_multigroup_
    just before submitting migrad, 

    This function does not need userpardicts and tools.translate 
    since the correspondence with Minuit parameters
    is given directly either by "function" or by "function_multi"
    '''
    from mujpy.tools.tools import multigroup_in_components, cstack, setkey
    
    model = dashboard['model_guess']  # guess surely exists
    # these are the only Minuit parameters [p[k] for k in range (nuser)]
    ntot = sum([len(model[k]['pardicts']) for k in range(len(model))])
    method_key = []
    pardicts = [pardict for component in model for pardict in component['pardicts']]
    mask_function_multi = multigroup_in_components(dashboard)
    # print('int2_multigroup_method_key tools debug: index function_multi {}\npardicts = {}'.format(mask_function_multi,pardicts))
#    if "userpardicts_guess" in dashboard.keys():
#        updicts =  dashboard["userpardicts_guess"]
#        print('int2_multigroup_method_key tools debug: userpardicts ')
#        for j,pd in enumerate(updicts):
#            print('{} {} = {}({}), {}, {} '.format(j,pd["name"],pd["value"],
#                                    pd["error"], pd["flag"],pd["limits"]))
    if sum(mask_function_multi):
        ngroups = len(pardicts[mask_function_multi.index(1)]["function_multi"])
    else:
        return []
    nint = -1 # initialize the index of the dashboard component parameters
    # p = [1.13,1.05,0.25,0.3,0.8,700,35,125,3.3,680,0.1] # debug delete
    # print('tools int2_multigroup_method_key debug: fake values k, p {}'.format([[k,par] for k,par in enumerate(p)]))
    bndmthd = {} # to avoid same name
    for j,component in enumerate(model):  # scan the model components
        name = component['name']
        keys = []
        bndmthd[name] = lambda x,*pars, name=name : cstack(the_model.__getattribute__(name),x,*pars)
        bndmthd[name].__doc__ = '"""'+name+'"""'
                            # this is the method to calculate a component, to set alpha, dalpha apart
        #print('\n\tools int2_multigroup_method_key debug: {}-th component name = {}'.format(j,bndmthd[name].__doc__))
        nint0 = nint
        for l in range(ngroups):
            key = []  
            nint = nint0
            for pardict in component['pardicts']: 
                nint += 1  # internal parameter index incremented always 
                if mask_function_multi[nint]>0:
#                    print('tools int2_multigroup_method_key debug: {}[{}] = {}'.format(pardict["name"],l,pardict["function_multi"][l])) 
                    key_as_lambda = set_key(pardict["function_multi"][l]) # NEW! calculates simple functions and speedup
                else:                
#                    print('tools int2_multigroup_method_key debug: {}[{}] = {}'.format(pardict["name"],l,pardict["function"])) 
                    key_as_lambda = set_key(pardict["function"][l]) # NEW! calculates simple functions and speedup
                # print('tools int2_multigroup_method_key debug: key_as_lambda(p) = {} **delete also p!'.format(key_as_lambda(p)))
                key.append(key_as_lambda) # the function key will be evaluated, key(p), inside mucomponents
            keys.append(key)
        #print('int2_method tools debug: appending {}-th bndmthd {} with {} groups x {} keys'.format(j,bndmthd[name].__doc__,len(keys),len(keys[0])))
        method_key.append([bndmthd[name],keys]) # vectorialized method, with keys 
        # keys = [[strp0g0, strp1g0,...],[strp0g1, strp1g1, ..],[strp0g2, strp1g2,...]..]
        # pars = [[p0g0, p1g0, ...],[p0g1, p1g1, ..],[p0g2, p1g2,...]..]
    return method_key

def int2_multirun_user_method_key(dashboard,the_model,nruns):
    '''
    input: 
        dashboard, the dashboard dict structure
        the_model is fit._the_model_ i.e. an instance of mumodel 
        nruns is the numer of runs in the suite
    output: a list of methods and a list of lists of keys, [[key,...,key],...,[key,..,key]]
    the internal list is same parameter, different runs
    the model components 
            for the use of mumodel._add_multirun_.
            method is a component function 
            accepting time and a list of parameters, e.g mumodel.bl(x,A,λ)
            key is string defining a lambda function that produces one method parameter for a specific run, 
            the list is for the same parameter over diffenet runs 
            keys is a list of lists for all the parameters (any flag) of the component
            the list of [binding,keys] is over the components of the model
    This list of [binding, keys] allows mumodel _add_multirun_ to use the minuit p list
    (n_globals global user values, followed by nruns replica of 
     n_locals local user values and a model specific number of local (~,!) component par values)
    to produce component-driven vectorized values, as many values in the vector as the runs
    In this way minuit fcn is a vector, one fcn per run,
    likewise asymm, asyme are vectors (see suite for multirun)
    and mumodel._chisquare_ cost function sums over individual runs for a unique global chisquare
    Invoked by the iMinuit initializing call
             self._the_model_._load_data_multirun_user_
    just before submitting migrad
    '''
    from mujpy.tools.tools import multigroup_in_components, cstack, translate_multirun, set_key
    from mujpy.tools.tools import get_functions_in
#            self._components_ is a list [[method,[key,...,key]],...,[method,[key,...,key]]], 
#                produced by int2_multirun_user_method_key() from mujpy.tools.tools
#                where method is an instantiation of a component, e.g. self.ml 
#                and value = eval(key) produces the parameter value
    model = dashboard['model_guess']  # guess surely exists, it is a list of component dicts, e.g. for mgbl 2 dicts
    method_key = []
    bndmthd = {} # to avoid same name
    n_locals =  [pardict["local"] for pardict in dashboard["userpardicts_guess"]].count(True)
    n_globals = len(dashboard["userpardicts_guess"])-n_locals
    kloc = n_globals+n_locals
    functions_in = get_functions_in(model,kloc-1)
    functions_out = translate_multirun(functions_in,n_locals,kloc,nruns)   
    

    # print('\n\ndebug tools int2_multirun_user_method_key functions_out = {}'.format(functions_out))
    for j,component in enumerate(model):  # scan the model components (as for the first run)
        name = component['name']
        keys = []
        # this method uses pars, a list of lists (runs) of parameter for this component, obtained by key(p) from minuit p
        bndmthd[name] = lambda x,*pars, name=name : cstack(the_model.__getattribute__(name),x,*pars)
        bndmthd[name].__doc__ = '"""'+name+'"""'
                            # no alpha in global multirun!
        # its pars are generated as a list of lists of the key_as_lambda functions
        for funcs in functions_out[j]: # funcs is a run, in the suite of runs
            key = []
            for func in funcs: # this is a parameter for this run, in the component parameters 
                #print('debug tools int2_multirun_user_method_key func = {}'.format(func))
                key_as_lambda = set_key(func) # NEW! calculates simple functions and speedup
                # function key will be evaluated as key(p) inside mucomponents
                key.append(key_as_lambda) # collect parameter key(s) of the component  
            keys.append(key) # create outer list adding component parameters for this run
        method_key.append([bndmthd[name],keys]) # vectorialized method, with its keys list of lists
        # appended to a list of [method,
        # print('debug tools int2_multirun_user_method_key: locals =\n{}'.format(globals()))
    return method_key

def int2_multirun_grad_method_key(dashboard,the_model,nruns):
    '''
    input: 
        dashboard, the dashboard dict structure
        the_model is fit._the_model_ i.e. an instance of global multirun mumodel 
        nruns is the numer of runs in the suite
    output
        minuit_ordered_grad_list 
         i.e. a list of lists [k,n,j,grad_bndmthd,dkey])], one for each minuit internal parameter p[m] 
                k n j are indices of run, component and parameter for which 
                dkndj_bndmthd is the <module> that computes the derivative of component n in run k with respect to parameter j 
                    (if par are the parameters for component n in run k, then dkndj_bndmthd(x,*par) calculates the derivative) 
                djdm is a <module> that computes the derivative of the user funct (e.g. "p[0]*p[21]') with respect to p[m]
                    (then djdm(p) is the value of the derivative)
        The products of these derivatives is sparse, i.e. non zero only for few values k,n,j
        The present method identifies each and every set of indices (k,n,j) for which the product
              gg_nj = dkndj_bndmthd(x,*par)*djdm(p) != 0       
        --------------------------- Usage
        To generate chisquare grad values for all minuit parameters, in order to optimize numpy miniut calculations;
        mucomponents _add_multirun_grad_ uses the 2D array      
                      gg = sum_n,j gg_nj 
        and multiplies it by the 2D array dcdf = 2(f-y)/e^2
        The m-th component of the chisquare gradient is sum(dcdf*gg,axis=None)
        -------------------------- General equation
        Assuming asymmetry data end errors y(k;i), ey(k;i), the expression for the chisquare gradient is
               sum_i,k {2(sum_n y_n(k;i,*par(k,n))-y(k;i))/ey(k;i)^2} * sum_n,j {partial y_n(k;i,*par(k,n)/partial par[k,n,j]} * {partial par[k,n,j]/partial p[m]}
          hereafter            dcdf                                   * sum_n,j           dkndj                                *               djdm        
    '''
    from mujpy.tools.tools import get_functions_in, diffunc, get_indices, get_number_minuit_internal
    from mujpy.tools.tools import multigroup_in_components, translate_multirun, set_key
    # firts generate dmethod_keys
    # dmethod_keys contains [[m_d,keys]...[m_d,keys]], as many as the model components (e.g. 2 for mgbl)
    # m_d is [method] if no derivative is required (e.g. bl) 
    # or [method, derivative_method] if derivative is required (e.g. mg)  
    # keys is [runkeys,...,runkeys] such that
    # par = [key(p) for key in runkeys] and method(x,*par) and derivative_method(self._x_,*par) produce the additive component for that run
    model = dashboard['model_guess']
    names = [component['name'] for component in model] 
    n_locals =  [pardict["local"] for pardict in dashboard["userpardicts_guess"]].count(True)
    n_globals = len(dashboard["userpardicts_guess"])-n_locals
    kloc = n_globals+n_locals

    functions_in = get_functions_in(model,kloc-1) # functions_in are the user func strings of the single-run model 
    functions_out = translate_multirun(functions_in,n_locals,kloc,nruns)  # user func strings of the multirun model
    minuit_ordered_grad_list = [[] for x in range(get_number_minuit_internal(nruns,n_globals,n_locals,model))] # this is the empty output container
    #print('debug tools int2_multirun_grad_method_key, minuit_grad_list = {}'.format(minuit_ordered_grad_list))
    for n_component, (component, component_name) in enumerate(zip(functions_out,names)): # the order is model components, runs, component parameter
        for k_run, run_component in enumerate(component):
            for j_parameter, func in enumerate(run_component):
                dfuncs, indices = diffunc(func)  # es func = 'p[0]*p[7]', dfuncs = ['p[7]','p[0]'] indices = [0,7]
                for dfunc,m_minuit_parameter, in zip(dfuncs,indices):               
                    grad_bndmthd = lambda x, *par, gname='_grad_'+component_name+'_'+str(j_parameter)+'_' : the_model.__getattribute__(gname)(x,*par)
                    grad_bndmthd.__doc__ = '"""'+'_grad_'+component_name+'_'+str(j_parameter)+'_"""'
                    # if e.g. component_name is 'bl'  methods must exist called _grad_bl_0_, _grad_bl_1_ ...  
                    # print('debug tools int2_multirun_grad_method_key, m_minuit_parameter = {}, k, n, j = {};{},{}'.format(m_minuit_parameter,k_run,n_component,j_parameter))
                    grad_list = minuit_ordered_grad_list[m_minuit_parameter]
                    grad_list.append([k_run,n_component,j_parameter,grad_bndmthd,set_key(dfunc)]) 
    return minuit_ordered_grad_list

def get_number_minuit_internal(nruns,n_globals,n_locals,model):
    k_mint = 0
    for j,component in enumerate(model):  # scan the model components (as for the first run)
        flags = [pardict["flag"] for pardict in component["pardicts"]] # these are the flags in the present component
        for k,flag in enumerate(flags): # as many flags as parameters in component
            if flag!="=":
                k_mint += 1
    return n_globals + nruns*(n_locals + k_mint)
    
def get_functions_in(model,kk):
    '''
    input 
        model = single-run model dashboard dict
        kk = kloc -1, is incremented at each free parameter of the model, so that it scans the internal minuit indices
             for these parameters (ignoring those determined by a user funct e.g. "p[0]*p[4]"
    output 
        functions_in = list of lists, one per component, of user functs, one per parameter, for the single-run model, 
                       all component parameters,  including "~" and "!", are translated to appropriate user funct
    '''
    functions_in = []
    for j,component in enumerate(model):  # scan the model components (as for the first run)
        flags = [pardict["flag"] for pardict in component["pardicts"]] # these are the flags in the present component
        function_in = [pardict["function"] for pardict in component["pardicts"]] # these are the original function (some are empty)
        for k,flag in enumerate(flags): # as many flags as parameters in component
            if flag!="=": # this parameter is among the minuit parameters
                kk += 1 # this is the minuit index of the current first run parameter
                # suppose mg with n_globals = 5 (0,1,2,3,4), n_locals = 1 (5)
                #   6 A = f(p[1],p[5]) 7 B = p[6] 8 φ = p[2] 9 σ = p[7] 
                # k     kk          n_equals
                # 0     -              1
                # 1   5+1+1-1=6        -
                # 2     -              2
                # 3   5+1+3-2=7        -
                function_in[k] = 'p['+str(kk)+']' # write a fake "function" to eval this parameter as 'p[kk]'
        functions_in.append(function_in)
    return functions_in
    
def set_key(string):   
    """
    input: the function string from the json or the mudash dashboard
         e.g. that written in the json file as 'function':'p[0]*(0.5+1/pi*arctan(p[2])'
              or typed into mudash text widget as 'p[0]*(0.5+1/pi*arctan(p[2])'
    output: key, a python method, such that in mumodel mucomponents _add_ the command 
            key(p) evaluates the formula 
            the evaluation knows simple numpy functions, see the import in string code, below
    """
    code = """
from numpy import cos, sin, tan, sinh, cosh, tanh, log, pi, exp, sqrt, real, abs, arctan
def foo():
"""
    string = '"lambda p: '+string+ '"' 
    string = "    key = eval('"+string+"')"
    # print('string ={}'.format(string))
    code = code + string + """
    return key
"""
    # print('code = {}'.format(code))
    exec(code,globals(),globals())
    key = eval(foo())
    return key   
        
#def fstack(npfunc,x,*pars):
#    '''
#    vectorialize npfunc
#    input: 
#        npfunc numpy function with input (x,*argv)
#        x time
#        *pars is a list of lists of parameters, 
#              list len is the output_function_array.shape[0]
#    output:
#        output_function_array
#            stacks vertically n replica of npfunc distributing parameters as in
#            (x, *argv[i]) for each i-th replica 
#    '''
#    # fstack reproduces the parameter input of a component according to         
#    # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
#    # where the outer list a replica of the same component method 
#    # either over several groups (multigroup) or over several runs (multirun)
#    # as of now this method does not work for the multirun multigroup userpar case (C2)

#    from numpy import vstack
#    for k,par in enumerate(pars):
#        if k:
#            # print('tools fstack debug: npfunc.__doc__: {}'.format(npfunc.__doc__))
#            f = vstack((f,npfunc(x,*par)))
#        else:
#            # print('tools fstack debug: k=0 npfunc.__doc__: {}'.format(npfunc.__doc__))
#            f = npfunc(x,*par)
#    return f
    
def cstack(npfunc,x,*pars):
    '''
    vectorialize npfunc
    input: 
        npfunc numpy function with input (x,*argv)
        x time
        *pars is a list of lists of parameters, 
              list len n is the output_function_array.shape[0]
    output:
        output_function_array
            stacks vertically n replica of npfunc distributing parameters as in
            (x, *argv[i]) for each i-th replica 
    '''
    # cstack reproduces the parameter input of a component according to         
    # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
    # where the outer list a replica of the same component method 
    # either over several groups (multigroup) or over several runs (multirun)
    # as of now this method does not work for the multirun multigroup userpar case (C2)

    from numpy import concatenate    
    return concatenate([npfunc(x,*par) for par in pars]).reshape(-1,x.shape[0])
    
def int2_calib_method_key(dashboard,the_model):
    '''
    NOT USED, remove
    input: the dashboard dict structure and the fit model 'alxx..' instance
           the actual model contains 'al' plus 'xx', ..
           the present method considers only the latter FOR PLOTTING ONLY 
           (USE int2_method for the actual calib fit)
    output: a list of methods for calib fits, in the order of the 'xx..' model components 
            (skipping al) for the use of mumodel._add_single_.
    Invoked by the iMinuit initializing call
             self._the_model_._load_data_, 
    just before submitting migrad, 
    self._the_model_ is an instance of mumodel 
     
    This function applies tools.translate to the parameter numbers in formulas
    since on the dash each parameter of each component gets an internal number,
    but alpha is popped and shared or formula-determined ('=') ones are not minuit parameters  
    '''
    from mujpy.tools.tools import translate

    model_guess = dashboard['model_guess']  # guess surely exists

    ntot = sum([len(model_guess[k]['pardicts']) for k in range(len(model_guess))])-1 # minus alpha
    lmin = [] # initialize the minuit parameter index of dashboard function indices 
    nint = -1 # initialize the number of internal parameters
    nmin = -1 # initialize the number of minuit parameters
    method_key = []
    function = [pardict['function'] for component in model_guess for pardict in component['pardicts']]
    for k in range(1,len(model_guess)):  # scan the model popping 'al' and its parameter 'alpha'
        name = model_guess[k]['name']
        # print('name = {}, model = {}'.format(name,self._the_model_))
        bndmthd = the_model.__getattribute__(name) 
        keys = []
        # isminuit = [] not used
        flag = [item['flag'] for item in model_guess[k]['pardicts']]
        for j,pardict in enumerate(model_guess[k]['pardicts']): 
            nint += 1  # internal parameter incremente always   
            if flag[j] == '=': #  function is written in terms of nint
                # nint must be translated into nmin 
                string = translate(nint,lmin,pardict['function']) # here is where lmin is used
                # translate substitutes lmin[n] where n is the index read in the function (e.g. p[3])
                keys.append(string) # the function will be eval-uated, eval(key) inside mucomponents
                # isminuit.append(False)
                lmin.append(0)
            else:# flag[j] == '~' or flag[j] == '!'
                nmin += 1
                keys.append('p['+str(nmin)+']')  # this also needs direct translation                      
                lmin.append(nmin) # 
                # isminuit.append(True)
        method_key.append([bndmthd,keys]) 
    return method_key

def int2_calib_multigroup_method_key(dashboard,the_model):
    '''
    NOT USED, remove
    input: the dashboard dict structure and the fit model 'alxx..' instance
           the actual model contains 'al' plus 'xx', ..
           the present method considers only the latter FOR PLOTTING ONLY 
           (USE int2_method for the actual calib fit)
    output: a list of methods for calib fits, in the order of the 'xx..' model components 
            (skipping al) for the use of mumodel._add_single_.
    Invoked by the iMinuit initializing call
             self._the_model_._load_data_, 
    just before submitting migrad, 
    self._the_model_ is an instance of mumodel 
     
    This function applies tools.translate to the parameter numbers in formulas
    since on the dash each parameter of each component gets an internal number,
    but alpha is popped and shared or formula-determined ('=') ones are not minuit parameters  
    '''
    from mujpy.tools.tools import translate
    print('int2_calib_multigroup_method_key tools debug: copy of non multigroup, adapt!')
    model_guess = dashboard['model_guess']  # guess surely exists

    ntot = sum([len(model_guess[k]['pardicts']) for k in range(len(model_guess))])-1 # minus alpha
    lmin = [] # initialize the minuit parameter index of dashboard function indices 
    nint = -1 # initialize the number of internal parameters
    nmin = -1 # initialize the number of minuit parameters
    method_key = []
    function = [pardict['function'] for component in model_guess for pardict in component['pardicts']]
    for k in range(1,len(model_guess)):  # scan the model popping 'al' and its parameter 'alpha'
        name = model_guess[k]['name']
        # print('name = {}, model = {}'.format(name,self._the_model_))
        bndmthd = the_model.__getattribute__(name) 
        keys = []
        # isminuit = [] not used
        flag = [item['flag'] for item in model_guess[k]['pardicts']]
        for j,pardict in enumerate(model_guess[k]['pardicts']): 
            nint += 1  # internal parameter incremente always   
            if flag[j] == '=': #  function is written in terms of nint
                # nint must be translated into nmin 
                string = translate(nint,lmin,pardict['function']) # here is where lmin is used
                # translate substitutes lmin[n] where n is the index read in the function (e.g. p[3])
                keys.append(string) # the function will be eval-uated, eval(key) inside mucomponents
                # isminuit.append(False)
                lmin.append(0)
            else:# flag[j] == '~' or flag[j] == '!'
                nmin += 1
                keys.append('p['+str(nmin)+']')  # this also needs direct translation                      
                lmin.append(nmin) # 
                # isminuit.append(True)
        method_key.append([bndmthd,keys]) 
    return method_key

def min2int(model_guess,values_in,errors_in):
    '''
    input:
        model_component from dashboard
        values_in Minuit.values
        errors_in Minuit.errors
    output: for all dashbord parameters
        names list of lists of parameter names
        values_out list of lists of their values
        errors_out list of lists of their errors
    reconstruct dashboard with Minuit best fit values and errors
    for print_components, compact fit summary 
    '''
    # 
    # initialize
    #
    from mujpy.tools.tools import translate

    names, values_out, p, errors_out, e = [], [], [], [], []
    nint = -1 # initialize
    nmin = -1
    lmin = []
    flag = [pardict['flag'] for component in model_guess for pardict in component['pardicts']]
#    flag = [pardict['flag'] for component in model_guess for pardict in component['pardicts']]
    for k,component in enumerate(model_guess):  # scan the model
        component_name = component['name']
        name, value, error = [], [], []
        label = model_guess[k]['label']
        
        for j,pardict in enumerate(model_guess[k]['pardicts']): # list of dictionaries, par is a dictionary
            nint += 1  # internal parameter incremented always
            if j==0:
                name.append('{}{}_{}'.format(component_name,pardict['name'],label))
            else:
                name.append('{}_{}'.format(pardict['name'],label))
            if flag[nint] != '=': #  skip functions, they are not new minuit parameter
                nmin += 1
                lmin.append(nmin)
                p.append(values_in[nmin]) # needed also by functions
                value.append(values_in[nmin])
                e.append(errors_in[nmin])
                error.append(errors_in[nmin]) # parvalue item is a string
            else: # functions, calculate as such
                # nint must be translated into nmin 
                string = translate(nint,lmin,pardict['function'])  
                p.append(eval(string))
                value.append(eval(string))
                e.append(eval(string.replace('p','e')))
                error.append(eval(string.replace('p','e')))
                lmin.append(0) # not needed
        names.append(name)
        values_out.append(value)
        errors_out.append(error)
    return names, values_out, errors_out # list of parameter values 

def min2int_multirun(dashboard,p,e,_the_runs_):
    '''
    input:
        dashboard;  userpardicts_guess and model_guess from 
            used only to retrieve "function" or "function_multi" 
            and "error_propagation_multi"
        p,e Minuit best fit parameter values and std
        _th_runs_ = list of run numbers in suite
    output: for all parameters
        names list of lists of parameter names
        pars list of lists ofparameter values
        epars list of lists of parameter errors
    used only in summary_multirun_global that prints name value(error) 
        one or more lines of global user parameters (the first list in the inner lists)
        one line per run local user parameters and local component parameters (the others)
    '''
    # 
    # initialize
    #
    from mujpy.tools.tools import multigroup_in_components

    names, pars, epars = [], [], []
    nameloc, npars, n_locals = [], -1, 0# inner list, components
    name, par, epar = [], [], []
    for k, pardict in enumerate(dashboard['userpardicts_guess']):
        if not pardict['local']:
            # name, par, epar are lists of globals
            name.append(pardict['name'])
            par.append(p[k])
            epar.append(e[k])
            npars += 1 
        else:
            # nameloc are bare names of locals 
            n_locals += 1
            nameloc.append(pardict['name'])
    # store the globals in name[0], pars[0], epars[0]
    names.append(name) 
    pars.append(par)
    epars.append(epar)
    model = dashboard['model_guess']
    for run in range(len(_the_runs_)):
        name, par, epar = [], [], []# inner list
        for k in range(n_locals):
            npars += 1
            # for brevity name appends a progressive index, not the run number as in minuit
            name.append(nameloc[k]+str(run)) # if run in _the_runs, use run.get_runNumber_int()
            par.append(p[npars])
            epar.append(e[npars])
        for component in model:  # scan the model components
            component_name = component['name']
            label = component['label']
            for j,pardict in enumerate(component['pardicts']): 
                if not pardict['flag']=='=': 
                    npars += 1  # internal parameter index incremented always
                    name.append('{}.{}_{}'.format(pardict['name'],label,str(run)))
                    par.append(p[npars]) 
                    epar.append(e[npars])
        names.append(name)
        pars.append(par)
        epars.append(epar) 
    return names, pars, epars  # list of lists of parameter names, values, errors

def min2int_multigroup(dashboard,p,e):
    '''
    input:
        userpardicts_guess from dashboard 
            (each dict corresponds to a Minuit parameter)
x\            used only to retrieve "function" 
            and "error_propagation_multi"
            p,e Minuit best fit parameter values and std
    output: for all parameters
        namesg list of lists of dashboard parameter names
        parsg list of lists of dashboard parameter values
        eparsg list of lists of dashboard parameter errors
    used only in summary_global
    '''
    # 
    # initialize
    #
    from mujpy.tools.tools import multigroup_in_components
    # print('min2int_multigroup in tools debug: dash {}'.format(dashboard))
    mask_function_multi = multigroup_in_components(dashboard)

    userpardicts = dashboard['userpardicts_guess']  
    e = [e[k] if pardict['flag']=='~' else 0 for k,pardict in enumerate(userpardicts)]
    # names = [pardict['name'] for pardict in userpardicts]

    model = dashboard['model_guess']
    pardicts = [pardict for component in model for pardict in component['pardicts']]
    ngroups = len(pardicts[mask_function_multi.index(1)]["function_multi"])
    nint = -1 # initialize
    namesg, parsg, eparsg = [], [], []
    for l in range(ngroups):
        nint0 = nint
        names, pars, epars = [], [], []
        for component in model:  # scan the model components
            component_name = component['name']
            label = component['label']
            # nint = nint0
            name, par, epar = [], [], [] # inner list, components
            for j,pardict in enumerate(component['pardicts']): 
                nint0 += 1  # internal parameter index incremented always 
                if j==0:
                    name.append('{}: {}_{}'.format(component_name,pardict['name'],label))
                else:
                    name.append('{}_{}'.format(pardict['name'],label))
                if mask_function_multi[nint0]:
                    par.append(eval(pardict["function_multi"][l])) 
                    try:
                        epar.append(eval(pardict["error_propagate_multi"][l]))
                    except:
                        epar.append(eval(pardict["function_multi"][l].replace('p','e')))
                else:                
                    par.append(eval(pardict["function"])) # the function will be eval-uated inside mucomponents
                    try:
                        epar.append(eval(pardict["error_propagate"]))
                    except:
                        epar.append(eval(pardict["function"].replace('p','e')))
                # print('tools min2int_multigroup debug: nint = {} name = {} par = {}, epar = {}'.format(nint0,name[-1],par[-1],epar[-1]))
            pars.append(par) # middel list, model
            names.append(name)
            epars.append(epar)
        namesg.append(names)
        parsg.append(pars)
        eparsg.append(epars) 
    return namesg, parsg, eparsg  # list of parameter values
    
def minglobal2sequential(p_out,p_in,method_keys,dashboard):
    '''
    translate global best fit results (values) into nruns equivalent sequential fits
    for plotting purposes: mufitplot(plot_range,the_fit) will access
    self.fit.lastfits and self.fit.dashboard_single if self.fit.C1 is True
    input:
       p_out is global best fit minuit parameters self.lastfit
       p_in is global fit minuit guess parameters values_in
       method_key is produced by int2_multirun_user_method_key
       dashboard is the fit global dashboard
       results toggles between plotting result of guess
    output:
       lastfits,  list of lists of best fits, in the style of multirun sequential single group (B1)
       dashboard_single (a dashboard that produces a single run best fit function for animated plots) 
    '''
    from copy import deepcopy    
    # values are: 
    #     first the user globals, 
    #     then the run replica: first the user locals, then the free model parameters 
    # for each run must reconstruct a simple model_guess, model_result dashboard (no userpardicts)
    # with its '=' flag parameters translated in '!' flag with value calculated by key (from function)
    lastfits = []
    userpars_g, userpars_r = dashboard["userpardicts_guess"],dashboard["userpardicts_result"]
    n_locals =  [pardict["local"] for pardict in userpars_g].count(True)
    n_globals = len(userpars_g)-n_locals
    # could also simply transfer p = self.fit.lastfit.values
#    for krun in range(len(values[1:])): # i.e. in range(n_runs)
    par = [] # list of single-run equivalent best fit parameters
    # scan method_keys = [[mthd,[[key,...,key],...,[key,,...,key]]],..,[mthd,[[key,...,key]...]]]
    #                   outer components, middle runs, inner component parameters
    dash = deepcopy(dashboard)
    dash.pop("userpardicts_guess")
    if "userpardicts_result" in dash.keys(): dash.pop("userpardicts_result")
    model_in, model_out = [],[]
    if "model_result" not in dash.keys():
        dash["model_result"] = dash["model_guess"] # make sure model_result exists in dash
    for jcomp,method_key in enumerate(method_keys): # components
        method, keys = method_key # component method and run-by-run list of keys
        model_guess_component = dash["model_guess"][jcomp]
        model_result_component = dash["model_result"][jcomp] 
        pars_in, pars_out = [], []
        for krun, runkeys in enumerate(keys): # single run                    
            par_in, par_out = [],[] # each run has its set of parameters
            for kpar,key in enumerate(runkeys):
                par_in.append(key(p_in))
                par_out.append(key(p_out))
                # transform all keys in "~" or "!"
                if model_guess_component["pardicts"][kpar]["flag"]=='=': # turn them in "!"
                    model_guess_component["pardicts"][kpar]["flag"] = '!'
                    model_result_component["pardicts"][kpar]["flag"] = '!' 
                if type(model_guess_component["pardicts"][kpar]["value"]) is list: 
                    model_guess_component["pardicts"][kpar]["value"] = model_guess_component["pardicts"][kpar]["value"][0]
                    model_result_component["pardicts"][kpar]["value"] = model_result_component["pardicts"][kpar]["value"][0]
                if type(model_guess_component["pardicts"][kpar]["error"]) is list: 
                    model_guess_component["pardicts"][kpar]["error"] = model_guess_component["pardicts"][kpar]["error"][0]
                    model_result_component["pardicts"][kpar]["error"] = model_result_component["pardicts"][kpar]["error"][0]
            pars_in.append(par_in) # list of run lists has its set of parameters
            pars_out.append(par_out) 
        dash["model_guess"][jcomp] = model_guess_component
        dash["model_result"][jcomp] = model_result_component
        model_in.append(pars_in)
        model_out.append(pars_out)
# model is in [model [run [component]]] nesting order
    pars_in, pars_out, p_in, p_out = [], [], [], []
    n_runs = len(model_in[0]) #
    n_components = len(model_in)      
    for component_in,component_out in zip(model_in,model_out): # component contains all runs
        for krun in range(n_runs):       
            pars_in.append(component_in[krun]) # 
            pars_out.append(component_out[krun])
    for krun in range(n_runs):
        par_in, par_out  = [], []
        for component in range(n_components):
            par_in += pars_in[krun+component] 
            par_out += pars_out[krun+component] 
        p_in.append(par_in)
        p_out.append(par_out)
    #print('debug tools minglobal2sequential: p_in = {}'.format(p_in))
       
# now the list is in the [run [model [component]]] nesting order
    return p_in,p_out, dash
            
def len_print_components(names,values,errors):
	'''
	input: for a component
		parameter names 
		parameter values 
		parameter errors 
	output:
	    max length of string to print, e.g.
	    "bl.A_fast 0.123(4) bl.λ_fast 12.3(4) bl.σ_fast 0(0)"
	'''
	from mujpy.tools.tools import value_error
	out = [' '.join([names[k],'=',value_error(values[k],errors[k])]) for k in range(len(names))]
	maxlen = len(max(out,key=len))
	return maxlen
    
def print_components(names,values,errors,maxlen):
	'''
	input: for a component
		parameter names 
		parameter values 
		parameter errors 
	output:
	    string to print, e.g.
	    "bl.A_fast 0.123(4) bl.λ_fast 12.3(4) bl.σ_fast 0(0)"
	'''
	from mujpy.tools.tools import value_error
	out = [' '.join([names[k],'=',value_error(values[k],errors[k])]) for k in range(len(names))]
	out = [out[k]+(maxlen-len(out[k]))*' ' for k in range(len(out))]
	return " ".join(out)
	
def len_print_components_multirun(names,values,errors):
	'''
	input: for a component
		parameter names 
		parameter values 
		parameter errors 
	output:
	    max length of string to print, e.g.
	    "bl.A_fast 0.123(4) bl.λ_fast 12.3(4) bl.σ_fast 0(0)"
	'''
	from mujpy.tools.tools import value_error
	outname = [' '+names[k] for k in range(len(names))]
	outval = [' '+value_error(values[k],errors[k]) for k in range(len(names))]
	maxlen = max(len(max(outname,key=len)),len(max(outval,key=len)))
	return maxlen
    
def print_components_multirun(names,values,errors,maxlen):
	'''
	input: for a component
		parameter names 
		parameter values 
		parameter errors 
	output:
	    strings to print, e.g.
	    "A.fast    λ.fast    σ.fast"
	    "0.123(4)  12.3(4)   0(0)"
	'''
	from mujpy.tools.tools import value_error
	outnam = [' '+names[k] for k in range(len(names))]
	outnam = [outnam[k]+(maxlen-len(outnam[k]))*' ' for k in range(len(outnam))]
	outval = [' '+value_error(values[k],errors[k]) for k in range(len(names))]
	outval = [outval[k]+(maxlen-len(outval[k]))*' ' for k in range(len(outval))]
	return "".join(outnam), "".join(outval)
	
def mixer(t,y,f0):
    '''
    mixer of a time-signal with a reference 
    input
        t time
        y the time-signal
        f0 frequency of the cosine reference
    output
        y_rrf = 2*y*cos(2*pi*f0*t)  
    t is 1d and y is 1-d, 2-d or 3-d but t.shape[0] == y.shape[-1]
    t is vstack-ed to be the same shape as y
    '''
    from mujpy.tools.tools import fft_filter
    from numpy import pi, cos, vstack, fft, delete
    ydim, tdim = len(y.shape), len(t.shape)
    # print('tools mixer debug 1: y t shape {}, {}'.format(y.shape,t.shape))
    if tdim == 1: # must replicate t to the same dimensions as y 
        if ydim ==2:
            for k in range(ydim):
                if k:
                    time = vstack((time,t))
                else:
                    time = t
            t = time
        elif ydim==3: # max is ydim = 3
            for j in range(len.shape[-1]):
                for k in len.shape[-2]:
                    if k:
                        time = vstack((time,t))
                    else:
                        time = t
                if j:
                    for l in len.shape[-1]:
                        tim = vstack((tim,time))
                    else:
                        tim = time
            t = tim 
    n = t.shape[-1] # apodize by zero padding to an even number
    yf = fft.irfft(fft_filter(t,fft.rfft(2*y*cos(2*pi*f0*t),n=n+1),f0),n=2*n)
    # now delete padded zeros 
    mindex = range(n,2*n)
    yf =delete(yf,mindex,-1)
    # print('tools mixer debug 3: yf shape {}'.format(yf.shape))
    return yf
    
def fft_filter(t,fy,f0):
    '''
    filter above 0.2*fy peak freq 
    works for 1-2 d
    '''
    from numpy import arange, mgrid, where
    # determine max frequency fmax
    leny = len(fy.shape)
    if leny == 1:
        dt = t[1]-t[0]
        # array f of fourier component indices (real fft, 0 to fmax)
        m = fy.shape
        f = arange(m) 
    elif leny == 2:
        dt = t[0,1]-t[0,0]
        # find peak in rfft below the rrf frequency f0
        # array f of fourier component indices (real fft, 0 to fmax)
        n,m = fy.shape
        _,f = mgrid[0:n,0:m] 
    else:
        dt = t[0,0,1]-t[0,0,0]
        l,n,m = fy.shape
        _,_,f = mgrid[0:l,0:n,0:m] 
                
    fmax = 1/2/dt
    mask = (f<=f0/fmax*m).astype(int)
    # find where fy has a peak, below the rrf frequency f0
    if leny == 1:
        npeak = where(abs(fy)==abs(mask*fy).max()).max()
    elif leny == 2:
        npeak = where(abs(fy)==abs(mask*fy).max())[1].max()
    else:
        npeak = where(fy==(mask*fy).max())[2].max()    
    mask = (f<=2*npeak).astype(int)
    # print('tools fft_filter debug 2: fy {},mask {} shape'.format(fy.shape,mask.shape))    
    return fy*mask
    
def model_name(dashboard):
    '''
    input the dashboard dictionary structure
    output the model name (e.g. 'mgbgbl') 
    '''    
    return ''.join([item for component in dashboard["model_guess"] for item in component["name"]])
    
def userpars(dashboard):
    '''
    checks if there are userpardicts in the fit dashboard
    alias of global type fit, of any kind (gg, gr, G)
    used by fit and plt switchyard
    '''
    return "userpardicts_guess" in dashboard

def userlocals(dashboard):
    '''
    input:
        full dashboard
    output:
        True is "userpardicts_local" in dashboard.keys 
    '''
    return "userpardicts_local" in dashboard    

def multigroup_in_components(dashboard):
    '''
    input full dashboard
    output mask list, 
        1 where "model_guess" 
        contains at least one component (dict) 
        whose "pardicts" (list) 
        contains a parameter dict 
        with at least one "function_multi":[string, string ..] key
        0 otherwise
    '''

    #print('multigroup_in_components tools debug: model_guess len {}'.format(len(dashboard["model_guess"])))
    #print('multigroup_in_components tools debug: pardicts len {}'.format(len(dashboard["model_guess"][0]['pardicts'])))
    #print('multigroup_in_components tools debug: pardict.keys len {}'.format(len(dashboard["model_guess"][0]['pardicts'][0].keys())))

    return ['function_multi' in pardict.keys() for component in dashboard["model_guess"]  for pardict in component["pardicts"]]                                
                            
    # contains 1 for all parameters that have "function_multi", 0 otherwise 
    # return [k for k,component in enumerate(component_function) if component>0]

def stringify_groups(groups):
    '''
    returns a unique string for many groups
    to use in json file name
    '''
    strgrp = []
    for group in groups: 
        fgroup, bgroup = group['forward'],group['backward']
        strgrp.append(fgroup.replace(',','_')+'-'+bgroup.replace(',','_'))
    return '_'.join(strgrp)

def modelstrip(name):
    '''
    strips numbers of external parameters at beginning of model name
    '''
    import re
    nexternals, ngroups = 0, 0
    # strip the name and extract number of external parameters
    try:
        nexternals = int('{}'.format(re.findall('^([0-9]+)',name)[0]))
        if nexternals>9:
            return []
        name = name[:-1]
    except:
        pass
#    try:
#        ngroups = int('{}'.format(re.findall('([0-9]+)$',name)[0]))
#        if ngroups>9:
#            return []
#        name = name[1:]
#    except:
#        pass
    return name, nexternals
            


##############
# MUGUI AUX ?
##############

def name_of_model(model_components,model):
    '''
    check if model_components list of dictionaries correstponds to model
    '''
    content = []
    for component in model_components:
        content.append(component["name"])
    return True if ''.join(content) == model else False

def create_model(model):
    '''
    create_model('daml') # adds e.g. the two component 'da' 'ml' model
    this method 
    does not check syntax (prechecked by checkvalidmodel)
       ? separates nexternals number from model name (e.g. '3mgml' -> 'mgml', 3)
       ?  starts switchyard for A1,A1,B1, B2, C1, C2 fits
    adds a model of components selected from the available_component tuple of  
    directories
    with zeroed values, stepbounds from available_components, flags set to '~' and zeros functions
    '''
    import string
    from mujpy.tools.tools import addcomponent, _available_components_
    # print('create_model: {}'.format(model))
    components = [model[i:i+2] for i in range(0, len(model), 2)]
    model_guess = [] # start from empty model
    for k,component_name in enumerate(components):
        component, emsg = addcomponent(component_name) # input a component name, output a component dictionary
        if component:
            model_guess.append(component) # list of dictionaries                
        # self.console('create model added {}'.format(component+label))
        else:
             return False, emsg

    return model_guess, '' # list of component dictionaries

def addcomponent(name):
    '''
    addcomponent('ml') # adds e.g. a mu precessing, lorentzian decay, component
    this method adds a component selected from _available_components_(), tuple of directories
    with zeroed values, error and limits from available_components, 
    flags set to '~' and zeros functions
    [plan also addgroupcomponents and addruncomponents (for A2, B2, C1, C2)]
    '''
    from copy import deepcopy
    from mujpy.tools.tools import _available_components_
    available_components =_available_components_() # creates list automagically from mucomponents
    component_names = [available_components[i]['name'] 
                            for i in range(len(available_components))]
    if name in component_names:
        k = component_names.index(name)
        npar = len(available_components[k]['pardicts']) # number of pars
        pars = deepcopy(available_components[k]['pardicts']) # list of dicts for 
        # parameters, {'name':'asymmetry','error':0.01,'limits':[0, 0]}

        # now remove parameter name degeneracy                   
        for j, par in enumerate(pars):
            pars[j]['name'] = par['name']
            if par['name']=='α':
                pars[j].update({'value':1.0}) # initilize
            elif par['name']=='A':
                pars[j].update({'value':0.1}) # initialize to not zero
            elif par['name']=='B':
                pars[j].update({'value':2.}) # initialize to TF20
            else:
                pars[j].update({'value':0}) # does not need initialization
            pars[j].update({'flag':'~'})
            pars[j].update({'function':''}) # adds these three keys to each pars dict
            pars[j]['error'] = par['error']
            pars[j]['limits'] = par['limits']                    
            # they serve to collect values in mugui
        # self.model_guess.append()
        return {'name':name,'label':'','pardicts':pars}, None # OK code, no message
    else:
        # self.console(
        error_msg = '\nWarning: '+name+' is not a known component. Not added.'
        return {}, error_msg # False error code, message


def chi2std(nu):
    '''
    computes 1 std for least square chi2
    '''
    import numpy as np
    from scipy.special import gammainc
    from scipy.stats import norm
    
    mm = round(nu/4)              
    hb = np.linspace(-mm,mm,2*mm+1)
    cc = gammainc((hb+nu)/2,nu/2) # see mulab: muchi2cdf(x,nu) = gammainc(x/2, nu/2);
    lc = 1+hb[min(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nu
    hc = 1+hb[max(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nu
    return lc, hc

def component(model,kin):
    '''
    returns the index of the component to which parameter k belongs in
    model = self.model_guess, in mugui, a list of complex dictionaries::
            [{'name':'da', 'pardicts':{'name':'lpha',...},
            {'name':'mg', 'pardicts':{       ...       }]
            
    kin is the index of a dashboard parameter (kint)
    '''
    from numpy import array, cumsum, argmax
    
    ncomp = len(model) # number of components in model
    npar = array([len(model[k]['pardicts']) for k in range(ncomp)]) # number of parameters of each component
    npars = cumsum(npar)
    return argmax(npars>kin)

#############
# GENERAL AUX
#############

def calib(dashboard):
    '''
    True if the first component is 'al'
    '''
    return dashboard['model_guess'][0]['name']=='al'
            
def derange(string,vmax,pack=1):
    '''
    derange(string,vmax,pack=1) 
    reads string 
    assumes it contains 2, 3, 4 or 5 csv or space separated values
    uses isinstance(vmax,float) to distinguish floats (fft) from integers (fit and plot) 

        5: start, stop, packe, last, packl       # for plot
        4: start, stop, last, packl              # for plot (packe is 1) 
        3: start, stop, pack
        2: start, stop (pack is added, pack default is 1)

    returns 2, 3, 4 or 5 floats or int, or 
    default values, 0,vmax,pack, if fails validity check (stop>start, bin <stop-start, last < vmax) 
    errmsg = '' in ok, a string indicates errors       
    '''
    
    # print('In derange, string = {}'.format(string))
    errmsg = ''
    x_range = string.split(',') # assume ',' is the separator
    if len(x_range)==1: # try ' ' as separator
        x_range = string.split(' ')
    if len(x_range)==1: # wrong syntax
        x_range = [vmax-vmax,vmax,pack] # default, int for int vmax, float for float vmax
        errmsg = 'no range'
    if not errmsg:
        try: # three items are they integers floats or misprints?
            if isinstance(vmax,float): # should be three floats
                x_range = [float(chan) for chan in x_range] # breaks if non digits in x_range 
            else: # should be three integers
                x_range = [int(chan) for chan in x_range] # breaks if non digits in x_range 
            if len(x_range)==2: # guarantees three items
                x_range.append(pack)
            if x_range[2]>(x_range[1]-x_range[0])//2: # True for fit_range[1]<fit_range[0]  or too large pack
                raise Exception
        except:
            x_range = [vmax-vmax,vmax,pack] # default
            errmsg = 'Syntax error, reset range to default. '
    # to re-compose a correct string use
    # string = ','.join([str(val) for val in x_range])
    # print('tools derange: x_range = {}'.format(x_range))
        
    return x_range, errmsg # a list of values (int or float as appropriate)
    
def derun(string):
    '''
    parses string, producing a list of runs; 
    expects comma separated items

    looks for 'l','l:m','l+n+m','l:m:-1' 
    where l, m, n are integers
    also more than one, comma separated 

    rejects all other characters

    returns a list of lists of integer
    '''
    import re
    s = []
    try:
    # substitute multiple consecutive spaces with ','
        string = re.sub("\s+", ",", string.strip())
    # systematic str(int(b[])) to check that b[] ARE integers
        for b in string.split(','): # csv
            kminus = b.find(':-1') # '-1' means reverse order
            kcolon = b.find(':') # ':' and '+' are mutually exclusive
            kplus = b.find('+')
            #print(kminus,kcolon,kplus)

            if kminus<0 and kcolon<0 and kplus<0: # single run
                int(b) # produces an Error if b is not an integer
                s.append([b]) # append single run string   
            else:
                if kminus>0 and kminus == kcolon:
                    return [], 'l:-1 is illegal'
                elif kplus>0:
                    # add files, append a list or run strings
                    ss = []
                    k0 = 0
                    while kplus>0: # str(int(b[]))
                        ss.append(int(b[k0:kplus])) 
                        k0 = kplus+1
                        kplus = b.find('+',k0)
                    ss.append(int(b[k0:]))
                    s.append([str(q) for q in ss])
                else:
                    # either kminus=-1 (just a range) or  kcolon<kminus, (range in reverse order)
                    # in both cases:
                    if kminus<0:
                        #print(int(b[:kcolon]),int(b[kcolon+1:]))
                        if int(b[:kcolon])>int(b[kcolon+1:]):
                            return [], 'l:m must have l<m'
                        for j in range(int(b[:kcolon]),int(b[kcolon+1:])+1):
                            s.append([str(j)]) # append single run strings
                    else:
                        ss = [] 
                        # # :-1 reverse order
                        if int(b[:kcolon])>int(b[kcolon+1:kminus]):
                            return ss, 'l:m:-1 must have l<m'
                        for j in range(int(b[:kcolon]),int(b[kcolon+1:kminus])+1):
                            ss.append([str(j)]) # append single run strings
                        ss = ss[::-1]
                        for sss in ss:
                            s.append(sss)
        return s, None
    except:
        return [], 'error to be debugged'
        
def run_shorthand(runstrings):
    '''
    write the runlist contained in runstrings (suite self.runs produced by derun)
        i.e. a list of lists, with separate run numbers in string format, the inner ones  to be added 
    in a compact string, with space separated notation
    e.g.
    '650:655,675,656:674' 
    '''
    # [[623],[624],[625],[626],[627,628,629], [631],[632],[633],[630]] -> 623:626 627+628+629 631:633 630
    runlists = [[int(run) for run in runstringlist] for runstringlist in runstrings]
    string = [[] for i in range(len(runlists))]
    index_runadds = [i for i in range(len(runlists)) if len(runlists[i])>1]
    index_runs = [i for i in range(len(runlists)) if len(runlists[i])==1]
    for j in index_runadds:
        string[j]='+'.join([str(k) for k in runlists[j]])
    k =  index_runs[0]
    string[k].append(runlists[k][0]) 

    up = None
    for j in range(1,len(string)):
        if j in index_runadds:
            k = j
        elif runlists[j - 1][0] + 1 == runlists[j][0]:
            if up != True:
                up = True
                k = j-1 
            string[k].append(runlists[j][0])
        elif runlists[j - 1][0] - 1 == runlists[j][0]:
            if up != False:
                up = False
                k = j-1
            string[k].append(runlists[j][0])
        else:
            up = None
            k = j
            string[k].append(runlists[j][0])
    ss = list(filter(([]).__ne__,string))
    for k,l in enumerate(ss):
        if isinstance(l,list):
            if len(l)==1:
                ss[k] = str(l[0])
            else:
                ss[k] = str(l[0])+':'+str(l[-1])
    s = ','.join(ss)
    return s

def findall(p, s):
    '''Yields all the positions of
    the pattern p in the string s.
    
    Used by translate.
    '''
    i = s.find(p)
    while i != -1:
        yield i
        i = s.find(p, i+1)

def find_nth(haystack, needle, n):
    '''
    Finds nth needle in haystack 

    Returns its first occurrence (0 if not present)

    Used by ?
    '''
    start = haystack.rfind(needle)
    while start >= 0 and n > 1:
        start = haystack.rfind(needle, 1, start-1)
        n -= 1
    return start
    
def get_datafilename(datafile,run):
    '''
    datafilename = template, e.g. '/fullpath/deltat_gps_tdc_0935.bin'
    run = string of run digits, e.g. '1001'
    returns '/fullpath/deltat_gps_tdc_1001.bin'
    '''
    
    import re
    dot_suffix = datafile[-4:]
    padded = re.match('.*?([0-9]+)$', datafile[:-4]).group(1) # run string of digits
    oldrun = str(int(padded)) # strip padding zeros
    datafileprefix = datafile[:datafile.find(oldrun)] # prefix up to original zero padding
    if len(run)-len(oldrun)>0:
        datafilename = datafileprefix[:len(oldrun)-len(run)]+run+dot_suffix
    elif len(run)-len(oldrun)==-1:
        datafilename = datafileprefix+'0'+run+dot_suffix
    elif len(run)-len(oldrun)==-2:
        datafilename = datafileprefix+'00'+run+dot_suffix
    elif len(run)-len(oldrun)==-3:
        datafilename = datafileprefix+'000'+run+dot_suffix
    else:
        datafilename = datafileprefix+run+dot_suffix
    return datafilename

def get_datafile_path_ext(datafile,run):
    '''
    datafilename = template, e.g. '/fullpath/deltat_gps_tdc_0935.bin'
    run = string of run digits, e.g. '1001'
    returns '/fullpath/deltat_gps_tdc_1001.bin'
    '''
    import os
    path = datafile[:datafile.rfind(os.path.sep)+1] # e.g. /afs/psi.ch/projec/bulkmusr/data/gps/d2022/tdc/', works in  WIN with '\' as separator
    fileprefix = datafile[datafile.rfind(os.path.sep)+1:datafile.rfind('.')]
    ext = datafile[datafile.rfind['.']+1-len(datafile)] # e.g. 'bin' or 'nxs'
    return path, fileprefix, ext

def get_grouping(groupcsv):
    """
    input
      groupcsv is a shorthand csv string, e.g. '1:3,5' or '1,3,5' etc.
      contained in self.suite.group[k]["forward] of self.suite.group[k]["backward"]
          (the k-th detector group of a multi group fit)
    output
     grouping is an np.array of indces, 0 based
    """
    import numpy as np

    # two shorthands: either a list, comma separated, such as 1,3,5,6 
    # or a pair of integers, separated by a colon, such as 1:3 = 1,2,3 
    # only one column is allowed, but 1, 3, 5 , 7:9 = 1, 3, 5, 7, 8, 9 
    # or 1:3,5,7 = 1,2,3,5,7  are also valid
    # no more complex nesting (3:5,5,8:10 is not allowed)
    #       get the shorthand from the gui Text 
    groupcsv = groupcsv.replace('.',',') # can only be a mistake: '.' means ','
    try:
        if groupcsv.find(':')==-1: # no colon, it's a pure csv
            grouping = np.array([int(ss) for ss in groupcsv.split(',')]) # read it
        else:  # colon found                 
            if groupcsv.find(',')==-1: # (no commas, only colon, must be n:m)
                nm = [int(w) for w in groupcsv.split(':')] # read n m
                grouping = np.array(list(range(nm[0],nm[1]+1))) # single counters
            else: # general case, mixed csv and colon
                p = groupcsv.split(':') # '1,2,3,4,6' '7,10,12,14' '16,20,23'
                ncolon = len(p)-1 
                grouping = np.array([])
                for k in range(ncolon):
                    q = p[k].split(',') # ['1' '2' '3' '4' '6']
                    if k>0:
                        last = int(q[0])
                        grouping = np.concatenate((grouping,np.array(list(range(first,last+1)))))
                        first = int(q[-1])
                        grouping = np.concatenate((grouping,np.array(list(int(w) for w in q[1:-1]))))
                    elif k==0:
                        first = int(q[-1])
                        grouping = np.concatenate((grouping,np.array(list(int(w) for w in q[:-1]))))
                q = p[-1].split(',') # '22','25'
                last = int(q[0])
                grouping = np.concatenate((grouping,np.array(list(range(first,last+1)))))
                grouping = np.concatenate((grouping,np.array(list(int(w) for w in q[1:]))))
        grouping -=1 # this is counter index, remove 1 for python 0-based indexing 
    except:
        grouping = np.array([-1]) # error flag
        
    return grouping
    
def get_group(grouping):
    '''
    reverse of get_grouping, 
    input 
        grouping is an np.array of indices of detectors , 0 based 
    output is 
        groupcsv shorthand as in self.group[k]["forward"} or self.group[k]["backward"}
          e.g '1:3,5' or '1,3,5' etc.
    '''
    import numpy as np
    # find sequences
    groups = []
    if grouping.size>1:
        grouping = np.sort(grouping)+1 # 1 base for csv
        gsequences = np.split(grouping, np.where(np.diff(grouping) != 1)[0]+1)
        for gsequence in gsequences:
            gstring = str(gsequence) if gsequence.size==1 else str(gsequence[0])+':'+str(gsequence[-1])
        groups.append(gstring)
    else:
        groups.append(str(grouping[0]))
    return ','.join(groups)
    
def getname(fullname):
    '''
    estracts parameter name from full parameter name (i.e. name + label)
    for the time being just the first letter
    '''
    return fullname[0]

def initialize_csv(Bstr, filespec, the_run ):
    '''
    writes beginning of csv row 
    with nrun T [T eT T eT] B 
    for ISIS [PSI]
    '''
    nrun = the_run.get_runNumber_int()
    # print('tools initialize_csv debug: nrun {}'.format(nrun))
    if filespec=='bin' or filespec=='mdu':
        TsTc, eTsTc = the_run.get_temperatures_vector(), the_run.get_devTemperatures_vector()
        n1,n2 = spec_prec(eTsTc[0]),spec_prec(eTsTc[1]) # calculates format specifier precision
        form = '{},{:.'
        form += '{}'.format(n1)
        form += 'f},{:.'
        form += '{}'.format(n1)
        form += 'f},{:.'
        form += '{}'.format(n2)
        form += 'f},{:.'
        form += '{}'.format(n2)
        form += 'f},{}' #".format(value,most_significant)'
        return form.format(nrun, TsTc[0],eTsTc[0],TsTc[1],eTsTc[1], Bstr[:Bstr.find('G')])
    elif filespec=='nxs':
        Ts = the_run.get_temperatures_vector()
        n1 = '1'       
        form = '{},{:.'
        form += '{}'.format(n1)
        form += 'f},{}' #".format(value)
        return form.format(nrun, Ts[0], Bstr[:Bstr.find('G')])

def minparam2_csv(dashboard,values_in,errors_in,multirun=0):
    '''
    transforms Minuit values Minuit errors in cvs format
    input:
        dashboard 
                 dashboard["model_guess"], for single group 
                 None, for multi group
                 dashboard, for single group multirun user (C1) 
        values_in, errors_in are  Minuit values Minuit errors
        if multirun is nruns !=0  (True) uses 
           min2int_multirun(dashboard,values_in,errors_in,multirun_nruns)
        else (multirun = 0 (False) uses 
           min2int(dashboard,values_in,errors_in)
    output:
        cvs partial row with parameters and errors for A1, A20 and B1, or A21
            list of partial rows (one per run) for C1
    '''
    from mujpy.tools.tools import min2int, min2int_multirun, spec_prec

    if multirun:
        _, values, errors = min2int_multirun(dashboard,values_in,errors_in,multirun)
        # must write rows with single run
        gvalues,gerrors = values[0],errors[0]
        rows = []
        for parvalues,parerrors in zip(values[1:],errors[1:]): #locals
            row = ''
            for parvalue,parerror in zip(parvalues,parerrors):
                n1 = spec_prec(parerror) # calculates format specifier precision
                form = ',{:.'
                form += '{}'.format(n1)
                form += 'f},{:.'
                form += '{}'.format(n1)
                form += 'f}'
                row += form.format(parvalue,parerror)
            for parvalue,parerror in zip(values[0],errors[0]): # globals replicated in every row
                n1 = spec_prec(parerror) # calculates format specifier precision
                form = ',{:.'
                form += '{}'.format(n1)
                form += 'f},{:.'
                form += '{}'.format(n1)
                form += 'f}'
                row += form.format(parvalue,parerror)
            rows.append(row) # these are as many rows as runs     
    else:
        (_, values, errors) = (min2int(dashboard,values_in,errors_in) if dashboard else
                                                        (None, [values_in], [errors_in]))    
        # from minuit parameters to component parameters
        # output is lists (components) of lists (parameters) 
        # else dashboard false (multigroup user) Minuit and user parameters coincide
        rows = '' # this is a single row, really
        for parvalues, parerrors in zip(values,errors): 
            for parvalue,parerror in zip(parvalues,parerrors):
                n1 = spec_prec(parerror) # calculates format specifier precision
                form = ',{:.'
                form += '{}'.format(n1)
                form += 'f},{:.'
                form += '{}'.format(n1)
                form += 'f}'
                rows += form.format(parvalue,parerror)
    return rows
    
def nextrun(datapath):
    '''
    assume datapath is path+fileprefix+runnumber+extension
    datafile is next run, runnumber incremented by one
    if datafile exists return next run, datafile
    else return runnumber and datapath
    '''
    import os
    from mujpy.tools.tools import muzeropad

    path, ext = os.path.splitext(datapath)
    lastchar = len(path)
    for c in reversed(path):
        try:
            int(c)
            lastchar -= 1
        except:
            break
    run = path[lastchar:]
    runnext = str(int(run)+1)
    datafile = path[:lastchar]+muzeropad(runnext)+ext
    run = runnext if os.path.exists(datafile) else run
    datafile = datafile if os.path.exists(datafile) else datapath                         
    return run, datafile

def thisrun(datapath):
    '''
    assume datapath is path+fileprefix+runnumber+extension
    datafile is present run
    if datafile exists returns path to datafile
    '''
    import os
    from mujpy.tools.tools import muzeropad

    path, ext = os.path.splitext(datapath)
    lastchar = len(path)
    for c in reversed(path):
        try:
            int(c)
            lastchar -= 1
        except:
            break
    run = path[lastchar:]
    datafile = path[:lastchar]+muzeropad(run)+ext
    return datafile

def prevrun(datapath):
    '''
    assume datapath is path+fileprefix+runnumber+extension
    datafile is prev run, runnumber decremented by one
    if datafile exists return prev run, datafile
    else return runnumber and datapath
    '''
    import os
    from mujpy.tools.tools import muzeropad

    path, ext = os.path.splitext(datapath)
    lastchar = len(path)
    for c in reversed(path):
        try:
            int(c)
            lastchar -= 1
        except:
            break
    run = path[lastchar:]
    runprev = str(int(run)-1)
    datafile = path[:lastchar]+muzeropad(runprev)+ext                            
    run = runprev if os.path.exists(datafile) else run
    datafile = datafile if os.path.exists(datafile) else datapath                         

    return run, datafile
    
def chi2_csv(chi2,lowchi2,hichi2,groups,offset):
    '''
    input:
        chi2, chi2-sdt, chi2+sdt, groups, offset (bins)
        groups is suite.groups and its len, 1 or more, identifies multigroup
    output:
        cvs partial row with these values and timestring
    '''
    from time import localtime, strftime
    
    echi = max(chi2-lowchi2,hichi2-chi2)
    n1 = spec_prec(echi) # calculates format specifier precision
    form = ',{:.'
    form += '{}'.format(n1)
    form += 'f},{:.'
    form += '{}'.format(n1)
    form += 'f},{:.'
    form += '{}'.format(n1)
    form += 'f}' # ' {} {}'
    row = form.format(chi2,chi2-lowchi2,hichi2-chi2)
    for group in groups:
        row += ',{}'.format(group["alpha"])
    row += ',{},{}'.format(offset,strftime("%d.%b.%H:%M:%S", localtime()))
    return row

def write_csv(header,row,the_run,file_csv,filespec,scan=None):
    '''
    input :
        header, the model specific csv header 
                to compare with that of the csv file
        row, the line to be added to the csv file
        the_run, run instance (first one for added runs)
        file_csv, full path/filename to csv file 
        filespec, 'bin', 'mdu' or 'nsx'
        scan, T, B or None
    output:
        two strings to write on console
    writes onto csv finding the right line
    writes a new file if csv does not exist or is incompatible (writes ~ version)
    '''
    from mujpy.tools.tools import get_title
    import os
    import re
    from datetime import datetime

    nrun = int(re.split(" |,|, ",row)[0])
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S") 
    if scan==None:  # order by nrun, first item in csv
        csv_index = 0
    elif scan=='T': # order by T, 4th or 2nd item in csv
        csv_index = 3 if filespec == 'bin' or filespec == 'mdu' else 1
    else:           # order by B, 6th or 4th item in csv
        csv_index = 5 if filespec == 'bin' or filespec == 'mdu' else 3
    rowvalue = float(re.split(" |,|, ",row)[csv_index]) # also nrun is transformed into float 

    if os.path.isfile(file_csv):
        try: # the file exists
            lineout = [] # is equivalent to False
            with open(file_csv,'r') as f_in:
                notexistent = True # the line does not exist
                for nline,line in enumerate(f_in.readlines()):
                    if nline==0:
                        if header!=line: # different headers, substitute present one
                            raise # exits this try 
                        else:
                            lineout.append(header)
                    elif float(re.split(" |,|, ",line)[csv_index]) < rowvalue: # append line first
                        lineout.append(line)
                    elif float(re.split(" |,|, ",line)[csv_index]) == rowvalue: # substitute an existing fit
                        lineout.append(row) # insert before last existing fit
                        notexistent = False
                    else: 
                        if notexistent:
                            lineout.append(row) # insert before last existing fit
                            notexistent = False
                        lineout.append(line) # insert all other existing fits
                if notexistent:
                    lineout.append(row) # append at the end
                    notexistent = False
            with open(file_csv,'w') as f_out:                 
                for line in lineout:
                    f_out.write(line)
            file_csv = file_csv[file_csv.rfind('/')+1:]
            return 'Run {}: {} ***'.format(nrun,
                   get_title(the_run)), '.  Log added to {}'.format(file_csv)

        except: # incompatible headers, save backup and write a new file
            os.rename(file_csv,file_csv+'~')
            with open(file_csv,'w') as f:
                f.write(header)
                f.write(row)
            file_csv = file_csv[file_csv.rfind('/')+1:]
            return 'Run {}: {} ***'.format(nrun,
                    get_title(the_run)),'.  Log in NEW {} [backup in {}]'.format(
                                                                         file_csv,
                                                                         file_csv+'~')
            
    else: # csv does not exist
        with open(file_csv,'w') as f:
            f.write(header)
            f.write(row)
        file_csv = file_csv[file_csv.rfind('/')+1:]
        return 'Run {}: {} ***'.format(nrun,
                        get_title(the_run)),'.  Log in NEW {}'.format(file_csv)

def get_title(run,notemp=False,nofield=False):
    '''
    form standard psi title
    '''
    title = [(run.get_sample()).rstrip()]
    title.append((run.get_orient()).rstrip())  
    if not notemp:
        tstr = run.get_temp()
        try:
            temp = float(tstr[:tstr.index('K')])
        except:
            temp = float(tstr)
        title.append('{:.1f}K'.format(temp))
    if not nofield:
        field = run.get_field()
        try:
            title.append('{:.0f}mT'.format(float(field[:field.index('G')])/10))
        except:
            title.append('{:.0f}mT'.format(float(field)/10))
    return ' '.join(title)    
    
def get_run_title(the_suite):
    '''
    output 
        list of run and title strings
            each run and group in the run replicates its run number + title
    used only in mufitplot (fit and fft  
    '''
    from mujpy.tools.tools import get_title
    run_title = []
    for run in the_suite._the_runs_:
        for kgroup in range(len(the_suite.grouping)):
                run_title.append(str(run[0].get_runNumber_int())+'-'+get_title(run[0]))
    return run_title
    
def get_nruns(the_suite):
    '''
    get nrun strings
    '''
    nruns = []
    print
    for k,run in enumerate(the_suite._the_runs_):
        nruns.append(str(run[0].get_runNumber_int()))
    return nruns


def get_run_number_from(path_filename,filespecs):
    '''
    strips number after filespecs[0] and before filespec[1]
    '''
    try:
        string =  path_filename.split(filespecs[0],1)[1]
        run = string.split('.'+filespecs[1],1)[0]
    except:
       run = '-1' 
    return str(int(run)) # to remove leading zeros

def muvalid(string):
    '''
    parse function 

    CHECK WITH MUCOMPONENT, THAT USES A DIFFERENT SCHEME

    accepted functions are RHS of agebraic expressions of parameters p[i], i=0...ntot  
    '''
    import re
    error_message = ''
    if string.strip() !='': # empty and blank strings are validated 
        pattern = re.compile(r"p\[(\d+)\]") # find all patterns p[*] where * is digits
        test = pattern.sub(r"a",string) # substitute "a" to "p[*]" in s
        #           strindices = pattern.findall(string)
        #           indices = [int(strindices[k]) for k in range(len(strindices))] # in internal parameter list
        #           mindices = ... # produce the equivalent minuit indices  
        try: 
            safetry(test) # should select only safe use (although such a thing does not exist!)
        except Exception as e:
            error_message = 'Function: {}. Tested: {}. Wrong or not allowed syntax: {}'.format(string,test,e)
    return error_message
    
def p2x(instring):
    '''
    replaces parameters e.g. p[2] with variable x2 in string
    returns substitude string and list of indices (ascii)
    '''
    import re
    patterna = re.compile(r"p\[(\d+)\]") # find all patterns p[*] where * is digits
    n = patterna.findall(instring) # all indices of parameters
    outstring = instring
    for k in n:
        strin = r"p\["+re.escape(k)+"\]"
        patternb = re.compile(strin)
        stri = r"x"+re.escape(k)  # variable
        outstring = patternb.sub(stri,outstring)
    return outstring, n
    
def errorpropagate(string,p,e):
    '''
    parse function in string 
    
    substitute p[n] with xn, with errors en
    calculate the partial derivative pdn = partial f/partial xn 
    return the sqrt of the sum of (pdn*en)**2
    '''
    from jax import grad
    import numpy as np
    funct,n = p2x(string) # from parameters p[n] to variables xn
    s = 'lambda '
    ss = ['x'+k+',' for k in n]
    args = ''.join(ss)[:-1]
    s = s + args + ': '+funct[1:] # removes the '='
    #  s = 'lambda xn,xm,... : expression of xn, xm, ...'
    f = eval(s) # defines a function of the parameters, called xn, xm, 
    variance = 0
    for k in n:
        exec('x'+k+'= p['+k+']')   # this assigns p[n] value to xn 
        d = grad(f,argnums=int(k)) # this is the derivative with respect to the k-th variable
        ss
        exec('variance += (d('+args+')*e['+k+'])**2')
    return np.sqrt(variance)
    
def group_shorthand(grouping):
    '''
    group_calib is the list of gorup dictionaries
    '''
    shorthand = []
    for group in grouping:
        fwd = '_'.join([str(s+1) for s in group['forward']])
        bkd = '_'.join([str(s+1) for s in group['backward']])
        shorthand.append(fwd+'-'+bkd)
    return '+'.join(shorthand)

def json_name(model,datafile,grouping,version,g=False):
    '''
    model is e.g. 'mlmg'
    datafile is e.g. '/afs/psi.ch/bulkmusr/data/gps/d2022/tdc/deltat_gps_tdc_1233.bin'
       must have a single '.'
    grp_calib is the list of dictionaries defining the groups
    g = True for global
    version is a label
    returns a unique name for the json dashboard file
    '''    
    from re import findall
    from mujpy.tools.tools import group_shorthand
    run = findall('[0-9]+',datafile)[-1]
    return model+'.'+run+'.'+group_shorthand(grouping)+'.'+version+'.json'
    
def muvaluid(string):
    '''
    Run suite fits: muvaluid returns True/False
    * checks the syntax for string function 
    corresponding to flag='l'. Meant for pars
    displaying large changes across the run suite,
    requiring different migrad start guesses::

    # string syntax: e.g. "0.2*3,2.*4,20."
    # means that for the first 3 runs value = 0.2,
    #            for the next 4 runs value = 2.0
    #            from the 8th run on value = 20.0

    '''
    try:
        value_times_list = string.split(',')
        last = value_times_list.pop()
        for value_times in value_times_list:
            value,times = value_times.split('*')
            dum, dum = float(value),int(times)
        dum = float(last)
        return True
    except:
        return False

def muvalue(lrun,string):
    '''
    Run suite fits: 

    muvalue returns the value 
    for the nint-th parameter of the lrun-th run
    according to string (corresponding flag='l').
    Large parameter change across the run suite
    requires different migrad start guesses.
    Probably broken!
    '''
    # string syntax: e.g. "0.2*3,2.*4,20."
    # means that for the first 3 runs value = 0.2,
    #            for the next 4 runs value = 2.0
    #            from the 8th run on value = 20.0

    value = []
    for value_times in string.split(','):
        try:  # if value_times contains a '*' 
            value,times = value_times.split('*') 
            for k in range(int(times)):
                value.append(float(value))
        except: # if value_times is a single value
            for k in range(len(value),lrun):
                value.append(float(value_times))
    # cannot work! doesn't check for syntax, can be broken; this returns a list that doesn't know about lrun
    return value[lrun]

def muzeropad(runs,nzeros=4):
    '''

    runs is a string containing the run number
    nzeros the number of digit chars in the filename
    PSI bin: nzeros=4
    ISIS nxs nzeros=8
    returns the runs string 
    with left zero padding to nzeros digits
    '''
    zeros='0'*nzeros
    if len(runs)<len(zeros):
        return zeros[:len(zeros)-len(runs)]+runs
    elif len(runs)==len(zeros):
        return runs

def path_file_dialog(path,spec):
    import tkinter
    from tkinter import filedialog
    import os
    tkinter.Tk().withdraw() # Close the root window
    spc, spcdef = '.'+spec,'*.'+spec
    in_path = filedialog.askopenfilename(initialdir = path,filetypes=((spc,spcdef),('all','*.*')))
    return in_path

def path_dialog(path,title):
    import tkinter
    from tkinter import filedialog
    import os
    tkinter.Tk().withdraw() # Close the root window
    in_path = filedialog.askdirectory(initialdir = path,title = title)
    
    return in_path

################
# PLOT METHODS #
################

def plot_parameters(nsub,labels,fig=None): 
    '''
    standard plot of fit parameters vs B,T (or X to be implemente)
    input
       nsub<6 is the number of subplots
       labels is a dict of labels, 
       e.g. {title:self.title, xlabel:'T [K]', ylabels: ['asym',r'$\lambda$',r'$\sigma$,...]}
       fig is the standard fig e.g self.fig_pars
       
    output 
       the ax array on which to plot 
       one dimensional (from top to bottom and again, for two columns)
       example 
         two asymmetry parameters are both plotfal=1 and are plotted in ax[0]
         a longitudinal lambda is plotflag=2 and is plotted in ax[1]
         ...
         a transverse sigma is plotflag=n and is plotted in ax[n-1]
         
    '''
    import matplotlib.pyplot as P
    nsubplots = nsub if nsub!=5 else 6 # nsub = 5 is plotted as 2x3 
    # select layout, 1 , 2 (1,2) , 3 (1,3) , 4 (2,2) or 6 (3,2)
    nrc = {
            '1':(1,[]),
            '2':(2,1),
            '3':(3,1),
            '4':(2,2),
            '5':(3,2),
            '6':(3,2)
            }
    figsize = {
                '1':(5,4),
                '2':(5,6),
                '3':(5,8),
                '4':(8,6),
                '5':(8,8),
                '6':(8,8)
                } 
    spaces = {
                '1':[],
                '2':{'hspace':0.05,'top':0.90,'bottom':0.09,'left':0.13,'right':0.97,'wspace':0.03},
                '3':{'hspace':0.05,'top':0.90,'bottom':0.09,'left':0.08,'right':0.97,'wspace':0.03},
                '4':{'hspace':0.,'top':0.90,'bottom':0.09,'left':0.08,'right':0.89,'wspace':0.02},
                '5':{'hspace':0.,'top':0.90,'bottom':0.09,'left':0.08,'right':0.89,'wspace':0.02},
                '6':{'hspace':0.,'top':0.90,'bottom':0.09,'left':0.08,'right':0.89,'wspace':0.02}
                }
    if fig: # has been set to a handle once
       fig.clf()
       if nrc[str(nsub)][1]: # not a single subplot
           fig,ax = P.subplots(nrc[str(nsub)][0],nrc[str(nsub)][1],
                               figsize=figsize[str(nsub)],sharex = 'col', 
                               num=fig.number) # existed, keep the same number
           fig.subplots_adjust(**spaces[str(nsub)]) # fine tune in dictionaries
       else: # single subplot
           fig,ax = P.subplots(nrc[str(nsub)][0],
                                figsize=figsize['1'],
                                num=fig.number) # existed, keep the same number
    else: # handle does not exist, make one
       if nrc[str(nsub)][1]: # not a single subplot
           fig,ax = P.subplots(nrc[str(nsub)][0],nrc[str(nsub)][1],
                               figsize=figsize[str(nsub)],sharex = 'col') # first creation
           fig.subplots_adjust(**spaces[str(nsub)]) # fine tune in dictionaries
       else: # single subplot
           fig,ax = P.subplots(nrc[str(nsub)][0],
                                figsize=figsize['1']) # first creation

    fig.canvas.manager.set_window_title('Fit parameters') # the title on the window bar
    fig.suptitle(labels['title']) # the sample title
    axout=[]
    axright = []
    if nsubplots>3: # two columns (nsubplots=6 for nsub=5)
        ax[-1,0].set_xlabel(labels['xlabel']) # set right xlabel
        ax[-1,1].set_xlabel(labels['xlabel']) # set left xlabel
        nrows = int(nsubplots/2) # (nsubplots=6 for nsub=5), 1, 2, 3
#        for k in range(0,nrows-1): 
#            ax[k,0].set_xticklabels([]) # no labels on all left xaxes but the last
#            ax[k,1].set_xticklabels([]) # no labels on all right xaxes but the last
        for k in range(nrows):
            axright.append(ax[k,1].twinx()) # creates replica with labels on right
            axright[k].set_ylabel(labels['ylabels'][nrows+k]) # right ylabels
            ax[k,0].set_ylabel(labels['ylabels'][k]) # left ylabels
            axright[k].tick_params(left=True,direction='in') # ticks in for right subplots
            ax[k,0].tick_params(top=True,right=True,direction='in') # ticks in for x axis, right subplots
            ax[k,1].tick_params(top=True,left=False,right=False,direction='in') # ticks in for x axis, right subplots
            ax[k,1].set_yticklabels([])
            axout.append(ax[k,0])    # first column
        for k in range(nrows):
            axout.append(axright[k])    # second column axout is a one dimensional list of axis   
    else: # one column
        ax[-1].set_xlabel(labels['xlabel']) # set xlabel
        for k in range(nsub-12): 
            ax[k].set_xticklabels([]) # no labels on all xaxes but the last
        for k in range(nsub):
            ylab = labels['ylabels'][k]
            if isinstance(ylab,str): # ylab = 1 for empty subplots
                ax[k].set_ylabel(ylab) # ylabels
                ax[k].tick_params(top=True,right=True,direction='in') # ticks in for right subplots
        axout = ax    # just one column
    return fig, axout


def set_bar(n,b):
    '''
    service to animate histograms
    e.g. in the fit tab

    extracted from matplotlib animate 
    histogram example
    '''
    from numpy import array, zeros, ones
    import matplotlib.path as path

    # get the corners of the rectangles for the histogram
    left = array(b[:-1])
    right = array(b[1:])
    bottom = zeros(len(left))
    top = bottom + n
    nrects = len(left)

    # here comes the tricky part -- we have to set up the vertex and path
    # codes arrays using moveto, lineto and closepoly

    # for each rect: 1 for the MOVETO, 3 for the LINETO, 1 for the
    # CLOSEPOLY; the vert for the closepoly is ignored but we still need
    # it to keep the codes aligned with the vertices
    nverts = nrects*(1 + 3 + 1)
    verts = zeros((nverts, 2))
    codes = ones(nverts, int) * path.Path.LINETO
    codes[0::5] = path.Path.MOVETO
    codes[4::5] = path.Path.CLOSEPOLY
    verts[0::5, 0] = left
    verts[0::5, 1] = bottom
    verts[1::5, 0] = left
    verts[1::5, 1] = top
    verts[2::5, 0] = right
    verts[2::5, 1] = top
    verts[3::5, 0] = right
    verts[3::5, 1] = bottom
    xlim = [left[0], right[-1]]
    return verts, codes, bottom, xlim

def set_fig(num,nrow,ncol,title,**kwargs): # unused? perhaps delete? check first 
    '''
    num is figure number (static, to keep the same window) 
    nrow, ncol number of subplots rows and columns
    kwargs is a dict of keys to pass to subplots as is
    initializes figures when they are first called 
    or after accidental killing
    '''
    import matplotlib.pyplot as P
    fig,ax = P.subplots(nrow, ncol, num = num, **kwargs)
    fig.canvas.manager.set_window_title(title)
    return fig, ax            
    
###############
# END OF PLOT #
###############

def rebin(x,y,strstp,pack,e=None):
    '''
    input:
        x is 1D intensive (time) 
        y [,e] are 1D, 2D or 3D intensive arrays to be rebinned
        pack > 1 is the rebinning factor, e.g it returns::
    
        xr = array([x[k*pack:k*(pack+1)].sum()/pack for k in range(int(floor((stop-start)/pack)))])
    
        strstp = [start,stop] is a list of slice indices 
       
        rebinning of x, y [,e] is done on the slice truncated to the approrpiate pack multiple, stopm
             x[start:stopm], y[start:stopm], [e[start:stopm]]      
    use either::

        xr,yr = rebin(x,y,strstp,pack)

    or::

       xr,yr,eyr = rebin(x,y,strstp,pack,ey) # the 5th is y error
    '''
    from numpy import floor, sqrt, zeros
    start,stop = strstp
#    print('tools rebin debug: start, stop, pack = {}, [], []'.format(start, stop, pack))
    m = int(floor((stop-start)/pack)) # length of rebinned xb
    mn = m*pack # length of x slice 
    xx =x[start:start+mn] # slice of the first 1d array
    xx = xx.reshape(m,pack) # temporaty 2d array
    xr = xx.sum(1)/pack # rebinned first ndarray
    if len(y.shape)==1:
        yb = zeros(m)
        yy = y[start:start+mn]  # slice row
        yy = yy.reshape(m,pack)  # temporaty 2d
        yr = yy.sum(1)/pack # rebinned row           
        if e is not None:
            ey = e[start:start+mn]   # slice row
            ey = ey.reshape(m,pack)  # temporaty 2d
            er = sqrt((ey**2).sum(1))/pack  # rebinned row - only good for ISIS 
    elif len(y.shape)==2:
        nruns = y.shape[0] # number of runs
        yr = zeros((nruns,m))
        if e is not None:
            er = zeros((nruns,m))
        for k in range(nruns): # each row is a run
            yy = y[k][start:start+mn]  # slice row
            yy = yy.reshape(m,pack)  # temporaty 2d
            yr[k] = yy.sum(1)/pack # rebinned row
            if e is not None:
                ey = e[k][start:start+mn]   # slice row
                ey = ey.reshape(m,pack)  # temporaty 2d
                er[k] = sqrt((ey**2).sum(1))/pack  # rebinned row        
    elif len(y.shape)==3:        
        ngroups,nruns = y.shape[0:2] # number of groups, runs
        yr = zeros((ngroups,nruns,m))
        
        if e is not None:
            er = zeros((ngroups,nruns,m))
        for k in range(ngroups): 
            for j in range(nruns):  
                yy = y[k][j][start:start+mn]  # slice row
                yy = yy.reshape(m,pack)  # temporaty 2d
                yr[k][j] = yy.sum(1)/pack # rebinned row
            if e is not None:
                ey = e[k][j][start:start+mn]   # slice row
                ey = ey.reshape(m,pack)  # temporaty 2d
                er[k][j] = sqrt((ey**2).sum(1))/pack  # rebinned row        
    if e is not None:
        return xr,yr,er
    else:
        return xr,yr

def rebin_decay(x,yf,yb,bf,bb,strstp,pack):
    '''
    input:
        x is 1D intensive (time)
        yf, yb 1D, 2D, 3D extensive arrays to be rebinned
        bf, bb are scalars or arrays (see musuite.single_for_back_counts and musuite.single_multigroup_for_back_counts)
        pack > 1 is the rebinning factor, e.g it returns::
    
        xr = array([x[k*pack:k*(pack+1)].sum()/pack for k in range(int(floor((stop-start)/pack)))])
        yr = array([y[k*pack:k*(pack+1)].sum() for k in range(int(floor((stop-start)/pack)))])
    
        strstp = [start,stop] is a list of slice indices 
       
        rebinning of x,y is done on the slice truncated to the approrpiate pack multiple, stopm
             x[start:stopm], y[start:stopm]        
    use::

        xr,yfr, ybr, bfr, bbr, yfmr, ybmr = rebin(x,yf,yb,bf,bb,yfm,ybm,strstp,pack)
    '''
    from numpy import floor, sqrt, exp, zeros, mean
    from mujpy.tools.tools import TauMu_mus

    start,stop = strstp
    m = int(floor((stop-start)/pack)) # length of rebinned xb
    mn = m*pack # length of x slice 
    xx =x[start:start+mn] # slice of the first 2D array
    xx = xx.reshape(m,pack) # temporaty 2d array
    xr = xx.sum(1)/pack # rebinned first ndarray
    bfr, bbr = bf*pack, bb*pack
    if len(yf.shape)==1:
        yfr = zeros(m)
        ybr = zeros(m) 
        yfr = yf[start:start+mn]  # slice row
        ybr = yb[start:start+mn]  # slice row
        yfr = yfr.reshape(m,pack)  # temporaty 2d
        ybr = ybr.reshape(m,pack)  # temporaty 2d
        yfr = yfr.sum(1) # rebinned row extensive          
        ybr = ybr.sum(1) # rebinned row extensive       
        yfmr, ybmr = mean((yfr-bfr)*exp(xr/TauMu_mus())), mean((ybr-bbr)*exp(xr/TauMu_mus()))   
    elif len(yf.shape)==2:
        nruns = yf.shape[0] # number of runs
        yfr = zeros((nruns,m))
        ybr = zeros((nruns,m)) 
        for k in range(nruns): # each row is a run, or a group
            yyf = yf[k][start:start+mn]  # slice row
            yyf = yyf.reshape(m,pack)  # temporaty 2d
            yfr[k] = yyf.sum(1) # rebinned row extesive
            yyb = yb[k][start:start+mn]  # slice row
            yyb = yyb.reshape(m,pack)  # temporaty 2d
            ybr[k] = yyb.sum(1) # rebinned row extesive
            bfr, bbr = bf[k]*pack, bb[k]*pack
            # print('tools,rebin_decay,debug: bfr {}, bbr {}'.format(bfr, bbr))
            yfmr, ybmr = mean((yfr[:][k]-bfr)*exp(xr/TauMu_mus())), mean((ybr[:][k]-bbr)*exp(xr/TauMu_mus()))
            
    elif len(yf.shape)==3:        # probably never used unless calib mode becomes a C2 case
        ngroups,nruns = yf.shape[0:2] # number of runs
        yfr = zeros((ngroups,nruns,m))
        ybr = zeros((nruns,m)) 
        for k in range(ngroups): 
            for j in range(nruns):  
                yyf = yf[k][j][start:start+mn]  # slice row
                yyf = yyf.reshape(m,pack)  # temporaty 2d
                yfr[k][j] = yyf.sum(1) # rebinned row extesive
                yyb = yb[k][j][start:start+mn]  # slice row
                yyb = yyb.reshape(m,pack)  # temporaty 2d
                ybr[k][j] = yyb.sum(1) # rebinned row extesive
                bfr, bbr = bf[k][j]*pack, bb[k][j]*pack
                yfmr, ybmr = mean((yfr[:][k][j]-bfr)*exp(xr/TauMu_mus())), mean((ybr[:][k][j]-bbr)*exp(xr/TauMu_mus()))
    return xr,yfr,ybr,bfr,bbr,yfmr,ybmr

def safetry(string):
    '''
    Used by muvalid
    '''
    from math import acos,asin,atan,atan2,ceil,cos,cosh,degrees,e,exp,floor,log,log10,pi,pow,radians,sin,sinh,sqrt,tan,tanh
    safe_list = ['a','acos', 'asin', 'atan', 'atan2', 'ceil', 'cos', 'cosh', 'degrees', 'e', 
                 'exp', 'floor', 'log', 'log10', 'pi', 'pow', 'radians', 'sin', 'sinh', 'sqrt', 'tan', 'tanh']
    # 	use the list to filter the local namespace
    a = 0.3
    safe_dict={}
    for k in safe_list:
        safe_dict[k]=locals().get(k)
    #    print(safe_dict[k])
    return eval(string,{"__builtins__":None},safe_dict)

def scanms(y,n):
    # produces guess for hifi t=0 bin, to be fed to a step fit function
    # check running average of (n bins,n skips,n bins) 
    # with two means m1,m2 and two variances s21,s22, against step pattern
    # compares m2-m1 with sqrt(s21+s22)
    from numpy import sqrt
    istart = []
    istop = []
    for k in range(y.shape[0]-n):
        m1,m2 = y[k:k+n].sum()/n, y[k+2*n:k+3*n].sum()/n
        s = sqrt(((y[k:k+n]-m1)**2).sum()/(n-1)+ ((y[k+2*n:k+3*n]-m2)**2).sum()/(n-1))
        if m2-m1>s:
            if not istart:
                istart = k+n
            elif not istop:
                istop = k+n
            elif istop == k+n-1:
                istop = k+n
        if istop and istart:
            if istop-istart == n:
                return istop
    return -1


def spec_prec(a):
    '''
    format specifier precision::

        0 for a > 1.0
        1 for 1.0 > a > 0.1
        2 for 0.1 > a > 0.01 etc.

    '''
    import numpy as np
    return int(abs(min(0.,np.floor(np.log10(a))))) 

def shorten(path,subpath):
    '''
    shortens path
    e.g. path, subpath = '/home/myname/myfolder', '/home/myname'
         shart = './myfolder' 
    '''
    short = path.split(subpath)
    if len(short)==2:
        short = '.'+short[1]
    return short

def exit_safe():
    '''
    opens an are you sure box?
    '''
    from tkinter.messagebox import askyesno
            
    answer = askyesno(title='Exit mujpy', message='Really quit?')
    return answer
    
def step(x,a,n,dn,b):
    from scipy.stats import norm
    # error function as step function for t=0 in HIFI
    return a+b*norm.cdf(x,n,dn)
       
def tlog_exists(path,run,ndigits):
    '''
    check if tlog exists under various known filenames types
    '''
    import os

    filename_psibulk = 'run_'+muzeropad(run,ndigits)+'.mon' # add definitions for e.g. filename_isis
    ok = os.path.exists(os.path.join(path,filename_psibulk)) # or os.path.exists(os.path.join(paths,filename_isis))
    return ok


def translate_multirun(functions_in,n_locals,kloc,nruns):
    '''
    functions_in  = [list of function strings], 
                    for the model components of a single-run model (obtained from get_functions_in)
                    where a "~","!" parameter dummy function has been redefined as 'p[k]'
                    and k is their dash index
    n_locals      = number of user_local parameters
    kloc          = index of first component first parameter 
                    in the single run model
                    (kloc-n_locals is the index of the first user_locals) 
    nruns         = number of runs in the suite
    functions_out = [list of lists of function strings, indices of component & parameter],
                    with translated minuit indices
                    outer list is components of the model, 
                    middle list is runs,
                    inner list is component parameter functions 
    used in int2_multirun_user_method_key and int2_multirun_grad_method_key
    '''
    # print('debug tools translate_multirun functions_in = {}'.format(functions_in))
    korig = kloc # minuit index index of first component first parameter in the single run model
    npar_run = n_locals # these will be the local parameters in each run, initialized to number of user_locals
    # extact this from functions: first scan a single run model to find how many parameters reference themeselves
    # the next loop is solely to determine npar_run = number of local parameters per run (korig is used but will be reset in the next loop)
    for funcs in functions_in: # list of functions for a component
        for func in funcs: # individual function for one parameter
            parloc = 'p['+str(korig)+']' #  original parameter
#            print('debug tools translate_multirun korig = {}, parloc = {}, func =  {}'.format(korig, parloc, func))
            if func.find(parloc)>=0: # if present this parameter references itself i.e. it is a "~" local parameter i.e. a free minuit parameter
             # (non-local are determines by user parameters)
#                print('debug tools translate_multirun found! korig = {}, func =  {}'.format(korig, func))
                korig += 1    # increment the minuit index of the single-run-model parameter
                npar_run +=1 # increment the number of local parameters per run
#    print('debug tools translate_multirun local parameters per run npar_run = {} minuit parameters =  {}'.format(npar_run, korig))

    fso = []
# next loop produces fso, appending a list funcs_out for each run, of lists func_out for each component, containing translated indices 
#                                                          from single-run-model to actual multirun model
    knew = kloc-n_locals # this is the index of the first user_local parameter
    for krun in range(nruns):
        funcs_out = []  # list of component lists of translated functions for a run
        korig = kloc # index for first run free parameter index in model, starts after user local replicas 
        knew += n_locals # minuit index includes run user local replicas, and is incremented at each run
        for funcs in functions_in: # funcs is list of functions for a component
            func_out = [] # list of translated functions for one component
            for k, func in enumerate(funcs): # individual function for parameter k
                parloc = 'p['+str(korig)+']' #  original parameter
                parnew = 'p['+str(knew)+']' # parameter for this run
                # print('debug tools translate_multirun knew = {}, parloc = {} parnew {}'.format(knew,parloc,parnew))
                if func.find(parloc)>=0: # if present
                    # print('debug tools translate_multirun func = {} becomes {}'.format(func,func.replace(parloc,parnew)))
                    func_out.append(func.replace(parloc,parnew)) # it is translated and appended
                    knew += 1 # increment the minuit index for the multirun model
                    korig +=1 # increment the single-run-model index
                else: # otherwise
                    func_out.append(func) # it is appended untranslated
                    # either way func_out[k] is appended  
                # now check for every parameter if contains user_local parameters                  
                for j in range(n_locals): # assign local index to user_local parameters
                    kuserorig = kloc-j-1 # index of user_local parameter for first run
                    parloc = 'p['+str(kuserorig)+']' # this is a user_local parameter
                    kusernew =  kloc+krun*npar_run-j-1 # index of user_local parameter for present run
                    parnew = 'p['+str(kusernew)+']' # this is the value for present run
                    if func_out[k].find(parloc)>=0: # if present
                        # print('debug tools translate_multirun func_out = {} becomes {}'.format(func_out[k],func_out[k].replace(parloc,parnew)))
                        func_out[k] = func_out[k].replace(parloc,parnew) # it is translated
                # at the end of this loop func.out is a list of func for the parameters of this component
            funcs_out.append(func_out) # adds this component to funcs_out, list of components in this run
        fso.append(funcs_out) # adds list of components for this run to list of runs
    # the next loop reshuffles fso to produce functions_out in the correct order model components, runs, component parameter
    functions_out = []
    ncomponents = len(fso[0]) # middle list is components
    for jcomp in range(ncomponents):
        frun = []
        for krun,run in enumerate(fso): # run is list of lists  for run krun
            frun.append(run[jcomp]) # run[jcomp] is the jcomp component functions for run krun
        # now frun is a list or runs for component jcomp    
        functions_out.append(frun) # now functions_out is a list of components, each a list or runs, each a list of parameter func
    return functions_out    
    
def get_indices(func):
    '''
    input 
      func is a user string function, e.g. 'p[0]*p[2]'
    output
      list of (string) indices found in the string, 
      in between 'p[' and ']', e.g. '0','2'
    '''
    from mujpy.tools.tools import findall
    return [func[i:j] for (i,j) in zip([k+1 for k in findall('[',func)],[l for l in findall(']',func)])]

def diffunc(func):
    '''
    input user function of the form 'p[0]*(1-p[1])'
          up to functions of three parameters (this could be easily extended)
    output a list of its derivatives [with respect to 'p[0]' and 'p[1]']
    and the list of their indices, [0,1]
    '''
    from mujpy.tools.tools import get_indices
    from sympy import symbols,diff,sympify,simplify
    from sympy import sin,cos,exp, sqrt,atan,pi
    # identify variables
    # first identify indices
    func = func.replace('abs','Abs').replace('arctan','atan')
#    indices =[func[i:j] for (i,j) in zip([k+1 for k in findall('[',func)],[l for l in findall(']',func)])]
    indices  = get_indices(func)
    ind = [int(k) for k in indices]
    if len(indices)==0: # no indices, func is the empty string
        return ['0'],ind
    elif len(indices)==1: # one index
        x = symbols('x')
        p0 = 'p['+indices[0]+']'
        fun = func.replace(p0,'x') # function of x
        f0 = str(diff(sympify(fun),x)).replace('x',p0).replace('Abs','abs').replace('atan','arctan')
        return [f0],ind
    elif len(indices)==2:  # two index
        x,y = symbols('x,y')
        p0,p1 = 'p['+indices[0]+']','p['+indices[1]+']'
        fun = func.replace(p0,'x').replace(p1,'y') # function of x,y
        f0 = str(sympify(diff(fun,x))).replace('x',p0).replace('y',p1).replace('Abs','abs').replace('atan','arctan')
        f1 = str(sympify(diff(fun,y))).replace('x',p0).replace('y',p1).replace('Abs','abs').replace('atan','arctan')
        return [f0,f1],ind
    elif len(indices)==3:   # three index
        x,y,z = symbols('x,y,z')
        p0,p1,p2 = 'p['+indices[0]+']','p['+indices[1]+']','p['+indices[2]+']'
        f0 = str(sympify(diff(fun,x))).replace('x',p0).replace('y',p1).replace('z',p2).replace('Abs','abs').replace('atan','arctan')
        f1 = str(sympify(diff(fun,y))).replace('x',p0).replace('y',p1).replace('z',p2).replace('Abs','abs').replace('atan','arctan')
        f2 = str(sympify(diff(fun,z))).replace('x',p0).replace('y',p1).replace('z',p2).replace('Abs','abs').replace('atan','arctan')
    # could be extended to four, five ...
        return [f0,f1,f2],ind
                
def translate(nint,lmin,function_in):
    '''
    input: 
        nint: dashbord index, 
        lmin: list of minuit indices replacement, one for each dashboard index, -1 is blanck
        function: single function string, of dashboard index nint, to be translated
    output: 
        function_out: single translated function
    Used in int2_method_key and min2int to replace parameter indices contained in function[nint] e.g.

    ::
 
       p[0]*2+p[3]

    by translate the internal parameter indices 0 and 3 (written according to the dashboard dict order)
    into the corresponding minuit parameter list indices, skipping shared parameters.

    e.g. if parameter 1 is shared with parameter 0, the minuit parameter index 3
    will be translated to 2  (skipping internal index 1)
    '''
    from copy import deepcopy
    from mujpy.tools.tools import findall
    # print(' nint = {}, lmin = {}\n{}'.format(nint,lmin,function_in))
    function_out = deepcopy(function_in)
    # search for integers between '[' and ']'
    start = [i+1 for i in  findall('[',function_out)]  
    # finds index of number after all occurencies of '['
    stop = [i for i in  findall(']',function_out)]
    # same for ']'
    nints = [function_out[i:j] for (i,j) in zip(start,stop)] 
    # this is a list of strings with the numbers to be replaced
    nmins = [lmin[int(function_out[i:j])] for (i,j) in zip(start,stop)]
    # replacements integers
    for lstr,m in zip(nints,nmins):
        function_out = function_out.replace(lstr,str(m))
    return function_out

def translate_nint(nint,lmin,function): # NOT USED any more?!!
    '''
    Used in int2_int and min2int to parse parameters contained in function[nint].value e.g.
    ::
 
       p[4]*2+p[7]

    and translate the internal parameter indices 4 and 7 (written according to the gui parameter list order)
    into the corresponding minuit parameter list indices, that skips shared and fixed parameters.

    e.g. if parameter 6 is shared with parameter 4 and parameter 2 is fixed, the minuit parameter indices
    will be 3 instead of 4 (skipping internal index 2) and 5 instead of 7 (skipping both 2 and 6)
    Returns lmin[nint]
    '''
    from mujpy.tools.tools import findall
    string = function[nint].value
    # search for integers between '[' and ']'
    start = [i+1 for i in  findall('[',string)]  
    # finds index of number after all occurencies of '['
    stop = [i for i in  findall(']',string)]
    # same for ']'
    nints = [string[i:j] for (i,j) in zip(start,stop)] 
    # this is a list of strings with the numbers
    nmins = [lmin[int(string[i:j])] for (i,j) in zip(start,stop)]
    return nmins

def value_error(value,error):
    '''
    value_error(v,e)
    returns a string of the format v(e) 
    '''
    from numpy import floor, log10, seterr
    eps = 1e-10 # minimum error
    if error>eps: # normal error
        exponent = int(floor(log10(error)))  
        most_significant = int(round(error/10**exponent))
        if most_significant>9:
            exponent += 1
            most_significant=1
        exponent = -exponent if exponent<0 else 0
        form = '"{:.'
        form += '{}'.format(exponent)
        form += 'f}({})".format(value,most_significant)'
    else:
        if abs(value)<eps:
            form = '"(0(0)"' # too small both
        else:
            form = '"{}(0)".format(value)' # too small error
    return eval(form)
    
    
def results():
    '''
    generate a notebook with some results
    '''
    import subprocess
    # write a python script
    script = '# Single Run Single Group Fit'
    script = script + '\n'#!/usr/bin/env python3'
    script = script + '\n# -*- coding: utf-8 -*-'
    script = script + '%matplotlib tk'
    script = script + '\n%cd /home/roberto.derenzi/git/mujpy/'
    script = script + '\nfrom mujpy.musuite import suite'	
    script = script + '\nimport json, re'
    script = script + '\nfrom os.path import isfile'
    script = script + '\nfrom mujpy.mufit import mufit'
    script = script + '\nfrom mujpy.mufitplot import mufitplot'
    script = script + "\njsonsuffix = '.json'\n"
    # notice: the new cell is produced by the \n at the end of the previous line followed by \n#
    script = script + '\n# Define log and data paths,   '
    script = script + '\n# detector grouping and its calibration  '
    script = script + '\nlogpath = "/home/roberto.derenzi/git/mujpy/log/"'
    script = script + '\ndatafile = "/home/roberto.derenzi/musrfit/MBT/gps/run_05_21/data/deltat_tdc_gps_0822.bin"'
    script = script + '\nrunlist = "822" # first run first'
    script = script + '\nmodelname = "mgml"'
    script = script + '\nversion = "1"'
    script = script + "\ngrp_calib = [{'forward':'3', 'backward':'4', 'alpha':1.13}]"
    script = script + "\ngroupcalibfile = '3-4.calib'"
    script = script + "\ninputsuitefile = 'input.suite'"
    script = script + "\ndashboard = modelname+'.'+re.search(r'\d+', runlist).group()+'.'+groupcalibfile[:groupcalibfile.index('.')]+'.'+version+jsonsuffix"
    script = script + '\nif not isfile(logpath+dashboard):' 
    script = script + "\n    print('Model definition dashboard file {} does not exist. Make one.'.format(logpath+dashboard))\n"
    script = script + "\n#  Can add 'scan': 'T' or 'B' for orderinng csv for increasing T, B, otherwise increasing nrun"
    script = script + "\ninput_suite = {'console':'print',"
    script = script + "\n               'datafile':datafile,"
    script = script + "\n               'logpath':logpath,"
    script = script + "\n               'runlist':runlist,"
    script = script + "\n               'groups calibration':groupcalibfile,"
    script = script + "\n               'offset':20"
    script = script + "\n              }  # 'console':logging, for output in Log Console, 'console':print, for output in notebook"
    script = script + '\nwith open(logpath+inputsuitefile,"w") as f:'
    script = script + "\n    json.dump(input_suite,f)"
    script = script + '\nwith open(logpath+groupcalibfile,"w") as f:'
    script = script + "\n    json.dump(grp_calib,f)"
    
    script = script + "\nthe_suite = suite(logpath+inputsuitefile,mplot=False) # the_suite implements the class suite according to input.suite\n"
    script = script + "\n# End of suite definition, this suiete is a single run"
    script = script + "\n# Now let's fit it according to the dashboard"
    script = script + "\nthe_fit = mufit(the_suite,logpath+dashboard)\n"
    script = script + "\n# Let's now plot the fit result"
    script = script + "\n# We plot the fit result, not the guess, over a single range"
    script = script + "\nfit_plot= mufitplot('0,20000,40',the_fit)#,guess=True) # '0,1000,4,24000,100' # two range plot"
    # save it in cache/notebook.py
    with open('/tmp/notebook.py',"w") as f:
        f.write(script)
    # compose notebook filename = 'xxyy.label.group.ipynb'
    bashCommand = 'p2j -o -t /home/roberto.derenzi/git/mujpy/getstarted/Delendo/xxyy.label.group.ipynb /tmp/notebook.py'
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    # issue os command 'p2j cache/notebook.py '+filename
    
def signif(x, p):
    '''
    write x with p significant digits
    '''
    from numpy import asarray,where,isfinite,abs,floor,log10,round
    x = asarray(x)
    x_positive = where(isfinite(x) & (x != 0), abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - floor(log10(x_positive)))
    return round(x * mags) / mags
        
