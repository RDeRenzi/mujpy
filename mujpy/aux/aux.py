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
    from mujpy.aux.aux import _available_components_, _errors_, _limits_
    
    available_components = _available_components_()
    dash_out = deepcopy(dash_in)   
    # these lists contain all parameter values in the dashboard, including their error steps and limits

    for k, component in enumerate(dash_out['model_guess']):
        steps = _errors_(component['name'],available_components)             
        limits = _limits_(component['name'],available_components) 
        for j,pardict in enumerate(component['pardicts']):
            dash_out['model_guess'][k]['pardicts'][j]['error'] = steps[j]                     
            dash_out['model_guess'][k]['pardicts'][j]['limits'] = limits[j]                         
    return dash_out
   
def _available_components_():
    '''
    returns a list of template dictionaries (one per fit component):

    Each dictionary contains 'name' and 'pardicts', 
    'pardicts' = list of parameter dictionaries, keys: 'name','error,'limits'

    ::  ({'name':'bl','pardicts':[{'name':'A','error':0.01,'limits'[0,0]},

                              {'name':'λ','error':0.01,'limits'[0,0]}}, 
                              ...)

    retreived magically from the mucomponents mumodel class.
    '''
    from mujpy.mucomponents.mucomponents import mumodel
    from iminuit import describe
    
    available_components = [] # generates the template of available components.
    for name in [module for module in dir(mumodel()) if module[0]!='_']: # magical extraction of component names
        pars = describe(mumodel.__dict__[name])[2:]            #  [12:] because the first two arguments are self, x
        _pars = [] 
        # print('pars are {}'.format(pars))
        for parname in pars:
        # parname, error, limits
        # In this template only
        #   {'name':'amplitude','error':0.01,'limits':[0, 0]}
        # parameter name will get a label later 
            error, limits = 0.01, [None, None] # defaults
            if parname == 'B' or parname == 'φ' or parname == 'Bd': error = 1.0
            if parname == 'β': error,limits = 0.05, [0, None]
            if parname == 'α': error,limits = 0.01, [0, None]
            # add here special cases for errors and limits, e.g. positive defined parameters
            _pars.append({'name':parname,'error':error,'limits':limits})
        available_components.append({'name':name,'pardicts':_pars})
    # [available_components[i]['name'] for i in range(len(available_components))] 
    # list of just mucomponents method names
    return available_components

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
    from mujpy.aux.aux import modelstrip
    
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

def _nparam(model):
    '''
    input: dashboard['model_guess']
    output: ntot, nmintot, nfree
    '''
    number_components = len(model)
    # print('_nparam aux debug: model {}'.format(model))
    ntot = sum([len(model[k]['pardicts']) 
                                 for k in range(number_components)]) # total number of component parameters
    flag = [pardict['flag'] for component in model for pardict in component['pardicts']]
    nmintot = ntot - sum([1 for k in range(ntot) if flag[k]=='=']) # ntot minus number of functions 
    nfree = nmintot - sum([1 for k in range(ntot) if flag[k]=='!']) # ntot minus number of fixed parameters 
    return ntot, nmintot, nfree
    
##################################################################
# int2min methods: generate guess values of minuit parameters
#  int2min : 
#  int2min_multigroup : assumes all parameters are in userpardicts
##################################################################

def int2min(model):
    '''
    input: 
        model 
            either dashboard["model_guess"] (after add_step_limits_to_model)
            or  dashboard["model_guess"] both lists of dicts
    output: a list of lists:  
        fitvalues: minuit parameter values, either guess of result
        fiterrors: their steps
        fitfixed: True/False for each
        fitlimits: [low, high] limits for each or [None,None]  
        parameter_name: name of parameter 'x_label' for each parameter

    To be used e.g. as:
            self.lastfit = Minuit(self._the_model_._chisquare_,
                             name=parameter_name,
                             *fitvalues)                                        
            self.lastfit.errors = fiterrors
            self.lastfit.limits = fitlimits
            self.lastfit.fixed = fitfixed
            self.lastfit.migrad()
    '''
    from mujpy.aux.aux import _nparam

    dum, ntot, dum  = _nparam(model)
    
    #####################################################
    # the following variables contain the same as input #
    # parameters to iMinuit, removing '='s (functions)  #
    #####################################################
                                                        
    fitval, fiterr, fitfix, fitlim = [], [], [], []           
    parameter_name = []

    nint = -1 # initialize
    for k in range(len(model)):  # scan the model components
        label = model[k]['label']
        for j, pardict in enumerate(model[k]['pardicts']):  # list of dictionaries
            nint += 1  # internal parameter incremented always   
            # print(pardict)
            if pardict['flag'] != '=': #  skip functions, only new minuit parameters
                fitval.append(float(pardict['value']))
                parameter_name.append(pardict['name']+'_'+label) 
                fiterr.append(float(pardict['error']))
                fitlim.append(pardict['limits'])
            if pardict['flag'] == '~':
                 fitfix.append(False)
            elif pardict['flag'] == '!':
                fitfix.append(True)
    # self.console('fitval = {}\nfiterr = {}\nfitfix = {}\nfitlim = {}\ncomp name = {},\npar name = {} '.format(fitval,fiterr,fitfix,fitlim,component_name,parameter_name)) 

    return fitval, fiterr, fitfix, fitlim, parameter_name

def int2min_multigroup(pardicts):
    '''
    input: 
        pardicts 
            either dashboard["userpardicts_guess"] if guess = True
            or  dashboard["userpardicts_result"] if gues = False
    output: a list of lists:  
        fitvalues: minuit parameter values, either guess of result
        fiterrors: their steps
        fitfixed: True/False for each
        fitlimits: [low, high] limits for each or [None,None]  
        parameter_name: name of parameter 'x_label' for each parameter
    this works for A2 single fit, multigroup with userpardicts parameters = Minuit parameters
    '''
    # pardicts = model["userpardicts_guess"]
    ntot  = len(pardicts)
    
    #####################################################
    # the following variables contain the same as input #
    # parameters to iMinuit, removing '='s (functions)  #
    #####################################################
                                                        
    fitval, fiterr, fitfix, fitlim = [], [], [], []           
    parameter_name = []
    for pardict in pardicts:  # scan the model components
        errstd = 'error' if 'error' in pardict.keys() else 'std'
        fitval.append(float(pardict['value']))
        parameter_name.append(pardict['name']) 
        fiterr.append(float(pardict[errstd]))
        if 'flag' in pardict.keys():
            if pardict['flag'] == '!':
                fitfix.append(True)
            elif pardict['flag'] == '~':
                fitfix.append(False)
            else:
                return False,_,_,_,_
            fitlim.append(pardict['limits'])
        # self.console('fitval = {}\nfiterr = {}\nfitfix = {}\nfitlim = {}\ncomp name = {},\npar name = {} '.format(fitval,fiterr,fitfix,fitlim,component_name,parameter_name)) 
    return fitval, fiterr, fitfix, fitlim, parameter_name
    
##################################
# method and key methods: provide component methods 
#                           and parameter key for eval(key) in _add_
#   int2_method_key :                single run single group 
#   int2_multigroup_method_key :     single run multi group


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
    This function applies aux.translate to the parameter numbers in formulas:
    dashboard "function" is written in terms of the internal parameter index,
    while Minuit parameter index skips shared or formula-determined ('=') parameters  
    '''
    from mujpy.aux.aux import translate

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
                key_as_lambda = eval('lambda p:'+string) # NEW! speedup
                keys.append(key_as_lambda) # the function key in keys will be evaluated, key(p), inside mucomponents
                # isminuit.append(False)
                lmin.append(0)
            else:# flag[j] == '~' or flag[j] == '!'
                nmin += 1
                key_as_lambda = eval('lambda p:'+'p['+str(nmin)+']') # NEW! speedup
                keys.append(key_as_lambda) # the function key in keys will be evaluated, key(p), inside mucomponents
                lmin.append(nmin) # 
                # isminuit.append(True)
        # print('int2_method aux debug: bndmthd = {}, keys = {}'.format(bndmthd,keys))
        method_key.append([bndmthd,keys]) 
    return method_key

def int2_multigroup_method_key(dashboard,the_model,guess=True):
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

    This function does not need userpardicts and aux.translate 
    since the correspondence with Minuit parameters
    is given directly either by "function" or by "function_multi"
    '''
    from mujpy.aux.aux import multigroup_in_components, fstack

    model = dashboard['model_guess']  # guess surely exists
    # these are the only Minuit parameters [p[k] for k in range (nuser)]
    ntot = sum([len(model[k]['pardicts']) for k in range(len(model))])
    method_key = []
    pardicts = [pardict for component in model for pardict in component['pardicts']]
    mask_function_multi = multigroup_in_components(dashboard)
    # print('int2_multigroup_method_key aux debug: index function_multi {}\npardicts = {}'.format(mask_function_multi,pardicts))
    if sum(mask_function_multi):
        ngroups = len(pardicts[mask_function_multi.index(1)]["function_multi"])
    else:
        return []
    nint = -1 # initialize the index of the dashboard component parameters
    # p = [1.13,1.05,0.25,0.3,0.8,700,35,125,3.3,680,0.1] # debug delete
    for component in model:  # scan the model components
        name = component['name']
        keys = []
        # print('name = {}, model = {}'.format(name,self._the_model_))
        bndmthd = lambda x,*pars : fstack(the_model.__getattribute__(name),x,*pars)
                            # this is the method to calculate a component, to set alpha, dalpha apart
        nint0 = nint
        for l in range(ngroups):
            key = []  
            nint = nint0
            for j,pardict in enumerate(component['pardicts']): 
                nint += 1  # internal parameter index incremented always 
                if mask_function_multi[nint]>0:
                    key_as_lambda = eval('lambda p:'+pardict["function_multi"][l]) # NEW! speedup
                else:                
                    key_as_lambda = eval('lambda p:'+pardict["function"]) # NEW! speedup
                # print('aux int2_multigroup_method_key debug: key_as_lambda(p) = {} **delete also p!'.format(key_as_lambda(p)))
                key.append(key_as_lambda) # the function key will be evaluated, key(p), inside mucomponents
            keys.append(key)
        # print('int2_method aux debug: bndmthd = {}, keys = {}'.format(bndmthd,keys))
        method_key.append([bndmthd,keys]) # vectorialized method, with keys 
        # keys = [[strp0g0, strp1g0,...],[strp0g1, strp1g1, ..],[strp0g2, strp1g2,...]..]
        # pars = [[p0g0, p1g0, ...],[p0g1, p1g1, ..],[p0g2, p1g2,...]..]
    return method_key
    
def fstack(npfunc,x,*pars):
    '''
    vectorialize npfunc
    input: 
        npfunc numpy function with input (x,*argv)
        x time
        *argv is variable number of lists of parameters, list len is the output_function_array.shape[0]
    output:
        output_function_array
            stacks vertically n replica of npfunc distributing parameters as in
            (x, *argv[i]) for each i-th replica 
    '''
    from numpy import vstack
    for k,par in enumerate(pars):
        if k:
            f = vstack((f,npfunc(x,*par)))
        else:
            f = npfunc(x,*par)
    return f
    
def int2_calib_method_key(dashboard,the_model):
    '''
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
     
    This function applies aux.translate to the parameter numbers in formulas
    since on the dash each parameter of each component gets an internal number,
    but alpha is popped and shared or formula-determined ('=') ones are not minuit parameters  
    '''
    from mujpy.aux.aux import translate

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
     
    This function applies aux.translate to the parameter numbers in formulas
    since on the dash each parameter of each component gets an internal number,
    but alpha is popped and shared or formula-determined ('=') ones are not minuit parameters  
    '''
    from mujpy.aux.aux import translate
    print('int2_calib_multigroup_method_key aux debug: copy of non multigroup, adapt!')
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
    from mujpy.aux.aux import translate

    names, values_out, p, errors_out, e = [], [], [], [], []
    nint = -1 # initialize
    nmin = -1
    lmin = []
    flag = [pardict['flag'] for component in model_guess for pardict in component['pardicts']]
    flag = [pardict['flag'] for component in model_guess for pardict in component['pardicts']]
    for k,component in enumerate(model_guess):  # scan the model
        component_name = component['name']
        name, value, error = [], [], []
        label = model_guess[k]['label']
        
        for j,pardict in enumerate(model_guess[k]['pardicts']): # list of dictionaries, par is a dictionary
            nint += 1  # internal parameter incremented always
            if j==0:
                name.append('{}: {}_{}'.format(component_name,pardict['name'],label))
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

def min2int_multigroup(dashboard,p,e):
    '''
    input:
        userpardicts_guess from dashboard (each dict corresponds to a Minuit parameter)
    output: for all parameters
        names list of lists of parameter names
    used only in summary_global
    '''
    # 
    # initialize
    #
    from mujpy.aux.aux import multigroup_in_components
    # print('min2int_multigroup in aux debug: dash {}'.format(dashboard))
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
                # print('aux min2int_multigroup debug: nint = {} name = {} par = {}, epar = {}'.format(nint0,name[-1],par[-1],epar[-1]))
            pars.append(par) # middel list, model
            names.append(name)
            epars.append(epar)
        namesg.append(names)
        parsg.append(pars)
        eparsg.append(epars) 
    return namesg, parsg, eparsg  # list of parameter values 

def print_components(names,values,errors):
	'''
	input: for a component
		parameter names 
		parameter values 
		parameter errors 
	output:
	    string to print, e.g.
	    "bl.A_fast 0.123(4) bl.λ_fast 12.3(4) bl.σ_fast 0(0)"
	'''
	from mujpy.aux.aux import value_error
	out = [' '.join([names[k],'=',value_error(values[k],errors[k])]) for k in range(len(names))]
	return " ".join(out)
	
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
    from mujpy.aux.aux import filter
    from numpy import pi, cos, vstack, fft, delete
    ydim, tdim = len(y.shape), len(t.shape)
    # print('aux mixer debug 1: y t shape {}, {}'.format(y.shape,t.shape))
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
    yf = fft.irfft(filter(t,fft.rfft(2*y*cos(2*pi*f0*t),n=n+1),f0),n=2*n)
    # now delete padded zeros 
    mindex = range(n,2*n)
    yf =delete(yf,mindex,-1)
    # print('aux mixer debug 3: yf shape {}'.format(yf.shape))
    return yf
    
def filter(t,fy,f0):
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
    # print('aux filter debug 2: fy {},mask {} shape'.format(fy.shape,mask.shape))    
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

    #print('multigroup_in_components aux debug: model_guess len {}'.format(len(dashboard["model_guess"])))
    #print('multigroup_in_components aux debug: pardicts len {}'.format(len(dashboard["model_guess"][0]['pardicts'])))
    #print('multigroup_in_components aux debug: pardict.keys len {}'.format(len(dashboard["model_guess"][0]['pardicts'][0].keys())))

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

def tilde_in_components(dashboard):
    '''
    input full dashboard
    output list of indices where  
        "model_guess" (list) 
        contains at least one component (dict) 
        whose "pardicts" (list) 
        contains a parameter dict 
        with at least one "flag":"~" key
    Indices refers to a single list 
    of dashboard component parameters
    The model name (e.g. mgbl) dictates the order
    '''
    component_flags = [pardict["flag"] for component in dashboard["model_guess"] 
                                       for pardict in component["pardicts"]]
    indices_tilde = [i for i, x in enumerate(component_flags) if x == "~"]

    return indices_tilde
    
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

def create_model(dashboard):
    '''
    create_model('daml') # adds e.g. the two component 'da' 'ml' model
    this method 
    does not check syntax (prechecked by checkvalidmodel)
    separates nexternals number from model name (e.g. '3mgml' -> 'mgml', 3)
    starts switchyard for A1,A1,B1, B2, C1, C2 fits
    adds a model of components selected from the available_component tuple of  
    directories
    with zeroed values, stepbounds from available_components, flags set to '~' and empty functions
    '''
    import string
    from mujpy.aux.aux import modelstrip, addcomponent


    components
    self.model_guess = [] # start from empty model
    for k,component in enumerate(components):
        label = string.ascii_lowercase[k] # was uppercase[k]
        if not addcomponent(component,dashboard):
            return False
        # self.console('create model added {}'.format(component+label))
    return True

def addcomponent(name,label):
    '''
    addcomponent('ml') # adds e.g. a mu precessing, lorentzian decay, component
    this method adds a component selected from _available_components_(), tuple of directories
    with zeroed values, error and limits from available_components, 
    flags set to '~' and empty functions
    [plan also addgroupcomponents and addruncomponents (for A2, B2, C1, C2)]
    '''
    from copy import deepcopy
    from mujpy.aux.aux import _available_components_
    available_components =_available_components_() # creates list automagically from mucomponents
    component_names = [available_components[i]['name'] 
                            for i in range(len(available_components))]
    if name in component_names:
        k = component_names.index(name)
        npar = len(available_components[k]['pardicts']) # number of pars
        pars = deepcopy(self.available_components[k]['pardicts']) # list of dicts for 
        # parameters, {'name':'asymmetry','error':0.01,'limits':[0, 0]}

        # now remove parameter name degeneracy                   
        for j, par in enumerate(pars):
            pars[j]['name'] = par['name']+label
            pars[j].update({'value':0.0})
            pars[j].update({'flag':'~'})
            pars[j].update({'function':''}) # adds these three keys to each pars dict
            # they serve to collect values in mugui
        # self.model_guess.append()
        return {'name':name,'pardicts':pars}, None # OK code, no message
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
        
def derange(string,vmax,int_or_float='int'):
    '''
    derange(string) 

    reads string 
    assuming 2, 3, 4 or 5 csv or space separated values, either int (default) or floats (specify type='float')::

        5: start, stop, packe, last, packl
        4: start, stop, last, packl (packe is 1)
        3: start, stop, pack
        2: start, stop

    returns 2, 3, 4 or 5 floats or int, or 
    two negative values, if fails either::

        validity check (stop>start and bin <stop-start) 
        check for sanity (last value less then a vmax)

    '''
#    print('In derange')
    try:  
        try:
            if int_or_float=='float':
                values = [float(x) for x in string.split(',')]
            else:
                values = [int(x) for x in string.split(',')]
        except:
            print(string)
        if len(values)==5: # start, stop, packe, last, packl
            if values[3]<values[1] or values[4]>values[3]-values[1] or values[2]>values[1]-values[0] or values[1]<values[0] or sum(n<0 for n in values)>0 or values[3]>vmax:
                
                return -5,-5
            return values[0],values[1],values[2],values[3],values[4] # start, stop, packe, last, packl
        elif len(values)==4: # start, stop, last, packl
            if values[2]<values[1] or values[3]>values[2]-values[1] or values[1]<values[0] or sum(n<0 for n in values)>0 or values[2]>vmax:
                return -4,-4
            return values[0],values[1],1,values[2],values[3] # start, stop, packe, last, packl
        elif len(values)==3: # start, stop, packe
            if values[2]>values[1] or values[1]<values[0] or sum(n<0 for n in values)>0 or values[1]>vmax:
                return -3,-3
            return values[0],values[1],values[2] # start, stop, pack
        elif len(values)==2: # start, stop
            if values[1]<values[0] or sum(n<0 for n in values)>0 or values[1]>vmax:
                return -2,-2
            return values[0],values[1] # start, stop
    except:
        # print(string,vmax,int_or_float)
        return -10,-10

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
    s = []
    try:
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

def get_grouping(groupcsv):
    """
    name = 'forward' or 'backward'

    * grouping(name) is an np.array wth detector indices
    * group.value[k] for k=0,1 is a shorthand csv like '1:3,5' or '1,3,5' etc.
    * index is present mugui.mainwindow.selected_index
    * out is mugui._output_ for error messages

    returns

    * grouping, group, index
         
    group and index are changed only in case of errors
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
    if filespec=='bin' or filespec=='mdu':
        TsTc, eTsTc = the_run.get_temperatures_vector(), the_run.get_devTemperatures_vector()
        n1,n2 = spec_prec(eTsTc[0]),spec_prec(eTsTc[1]) # calculates format specifier precision
        form = '{} {:.'
        form += '{}'.format(n1)
        form += 'f}  {:.'
        form += '{}'.format(n1)
        form += 'f}  {:.'
        form += '{}'.format(n2)
        form += 'f}  {:.'
        form += '{}'.format(n2)
        form += 'f} {}' #".format(value,most_significant)'
        return form.format(nrun, TsTc[0],eTsTc[0],TsTc[1],eTsTc[1], Bstr[:Bstr.find('G')])
    elif filespec=='nxs':
        Ts = the_run.get_temperatures_vector()
        n1 = '1'       
        form = '{} {:.'
        form += '{}'.format(n1)
        form += 'f} {}' #".format(value)
        return form.format(nrun, Ts[0], Bstr[:Bstr.find('G')])

def minparam2_csv(dashboard,values_in,errors_in):
    '''
    input:
        dashboard dashboard[model_guess"], for single group or None, for multi group
        multigroup True/False
        Minuit values
        Minuit errors
    output:
        cvs partial row with parameters and errors
    '''
    from mujpy.aux.aux import min2int, spec_prec

    # Minuit and user parameters coincide
    (_, values, errors) = (min2int(dashboard,values_in,errors_in) if dashboard else
                                                        (None, [values_in], [errors_in]))
        
        # from minuit parameters to component parameters
    # output is lists (components) of lists (parameters) 
    row = ''
    for parvalues, parerrors in zip(values,errors): 
        for parvalue,parerror in zip(parvalues,parerrors):
            n1 = spec_prec(parerror) # calculates format specifier precision
            form = ' {:.'
            form += '{}'.format(n1)
            form += 'f}  {:.'
            form += '{}'.format(n1)
            form += 'f}'
            row += form.format(parvalue,parerror)
    return row

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
    form = ' {:.'
    form += '{}'.format(n1)
    form += 'f}  {:.'
    form += '{}'.format(n1)
    form += 'f}  {:.'
    form += '{}'.format(n1)
    form += 'f}' # ' {} {}'
    row = form.format(chi2,chi2-lowchi2,hichi2-chi2)
    for group in groups:
        row += ' {}'.format(group["alpha"])
    row += ' {} {}'.format(offset,strftime("%d.%b.%H:%M:%S", localtime()))
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
    from mujpy.aux.aux import get_title
    import os
    from datetime import datetime

    nrun = int(row.split(" ")[0])
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S") 
    if scan==None:  # order by nrun, first item in csv
        csv_index = 0
    elif scan=='T': # order by T, 4th or 2nd item in csv
        csv_index = 3 if filespec == 'bin' or filespec == 'mdu' else 1
    else:           # order by B, 6th or 4th item in csv
        csv_index = 5 if filespec == 'bin' or filespec == 'mdu' else 3
    rowvalue = float(row.split(" ")[csv_index]) # also nrun is transformed into float 

    if os.path.isfile(file_csv):
        try: # the file exists
            lineout = [] # is equivalent to False
            with open(file_csv,'r') as f_in:
                notexistent = True # the line dows not exist
                for nline,line in enumerate(f_in.readlines()):
                    if nline==0:
                        if header!=line: # different headers, substitute present one
                            raise # exits this try 
                        else:
                            lineout.append(header)
                    elif float(line.split(" ")[csv_index]) < rowvalue: # append line first
                        lineout.append(line)
                    elif float(line.split(" ")[csv_index]) == rowvalue: # substitute an existing fit
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
        temp = float(tstr[:tstr.index('K')-1])
        title.append('{:.1f}K'.format(temp))
    if not nofield:
        title.append('{:.0f}mT'.format(float(run.get_field())/10))
    return ' '.join(title)    

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
    here = os.getcwd()
    os.chdir(path)
    tkinter.Tk().withdraw() # Close the root window
    spc, spcdef = '.'+spec,'*.'+spec
    in_path = filedialog.askopenfilename(initialdir = path,filetypes=((spc,spcdef),('all','*.*')))
    os.chdir(here)
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

def plotile(x,xdim=0,offset=0):
    '''
    Produces a tiled plot, in the sense of np.tile e.g.

    ::

        x.shape = (1,1000) 
        y.shape = (4,1000)
        xt = plotile(x,4)
        yt = plotile(y,offset=0.1) 

    '''
    # x is an array(x.shape[0],x.shape[1])
    # xoffset is a step offset
    # xdim = x.shape[0] if xdim == 0 else xdim
    # each row is shifted by xoffset*n, where n is the index of the row  
    # 
    # 
    from copy import deepcopy
    from numpy import tile, arange
    xt = deepcopy(x)
    if xdim != 0: # x is a 1D array, must be tiled to xdim
        xt = tile(xt,(int(xdim),1))
    if offset != 0:
        xt += tile(offset*arange(xt.shape[0]),(x.shape[1],1)).transpose()
    return xt

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
    from numpy import floor, sqrt, empty

    start,stop = strstp
    m = int(floor((stop-start)/pack)) # length of rebinned xb
    mn = m*pack # length of x slice 
    xx =x[start:start+mn] # slice of the first 1d array
    xx = xx.reshape(m,pack) # temporaty 2d array
    xr = xx.sum(1)/pack # rebinned first ndarray
    if len(y.shape)==1:
        yb = empty(m)
        yy = y[start:start+mn]  # slice row
        yy = yy.reshape(m,pack)  # temporaty 2d
        yr = yy.sum(1)/pack # rebinned row           
        if e is not None:
            ey = e[start:start+mn]   # slice row
            ey = ey.reshape(m,pack)  # temporaty 2d
            er = sqrt((ey**2).sum(1))/pack  # rebinned row - only good for ISIS 
    elif len(y.shape)==2:
        nruns = y.shape[0] # number of runs
        yr = empty((nruns,m))
        if e is not None:
            er = empty((nruns,m))
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
        yr = empty((ngroups,nruns,m))
        if e is not None:
            er = empty((ngroups,nruns,m))
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
    from numpy import floor, sqrt, exp, empty, mean
    from mujpy.aux.aux import TauMu_mus

    start,stop = strstp
    m = int(floor((stop-start)/pack)) # length of rebinned xb
    mn = m*pack # length of x slice 
    xx =x[start:start+mn] # slice of the first 2D array
    xx = xx.reshape(m,pack) # temporaty 2d array
    xr = xx.sum(1)/pack # rebinned first ndarray
    bfr, bbr = bf*pack, bb*pack
    if len(yf.shape)==1:
        yfr = empty(m)
        ybr = empty(m) 
        yfr = yf[start:start+mn]  # slice row
        ybr = yb[start:start+mn]  # slice row
        yfr = yfr.reshape(m,pack)  # temporaty 2d
        ybr = ybr.reshape(m,pack)  # temporaty 2d
        yfr = yfr.sum(1) # rebinned row extensive          
        ybr = ybr.sum(1) # rebinned row extensive       
        yfmr, ybmr = mean((yfr-bfr)*exp(xr/TauMu_mus())), mean((ybr-bbr)*exp(xr/TauMu_mus()))   
    elif len(yf.shape)==2:
        nruns = yf.shape[0] # number of runs
        yfr = empty((nruns,m))
        ybr = empty((nruns,m)) 
        for k in range(nruns): # each row is a run, or a group
            yyf = yf[k][start:start+mn]  # slice row
            yyf = yyf.reshape(m,pack)  # temporaty 2d
            yfr[k] = yyf.sum(1) # rebinned row extesive
            yyb = yb[k][start:start+mn]  # slice row
            yyb = yyb.reshape(m,pack)  # temporaty 2d
            ybr[k] = yyb.sum(1) # rebinned row extesive
            bfr, bbr = bf[k]*pack, bb[k]*pack
            # print('aux,rebin_decay,debug: bfr {}, bbr {}'.format(bfr, bbr))
            yfmr, ybmr = mean((yfr[:][k]-bfr)*exp(xr/TauMu_mus())), mean((ybr[:][k]-bbr)*exp(xr/TauMu_mus()))
            
    elif len(yf.shape)==3:        # probably never used unless calib mode becomes a C2 case
        ngroups,nruns = yf.shape[0:2] # number of runs
        yfr = empty((ngroups,nruns,m))
        ybr = empty((nruns,m)) 
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
