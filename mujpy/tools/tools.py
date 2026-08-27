"""
########################
# FFT AUTO PHASE METHODS
########################
"""
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

"""
##############
# MUSUITE AUX
##############
"""

def short_path(path,startup_path):
    '''
    try to shorten path keeping main log info

    input path, startup_path
        tries to remove or shorted common beginning
        used in musuite.py mufit.py 
    '''

    from pathlib import Path
    home = Path.home().as_posix()
    path = path[2:] if path[0:2]=='..' else path
    if path == startup_path:
        return './'
    else:
        # print('tools short_path path {}\n startup_path {}'.format(path,startup_path))
        index = next((i for i, (char1, char2) in enumerate(zip(path, startup_path)) if char1 != char2), None)

        index = len(startup_path) if index == None else index
        k = index -1 if path[index-1] == '/' else index
        # print('index {}, k {}'.format(index,k))
        short_path = '~/'+path[index:] if path[0:k] == home else '.'+path[k:]
        return short_path

"""
##############
# MUFIT AUX
##############
"""

def TauMu_mus():
    """
    muon mean lifetime in microsecond

    from Particle Data Group 2017 
    (not present in scipy.constants)
    """

    return 2.1969811 
    
def _errors_(component):
    """
    suggested steps for each mumodel component parameter

    inputs: one legal mucomponent name contained in _available_components_()
    output: a list of errors (steps), one for each parameter of this component
        used in tools.add_step_limits
    """
    from mujpy.tools.tools import _available_components_ 
    available_components = _available_components_()
    #print(component,available_components)
    k = [item['name'] for item in available_components].index(component)
    return [pardict["error"] for pardict in available_components[k]['pardicts']] 

def _limits_(component):
    """
    suggested limits for each mumodel component parameter

    inputs: one legal mucomponent name contained in _available_components_()
    output: a list of lists of limits (low, high), one for each parameter of this component
        used in tools.add_step_limits
    """
    from mujpy.tools.tools import _available_components_ 
    available_components = _available_components_()
    name = [item['name'] for item in available_components].index(component)
    return [pardict["limits"] for pardict in available_components[name]['pardicts']] 

def _pospar_(component):
    """
    pospar, True or False for each mumodel component parameter

    inputs: one legal mucomponent name contained in _available_components_()
    output: a list of True or False for each parameter of this component
        used in tools.add_step_limits
    """
    from mujpy.tools.tools import _available_components_ 
    available_components = _available_components_()
    #print(component,available_components)
    name = [item['name'] for item in available_components].index(component)
    pardicts = []
    for pardict in available_components[name]['pardicts']:
        p = True if 'pospar' in pardict else False
        pardicts.append(p)
    return pardicts

def add_step_limits(model_in):
    """
    add error, limits [and pospar] to pardicts of model_in, a  list of components

    input:  dashboard['model_guess'] list, or any equivalent list
            its pardicts already contain keys 'name','value','flag','function'
            (sequential fits A1, A20, B1, B20 assumed)
    output: model is a deepcopy including 'error', 'limits'[, 'positive_parity']
            if model_in passes invalid_err_lim test.
        used in mudashed
    """

    from copy import deepcopy
    from mujpy.tools.tools import _errors_, _limits_,_pospar_,invalid_err_lim
  
    model = deepcopy(model_in)   
    # these lists contain all parameter values in the dashboard, including their error steps and limits

    for component in model:
        #values = _values_(component['name']) 
        steps = _errors_(component['name'])             
        limits = _limits_(component['name'])
        pospar = _pospar_(component['name'])
        for j,pardict in enumerate(component['pardicts']):
            pardict['error'] = steps[j]                  
            pardict['limits'] = limits[j]
            if pospar[j]: # only True is added
                pardict['positive_parity'] = pospar[j]
    return model
   
def _available_components_():
    """
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
        used in mufit and mudashed
    """

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
            if parname in positive_parity:
                _pars.append({'name':parname,'error':error,'limits':limits,'positive_parity':True})
            else:
                _pars.append({'name':parname,'error':error,'limits':limits})
        available_components.append({'name':name,'pardicts':_pars,'tip':tip})
    # [available_components[i]['name'] for i in range(len(available_components))] 
    # list of just mucomponents method names
    return available_components

def check_function(dashboard,groups):
    """
    checks function syntax for max 100 parameters [must be extended]

        input
            dashboard (full dict)
            groups (suite.groups)
        output 
            kc, kp, kg indices of failing component, parameter, group 
    checks that
        len of any 'function_multi' list is same as len(groups)
        only model contains either 'function' or 'function_multi' (pops the empty one)
        each string is executable
        called in mufit._dash_load_
    """

    from numpy import linspace, cos, sin, tan, sinh, cosh, tanh, log, pi, exp, sqrt, real, abs, arctan
    model = dashboard["model_guess"]
    model_pardicts = [pardict for component in model for pardict in component['pardicts'] if pardict['flag'] in ['~','!']]
    n_par = len(dashboard['globpardicts_guess']) if 'globpardicts_guess' in dashboard else len(model_pardicts)
    #           min_pars used in 'function'/'function_multi' of model are either len(dashboard['globpardicts_guess']
    #                                                                     or len model_pardicts
    p = linspace(.1,2.0,n_par) # mock p values for evaluated 'function'/'function_multi' strings

    # local fits min_pars coincide with non '=' model parameters 
    if 'globpardicts_guess' in dashboard:
    # global fits minuit parameters coincide with globpardicts, 
        for kp, pardict in enumerate(dashboard['globpardicts_guess']):
            if pardict['flag'] == '=': 
                print("++++ Wrong json syntax, only '~','!','#' flags allowed in 'globpardicts_guess'")
                return -1, kp, 0

        for kc, component in enumerate(model):
            for kp,pardict in enumerate(component['pardicts']):
                if pardict['flag'] != '=':
                    print('++++ Wrong json syntax, only '+"'='"+' flag allowed in '+"'model_guess'")
                    return kc, kp, 0
                if 'function' in pardict.keys() and len(pardict['function']):
                    kgroup = 0 # default
                    if 'function_multi' in pardict.keys() and len(pardict['function_multi'][0]):
                        print("++++ Wrong json syntax, 'function','function_multi' mutually exclusive in 'model_guess'")
                        return kc, kp, kgroup
                    if 'function_multi' in pardict.keys(): # zero len
                        pardict.pop('function_multi')
                    string = pardict['function']
                    # print('user function {}, kc,kp,kg={},{},{}'.format(string, kc,kp,kgroup))
                    eval(string)
                elif 'function_multi' in pardict.keys() and len(pardict['function_multi'][0]):
                    if 'function' in pardict.keys(): # zero len
                        pardict.pop('function')
                    if len(pardict['function_multi']) != len(groups):
                        print("++++ Wrong json syntax, need as many 'function_multi' strings as groups in 'model_guess'")
                        return kc, kp, len(groups)
                    for kgroup in range(len(groups)): # loop groups
                        string = pardict['function_multi'][kgroup]
                        # print('function multi {}, kc,kp,kg={},{},{}'.format(string, kc,kp,kgroup)) 
                        eval(string)
                else:
                    print("++++ Wrong json syntax, missing 'function' or 'function_multi' strings in 'model_guess'")
                    return kc, kp, 0

    else: # local fit
        kgroup = 0 # default
        for kc,component in enumerate(model):
            for kp,pardict in enumerate(component['pardicts']):
                if pardict['flag'] == '=':
                    if 'function' in pardict.keys() and len(pardict['function']):
                   # print('function {}, kc,kp,kg={},{},{}'.format(string, kc,kp,kgroup))
                        eval(pardict['function'])
                    else:
                        print("++++ Wrong json syntax, missing or empty 'function' in 'model_guess'")
                        return kc, kp, 0
    return -1,-1,-1 # none has failed the syntax check

def get_gtotals(suite):
    """
    calculates the grand totals and group totals for multi run multi group 

    input is self.suite of class musuite
    returns totalcounts (list of run total str)
            groupcounts (list of lists group total str) 
            nsbin, maxbin (str)
        used in mudashed
    """

    from numpy import array, concatenate
    # called only by self.suite after having loaded a run or a run suite

    ###################
    # grouping set 
    # suite.grouping['forward'] and suite.grouping['backward'] are np.arrays of integers
    # initialize totals
    ###################
    grc = []
    for k,grpdict in enumerate(suite.grouping):
        if not k: # k is 0
            gr = concatenate((grpdict['forward'],grpdict['backward']))
            grc.append(gr)
        else:
            g = concatenate((grpdict['forward'],grpdict['backward']))
            grc.append(g)
            gr = concatenate((gr,g))

    ts,gs =  [],[]
    n1 = suite.offset+suite.nt0[0]
    for k,runs in enumerate(suite._the_runs_):
        for j,run in enumerate(runs): # add values for runs to add
            tsum = 0
            ggs = []
            gts = []
            for j,group in enumerate(grc):
                gggs =[]
                for counter in group:
                    gsum = 0
                    if suite.datafile[-3:]=='bin' or suite.datafile[-3:]=='mdu' or suite.datafile[-4.:]=='root':
                        n1 = suite.offset+suite.nt0[counter] 
                    histo = array(run.get_histo_vector(counter,1)).sum() 
                    gsum += histo
                if j==0:
                    gsum0 = gsum
                gggs.append('{:.2f}'.format(gsum0/1e6)+'Mev'+', {:.2f}'.format(gsum/1e6)+'Mev')
            ggs.append(gggs)
            for counter in range(run.get_numberHisto_int()):
                if suite.datafile[-3:]=='bin' or suite.datafile[-3:]=='mdu' or suite.datafile[-4.:]=='root':
                    n1 = suite.offset+suite.nt0[counter] 
                histo = array(run.get_histo_vector(counter,1)).sum() 
                tsum += histo
            gts.append('{:.2f}'.format(tsum/1e6)+'Mev')
        gs.append(ggs)
        ts.append(gts)
        # print('In get totals inside loop,k {}, runs {}'.format(k,runs))
    return ts, gs, '{:.3}'.format(suite._the_runs_[0][0].get_binWidth_ns()), str(suite.histoLength)
   
###############################################################################
# mufit tools:
# int2min methods: generate guess values, errors, fixed, limits, names, pospar
#                  for minuit parameters
#   int2min :               value, errors etc for Minuit, model_guess directs the fit
#   int2min_global :        same when globpardicts are Minuit parameters
#   int2fft :               needs refurbishing
#   int2_method_key:        generates components for mumodel _load_ fit function
#   int2_global_method_key  same for globpardicts driven cases
#   min2int:                from minuit parameters to model_guess parameters
#   min2int_global:         same for globpardicts driven case
###############################################################################
 
def _nparam(model):
    """
    return numbers of internal, minuit and free ('~'] parameters

    input: dashboard['model_guess']
    output: ntot, nmintot, nfree
        used in tools.int2min, used in mufit
    """

    number_components = len(model)
    ntot = sum([len(model[k]['pardicts']) 
                                 for k in range(number_components)]) # total number of component parameters
    flag = [pardict['flag'] for component in model for pardict in component['pardicts']]
    nmintot = ntot - sum([1 for k in range(ntot) if flag[k]=='=']) # ntot minus number of functions 
    nfree = nmintot - sum([1 for k in range(ntot) if flag[k]=='!']) # ntot minus number of fixed parameters 
    return ntot, nmintot, nfree
    

def int2min(dashboard,runs,guess=True):
    """
    returns values, errors, fixed, limits, names, pospar to initialize Minuit A1 A20 B1 B20 & calib
    
    input: 
        dashboard (json file)
        runs,              not used, compatibility with other int2min_global
        guess = True  dashboard["model_guess"], default, lists of dicts
                False dashboard["model_results"] 
    output: a list of lists:  
        values: minuit parameter values, either guess of result
        errors: their steps/stds
        fixed: True/False for each
        limits: [low, high] limits for either or [None,None]  
        names: name of parameter 'x_label' for each parameter
        pospar: parameter for which component is positive parity, eg s in e^{-(s*t)^2/2}
        used in mufit mufitplot
    """

    # no globpardicts, Minuit performs a 1d asymmetry fit 
    from mujpy.tools.tools import _nparam

    model = dashboard['model_guess'] if guess or 'model_results' not in dashboard.keys() else dashboard['model_results']
    dum, ntot, dum  = _nparam(model)
    
    positive_parity = ['Δ','σ']                                                    
    #####################################################
    # the following variables contain the same as input #
    # parameters to iMinuit, removing '='s (functions)  #
    #####################################################
    
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

    return val, err, fix, lim, names, pospar

def int2min_global(dashboard,runs,guess=True):
    """
    returns values, errors, fixed, limits, names, pospar to initialize Minuit C1 C2 A21 B21 & calib
     
    input:
        dashboard: json loaded dict, contains "globpardicts" 
        runs: list of lists of run numbers (strings)
    output:
        val, err, fix, lim, nam, pospar: ordered lists of minuit parameters
            N globpardicts, in internal order
               ("_"+runs[0] added to the M "#" parameter name) 
            replicas of groups of the M "#" parameters, in internal precedence, 
                for runs[1:] with "_"+runs[k] added to name
    total: min_par = N-M+len(runs)*M minuit parameters (not an output)
    Works for any global fit, both plain and calib, identified by the presence of globpardicts
        multirun global fits C1 C2 contain hash "#" user parameter flag(s)
        used in mufit
    """

    # C1 C2 have hashes, A21 and B21 do not, B21 is ok if runs is sliced to one (at a time) 

    # positive_parity = ['Δ','σ']  # not used
    val, err, fix, lim, nam, pospar = [], [], [], [], [], []           

    # first scan globpardicts and model "#" parameters
    pardicts = dashboard["globpardicts_guess"] if guess or not "globpardicts_result" in dashboard.keys() else dashboard["globpardicts_result"]

    # model = dashboard["model_guess"] if guess or not "model_result" in dashboard.keys() else dashboard["model_result"]

    # scan pardicts and add their val, err, fix, lim, nam, pospar to the list
    hash_parameter_id = []
    min_par = len(pardicts)
    for k_internal,pardict in enumerate(pardicts): # k_internal is the internal order (same of dash) 
        if 'flag' not in pardict.keys():
            pardict['flag'] = '~'
        # set defaults if not set and check lists
        if type(pardict['value']) is list:
            first_value = pardict['value'][0]
            if len(pardict['value']) < len(runs):
                last_value = pardict['value'][-1]
                for j in range(len(pardict['value']),len(runs)):
                    pardict['value'].append(last_value)
        else:
            first_value = pardict['value']
        val.append(first_value)
        if 'limits' not in pardict.keys():
            pardict['limits'] = [ None, None ]
        if 'error' not in pardict.keys():
            pardict['error'] = first_value/20. if first_value else 0.1
        errstd = 'error' if guess else 'std'
        if type(pardict[errstd]) is list:
            first_error = pardict[errstd][0]
            if len(pardict[errstd]) < len(runs):
                last_error = pardict[errstd][-1]
                for j in range(len(pardict[errstd]),len(runs)):
                     pardict[errstd].append(last_error)
        else:
            first_error = pardict[errstd]
        if 'positive_parity' in pardict.keys(): 
            if pardict['positive_parity']:
                pospar.append(k_internal)
                pardict['limits'][0] = 0. 
            # print('debug tools int2min_multirun: pospar {} lim({}) = {}'.format(pardict["name"],k, pardict['limits']))
                    
        err.append(pardict[errstd])
        lim.append(pardict['limits'])
        if pardict['flag'] in ['~','!']:
            nam.append(pardict['name'])
            fix.append(False)
            if pardict['flag'] == '!':
                fix.append(True)
        elif pardict['flag'] == '#': # for A21, B21 fits this is never the case and leaves hash_parameter_id empty
            nam.append(pardict['name']+'_'+runs[0][0])
            hash_parameter_id.append(k_internal)
            fix.append(False)
        else: # 'flag' is an unidentified character
            return False,_,_,_,_,_,_     

# now scan the other runs
    for krun,run in enumerate(runs[1:]): # run is a list of run number strings
        # the loop is skipped for A21, B21

        krun += 1 # first '#' minuit parameter is already allocated with userpars
        for k_internal in hash_parameter_id:
            if type(pardicts[k_internal]['value']) is list:
                value = pardicts[k_internal]["value"][krun]
            else:
                value = pardicts[k_internal]["value"]
            # either list of guess values, one per run, or single guess value 
            name = pardicts[k_internal]["name"]+'_'+run[0]
            keys = pardicts[k_internal].keys()
            err_std = 'error' if guess or 'std' not in keys else 'std'
            if type(pardicts[k_internal][err_std]) is list:
                errstd = pardicts[k_internal][err_std][krun]
            else:
                errstd = pardicts[k_internal][err_std]
            limits = pardicts[k_internal]['limits']
            fix.append(False) # can only be not-fixed 
            if 'positive_parity' in keys: 
                pospar.append(min_par)
            val.append(value)
            nam.append(name)
            err.append(errstd)
            lim.append(limits)
            
            min_par += 1 # increment after, final min_par is number of minuit parametersi

#    print('debug tools int2min multirun multigroup: names {}, values {}'.format(nam[0],val[0]))
#    for k in range(1,len(val)):
#        print('                                       {},        {}'.format(nam[k],val[k]))
#    print('int2_multirun_multigroup len(val) {}, len(err) {}, len(fix) {}, len(nam) {}'.format(len(val), len(err), len(fix), len(nam)))
    return val, err, fix, lim, nam, pospar # all simple lists of sequential parameters, minuit order 
   
##################################
# method and key methods: provide component methods 
#                           and parameter key for eval(key) in _add_
#   int2_method_key :                single run single group 
#   int2_global_method_key :         single run multi group
#   int2_multirun_grad_method_key :  same with grad

def int2_method_key(dashboard,the_model,nruns):
    """
    returns A1, A20, B1, B20 & calib [[method [keys ...]], ...] for mumodel._load_

    input: 
       dashboard, the dashboard dict structure 
       the_model,  a fit model instance (not necessarily loaded)
       nruns,        unused here, for homogeneity with int2_global_model_key
    output: 
       a list of lists, the inner lists contain each
         method,  a mumodel component method, in the order of the model components
                   for the use of mumodel._add_.
         keys,   a list of as many lambda functions as the parameters of teh component
                 hard coding the translated "function" string for fast evaluation.
    This function applies tools.translate to the parameter numbers in formulas:
    dashboard "function" is written in terms of the internal parameter index,
    while Minuit parameter index skips shared or formula-determined ('=') parameters  
        used in mufit mufitplot
    """

    from mujpy.tools.tools import translate, set_key

    model_guess = dashboard['model_guess']  # guess surely exists and model is the same

    ntot = sum([len(model_guess[k]['pardicts']) for k in range(len(model_guess))]) # max nint >= max nmin
    lmin = [] # initialize the minuit parameter index of dashboard function indices 
    nint = -1 # initialize the number of internal parameters
    nmin = -1 # initialize the number of minuit parameters
    method_key = []
    function = [pardict['function'] for component in model_guess for pardict in component['pardicts']]
    for k in range(len(model_guess)):  # scan the model
        name = model_guess[k]['name']
        # print('name = {}, model = {}'.format(name,self._the_model_))
        is_al = name=='al'
        bndmthd = None if is_al else the_model.__getattribute__(name) 

        keys = [] # keys are appended also for 'al' (never used)
        flag = [item['flag'] for item in model_guess[k]['pardicts']]
        for j,pardict in enumerate(model_guess[k]['pardicts']): 
            nint += 1  # internal parameter incremente always   
            if flag[j] == '=': #  function is written in terms of nint
                # nint must be translated into nmin 
                string = translate(lmin,pardict['function']) # here is where lmin is used
                if not string: # there was a exception in translate, stopping the fit
                    return False
                # translate substitutes lmin[n] where n is the index read in the function (e.g. p[3])
                # n is the parameter index but previous indices may also be "=", i.e. not in the minuit indices
                #print('debug int2_method_key: key = {}, p =  {}'.format(
                key_as_lambda = set_key(string) # NEW! calculates simple functions and speedup
                keys.append(key_as_lambda) # the function key in keys will be evaluated, key(p), inside mucomponents
                lmin.append(ntot) # an illegal index, cannot exceed ntot-1
            else:# flag[j] == '~' or flag[j] == '!'
                nmin += 1
                key_as_lambda = set_key('p['+str(nmin)+']') # NEW! calculates simple functions and speedup
                keys.append(key_as_lambda) # the function key in keys will be evaluated, key(p), inside mucomponents
                lmin.append(nmin) # 
                # isminuit.append(True)
        method_key.append([bndmthd,keys]) 
    return method_key

def int2_global_method_key(dashboard,the_model,nruns):
    """
    returns A21, B21, C1, C2 & calib [[method [keys ...]], ...] for mumodel._load_

    input: 
        dashboard, the  full dashboard dict structure
        the_model, i.e. fit._the_model_, is an instance of mumodel
        nruns = len(runs) in fit slice (not from suite!), ngroups = from dashboard "function_multi")
    output: method_key = [method, keys] for all mumodel._add_ versions
                method is nfunc(t,*par) with pars list: 
                2d vector nfunc(t,*par) by cstack, pars list of lists for A21, B21, C1 fits
                3d vector nfunc(t,*par) by ccstack, pars list of lists of lists for C2 fits
                (calib included)
                :: 
                     e.g mumodel.bl(t,A,λ) for runs [['252'],['253'],['254']], 2-group fit 
                         with userpars A1, A2 '~', and  λ '#'
                the 3d method arguments are
                (t,[[[A1,λ_252],[A2,λ_252]],[[A1,λ_253],[A2,λ_253]],[[A1, λ_254],[A2, λ_254]])
                    run0 gr0        gr1     run1 gr0       gr1      run3 gr0        gr1        
                encoded in a lambda function for each key string 
                the run-list of group-lists of strings from "function", "function_multi" 
                with run-specific minuit parameters for λ '#'
        calls relocate_keys to sort the lambda functions
        used in mufit
    """

    from mujpy.tools.tools import function_multi_in_components, ccstack, cstack, nostack
    from mujpy.tools.tools import relocate_keys
    #from numpy import array, array_equal

    ################################################################
    # New unified global scheme: A21, B21, C1, C2 + calib versions #
    #                                                              #
    #    minuit parameters are globpardicts (guess or _results)    #
    #         flag =  '~' global                                   #
    #                 '!' fixed                                    #
    #                 '#' local, replicated for each run           #
    #    model parameters all have '=' flag,                       #
    #                               either function or f_multi     #
    #      (when k_int, k'_int are '#' local userpars              #
    #        strings can only be 'p[k_in]'  (function)             #
    #        or ['p[k_int]','p[k'_int]',...] (function_multi)      #
    #      more complex functions only refer to global userpars)   #
    ################################################################

    model = dashboard['model_guess'] # same model as for '_results'
    mask_function_multi = function_multi_in_components(dashboard)# True for function_multi, False otherwise
    pardicts = [pardict for component in model for pardict in component['pardicts']]
    ngroups = len(pardicts[mask_function_multi.index(True)]['function_multi']) if sum(mask_function_multi) else 1 # counts groups
    # dofit_ must check that the global multigroup fits do have at least one 'function_multi'
    stack = ccstack if nruns > 1 and ngroups > 1 else cstack if nruns > 1 or ngroups > 1 else nostack 
    bndmthd = {} # to load the [name] key
    method_key = []
    relkeys = relocate_keys(dashboard,nruns,ngroups)
    for kc,component in enumerate(model):  # scan model components  
        name = component['name'] # must mnatch those in mucomponents
        is_al = name=='al'
        bndmthd[name] = None if is_al else lambda x,*pars, name=name : stack(the_model.__getattribute__(name),x,nruns,ngroups,*pars)
        # A21, B21, C1, C2
        if nruns >1 and ngroups > 1: # C2 needs a run list of group lists of parameter lists for component kc
            # keys = [[keygroup[kc]] for keyrun in relkeys for keygroup in keyrun]
            keys = [[keygroup[kc] for keygroup in keyrun] for keyrun in relkeys]
        #elif ngroups > 1: # A21 or B21, same fit, need a group list of parameter lists for component kc
        #    keys = [keygroup[kc] for keygroup in relreys]
        #elif nruns > 1: # C1 needs a run list of parameter lists for component kc
        #    keys = [keyrun[kc] for keyrun in relkeys]
        else:
            keys = [k[kc] for k in relkeys] # k stands for keygroup (A21 or B21) or keyrun (C1)
        method_key.append([bndmthd[name],keys]) # vectorialized method, with keys 
    return method_key

def relocate_keys(dashboard,nruns,ngroups):
    """
    activates keys for function or function_multi strings in global fits

    input
        dashboard
        nruns       number of runs in fit slice
        ngroups     number of groups in suite
    output
        keys        list of ... lists (list_depth 2 for A21, B21, C1, or 3, for C2)
        called by tools.int2_global_method_key
    """

    from mujpy.tools.tools import  translate, get_indices, set_key

    model = dashboard['model_guess'] # same model as for '_results'
    globpardicts = dashboard['globpardicts_guess'] # same as '_results'
    mask_hash =  [pardict['flag'] == '#' for pardict in globpardicts] # True for flag = "#", they also must have either function of function_multi strings
    hash_indices = [k for k,m in enumerate(mask_hash) if m] # e.g. 3,5 if userpar[3,5] are hash user parameters]
    hash_strings = [str(k) for k in hash_indices]
    #print('mufit int2_global_method_key hash_indices = {}'.format(hash_indices))
    khash = len(globpardicts) # index for the next minuit parameter, the first nhash being userparameters
    lmin = list(range(khash)) # initially just the indices of the userparameter


    if nruns >1 and ngroups > 1: # C2
        keys = []
        for krun in range(nruns):
            keygroup = []
            store = []
            for kgroup in range(ngroups):
                keycomp = []
                #print('mufit int2_global_method_key -------------- krun kgroup {} {}'.format(krun,kgroup))
                for kc,component in enumerate(model):  # scan model components  
                    # nintcomp = nint # increments nintcomp for next component
                    key = []
                    for pardict in component['pardicts']:
                        string = pardict["function"] if "function" in pardict.keys() else pardict["function_multi"][kgroup]
                        # string is either 'p[k1]' or e.g. '(1-p[k1)*p[k2])]
                        ind_strings =  get_indices(string)
                        #print('string is {}'.format(string))
                        if krun!=0: # skip krun=0, already in userpars, doesn't need translation
                            ind = [int(k) for k in ind_strings] # substitute the indices of the first run
                            for k in ind: 
                                if k in hash_indices: # a hash parameter
                                    #print('mufit int2_g k, lmin[k], khash] {},{},{}'.format(k, lmin[k],khash)) 
                                    if not k in store: 
                                        lmin[k] = khash # update lmin for this translation
                                        store.append(k)
                                        khash += 1 # increment after encoding, finishes with khash = number of minuit parameters
                                    string = translate(lmin,string) # substitutes nhash, encoded in lmin to k_internal 
                                    #print('translated string is {}'.format(string))
                                    if not string: # there was a exception in translate, stopping the fit
                                        return False
                        key_as_lambda = set_key(string) # key(p), non hash is encoded as is
                        key.append(key_as_lambda) # key are for all parameters of component kc
                    keycomp.append(key) # over components
                keygroup.append(keycomp)
            keys.append(keygroup)
    elif ngroups>1: # A21 or B21 no hash parameters
        keys = []
        for kgroup in range(ngroups):
            keycomp = []
            #store = []
            for kc,component in enumerate(model):  # scan model components  
                # nintcomp = nint # increments nintcomp for next component
                key = []
                for pardict in component['pardicts']:
                    string = pardict["function"] if "function" in pardict.keys() else pardict["function_multi"][kgroup]
                    key_as_lambda = set_key(string) # key(p), non hash is encoded as is
                    key.append(key_as_lambda) # for all parameters of component kc
                keycomp.append(key) # over components
            keys.append(keycomp) # over runs
    elif nruns > 1: # C1, pure hash translation
        keys = []
        for krun in range(nruns):
            keycomp = []
            store = []
            for component in model:  # scan model components  
                # nintcomp = nint # increments nintcomp for next component
                key = []
                for pardict in component['pardicts']:
                    string = pardict["function"]
                    # string is either 'p[k1]' or e.g. '(1-p[k1)*p[k2])]
                    ind_strings =  get_indices(string)
                    #print('string is {}'.format(string))
                    if krun!=0: # skip krun=0, already in userpars, doesn't need translation
                        ind = [int(k) for k in ind_strings] # substitute the indices of the first run
                        #print('mufit int2_global_method_key ind = {}'.format(ind))
                        for k in ind: 
                            if k in hash_indices: # a hash parameter
                                #print('mufit int2_g k, lmin[k], khash] {},{},{}'.format(k, lmin[k],khash)) 
                                if not k in store: 
                                    lmin[k] = khash # update lmin for this translation
                                    store.append(k)
                                    khash += 1 # increment after encoding, finishes with khash = number of minuit parameters
                                string = translate(lmin,string) # substitutes nhash, encoded in lmin to k_internal 
                                if not string: # there was a exception in translate, stopping the fit
                                    return False
                    #print('tools int2_global_method_keys krun {} string {}'.format(krun,string))
                    key_as_lambda = set_key(string) # key(p). non hash is encoded as is
                    key.append(key_as_lambda) # for all parameters of the component
                keycomp.append(key) # for all components of the model
            keys.append(keycomp) # for all runs 
    return keys

def set_alpha(pars,nrun,ngroup):
    """
    returns np.array of alpha values from Minuit pars to multiply mumodel._yb_

    input pars, minuit parameters of an al.... multi(run, group, both)  fit
          nrun, number of runs in slice
          ngoup, number of groups in slice
        used in mucomponents
    """

    from numpy import array, kron, ones, newaxis
    #alpha = pars[0:ngroup] # expects an many al parameters as the number of detector groups
    # is converted into numy array that allows direct numpy multiplication with the backward yb array 
    # to produce an array of the same shape as yb
    #return kron(ones((nrun,1)),kron(alpha,ones((1,1))))[:,:,newaxis]
    #print('tools set_alpha pars =  {}'.format(pars))
    if (nrun,ngroup) == (1,1): # 1d
        return array([pars[0]])
    elif nrun == 1: # 2d
        return array([[pars[k]] for k in range(ngroup)])
    elif ngroup == 1: # 2d
        return array([[pars[0]] for k in range(nrun)])
    else: # 3d
        alpha = array([[[pars[k]] for k in range(ngroup)] for kr in range(nrun)])
        return alpha 

def set_key(string):   
    """
    coding p=[0.1,0.2,0.3], key = set_key('p[2]'), key(p) produces 0.3

    input: the function string from the json or the mudash dashboard
         e.g. that written in the json file as 'function':'p[0]*(0.5+1/pi*arctan(p[2])'
              or typed into mudash text widget as 'p[0]*(0.5+1/pi*arctan(p[2])',
              (for global fits: indices are already translated!)  
    output: key, a python method, such that in mumodel mucomponents _add_ the command 
            key(p) evaluates the formula 
            the evaluation knows simple numpy functions, see the import in string code, below
    COMMENTS: maybe not the simplest, but it works
    try the following in ipython3
        from mujpy.tools.tools import set_key
        p = [1,2,3,4]
        key = set_key('p[2]')
        key(p)
    Out[] 3
        used in tools relocate_keys, int2_multirun_grad_method_key, int2_method_keys
    """

#    print('debug tools set_key: string = {}'.format(string))
    code = """
from numpy import cos, sin, tan, sinh, cosh, tanh, log, pi, exp, sqrt, real, abs, arctan
def foo():
"""
    string = "    key = eval('"+'lambda p: '+string +"')"
    # print('string ={}'.format(string))
    code = code + string + """
    return key
"""
    # print('code = {}'.format(code))
    exec(code,globals(),globals())   # foo is defined by executing code
    return eval('foo()') # key = set_key('p[2') when p = [10,11,12,13], key(p) returns 12  

def int2_multirun_grad_method_key(dashboard,the_model,nruns):
    """
    methods_keys for analytic gradients [not refurbished]

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
    """

    from mujpy.tools.tools import get_functions_in, diffunc, get_indices, get_number_minuit_internal
    from mujpy.tools.tools import translate_multirun, set_key 
    # first generate dmethod_keys
    # dmethod_keys contains [[m_d,keys]...[m_d,keys]], as many as the model components (e.g. 2 for mgbl)
    # m_d is [method] if no derivative is required (e.g. bl) 
    # or [method, derivative_method] if derivative is required (e.g. mg)  
    # keys is [runkeys,...,runkeys] such that
    # par = [key(p) for key in runkeys] and method(x,*par) and derivative_method(self._x_,*par) produce the additive component for that run
    model = dashboard['model_guess']
    names = [component['name'] for component in model] 
    n_locals =  [pardict["local"] for pardict in dashboard["globpardicts_guess"]].count(True)
    n_globals = len(dashboard["globpardicts_guess"])-n_locals
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
       
def nostack(npfunc,x,nruns,ngroup,*pars):
    """
    mock stacking
    x dummy
    nruns duummy
    ngroups dummy
    pars dummy
          no need to pass args, mucomponents expect (x,*pars)
    output is npfunc
        used by tools.int2_global_method_key
    """
    return npfunc

def cstack(npfunc,x,nruns,ngroups,*pars):
    """
    vectorialize npfunc
    input: 
        npfunc numpy function with input (x,*argv)
        x time
        nruns 1 or >1
        ngroups 1 or >1 (for compatibility with ccstack)
        *pars is a list of lists of parameters, 
              list len n is the output_function_array.shape[0]
    output:
        output_function_array
            stacks vertically n replica of npfunc distributing parameters as in
            (x, *argv[i]) for each i-th replica 
        used by tools.int2_global_method_key mufitplot
    """
    # cstack reproduces the parameter input of a component according to         
    # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
    # where the outer list a replica of the same component method 
    # either over several groups (multigroup) or over several runs (multirun)
    # as of now this method does not work for the multirun multigroup userpar case (C2)

    from numpy import array #concatenate #try with A21 before replacing
    # print('debug tools.cstack: npfunc = {} pars = {}'.format(npfunc,pars))
    # reshape(-1,x-shape[0]) makes as many rows as necessary, each with x.shape[0] columns   
    #return concatenate([npfunc(x,*par) for par in pars]).reshape(max(nruns,ngroups),x.shape[0])
    if nruns == 1:
        return array([npfunc(x,*par) for par in pars])
    else:
        return array([[npfunc(x,*par)] for par in pars])
 
def ccstack(npfunc,x,nruns,ngroups,*pars):
    """
    vectorialize npfunc
    input: 
        npfunc numpy function with input (x,*pars)
        x time
        *pars is a list of lists of lists of parameters, 
    output:
        output_function_array
            stacks vertically n replica of npfunc distributing parameters as in
            (x, *argv[i]) for each i-th replica 
            output_function_array(x,*argv) produces an array of shape (nruns,ngroups,x.shape[0])
        used by tools.int2_global_method_key
    """
    # cstack reproduces the parameter input of a component according to         
    # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
    # where the outer list a replica of the same component method 
    # either over several groups (multigroup) or over several runs (multirun)
    # as of now this method does not work for the multirun multigroup userpar case (C2)

    from numpy import array
    return array([[npfunc(x,*par) for par in parg] for parg in pars])

def min2int(dashboard,values_in,errors_in,krun,kgroup,nruns,ngroups):
    """
    returns names, values, stds, shared list of parameter lists for summary print_components

    input:
        dashboard
        values_in Minuit.values
        errors_in Minuit.errors
        krun,kgroup:    index of run and group in suite._the_runs_
        nruns, ngroups: slice dimensions
              last four not used, homogeneity with the _global version
    output: for all model_guess parameters
        names      component list of parameter lists of names
        values_out component list of parameter lists of values 
        errors_out component list of parameter lists of stds
        shared     component list of parameter list of 'flag'=='=' 
    ready for print_components in summary
    works for: A1 and its calib
               A20, B1, B20 and their calib, iterated in mufit.dofit_
        used by mufit
    """

    from mujpy.tools.tools import translate
    # 
    # initialize
    #
    model_guess = dashboard['model_guess']
    ntot = sum([len(model_guess[k]['pardicts']) for k in range(len(model_guess))]) 
    # total number of internal parameters >= number of Minuit parameters
    
    names, values_out, p, errors_out, e, shared = [], [], [], [], [], []
    nint = -1 # initialize
    nmin = -1
    lmin = []
    flag = [pardict['flag'] for component in model_guess for pardict in component['pardicts']]
    for k,component in enumerate(model_guess):  # scan the model
        component_name = component['name']
        name, value, error, share = [], [], [], []
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
                # print('diretto p_in[{}] = {} -> {}'.format(nmin,p[-1],value[-1]))
                e.append(errors_in[nmin])
                error.append(errors_in[nmin]) # parvalue item is a string
                share.append(False)
            else: # functions, calculate as such
                # nint must be translated into nmin 
                string = translate(lmin,pardict['function'])  
                p.append(eval(string))
                value.append(eval(string))
                # print('{} shared p_in = {} -> {}'.format(string,p[-1],value[-1]))
                e.append(eval(string.replace('p','e')))
                error.append(eval(string.replace('p','e'))) # this is where e is used
                share.append(True)
                lmin.append(ntot) # an illegal index, was working with lmin.append(0)
        names.append(name)
        values_out.append(value)
        errors_out.append(error)
        shared.append(share)
    return names, values_out, errors_out, shared # list of parameter values 
 
def min2int_global(dashboard,p,e,krun,kgroup,nruns,ngroups):
    """
    names, values, stds, nonhash list of parameter lists for summary_global print_components

    input:
        dashboard
        p               Minuit.values
        e               Minuit.errors
        krun,kgroup:    indices of run, group in suite._the_runs_
        nruns, ngroups  slice dimensions
    output: for all the dashboard['model_guess'] parameters
        names      component list of parameter lists of names
        values_out component list of parameter lists of their values 
        errors_out component list of parameter lists of their stds
        nonhash    component list of parameter list True if 'flag'!='#'
    ready for for print_components in summary_global
    iterated in mufit.dofit_ for A21 B21 C1 C2 & calib
        used by mufit
    """

    from mujpy.tools.tools import function_multi_in_components, error_propagation, get_indices
    from numpy import cos, sin, tan, sinh, cosh, tanh, log, pi, exp, sqrt, real, abs, arctan

    globpardicts = dashboard['globpardicts_guess'] # values are taken from p,e, anyway
    model = dashboard['model_guess']
    pardicts = [pardict for component in model for pardict in component['pardicts']]
    mask_function_multi = ['function_multi' in pardict.keys() for component in model for pardict in component['pardicts']]
    mask_hash =  [pardict['flag'] == '#' for pardict in globpardicts]  # works also for A21, B21, no hashes, masl_hash = []
    lmin = list(range(len(globpardicts))) # initialized to [0,1,2,3 ...]#)
    hash_indices = [k for k,m in enumerate(mask_hash) if m]
    new_hash_indices = [k for k,m in enumerate(mask_hash) if m] 
    if mask_hash:
        for k in range(len(hash_indices),len(p)):
            new_hash_indices.append(k) # builds indices of run specific Minuit parameters
    #print('tools min2int_global new_hash_ind {}'.format(new_hash_indices))
    names, values_out, stds_out, nonhash = [],[],[],[]
    khash = len(globpardicts)
    #print('tools min2int_global nruns ngroups krun kgroup {} {} {} {}'.format(nruns,ngroups,krun,kgroup))
    for kr in range(nruns):
        store = hash_indices if kr==0 else []
        for kg in  range(ngroups):
            for component in model:  # scan model components  
                component_name = component['name']
                label = component['label']
                namep, valp, stdp, nonhp = [],[],[],[]
                for j,pardict in enumerate(component['pardicts']): 
                    if j==0: # first parameter of this component
                        namep.append('{}: {}_{}'.format(component_name,pardict['name'],label))
                    else:
                        namep.append('{}_{}'.format(pardict['name'],label))
                    string = pardict["function_multi"][kgroup] if 'function_multi' in pardict.keys() else pardict["function"]
                    ind_strings = get_indices(string) # string is 'p[k_internal]'
                    ind = [int(k) for k in ind_strings] # substitute the function/function_multi indices for hash parameters 
                    for k in ind: 
                        if k in hash_indices: # a hash parameter
                            if not k in store:
                                #print('tools min2int_global lmin k khash {} {} {}'.format(lmin,k,khash))
                                lmin[k] = new_hash_indices[khash] # update lmin for this translation
                                store.append(k)
                                khash += 1
                    ishash = True if any([k in hash_indices for k in ind]) else False
                    string = translate(lmin,string) if ishash else string
                    if kr==krun and kg==kgroup:
                        if not ishash:
                            nonhp.append(True) 
                        else:
                            nonhp.append(False)# does not appènd to the ur hash parameter (krun=0)
                        key = set_key(string)
                        val = key(p)
                        err = eval(error_propagation(string))
                        valp.append(val) #csv_format.format(val)) # formatted value
                        stdp.append(err) #csv_format.format(err)) # formatted error
                if kr==krun and kg==kgroup:
                    #print('tools min2int_global kr = {}, kg = {}'.format(kr,kg))
                    names.append(namep)
                    values_out.append(valp) # all groups in the same run list
                    stds_out.append(stdp) # formats the error
                    nonhash.append(nonhp)
                # always a component list of parameter lists, list_depth = 2

    return names, values_out, stds_out, nonhash

def min2int_names(dashboard):
    """
    like min2int, produces only names

    input: dashboard
    output: names, component list of parameter lists, for csv
            name of parameter[0] contains component name, 
            all names end with '_'+label, eg mgA_0, B_0, φ_0, σ_0 
    called by mufit.prepare_csv_row and prepare_usercsv_row
        used by mufit
    """

    model_guess = dashboard['model_guess']
    names = []
    nint = -1 # initialize
    nmin = -1
    for k,component in enumerate(model_guess):  # scan the model
        component_name = component['name']
        name = []
        label = model_guess[k]['label']
        for j,pardict in enumerate(model_guess[k]['pardicts']): # list of dictionaries, par is a dictionary
            nint += 1  # internal parameter incremented always
            if j==0:
                name.append('{}{}_{}'.format(component_name,pardict['name'],label))
            else:
                name.append('{}_{}'.format(pardict['name'],label))
        names.append(name)
    return names # list of lists of parameter names

def chunk(seq,size):
    """
    use as    for v in chuck(a,n):  to extract n items at a time from a
        used in mufit
    """

    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def list_depth(L):
    """
    returns max depth of a list of ... lists
        used in mufit and mucomponents
    """

    if isinstance(L, list):
        return 1 + max(list_depth(item) for item in L)
    else:
        return 0 
            
def print_components(names,values,errors,shared,global_fit=False):
    """
    returns string with all name = value(std) formatted parameters of a model component

    input: for a component
    	parameter names 
    	parameter values 
    	parameter errors 
        parameter shared (boolean) / nonhash (boolean)
                      is 'flag' == '='     not global_fit
                         'flag' != '#'     global_fit  
        global_fit = False (default), boolean
    output:
        string to print
        "bl.A_fast 0.123(4) bl.λ_fast 12.3(4)S bl.σ_fast 0(0)"
        adds 'G' (global_fit) or 'S' to error of global or shared parameters
        used in mufit
    """

    from mujpy.tools.tools import value_error
    ch = '' if global_fit else 's'
    val_err = [value_error(values[k],errors[k])+ch if shared[k] else value_error(values[k],errors[k]) for k in range(len(shared))]
    out = [' '.join([names[k],'=',val_err[k]]) for k in range(len(names))]
    maxlen = len(max(out,key=len))
    return " ".join([out[k]+(maxlen-len(out[k]))*' ' for k in range(len(out))])

def print_csv_components(values,errors):
    """
    returns string with csv formatted value, std, per all parameters of all model components

    input: for a component
        parameter values
        parameter errors
    output:
        string to print, see print_components for this e.g.
        "0.123,0.004,12.3,0.4,0,0,"
        notice! ends with ","
        used in mufit
    """

    from mujpy.tools.tools import value_error_csv
    return ''.join(value_error_csv(values[k],errors[k]) for k in range(len(values)))

def mixer(t,y,f0):
    """
    returns rotating frame version of input t at frequency f0 [MHz]

    mixer of a time-signal with a reference 
    input
        t time
        y the time-signal
        f0 frequency of the cosine reference
    output
        y_rrf = 2*y*cos(2*pi*f0*t)  
    t is 1d and y is 1-d, 2-d or 3-d but t.shape[0] == y.shape[-1]
    t is vstack-ed to be the same shape as y
        used in mufitplot
    """

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
    """
    returns fy filtered above the 0.2*fy peak freq 

    works for 1-2 d
        used in tools.mixer 
    """

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

def function_multi_in_components(dashboard):
    """
    list of True/False for model_guess components pardicts if they contain the 'function_multi' key

    input full dashboard
    output mask list, 
        len is sum of number of parameters in components of "model_guess"
        True where parameter pardict contains "function_multi" key
        False otherwise
        used inn mufitplot, tools.min2int_global, tools.int2_global_method_key
    """

    return ['function_multi' in pardict.keys() for component in dashboard["model_guess"]  for pardict in component["pardicts"]]                                
def stringify_group(group):
    """
    returns a unique string for a group

        uses in tools.stringify_groups and in mufit
    """

    fgroup, bgroup = group['forward'],group['backward']
    return fgroup.replace(',','_')+'-'+bgroup.replace(',','_')

def stringify_groups(groups,joinch='_'):
    """
    returns a unique string for many groups

        used in mufit
    """

    from mujpy.tools.tools import stringify_group
    strgrp = []
    for group in groups: 
        strgrp.append(stringify_group(group))
    return joinch.join(strgrp)

def rshp(y):
    '''
    reshape numpy array from suite to plot shape

    suite array are 3d 2d or 1d 
    plot must be 1d, for static plot
                 2d, for anim
    returns a squeezed version and 3d is reshaped to 2d
        used in mufitplot and mucomponents
    '''

    if len(y.shape)==3:
        n0,n1 = (y.shape[0]*y.shape[1], y.shape[2]) 
        y_out = y.reshape((n0,n1)).squeeze()
    else:
        y_out = y.squeeze()
    return y_out

##############
# MUDASHED AUX need pruning deprecated methods
##############

def glob2widgets(kp,pardict,flags,keylen):
    """
    unwraps pardict in widgets

    order is 'p[k]' 'name' 'value' 'flag' 'error' 'limits' 'par>0'
        used in mudashed
    """

    from ipywidgets.widgets import HBox, Label, Combobox, FloatText, Dropdown, Text, Checkbox, Layout
                #
    return HBox([
    Label(value=str(kp),layout=Layout(width=keylen[0])), 
    Combobox(options=['α','λ','σ','φ','Δ','β','θ','δ','ν','τ'], #placeholder='ty+sel',
             value=pardict['name'],layout=Layout(width=keylen[1])), 
    FloatText(value=pardict['value'],layout=Layout(width=keylen[2])), 
    Dropdown(options = flags,value=pardict['flag'],layout=Layout(width=keylen[3])), 
    FloatText(value=pardict['error'],layout=Layout(width=keylen[4])),
    Text(value=tup2str(pardict['limits']),tooltip='e.g. 0,1\nor 0,None\nor None,None',layout=Layout(width=keylen[5])),
    Checkbox(value=pardict['positive_parity'],indent=False,layout=Layout(width=keylen[6]))]
    )

def comp2widgets(component,k):
    '''
    returns mudashed HBox of widgets for component dict, index k
    
        Text str(k):name, disabled = True
        Label tag
        Text tag [default str(k)]
        Checkbox FFT value = True 
        Label FFT 
        used in mudashed
    '''

    from ipywidgets.widgets import HBox, Text, Label, Checkbox, Layout
    width = ['15%','8%','5%','12%','7%','4%'] # total 36%
    widgets = [Label(value='component',layout=Layout(width=width[0])),
               Text(value=str(k)+':'+component['name'],disabled=True,layout=Layout(width=width[1])),
               Label(value='tag',layout=Layout(width=width[2])),
               Text(value=str(k),tooltip='par name tag',layout=Layout(width=width[3])),
               Label(value='FFT',layout=Layout(width=width[4])),
               Checkbox(value=True,tooltip='subtract/n(residues FFT)',indent=False,layout=Layout(width=width[5]))]
    return HBox(widgets)

def par2labels(pardict):
    '''
    returns mudashed HBox of Labels for parameters

    see par2widgets
        used in mudashed
    '''

    from ipywidgets.widgets import HBox, Label,  Layout
    width = ['5%','8%','15%','10%','30%'] # total 48%
    widgets = [Label(value='p[k]',layout=Layout(width=width[0])),
               Label(value='Name',layout=Layout(width=width[1])),
               Label(value='Value',layout=Layout(width=width[2])),
               Label(value='Flag',layout=Layout(width=width[3])),
               Label(value='Function',layout=Layout(width=width[4]))]
    return HBox(widgets)

def par2widgets(pardict,k,glob=False):
    '''
    returns mudashed HBox of widgets for pardict, index k, glob selects flag options

        Label str(k)
        Text name disabled = True
        Text value
        Dropdown Flag 
        Text Function
        used in mudashed
    '''

    from ipywidgets.widgets import HBox, Text, FloatText, Label, Dropdown, Layout
    width = ['5%','8%','15%','10%','30%'] # total 48%
    options = ['='] if glob else ['~','!','=']
    value = 0 if 'value' not in pardict else pardict['value']
    if 'function_multi' in pardict:
        function = ';'.join(pardict['function_multi'])
    elif 'function' in pardict:
        function = pardict['function']
    else:
        function = ''
    widgets = [Label(value=str(k),layout=Layout(width=width[0])),
               Text(value=pardict['name'],disabled=True,layout=Layout(width=width[1])),
               FloatText(value=value,tooltip='a float',layout=Layout(width=width[2])),
               Dropdown(options=options,value=pardict['flag'],layout=Layout(width=width[3])),
               Text(value=function,layout=Layout(width=width[4]))]
    return HBox(widgets)

def read_pardict_from_widgets(kids,kmax):
    """
    reads from kids, returns name, value, flag, function (or ValueError)

        kids is an HBox.children of widgets: here accesses .value at index 1,2,3,4(see par2widgets)
        global fits ignore value (nothing ad hoc here)
        strictly allowed flags already imposed by par2widgets Dropdown options
        muvalid checks math syntax
        verify 0 < indices < kmax (number of globals or present internal index k)  
        value checks done by add_step_limits
            name, None, flag, function or function_multi
                           list of 1 or more strings 
        else
            name, value, flag, function
            (errors, limits and pospar are dealt by adds_step_limits                            
                 
        used by mudashed.dump_dashed
    """
    from mujpy.tools.tools import muvalid
    import re

    na, va, fl, fum = [widg.value for widg in kids.children[1:]]
# [skip k] name value flag function or function_multi
    fu = re.split("[;,]",fum)
    for f in fu: 
    # check function syntax
        invalid = False
        if f: invalid = muvalid(f) # avoid check on empty function
        if invalid:
            return invalid # a string
        # check indices
        invalid = not all([int(k)<kmax and int(k)>=0 for k in get_indices(f)])
        if invalid:
            return 'index out of bounds'
    #if not glob: # check invalid 
    pardict = {'name':na,'flag':fl}
    if len(fu)==1:
        pardict['function']=fu[0]
    else:
        pardict['function_multi'] = fu 
    return pardict # a dict

def tup2str(tup):
    """
    translates a tuple of floats to a csv string (for ipwidgets Text)

    (float includes None) 
        used in tools.glob2widgets
    """
    return ','.join(str(tup[i]) for i in range(len(tup)))

def invalid_err_lim(value,error,limits):
    """
    returns a specific ValueError if parameters are inconsistent
        used in mudashed
    """
    if value !=0 and error > abs(value): 
        invalid = 'abs(value)<=error' 
    elif limits[0] and value < limits[0]:
        invalid = 'value < limits[0]'
    elif limits[1] and value > limits[1]:
        invalid = 'value > limits[1]'
    else:
        invalid = '' # False, i.e valid

def validmodel(model):
    """
    checks validity of model name, e.g. "almlmg"
        used in mudashed
    """

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

def chi2std(nu):
    """
    computes 1 std for least square chi2
        used in mufit
    """
    import numpy as np
    from scipy.special import gammainc
    from scipy.stats import norm
    
    mm = round(nu/4)              
    hb = np.linspace(-mm,mm,2*mm+1)
    cc = gammainc((hb+nu)/2,nu/2) # see mulab: muchi2cdf(x,nu) = gammainc(x/2, nu/2);
    lc = 1+hb[min(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nu
    hc = 1+hb[max(list(np.where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nu
    return lc, hc

#############
# GENERAL AUX
#############

def derange(string,vmax,pack=1):
    """
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
        used in mufit and mufitplot
    """
    
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
    """
    parses string, producing a list of runs; 
    expects comma separated items

    looks for 'l','l:m','l+n+m','l:m:-1' 
    where l, m, n are integers
    also more than one, comma separated 

    rejects all other characters

    returns a list of lists of integer
        used in musuite and mudashed
    """
    import re
    from copy import deepcopy as cp
    s = []
    try:
    # substitute ',' followed by spaces by plain ','
        string_in = cp(string)
        string = re.sub(r",\s+", ",", string.strip())
    # substitute multiple consecutive spaces with ',' in strings separated by spaces only
        string = re.sub(r"\s+", ",", string) #
    # systematic str(int(b[])) to check that b[] ARE integers
        for b in string.split(','): # csv
            kminus = b.find(':-1') # '-1' means reverse order
            kcolon = b.find(':') # ':' and '+' are mutually exclusive
            kplus = b.find('+')
            #print(kminus,kcolon,kplus)

            if kminus<0 and kcolon<0 and kplus<0: # single run, no run addition
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
        return [], 'not a valid string: {}'.format(string_in)

def findall(p, s):
    """Yields (provides an iterator for) 
    all the positions of the pattern p in the string s.
    
    usage:
        for i in findall('x','xaxxa, che xifo!'):
            print(i) 
        used by translate and other tools methods (but see also re.findall)
    """
    i = s.find(p)
    while i != -1:
        yield i
        i = s.find(p, i+1)
   
def get_datafilename(datafile,run):
    """
    datafilename = template, e.g. '/fullpath/deltat_gps_tdc_0935.bin'
    run = string of run digits, e.g. '1001'
    returns '/fullpath/deltat_gps_tdc_1001.bin'
        used in musuite
    """
    
    import re
    datafile_suffix = datafile[-5:].split('.')[1]
    dot = '.'
    dot_suffix = dot+datafile_suffix
    datafile_nosuffix = datafile.split(dot_suffix)[0]
    # print(datafile_nosuffix)    
    # assuming that the datafile_nosuffix is path+datafile_header+run_number e.g. path+'deltat_tdc_gps_3245'
    padded = re.match('.*?([0-9]+)$', datafile_nosuffix).group(1) # run string of digits, including padding zeros
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
    # print('datafilename in tools: {}'.format(datafilename))
    return datafilename

def get_grouping(groupcsv):
    """
    input
      groupcsv is a shorthand csv string, e.g. '1:3,5' or '1,3,5' etc.
      contained in self.suite.group[k]["forward] of self.suite.group[k]["backward"]
          (the k-th detector group of a multi group fit)
    output
     grouping is an np.array of indces, 0 based
        used in musuite
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
                grouping = np.concatenate((grouping,np.array(list(int(w) for w in q[1:])))).astype(int)

        grouping -=1 # this is counter index, remove 1 for python 0-based indexing 
    except:
        grouping = np.array([-1]) # error flag
        
    return grouping

def init_csv_row(filespec, the_run, group = False):
    """
    writes beginning of csv row with nrun T [T eT T eT] B

    filespec    data mime to identify facility
    the_run     this run instance
    [group]     stringify(group)
    updated for ISIS [PSI] 
            *** must update for root datasets
        used in mufit
    """

    nrun = str(the_run.get_runNumber_int())
    #print('tools init_csv_row nrun {}'.format(nrun))
    Bstr = the_run.get_field()
    if filespec=='bin' or filespec=='mdu':
        TsTc, eTsTc = the_run.get_temperatures_vector(), the_run.get_devTemperatures_vector()
        n1,n2 = spec_prec(eTsTc[0]),spec_prec(eTsTc[1]) # calculates format specifier precision
        form = '{},{:.'+'{}'.format(n1)+'f},{:.'+'{}'.format(n1)
        form += 'f},{:.'+'{}'.format(n2)+'f},{:.'+'{}'.format(n2)+'f},{:.0f},'
        if group:
            form += '{},'
            return form.format(nrun, TsTc[0],eTsTc[0],TsTc[1],eTsTc[1], float(Bstr[:Bstr.find('G')]),group)
        else:
            return form.format(nrun, TsTc[0],eTsTc[0],TsTc[1],eTsTc[1], float(Bstr[:Bstr.find('G')]))
    elif filespec=='nxs':
        Ts = the_run.get_temperatures_vector()
        n1 = '1'       
        form = '{},{:.'+'{}'.format(n1)+'f},{},' 
        if group:
            form += '{},'
            return form.format(nrun, Ts[0], Bstr[:Bstr.find('G')],group)
        else:
            return form.format(nrun, Ts[0], Bstr[:Bstr.find('G')])
   
def chi2_csv(chi2,lowchi2,hichi2,alpha,offset):
    """
    input:
        chi2, chi2-sdt, chi2+sdt, alpha, offset (bins)
    output:
        cvs partial row with these values and timestring
        used by mufit
    """
    from time import localtime, strftime
    
    #echi = min(chi2-lowchi2,hichi2-chi2) # was max
    n1=3   #n1 = spec_prec(echi) # calculates format specifier precision
    form = '{:.'+'{}'.format(n1)+'f},{:.'+'{}'.format(n1)+'f}'
    form += ',{:.'+'{}'.format(n1)+'f}'
    row = form.format(chi2,lowchi2,hichi2)
    if isinstance(alpha,list):
        for a in alpha:
            row += ',{:.4f}'.format(a)
    else:
        row += ',{:.4f}'.format(alpha)
    row += ',{},{}'.format(offset,strftime("%d.%b.%H:%M:%S", localtime()))
    return row

def write_csv(header,row,the_run,file_csv,filespec,scan=None):
    """
    input :
        header, the model specific csv header 
                to compare with that of the csv file
        row, the line to be added to the csv file
        the_run, run instance (first one for added runs)
        file_csv, full path/filename to csv file 
        filespec, 'bin', 'mdu' or 'nsx'
       scan: T, B or None
    output:
        two strings to write on console
    writes onto csv finding the right line
    writes a new file if csv does not exist or is incompatible (writes ~ version)
        used by mufit
    """
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
                notexistent = True # assume row yet non-existent fit
                for nline,line in enumerate(f_in.readlines()):
                    if nline==0:
                        if header!=line: # different headers, use more recent
#                           warn = '**** write_csv\n{}\n{}\n'.format(header,line)
                            raise # exits this try 
                        else:
                            lineout.append(header)
                            checkgroup =[x for x, g in enumerate(re.split(" |,|, ",header)) if g.find('group')+1]
                    elif float(re.split(" |,|, ",line)[csv_index]) < rowvalue: # reinsert older fit
                        lineout.append(line)
                    elif float(re.split(" |,|, ",line)[csv_index]) == rowvalue:
                        if checkgroup: 
                            if re.split(" |,|, ",line)[checkgroup[0]]== re.split(" |,|, ",row)[checkgroup[0]]:
                                lineout.append(row) # substitute existing fit, multigroup sequential
                                notexistent = False
                            else:
                                lineout.append(line) # reinsert older fit, multigroup sequential
                        else:
                            lineout.append(row) # substitute existing fit, singlegroup
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
            strgrp0 = 'minuit parameters' if file_csv[0]=='U' else '--'  
            strgrp = re.split(" |,|, ",row)[checkgroup[0]] if checkgroup else strgrp0
            return 'Run {}: {} ***'.format(nrun,
                                           get_title(the_run)), '{} {}:  row added to {}'.format(nrun,strgrp,file_csv)

        except Exception as exc: # incompatible headers, save backup and write a new file
            #print('write_csv excetption: {}'.format(exc))
            os.rename(file_csv,file_csv+'~')
            with open(file_csv,'w') as f:
                f.write(header)
                f.write(row)
            file_csv = file_csv[file_csv.rfind('/')+1:]
            return 'Run {}: {}'.format(nrun,
                    get_title(the_run)),'.  NEW file {} [backup in {}]'.format(
                                                                         file_csv,
                                                                         file_csv+'~')
            
    else: # csv does not exist
        print('file {} not found'.format(file_csv))
        with open(file_csv,'w') as f:
            f.write(header)
            f.write(row)
        file_csv = file_csv[file_csv.rfind('/')+1:]
        return 'Run {}: {} ***'.format(nrun,
                        get_title(the_run)),'.  Log in NEW {}'.format(file_csv)

def get_title(run,notemp=False,nofield=False):
    """
    form standard psi title
        used in tools mufit musuite mudashed
    """
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
    """
    output 
        list of run and title strings
            each run and group in the run replicates its run number + title
        used only in mufitplot (fit and fft)  
    """
    from mujpy.tools.tools import get_title
    run_title = []
    for run in the_suite._the_runs_:
        for kgroup in range(len(the_suite.grouping)):
                run_title.append(str(run[0].get_runNumber_int())+'-'+get_title(run[0]))
    return run_title
    
def get_nruns(the_suite):
    """
    get nrun strings
        used in mufitplot
    """
    nruns = []
    for k,run in enumerate(the_suite._the_runs_):
        nruns.append(str(run[0].get_runNumber_int()))
    return nruns

def muvalid(string):
    """
    parse function 

    CHECK WITH MUCOMPONENT, THAT USES A DIFFERENT SCHEME

    accepted functions are RHS of agebraic expressions of parameters p[i], i=0...ntot  
        used in tools.read_pardict_from_widgets
    """
    import re
    error_message = ''
    if string.strip() !='': # empty and blank strings are validated 
        pattern = re.compile(r"p\[(\d+)\]") # find all patterns p[*] where * is digits
        test = pattern.sub(r"a",string) # substitute "a" to "p[*]" in s
        try: 
            safetry(test) # should select only safe use (although such a thing does not exist!)
        except Exception as e:
            error_message = 'Function: {}. Wrong or not allowed syntax: {}'.format(string,e)
    return error_message

def muzeropad(runs,nzeros=4):
    """

    runs is a string containing the run number
    nzeros the number of digit chars in the filename
    PSI bin: nzeros=4
    ISIS nxs nzeros=8
    returns the runs string 
    with left zero padding to nzeros digits
    """
    zeros='0'*nzeros
    if len(runs)<len(zeros):
        return zeros[:len(zeros)-len(runs)]+runs
    elif len(runs)==len(zeros):
        return runs

def path_file_dialog(path,spec):
    """
    launch tkinter filedialog in path, spec is filename after dot
        used in mudashed
    """

    import tkinter
    from tkinter import filedialog
    import os
    tkinter.Tk().withdraw() # Close the root window
    spc, spcdef = '.'+spec,'*.'+spec
    in_path = filedialog.askopenfilename(initialdir = path,filetypes=((spc,spcdef),('all','*.*')))
    return in_path

def rebin(x,y,strstp,pack,e=None):
    """
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
       used in mufitplot mucomponents
    """

    from numpy import floor, sqrt, zeros, where
    # yy is a slice of  y, reshape acts as if yy where a deepcopy, y is not modified
    start,stop = strstp
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
        nruns_or_groups = y.shape[0] # number of runs/groups
        yr = zeros((nruns_or_groups,m))
        if e is not None:
            er = zeros((nruns_or_groups,m))
        for k in range(nruns_or_groups): # each row is a run
            yy = y[k][start:start+mn]  # slice row
            yy = yy.reshape(m,pack)  # temporaty 2d
            yr[k] = yy.sum(1)/pack # rebinned row
            if e is not None:
                ey = e[k][start:start+mn]   # slice row
                ey = ey.reshape(m,pack)  # temporaty 2d
                er[k] = sqrt((ey**2).sum(1))/pack  # rebinned row        
    elif len(y.shape)==3:        
        nruns, ngroups, _ = y.shape # number of groups, runs
        yr = zeros((nruns,ngroups,m))
 
        if e is not None:
            er = zeros((nruns,ngroups,m))
        for krun in range(nruns): 
            for jgroup in range(ngroups):  
                yy = y[krun][jgroup][start:start+mn]  # slice 
                yy = yy.reshape(m,pack)  # temporaty 2d
                yr[krun][jgroup] = yy.sum(1)/pack # ebinned row
                if e is not None:
                    ey = e[krun][jgroup][start:start+mn]   # slice row
                    ey = ey.reshape(m,pack)  # temporaty 2d
                    er[krun][jgroup] = sqrt((ey**2).sum(1))/pack  # rebinned row
#                if list(where(er[krun][jgroup]==0)[0]):
#                    print('rebin: zero er!!!')
#                else:
#                    print(where(er[krun][jgroup]==0))
    if e is not None:
        return xr,yr,er
    else:
        return xr,yr

def safetry(string):
    """
        used by muvalid(used by read_pardict_from_widgets)
    """
    from math import acos,asin,atan,atan2,ceil,cos,cosh,degrees,e,exp,floor,log,log10,pi,radians,sin,sinh,sqrt,tan,tanh
    safe_list = ['a','acos', 'asin', 'atan', 'atan2', 'ceil', 'cos', 'cosh', 'degrees', 'e', 
                 'exp', 'floor', 'log', 'log10', 'pi', 'pow', 'radians', 'sin', 'sinh', 'sqrt', 'tan', 'tanh']
    # 	use the list to filter the local namespace
    a = 0.3
    safe_dict={}
    for k in safe_list:
        safe_dict[k]=locals().get(k)
    return eval(string,{"__builtins__":None},safe_dict)

def scanms(y,n):
    """
    produces guess for hifi t=0 bin, to be fed to a step fit function
        used by musuite 
    """
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
    """
    format specifier precision::

        0 for a > 1.0
        1 for 1.0 > a > 0.1
        2 for 0.1 > a > 0.01 etc.
        used by mufit
    """
    import numpy as np
    return int(abs(min(0.,np.floor(np.log10(abs(a)))))) 
   
def step(x,a,n,dn,b):
    """
    step function from norm.cdf
        used in musuite
    """
    from scipy.stats import norm
    # error function as step function for t=0 in HIFI
    return a+b*norm.cdf(x,n,dn)

def get_indices(func):
    """
    input 
      func is a user string function, e.g. 'p[0]*p[2]'
    output
      list of (string) indices found in the string, 
      in between 'p[' and ']', e.g. '0','2'
        used by tools int2_method_key, int2_global_method_key, diffunc (grad), relocate_keys, read_pardict_from_widgets
    """
    from mujpy.tools.tools import findall
    return [func[i:j] for (i,j) in zip([k+1 for k in findall('[',func)],[l for l in findall(']',func)])]

def diffunc(func):
    """
    input user function of the form 'p[0]*(1-p[1])'
          up to functions of three parameters (this could be easily extended)
    output a list of its derivatives [with respect to 'p[0]' and 'p[1]']
    and the list of their indices, [0,1]
        grad method used in tools
    """
    from mujpy.tools.tools import get_indices
    from sympy import symbols,diff,sympify,simplify
    from sympy import sin,cos,exp, sqrt,atan,pi
    # identify variables
    # first identify indices
    func = func.replace('abs','Abs').replace('arctan','atan')
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
                
def translate(lmin,function_in):
    """
    input: 
        lmin: list of minuit indices replacement, one for each dashboard index, -1 is blank
        function_in: single function string, of dashboard index nint, to be translated
    output: 
        function_out: single translated function
    ::
 
       translate([0,0,1,2],'p[0]*2+p[3]') yields 'p[0]+2*p[2]'

    e.g. if parameter 1 is shared with parameter 0, the minuit parameter index 3
    will be translated to 2  
        used in int2_method_key and min2int to replace parameter indices contained in function[nint] e.g.
    """
    from copy import deepcopy
    from mujpy.tools.tools import findall
    # print(' nint = {}, lmin = {}\n{}'.format(nint,lmin,function_in))
    function_out = deepcopy(function_in)
    # search for integers between '[' and ']'
    start = [i+1 for i in findall('[',function_out)]  
    # finds index of number after all occurencies of '['
    stop = [i for i in findall(']',function_out)]
    # same for ']'
    nints = [function_out[i:j] for (i,j) in zip(start,stop)] 
    # this is a list of strings with the numbers to be replaced
    try: 
        nmins = [lmin[int(function_out[i:j])] for (i,j) in zip(start,stop)]
    # replacements integers
        for lstr,m in zip(nints,nmins):
            function_out = function_out.replace(lstr,str(m))    
        return function_out
    except Exception as err: # supposes that err is list index out of range
        print('Exception: {}'.format(err))
        print("If Exception is 'list index out of range' then probably")
        print("**** One or more model parameters with flag  '='")
        print("     point to other model parameters with an '=' flag")
        print("**** Check function syntax and model layout")
        return False

def value_error(value,error):
    """
    value_error(v,e)
    returns a string of the format v(e) 
        used in mufit
    """
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
    
def value_error_csv(value,error):
    """
    value_error_csv(v,e)
    returns a string of the format v, e, with the correct precision
        used in tools.print_csv_components, used in mufit
    """
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
        form += 'f},{:.'
        form += '{}'.format(exponent)
        form += 'f},".format(value,error)'
    else:
        if abs(value)<eps:
            form = '"0,0,"' # too small both
        else:
            form = '"{},0,".format(value)' # too small error
    return eval(form)

def version_flag(mufit_method):
    """
    generate a label with the fit type A1-C2
    to be included in the log and csv filename version 
        used in mufit
    """
    if mufit_method.A1() or mufit_method.A1_calib():
        return "A1"
    elif mufit_method.A20() or mufit_method.A20_calib():
        return "A20"
    elif mufit_method.A21() or mufit_method.A21_calib():
        return "A21"
    elif mufit_method.B1() or mufit_method.B1_calib():
        return "B1"
    elif mufit_method.B20() or mufit_method.B20_calib():
        return "B20"
    elif mufit_method.B21() or mufit_method.B21_calib():
        return "B21"
    elif mufit_method.C1() or mufit_method.C1_calib():
        return "C1"
    elif mufit_method.C2() or mufit_method.C2_calib():
        return "C2"
    else:
        return "_"

def transliterate(string_in):
    """
    input: function or function_multi string_in, e.g. '(1-p[2])*p[3]' 
        searches 'p[d]' and substitutes progressively with 'x', 'y', 'z', 'u', 'v', 'w'
    output: string_out '(1-x)*y' and list indices = ['2','3'] of unique indices (strings)
        used by tools error_propagation with sympy
    """
    import re
    symbols = ['x','y','z','u','v','w']
    indices = re.findall(r'p\[(\d+)\]',string_in)
    indices = list(dict(zip(indices,indices)))
    nind = len(indices)
    for k,x in enumerate(symbols[:nind]):
        string_in = re.sub(r'p\['+indices[k]+r'\]',x,string_in)
    return string_in, indices

def litertransate(string_in,indices):
    """
    input:  string error_propagation formula in terms of x, y, z, ...
            indices, list ['i',...] of 'p[i]' to be substituted to the above
    output: substituted string
        used by tools error_propagation with sympy
    """
    import re
    symbols = ['x','y','z','u','v','w']
    nind = len(indices) 
    for k,x in enumerate(symbols[:nind]):
        string_in = re.sub(x,'p['+indices[k]+']',string_in)
    return string_in

def error_propagation(string_in):
    """
    input: string from 'function' of 'function_multi' dictionary, e.g. '(1-p[2])*p[3]'
    output: string for error propagation, e.g. "sqrt((p[3]*e[2])**2+((1-p[2])*e[2])**2)"
            calculated with transliterate, sympy, litertransate
        used by  min2int_global, used by mufit
    """
    from mujpy.tools.tools import transliterate, litertransate
    from sympy import sympify, diff, symbols
    x,y,z,u,v,w = symbols('x y z u v w')
    string, indices = transliterate(string_in)
    expression = sympify(string)
    nind = len(indices)
    variables = [x,y,z,u,v,w]
    string_out = 'sqrt('
    for k,x in enumerate(variables[:nind]):
        dx = str(diff(expression,x))
        # string_out += '('+re.sub(r'p\['+indices[k]+r'\]',x,string_in)+'*e['+str(indices[k])+'])**2'
        string_out += '(('+litertransate(dx,indices)+')*e['+str(indices[k])+'])**2'
        string_out += '+'
    string_out = string_out[:-1]+')'
    return string_out

def limits(string):
    """
    translates string 'x,y' to list of floats or None 

        used by mudashed.dump_dashed
    """
    nones = string.count('None')
    return  [None, None] if nones == 2 else [None, float(string.split(',')[1])] if nones == 1 and string.index('None')== 0 else [float(string.split(',')[0]),None] if nones ==1 else [float(s) for s in string.split(',')]

"""
 REMEMBER: TOOLS METHODS DO NOT NEED TO IMPORT OTHER TOOLS METHODS!
     MAY REMOVE ALL from tools.tools import ...
     some are already
     leaving import may be useful to grep if used by other tools methods
"""
