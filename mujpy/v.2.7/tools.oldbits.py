    
#def _available_gradients_(component):
#    """
#    returns True if the component has an analytic gradient
#    i.e. for component name xx in the mucomponents mumodel class, 
#    a method _grad_xx_ in the same class.
#    """
#    from mujpy.mucomponents.mucomponents import mumodel
#    
#    methods_with_grad = [module[6:8] for module in dir(mumodel()) if module[0:6]=='_grad_']: # magical extraction of component names
#    return component in methods_with_grad
    
def get_fit_range(string):
    """
    unused? 
    transforms a valid string for fit_range into a list of integers
    """

    fit_range = []
    for chan in string.split(','):
        fit_range.append(int(chan))
    return fit_range

def get_totals(suite):
    """
    calculates the grand totals and group totals for a single run 

    input is self.suite of class musuite
    returns strings totalcounts groupcounts nsbin maxbin
    """

    from numpy import array, concatenate
    # called only by self.suite after having loaded a run or a run suite

    ###################
    # grouping set 
    # suite.grouping['forward'] and suite.grouping['backward'] are np.arrays of integers
    # initialize totals
    ###################
    
    for k,grpdict in enumerate(suite.grouping):
        if not k: # k is 0
            gr = concatenate((grpdict['forward'],grpdict['backward']))
        else:
            gr = concatenate((gr,concatenate((grpdict['forward'],grpdict['backward']))))
    ts,gs =  [],[]

    for k,runs in enumerate(suite._the_runs_):
        tsum, gsum = 0, 0
        for j,run in enumerate(runs): # add values for runs to add
            n1 = suite.offset+suite.nt0[0]
            for counter in range(run.get_numberHisto_int()):
                if suite.datafile[-3:]=='bin' or suite.datafile[-3:]=='mdu' or suite.datafile[-4.:]=='root':
                    n1 = suite.offset+suite.nt0[counter] 
                histo = array(run.get_histo_vector(counter,1)).sum() 
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

def int2fft(model):
    """
    deprecate    retrieves which components to subtract for residues
    input: 
        model 
            dashboard["model_guess"] 
    output: 
        fft_subtract: a list of boolean values, one per model component
            fft flag True, component subtracted in residues 
    """
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
    
def model_name(dashboard):
    """
    returns the conventional model name (e.g. 'mgbgbl')

    input the dashboard dictionary structure
        used in tools.plot
    """

    return ''.join([item for component in dashboard["model_guess"] for item in component["name"]])

def name_of_model(model_components,model):
    """
    mudash check if model_components list of dictionaries correstponds to model
    """

    content = []
    for component in model_components:
        content.append(component["name"])
    return True if ''.join(content) == model else False

def create_model(model):
    """
    creates the model_guess section of a dashboard, based on the model string (e.g. 'almg')
    
    prechecked by validmodel
    invokes addcomponent
    """

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
    """
    adds a component dict to the model_guess list

    with keys 'name' and a skeleton 'pardicts' from _available_components_ 
    """
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
        # does notb happen if checked by validmodel
        error_msg = '\nWarning: '+name+' is not a known component. Not added.'
        return {}, error_msg # False error code, message

def component(model,kin):
    """
    returns the index of the component to which parameter k belongs in
    model = self.model_guess, in mugui, a list of complex dictionaries::
            [{'name':'da', 'pardicts':{'name':'lpha',...},
            {'name':'mg', 'pardicts':{       ...       }]
            
    kin is the index of a dashboard parameter (kint)
    """
    from numpy import array, cumsum, argmax
    
    ncomp = len(model) # number of components in model
    npar = array([len(model[k]['pardicts']) for k in range(ncomp)]) # number of parameters of each component
    npars = cumsum(npar)
    return argmax(npars>kin)

def calib(dashboard):
    """
    True if the first component is 'al'
    """
    return dashboard['model_guess'][0]['name']=='al'
         
def run_shorthand(runstrings):
    """
    write the runlist contained in runstrings (suite self.runs produced by derun)
        i.e. a list of lists, with separate run numbers in string format, the inner ones  to be added 
    in a compact string, with space separated notation
    e.g.
    '650:655,675,656:674' 
    """
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

def find_nth(haystack, needle, n):
    """
    Finds nth needle in haystack 

    Returns its first occurrence (0 if not present)

    Used by ?
    """
    start = haystack.rfind(needle)
    while start >= 0 and n > 1:
        start = haystack.rfind(needle, 1, start-1)
        n -= 1
    return start
 
def get_datafile_path_ext(datafile,run):
    """
    datafilename = template, e.g. '/fullpath/deltat_gps_tdc_0935.bin'
    run = string of run digits, e.g. '1001'
    returns '/fullpath/deltat_gps_tdc_1001.bin'
    """
    import os
    path = datafile[:datafile.rfind(os.path.sep)+1] # e.g. /afs/psi.ch/projec/bulkmusr/data/gps/d2022/tdc/', works in  WIN with '\' as separator
    fileprefix = datafile[datafile.rfind(os.path.sep)+1:datafile.rfind('.')]
    ext = datafile[datafile.rfind['.']+1-len(datafile)] # e.g. 'bin' or 'nxs'
    return path, fileprefix, ext
    
def get_group(grouping):
    """
    reverse of get_grouping, 
    input 
        grouping is an np.array of indices of detectors , 0 based 
    output is 
        groupcsv shorthand as in self.group[k]["forward"} or self.group[k]["backward"}
          e.g '1:3,5' or '1,3,5' etc.
    """
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
    """
    estracts parameter name from full parameter name (i.e. name + label)
    for the time being just the first letter
    """

    return fullname[0]
   
def nextrun(datapath):
    """
    assume datapath is path+fileprefix+runnumber+extension
    datafile is next run, runnumber incremented by one
    if datafile exists return next run, datafile
    else return runnumber and datapath
    """
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
    """
    assume datapath is path+fileprefix+runnumber+extension
    datafile is present run
    if datafile exists returns path to datafile
    """
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
    """
    assume datapath is path+fileprefix+runnumber+extension
    datafile is prev run, runnumber decremented by one
    if datafile exists return prev run, datafile
    else return runnumber and datapath
    """
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
 
def get_run_number_from(path_filename,filespecs):
    """
    strips number after filespecs[0] and before filespec[1]
    """
    try:
        string =  path_filename.split(filespecs[0],1)[1]
        run = string.split('.'+filespecs[1],1)[0]
    except:
       run = '-1' 
    return str(int(run)) # to remove leading zeros
    
def p2x(instring):
    """
    replaces parameters e.g. p[2] with variable x2 in string
    returns substitude string and list of indices (ascii)
    """
    import re
    patterna = re.compile(r"p\[(\d+)\]") # find all patterns p[*] where * is digits
    n = patterna.findall(instring) # all indices of parameters
    outstring = instring
    for k in n:
        strin = r"p\["+re.escape(k)+r"\]"
        patternb = re.compile(strin)
        stri = r"x"+re.escape(k)  # variable
        outstring = patternb.sub(stri,outstring)
    return outstring, n
    
def errorpropagate(string,p,e):
    """
    parse function in string 
    
    substitute p[n] with xn, with errors en
    calculate the partial derivative pdn = partial f/partial xn 
    return the sqrt of the sum of (pdn*en)**2
    """
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
    """
    group_calib is the list of group dictionaries
        used in json_name
    """
    shorthand = []
    for group in grouping:
        fwd = '_'.join([str(s+1) for s in group['forward']])
        bkd = '_'.join([str(s+1) for s in group['backward']])
        shorthand.append(fwd+'-'+bkd)
    return '+'.join(shorthand)

def json_name(model,datafile,grouping,version,g=False):
    """
    model is e.g. 'mlmg'
    datafile is e.g. '/afs/psi.ch/bulkmusr/data/gps/d2022/tdc/deltat_gps_tdc_1233.bin'
       must have a single '.'
    grp_calib is the list of dictionaries defining the groups
    g = True for global
    version is a label
    returns a unique name for the json dashboard file
    """    

    from re import findall
    from mujpy.tools.tools import group_shorthand
    run = findall('[0-9]+',datafile)[-1]
    return model+'.'+run+'.'+group_shorthand(grouping)+'.'+version+'.json'
    
def muvaluid(string):
    """
    deprecated, now lists are dealt with in json

    Run suite fits: muvaluid returns True/False
    * checks the syntax for string function 
    corresponding to flag='l'. Meant for pars
    displaying large changes across the run suite,
    requiring different migrad start guesses::

    # string syntax: e.g. "0.2*3,2.*4,20."
    # means that for the first 3 runs value = 0.2,
    #            for the next 4 runs value = 2.0
    #            from the 8th run on value = 20.0

    """
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
    """
    Run suite fits: 

    muvalue returns the value 
    for the nint-th parameter of the lrun-th run
    according to string (corresponding flag='l').
    Large parameter change across the run suite
    requires different migrad start guesses.
    Probably broken!
    """
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

def path_dialog(path,title):
    import tkinter
    from tkinter import filedialog
    import os
    tkinter.Tk().withdraw() # Close the root window
    in_path = filedialog.askdirectory(initialdir = path,title = title)
    
    return in_path

################
# PLOT METHODS #
#   see duplicate in tools.plot
################

def plot_parameters(nsub,labels,fig=None): 
    """
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
         used in v.1
    """
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
    """
    service to animate histograms
    e.g. in the fit tab

    extracted from matplotlib animate 
    histogram example
    """
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
    """
    num is figure number (static, to keep the same window) 
    nrow, ncol number of subplots rows and columns
    kwargs is a dict of keys to pass to subplots as is
    initializes figures when they are first called 
    or after accidental killing
    """
    import matplotlib.pyplot as P
    fig,ax = P.subplots(nrow, ncol, num = num, **kwargs)
    fig.canvas.manager.set_window_title(title)
    return fig, ax            
    
###############
# END OF PLOT #
###############

def slice(y,e,krun,kgroup):
    """
        used in v.2.7 mucomponents
    """
    if (krun,kgroup)==(-1,-1): # to be used for A1
        return y, e
    elif krun == -1:
        if len(y.shape) <= 2:
            return y,e
        else:
            return y[:][kgroup][:],e[:][kgroup][:]
    elif kgroup == -1:
        if len(y.shape) <= 2:
            return y,e
        else:
            return y[krun][:][:],e[krun][:][:]
    else:
        if len(y.shape) == 1:
            return y,e
        elif len(y.shape) == 2:
            if krun == 0:
                return y[kgroup][:],e[kgroup][:]
            else: # must be kgroup = 0 otherwise len(y.shape) is 3
                return y[krun][:],e[krun][:]
        else:
            return y[krun][kgroup][:],e[krun][kgroup][:]

def shorten(path,subpath):
    """
    shortens path

    e.g. path, subpath = '/home/myname/myfolder', '/home/myname'
         short = './myfolder' 
        used in v.1
    """

    short = path.split(subpath)
    if len(short)==2:
        short = '.'+short[1]
    return short

def exit_safe():
    """
    opens an are you sure box?
        used in v.1
    """
    from tkinter.messagebox import askyesno
            
    answer = askyesno(title='Exit mujpy', message='Really quit?')
    return answer
        
def tlog_exists(path,run,ndigits):
    """
    check if tlog exists under various known filenames types
        used in v.1
    """
    import os

    filename_psibulk = 'run_'+muzeropad(run,ndigits)+'.mon' # add definitions for e.g. filename_isis
    ok = os.path.exists(os.path.join(path,filename_psibulk)) # or os.path.exists(os.path.join(paths,filename_isis))
    return ok

def translate_nint(nint,lmin,function): # NOT USED any more?!!
    """
    Used in int2_int and min2int to parse parameters contained in function[nint].value e.g.
    ::
 
       p[4]*2+p[7]

    and translate the internal parameter indices 4 and 7 (written according to the gui parameter list order)
    into the corresponding minuit parameter list indices, that skips shared and fixed parameters.

    e.g. if parameter 6 is shared with parameter 4 and parameter 2 is fixed, the minuit parameter indices
    will be 3 instead of 4 (skipping internal index 2) and 5 instead of 7 (skipping both 2 and 6)
    Returns lmin[nint]
    """
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
  
def hash_translate(kin,hash_id,krun):
    """
    is this used?
    translates indices in functions according to 2026 fits
    """
    kout = []
    for k in kin:
        if int(k) in hash_id:
            kout.append(str(int(k)+len(hash_id)*krun))
        else:
            kout.append(k)
    return kout
    
def results():
    """
    generate a notebook with some results
        never used ?
    """
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
    script = script + "\ninputsuitefile = 'input.suite'\n"
    script = script + r"dashboard = modelname+'.'+re.search(r'\d+', runlist).group()+'.'+groupcalibfile[:groupcalibfile.index('.')]+'.'+version+jsonsuffix"
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
    script = script + "\n# End of suite definition, this suite is a single run"
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
    """
    write x with p significant digits
        used in v.2.7
    """
    from numpy import asarray,where,isfinite,abs,floor,log10,round
    x = asarray(x)
    x_positive = where(isfinite(x) & (x != 0), abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - floor(log10(x_positive)))
    return round(x * mags) / mags

def int2_calib_method_key(dashboard,the_model):
    """
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
    """
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
                string = translate(lmin,pardict['function']) # here is where lmin is used
                # translate substitutes lmin[n] where n is the index read in the function (e.g. p[3])
                keys.append(string) # the function will be eval-uated, eval(key) inside mucomponents
                # isminuit.append(False)
                lmin.append(ntot+1) # an illegal index, was working with lmin.append(0)
            else:# flag[j] == '~' or flag[j] == '!'
                nmin += 1
                keys.append('p['+str(nmin)+']')  # this also needs direct translation                      
                lmin.append(nmin) # 
                # isminuit.append(True)
        method_key.append([bndmthd,keys]) 
    return method_key


def int2min_multirun(dashboard,runs,guess=True):
    """
    deprecated
    input: 
        dashboard
        runs not used, for compatibility
        guess = True (defaut), False for results in nmufitplot
            either dashboard["globpardicts_guess"] if guess = True
            or  dashboard["globpardicts_result"] if guess = False
    output: a list of lists:  
        values: minuit parameter values, either guess of result
        errors: their steps
        fixed: True/False for each
        limits: [low, high] limits for each or [None,None]  
        name: name of parameter 'x_label' for each parameter
        pospar: parameter for which component is positive parity, eg s in e^{-(s*t)^2/2}
    """

    pardicts = dashboard["globpardicts_result"] if "globpardicts_result" in dashboard.keys() else dashboard["globpardicts_guess"]    
    model = dashboard["model_result"] if "model_result" in dashboard.keys() else dashboard["model_guess"]

    positive_parity = ['Δ','σ']                                                    
    #####################################################
    # the following variables contain the same as input #
    # parameters to iMinuit, removing '='s (functions)  #
    #####################################################
                                                        
    val, err, fix, lim = [], [], [], []           
    name = []
    pospar = [] # contains index of positive parity parameters, to rerun with no limits
    hash_parameter_id = []
    # first scan the global and local user parameters 
    
    for k,pardict in enumerate(pardicts):  # scan the model components
        if 'positive_parity' in pardict.keys(): 
            if pardict["positive_parity"]:
                posparc.append(k)
                pardict['limits'][0] = 0. 
            # print('debug tools int2min_multirun: pospar {} lim({}) = {}'.format(pardict["name"],k, pardict['limits']))
        errstd = 'error' if 'error' in pardict.keys() else 'std'
        val.append(pardict['value'])
        name.append(pardict['name'])
        err.append(pardict[errstd])
        if 'limits' in pardict.keys():
            lim.append(pardict['limits'])
        if 'flag' in pardict.keys():
            if pardict['flag'] == '!':
                    fix.append(True)
                elif pardict['flag'] == '~':
                    fix.append(False)
                elif pardict['flag'] == '#':
                    fix.append(False)
 
        else: # no 'flag' assume '~'
            fix.append(False)
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

    print('debug tools int2min_multirun: pospar_local {}\nuser_local = {}'.format(pospar_loc,user_local))
    print('debug tools int2min: runs {}'.format([runs[k][0] for k in range(len(runs))])) 

def int2min_multigroup(dashboard,runs,guess=True):
    """
    Deprecated

    input: 
        dashboard
        runs not used, for compatibility
        guess = True (defaut), False for results in nmufitplot
            either dashboard["globpardicts_guess"] if guess = True
            or  dashboard["globpardicts_result"] if guess = False
    output: a list of lists:  
        values: minuit parameter values, either guess of result
        errors: their steps
        fixed: True/False for each
        limits: [low, high] limits for each or [None,None]  
        name: name of parameter 'x_label' for each parameter
        pospar: parameter for which component is positive parity, eg s in e^{-(s*t)^2/2}
    """

    pardicts = dashboard['globpardicts_guess'] if guess or 'globpardicts_results' not in dashboard.keys() else dashboard['globpardicts_results']
    
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
        if 'flag' in pardict.keys():
            if pardict['flag'] == '!':
                fix.append(True)
            elif pardict['flag'] == '~':
                fix.append(False)
            else:
                return False,_,_,_,_,_,_
        # self.console('val = {}\nerr = {}\nfix = {}\nlim = {}\python list with more repeated valuesncomp name = {},\npar name = {} '.format(val,err,fix,lim,name)) 
    return val, err, fix, lim, name, pospar

def int2_multirun_user_method_key(dashboard,the_model,nruns):
    """
    deprecated
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
    """
    from mujpy.tools.tools import cstack, translate_multirun, set_key#, function_multi_in_components
    from mujpy.tools.tools import get_functions_in
#            self._components_ is a list [[method,[key,...,key]],...,[method,[key,...,key]]], 
#                produced by int2_multirun_user_method_key() from mujpy.tools.tools
#                where method is an instantiation of a component, e.g. self.ml 
#                and value = eval(key) produces the parameter value
    model = dashboard['model_guess']  # guess surely exists, it is a list of component dicts, e.g. for mgbl 2 dicts
    method_key = []
    bndmthd = {} # to avoid same name
    n_locals =  [pardict["local"] for pardict in dashboard["globpardicts_guess"]].count(True)
    n_globals = len(dashboard["globpardicts_guess"])-n_locals
    kloc = n_globals+n_locals
    functions_in = get_functions_in(model,kloc-1)
    functions_out = translate_multirun(functions_in,n_locals,kloc,nruns)   
    

    # print('\n\ndebug tools int2_multirun_user_method_key functions_out = {}'.format(functions_out))
    for j,component in enumerate(model):  # scan the model components (as for the first run)
        name = component['name']
        keys = []
        # this method uses pars, a list of lists (runs) of parameter for this component, obtained by key(p) from minuit p
        
        bndmthd[name] = lambda x,*pars, name=name : cstack(the_model.__getattribute__(name),x,nruns,1,*pars)
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


def get_number_minuit_internal(nruns,n_globals,n_locals,model):
    """
    Deprecated

    uses old n_locals
    """
    k_mint = 0
    for j,component in enumerate(model):  # scan the model components (as for the first run)
        flags = [pardict["flag"] for pardict in component["pardicts"]] # these are the flags in the present component
        for k,flag in enumerate(flags): # as many flags as parameters in component
            if flag!="=":
                k_mint += 1
    return n_globals + nruns*(n_locals + k_mint)
    
def get_functions_in(model,kk):
    """
    Deprecated

    input 
        model = single-run model dashboard dict
        kk = kloc -1, is incremented at each free parameter of the model, so that it scans the internal minuit indices
             for these parameters (ignoring those determined by a user funct e.g. "p[0]*p[4]"
    output 
        functions_in = list of lists, one per component, of user functs, one per parameter, for the single-run model, 
                       all component parameters,  including "~" and "!", are translated to appropriate user funct
    """
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


def checkvalidmodel(name,component_names):
    """
    Deprecated

    checkvalidmodel(name) checks that name   
    ::      A1, B1: 2*component string of valid component names, e.g.
                        'mgmgbl'                
                        following not valid anymore
    ::      or A2, B2: same, ending with 1 digit, number of groups (max 9 groups), 
                        'mgmgml' (2 groups)
    ::      or C1: same, beginning with 1 digit, number of external minuit parameters (max 9)
                        '3mgml' (3 external parameters e.g. A, f, phi)
    ::      or C2: same, both previous options
                        '3mgml2' (3 external parameters, 2 groups)  
    """

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

def int2_calib_multigroup_method_key(dashboard,the_model):
    """
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
    """
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
                string = translate(lmin,pardict['function']) # here is where lmin is used
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

def min2int_multirun(dashboard,p,e,_the_runs_):
    """
    Deprecated

    input:
        dashboard;  globpardicts_guess and model_guess from 
            used only to retrieve "function" or "function_multi" 
            and "error_propagation_multi"
        p,e Minuit best fit parameter values and std
        _the_runs_ = list of run numbers in suite
    output: for all parameters
        names list of lists of parameter names
        pars list of lists ofparameter values
        epars list of lists of parameter errors
    used only in summary_multirun_global that prints name value(error) 
        one or more lines of global user parameters (the first list in the inner lists)
        one line per run local user parameters and local component parameters (the others)
    """

    # 
    # initialize
    #
    names, pars, epars = [], [], []
    nameloc, npars, n_locals = [], -1, 0# inner list, components
    name, par, epar = [], [], []
    for k, pardict in enumerate(dashboard['globpardicts_guess']):
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
    for run in range(nruns):
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
    """
    Deprecated

    input:
        dashboard:  full dashboard 
                    to retrieve "function" and "error_propagation_multi"
        p,e:        Minuit best fit parameter values and std
    output: for all groups, all components, all parameters  
        namesg:     list of lists of list of dashboard parameter names
        parsg:      list of lists of list of dashboard parameter values
        eparsg:     list of lists of list of dashboard parameter errors
    used only in summary_global
    e.g. bgbl for 2 groups yields namesg = [[['bgA0','σ0'],['bgA1','λ1']],[['bgA0','σ0'],['bgA1','λ1']]]
                                                   first group                   second group     
    """

    # 
    # initialize
    #
    from mujpy.tools.tools import function_multi_in_components
    from numpy import cos, sin, tan, sinh, cosh, tanh, log, pi, exp, sqrt, real, abs, arctan

    for k in range(len(p)):
        mask_function_multi = function_multi_in_components(dashboard)
    globpardicts = dashboard['globpardicts_guess']  
    e = [e[k] if pardict['flag']=='~' else 0 for k,pardict in enumerate(globpardicts)] # redundant? e=0 for flag "!"
    # names = [pardict['name'] for pardict in globpardicts]

    model = dashboard['model_guess']
    pardicts = [pardict for component in model for pardict in component['pardicts']]
    ngroups = len(pardicts[mask_function_multi.index(1)]["function_multi"]) # liist.index(1) is the index of the first occurrence
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
                        # print('excepted')
                        epar.append(eval(pardict["function_multi"][l].replace('p','e')))
                else:                
                    par.append(eval(pardict["function"])) # the function will be eval-uated inside mucomponents
                    try:
                        epar.append(eval(pardict["error_propagate"]))
                    except:
                        epar.append(eval(pardict["function"].replace('p','e')))
            pars.append(par) # middel list, model
            names.append(name)
            epars.append(epar)
        namesg.append(names)
        parsg.append(pars)
        eparsg.append(epars) 
    return namesg, parsg, eparsg  # list of list of list of parameter names, values, stds

def len_print_components_multirun(names,values,errors):
	"""
    Deprecated

	input: for a component
		parameter names 
		parameter values 
		parameter errors 
	output:
	    max length of string to print, e.g.
	    "bl.A_fast 0.123(4) bl.λ_fast 12.3(4) bl.σ_fast 0(0)"
	"""

	from mujpy.tools.tools import value_error
	outname = [' '+names[k] for k in range(len(names))]
	outval = [' '+value_error(values[k],errors[k]) for k in range(len(names))]
	maxlen = max(len(max(outname,key=len)),len(max(outval,key=len)))
	return maxlen
    
def print_components_multirun(names,values,errors,maxlen):
	"""
    Deprecated

	input: for a component
		parameter names 
		parameter values 
		parameter errors 
	output:
	    strings to print, e.g.
	    "A.fast    λ.fast    σ.fast"
	    "0.123(4)  12.3(4)   0(0)"
	"""

	from mujpy.tools.tools import value_error
	outnam = [' '+names[k] for k in range(len(names))]
	outnam = [outnam[k]+(maxlen-len(outnam[k]))*' ' for k in range(len(outnam))]
	outval = [' '+value_error(values[k],errors[k]) for k in range(len(names))]
	outval = [outval[k]+(maxlen-len(outval[k]))*' ' for k in range(len(outval))]
	return "".join(outnam), "".join(outval)
	
def modelstrip(name):
    """
    Deprecated

    strips numbers of external parameters at beginning of model name
    """

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

#def globpars(dashboard):
#    """
#    deprecated, not used
#    True if there are globpardicts in the fit dashboard (== mufit.global_fit)
#
#    alias of global type fit A21, B21, C1, C2 
#    """
#
#    return "globpardicts_guess" in dashboard

def translate_multirun(functions_in,n_locals,kloc,nruns):
    """
    deprecated
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
    used in int2_multigroup_method_keyrun_user_method_key and int2_multirun_grad_method_key
    """
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

def minparam2_csv(dashboard,values_in,errors_in,multirun=0):
    """
    deprecated
    not needed go the common way: init_csv_row prepare_csv_row write_csv for all
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
    """
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
                form = ',{:.'+'{}'.format(n1)+'f},{:.'+'{}'.format(n1)+'f}'
                row += form.format(parvalue,parerror)
            for parvalue,parerror in zip(values[0],errors[0]): # globals replicated in every row
                n1 = spec_prec(parerror) # calculates format specifier precision
                form = ',{:.'+'{}'.format(n1)+'f},{:.'+'{}'.format(n1)+'f}'
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
                form = ',{:.'+'{}'.format(n1)+'f},{:.'+'{}'.format(n1)+'f}'
                rows += form.format(parvalue,parerror)
    return rows

    def fstack(f,n,val):
    """
    stack one layer of mumodel functions of any dimension
    """
    from numpy import vstack

        for k in range(1,n):
            f = vstack(f,[self._add_(x,*val[k])])
        return f

 def path_file_dialog(path,spec,root=None):
    """
    totally deprecated, breaks under windows 
    launch tkinter filedialog in path, spec is filename after dot
        used in mudashed
    """

    #from ipywidgets.widgets import FileUpload
    from tkinter import filedialog, Tk
    import os


    # out Output in tab dialogs 
    # observe is 
    # def on_file_upload(c):
    #     with out;
    #     out.clear_output()
    #     if not uploader.value:
    #         return
    #     uploaded_file = uploader.value[0]
    #     file_name = uploaded_file['name']
    #     file_content = uploaded_file['content']

    # uploader = FileUpload(accept=spec, multiple=False) # in tools to be able to select different specs
    #
    try:
        root.deiconify()
    except:
        root = Tk() # Close the root window
        root.geometry("+400+10")
    spc, spcdef = '.'+spec,'*.'+spec
    in_path = filedialog.askopenfilename(initialdir = path, filetypes=((spc,spcdef),('all','*.*')))
    in_path = '' if in_path == () else in_path
    root.withdraw()
    return in_path,root

def make_links(test):
    """
    deprecated
    if getcwd() is writeable, ln -a groups to tests/group and, if test, generate data, fit directories, 
    """

    from mujpy import __file__ as MuJPyName
    from os import getcwd, symlink, access, W_OK, remove, mkdir, listdir, rmdir
    from os.path import join, dirname, isdir, islink, isfile
    from mujpy.tools.tools import can_symlink
    from shutil import copyfile as cp

    startuppath = getcwd()
    writeable = access(startuppath, W_OK)
    ln_cp = symlink if can_symlink() else cp
    if writeable:
        # duplicate grp locally
        grp_dir = join(startuppath,"groups")
        src_dir = join(join(dirname(MuJPyName),"tests"),"groups")
        if not isdir(grp_dir): 
            ln_cp(src_dir,grp_dir)

        if test:
            test = test.upper()
            data_dir = join(startuppath,"data")
            data = "data_gps" if test == 'GPS' else "data_root" if test == 'LEM' else "data_nexus"
            mujpy_data_dir = join(join(dirname(MuJPyName),"tests"),data)
            fit = "fit_gps" if test == 'GPS' else "fit_root" if test == 'LEM' else "fit_nexus"
            fit_dir = join(startuppath,"fit")
            mujpy_fit_dir = join(join(dirname(MuJPyName),"tests"),fit)
            # on Windows data_dir isdir even when copying test data
            if not islink(data_dir) and isdir(data_dir): # a real data directory exists
                return test,None,None
            else:
                if islink(data_dir): 
                    remove(data_dir) # stale link, remove
                if isdir(fit_dir): 
                    for file in listdir(fit_dir):
                        pathfile = join(fit_dir,file)
                        if islink(pathfile) or isfile(pathfile): 
                            remove(pathfile) # stale jsons
                        elif isdir(pathfile):
                            for fil in listdir(pathfile):
                                pathfil = join(pathfile,fil) 
                                remove(pathfil)
                            rmdir(pathfile)
                else: 
                    mkdir(fit_dir)
                if can_symlink():
                    symlink(mujpy_data_dir,data_dir)# ln -s mujpy_data_dir in data_dir
                else: 
                    mkdir(data_dir)
                    for file in listdir(mujpy_data_dir):
                        cp(join(mujpy_data_dir,file),join(data_dir,file))

                for file in listdir(mujpy_fit_dir):
                    pathfile = join(mujpy_fit_dir,file) 
                    if isfile(pathfile): ln_cp(pathfile,join(fit_dir,file)) # ln -s file in fit_dir
    else:
        data_dir = True # test True, False means abort 
    return test,data_dir,writeable # allow check writeable

def tk_choose(text,title,options,root=None):
    from tkinter import Tk, Label, Button, Radiobutton, IntVar
    #    ^ Use capital T here if using Python 2.7
    try:
        root.deiconify()
    except:
        root = Tk() # Close the root window
        root.geometry("+400+10")
    root.title(title)
    Label(root, text=text).pack()
    Button(text="Submit", command=root.destroy).pack()
    v = IntVar()
    for i, option in enumerate(options):
        Radiobutton(root, text=option, variable=v, value=i, command=root.destroy).pack(anchor="w")
    root.mainloop()
    if v.get() == 0: return None
    return options[v.get()]

def tk_error(text,title,root=None):
    """
    popup warning for generic typo
    """

    from tkinter import Tk, Label, Button #messagebox as mb
    try:
        root.deiconify()
    except:
        root = Tk() # Close the root window
        root.geometry("+400+10")
    root.title(title)
    label = Label(root, text = text)
    label.pack()
    button = Button(root, text='OK', width=25, command=root.destroy)
    button.pack()
    root.geometry('600x100+400+10')
    root.mainloop()
    return root

def group_syntax(text,root=None):
    """
    popup warning for group syntax
    """

    from tkinter import Tk, Label, Button #messagebox as mb
    try:
        root.deiconify()
    except:
        root = Tk() # Close the root window
        root.geometry("+400+10") 
    root.title("Watch group syntax")
    label = Label(root, text = text)
    label.pack()
    button = Button(root, text='OK', width=25, command=root.destroy)
    button.pack()
    root.geometry('400x100+400+0')
    root.mainloop()
    return root

def savetests():
    """
    deprecated
    save tests (tests.py, almgml.822.3-4.1_fit.py etc.) to local path

    who invokes this?
    """

    from os import getcwd, symlink, listdir
    from os.path import join, isfile, dirname
    from test.support.os_helper import can_symlink
    from shutil import copyfile as cp
    from mujpy import __file__ as MuJPyName

    ln_cp = symlink if can_symlink() else cp
    here = getcwd()
    test = join(dirname(MuJPyName),'tests')
    for fil in listdir(test):
        file = join(test,fil)
        print('file {}'.format(file))
        if isfile(file) and fil[-3:]=='.py': 
            ln_cp(file,join(here,fil))
            print('ln -s {} ./'.format(file))


