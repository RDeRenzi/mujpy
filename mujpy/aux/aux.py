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
# MUGUI AUX
##############

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
    model = self.model_components, in mugui, a list of complex dictionaries
        [{'name':'da', 'pars':{'name':'lpha',...},
         {'name':'mg', 'pars':{       ...       }]
    kin is the index of a dashboard parameter (kint)
    '''
    from numpy import array, cumsum, argmax
    
    ncomp = len(model) # number of components in model
    npar = array([len(model[k]['pars']) for k in range(ncomp)]) # number of parameters of each component
    npars = cumsum(npar)
    return argmax(npars>kin)
        
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

    looks for 'l','l:m','l+n+m' 
    where l, m, n are integers

    rejects all other characters

    returns a list of lists of integer
    '''
    s = []
    try:
    # systematic str(int(b[])) to check that b[] ARE integers
        for b in string.split(','): # csv
            kcolon = b.find(':') # ':' and '+' are mutually exclusive
            kminus = b.find('-') # '-' is 'equivalent to ':' but only one of the two is admitted
            kplus = b.find('+')
            if kcolon>0: # value might be a range
                if int(b[:kcolon])>int(b[kcolon+1:]):
                    errmsg='typo!'
                for j in range(int(b[:kcolon]),int(b[kcolon+1:])+1):
                    s.append([str(j)]) # strings
            elif kminus>0:
                if int(b[:kminus])>int(b[kminus+1:]):
                    errmsg='typo!'
                for j in range(int(b[:kminus]),int(b[kminus+1:])+1):
                    s.append([str(j)]) # strings
            elif kplus>0:
                ss = []
                k0 = 0
                while kplus>0: # str(int(b[]))
                    ss.append(int(b[k0:kplus])) 
                    k0 = kplus+1
                    kplus = b.find('+',k0)
                ss.append(int(b[k0:]))
                s.append([str(q) for q in ss])
            else:# value should be an integer
                int(b) # produces an Error if b is not an integer
                s.append([b]) # added as a string
    except Exception as e:
        s = []
        errmsg = e
    if 'errmsg' not in locals():
        errmsg = None
    return s, errmsg

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

def get_title(run,notemp=False,nofield=False):
    '''
    form standard psi title
    '''
    if notemp:
        return '{} {} {}'.format(run.get_sample(),run.get_field(),run.get_orient())
    elif nofield:
        return '{} {} {}'.format(run.get_sample(),run.get_orient(),run.get_temp())
    return '{} {} {} {}'.format(run.get_sample(),run.get_field(),run.get_orient(),run.get_temp())    

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

def muzeropad(runs):
    '''

    runs is a string containing the run number
    Utility of the suite tab, not a method.

    Future::

    1) determine how many leading zeros for padding
       read a file from data dir
       check number of digits before '.'
       count number of digits in run
       zero pad

    now::

    0) assumes no more than 4 digits 
       left pads zeros to 4 digits
       e.g. 32 becomes 0032
    
    '''
    zeros='0000'
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
    in_path = filedialog.askopenfilename(filetypes=((spc,spcdef),('all','*.*')))
    os.chdir(here)
    return in_path


def plot_parameters(nsub,labels,fig=None): 
    '''
    standard plot of fit parameters vs B,T (or X to be implemente)
    input
       nsub<6 is the number of subplots
       labels is a dict of labels, e.g. {title:self.title, xlabel:'T [K]', ylabels: ['asym',r'$\lambda$',r'$\sigma$,...]}
       fig is the standard fig e.g self.fig_pars
    output the ax array on which to plot 
       one dimensional, from top to bottom 
                        and again for two columns
    e.g. two asymmetry parameters are both plotfal=1 and are plotted in ax[0]
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

    fig.canvas.set_window_title('Fit parameters') # the title on the window bar
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

def rebin(x,y,strstp,pack,e=None):
    '''
    * x,y[,e] are 2D arrays to be rebinned
    * pack is the rebinning factor, e.g it returns::

        xb = array([x[0:pack].sum()/pack])

    * strstp = [start,stop] is a list of indices

        rebinning is done on slices of x,y[,e]
        such as x[start:stop]

    use either::

        xb,yb = rebin(x,y,strstp,pack)

    or::

       xb,yb,eyb = rebin(x,y,strstp,pack,ey) # the 5th is y error

    Works also with pack = 1

    Works for 1d array x and 2D ndarrays y, ey returning 2D arrays
    xb, yb, eb, e.g.::
 
         xb.shape = (1,25000), 
         yb.shape = (nruns,25000) 

    '''
    from numpy import floor, sqrt,empty, array
    start,stop = strstp
    if pack==1: # no rebinning, just a slice of 2D arrays
        xx = x[:,start:stop] 
        yy = y[:,start:stop]
        # print('in rebin: shape xx {}, yy {}'.format(xx.shape,yy.shape)) 
        if e is None:
            return xx,yy
        else:
            ee = e[:,start:stop]
            return xx,yy,ee
    else:
        m = int(floor((stop-start)/pack)) # length of rebinned xb
        mn = m*pack # length of x slice 
        xx =x[:,start:start+mn] # slice of the first 2D array
        xx = xx.reshape(m,pack) # temporaty 2d array
        xb = array([xx.sum(1)/pack]) # rebinned first ndarray
        nruns = y.shape[0] # number of runs
        yb = empty((nruns,m))
        if e is not None:
            eb = empty((nruns,m))
        for k in range(nruns): # each row is a run
            yy = y[k][start:start+mn]  # slice row
            yy = yy.reshape(m,pack)  # temporaty 2d
            yb[k] = yy.sum(1)/pack # rebinned row
            if e is not None:
                ey = e[k][start:start+mn]   # slSice row
                ey = ey.reshape(m,pack)  # temporaty 2d
                eb[k] = sqrt((ey**2).sum(1))/pack  # rebinned row
        if e is not None:
            return xb,yb,eb
        else:
            return xb,yb

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

def spec_prec(a):
    '''
    format specifier precision::

        0 for a > 1.0
        1 for 1.0 > a > 0.1
        2 for 0.1 > a > 0.01 etc.

    '''
    import numpy as np
    return int(abs(min(0.,np.floor(np.log10(a))))) 

def tlog_exists(path,run):
    '''
    check if tlog exists under various known filenames types
    '''
    import os

    filename_psibulk = 'run_'+muzeropad(run)+'.mon' # add definitions for e.g. filename_isis
    ok = os.path.exists(os.path.join(path,filename_psibulk)) # or os.path.exists(os.path.join(paths,filename_isis))
    return ok

def translate(nint,lmin,function):
    '''
    Used in int2_int and min2int to parse parameters contained in function[nint].value e.g.

    ::
 
       p[4]*2+p[7]

    and translate the internal parameter indices 4 and 7 (written according to the gui parameter list order)
    into the corresponding minuit parameter list indices, that skips shared and fixed parameters.

    e.g. if parameter 6 is shared with parameter 4 and parameter 2 is fixed, the minuit parameter indices
    will be 3 instead of 4 (skipping internal index 2) and 5 instead of 7 (skipping both 2 and 6)
    Returns executable formula
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
    for lstr,m in zip(nints,nmins):
        string = string.replace(lstr,str(m))
    return string

def translate_nint(nint,lmin,function):
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
    else: # too small error
        form = 'print("{} +- {:.1e}".format(value,error))'
    return eval(form)
