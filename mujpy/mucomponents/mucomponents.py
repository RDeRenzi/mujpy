from numpy import cos, sin, pi, exp, sqrt, log, real, nan_to_num, inf, ceil, linspace, zeros, empty, ones, hstack, fft, sum, zeros_like, abs, array, where, arctan, pi
from scipy.special import dawsn,erf, j0, j1
from scipy.constants import physical_constants as C
from iminuit.util import make_func_code

class mumodel(object):
    '''
    Defines the components of the fitting model. Provides a chi_square function for Minuit.
    Fields in mT after substitution of self._gamma_mu_ with self._gamma_Mu_MHzper_mT
    '''
    def __init__(self):
        ''' 
        Defines few constants and _help_ dictionary
        '''
        self._radeg_ = pi/180.
        self._gamma_Mu_MHzperT = 3.183345142*C['proton gyromag. ratio over 2 pi'][0]  # numbers are from Particle Data Group 2017
        self._gamma_mu_ = 135.538817
        self._gamma_Mu_MHzper_mT = self._gamma_Mu_MHzperT*1e-3
        self._help_ = {'bl':r'Lorentz decay: $A\exp(-\lambda\,t)$',
                     'bg':r'Gauss decay: $A\exp(-0.5(\sigma\,t)^2)$',
                     'bs':r'stretched decay: $A\exp(-0.5(\Lambda\,t)^\beta)$',
                     'ba':r'Lorentz and Gauss decay: $A\exp(-\lambda\,t)\exp(-0.5(\sigma\,t)^2)$',
                     'da':r'Linearized dalpha correction: $f = \frac{2f_0(1+\alpha/\mbox{dalpha})-1}{1-f_0+2\alpha/dalpha}$',
                     'ml':r'Lorentz decay: $A\cos[2\pi(\gamma_\mu B\, t +\phi/360)]\exp(-\lambda\,t)$',
                     'mg':r'Lorentz and Gauss decay: $A\cos[2\pi(\gamma_\mu B\, t +\phi/360)]\exp(-\lambda\,t)\exp(-0.5(\sigma\,t)^2)$',
                     'mu':r'Gauss decay: $A\cos[2\pi(\gamma_\mu B\, t +\phi/360)]\exp(-0.5(\sigma\,t)^2)$',
                     'ms':r'Gauss decay: $A\cos[2\pi(\gamma_\mu B\, t +\phi/360)]\exp(-(\Lambda\,t)^\beta)$',
                     'jl':r'Lorentz Bessel: $Aj_0[2\pi(\gamma_\mu B\, t +\phi/360)]\exp(-\lambda\,t)$',
                     'jg':r'Gauss Bessel: $A j_0[2\pi(\gamma_\mu B\, t +\phi/360)]\exp(-0.5(\sigma\,t)^2)$',
                     'js':r'Gauss Bessel: $A j_0[2\pi(\gamma_\mu B\, t +\phi/360)]\exp(-(\Lambda\,t)^\beta)$',
                     'fm':r'FMuF: $A/6[3+\cos 2*\pi\gamma_\mu\mbox{dipfield}\sqrt{3}\, t + \
               (1-1/\sqrt{3})\cos \pi\gamma_\mu\mbox{dipfield}(3-\sqrt{3})\,t + \
               (1+1/\sqrt{3})\cos\pi\gamma_\mu\mbox{dipfield}(3+\sqrt{3})\,t ]\exp(-\mbox{Lor_rate}\,t)$', 
                     'kg':r'Gauss Kubo-Toyabe: static and dynamic, in zero or longitudinal field by G. Allodi [Phys Scr 89, 115201]',
                     'kl':r'Lorentz Kubo-Toyabe: static, in zero or longitudinal field by G. Allodi [Phys Scr 89, 115201]',
                     'kd':r'Lorentz Kubo-Toyabe: static, in zero field, multiplied by Lorentz decay, by G. Allodi [Phys Scr 89, 115201]'}
        self._axis_ = None # for self._chisquare_ when set 0,1 sums only along that axis
        self._n0truecomponents_ = 0 # initialize to zero, calib fits set it to 1 to extract alpha 
        # ---- end generic __init__
    
###############################################################################
# General structure:
# _load_... loads data and defines which _add_... method
# _add_... distributes Minuit parameters to components, adds components into model function
# organized by mufit dofit_...which invokes tools methods:
# int2min_... to pass guess values as
#        int2min val err fix lim names pospar are lists
#        int2min_multigroup >> >> 
#        int2min_userpar >> are 
# int2_..._method_key to prepare
#       bndmthd methods that take k = key_as_lambda functions 
#       such that          bndmthd(x,*eval('k(p)'))  is the component value
# fstack to vectorize muligroup, multirun to produce ngroups, nruns component values
#
# calib is dealt differently in fit where 
#                    alpha is a parameter and asymmetry,errors must be recalculated
#                    and in plot where alpha is used to calculate asymmetry, error data
# 
# for data array shape ..............................   see suite asymmetry
#    single asymmetry is an array of shape (nbins,)               _single
#    N multigroup or N multirun are arrays of shape (N,nbins)     _multirun or _multigroup
#    M multigroup and N multirun are arrays of shape (N,M,nbins)  _multirun_multigroup
# fft has its own
#####################################################################################################
# load methods:
#  _load_                   : single run, single group 
#  _load_calib_             : single run single group with recalculation of asymmetry and errors
#                             _add_ = _add_calib_, Minuit must call _chisquare_calib_
#  _load_multigroup_        : single run multigroup global (a unique chisquare)
#  _load_calib_multigroup_  : single run multigroup global with recalculation of asymmetry and errors
#                             _add_ = _add_calib_multigroup_, Minuit calls _chisquare_calib_  
#  _load_multirun_user      : multi run single group global (a unique chisquare)
#  _load_multirun_grad      : multi run single group global (a unique chisquare)
#####################################################################################################
# different calib asymmetry method: must be on the fly during Minuit, i.e. it must rebin in here
# this interferes in mufitplot, where asymm is already rebinned. Requires two distinct _add_ methods
#
# mufitplot.chi calls _add_norebin_ [= _add_, normal fits] doesn't rebin when _add_ is _add_calib
#####################################################################################################

    def _load_(self,x,y,components,e=1):
        '''
        input: 
            x, y, e are numpy arrays, y, e are 1d 
            e = 1 or missing yields unitary errors 
            components is a list [[method,[key,...,key]],...,[method,[key,...,key]]], 
                produced by int2_method() from mujpy.tools.tools
                where method is an instantiation of a component, e.g. self.ml 
                and value = eval(key) produces the parameter value
        '''
        self._x_ = x
        self._y_ = y       
        self._components_ = components
        self._da_index_ = []
        self._add_ = self._add_single_
        self._add_norebin_ = self._add_single_
        self._ntruecomponents_ = len(components)
        self._ndata = y.size

        self._include_all_() # to rcover from possible fft mode
        try:
            if isinstance(e,int):
                self._e_ = ones((y.shape))
            else:
                if len(y.shape)>1:
                    if e.shape!=y.shape or x.shape[0]!=e.shape[-1]:
                        raise ValueError('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                       y.shape,
                                                                                          e.shape))
                    self._e_ = e     
                elif e.shape!=y.shape or e.shape[0]!=x.shape[0]:
                        raise ValueError('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                       y.shape,
                                                                                       e.shape))          
                else:
                    self._e_ = e
        except ValueError as err:
            return False, err       
        return True, '' # no error

    def _load_calib_(self,x,yf,yb,eyf,eyb,returntup,components): 
        '''
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
            self._add_ = self._add_calib_multigroup_ returntup
            no fft of residues
        '''
        from mujpy.tools.tools import rebin
        self._t_ = x  # pristine time array, unbinned
        self._yf_ = yf 
        self._yb_ = yb
        self._eyf_ = eyf 
        self._eyb_ = eyb
        self._start, self._stop, self._pack = returntup
        self._x_,_,_ = rebin(x,yf,[self._start,self._stop],self._pack,e=eyf)
        self._components_ = components
        self._add_ = self._add_calib_
        self._add_norebin_ = self._add_calib_norebin_
        self._ntruecomponents_ = len(components)
        self._n0truecomponents_ = 1 # used when first component is 'al'
        self._ndata = yf.size

    def _load_multigroup_(self,x,y,components,e=1): 

        '''
        input: 
            x, y, e are numpy arrays, y, e are 2d 
            e = 1 or missing yields unitary errors 
            components is a list [[method,[[key,...],...,[key,...]]],...,[method,[[key...],...,[key,...]]]], 
                produced by int2_method() from mujpy.tools.tools
                where method is an instantiation of a component, e.g. self.ml 
                and value = eval(key) produces the parameter value (the inner key list is for different groups)
            _add_multigroup_ must produce a 2d function f.shape(ngroup,nbins)
            therefore _components_ must be a np.vectorize of ngroup copies of the method 
            method da not allowed here, no need for alpha
            no fft of residues
        '''
        self._x_ = x
        self._y_ = y
        self._components_ = components
        self._ntruecomponents_ = len(components)
        self._add_ = self._add_multigroup_
        self._add_norebin_ = self._add_multigroup_
        self._ndata = y.size

        try:
            if isinstance(e,int):
                self._e_ = ones((y.shape))
            else:
                if len(y.shape)>1:
                    # print('_load_multigroup_ mucomponents debug: x,y,e not e=1')
                    if e.shape!=y.shape or x.shape[0]!=y.shape[-1]:
                        raise ValueError('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                       y.shape,
                                                                                       e.shape))          
                elif e.shape!=y.shape or e.shape[0]!=x.shape[0]:
                    raise ValueError('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                           y.shape,
                                                                                           e.shape))          
            self._e_ = e
        except ValueError as e:
            return False, e       
        return True, '' # no error

    def _load_calib_multigroup_(self,x,yf,yb,eyf,eyb,returntup,components,e=True): 
        '''
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
        '''
        from mujpy.tools.tools import rebin
        # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
        # self._ntruecomponents_ = number of components apart from dalpha 
        self._t_ = x # must create a protected pristine array, self._x_ is recalculated fromn it at each call 
        self._yf_ = yf 
        self._yb_ = yb
        self._eyf_ = eyf
        self._eyb_ = eyb
        self._start, self._stop, self._pack = returntup
        self._x_,_,_ = rebin(x,yf,[self._start,self._stop],self._pack,e=eyf)
        self._components_ = components
 
        self._add_ = self._add_calib_multigroup_
        self._add_norebin_ = self._add_calib_multigroup_norebin_
        self._n0truecomponents_ = 1
        self._ntruecomponents_ = len(components)
        self._ndata = yf.size
        # print('debug mucomponents._load_calib_multigroup_: components = {}'.format(self._components_[self._n0truecomponents_:]))
#        for k, val in enumerate(components):
#            self._ntruecomponents_ += 1
#            # print('mucomponents _load_calib_multigroup_ debug: val = {}'.format(val))
#            self._components_.append(val) # store again [method, [key,...,key]], ismin
        if e:
            self._e_ = eyf-eyf # zeros 
        return True, '' # no error

        
    def _load_multirun_user_(self,x,y,components,e=1):
        '''
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
        '''
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

    def _load_multirun_grad_(self,minuit_ordered_grad_list):
        '''
        minuit_ordered_grad_list is
        self._glocals_ to loop 
            over minuit parameters and calculate gradient components, 
                over runs, whose number n = y.shape[0], for globals
                on single run, for locals
                    over component parameters that contain that minuit parameter, i
                    ncluding the derivative of the user function.
        constructed in mujpy.tools.tools  int2_multirun_grad_method_key
        '''
        self._minuit_grad_list_ = minuit_ordered_grad_list
#        for m,grad_list in enumerate(minuit_ordered_grad_list): 
#            for [k,n,j,_,_] in grad_list:
#                print('debug mucomponents _load_multirun_grad_: min = {}, k,n,j = {},{},{}'.format(m,k,n,j))

#####################################################################################################
# add methods:
# the 'time' array is x, not self._x_ because they are invoked by plot with different t vectors
#  _add_single_           : single run, single group 
#  _add_calib_            : single run single group with recalculation of asymmetry and errors
#                           first _load_calib_single_data_, Minuit must call _chisquare_calib_
#  _add_calib_norebin_    : single run single group asymmetry and errors norebin, for mufitplot
#  _add_multigroup_       : single run multigroup global (a unique chisquare)
#  _add_calib_multigroup_ : single run multigroup global with recalculation of asymmetry and errors
#                           first _load_calib_multigroup_data_, FMinuit must call _chisquare_calib_  
#  _add_calib_multigroup_norebin_ 
#                         : single run multigroup global asymmetry and errors norebin, for mufitplot
#  _add_multirun_         : multi run single group global (a unique chisquare)
#  _add_fft_              : single run, single group for partial residues
#####################################################################################################

    def _asymmetry_(self,alpha):
        '''
        input alpha
        output asymm, asyme
        '''
        from numpy import where
        denominator = (self._yf_ + alpha*self._yb_)
        i = where(denominator==0)
        denominator[i] = 1 
        nominator = (self._yf_ - alpha*self._yb_)
        y = nominator/denominator 
        e = 2*alpha/denominator**2*sqrt((self._yb_*self._eyf_)**2 + (self._yf_*self._eyb_)**2) 
        e[where(e==0)] = 1  
        return y,e

    def _add_single_(self,x,*argv): 
        '''
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    val1,val2,val3,val4,val5,val6, ... at this iteration 
                    argv is a list of values [val1,val2,val3,val4,val5,val6, ...]
        _add_single_ DISTRIBUTES THESE PARAMETER VALUES::
              asymmetry fit with fixed alpha
              order driven by model e.g. blml
        '''      
        f = zeros_like(x)  # initialize a 1D array
        p = argv
        for j in range(self._n0truecomponents_,self._ntruecomponents_): 
            component = self._components_[j][0]
            keys = self._components_[j][1] 
            pars = [key(p) for key in keys] 
            f += component(x,*pars)  
        return f

    def _add_calib_(self,x,*argv):

        '''
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    alpha,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [alpha,val1,val2,val3,val4,val5, ...]
        _add_calib_ DISTRIBUTES THESE PARAMETER VALUES::
              asymmetry fit with fitted alpha
              order driven by model e.g. alml
        '''      
        from mujpy.tools.tools import rebin
        alpha = argv[0]
        y,e = self._asymmetry_(alpha)
        t,self._y_,self._e_ = rebin(x,y,[self._start,self._stop],self._pack,e=e)
        return self._add_single_(t,*argv) 

    def _add_calib_norebin_(self,x,*argv):

        '''
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    alpha,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [alpha,val1,val2,val3,val4,val5, ...]
        _add_calib_ DISTRIBUTES THESE PARAMETER VALUES::
              asymmetry fit with fitted alpha, norebin, for mufitplot.chi
              order driven by model e.g. alml
        '''      
        from mujpy.tools.tools import rebin
        alpha = argv[0]
        y,e = self._asymmetry_(alpha)
        return self._add_single_(x,*argv) 

    def _add_multigroup_(self,x,*argv):   
        '''
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    val0,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [val0,val1,val2,val3,val4,val5, ...]
        _add_multigroup_ DISTRIBUTES THESE PARAMETER VALUES::
              asymmetry fit with fixed alpha
              order driven by 
              first global user parameters 
              then local run parameters, first user local and then "~","!" pars in model, e.g. mgbl 
        method is vectorized as a vstack over groups, whose number n = y.shape[0]
        and produce a n-valued np.array function f, f[k,:] for y[k,:],e[k,:] 
        '''      
        from numpy import savez,set_printoptions
        set_printoptions(precision = 2)
        f = zeros((self._y_.shape[0],x.shape[0]))  # initialize a 2D array shape (groups,bins)   
        p = argv 
        
        for method, keys in self._components_[self._n0truecomponents_:]: # works also for calib [self._n0truecomponents_ skips alphas]
            pars = [[key(p) for key in groups_key] for groups_key in keys]
            f += method(x,*pars)  
        return f
          
    def _add_multirun_(self,x,*argv):   
        '''
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    val0,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [val0,val1,val2,val3,val4,val5, ...]
        _add_multirun_ DISTRIBUTES THESE PARAMETER VALUES::
              asymmetry fit with fixed alpha
              order driven by 
              first global user parameters 
              then local run parameters, first user local and then "~","!" pars in model, e.g. mgbl 
        method is vectorized as a vstack over runs, whose number n = y.shape[0]
        and produce a n-valued np.array function f, f[k,:] for y[k,:],e[k,:] 
        '''    
        f = zeros((self._y_.shape[0],x.shape[0]))  # initialize a 2D array shape (groups,bins)   
        p = argv 
        for method, keys in self._components_: #
            pars = [[key(p) for key in groups_key] for groups_key in keys]
            f += method(x,*pars) 
        return f

    def _add_calib_multigroup_(self,x,*argv):   
        '''
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    val0,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [val0,val1,val2,val3,val4,val5, ...]
        _add_calib_multigroup_ REBINS ASYMMETRY AND DISTRIBUTES THESE PARAMETER VALUES::
              asymmetry fit with fixed alpha
              order driven by model e.g. mgbl
        must loop over groups, whose number n = y.shape[0]
        and produce a n-valued np.array function f, f[k] for y[k],e[k] 
        '''      
        from mujpy.tools.tools import rebin
        alpha = array([[argv[0]],[argv[1]]])
        y,e = self._asymmetry_(alpha)
        t,self._y_,self._e_ = rebin(x,y,[self._start,self._stop],self._pack,e=e)
        return self._add_multigroup_(t,*argv)

    def _add_calib_multigroup_norebin_(self,x,*argv):   
        '''
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    val0,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [val0,val1,val2,val3,val4,val5, ...]
        _add_calib_multigroup_ calculates ASYMMETRY AND DISTRIBUTES THESE PARAMETER VALUES::
              asymmetry fit with fixed alpha, norebin, for mufitplot.chi
              order driven by model e.g. mgbl
        must loop over groups, whose number n = y.shape[0]
        and produce a n-valued np.array function f, f[k] for y[k],e[k] 
        '''      
        from mujpy.tools.tools import rebin
        p = argv  # way of passing minuit parameters to key(p) to generate pars
        alpha = array([[argv[0]],[argv[1]]])
        y,e = self._asymmetry_(alpha)
        return self._add_multigroup_(x,*argv)

    def _add_multirun_grad_(self,*argv):
        '''
        input:
            *argv   passed as a variable number of parameter values 
                    val0,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [val0,val1,val2,val3,val4,val5, ...]
           (here independent variable is self._x_ time array by default, no plots)
        output:
            np.array(grad) 
              whose m-th value is the chisquare gradient 
              with respect to internal minuit parameter p[m]
        requires previous calls to _load_data_multirun_grad_  to define self._minuit_grad_list_
                       both constructed in int2_multirun_grad_method_key() from mujpy.tools.tools
        self._gradients_ used to store f values and g derivatives of all components, to minimize numpy array calculations;
        self._glocals_ to loop 
            over minuit parameters and calculate gradient components, 
                over runs, whose number n = y.shape[0], for globals
                on single run, for locals
                    over component parameters that contain that minuit parameter, i
                    including the derivative of the user function.
        _add_multirun_grad_ (user is implicit, no grads in sequential fit) 
        '''    
        p = argv 
        grad = zeros(len(p))
        dcdf = 2*(self._add_(self._x_,*argv) - self._y_)/self._e_**2
        # print('debug grad: shape dcdf = {}'.format(dcdf.shape))
        pars = [[[key(p) for key in component_run]  for component_run in component_runs] for _,component_runs in self._components_]
        for m, grad_list in enumerate(self._minuit_grad_list_):
            gg = zeros((self._y_.shape[0],self._x_.shape[0]))
            for [k,n,j,dkndj,djdm] in grad_list:
                par = pars[n][k]
                dk,dj = dkndj(self._x_,*par),djdm(p)
                #print ('debug grad: m= {}, k = {}, n = {}, j = {}, par = {}'.format(m,k,n,j,    par))
                #print ('debug: dkddj shape = {}, djdm = {}'.format(dk.shape,dj))
                gg[k] +=  dk*dj
            grad[m] = sum(dcdf*gg,axis=None)           
        return grad     

    def _add_fft_(self,x,y,*argv,calib = False):
        '''
        input:
            x time array, 1d
            *argv Minuit parameters
          Called only externally. Produces f for the fft of residues in mufitplot::
          
            y - f
          
          Components can be selectively added to f
          i.e. subtracted in the residues::
   
            f += method(x,*pars) if self._include_components[j] else 0. 
            
          For the time being only single run single group (f is 1d)
        '''
        f = zeros_like(y)  # initialize a 1D array
        for k,p in enumerate(argv):
            # print('mucomponent _add_fft_ debug: p {}'.format(p[-1]))
            ntruecomp = 1 if calib else 0
            # print('add_single mucomponents debug: p = {}'.format(p))
            for j in range(ntruecomp,self._ntruecomponents_): # all components in model excluding da
                component = self._components_[j][0]
                keys = self._components_[j][1] 
                # print('mucomponent _add_fft_ debug: component {}, keys {}'.format(component,keys))
                # print('add_fft_ mucomponents debug: keys = {}'.format(keys))
#                for key in keys:
#                    print('mucomponent add_fft_ debug: key(p) = {}'.format(key(p)))
                pars = [key(p) for key in keys] # NEW! spedup, evaluates p[1], p[2] etc.
                # print('y:{},x:{},f:[]'.format(self._y_.shape,x.shape,f.shape))
                # print('pars = {}'.format(pars))
                # print('f.shape = {}, zeros.shape = {}'.format(f.shape,zeros_like(x).shape))
                f[k,:] += component(x,*pars) if self._include_components[j] else 0. # new 2.0
                                     # must contain x, for plot x != self._x_
            # remember *p.comp means 'pass as many arguments as required by component, exausting the list p_comp'
#        if self._da_index_:  # linearized correction 
#            dalpha = p[self._da_index_-1]
#            dada = dalpha/self._alpha_
#            f = ((2.+dada)*f-dada)/((2.+dada)-dada*f) if self._include_da else f
        return f

    def _fft_init(self,include_components):
        '''
        input:
          include_components 
                True to subtract in residues  asymm - f
        used to generate a partial residue for FFT
        '''
        self._include_components = include_components 

    def _include_all_(self):
        '''
        reset to normal fit mode (initially of after fft)
        '''
        self._include_components = [True]*self._ntruecomponents_

    def al(self,x,α):
        '''
        alpha calibration 
        x [mus], α
        x dummy, for compatibility
        '''
        # empty method  (could remove x from argument list ?)
        # print('al = {}'.format(α))
        return []
        al.func_code = make_func_code(["α"])                
           
    def bl(self,x,A,λ): 
        '''
        Lorentzian decay, A*exp(-x*λ)
        x [mus], A, λ [mus-1]
        '''
        # x need not be self.x (e.g. in plot)
        # λ = -87. if λ < -87. else λ
        return A*exp(-x*λ)
        bl.func_code = make_func_code(["A","λ"])

    def _grad_bl_0_(self,x,A,λ): 
        '''
        derivative of bl with respect to A in terms of self.bl
        x [mus], A, λ [mus-1]  
        '''
        return self.bl(x,A,λ)/A

    def _grad_bl_1_(self,x,A,λ): 
        '''
        derivative of bl with respect to λ in terms of self.bl
        x [mus], A, λ [mus-1]  
        '''
        return -x*self.bl(x,A,λ)

    def bg(self,x,A,σ): 
        '''
        Gaussian decay, A*exp(-0.5*(x*σ)**2)
        x [mus], A, σ [mus-1] (positive parity)
        '''
        # x need not be self.x (e.g. in plot)        
        return A*exp(-0.5*(x*σ)**2)
        bg.func_code = make_func_code(["A","σ"])

    def _grad_bg_0_(self,x,A,σ): 
        '''
        derivative of bg with respect to A in terms of self.bg
        x [mus], A, σ [mus-1]  
        '''
        return self.bg(x,A,σ)/A
        
    def _grad_bg_1_(self,x,A,σ): 
        '''
        derivative of bg with respect to σ in terms of self.bg
        x [mus], A, σ [mus-1]  
        '''
        return -x**2*σ*self.bg(x,A,σ)

    def ba(self,x,A,λ,σ): 
        '''
        Lorentzian times Gaussian decay, A*exp(-x*λ)*exp(-0.5*(x*σ)**2)
        x [mus], A, λ [mus-1], σ [mus-1] (positive parity)
        '''
        # x need not be self.x (e.g. in plot)
        return A*exp(-x*λ)*exp(-0.5*(x*σ)**2)
        ba.func_code = make_func_code(["A","λ","σ"])

    def _grad_ba_0_(self,x,A,λ,σ): 
        '''
        derivative of ba with respect to A in terms of self.ba
        x [mus], A, σ [mus-1]  
        '''
        return self.ba(x,A,λ,σ)/A

    def _grad_ba_1_(self,x,A,λ,σ): 
        '''
        derivative of ba with respect to λ in terms of self.ba
        x [mus], A, σ [mus-1]  
        '''
        return -x*self.ba(x,A,λ,σ)

    def _grad_ba_2_(self,x,A,λ,σ): 
        '''
        derivative of ba with respect to σ in terms of self.ba
        x [mus], A, σ [mus-1]  
        '''
        return -x**2*σ*self.ba(x,A,λ,σ)

    def bs(self,x,A,Λ,β): 
        '''
        stretched decay A*exp(-(x*Λ)**β), 
        x [mus], A, Λ [mus-1] (>0), β (>0)
        '''
        # x need not be self.x (e.g. in plot)
        return A*exp(-(x*Λ)**β)
        bs.func_code = make_func_code(["A","Λ","β"])

    def _grad_bs_0_(self,x,A,Λ,β): 
        '''
        derivative of bs with respect to A in terms of self.bs
        x [mus], A, Λ [mus-1] (>0), β (>0)  
        '''
        return self.bs(x,A,Λ,β)/A

    def _grad_bs_1_(self,x,A,Λ,β): 
        '''
        derivative of bs with respect to Λ in terms of self.bs
        x [mus], A, Λ [mus-1] (>0), β (>0)  
        '''
        return -β/Λ*(Λ*x)**β*self.bs(x,A,Λ,β)

    def _grad_bs_3_(self,x,A,Λ,β): 
        '''
        derivative of bs with respect to β in terms of self.bs
        x [mus], A, Λ [mus-1] (>0), β (>0)  
        '''
        return -log(Λ*x)*(Λ*x)**β*self.bs(x,A,Λ,β)

    def ml(self,x,A,B,φ,λ): 
        '''
        precession A cos(2 pi _gamma_Mu_MHzper_mT B x+φ _radeg_) times Lorentzian decay, 
        x [mus], A, B [mT], φ [deg], λ [mus-1]
        '''
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)
        ml.func_code = make_func_code(["A","B","φ","λ"])

    def _derivative_ml_(self,x,A,B,φ,λ): 
        '''
        derivative of mlwith respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times Lorentzian decay
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        '''
        return -A*sin(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)

    def _grad_ml_0_(self,x,A,B,φ,λ): 
        '''
        derivative of ml with respect to A in terms of self.ml and self._derivative_ml_
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        '''
        return self.ml(x,A,B,φ,λ)/A

    def _grad_ml_1_(self,x,A,B,φ,λ): 
        '''
        derivative of ml with respect to B in terms of self.ml and self._derivative_ml_
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        '''
        return -2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_ml_(x,A,B,φ,λ)

    def _grad_ml_2_(self,x,A,B,φ,λ): 
        '''
        derivative of ml with respect to φ in terms of self.ml and self._derivative_ml_
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        '''
        return -self._radeg_*self._derivative_ml_(x,A,B,φ,λ)

    def _grad_ml_3_(self,x,A,B,φ,λ): 
        '''
        derivative of ml with respect to λ in terms of self.ml and self._derivative_ml_
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        '''
        return -x*self.ml(x,A,B,φ,λ)

    def mg(self,x,A,B,φ,σ): 
        '''
        precession A cos(2 pi _gamma_Mu_MHzper_mT B x+φ _radeg_) times Gaussian decay, 
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  (positive parity)
        '''
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-0.5*(x*σ)**2)
        mg.func_code = make_func_code(["A","B","φ","σ"])
        
    def _derivative_mg_(self,x,A,B,φ,σ): 
        '''
        derivative of mg with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times Gaussian decay,
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        '''
        return -A*sin(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-0.5*(x*σ)**2)

    def _grad_mg_0_(self,x,A,B,φ,σ): 
        '''
        derivative of mg with respect to A in terms of self.mg and self._derivative_mg_
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        '''
        return self.mg(x,A,B,φ,σ)/A
        
    def _grad_mg_1_(self,x,A,B,φ,σ): 
        '''
        derivative of mg with respect to B in terms of self.mg and self._derivative_mg_
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        '''
        return 2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_mg_(x,A,B,φ,σ)
        
    def _grad_mg_2_(self,x,A,B,φ,σ): 
        '''
        derivative of mg with respect to φ in terms of self.mg and self._derivative_mg_
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        '''
        return -self._radeg_*self._derivative_mg_(x,A,B,φ,σ)
        
    def _grad_mg_3_(self,x,A,B,φ,σ): 
        '''
        derivative of mg with respect to σ in terms of self.mg and self._derivative_mg_
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        '''
        return -x**2*σ*self.mg(x,A,B,φ,σ)
        
    def mu(self,x,A,B,φ,λ,σ): 
        '''
        precession A cos(2 pi _gamma_Mu_MHzper_mT B x+φ _radeg_) times Gaussian times Lorentzian decays, 
        x [mus], A, B [mT], φ [degrees], 
        λ [mus-1], σ [mus-1]  (positive parity)
        '''
        # x need not be self.x (e.g. in plot)
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)*exp(-0.5*(x*σ)**2)
        mu.func_code = make_func_code(["A","B","φ","λ","σ"])

    def _derivative_mu_(self,x,A,B,φ,λ,σ): 
        '''
        derivative of mu with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times Lorentzian times Gaussian decay,
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        '''
        return -A*sin(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)*exp(-0.5*(x*σ)**2)

    def _grad_mu_0_(self,x,A,B,φ,λ,σ): 
        '''
        derivative of mu with respect to A in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        '''
        return self.mu(x,A,B,φ,λ,σ)/A

    def _grad_mu_1_(self,x,A,B,φ,λ,σ): 
        '''
        derivative of mu with respect to B in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        '''
        return -2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_mu_(x,A,B,φ,λ,σ)

    def _grad_mu_2_(self,x,A,B,φ,λ,σ): 
        '''
        derivative of mu with respect to φ in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        '''
        return -self._radeg_*self._derivative_mu_(x,A,B,φ,λ,σ)

    def _grad_mu_3_(self,x,A,B,φ,λ,σ): 
        '''
        derivative of mu with respect to λ in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        '''
        return -x*self.mu(x,A,B,φ,λ,σ)

    def _grad_mu_4_(self,x,A,B,φ,λ,σ): 
        '''
        derivative of mu with respect to σ in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        '''
        return -x**2*σ*self.mu(x,A,B,φ,λ,σ)

    def ms(self,x,A,B,φ,Λ,β): 
        '''
        precession A cos(2 pi _gamma_Mu_MHzper_mT B x+φ _radeg_) times stretched decay, 
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        '''
        # x need not be self.x (e.g. in plot)
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-(x*Λ)**β)
        ms.func_code = make_func_code(["A","B","φ","Λ","β"])

    def _derivative_ms_(self,x,A,B,φ,Λ,β): 
        '''
        derivative of ms with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times stretched decay
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        '''
        return -A*sin(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-(x*Λ)**β)

    def _grad_ms_0_(self,x,A,B,φ,Λ,β): 
        '''
        derivative of ms with respect to A in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0) 
        '''
        return self.ms(x,A,B,φ,Λ,β)/A

    def _grad_ms_1_(self,x,A,B,φ,Λ,β): 
        '''
        derivative of ms with respect to B in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0) 
        '''
        return -2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_ms_(x,A,B,φ,Λ,β)

    def _grad_ms_2_(self,x,A,B,φ,Λ,β): 
        '''
        derivative of ms with respect to φ in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        '''
        return -self._radeg_*self._derivative_ms_(x,A,B,φ,Λ,β)

    def _grad_ms_3_(self,x,A,B,φ,Λ,β): 
        '''
        derivative of ms with respect to Λ in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        '''
        return -β/Λ*(Λ*x)**β*self.ms(x,A,B,φ,Λ,β)

    def _grad_ms_4_(self,x,A,B,φ,Λ,β): 
        '''
        derivative of ms with respect to β in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        '''
        return -log(Λ*x)*(Λ*x)**β*self.ms(x,A,B,φ,Λ,β)


    def fm(self,x,A,B,λ):
        '''
        FmuF (powder average)
        according to Book  
        x [mus], A, B [mT], λ [mus-1]
        B is Bdip
        '''
        # x need not be self.x (e.g. in plot)
        return A/6.0*(1.+cos(2*pi*self._gamma_Mu_MHzper_mT*B*x)+
               2.*(cos(pi*self._gamma_Mu_MHzper_mT*B*x)+
                   cos(3*pi*self._gamma_Mu_MHzper_mT*B*x) ))*exp(-x*λ)
        fm.func_code = make_func_code(["A","B","λ"])

    def jl(self,x,A,B,φ,λ): 
        '''
        Bessel j0 precession times Lorentzian decay, 
        x [mus], A, B [mT], φ [degrees], λ [mus-1]
        '''
        # x need not be self.x (e.g. in plot)
        return A*j0(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)
        jl.func_code = make_func_code(["A","B","φ","λ"])

    def _derivative_jl_(self,x,A,B,φ,λ): 
        '''
        derivative of jl with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        - A J1(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times Lorentzian decay,
        x [mus], A, B [mT], φ [degrees], λ [mus-1]
        '''
        return -A*j1(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)

    def _grad_jl_0_(self,x,A,B,φ,λ): 
        '''
        derivative of jl with respect to A in terms of self.ml and self._derivative_ml_
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        '''
        return self.jl(x,A,B,φ,λ)/A

    def _grad_jl_1_(self,x,A,B,φ,λ): 
        '''
        derivative of jl with respect to B in terms of self.ml and self._derivative_ml_
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        '''
        return -2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_jl_(x,A,B,φ,λ)

    def _grad_jl_2_(self,x,A,B,φ,λ): 
        '''
        derivative of jl with respect to φ in terms of self.ml and self._derivative_ml_
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        '''
        return -self._radeg_*self._derivative_jl_(x,A,B,φ,λ)

    def _grad_ml_3_(self,x,A,B,φ,λ): 
        '''
        derivative of ml with respect to λ in terms of self.ml and self._derivative_ml_
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        '''
        return -x*self.jl(x,A,B,φ,λ)


    def jg(self,x,A,B,φ,σ): 
        '''
        Bessel j0 precession times Gaussian decay, 
        x [mus], A, B [mT], φ [degrees], σ [mus-1] (positive parity)
        '''
        # x need not be self.x (e.g. in plot)
        return A*j0(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-0.5*(x*σ)**2)
        jg.func_code = make_func_code(["A","B","φ","σ"])

    def _derivative_jg_(self,x,A,B,φ,σ): 
        '''
        derivative of jg with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times Gaussian decay,
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        '''
        return -A*j1(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-0.5*(x*σ)**2)

    def _grad_jg_0_(self,x,A,B,φ,σ): 
        '''
        derivative of jg with respect to A in terms of self.mg and self._derivative_mg_
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        '''
        return self.jg(x,A,B,φ,σ)/A
        
    def _grad_jg_1_(self,x,A,B,φ,σ): 
        '''
        derivative of jg with respect to B in terms of self.mg and self._derivative_mg_
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        '''
        return 2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_jg_(x,A,B,φ,σ)
        
    def _grad_jg_2_(self,x,A,B,φ,σ): 
        '''
        derivative of jg with respect to φ in terms of self.mg and self._derivative_mg_
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        '''
        return -self._radeg_*self._derivative_jg_(x,A,B,φ,σ)
        
    def _grad_jg_3_(self,x,A,B,φ,σ): 
        '''
        derivative of jg with respect to σ in terms of self.mg and self._derivative_mg_
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        '''
        return -x**2*σ*self.jg(x,A,B,φ,σ)
        
    def j0(self,x,A,B,φ,λ,σ): 
        '''
        precession A j1(2 pi _gamma_Mu_MHzper_mT B x+φ _radeg_) times Gaussian times Lorentzian decays, 
        x [mus], A, B [mT], φ [degrees], 
        λ [mus-1], σ [mus-1]  (positive parity)
        '''
        # x need not be self.x (e.g. in plot)
        return A*j0(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)*exp(-0.5*(x*σ)**2)
        mu.func_code = make_func_code(["A","B","φ","λ","σ"])

    def _derivative_j0_(self,x,A,B,φ,λ,σ): 
        '''
        derivative of j0 with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times Lorentzian times Gaussian decay,
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        '''
        return -A*j1(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)*exp(-0.5*(x*σ)**2)

    def _grad_j0_0_(self,x,A,B,φ,λ,σ): 
        '''
        derivative of j0 with respect to A in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        '''
        return self.j0(x,A,B,φ,λ,σ)/A

    def _grad_j0_1_(self,x,A,B,φ,λ,σ): 
        '''
        derivative of j0 with respect to B in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        '''
        return -2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_j0_(x,A,B,φ,λ,σ)

    def _grad_j0_2_(self,x,A,B,φ,λ,σ): 
        '''
        derivative of j0 with respect to φ in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        '''
        return -self._radeg_*self._derivative_j0_(x,A,B,φ,λ,σ)

    def _grad_j0_3_(self,x,A,B,φ,λ,σ): 
        '''
        derivative of j0 with respect to λ in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        '''
        return -x*self.j0(x,A,B,φ,λ,σ)

    def _grad_j0_4_(self,x,A,B,φ,λ,σ): 
        '''
        derivative of j0 with respect to σ in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        '''
        return -x**2*σ*self.j0(x,A,B,φ,λ,σ)


    def js(self,x,A,B,φ,Λ,β): 
        '''
        Bessel j0 precession times stretched decay, 
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        '''
        # x need not be self.x (e.g. in plot)
        return A*j0(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-(x*Λ)**β)
        js.func_code = make_func_code(["A","B","φ","Λ","β"])

    def _derivative_js_(self,x,A,B,φ,Λ,β): 
        '''
        derivative of js with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times stretched decay
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        '''
        return -A*j1(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-(x*Λ)**β)

    def _grad_js_0_(self,x,A,B,φ,Λ,β): 
        '''
        derivative of js with respect to A in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0) 
        '''
        return self.js(x,A,B,φ,Λ,β)/A

    def _grad_js_1_(self,x,A,B,φ,Λ,β): 
        '''
        derivative of js with respect to B in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0) 
        '''
        return -2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_js_(x,A,B,φ,Λ,β)

    def _grad_js_2_(self,x,A,B,φ,Λ,β): 
        '''
        derivative of js with respect to φ in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        '''
        return -self._radeg_*self._derivative_js_(x,A,B,φ,Λ,β)

    def _grad_js_3_(self,x,A,B,φ,Λ,β): 
        '''
        derivative of ms with respect to Λ in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        '''
        return -β/Λ*(Λ*x)**β*self.js(x,A,B,φ,Λ,β)

    def _grad_js_4_(self,x,A,B,φ,Λ,β): 
        '''
        derivative of js with respect to β in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        '''
        return -log(Λ*x)*(Λ*x)**β*self.js(x,A,B,φ,Λ,β)
        
# kubo toyabe and fm gradients not implemented

    def _kg(self,t,w,Δ):
        '''
        auxiliary component for a static Gaussian Kubo Toyabe in longitudinal field, 
        t [mus], w [mus-1], Δ [mus-1], 
        w = 2*pi*gamma_mu*L_field
        The first derivative of dawsn(x) is 1-2*x*dawsn(x)
        '''
        # note that t can be different from self._x_
        Dt = Δ*t
        DDtt = Dt**2
        DD = Δ**2
        sqr2 = sqrt(2)
        argf = w/(sqr2*Δ)
        fdc = dawsn(argf)
        wt = w*t
        if (w!=0): # non-vanishing Longitudinal Field
            Aa = real(exp(-0.5*DDtt + 1j*wt)*dawsn(-argf - 1j*Dt/sqr2) )
            Aa[Aa == inf] = 0 # bi-empirical fix
            nan_to_num(Aa,copy=False) # empirical fix 
            A=sqr2*(Aa + fdc)
            f = 1. - 2.*DD/w**2*(1-exp(-.5*DDtt)*cos(wt)) + 2.*(Δ/w)**3*A
        else:
            f = (1. + 2.*(1-DDtt)*exp(-.5*DDtt))/3.
        return f

    def _kl(self,t,w,Δ):
        '''
        static Lorentzian Kubo Toyabe in longitudinal field, 
        t [mus], w [mus-1], Δ [mus-1], 
        w = 2*pi*gamma_mu*L_field
        '''
        # note that t can be different from self._x_
        Dt = Δ*t
        wt = w*t
        dt = t[1]-t[0]
        Dtt = Δ*t[1:] # eliminate first point when singular at t=0
        wtt = w*t[1:] # eliminate first point when singular at t=0
        if w*Δ: # non-vanishing Longitudinal Field
            if abs(w/Δ)<2e-9:
                f = (1. + 2.*(1-Dt)*exp(-Dt))/3.
            else:
                
                if t[0]: # singularity at t=0
                    c = Δ/wtt**2.*(1+Dtt) 
                    f =append(-2/3*Δ, exp(-Dtt)*(sin(wtt)/wtt*(c-Δ)-c*cos(wtt))) # put back first point
                else: # no singularities
                    c = Δ/wt**2.*(1+Dt)
                    f = exp(-Dt)*(sin(wt)/wt*(c-Δ)-c*cos(wt))
                f = 2*cumsum(f*dt)+1 # simplified integral, accuracy < 1e-3;
        else:
            f = (1. + 2.*(1-Dt)*exp(-Dt))/3.
        return f

    def _kgdyn(self,x,w,Δ,ν,*argv):
        ''' 
        auxiliary dynamization of Gaussian Kubo Toyabe 
        by G. Allodi 
        N: number of sampling points;
        dt: time interval per bin [i.e. time base is t = dt*(0:N-1)]
        w [mus-1], Δ [mus-1], ν [MHz] 
        (longitudinal field freq, Gaussian distribution, scattering frequency 
        % alphaN: [optional argument] weighting coefficient alpha times N. Default=10 
        '''
        alphaN = 10. if not argv else argv[0] # default is 10.
        dt = x[1]-x[0]
        N = x.shape[0] + int(ceil(x[0]/dt)) # for function to include t=0
        Npad = N * 2 # number of total time points, includes as many zeros
        t = dt*linspace(0.,Npad-1,Npad)
        expwei = exp(-(alphaN/(N*dt))*t)

        gg = self._kg(t,w,Δ)*(t < dt*N)  #  padded_KT, here t is not x 
        # gg = 1/3*(1 + 2*(1 - s^2*tt.^2).*exp(-(.5*s^2)*tt.^2)) % 

        ff = fft.fft(gg*expwei*exp(-ν*t)) # fft(padded_KT*exp(-jump_rate*t))
        FF = exp(-ν*dt)*ff/(1.-(1.-exp(-ν*dt))*ff) # (1-jump_rate*dt*ff)  

        dkt = real(fft.ifft(FF))/expwei  # ifft
        dkt = dkt[0:N] # /dkt(1) 

        #if (nargout > 1),
        #   t = t[0:intN-1]
        return dkt
         
    def kg(self,x,A,BL,Δ,ν):
        '''
        Gauss Kubo Toyabe in (fixed) long field, static or dynamic
        x [mus], A, BL [mT], Δ [mus-1] (positive parity), ν (MHz)
        '''
        # x need not be self.x (e.g. in plot)
        N = x.shape[0]
        w = 2*pi*BL*self._gamma_Mu_MHzper_mT
        if ν==0: # static 
           f = self._kg(x,w,Δ) # normalized to 1. In this case t = x
        else :            # dynamic
           # P=[w Δ];
 
           f = self._kgdyn(x,w,Δ,ν)
# function generated from t=0, shift result nshift=data(1,1)/dt bins backward
           dt = x[1]-x[0]
           nshift = x[0]/dt
           Ns = N + ceil(nshift)
           if Ns%2: # odd
               Np = Ns//2
               Nm = -Np
           else: # even
               Np = Ns//2-1
               Nm = -Ns//2
           n = hstack((inspace(0,Np,Np+1),linspace(Nm,-1.,-Nm)))
           f = fft.ifft(fft.fft(f)*exp(nshift*1j*2*pi*n/Ns)) # shift back
        # multiply by amplitude
        f = A*real(f[0:N])
        return f
        kg.func_code = make_func_code(["A","BL","Δ","ν"])

    def kl(self,x,A,BL,Γ):
        '''
        Lorent Kubo Toyabe in (fixed) long field, static 
        x [mus], A, BL [mT], Γ [mus-1] 
        '''
        # x need not be self.x (e.g. in plot)
        # (dynamic makes no sense)
        w = 2*pi*BL*self._gamma_Mu_MHzper_mT
        return A*self._kl(x,w,Γ)
        kl.func_code = make_func_code(["A","BL","Γ"])

    def kd(self,x,A,Δ,λ):
        '''
        Gauss Kubo Toyabe static times Lorentz decay
        x [mus], A, B [T], Δ [mus-1], ν (MHz)
        '''
        # x need not be self.x (e.g. in plot)
        return A*self._kg(x,0,Δ)*exp(-x*λ)
        kd.func_code = make_func_code(["A","Δ","λ"])
        #kd.limits = [[None,None],[0.,None],[None,None]]
        #kd.error = [0.002,0.05,0.05]

    def ks(self,x,A,Δ,Λ,β):
        '''
        Gauss Kubo Toyabe times stretched decay
        x [mus], A, B [T], Δ [mus-1], Λ [mus-1] (>0), β (>0)
        '''
        # x need not be self.x (e.g. in plot)
        return A*self._kg(x,0,Δ)*exp(-(x*Λ)**β)
        ks.func_code = make_func_code(["A","Δ","Λ","β"])
        #kd.limits = [[None,None],[0.,None],[None,None]]
        #kd.error = [0.002,0.05,0.05]

    def _chisquare_(self,*argv):
        '''
        Signature provided at Minuit invocation by 
        optional argument forced_parameters=parnames
        where parnames is a tuple of parameter names::

           e.g. parnames = ('asym','field','phase','rate') 

        Works also for global fits, 
        where sum (...,axis=None) yields the sum over all indices.

        Provides partial chisquares over individual runs or groups if self._axis_ = 1 
        None is default and sum is over all indices::
        ''' 

        # print('_chisquare_ mucomponents debug: {} {} {}'.format(self._x_.shape,self._y_.shape,self._e_.shape))
        from numpy import abs # finfo, where, array, 
        # Mepsi = finfo('d').max/10.
        num = abs(self._add_(self._x_,*argv) - self._y_)
#        self._ndata = self._y_.size # number of data points
        normsquaredev = (num/self._e_)**2
        # divergence = normsquaredev>Mepsi
#        if divergence.any():
#            print('Warning: big numbers in chisquare {}'.format(normsquaredev[divergence]))
        return sum(normsquaredev,axis=self._axis_ )

    def _chisquare_calib_(self,*argv):
        '''
        Signature provided at Minuit invocation by 
        optional argument forced_parameters=parnames
        where parnames is a tuple of parameter names::

           e.g. parnames = ('asym','field','phase','rate') 

        Works also for global fits, 
        where sum (...,axis=None) yields the sum over all indices.

        Provides partial chisquares over individual runs or groups if self._axis_ = 1 
        None is default and sum is over all indices::
        '''
        # print('_chisquare_ mucomponents debug: {} {} {}'.format(self._x_.shape,self._y_.shape,self._e_.shape))
        from numpy import abs 
        # Mepsi = finfo('d').max/10.
        num = abs(self._add_(self._t_,*argv) - self._y_) # self._t_ is transformed in self._x_ by _add_calib_... 
                                                         # which also generates self._y_ self._e_ on the spot
#        self._ndata = self._y_.size
        normsquaredev = (num/self._e_)**2
        # divergence = normsquaredev>Mepsi
#        if divergence.any():
#            print('Warning: big numbers in chisquare {}'.format(normsquaredev[divergence]))
        return sum(normsquaredev,axis=self._axis_ )

    def _grad_chisquare_(self,*argv):
        '''
        option for global multirun fits, 
        where sum (...,axis=None) yields the sum over all indices.

        Provides gradient of chisquare with respect to p, i.e. along the i-th parameter p_i it is
        sum_j 2[y(t_j,p)-y_ej)]/e_j^2 sum_k d y_k(t_j,p) /dp_i
        where j are bins, k are components in the model
        
        The first factor is common to a all grad components
        The second factor mus be selected. 
        y_k may not depend on p_i, dy_k/dp_i = 0
        And if y_k' does depend, p_i will be its l-th parameter, and we must use the l-th component of its gradient
        Need 
        
        ''' 
        # print('_chisquare_ mucomponents debug: {} {} {}'.format(self._x_.shape,self._y_.shape,self._e_.shape))
        from numpy import abs # finfo, where, array, 
        # Mepsi = finfo('d').max/10.
#        self._ndata = self._y_.size
        num = abs(self._add_(self._x_,*argv) - self._y_)
        normsquaredev = (num/self._e_)**2
        # divergence = normsquaredev>Mepsi
#        if divergence.any():
#            print('Warning: big numbers in chisquare {}'.format(normsquaredev[divergence]))
        return sum(normsquaredev,axis=self._axis_ )

    def _chisquare_single_(self,*argv,k=0,l=None):
        '''
        inputs:
            argv ar single run single group fit parameters
            k[, l] are indices of _y_ and _e_ multidimensional arrays
        Used only in mufitplot (is it still?)
        Provides partial chisquares over individual runs or groups
        ''' 
       # print('_chisquare_ mucomponents debug: {} {} {}'.format(self._x_.shape,self._y_.shape,self._e_.shape))
#        self._ndata = self._y_[k,:].size
        if l is None:
            
            return sum(  ( (self._add_(self._x_,*argv) - self._y_[k,:]) /self._e_[k,:])**2 )
        else:
            return sum(  ( (self._add_(self._x_,*argv) - self._y_[k,l,:]) /self._e_[k,l,:])**2 )
            
    from iminuit import Minuit as _M
    _chisquare_.errordef = _M.LEAST_SQUARES

