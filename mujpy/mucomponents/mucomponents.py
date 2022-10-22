from numpy import cos, pi, exp, sqrt, real, nan_to_num, inf, ceil, linspace, zeros, empty, ones, hstack, fft, sum, zeros_like, abs, array, where
from scipy.special import dawsn,erf, j0
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
        self._gamma_mu_ = 135.5
        self._gamma_Mu_MHzper_mT = self._gamma_Mu_MHzperT*1e-3
        self._help_ = {'bl':r'Lorentz decay: $\mbox{asymmetry}\exp(-\mbox{Lor_rate}\,t)$',
                     'bg':r'Gauss decay: $\mbox{asymmetry}\exp(-0.5(\mbox{Gau_rate}\,t)^2)$',
                     'bs':r'Gauss decay: $\mbox{asymmetry}\exp(-0.5(\mbox{rate}\,t)^\beta)$',
                     'da':r'Linearized dalpha correction: $f = \frac{2f_0(1+\alpha/\mbox{dalpha})-1}{1-f_0+2\alpha/dalpha}$',
                     'mg':r'Gauss decay: $\mbox{asymmetry}\cos[2\pi(\gamma_\mu \mbox{field}\, t +\mbox{phase}/360)]\exp(-0.5(\mbox{Gau_rate}\,t)^2)$',
                     'ml':r'Gauss decay: $\mbox{asymmetry}\cos[2\pi(\gamma_\mu \mbox{field}\, t +\mbox{phase}/360)]\exp(-\mbox{Lor_rate}\,t)$',
                     'ms':r'Gauss decay: $\mbox{asymmetry}\cos[2\pi(\gamma_\mu \mbox{field}\, t +\mbox{phase}/360)]\exp(-(\mbox{rate}\,t)^\beta)$',
                     'jg':r'Gauss Bessel: $\mbox{asymmetry} j_0[2\pi(\gamma_\mu \mbox{field}\, t +\mbox{phase}/360)]\exp(-0.5(\mbox{Lor_rate}\,t)^2)$',
                     'jl':r'Lorentz Bessel: $\mbox{asymmetry}j_0[2\pi(\gamma_\mu \mbox{field}\, t +\mbox{phase}/360)]\exp(-0.5(\mbox{Lor_rate}\,t)^2)$',
                     'fm':r'FMuF: $\mbox{asymmetry}/6[3+\cos 2*\pi\gamma_\mu\mbox{dipfield}\sqrt{3}\, t + \
               (1-1/\sqrt{3})\cos \pi\gamma_\mu\mbox{dipfield}(3-\sqrt{3})\,t + \
               (1+1/\sqrt{3})\cos\pi\gamma_\mu\mbox{dipfield}(3+\sqrt{3})\,t ]\exp(-\mbox{Lor_rate}\,t)$', 
                     'kg':r'Gauss Kubo-Toyabe: static and dynamic, in zero or longitudinal field by G. Allodi [Phys Scr 89, 115201]',
                     'kl':r'Lorentz Kubo-Toyabe: static, in zero or longitudinal field by G. Allodi [Phys Scr 89, 115201]',
                     'kd':r'Lorentz Kubo-Toyabe: static, in zero field, multiplied by Lorentz decay, by G. Allodi [Phys Scr 89, 115201]'}
        self._axis_ = None # for self._chisquare_ when set 0,1 sums only along that axis
    # ---- end generic __init__
    
###############################################################################
# General structure:
# _load_... loads data and defines which _add_... method
# _add_... distributes Minuit parameters to components, adds components into model function
# calib is dealt differently in fit where 
#                    alpha is a parameter and asymmetry,errors must be recalculated
#                    and in plot where alpha is used to caclulate asymmetry, error data
# fft has its own

    def _load_data_(self,x,y,_components_,e=1):
        '''
        input: 
            x, y, e are numpy arrays, y, e are 1d 
            e = 1 or missing yields unitary errors 
            _components_ is a list [[method,[key,...,key]],...,[method,[key,...,key]]], 
                produced by int2_method() from mujpy.aux.aux
                where method is an instantiation of a component, e.g. self.ml 
                and value = eval(key) produces the parameter value
        '''
        # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
        self._x_ = x
        self._y_ = y        # self._global_ = True if _nglobals_ is not None else False
        self._components_ = []
        self._da_index_ = []
        self._add_ = self._add_single_
        self._ntruecomponents_ = 0
        self._n0truecomponents_ = 0
        

        for k, val in enumerate(_components_):
            if val[0]: # val[0] is directly the method for all components but dalpha
                # print('mucomponents _load_data_ debug: keys = {}'.format(val[1]))
                self._ntruecomponents_ += 1
                self._components_.append(val) # store again [method, [key,...,key]], ismin
            else:  # when the method is da  (val[0] was set to [], i.e. False)
                npar = sum([len(comp[1]) for comp in _components_])
                p = range(npar)
                self._da_index_ = 1+val[1][0](p) # position in minuit parameter list +1 to pass logical test
                # print('_da_index_ = {}'.format(self._da_index_-1))
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
                elif e.shape!=y.shape or e.shape[0]!=x.shape[0]:
                        raise ValueError('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                       y.shape,
                                                                                       e.shape))          
                else:
                    self._e_ = e
                    # print('mucomponents _load_data_ debug: {}'.format(self._e_))
        except ValueError as err:
            return False, err       
        return True, '' # no error

    def _load_calib_single_data_(self,x,yf,yb,bf,bb,yfm,ybm,components): 
        '''
        fit with alpha as free parameter
        input:
           x suite.time, 1d of shape nbins
             see suite.time
           yf,yb 1d array of shape nbins
           bf,bb,yfm,yfb are scalars (see suite.single_for_back_counts)
        components is a list [[method,[key,...,key]],...,[method,[key,...,key]]  where
        arg _components_ is the output of int2_int() from mujpy.aux.aux
        for calib the first method is empty (loke da) and its single parameter is alpha
        '''
        self._x_ = x 
        self._yf_ = yf 
        self._yb_ = yb
        self._bf_ = bf
        self._bb_ = bb 
        self._yfm_ = yfm
        self._ybm_ = ybm
        self._components_ = components
        self._add_ = self._add_calib_single_
        # print('_load_calib_single_data_ in mucomponents debug: self._add_ = {}'.format(self._add_))
        self._ntruecomponents_ = 0
        self._n0truecomponents_ = 1

        for k, val in enumerate(components):
            #if val[0]: # val[0] is directly the method for all components but al
            # however al is taken care of by range(1,self._ntruecomponents_) in _add_calib_single_
            self._ntruecomponents_ += 1 # total number of components including 'al'
        # print('_load_calib_single_data_ mucomponents debug: self._ntruecomponents_ = {}'.format(self._ntruecomponents_))

    def _load_data_calib_(self,x,y,components,e=1): 

        '''
        for normal asymmetry plot of a calib best fit (or guess)
        input: 
            x, y, e are numpy arrays, y, e are 1d 
            e = 1 or missing yields unitary errors 
            _components_ is a list [[method,[key,...,key]],...,[method,[key,...,key]]], 
                produced by int2_method() from mujpy.aux.aux
                where method is an instantiation of a component, e.g. self.ml 
                and value = eval(key) produces the parameter value
        the argument _components_ is stored in self._components_
        '''
        # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
        # self._ntruecomponents_ = number of components apart from dalpha 
        self._x_ = x
        self._y_ = y        # self._global_ = True if _nglobals_ is not None else False
        self._components_ = components
        self._add_ = self._add_single_ # calib_ #  always preceed by _load_data_calib_
        self._ntruecomponents_ = 0
        self._n0truecomponents_ = 1

        for k, val in enumerate(components):
                self._ntruecomponents_ += 1
        self._include_all_() # to recover from possible fft mode
        try:
            if isinstance(e,int):
                self._e_ = ones((y.shape))
            else:
                if len(y.shape)>1:
                    if e.shape!=y.shape or x.shape[0]!=e.shape[-1]:
                        raise ValueError('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                       y.shape,
                                                                                       e.shape))          
                elif e.shape!=y.shape or e.shape[0]!=x.shape[0]:
                        raise ValueError('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                       y.shape,
                                                                                       e.shape))          
                else:
                    self._e_ = e
        except ValueError as e:
            return False, e       
        return True, '' # no error

    def _load_data_multigroup_(self,x,y,components,e=1): 

        '''
        input: 
            x, y, e are numpy arrays, y, e are 2d 
            e = 1 or missing yields unitary errors 
            components is a list [[method,[[key,...],...,[key,...]]],...,[method,[[key...],...,[key,...]]]], 
                produced by int2_method() from mujpy.aux.aux
                where method is an instantiation of a component, e.g. self.ml 
                and value = eval(key) produces the parameter value (the inner key list is for different groups)
            _add_multigroup_ must produce an 2d function f.shape(ngroup,nbins)
            therefore _components_ must be a np.vectorize of ngroup copies of the method 
            method da not allowed here, no need for alpha
            no fft of residues
        '''
        # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
        # self._ntruecomponents_ = number of components apart from dalpha 
        self._x_ = x
        self._y_ = y        # self._global_ = True if _nglobals_ is not None else False
        self._components_ = components
        self._ntruecomponents_ = len(components)
        self._add_ = self._add_multigroup_
        self._n0truecomponents_ = 0
        # print('mucomponents _load_data_multigroup_ debug: {} components, n of parameters/component: {}'.format(len(components),[len(par) for group in components for par in group[1]]))
        try:
            if isinstance(e,int):
                self._e_ = ones((y.shape))
            else:
                if len(y.shape)>1:
                    # print('_load_data_multigroup_ mucomponents debug: x,y,e not e=1')
                    if e.shape!=y.shape or x.shape[0]!=y.shape[-1]:
                        # print('_load_data_multigroup_ mucomponents debug: x,y,e different shape[0]>1')
                        raise ValueError('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                       y.shape,
                                                                                       e.shape))          
                elif e.shape!=y.shape or e.shape[0]!=x.shape[0]:
                    # print('_load_data_multigroup_ mucomponents debug: x,y,e different shape[0]=1')
                    raise ValueError('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                           y.shape,
                                                                                           e.shape))          
            # print('_load_data_multigroup_ mucomponents debug: defining self._e_')
            self._e_ = e
            # print('mucomponents _load_data_multigroup_ debug: self._x_ {}, self._y_ {}, self._e_  {}shape'.format(self._x_,self._y_,self._e_) 
        except ValueError as e:
            return False, e       
        return True, '' # no error

    def _load_data_calib_multigroup_(self,x,yf,yb,bf,bb,yfm,ybm,components,e=1): 

        '''
        input: 
            x, yf, yb  are numpy arrays, yf, yb are 2d
            bf, bb, yfm, ybm are backgrounnd f and b, and its exponential average, f and b 
            e = 1, since errors are calculated (should be simplified)
            components is a list [[method,[[key,...],...,[key,...]]],...,[method,[[key...],...,[key,...]]]], 
                produced by int2_method() from mujpy.aux.aux
                where method is an instantiation of a component, e.g. self.ml 
                and value = eval(key) produces the parameter value (the inner key list is for different groups)
            _add_multigroup_ must produce an 2d function f.shape(ngroup,nbins)
            therefore _components_ must be a np.vectorize of ngroup copies of the method 
            method da not allowed here, no need for alpha
            no fft of residues
        '''
        # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
        # self._ntruecomponents_ = number of components apart from dalpha 
        self._x_ = x
        self._yf_ = yf 
        self._yb_ = yb
        self._bf_ = bf
        self._bb_ = bb 
        self._yfm_ = yfm
        self._ybm_ = ybm
        self._components_ = []
        self._add_ = self._add_calib_multigroup_
        self._ntruecomponents_ = 0
        self._n0truecomponents_ = 1

        for k, val in enumerate(components):
            self._ntruecomponents_ += 1
            # print('mucomponents _load_data_calib_multigroup_ debug: val = {}'.format(val))
            self._components_.append(val) # store again [method, [key,...,key]], ismin
        try:
            if isinstance(e,int):
                self._e_ = ones((yf.shape))
            else:
                if len(yf.shape)>1:
                    # print('_load_data_multigroup_ mucomponents debug: x,y,e not e=1')
                    if e.shape!=yf.shape or x.shape[0]!=yf.shape[-1]:
                        print('_load_data_multigroup_ mucomponents debug: x,y,e different shape[0]>1')
                        raise ValueError('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                       y.shape,
                                                                                       e.shape))          
                elif e.shape!=yf.shape or e.shape[0]!=x.shape[0]:
                    # print('_load_data_multigroup_ mucomponents debug: x,y,e different shape[0]=1')
                    raise ValueError('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                           y.shape,
                                                                                           e.shape))          
            # print('_load_data_multigroup_ mucomponents debug: defining self._e_')
            self._e_ = e
        except ValueError as e:
            return False, e       
        return True, '' # no error

    def _load_data_multigroup_calib_(self,x,y,components,e=1): 

        '''
        for normal asymmetry plot of a calib best fit (or guess)
        input: 
            x, y, e are numpy arrays, y, e are 2d 
            e = 1 or missing yields unitary errors 
            components is a list [[method,[[key,...],...,[key,...]]],...,[method,[[key...],...,[key,...]]]], 
                produced by int2_method() from mujpy.aux.aux
                where method is an instantiation of a component, e.g. self.ml 
                and value = eval(key) produces the parameter value (the inner key list is for different groups)
            _add_multigroup_ must produce an 2d function f.shape(ngroup,nbins)
            therefore _components_ must be a np.vectorize of ngroup copies of the method 
            method da not allowed here, no need for alpha
            no fft of residues
        '''
        self._x_ = x
        self._y_ = y        # self._global_ = True if _nglobals_ is not None else False
        self._components_ = components
        self._add_ = self._add_multigroup_ # single_multigroup_calib_ 
        self.ntruecomp = 1
        self._ntruecomponents_ = 0
        self._n0truecomponents_ = 1

        for k, val in enumerate(components):
                self._ntruecomponents_ += 1
        try:
            if isinstance(e,int):
                self._e_ = ones((y.shape))
                # print('mucomponents _load_data_multigroup_calib_ debug, self._e_ = uno'.format(self._e_))    
            else:
                if len(y.shape)>1:
                    if e.shape!=y.shape or x.shape[0]!=e.shape[-1]:
                        raise ValueError('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                       y.shape,
                                                                                       e.shape))          
                    else:
                        self._e_ = e
                        # print('mucomponents _load_data_multigroup_calib_ debug, self._e_ = e')
                elif e.shape!=y.shape or e.shape[0]!=x.shape[0]:
                        raise ValueError('x, y, e have different lengths, {},{},{}'.format(x.shape,
                                                                                       y.shape,
                                                                                       e.shape))          
        except ValueError as e:
            return False, e 
        return True, '' # no error

#####################################################################################################
# add methods:
# the 'time' array is x, not self._x_ because they are invoked by plot with different t vectors
#  _add_single_ :           single run, single group (still has da)
#  _add_calib_single_ :     single run single group with recalculation of asymmetry and errors
#  _add_single_calib_ :     single run single group for calib plot (normal asymmetry)
#                                  use _add_single_ after _load_data_calib_
#  _add_multigroup_ :       single run multigroup global (for single chisquare)
#  _add_calib_multigroup_ : single run multigroup global with recalculation of asymmetry and errors
#  _add_fft_ : single run, single group for partial residues
#####################################################################################################

    def _add_single_(self,x,*argv): 
        '''
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    val1,val2,val3,val4,val5,val6, ... at this iteration 
                    argv is a list of values [val1,val2,val3,val4,val5,val6, ...]

        _add_ DISTRIBUTES THESE PARAMETER VALUES::

              asymmetry fit with fixed alpha
              order driven by model e.g. blml

        Removed "da" forever - use calib instead

        Removed FFT, see _add_fft_ instead

        '''      

        f = zeros_like(x)  # initialize a 1D array
        p = argv 
        # print('add_single mucomponents debug: p = {}'.format(p))
        for j in range(self._ntruecomponents_): # all components in model excluding da
            component = self._components_[j][0]
            keys = self._components_[j][1] 
            # print('add_single mucomponents debug: keys = {}'.format(keys))
            pars = [key(p) for key in keys] # NEW! spedup, evaluates p[1], p[2] etc.
            # print('y:{},x:{},f:[]'.format(self._y_.shape,x.shape,f.shape))
            # print('pars = {}'.format(pars))
            # print('f.shape = {}, zeros.shape = {}'.format(f.shape,zeros_like(x).shape))
            f += component(x,*pars)  # if self._include_components[j] else 0. # removed 2.0
                                     # must contain x, for plot x != self._x_
            # remember *p.comp means 'pass as many arguments as required by component, exausting the list p_comp'
#        if self._da_index_:  # linearized correction 
#            dalpha = p[self._da_index_-1]
#            dada = dalpha/self._alpha_
#            f = ((2.+dada)*f-dada)/((2.+dada)-dada*f) if self._include_da else f
        return f
        
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

    def _add_calib_single_(self,x,*argv):
        '''
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    alpha,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [alpha,val1,val2,val3,val4,val5, ...]

        _add_calib_single_ DISTRIBUTES THESE PARAMETER VALUES::

              asymmetry fit with fitted alpha
              order driven by model e.g. alml

        NO FFT mode, no check on self._include_components
        '''      
        from numpy import where
        from mujpy.aux.aux import TauMu_mus
                
        f = zeros_like(x)  # initialize a 1D array
        p = argv 
        alpha = p[0]
        #print('_add_calib_single_ debug alpha = {}, p = {}'.format(alpha,p))

        # compute asymmetry and error (needed only by fit, for plot it's a small overhead)
        denominator = (self._yfm_ + alpha*self._ybm_)*exp(-x/TauMu_mus()) # f+b normalization count
        self._y_ = (self._yf_ - alpha*self._yb_ - (self._bf_ - alpha*self._bb_)) / denominator 
        errexp = sqrt(self._yf_ + alpha**2*self._yb_) # equivalent f+b counts
        errexp[where(errexp==0)] = 1  #   set to 1 the minimum error for zero equivalent f+b counts
        self._e_ = errexp / denominator 
                
        for j in range(self._n0truecomponents_,self._ntruecomponents_): # all components in model, excluding alpha
            method = self._components_[j][0]
            keys = self._components_[j][1] 
            pars = [key(p) for key in keys] # NEW! spedup, evaluates p[1], p[2] etc.
            # print('_add_calib_single_ debug y:{},x:{},f:[]'.format(self._y_.shape,x.shape,f.shape))
            # print('_add_calib_single_ debug pars = {}'.format(pars))
            # print('_add_calib_single_ debug f.shape = {}, zeros.shape = {}'.format(f.shape,zeros_like(x).shape))
            f += method(x,*pars)  # must contain x, for plot x != self._x_
            # remember *p.comp means 'pass as many arguments as required by component, exausting the list p_comp'

        return f     

    def _add_multigroup_(self,x,*argv):   
        '''
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    val0,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [val0,val1,val2,val3,val4,val5, ...]

        _add_ DISTRIBUTES THESE PARAMETER VALUES::

              asymmetry fit with fixed alpha
              order driven by model e.g. mgbl
        must loop over groups, whose number n = y.shape[0]
        and produce a n-valued np.array function f, f[k] for y[k],e[k] 
        '''      

        f = zeros((self._y_.shape[0],x.shape[0]))  # initialize a 2D array shape (groups,bins)   
        p = argv 
        
        # self._component_ contains [bndkeys,...,bndkeys], as many as the model components (e.g. 2 for mgbl)
        # bndkeys is [method, [keys_1,keys_2]] if there are 2 groups, keys_i is a list of keys for group i=1,2   
        # such that method(x,*par_i),  produce the additive function component for group i
        # and par_i[k] = eval(keys_i[k])
#        j = -1
#        for method, keys in self._components_:# all components in model 
#            j +=1
#            pars = [[key(p) for key in groups_key] for groups_key in keys]
            # print('mucomponents _add_multigroup_ debug: component {}-th: {}\npars {}'.format(j,method.__doc__,pars))
        j = -1
        for method, keys in self._components_:
            j += 1
        #j in range(self._n0truecomponents_,self._ntruecomponents_): # all components in model excluding da
            #component = self._components_[j][0]
            #keys = self._components_[j][1] # = [keys_1,keys_2,...]
            
            
            # keys = [[p0g0, p0g1,...],[p1g0, p1g1, ..],[p2g0, p2,g1,...]..]
            # print('add_multigroup mucomponents debug: key = {}'.format(keys))
            pars = [[key(p) for key in groups_key] for groups_key in keys]
            # print('mucomponents _add_multigroup_ debug:component {}-th: {}\npars {}'.format(j,method.__doc__,pars))
            f += method(x,*pars)
            # print('mucomponents _add_multigroup_ debug: f[:,0] {}'.format(f[:,0]))
            # f += component(x,*pars) # x is 1d, component is vstacked, with shape (groups,bins) 
            # remember *p.comp means 'pass as many arguments as required by component, exausting the list p_comp'

            # print('add_multigroup mucomponents debug: pars = {}'.format(pars))
            # pars = [[eval(key) for key in groups_key] for groups_key in keys]
            # print('add_multigroup mucomponents debug: y:{},x:{},f:[]'.format(self._y_.shape,x.shape,f.shape))
            # print('add_multigroup mucomponents debug: pars = {}'.format(pars))
            # print('add_multigroup mucomponents debug: f.shape = {}, zeros.shape = {}'.format(
            #                                                         f.shape,zeros_like(x).shape))
#            warn = array(where(abs(f)>1))
#            if warn.size:
#                print('Warning, model is getting too big in {}'.format(warn))
        return f     

    def _add_calib_multigroup_(self,x,*argv):   
        '''
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    val0,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [val0,val1,val2,val3,val4,val5, ...]

        _add_ DISTRIBUTES THESE PARAMETER VALUES::

              asymmetry fit with fixed alpha
              order driven by model e.g. mgbl
        must loop over groups, whose number n = y.shape[0]
        and produce a n-valued np.array function f, f[k] for y[k],e[k] 
        '''      
        from mujpy.aux.aux import TauMu_mus
        from numpy import where,sqrt,exp,array

        f = zeros((self._yf_.shape[0],x.shape[0]))  # initialize a 2D array        
        p = argv 
        # print('mucomponents _add_calib_multigroup_ debug: Minuit p = {}'.format(p))
        alpha = []
        for group in self._components_[0][1]:
            # print('mucomponents _add_calib_multigroup_ debug: group = {}'.format(group))
            key = group[0]
            # print('mucomponents _add_calib_multigroup_ debug: component p = {}'.format(key(p)))
            alpha.append([key(p)])
        alpha= array(alpha)
        # print('mucomponents _add_calib_multigroup_ debug: alpha = {}'.format(alpha))
        #alpha = alpha.transpose() # shape is (ngroups,1)
        # compute asymmetry and error (needed only by fit, for plot it's a small overhead)
        # can multiply 2-d np.arrays a*A*b if a.shape,A.shape,b.shape = ((1, n), (m, n), (m, 1))
        # caution: self._yf_ self._yb_ are (ngroups,nbins), x is (1,nbins) and alpha is (ngropus,1), hence
        #          alpha multiplies from the right 
        #          x functions multipy from the left
        
        denfactorleft = exp(-x/TauMu_mus())
        denfactorright = self._yfm_ + self._ybm_*alpha
        denominator = denfactorleft*denfactorright # f+b normalization count
        self._y_ = (self._yf_ - self._yb_*alpha - (self._bf_ - self._bb_*alpha)) / denominator 
        errexp = sqrt(self._yf_ + self._yb_*alpha**2) # equivalent f+b counts
        errexp[where(errexp==0)] = 1  #   set to 1 the minimum error for zero equivalent f+b counts
        self._e_ = errexp / denominator 
        
        # self._component_ contains [bndkeys,...,bndkeys], as many as the model components (e.g. 2 for mgbl)
        # bndkeys is [method, [keys_1,keys_2]] if there are 2 groups, keys_i is a list of keys for group i=1,2   
        # such that method(x,*par_i),  produce the additive function component for group i
        # and par_i[k] = eval(keys_i[k])   
        for j in range(self._n0truecomponents_,self._ntruecomponents_): # all components in model excluding "al", which must always be the first
            component = self._components_[j][0]
            keys = self._components_[j][1] # = [keys_1,keys_2,...]
            # keys = [[p0g0, p0g1,...],[p1g0, p1g1, ..],[p2g0, p2,g1,...]..]
            # print('add_multigroup mucomponents debug: key = {}'.format(keys))
            pars = [[key(p) for key in groups_key] for groups_key in keys]# NEW! spedup, evaluates p[1], p[2] etc.
            f += component(x,*pars)  # must contain x, 
                                                 # for plot x != self._x_
            # remember *p.comp means 'pass as many arguments as required by component, exausting the list p_comp'

            # print('add_multigroup mucomponents debug: pars = {}'.format(pars))
            # pars = [[eval(key) for key in groups_key] for groups_key in keys]
            # print('add_multigroup mucomponents debug: y:{},x:{},f:[]'.format(self._y_.shape,x.shape,f.shape))
            # print('add_multigroup mucomponents debug: pars = {}'.format(pars))
            # print('add_multigroup mucomponents debug: f.shape = {}, zeros.shape = {}'.format(
            #                                                         f.shape,zeros_like(x).shape))
        return f     

    def _add_single_calib_(self,x,*argv):
        '''
        use _add_single_ instead, after _load_data_calib_
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    alpha,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [alpha,val1,val2,val3,val4,val5, ...]

        _add_single_calib_ DISTRIBUTES THESE PARAMETER VALUES for plots::

              asymmetry fit with fitted alpha
              version for plotting calib fits as normal asymmetry fits
              order driven by model e.g. alml
        NO FFT mode, no check on self._include_components
        '''      
#        from numpy import where
#        from mujpy.aux.aux import TauMu_mus
#                
#        f = zeros_like(x)  # initialize a 1D array
#        p = argv 
#                
#        for j in range(self._n0truecomponents_,self._ntruecomponents_): 
#            method = self._components_[j][0]
#            keys = self._components_[j][1] 
#            pars = [key(p) for key in keys] 
#            f += method(x,*pars)  
#        return f   

    def _add_single_multigroup_calib_(self,x,*argv):
        '''
        use instead _add_multigroup_ after _load_data_multigroup_calib_ 
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    alpha,val1,val2,val3,val4,val5, ... at this iteration 
                    argv is a list of values [alpha,val1,val2,val3,val4,val5, ...]

        _add_single_calib_ DISTRIBUTES THESE PARAMETER VALUES for plots::

              asymmetry fit with fitted alpha
              version for plotting calib fits as normal asymmetry fits
              order driven by model e.g. alml
        NO FFT mode, no check on self._include_components
        '''
#        from numpy import where
#        from mujpy.aux.aux import TauMu_mus
#                
#        f = zeros((self._y_.shape[0],x.shape[0]))  # initialize a 1D array
#        p = argv 
#        # alpha = p[0]
#        # print('_add_single_multigroup_calib_ debug alpha = {}, p = {}'.format(alpha,p))
#                
#        for j in range(self._n0truecomponents_,self._ntruecomponents_): # all components in model excluding "al", which must always be the first
#            component = self._components_[j][0]
#            keys = self._components_[j][1] # = [keys_1,keys_2,...]
#            # keys = [[p0g0, p0g1,...],[p1g0, p1g1, ..],[p2g0, p2,g1,...]..]
#            # print('add_multigroup mucomponents debug: key = {}'.format(keys))
#            pars = [[key(p) for key in groups_key] for groups_key in keys]# NEW! spedup, evaluates p[1], p[2] etc.
#            f += component(x,*pars)  # must contain x, 
#                                                 # for plot x != self._x_
#            # remember *p.comp means 'pass as many arguments as required by component, exausting the list p_comp'

#            # print('add_multigroup mucomponents debug: pars = {}'.format(pars))
#            # pars = [[eval(key) for key in groups_key] for groups_key in keys]
#            # print('add_multigroup mucomponents debug: y:{},x:{},f:[]'.format(self._y_.shape,x.shape,f.shape))
#            # print('add_multigroup mucomponents debug: pars = {}'.format(pars))
#            # print('add_multigroup mucomponents debug: f.shape = {}, zeros.shape = {}'.format(
#            #                                                         f.shape,zeros_like(x).shape))
#        return f     


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
        fit component for alpha calibration 
        x [mus], dα
        x for compatibility, here it is dummy anyway
        '''
        # empty method  (could remove x from argument list ?)
        # print('al = {}'.format(α))
        return []
        al.func_code = make_func_code(["α"])                
           
    def bl(self,x,A,λ): 
        '''
        fit component for a Lorentzian decay, 
        x [mus], A, λ [mus-1]
        x need not be self.x (e.g. in plot)
        '''
        # λ = -87. if λ < -87. else λ
        return A*exp(-x*λ)
        bl.func_code = make_func_code(["A","λ"])

    def bg(self,x,A,σ): 
        '''
        fit component for a Gaussian decay, 
        x [mus], A, σ [mus-1] (positive parity)
        x need not be self.x (e.g. in plot)
        '''
        return A*exp(-0.5*(x*σ)**2)
        bg.func_code = make_func_code(["A","σ"])

    def bs(self,x,A,Λ,β): 
        '''
        fit component for a stretched decay, 
        x [mus], A, Λ [mus-1] (>0), β (>0)
        x need not be self.x (e.g. in plot)
        '''
        return A*exp(-(x*Λ)**β)
        bs.func_code = make_func_code(["A","Λ","β"])

    def da(self,x,dα):
        '''
        fit component for linearized alpha correction
        x [mus], dα
        x for compatibility, here it is dummy anyway
        '''
        # the returned value will not be used, correction in _add_
        # print('dα = {}'.format(dα))
        return [] # zeros(self.x.shape[0]) # alternatively, return [], remove x from argument list
        da.func_code = make_func_code(["dα"])


    def ml(self,x,A,B,φ,λ): 
        '''
        fit component for a precessing muon with Lorentzian decay, 
        x [mus], A, B [mT], φ [degrees], λ [mus-1]
        x need not be self.x (e.g. in plot)
        '''
#        import warnings
#        warnings.filterwarnings("error")
        # print('a={}, B={}, ph={}, lb={}'.format(asymmetry,field,phase,Lor_rate))
#        try:
        # λ = -87. if λ < -87. else λ
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)
#        except RuntimeWarning:
#            print('debug: λ = {}'.format(λ))
#            return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)      
        ml.func_code = make_func_code(["A","B","φ","λ"])

    def mg(self,x,A,B,φ,σ): 
        '''
        fit component for a precessing muon with Gaussian decay, 
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  (positive parity)
        x need not be self.x (e.g. in plot)
        '''
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-0.5*(x*σ)**2)
        mg.func_code = make_func_code(["A","B","φ","σ"])

    def md(self,x,A,B,φ,λ,σ): 
        '''
        fit component for a precessing muon with Gaussian and Lorentzian independent decaya, 
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  (positive parity)
        x need not be self.x (e.g. in plot)
        '''
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)*exp(-0.5*(x*σ)**2)
        mg.func_code = make_func_code(["A","B","φ","λ","σ"])

    def ms(self,x,A,B,φ,Λ,β): 
        '''
        fit component for a precessing muon with stretched decay, 
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        x need not be self.x (e.g. in plot)
        '''
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-(x*Λ)**β)
        ms.func_code = make_func_code(["A","B","φ","Λ","β"])

    def fm(self,x,A,B,λ):
        '''
        fit component for FmuF (powder average)
        x [mus], A, B [mT], λ [mus-1]
        x need not be self.x (e.g. in plot)
        '''
        return A/6.0*( 3.+cos(2*pi*self._gamma_Mu_MHzper_mT*B*sqrt(3.)*x)+
               (1.-1./sqrt(3.))*cos(pi*self._gamma_Mu_MHzper_mT*Bd*(3.-sqrt(3.))*x)+
               (1.+1./sqrt(3.))*cos(pi*self._gamma_Mu_MHzper_mT*Bd*(3.+sqrt(3.))*x) )*exp(-x*λ)
        fm.func_code = make_func_code(["A","B","λ"])

    def jl(self,x,A,B,φ,λ): 
        '''
        fit component for a Bessel j0 precessing muon with Lorentzian decay, 
        x [mus], A, B [mT], φ [degrees], λ [mus-1]
        x need not be self.x (e.g. in plot)
        '''
        return A*j0(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)
        jl.func_code = make_func_code(["A","B","φ","λ"])

    def jg(self,x,A,B,φ,σ): 
        '''
        fit component for a Bessel j0 precessing muon with Gaussian decay, 
        x [mus], A, B [mT], φ [degrees], σ [mus-1] (positive parity)
        x need not be self.x (e.g. in plot)
        '''
        return A*j0(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-0.5*(x*σ)**2)
        jg.func_code = make_func_code(["A","B","φ","σ"])

    def js(self,x,A,B,φ,Λ,β): 
        '''
        fit component for a Bessel j0 precessing muon with Gaussian decay, 
        x [mus], A, B [mT], φ [degrees], σ [mus-1]
        x need not be self.x (e.g. in plot)
        '''
        return A*j0(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-(x*Λ)**β)
        jg.func_code = make_func_code(["A","B","φ","Λ","β"])

    def _kg(self,t,w,Δ):
        '''
        auxiliary component for a static Gaussian Kubo Toyabe in longitudinal field, 
        t [mus], w [mus-1], Δ [mus-1], note that t can be different from self._x_
        w = 2*pi*gamma_mu*L_field
        '''
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
        t [mus], w [mus-1], Δ [mus-1], note that t can be different from self._x_
        w = 2*pi*gamma_mu*L_field
        '''
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
        Gaussian Kubo Toyabe in (fixed) longitudinal field, static or dynamic
        x [mus], A, BL [mT], Δ [mus-1] (positive parity), ν (MHz)
        x need not be self.x (e.g. in plot)
        '''
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
        Lorentzian Kubo Toyabe in (fixed) longitudinal field, static (dynamic makes no sense)
        x [mus], A, BL [mT], Γ [mus-1] 
        x need not be self.x (e.g. in plot)
        '''
        w = 2*pi*BL*self._gamma_Mu_MHzper_mT
        return A*self._kl(x,w,Γ)
        kl.func_code = make_func_code(["A","BL","Γ"])

    def kd(self,x,A,Δ,λ):
        '''
        static Gaussian Kubo Toyabe times an (independent) exponential decay
        x [mus], A, B [T], Δ [mus-1], ν (MHz)
        x need not be self.x (e.g. in plot)
        '''
        return A*self._kg(x,0,Δ)*exp(-x*λ)
        kd.func_code = make_func_code(["A","Δ","λ"])
        #kd.limits = [[None,None],[0.,None],[None,None]]
        #kd.error = [0.002,0.05,0.05]

    def ks(self,x,A,Δ,Λ,β):
        '''
        static Gaussian Kubo Toyabe times an (independent) stretched exponential decay
        x [mus], A, B [T], Δ [mus-1], ν (MHz)
        x need not be self.x (e.g. in plot)
        '''
        return A*self._kg(x,0,Δ)*exp(-(x*Λ)**β)
        kd.func_code = make_func_code(["A","Δ","Λ","β"])
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
        from numpy import finfo, where, array, abs
        Mepsi = finfo('d').max/10.
        num = abs(self._add_(self._x_,*argv) - self._y_)
        normsquaredev = (num/self._e_)**2
        divergence = normsquaredev>Mepsi
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
        if l is None:
            
            return sum(  ( (self._add_(self._x_,*argv) - self._y_[k,:]) /self._e_[k,:])**2 )
        else:
            return sum(  ( (self._add_(self._x_,*argv) - self._y_[k,l,:]) /self._e_[k,l,:])**2 )
            
    from iminuit import Minuit as _M
    _chisquare_.errordef = _M.LEAST_SQUARES

