from numpy import cos, pi, exp, sqrt, real, nan_to_num, inf, ceil, linspace, zeros, empty, ones, hstack, fft, sum, zeros_like
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
                     'kg':r'Gauss Kubo-Toyabe: static and dynamic, in zero or longitudinal field by G. Allodi [Phys Scr 89, 115201]'}
        self._axis_ = None # for self._chisquare_ when set 0,1 sums only along that axis
        # self._alpha_ =  [] to be deleted
    # ---- end generic __init__

    def _load_data_(self,x,y,_components_,_alpha_,e=1):
        '''
        input: 
            x, y, e are numpy arrays, y, e are 1d 
            e = 1 or missing yields unitary errors 
            _components_ is a list [[method,[key,...,key]],...,[method,[key,...,key]]], 
                produced by int2_method() from mujpy.aux.aux
                where method is an instantiation of a component, e.g. self.ml 
                and value = eval(key) produces the parameter value
            _alpha_ is the fixed value used to produce asymmetry y and needed only by method da
        the argument _components_ is stored in self._components_
        '''
        # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
        # self._alpha_, self._da_index_ (index of dalpha or [])
        # self._ntruecomponents_ = number of components apart from dalpha 
        self._x_ = x
        self._y_ = y        # self._global_ = True if _nglobals_ is not None else False
        self._alpha_ = _alpha_
        self._components_ = []
        self._da_index_ = []
        self._add_ = self._add_single_
        self._ntruecomponents_ = 0

        for k, val in enumerate(_components_):
            if val[0]: # val[0] is directly the method for all components but dalpha
                self._ntruecomponents_ += 1
                self._components_.append(val) # store again [method, [key,...,key]], ismin
            else:  # when the method is da  (val[0] was set to [], i.e. False)
                self._da_index_ = 1+int(val[1][0][2:val[1][0].find(']')]) # position in minuit parameter list +1 to pass logical test
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
        for k, val in enumerate(components):
            #if val[0]: # val[0] is directly the method for all components but al
            # however al is taken care of by range(1,self._ntruecomponents_) in _add_calib_single_
            self._ntruecomponents_ += 1 # total number of components including 'al'
        # print('_load_calib_single_data_ mucomponents debug: self._ntruecomponents_ = {}'.format(self._ntruecomponents_))

    def _load_data_calib_(self,x,y,components,_alpha_,e=1): 

        '''
        for normal asymmetry plot of a calib best fit (or guess)
        input: 
            x, y, e are numpy arrays, y, e are 1d 
            e = 1 or missing yields unitary errors 
            _components_ is a list [[method,[key,...,key]],...,[method,[key,...,key]]], 
                produced by int2_method() from mujpy.aux.aux
                where method is an instantiation of a component, e.g. self.ml 
                and value = eval(key) produces the parameter value
            _alpha_ is the fixed value used to produce asymmetry y and needed only by method da
        the argument _components_ is stored in self._components_
        '''
        # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
        # self._alpha_, self._da_index_ (index of dalpha or [])
        # self._ntruecomponents_ = number of components apart from dalpha 
        self._x_ = x
        self._y_ = y        # self._global_ = True if _nglobals_ is not None else False
        self._alpha_ = _alpha_
        self._components_ = components
        self._add_ = self._add_single_calib_
        self._ntruecomponents_ = 0

        for k, val in enumerate(components):
                self._ntruecomponents_ += 1
        self._include_all_() # to recover from possible fft mode
        try:
            if isinstance(e,int):
                self._e_ = ones((y.shape))
            else:
                if len(y.shape)>1:
                    if e.shape!=y.shape or x.shape[0]!=eshape[-1]:
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
        # self._alpha_, self._da_index_ (index of dalpha or [])
        # self._ntruecomponents_ = number of components apart from dalpha 
        self._x_ = x
        self._y_ = y        # self._global_ = True if _nglobals_ is not None else False
        self._components_ = []
        self._add_ = self._add_multigroup_
        self._ntruecomponents_ = 0

        for k, val in enumerate(components):
            self._ntruecomponents_ += 1
            self._components_.append(val) # store again [method, [key,...,key]], ismin
        try:
            if isinstance(e,int):
                self._e_ = ones((y.shape))
            else:
                if len(y.shape)>1:
                    # print('_load_data_multigroup_ mucomponents debug: x,y,e not e=1')
                    if e.shape!=y.shape or x.shape[0]!=y.shape[-1]:
                        print('_load_data_multigroup_ mucomponents debug: x,y,e different shape[0]>1')
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
        except ValueError as e:
            return False, e       
        return True, '' # no error

    def _add_single_(self,x,*argv):   ## must be (self,x,*argv): it can also be invoked by plot over x != self._x_
        '''
         input: 
            x       time array
            *argv   passed as a variable number of parameter values 
                    val1,val2,val3,val4,val5,val6, ... at this iteration 
                    argv is a list of values [val1,val2,val3,val4,val5,val6, ...]

        _add_ DISTRIBUTES THESE PARAMETER VALUES::

              asymmetry fit with fixed alpha
              order driven by model e.g. blml

        UNIFIED WITH FFT, where, for each component f adds it::
   
            if self._fft_include_components[j] else 0. 
            if if self._fft_include_da else f

        FINAL COMMMENT: eval implies both a python time overhead and a security breach
                        is there a way to avoid it, implementing free parameter functions?    
        '''      

        f = zeros_like(x)  # initialize a 1D array
        p = argv 
        # print('add_single mucomponents debug: p = {}'.format(p))
        for j in range(self._ntruecomponents_): # all components in model excluding da
            component = self._components_[j][0]
            keys = self._components_[j][1] 
            # print('add_single mucomponents debug: keys = {}'.format(keys))
            pars = []
            for l in range(len(keys)):
                pars.append(eval(keys[l])) # typically evaluates p[1], p[2] etc.
            # print('y:{},x:{},f:[]'.format(self._y_.shape,x.shape,f.shape))
            # print('pars = {}'.format(pars))
            # print('f.shape = {}, zeros.shape = {}'.format(f.shape,zeros_like(x).shape))
            f += component(x,*pars) if self._include_components[j] else 0. # must contain x, 
                                                                           # for plot x != self._x_
            # remember *p.comp means 'pass as many arguments as required by component, exausting the list p_comp'
        if self._da_index_:  # linearized correction 
            dalpha = p[self._da_index_-1]
            dada = dalpha/self._alpha_
            f = ((2.+dada)*f-dada)/((2.+dada)-dada*f) if self._include_da else f
        return f     

    def _add_single_calib_(self,x,*argv):
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
                
        for j in range(1,self._ntruecomponents_): # all components in model, excluding alpha
            method = self._components_[j][0]
            keys = self._components_[j][1] 
            pars = []
            for l in range(len(keys)):
                pars.append(eval(keys[l])) # typically evaluates p[1], p[2] etc.
            # print('_add_calib_single_ debug y:{},x:{},f:[]'.format(self._y_.shape,x.shape,f.shape))
            # print('_add_single_calib_ debug pars = {} for component {}/{}'.format(pars,j+1,self._ntruecomponents_))
            # print('_add_calib_single_ debug f.shape = {}, zeros.shape = {}'.format(f.shape,zeros_like(x).shape))
            f += method(x,*pars)  # must contain x, for plot x != self._x_
            # remember *p.comp means 'pass as many arguments as required by component, exausting the list p_comp'

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
                
        for j in range(1,self._ntruecomponents_): # all components in model, excluding alpha
            method = self._components_[j][0]
            keys = self._components_[j][1] 
            pars = []
            for l in range(len(keys)):
                pars.append(eval(keys[l])) # typically evaluates p[1], p[2] etc.
            # print('_add_calib_single_ debug y:{},x:{},f:[]'.format(self._y_.shape,x.shape,f.shape))
            # print('_add_calib_single_ debug pars = {}'.format(pars))
            # print('_add_calib_single_ debug f.shape = {}, zeros.shape = {}'.format(f.shape,zeros_like(x).shape))
            f += method(x,*pars)  # must contain x, for plot x != self._x_
            # remember *p.comp means 'pass as many arguments as required by component, exausting the list p_comp'

        return f     

    def _add_multigroup_(self,x,*argv):   ## must be (self,x,*argv): it can also be invoked by plot over x != self._x_
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

        f = zeros((self._y_.shape[0],x.shape[0]))  # initialize a 2D array
        
        p = argv 
        
        # self._component_ contains [bndkeys,...,bndkeys], as many as the model components (e.g. 2 for mgbl)
        # bndkeys is [method, [keys_1,keys_2]] if there are 2 groups, keys_i is a list of keys for group i=1,2   
        # such that method(x,*par_i),  produce the additive function component for group i
        # and par_i[k] = eval(keys_i[k])   
        for j in range(self._ntruecomponents_): # all components in model excluding da
            component = self._components_[j][0]
            keys = self._components_[j][1] # = [keys_1,keys_2,...]
            # keys = [[p0g0, p0g1,...],[p1g0, p1g1, ..],[p2g0, p2,g1,...]..]
            # print('add_multigroup mucomponents debug: key = {}'.format(keys))
            pars = []
            for groups_key in keys:
                pargroup =[]
                for key in groups_key:
                    pargroup.append(eval(key))
                pars.append(pargroup)
            # print('add_multigroup mucomponents debug: pars = {}'.format(pars))
            # pars = [[eval(key) for key in groups_key] for groups_key in keys]
            # print('add_multigroup mucomponents debug: y:{},x:{},f:[]'.format(self._y_.shape,x.shape,f.shape))
            # print('add_multigroup mucomponents debug: pars = {}'.format(pars))
            # print('add_multigroup mucomponents debug: f.shape = {}, zeros.shape = {}'.format(f.shape,zeros_like(x).shape))
            f += component(x,*pars)  # must contain x, 
                                                 # for plot x != self._x_
            # remember *p.comp means 'pass as many arguments as required by component, exausting the list p_comp'
        return f     

    def _fft_init(self,include_components,include_da=True):
        '''
        saves the string of component flags used to generate a partial residue for FFT
        '''
        self._include_components = include_components # True means include, False, do not include 
                                              # in function f, for partial residues = asymm - f
        self._include_da = include_da # True means "da is a component" and "include it" in function f,
                                              # for partial residues = asymm - f, False, all other combinations

    def _include_all_(self):
        '''
        reset to normal fit mode (initially of after fft)
        '''
        self._include_components = [True]*self._ntruecomponents_
        self._include_da = True

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
        return A*exp(x*λ)
        bl.func_code = make_func_code(["A","λ"])

    def bg(self,x,A,σ): 
        '''
        fit component for a Gaussian decay, 
        x [mus], A, σ [mus-1]
        x need not be self.x (e.g. in plot)
        '''
        return A*exp(-0.5*(x*σ)**2)
        bg.func_code = make_func_code(["A","σ"])

    def bs(self,x,A,λ,β): 
        '''
        fit component for a stretched decay, 
        x [mus], A, λ [mus-1], β (>0)
        x need not be self.x (e.g. in plot)
        '''
        return A*exp(-(x*λ)**β)
        bs.func_code = make_func_code(["A","λ","β"])

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
        # print('a={}, B={}, ph={}, lb={}'.format(asymmetry,field,phase,Lor_rate))
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)
        ml.func_code = make_func_code(["A","B","φ","λ"])

    def mg(self,x,A,B,φ,σ): 
        '''
        fit component for a precessing muon with Gaussian decay, 
        x [mus], A, B [mT], φ [degrees], σ [mus-1]
        x need not be self.x (e.g. in plot)
        '''
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-0.5*(x*σ)**2)
        mg.func_code = make_func_code(["A","B","φ","σ"])

    def ms(self,x,A,B,φ,λ,β): 
        '''
        fit component for a precessing muon with stretched decay, 
        x [mus], A, B [mT], φ [degrees], λ [mus-1], β (>0)
        x need not be self.x (e.g. in plot)
        '''
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-(x*λ)**beta)
        ms.func_code = make_func_code(["A","B","φ","λ","β"])

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
        x [mus], A, B [mT], φ [degrees], σ [mus-1]
        x need not be self.x (e.g. in plot)
        '''
        return A*j0(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-0.5*(x*σ)**2)
        jg.func_code = make_func_code(["A","B","φ","σ"])

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

    def _kgdyn(self,x,w,Δ,ν,*argv):
        ''' 
        auxiliary dynamization of Gaussian Kubo Toyabe 
        by G. Allodi 
        N: number of sampling points;
        dt: time interval per bin [i.e. time base is t = dt*(0:N-1)]
        w [mus-1], Δ [mus-1], ν [MHz] 
        (longitudinal field freq, dGaussian distribution, scattering frequency 
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
         
    def kg(self,x,A,B,Δ,ν):
        '''
        Gaussian Kubo Toyabe in longitudinal field, static or dynamic
        x [mus], A, B [T], Δ [mus-1], ν (MHz)
        x need not be self.x (e.g. in plot)
        '''
        N = x.shape[0]
        w = 2*pi*B*self._gamma_Mu_MHzper_mT
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
        kg.func_code = make_func_code(["A","B","Δ","ν"])

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
        return sum(  ( (self._add_(self._x_,*argv) - self._y_) /self._e_)**2 ,axis=self._axis_ )
        
    from iminuit import Minuit as _M
    _chisquare_.errordef = _M.LEAST_SQUARES

