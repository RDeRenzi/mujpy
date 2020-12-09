from numpy import cos, pi, exp, sqrt, real, nan_to_num, inf, ceil, linspace, zeros, empty, ones, hstack, fft, sum, zeros_like
from scipy.special import dawsn,erf, j0
from scipy.constants import physical_constants as C

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
        # self._alpha_ =  [] to be deleted
    # ---- end generic __init__

    def _add_(self,x,*argv):
        '''
        e.g. a blmg global model

        * argv will be a list parameter values (val1,val2.val3,val4,val5,val6, ...) at this iteration 

        _add_ DISTRIBUTES THESE PARAMETER VALUES::

              SINGLE OR SUITE OF RUNS WITH self._global_==False
              order driven by model e.g. blml
              see mugui int2_int, _add_ replicates the same loops
                loop over model components (bl & ml)
                    loop over pars (component parameters, two and four, respectively)
                        eval(pars) # pars is a string set by int2_int as keys

        use to plot as follows::
 
          j = -1
          yoffset = 0.05
          for y,e in zip(self.asymm,self.asyme)
              j += 1
              plt.errorbar(x,y+j*yoffset,yerr=e)
              plt.plotcompare two branches(x,mumodel()._add_global(x,val1,val2.val3,val4,val5,val6,run=j)+j*yoffset)

        UNIFIED WITH FFT, where, for each component f adds it::
   
            if self._fft_include_components[j] else 0. 
            if if self._fft_include_da else f

              WITH self._global_==True     
              order driven by model, e.g. blml, 
              global, global pars, local constants, local pars (see int2_int)
                 loop over global
                 loop over model components (bl & ml)  
                     loop over global pars
                         eval(pars) 
                 loop over runs
                     loop over local pars
                     loop over model components (bl & ml)
                         loop over local pars
                             eval(pars)  

        FINAL COMMMENT: eval implies both a python time overhead and a security break
                        is there a way to avoid it, implementing free parameter functions?    
        '''      

        # _load_data_ loads also: 
        # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
        # self._alpha_, self._da_index_ (index of dalpha or [])
        # self._nglobals_ = number of globals
        # self._ntruecomponents_ = number of components apart from dalpha 

        # With nglobal global parameters and nlocal local parameters 
        # there must be nruns*nlocal+nglobal minuit parameters
        # Additionally with nfixed local constants 
        # self.locals contains nruns*nfixed additional constants 
        nfixed = 0
        # print('x[0] = {}, x[-1] = {}, nbins = {}'.format(x[0],x[-1],x.shape[0]))

        if self._global_:
            if self._locals_:
                nfixed = self._locals_.shape[1]*self._locals_.shape[0] # number of local constants
            f = zeros((self._y_.shape[0],x.shape[0])) # initialize a 2D array
            for run in range(self._y_.shape[0]): # range(nruns)
                # argv is a list
                p = empty(nfixed+len(argv)) # internally, also locals are parameters
                if self._locals_:
                    p[:nfixed] = self._locals_[run][:] # locals for this run 
                # nothing, if run <1 nfixed = 0 and self._locals_ = None
                p[nfixed:nfixed+self._nglobals_] = argv[:self._nglobals_] # store global minuit parameter
                kp, ka = nfixed+self._nglobals_, self._nglobals_  # new version
                for j in range(self._ntruecomponents_): # all components in model excluding da
                    component = self._components_[j][0]
                    pars = self._components_[j][1]
                    ismin =  self._components_[j][2]
                    p_comp = []
                    for l in range(len(pars)):
                        p_comp.append(eval(pars[l]))
                        if ismin[l]:     # new version
                            p[kp] = argv[ka]
                            kp += 1
                            ka += 1
                    f[run,:] += component(x,*p_comp) if self._include_components[j] else 0.
                if self._da_index_:  # linearized correction 
                    dalpha = p[self._da_index_]
                    dada = dalpha/self._alpha_
                    f[run,:] = ((2.+dada)*f-dada)/((2.+dada)-dada*f)  if self._include_da else f
        else: # non global, single run fit, possibly one of a suite
            if self._locals_:
                nfixed = self._locals_.shape[0] # number of local constants per run
            f = zeros_like(x)  # initialize a 1D array
            p = empty(nfixed+len(argv)) # internally, also locals are parameters
            if self._locals_:
                p[:nfixed] = self._locals_ # locals for this run 
            p[nfixed:nfixed+self._nglobals_] = argv[:self._nglobals_] # store global minuit parameters
            p[(nfixed+self._nglobals_):] = argv[(self._nglobals_):]
            for j in range(self._ntruecomponents_): # all components in model excluding da
                component = self._components_[j][0]
                pars = self._components_[j][1] 
                p_comp = []
                for l in range(len(pars)):
                    p_comp.append(eval(pars[l]))
                # print('y:{},x:{},f:[]'.format(self._y_.shape,x.shape,f.shape))
                # print('f.shape = {}, zeros.shape = {}'.format(f.shape,zeros_like(x).shape))
                f += component(x,*p_comp) if self._include_components[j] else 0.
            if self._da_index_:  # linearized correction 
                dalpha = p[self._da_index_-1]
                dada = dalpha/self._alpha_
                f = ((2.+dada)*f-dada)/((2.+dada)-dada*f) if self._include_da else f
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

    def _load_data_(self,x,y,_components_,_alpha_,e=1,_nglobals_=None,_locals_=None):
        ''' ,
        use as        _load_data_(x,y,_components_,_alpha_,e,_nglobals_,_locals_,_ntruecomponents_)
        * x, y, e are numpy arrays
        * e = 1 yields unitary errors 
        * _components_ is a list [[method,[key,...,key]],...,[method,[key,...,key]]], where
        ** method is an instantiation of a component, e.g. self.ml 
        ** value = eval(key) produces the parameter value
        * _alpha_ is ditto
        * _nglobals_ = index of global parameter in iminuit parameter list (only for global fits)
        * _locals_ = values of local constants (e.g. B, T extrated form data file by musr2py), a 2D array for global fits

        Strategy to accommodate single runs, multi run suites and global fits:
        for global fit y is a 2D array and self._global_ = True
        for single and multi fit y is a 1D array 
        '''
        # self._components_ = [[method,[key,...,key]],...,[method,[key,...,key]]], and eval(key) produces the parmeter value
        # self._alpha_, self._da_index_ (index of dalpha or [])
        # self._nglobals_ = number of globals
        # self._locals_ = np.array nruns x nlocals local values 
        # self._ntruecomponents_ = number of components apart from dalpha 
        # local, y is a simple vector
        self._x_ = x
        self._y_ = y
        self._global_ = True if _nglobals_ is not None else False
        self._alpha_ = _alpha_
        self._components_ = []
        self._da_index_ = []
        self._ntruecomponents_ = 0
        for k, val in enumerate(_components_):
            if val[0]: # val[0] is directly the method for all components but dalpha
                self._ntruecomponents_ += 1
                self._components_.append(val) # store again [method, [key,...,key]] # also val[3] = isminuit, in new version
            else:  # when the method is da  (val[0] was set to [], i.e. False)
                self._da_index_ = 1+int(val[1][0][2:val[1][0].find(']')]) # position in minuit parameter list +1 to pass logical test
                # print('_da_index_ = {}'.format(self._da_index_-1))
        self._include_all_()
        if _nglobals_:
            self._nglobals_ = _nglobals_
        else:
            self._nglobals_ = 0 # to be used as index
        self._locals_ = _locals_
        if self._global_:
            # global
            try:
                # (should be y.shape[1]=x.shape[0])
                if y.shape[1]!=x.shape[0]: # not global, error!
                    print('mumodel._load_data: x, y have different lengths')
                    return -1 # error
                if isinstance(e,int):
                    self._e_ = ones((y.shape[0],x.shape[0]))
                else:
                    if e.shape[1]!=x.shape[1]:
                        raise ValueError('x, e have different lengths, {},{}'.format(x.shape,e.shape))           
                    else:
                        self._e_ = e
            except: # y.shape[1] does not exist
                print('mumodel._load_data: this is a single run!')
                return -1 # error
        else:
            # local
            if isinstance(e,int):
                self._e_ = ones((x.shape[0]))
            else:
                if e.shape[0]!=x.shape[0]:
                    raise ValueError('x, e have different lengths, {},{}'.format(x.shape,e.shape))          
                else:
                    self._e_ = e
        return 0 # no error

    def bl(self,x,asymmetry,Lor_rate): 
        '''
        fit component for a Lorentzian decay, 
        x [mus], asymmetry, Lor_rate [mus-1]
        '''
        return asymmetry*exp(-x*Lor_rate)

    def bg(self,x,asymmetry,Gau_rate): 
        '''
        fit component for a Gaussian decay, 
        x [mus], asymmetry, Gau_rate [mus-1]
        '''
        return asymmetry*exp(-0.5*(x*Gau_rate)**2)

    def bs(self,x,asymmetry,rate,beta): 
        '''
        fit component for a stretched decay, 
        x [mus], asymmetry, rate [mus-1], beta (>0)
        '''
        return asymmetry*exp(-(x*rate)**beta)

    def da(self,x,dalpha):
        '''
        fit component for linearized alpha correction
        x [mus], dalpha
        '''
        # the returned value will not be used, correction in _add_
        # print('dalpha = {}'.format(dalpha))
        return zeros(x.shape[0])

    def ml(self,x,asymmetry,field,phase,Lor_rate): 
        '''
        fit component for a precessing muon with Lorentzian decay, 
        x [mus], asymmetry, field [T], phase [degrees], Lor_rate [mus-1]
        '''
        # print('a={}, B={}, ph={}, lb={}'.format(asymmetry,field,phase,Lor_rate))
        return asymmetry*cos(2*pi*self._gamma_Mu_MHzper_mT*field*x+phase*self._radeg_)*exp(-x*Lor_rate)

    def mg(self,x,asymmetry,field,phase,Gau_rate): 
        '''
        fit component for a precessing muon with Gaussian decay, 
        x [mus], asymmetry, field [T], phase [degrees], Gau_rate [mus-1]
        '''
        return asymmetry*cos(2*pi*self._gamma_Mu_MHzper_mT*field*x+phase*self._radeg_)*exp(-0.5*(x*Gau_rate)**2)

    def ms(self,x,asymmetry,field,phase,rate,beta): 
        '''
        fit component for a precessing muon with stretched decay, 
        x [mus], asymmetry, field [T], phase [degrees], rate [mus-1], beta (>0)
        '''
        return asymmetry*cos(2*pi*self._gamma_Mu_MHzper_mT*field*x+phase*self._radeg_)*exp(-(x*rate)**beta)

    def fm(self,x,asymmetry,dipfield,Lor_rate):
        '''
        fit component for FmuF (powder average)
        '''
        return asymmetry/6.0*( 3.+cos(2*pi*self._gamma_Mu_MHzper_mT*dipfield*sqrt(3.)*x)+
               (1.-1./sqrt(3.))*cos(pi*self._gamma_Mu_MHzper_mT*dipfield*(3.-sqrt(3.))*x)+
               (1.+1./sqrt(3.))*cos(pi*self._gamma_Mu_MHzper_mT*dipfield*(3.+sqrt(3.))*x) )*exp(-x*Lor_rate)

    def jl(self,x,asymmetry,field,phase,Lor_rate): 
        '''
        fit component for a Bessel j0 precessing muon with Lorentzian decay, 
        x [mus], asymmetry, field [T], phase [degrees], Lor_rate [mus-1]
        '''
        return asymmetry*j0(2*pi*self._gamma_Mu_MHzper_mT*field*x+phase*self._radeg_)*exp(-x*Lor_rate)

    def jg(self,x,asymmetry,field,phase,Gau_rate): 
        '''
        fit component for a Bessel j0 precessing muon with Lorentzian decay, 
        x [mus], asymmetry, field [T], phase [degrees], Lor_rate [mus-1]
        '''
        return asymmetry*j0(2*pi*self._gamma_Mu_MHzper_mT*field*x+phase*self._radeg_)*exp(-0.5*(x*Gau_rate)**2)

    def _kg(self,x,w,Gau_delta):
        '''
        auxiliary component for a static Gaussian Kubo Toyabe in longitudinal field, 
        x [mus], w [mus-1], Gau_delta [mus-1]
        w = 2*pi*gamma_mu*L_field
        '''
        Dt = Gau_delta*x
        DDtt = Dt**2
        DD = Gau_delta**2
        sqr2 = sqrt(2)
        argf = w/(sqr2*Gau_delta)
        fdc = dawsn(argf)
        wx = w*x
        if (w!=0): # non-vanishing Longitudinal Field
            Aa = real(exp(-0.5*DDtt + 1j*wx)*dawsn(-argf - 1j*Dt/sqr2) )
            Aa[Aa == inf] = 0 # bi-empirical fix
            nan_to_num(Aa,copy=False) # empirical fix 
            A=sqr2*(Aa + fdc)
            f = 1. - 2.*DD/w**2*(1-exp(-.5*DDtt)*cos(wx)) + 2.*(Gau_delta/w)**3*A
        else:
            f = (1. + 2.*(1-DDtt)*exp(-.5*DDtt))/3.
        return f

    def _kgdyn(self,x,w,Gau_delta,jump_rate,*argv):
        ''' 
        auxiliary dynamization of Gaussian Kubo Toyabe 
        by G. Allodi 
        N: number of sampling points;
        dt: time interval per bin [i.e. time base is t = dt*(0:N-1)]
        w [mus-1], Gau_delta [mus-1], jump_rate [MHz] 
        (longitudinal field freq, dGaussian distribution, scattering frequency 
        % alphaN: [optional argument] weighting coefficient alpha times N. Default=10 
        '''
        alphaN = 10. if not argv else argv[0] # default is 10.
        dt = x[1]-x[0]
        N = x.shape[0] + int(ceil(x[0]/dt)) # for function to include t=0
        Npad = N * 2 # number of total time points, includes as many zeros
        t = dt*linspace(0.,Npad-1,Npad)
        expwei = exp(-(alphaN/(N*dt))*t)

        gg = self._kg(t,w,Gau_delta)*(t < dt*N)  #  padded_KT
        # gg = 1/3*(1 + 2*(1 - s^2*tt.^2).*exp(-(.5*s^2)*tt.^2)) % 

        ff = fft.fft(gg*expwei*exp(-jump_rate*t)) # fft(padded_KT*exp(-jump_rate*t))
        FF = exp(-jump_rate*dt)*ff/(1.-(1.-exp(-jump_rate*dt))*ff) # (1-jump_rate*dt*ff)  

        dkt = real(fft.ifft(FF))/expwei  # ifft
        dkt = dkt[0:N] # /dkt(1) 

        #if (nargout > 1),
        #   t = t[0:intN-1]
        return dkt
         
    def kg(self,x,asymmetry,L_field,Gau_delta,jump_rate):
        '''
        Gaussian Kubo Toyabe in longitudinal field, static or dynamic
        x [mus], asymmetry, L_field [T], Gau_delta [mus-1], jump_rate (MHz)
        '''
        N = x.shape[0]
        w = 2*pi*L_field*self._gamma_Mu_MHzper_mT
        if jump_rate==0: # static 
           f = self._kg(x,w,Gau_delta) # normalized to 1.
        else :            # dynamic
           # P=[w Gau_delta];
 
           f = self._kgdyn(x,w,Gau_delta,jump_rate)
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
        f = asymmetry*real(f[0:N])
        return f

    def _chisquare_(self,*argv,axis=None):
        '''
        Signature provided at Minuit invocation by 
        optional argument forced_parameters=parnames
        where parnames is a tuple of parameter names::

           e.g. parnames = ('asym','field','phase','rate') 

        Works also for global fits, 
        where sum (...,axis=None) yields the sum over all indices.

        Provides partial chisquares over individual runs if invoked as::

            self._chisquare_(*argv,axis=1)

        '''
        return sum(  ( (self._add_(self._x_,*argv) - self._y_) /self._e_)**2 ,axis=axis )

