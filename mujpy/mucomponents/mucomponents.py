from numpy import cos, sin, pi, exp, sqrt, log, real, nan_to_num, inf, ceil, linspace, zeros, empty, ones, hstack, fft, sum, zeros_like, abs, array, where, arctan, pi, any, set_printoptions
from scipy.special import dawsn,erf, j0, j1
from scipy.constants import physical_constants as C
from iminuit.util import make_func_code

class mumodel(object):
    """
    Defines the components of the fitting model. Provides a chi_square function for Minuit.

    #################################################################################################
    # mumodel knows it all, for homogeneous fit and plot
    #                           fit types: A1, A20, A21, B1, B20, B21, C1, C2 (see mufit)
    #                           calib fits: alpha is a parameter, asymm, asyme must be recalculated
    #----------- Methods --------------------------------------------------------------------------
    # _load_(suite, methods_keys)       unique entry point for any mufit specific model
    #                                   selects a suite data slice 
    #                                   overloads _add_ as _add_plain or _add_calib_ 
    # _rebin_(tup)                      mufit data frontend, rebins the slice
    #                                   calls _set_alpha_, _asymmetry_ if _calib else suite.asymmetry_slice in _load_
    # _rebin_plot_(tup,lastfits)        mufitplot data frontend, rebins the entire suite
    #                                   allows for suite > slice (sequential fits A20, B1, B20, B21)
    #                                   calls _set_alpha_plot, _asymmetry_plot_ if _calib else suite.slice_asymmetry(-1,-1)
    #---------- General structure ------------------------------------------------------------------
    # _fit_types                        overload mufit.fit_types, called in _alpha_plot_
    # _add_                             overloaded method to distributes Minuit parameters to fit components and
    #                                                     to add components into model function
    #                                   called by mufit dofit_ which invokes tools methods:
    #       int2min_...                       to pass val err fix lim names pospar as plain lists
    #           int2min for A1, A20, B1, B20 fits
    #           int2min_global for A21, B21, C1, C2 fits
    #       int2_..._method_key             to prepare [bndmthd [...[k,..]..]]
    #                                       each bndmthd method taking k = key_as_lambda functions for each marameter
    #                                       such that       pars = [eval('k(p)'), ...]     
    #                                                       bndmthd(x,*pars)  is the component value
    # _add_plot_                        overloaded method to distributes lastfits Minuit parameters to fit components and
    #                                                     to add components into model function
    #                                                     also for suite > fit slices (sequential fits A20, B1, B20, B21)
    # ------------------------------------------------------------------------------------------------------
    # Data, alphas and model functions share a (squeezed) numpy array structure: array([[[time bins], ...  ], ... ])
    #                                           1d, 2d or 3d                                   bins   groups  runs
    #                                           e.g. fit A1 has a 1d suite=slice, fit A20 2d suite, 1d slice, fit C2 3d suite=slice
    #       parameters replicate the structure as lists, hence keys in methods_keys do as well
    #       tools cstack or ccstack vectorize methods
    #       _set_alpha and _alpha_plot_ vectorize alphas
    #       _add_ and _add_plot: vectorizes 
    ####################################################################################################
    # methods not starting with _ are fit components, 
    # tools method _available_components_ automagically creates a dictionary ready for the dashboard.json but for the key 'value'
    #       mufit collects several methods in a model identified as 
    #                       ''.join(method.__getattribute__(name)) for method in model
    #
    # _chisquare_            Minuit cost function at self._x_ through self._add_, self_y_, self._e_
    # _chisquare_calib_      Minuit cost function at self._t_ through self._add_, self_y_, self._e_
    #                                                         recomputed by self._add_ -> calib
    ####################################################################################################
    """

    def __init__(self):
        """ 
        Defines few constants and _help_ dictionary
        """

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
        # ---- end generic __init__

    def _load_(self,suite,components,krun,kgroup,fit_types):
        """
        mufit entry, initialized a mumodel instance for best fit and for plot

        works for A1, A20, A21, B1, B20, B21, C1, C2 and calib
        input: 
            suite,          the fit suite instance (for mufitplot)
            components,     aka methods_keys 
                            a structure of methods and their list of lists ... of key
                               [[method,[key...]],...,[method,[key...]]], 
                            produced mujpy.tools.tools int2_[global_]method()
                            such that key(p) is a parameter for method component
            krun, kgroup    determines fit slice <= suite (-1 for the suite span)
            fit_types       A1...C2 = fit_types() passed to _alpha_plot_
            detects calib fits, 'alxxyy', since the 'al' method is None
        """

        from mujpy.tools.tools import rebin
        from mujpy.tools.tools import list_depth

        # no wrappers, mufitplot calls _asymmetry_plot_
        # self._asymmetry_ = self._asymmetry_fit_ #
        
        self._fit_types = fit_types # for _alpha_plot_
        #A1,A20,A21,B1,B20,B1,C1,C2 = fit_types()
        self._calib = not components[0][0] # the first component method, None for 'al' fits
        # preload data slice
        self._t_ = suite.time
        nruns, ngroups = len(suite.runs), len(suite.groups) # suite dimensions

        # selects fit slice dimensions 
        if self._calib:
            self._yf_, self._yb_, self._eyf_, self._eyb_ = suite.slice_for_back_counts(krun,kgroup)
            #print('mumodel _load_ krun, kgroup {}, {}'.format(krun,kgroup))
            #self._add_ = self._add_calib_ # for _add_plot_ [and for _chisquare_calib_]
        else: 
            self._a_,self._ea_  = suite.asymmetry_slice(krun,kgroup)
            #self._add_ = self._add_plain_ # for _add_plot_ [and for _chisquare_]
        #self._the_run_ = suite._the_runs_[krun][0] # useless, in suites it ends as the last loaded  
        self.suite = suite
        self._nruns_, self._ngroups_ = nruns, ngroups
        self._ntruecomponents_ = len(components) # remove? used in fft
        self._n0component = 1 if self._calib else 0

        # this blockdefined self._lambdakeys_ and self._keys for use in self._add_
        key_depth = list_depth(components) # 3,4,5, respectively, for single, multirg, or multir_multig
        # print('mumodel _load_ key_depth {}'.format(key_depth))
        # key_depth: 3 (A1, A20, B1, B20),   
        #            4 (A21, B21, C1)
        # _add_ loops for mthd, keys in zip(self._methods_,self._keys_): 
        #                 pars = keys(p,keys)
        #                 f += methd(t,*pars)
        #print('mumodel _load_ list_depth(components) {}'.format(key_depth))
        self._methods_ = [mthd_keys[0] for mthd_keys in components] # extracts methods 
        self._keys_ = [mthd_keys[1] for mthd_keys in components] # extracts corresponding keys
        self._lambdakeys_ = []
        for component in components:
            if key_depth==3: # A1 A20 B1 B20 
                lkeys = lambda p,kys: [key(p) for key in kys]  
            elif key_depth==4: # B21 C1 
                lkeys = lambda p,kys: [[key(p) for key in rgkey] for rgkey in kys] 
            else: # C2
                lkeys = lambda p,kys: [[[key(p) for key in gkey] for gkey in rkey] for rkey in kys] # key_depth==4
            self._lambdakeys_.append(lkeys)
        #p = [0.556, 0.267, 0.249, 10.02, 35.0, 0.05, 0.8]
        #pars = [key(p) for key in self._keys[0][0]]
        #print('mumodel _load_ pars {}'.format(pars))
        #print('mumodel _load_ f = {}'.format(self._methods_[0](self._t_,*pars)))
            # so that in add_plain_ 
            # for component,lambdapars,keys in (zip(self._methods_[n0:],self._lambdakeys_[n0:],self._keys_[n0:])):
            #     pars = lambdapars(p,keys) 
            #     f += component(x,*pars)
        
    def _set_alpha_(self,argv):
        """
        set_alpha for mufit calib fits, for the suite slice
        
        selected by _load_ into _yf_, _yb_, _eyf_, _eyb_ 
        """
        
        from mujpy.tools.tools import set_alpha
        self._alpha = set_alpha(argv,self._nruns_,self._ngroups_)

    def _alpha_plot_(self,lastfits):
        """
        returns alpha for mufitplot calib, for the entire suite
        
        lastfits is a list of lists of values (results of sequential fits)
        if suite == fit slice it also sets self._alpha
        """
        from mujpy.tools.tools import set_alpha
        from numpy import array
        
        A1,A20,A21,B1,B20,B21,C1,C2 = self._fit_types()

        if A1 or A21 or C1 or C2: # suite == slice
            fit = lastfits[0] # if A1 else lastfits
            self._set_alpha_(fit)
            if C2:
                nr,ng = self._alpha.shape[0:2]
                self._alpha = self._alpha.reshape((nr*ng,-1))
            return self._alpha
        elif A20 or B1 or B21: # suite > slice
            if B21: # B21
                return array([[lfit[k]] for lfit in lastfits for k in range(self._ngroups_) ]) # one lfit per run, ngroup alphas
            else: # A20 or B1, lastfits are per groups or per runs, respectively
                return array([[lfit[0]] for lfit in lastfits]) # A20, B1: each lfit ha only one alpha, values[0] 
        else: # B20
            return array([[lastfits[krun*self._ngroups_+kgroup][0]] for kgroup in range(self._ngroups_) for krun in range(self._nruns_)]) 
            # B20 lfit lists ngroup fits run 0, ... ngroup fits run nrun-1, each one alpha
        # print('********* entered mucomponents _reload_ self._x_.shape,self._y_.shape {},{}'.format(self._x_.shape,self._y_.shape))

    def _rebin_(self,tup,*argv):
        """
        data frontend for mufit.dofit_ defines fit range from the unbinned _load_ slice
        
        (self. is implied in front of _ variables)
        input:
            tup = (start, stop, pack)
            *argv      = Minuit parameters, only for calib _set_alpha_
            if self._calib
                self._alpha = set_alpha(argv,self._nruns_,self._ngroups_)
                a, ea = _asymmetry_() unbinned, uses
                _yf_, _yb_, _eyf_, _eyb_ unbinned _load_ slice 
            else
                _a_, _ea_  unbinned _load_ slice
        output:
            rebinned _x_, _y_, _e_   _load_ slice
        """
        
        from mujpy.tools.tools import rebin, set_alpha

        self._include_all_() # to rcover from possible fft mode
        self._start_, self._stop_, self._pack_ = tup
        self._slice_runs_groups_()
        if self._calib:
            self._x_, self._yf_, self._eyf_ = rebin (self._t_,self._yf_,[self._start_,self._stop_],self._pack_,e=self._eyf_)
            self._x_, self._yb_, self._eyb_ = rebin (self._t_,self._yb_,[self._start_,self._stop_],self._pack_,e=self._eyb_)
            self._alpha = set_alpha(argv,self._slice_nruns_,self._slice_ngroups_)
            self._y_,self._e_ = self._asymmetry_()
            #print('mufit _rebin_ self._y_.shape {}'.format(self._y_.shape))
            self._radd_plain_(self._x_,*argv)
        else:
            self._x_,self._y_,self._e_ = rebin(self._t_,self._a_,[self._start_,self._stop_],self._pack_,e=self._ea_)

    def _slice_runs_groups_(self):
        """
        saves mufit self._slice_nruns_ self._slice_ngroups_

        the mumodel slice equivalents of suite self._nruns_ self._ngroups_ 
        """
        A1,A20,A21,B1,B20,B21,C1,C2 = self._fit_types()
        if A1 or A20 or B1 or B20:
            self._slice_nruns_, self._slice_ngroups_ = 1, 1 # single run single group fits
        elif A21 or B21:
            self._slice_nruns_, self._slice_ngroups_ = 1, self._ngroups_ # single run multi group fits
        elif C1: 
            self._slice_nruns_, self._slice_ngroups_ = self._nruns_, 1 # single run single group fits
        else: # C2
            self._slice_nruns_, self._slice_ngroups_ = self._nruns_, self._ngroups_ # single run single group fits
    

    def _rebin_plot_(self,tup,lastfits):
        """
        data frontend for muplotfit, used as follows
        
        input:
            tup = start, stop, pack
            lastfits from Minuit (take also care of guess=True)
        calls
            _asymmetry_plot_ which calls _alpha_plot_, if _calib
            suite.asymmetry_slice(-1,-1)  otherwise       
        directly loads self._x_, self._y_, self._e_  for the suite 
        """
        from mujpy.tools.tools import rebin, rshp

        start, stop, pack = tup
        if self._calib:
            yf, yb, eyf, eyb =  self.suite.slice_for_back_counts(-1,-1)
            #x, self._yf_, self._eyf_ = rebin (self._t_,yf,[self._start_,self._stop_],self._pack_,e=eyf)
            x, yf, eyf = rebin (self._t_,yf,[start,stop],pack,e=eyf)
            #print('mumodel _rebin_plot_ self._yf_.shape {}'.format(self._yf_.shape))
            self._yf_, self._eyf_ = rshp(yf),rshp(eyf)
            #x, self._yb_, self._eyb_ = rebin (self._t_,yb,[self._start_,self._stop_],self._pack_,e=eyb)
            x, yb, eyb = rebin (self._t_,yb,[start,stop],pack,e=eyb)
            self._yb_, self._eyb_ = rshp(yb),rshp(eyb)
            y,e = self._asymmetry_plot_(lastfits) # entire suite, uses best fit alphas
        else:
            a, ea = self.suite.asymmetry_slice(-1,-1) # entire suite, uses suite.groupings alphas
            x, y ,e = rebin(self._t_,a,[start,stop],pack,e=ea)
        self._x_,self._y_,self._e_ = rshp(x), rshp(y), rshp(e) 
        #print('mumodel _rebin_plot_ self._y_.shape {}'.format(self._y_.shape))
        #print('mumodel _rebin_plot_ self._x_.shape {}, self._y_.shape {}'.format(self._x_.shape,self._y_.shape))

    def _asymmetry_(self):
        """
        mufit.dofit_ calib asymmetry, called by model._rebin_calib_
        
            self._alpha is preset by _set_alpha_ in _add_calib_
        output unbinned asymm, asyme, slice according to _load_
            from unbinned slice self._yf_ self._yb_ self._eyf_ self._eyb_
        """

        from numpy import where

        # print('mumodel _asymmetry_ self._alpha {}, self._alpha.shape {}'.format(self._alpha,self._alpha.shape))
        denominator = (self._yf_ + self._alpha*self._yb_)
        i = where(denominator==0)
        denominator[i] = 1 
        nominator = (self._yf_ - self._alpha*self._yb_)
        y = nominator/denominator 
        e = 2*self._alpha/denominator**2*sqrt((self._yb_*self._eyf_)**2 + (self._yf_*self._eyb_)**2) 
        e[where(e==0)] = 1  
        return y,e

    def _asymmetry_plot_(self,lastfits):
        """
        mufitplot.chi calib _alpha_plot_ & asymmetry, called by model._rebin_plot_ 

        input
            lastfits     passed to _alpha_plot_ to produce alpha for the suite
        output
            asymm
            asyme        unbinned, for the entire suite
            from rebinned yf yb eyf eyb by _rebin_plot_
        """

        from numpy import where

        #print('_asymmetry_plot_ lastfits {}'.format(lastfits))
        alpha = self._alpha_plot_(lastfits)
        # print('_asymmetry_plot_ len(alpha) {} alpha {}'.format(len(alpha),alpha))
        denominator = (self._yf_ + alpha*self._yb_) # already rebinned
        i = where(denominator==0)
        denominator[i] = 1 
        nominator = (self._yf_ - alpha*self._yb_)
        y = nominator/denominator 
        e = 2*alpha/denominator**2*sqrt((self._yb_*self._eyf_)**2 + (self._yf_*self._eyb_)**2) 
        e[where(e==0)] = 1  
        return y,e

#####################################################################################
# add methods:                                                                      #
# 'time' array is x, not self._x_, plot invokes _add_... with different t vectors   #
#  _add_plain_            : any slice dimension,                                    #
#  _add_calib_            : any slice dimension,                                    #
#                           recalculates asymmetry with minuit alpha parameters     #
#                           first _load_calib_, Minuit calls _chisquare_calib_      #
# uadd_plot_              : entire suite for mufitplot                              #
#  _add_fft_              : single run, single group for partial residues           #
#####################################################################################

    def _radd_plain_(self,x,*argv):
        """
        mock _add_plain_ (see) for debug print purposes

        input
            x       time array
            argv   Minuit parameters
        """
        n0 = self._n0component
        f = zeros((self._slice_nruns_,self._slice_ngroups_,x.shape[0]))
        #print('mumodel _add_plain_ f.shape= {}'.format(f.shape))
        p = argv
        k = 0
        for component,lambdapars,keys in zip(self._methods_[n0:],self._lambdakeys_[n0:],self._keys_[n0:]):
            pars = lambdapars(p,keys)
            #print('mumodel _radd_plain_ k {}'.format(k))
            #for pa in pars:
            #    for par in pa:
            #        print('mumodel _radd_plain_ par {:.4f},{:.4f},{:.4f},{:.4f}'.format(par[0],par[1],par[2],par[3]))
            f += component(x,*pars)
            k += 1
        return f.squeeze()

    def _add_plain_(self,x,*argv):
        """
        overloads _add_ in cost function _chisquare_

        input:
            x       time array
            *argv   passed as a variable number of parameter values 
                    val1,val2,val3,val4,val5,val6, ... at this iteration 
                    argv is a list of values [val1,val2,val3,val4,val5,val6, ...]
        unique method for all non calib fits slice dimensions
              asymmetry fit with fixed alpha
        """
        
        set_printoptions(precision=2)
        n0 = self._n0component
        f = zeros((self._slice_nruns_,self._slice_ngroups_,x.shape[0]))
        #print('mumodel _add_plain_ f.shape= {}'.format(f.shape))
        p = argv
        for component,lambdapars,keys in (zip(self._methods_[n0:],self._lambdakeys_[n0:],self._keys_[n0:])):
            pars = lambdapars(p,keys)
            #print('mumodel _add_plain_ pars = {}'.format(pars))
            f += component(x,*pars)
        return f.squeeze()

    def _add_calib_(self,x,*argv):

        """
        calls _rebin_calib_ and wraps to _add_plain_

        overloads _add_ in cost function _chisquare_calib_ 
        """

        from mujpy.tools.tools import rebin, set_alpha
        #print('mumodel _add_calib_ argv = {}'.format(argv))
        self._alpha = set_alpha(argv,self._slice_nruns_,self._slice_ngroups_)
        #if not self._debug:
        #   print('mumodel _add_calib_ self._alpha = {}'.format(self._alpha))
        #   self._debug = 1
        self._y_,self._e_ = self._asymmetry_()
        #self._x_,self._y_,self._e_ = rebin(self._t_,y,[self._start_,self._stop_],self._pack_,e=e)
        return self._add_plain_(x,*argv) 

    def _add_plot_(self,x,lastfits):
        """
        produces suite f in mufitplot.chi for sequential fits A20, B1, B20, B21 (other fits call _add_)

        input:
            x is rebinned time vector
            lastfits is mufit.lastfits, a list of mufit.lastfit instances
                pars are mufit.lastfit.values (takes also care of guess)
            self._add_(x,*argv) returns the function used in the fit
                                            A20, B1, B20 a 1d array
                                            B21 a 2d array
            _add_plot_ vstacks it according to the full suite
                                            A20             _ngroups_
                                            B1, B21 _nruns_
                                            B20     _nruns_,_ngroups_
        """

        from numpy import array # , set_printoptions
        from mujpy.tools.tools import list_depth,rshp
        # distinguishes A20, B1, B21, B20 
        # by _nruns_, _ngroups_ and list_depth(_add_(x,*argv))
        mufit_f = self._add_plain_(x,*lastfits[0]) # this lastfit surely exists, this is the function used in mufit
        if isinstance(mufit_f, list): # never a list!
            f_depth = list_depth(mufit_f) # all fits have same f_depthi
            print('mumodel _add_plot_ f_depth = {}'.format(f_depth))
        else:
            f_depth = len((rshp(mufit_f)).shape) # this is always 2 at most
            #set_printoptions(precision = 2)
            #print('mumodel _add_plain_plot_ f {}'.format(mufit_f))
        if f_depth == 1: # A20 and B20 
            val = [values for values in lastfits] # [0,1]
            f = array([[self._add_plain_(x,*par)] for par in val]) # self._add_plain_ is a 1d np.array, f and self.model._y_ are 2d
            #print('mumodel _add_plain_plot_ val = {}'.format(val))
            #print('mumodel _add_plain_plot_ f = {}'.format(f))
        elif f_depth == 2 and self._ngroups_ == 1: # B1
            f = array([self._add_plain_(x,*par) for par in lastfits])
        elif f_depth == 2 and self._ngroups_ > 1: # B21 and C1
            if self._nruns_ == self._slice_nruns_: # C2
                f = rshp(mufit_f)
            else: # B21
                f = array([rshp(self._add_plain_(x,*par))[k,:] for par in lastfits for k in range(self._ngroups_)])
        # self._add_ = self._add_plain_ # self._add_calib_ if self._calib else
        # self._add_ = self._add_plain_ # self._add_calib_ if self._calib else 
        #else: # B20 and C2
        #    if len(lastfits) > 1: # B20
        #        values = [[lastfits[krun*self._ngroups_+kgroup] for kgroup in range(self._groups_)] for krun in range(self._nruns_)]
        #        f = array([[[self._add_plain_(x,*par)] for val in values] for par in val]) # run 0 
        #    else: # C2
        #        nr,ng,nb = f.shape
        #        f = mufit_f.reshape((nr*ng,nb))
        #        print('mumodel _add_  reshaped f.shape {}'.format(f.shape))
        return f

    def _load_multirun_grad_(self,minuit_ordered_grad_list):
        """
        probably broken in 2026, grads are slower that numerical ones
        
        maybe retry timing simplest fits
        minuit_ordered_grad_list is
        self._glocals_ to loop (broken) 
            over minuit parameters and calculate gradient components, 
                over runs, whose number n = y.shape[0], for globals
                on single run, for locals
                    over component parameters that contain that minuit parameter, i
                    including the derivative of the user function.
        constructed in mujpy.tools.tools  int2_multirun_grad_method_key
        """

        self._minuit_grad_list_ = minuit_ordered_grad_list


    def _add_multirun_grad_(self,*argv):
        """
        Deprecated

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
        _add_multirun_grad_ (glob is implicit, no grads in sequential fit) 
        """

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
        """
        Produces f for the fft of residues in mufitplot::

        input:
            x time array, 1d
            *argv Minuit parameters
          
            y - f
          
          Components can be selectively added to f
          i.e. subtracted in the residues::
   
            f += method(x,*pars) if self._include_components[j] else 0. 
            
          For the time being only single run single group (f is 1d)
        """

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
        """
        Generates partial residues for FFT

        input:
          include_components 
                True to subtract in residues  asymm - f
        """

        self._include_components = include_components 

    def _include_all_(self):
        """
        mudash reset to normal fit mode (initially of after fft)
        """

        self._include_components = [True]*self._ntruecomponents_

    def al(self,x,α):
        """
        alpha calibration

        x [mus], α
        x dummy, for compatibility
        """
        
        # empty method  (could remove x from argument list ?)
        # print('al = {}'.format(α))
        return []
        al.func_code = make_func_code(["α"])                
           
    def bl(self,x,A,λ): 
        """
        Lorentzian decay, A*exp(-x*λ)
        
        x [mus], A, λ [mus-1]
        """
        
        # x need not be self.x (e.g. in plot)
        # λ = -87. if λ < -87. else λ
        return A*exp(-x*λ)
        bl.func_code = make_func_code(["A","λ"])

    def _grad_bl_0_(self,x,A,λ): 
        """
        derivative of bl with respect to A in terms of self.bl
        
        x [mus], A, λ [mus-1]  
        """
        
        return self.bl(x,A,λ)/A

    def _grad_bl_1_(self,x,A,λ): 
        """
        derivative of bl with respect to λ in terms of self.bl
        
        x [mus], A, λ [mus-1]  
        """
        
        return -x*self.bl(x,A,λ)

    def bg(self,x,A,σ): 
        """
        Gaussian decay, A*exp(-0.5*(x*σ)**2)
        
        x [mus], A, σ [mus-1] (positive parity)
        """
        
        # x need not be self.x (e.g. in plot)        
        return A*exp(-0.5*(x*σ)**2)
        bg.func_code = make_func_code(["A","σ"])

    def _grad_bg_0_(self,x,A,σ): 
        """
        
        derivative of bg with respect to A in terms of self.bg
        x [mus], A, σ [mus-1]  
        """
        
        return self.bg(x,A,σ)/A
        
    def _grad_bg_1_(self,x,A,σ): 
        """
        derivative of bg with respect to σ in terms of self.bg
        
        x [mus], A, σ [mus-1]  
        """

        return -x**2*σ*self.bg(x,A,σ)

    def ba(self,x,A,λ,σ): 
        """
        Lorentzian times Gaussian decay, A*exp(-x*λ)*exp(-0.5*(x*σ)**2)
        
        x [mus], A, λ [mus-1], σ [mus-1] (positive parity)
        """
        
        # x need not be self.x (e.g. in plot)
        return A*exp(-x*λ)*exp(-0.5*(x*σ)**2)
        ba.func_code = make_func_code(["A","λ","σ"])

    def _grad_ba_0_(self,x,A,λ,σ): 
        """
        derivative of ba with respect to A in terms of self.ba
        
        x [mus], A, σ [mus-1]  
        """
        
        return self.ba(x,A,λ,σ)/A

    def _grad_ba_1_(self,x,A,λ,σ): 
        """
        derivative of ba with respect to λ in terms of self.ba
        
        x [mus], A, σ [mus-1]  
        """
        
        return -x*self.ba(x,A,λ,σ)

    def _grad_ba_2_(self,x,A,λ,σ): 
        """
        derivative of ba with respect to σ in terms of self.ba
        
        x [mus], A, σ [mus-1]  
        """
        
        return -x**2*σ*self.ba(x,A,λ,σ)

    def bs(self,x,A,Λ,β): 
        """
        stretched decay A*exp(-(x*Λ)**β), 
        
        x [mus], A, Λ [mus-1] (>0), β (>0)
        """
        
        # x need not be self.x (e.g. in plot)
        return A*exp(-(x*Λ)**β)
        bs.func_code = make_func_code(["A","Λ","β"])

    def _grad_bs_0_(self,x,A,Λ,β): 
        """
        derivative of bs with respect to A in terms of self.bs
        
        x [mus], A, Λ [mus-1] (>0), β (>0)  
        """
        
        return self.bs(x,A,Λ,β)/A

    def _grad_bs_1_(self,x,A,Λ,β): 
        """
        derivative of bs with respect to Λ in terms of self.bs
        
        x [mus], A, Λ [mus-1] (>0), β (>0)  
        """
        
        return -β/Λ*(Λ*x)**β*self.bs(x,A,Λ,β)

    def _grad_bs_3_(self,x,A,Λ,β): 
        """
        derivative of bs with respect to β in terms of self.bs
        
        x [mus], A, Λ [mus-1] (>0), β (>0)  
        """
        
        return -log(Λ*x)*(Λ*x)**β*self.bs(x,A,Λ,β)

    def ml(self,x,A,B,φ,λ): 
        """
        precession A cos(2 pi _gamma_Mu_MHzper_mT B x+φ _radeg_) times Lorentzian decay, 
        
        x [mus], A, B [mT], φ [deg], λ [mus-1]
        """
        
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)
        ml.func_code = make_func_code(["A","B","φ","λ"])

    def _derivative_ml_(self,x,A,B,φ,λ): 
        """
        derivative of mlwith respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times Lorentzian decay
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        """
        
        return -A*sin(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)

    def _grad_ml_0_(self,x,A,B,φ,λ): 
        """
        derivative of ml with respect to A in terms of self.ml and self._derivative_ml_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        """
        
        return self.ml(x,A,B,φ,λ)/A

    def _grad_ml_1_(self,x,A,B,φ,λ): 
        """
        derivative of ml with respect to B in terms of self.ml and self._derivative_ml_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        """
        
        return -2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_ml_(x,A,B,φ,λ)

    def _grad_ml_2_(self,x,A,B,φ,λ): 
        """
        derivative of ml with respect to φ in terms of self.ml and self._derivative_ml_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        """
        
        return -self._radeg_*self._derivative_ml_(x,A,B,φ,λ)

    def _grad_ml_3_(self,x,A,B,φ,λ): 
        """
        derivative of ml with respect to λ in terms of self.ml and self._derivative_ml_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        """
        
        return -x*self.ml(x,A,B,φ,λ)

    def mg(self,x,A,B,φ,σ): 
        """
        precession A cos(2 pi _gamma_Mu_MHzper_mT B x+φ _radeg_) times Gaussian decay, 
        
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  (positive parity)
        """
        
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-0.5*(x*σ)**2)
        mg.func_code = make_func_code(["A","B","φ","σ"])
        
    def _derivative_mg_(self,x,A,B,φ,σ): 
        """
        derivative of mg with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times Gaussian decay,
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        """
        
        return -A*sin(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-0.5*(x*σ)**2)

    def _grad_mg_0_(self,x,A,B,φ,σ): 
        """
        derivative of mg with respect to A in terms of self.mg and self._derivative_mg_
        
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        """
        
        return self.mg(x,A,B,φ,σ)/A
        
    def _grad_mg_1_(self,x,A,B,φ,σ): 
        """
        derivative of mg with respect to B in terms of self.mg and self._derivative_mg_
        
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        """
        
        return 2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_mg_(x,A,B,φ,σ)
        
    def _grad_mg_2_(self,x,A,B,φ,σ): 
        """
        derivative of mg with respect to φ in terms of self.mg and self._derivative_mg_
        
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        """
        
        return -self._radeg_*self._derivative_mg_(x,A,B,φ,σ)
        
    def _grad_mg_3_(self,x,A,B,φ,σ): 
        """
        derivative of mg with respect to σ in terms of self.mg and self._derivative_mg_
        
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        """
        
        return -x**2*σ*self.mg(x,A,B,φ,σ)
        
    def mu(self,x,A,B,φ,λ,σ): 
        """
        precession A cos(2 pi _gamma_Mu_MHzper_mT B x+φ _radeg_) times Gaussian times Lorentzian decays, 
        
        x [mus], A, B [mT], φ [degrees], 
        λ [mus-1], σ [mus-1]  (positive parity)
        """
        
        # x need not be self.x (e.g. in plot)
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)*exp(-0.5*(x*σ)**2)
        mu.func_code = make_func_code(["A","B","φ","λ","σ"])

    def _derivative_mu_(self,x,A,B,φ,λ,σ): 
        """
        
        derivative of mu with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times Lorentzian times Gaussian decay,
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        """
        
        return -A*sin(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)*exp(-0.5*(x*σ)**2)

    def _grad_mu_0_(self,x,A,B,φ,λ,σ): 
        """
        derivative of mu with respect to A in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        """
        
        return self.mu(x,A,B,φ,λ,σ)/A

    def _grad_mu_1_(self,x,A,B,φ,λ,σ): 
        """
        
        derivative of mu with respect to B in terms of self.mu and self._derivative_mu_
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        """
        
        return -2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_mu_(x,A,B,φ,λ,σ)

    def _grad_mu_2_(self,x,A,B,φ,λ,σ): 
        """
        derivative of mu with respect to φ in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        """
        
        return -self._radeg_*self._derivative_mu_(x,A,B,φ,λ,σ)

    def _grad_mu_3_(self,x,A,B,φ,λ,σ): 
        """
        derivative of mu with respect to λ in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        """
        
        return -x*self.mu(x,A,B,φ,λ,σ)

    def _grad_mu_4_(self,x,A,B,φ,λ,σ): 
        """
        derivative of mu with respect to σ in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        """
        
        return -x**2*σ*self.mu(x,A,B,φ,λ,σ)

    def ms(self,x,A,B,φ,Λ,β): 
        """
        precession A cos(2 pi _gamma_Mu_MHzper_mT B x+φ _radeg_) times stretched decay, 
        
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        """
        
        # x need not be self.x (e.g. in plot)
        return A*cos(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-(x*Λ)**β)
        ms.func_code = make_func_code(["A","B","φ","Λ","β"])

    def _derivative_ms_(self,x,A,B,φ,Λ,β): 
        """
        derivative of ms with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times stretched decay
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        """
        
        return -A*sin(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-(x*Λ)**β)

    def _grad_ms_0_(self,x,A,B,φ,Λ,β): 
        """
        derivative of ms with respect to A in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0) 
        """
        
        return self.ms(x,A,B,φ,Λ,β)/A

    def _grad_ms_1_(self,x,A,B,φ,Λ,β): 
        """
        derivative of ms with respect to B in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0) 
        """
        
        return -2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_ms_(x,A,B,φ,Λ,β)

    def _grad_ms_2_(self,x,A,B,φ,Λ,β): 
        """
        derivative of ms with respect to φ in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        """
        
        return -self._radeg_*self._derivative_ms_(x,A,B,φ,Λ,β)

    def _grad_ms_3_(self,x,A,B,φ,Λ,β): 
        """
        derivative of ms with respect to Λ in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        """
        
        return -β/Λ*(Λ*x)**β*self.ms(x,A,B,φ,Λ,β)

    def _grad_ms_4_(self,x,A,B,φ,Λ,β): 
        """
        derivative of ms with respect to β in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        """
        
        return -log(Λ*x)*(Λ*x)**β*self.ms(x,A,B,φ,Λ,β)


    def fm(self,x,A,B,λ):
        """
        FmuF (powder average)
        
        according to Book  
        x [mus], A, B [mT], λ [mus-1]
        B is Bdip
        """
        
        # x need not be self.x (e.g. in plot)
        return A/6.0*(1.+cos(2*pi*self._gamma_Mu_MHzper_mT*B*x)+
               2.*(cos(pi*self._gamma_Mu_MHzper_mT*B*x)+
                   cos(3*pi*self._gamma_Mu_MHzper_mT*B*x) ))*exp(-x*λ)
        fm.func_code = make_func_code(["A","B","λ"])

    def jl(self,x,A,B,φ,λ): 
        """
        Bessel j0 precession times Lorentzian decay, 
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1]
        """
        # x need not be self.x (e.g. in plot)
        
        return A*j0(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)
        jl.func_code = make_func_code(["A","B","φ","λ"])

    def _derivative_jl_(self,x,A,B,φ,λ): 
        """
        derivative of jl with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        
        - A J1(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times Lorentzian decay,
        x [mus], A, B [mT], φ [degrees], λ [mus-1]
        """
        
        return -A*j1(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)

    def _grad_jl_0_(self,x,A,B,φ,λ): 
        """
        derivative of jl with respect to A in terms of self.ml and self._derivative_ml_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        """
        
        return self.jl(x,A,B,φ,λ)/A

    def _grad_jl_1_(self,x,A,B,φ,λ): 
        """
        derivative of jl with respect to B in terms of self.ml and self._derivative_ml_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        """
        
        return -2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_jl_(x,A,B,φ,λ)

    def _grad_jl_2_(self,x,A,B,φ,λ): 
        """
        derivative of jl with respect to φ in terms of self.ml and self._derivative_ml_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        """
        
        return -self._radeg_*self._derivative_jl_(x,A,B,φ,λ)

    def _grad_ml_3_(self,x,A,B,φ,λ): 
        """
        derivative of ml with respect to λ in terms of self.ml and self._derivative_ml_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1]  
        """
        
        return -x*self.jl(x,A,B,φ,λ)


    def jg(self,x,A,B,φ,σ): 
        """
        Bessel j0 precession times Gaussian decay, 
        
        x [mus], A, B [mT], φ [degrees], σ [mus-1] (positive parity)
        """
        
        # x need not be self.x (e.g. in plot)
        
        return A*j0(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-0.5*(x*σ)**2)
        jg.func_code = make_func_code(["A","B","φ","σ"])

    def _derivative_jg_(self,x,A,B,φ,σ): 
        """
        derivative of jg with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times Gaussian decay,
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        """
        
        return -A*j1(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-0.5*(x*σ)**2)

    def _grad_jg_0_(self,x,A,B,φ,σ): 
        """
        derivative of jg with respect to A in terms of self.mg and self._derivative_mg_
        
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        """
        
        return self.jg(x,A,B,φ,σ)/A
        
    def _grad_jg_1_(self,x,A,B,φ,σ): 
        """
        derivative of jg with respect to B in terms of self.mg and self._derivative_mg_
        
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        """
        
        return 2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_jg_(x,A,B,φ,σ)
        
    def _grad_jg_2_(self,x,A,B,φ,σ): 
        """
        derivative of jg with respect to φ in terms of self.mg and self._derivative_mg_
        
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        """
        
        return -self._radeg_*self._derivative_jg_(x,A,B,φ,σ)
        
    def _grad_jg_3_(self,x,A,B,φ,σ): 
        """
        derivative of jg with respect to σ in terms of self.mg and self._derivative_mg_
        
        x [mus], A, B [mT], φ [degrees], σ [mus-1]  
        """
        
        return -x**2*σ*self.jg(x,A,B,φ,σ)
        
    def j0(self,x,A,B,φ,λ,σ): 
        """
        precession A j1(2 pi _gamma_Mu_MHzper_mT B x+φ _radeg_) times Gaussian times Lorentzian decays, 
        
        x [mus], A, B [mT], φ [degrees], 
        λ [mus-1], σ [mus-1]  (positive parity)
        """
        
        # x need not be self.x (e.g. in plot)
        return A*j0(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)*exp(-0.5*(x*σ)**2)
        mu.func_code = make_func_code(["A","B","φ","λ","σ"])

    def _derivative_j0_(self,x,A,B,φ,λ,σ): 
        """
        derivative of j0 with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times Lorentzian times Gaussian decay,
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        """
        
        return -A*j1(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-x*λ)*exp(-0.5*(x*σ)**2)

    def _grad_j0_0_(self,x,A,B,φ,λ,σ): 
        """
        derivative of j0 with respect to A in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        """
        
        return self.j0(x,A,B,φ,λ,σ)/A

    def _grad_j0_1_(self,x,A,B,φ,λ,σ): 
        """
        derivative of j0 with respect to B in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        """
        
        return -2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_j0_(x,A,B,φ,λ,σ)

    def _grad_j0_2_(self,x,A,B,φ,λ,σ): 
        """
        derivative of j0 with respect to φ in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        """
        
        return -self._radeg_*self._derivative_j0_(x,A,B,φ,λ,σ)

    def _grad_j0_3_(self,x,A,B,φ,λ,σ): 
        """
        derivative of j0 with respect to λ in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        """
        
        return -x*self.j0(x,A,B,φ,λ,σ)

    def _grad_j0_4_(self,x,A,B,φ,λ,σ): 
        """
        derivative of j0 with respect to σ in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], λ [mus-1], σ [mus-1]  
        """
        
        return -x**2*σ*self.j0(x,A,B,φ,λ,σ)


    def js(self,x,A,B,φ,Λ,β): 
        """
        Bessel j0 precession times stretched decay, 
        
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        """
        
        # x need not be self.x (e.g. in plot)
        return A*j0(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-(x*Λ)**β)
        js.func_code = make_func_code(["A","B","φ","Λ","β"])

    def _derivative_js_(self,x,A,B,φ,Λ,β): 
        """
        derivative of js with respect to total phase alpha =  2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_,
        
        - A sin(2 pi _gamma_Mu_MHzper_mT B x + φ _radeg_) times stretched decay
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        """
        
        return -A*j1(2*pi*self._gamma_Mu_MHzper_mT*B*x+φ*self._radeg_)*exp(-(x*Λ)**β)

    def _grad_js_0_(self,x,A,B,φ,Λ,β): 
        """
        derivative of js with respect to A in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0) 
        """
        
        return self.js(x,A,B,φ,Λ,β)/A

    def _grad_js_1_(self,x,A,B,φ,Λ,β): 
        """
        derivative of js with respect to B in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0) 
        """
        
        return -2*pi*self._gamma_Mu_MHzper_mT*x*self._derivative_js_(x,A,B,φ,Λ,β)

    def _grad_js_2_(self,x,A,B,φ,Λ,β): 
        """
        derivative of js with respect to φ in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        """
        
        return -self._radeg_*self._derivative_js_(x,A,B,φ,Λ,β)

    def _grad_js_3_(self,x,A,B,φ,Λ,β): 
        """
        derivative of ms with respect to Λ in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        """
        
        return -β/Λ*(Λ*x)**β*self.js(x,A,B,φ,Λ,β)

    def _grad_js_4_(self,x,A,B,φ,Λ,β): 
        """
        derivative of js with respect to β in terms of self.mu and self._derivative_mu_
        
        x [mus], A, B [mT], φ [degrees], Λ [mus-1] (>0), β (>0)
        """
        
        return -log(Λ*x)*(Λ*x)**β*self.js(x,A,B,φ,Λ,β)
        
# kubo toyabe and fm gradients not implemented

    def _kg(self,t,w,Δ):
        """
        auxiliary component for a static Gaussian Kubo Toyabe in longitudinal field, 
        
        t [mus], w [mus-1], Δ [mus-1], 
        w = 2*pi*gamma_mu*L_field
        The first derivative of dawsn(x) is 1-2*x*dawsn(x)
        """
        
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
        """
        static Lorentzian Kubo Toyabe in longitudinal field, 
        
        t [mus], w [mus-1], Δ [mus-1], 
        w = 2*pi*gamma_mu*L_field
        """
        
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
        """ 
        auxiliary dynamization of Gaussian Kubo Toyabe by G. Allodi 
        
        N: number of sampling points;
        dt: time interval per bin [i.e. time base is t = dt*(0:N-1)]
        w [mus-1], Δ [mus-1], ν [MHz] 
        (longitudinal field freq, Gaussian distribution, scattering frequency 
        % alphaN: [optional argument] weighting coefficient alpha times N. Default=10 
        """
        
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
        """
        Gauss Kubo Toyabe in (fixed) long field, static or dynamic
        
        x [mus], A, BL [mT], Δ [mus-1] (positive parity), ν (MHz)
        """
        
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
           # WARNING! was is inspace)0,Np,,Np+1= but what is inspace?
           n = hstack((linspace(0,Np,Np+1),linspace(Nm,-1.,-Nm))) # 
           f = fft.ifft(fft.fft(f)*exp(nshift*1j*2*pi*n/Ns)) # shift back
        # multiply by amplitude
        f = A*real(f[0:N])
        return f
        kg.func_code = make_func_code(["A","BL","Δ","ν"])

    def kl(self,x,A,BL,Γ):
        """
        Lorent Kubo Toyabe in (fixed) long field, static 
        
        x [mus], A, BL [mT], Γ [mus-1] 
        """
        
        # x need not be self.x (e.g. in plot)
        # (dynamic makes no sense)
        w = 2*pi*BL*self._gamma_Mu_MHzper_mT
        return A*self._kl(x,w,Γ)
        kl.func_code = make_func_code(["A","BL","Γ"])

    def kd(self,x,A,Δ,λ):
        """
        Gauss Kubo Toyabe static times Lorentz decay
        
        x [mus], A, B [T], Δ [mus-1], ν (MHz)
        """
        
        # x need not be self.x (e.g. in plot)
        return A*self._kg(x,0,Δ)*exp(-x*λ)
        kd.func_code = make_func_code(["A","Δ","λ"])
        #kd.limits = [[None,None],[0.,None],[None,None]]
        #kd.error = [0.002,0.05,0.05]

    def ks(self,x,A,Δ,Λ,β):
        """
        Gauss Kubo Toyabe times stretched decay
        
        x [mus], A, B [T], Δ [mus-1], Λ [mus-1] (>0), β (>0)
        """
        
        # x need not be self.x (e.g. in plot)
        return A*self._kg(x,0,Δ)*exp(-(x*Λ)**β)
        ks.func_code = make_func_code(["A","Δ","Λ","β"])
        #kd.limits = [[None,None],[0.,None],[None,None]]
        #kd.error = [0.002,0.05,0.05]

    def _chisquare_(self,*argv,debug=False):
        """
        Provides chisquares, partial over individual runs or groups if self._axis_ = 1 
        
        None is default and sum is over all indices::
        Signature provided at Minuit invocation by 
        optional argument forced_parameters=parnames
        where parnames is a tuple of parameter names::

           e.g. parnames = ('asym','field','phase','rate') 

        Works also for global fits, 
        where sum (...,axis=None) yields the sum over all indices.
        """ 

        # Mepsi = finfo('d').max/10.
        num = self._add_plain_(self._x_,*argv)- self._y_
        squaredev = (num/self._e_)**2
        if debug:
            return sum(squaredev), num, self._e_
        else:
            return sum(squaredev,axis=self._axis_ )

    def _chisquare_calib_(self,*argv):
        """
        Provides chisquares calib, partial over individual runs or groups if self._axis_ = 1 

        None is default and sum is over all indices::
        Signature provided at Minuit invocation by 
        optional argument forced_parameters=parnames
        where parnames is a tuple of parameter names::

           e.g. parnames = ('asym','field','phase','rate') 

        Works also for global fits, 
        where sum (...,axis=None) yields the sum over all indices.
        """

        # print('_chisquare_ mucomponents debug: {} {} {}'.format(self._x_.shape,self._y_.shape,self._e_.shape))
        from numpy import abs 
        # Mepsi = finfo('d').max/10.
        num = abs(self._add_calib_(self._x_,*argv) - self._y_) # self._t_ is set to suite.time by _load_calib_ 
                                                         # _add_calib_ generates self._x_ self._y_ self._e_  
#        self._ndata = self._y_.size
        normsquaredev = (num/self._e_)**2
        # divergence = normsquaredev>Mepsi
#        if divergence.any():
#            print('Warning: big numbers in chisquare {}'.format(normsquaredev[divergence]))
        return sum(normsquaredev,axis=self._axis_ )

    def _grad_chisquare_(self,*argv):
        """
        option for global multirun fits, sum (...,axis=None) yields the sum over all indices.

        Provides gradient of chisquare with respect to p, i.e. along the i-th parameter p_i it is
        sum_j 2[y(t_j,p)-y_ej)]/e_j^2 sum_k d y_k(t_j,p) /dp_i
        where j are bins, k are components in the model
        
        The first factor is common to a all grad components
        The second factor mus be selected. 
        y_k may not depend on p_i, dy_k/dp_i = 0
        And if y_k' does depend, p_i will be its l-th parameter, and we must use the l-th component of its gradient
        """ 

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

