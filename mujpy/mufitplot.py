class mufitplot(object):
    """
    class to plot fit guess/results and fft the residues

    input 
       plot_range: 2 to 5 values
          None (skip time plot)
          start stop
          start stop pack
          start stopearly packearly stop packlate
       the_fit: the mufit instance, 
          must be run before 
       guess: True, False (default)
       rotating_frame_frequencyMHZ: 0.0 default
       fit_range: 1 or 3 values
          None default
          start, stop, sigma (LB) in MHz  
    mufitplot class 
    produces the fitplot
    multiple runs sequential produce anim
    can plot in the rotating frame
    optional fft of residues toggles
    fig handle is now handled automatically 
    """

# A bit of discussion to be removed when done
# plot needs to access a version of mumodel, 
# that already stored in the_fit,
# with different ranges for rebinned data, 
# denser function plot, original data
# mumodel() is needed for generating model function
# and calculating chi2.
# everything can be obtained from the_fit, including dashboard initial guess
# Maybe make another version that does multiplots from the dashboards
#
#
#  the function that performs the plot is mujpy.tools.plot import set_single_fit, self_sequence_fit 
#        invoked by self.plot_run in 
#        plot_singlerun or plot_singlerun_multigroup_sequential or plot_singlerun_multigroup_userpar

    def __init__(self,plot_range,the_fit,guess=False,rotating_frame_frequencyMHz = 0.0,fft_range = None,real=True,fig_fit=None,fig_fft=None):
        """
        init args are plot_range and a mufit instance [options]

        input:
            plot_range a csv string, start, stop [,pack [start_r, stop_r]]  
            the_fit mufit instance 
                    self.dashboard = the_fit.dashboard best fit json file
                    mufit(...,initialize_only=True) for dry Minuit to plot guess
            [guess = True for initial guess plot, default False for Minuit results]
            [rotating_frame_frequencyMHz, ditto if !=0.0 (default)]
            [fft_range ditto in !None (default)]
            [real = False for fft power]
            finally calls self.fig = set_figure(self.fig, 
                                    self.model_name(), class method ok
                                    early_late, True/False, 3 or 5 integers in plot_range
                                    data, [t, y, ey, f_res, tf, f, dy_fit], data from the_model, f, f_res from chi
                                    group, [fgroup, bgroup, alpha] trivial list over groups
                                    run_title, from tools get_title
                                    chi_dof, [nu_fit,chi_fit] from chi
                                    data_late, [t_late,y_late,ey_late,f_late_res,tfl,fl,dy_fit_early,dy_fit_late] as data
                                    chi_dof_late,rrf=self.rotating_frame_frequencyMHz)
            fig_fit fig_handle for handle reuse in mudash]
        """

#        from mujpy.mucomponents.mucomponents import mumodel
        
# reassigned, for backward compatibility  
        # print('mufitplot __init__ debug: inside')  
        from mujpy._version import __version__
        self.__version__ = __version__
        self.fit_types = the_fit.fit_types
        self.calib = the_fit.calib()
        self.global_fit = the_fit.global_fit()
        self.model = the_fit.the_model
        self.suite = the_fit.suite
        self.log = the_fit.log
        self.lastfit = the_fit.lastfit
        self.lastfits = the_fit.lastfits
        #print('MUFITPLOT __INIT__')
        #for fit in self.lastfits:
        #    print(fit)
        #self.add = the_fit.the_model._add_plot_

#        # can do without pars class argument by
#        A1,A20,A21,B1,B20,B21,C1,C2 = self.fit_types()
#        if A1 or A21 or C1 or C2: # single lastfits  
#            pars = self.lastfit
#        else: # A20 or B1 or B20 or B21, multiple lastfits
#            pars = []
#            for lastfit in self.lastfits:
#                values = lastfit.values                
#                pars.append(values) # list of list
        if not self.suite.loadfirst: 
            self.log('Sorry, no access to data... quitting mufitplot')
            ok = False
        else:
            ok = True 
            self.dashboard = the_fit.dashboard 
            # print(self.dashboard)
        # print('mufitplot __init__ debug: the_fit.C1() =\n{}'.format(self.fit.C1()))
        self.fig = fig_fit
        self.fig_fft = fig_fft

        if ok: 
            self.rotating_frame_frequencyMHz = rotating_frame_frequencyMHz
            self.guess = guess

# used in choosefftplot, non really needed
            if "userpardicts_guess" in self.dashboard.keys():
                self.jsonmodel = "global results"
            self.jsonmodel = "model_guess" 
            if not self.guess:
                if "userpardicts_result" in self.dashboard.keys():   
                    self.jsonmodel = "userpardicts_result"
                elif "model_result" in self.dashboard:
                    self.jsonmodel = "model_result"  
           
            if plot_range is not None:
#                self.chooseplot(plot_range)
                self.plot_run(plot_range)
                if fft_range is not None: # fft always plots fit as well
                    self.choosefftplot(fft_range,real)
            else:
                if fft_range is None:
                    self.log('Exiting mufitplot: nothing to do.')    
         
    def plot_run(self,plot_range):
        """
        recovers data and statistics from mufit.mumodel and plots fit result/guess
        
        input :
            plot_range
                (start,stop)
                (start,stop,pack)
                (start_early,stop_early,pack_early,stop_late,pack_late)
        calls
            self.chi, which calls self.reload((start,stop,pack))
                            uploads rebinned data slice in self.model._x_, self.model._y_, self.model._e_
                            produces ndof, f, chi for the chosen slice  (through model)
        finally calls 
            set_single_fit or set_sequence_fit, from tools.plot
                                    (produce actual figures)
            Draws either static or anim figures using tools.plot functions

            Data are rebinned with reload, using mumodel methods.
            Rot frame: self.rrf_y, self.rrf_e are unbinned mixed data/errors
            Here they are just rebinned and their fit with same pack is mixed
        """

        from mujpy.tools.tools import derange, rebin, rshp
        from mujpy.tools.tools import get_run_title,  mixer
        from mujpy.tools.plot import set_single_fit, set_sequence_fit
        from numpy import ones, tile 

        ############################àààààààààààààààààààààà###
        # standard fit tools.plot has four main subplots    #
        #     plot_run passes the following lists           #
        #   data = [t, y, ey, f_res, tf, f, dy_fit]         #
        #   chi_dof = [nu_fit, chi_fit]                     #
        #   data_late = [t_late,y_late,ey_late,             #
        #                f_late_res,tfl,fl,                 #
        #                dy_fit_early,dy_fit_late]          #
        #   chi_dof_late = [nu_fit_early, nu_fit_late,      #
        #                   chi_fit_early, chi_fit_late,    #
        #                   w_early, W_late]                #
        #------------   subplots description    ----------  #
        # 1) data and fit [split into early and late]       #
        #                [if early_late = True]a            #
        # 1a)   data     errorbar(t,y,yerr=ey)              #
        #       fit      plot(tf,f)                         #    
        # 1b)   data errorbar(t_late,y_late,yerr=ey_late)   #
        #       fit      plot(tfl,fl)                       #
        # 2) residues, below data [also early_late]         #
        # 2a)   errorbar(t,y - fres,yerr=ey)                #
        # 2b)   errorbar(t_late,y_late-f_late_res)          #
        # 3) residues histogram                             #
        #       histogram(dy_fit[i])                        #
        #      [histogram(dy_fit_early[i],weights=w_early)  #
        #       histogram(dy_fit_late[i],weights=w_late)]   #
        # 4) info: [partial] chisquares, dof                #
        #               chi_fit and nu                      #
        #####################################################

        run_title = get_run_title(self.suite) # always a list, even for single
        #A1,A20,A21,B1,B20,B21,C1,C2 = self.fit_types()
        #string = 'global ' if A21 or B21 or C1  or C2 else ''
        string = 'global ' if self.global_fit else ''
        if self.guess:
            run_title = [title + ": "+string+"guess values" for title in run_title]
        else:
            run_title = [title + ": "+string+"fit results" for title in run_title]

# self.chi is the main method below, for a given tuple
# returns ndof fit function, chi2 
# reloads self.model._x_, self.model._y_, self.model._e_ -> local variables
#
# fittup: t_fit,y_fit_e_fit,f_fit for dy_fit=(y_fit-f_fit)/ey_fit, for chi2 histograms
# if early_late 
# lattup: t_late,y_late,ey_late,f_late_res  -> tile(t_late) RF(y_late,ey_late,f_late_res)  if freq > 0
# fltup: tfl, fl (fit plot late) -> tile(tfl) RF(fl) if freq > 0,
# flattup:  _,y_fit_late,ey_fit_late,f_fit_late for dy_fit_late, for chi2 histograms
# else
# plottup: t,y,ey,f_res -> tile(t) RF(y,e,f) if freq > 0
# plotftup: tf, f (fit plot [early]) -> tile(tf) RF(f) if freq > 0

# before reloading at the end of mufit, self.model._a_ self.model_ea_ are the unbinned asymm, asyme
        t = self.suite.time
        if self.rotating_frame_frequencyMHz: 
        # preload unbinned rotating frame data in self.rrf_y self.rrf_e
        # all other rrf frames obtained by rebinning rrf_y, rrf_e
            time = tile(t,self.model._a_.shape[0:]+(1,)) if len(self.model._a_.shape)>len(self.model._t_.shape) else t
            # tiles self.model._t_ to the shape of self.model._a_
            self.rrf_y = mixer(time,self.model._a_,self.rotating_frame_frequencyMHz)
            self.rrf_e = mixer(time,self.model._ea_,self.rotating_frame_frequencyMHz)
            # self.rrf_y, self.rrf_e are rotating frame unbinned data, stds
        fittup,ermsg = derange(self.dashboard["fit_range"],self.suite.histoLength)# range as tuple
        fittup = fittup if len(fittup) == 3 else (fittup[0],fittup[1],1) # if pack is missing
        fit_start,fit_stop,fit_pack = fittup[0:3] # only one range, trailing in longer fittup ignored
        
        # t_fit,y_fit,ey_fit,f_fit
        nu_fit, f_fit, chi_fit = self.chi(fittup) # fit ndof, f, chi2 & rebins self.model._x_, self.model._y_, self.model._e_
        t_fit,y_fit,ey_fit = rshp(self.model._x_), rshp(self.model._y_), rshp(self.model._e_) # for dy_fit = (y_fit-y_fit)/ey_fit

        if self.rotating_frame_frequencyMHz:
            time = tile(t_fit,f_fit.shape[0:]+(1,)) if len(self.model._a_.shape)>len(self.model._t_.shape) else t_fit
            # tiles
            _, y_fit, ey_fit = rebin(self.suite.time,self.rrf_y,[fittup[0],fittup[1]],fittup[2],e=self.rrf_e)  # y_fit ey_fit in rot frame
            f_fit = mixer(time,f_fit,self.rotating_frame_frequencyMHz)                    # fit in rot frame

        plottup,ermsg = derange(plot_range,self.suite.histoLength)

        if len(plottup)==5: 
        #########################################
        # double range plot                     #
        #########################################
            early_late = True
            start, stop, pack, last, packlate = plottup # required by rebin rrf late
        
            # t_late,y_late,ey_late,f_late_res
            lattup = (stop,last,packlate)
            _, f_late_res, _  = self.chi(lattup) # load rebinned rshp(self.model._x_), rshp(self.model._y_), rshp(self.model._e_)
            # f_late_res for calculating residues on data points 
            t_late,y_late,ey_late = rshp(self.model._x_), rshp(self.model._y_), rshp(self.model._e_) # rebinned late data slice         

            # tfl,fl for late plot curve
            packfit = int(packlate/2)
            fltup = (stop, last, packfit)
            _,fl,_ = self.chi(fltup)
            tfl = rshp(self.model._x_)
            
            #plotfit = (start, stop, packfit)
            #_, fl, _ = self.chi(plotfit) # fl for plotting function  
            if self.rotating_frame_frequencyMHz: # (t_late,y_late,f_late_res,ey_late) the last three must be transformed to rrf
                if len(self.model._a_.shape)>len(self.model._t_.shape):
                    time = tile(t_late,f_fit.shape[0:]+(1,)) # tiles
                    timf = tile(tfl,f_fit.shape[0:]+(1,)) # tiles
                else:
                    time, timf = t_late, tfl
                _,y_late,ey_late = rebin(self.suite.time,self.rrf_y,[start,stop],packlate,e=self.rrf_e)
                f_late_res = mixer(time, f_late_res, self.rotating_frame_frequencyMHz) 
                fl = mixer(timf, fl, self.rotating_frame_frequencyMHz) 
            
            fit_late_start = int(stop/packlate*fit_pack) # divides fit_range in early and late
            fitlattup = (fit_late_start,fit_stop,fit_pack) # for dy_fit_late histogram

            # redo the same with the original fit binning fit_pack, only for histo and chi2
            nu_fit_late, f_fit_late, chi_fit_late = chi(fitlattup) # chi(t_fit_late,y_fit_late,ey_fit_late,pars)
            _,y_fit_late,ey_fit_late = rshp(self.model._x_),rshp(self.model._y_),rshp(self.model._e_)        

#            nu_fit_early, f_fit_early, chi_fit_early = chi(plottup[:3]) #chi(t_fit_early,y_fit_early,ey_fit_early,pars)
#            t_fit_early,y_fit_early,ey_fit_early = rshp(self.model._x_), rshp(self.model._y_), rshp(self.model._e_) 
#            
#            # single slices as in fit 
#            if self.rotating_frame_frequencyMHz:
#                _,y_fit_early,ey_fit_early = rebin(self.suite.time,self.rrf_y,[plottup[0],plottup[1]],plottup[2],e=self.rrf_e)
#                f_fit_early = mixer(t_fit_early, f_fit_early,self.rotating_frame_frequencyMHz) 
        else:
        ###################
        # single range plot
        ###################
            early_late = False
            pack = 1
            chi_fit_late, nu_fit_late, chi_fit_early, nu_fit_early = None, None, None, None
            if len(plottup)==2: # plot start stop
                start, stop = plottup
                plottup = (plottup[0],plottup[1],1)
            else:
                start, stop, pack = plottup[:3]
#        self.log('plot_range= {}'.format(plot_range))
#        self.log('start, stop, pack = {},{},{}'.format(start, stop, pack))

        # the following serve both for full range and for early range
        _,f_res, _ = self.chi(plottup[:3]) # chi(t,None,None,pars)
        t,y,ey = rshp(self.model._x_), rshp(self.model._y_), rshp(self.model._e_) 

        packfit =  int(pack/2)
        plotftup = (start, stop,packfit)
        _,f,_ = self.chi(plotftup)
        tf,_,_ = rshp(self.model._x_), rshp(self.model._y_), rshp(self.model._e_)# rebin(self.suite.time,asymm,[start,stop],packfit,e=asyme)
#        else:
        # print('mufitplot plot_run debug: min max f_res {}, min max f {}'.format([f_res.min(),f_res.max()],[f.min(),f.max()])) 
        if self.rotating_frame_frequencyMHz: # (t,y,f_res,ey) the last three must be transformed to rrf
                                             # (tf,f) non rebinned, the last must be transformed to rrf
            if len(self.model._a_.shape)>len(self.model._t_.shape):
                time = tile(t,y.shape[0:]+(1,)) # tiles
                timf = tile(tf,f.shape[0:]+(1,)) # tiles
            else:
                time, timf = t_late, tfl
            _,y,ey = rebin(self.suite.time,self.rrf_y,[start,stop],pack,e=self.rrf_e)
            f_res = mixer(time, f_res,self.rotating_frame_frequencyMHz) 
            f = mixer(timf, f,self.rotating_frame_frequencyMHz) 
        # assume self.suite.single() 
        # prepare figure  
        fgroup = [self.suite.groups[k]['forward'] for k in range(len(self.suite.groups))]
        bgroup = [self.suite.groups[k]['backward'] for k in range(len(self.suite.groups))]
        alpha = [self.suite.groups[k]['alpha'] for k in range(len(self.suite.groups))]
        group = [fgroup, bgroup, alpha]
        dy_fit = (y_fit-f_fit)/ey_fit

        ########################################
        # data passed to set_figure_fit: 
        # exp data, fcn for residue plot  (t,y,ey,f_res)
        # time and fcn for model plot (tf,f)
        data = [t, y, ey, f_res, tf, f, dy_fit]

        chi_dof =[nu_fit,chi_fit]

        if self.suite.single() and not self.suite.multi_groups():
        # single draws a static plot
            set_figure_fit = set_single_fit
        else:
        # many runs are dealt by animation
            set_figure_fit = set_sequence_fit
        # both plot data, fit, residues, chisquare stat histo, info
        #       by two separate methods, from mujpy.tools.tools.plot
        #       The static method expects 1d y, ey, f, dy arrays
        #       The anim method expects 2d y, ey, f, dy arrays
        #       Reshape in 

        if early_late:
            dy_fit_early = (y_fit_early-f_fit_early)/ey_fit_early
            dy_fit_late = (y_fit_late-f_fit_late)/ey_fit_late            
            w_early  = nu_fit*ones(t_fit_early.shape[0])/nu_fit_early
            w_late = nu_fit*ones(t_fit_late.shape[0])/nu_fit_late
            data_late = [t_late,y_late,ey_late,f_late_res,tfl,fl,dy_fit_early,dy_fit_late]
            chi_dof_late = [nu_fit_early,nu_fit_late,
                            chi_fit_early,chi_fit_late,
                            w_early,w_late]
        else:
            data_late, chi_dof_late = None, None
        # self.log('Debug-mufitplot: set_figure  = {}'.format(set_figure_fit))
        self.fig = set_figure_fit(  self.fig,self.model_name(),
                                    early_late,
                                    data,
                                    group,
                                    run_title,
                                    chi_dof,
                                    data_late,
                                    chi_dof_late,rrf=self.rotating_frame_frequencyMHz)#,canvas=self.suite.canvas)
        return True, ''
    
    def chi(self,tup):
        """
        reloads data and statistics through mufit.mumodel

        input:
            tup = start, stop, pack
        output:
            ndof number of degrees of freedom [*approximated] for each 1d slice of self.model._y_
            f fit function, same dim as self.model._y_
            chi2_r [*approximate] reduced chisquare of each 1d slice 
        valid for any fit, A1, A20, A21*, B1, B20, B21*, C1*, C2*
                                                  and their calib
        """
        
        from mujpy.tools.tools import rshp

        self.model._rebin_plot_(tup,self.lastfits) # reloads start stop pack into mumodel and reloads _x_, _y_, _e_ which are linked 
        #print('mufitplot chi after _rebin_plot_ self.model._y_.shape {}'.format(self.model._y_.shape))
        #   in mufitplot as self.model._x_, self.model._y_, self.model._e_
        # now also calib is rebinned
        f = rshp(self.model._add_plot_(self.model._x_,self.lastfits)) # self.model._add_plain/_add_calib_
        #print('mufitplot chi after _add_plot_ self.model._y_.shape {}'.format(self.model._y_.shape))
        #print('chi f.shape = {}'.format(f.shape))
        A1,A20,A21,B1,B20,B21,C1,C2 = self.fit_types()
        n_runs, n_groups = self.suite.nruns, len(self.suite.groups) # suite dimensions
        if A1: # 1d yin and f, par is a simple list
                    #self.model._x_,self.model._y_,self.model._e_ = self.fit._x_, self.fit._y_, self.fit._e_
        #print('mufitplot chi f.shape {}, self.model._x_.shape {}, self.model._y_.shape {}'.format(f.shape,self.model._x_.shape,self.model._y_.shape))
        #self.model._x_ = self.model._x_
            #self.model._y_ = self.model._y_
            #self.model._e_ = self.model._e_
#            print('f.shape {}, self.model._x_.shape {} elf.y.shape {} self.model._e_.shape {}'.format(f.shape,self.model._x_.shape,self.model._y_.shape,self.model._e_.shape))
#            self.model._x_ = self.model._x_.squeeze()
            #self.model._y_ = self.model._y_.squeeze()
            #self.model._e_ = self.model._e_.squeeze()
            ndof = self.model._y_.size - self.lastfit.nfit # degrees of freedom
            # chi2 = self.fit.lastfit.fval/ndof
            chi2_r = (((self.model._y_-f)/self.model._e_)**2).sum()/ndof # exact
        elif A21: # 2d yin and f, par is a list
            n_groups = f.shape[0]
            ndof = self.model._y_.shape[-1] - round(self.lastfit.nfit/n_groups) # approximate
            chi2_r = []
            for kgroup in range(n_groups):
                chi2 = (((self.model._y_[kgroup,:]-f[kgroup,:])/self.model._e_[kgroup,:])**2).sum()/ndof
                chi2_r.append(chi2) # list of n_groups approximate chi2_r
        elif A20 or B1: # suite 2d yin and f, par is a list of lists
        #self.model._x_,self.model._y_,self.model._e_ = self.fit._x_, self.fit._y_, self.fit._e_
        #print('mufitplot chi f.shape {}, self.model._x_.shape {}, self.model._y_.shape {}'.format(f.shape,self.model._x_.shape,self.model._y_.shape))           n_groups = f.shape[0]
            # print('mufitplot chi self.model._y_.shape {}, f.shape {}'.format(self.model._y_.shape, f.shape))
            nfit = self.lastfit.nfit if A20 or B1 else round(self.lastfit.nfit/n_runs)
            ndof = self.model._y_.shape[-1] - nfit # exact, all equal
            chi2_r = []
            #print('mufitplot chi len(self.lastfits),self.model._y_.shape {}, {}'.format(len(self.lastfits),self.model._y_.shape))
            for k in range(len(self.lastfits)): # either ngroups (A20) or nruns (B1) or 1 (C1)
                chi2 = (((self.model._y_[k,:]-f[k,:])/self.model._e_[k,:])**2).sum()/ndof
                chi2_r.append(chi2) # list of n_groups/n_runs exact chi2_r
        else: # if B20 or B21 or C1 or C2: # 3d yin and f, par is a list of lists 
            nfit = self.lastfit.nfit if B20 else round(self.lastfit.nfit/n_groups) if B21 else round(self.lastfit.nfit/(n_runs*n_groups))
            ndof = self.model._y_.shape[-1] - nfit
            chi2_r = []
            for krun in range(n_runs):
                #chi2 = []
                for kgroup in range(n_groups):
                    #chisquare =  (((self.model._y_[krun,kgroup,:]-f[krun,kgroup,:])/self.model._e_[krun,kgroup,:])**2).sum()/ndof if C2 else (
                    #        ((self.model._y_[krun*n_groups+kgroup,:]-f[krun*n_groups+kgroup,:])/self.model._e_[krun*n_groups+kgroup,:])**2).sum()/ndof
                    chisquare = (((self.model._y_[krun*n_groups+kgroup,:]-f[krun*n_groups+kgroup,:])/self.model._e_[krun*n_groups+kgroup,:])**2).sum()/ndof
                    chi2_r.append(chisquare)
                #chi2_r.append(chi2) 
        #print('mufitplot plot_run n_runs {}, len(chi2_r) {}'.format(n_runs,len(chi2_r)))
        return ndof, f, chi2_r
    
    def model_name(self):
        """
        returns model name (e.g. 'mlbg')

        used by plot_run
        """

        model = self.dashboard["model_guess"]
        return ''.join([component["name"] for component in model]) 
 
    def choosefftplot(self,fft_range,real):    
        """
        either amplitude or power fft of residues (refurbish!)

        distinguishes  
          single-multi run
          single-multi group
          sequential-global  
        """

        from mujpy.tools.tools import function_multi_in_components, userlocals
        from mujpy.tools.tools import get_nruns, int2min, int2min_multigroup
        from numpy import array

        # single - multi run sequential  A1 B1 both produce list of pars 
        # single - multi group sequential A1 A20 both produce list of pars 
        # multi run sequential multi group global B2 produces list of pars
        # multi group global A21 produces par
        # multi run global C1 produces par
        # multi run multi group global C2 produces pari

        # must simplify, pars are already in self.pars, asymm, asyme are obtained with 
        if self.suite.single():
            if self.suite.multi_groups(): # A2
                # print('mufitplot choosefftplot debug: A2')
                self.log('Multigroup fft animation: toggle pause/resume by clicking on the plot')
                if sum(function_multi_in_components(self.dashboard)): # A21 single chi2
                    #userpars = "userpardicts_guess" if self.guess else "userpardicts_results"
                    #pardicts = self.dashboard[userpars]
                    pars,_,_,_,_,_ = int2min_multigroup(self.dashboard,self.suite.runs,guess=self.guess)
                    # ok = self.single_fft_plot_multi_global(fft_range,pars,real)
                else: # A20 as many chi2 as groups now results are saved in single group dashboards
                    if self.jsonmodel=="model_result":
                        pars = [array(lastfit.values) for lastfit in self.lastfits]
                        # print('mufitplot choosefftplot debug: A20 fit pars = {}'.format(pars))
                    else:
                        par,_,_,_,_,_ = int2min(self.dashboard,self.suite.runs)                    
                        pars = [array(par) for k in range(len(self.suite.groups))]
                        # print('mufitplot choosefftplot debug: A20 guess pars = {}'.format(pars))
                    # ok = self.single_fft_plot_multi_sequential(fft_range,pars,real)
                asymm, asyme = self.suite.asymmetry_multigroup()
            else:  # A1 simple single plot
                pars,_,_,_,_,_ = int2min(self.dashboard,self.suite.runs)
                # ok = self.single_fft_plot(fft_range,pars,real)
                asymm, asyme = self.suite.asymmetry_single(self.suite._the_runs_[0],0)
        else: 
            if 'userpardicts_guess' in self.dashboard:
                if not self.suite.multi_groups(): # C1, userpardicts no multigroup
                    self.log('Multirun fft animation: toggle pause/resume by clicking on the plot')
                    pass# not yet
                else: 
                    self.log('Multi group and run fft animation: toggle pause/resume by clicking on the plot')
                    if userlocals(self.dashboard): # C2
                        pass# not yet
                    else: # B21 userpardicts, function_multi_in_components no tilde_incomponents
                        pars = [] 
                        # not yet
#                        for run in get_nruns(self.suite):
#                            par = 
#                            pars.append(par)
            else:
                pars = [] 
                if self.suite.multi_groups(): 
                    for run in get_nruns(self.suite): # B1 no multigroup
                        for group in self.suite.groups: # or B20 no userpardicts multigroup
                            par,_,_,_,_,_ = int2min(self.dashboard,self.suite.runs)
                            pars.append(par)
                            asymm, asyme = self.suite.asymmetry_multirun_multigroup()
                else:
                    for run in get_nruns(self.suite): # B1 no multigroup
                        par,_,_,_,_,_ = int2min(self.dashboard,self.suite.runs)
                        pars.append(par)
                    asymm, asyme = self.suite.asymmetry_multirun(0)
                # ok = self.sequential_fft_plot(fft_rangepars,real)

        # print('mufitplot choosefftplot debug: pars = {}'.format(pars))
        self.plot_fft(fft_range,pars,asymm,asyme,real)
 
    def plot_fft(self,fft_range,pars,asymm,asyme,real):
        """
        plots fft of residues 

        input:
            fft_range
                start,stop,sigma MHz
            pars list (single, calib) or list of lists (sequence)
            asymm, asyme  1d (single, calib) or 2d (sequence)
        uses data as dictated by self.dashboard["fit_range"]
        
        calls either 
            self.chi_fft (only for model function) 
            set_single_fft or set_sequence_fft, from tools.plot
            first version single only
                                    (produces actual figures
                                     using tools.plot functions)
        """

        from mujpy.tools.tools import derange, rebin, autops, ps
        from mujpy.tools.tools import get_run_title, userpars
        from mujpy.tools.plot import set_figure_fft
        from copy import deepcopy
        from numpy import ones, exp, linspace, sqrt, mean, fft
        from numpy import hstack, linspace, zeros, mgrid

        # print('plot_fft muplotfit debug: pars = {}'.format(pars))

        run_title = get_run_title(self.suite)    # always a list, even for single 
        string = 'global ' if userpars(self.dashboard) else ''  
        if len(self.suite.groups)>1:
            strgrps = [groups['forward']+'-'+groups['backward'] for groups in self.suite.groups]
        else:
            strgrps = self.suite.groups['forward']+'-'+self.suite.groups['backward']
        if self.guess:
            if len(run_title)==len(self.suite.groups):
                run_title = [runtitle + " ("+string+"guess) group" + strgrp for runtitle,strgrp in zip(run_title,strgrps)]
            else:
                run_title = [runtitle + " ("+string+"guess)" for runtitle in run_title]
        else:
            if len(run_title)==len(self.suite.groups):
                run_title = [runtitle + " ("+string+"fit) group" + strgrp for runtitle,strgrp in zip(run_title,strgrps)]
            else:
                run_title = [runtitle + " ("+string+"fit)" for runtitle in run_title]
        #############################
        # rebinning of data as in fit 
        # this works for single 
        # and for sequential
        #############################
        fittup,ermsg = derange(self.dashboard["fit_range"],self.suite.histoLength) # range as tuple
        fit_pack = 1
        if len(fittup)==3: # plot start stop pack
            fit_start, fit_stop, fit_pack = fittup[0], fittup[1], fittup[2]
        elif len(fittup)==2: # plot start stop
            fit_start, fit_stop = fittup[0], fittup[1]

        t_fit,y_fit,ey_fit = rebin(self.suite.time,asymm,[fit_start,fit_stop],fit_pack,e=asyme)

        f_fit_res, f_fit = self.chi_fft(t_fit,y_fit,*pars) # returns both partial and full model
        
        fgroup = [self.suite.groups[k]['forward'] for k in range(len(self.suite.groups))]
        bgroup = [self.suite.groups[k]['backward'] for k in range(len(self.suite.groups))]
        alpha = [self.suite.groups[k]['alpha'] for k in range(len(self.suite.groups))]
        group = [fgroup, bgroup, alpha]
        
        dt = t_fit[1]-t_fit[0]  # time bin
        fmax = 0.5/dt  # max frequancy available
        l = (fit_stop-fit_start)//fit_pack # floor division, dimension of data
        df = 1/(dt*l)
        n = 2*l # not a power of 2, but surely even
        nf = hstack((linspace(0,l,l+1,dtype=int), linspace(-l+1,-1,l-2,dtype=int)))
        dfa = 1/n/dt         # digital frequency resolution
        f = nf*dfa  # all frequencies, l+1 >=0 followed by l-1 <0

        fstart, fstop, fsigma = float(fft_range[0]), float(fft_range[1]), float(fft_range[2])
        start, stop = int(round(fstart/dfa)), int(round(fstop/dfa))
        f = deepcopy(f[start:stop]) # selected slice

# asymmetry has 1,2,3 dimensions to distinguish run from group
# here they are on same par: the data are reshaped into a 2 d array
        if len(asymm.shape)==1:
            y = zeros(n) # for data zero padded to n
            yf = zeros(n) # fit function for rephasing, zero padded to n
            ey = zeros(n)
            y[0:l] = y_fit[fit_start:fit_stop] - f_fit_res # zero padded partial residues
            yf[0:l] = f_fit
            ey[0:l] = ey_fit[fit_start:fit_stop] #  slice of time stds
            t = dt*linspace(0,n-1,n)
        elif len(asymm.shape)==2:
            # print('mufitplot plot_fft debug: nruns/ngroups {}, nbins {}'.format(asymm.shape[0],y_fit.shape[-1]))
            y = zeros((y_fit.shape[0],n)) # for data zero padded to n
            yf = zeros((y_fit.shape[0],n)) # fit function for rephasing, zero padded to n
            ey = zeros((ey_fit.shape[0],n))
            # print('mufitplot plot_fft debug: shapes y {}, asymm {}, f_fit_res {}'.format(y[:,0:l].shape,y_fit[:,fit_start:fit_stop].shape,f_fit_res.shape))
            y[:,0:l] = y_fit[:,fit_start:fit_stop]- f_fit_res # zero padded partial residues
            yf[:,0:l] = f_fit
            ey[:,0:l] = ey_fit[:,fit_start:fit_stop] #  slice of time stds
            _,t = dt*mgrid[0:asymm.shape[0],0:n]
        else: # 3 is reshaped to 2d
            nruns,ngroups = asymm.shape[0], asymm.shape[1]
            y = zeros((nruns*ngroups,n)) # for data zero padded to n
            yf = zeros((nruns*ngroups,n)) # fit function for rephasing, zero padded to n
            ey = zeros((nruns*ngroups,n))
            y[:,:,0:l] = y_fit[:,:,fit_start:fit_stop]- f_fit_res # zero padded partial residues
            yf[:,:,0:l] = f_fit # zero padded full fit function
            ey[:,:,0:l] = ey_fit[:,:,fit_start:fit_stop] #  slice of time stds
            y = y.reshape(nruns*ngroups,n)  # reshaped to 2d
            yf = yf.reshape(nruns*ngroups,n)  # reshaped to 2d
            ey = ey.reshape(nruns*ngroups,n)  # reshaped to 2d
            _,t = dt*mgrid[0:nruns*ngroups,0:n] # time for zero padded data
        # print('mufitplot plot_fft debug: t.min t.max = {},{}'.format(t.min(),t.max()))
        filter_apo = exp(-(t*fsigma)**3) # hypergaussian filter mask
                                           # is applied as if first good bin were t=0
        filter_apo = filter_apo/filter_apo.sum(axis=1)[0]/dt # approx normalization
        dfa = 1/n/dt         # digital frequency resolution
        y *= filter_apo # zero padded, filtered partial residues
        yf *= filter_apo # zero padded, filtered full fit function
        ey *= filter_apo # filteres 
        fft_e = (sqrt(mean(fft.fft(ey)**2,axis=-1))).real#.reshape(ey.shape[0],1)*ones(y.shape) # one per asymmetry
        fft_amplitude = fft.fft(y)  # amplitudes (complex), matrix with rows fft of each run
        fftf_amplitude = fft.fft(yf)  # amplitudes (complex), same for fit function
        if real:
            ########################
            # REAL PART
            # APPLY PHASE CORRECTION
            # try acme
            ########################
            if len(y_fit.shape)==1:
                fftf_amplitude[start:stop], p0, p1, out = autops(fftf_amplitude[start:stop],'acme') # fix phase on theory 
                out = out[:out.index('\n')]
                self.log('Autophase: '+out)
                fft_amplitude[start:stop] = ps(fft_amplitude[start:stop], p0=p0 , p1=p1).real 
                ap = deepcopy(fft_amplitude[start:stop].real)
                apf = deepcopy(fftf_amplitude[start:stop].real)            
            else:
                for k in range(fftf_amplitude.shape[0]):
                    fftf_amplitude[k,start:stop], p0, p1, out = autops(
                                              fftf_amplitude[k,start:stop],'acme') # fix phase on theory 
                    out = out[:out.index('\n')]
                    self.log('Autophase {}: '.format(k)+out)
                    # fft_amplitude[k,start:stop], p0, p1, out = autops(fft_amplitude[k,start:stop],'acme') # fix phase on theory 
                    fft_amplitude[k,start:stop] = ps(fft_amplitude[k,start:stop],p0=p0,p1=p1).real # same correction as theory
                ap = deepcopy(fft_amplitude[:,start:stop].real)
                apf = deepcopy(fftf_amplitude[:,start:stop].real)            
            ylabel = 'FFT Real part [arb. units]'
        else:
            ##################
            # POWER
            ##################
            if len(asymm.shape)==1:
                ap = fft_amplitude.real[start:stop]**2+fft_amplitude.imag[start:stop]**2
                apf = fftf_amplitude.real[start:stop]**2+fftf_amplitude.imag[start:stop]**2
            else:
                ap = fft_amplitude.real[:,start:stop]**2+fft_amplitude.imag[:,start:stop]**2
                apf = fftf_amplitude.real[:,start:stop]**2+fftf_amplitude.imag[:,start:stop]**2
            ylabel = 'FFT Power [arb. units]'
        self.fig_fft = set_figure_fft(self.fig_fft,self.model_name(),ylabel,
                                    f,
                                    ap,
                                    apf,
                                    fft_e,
                                    group,
                                    run_title)
#        self.the_model._include_all_() # usual.the_model.mode: all components included
        
    def chi_fft(self,t_fit,y_fit,*pars):
        """
        returns f for the fft of residues

        for the moment works only for single or sequential fits
        """

        from mujpy.tools.tools import int2_method_key, cstack
        from numpy import vstack
        the_model = self.the_model
        # can be done without load, use new self._load
        method_keys = int2_method_key(self.dashboard,self.the_model)
        self.the_model._load_(t_fit,y_fit,method_keys) # wrong must pass rangetup to self.the_model
        fft_include_components = []      
        for component in self.dashboard["model_guess"]:
            if "fft" in component.keys():
                fft_include_components.append(component["fft"]) # conditionally subtract component
            else:
                fft_include_components.append('False') # fft of data if fft is not in model dict
        self.the_model._fft_init(fft_include_components) # sets 
        # distinguish list from list of lists 
        # this applies to A20 single run multigroup sequential fits
        # print('chi_fft mufitplot debug: pars = {}'.format(pars))
        fres = self.the_model._add_fft_(t_fit,y_fit,*pars)
        if isinstance(pars,tuple):
            # print('mufitplot chi_fft debug: pars is tuple')
            f = cstack(t_fit,*pars)
        else: 
            # print('mufitplot chi_fft debug: pars is not tuple')
            self.the_model._add_(t_fit,*pars)
        # debug
#        from matplotlib.pyplot import draw, subplots 
#        fig,ax = subplots()
#        for k in range(f.shape[0]):
#            ax.plot(t_fit,f[k,:])
#            ax.plot(t_fit,fres[k,:],alpha=0.5)
#        draw()    
        return fres, f

##############
# Deprecated
##############
                              
    def reload(self,tup):
        """
        Deprecated
        reloads data through mufit.mumodel and mufit.suite

            utilizing self.pars, self.dashboard, self.suite
        into self.x, self.y, self.e 
        that overloads mucomponents 
            self._x_, self._y_, self._e_, required by the cost function
        """

        from mujpy.tools.tools import set_alpha, rebin, rshp
        if self.calib: 
            #self.fit._alpha = set_alpha(self.pars,len(self.suite.runs),len(self.suite.groups))
            self.model._reload_calib_(tup)
        else:
            self.model._reload_(tup) # 
        self.x, self.y, self.e = rshp(self.model._x_), reshaplot(self.model._y_), reshaplot(self.model._e_)
 
#    def single_chi(self):
#        """
#        unused
#        output:
#            True if chi_1 is required (single cost function in global mufitplot: A1,A21,C1,C2 adn their calib)
#            False if chi_2 is required (multi cost finctions in sequential mufitplots: A20,B1,B20,B21 and their calib)
#        """
#        #       multi_groups suite.single userpars function_multi_in_components userlocals  single_cost_function
#        # A1     False          True        False       False               False        *
#        # A20    True           True        False       False               False
#        # A21    True           True        True        True                False        *
#        # B1     False          False       False       False               False
#        # B20    True           False       False       False               False
#        # B21    True           False       True        True                False
#        # C1     False          False       True        False               True         *
#        # C2     True           False       True        True                True         *
#        #  True : A1 1d, A21 2d multigroup userpar, C1 2d userpar, C2 3d multgrup userpar
#        #  False: B1 2d, A20, B20 2d multigroup, B21 2d multigroup userpar
##        from mujpy.tools.tools import userlocals, userpars
##        return (userlocals(self.dashboard) or 
##                   (self.suite.single() and 
##                        (not self.suite.multi_groups() or userpars(self.dashboard)))) 
#        return self.fit.A1() or self.fit.A21() or self.fit.A1_calib() or self.fit.A21_calib() or self.fit.C1() or self.fit.C1_calib() or self.fit.C2() or self.fit.C2_calib()
 
#    def chooseplot(self,plot_range):
#        """
#        deprecated
#            switch for single (A1), sequential (B1),
#            ... 
#        #       multi_groups suite.single userpars function_multi_in_components userlocals
#        # A1     False          True        False       False               False *
#        # A20    True           True        False       False               False
#        # A21    True           True        True        True                False *
#        # B1     False          False       False       False               False
#        # B20    True           False       False       False               False does not exist
#        # B21    True           False       True        True                False
#        # C1     False          False       True        False               True  *
#        # C2     True           False       True        True                True  *
#        """
#        from mujpy.tools.tools import function_multi_in_components, userpars
#
#        A1, A20, A21, B1, B20, B21, C1, C2 = self.fit_types()
#        if A1:
#            ok, msg = self.plot_singlerun(plot_range)
#        elif A20 or A21:
#            self.log('Multigroup fit animation: toggle pause/resume by clicking on the plot')
#            ok, msg = self.plot_singlerun_multigroup(plot_range)
#        else:
#            self.log('Multirun fit animation: toggle pause/resume by clicking on the plot')
#            if B1():
#                ok, msg = self.plot_multirun_singlegroup_sequential(plot_range)
#            elif B21():
#                ok, msg = self.plot_multirun_multigroup_userpar(plot_range)
#            elif C1():
#                self.log('No C1 yet. Exiting mufitplot without a plot')
#                ok, msg = False, 'C1'
#            else: #  C2()
#                ok, msg = self.plot_multirun_multigroup_userpar(plot_range)
#       
#        if not ok:
#            self.log('Exiting mufitplot without a plot: '+msg)
#    
#    def plot_singlerun(self,plot_range):
#        """
#        input plot_range a csv string, start, stop [,pack [start_r, stop_r]]
#              obtains asymm, asyme = self.the_model._asymmetry_(alpha)
#        A1
#        """
#        from mujpy.tools.tools import int2min, mixer, calib
#        from numpy import cos, pi
##        self.log('muplotfit plot_singlerun debug')
#        kgroup = 0 # default single group
#        pars = self.lastfit.values
#        #pars,_,_,_,_,_ = int2min(self.dashboard,self.suite.runs)
#        if calib(self.dashboard): # A1 calib
#            self.suite.grouping[kgroup]['alpha'] = pars[0] 
#            self.suite.groups[kgroup]['alpha'] = pars[0]
#            asymm, asyme = self.the_model._asymmetry_(pars[0])
#        else: # A1
#            asymm, asyme = self.suite.asymmetry_single(self.suite._the_runs_[0],0)
#        if self.rotating_frame_frequencyMHz:
#            self.rrf_asymm = mixer(self.suite.time,asymm,self.rotating_frame_frequencyMHz)
#            self.rrf_asyme = mixer(self.suite.time,asyme,self.rotating_frame_frequencyMHz)
#        #self.log('mufitplot: Inside single plot; debug mode')
#        return self.plot_run(plot_range,pars,asymm,asyme)        
#
#    def plot_singlerun_multigroup(self,plot_range):
#        """
#        remove!
#        method chosen to standardize, all done
#        input plot_range, passed to self.plot_run together with
#            pars, list of one list of fit parameter values
#            asymm, asyme 2d # A20, A20_calib, A21, A21 calib
#        calib retrieves alpha as [par[0],par[1],...] and
#              obtains asymm, asyme = self.the_model._asymmetry_(alpha)
#        A20 A21
#        """
#        from mujpy.tools.tools import mixer, calib, int2min_multigroup, function_multi_in_components, int2_method_key, int2_multigroup_method_key, set_alpha, rebin
#        from numpy import vstack, array
#        
        

#    def plot_multirun_multigroup_userpar(self,plot_range):
#        """
#         inputs: 
#            plot_range
#        # C2 fit 
#        # asymm, asyme are 3d
#        """
#        from mujpy.tools.tools import int2min_multirun_multigroup, mixer, calib
#        from numpy import cos, pi, vstack
##        userpardicts = (self.dashboard["userpardicts_guess"] if self.guess else 
##                        self.dashboard["userpardicts_result"])
##        pardict = self.dashboard["model_guess"][0]["pardicts"][0]
#        pars,_,_,_,_,_ = int2min_multirun_multigroup(self.dashboard,self.suite.runs)
## no calib mode!
##        p = pars
##        if calib(self.dashboard):
##            for kgroup,group in enumerate(self.suite.groups):
##                group['alpha'] = eval(pardict["function_multi"][kgroup])
##                self.suite.groups[kgroup]["alpha"] = eval(pardict["function_multi"][kgroup])
#        asymm, asyme = self.suite.asymmetry_multirun_multigroup()
#        if self.rotating_frame_frequencyMHz:
#            for krun in range(asymm.shapcalibe[0]):
#                for kgroup in range(asymm.shape[1]):
#                    if not kgroup: # kgroup=0
#                        tim = self.suite.time
#                    else:
#                        tim = vstack((time,self.suite.time))
#                if not krun: # krun=0
#                    time = array([tim])
#                else:
#                    time = vstack(time,array([tim]))                      
#            self.rrf_asymm = mixer(time,asymm,self.rotating_frame_frequencyMHz)
#            self.rrf_asyme = mixer(time,asyme,self.rotating_frame_frequencyMHz)
#        #self.log('mufitplot: Inside single plot; debug mode')
#        return self.plot_run(plot_range,pars,asymm,asyme)        
#
#
#    def plot_multirun_singlegroup_sequential(self, plot_range):
#        """
#        input plot_range, passed to self.plot_run together with
#            pars, list of lists  of fit parameter values
#            asymm, asyme 2d
#            # B1 fit
#        """
#        from mujpy.tools.tools import mixer
#        from numpy import vstack
#        
#        kgroup = 0 # default single 
#        # dashboard must become a suite thing: each sequential fit has its own
#        pars = []
#        # self.log('mufitplot: Inside sequential plot; debug mode')
#        for k,lastfit in enumerate(self.fit.lastfits): # loops over runs
#            values = lastfit.values   # list of values for a run
#            
#            #print('plot_run muplotfit debug: pars = {}'.format(pars)) run {},  values = {}'.format(run[0].get_runNumber_int(),values))
#            pars.append(values) # list of lists one per run
#        asymm, asyme = self.suite.asymmetry_multirun(kgroup) # 
#        if self.rotating_frame_frequencyMHz:
#            for k in range(asymm.shape[0]):
#                if not k:
#                    time = self.suite.time
#                else:
#                    time = vstack((time,self.suite.time))
#            self.rrf_asymm = mixer(time,asymm,self.rotating_frame_frequencyMHz)
#            self.rrf_asyme = mixer(time,asyme,self.rotating_frame_frequencyMHz)
#        return self.plot_run(plot_range,pars,asymm,asyme)        
#  
#    def plot_multirun_singlegroup_user(self,plot_range):  
#        """
#        input plot_range, passed to self.plot_run together with
#            pars, list of lists of fit parameter values
#            asymm, asyme 2d
#            # C1 fit (presented to plot_run like A21 singlerun_multigroup_userpar)
#            the difference is that A21 all Minuit pars are userpardicts, 
#            while here there is a mixed situation
#        """
#        from mujpy.tools.tools import mixer
#        from numpy import vstack
#        
#        kgroup = 0 # default
#        # dashboard must become a suite thing: each sequential fit has its own
#        pars = list(self.fit.lastfit.values)  # self.lastfit.values is the damn ValueView
#        # self.log('mufitplot: Inside sequential plot; debug mode')
#        # special arrangement for C1 fits plotted here as B1 fits
#        # self.fit.lastfits[0] is guess, self.fit.lastfits[1] is result
#        # both are Minuit instances with data loaded as multirun sequential
#        asymm, asyme = self.suite.asymmetry_multirun(kgroup) # 
#        if self.rotating_frame_frequencyMHz:
#            for k in range(asymm.shape[0]):
#                if not k:
#                    time = self.suite.time
#                else:
#                    time = vstack((time,self.suite.time))
#            self.rrf_asymm = mixer(time,asymm,self.rotating_frame_frequencyMHz)
#            self.rrf_asyme = mixer(time,asyme,self.rotating_frame_frequencyMHz)
#        return self.plot_run(plot_range,pars,asymm,asyme)        
#    def chi_1(self,t,yin,eyin,pars):
#        """
#        deprecated
#        input:
#            t, yin, eyin are time, asymmetry, error
#            pars fit parameter values (single run or calib) 
#               (fitvalues obtained by  appropriate int2min)
#            kgroup index of grouping assumed 0
#        output:
#            number_dof, fit function, chi2_r
#            chi2_r is the [approx] chi_2 of each 1d asymmetry
#        valid for any fit that produces a single chi: A1, C1, C2
#                                                 and their calib
#        """
#        #from mujpy.tools.tools import function_multi_in_components, int2_multigroup_method_key, int2_multirun_user_method_key
#        
#        if yin is None:
#            ndof = None
#            chi2 = None
#        else: # yin and data are not None 
#              # lastfit.nfit is n. of Minuit parameters
#
##            if sum(function_multi_in_components(self.dashboard)):
##                if calib(self.dashboard):
##                # must be handles as a standard global multigroup plot
##                    methods_keys = methods_keys[1:]
##                    # ok, msg = self.the_model._load_calib_multigroup_(t,yf,yb,eyf,eyb,returntup,methods_keys) # use with full methods_keys
##                    # rebinning is done in mucomponents _add_calib_multigroup_
##                #else:
##                ok, msg = self.the_model._load_multigroup_(t,
##                                                                  yin,
##                                                                  methods_keys,
##                                                                  e=eyin)
##            elif self.fit.C1():
##                ok, msg = self.the_model._load_multirun_user_(t,
##                                                                      yin,
##                                                                      methods_keys,
##                                                                      e=eyin)
##            elif calib(self.dashboard):
##                methods_keys = methods_keys[1:]
##                ok, msg = self.the_model._load_calib_(t,
##                                                             yin,
##                                                             methods_keys,
##                                                             e=eyin) 
##            else:
##                ok, msg = self.the_model._load_(t,
##                                                       yin,
##                                                       methods_keys,
##                                                       e=eyin)             
##            if ok:
##                # print('mufitplot chi_1 debug: pars = {}'.format(pars))
##                # self.debug(pars) # works only for global fits
##                f = self.the_model._add_(t,*pars)
##                chi2 = self.the_model._chisquare_(*pars)/ndof # chi2 in plot
##                if sum(function_multi_in_components(self.dashboard)) or self.fit.C1():
##                    # global cases: chi2 is computed on the fly 
##                    # over true plot data, loaded above by the appropriate _load_...
##                    self.the_model._axis_ = 1 # this produces separate fcn values per run/grouping
##                    ndof = int(ndof/yin.shape[0])  # number of dof per run or grouping = total number of dof / number of runs or groupings 
##                    chi2 = self.the_model._chisquare_(*pars)/ndof # chi2 in plot
##                    self.the_model._axis_ = 1
##                    # print('debug mufitplot chi_1: ndof = {}'.format(ndof))
##                    # for k in range(yin.shape[0]):
##                    #     print('* debug chi2({}) = {}'.format(k,chi2[k]))
##            else: 
##                self.log(msg)
##                return None, None, None            
#        return ndof,f,chi2
#        
#    def chi_2(self,t,yin,eyin,pars):
#        """
#        input:
#            t, yin, eyin are time 1d, asymmetry, error 2 or 3d
#            pars list of lists [of lists]  of fit parameter values for sequence 
#               (fitvalues as obtained by int2min)
#            nruns, ngroups numbers to distinguish 2d and 3d
#            otherwise sequential or global groups 
#        output:
#            number_dof, fit function (2d), chi2_r (list)
#        deals also with cases where only f is needed (plot function)
#        #  B1 2d (runs) sequential, pars is a list, uses _load_ (2d yin nruns lastfits) 
#        #  A20 2d (groups) sequential, pars is a list, uses _load_ (2d yin ngroup lastfits)
#        # to be implementedint2_multigroup_method_key
#        #  B20 3d (group,runs) sequential, pars is a list, uses _load_multigroup_ (3d yin nruns x ngroups lastfits)
#        #  B21 3d multigroup userpar, to be done (3d yin nruns lastfits) only cas in which len(yin.shape)>depth(lastfits)
#        """
#        from mujpy.tools.tools import _nparam, int2_method_key, calib, cstack, ccstack
# 
#        depth = lambda L: isinstance(L, list) and max(map(depth, L))+1 
#        nruns, ngroups = len(self.suite.runs), len(self.suite.groups)
#        stack = ccstack if nruns > 1 and ngroups > 1  else cstack
#        f = stack(self.the_model._add_norebin_,t,nruns,ngroups,*pars) # stack vstacks the calculated function to match vstacked asymmetry yin 
#        chi2 = []
#        ndof = []
#        if yin is not None: # data are not None
#            
#            if depth(self.fit.lastfits)== 3: # B20, not used yet
#                ndof = yin.shape[2] - self.fit.lastfits[0][0].nfit 
#                for krun in range(len(self.fit.lastfits)):
#                    chi = []
#                    for kgroup in range(len(self.fit.lastfits[krun])):
#                        chi.append((((yin[krun,kgroup,:]-f[krun,kgroup,:])/eyin[krun,kgroup,:])**2).sum()/ndof)
#                    chi2.append(chi)
#            elif len(yin.shape)==3: # B21, lastfits is a list of lists, not used yet
#                ndof = yin.shape[2] - self.fit.lastfits[0][0].nfit 
#                for krun in range(len(self.fit.lastfits)):
#                    chi2.append((((yin[krun,:,:]-f[krun,:,:])/eyin[krun,:,:])**2).sum()/ndof)
#            else: # B1, A20, 2d data, only one in use
#                ndof = yin.shape[1] - self.fit.lastfits[0].nfit 
#                for k in range(len(self.lastfits)):
#                    chi2.append((((yin[k,:]-f[k,:])/eyin[k,:])**2).sum()/ndof)
#         
##        _,_, freepars = _nparam(self.dashboard[self.model])
##        ndof = len(t) - freepars # degrees of freedom in plot
###        print('debug mufitplot chi_2: freepars = {} nu = {}'.format(freepars,nu))
##        if yin is not None: # data are not None
##            mthdk = int2_method_key(self.dashboard,self.the_model.
##            if calib(self.dashboard):
##                mthdk = mthdk[1:]
##            ok, msg = self.the_model._load_(t,yin,mthdk,e=eyin) 
##            if ok:
##                f = cstack(self.the_model._add_,t,*pars)
##                # print('mufitplot chi_2 debug: pars {}\nf.shape = {}'.format(pars,f.shape))
##            else:
##                self.log(msg)
##                return None, None, None
##            chi2 = [self.lastfits[k].fval/ndof for k in range(len(self.lastfits))]
###            for k in range(len(lastfits)): # works only for 2d cases
###                chi2.append(lastfits[k]) # chi2 in plot
##                # print('mufitplot chi_2 debug: k {}, chi2 = {}'.format(k,chi2[-1]))
##        else: # only f is really needed, this can be called only after a call with yin = None
##            f = cstack(self.the_model._add_,t,*pars)
##            chi2 = [None for k in range(len(pars))]
#        return ndof,f,chi2
                
#    def fstack(self,t,*pars):
#        """
#        deprecate
#        done by model._add_
#        """
#        from numpy import vstack
#        # print('fstack mufitplot debug: pars = {}'.format(pars))
#        for k,par in enumerate(pars):
#            # print('fstack mufitplot debug: par = {}'.format(par))
#            fin = self.the_model._add_(t,*par)
#            if k==0: # f for histogram
#               f = fin
#            else:
#               f = vstack((f,fin))
#        return f

                   
    def debug(self,pars):
        # can be deleted after debugging stage of mujpy 2.0
        from mujpy.tools.tools import min2int_multigroup
        names, ps, _ = min2int_multigroup(self.dashboard,pars,pars)
        print('muplotfit debug: internal dashboard parameters {}'.format([[name,p] for name,p in  zip(names,ps)]))
