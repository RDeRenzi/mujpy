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

###
# fitplot
###    

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

def set_fig(num,nrow,ncol,title,**kwargs):  # NOT CLEAR WHERE IT IS USED
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

def set_single_fit(fig,model,early_late,data,group,run_title,chi_dof,data_late,chi_dof_late,rrf=0):
    '''
    input:
        
        flags = [early_late, anim] True/False
        data = [t,y,ey,f_res,tf,f,dy_fit] 
            used to errorbar(t,y,yerr=ey), plot(tf,f), plot(t,y-f_res), chi(dy_fit)
        if early_late:
            data_late = [t_late,y_late,ey_late,f_res_late,tfl,fl,dy_fit_early,dy_fit_late]
        group = [fgroup,bgroup,,alpha]
        run_title, to print on figure
        chi_dof = [nu,chi], number of dof and chi2
                chi are scalar in single
        if early_late:
            chi_dof_late = [nu_e,nu_l,chi_e,chi_l,w_e,w_l], 
            w_e, w_l = nu_fit/nu_e\l*ones(dy_fit_early\late.shape[0])
    output: 
        fig, to reuse fig window
    recovers or creates figure subplots 
    for single fit or anim
    '''
    import matplotlib.pyplot as P
    from matplotlib.pyplot import rcParams, suptitle
    from numpy import hstack   
    from mujpy.tools.plot import errorb, plot_fit, decorate_data, decorate_data_late
    from mujpy.tools.plot import plot_res, decorate_res, decorate_res_late, plot_txt, plot_chi2
    from mujpy.tools.plot import draw
        
    font = {'family':'Ubuntu','size':10}
    P.rc('font', **font)
    prop_cycle = rcParams['axes.prop_cycle']
    color = prop_cycle.by_key()['color']
    # data and residues, data_early and residue early: color[0] 
    # fit, both early and late: color[1]
    # data_late and residue_late: color [2]

    t, y, ey, f_res, tf, f, dy_fit = data
    fgroup, bgroup, alpha = group
    nu_fit, chi_fit = chi_dof  
      
    if early_late:
        ncols, width_ratios = 3,[2,2,1]
    else:
        ncols, width_ratios = 2,[4,1]

    try:    
        fig.clf()
        fig,ax = P.subplots(2,ncols,sharex = 'col', 
                     gridspec_kw = {'height_ratios':[3,1],'width_ratios':width_ratios},
                                    num=fig.number)
        fig.subplots_adjust(hspace=0.05,top=0.90,bottom=0.12,right=0.97,wspace=0.03)
    except: # handle does not exist, make one
        fig,ax = P.subplots(2,ncols,figsize=(6,4),sharex = 'col',
                     gridspec_kw = {'height_ratios':[3, 1],'width_ratios':width_ratios})
        fig.canvas.manager.set_window_title('Fit')
        fig.subplots_adjust(hspace=0.05,top=0.90,bottom=0.12,right=0.97,wspace=0.03)

    
    if early_late: 
        t_late, y_late, ey_late, f_late_res, tfl, fl, dy_fit_early, dy_fit_late = data_late
        nu_fit_early,nu_fit_late, chi_fit_early,chi_fit_late,w_early,w_late = chi_dof_late
        [[ax_early,ax_late,ax_txt],[ax_early_res,ax_late_res,ax_chi]] = ax
        ytot, rtot = hstack((y,y_late)), hstack((y-f_res,y_late-f_late_res))
        ym,yM,rm,rM = ytot.min()-0.05,ytot.max()+0.01,rtot.min()-0.005,rtot.max()+0.005
        _,_,_ = errorb(ax_early,t,y,ey,color[0])
        _ = plot_fit(ax_early,tf,f,color[1])
        _ = plot_res(ax_early_res,t,y-f_res,color[0])
        _,_,_ = errorb(ax_late,t_late,y_late,ey_late,color[2])
        _ = plot_fit(ax_late,tfl,fl,color[1])
        _ = plot_res(ax_late_res,t_late,y_late-f_late_res,color[2])
        decorate_data(ax_early,t,ym,yM)
        _,_,_,_ = decorate_res(ax_early_res,t,ey,rm,rM)
        decorate_data_late(ax_late,t_late,ym,yM)
        _,_,_,_ = decorate_res_late(ax_late_res,t_late,ey_late,rm,rM)
    else:
        [[ax_fit,ax_txt],[ax_res,ax_chi]] = ax 
        ym,yM,rm,rM = y.min()-0.05,y.max()+0.01,(y-f_res).min()-0.005,(y-f_res).max()+0.005
        nu_fit_early,nu_fit_late,dy_fit_early,w_early = None, None, None, None
        dy_fit_late,w_late,chi_fit_early,chi_fit_late = None, None, None, None
        _,_,_ = errorb(ax_fit,t,y,ey,color[0])
        _ = plot_fit(ax_fit,tf,f,color[1])
        _ = plot_res(ax_res,t,y-f_res,color[2])
        decorate_data(ax_fit,t,ym,yM)
        _,_,_,_ = decorate_res(ax_res,t,ey,rm,rM)
        
    # leftmost text (top) and error distribution plot (bottom) - chi are scalar in single
    #(ax,nu_fit,nu_early,nu_late,chi_fit,chi_early,chi_late,fgroup,bgroup,alpha,ylim)
    _,_,_,_,_,_,_,_,_ = plot_txt(ax_txt,model,
                                 nu_fit,
                                 nu_fit_early,
                                 nu_fit_late,
                                 chi_fit,
                                 chi_fit_early,
                                 chi_fit_late,
                    fgroup[0],bgroup[0],alpha[0],[-5,ym])
    # (ax,dy_fit,nu_fit,dy_early,w_early,dy_late,w_late)
    _,_,_,_,_,_,_,_ = plot_chi2(ax_chi,
                                      dy_fit,
                                      nu_fit,
                                      dy_fit_early,
                                      w_early,
                                      dy_fit_late,
                                      w_late)   
    if rrf:
        string = r'$\nu_R=$'+'{:.1f}MHz   '.format(rrf)
        suptitle(string+run_title,x=-0.125,ha='right')    
    else:
        suptitle(run_title)                           
    draw(fig)
    return fig

def set_sequence_fit(fig,model,early_late,data,group,run_title,chi_dof,data_late,chi_dof_late,rrf=0):
    '''
    input:
        
        flags = [early_late, anim] True/False
        data = [t,y,ey,f_res,tf,f,dy_fit] 
            used to errorbar(t,y,yerr=ey), plot(tf,f), plot(t,y-f_res), chi(dy_fit)
        if early_late:
            data_late = [t_late,y_late,ey_late,f_res_late,tfl,fl,dy_fit_early,dy_fit_late]
        group = [fgroup,bgroup,,alpha]
        run_title, to print on figure
        chi_dof = [nu,chi], number of dof and chi2
            chi are lists in sequence
        if early_late:
            chi_dof_late = [nu_e,nu_l,chi_e,chi_l,w_e,w_l], 
            w_e, w_l = nu_fit/nu_e\l*ones(dy_fit_early\late.shape[0])
    output: 
        fig, to reuse fig window
    recovers or creates figure subplots 
    for single fit or anim
    '''
    import matplotlib.pyplot as P
    from matplotlib.pyplot import rcParams

    import matplotlib.animation as animation
    from numpy import hstack, vstack   
    from mujpy.tools.plot import errorb, plot_fit, decorate_data, decorate_data_late
    from mujpy.tools.plot import plot_res, decorate_res, decorate_res_late, plot_txt, plot_chi2
    from mujpy.tools.plot import draw
    from datetime import datetime

    def animate_fit(i): 
        '''
        anim function
        update errorbar data, fit, residues and their color,
               chisquares, their histograms 

        '''
        from  numpy import histogram, array
        # early_late == False if data_late == None
        #early_late = False if data_late == None else False

        global ax_0
        
        #print('animate_fit: plot debug: Hey, I am here!, i = {}'.format(i))
        if rrf:
            supti.set_text(stringrrf+run_title[i]) 
        else:
            supti.set_text(run_title[i]) 
        line.set_ydata(y[i]) # begin errorbar
        segs = [array([[q,w-a],[q,w+a]]) for q,w,a in zip(t,y[i],ey[i])]
        ye[0].set_segments(segs) 

        # update fit plot 
        fline.set_ydata(f[i]) # fit

        # update residues
        res.set_ydata(y[i]-f_res[i]) # residues
        linesp.set_ydata(ey[i])        
        linesm.set_ydata(-ey[i])        
        line2sp.set_ydata(2*ey[i])        
        line2sm.set_ydata(-2*ey[i])        

        # update chi2 histogram, entire fit   
        nhistf,_ = histogram(dy_fit[i],xbin)
        top = bottom + nhistf
        vertf[1::5, 1] = top
        vertf[2::5, 1] = top

        # update text
        if len(fgroup)==1:
            string1 = 'F-B: {} - {}\n'.format(fgroup[0],bgroup[0])
            string1 += r'$\alpha=$ {:.4f}'.format(alpha[0])
        else:
            string1 = 'F-B: {} - {}\n'.format(fgroup[i],bgroup[i])
            string1 += r'$\alpha=$ {:.4f}'.format(alpha[i])
        string1 += '\n$\chi^2_f=$ {:.4f}\n ({:.2f}-{:.2f})\n{} dof'.format(chi_fit[i],lc,hc,nu_fit) 
        text1.set_text(string1)

        if early_late:
        # update late replica
            # begin errorbar
            linel.set_ydata(y_late[i]) 
            segs = [array([[q,w-a],[q,w+a]]) for q,w,a in zip(t_late,y_late[i],ey_late[i])]
            yel[0].set_segments(segs) 

            # fit
            flinel.set_ydata(fl[i]) 
            
            # residues
            resl.set_ydata(y_late[i]-f_late_res[i]) 
            # std, 2std lines
            linespl.set_ydata(ey_late[i])        
            linesml.set_ydata(-ey_late[i])        
            line2spl.set_ydata(2*ey_late[i])        
            line2sml.set_ydata(-2*ey_late[i])        

            # chi2 histograms
            nhiste,_ = histogram(dy_fit_early[i],xbin,weights=w_early)
            top = bottom + nhiste
            verte[1::5, 1] = top
            verte[2::5, 1] = top
            nhistl,_ = histogram(dy_fit_late[i],xbin,weights=w_late)
            top = bottom + nhistl
            vertl[1::5, 1] = top
            vertl[2::5, 1] = top
            
            #nufitplot.set_ydata(nufit[i]*yh)
            string2 = '$\chi^2_e=$ {:.4f}\n({:.2f}-{:.2f})'.format(chi_fit_early[i],lce,hce)
            string3 = '$\chi^2_l=$ {:.4f}\n({:.2f}-{:.2f})'.format(chi_fit_late[i],lcl,hcl)
            text2.set_text(string2)
            text3.set_text(string3)

            return (line, ye, fline, res, linel, yel, flinel, resl, linesp, linesm, line2sp,       
                line2sm, linespl, linesml, line2spl, line2sml, vertf, vertl, verte, text1,text2, text3)
        return line, ye, fline, res, linesp, linesm, line2sp, line2sm, vertf, text1
              
    def init_animate_fit():
        '''
        anim init function
        blitting (see wikipedia)
        to give a clean slate 

        '''
        from  numpy import histogram, array
        global ax_0

        # print('init_animate_fit: plot debug: Hey, I am here!')
        if rrf:
            supti.set_text(stringrrf+run_title[0]) 
        else:
            supti.set_text(run_title[0]) 
        line.set_ydata(y[0]) # begin errorbar
        segs = [array([[q,w-a],[q,w+a]]) for q,w,a in zip(t,y[0],ey[0])]
        ye[0].set_segments(segs)
        fline.set_ydata(f[0]) # fit
        res.set_ydata(y[0]-f_res[0]) # residues
        linesp.set_ydata(ey[0])        
        linesm.set_ydata(-ey[0])        
        line2sp.set_ydata(2*ey[0])        
        line2sm.set_ydata(-2*ey[0])        

        # update chi2 histogram, entire fit   
        nhistf,_ = histogram(dy_fit[0],xbin)
        top = bottom + nhistf
        vertf[1::5, 1] = top
        vertf[2::5, 1] = top

        # update text
        string1 = 'F-B: {} - {}\n'.format(fgroup[0],bgroup[0])
        string1 += r'$\alpha=$ {:.4f}'.format(alpha[0])
        string1 += '\n$\chi^2_f=$ {:.4f}\n ({:.2f}-{:.2f})\n{} dof'.format(chi_fit[0],lc,hc,nu_fit) 
        text1.set_text(string1)

        if early_late:
        # update late replica
            # begin errorbar
            linel.set_ydata(y_late[0]) 
            segs = [array([[q,w-a],[q,w+a]]) for q,w,a in zip(t_late,y_late[0],ey_late[0])]
            yel[0].set_segments(segs) 

            # fit
            flinel.set_ydata(fl[0]) 
            
            # residues
            resl.set_ydata(y_late[0]-f_late_res[0]) 
            # std, 2std lines

            # std, 2std lines
            linespl.set_ydata(ey_late[0])        
            linesml.set_ydata(-ey_late[0])        
            line2spl.set_ydata(2*ey_late[0])        
            line2sml.set_ydata(-2*ey_late[0])        
            # chi2 histograms
            nhiste,_ = histogram(dy_fit_early[0],xbin,weights=w_early)
            top = bottom + nhiste
            verte[1::5, 1] = top
            verte[2::5, 1] = top
            nhistl,_ = histogram(dy_fit_late[0],xbin,weights=w_late)
            top = bottom + nhistl
            vertl[1::5, 1] = top
            vertl[2::5, 1] = top
            
            #nufitplot.set_ydata(nufit[i]*yh)
            string2 = '$\chi^2_e=$ {:.4f}\n({:.2f}-{:.2f})'.format(chi_fit_early[0],lce,hce)
            string3 = '$\chi^2_l=$ {:.4f}\n({:.2f}-{:.2f})'.format(chi_fit_late[0],lcl,hcl)
            text2.set_text(string2)
            text3.set_text(string3)
            return (line, ye, fline, res, linel, yel, flinel, resl, linesp, linesm, line2sp,
            line2sm, linespl, linesml, line2spl, line2sml, vertf, vertl, verte, text1,text2, text3)
        return line, ye, fline, res, linesp, linesm, line2sp, line2sm, vertf, text1


    global paused
    global anim_fit,ax_0

    paused = False
    def toggle_pause(*args, **kwargs):
        global paused
        global anim_fit
        if paused:
            anim_fit.event_source.start() # matplotlib.__version__ >= 3.4 
            # animation.event_source.start()
        else:
            anim_fit.event_source.stop() # if matplotlib.__version__ >= 3.4 
            # animation.event_source.stop()
        paused = not paused    

## set_sequence_fit begins here
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")             
    font = {'family':'Ubuntu','size':10}
    P.rc('font', **font)
    prop_cycle = rcParams['axes.prop_cycle']
    color = prop_cycle.by_key()['color']
    anim_delay = 1000.
#        # unpack data, chi_dof
#        # data packed by mufitplot: 
#        # chi_dof packed by mufitplot: 

    t, y, ey, f_res, tf, f, dy_fit = data
#    print('tools plot debug: max, min t {}, y {}, ey {}, f_res {}, tf {}, f {}, dy_fit {}'.format(
#                                                                                 [t.min(),t.max()],
#                                                                                 [y.min(),y.max()],
#                                                                               [ey.min(),ey.max()],
#                                                                         [f_res.min(),f_res.max()],
#                                                                               [tf.min(),tf.max()],
#                                                                                 [f.min(),t.max()],
#                                                                      [dy_fit.min(),dy_fit.max()]))
    fgroup, bgroup, alpha = group
    nu_fit, chi_fit = chi_dof  
    n = len(chi_fit)
      
    if early_late:
        ncols, width_ratios = 3,[2,2,1]
    else:
        ncols, width_ratios = 2,[4,1]

    try:    
        fig.clf()
        fig,ax = P.subplots(2,ncols,sharex = 'col', 
                     gridspec_kw = {'height_ratios':[3,1],'width_ratios':width_ratios},
                                    num=fig.number)
        fig.subplots_adjust(hspace=0.05,top=0.90,bottom=0.12,right=0.97,wspace=0.03)
    except: # handle does not exist, make one
        fig,ax = P.subplots(2,ncols,figsize=(6,4),sharex = 'col',
                     gridspec_kw = {'height_ratios':[3, 1],'width_ratios':width_ratios})
        fig.canvas.manager.set_window_title('Fit')
        fig.subplots_adjust(hspace=0.05,top=0.90,bottom=0.12,right=0.97,wspace=0.03)

    if early_late:   
        t_late, y_late, ey_late, f_late_res, tfl, fl, dy_fit_early, dy_fit_late = data_late
        nu_fit_early,nu_fit_late, chi_fit_early,chi_fit_late,w_early,w_late = chi_dof_late
        [[ax_early,ax_late,ax_txt],[ax_early_res,ax_late_res,ax_chi]] = ax
        #print('debug-plot: shape dy_fit_late {}, w_late {}'.format(dy_fit_late.shape, w_late.shape))

        for k in range(n):  # ytot rtot for gettimg absolute max min of data and residue plots
            yt, rt = hstack((y[k],y_late[k])), hstack((y[k]-f_res[k],y_late[k]-f_late_res[k]))
            if k==0:
                ytot, rtot = yt, rt
            else:
                ytot, rtot = vstack((ytot,yt)), vstack((rtot,rt))            
        ym,yM,rm,rM = ytot.min()-0.05,ytot.max()+0.01,rtot.min()-0.005,rtot.max()+0.005
        xtgl,ytgl = 1.35*t.min()-0.35*t.max(), 1.65*rm-0.65*rM
        ax_early_res.text(xtgl,ytgl,'Click to toggle\npause/resume',fontsize='small')
        ax_early.text(xtgl,yM+0.01,dt_string,fontsize='small')
 
 
        line, xe, ye = errorb(ax_early,t,y[0],ey[0],color[0])
        fline = plot_fit(ax_early,tf,f[0],color[1])
        res = plot_res(ax_early_res,t,y[0]-f_res[0],color[0])
        linel, xel, yel = errorb(ax_late,t_late,y_late[0],ey_late[0],color[2])
        flinel = plot_fit(ax_late,tfl,fl[0],color[1])
        resl = plot_res(ax_late_res,t_late,y_late[0]-f_late_res[0],color[2])
        decorate_data(ax_early,t,ym,yM)
        linesp,linesm,line2sp,line2sm = decorate_res(ax_early_res,t,ey[0],rm,rM)
        decorate_data_late(ax_late,t_late,ym,yM)
        linespl,linesml,line2spl,line2sml = decorate_res_late(ax_late_res,t_late,ey_late[0],rm,rM)
        
        # if multigroup_in_components (here unknown) nufit=0.5*nufit, find a hack
        # nufit 

        vertf, nhistf, bottom, xbin, vertl, nhistl, verte, nhiste = plot_chi2(
                                                            ax_chi,dy_fit[0],nu_fit,
                                                            dy_fit_early[0],w_early,
                                                            dy_fit_late[0],w_late)
            #(ax,nu_fit,nu_early,nu_late,chi_fit,chi_early,chi_late,fgroup,bgroup,alpha,ylim)
        text1, lc, hc, text2, lce, hce, text3, lcl, hcl  = plot_txt(
                                                                    ax_txt,model,
                                                                    nu_fit,
                                                                    nu_fit_early,
                                                                    nu_fit_late,
                                                                    chi_fit[0],
                                                                    chi_fit_early[0],
                                                                    chi_fit_late[0],
                                                                    fgroup[0],
                                                                    bgroup[0],
                                                                    alpha[0],
                                                                    [-5,ym])
#        chi_dof = chi_dof + [lc,hc]
#        chi_dof_lhc = chi_dof_late + [lce,hce,lcl,hcl]
# 
#        anim_handles = [line,xe,ye,fline,res,
#                        linel,xel,yel,flinel,resl,
#                        linesp,linesm,line2sp,line2sm,
#                        linespl,linesml,line2spl,line2sml,
#                        vertp, codep, bottom, xbin, 
#                        vertl, codel, 
#                        verte, codee,
#                        text1, text2, text3]
    else:
        [[ax_fit,ax_txt],[ax_res,ax_chi]] = ax 
        ym,yM,rm,rM = y.min()-0.05,y.max()+0.01,(y-f_res).min()-0.005,(y-f_res).max()+0.005
        xtgl,ytgl = 1.125*t.min()-0.125*t.max(), 1.5*rm-0.5*rM
        ax_res.text(xtgl,ytgl,'Click to toggle pause/resume',fontsize='medium')
        ax_fit.text(xtgl,yM,dt_string,fontsize='small')
        nu_fit_early,nu_fit_late,dy_early = n*[None], n*[None], n*[None],
        w_early,dy_late,w_late = n*[None], n*[None], n*[None]


        line, xe, ye = errorb(ax_fit,t,y[0],ey[0],color[0])
        fline = plot_fit(ax_fit,tf,f[0],color[1])
        res = plot_res(ax_res,t,y[0]-f_res[0],color[0])
        decorate_data(ax_fit,t,ym,yM)
        linesp,linesm,line2sp,line2sm = decorate_res(ax_res,t,ey[0],rm,rM)
        #vertp, codep, bottom, xbin, patchfit, vertl, codel, verte, codee, patchlate, patchearly
        vertf, nistf, bottom,xbin,_,_,_,_ = plot_chi2(ax_chi,
                                                 dy_fit[0],nu_fit,
                                                 None,None,
                                                 None,None)

        # (ax,nu_fit,nu_early,nu_late,chi_fit,chi_early,chi_late,fgroup,bgroup,alpha,ylim)
        text1, lc, hc,_,_,_,_,_,_ = plot_txt(ax_txt,model,nu_fit,None,None,
                                    chi_fit[0],None,None,
                                    fgroup[0],bgroup[0],alpha[0],[-5,ym])

    ax_0 = ax_early if early_late else ax_fit
    title = run_title
    # print('set_sequence plot debug: next launcing animations')
    anim_fit = animation.FuncAnimation(fig,animate_fit, 
                                            frames=range(len(run_title)),
                                            init_func=init_animate_fit,
                                            interval=anim_delay,
                                            repeat=True,
                                            blit=False)

    fig.canvas.mpl_connect('button_press_event', toggle_pause)
    if rrf:
        stringrrf = r'$\nu_R=$'+'{:.1f}MHz   '.format(rrf)
        supti = P.suptitle(stringrrf+run_title[0])#,x=-0.125,ha='right')    
    else:
        supti = P.suptitle(run_title[0]) # 'large','small'
    draw(fig)
    return fig

def set_figure_fft(fig_fft,model_name,ylabel,f,ap,apf,ep,group,run_title):
    '''
    input:
        fig_fft figure handle
        model_name e.g. "mgml"
        f frequency slice
        a data fft real part slice
        af partial model fft real part slice
        e, scalar, std on a, af
    output the figure (anim for 2,3d a, af)
    '''
    from numpy import arange, array
    from mujpy.tools.tools import autops, ps, _ps_acme_score, _ps_peak_minima_score
    from copy import deepcopy
    import matplotlib.pyplot as P
    from matplotlib.path import Path
    import matplotlib.patches as patches
    import matplotlib.animation as animation

    ###################
    # PYPLOT ANIMATIONS
    ###################
    def animate_fft(i):
        '''
        anim function
        update fft data, fit fft and their color 

        '''
        # color = next(ax_fft._get_lines.prop_cycler)['color']
        ax_fft.set_title(run_title[i])
        marks.set_ydata(ap[i,:])
        segs = [array([[q,w-a],[q,w+a]]) for q,w,a in zip(f,ap[i],ap[0]-ap[0]+ep[0])]
        ye[0].set_segments(segs)
        fline.set_ydata(apf[i,:])
        return marks, ye, fline


    def init_animate_fft():
        '''
        anim init function
        blitting (see wikipedia)
        to give a clean slate 

        '''
        ax_fft.set_title(run_title[0])
        marks.set_ydata(ap[0,:])
        segs = [array([[q,w-a],[q,w+a]]) for q,w,a in zip(f,ap[0],ap[0]-ap[0]+ep[0])]
        # print('plot init_animate_fft debug: ye = {}'.format(ye))
        ye[0].set_segments(segs)
        fline.set_ydata(apf[0,:])
        return marks, ye, fline 

    global paused_fft
    global anim_fft

    paused_fft = False
    def toggle_pause_fft(*args, **kwargs):
        global paused_fft
        global anim_fft
        if paused_fft:
            anim_fft.event_source.start() # matplotlib.__version__ >= 3.4 
            # animation.event_source.start()
        else:
            anim_fft.event_source.stop() # if matplotlib.__version__ >= 3.4 
            # animation.event_source.stop()
        paused_fft = not paused_fft   
    
    ########################
    # build or recall Figure
    ########################
    font = {'family':'Ubuntu','size':10}
    P.rc('font', **font)
    anim_delay = 1000.

    if fig_fft: # has been set to a handle once
        fig_fft.clf()
        fig_fft,ax_fft = P.subplots(num=self.fig_fft.number)
    else: # handle does not exist, make one
        fig_fft,ax_fft = P.subplots(figsize=(6,4))
        fig_fft.canvas.manager.set_window_title('FFT')
    ax_fft.set_xlabel('Frequency [MHz]')
    ax_fft.set_title(run_title[0]+' - '+model_name)
    xm, xM = f.min(),f.max()                    
    ax_fft.set_xlim(xm,xM)


    yM = 1.02*max(ap.max(),apf.max())
    # print('plot set_figure_fft debug: ap.min()={}, apf.min()={}'.format(ap.min(),apf.min())) 
    ym = min(0,1.02*ap.min(),1.02*apf.min())
    if len(ap.shape)==1: # single 
        ax_fft.errorbar(f,ap,'o',yerr=ep,ms=2,alpha=0.8)
        ax_fft.plot(f,apf,'-',lw=1,alpha=0.8)
    ##########################
    # animation initial plot #
    ##########################
    else: # multigroup and/or multirun
        # print('plot set_figure_fft debug: ep = {}'.format(ep[0]))

        marks, xe, ye, = ax_fft.errorbar(f,ap[0,:],fmt='o',yerr=ep[0],ms=2,alpha=0.8)
        fline, = ax_fft.plot(f,apf[0,:],'-',lw=1,alpha=0.8)
        #######
        # anim
        #######
    ax_fft.set_ylim(ym,yM)
    
    ax_fft.text(0.75*xM,ym-0.125*(yM-ym),'Click to pause/resume',fontsize='medium')
#    ydat = 1.01*yM
#    ax_early.text(,ydat,dt_string,fontsize='small')
    if len(ap.shape)>1:
        anim_fft = animation.FuncAnimation(fig_fft, animate_fft, 
                                                frames=range(len(run_title)),
                                                init_func=init_animate_fft,
                                                interval=anim_delay,
                                                repeat=True,
                                                blit=False)
    fig_fft.canvas.mpl_connect('button_press_event', toggle_pause_fft)

#        # print('f.shape = {}, ap.shape = {}'.format(f.shape,ap.shape))  
#        ax_fft.plot(f[k],ap[k],'o',ms=2,alpha=0.5,color=color[k]) # f, ap, apf are     plotiled!
#        ax_fft.plot(f[k],apf[k],'-',lw=1,alpha=0.5,color=color[k])
#        ax_fft.fill_between([f[0,0],f[0,-1]],[k*yoffset,k*yoffset],[k*yoffset+fft_e[k],k*yoffset+fft_e[k]],facecolor=color[k],alpha=0.2)
        ###################
        # errors, alpha_version for single
        ################### 
#                    if self._single_: 

    ax_fft.set_ylabel(ylabel)

    fig_fft.canvas.manager.window.tkraise()
    P.draw()


def errorb(ax,t,y,ey,color):
    '''
    input:
        ax = axis handle, 
        t, y, ey = time data and error
        color = symbol and errorbar color 
    output:
        line, xe, ye = three handles for the animations
    Works for early, late, whole, depending on input 
    To ignore handles for single run: 
            dum, = errorb(...)
    Draws errorbar of data
    '''
    line, xe, ye, = ax.errorbar(t,y,yerr=ey,fmt='o',elinewidth=1.0,
                                    ecolor=color,mec=color,mfc=color,
                                    ms=2.0,zorder=0, alpha=0.5) # errorbar of data
    return line, xe, ye

def plot_fit(ax,t,f,color):
    '''
    input:
        ax, axis handle
        t, f, time fit
        color = line color 
    output:
        fline = fit line handle for animation
    Ignore handle for single run
        dum = plot_fit(...)
    Draws overlayed fit     
    '''
    from matplotlib.pyplot import rcParams

    fline, = ax.plot(t,f,'-',lw=1.5,alpha=1,zorder=2,color=color) # fit
    return fline
     
def decorate_data(ax,t,ym,yM):
    '''
    input:
        ax, axis handle
         
    sets limits, writes x_label, y_label 
    '''
    ax.set_xlim(0,t.max())
    ax.set_ylim(ym,yM)
    ax.set_ylabel('Asymmetry')

def decorate_data_late(ax,t,ym,yM):
    '''
    input:
        ax, axis handle
         
    sets limits,writes x_label
    '''
    from matplotlib import ticker
    ax.set_xlim(t[0],t.max())
    ax.set_ylim(ym,yM)
    ax.yaxis.set_major_formatter(ticker.NullFormatter())# 

def plot_res(ax,t,dy,color): 
    '''
    input:
        ax, axis handle
        t, dy, time residues
        color_index = 0, 2 for early (whole) or late
    output:
        res = residue line handle for animation
    Ignore handle for single run
        dum = plot_res(...)
    Draws residues         
    '''
    res, = ax. plot(t,dy,'-',lw=1.0,alpha=0.8,zorder=2,color=color) # residues 
    return res

def decorate_res(ax,t,ey,rm,rM):
    '''
    input:
        ax, axis handle
         
    sets limits, writes x_label, y_label, 1 & 2 std lines 
    '''
    from matplotlib.pyplot import rcParams

    prop_cycle = rcParams['axes.prop_cycle']
    color = prop_cycle.by_key()['color']
    linesp, = ax. plot(t,ey,'-',lw=1.0,alpha=0.3,zorder=2,color=color[1]) # residues
    linesm, = ax. plot(t,-ey,'-',lw=1.0,alpha=0.3,zorder=2,color=color[1]) # residues
    line2sp, = ax. plot(t,2*ey,'-',lw=1.0,alpha=0.3,zorder=2,color=color[2]) # residues
    line2sm, = ax. plot(t,-2*ey,'-',lw=1.0,alpha=0.3,zorder=2,color=color[2]) # residues
    ax.plot([0,t.max()],[0.,0.],'k-',lw=0.5,alpha=0.3,zorder=0) # zero line
    ax.set_xlim(0,t.max())
    ax.set_ylim(rm,rM)
    ax.set_ylabel('Residues')
    ax.set_xlabel(r'Time [$\mu$s]')
    return linesp, linesm, line2sp, line2sm

def decorate_res_late(ax,t,ey,rm,rM):
    '''
    input:
        ax, axis handle
         
    sets limits, writes x_label, 1 & 2 std lines 
    '''    
    from matplotlib import ticker
    from matplotlib.pyplot import rcParams

    prop_cycle = rcParams['axes.prop_cycle']
    color = prop_cycle.by_key()['color']
    linespl, = ax. plot(t,ey,'-',lw=1.0,alpha=0.3,zorder=2,color=color[1]) # residues
    linesml, = ax. plot(t,-ey,'-',lw=1.0,alpha=0.3,zorder=2,color=color[1]) # residues
    line2spl, = ax. plot(t,2*ey,'-',lw=1.0,alpha=0.3,zorder=2,color=color[2]) # residues
    line2sml, = ax. plot(t,-2*ey,'-',lw=1.0,alpha=0.3,zorder=2,color=color[2]) # residues
    ax.plot([t[0],t.max()],[0.,0.],'k-',lw=0.5,alpha=0.3,zorder=0) # zero line
    ax.set_xlim(t[0],t.max())
    ax.set_ylim(rm,rM)
    ax.yaxis.set_major_formatter(ticker.NullFormatter())# 
    ax.set_xlabel(r'Time [$\mu$s]')
    return linespl, linesml, line2spl, line2sml

def plot_txt(ax,model,nu_fit,nu_early,nu_late,chi_fit,chi_early,chi_late,fgroup,bgroup,alpha,ylim):
    '''
    write stuff on top right
    input: 
        ax
        nu_fit, nu_early, nu_late: dofs
        chi_fit, chi_early, chi_late: reduced chi2
        fgroup, bgroup, alpha: fit group
        ylim: data ylim 
    background white for fit, 
               color[0] for data, data_early
               color[2] for data_late 
    '''
    from scipy.special import gammainc 
    from scipy.stats import norm, chi2
    from numpy import linspace, where
    from matplotlib.pyplot import rcParams

    prop_cycle = rcParams['axes.prop_cycle']
    color = prop_cycle.by_key()['color']
    dylim = ylim[1]-ylim[0]
    dyl = dylim/20
    transalpha = 0.5
    pad = 3
    xtxt = -4.5

    mm = round(nu_fit/4)
    hb = linspace(-mm,mm,2*mm+1)
    cc = gammainc((hb+nu_fit)/2,nu_fit/2) # muchi2cdf(x,nu) = gammainc(x/2, nu/2);
    lc = 1+hb[min(list(where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nu_fit
    hc = 1+hb[max(list(where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nu_fit
    
    string1 = 'F-B: {} - {}\n'.format(fgroup,bgroup)
    string1 += r'$\alpha=$ {:.4f}'.format(alpha)
    string1 += '\n$\chi^2_r=$ {:.3f}\n ({:.2f}-{:.2f})\n{} dof'.format(chi_fit,lc,hc,nu_fit)
    text1 = ax.text(xtxt,ylim[0]+0.55*dylim,string1,bbox={'facecolor': 'white', 'pad': pad}) 
    ax.axis('off')
    ax.set_ylim(ylim)

    # print('plot_txt in tools.plot debug: chi_late = {}'.format(chi_late))
    if chi_late is not None: 
        mm = round(nu_early/4)
        hb = linspace(-mm,mm,2*mm+1)
        cc = gammainc((hb+nu_early)/2,nu_early/2) # muchi2cdf(x,nu) = gammainc(x/2, nu/2);
        lce = 1+hb[min(list(where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nu_early
        hce = 1+hb[max(list(where((cc<norm.cdf(1))&(cc>norm.cdf(-1))))[0])]/nu_early
        string2 = '$\chi^2_e$ = {:.3f}\n({:.2f}-{:.2f})'.format(chi_early,lce,hce) 
        mml = round(nu_late/4)
        hbl = linspace(-mml,mml,2*mml+1)
        ccl = gammainc((hbl+nu_late)/2,nu_late/2) # muchi2cdf(x,nu) = gammainc(x/2, nu/2);
        lcl = 1+hbl[min(list(where((ccl<norm.cdf(1))&(ccl>norm.cdf(-1))))[0])]/nu_late
        hcl = 1+hbl[max(list(where((ccl<norm.cdf(1))&(ccl>norm.cdf(-1))))[0])]/nu_late
        string3 = '$\chi^2_l$ = {:.3f}\n({:.2f}-{:.2f})'.format(chi_late,lcl,hcl) 
        text2 = ax.text(xtxt,ylim[0]+0.28*dylim,string2,
                        bbox={'facecolor': color[0], 'alpha': transalpha, 'pad': pad})
        dyl = dylim/50
        text3 = ax.text(xtxt,ylim[0]+dyl,string3,bbox={'facecolor': color[2], 'alpha': transalpha, 'pad': pad})
        text4 = ax.text(xtxt,ylim[0]+dylim,model+' fit')
        return text1, lc, hc, text2, lce, hce, text3, lcl, hcl

    return text1, lc, hc, None, None, None, None, None, None
    

def plot_chi2(ax,dy_fit,nu_fit,dy_early,w_early,dy_late,w_late):
    '''
    input:
        ax, axis handle
        dy_fit, normalized deviation for fit, on fi-range
        dy, normalized deviation for y, plot (full range or early)
        w_plot, n_dof(fit)/n_dof(data_plot)*data_plot.shape[0]
        dy_late, normalized deviation for y_late, plot (late)  or None
        w_late, n_dof(fit)/n_dof(data_plot)*data_plot.shape[0] or None
    histogram of chi2 distribution
    '''
    from numpy import linspace, ndarray, histogram
    from scipy.stats import norm
    from matplotlib import ticker
    from matplotlib.pyplot import rcParams
    import matplotlib.path as path
    import matplotlib.patches as patches

    prop_cycle = rcParams['axes.prop_cycle']
    color = prop_cycle.by_key()['color']
    alpha = 0.5
    ########################
    # chi2 distribution: fit
    ########################
    xbin = linspace(-5.5,5.5,12)
    # self.ax_fit[(1,-1)].set_ylim(0, 1.15*nhist.max())
    #nhist = ax.hist(dy_fit[0],xbin,rwidth=0.9,fc='w',ec='k',lw=0.7)
    nhistf,_ = histogram(dy_fit,xbin) # fc, lw, alpha set in patches
    vertf, codef, bottom, xlim = set_bar(nhistf,xbin) # fc, lw, alpha set in patches
    barpathf = path.Path(vertf, codef)
    patchfit = patches.PathPatch(barpathf,ec='k', facecolor='w', alpha=alpha,lw=0.7,zorder=2)
    ax.add_patch(patchfit)  #hist((yfit-ffit)/eyfit,xbin,rwidth=0.9,fc='w',ec='k',lw=0.7)
    ax.set_xlim(xlim[0],xlim[1])
    
    #########################################
    # chi2 distribution: plots, scaled to fit
    #########################################
    if isinstance(dy_late,ndarray):
        nhistl,_ = histogram(dy_late,xbin,weights=w_late)
        vertl, codel,_,_ = set_bar(nhistl,xbin) # fc, lw, alpha set in patches
        barpathl = path.Path(vertl, codel)
        patchlate = patches.PathPatch(barpathl, ec=color[2], facecolor=color[2], alpha=alpha,lw=0.7)
        ax.add_patch(patchlate)  #hist((yfit-ffit)/eyfit,xbin,rwidth=0.9,fc='w',ec='k',lw=0.7)
        nhiste,_ = histogram(dy_early,xbin,weights=w_early)
        # ],weights=w_early[0],rwidth=0.9,fc=color[0],alpha=alpha)
        verte, codee,_,_ = set_bar(nhiste,xbin) # fc, lw, alpha set in patches
        barpathe = path.Path(verte, codee)
        patchearly = patches.PathPatch(barpathe, ec=color[0], facecolor=color[0], alpha=alpha,lw=0.7)
        ax.add_patch(patchearly)  #hist((yfit-ffit)/eyfit,xbin,rwidth=0.9,fc='w',ec='k',lw=0.7)
    ###############################
    # chi2 dist theo curve & labels 
    ###############################
    xh = linspace(-5.5,5.5,23)
    yh = norm.cdf(xh+1)-norm.cdf(xh)
    ax.plot(xh+0.5,nu_fit*yh,'-',color=color[1])
    ax.set_xlabel("$\sigma$")
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.set_xlim([-5.5, 5.5])
    if isinstance(dy_late,ndarray):
        return vertf, nhistf, bottom, xbin, vertl, nhistl, verte, nhiste
    else:
        return vertf, nhistf, bottom, xbin, None, None, None, None

def draw(fig):
    '''
    '''
    import matplotlib.pyplot as P
    cfm = P.get_current_fig_manager()
    if hasattr(cfm.window,'attributes'):
        cfm.window.attributes('-topmost', True) # Tkinter backend
        cfm.window.attributes('-topmost', False)
    #cfm.window.activateWindow()
    #cfm.window.raise_()
    #fig.canvas.manager.window.tkraise()# fig.canvas.manager.window.raise_()
    P.draw()
    
