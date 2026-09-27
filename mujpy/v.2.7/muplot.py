class multiplot(object):
    '''
    plot class (let's see)
    '''
    def __init__(self,time,asymm,title,nscan,histoLength):
        '''
        input: if suite is the multiple run instance
            time, asymm - 1d and 2d numpy arrays 
                  e.g. from rebin(suite.time,suite.asymmetry_multirun(),(0,20000),100) 
            title  - list
                  e.g. [get_title(run[0]) for run in suite._the_runs_]
            nscan - list
                  e.g. [run[0].get_runNumber_int() for run in suite._the_runs_]
                       [groups = [grp['forward']+'-'+grp['backward'] for grp in the_suite.groups]
            histoLength - max length of each asymmetry array
                  e.g.  musuite.histoLength
        method:
            multiplot.display(anim=True) 1s sequence of run plots, paused by muose click
        '''
        from numpy import array
        self.time = array([time])
        self.asymm = asymm
        self.title = title
        self.scan = nscan
        self.histoLength = histoLength
        self.multi_offset = 0.1
    
    def display(self,groups = False, anim = False, anim_delay = 1000):       
        '''
        input: 
        produces plot
        multiplot_range = self.multiplot_range
        anim = True, False
        anim_delay = delay between frames (ms)
        output: MULTIPLOT display:
        anim_multiplot to be .paused() and .resumed()  by toggle_pause
        '''
        import matplotlib.pyplot as P
        from numpy import array
        from mujpy.aux.aux import set_fig, derange, rebin #, animate_multiplot, init_animate_multiplot
        import matplotlib.animation as animation

        ###################
        # PYPLOT ANIMATIONS
        ###################
        def animate_multiplot(i):
            '''
            anim function
            update multiplot data and its color 
            '''
            line.set_ydata(self.asymm[i])
            line.set_color(color[i])
            self.ax.set_title(str(self.scan[i])+': '+ self.title[i])
            return line, 


        def init_animate_multiplot():
            '''
            anim init function
            to give a clean slate 
            '''
            line.set_ydata(self.asymm[0])
            line.set_color(color[0])
            self.ax.set_title(str(self.scan[0])+': '+ self.title[0])
            return line, 


        def toggle_pause(*args, **kwargs):
            if self.paused:
                self.anim_multiplot.resume() # event_source.start() # matplotlib.__version__ >= 3.4 
                # animation.event_source.start()
            else:
                self.anim_multiplot.pause() #event_source.stop() # if matplotlib.__version__ >= 3.4 
                # animation.event_source.stop()
            self.paused = not self.paused   
        dpi = 100.
        if len(self.asymm.shape)==1:
            anim = False # make sure
            nscans, nbins = 1, self.asymm.shape[0]
        else:
            nscans,nbins = self.asymm.shape

        #print('start, stop, pack = {},{},{}'.format(start,stop,pack))
        #print('shape time {}, asymm {}'.format(time.shape,asymm.shape))
        y = 4. # normal y size in inches
        x = 6. # normal x size in inches
        my = 12. # try not to go beyond 12 run plots

        ##############################
        #  set figure, axes 
        ##############################
        kwargs = {'figsize':(x,y),'dpi':100.}  
       
        fig, self.ax = set_fig(1,1,1,'Multiplot',**kwargs)

        screen_x, screen_y = P.get_current_fig_manager().window.wm_maxsize() # screen size in pixels
        y_maxinch = float(screen_y)/float(fig.dpi) # maximum y size in inches

        ########## note that "inches" are conventional, since they depend on the display pitch  
        # print('your display is y_maxinch = {:.2f} inches'.format(y_maxinch))
        ########## XPS 13 is 10.5 "inches" high @160 ppi (cfr. conventional fig.dpi = 100)
        bars = 1. # overhead y size(inches) for three bars (tools, window and icons)
        dy = 0. if anim else (y_maxinch-y-1)/my   # extra y size per run plot
        y = y + nscans*dy if nscans < 12 else y + 12*dy # size, does not dilate for anim 
        # fig.set_size_inches(x,y, forward=True)

        ##########################
        #  plot data and fit curve
        ##########################
        color = []
        color.append(next(self.ax._get_lines.prop_cycler)['color'])
        for run in range(1,nscans):
            color.append(next(self.ax._get_lines.prop_cycler)['color'])
        anim_multiplot = []
        if anim:
        #############
        # animation
        #############
            ##############
            # initial plot
            ##############
            ylow, yhigh = self.asymm.min()*1.02, self.asymm.max()*1.02
            line, = self.ax.plot(self.time[0],self.asymm[0],'o-',ms=2,lw=0.5,color=color[0],alpha=0.5,zorder=1)
            self.ax.set_title(str(self.scan[0])+': '+self.title[0])
            self.ax.plot([self.time[0][0],self.time[0][-1]],[0,0],'k-',lw=0.5,alpha=0.3)
            self.ax.set_xlim(self.time[0][0],self.time[0][-1])
            self.ax.set_ylim(ylow,yhigh)
            self.ax.set_ylabel('Asymmetry')
            self.ax.set_xlabel(r'time [$\mu$s]')
            #######
            # anim
            #######
            self.anim_multiplot = animation.FuncAnimation(fig, animate_multiplot, nscans,       
                                                     init_func=init_animate_multiplot,
                                                     interval=anim_delay, blit=False)
            self.paused = False
            fig.canvas.mpl_connect('button_press_event', toggle_pause)
            P.suptitle('Click to toggle pause/resume',fontsize='small')            

        ###############################
        # tiles with offset
        ###############################
        else: 
            aoffset = self.asymm.max()*self.multi_offset*array([[run] for run in range(nscans)])
            self.asymm = self.asymm + aoffset # exploits numpy broadcasting
            ylow,yhigh = min([0,self.asymm.min()+0.01]),self.asymm.max()+0.01
            if nscans>1:
                for run in range(nscans):
                    self.ax.plot(self.time[0],self.asymm[run],'o-',
                                lw=0.5,ms=2,alpha=0.5,color=color[run],zorder=1)
                    self.ax.plot([self.time[0][0],self.time[0][-1]],
                                           [aoffset[run],aoffset[run]],'k-',lw=0.5,alpha=0.3,zorder=0)
                    self.ax.text(self.time[-1]*1.025,aoffset[run],self.run[run])
                self.ax.set_title(self.title[run])
            else:
                self.ax.plot(self.time,self.asymm,'o-',lw=0.5,ms=2,alpha=0.5,color=color[0],zorder=1)
                self.ax.set_title(self.title)            
            self.ax.set_xlim(self.time[0][0],self.time[0][-1]*9./8.)
            self.ax.set_ylim(ylow,yhigh)
            # print('axis = [{},{},{},{}]'.format(time[0,0],time[0,-1]*9./8.,ylow,yhigh))
            self.ax.set_ylabel('Asymmetry')
            self.ax.set_xlabel(r'time [$\mu$s]')
            # self.fig_multiplot.tight_layout()
        fig.canvas.manager.window.tkraise()
        P.draw()
        return anim_multiplot

