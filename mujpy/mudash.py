class dash(object):
    '''
    A graphic handler class for the mujpy dashboard
    python3 only.
    Runs in jupyter notebook inside jupyter lab (OS independent).
    Based on ipywidgets, jupyterlab tab
    The main method designs a tabbed gui to be run in a Jupyter notebook. 
    Public methods correspond to the main tabs (fft, fit, setup, suite, ...). 
    
    **** Output *****
        working solution in mujpy/getstarted/Delendo/MuDosh.ipynb
        __init__ calls self.dash_init which defines board with html style overflow: auto;, display: flex; flex-direction: column-reverse;
        console, gui, fit, group, about are private methods (no self) of self.dash
      
      later (eg on_load_nxt) both end up TWICE! in log console, despite same sys.stdout & err and same __stdout__ & __stderr__

    Relies on the mujpy classes
       musuite
       [muplot]
       mufit
       mufitplot
       mumodel in mucomponents/mucomponents
       muedge        >>
       muprompt      >>
    and on the library 
       aux in aux/aux
       
    function __init__(self)

      * initiates an instance and a few attributes, 
      * launches the gui.
      * Use as follows::

             from mugui import mudash
             dashboard = mudash() # instance is dashboard

    '''
##########################
# INIT
##########################
    def __init__(self, suite, model = 'alml', json_file='', version = '1', fit_range = [0,20000,4], load=False, userpar = []):
        '''
        * Launches the gui, initiates a few attributes, 
        * pass empty model, e.g. model = 'mgmgbl' for a two site polycrystal magnet, 
        *            [default 'alml' for grp calibration with one muon component]
        * pass best fit model with json = filename [default json = '']
        * fit_range is a list of integer values, passed on to fit
        * version is a label, fit range a list, 
        * load = True requests gui to load best-fit [default load = False] 
        * userpar = list of dictionaries, keys 'name', 'value', 'std', for global fits, [default is [], no globals]
        *
        * Use a template notebook, basically as follows::

             from mujpy.mudash import mudash as MG
             from mujpy.musuite import suite
             the_suite = suite(,mplot=False)
             MuJPy = MG(mysuite,model='mgml') # instance is MuJPy

        '''
        # print('__init__ now ...')
#        from mujpy import __file__ as MuJPyName
        from mujpy.aux.aux import validmodel
        import os
        import numpy as np
        from IPython.display import display# ,  HTML
       
# some initializations
        try:
            suite.single()
        except:
            print('no valid suite! aborting ...')
            return 
        self.suite = suite # this is an instance of suite that must be created beforehand
        
        # inputs to methods with sanity checks
        if not isinstance(model,str): model ='alml' # default
        self.model = model if validmodel(model) else 'alml' # this is a model name of type 'mgml', with default 'alml' calibration fit
        self.model_components = [] # initialized, will contain model_guess or model_results 
        if not isinstance(version,str): version ='1'          
        self.version = version if version else '1' # a free label for this fit, default '1'
        self.load = load # if load then a file selection gui for .json files is presented
        self.dashboard_file = json_file # this is an optional json dashboard in the path

        self.__logopath__ = os.path.join(self.suite.__path__,"logo") # logo path
        ############################################################################
        # other paths already set in suite
        #       self.suite.__startuppath__ is  working directory
        #       self.suite.__templatepath__ is  mujpy root directory + '/templates/
        #       self.suite.__fitpath__ is  working directory + '/fit/
        #       self.suite.__csvpath__ is  working directory + '/csv/
        #       self.suite.__cachepath__  is  working directory + '/cache/ (for copy of summary)
        ############################################################################

        self.mujpy_width = '880px'   
        ############################    
        # output_width = '660px'   #
        # log_height = '600px'     #
        ############################
        self.button_color = 'lightblue'
        self.border_color = '2px solid dodgerblue'

        self.buffer = []  # per board Output
        self.the_fit = []
        self.plot_range = fit_range
        self.fig_fit = None
        self.fig_fft = None
        for par in userpar:
            if 'name' not in par:
                print('*** Global parameters in userpar list must contain at least their names')
                print("    suggested e.g. [{'name':'B[mT]},{'name':'λ[µs-1]'}]")
                print('aborting ...') 
                return
        self.userpar = userpar
        self.dash_init(fit_range)
        
    def dash_init(self,fit_range):
        from ipywidgets import Output, HTML, HBox, VBox, Tab,Layout
        from IPython.display import display        
        # print('**** dash_init Begins ****')
        ##########################################
        # This actually produces the gui interface
        # with tabs: fit group about
        ##########################################        
        board = Output()
        board.add_class('tail-f')
        # html_style = HTML('<style>.tail-f{height: 420px; width: 550px; border:2px solid limegreen; overflow:auto; display:flex; flex-direction:column;color-scheme:light}</style>')
        html_style = HTML('<style>.tail-f {height: 420px;'\
                                         'width: 650px;'\
                                         'border: 3px solid limegreen;'\
                                         'background-color: #ddd;'\
                                         'color: #ffaaaa !important;'\
                                         'overflow:auto;'\
#                                         'display:flex;'\
#                                         'flex-direction:column-reverse;}'\
                                         '</style>')
        tabs_content = ['fit ','group','about ']  # to add a tab: add its name as tabs_content[n], 
        # add a module def name(mainwindow): returning mainwindow.children[n].children
        tabs = [VBox(description=name,layout=Layout(border=self.border_color)) for name in tabs_content]
        self.mainwindow = Tab(children = tabs,layout=Layout(width=self.mujpy_width,border=self.border_color))# 
        self.mainwindow.selected_index = 0 # to stipulate that the first display is on tab 0, fit
        for i,tab_title in enumerate(tabs_content):
            self.mainwindow.set_title(i, tab_title)

        def console(string):
            self.buffer.append(string)
            reversebuf = list(reversed(self.buffer))
            board.clear_output()
            with board:
               for line in reversebuf:
                    print(line)
        def erlog(string):
            print(string)
        self.erlog = erlog 
        self.log = console
        self.suite.console = console # redirects suite output to log Output widget
        self.log('********************** LOG ABOVE HERE FILLS UPWARDS **********************') 
                # 12345678901234567890123456789012345678901234567890123456789012345678901234
                # 0        1         2         3         4         5         6         7
        # now tabs are defined as methods: about, gui, group, fit 
        # dash_init continues at the bottom, after all these methods, and invokes them
##########################
# ABOUT
##########################
        def about():
            '''
            about tab:

            - a few infos (version and authors)

            '''
            from ipywidgets import Textarea, Layout, HTML, VBox

            _version = 'MuJPy          version '+'2.0' # increment while progressing
            _authors = '\n\n  Authors: Roberto De Renzi, Pietro Bonfà (*)'
            _blahblah = ('\n\n A Python MuSR data analysis tool'+
                         '\n  designed for Jupyterlab.'+
                         '\n  Released under the MIT licence')
            _html = ('\n  See the <a href="https://mujpy.readthedocs.io/en/latest/Tutorial.html"target="_blank">ReadTheDocs Tutorial</a>')
            _pronounce = ('\n  Pronounce it as mug + pie')
            _additional_credits_ = ('\n ---------------------\n (*) dynamic Kubo-Toyabe algorithm by G. Allodi\n     MuSR_td_PSI by A. Amato and A.-R. Raselli \n     acme algorithm code from NMRglue, by Jonathan J. Helmus')
            _about_text = _version+_blahblah+_pronounce+_authors+_additional_credits_
            _about_area = Textarea(value=_about_text,
                                       description='MuJPy',
                                       layout=Layout(width='98%',height='250px'),
                                       disabled=True)
            _about_html = HTML(value=_html,
                                       description='Docs',
                                       layout=Layout(width='98%',height='250px'),
                                       disabled=True)
            # now collect the handles of the three horizontal frames to the main fit window (see self.tabs_contents for index)
            abouttab = [VBox([_about_area,_about_html])]  # add the list of widget handles as the third tab, fit
            return abouttab


##########################
# GUI
##########################
        def gui():
            '''
            Main gui layout. Based on ipywidgets.
            it designs 
            
            * an external frame
            
            * the logo and title header
            
            * the tab structure.

            It requires an instance of suite being passed to it
            Thus the jupyter notebook first cells (prologue) determine:
            
                datafile (the prototype, also selecting data and log dirs:
                    data dir wherever, 
                    log dir is ./log from the folder where jupyterlab is launched/run) 
                group_calib dictionary
                input_suite dictionary
                global dictionary 
                model = model
            and then create 
                a suite object, 
                a mudash object and passing the first to the second        
            so that 
                single or multi group
                single or multi run 
                and global flags determine the layout of gui    
                (e.g load next previous appears only for single run
                 add is dealt within suite)
                 
            mudash shows  

                run number, title, 
                total counts, group counts, ns/bin
                comment, start stop date, next run
                for the run or for the prototype run of the suite 

            At the end (Araba.Phoenix) the gui method redefines self.gui 
            as a Vbox named 'whole', containing the entire gui structure                      
            '''
            def on_load_nxt(b):
                '''
                load next run (if it exists)
                '''
                from time import sleep
                from mujpy.musuite import suite 
                from mujpy.aux.aux import nextrun, get_title, get_totals, get_datafilename
                # print('self.nrun[0] = {}'.format(self.nrun[0]))

                if text_next_label.value:
                    grp_calib = self.suite.groups
                    offset = self.suite.offset
                    runnext, datafile = nextrun(self.suite.datafile)
                    # self.log('Next file is {}'.format(datafile))
                    self.suite = suite(datafile, text_next_label.value , 
                                       grp_calib , offset, self.suite.__startuppath__,
                                       console='self.dash.log',dash=self)
                    # update the Gui Header
                    text_run.value = self.suite.runs[0][0]
                    text_title.value = get_title(self.suite._the_runs_[0][0])
                    text_comment.value = self.suite._the_runs_[0][0].get_comment()
                    # self.log('file is {}'.format(self.suite.datafile))
                    text_filepath.value = self.suite.datafile
                    text_totalcounts.value, text_groupcounts.value, text_maxbin.value, text_nsbin.value = get_totals(self.suite)
                    # self.log('Next = {}'.format(str(int(self.suite.runs[0][0])+1)))
                    if os.path.exists(get_datafilename(self.suite.datafile,str(int(self.suite.runs[0][0])+1))): # next file
                        text_next_label.value = str(int(self.suite.runs[0][0])+1)
                        text_next_label.background_color = 'white'
                        text_start_date.value = ' '.join(self.suite._the_runs_[0][0].get_timeStart_vector())
                        text_stop_date.value = ' '.join(self.suite._the_runs_[0][0].get_timeStop_vector())
                    else:
                        text_next_label.value = ''
                        text_next_label.background_color = "mistyrose"
    #                runnext, datafile = nextrun(self.suite.datafile) 
    #                text_next_label.value = runnext if runnext != run else ''
    #                text_next_label.background_color = 'white' if  runnext != run else 'mistyrose'
                # print('self.nrun[0] = {}'.format(self.nrun[0]))
                else:
                    text_next_label.background_color = 'red'
                    sleep(0.2)
                    text_next_label.background_color = 'mistyrose'
                    

            def on_load_datafile(b):
                '''
                load from gui
                '''
                from mujpy.musuite import suite 
                from mujpy.aux.aux import nextrun, path_file_dialog,  get_title, get_datafilename, get_totals
                import os
                import re
                # from datetime import datetime
                # print('self.nrun[0] = {}'.format(self.nrun[0]))
                grp_calib = self.suite.groups
                offset = self.suite.offset  
                path, filename  = os.path.split(self.suite.datafile)          
                datafile =  path_file_dialog(path,filename[filename.rfind('.')+1:])
                # self.log('debug mudash on_load_datafile File {}'.format(datafile))
                run = str(int(next(re.finditer(r'\d+$', os.path.splitext(datafile)[0])).group(0)))
                try:    
                    # now = datetime.now().time() # time object
                    # self.log("now = {}".format(now))                   
                    self.suite = suite(datafile, run, grp_calib , offset,
                                        self.suite.__startuppath__,console='self.dash.log',dash=self)
                    # update the Gui Header
                    text_run.value = self.suite.runs[0][0]
                    text_title.value = get_title(self.suite._the_runs_[0][0])
                    text_comment.value = self.suite._the_runs_[0][0].get_comment()
                    text_filepath.value = self.suite.datafile
                    text_totalcounts.value, text_groupcounts.value, text_maxbin.value, text_nsbin.value = get_totals(self.suite)
                    runnext, datafile = nextrun(self.suite.datafile) 
                    text_next_label.value = runnext if runnext != run else ''
                    text_next_label.background_color = 'white' if  runnext != run else 'mistyrose'
                    text_start_date.value = ' '.join(self.suite._the_runs_[0][0].get_timeStart_vector())
                    text_stop_date.value = ' '.join(self.suite._the_runs_[0][0].get_timeStop_vector())
                except Exception as e:
                    self.erlog(e)
                    self.erlog('Run {} does not exist'.format(run))
                # print('self.nrun[0] = {}'.format(self.nrun[0]))
                    
            def on_load_prv(b):
                '''
                load previous run (if it exists)
                '''
                from time import sleep
                from mujpy.musuite import suite 
                from mujpy.aux.aux import prevrun, get_title, get_totals
                runprev, datafile = prevrun(self.suite.datafile)                  
                # print('self.nrun[0] = {}'.format(self.nrun[0]))
                grp_calib = self.suite.groups
                offset = self.suite.offset
                if runprev!=str(self.suite.runs[0][0]):        
                    # self.log('Next = {}'.format(str(int(self.suite.runs[0][0])+1)))
                    text_next_label.value = str(self.suite.runs[0][0])
                    self.suite = suite(datafile, runprev, grp_calib, offset,
                                       self.suite.__startuppath__,console='self.dash.log',dash=self)
                    text_next_label.background_color = 'white'                
                    text_run.value = self.suite.runs[0][0]
                    text_title.value = get_title(self.suite._the_runs_[0][0])
                    text_comment.value = self.suite._the_runs_[0][0].get_comment()
                    # self.log('file is {}'.format(self.suite.datafile))
                    text_totalcounts.value, text_groupcounts.value, text_maxbin.value, text_nsbin.value = get_totals(self.suite)
                    text_filepath.value = self.suite.datafile
                    text_start_date.value = ' '.join(self.suite._the_runs_[0][0].get_timeStart_vector())
                    text_stop_date.value = ' '.join(self.suite._the_runs_[0][0].get_timeStop_vector())
                    # update the Gui Header
                else: # runprev does not exist
                    self.log('Run {} does not exist'.format(runprev))
                    text_title.background_color = 'mistyrose'  
                    sleep(0.2)
                    text_title.background_color = 'white'  
                    
                # print('self.nrun[0] = {}'.format(self.nrun[0]))

                                           
            from ipywidgets import Image, Text, Layout, HBox, Output, VBox, Tab
            from ipywidgets import Button, Label, Box, HTML
            from mujpy.aux.aux import get_totals, get_title, get_datafilename   
            import os                      

            file = open(os.path.join(self.__logopath__,"logo.png"), "rb")
            image = file.read()
            logo = Box([Image(value=image)],
                        layout=Layout(width='114px',
                        height='100px'))#,object_fit='none',margin='0 0 0 0'))
            
            
            width = '100%'
    #----------------- firsttrow
            text_run = Text(description='# ',
                                     description_tooltip='Run number',
                                     value=self.suite.runs[0][0],
                                     layout=Layout(width='8%'),
                                     disabled=True)                        # 8%
            text_run.style.description_width = '15%'

            text_title = Text(description_tooltip ='Title', 
                              value=get_title(self.suite._the_runs_[0][0]),
                              layout=Layout(width='24%'),
                              disabled=True)                               # 32%
            text_title.style.description_width = '20%'
            text_comment = Text(value = self.suite._the_runs_[0][0].get_comment(),
                                description_tooltip='Comment',
                                layout=Layout(width='38%'),
                                disabled=True)                              # 70%
            ns, gs , nsbin, maxbin= get_totals(self.suite)
            text_totalcounts = Text(value=ns,
                                    description='Tot ',
                                    description_tooltip='Total e+ counts',
                                    layout=Layout(width='13%'),
                                    disabled=True)                         # 82%
                                    # events                               # 84%
            text_totalcounts.style.description_width = '18%' 
            text_groupcounts = Text(value=gs,
                                    description='Grp ',
                                    description_tooltip='Group e+ counts',
                                    layout=Layout(width='13%'),
                                    disabled=True)                         # 96%
                                    # events                               # 100%
            text_groupcounts.style.description_width = '20%'
            event_units = Label('ev ',layout={'width':'2%','height':'16pt'})   # 2%                      

     #----------------- secondtrow
            # text_comment.style.description_width = '10%'
            spacer = Label(' ',
                                 layout={'width':'2%','height':'16pt'})  
            text_filepath = Text(value=self.suite.datafile,
                                 layout={'width':'70%','height':'16pt'})   # 72%
            text_maxbin = Text(description='Max',
                               value = maxbin,
                               description_tooltip='bins per histogram',
                               layout=Layout(width='12%'),
                               disabled=True)                               # 76%
            text_maxbin.style.description_width = '25%'                           

            maxbin_units = Text(value =' bins',
                                layout=Layout(width='6%'),
                                disabled=True)                              # 82%

            text_nsbin = Text(value = nsbin ,
                              description_tooltip='ns/bin',
                              layout=Layout(width='10%'),
                              disabled=True)                                # 92%
            nsbin_units = Text(value = ' ns/bin',
                               layout=Layout(width='8%'),
                               disabled=True)                               # 100%

            # text_nsbin.style.description_width = '40%'


    #----------------- thirdrow
            text_start_date = Text(description='Start ',
                                         value=' '.join(self.suite._the_runs_[0][0].get_timeStart_vector()) ,
                                         layout=Layout(width='25%'),
                                         disabled=True)                    # 25%
            text_start_date.style.description_width = '20%'
    #        exit_button = Button(description='Exit',
    #                           tooltip='cleanly exit mujpy',
    #                           layout=Layout(width='8%'))                   # 36%
    #        exit_button.on_click(on_exit)                   
    #        exit_button.style.button_color = self.button_color
            spacer = Label(description='     ',layout=Layout(width='2%'))   # 27%
            text_stop_date = Text(description='Stop ',
                                         value=' '.join(self.suite._the_runs_[0][0].get_timeStop_vector()) ,
                                         layout=Layout(width='25%'),
                                         disabled=True)                    # 52%
            text_stop_date.style.description_width = '20%'
            
    #        Ap_button = Button(description='<Add',
    #                           tooltip='Add previous run\n(refers to Last+\nor to Next#)',
    #                           layout=Layout(width='7%'))                   # 
    #        Ap_button.on_click(on_add_prv)
    #        Ap_button.style.button_color = 'darkgray'#self.button_color
    #        
    #        An_button = Button(description='Add>',
    #                           tooltip='Add next run\n(refers to Last+\nor to Next#)',
    #                           layout=Layout(width='7%'))                   # 
    #        An_button.on_click(on_add_nxt)
    #        An_button.style.button_color = 'darkgray'#self.button_color

    #        last_add = Text(description='Last +',
    #                        description_tooltip='Last added run',
    #                        disabled=True,
    #                        layout=Layout(width='18%'))                     # 
    #        last_add.style.description_width = '27%'

            # spacer                                                        # 54%
            Lp_button = Button(description='<',
                               tooltip='Loads previous run\n(if it exists)',
                               layout=Layout(width='6%'))                   # 60%
            Lp_button.on_click(on_load_prv)
            Lp_button.style.button_color = self.button_color
            Ld_button = Button(description='Load data',
                               tooltip='Opens a run file GUI\nin the data path (see setup tab)',
                               layout=Layout(width='10%'))                   # 70% 
            Ld_button.style.button_color = self.button_color
            Ld_button.on_click(on_load_datafile)

            Ln_button = Button(description='>',
                               tooltip='Loads next # run\n(if it exists)',
                               layout=Layout(width='6%'))                   # 76%
            Ln_button.on_click(on_load_nxt)
            Ln_button.style.button_color = self.button_color
            if os.path.exists(get_datafilename(self.suite.datafile,str(int(self.suite.runs[0][0])+1))): # next file
                value = str(int(self.suite.runs[0][0])+1)
                color = 'white'
            else:
                value = ''
                color = "mistyrose"
            text_next_label = Text(description='Next #',
                              value=value,
                              background_color = color,
                              description_tooltip='Next run to load',
                              disabled=True,
                              layout=Layout(width='14%'))                   # 100%
            text_next_label.style.description_width = '50%'


            firstrow = HBox(layout=Layout(width=width))
            firstrow.children = [text_run   , text_title, 
                                 text_comment,
                                 text_totalcounts, event_units,
                                 text_groupcounts, event_units]
            secondrow = HBox(layout=Layout(width=width))
            secondrow.children = [spacer,text_filepath, 
                                  text_maxbin, maxbin_units, 
                                  text_nsbin, nsbin_units]
            thirdrow = HBox(layout=Layout(width=width))
            thirdrow.children = [text_start_date, spacer, 
                                 text_stop_date, spacer]
            if self.suite.single():
                 thirdrow.children = [text_start_date, spacer, 
                                 text_stop_date, spacer, spacer,spacer,
                                 Lp_button, Ld_button, Ln_button, 
                                 text_next_label]

            titlewindow = VBox(layout=Layout(width='100%')) 
            titlewindow.children = [firstrow, secondrow, thirdrow] 
            
            titlelogowindow = HBox(layout=Layout(width=self.mujpy_width))
            titlelogowindow.children = [logo, titlewindow]

            # main layout: tabs
                          

            return titlelogowindow
                        	
##########################
# GROUP TAB
##########################
        def group():
            '''
            gui version of grouping
            passed to suite via input (shorthand) grp_calib so that suite.groups = grp_calib 
            converted to lists by get_grouping (aux) in suite.grouping
            and passed to mudash via suite
            '''
        ######### here starts the group method of mudash ####################################################
            from ipywidgets import Text, Layout, Button, HBox, FloatText, \
                                   Checkbox, VBox, Dropdown, ToggleButton, Label
            from mujpy.aux.aux import _available_components_
            textheight = '23px'
            labelheight = '23px'
            buttonheight = '48px'
            bwidth = '92%'
            # self.suite.grouping is a list of dictionaries with the two np.arrays of detetctor indices, forw and backw, plus an alpha value
            rows = [] # empty VBox content
            groups = [] # empty HBox content
            # group_index = k where k is the index of groups
            # no need for observe, the input values for group k are already stored as 
            # groups[3*k]   0 3 6 ... :  forward  (np.array of indices)
            # groups[3*k+1] 1 4 7 ... :  backward (np.array of indices)
            # groups[3*k+2] 2 5 8 ... :  alpha    (float)
            #
            for k,grp in enumerate(self.suite.grouping): # here groups are lists
                # group_index += 1            
                groups.append(Text(value=','.join(str(e+1) for e in grp['forward']), # convert np,array([1,2,3]) into '1,2,3'
                                   description='forward', # value is a np.array ofintegers
                                   layout=Layout(width='17%',height=textheight),
                                   description_tooltip='forward group\n1:4 or 1,2,3,4\nor 1,3:5',
                                   continuous_update=False))                                  # 17%
                groups[k%2*3].style.description_width='40%'
                #groups[k%2*3].observe(on_forward,'value')
                groups.append(Text(value=','.join(str(e+1) for e in grp['backward']),
                                   description='backward',
                                   layout=Layout(width='20%',height=textheight),
                                   description_tooltip='backward group\n6:9 or 6,7,8,9\nor 6,8:10',
                                   continuous_update=False))                                  # 37%
                groups[k%2*3+1].style.description_width='40%'
                #groups[k%2*3+1].observe(on_backward,'value')
                groups.append(FloatText(value=grp['alpha'],
                                        description='α',    # value is a float
                                        layout=Layout(width='12%',height=textheight),
                                        description_tooltip='α=Nf/Nb\nratio of forward to backward\nunpolarized count rates',
                                        continuous_update=False))                                  # 49%
                groups[k%2*3+2].style.description_width='20%'
                #groups[k%2*3+2].observe(on_alpha,'value')
                if not (k+1)%2:
                    rows.append(HBox(groups,layout=Layout(width='100%',height='100%')))
                    groups = []
            if groups:
                rows.append(HBox(groups,layout=Layout(width='100%',height='100%')))
            grouptab = [VBox(rows,layout=Layout(border = self.border_color,width='100%',height='100%'))] 
            return grouptab
    
##########################
# FIT
##########################
        def fit(fit_range=None,init=False): # '' is reset to a 'alml' default before create_model(self.model)
                                       #  fit() produces a different layout
                                       #  if self.model_components exists and coincides with 'mgmgml' values are not erased
            '''
            input fit_range is a list of three integers, start stop, pack
            fit tab of mudash, assumes 
               an instance of suite has been run, taking care of
                  t0 = 0
                  forw and backw groups and alpha when needed
                  self.suite.offset, an integer 
               model name and version label  
            sets 
            - fit and plot ranges
            - dashboard according to 
                - single/multi run
                - single/multi group
                - global/sequential          
            - displays: model name
            - fit, plot buttons
            - edit components of model
            - to select parameters value, fix, function, fft subtract etc
            - to produce the dashboard json file that is the input of mufit

            ::

            # the calculation is performed in independent class mumodel, in library mucomponents
            # the methods are "inherited" by mugui 
            # via the reference instance self._the_model_, initialized in steps: 
            #     __init__ share initial attributes (constants) 
            #     _available_components_ automagical list of mucomponents 
            #     clear_asymmetry: includes reset check when suite is implemented
            #     create_model: lay out self._the_model_
            #     delete_model: for a clean start
            #     functions use eval, evil but needed, checked by muvalid, 
            #     int2_int to pass instances of mumodel methods
            #     iminuit2 requires them to be formatted by int2min as 
            #              numpy arrays fitval, fiterr, fitfix, fitlim, stored here as
            #              fitarg={"fitval":fitval,'"fiterr":fiterr, "fitfix":fitfix, "fitlim":fitlim}
            #     help  
            #     load 
            #     save_fit/load_ft save results in mujpy format (dill)
            #     write_csv produces a qtiplot/origin loadable summary
            # 
            # Six fit types: 
            #       single A1, multi single A2, 
            #       suite B1, multi suite B2, 
            #       suite C1, multi suite C2.
            # Multi optimizes a single chi2 over multiple groups of detectors
            # Single is the standard fit of one time-differential asymmetry, A1
            # Multi single optimizes a single chi2 over multiple groups, A2.
            # Suite either iterates a single fit over several runs, B1,
            #       or optimizes a single chi2 over the same runs, C1.
            # Multi suite either iterates a single fit over several runs with multiple groups, B2,
            #             or optimizes a single chi2 over the same runs, C2.
            # with global (one distinct value for all asymmetries, all runs and all groups), 
            #      glocal (one distinct value per group for all runs)
            #      or rlocal (one distinct value for all groups and each run) a parameters
            # another distinction is between 'internal' parameters from mumodel components
            # and 'external' parameters, both those extracted from each run header, such as T,B,theta
            # and those defined by the model, 
            # such as e.g. total asymmetry, fractions, T=0 order parameters, Tc etc
            '''


            from mujpy.mucomponents.mucomponents import mumodel
            import numpy as np
            
            def load_fit(b):
                '''
                loads fit values such that the same fit can be reproduced on the same data
                '''
                from mujpy.aux.aux import path_file_dialog
                import json
                if dropdwn_empty_tmplate.value==2:            
                    dashboard_file = path_file_dialog(self.suite.__templatepath__,'json') # returns the full path and filename
                else:
                    dashboard_file = path_file_dialog(self.suite.__fitpath__,'json')
                if dashboard_file:
                    self.dashboard_file = dashboard_file
                    load_dashboard()
                    

            def load_dashboard():
                '''
                from self.dashboard_file, already verifies
                '''
                import json
                try:        
                    with open(self.dashboard_file,"r") as f:
                        fit_json = f.read()
                        fit_dict = json.loads(fit_json)

                    self.userpar = []  # must be cleared in case neither of userpardicts_results of _guess exist
                    self.version = fit_dict['version']
                    fit_range = json.loads('['+fit_dict['fit_range']+']') # is a list
                    # self.log('debug mudash load_dashboard: fit_range = {}'.format(fit_range))
                    if 'userpardicts_result' in fit_dict:
                        self.userpar = fit_dict['userpardicts_result']
                    elif 'userpardicts_guess' in fit_dict:
                        self.userpar = fit_dict['userpardicts_guess']
                    if 'model_result' in fit_dict:
                        self.model_components = fit_dict['model_result']
                    else:                             # 'model_guess' always exists
                        self.model_components = fit_dict['model_guess']
                    if 'grp_calib' in fit_dict:
                        self.suite.groups = fit_dict['grp_calib']
                        self.suite.store_groups()
                    self.model = ''.join([component["name"] for component in self.model_components])
                    # self.log('debug mudash load_dashboard: self.model = {}'.format(self.model))
                    fit(fit_range=fit_range) # re-initialize the tab with a new model, self.model_components exists so gui values are taken from there
                    # self.log(' group( mudash load_dashboard: loaded from {}'.format(self.dashboard_file))
                    
                except Exception as e:
                    self.erlog('Problems with reading {} file\nException: {}'.format(self.dashboard_file,e))

            def on_fit_request(b):
                '''
                This is the entry to iminuit and fitplot 
                triggered by Fit Button
                retrieve data from the gui dashboard 
                copy latest groupings from the group tab to self.suite
                and produce the json dashboard file            
                '''
                from mujpy.aux.aux import json_name
                import json
                from numpy import array
                from mujpy.mufit import mufit
                from mujpy.mufitplot import mufitplot
    #            from subnotebook import Return, default_values
                
                # copy groupings from group tab
                # group_index = k where k is the index of groups
                # no need for observe, the input values for group k are already stored as 
                # groups[3*k]   0 3 6 ... :  forward  (np.array of indices, as requested of self.suite.grouping)
                # groups[3*k+1] 1 4 7 ... :  backward (np.array of indices, as requested of self.suite.grouping)
                # groups[3*k+2] 2 5 8 ... :  alpha    (float)
                
                initialize_only = False
                try:
                    if b==False:
                        initialize_only = True                        
                except:
                    initialize_only = False
                    
                #self.log('debug mudash on_fit_request\nplot {}'.format(initialize_only))
                group_index = -1 
                # determine if number of groups was changed
                rows = self.mainwindow.children[1].children # the granddaughter
                group_number = sum([len(groups.children[0].children)//3 for groups in rows])
                if group_number > len(self.suite.grouping):
                    for grp in range(len(self.suite.grouping),group_number):
                        self.suite.grouping.append({}) 
                elif group_number < len(self.suite.grouping):
                    for grp in range(group_number,len(self.suite.grouping)):
                        self.suite.grouping.pop()
                for k,row in enumerate(rows): 
                    for j in range(0,len(row.children),3):  # row contains once or twice the following widgets: forward Text, backward Text, alpha FloatText
                        group_index += 1
                        forward = array([int(s)-1 for s in row.children[0].children[j].value.split(',')])  
                        backward = array([int(s)-1 for s in row.children[0].children[j+1].value.split(',')])  
                        alpha = float(row.children[0].children[j+2].value)                   
                        self.suite.grouping[group_index]['forward'] = forward # convert '1,2,3' into array([1,2,3])
                        self.suite.grouping[group_index]['backward'] = backward
                        self.suite.grouping[group_index]['alpha'] = alpha
                        self.suite.groups[group_index]['alpha'] = alpha
                if self.userpar:
                    # fill self.userpar from list_text_userparvalue[k], list_text_userparstd.value
                    for k,par in enumerate(self.userpar):
                        par['value'] = float(ist_text_userparvalue[k].value)
                        par['flag'] = '~'
                        par['error'] = float(list_text_userparstd[k].value)
                        lim0 = None if list_text_userparlim0[k].value == 'None' else float(list_text_userparlim0[k].value)
                        lim1 = None if list_text_userparlim1[k].value == 'None' else float(list_text_userparlim1[k].value)
                        par['limits'] = [lim0, lim1]
                model_guess = []
                nint = len(self.userpar)-1
                for k,comp in enumerate(self.model_components):  # scan the model, already has component "name" and a basic "pardicts"
                # comp is a dictionary for the k-th component of the model, according to addcomponent
                # here we add keys to it so that it can generate the json dashboard
                # build dash["model_guess"], a list of dictionaries with keys "name", "label", "pardicts", for the json file
                    comp["label"] = list_componentlabel[k].value                
                    pardicts = []
                    for j,pardict in enumerate(self.model_components[k]['pardicts']): # update and complete the pardict for each parameter 
                        nint += 1
                        pardict["value"] = float(list_parvalue[nint].value) # is a float
                        pardict["flag"] = list_flag[nint].value
                        fun = list_function[nint].value.split(";") # already checked for sanity
                        if len(fun)>1:
                            pardict["function_multi"] = fun
                            pardict["function"] = ""
                        else:
                            pardict["function"] = fun[0]
                            pardict.pop("function_multi",None)
                        # pardict["error"] is already set in dict_model_component
                        # deactivate error_propagate_multi, that was devised in the hand-written json dashboard 
                        # Error propagation is dealt with numerically by jax grad, as in mujpy.aux.aux errorpropagate, at the best fit
                        pardicts.append(pardict) 
                    comp["pardicts"] = pardicts # list of dictionaries
                    model_guess.append(comp)  
                dash = {}
                dash['version'] = self.version
                dash["fit_range"] = text_fit_range.value # is a string
                dash["offset"] = self.suite.offset
                if self.userpar:
                    dash["userpardicts_guess"] = self.userpar
                dash["model_guess"] = model_guess
#                self.log('  on_fit_request: self.model = {}'.format(self.model))                            
                # dashboard must be stored in a json file
                self.dashboard_file = self.suite.__templatepath__+json_name(self.model,self.suite.datafile,self.suite.grouping,self.version,g=self.userpar)
                with open(self.dashboard_file,"w") as f:
                    json.dump(dash,f,indent=4)
                    
                # initialize_only = True means only plot, False means do the fit first
                #                   plot_type.value = 1 Fit, =2 Guess
                if initialize_only and plot_type.value==1 and self.the_fit: # fit exists, just plot it fit 
                    the_fitplot = mufitplot(text_plot_range.value,self.the_fit,fig_fit=self.fig_fit)
                elif initialize_only and plot_type.value==2: # plot guess
                    the_fit = mufit(self.suite,self.dashboard_file,initialize_only=initialize_only,dash=self)
                    the_fitplot = mufitplot(text_plot_range.value,the_fit,fig_fit=self.fig_fit)
                else: # do fit
                    self.the_fit = mufit(self.suite,self.dashboard_file,dash=self)
                    the_fitplot = mufitplot(text_plot_range.value,self.the_fit,fig_fit=self.fig_fit)                
                self.fig_fit = the_fitplot.fig
                self.fig_fft = the_fitplot.fig_fft
                
            def on_flag_changed(change):
                '''
                observe response of fit tab widgets:
                set disabled on corresponding function (True if flag=='!' or '~', False if flag=='=') 
                '''
                dscr = change['owner'].description # description is three chars ('val','fun','flg') followed by an integer nint
                                                   # iterable in range(ntot), total number of internal parameters
                n = int(dscr[4:]) # description='func'+str(nint), skip 'func'
                list_function[n].disabled=False if change['new']=='=' else True


            def on_function_changed(change):
                '''
                observe response of fit tab widgets:
                check for validity of function syntax

                '''
                from mujpy.aux.aux import muvalid

                dscr = change['owner'].description # description is three chars ('val','fun','flg') followed by an integer nint
                                                   # iterable in range(ntot), total number of internal parameters
                n = int(dscr[4:]) # description='func'+str(nint), skip 'func'
                funct = change['new'].split(";") # in case it's a function_multi (global multi group fit)
                for k,fun in enumerate(funct):
                    error_message = muvalid(change['new'])
                    if error_message:
                        list_function[n].value = ''
                        self.log(error_message+fun+' passed in '+change['new']) 

            def on_parvalue_changed(change):
                '''
                observe response of fit tab widgets:
                check for validity of function syntax
                '''
                # self.log('**** after first parvalue change ****') 
                dscr = change['owner'].description # description is three chars ('val','fun','flg') followed by an integer nint
                                                   # iterable in range(ntot), total number of internal parameters
                n = int(dscr[5:]) # description='value'+str(nint), skip 'func'
                try:
                    float(list_parvalue[n].value)
                    list_parvalue[n].background_color = "white"
                except:
                    list_parvalue[n].value = '0.0' 
                    list_parvalue[n].background_color = "mistyrose"

            def on_offset(change):
                '''
                store new version
                '''
                self.suite.offset = change['owner'].value
                
            def on_plot_request(b):
                '''
                plot wrapper
                '''
                on_fit_request(False)
                
            def on_range(change):
                '''
                observe response of FIT, PLOT range widgets:
                check for validity of function syntax
                '''
                from mujpy.aux.aux import derange
                import time

                fit_or_plot = change['owner'].description[0] # description is a long sentence starting with 'fit range' or 'plot range'
                if fit_or_plot=='f':
                    name = 'fit'
                else:
                    name = 'plot'
                x_range,errmsg = derange(change['owner'].value,self.suite.histoLength) 
                # print('fit or plot = {}, range = {},{},{}, errmsg = {}'.format(name,x_range[0],x_range[1],x_range[2],errmsg))
                if errmsg: #  contains and error message
                                       # go back to defaults
                    if name == 'fit':
                        fit_range = [0,self.suite.histoLength,1] # default
                        string = ','.join([str(val) for val in fit_range])         
                        text_fit_range.background_color = "mistyrose"
                        time.sleep(3)
                        text_fit_range.value = string
                    else:
                        self.plot_range = [0,self.suite.histoLength,1] # default
                        string = ','.join([str(val) for val in plot_range])         
                        text_plot_range.value = string
                        text_plot_range.background_color = "mistyrose"
                else:  # no errors
                    self.plot_range = x_range
                    string = ','.join([str(val) for val in x_range])         
                    if name == 'fit':
                        text_fit_range.background_color = "white"
                        text_fit_range.value = string        
                    else:
                        text_plot_range.background_color = "white"
                        text_plot_range.value=string
       
            def on_add_userpar(b):
                '''
                remove the clicked user parameter
                '''
                import tkinter as tk
                from json import loads as str2lst
                def create_userpar():
                    self.userpar.append({'name':e1.get()})
                    lst = tk.grid_slaves()
                    for l in lst:
                        l.destroy()
                    fit_range = str2lst('['+text_fit_range.value+']')
                    fit(fit_range=fit_range)  # keeps existing text_fit_range.value
                master = tk.Tk()
                master.title('mujpy user parameter')
                tk.Label(master, text="Enter new parameter name").grid(row=0)
                e1 = tk.Entry(master)
                e1.grid(row=0, column=1)
                tk.Button(master, 
                          text='Cancel', 
                          command=master.quit).grid(row=1, 
                                        column=0, 
                                        sticky=tk.W, 
                                        pady=4)
                tk.Button(master, 
                          text='Create', command=create_userpar).grid(row=1, 
                                                           column=1, 
                                                           sticky=tk.W, 
                                                           pady=4)            
                tk.mainloop()
                            
                
            def on_remove_userpar(change):
                '''
                remove the clicked user parameter
                '''
                import tkinter as tk
                from tkinter.messagebox import askyesno
                from json import loads as str2lst
                k = -change['owner'].value
                # are you sure?
                answer = askyesno(message='Delete parameter {}?'.format(str(k)+' - '+list_userparname[k]))
                if answer:
                    self.userpar.pop(k)
                    list_userparname.pop(k) 
                    list_text_userparvalue.pop(k) 
                    list_text_userparstd.pop(k)
                    list_text_userparlim0.pop(k)
                    list_text_userparlim1.pop(k)
                #version = change.value
                fit_range = str2lst('['+text_fit_range.value+']')
                fit(fit_range=fit_range)  # keeps existing text_fit_range.value
                
            def on_add_component(b):
                '''
                add a component
                '''
                from tkinter import Tk, StringVar, OptionMenu, Button, Label
                from json import loads as str2lst
                # Create object
                def add():
                    text = clicked.get()
                    components.append(text)
                    label.config(text=components)
                    # Dropdown menu options
                def save():
                    add()
                    if components:
                        self.model = self.model+''.join(components)
                        self.model_components = [] # initialized, will contain model_guess or model_results 
                        # model_handle.value=self.model
                        self.mainwindow.set_title(0, 'fit '+self.model)
                        fit_range = str2lst('['+text_fit_range.value+']')
                        fit(fit_range=fit_range) # keeps existing text_fit_range.value
                        quit()
                def quit():
                    mucomponent_chooser.destroy()
                    # Dropdown menu options
                components = []
                mucomponent_chooser = Tk()
                mucomponent_chooser.geometry( "200x100" )
            
                options = [ c["name"] for c in list_available_components ]
      
                # datatype of menu text
                clicked = StringVar()
                # initial menu text
                clicked.set( "ml" )
                # Create Dropdown menu
                drop = OptionMenu( mucomponent_chooser, clicked , *options )
                drop.grid(row=0)
                # Create button, it will change label text
                button = Button( mucomponent_chooser, text = "Add" , command = add ).grid(row=1, column=0)
                button = Button( mucomponent_chooser, text = "Quit" , command = quit ).grid(row=1, column=2)
                button = Button( mucomponent_chooser, text = "Add&Save" , command = save ).grid(row=1, column=1)
                # Create Label
                label = Label( mucomponent_chooser, text = "Chosen: " )
                label.grid(row=2,column=1)
      
                # Execute tkinter
                mucomponent_chooser.mainloop()    
            
            
            def on_remove_component(b):
                '''
                remove this component from the model
                need to identify which
                '''
                self.log('Not yet. Wait to see layout.')
                # select component to delete
                # confirm deletion
                # delete component from self.model and from self.model_components 
                # rerun gui cell without cancelling values of surviving components
                fit_range = str2lst('['+text_fit_range.value+']')
                # fit(fit_range=fit_range) # keeps existing text_fit_range.value
                    
                
            def on_version(change):
                '''
                store new version
                '''
                self.version = change['owner'].value
            def on_model(change):
                '''
                check syntax
                store new name
                araba phoenix rebirth of gui 
                '''
                from mujpy.aux.aux import validmodel, create_model
                from json import loads as str2lst
                model = (change['owner'].value).strip() # removes accidental lead & trail blanks
                # self.log('{} is a new model'.format(model))
                if model: # when model is set to '', below on_model is alerted again but nothing happens
                    if validmodel(model):
                        # self.log('{} is a valid model'.format(model))
                        self.model = model
                        self.model_components,emsg = create_model(model)
                        self.mainwindow.set_title(0, 'fit '+self.model)
                        fit_range = str2lst('['+text_fit_range.value+']')
                        fit(fit_range=fit_range) # keeps existing text_fit_range.value
                        model_handle.value = ''
                        model_handle.background_color = "white"
                    else:
                        self.log('{} is not a valid model'.format(model))
                        model_handle.value = ''
                        model_handle.background_color = "mistyrose"
                    
       
######### here starts the fit method of mudash ####################################################
            from ipywidgets import Text, IntText, Layout, Button, HBox, FloatText, \
                                   Checkbox, VBox, Dropdown, ToggleButton, Label
            from mujpy.aux.aux import _available_components_
            from mujpy.aux.aux import path_file_dialog, create_model, addcomponent, name_of_model
            import os
            
            list_available_components =_available_components_() # creates list automagically from mucomponents

            # list of just mucomponents method name, used just in addcomponent

            #------------------------- oneframe

            # step 1) use fit() and self.model_components 
            # if self.dashboard_file is an existing model/*.json file or fit/*_fit.json file load it, must work only first time
            #      make sure that when reloading models self.dashboard_file is cleared
            #   two options, from model/*.json (default) or from fit/*_fit.json if tick besides load button
            #   creates eempty or filled self.model_components and goes to step 2)
            # if init is False skip following if True then
            #   if self.model_components exists -> check that agrees with self.model and go directly to step 2) (use it) 
            #                    includes update (move fit results to self.model_components)
            #                    add and remove component (must edit coherently self.model and self.model_components consistently)
            #                    if self.model_components does not exist or empty
            #   create empty self.model_components with create_model(self.model)
            #   no self.model, fall back to 'alml'
            # 
            # 
            # empty model, editing models (add, remove components), loading templates, loading fits (perhaps a tick?)
            # check that the following does actually restart properly the gui in all instances: 

            if init: # either load dash or fall back to defaults
                if self.model_components:  # this is top priority, fit() recalls gui generated mods on model 
                    # self.log('debug self.model_components in locals()')
                    if not name_of_model(self.model_components,self.model): # if model and model_components do not agree create an empty self.model
                        model_components, emsg = create_model(self.model) # try 
                        if model_components: # is not empty
                            self.model_components = model_components
                        else: # del.model was not valid
                            self.log(emsg) # print error message from the attempt
                            self.log('Tring to create model {}. Sorry, try again, fall back to model alml'.format(self.model))
                            self.model='alml' # fall back
                            self.model_components, emsg = create_model(self.model) # produces self.model_components (ignores labels)
                    else: # self.model does match self.model_components
                        # self.log('self.model_components matches self.model')
                        pass # self.model_components is ok
                        
                elif os.path.isfile(self.dashboard_file): # self.model_components does not exist, startup
                    load_dashboard()   # self.dashboard_file exists as such 
                elif os.path.isfile(os.path.join(self.suite.__templatepath__,self.dashboard_file)): #  or in ./model/
                    self.dashboard_file = os.path.join(self.suite.__templatepath__,self.dashboard_file)
                    load_dashboard()
                elif os.path.isfile(os.path.join(self.suite.__fitpath__,self.dashboard_file)): # or in ./fit/
                    self.dashboard_file = os.path.join(self.suite.__fitpath__,self.dashboard_file)
                    load_dashboard()
                else: # no self.dashboard_file, create an empty model from self.model, passed at startup
                    self.model_components, emsg = create_model(self.model) # produces self.model_components (ignores labels)
                  
            textheight = '23px'
            labelheight = '23px'
            buttonheight = '48px'
            bwidth = '92%'

            fit_button = Button(description='Fit',
                                 tooltip='Execute the Minuit fit',
                                 layout=Layout(width='6%',height=buttonheight))              # 6%
            fit_button.style.button_color = self.button_color
            fit_button.on_click(on_fit_request)

            offset_width = '6%'
            offset = IntText(value=self.suite.offset,
                             layout=Layout(width=bwidth,height=textheight),
                             description_tooltip='First good bin\n(number of bins skipped\nafter prompt peak)',
                             continuous_update=False) # offset, is an integer
                                                                                              # 12%
            offset_description = Label('  offset',layout=Layout(width=bwidth,height=labelheight))
            # initialized to 50, only input is from an IntText, integer value, or saved and reloaded from mujpy_setup.pkl
            offset.observe(on_offset)

            default_fit_range = [0,self.suite.histoLength,1]
            # fit_range may be None if not initialized, a list of integers otherwise
            # convert to string or check and convert
            #self.log('mudash fit debug: fit_range = {}'.format(fit_range))  
            if not fit_range: # first call not initialized go to default
                fit_range = default_fit_range
                string = ','.join([str(val) for val in fit_range])
            else: # a list of integers, first call, fit_range passed from __init__, check its sanity
                if all([isinstance(x,int) for x in fit_range]) and fit_range[2]<=(fit_range[1]-fit_range[0]):
                    string = ','.join([str(val) for val in fit_range]) # must be converted to string
                else: # wrong syntax, go to defaults
                    self.log('Error in range: {}. Back to default {}'.format(
                      ','.join([str(val) for val in fit_range]),','.join([str(val) for val in default_fit_range])))
                    fit_range = default_fit_range
                    string = ','.join([str(val) for val in fit_range])
            # self.log('mudash fit debug: string = {}'.format(string))          
            fit_range_width = '11%'              
            text_fit_range = Text(value=string,description='f',
                                  layout=Layout(width=bwidth,height=textheight),
                                  description_tooltip='start,stop[,pack]\ne.g. 0,20000,10\ni.e. starts from offset\nuse 20000 bins\n pack them 10 by  10',
                                  continuous_update=False)                                  # 23%
            text_fit_range.style.description_width='0%'
            fit_range_description = Label('     fit range',layout=Layout(width=bwidth,height=labelheight))
            text_fit_range.observe(on_range,'value')

            plot_button = Button (description='Plot',
                                  tooltip='Generate a plot',
                                  layout=Layout(width='6%',height=buttonheight))             # 29%
            plot_button.style.button_color = self.button_color
            plot_button.on_click(on_plot_request)

            plot_range_width = '18%'
            string = ','.join([str(val) for val in self.plot_range])
            text_plot_range = Text(value=string,
                                   description='p',
                                   description_tooltip='start,stop[,pack][,last,pack]\ne.g 0,20000,10 \nor 0,2000,10,20000,100 \npack 10 up to bin 2000\npack 100 from bin 2001 to bin 20000',
                                   layout=Layout(width=bwidth,height=textheight),
                                   continuous_update=False)                        # 47%
            text_plot_range.style.description_width='0%'
            plot_range_description = Label('      plot range',layout=Layout(width=bwidth,height=labelheight))
            text_plot_range.observe(on_range,'value')
            self.mainwindow.set_title(0, 'fit '+self.model)
            
            type_width = '11%'
            plot_type = Dropdown(options=[('Fit',1),
                                                    ('Guess',2)],
                                                    value=1,
                                                    description_tooltip='Plot best fit\n(if available)\nor initial guess',
                                                    layout=Layout(width=bwidth,height=textheight)) 
                                                                                             # 58%
            type_description = Label('fit/guess',layout=Layout(width=bwidth,height=labelheight),
                                                    description_tooltip='Plot best fit\n(if available)\nor initial guess')

            addcomponent_button = Button(description='+\ncomponent',
                                         layout=Layout(width='8%',height=buttonheight),
                                         tooltip='Add new component')              # 66%
            addcomponent_button.on_click(on_add_component)
            addcomponent_button.style.button_color = self.button_color
            
            model_width = '8%'
            model_handle = Text(value = ' ',
                                layout=Layout(width=bwidth,height=textheight),
                                description_tooltip='Acronyim for new model\n from scratch',
                                continuous_update=False)                          #  74%
            model_handle.observe(on_model,'value')
            model_description = Label('new model',layout=Layout(width=bwidth,height=labelheight))
            
            version_width = '9%'
            version = Text(value=self.version,
                           description_tooltip='String to distinguish among model output files',
                           layout=Layout(width=bwidth,height=textheight))   # 83%
            version.observe(on_version,'value')
            version_description = Label('version',layout=Layout(width=bwidth,height=labelheight))
                
            loadbutton = Button(description='Load',
                                tooltip='Opens GUI to load\none of the existing\nfit templates',
                                layout=Layout(width='6%',height=buttonheight))      #  89%
            loadbutton.style.button_color = self.button_color
            loadbutton.on_click(load_fit)

            #        update_button = Button (description='Update',
            #                                tooltip='Update parameter starting guess\nfrom latest fit\n(must have fitted this model once).',
            #                                layout=Layout(width='7.3%'))                     # %
            #        update_button.style.button_color = self.button_color
            #        update_button.on_click(on_update)

            template_width = '11%'
            dropdwn_empty_tmplate = Dropdown(options=[('Best Fit',1),
                                                    ('Template',2)],
                                                    value=1,
                                                    description_tooltip='Load either\ngeneric template or\nBest fit for given run',
                                                    layout=Layout(width=bwidth,height=textheight))      # 100%
            template_description = Label('model',layout=Layout(width=bwidth,height=labelheight),description_tooltip='Load either\ngeneric template or\nBest fit for given run')
       
    #-------------------- fill model template into input widgets, two columns



            nint = -1 # internal parameter count, each widget its unique name

            list_parname, list_parvalue, list_flag, list_function = [],[],[],[]            

            # define userpar in a Jupyter notebook, e.g. as
            # usepar = [{'name':'B[mT]'},{'name':'A21'},{'name':'A34'},{'name':'φ21},{'name':'φ34'}]
            # before invoking
            # mudash(usepar=usepar)

            if self.userpar:
                list_userparname, list_text_userparvalue, list_text_userparstd, list_text_userparlim0, list_text_userparlim1, list_button_userparm = [] , [], [], [], [] # lists of handles, index runs according to internal parameter count nint
                
                nuserpar = len(self.userpar)
                s_n,s_nam,s_val,s_std,s_lim,s_add ='P[#]','Name','Value','Std','Limits','+'+str(nuserpar+1)
                userhead = HBox([Label(s_n,layout={'width':'3%','height':'16pt'},description_tooltip='number to be used in Functions.'),
                                Label(s_nam,layout={'width':'6%','height':'16pt'}),
                                Label(s_val,layout={'width':'9%','height':'16pt'},description_tooltip='initial guess.'),
                                Label(s_std,layout={'width':'4%','height':'16pt'},description_tooltip='initial step'),
                                Label(s_lim,layout={'width':'8%','height':'16pt'},description_tooltip='e.g. 0, None, None is no limit'),
                                Button(s_add,layout={'width':'3%'})
                                ])
                userhead[5].style.button_color = self.button_color
                userhead[5].on_click(on_add_userpar)
                one_list, two_list, three_list = [],[],[]
                
                nuserpar = len(self.userpar)
                for k in range(nuserpar):
                    nint += 1      # all parameters are internal parameters, first is pythonically zero 
                        
                    nintlabel_handle = Text(value=str(nint),
                                            layout=Layout(width='3%'),
                                            disabled=True)                               # 3%
                    value = self.userpar[k]['name'] # this one must exist, checked at _init__
                    list_userparname.append( 
                                        Text(value=value,
                                          layout=Layout(width='6%'),
                                          disabled=True))                                 # 9%
                    # self.log('\n{} - comp {} par {} appended'.format(nint,self.model_components[k]['name'],name))

                    value = '' if not 'value' in self.userpar[k] else self.userpar[k]['value']
                    list_text_userparvalue.append(
                                      Text(value=value,
                                      layout=Layout(width='9%'),
                                      continuous_update=False))                          # 18%
                    # parvalue handle must be unique and stored at position nint, it will provide the initial guess for the fit

                    value = '' if not 'std' in self.userpar[k] else self.userpar[k]['std']
                    list_text_userparstd.append( 
                                     Text(value=value,
                                     layout=Layout(width='4%'),
                                      continuous_update=False))                          # 22%

                    value = 'None' if not 'limits' in self.userpar[k] else self.userpar[k]['limits'][0]
                    list_text_userparlim0.append( 
                                     Text(value=value,
                                     layout=Layout(width='4%'),
                                      continuous_update=False))                          # 26%

                    value = 'None' if not 'limits' in self.userpar[k] else self.userpar[k]['limits'][2]
                    list_text_userparlim0.append( 
                                     FloatText(value=value,
                                     layout=Layout(width='4%'),
                                      continuous_update=False))                          # 30%
                    value = -k 
                    list_button_userparm.append( 
                                     Button(value=value,
                                     layout=Layout(width='2%'),
                                     tooltip='Click to remove this parameter'))                          # 32%
                    list_button_userparm[k].on_click(on_remove_userpar)
                    list_button_userparm[k].style.button_color = self.button_color


                    par_handle = HBox([nintlabel_handle, list_userparname[nint], list_text_userparvalue[nint], list_text_userparstd[nint]],layout=Layout(width='100%'))
                    if k%3==0: 
                        if k//3==0:
                            one_list.append(userhead) 
                        one_list.append(par_handle) # append it to the left if k even
                    elif k%3==1: 
                        if k//4==0:
                            two_list.append(userhead) 
                        two_list.append(par_handle) # or to the right if k odd  
                    elif k%3==2: 
                        if k//4==0:
                            three_list.append(userhead) 
                        three_list.append(par_handle) # or to the right if k odd  

                if nuserpar//3: # 3 or more userpar, all lists are populated
                    userframe_handle = HBox([one_list,two_list,three_list],layout={'width':width})
                elif nuserpar==2: # only 2 userpar
                    userframe_handle = HBox([one_list,two_list],layout={'width':width/4})
                elif nuserpar==1: # Only 1 userpar
                    userframe_handle = HBox(one_list,layout={'width':width/4})

            #-------------------- create the model template
            # nint = len(self.userpar)-1
            for k in range(len(self.userpar)): # if len(self.userpar)>0 keeps parameter k -> list_parname[k] etc.
                list_parname.append([])
                list_parvalue.append([])
                list_flag.append([])
                list_function.append([])
                # append an empy list for each so that nint is correct the number also for these handles
                 
            s_n,s_nam,s_val,s_flag,s_func,s_plot ='p[#]','Name','Value','Flag','Function','Plot'
            add_component_button = Button(layout={"width":'8%','height':'16pt'},description="+",description_tooltip='adds a component to the model')
            add_component_button.style.button_color = self.button_color
            add_component_button.on_click(on_add_component)
            dashheadplus = HBox([
                            Label(s_n,layout={'width':'8%','height':'16pt'},description_tooltip='Number to be used in Functions.'),
                            Label(s_nam,layout={'width':'20%','height':'16pt'}),
                            Label(s_val,layout={'width':'20%','height':'16pt'},description_tooltip='initial guess.'),
                            Label(s_func,layout={'width':'44%','height':'16pt'},description_tooltip='~ is free\n! is fixed,\n= activates function.\np[0] replaces by par 0 (use only previous numbers).\nSimple algebra is allowed, e.g.\n0.5*p[0]+0.5*p[4].')
                            ])
            dashhead  = HBox([Label(s_n,layout={'width':'16%','height':'16pt'},description_tooltip='Number to be used in Functions.'),
                            Label(s_nam,layout={'width':'20%','height':'16pt'}),
                            Label(s_val,layout={'width':'20%','height':'16pt'},description_tooltip='initial guess.'),
                            Label(s_func,layout={'width':'44%','height':'16pt'},description_tooltip='~ is free\n! is fixed,\n= activates function.\np[0] replaces by par 0 (use only previous numbers).\nSimple algebra is allowed, e.g.\n0.5*p[0]+0.5*p[4].')
                            ])
            leftframe_list, rightframe_list = [],[]
      
            list_fftcheck,list_componentlabel,list_remove_component = [],[],[]        
            nleft,nright = 0,0
            # self.log('n of components {}'.format(len(self.model_components)))
            for j in range(len(self.model_components)):  # scan the model
                list_fftcheck.append(Checkbox(description='FFT',
                                              description_tooltip='uncheck for showing this component\nin the FFT of residues',
                                              value=True,
                                              layout=Layout(width='16%')))              
                list_componentlabel.append(Text(description='tag',
                                              description_tooltip='remove minuit parameter name degeneracy\nuse descriptive tags\ndefault is component number',
                                              value=str(j),
                                              layout=Layout(width='30%')))              
                list_componentlabel[j].style.description_width='20%'
                list_remove_component.append(Button(layout={"width":'10%'},description="del",tooltip='Remove this component from the model'))
                list_remove_component[j].style.button_color = self.button_color
                list_remove_component[j].on_click(on_remove_component)
            
                list_fftcheck[j].style.description_width='2%'
                header = HBox([ Text(value=self.model_components[j]['name'],
                                     disabled=True,
                                     layout=Layout(width='8%')), 
                                list_componentlabel[j],                      #  8 %
                                list_fftcheck[j],
                                list_remove_component[j]]) # component header HBox              # 24 %
                                                   # composed of the name (e.g. 'da'), the label, the FFT flag, the trash button
                                                   # fft will be applied to a 'residue' where only checked components
                                                   # are subtracted
                                                   
                if j%2==0:                         # and ...
                    leftframe_list.append(header)  # append it to the left if k even
                    if j==0:
                          leftframe_list.append(dashheadplus)
                else:
                    rightframe_list.append(header) # or to the right if k odd 
                    if j==1:
                          rightframe_list.append(dashhead)
      
                                                   # list of HBoxes, headers and pars
                nleftorright = 0                                   
                for k in range(len(self.model_components[j]['pardicts'])): # make a new par for each parameter 
                                                                                    # and append it to component_frame_content
                    nint += 1      # all parameters are internal parameters, first is pythonically zero 
                    nleftorright += 1
                    nintlabel_handle = Text(value=str(nint),
                                            layout=Layout(width='7%'),
                                            disabled=True)                               # 7%
                    name = self.model_components[j]['pardicts'][k]['name']
                    baloon = ''
                    rates = ['Δ','Γ','λ','σ']
                    fields = ['B','BL']
                    if name in fields:
                        baloon = '[mT]'
                    elif name in rates:
                        baloon = '[mus-1]'
                    elif 'φ' in name:
                        baloon = '[deg]'
                    elif 'ν' in name:
                        baloon = '[MHz]'
                    list_parname.append( 
                                        Text(value=name,
                                          description_tooltip=baloon,
                                          layout=Layout(width='16%'),
                                          disabled=True))                                 # 23%
                    # self.log('\n{} - comp {} par {} appended'.format(nint,self.model_components[k]['name'],name))

                    list_parvalue.append(
                                      Text(value='{:.4f}'.format(self.model_components[j]['pardicts'][k]['value']),
                                      layout=Layout(width='15%'),
                                      description='value'+str(nint),
                                      continuous_update=False))                          # 38%
                    list_parvalue[nint].style.description_width='0%' # hidden, used in on_parvalue_changed
                    try:
                        list_parvalue[nint].value = _parvalue[nint]
                    except:
                        pass
                    # parvalue handle must be unique and stored at position nint, it will provide the initial guess for the fit

                    list_flag.append(Dropdown(options=['~','!','='], 
                                     value=self.model_components[j]['pardicts'][k]['flag'],
                                     layout=Layout(width='11%'),
                                     description_tooltip='e.g.\n=p[2]/2\nfor single group\n=0.5*p[2];=p[3]/2\nfor a global two-group fit',
                                     description='flag'+str(nint)))                      # 49%
                    list_flag[nint].style.description_width='0%'
                    try:
                        list_flag[nint].value = _flag[nint]
                    except:
                        pass
                     # flag handle must be unique and stored at position nint, it will provide (eventually) the nonlinear relation to be evaluated

                     
                    # function works also for multi group global fits, which require different definitions for different groups
                    # the separator between two group definition is ';' so that the syntax for two groups is e.g.
                    # =p[2];=p[3]    
                                 
                    list_function.append(Text(value=self.model_components[j]['pardicts'][k]['function'],
                                         layout=Layout(width='36%'),
                                         description='func'+str(nint),
                                         description_tooltip='e.g.\n=p[2]/2\nfor single group\n=0.5*p[2];=p[3]/2\nfor a global two-group fit',
                                         continuous_update=False))                       # 85%

                    # function handle must be unique and stored at position nint, it will provide (eventually) the nonlinear relation 

                    fdis = False if self.model_components[j]['pardicts'][k]['flag']=='=' else True 
                    list_function[nint].disabled = fdis # enabled only if flag='='
                    list_function[nint].style.description_width='0%'
                    # now put this set of parameter widgets for the new parameter inside a parameter HBox
                    par_handle = HBox([nintlabel_handle, list_parname[nint], list_parvalue[nint], list_flag[nint], list_function[nint]],layout=Layout(width='100%'))
                               # handle to an HBox of a list of handles; notice that parvalue, flag and function are lists of handles
                    
                    # now make value flag and function active 
                    list_parvalue[nint].observe(on_parvalue_changed,'value')
                    list_flag[nint].observe(on_flag_changed,'value') # when flag[nint] is modified, function[nint] is z(de)activated
                    list_function[nint].observe(on_function_changed,'value') # when function[nint] is modified, it is validated
     
                    if j%2==0:                                         # and ...
                        leftframe_list.append(par_handle) # append it to the left if k even
                        nleft += nleftorright
                    else:
                        rightframe_list.append(par_handle) # or to the right if k odd  
                        nright += nleftorright 
            

            width = '99.8%'
            widthhalf = '100%'
            leftframe_handle = VBox(layout=Layout(width=widthhalf),
                                    children=leftframe_list)#,layout=Layout(width='100%')
            rightframe_handle = VBox(layout=Layout(width=widthhalf),
                                     children=rightframe_list)# ,layout=Layout(width='100%')

    #------------------ include frames in boxes        
            oneframe_handle = HBox([fit_button,
                                    VBox([offset,offset_description],
                                         layout=Layout(width=offset_width)),
                                    VBox([text_fit_range,fit_range_description],
                                         layout=Layout(width=fit_range_width)),
                                    plot_button,
                                    VBox([text_plot_range,plot_range_description],
                                         layout=Layout(width=plot_range_width)),
                                    VBox([plot_type,type_description],
                                         layout=Layout(width=type_width)),
                                    addcomponent_button,
                                    VBox([version,version_description],
                                         layout=Layout(width=version_width)),
                                    VBox([model_handle,model_description],
                                         layout=Layout(width=model_width)),
                                    loadbutton,
                                    VBox([dropdwn_empty_tmplate,template_description],
                                         layout=Layout(width=template_width))],
                                    layout=Layout(width=width)
                                    )

            bottomframe_handle = HBox([leftframe_handle,rightframe_handle],
                                       layout=Layout(width=width,border='1px solid dodgerblue')
                                      )

            # now collect the handles of the three horizontal frames to the main fit window 
            if self.userpar:
                fittab = [VBox([oneframe_handle,
                                userframe_handle,
                                bottomframe_handle],layout=Layout(width=width))]
            else: 
                fittab = [VBox([oneframe_handle,
                                bottomframe_handle],layout=Layout(width=width))]

            self.mainwindow.children[0].children = fittab

 # HERE CONTINUES dash_init   
        # now popolate self.mainwindow tabs    
        titlelogowindow = gui() # define top gui
        
        fit(fit_range=fit_range,init=True) # this is the fit tab, fit_range is passed from __init__
        self.mainwindow.children[0].layout = Layout(border = self.border_color,width='100%')        
        
        grouptab = group()          #   
        self.mainwindow.children[1].children = grouptab

        abouttab = about() # fine as is (fine tune later)
        self.mainwindow.children[2].children = abouttab
        self.mainwindow.children[2].layout = Layout(border = self.border_color,width='100%')
        # Araba.Phoenix: gui starts as method and it redefines itself as a widget
        gui_whole = VBox([titlelogowindow, self.mainwindow],
                         description='dash',
                         layout=Layout(width='auto'))
        dash = HBox([gui_whole,html_style,board],layout=Layout(border = '2px',width='auto'))

        display(dash) # prints dash and html-styled board, side by side

