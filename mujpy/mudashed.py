class dashed(object):
    '''
    barebone ipywidgets fit editor for jupiter nb, produces dashboard.json files

    a subset of a full fledged mudash GUI interface
        changing MN (model) or NG (number of globals) require rerun
        cannot delete components or single global parameters
    '''


##########################
# INIT
##########################
    def __init__(self,facility='PSI',sidecar=False,test=None):
        '''
            Launches the simple gui, that requires only an instance of mujpy.musuite suite 
        '''
 
        from mujpy._version import __version__
        from mujpy.tools.tools import make_links

        self.__version__ = __version__
        self.facility = facility
        self.sidecar = sidecar
        self.mudashed_width = '900px'   
        self.output_width = '900px'
        self.textheight = '23px'
        self.labelheight = '23px'
        self.buttonheight = '48px'

        self.suite_button_color = '#d6cbbf'
        self.command_button_color = '#97b3ae'
        self.global_button_color = '#c29c94'
        self.model_button_color = '#c0b1ab'
    

        # initialize dashboard, a dictionary
        # self.log = self.suite.console
        self.root = None # initialize for tkinter
        self.fig_fit = None
        self.fig_fft = None
        self.loading_dash = False # ensures widgets are built for empty model, unless loaded from file
        self.test, self.data_dir, writeable_folder  = make_links(test) # links groups, plus data and fit if test
        if writeable_folder and self.data_dir: # normal and test mode
            self.board()
        elif not Self.data_dir: # test attempt on a folder with a permanent data dir
            print('please start test or demos from an empty folder, e.g. tmp')
        else: # start outside %HOME
            print('please start from a writeable folder')

    def log(self,string):
        """
        written output
        """

        with self.board_box: # just an Output
            # THIS MUST BE print, NOT self.log!
            print(string) 

    def command2dash(self):
        """
        activates second row self.command_1 of mudashed command_box

            Button Fit 8%
            Label FR 5%
            Text FR range 12%
            Label VR 5%
            Text VR tag 8%
            Button Plot 8%
            Dropdown Fit/Guess 10%
            Label PR 5%
            Text PR range 16%
            Label RF 5%
            FloatText RF frequency 8%
            Button FFT 8%
            Label nuR 5%
            Text nuR range 10%
            Label LB 5%
            FloatText LB 8%
            Tot 25% + 24% + 48% = 97% 
        """
        # begins below def build_dashed
        def on_Fit(b):
            # build_dashed()
            # the_fit = mufit(self.suite,dashboard_file)
            # with self.figure_box:
            #     the_plot = mufitplot(plot_range,the_fit)

            import json
            from mujpy.mufit import mufit
            from mujpy.mufitplot import mufitplot
            if self.build_dashed(): # creates self.dashboard and returns True if no validation raise occurred
                dashboard_file = self.suite.__fitpath__+'dashed.json'
                with open(dashboard_file,'w') as f:
                    json.dump(self.dashboard,f) # mufit wants to read this from a file
                self.board_box.clear_output()
                the_fit = mufit(self.suite,dashboard_file,dash_log = self.log) # writes text to board_box
                #self.figure_box.clear_output()
                plot_range = self.command_1.children[8].value
                rotfreq = self.command_1.children[9].value
                fft_range = self.command_1.children[12].value # not yet in use 
                lb = self.command_1.children[13].value # not yet in use 
                self.tab.selected_index = 2
                the_plot = mufitplot(plot_range, the_fit, rotating_frame_frequencyMHz = rotfreq, plot_out = self.figure_box, fig_fit = self.fig_fit) # plots in self.figure_box
                self.fig_fit = the_plot.fig
                

        def on_Plot(b):
            """
            invokes mufitplot from mudashed

            """
            
            # read dashed widget values (including guess, rotfreq
            # json.dump again
            # self.figure_box.clear_output()
            # if guess: mufit(plot_range,dashboard_file, no_fit = not guess,out=self.figure_box)
            # mufitplot(plot_range, guess = guess, rotating_frame_frequencyMHz = rotfreq, plot_out = self.figure_box)
            import json
            from mujpy.mufit import mufit
            from mujpy.mufitplot import mufitplot
            if self.build_dashed(): # creates self.dashboard and returns True if no validation raise occurred
                dashboard_file = self.suite.__fitpath__+'dashed.json'
                self.log('dumping {}'.format(dashboard_file))
                with open(dashboard_file,'w') as f:
                    json.dump(self.dashboard,f) # mufit wants to read this from a file
                guess = self.command_1.children[6].value=='Guess'
                self.the_fit = mufit(self.suite,dashboard_file,no_fit = guess, dash_log = self.log) # writes text to board_box
                #self.figure_box.clear_output()
                plot_range = self.command_1.children[8].value
                rotfreq = self.command_1.children[9].value
                the_plot = mufitplot(plot_range, 
                                     self.the_fit, 
                                     rotating_frame_frequencyMHz = rotfreq, 
                                     plot_out = self.figure_box, 
                                     fig_fit = self.fig_fit) # plots in self.figure_box
                self.fig_fit = the_plot.fig


        def on_FFT():
            # later
            self.log('Nothing yet!')

        from ipywidgets.widgets import Layout, Button, Label, Text, Dropdown, FloatText, HTML
        command1_width = ['10%','5%','10%','5%','6%','12%','10%','5%','16%','10%','8%','5%','8%','10%']
        fit_range = self.dashboard['fit_range'] if 'dashboard' in self.__dir__() else '0,20000,40' if self.suite._the_facility_ == 'PSI' else '0,1000,1'
        version = self.dashboard['version'] if 'dashboard' in self.__dir__() else '1'

        buttonfit = Button(description='Fit',layout=Layout(width=command1_width[0]))
        buttonfit.on_click(on_Fit)
        buttonfit.style.button_color = self.command_button_color
        buttonplot = Button(description='Plot',layout=Layout(width=command1_width[5]))
        buttonplot.on_click(on_Plot)
        buttonplot.style.button_color = self.command_button_color
        buttonfft = Button(description='FFT',layout=Layout(width=command1_width[10]))
        buttonfft.on_click(on_FFT)
        buttonfft.style.button_color = self.command_button_color
        RF_float = FloatText(value=0,
                             description = 'RF', 
                             tooltip='Rot Frame\nFreq [MHz]',
                             layout=Layout(width=command1_width[9]))
        RF_float.style.description_width = '25%'
        LB_float = FloatText(value='0.3',
                             description='LB',
                             tooltip='μs-1\nfft filter',
                             layout=Layout(width=command1_width[13]))
        LB_float.style.description_width = '25%'
        widgets = [buttonfit,
                   Label(value='FR',layout=Layout(width=command1_width[1])),
                   Text(value=fit_range,tooltip='fit start,stop,pack',layout=Layout(width=command1_width[2])),
                   Label(value='VS',layout=Layout(width=command1_width[3])),
                   Text(value=version,tooltip='version label',layout=Layout(width=command1_width[4])),
                   buttonplot,
                   Dropdown(options=['Fit','Guess'],value='Fit',layout=Layout(width=command1_width[6])),
                   Label(value='PR',layout=Layout(width=command1_width[7])),
                   Text(value=fit_range,tooltip='srt,stp,pack\nsrt,stp0,pack0,stp,pack',layout=Layout(width=command1_width[8])),
                   RF_float,
                   buttonfft,
                   Label(value='νR',layout=Layout(width=command1_width[11])),
                   Text(value='0,50.',tooltip='start, stop\n[MHz]',layout=Layout(width=command1_width[12])),
                   LB_float 
                    ]
#        with self.board_box:
#            self.log('command2dash self.command_1.children {}'.format(widgets))
 
        self.command_1.children = widgets
           
    def json2dash(self):
        '''
        builds second stage widgets & value from self.dashboard - LD and LF - or from NG_int.value and self.MN_text_value

        assumes self.dashboard is loaded before calling json2dash and valid (check!)
        adds rows beyond first one, according to json content of file
            ['model_guess'] and ['globpardicts_guess']
        '''

        from ipywidgets.widgets import HBox, VBox, Label, HTML, Text, Layout
        from mujpy.tools.tools import glob2widgets, comp2widgets, par2widgets, par2labels 
        from mujpy.tools.tools import _available_components_, add_step_limits
        
        hashed = False
        hspacer = Label(' ',layout={'width':'42%','height':'16pt'})
        glotitle = Text(value='global parameters',disabled=True,layout={'width':'14%','height':'16pt'})
        gbt_style = "<style>.gbt_input input { background-color:#FFDAB9 !important; }</style>"
        glotitle.add_class('gbt_input')
        modtitle = Text(value='model parameters',disabled=True,layout={'width':'14%','height':'16pt'}) 
        mod_style = "<style>.mod_input input { background-color:#DDC1B0 !important; }</style>"
        modtitle.add_class('mod_input')
        model_title = HBox([hspacer,HTML(mod_style),modtitle,hspacer])
        if 'dashboard' in self.__dir__(): # 
            if 'globpardicts_guess' in self.dashboard:
                pardicts = self.dashboard['globpardicts_guess']
                self.command_0.children[0].value='global fit'
                self.NG_int.value = len(pardicts)
                hashed = '#' in [pardict['flag'] for pardict in pardicts]
        else:
            if self.NG_int.value > 0:
                pardicts = []
                for kp in range(self.NG_int.value):
                    pardicts.append({'name':'','value':0.0,'flag':'~','error':0.001,'limits':'None,None','positive_parity':False})

        # skip this if fit is sequential
        glob = False
        if self.NG_int.value > 0: # 'globpardicts_guess' in self.dashboard: 
            glob = True
            global_title = HBox([hspacer,HTML(gbt_style),glotitle,hspacer])
            flags = ['~','!','#'] if hashed or 'dashboard' not in self.__dir__() else ['~','!']
            # unwraps unique pardict keys into labels
            keys = list(max(pardicts,key=len).keys()) # list of all unique pardicts keys
            keys.insert(0,'p[k]') # order is 'p[k]' 'name' 'value' 'flag' 'error' 'limits' 'par>0'
            if 'positive_parity' in keys:
                keys[keys.index('positive_parity')]='par>0'
                # p[k] name value flag error limits par>0
            keylen = ['5%','12%','18%','12%','16%','20%','9%']
            labels = [Label(value=key,layout=Layout(width=klen)) for key,klen in zip(keys,keylen)]
            
            n_columns = self.NG_int.value//2+self.NG_int.value%2
            
            # now add one row of widgets per pardict
            left_column, right_column = [HBox(labels)],[HBox(labels)] # first row of labels
            for kp,pardict in enumerate(pardicts[:n_columns]):
                if 'positive_parity' not in pardict: pardict['positive_parity']=False
                left_column.append(glob2widgets(kp,pardict,flags,keylen))
            for kp,pardict in enumerate(pardicts[n_columns:]):
                if 'positive_parity' not in pardict: pardict['positive_parity']=False
                right_column.append(glob2widgets(kp+n_columns,pardict,flags,keylen))
            self.global_box.children = [global_title,HBox([VBox(left_column),VBox(right_column)])]
            # global lists k name value flag error limits pospar
            # always add model, = self.dashboard or from _available_components_, + label + value + error + limits  
        if 'dashboard' in self.__dir__():
            model = self.dashboard['model_guess'] # list of components
            model_name = ''.join([component['name'] for component in model])
            self.MN_text.value = model_name 
        # MN_text calls on_MN but  empty global_box and model_bosonly if self.dashboard does not exist  
        else:
            components = [self.MN_text.value[i:i + 2] for i in range(0, len(self.MN_text.value), 2)]
            model = [av for component in components for av in _available_components_() if component == av['name']] 
            for kc,compdict in enumerate(model): # kc is component index and npar its number of pardicts
                # add label
                compdict['label'] = str(kc)
                for pardict in compdict['pardicts']:
                    # add value
                   pardict['value'] = 0.0
                   pardict['flag'] = '=' if glob else '~'
        if self.NG_int.value==0:
            model = add_step_limits(model) # adds errors, limits, pospar 

        # now transfer model from either sources to the same dashed
        kp = self.NG_int.value
        left_column, right_column = [],[]
        labels = par2labels(model[0]['pardicts'][0]) # returna an HBox of labels for par2widges
            #for kc in range(0,len(model),2):
        for kc,compdict in enumerate(model):
            # unpack component [k]
            column = right_column if kc%2 else left_column
            column.append(comp2widgets(compdict,kc)) # returns an Hbox of component widgets
            for j,pardict in enumerate(compdict['pardicts']):
                if j==0: # first parameter, self.log labels
                    column.append(labels)
                kp += 1
                #self.log('json2dash pardict["flag"] = {}'.format(pardict['flag']))
                column.append(par2widgets(pardict,kp-1,glob=glob)) # returns am HBox of parameter widgets

        layout_column = Layout(width='50%')
        columns = [model_title,HBox([VBox(left_column,layout=layout_column), VBox(right_column,layout=layout_column)])]

        self.model_box.children = columns

    def build_dashed(self):
        """
        reads dashed widget values and builds self.dashboard

            validates
            - values with errors, limits, invalid_err_lim
            uses read_pardict_from_widgets which also validates
            - function math syntax, muvalid
            - 0 <= get_indices < kmax = NG_int.value for global fits and mudashed index of this parameter
        """

        from mujpy.tools.tools import limits, _available_components_, read_pardict_from_widgets
        from mujpy.tools.tools import invalid_err_lim, add_step_limits

        self.dashboard = {}
        kids = self.command_1.children
        self.dashboard['version'] = kids[4].value # (string)
        self.dashboard['fit_range'] = kids[2].value # (string)
        self.dashboard['offset'] = str(self.suite.offset)
        glob = self.global_box.children # empty list is false
        if glob: 
            pardicts = []
            #                    kid[0]             kid[1]     kidd0  kiddd0...n   Kidd1 kiddd0...m
            # global_box = VBox([HBox([hsp,gt,hsp]),HBox([VBox([HBox(leftparwidgs),HBox(rightparwidgs)])])])
            n_col_left, n_col_right = len(self.global_box.children[1].children[0].children),len(self.global_box.children[1].children[1].children)
            kmax = n_col_left+n_col_right-2 # number of parameters aka NG_int.value (n_col includes labels)
            # needed for ki, i.e. for read_pardict_from_widgets
            for k,parwidgs in enumerate(self.global_box.children[1].children[0].children): # left column
                pardict = {}
                if k: # skips k=0 labels
                    invalid = invalid_err_lim(parwidgs.children[2].value,
                                              parwidgs.children[4].value,
                                              limits(parwidgs.children[5].value)) 
                    if invalid: 
                        raise ValueError('global parameter str(k-1):'+invalid)
                    pardict['name'] = parwidgs.children[1].value # string
                    pardict['value'] = parwidgs.children[2].value # float
                    # for sequentila fits may be a list, i.e. a string as w2idgets value
                    pardict['flag'] = parwidgs.children[3].value # string
                    pardict['error'] = parwidgs.children[4].value # float
                    pardict['limits'] = limits(parwidgs.children[5].value) # translates list of csv string to values
                    if parwidgs.children[6].value: pardict['positive_parity'] = parwidgs.children[6].value
                    pardicts.append(pardict)
            for k,parwidgs in enumerate(self.global_box.children[1].children[1].children): # right column
                pardict = {}
                if k: # skips k=0 labels
                    invalid = invalid_err_lim(parwidgs.children[2].value,
                                              parwidgs.children[4].value,
                                              limits(parwidgs.children[5].value)) 
                    if invalid: 
                        raise ValueError('global parameter str(k-1):'+invalid)
                    pardict['name'] = parwidgs.children[1].value # string
                    pardict['value'] = parwidgs.children[2].value # float
                    # for sequentila fits may be a list, i.e. a string as w2idgets value
                    pardict['flag'] = parwidgs.children[3].value # string
                    pardict['error'] = parwidgs.children[4].value # float
                    pardict['limits'] = limits(parwidgs.children[5].value) # translates list of csv string to values
                    if parwidgs.children[6].value: pardict['positive_parity'] = parwidgs.children[6].value
                    pardicts.append(pardict)
            self.dashboard['globpardicts_guess'] = pardicts
        model = []
        components = [self.MN_text.value[i:i + 2] for i in range(0, len(self.MN_text.value), 2)]
        avc = _available_components_()
        nparam = [len(av['pardicts']) for component in components for av in avc if component == av['name']] 
        left, right = self.model_box.children[1].children[0], self.model_box.children[1].children[1] 
        # list of two VBox, to be read, left and right columns
        rowleft, rowright, ki = 0, 0, 0

        #self.log('debug build_dashed start pardicts')
        for kc,npar in enumerate(nparam): # ks is component index and npar its number of parameters
            # skip (component+legend) and read npar rows
            # writes    name,  None, flag, [function] or function_multi if glob
            #           name, value, flag, [function] otherwise
            compdict = {'name':components[kc],'label':str(kc)}
            pardicts = [] # list of pardict
            if kc%2: # odd, right
                # see json2dash, component widgets:
                rowright += 2 # component and legend rows  
                for kp in range(npar):
                    kmax = kmax if glob else ki
                    pardict = read_pardict_from_widgets(right.children[rowright],kmax) # from widgets
                    if not isinstance(pardict,dict): self.log('debug right build_dash pardict {}'.format(pardict)) # is an erro message from read_pardict_from_widgets
                    pardicts.append(pardict) # returns a pardict
                    rowright += 1 # incremented only for odd kc
                    ki +=1 # internal dashed parameter index incremented always
            else: # even, 0, 2, ... left
                rowleft += 2 # component and legend rows
                for kp in range(npar):
                    kmax = kmax if glob else ki
                    pardict = read_pardict_from_widgets(left.children[rowleft],kmax)
                    if not isinstance(pardict,dict): self.log('debug left build_dash pardict {}'.format(pardict)) # is an erro message from read_pardict_from_widgets
                    pardicts.append(pardict) # returns a pardict
                    rowleft += 1 # incremented only for even kc
                    ki +=1 # internal dashed parameter index incremented always
            compdict['pardicts'] = pardicts
            model.append(compdict)
        if not glob:
            self.dashboard['model_guess'] = add_step_limits(model) # adds errors, limits, pospar 
        else:
            self.dashboard['model_guess'] = model
        return True
        #except ValueError as e:
        #    self.log('build_dashed ValueError: {}'.format(e))
        #    return False


    def board(self):
        '''
        entry, draws the gui editor in 3 stages, suite, model selection, editor

            for suite input and information
            for command actions (Fit,Plot,FFT,Ranges ...)
            for 'model_guess' list of components and their pardicts
            for 'globpardicts_guess' list of global parameters
        '''

        from ipywidgets.widgets import Output, ToggleButtons, Button, Label, Layout, Text, IntText
        from ipywidgets.widgets import Dropdown, FloatText, HBox, VBox, HTML, Box, Image, Textarea
        from ipywidgets.widgets import Tab
        from mujpy.musuite import suite
        from mujpy import __file__ as MuJPyName
        from mujpy._version import __version_tuple__ as version_tup
        from mujpy.tools.tools import _available_components_
        from datetime import datetime
        from sidecar import Sidecar
        import os

        ##################################################################################################
        # Use from scratch
        # Three stages and four boxes, suite, command_, global_ and method_ box
        # first suite box, musuite is instantiated by a run list change [press Enter]
        # second first row of command box,
        # ToggleButtons select sequential/global fit 
        # IntText number of global parameter(NG) [default 0, disabled=True], update before 
        #   - ToggleButton.on_click toggles NG=1, disable=False
        # MN_text select model acronym [press Enter], LL load last fit/dashed.json, LF selects fit/*.json  
        #   - MN, LL, LF triggers third where
        #     self.command2dash() self.json2dash() build model_box, global_box children for all three cases
        # [widgets with tooltips, some actions preceeded by two letter label]
        ###################################################################################################

        def on_fit_type(change):
            '''
            toggle NG disabled False/True and set self.NG_int.value = 1,0

            '''

            if fit_type.value == 'sequential fit':
                self.NG_int.value = 0
                self.NG_int.disabled = True
            else:
                #with self.board_box:
                #    self.log('fit_type global fit')
                self.NG_int.value = 1
                self.NG_int.disabled = False

        def on_LL(change):
            '''
            Load Last dashed.json if it exists
            '''

            import json
            import os

            self.loading_dash = True # to allow json2dash to write MN_text
            file_json = self.suite.__fitpath__+'dashed.json'
            if os.path.isfile(file_json):
                with open(file_json,'r') as f:
                    self.dashboard = json.load(f) # copies json dict to self.dashboard
                self.command2dash()
                self.json2dash() # builds widgets for this model
            else:
                self.log('>>>>>>>>>>>>>>>> file dashed.json not found')
            self.loading_dash = False # builds widgets for this model
 
        def on_LF(change):
            '''
            Choose fit model from ./fit/ folder to load
            '''

            from mujpy.tools.tools import path_file_dialog 
            import json
            import os

            self.loading_dash = True # to allow json2dash to write MN_text
            file_json, self.root = path_file_dialog(self.suite.__fitpath__,'json',root = self.root)
            # self.log('Trying to load {} ...'.format(file_json))
            if os.path.isfile(file_json):
                if file_json[-4:]=='json':
                    with open(file_json,'r') as f:
                        self.dashboard = json.load(f) # copies json dict to self.dashboard
                    self.command2dash()
                    self.json2dash() # builds widgets for this model
            else:
                self.log('no valid json file was selected {}'.format(file_json))
            self.loading_dash = False # to allow dash to change MN_text
    
        def on_MN(change):
            '''
            if valid model, adds widgets to regenerate the second stage of mudashed (complete editor) 
            '''
    
            # command box, global box: VBox of rows (HBox)
            # model_box: HBox of two columns (VBox) of components (VBox) of rows, pardicts (HBox) of component widgets
            from mujpy.tools.tools import validmodel 
            from json import loads as str2lst
            from ipywidgets import Text, IntText, Layout, Button, HBox,  \
                                   VBox, ToggleButtons, Label, FloatText
            from tkinter.messagebox import askyesno, showerror
            #self.log('change["new"] is {}'.format(change['new']))

            if not self.loading_dash: # to allow json2dash to set model_name without interference 

                if self.NG_int.value==1: # suspicious!
                    no = not askyesno(title='Check!', message="NG=1 global parameter\ndo you need more?")
        
                model = self.MN_text.value.strip() # removes accidental lead & trail blanks
        
                no = True               
                if model and no: # when model is set to '', below on_model is alerted again but nothing happens

                    if validmodel(model):
#                        with self.board_box:
#                            self.log('MN {} is valid'.format(self.MN_text.value))                   # self.log('{} is a valid model'.format(model))
                        # model_box stricly needs to be self since on_MN has no return
                        self.command2dash()
                        self.json2dash()  # in this case self.dashboard is not loaded from a json file 
                                            # an model variable equivalent to self.dashboard is created 
                                            # model is used to build the second stage gui
                        #with self.board_box:
                            #self.log('self.command_1.children {}'.format(self.command_1.children))
                            #self.log('self.global_box.children {}'.format(self.global_box.children))
 
                        #display(panels)
                    else:
                        showerror(title='Wrong model syntax',message='{} not made of valid components!'.format(model))
                        self.MN_text.value = ''
                else:
                    self.MN_text.value = ''

        def on_RL(change):
            """
            start suite from run list input, checks path file exists

            beware: as of ipywidgets v. 8.1.5 this continuous_update=False is a bit of a mess
                    change['new'] is initially a dict instead of a value
                    and Enter triggers a double call.
                    proceeds only id change['new'] is not a dict
            """

            import os
            from mujpy.musuite import suite
            from mujpy.tools.tools import derun, get_title, get_gtotals, get_grouping, group_syntax
            from numpy import all


            runlist  = change['new']
            if not isinstance(runlist,dict):
                datafile = self.suite_box.children[2].children[2].value  # path Text value
                #self.log('runlist = {}'.format(runlist))

                try:
                    grp = self.suite_box.children[1].children[1].value
                    forward, backward = grp.split('-')
                    if all(get_grouping(forward)>=0) and all(get_grouping(backward)>=0):
                        grp_calib = [{'forward':forward, 
                                  'backward':backward, 
                                  'alpha':float(self.suite_box.children[1].children[2].value)}]
                    else:
                        raise NameError('No Group0')
                except ValueError as e:
                    f,b = get_grouping(forward), get_grouping(backward)
                    if isinstance(f,str):
                        e = f
                        if isinstance(b,str): e += ';'+b
                    elif isinstance(b,str): e = b
                    #self.log('Exception {}'.format(e))
                    #self.log('group syntax error: {}'.format(grp))
                    text = 'Exception {}'.format(e)
                    text += '\nGroup0 syntax error: {}'.format(grp)
                    self.root = group_syntax(text,root=self.root)
                    #self.log('Group0 group_syntax return a self.root = {}'.format(self.root))
                    return
                grp = self.suite_box.children[1].children[4].value
                if grp:
                    try:
                        alph = self.suite_box.children[1].children[5].value 
                        groups,alphas = grp.split(';'),alph
                        for group,alpha in zip(groups,alphas):
                            forward1, backward1 = group.split('-')
                            if all(get_grouping(forward1)>=0) and all(get_grouping(backward1)>=0):
                                grp_calib.append({'forward':forward1, 
                                      'backward':backward1, 
                                       'alpha':float(alpha)})
                    except ValueError as e:
                        f,b = get_grouping(forward), get_grouping(backward)
                        if isinstance(f,str):
                            e = f
                            if isinstance(b,str): e += ';'+b
                        elif isinstance(b,str): e = b
                        #self.log('Exception {}'.format(e))
                        #self.log('other group syntax error: {}, alpha {}'.format(grp,alph))
                        #self.tab.selected_index = 2
                        text = 'Exception {}'.format(e)
                        text += '\nGroups ... syntax error: {}'.format(grp)
                        self.root = group_syntax(text, root=self.root)
                        #self.log('Groups ... group_syntax return a self.root = {}'.format(self.root))
                        return

                offset = self.suite_box.children[1].children[8].value
                if os.path.isfile(datafile):
                    if runlist:
                        self.suite = suite(datafile , runlist , grp_calib , offset , 'CettoLaqualunque',console=self.log) #startuppath is set in suite
                        # self.log info
                        if self.suite.loadfirst: # suite loaded the data
                            starttime_options = [' '.join(self.suite._the_runs_[k][0].get_timeStart_vector()) for k in range(self.suite.nruns)]
                            starttime_options.insert(0,'Run start times')
                            self.suite_box.children[0].children[0].options = starttime_options
                            self.suite_box.children[0].children[0].value = starttime_options[1]

                            stoptime_options = [' '.join(self.suite._the_runs_[k][0].get_timeStop_vector()) for k in range(self.suite.nruns)]
                            stoptime_options.insert(0,'Run stop times')
                            #self.log('nruns {}, SD_options = {}'.format(self.suite.nruns,starttime_options))
                            self.suite_box.children[0].children[1].options = stoptime_options
                            self.suite_box.children[0].children[1].value = stoptime_options[1]

                            title_options = [get_title(self.suite._the_runs_[k][0]) for k in range(self.suite.nruns)]
                            title_options.insert(0,'Titles')
                            self.suite_box.children[0].children[2].options = title_options
                            self.suite_box.children[0].children[2].value = title_options[1] 

                            comment_options = [self.suite._the_runs_[k][0].get_comment() for k in range(self.suite.nruns)]
                            comment_options.insert(0,'Comments')
                            self.suite_box.children[0].children[3].options = comment_options 
                            self.suite_box.children[0].children[3].value = comment_options[1]

                            totalcounts, groupcounts, nsbin, maxbin = get_gtotals(self.suite)
                            self.suite_box.children[0].children[4].value = nsbin
                            self.suite_box.children[0].children[5].value = maxbin 
                            goptions = ['Group counts']
                            toptions = ['Total counts']
                            runs,e = derun(runlist)
                            for runadd,groupcount,totalcount in zip(runs,groupcounts,totalcounts):
                                run = ','.join([run for run in runadd])
                                for k,gc in enumerate(groupcount):
                                    counts = ': '+gc
                                    goptions.append(run+'.'+str(k)+counts)
                                toptions.append(run+': '+totalcount[0])
                            
                            self.suite_box.children[1].children[6].options = goptions
                            self.suite_box.children[1].children[7].options = toptions
                            self.command_box.children[0].children = self.command_0.children
                        else:
                            tk_error('No runs luaded, runlist {}?'.format(runlist),'suite error',root=self.root)
                            self.root = None # destroyed by tk_error
                    else:
                        self.log('Please specify runlist')
                else:
                    self.log('File {} not found'.format(datafile))
                    self.log('paths must be either in startup path and below or absolute')
            else:
                self.log('debug, on_RL, change["new"] was still a dict instead of the new text value')

        def on_LG(b):
            """
            group file load, tkinter 

            """

            import os
            from mujpy.tools.tools import path_file_dialog

            startpath = os.getcwd()
            grouppath = startpath+os.path.sep+'groups'+os.path.sep
            if os.path.exists(grouppath):
                groupfile,self.root = path_file_dialog(grouppath,'grp', root=self.root)
            else:
                self.log('Folder {} does not exist'.format(datapath))
            for kg, group in enumerate(groupfile):
                    self.suite_box.children[1].children[1].value = groupfile 

        def on_RL_button(b):
            """
            simulate RL text change
            """
            
            runlist = self.suite_box.children[2].children[5].value
            #self.suite_box.children[2].children[5].value =  ''
            self.suite_box.children[2].children[5].value = runlist 

        def on_DL(b):
            """
            data file load, tkinter 
            """

            import os
            from mujpy.tools.tools import path_file_dialog

            startpath = os.getcwd()
            datapath = startpath+os.path.sep+'data'+os.path.sep
            if os.path.exists(datapath):
                datafile, self.root = path_file_dialog(datapath,'*', root=self.root)
            else:
                self.log('Folder {} does not exist'.format(datapath))
            if datafile:
                self.suite_box.children[2].children[2].value = datafile 

        def on_fetch(change):
            """
            fetch PSI data
            """

            """
            self.fetch_box = VBox([HBox([Dropdown(options=areas,
                                        description='instrument'
                                        value='GPS',
                                        layout=Layout(width='20%')),
                                Dropdown(options=years,
                                        description='year'
                                        value=current_yr,
                                        layout=Layout(width='20%')),
                                IntText(value = 1,
                                        description = 'Start run #',
                                        layout=Layout(width='20%')),
                                IntText(value = 2
                                        description = 'Stop run #',
                                        layout=Layout(width='20%')),
                                        fetch_button         
                                    ]),
                          Textarea(value='',disabled=True,layout=Layout(width='900px',height='160px'))])
            """
            from mujpy.tools.tools import fetch_PSI_data
            from os.path import split
            from os import isdir
            area = self.fetch_box.children[0].children[0].value
            year = self.fetch_box.children[0].children[1].value
            run_start = self.fetch_box.children[0].children[2].value
            run_stop = self.fetch_box.children[0].children[3].value
            datapath = self.suite_box.children[2].children[2].value
            if datapath:
                datapath = datapath if isdir(datapath) else split(datapath)[0]
                error = fetch_PSI_data(year,area,run_start,run_stop,datapath)
                if error:
                    self.fetch_box.children[1].value += 'Error searching PSI database\n {}'.format(error)
                else:
                    self.fetch_box.children[1].value = 'Loaded {} file(s) in data/ path'.format(run_stop-run_start+1)
            else:
                self.fetch_box.children[1].value = 'No data path present, write one in the Fit tab'

                        
                

##########################
# initiate gui first stage
##########################   
        """
        tabs Fit Log Help About
        info [SD] Text^ [PD] Text^ [TT] Text^ [CM] Text^ [NS] Text^ [MB] Text^ 
             [3]  15    [3]  15         20         20         10         10    tot 90
             [GR] Text α FloatText [GR] Text α FloatText
                  16   3  16            16   3  16        tot 70
             [#]  Path Text DLButton Run List Text* OF IntText  
             [8]  5    24   10       10       16    3  10       tot 78

        self.suite_box [^ disabled, Buttons all on_click, * observe]
        """


        info_width = ['16%','16%','27%','27%','6%','8%']
        SD_drop = Dropdown(options = ['Run start times'], value='Run start times',
                        layout = Layout(width=info_width[0]),
                        tooltip = 'run start times',
                           )# this is not disabled = True, but has no observe
        SD_drop.add_class("custom-grey-dropdown") #  and is color gray
        PD_drop = Dropdown(options = ['Run stop times'], value='Run stop times',
                        layout = Layout(width=info_width[1]),
                        tooltip = 'run stop time')
        PD_drop.add_class("custom-grey-dropdown") #  and is color gray
        TL_drop = Dropdown(options = ['Titles'], value='Titles',
                        layout = Layout(width=info_width[2]),
                        tooltip = 'titles')
        TL_drop.add_class("custom-grey-dropdown") #  and is color gray
        CM_drop = Dropdown(options = ['Comments'], value='Comments',
                        layout = Layout(width=info_width[3]),
                        tooltip = 'comments')
        CM_drop.add_class("custom-grey-dropdown") #  and is color gray
        NS_text = Text(value='',
                        layout = Layout(width=info_width[4]),
                        tooltip = 'ns/bin',
                        disabled = True)
        MB_text = Text(value='0',
                          layout=Layout(width=info_width[5]),
                          tooltip='max bins',disabled = True)

        suite_info = HBox([
                       SD_drop,
                       PD_drop,
                       TL_drop,
                       CM_drop,
                       NS_text,
                       MB_text
                       ])

        groups_width=['7%','15%','7%','18%','14%','10%']
        GR_label = [Label('Group 0:',layout=Layout(width=groups_width[0])),
                    Label('Groups ...:',layout=Layout(width=groups_width[0]))]
        GR_text =  [Text(value='3-4',layout=Layout(width=groups_width[1]),tooltip = '3-4\n2,3-4,1\n[fwd-bwd]'),
                    Text(value='',placeholder='option, or LG',layout=Layout(width=groups_width[1]),tooltip = 'fw1-bw1;fw2-bw2')]
        alpha_text = [Text(value='1.0',layout=Layout(width=groups_width[2]),tooltip = 'group α'),
                      Text(value='1.0',layout=Layout(width=groups_width[2]),tooltip = 'group α')]
        GT_dropdown = Dropdown(value = 'Group counts',
                               options = ['Group counts','Run: 0, 0'],
                               #disabled = True,
                               layout = Layout(width=groups_width[3]))
        TO_dropdown = Dropdown(value = 'Total counts',
                               options = ['Total counts','Run: 0'],
                               layout = Layout(width=groups_width[4]))
        value = 20 if self.facility == 'PSI' else 7
        OF_inttext = IntText(value = value,
                             description = 'OF',
                             layout = Layout(width=groups_width[5]),
                             tooltip = 'first good bin')
        OF_inttext.style.description_width = '30%'


        #GT_text = [Text(value='0',layout=Layout(width=groups_width[3]),tooltip='Group total',disabled = True),
        #           Text(value='0',layout=Layout(width=groups_width[7]),tooltip='Group total',disabled = True)]

        suite_groups = HBox([
                            GR_label[0],
                            GR_text[0],
                            alpha_text[0],
                            GR_label[1],
                            GR_text[1],
                            alpha_text[1],
                            GT_dropdown,
                            TO_dropdown,
                            OF_inttext
                            ])
 
        runs_width = ['7%','4%','54%','5%','16%']
        LG_button = Button(description='LG',
                           tooltip = 'Groups file selection',
                           layout = Layout(width=runs_width[0]))
        LG_button.on_click(on_LG)
        LG_button.style.button_color = self.suite_button_color
        DL_button = Button(description='DL',
                           tooltip = 'Data file selection',
                           layout = Layout(width=runs_width[0]))
        DL_button.on_click(on_DL)
        DL_button.style.button_color = self.suite_button_color
        #if self.test:
        #    RL_value = '822' if self.test=='GPS' else '3561' if self.test=='LEM' else '126645'
        #else:
        RL_value = ''
        RL_text = Text(value=RL_value,
                       placeholder='run #s, Enter or RL',
                      tooltip='e.g. 822\nor 822,823:827:-1',
                      layout=Layout(width=runs_width[4]),
                      continuous_update = False)
        RL_text.observe(on_RL,names='value')
        RL_button = Button(description='RL',
                           tooltip = 'Load run list',
                           layout = Layout(width=runs_width[0]))
        RL_button.on_click(on_RL_button)
        RL_button.style.button_color = self.suite_button_color

        path_value = self.data_dir if self.test else ''
        suite_runs = HBox([
                        LG_button,
                        Label(value='path',
                            layout=Layout(width=runs_width[1])),
                        Text(value = path_value,
                             placeholder = 'to proto-data-file (enables DL)', 
                            tooltip = 'path to proto-run/n(from start path\nor absolute',
                            layout = Layout(width=runs_width[2]),
                            continuous_update = False),
                        DL_button,
                        Label(value='run list',
                              layout=Layout(width=runs_width[3])),
                        RL_text,
                        #Label(value='OF',
                        #      layout = Layout(width=runs_width[5]),tooltip = 'first good bin'),
                        RL_button
                        ])
                            
        self.suite_box = VBox([suite_info,suite_groups,suite_runs],layout=Layout(width=self.mudashed_width,border='1.5px solid DarkGoldenrod'))
        
        custom_css = """
        <style>
            .jp-OutputArea-output pre { white-space: pre !important; }
            .container { width:100% !important; }
            /* background and text color for unselected ToggleButtons */
            .widget-toggle-button {
                background-color: #b5d6d0;
                color: #888888 ;
            }

            /* background and text color for SELECTED ToggleButtons (active) */
            .widget-toggle-button.mod-active {
                background-color: #97b3ae;
                color: #000000;
                border-color: #4b5957 ;
            }
            /* SD_ PD_ TL_ CM_ dropdown in nsuite_info */
            .custom-grey-dropdown select option {
                color: grey !important;
            }
        </style>

        """
        
        # 2. Inietta il CSS nel notebook tramite un widget HTML
        css_widget = HTML(value=custom_css)
        command_width = ['38%','21%','11%','14%','8%','8%']
        self.figure_box = Output(layout=Layout(width='100%',height='410px'))# width='900px'
        self.board_box = Output(layout=Layout(width='100%',height='650px',overflow_y='auto'))
        fit_type = ToggleButtons(options = ['sequential fit','global fit'],
                                 value = 'sequential fit',
                                 tooltips = ['A1 A20 B1 B20\nsingle asymmetry fit','A21 B21 C1 C2\nmulti asymmetries fit'],
                                 layout = Layout(width=command_width[0]))
        fit_type.observe(on_fit_type)
        fit_type.style.description_width='0%'

        self.NG_int = IntText(value = 0,description='global parameters',
                        tooltip = 'parameters number',
                        layout = Layout(width=command_width[1],height=self.textheight),
                        disabled = True)
        self.NG_int.style.description_width = '58%'

        MN_label = Label(value='model acronym',layout=Layout(width=command_width[2]))
        self.MN_text = Text(value = '',
                            placeholder = 'xx, Enter to lauch',
                            tooltip = 'e.g. mg\n   almgml',
                            layout = Layout(width=command_width[3],height=self.textheight),
                            continuous_update=False) # requires CR
        self.MN_text.observe(on_MN)

        LL_button = Button(description = 'LL',
                           tooltip = 'Load last model\nif exists',
                           layout = Layout(width=command_width[4]))
        LL_button.on_click(on_LL)
        LL_button.style.button_color = self.command_button_color
        
        LF_button = Button(description = 'LF',
                           tooltip = 'Load fit file',
                           layout = Layout(width=command_width[5]))
        LF_button.on_click(on_LF)
        LF_button.style.button_color = self.command_button_color

        layout = Layout(width=self.mudashed_width,border='1px solid CadetBlue')
        board_width='930px'
        command_0 = HBox([])
        self.command_0 = HBox([fit_type,self.NG_int,MN_label,self.MN_text,LL_button,LF_button],layout={'width':self.mudashed_width})
        self.command_1 = HBox([],layout={'width':self.mudashed_width})
        self.command_box = VBox([command_0,self.command_1],layout=layout)
        #self.command_box.add_class("command_box_style",layout=layout)
        self.global_box = VBox([],layout=Layout(width=self.mudashed_width,border='1px solid Coral'))
        #self.global_box.add_class("global_box_style",layout=layout)
        self.model_box = VBox([],layout=Layout(width=self.mudashed_width,border='1px solid RosyBrown')) # filled with neft and right columns
        #self.model_box.add_class("model_box_style",layout=layout)
        dash = VBox([self.suite_box,
                     self.command_box,
                     self.global_box,
                     self.model_box],
                    layout={'width':'100%','border':self.model_button_color}) # 'width':board_width
        #panels = HBox([dash,VBox([self.figure_box,self.board_box])],layout={'width':'100%'})
        #now = datetime.now()
        #dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        areas = ['LEM','GPS','LTF','VMS','Dolly','GPD','HAL','FLAME']
        current_yr = datetime.today().strftime("%Y")
        years = [str(yr) for yr in range(2003,int(current_yr)+1)]
        fetch_button = Button(description = 'Fetch from PSI',
                                        layout=Layout(width='20%'))
        fetch_button.on_click(on_fetch)
        fetch_button.style.button_color = self.suite_button_color
        self.fetch_box = VBox([HBox([Dropdown(options=areas,
                                       description='instrument',
                                        value='GPS',
                                        layout=Layout(width='20%')),
                                Dropdown(options=years,
                                        description='year',
                                        value=current_yr,
                                        layout=Layout(width='20%')),
                                IntText(value = 1,
                                        description = 'start run',
                                        layout=Layout(width='20%')),
                                IntText(value = 2,
                                        description = 'stop run',
                                        layout=Layout(width='20%')),
                                        fetch_button         
                                    ]),
                          Textarea(value='Here you may fetch PSI data files.\nAlready have them? Click on Fit tab and selected data/ path',disabled=True,layout=Layout(width='900px',height='160px'))])
        help_box = Textarea(
                            disabled=True, # Prevents users from editing the text
                            layout=Layout(width='900px',height='660px'))
        help_text = 'Available components,  by unique two letters, and their Minuit parameters, (x is a time array [μs])'
        for c in _available_components_():
            help_text += '\n{}: {}'.format(c['name'],c['tip'].replace('\n    ',' ',2)) 
        help_box.value = help_text
        logo_file = open(os.path.join(os.path.join(os.path.dirname(MuJPyName),"logo"),"logo.png"), "rb")
        logo_image = logo_file.read()
        logo = Box([Image(value=logo_image)],layout=Layout(width='114px',height='100px'))
        about_text = "mujpy        "+'v'+'.'.join([str(version_tup[k]) for k in range(3)])
        about_text += "\npython μSR data analysis"
        about_text += "\nby R. De Renzi 2017-2026"
        about_text += "\n_________________________________________________"
        about_text += "\ncontributors, direct and indirect"
        about_text += "\nmusr2py: P. Bonfà (wrapper), A. Amato, A. Raselli"
        about_text += "\ndynamical KT: G. Allodi"
        about_text += "\nideas stolen from: A. Suter (musrfit)"
        about_text += "\ncgi-bin fetch and ideas: Z. Salman"
        about_text += "\nroot by uproot (muroot2py wrapper)"
        about_text += "\nnexus by nexusformat /muisis2py wrapper)"

        about = HBox([logo,Textarea(
            value = about_text,
            disabled=True, # Prevents users from editing the text
            layout=Layout(width='786px',height='160px')           #,height='250px' # Height constraint triggers the scrollbar
            )],layout=Layout(width='900px'))       
        if self.sidecar:
            display(css_widget,dash)
            SC = Sidecar(title='Mudashed log: {}'.format(datetime.now().strftime("%d/%m/%Y %H:%M")))
            with SC:
                display(self.board_box)
        else: 
            self.tab = Tab([dash,self.fetch_box,self.board_box,help_box,about],lyout=Layout(width='920px'))
            self.tab.titles = ['Fit','Fetch data','Log','Help','About']
            self.tab.selected_index = 0
            display(css_widget,self.tab)
        # Button( icon = 'fa-trash' #, <i class="fa-thin fa-trash"></i>
        #https://stackoverflow.com/questions/60116974/what-is-the-icon-argument-for-ipywidgets-button
