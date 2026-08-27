'''
Version history
version  3.0    new skimmed version of class musuite, mufit, mufitplot, mucomponents, 
                libraries tools.tools tools.plot 
                new Jupyter class mudashed.py, ipywidgets interface 
                with Output boxes for log and matplotlib animations
__________________________________________________________________________
version  2.9   rewritten mufit and mufitplot to support all eight fit types
               single run           A1, A20, A21 (single group, multi group sequential, multigroup global)
               multirun sequential  B1, B20, B21 (>>          , >>                    , >>)
               multirun global      C1, C2       (global single group, global multigroup)  
               each with its alpha calibration version
Version  2.8   new classes muroot2py.muroot2py muroot2py based on uproot
               and muisis2py.muisis2py munxs2py based on nexusformat.
               debugged A1 A1_calib A20 A20_calib A21 A21_calib B1 modes
               new example/test.py, individual python3 command scripts
               and models in example/fit/....1_fit.json
               mudash is not adequately implemented yet
________________________________________________________________________               
--- 11/2024 transferred to  nex pc34
Version  2.7.3.3a mistake in nexus, is this the last?
Version  2.7.3.2 FIxin another silly mistake in docs for the distro (arghhh)
Version  2.7.3.1 FIxin silly mistake in docs for the distro
Version  2.7.3 mudash works decently for non global, 
                                         global single run multi group and 
                                         global multirun single group fits
Version  2.7.2 introduces grad in minuit, analytical gradients for mg for testing purposes.
               gradient works for all components but Kuto Toyabe only in single group multirun global fit
               works but has two drawbacks: it is slower and it requires a guess much closer to minimum
               From this version this cdocstring is reversed (top is newest) 
Version  2.7, introduces fit C1 (multirun singlegroup global) working form the notebook
         2.7.1 introduces officially versioning (this ur-class, to be inherited by all)
Version  2.0, separates the following functions
                       class musuite data input
                       class mufit performs the fit, based on a specific json input file
                       class mufitplot produces the fit plots
                       class mudash is the gui that produces the json file, runs the fit etc.
                       up to here version is not officially defined
                       but subsequent versions included fits A1, A20 (sequential multigroup), 
                       then fit A21 (global multigroup),
                       then fits A1_calib. A2_0_calib, A21_calib (include alpha fit),
                       then fits B1, B20 (sequential multirun and multigroup)
                       then fit B21 (sequential multirun, global multigroup)
                       initially run by a notebook and later by the gui
Versions 1.x, not documented, based on mujpy.mugui, enclusing both gui and functions
              includes libraries aux.aux and aux.plot, utilities
                       class musr2py and muisis2py for musr data input
                       class mucomponents for fit functions
                       class muedge and muprompt for t=0 determination at isis and psi bulkmusr  

'''
__version__ = "2.8"
