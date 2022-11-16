'''
Version history
Versions 1.x, not documented, based on mujpy.mugui, enclusing both gui and functions
              includes libraries aux.aux and aux.plot, utilities
                       class musr2py and muisis2py for musr data input
                       class mucomponents for fit functions
                       class muedge and muprompt for t=0 determination at isis and psi bulkmusr  
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
Version  2.7, introduces fit C1 (multirun singlegroup global) working form the notebook
         2.7.1 introduces officially versioning (this ur-class, to be inherited by all)
'''
__version__ = "2.7.1"
