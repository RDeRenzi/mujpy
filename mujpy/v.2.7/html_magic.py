from os import getcwd
startuppath = getcwd()
from mujpy.musuite import suite
from mujpy.mudashed import dashed

#get_ipython().run_cell_magic('html', '',                              '<style> \n.lbl_bg{\n    width:auto;\n    background-color: yellow;}\n.box_style{\n    width:40%;\n    border : 2px solid red;\n    height: auto;\n    background-color:black;\n}\n                            </style>\n\nfrom ipywidgets import widgets\n\n\nlbl  = widgets.Label(value=  \'Test\' )\n# lbl  = widgets.HTMLMath(value=  \'Test\' ) # Alternate way using HTMLMath\nlbl.add_class(\'lbl_bg\')\nhBox = widgets.HBox([lbl],layout=widgets.Layout(justify_content= \'flex-end\'))\nhBox.add_class("box_style")\ndisplay(hBox)\n')

datafile = '/home/roberto.derenzi/musrfit/MBT/gps/run_05_21/data/deltat_tdc_gps_0822.bin'
runlist = '822' # first run first
offset = '20'
grp_calib = [{'forward':'3', 'backward':'4', 'alpha':1.13},{'forward':'2', 'backward':'1', 'alpha':1.13}]
#
the_suite = suite(datafile, runlist , grp_calib , offset, startuppath)
the_dash = dashed(the_suite)

