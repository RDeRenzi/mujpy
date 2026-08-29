# For ipython only. Choose your style: uncomment one or more of the following magics
# %matplotlib tk 
# for TkAgg backend (in ipython it requires tkinter installed system-wide
# or jupyter-lab must be running
# %matplotlib widget # for output in external window without need to P.show(), useful in ipython
from os import getcwd #, chdir
import matplotlib.pyplot as P
startuppath = getcwd() # before changing to git mujpy (if mujpy installed this is not needed)
# chdir('/home/roberto.derenzi/git/mujpy/')
# print(getcwd())
from mujpy.musuite import suite
from mujpy.mufit import mufit
from mujpy.mufitplot import mufitplot
datafile = 'data/deltat_tdc_gps_0822.bin'
runlist = '822' # first run first
offset = '20'
grp_calib = [{'forward':'3', 'backward':'4', 'alpha':1.13},{'forward':'2', 'backward':'1', 'alpha':1.13}]
#
the_suite = suite(datafile, runlist , grp_calib , offset, startuppath)
dashboard_file = startuppath+'/fit/almgml.822.3-4+2-1.1_fit.json' # must exist
the_fit = mufit(the_suite,dashboard_file)#,no_fit=True)
print('>>>>>>>>>>>>>> close (x) the figure to finish')
plot_range = '0,7000,40,25000,200'
the_fitplot = mufitplot(plot_range,the_fit)
P.show()
