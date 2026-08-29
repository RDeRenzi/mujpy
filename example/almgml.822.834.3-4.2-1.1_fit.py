# For ipython only. Choose your style: uncomment one or more of the following magics
# %matplotlib tk 
# for TkAgg backend (in ipython it requires tkinter installed system-wide
# or jupyter-lab must be running
# %matplotlib widget # for output in external window without need to P.show(), useful in ipython
from os import getcwd, chdir
import matplotlib.pyplot as P
startuppath = getcwd() # before changing to git mujpy (if mujpy installed this is not needed)
# chdir('/home/roberto.derenzi/git/mujpy/')
# print(getcwd())
from mujpy.musuite import suite
from mujpy.mufit import mufit
from mujpy.mufitplot import mufitplot
datafile = 'data/deltat_tdc_gps_0822.bin'
runlist = '822,827:834:-1' # first run first
offset = '20'
grp_calib = [{'forward':'3', 'backward':'4', 'alpha':1.13},{'forward':'2', 'backward':'1', 'alpha':1.135}]
#
the_suite = suite(datafile, runlist , grp_calib , offset, startuppath)
dashboard_file = startuppath+'/fit/almgml.822.3-4.1_fit.json' # must exist
the_fit = mufit(the_suite,dashboard_file,scan='T')
plot_range = '0,20000,40'
print('>>>>>>>>>>>>>> close (x) the figure to finish')
the_fitplot = mufitplot(plot_range,the_fit)
P.show()
