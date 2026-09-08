from os import getcwd #, chdir
import matplotlib.pyplot as P
startuppath = getcwd()
from mujpy.tools.tools import make_links
from mujpy.musuite import suite
from mujpy.mufit import mufit
from mujpy.mufitplot import mufitplot
make_links('GPS')
datafile = 'data/deltat_tdc_gps_0822.bin'
runlist = '822,827:834:-1' # first run first
offset = '20'
grp_calib = [{'forward':'3', 'backward':'4', 'alpha':1.13}]
#
the_suite = suite(datafile, runlist , grp_calib , offset, startuppath)
dashboard_file = startuppath+'/fit/mgml.822.3-4.C1.1_fit.json' # must exist
the_fit = mufit(the_suite,dashboard_file,scan='T' )#,no_fit=True)#,verbose=True, scan = 'T','B','[',)
print('>>>>>>>>>>>>>> close (x) the window to finish')
plot_range = '0,20000,40'
the_fitplot = mufitplot(plot_range,the_fit,guess=True)
P.show()
