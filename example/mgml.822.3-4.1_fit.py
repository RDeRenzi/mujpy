# For ipython only. Choose your style: uncomment one or more of the following magics
from os import getcwd, chdir
startuppath = getcwd() # before changing to git mujpy (if mujpy installed this is not needed)
from mujpy.musuite import suite
from mujpy.mufit import mufit
from mujpy.mufitplot import mufitplot
from matplotlib.pyplot import show
datafile = '../data/deltat_tdc_gps_0822.bin'
runlist = '822' # first run first
offset = '20'
grp_calib = [{'forward':'3', 'backward':'4', 'alpha':1.13}]
#
the_suite = suite(datafile, runlist , grp_calib , offset, startuppath)
dashboard_file = startuppath+'/fit/mgml.822.3-4.1_fit.json' # must exist
the_fit = mufit(the_suite,dashboard_file) #,no_fit=True)
plot_range = '0,20000,40'
print('>>>>>>>>>>>>>> close (x) the figure to finish')
the_fitplot = mufitplot(plot_range,the_fit)
show()
