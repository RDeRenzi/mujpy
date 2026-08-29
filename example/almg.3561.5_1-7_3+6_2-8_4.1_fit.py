# try out root datafiles 
from os import getcwd #, chdir
import matplotlib.pyplot as P
startuppath = getcwd() # before changing to git mujpy (if mujpy installed this is not needed)
# chdir('/home/roberto.derenzi/git/mujpy/')
# print(getcwd())
from mujpy.musuite import suite
from mujpy.mufit import mufit
from mujpy.mufitplot import mufitplot
datafile = 'data/root/lem24_his_3561.root'
runlist = '3561' # first run first
offset = '20'
grp_calib = [{'forward':'5,1', 'backward':'7,3', 'alpha':1.05},{'forward':'6,2', 'backward':'8,4', 'alpha':1.07}]
#
the_suite = suite(datafile, runlist , grp_calib , offset, startuppath)
dashboard_file = startuppath+'/fit/almg.3561.5_1-7_3+6_2-8_4.1_fit.json' # must exist
the_fit = mufit(the_suite,dashboard_file)#,no_fit=True)
print('>>>>>>>>>>>>>> close (x) the figure to finish')
plot_range = '0,25000,200'
the_fitplot = mufitplot(plot_range,the_fit)
P.show()
