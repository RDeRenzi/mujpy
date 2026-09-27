from os import getcwd #, chdir
startuppath = getcwd()
from mujpy.musuite import suite
from mujpy.tools.tools import get_title
datafile = '../data/deltat_tdc_gps_0822.bin'
runlist = '822,823:834:-1' # first run first
offset = '20'
grp_calib = [{'forward':'3', 'backward':'4', 'alpha':1.13}]
the_suite = suite(datafile, runlist , grp_calib , offset, startuppath)
for runs in the_suite._the_runs_:
    run = runs[0]
    print('{}: {}  from {} to  {} Tot {} Mev'.format(run.get_runNumber_int(),get_title(run),run.get_timeStart_vector(),run.get_timeStop_vector(),sum(run.get_eventsHisto_vector())))
