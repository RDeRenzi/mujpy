from mujpy.musuite import suite
from mujpy.mufit import mufit
from mujpy.mufitplot import mufitplot
from os import getcwd
from os.path import join
import matplotlib.pyplot as P
from importlib import resources
class tests():
    """
    collects tests of selected mujpy package data sets

    use as 
    """
    def __init__(self,test):

        self.test = test.upper()
        self.startuppath = getcwd()
        # self.do_test() # leave this out 
        self.data = 'data_'+test.lower()
        self.fit = 'fit_'+test.lower()
        self.groups = 'groups'

    def run_test(self):
        """
        generic run method
        """

        from mujpy.musuite import suite
        from mujpy.mufit import mufit
        from mujpy.mufitplot import mufitplot

        try:
            the_suite = suite(self.datafile, 
                              self.runlist , 
                              self.grp_calib , 
                              self.offset, 
                              self.startuppath)
            the_fit = mufit(the_suite,
                            self.dashboard_file)
            print('>>>>>>>>>>>>>> close (x) the figure to finish')
            mufitplot(self.plot_range,
                      the_fit)
            P.show()
            return ''
        except ValueError as e:
            return e

    def do_tests(self):
        """
        call different tests
        """

        if self.test == 'GPS':
        
            self.offset = '20'
            self.plot_range = '0,20000,40'
            self.datafile = resources.files("mujpy.tests").joinpath(self.data).joinpath('deltat_tdc_gps_0822.bin')
            fitdir = resources.files("mujpy.tests").joinpath(self.fit)

            models = ['mgml','almgml']
            for model in models:
            # A1 fit
                print('Demo of {} fit: A1, single run single group'.format(model)) 
                ok = self.runlist = '822'
                self.grp_calib = [{'forward':'3', 'backward':'4', 'alpha':1.13}]
                self.dashboard_file = fitdir.joinpath(model+'.822.3-4.1_fit.json')
                ok = self.run_test()
            # B1 fit
                print('Demo of {} fit: B1, sequential runs single group'.format(model)) 
                ok = self.runlist = '822,833,831,829,827'
                ok = self.run_test()
            # A20 fit
                print('Demo of {} fit: A20, single run sequential groups'.format(model)) 
                ok = self.runlist = '822'
                self.grp_calib = [{'forward':'3', 'backward':'4', 'alpha':1.13},{'forward':'2', 'backward':'1', 'alpha':1.13}]
                ok = self.run_test()
            # B20 fit
                print('Demo of {} fit: B20, sequential runs sequential groups'.format(model)) 
                ok = self.runlist = '822,833,831,829,827'
                ok = self.run_test()           
            # A21 fit
                print('Demo of {} fit: A21, single run global groups'.format(model))
                ok = self.runlist = '822'
                self.dashboard_file = fitdir.joinpath(model+'.822.3-4+2-1.1_fit.json')
                ok = self.run_test()
            # B21 fit
                print('Demo of {} fit: B21, sequential runs global groups'.format(model)) 
                ok = self.runlist = '822,833,831,829,827'
                ok = self.run_test()
            # C1 fit
                print('Demo of {} fit: C1, global runs single group'.format(model)) 
                self.grp_calib = [{'forward':'3', 'backward':'4', 'alpha':1.13}]
                self.dashboard_file = fitdir.joinpath(model+'.822.3-4+2-1.C1.1_fit.json')
                ok = self.run_test()
            # C2 fit
                print('Demo of {} fit: C2, global runs and groups'.format(model)) 
                self.grp_calib = [{'forward':'3', 'backward':'4', 'alpha':1.13},{'forward':'2', 'backward':'1', 'alpha':1.13}]
                self.dashboard_file = fitdir.joinpath(model+'.822.3-4+2-1.C2.1_fit.json')
                ok = self.run_test()

if __name__ == "__main__":
    the_test =unittest('gps')
