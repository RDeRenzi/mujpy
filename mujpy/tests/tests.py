import sys
import time
import traceback

class ParameterizedRunner:
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.results = []

    def run_test_case(self, index, test_dict, test_logic_func):
        """Runs the core test logic with a specific model, fit_type, runlist, grp combination."""

        model = test_dict["model"]
        fit_type = test_dict["fit_type"]
        case_label = f"Case {index:02d} [{model} | {fit_type}]"
        print(f"Running {case_label:.<50}", flush=True ,end ="")
        
        start_time = time.time()
        
        # Execute the test logic which internally handles try/except
        # and returns (is_success, status_message)
        try:
            is_success, message = test_logic_func(test_dict)
            duration = (time.time() - start_time) * 1000
            
            if is_success:
                print(f" [PASS] ({duration:.1f}ms)")
                self.passed += 1
                self.results.append((case_label, "SUCCESS", message, None))
            else:
                print(f" [FAIL] ({duration:.1f}ms)")
                self.failed += 1
                self.results.append((case_label, "FAILED", message, None))
                
        except Exception as e:
            # Captures unexpected framework crashes outside your try/except block
            duration = (time.time() - start_time) * 1000
            print(f" [ERROR] ({duration:.1f}ms)")
            self.failed += 1
            self.results.append((case_label, "CRASHED", f"Unexpected {type(e).__name__}: {e}", traceback.format_exc()))

    def report(self):
        """Prints a scannable summary report and exits with correct status codes."""
        print("\n" + "="*60)
        print("TEST SUMMARY")
        print("="*60)
        print(f"Total executed: {self.passed + self.failed}")
        print(f"Passed:         {self.passed}")
        print(f"Failed/Errors:  {self.failed}")
        
        print("\n" + "-"*60)
        print("DETAILED MESSAGES")
        print("-"*60)
        for label, status, msg, tb in self.results:
            icon = "🟢" if status == "SUCCESS" else "❌"
            print(f"{icon} {label} -> Status: {status} | Msg: {msg}")
            if tb:
                print(f"   Traceback:\n" + "\n".join(f"   {line}" for line in tb.splitlines()[-3:]))
        
        # Exit behavior for CI/CD or automation scripting
        if self.failed > 0:
            sys.exit(1)
        sys.exit(0)

    def test_matrix(self):
        """ setup the test matrix """

        tst_matrix = [
                {    "fit_type": "A1: single run, single group", 
                        "runlist":"822", 
                        "grp":[{'forward':'3', 'backward':'4', 'alpha':1.13}],
                        "json":'.822.3-4.1_fit.json'
                 },
                {    "fit_type": "A20: single run, sequential groups",
                        "runlist":"822", 
                        "grp":[{'forward':'3', 'backward':'4', 'alpha':1.13},{'forward':'2', 'backward':'1', 'alpha':1.13}],
                        "json":'.822.3-4.1_fit.json'
                 },
                {  "fit_type": "B1: sequential runs, single group",
                        "runlist":"822,833,831,829,827", 
                        "grp":[{'forward':'3', 'backward':'4', 'alpha':1.13}],
                        "json":'.822.3-4.1_fit.json'
                 },
                {  "fit_type": "B20: sequential runs and groups",
                        "runlist":"822,833,831,829,827",
                        "grp":[{'forward':'3', 'backward':'4', 'alpha':1.13},{'forward':'2', 'backward':'1', 'alpha':1.13}],
                        "json":'.822.3-4.1_fit.json'
                 },
                {      "fit_type": "A21: single run, global groups",
                        "runlist":"822", 
                        "grp":[{'forward':'3', 'backward':'4', 'alpha':1.13},{'forward':'2', 'backward':'1', 'alpha':1.13}],
                        "json":'.822.3-4+2-1.1_fit.json'
                 },
                {      "fit_type": "B21: sequential runs, global groups",
                        "runlist":"822,833,831,829,827", 
                        "grp":[{'forward':'3', 'backward':'4', 'alpha':1.13},{'forward':'2', 'backward':'1', 'alpha':1.13}],
                        "json":'.822.3-4+2-1.1_fit.json'
                 },
                {    "fit_type": "C1: global runs, single group",
                        "runlist":"822,833,831,829,827", 
                        "grp":[{'forward':'3', 'backward':'4', 'alpha':1.13}],
                        "json":'.822.3-4.C1.1_fit.json'
                 },
                {    "fit_type": "C2: global runs and groups",
                        "runlist":"822,833,831,829,827", 
                        "grp":[{'forward':'3', 'backward':'4', 'alpha':1.13},{'forward':'2', 'backward':'1', 'alpha':1.13}],
                        "json":'.822.3-4+2-1.C2.1_fit.json'}
                    ]
        test_matrix = []
        for model in ["mgml","almgml"]:
            for test_dict in tst_matrix:
                new_dict = test_dict.copy()
                new_dict["model"] = model
                test_matrix.append(new_dict)
        return test_matrix
# -------------------------------------------------------------
# Your Core Test Logic Function
# -------------------------------------------------------------
    def execute_model_evaluation(self,test_dict):
        """
        Actual test execution
    
        Returns a tuple: (bool_success, string_message)
        """
    
        from mujpy.musuite import suite
        from mujpy.mufit import mufit
        from mujpy.mufitplot import mufitplot
        from os import getcwd
        from os.path import join
        import matplotlib.pyplot as P
        from importlib import resources
        startuppath = getcwd()
        runlist = test_dict["runlist"]
        grp = test_dict["grp"]
        model = test_dict["model"]
        fit_type = test_dict["fit_type"]
        json = test_dict["json"]
        datafile = join(resources.files("mujpy.tests").joinpath('data_gps'),'deltat_tdc_gps_0822.bin')
        dashboard_file = join(resources.files("mujpy.tests").joinpath('fit_gps'),model+json)
        offset = "20"
        plot_range = "0,20000,40"
        try:
            the_suite = suite(datafile, 
                              runlist, 
                              grp, 
                              offset, 
                              startuppath)
            the_fit = mufit(the_suite,
                            dashboard_file)
            print('>>>>>>>>>>>>>> close (x) the figure to proceed')
            mufitplot(plot_range,
                      the_fit)
            P.show()      
    
            return True, fit_type
                    
        except Exception as e:
            # If your inner script has its own try/except capturing logic, 
            # format it into the 'failed' message here:
            return False, f"Caught Exception: {str(e)}"

    def test_single(self,index=0):
        """
        runs single tests, index=i with 0<=i>=15
        """

        test_matrix = self.test_matrix()
        print("Starting single Test Run") #.format([t['model']+t['fit_type'] for t in test_matrix]))
        if isinstance(index,int) and index>=0 and index<=len(test_matrix)-1:
            self.run_test_case(
                index = index,
                test_dict = test_matrix[index],
                test_logic_func=self.execute_model_evaluation
                )
            self.report()
        else:
            print('Test index {} is out of range (0<= index <={})'.format(index,len(test_matrix)))

    def list_tests(self):
        """list available tests"""

        test_matrix = self.test_matrix()
        print("List of tests")
        for k,t in enumerate(test_matrix):
            print('{} - {} fit: {}'.format(k,t['model'],t['fit_type']))

# -------------------------------------------------------------
# Execution Block
# -------------------------------------------------------------
if __name__ == "__main__":
    # Define your 16 comprehensive variations transparently

    runner = ParameterizedRunner()
    test_matrix = runner.test_matrix()
    print(f"Starting Parameterized Matrix Test Run ({len(test_matrix)} cases)...\n")
    
    for index, case in enumerate(test_matrix, start=1):
        runner.run_test_case(
            index=index,
            test_dict = case,
            test_logic_func=runner.execute_model_evaluation
        )
        
    runner.report()

