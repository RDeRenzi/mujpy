import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
from ctypes import cdll
from ctypes import c_int, c_char_p, c_void_p, c_double
from numpy.ctypeslib import ndpointer

#######################################
# see MuSR_td_PSI_bin.cpp by A. Amato #
# returns 0 for successful reading    #
#######################################

c_double_p = ctypes.POINTER(ctypes.c_double)  # universal double pointer
c_int_p = ctypes.POINTER(ctypes.c_int)        # universal int pointer
                                              # used to pass output results together
 					      # with c_int, c_char_p, c_void_p, c_double

### WORKAROUND, NOT VERY GOOD, CONSULT PIETRO: would not work unless MuSR2py.so si compiled and installed in this particular directory
lib = cdll.LoadLibrary('/usr/local/lib/MuSR2py.so')        # defines the c wrapper library # add smthng like
                                            #$ LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/roberto.derenzi/git/mujpy/mujpy/musr2py/
                                            #                                    path where MuSR2py is
                                            #$ export LD_LIBRARY_PATH 
                                            #
class musr2py(object):                        # defines the python class
  def __init__(self):
    """
    typical usage:: 

      musr2py.read("mydata.bin")  # this is the run data file
      n = musr2py.get_numberTemperature_int()
      # n is now an integer with the number of elements of the PSI monitor values

    """
    
    # set types for all functions
    
    #init                                       defines initializaton method
    lib.init.restype = c_void_p
    
    # read                                      defines method for reading a file
    lib.read.restype = c_int                    # out
    lib.read.argtypes = [c_void_p, c_char_p]    # in

                         
    #int get_histoLength_bin(void *ptr)
    lib.get_histoLength_bin.restype = c_int    # out
    lib.get_histoLength_bin.argtypes = [c_void_p]   # in
    
    #int get_numberHisto_int(void *ptr)
    lib.get_numberHisto_int.restype = c_int    # out
    lib.get_numberHisto_int.argtypes = [c_void_p]   # in
    
    #int get_histo_array_int(void *ptr, int histogram, int * array, int size)
    lib.get_histo_array_int.restype = c_int    # out
    lib.get_histo_array_int.argtypes = [c_void_p, c_int, ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), c_int] # in

    #int get_binWidth_ns(void *ptr)
    lib.get_binWidth_ns.restype = c_double # out
    lib.get_binWidth_ns.argtypes = [c_void_p] # in
    
    #int get_t0_double(void *ptr, int histogram)
    lib.get_t0_double.restype = c_double # out
    lib.get_t0_double.argtypes = [c_void_p, c_int] # in
    
    #int get_sample(void *ptr)
    lib.get_sample.restype = c_int # out ?
    lib.get_sample.argtypes = [c_void_p, c_char_p] # in
    
    #int get_field(void *ptr)
    lib.get_field.restype = c_int # out ?
    lib.get_field.argtypes = [c_void_p, c_char_p] # in
    
    #int get_orient(void *ptr)
    lib.get_orient.restype = c_int # out ?
    lib.get_orient.argtypes = [c_void_p, c_char_p] # in
    
    #int get_temp(void *ptr)
    lib.get_temp.restype = c_int # out ?
    lib.get_temp.argtypes = [c_void_p, c_char_p] # in
    
    #int get_comment(void *ptr)
    lib.get_comment.restype = c_int # out ?
    lib.get_comment.argtypes = [c_void_p, c_char_p] # in
    
    #int get_timeStart_vector(void *ptr)
    lib.get_timeStart_vector.restype = c_int    # out
    lib.get_timeStart_vector.argtypes = [c_void_p, c_char_p] # in

    #int get_timeStop_vector(void *ptr)
    lib.get_timeStop_vector.restype = c_int    # out
    lib.get_timeStop_vector.argtypes = [c_void_p, c_char_p] # in

    #int get_numberTemperature_int(void *ptr)
    lib.get_numberTemperature_int.restype = c_int    # out
    lib.get_numberTemperature_int.argtypes = [c_void_p]   # in

    #int get_temperatures_vector(void *ptr, double * array)
    lib.get_temperatures_vector.restype = c_int    # out
    lib.get_temperatures_vector.argtypes = [c_void_p, ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")] # in

    #int get_devTemperatures_vector(void *ptr, double * array)
    lib.get_devTemperatures_vector.restype = c_int    # out
    lib.get_devTemperatures_vector.argtypes = [c_void_p, ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")] # in

    #int get_runNumber_int(void *ptr)
    lib.get_runNumber_int.restype = c_int    # out
    lib.get_runNumber_int.argtypes = [c_void_p]   # in
    
    #int get_eventsHisto_vector(void *ptr, int * array)
    lib.get_eventsHisto_vector.restype = c_int # out
    lib.get_eventsHisto_vector.argtypes = [c_void_p, ndpointer(ctypes.c_long, flags="C_CONTIGUOUS")] # in notice long

    # destroy                              defines method for deallocation
    lib.destroy.argtypes = [c_void_p]
    
    
    self.obj = lib.init()

    
      
  def __exit__(self):
    lib.destroy(self.obj)
  
  def read(self, filename):
    """
    usage::
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
    """
    fnamePtr = ctypes.create_string_buffer(filename.encode('utf8'),120) # pointer to string  and forced conversion to utf8
    return lib.read(self.obj, fnamePtr)
    
    
  def get_numberHisto_int(self):
    """
    usage::
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     n = run.get_numberHisto_int()
     # n is now an integer with the number of histrograms of this run
    """
    return lib.get_numberHisto_int(self.obj)
     
  def get_histoLength_bin(self):
    """
    usage::
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     n = run.get_histoLength_bin()
     # n is now an integer with the number of bins in the histrograms of this run

    """
    return lib.get_histoLength_bin(self.obj)
    
  def get_histo_array_int(self, histogram):
    """
    usage::
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     h = run.get_histo_array_int(2)
     # h is now a numpy array of integers with the counts of the 3rd histogram

    (remember python indices starts from 0)
    """
    size = self.get_histoLength_bin()
    dataPtr = np.empty(size, dtype=np.int32)
    ret = lib.get_histo_array_int(self.obj, histogram, dataPtr, dataPtr.size)
    if ret == 0:
      return(dataPtr)
  
  def get_numberTemperature_int(self):
    """
    usage:: 
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     n = run.get_numberTemperature_int()
     # n is now an integer with the number of elements of the PSI monitor values

    """
    return lib.get_numberTemperature_int(self.obj)

  def get_temperatures_vector(self):
    """
    usage:: 
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     T = run.get_temperatures_vector()
     # T is now a numpy array of PSI monitor values

    """

    size = self.get_numberTemperature_int() # determine how many elements in vector 
    temperaturePtr = np.empty(size, dtype=np.float64)  # allocate an empty numpy array
                                                       # of float 64 of length size
                                                       # this is the pointer
    ret = lib.get_temperatures_vector(self.obj, temperaturePtr) # call the wrapper
    if ret == 0:
      return(temperaturePtr)  # return the pointer to the numpy array
                              # use T = mp2.get_temperatures_vector()
    else:
      print("Error in MuSR_td_PSI_bin, returning empty temperature numpy array")
      return(temperaturePtr)  
  
  def get_devTemperatures_vector(self):
    """
    usage::
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     eT = run.get_devTemperatures_vector()
     # eT is now a numpy array of PSI monitor std values

    """

    size = self.get_numberTemperature_int() # determine how many elements in vector 
    etemperaturePtr = np.empty(size, dtype=np.float64) # allocate an empty numpy array
                                                       # of float 64 of length size
                                                       # this is the pointer
    ret = lib.get_devTemperatures_vector(self.obj, etemperaturePtr) # call the wrapper
    if ret == 0:
      return(etemperaturePtr)  # return the pointer to the numpy array
                              # use T = mp2.get_temperatures_vector()
    else:
      print("Error in MuSR_td_PSI_bin, returning empty temperature std numpy array")
      return(etemperaturePtr)  
  
  def get_binWidth_ns(self):
    """
    usage:: 
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     dt = run.get_binWidth_ns()
     # dt is now a float with the time resolution in ns

    """
    return lib.get_binWidth_ns(self.obj)
  
  def get_t0_double(self, histogram):
    """
    usage::
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     t0 = run.get_t0_double(0)
     # t0 is now a float with the time of t0 in ns

    """
    return lib.get_t0_double(self.obj, histogram)

  def get_sample(self):
    """
    usage:: 
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     samplename = run.get_sample()
     # samplename is now a string with the sample name

    """
    samplePtr = ctypes.create_string_buffer(60) # pointer to the string
    ret = lib.get_sample(self.obj, samplePtr)
    if ret == 0:
      return samplePtr.value.decode('utf8')

  def get_field(self):
    """
    usage::
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     fieldstr = run.get_field()
     # fieldstr is now a string with the field value 

    """
    fieldPtr = ctypes.create_string_buffer(60) # pointer to the string
    ret = lib.get_field(self.obj, fieldPtr)
    if ret == 0:
      return fieldPtr.value.decode('utf8')

  def get_orient(self):
    """
    usage::
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     orientstr = run.get_orient()
     # orientstr is now a string with the sample orientation

    """
    orientPtr = ctypes.create_string_buffer(60) # pointer to the string
    ret = lib.get_orient(self.obj, orientPtr)
    if ret == 0:
      return orientPtr.value.decode('utf8')

  def get_temp(self):
    """
    usage::
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     tempstr = run.get_temp()
     # tempstr is now a string with the nominal temperature

    """
    tempPtr = ctypes.create_string_buffer(60) # pointer to the string
    ret = lib.get_temp(self.obj, tempPtr)
    if ret == 0:
      return tempPtr.value.decode('utf8')

  def get_comment(self):
    """
    usage: 
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     comment = run.get_comment()
     # comment is now a string with the run comment
    """
    commentPtr = ctypes.create_string_buffer(62) # pointer to the string
    ret = lib.get_comment(self.obj, commentPtr)
    if ret == 0:
      return commentPtr.value.decode('utf8')

  def get_eventsHisto_vector(self):
    """
    usage::
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     eventsHisto = run.get_eventsHisto_vector()
     # eventsHisto is a numpy array of integers containing the number of events per histo

    """
    size = self.get_numberHisto_int() # determine how many elements in vector 
    eventsHistoPtr = np.empty(size, dtype=np.int64) # allocate an empty numpy array
                                                       # of int 64 of length size 
                                                       # this is the pointer 
    ret = lib.get_eventsHisto_vector(self.obj, eventsHistoPtr) # call the wrapper
    if ret == 0:
      return(eventsHistoPtr)  # return the pointer to the numpy array
                              # use events = musr2py.get_eventsHisto_vector()
    else:
      print("Error in MuSR_td_PSI_bin, returning empty number of events  per histo std numpy array")
      return(eventsHistoPtr)  

  def get_runNumber_int(self):
    """
    usage: 
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     nrun = run.get_runNumber_int()
     # n is now an integer with the number of this run
    """
    return lib.get_runNumber_int(self.obj)

  def get_timeStart_vector(self):
    """
    usage::
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     timeStart = run.get_timeStart_vector()
     # timeStart is the start run hh:mm:ss dd:mm:yyy string 

    """
    timedatePtr = ctypes.create_string_buffer(18)
    ret = lib.get_timeStart_vector(self.obj, timedatePtr)
    if ret == 0:
      return timedatePtr.value.decode('utf8')

  def get_timeStop_vector(self):
    """
    usage::
     from musr2py import musr2py as psiload
     run = psiload()  # this is the run data file
     run.read("mydata.bin")
     timeStop = run.get_timeStop_vector()
     # timeStop is the stop run hh:mm:ss dd:mm:yyy string

    """
    timedatePtr = ctypes.create_string_buffer(18)
    ret = lib.get_timeStop_vector(self.obj, timedatePtr)
    if ret == 0:
      return timedatePtr.value.decode('utf8')

#  def test(self):
if __name__ == '__main__':
    """
    usage::
     write a script with
     from musr2py import musr2py
     # reads many things from a PSI data file test.bin

    """
    m2p = musr2py()
    m2p.read("test.bin")
    print(str(m2p.get_numberHisto_int())+' histograms in this run')
    print(m2p.get_histo_array_int(0))
    print(m2p.get_histo_array_int(1))
    print(m2p.get_histo_array_int(2))
    print(m2p.get_histo_array_int(3))
    print(str(m2p.get_binWidth_ns())+' ns/bin')
    print(str(m2p.get_numberTemperature_int())+' recorded monitors')
    np.set_printoptions(precision=3)
    print(str(m2p.get_temperatures_vector())+' +- ') 
    print(str(m2p.get_devTemperatures_vector())+' K')
    print(m2p.get_sample()+' '+m2p.get_temp()+' '+m2p.get_field()+' G  '+m2p.get_temp())
    print('Comment: '+m2p.get_comment())
    print('Start run '+m2p.get_timeStart_vector())
    print('Stop run '+m2p.get_timeStop_vector())
    print('Number of events per histo  = {}'.format(m2p.get_eventsHisto_vector()))
