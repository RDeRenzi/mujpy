from munxs.tree import NeXusTree

class muisis2py(NeXusTree):                        # defines the python class
    
  def get_numberHisto_int(self):
    """
    usage::
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      n = run.get_numberHisto_int()
      # n is now an integer with the number of histrograms of this run

    """
    self.opengrouppath(b'/run/histogram_data_1')
    return self.readpath(b'counts').shape[0]
    
  def get_histoLength_bin(self):
    """
    usage::
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      n = run.get_histoLength_bin()
      # n is now an integer with the number of bins in the histrograms of this run

    """
    self.opengrouppath(b'/run/histogram_data_1')
    return self.readpath(b'counts').shape[1]
    
  def get_histo_array_int(self, histogram):
    """
    usage::
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      h = run.get_histo_array_int(2)
      # h is now a numpy array of integers with the counts of the 3rd histogram

    (remember python indices starts from 0)
    """
    self.opengrouppath(b'/run/histogram_data_1')
    return self.readpath(b'counts')[histogram,:]
    
  def get_numberTemperature_int(self):
    """
    usage:: 
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      n = run.get_numberTemperature_int()
      # n is now an integer with the number of elements of the ISIS temperature log values

    """
    self.opengrouppath(b'/run/temperature_log_1')
    return self.readpath(b'values').shape[0]
    
  def get_temperatures_vector(self):
    """
    usage:: 
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      T = run.get_temperatures_vector()
      # T is now a numpy array of ISIS temperature log values
    """
    self.opengrouppath(b'/run/temperature_log_1')
    return self.readpath(b'values')
    
  def get_binWidth_ns(self):
    """
    usage:: 
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      dt = run.get_binWidth_ns()
      # dt is now a float with the time resolution in ns

    """
    self.opengrouppath(b'/run/histogram_data_1')
    return float(self.readpath(b'resolution'))/1000.
  
  def get_t0_double(self, histogram):
    """
    usage::
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      t0 = run.get_t0_double()
      # t0 is now a float with the time of t0 in ns

    """
    self.opengrouppath(b'/run/histogram_data_1')
    return self.readpath(b'time_zero')

  def get_sample(self):
    """
    usage:: 
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      samplename = run.get_sample()
      # samplename is now a string with the sample name

    """
    self.opengrouppath(b'/run/sample')
    return self.readpath(b'name').decode('utf-8')

  def get_field(self):
    """
    usage::
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      fieldstr = run.get_field()
      # fieldstr is now a string with the field value 

    """
    self.opengrouppath(b'/run/sample')
    return str(self.readpath(b'magnetic_field'))

  def get_orient(self):
    """
    usage::
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      orientstr = run.get_orient()
      # orientstr is now a string with the sample orientation

    """
    self.opengrouppath(b'/run/sample')
    return self.readpath(b'shape').decode('utf-8')

  def get_temp(self):
    """
    usage::
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      tempstr = run.get_temp()
      # tempstr is now a string with the nominal temperature

    """
    self.opengrouppath(b'/run/sample')
    return str(self.readpath(b'temperature'))

  def get_comment(self):
    """
    usage: 
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      comment = run.get_comment()
      # comment is now a string with the run comment
    """
    self.opengrouppath(b'/run/instrument')
    comment = self.readpath(b'name').decode('utf-8')
    comment += ' '
    self.opengrouppath(b'/run/sample')
    comment += self.readpath(b'environment').decode('utf-8')
    comment += ' '
    self.opengrouppath(b'/run/sample')
    comment += self.readpath(b'magnetic_field_state').decode('utf-8')
    return comment

  def get_eventsHisto_vector(self):
    """
    usage::
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      eventsHisto = run.get_eventsHisto_vector()
      # eventsHisto is a numpy array of integers containing the number of events per histo

    """
    self.opengrouppath(b'/run/histogram_data_1')
    return self.readpath(b'counts').sum(axis=1)

  def get_runNumber_int(self):
    """
    usage: 
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      nrun = run.get_runNumber_int()
      # n is now an integer with the number of this run
    """
    self.opengrouppath(b'/run')
    return self.readpath(b'number')

  def get_timeStart_vector(self):
    """
    usage::
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      timeStart = run.get_timeStart_vector()
      # timeStart is the start run yyyy-mm-ddThh:mm:ss string 

    """
    self.opengrouppath(b'/run')  
    return self.readpath(b'start_time').decode('utf-8')
      
  def get_timeStop_vector(self):
    """
    usage::
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      timeStop = run.get_timeStop_vector()
      # timeStop is the stop run yyyy-mm-ddThh:mm:ss string

    """
    self.opengrouppath(b'/run')
    return self.readpath(b'stop_time').decode('utf-8')

if __name__ == '__main__':
    """
    usage::
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      # reads many things from an ISIS data file EMU00020882.nxs

    """
    import numpy as np
    m2p = muisis2py("../../data/nexus/EMU00020882.nxs",'r')
    print(str(m2p.get_numberHisto_int())+' histograms in this run')
    print(m2p.get_histo_array_int(0))
    print(m2p.get_histo_array_int(1))
    print(m2p.get_histo_array_int(2))
    print(m2p.get_histo_array_int(3))
    print(str(m2p.get_binWidth_ns())+' ns/bin')
    print(str(m2p.get_numberTemperature_int())+' recorded temperature logs')
    np.set_printoptions(precision=3)
    print(str(m2p.get_temperatures_vector())+' K') 
    print(m2p.get_sample()+' '+m2p.get_temp()+' '+m2p.get_field()+' G  ')
    print('Comment: '+m2p.get_comment())
    print('Start run '+m2p.get_timeStart_vector())
    print('Stop run '+m2p.get_timeStop_vector())
    print('Number of events per histo  = {}'.format(m2p.get_eventsHisto_vector()))


