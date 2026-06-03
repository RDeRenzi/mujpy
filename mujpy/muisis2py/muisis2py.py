from nexusformat.nexus import *
import numpy as np

class munxs2py():                        # defines the python class

    def __init__(self,path_filename):
        self.readingOOK = self.read(path_filename)  # reads data

    def read(self,path_filename):
        with nxload(path_filename) as run:
            try:
               	self.run = run
               	read_ok = True
            except:
                read_ok = False
        return read_ok

    def _rbn(self,hist,binning):
        return hist[:(hist.size//binning)*binning].reshape(hist.size//binning,binning).sum(axis=1)

    def Filename(self):
        return str(self.run.file_name).split('\\')[-1]
        
    def readingOK(self):
        return self.readingOOK

    def Show(self):
        k1 = len(self.run["RunHeader;1"].members['fFolders'][:])
        for k in range(k1):
            k2 = len(self.run["RunHeader;1"].members['fFolders'][k][:])
            for l in range(k2):
                print('{}'.format(self.run["RunHeader;1"].members['fFolders'][k][l][:]))
                
    def get_proposal_number(self):
        return 'RB{}'.format(self.run.raw_data_1.experiment_identifier)
        
    def get_proposer(self):
        return str(self.run.raw_data_1.user_1.name)
 
    def get_facility(self): 
        return str(self.run.raw_data_1.instrument.source.name)

    def get_numberHisto_int(self):
        """
        usage::
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          n = run.get_numberHisto_int()
          # n is now an integer with the number of histrograms of this run
        """
        return np.array(self.run.raw_data_1.detector_1.counts).shape[1]
    
    def get_histoLength_bin(self):
        """
        usage::
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          n = run.get_histoLength_bin()
          # n is now an integer with the number of bins in the histrograms of this run

        """
        return np.array(self.run.raw_data_1.detector_1.counts).shape[2]
    
    
    def get_RedGreen_offsets(self):
        """
        usage::
          run = muld(path2file)  # this is the class instance initiated on the data file
          period = run.get_RedGreen_offsets() # offsets is now a list of integer imdices for the first dimension of the histo vector
          # check on a RedGreenm mode data set
       """
        return [int(self.run.raw_data_1.detector_1.period_index)-1] 
    
    def get_binWidth_us(self):
        return float(self.run.raw_data_1.detector_1.raw_time[1])
    
    def get_binWidth_ns(self):
        return float(self.run.raw_data_1.detector_1.raw_time[1]*1000)
    
    def get_histo_int(self,histogram,nbin):
        """
        usage::
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          # histogram in range(self.get_numberHisto_int())
          h = run.get_histo_int(3,100)
          # h is now the integer count of bin 100 in the 3rd histogram (PSI convention), python index 2
        (remember python indices starts from 0)
        """
        return int(self.run.raw_data_1.detector_1.counts[0,histogram,nbin])

    def get_histo_vector(self,histogram,binning=1):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          # histogram in range(self.get_numberHisto_int())
          h = run.get_histo_vector(offset)
          # h is a column of numpy row array (not list!) of integers containing the number of events per histo for the histos with a given offset
        """
        if binning == 1:
            return np.array(self.run.raw_data_1.detector_1.counts[0,histogram,:])
        else:
            return self._rbn(np.array(self.run.raw_data_1.detector_1.counts[0,histogram,:]),binning)
    
    def get_histo_vector_no0(self,histogram,binning=1):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          # histogram in range(self.get_numberHisto_int())
          h = run.get_histo_vector_no0(offset)
          # h is a column of numpy row array (not list!) of integers containing the number of events per histo for the histos with zeros replaced by 0.1
        """
        hist = np.array(self.run.raw_data_1.detector_1.counts[0,histogram,:])
        if binning != 1:
            hist = self._rbn(hist,binning)
        return np.where(hist != 0, hist, 0.1)
    
    def get_histo_fromt0_vector(self,histogram,binning=1,offset=0):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          # histogram in range(self.get_numberHisto_int())
          h = run.get_histo_fromt0_vector(offset)
          # h is a column of numpy row array (not list!) of integers containing the number of events per histo for the histos with a given offset after t0 
        """
        hist = np.array(self.run.raw_data_1.detector_1.counts[0,histogram,:])[self.get_t0_int(histogram)+offset:]
        if binning == 1:
            return hist
        else:
            return self._rbn(hist,binning)

    def get_numberTemperature_int(self):
        """
        usage:: 
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          n = run.get_numberTemperature_int()
          # n is now an integer with the number of elements of the ISIS temperature log values

        """
        return 1
    
    def get_temperatures_vector(self):
        """
        usage:: 
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          T = run.get_temperatures_vector()
          # T is now a numpy array of ISIS temperature log values
        """
        return np.array(self.run.raw_data_1.selog.Temp_Sample.value_log.value)
    
    def get_binWidth_ns(self):
        """
        usage:: 
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          dt = run.get_binWidth_ns()
          # dt is now a float with the time resolution in ns

        """
        return float(self.run.raw_data_1.instrument.detector_1.resolution)/1000.
  
    def get_binWidth_us(self):
        """
        usage:: 
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          dt = run.get_binWidth_ns()
          # dt is now a float with the time resolution in ns

        """
        return float(self.run.raw_data_1.instrument.detector_1.resolution)/1.e6
  
    def get_t0_int(self, histogram=0):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          t0 = run.get_t0_double()
          # t0 is now a float with the t0 bin 
        """
        return round(self.get_t0_double(histogram))
        
    def get_t0_double(self, histogram=0):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          t0 = run.get_t0_double()
          # t0 is now a float with the t0 bin 
        """
        return float(self.run.raw_data_1.detector_1.time_zero)/float(self.run.raw_data_1.instrument.detector_1.resolution)*1e6


    def get_sample(self):
        """
        usage:: 
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          samplename = run.get_sample()
          # samplename is now a string with the sample name

        """
        return str(self.run.raw_data_1.sample.name)

    def get_field(self):
        """
        usage::
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          fieldstr = run.get_field()
          # fieldstr is now a string with the field value 

        """
        return '{} {}'.format(str(self.run.raw_data_1.sample.magnetic_field),self.run.raw_data_1.sample.magnetic_field.units[0])

    def get_orient(self):
        """
        usage::
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          orientstr = run.get_orient()
          # orientstr is now a string with the sample orientation

        """
        return str(self.run.raw_data_1.sample.shape)

    def get_temp(self):
        """
        usage::
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          tempstr = run.get_temp()
          # tempstr is now a string with the nominal temperature

        """
        return '{} {}'.format(str(self.run.raw_data_1.sample.temperature),self.run.raw_data_1.sample.temperature.units[0])

    def get_comment(self):
        """
        usage: 
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          comment = run.get_comment()
          # comment is now a string with the run comment
        """
        return str(self.run.raw_data_1.notes)

    def get_eventsHisto_vector(self):
        """
        usage::
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          eventsHisto = run.get_eventsHisto_vector()
          # eventsHisto is a numpy array of integers containing the number of events per histo

        """
        return np.array(self.run.raw_data_1.detector_1.counts)[0].sum(axis=1)

    def get_runNumber_int(self):
        """
        usage: 
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          nrun = run.get_runNumber_int()
          # n is now an integer with the number of this run
        """
        return int(self.run.raw_data_1.run_number)

#  def get_timeTemperature_vector(self):
#    """
#    usage::
#      from muisis2py import muisis2py as muld
#      path2file = 'path and filename'
#      run = muld(path2file,'r')  # this is the run data nexus file
#      t = run.get_timeTemperature_vector()
#      # timeStart is the start run yyyy-mm-ddThh:mm:ss string 
#      
#    """
#    from datetime import datetime as DT
#    from numpy import linspace
#    t_format = "%Y-%m-%dT%H:%M:%S"
#    delta = DT.strptime(self.get_timeStop_vector(),t_format) - DT.strptime(self.get_timeStart_vector(),t_format)
#    nlogs = self.get_numberTemperature_int()
#    return linspace(0,delta.total_seconds(),nlogs)

    def get_timeStart_vector(self):
        """
        usage::
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          timeStart = run.get_timeStart_vector()
          # timeStart is the start run yyyy-mm-ddThh:mm:ss string 

        """
        return str(self.run.raw_data_1.start_time).split('T')
      
    def get_timeStop_vector(self):
        """
        usage::
          from muisis2py import muisis2py as muld
          path2file = 'path and filename'
          run = muld(path2file,'r')  # this is the run data nexus file
          timeStop = run.get_timeStop_vector()
          # timeStop is the stop run yyyy-mm-ddThh:mm:ss string

        """
        return str(self.run.raw_data_1.end_time).split('T')

#from munxs.tree import NeXusTree

#class muisis2py(NeXusTree):                        # defines the python class
class muisis2py(): # broken since 
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
    return str(self.readpath(b'magnetic_field'))+'G '

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
    return str(self.readpath(b'temperature'))+'K '

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

  def get_timeTemperature_vector(self):
    """
    usage::
      from muisis2py import muisis2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      t = run.get_timeTemperature_vector()
      # timeStart is the start run yyyy-mm-ddThh:mm:ss string 
      
    """
    from datetime import datetime as DT
    from numpy import linspace
    t_format = "%Y-%m-%dT%H:%M:%S"
    delta = DT.strptime(self.get_timeStop_vector(),t_format) - DT.strptime(self.get_timeStart_vector(),t_format)
    nlogs = self.get_numberTemperature_int()
    return linspace(0,delta.total_seconds(),nlogs)

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
      from muisis2py import munxs2py as muld
      path2file = 'path and filename'
      run = muld(path2file,'r')  # this is the run data nexus file
      # reads many things from an ISIS data file EMU00020882.nxs

    """
    import numpy as np
    m2p = munxs2py("../../data/nexus/EMU00020882.nxs",'r')
    print(str(m2p.get_numberHisto_int())+' histograms in this run')
    print(m2p.get_histo_array_int(0))
    print(m2p.get_histo_array_int(1))
    print(m2p.get_histo_array_int(2))
    print(m2p.get_histo_array_int(3))
    print(str(m2p.get_binWidth_ns())+' ns/bin')
    print(str(m2p.get_numberTemperature_int())+' recorded temperature logs')
    np.set_printoptions(precision=3)
    print(m2p.get_temperatures_vector()) 
    print(m2p.get_sample()+' '+m2p.get_temp()+' K'+m2p.get_field()+' G  ')
    print('Comment: '+m2p.get_comment())
    print('Start run '+m2p.get_timeStart_vector())
    print('Stop run '+m2p.get_timeStop_vector())
    print('Number of events per histo  = {}'.format(m2p.get_eventsHisto_vector()))


