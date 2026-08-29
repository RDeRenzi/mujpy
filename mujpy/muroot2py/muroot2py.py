import uproot
import numpy as np

class muroot2py():                        # defines the python class

    """
    the user method names are standardized, when possible, to those of MuSR_td_PSI_bin() by Amato Raselli
    usage::
      from muroot2py import muroot2py as muld
      path2file = 'path and filename'
      run = muld()  # this is the method instance
     """
    def __init__(self, path_filename):
        self.readingOOK = self.read(path_filename)  # reads data
        self._header_dictionary() # creates header dictionary self._header_dictionary
        self._histo_dictionary()  # creates histo dictionary self._histo_dictionary
                                  # used by all get_ calls
    def read(self,path_filename):
        """
        usage::
          run = muld(path2file)  # this is the method instance 
          initiated on the data file
          returns True (success) or False (insuccess)
        """    
        with uproot.open(path_filename) as run:
            try:
               	self.run = run
               	read_ok = True
            except:
                read_ok = False
        return read_ok

# internal methods

    def _runheader(self):
        return self.run["RunHeader"].members["fFolders"]

    def _extract(self,num,key):	
        return str(self._runheader()[0][num]).partition(key)[2].partition(' -@')[0]

    def _temperature_error_float(self):
        return float(self._header_dictionary['Sample Temperature string'].partition(' +- ')[2].partition(' ')[0])
        
    def _header_key(self,k):
        return str(self._runheader()[0][k]).split(':')[0].split('- ')[1]
        
    def _header_type(self,k):
        return int(str(self._runheader()[0][k])[-1])
        
    def _header_value(self,k):
        header = str(self._runheader()[0][k])
        muroot_type = int(header[-1])
        value = header.partition(self._header_key(k)+': ')[2].partition(' -@')[0]
        if muroot_type == 0: # string
            pass
        elif muroot_type == 1: # integer
            value = int(value)
        elif muroot_type == 2: # float
            value = float(value)
        elif muroot_type == 3: # physical quantity
            value = value.split()[0]
            value = float(value) if '.' in value else int(value)
        return value

    def _header_dictionary(self): 
        self._header_dictionary = {}
        for k in range(len(self._runheader()[0][:])):
            self._header_dictionary[self._header_key(k)] = self._header_value(k)
            if self._header_key(k) == 'Sample Temperature' or self._header_key(k) == 'Sample Magnetic Field':
                header = str(self._runheader()[0][k])
                self._header_dictionary[self._header_key(k)+' string'] = header.partition(self._header_key(k)+': ')[2].partition(' -@')[0]

    def _histo_key(self,histogram,k):
        return str(self._runheader()[1][histogram][k]).split(':')[0].split('- ')[1]
        
    def _histo_type(self,histogram,k):
        return int(str(self._runheader()[1][histogram][k])[-1])
        
    def _histo_value(self,histogram,k):
        header = str(self._runheader()[1][histogram][k])
        muroot_type = int(header[-1])
        value = header.partition(self._histo_key(histogram,k)+': ')[2].partition(' -@')[0]
        if muroot_type == 0: # string
            pass
        elif muroot_type == 1: # integer
            value = value
        elif muroot_type == 2: # float
            value = float(value)
        elif muroot_type == 3: # physical quantity
            value = value.split()[0]
            value = float(value) if '.' in value else int(value)
        return value

    def _histo_dictionary(self): 
        self._histo_dictionary = {}       
        for histo in range(len(self._runheader()[1][:])):
            for k in range(len(self._runheader()[1][histo][:])):
                self._histo_dictionary[self._histo_key(histo,k)+' '+str(histo)] = self._histo_value(histo,k)

    def _histos(self):
        return self.run["histos"].members['fFolders'][0].members['fFolders'] 

    def _plots(self):
        if self._header_dictionary['Instrument']=='LEM':
            return self.run["histos"].members['fFolders'][5].members['fFolders'] 
        else:
            return self.run["histos"].members['fFolders'][1].members['fFolders'] 
    
    def _rbn(self,hist,binning):
        return hist[:(hist.size//binning)*binning].reshape(hist.size//binning,binning).sum(axis=1)
        
# user methods

    def Filename(self):
        return self._header_dictionary['File Name']
        
    def readingOK(self):
        return self.readingOOK

    def Show(self):
        k1 = len(self.run["RunHeader;1"].members['fFolders'][:])
        for k in range(k1):
            k2 = len(self.run["RunHeader;1"].members['fFolders'][k][:])
            for l in range(k2):
                print('{}'.format(self.run["RunHeader;1"].members['fFolders'][k][l][:]))
                
    def get_binWidth_ns(self):
        """
        usage:: 
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          dt = run.get_binWidth_ns()
          # dt is now a float with the time resolution in ns

        """
        return self._header_dictionary['Time Resolution']

    def get_binWidth_us(self):
        """
        usage:: 
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          dt = run.get_binWidth_ns()
          # dt is now a float with the time resolution in ns

        """
        return self._header_dictionary['Time Resolution']/1000.

    def get_histoLength_bin(self):
        """
        usage::
          run = muld(path2file)  # this is the class instance initiated on the data file
          l = run.get_histoLength_bin() # l is now an integer with the number of bins in PSI detector histo
        """
        return int(self._histo_dictionary['Histo Length 0'])

    def get_numberHisto_int(self):
        """
        usage::
          run = muld(path2file)  # this is the class instance initiated on the data file
          l = run.get_numberHisto_int(histo) # l is now an integer with the number of detectors
       """
        return self._header_dictionary['No of Histos']
    
    def get_RedGreen_offsets(self):
        """
        usage::
          run = muld(path2file)  # this is the class instance initiated on the data file
          offsets = run.get_RedGreen_offsets() # offsets is now a list of integer offsets for PSI RedGreen mode
          # beware, histograms are contuiguous, the offset is in the PSI Histo Number
       """
        return [int(k) for k in self._header_dictionary['RedGreen Offsets'].split(';')]
    
    def get_histo_int(self, histogram,nbin):
        """
        usage::
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          # histogram in range(self.get_numberHisto_int())
          h = run.get_histo_int(3,100)
          # h is now the integer count of bin 100 in the 3rd histogram (PSI convention), python index 2
        (remember python indices starts from 0)
        """
        return int(self._histos()[histogram].counts()[nbin])

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
            return np.array(self._histos()[histogram].counts())
        else:
            return self._rbn(np.array(self._histos()[histogram].counts()),binning)
    
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
        hist = np.array(self._histos()[histogram].counts())
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
        hist = np.array(self._histos()[histogram].counts())[self.get_t0_int(histogram)+offset:]
        if binning == 1:
            return hist
        else:
            return self._rbn(hist,binning)

    def get_histo_fromt0_minus_bckgrd_vector(self,histogram,first_bkg,last_bkg,binning=1,offset=0):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          h = run.get_histo_vector(histogram,offset)
          # histogram in range(self.get_numberHisto_int())
          # h is a column of numpy row array (not list!) of integers containing the number of events per histo for the histos with a given offset after t0
          # with background subrtaction as estimated between first_bkg and last_bkg included
        """
        hist = np.array(self._histos()[histogram].counts())[self.get_t0_int(histogram)+offset:]
        background = np.array(self._histos()[histogram].counts())[first_bkg:last_bkg+1].sum()/(last_bkg+1-first_bkg)
        if binning == 1:
            return hist
        else:
            return self._rbn(hist,binning)-background   
    
    def get_instrument(self):
        return self._header_dictionary['Instrument']

    def get_title(self):
        return self._header_dictionary['Run Title']

    def get_runNumber_int(self):
        return self._header_dictionary['Run Number']

    def get_beamline(self):
        '''
        '''
        return self._header_dictionary['Muon Source']

    def get_transport_energy(self):
        if self._header_dictionary['Instrument']=='LEM':
            return self._header_dictionary['Moderator HV']

    def get_sample_HV(self):
        if self._header_dictionary['Instrument']=='LEM':
           return self._header_dictionary['Sample HV']

    def get_implantation_energy(self):
        if self._header_dictionary['Instrument']=='LEM':
            return self._header_dictionary['Implantation Energy']

    def get_spin_angle(self):
        if self._header_dictionary['Instrument']=='LEM':
            return self._header_dictionary['Muon Spin Angle']

    def get_proposal_number(self):
        return self._header_dictionary['Proposal Number']
        
    def get_proposer(self):
        if self._header_dictionary['Instrument']=='LEM':
            return self._header_dictionary['Main Proposer']
 
    def get_comment(self):
        """
        usage: 
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          run.red(path2file)
          comment = run.get_comment()
          # comment is now a string with the run comment
        """
        return self._header_dictionary['Comment']

    def get_temp(self):
        """
        usage:: 
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          T = run.get_temp()
          # T is now a string with error
        """
        return self._header_dictionary['Sample Temperature string'].split(' +- ')[0]+'K'

    def get_temperature(self):
        """
        usage:: 
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          T = run.get_temperature()
          # T is now a float ?
        """
        return self._header_dictionary['Sample Temperature']

    def get_t0_int(self, histogram):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          t0 = run.get_t0_double()
          # t0 is now a float with the t0 bin 
        """
        return int(self._histo_dictionary['Time Zero Bin '+str(histogram)])

    def get_t0_double(self, histogram):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          t0 = run.get_t0_double()
          # t0 is now a float with the t0 bin 
        """
        return self._histo_dictionary['Time Zero Bin '+str(histogram)]

    def get_sample(self):
        """
        usage:: 
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          samplename = run.get_sample()
          # samplename is now a string with the sample name

        """
        return self._header_dictionary['Sample Name']

    def get_orient(self):
        """
        usage:: 
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          orient = run.get_orient()
          # orient is now a string with the sample orientation info

        """
        if 'Sample Orientation' in self._header_dictionary.keys():
            return self._header_dictionary['Sample Orientation']
        else:
            return 'n/a'

    def get_field_str(self):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          fieldstr = run.get_field()
          # fieldstr is now a string with the field value 
        """
        return self._header_dictionary['Sample Magnetic Field string']

    def get_field(self):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          fieldstr = run.get_field()
          # fieldstr is now a float field value 
        """
        return self._header_dictionary['Sample Magnetic Field']

    def get_timeTemperature_vector(self):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          t,T = run.get_timeTemperature_vector()
          # a T plot is obtained as
          from matplotlib.pyplot import subplots,show
          from matplotlib import dates, units
          import datetime as DT
          converter = dates.ConciseDateConverter()
          units.registry[np.datetime64] = converter
          units.registry[DT.date] = converter
          units.registry[DT.datetime] = converter
          f,a = subplots()
          a.plot(t,T)
          show()          
        """
        from datetime import datetime as DT
        from numpy import arange
        t_format = "%Y-%m-%d %H:%M:%S"
        delta = (DT.strptime(self.get_timeStop(),t_format) - DT.strptime(self.get_timeStart(),t_format))/np.timedelta64(1,'s')
        T = self._plots()[3].counts() if self.get_instrument()=='LEM' else self._plots()[1].counts()
        ndelta = int(delta/T.size)
        t = arange(self.get_timeStart(),self.get_timeStop(),dtype='datetime64['+str(ndelta)+'s]')
        n = t.size-T.size
        if n>0:
            t = t[:-n]
        else:
            T = T[:n]
        return t,T
        
    def get_nameHisto(self,histogram):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          # histogram in range(self.get_numberHisto_int()*len(self.get_RedGreen_offsets())-1
          name = run.get_nameHisto(histogram)
          # name is now a string descriptor of the detector
        """
        if histogram<0 or histogram>=self.get_numberHisto_int()*len(self.get_RedGreen_offsets()):
            return 'Histogram {} does not exist. Max histogram is {}'.format(histogram,self.get_numberHisto_int()*len(self.get_RedGreen_offsets())-1)
        else:
            return self._histo_dictionary['Name '+str(histogram)]
   
    def get_histoNames_vector(self):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          # histogram in range(self.get_numberHisto_int()*len(self.get_RedGreen_offsets())-1
          name = run.get_histoNames_vector()
          # name is now a list of string descriptors of the detectors
        """
        return [self._histo_dictionary['Name '+str(k)] for k in range(self.get_numberHisto_int()*len(self.get_RedGreen_offsets()))]
   
    def get_numberTemperature_int(self):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          num = run.get_numberTemperature_int()
          # num is now the number of temperature sensors, 1 for root files and 4 for bin files
          # use self.Show() lines 127 ... for more info on sensors
        """
        return 1
    
    def get_temperatures_vector(self):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          T = run.get_temperatures_vector()
          # T is now a list of temperature sensor values, 1 for root files and 4 for bin files (only first 2 non zero)
        """
        return [self._header_dictionary['Sample Temperature']]
    
    def get_devTemperatures_vector(self):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          
          T = run.get_devTemperatures_vector()
          # T is now a list of temperature sensor std values, 1 for root files and 4 for bin files (only first 2 non zero)
        """
        return [self._temperature_error_float()]
    
    def get_timeStart_vector(self):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          timeStart = run.get_timeStart_vector()
          # timeStart is the start run yyyy-mm-dd hh:mm:ss string 

        """
        return self._header_dictionary["Run Start Time"].strip().split()
      
    def get_timeStop_vector(self):
        """
        usage::
          from muroot2py import muroot2py as muld
          path2file = 'path and filename'
          run = muld(path2file)  # this is the class instance initiated on the data file
          timeStop = run.get_timeStop_vector()
          # timeStop is the stop run yyyy-mm-dd hh:mm:ss string

        """
        return self._header_dictionary["Run Stop Time"].strip().split()

# commented methods below work for lem only
#    def _filename(self):
#        return str(self._extract(6,'File Name: '))

#    def _runtitle(self):
#        return str(self._extract(7,'Run Title: '))

#    def _runnumber(self):
#        return str(self._extract(8,'Run Number: '))

#    def _starttime(self):
#        return str(self._extract(9,'Start Time: '))

#    def _stoptime(self):
#        return str(self._extract(10,'Stop Time: '))

#    def _instrument(self):
#        return str(self._extract(13,'Instrument: '))

#    def _setup(self):
#        return str(self._extract(17,'Setup: '))

#    def _sample(self):
#        return str(self._extract(19,'Sample Name: '))

#    def _temperature_str(self):
#        return str(self._extract(20,'Sample Temperature: '))

#    def _temperature_float(self):
#        return float(self._temperature_str().partition(' +- ')[0])

#    def _temperature_unit_str(self):
#        return self._temperature_str().partition(' +- ')[2].partition(' ')[2]

#    def _magfield_str(self):
#        return str(self._extract(21,'Sample Magnetic Field: '))

#    def _magfield_float(self):
#        return float(self._magfield_str().partition(' +- ')[0])

#    def _magfield_error_float(self):
#        return float(self._magfield_str().partition(' +- ')[2].partition(' ')[0])

#    def _magfield_unit_str(self):
#        return self._magfield_str().partition(' +- ')[2].partition(' ')[2]

#    def _nhisto_int(self):
#        return int(str(self._extract(22,'No of Histos: ')).partition[2].partition(' -@'))

#    def _tres_nsec(self):
#        return float(str(self._extract(23,'Time Resolution: ')).partition(' ns')[0])

#    def _redgreen_offset(self):
#        return np.array([int(x) for x in str(self._extract(24,'RedGreen Offsets: ')).partition('RedGreen Offsets: ')[2].partition(' -@')[0].split(';')])

#    def _trasport_energy_keV(self):
#        return float(self._extract(27,'Moderator HV: ').partition(' keV')[0])

#    def _sample_keV(self):
#        return float(self._extract(28,'Sample HV: ').partition(' keV')[0])

#    def _implant_energy_keV(self):
#        return float(self._extract(29,'Implantation Energy: ').partition(' keV')[0])

#    def _spin_angle_deg(self):
#        return float(self._extract(30,'Muon Spin Angle: ').partition(' degree')[0])
       
#    def _detector_label(self,histogram):
#        return str(self._runheader()[1][histogram][0]).partition('Name: ')[2].partition(' -@')[0]

#    def _detector_number(self,histogram):
#        return int(str(self._runheader()[1][histogram][1]).partition('Number: ')[2].partition(' -@')[0])

#    def _detector_t0(self,histogram):

#    def _proposal_number(self): 
#    def _main_proposer(self): # only lem 
       
if __name__ == '__main__':
    """
    usage::
      from muroot2py import muroot2py as muld
      path2file = 'path and filename'
      run = muld()  # this is the instance
      
      # reads many things from a PSI root data file 

    """
    m2p = muroot2py()
    m2p.read("../../data/root/lem24_his_0662.root")
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

