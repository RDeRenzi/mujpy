from uproot import ReadOnlyFile
from uproot.dynamic import Model_TFolder_v1
from uproot.reading import ReadOnlyDirectory
import numpy as np

class muroot2py():                        # defines the python class
    """
    the user method names are standardized, when possible, to those of MuSR_td_PSI_bin() by Amato Raselli
    
    .. code::

      from mujpy.muroot2py import muroot2py as muld
      run(path_and_filename) # this is the loaded instance 

        reads: 
        lem 2006- lem06_his_nnnn.root
        ins 2023-202x deltat_tdc_ins_yyyy_nnnn.root where ins:x = dolly:4, vms:-, hifi:-
        gps 2023- deltatbc_tdc_gps_yyyy_nnnn.root  <-^ all these are TFolders type 
        flame 26- flameyy_his_nnnn.root            <-- this is TDirectory tye

    .. temporary methods description::

        self._runheader()       returns header (both types) in root format
        self._oldroot()         True for Tfolders, False for TDirectory (flame26)
        self._info_header()     self._runheader() as list of strings
        self._extract(k)        the k-th header string
        self._header_key(k)     extracts the key from the k-th string 
        self._header_value(k)   extracts the value from the k-th string so that
                                    self._header_dict[self._header_key(k)]
                                                         is self._header_value(k)
        self
        self._histos()          returns histos (both types) root tree

    header type 5 is ReadGreen Offset: 0 -@, treated as -@1
    header type 3 Physical quantity treated as str up to an optional ";"
    """

    def __init__(self, path_filename):
        self.readingOOK = self.read(path_filename)  # reads data
        self._make_header_dict() # creates header dictionary self._header_dict
        self._make_histo_dicts()  # creates histo dictionary self._histo_dicts
                                  # used by all get_ calls
    def read(self,path_filename):
        """
        read
        
        .. code::

          run = muld(path2file)  # this is the loaded method instance 
          initiated on the data file
          returns True (success) or False (insuccess)
        """    
        try:
            self.run = ReadOnlyFile(path_filename).root_directory # universal, Tfolder_v1 and TDirectory 
            read_ok = True
        except:
            read_ok = False
        return read_ok

# internal methods

    def _runheader(self):
        """header (v1 and TDictionary)"""

        if self._oldroot():
            header = self.run['RunHeader'].members["fFolders"][0]
        else: # ReadOnlyDirectory
            header = self.run["RunHeader/RunInfo"].items()
        return header 

    def _oldroot(self):
        """ self.run.is_instance('TFolder')"""

        # _newroot = 'classnames' in self.run['RunHeader'].__dir__() 
        return 'members' in self.run['RunHeader'].__dir__()

    def _info_header(self):
        """list of strings extractor (v1 and TDictionary)"""

        return [str(x) for x in self._runheader()] if self._oldroot() else [str(x) for x,_ in self._runheader()] 

    def _extract(self,k):
        """single header string extractor (v1 and TDictionary)"""

        return self._info_header()[k]

    def _header_key(self,k):
        """key for uproot to header dict(v1 and TDictionary)"""

        # string = str(self._runheader()[0][k])
        return self._extract(k).split(':')[0].split('- ')[1]
        
    def _header_type(self,k):
        """deprecated"""

        return int(self._extract(k)[-1])
        
    def _header_value(self,k):
        """value for uproot to header dict(v1 and TDictionary)"""
        jtyp = self._extract(k).split(': ')[-1].split(' -@')[1].split(';')[0]
        val = self._extract(k).split(': ')[1].split(' -@')[0]
        return val if jtyp in ['0','3','5'] else float(val) if (jtyp=='2' or '.' in val) else int(val)

    def _make_header_dict(self):
        """uproot to header dict (v1 and TDictionary)"""

        self._header_dict = {}
        for k in range(len(self._info_header())):
            self._header_dict[self._header_key(k)] = self._header_value(k)
            #if self._header_key(k) == 'Sample Temperature' or self._header_key(k) == 'Sample Magnetic Field':
            #    header = str(self._runheader()[0][k])
            #    self._header_dict[self._header_key(k)+' string'] = header.partition(self._header_key(k)+': ')[2].partition(' -@')[0]

    def _histo_info_list(self,k_histo):
        """list of strings for histo name, number, length, t0 bin, first good, last good, [AddedHistos]"""

        if self._oldroot():
            info = self.run['RunHeader'].members["fFolders"][1][k_histo]
            info_list = [str(x) for x in info]
        else: # ReadOnlyDirectory
            histo = "RunHeader/DetectorInfo/Detector00"+str(k_histo+1)
            info = self.run[histo].items()
            info_list = [str(x) for x,_ in info]
        return info_list

    def _histo_key(self,k_histo,k):
        """key for uproot to header dict(v1 and TDictionary)"""

        # string = str(self._runheader()[0][k])
        return self._histo_info_list(k_histo)[k].split(':')[0].split('- ')[1]

    def _histo_values(self,k_histo):
        """values for uproot to histo _dicts (v1 and TDictionary)"""
        vals = []
        for string in self._histo_info_list(k_histo):
            val = string.split(': ')[-1].split(' -@')
            jtyp, val = val[1].split(';')[0], val[0].split(';')[0]
            val = val if jtyp=='0' or jtyp=='3' else float(val) if (jtyp=='2' or '.' in val) else int(val)
            vals.append(val)
        return vals

    def _make_histo_dicts(self):
        """uproot to _histo_dicts list of dict"""

        number_histos = self._header_dict['No of Histos']
        self._histo_dicts = [
                {self._histo_key(k_histo,k):val for k,val in enumerate(self._histo_values(k_histo))} for 
                k_histo in range(number_histos)]

    def _histo_number(histo_number):
        """from Histo Number, e.g. lem 21, to 0-based k in k-th TH1F, e.g. lem 8"""

        number_histos = self._header_dict('No of Histos')
        offsets = self.get_RedGreen_offsets().split('; ')
        index = [histo_number<int(o) for o in offsets].index(True)-1
        offset = int(offsets[index])
        return histo_number%offset-1+histo_number//offset*number_histos

    def _histo(self,k):
        """np.array of counts from k-th TH1F vector"""

        if self._oldroot():
            histo= self.run['histos'].members["fFolders"][0].members['fFolders'][k]
        else: # ReadOnlyDirectory
            histo = self.run["histos/DecayAnaModule/hDecay00"+str(k+1)]
        # uproots arrays ARE numpy arrays
        return histo.counts()

    def _plots(self):
        """deprecated?"""

        if self._header_dict['Instrument']=='LEM':
            return self.run["histos"].members['fFolders'][5].members['fFolders'] 
        else:
            return self.run["histos"].members['fFolders'][1].members['fFolders'] 
    
    def _rbn(self,hist,binning):
        """assumes hist is a numpy array, as from self._histo(k)""" 
        return hist[:(hist.size//binning)*binning].reshape(hist.size//binning,binning).sum(axis=1)
        
# user methods

    def Filename(self):
        """header filename"""

        return self._header_dict['File Name']
        
    def readingOK(self):
        """True if instance loaded the data file"""

        return self.readingOOK

    def Show(self):
        """print RunHeader, fFolders"""

        for x in self._info_header(): print(' '+x)
                
    def get_binWidth_ns(self):
        """bin width in ns"""
        
        return float(self._header_dict['Time Resolution'].split(' ns')[0])

    def get_binWidth_us(self):
        """bin width in us"""

        return self.get_binWidth_ns/1000.

    def get_histoLength_bin(self):
        """histogram0 length in bins"""

        return int(self._histo_dicts[0]['Histo Length'])

    def get_numberHisto_int(self):
        """number of histograms"""

        return self._header_dict['No of Histos']
    
    def get_RedGreen_offsets(self):
        """list of offsets (Red-Green mode)"""

        str_offsets = self._header_dict['RedGreen Offsets'].split('; ')
        #for line in self._info_header(): 
        #    if 'RedGreen Offsets' in line: 
        #        str_offsets = line.split(': ')[1].split(' -@')[0].split('; ')
        #        break
        return [int(x) for x in str_offsets]
    
    def get_histo_int(self,k_histo,k):
        """counts in k_histo-th detector histo(k)"""

        return int(self._histo(k_histo)[k])

    def get_histo_vector(self,k_histo,binning=1):
        """np.array of counts for k_histo-th detector""" 

        return self._histo(k_histo) if binning==1 else self._rbn(self,_histo(k_histo),binning)
    
    def get_histo_vector_no0(self,k_histo,binning=1):
        """np.array of counts for k_histo-th detector (0.1 for 0)"""

        hist = self.get_histo_vector(k_histo,binning=binning)
        return np.where(hist != 0, hist, 0.1)
    
    def get_histo_fromt0_vector(self,k_histo,binning=1,offset=0):
        """np.array of counts for k_histo-th detector from t0 bin"""

        return self.get_histo_vector(k_histo,binning=binning)[self.get_t0_int(k_histo)+offset:]

    def get_histo_fromt0_minus_bckgrd_vector(self,k_histo,first_bkg,last_bkg,binning=1,offset=0):
        """np.array of counts for k_histo-th detectori from t0 bin with background subtraction from t<t0 mean"""

        hist = self.get_vector(k_histo)
        background = hist[first_bkg:last_bkg+1].sum()/(last_bkg+1-first_bkg)
        hist = self.get_histo_fromt0_vector(k_histo,binning=binning,offset=offset)
        return hist-background*binning   
    
    def get_instrument(self):
        """instrument from run header"""

        return self._header_dict['Instrument'].strip()

    def get_title(self):
        """title from run header"""

        return self._header_dict['Run Title']

    def get_runNumber_int(self):
        """integer run number"""

        return int(self._header_dict['Run Number'])

    def get_beamline(self):
        """beamline from run header"""

        return self._header_dict['Muon Source']

    def get_transport_energy(self):
        """transport energy from header (LEM)"""

        if self._header_dict['Instrument']=='LEM':
            return self._header_dict['Moderator HV']
        else: return "" 

    def get_sample_HV(self):
        """sample HV from header (LEM)"""

        if self._header_dict['Instrument']=='LEM':
           return self._header_dict['Sample HV']
        else: return "" 

    def get_implantation_energy(self):
        """
        Implantation Energy from header (LEM)
        """

        if self._header_dict['Instrument']=='LEM':
            return self._header_dict['Implantation Energy']
        else: return "" 

    def get_spin_angle(self):
        """spin rotator angle"""

        return self._header_dict['Muon Spin Angle']

    def get_proposal_number(self):
        """Proposal number"""

        return self._header_dict['Proposal Number']
        
    def get_proposer(self):
        """PI (LEM)"""

        if inst=='LEM' or not self._oldroot:
            return self._header_dict['Main Proposer']
        else: return "" 
 
    def get_comment(self):
        """comment"""

        return self._header_dict['Comment']

    def get_temp(self):
        """temperature without error (string)"""

        return self._header_dict['Sample Temperature'].split(' +-')[0]

    def get_temperature(self):
        """temperature (float)"""

        return float(self._header_dict['Sample Temperature'].split(' +-')[0])

    def get_t0_int(self, k_histo):
        """t0 bin from run histo k-histo header, IF you trust it"""  

        return int(self._histo_dicts[k_histo]['Time Zero Bin'])

    def get_t0_double(self, k_histo):
        """t0 float as float"""

        return float(self._histo_dicts['Time Zero Bin '+str(histogram)])

    def get_sample(self):
        """sample"""

        return self._header_dict['Sample Name']

    def get_orient(self):
        """orientation"""

        if 'Sample Orientation' in self._header_dict.keys():
            return self._header_dict['Sample Orientation']
        else:
            return ''

    def get_field_str(self):
        """field value string"""

        return self._header_dict['Sample Magnetic Field'].split(' G')[0].split(' G')[0].split(' +-')[0]

    def get_field(self):
        """field float"""

        return float(self._header_dict['Sample Magnetic Field'].split(' G')[0].split(' +-')[0])

    def get_timeTemperature_vector(self):
        """temperature time log 

        .. code::

          from mujpy.muroot2py import muroot2py as muld
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
        
    def get_nameHisto(self,k_histo):
        """k_histo-th histogram descriptor"""

        return self._histo_dicts['Name '+str(k_histo)]
   
    def get_histoNames_vector(self):
        """list of histogram descriptors"""

        return [self.get_nameHist(k) for k in range(self.get_numberHisto_int()*len(self.get_RedGreen_offsets()))]
   
    def get_numberTemperature_int(self):
        """number of T sensors (fixed to 1 arbitrarily)"""

        return 1
    
    def get_temperatures_vector(self):
        """list of T sensor values"""

        return [float(self._header_dict['Sample Temperature'].split(' +-')[0])]
    
    def get_devTemperatures_vector(self):
        """list of T sensor stds"""

        tempstr = self._header_dict['Sample Temperature']
        eT = 0 if '+-' not in tempstr else float(tempstr.split(' +- ')[1].split(' K')[0])
        return [eT]
    

    def get_timeStart_vector(self):
        """run start time"""

        return self._header_dict["Run Start Time"].strip().split()
      
    def get_timeStop_vector(self):
        """run stop time"""

        return self._header_dict["Run Stop Time"].strip().split()
