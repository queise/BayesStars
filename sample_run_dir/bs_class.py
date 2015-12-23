import time
import subprocess
import math as m
from bs_func2 import *
from bs_obsdata import *

class cl_stmod(object):

    def __init__(self, mass, metal, age):
        self.mass  = mass
        self.metal = metal
        self.age   = age
        self.modelname = ('M%.3f-Z%03d-A%04.0f') % (mass, int(self.metal*10000.), age)
        self.eosopal = ('eos_opal_%03d.bin') % (int(round(self.metal,3)*10000.))
        self.Y0 = 0.24 + 3*self.metal
        self.X0 = 1. - self.metal - self.Y0
        self.path2exe = '' # the executables dmp2k.out and SepCalc must be here 
        self.localdir = '/data/BayesStars/'  # '/local/tmp/jordi/allst/'
        self.path2out = self.localdir + 'outputs/'
        self.path2inp = self.localdir + 'inputs/'
        self.path2def = '/data/BayesStars/default_input_files/'
#        self.file_don = self.path2exe + self.modelname + '.don' # not used, it always read in self.path2exe + 'sun.don'
        self.file_red = self.path2inp + self.modelname + '.redistrb.in'
        self.file_adi = self.path2inp + self.modelname + '.adipls.c.in'
        self.file_fin = self.path2out + self.modelname + '.fin'
        self.file_osc = self.path2out + self.modelname + '.osc'
        self.file_gsm = self.path2out + self.modelname + '.gsm'
        self.file_fre = self.path2out + self.modelname + '.freqs'
        self.file_tot = 'Summary.dat'
        self.file_run = self.path2out + self.modelname + '.run'
        self.OK = True
        self.Chi = -np.inf
        self.Like = -np.inf
        self.logg = -np.inf ; self.teff = -np.inf ; self.zsx = -np.inf ; self.m_ls012 = -np.inf ; self.m_ss02 = -np.inf
        self.model_chars_loaded = False

        # observational data (see bs_obsdata.py):
        self.obs_data = load_observed_data()

        if not os.path.exists(self.path2out):
            os.makedirs(self.path2out)
        if not os.path.exists(self.path2inp):
            os.makedirs(self.path2inp)


    def run_st_evolution(self):

        # Runs dmp2k code:
        print '--- Running dmp2k.out...'

        with open(self.file_run, 'w') as f:
            f.write('this model is running...\n mass  metal  X0  Y0  age\n %.3e  %.3e  %.3e  %.3e  %.3e \n' %
                    (self.mass, self.metal, self.X0, self.Y0, self.age))

        subprocess.Popen("./dmp2k.out", stdin=subprocess.PIPE, shell=True).communicate(
            input='2\n'+'n\n'+self.path2def+'m010.zams\n'+self.modelname+'\n'+
            self.path2exe+'\n'+str(self.mass)+'\n'+str(self.X0)+'\n'+str(self.Y0)+'\n'+
            self.eosopal+'\n'+str(self.age)+'\n')[0]
        print '*** DMP run ended on '+time.strftime("%d %b %Y %H:%M:%S",time.localtime())+'.'

#       os.remove(self.file_run)

	    # Checks if the code ran without errors:
        if not os.path.isfile(self.path2out + self.modelname + '.ok'):
            self.OK = False
        else:
            self.load_model_chars() # loads chars of stellar model from modelname.fin

        # Cleans useless files:
        for extension in ['.lis', '.rad']: # .rcz, .age
            try:
                os.remove(self.path2out + self.modelname + extension)
            except OSError:
                pass


    def run_oscillations(self):
        # Prepares input file *.redistrb.in:
        shutil.copyfile(self.path2def + 'redistrb.in', self.file_red)
        replace_infile(self.file_red, 'MODELNAME', self.path2out + self.modelname)
        
        # Converts .osc.for ascii file to binary .osc:
        subprocess.call("form-amdl.d 2 " + self.file_osc + ".for " + self.file_osc, shell=True)
        print 'conversion .osc done'
        
        # Redistributes the mesh (creates p4 file)
        subprocess.call("redistrb.c.d " + self.file_red, shell=True)
        print 'redistribution done'
    
        # Prepares adipls.c.in file:
        shutil.copyfile(self.path2def + 'adipls.c.in', self.file_adi)
        replace_infile(self.file_adi, 'MODELNAME', self.path2out + self.modelname)
        replace_infile(self.file_adi, 'DEFINPUTSDIR', self.path2def[:-1])
        
        # Calculates the oscillation frequencies:
        subprocess.call("adipls.c.d " + self.file_adi, shell=True)
        
        # Checks if the code ran without errors:
        if not is_non_zero_file(self.file_gsm):
            self.OK = False
            return

        # Converts .gsm binary to ascii .gsm.for, saves output to .freqs file and cleans file:
        f = open(self.file_fre, 'w') ; ftemp = open(self.file_fre + ".tmp", 'w')
        subprocess.call("form-agsm.d 1 " + self.file_gsm + " " + self.file_gsm + ".for", shell=True, stdout=f)
        subprocess.call("grep -A1000 ' 1     0' " + self.file_fre, shell=True, stdout=ftemp)
        f.close() ; ftemp.close()
        shutil.move(self.file_fre + ".tmp", self.file_fre)

        # Cleans useless files:
        for extension in ['.gsm', '.gsm.for', '.osc', '.p4']:
            try:
                os.remove(self.path2out + self.modelname + extension)
            except OSError:
                pass


    def calc_seismicparameters(self):

        # load frequencies of stellar model:
        self.freq = load_model_freqs(self)

        # calculate large and small separations:
        self.m_ls012, self.m_ss02 = calc_separations(self)

        # calculate ratios 010 and slope within range:
        self.sl_r010 = calc_ratios_and_slopes(self)


    def calc_seismicparameters_with_SepCalc(self):
        """
        *****************************************************
        ** NOT USED IN PRESENT CONFIGURATION OF BayesStars **
        *****************************************************
        The seismic parameters and the chi-squareds are calculated by the fortran program SepCalc.
        Program available at:
        https://github.com/queise/SepCalc
        """
        # Prepares input.dat file:
        shutil.copyfile(self.path2def + 'input.dat', self.path2exe + 'input.dat')
        replace_infile(self.path2exe + 'input.dat', 'MODELNAME', self.path2out + self.modelname)
        replace_infile(self.path2exe + 'input.dat', 'STARNAME', self.path2def + 'Dushera')
        
        # Runs fortan code SepCalc:
        subprocess.call("./SepCalc", shell=True)
        
        # Reads relevant results:
        with open(self.path2out + self.modelname + '.Chi2') as f:
            self.Chi, self.Factor = [float(x) for x in f.readline().split()]


    def calc_Likelihood(self):

        # Calculates Chi and the Factor that multiplies in the Likelihood:
        self.calc_Chi_and_Factor()

        # Calculates Likelihood:
        self.Like = m.log(self.Factor) - (self.Chi/2)
        if m.isnan(self.Like):
            self.Like = -np.inf


    def load_model_chars(self):
        # The results from the stellar modelling are loaded from the file self.file_fin created by the dmp2k code.
        # If you want to pass different parameters to BayesStars, edit and recompile dmp2k.
        with open(self.file_fin, 'r') as f:
            age, rad, lum, self.teff, self.zsx, self.logg, rcz =  [float(x) for x in f.readline().strip().split(' ') if x]
        self.model_chars_loaded = True


    def calc_Chi_and_Factor(self):

        if not self.model_chars_loaded: self.load_model_chars()

        num_params_Chi = 5

        # non-seismic parameters
        obs_logg, obs_std_logg, obs_teff, obs_std_teff, obs_ZXs, obs_std_ZXs = self.obs_data[:6]
        self.Chi = (( obs_logg - self.logg ) / obs_std_logg )**2 + (( obs_teff - self.teff ) / obs_std_teff )**2 + \
                   (( obs_ZXs - self.zsx ) / obs_std_ZXs )**2

        # seismic parameters:
        obs_large_sep, obs_std_large_sep, obs_small_sep, obs_std_small_sep = self.obs_data[8:12]
        self.Chi += (( obs_large_sep - self.m_ls012) / obs_std_large_sep )**2 + \
                    (( obs_small_sep - self.m_ss02) / obs_std_small_sep )**2

        self.Chi = m.sqrt(self.Chi) / num_params_Chi

        # factor for Likelihood
        self.Factor = ((2.*m.pi)**(-num_params_Chi/2.)) / (obs_std_logg * obs_std_teff * obs_std_ZXs *
                                                           obs_std_large_sep * obs_std_small_sep)


    def print_results(self, Errflag=''):
        # Print to the output file self.file_tot the results from the modelling (seismic and non-seismic) as well as
        # the Chi and Likelihood:
        with open(self.file_tot, 'a') as f:
            f.write('%.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %s \n' % (self.mass, self.metal,
                self.age, self.logg, self.teff, self.zsx, self.m_ls012, self.m_ss02, self.Chi, self.Like, Errflag))