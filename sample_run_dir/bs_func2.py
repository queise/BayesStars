import os
import shutil
import numpy as np
from scipy import stats

###########################
#    SEISMIC FUNCTIONS
###########################

def load_model_freqs(self):
    """
    Loads the oscillation frequencies of a stellar model from the file self.modelname+'.freqs', filtering out the modes
    that are not within the range of observed frequencies (f_min and f_max defined in bs_obsdata.py)
    This file is created by the Adipls code.

    Inputs
    ------
    self:   is an object of the class cl_stmod (see bs_class.py)

    Returns
    -------
    self.freq: a dictionary which keys are the angular and radial order of the oscillations modes.
               Ex: for the l=0, n=13 mode with freq=0.512 mHz --> self.freq[0][13] = 0.512

    The function also stores in the class cl_stmod other information that is used afterwards when the separations and
    ratios are calculated, like self.l_n, self.l_max, self.n_max, self.n_min.
    """
    print self.modelname
    f_min, f_max = self.obs_data[6:8]

    # we load the observed frequencies:
    self.l_n = np.loadtxt(self.file_fre, usecols = (1,2), dtype=int) #[('f0',int),('f1',int),('f2',float),('f3',float)])
    freqs = np.loadtxt(self.file_fre, usecols = (6,))
    # radial orders (n) min and max for each and degree (l)
    self.l_max = min(np.max([row[0] for row in self.l_n]),3) # maximum l=3
    self.n_max = np.zeros(self.l_max+1, dtype=int)
    self.n_min = np.zeros(self.l_max+1, dtype=int)
    for l in xrange(self.l_max+1):
        all_ns_of_l = [n for ls,n in self.l_n if ls==l]
        self.n_max[l] = np.max(all_ns_of_l)
        self.n_min[l] = np.min(all_ns_of_l)
        #print ("l=%i \t nmin=%i \t nmax=%i " % (l, self.n_min[l], self.n_max[l]))
    # creating dictionaries for freq[l][n] and err[l][n] and declaring others:
    self.freq = {} ; self.largesep = {}
    for l in xrange(self.l_max+1):
        self.freq[l] = {}
        self.largesep[l] = {}
    for (l,n),f in zip(self.l_n,freqs):
        # store freq only if within l=0 and l=lmax:
        if l>self.l_max: break
        # store freq only if freq within observed frequency ranges (defined in bs_obsdata):
        if f>=f_min and f<=f_max:
            self.freq[l][n] = f
        #print ("l=%i \t n=%i \t %.2f" % (l,n,self.freq[l][n]))

    return self.freq


def calc_separations(self):
    """
    Calculates the mean large (for l=0,1,2) and small (for l=0,2) separations for f_min<f<f_max (already filtered
    in self.freq).

    Inputs
    ------
    self:   is an object of the class cl_stmod (see bs_class.py)

    Returns
    --------
    self.m_ls012
    self.m_ss02
    """

    # large separations and mean for l=0,1,2:
    for l,n in self.l_n:
        if l<=self.l_max and n<self.n_max[l] and (n in self.freq[l]) and ((n+1) in self.freq[l]):
            self.largesep[l][n] = self.freq[l][n+1] - self.freq[l][n]
    m_ls = {}
    for l in xrange(3):
        m_ls[l] = np.mean(self.largesep[l].values())
        #print ("m_ls[%i] = %.4f" % (l,m_ls[l]))
    self.m_ls012 = np.mean(m_ls.values())
    print 'm_ls012:',self.m_ls012

    # small separation l=0,2:
    smallsep02 = {}
    for n in xrange(max(self.n_min[0],self.n_min[2]+1), min(self.n_max[0],self.n_max[2]+1)+1):
        if (n in self.freq[0]) and ((n-1) in self.freq[2]):
            smallsep02[n] = self.freq[0][n] - self.freq[2][n-1]
    self.m_ss02 = np.mean(smallsep02.values())
    print 'm_ss02:',self.m_ss02

    return self.m_ls012, self.m_ss02


def calc_ratios_and_slopes(self):
    """
    Calculates the ratios r01, r10, builds r010, and calculates its slope within a frequency range defined by f_min_sl
    and f_max_sl (set in bs_obsdata.py).

    Inputs
    ------
    self:   is an object of the class cl_stmod (see bs_class.py)

    Returns
    -------
    self.sl_r010:   slope of r010 within freq range
    """
    # r01
    r01 = [] ; f_r01 = []
    for n in xrange(max(self.n_min[0]+1,self.n_min[1]+1), min(self.n_max[0]-1,self.n_max[1])+1):
        if (n in self.freq[0]) and (n-1 in self.freq[0]) and (n-1 in self.freq[1]) and \
            (n in self.freq[1]) and (n+1 in self.freq[0]):
            f_r01.append(self.freq[0][n])
            r01.append((self.freq[0][n-1] - 4*self.freq[1][n-1] + 6*self.freq[0][n] - 4*self.freq[1][n] +
                        self.freq[0][n+1]) / (8*(self.freq[1][n]-self.freq[1][n-1])))
    self.m_r01 = np.mean(r01) # .values())
    print 'mean r01:',self.m_r01

    # r10
    r10 = [] ; f_r10 = []
    for n in xrange(max(self.n_min[0],self.n_min[1]+1), min(self.n_max[0]-1,self.n_max[1]-1)+1):
        if (n in self.freq[1]) and (n-1 in self.freq[1]) and (n in self.freq[0]) and \
            (n+1 in self.freq[0]) and (n+1 in self.freq[1]):
            f_r10.append(self.freq[1][n])
            r10.append(-(self.freq[1][n-1] - 4*self.freq[0][n] + 6*self.freq[1][n] - 4*self.freq[0][n+1] +
                         self.freq[1][n+1]) / (8*(self.freq[0][n+1]-self.freq[0][n])))
    self.m_r10 = np.mean(r10) #.values())
    print 'mean r10:',self.m_r10

    # r010
    r010 = [] ; f_r010 = []
    for i in xrange(max(len(r01),len(r10))):
        if i < len(r01):
            f_r010.append(f_r01[i])
            r010.append(r01[i])
        if i < len(r10):
            f_r010.append(f_r10[i])
            r010.append(r10[i])
    # average and slope in particular range:
    f_min_sl, f_max_sl = self.obs_data[12:14]
    print 'f_min_sl, f_max_sl:',f_min_sl, f_max_sl
    f_r010_r = [f_r010[i] for i in xrange(len(r010)) if f_r010[i]>f_min_sl and f_r010[i]<f_max_sl]
    r010_r = [r010[i] for i in xrange(len(r010)) if f_r010[i]>f_min_sl and f_r010[i]<f_max_sl]
    self.m_r010 = np.mean(r010_r)
    print 'm_r010 = ',self.m_r010
    # slope
    self.sl_r010 = stats.linregress(f_r010_r,r010_r)[0]
    print 'sl_r010 = ',self.sl_r010

    return self.sl_r010


###########################
#    OTHER FUNCTIONS
###########################

def replace_infile(file_path, text2search, text2replace):
    temp_path = file_path+'.tmp'
    new_file = open(temp_path,'w')
    old_file = open(file_path)
    for line in old_file:
        new_file.write(line.replace(text2search, text2replace))
    new_file.close()
    old_file.close()
    os.remove(file_path)
    shutil.move(temp_path, file_path)

def is_non_zero_file(fpath):
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False

