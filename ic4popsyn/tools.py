import numpy as np
# mandatory for GetFromUrl 
import subprocess
import os 
import tarfile

class GetFromUrl:
    
    def __init__(self, name="mobse"):
        if name.lower() == 'mobse':
            self.version = "v1.0" #It will be active in the next future
            self.sourcename =  'source-code-'+self.version
            self.url = 'https://gitlab.com/mobse/source-code/-/archive/v1.0/source-code-v1.0.tar.gz'
            self.name = name
        if name.lower() == 'sevn':
            print('WARNIG: SEVN is still not publicly available.')
            self.version = "v1." #It will be active in the next future
            self.sourcename = self.version
            self.url = 'https://github.com/GiacobboNicola/IC4popsyn/archive/refs/tags/v0.1.tar.gz'
            self.name = name


    def directory(self):
        maindir = os.path.abspath(os.path.dirname(__file__))
        subdir = os.path.join(maindir, self.name)
        return subdir
    
    def tar_from_url(self):
        subprocess.run(["curl","-L","-O",self.url])
        #subprocess.run(["wget",self.url])

    def rename_dir(self):
        path_to_save = ''
        if os.path.exists(path_to_save+self.name):
            print(f"WARNING: {self.name} already exists! We are renaming it!")
            counter = 0
            while os.path.exists(path_to_save+self.name+'_n{0}'.format(counter)):
                counter += 1
                if counter > 2: 
                    print("    |----> Be carefull! Too many copy of the folder.")
            os.rename(path_to_save+self.name, path_to_save+self.name+'_n{0}'.format(counter))
        os.rename(self.sourcename, self.name)
        
    def main(self):
        print("-----------------------------------------------------------------------------------")
        print("downloading version", self.version) 
        print("from", self.url) 
        print("to", self.directory())
        print("-----------------------------------------------------------------------------------")

        self.tar_from_url()

        tar = tarfile.open(self.sourcename+'.tar.gz', 'r:gz')
        os.remove(self.sourcename+'.tar.gz')

        tar.extractall()
        tar.close()

        self.rename_dir()

        print("download completed")
        print("-----------------------------------------------------------------------------------\n")


def a2p(sep, m1, m2):
    """
    It computes the period (day) given m1 (Msun), m2 (Msun) and the sep (Rsun).
    """
    yeardy=365.24
    AURsun=214.95
    p = ((sep/AURsun)**3./(m1 + m2))**(0.5)
    return p, p*yeardy

def p2a(p, m1, m2):
    """
    It computes the separation (Rsun) given m1 (Msun), m2 (Msun) and p (days).
    """        
    yeardy=365.24
    AURsun=214.95
    p = p/yeardy
    a = AURsun*(p*p*(m1 + m2))**(1./3.)
    return a 
    
def IMF(number_of_stars, mass_bouders=[0.1,0.5,150], alphas=[-1.3,-2.3]):
    """
    It samples stellar masses from a broken-power-law.
    Input:  
        Number of stars
        mass_bouders=[0.1,0.5,150]
        alphas=[-1.3,-2.3]
    Default: 
        f(m) = m^alpha (Kroupa01)
        0.1 <= m <= 0.5   alpha = -1.3
        0.5 <= m <= 150   alpha = -2.3
    """
    random = np.random.random

    breaks = np.array(mass_bouders)
    slopes = np.array(alphas)

    bins = len(slopes)

    fraction_per_bin, stars_per_bin = [], []

    # it computes the fractions of stars per each mass range definded by the breaks 
    for i,slope in enumerate(slopes):
        if slope == -1:
            factor = np.log(breaks[i+1] / breaks[i])
        else:
            factor = (breaks[i+1]**(slope+1) - breaks[i]**(slope+1)) / (slope+1)

        j=0
        for j in range(bins - i -1):
            factor *= breaks[-j-2]**(slopes[-j-1] - slopes[-j-2])
        stars_per_bin.append(factor)

    total = sum(stars_per_bin,0.0)
    fraction_per_bin = np.array(stars_per_bin)/total

    cumulative_fractions = np.array([sum(fraction_per_bin[:i]) for i in range(bins+1)])
    cumulative_fractions[-1] = 1.0 # to prevent problem with the rounding
    kfacs = pow(breaks[1:] / breaks[:-1], slopes + 1.0) - 1.0 # = (slope + 1) / K 
    inv_slopes_plus1 = np.array([np.inf if slope==-1 else (1.0 / (slope + 1.0)) for slope in slopes])

    # smart way to compute masses (see https://amusecode.github.io/ and doc/imf.pdf for details) 
    random_numbers = random(number_of_stars) #generate number_of_stars ramdom values in [0.1] 
    # for each ran. num. we save in which 'fraction' it belongs
    indices = np.searchsorted(cumulative_fractions[:-1], random_numbers, 'right') - 1
    # to have the ralative 'position' of the ran. num. inside the mass range (e.i. between the breaks)
    scaled = ((random_numbers - cumulative_fractions[indices]) /
        (cumulative_fractions[indices+1] - cumulative_fractions[indices])) 
    result = np.empty_like(random_numbers)
    zerodiv = slopes[indices]==-1
    normal = np.logical_not(zerodiv)
    result[zerodiv] = pow(breaks[1:][indices[zerodiv]] / breaks[:-1][indices[zerodiv]], scaled[zerodiv])
    result[normal] = pow(1.0 + kfacs[indices[normal]] * scaled[normal], inv_slopes_plus1[indices[normal]])

    return breaks[:-1][indices] * result

def power_law(N, low, up, eta):
    """
    It samples elements from a power-law distribution.
    Input:  N = number of systems
            low = lower limit
            up = upper limit
            eta = exponent of the power law
    """
    X = np.random.random(N)
    Y = pow((pow(up,1.+eta) - pow(low,1.+eta)) * X + \
            pow(low,1.+eta),1./(1.+eta))
    return Y

def vec_power_law(low, up, eta):
    _vec_power_law = np.vectorize(power_law)
    return _vec_power_law(1,low,up,eta)

def eccvsP(P):
    """
    Eccentricity as a function of the period (eq. 3 in M&DS2017).
    P must be in days.
    """
    return 1 - (P / 2.)**(-2./3.) # eq.3 in M&DS17