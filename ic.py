import argparse as arg
import numpy as np
from numpy.lib.shape_base import hsplit
import pandas as pd
import progressbar as pb

class Population:
    def __init__(self, number_of_binaries, **kwards):

        self.Nbin = number_of_binaries
        self.Nstars = int(2.*self.Nbin)
        # Update paramenters needed for draw distributions
        self.setup_params(self, kwards)
        # Upload the method used to save inputs for our pop-syn
        self.save_mobse_input = self._save_mobse_input
        self.save_sevn_input = self._save_sevn_input
        # Check for IC model
        self.model = None
        if 'model' in kwards:
            self.model = kwards['model']

        # Build dataframe with minimal info
        self.binaries = pd.DataFrame(columns=['m1','m2','p','ecc','a'])

        # Select model (if specified)
        if self.model.lower() == 'sana12':
            print('Building a population of binaries based on Sana+2012 and Kroupa2001')
            self.Sana_etal12(self)
        elif self.model.lower() == 'M&DS17':
            print('Be patiente, we will implement this model soon')
        else:
            print('No model selected (binaries is empty)! Build your own model.')
        
    @classmethod
    def setup_params(cls, self, kwards):
        params = {
            'alphas':[-1.3,-2.3],
            'mass_range':[0.1,0.5,150],
            'P_min': 0.15,
            'P_max': 5.5,
            'P_slope': -0.55,
            'e_min': 0.0,
            'e_max': 0.9999,
            'e_slope': -0.45,
            'q_min': 0.1,
            'q_max': 1.0,
            'q_slope': -0.1,
            'mass_min': 0.1
        }

        for key in params.keys():
            if key in kwards:
                setattr(self, key, kwards[key])
            else:
                setattr(self, key, params[key])

    @classmethod    
    def Sana_etal12(cls, self):
        # The beginning of the main function
        pbar = pb.ProgressBar().start()
        self.binaries['m1'] = self.IMF(self.Nbin, self.mass_range, self.alphas)
        pbar.update((1/6)*100)
        q = self.power_law(self.Nbin, self.q_min, self.q_max, self.q_slope)
        pbar.update((2/6)*100)
        secondary = self.binaries['m1'] * q
        self.binaries['m2'] = np.where(secondary < self.mass_min, self.mass_min, secondary)
        pbar.update((3/6)*100)
        self.binaries['p'] = pow(10.,self.power_law(self.Nbin, self.P_min, self.P_max, self.P_slope))
        pbar.update((4/6)*100)
        self.binaries['ecc'] = self.power_law(self.Nbin, self.e_min, self.e_max, self.e_slope)
        pbar.update((5/6)*100)
        self.binaries['a'] = self.p2a(self.binaries['p'], self.binaries['m1'], self.binaries['m2'])        
        pbar.finish()   
        # The end
    
    # To use pandas dataframe use to define the object
    def _save_mobse_input(self, name, met, tmax, backup=10):
        """
        Input: 
            fname = name of the output file
            met = metallicity (string)
            tmax = max time (float)
        """
        Z = np.full(self.Nbin,float(met))
        t = np.full(self.Nbin,tmax)
        ind = np.arange(1,self.Nbin+1,1)
        m1, m2, p, ecc = np.hsplit(self.binaries[['m1','m2','p','ecc']].to_numpy(),4)
        np.savetxt(name+"_"+met+".in", 
                np.c_[ind, m1, m2, p, ecc, Z, t], 
                fmt=('%i %4.4f %4.4f %10.4f %1.4f %1.4f %6.1f'),
                header=str(self.Nbin-backup), delimiter=' ', comments='')

    def _save_sevn_input(self, name, z1, z2, o1, o2, tend, tstart,\
                        dt, sn1, sn2, dtout):
        """
        Input: 
            z1, z2 = metallicity of the stars 1 and 2 (arrays)
            o1, o2 = spin of the stars 1 and 2 [1/yr] (arrays)
            tend = ending time of the simulation [Myr] (float)
            tstart = starting time of the simulation [Myr] (float)
            dt = time step of the integration [Myr] (not used)
            sn1, sn2 = supernova explosion mechanism of the stars (string: _delayed_ / _rapid_ / _startrack_)
            dtout = time step for printing outputs [Myr] (float)
        """
        np.savetxt(name+".in", 
            np.c_[z1, z2], 
            fmt=('%1.4f %1.4f'),
            delimiter=' ', comments='')

    @staticmethod
    def save_mobse_input(name, m1, m2, p, ecc, met, tmax, backup=10):
        """
        Input: 
            fname = name of the output file
            m1 = primary masses (array)
            m2 = secondary masses (array)
            p = periods (array)
            ecc = eccentricity (array)
            met = metallicity (string)
            tmax = max time (float)
        """
        Nbin = len(m1)
        Z = np.full(Nbin,float(met))
        t = np.full(Nbin,tmax)
        ind = np.arange(1,Nbin+1,1)
        np.savetxt(name+"_"+met+".in", 
                np.c_[ind, m1, m2, p, ecc, Z, t], 
                fmt=('%i %4.4f %4.4f %10.4f %1.4f %1.4f %6.1f'),
                header=str(Nbin-backup), delimiter=' ', comments='')

    @staticmethod
    def save_sevn_input(name, m1, m2, z1, z2, o1, o2, a, e, tend, tstart,\
                        dt, sn1, sn2, dtout):
        """
        Input: 
            fname = name of the output (string)
            m1, m2 = mass of the stars 1 and 2 [Msun] (arrays)
            z1, z2 = metallicity of the stars 1 and 2 (arrays)
            o1, o2 = spin of the stars 1 and 2 [1/yr] (arrays)
            a = binary separation [Rsun] (array)
            e = binary eccentricity (array)
            tend = ending time of the simulation [Myr] (float)
            tstart = starting time of the simulation [Myr] (float)
            dt = time step of the integration [Myr] (not used)
            sn1, sn2 = supernova explosion mechanism of the stars (string: _delayed_ / _rapid_ / _startrack_)
            dtout = time step for printing outputs [Myr] (float)
        """
        np.savetxt(name+".in", 
            np.c_[m1, m2], 
            fmt=('%4.4f %4.4f'),
            delimiter=' ', comments='')
            
    @staticmethod
    def a2p(sep, m1, m2):
        """
        Compute the period (day) given m1 (Msun), m2 (Msun) and the sep (Rsun)
        """
        yeardy=365.24
        AURsun=214.95
        p = ((sep/AURsun)**3./(m1 + m2))**(0.5)
        return p, p*yeardy

    @staticmethod
    def p2a(p, m1, m2):
        """
        Compute the separation (Rsun) given m1 (Msun), m2 (Msun) and p (days)
        """        
        yeardy=365.24
        AURsun=214.95
        p = p/yeardy
        a = AURsun*(p*p*(m1 + m2))**(1./3.)
        return a 
        
    @staticmethod
    def IMF(number_of_stars, mass_bouders=[0.1,0.5,150], alphas=[-1.3,-2.3]):
        """
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

    @staticmethod
    def power_law(N, low, up, eta):
        """
        Input:  N = number of systems
                low = lower limit
                up = upper limit
                eta = exponent of the power law
        """
        X = np.random.random(N)
        Y = pow((pow(up,1.+eta) - pow(low,1.+eta)) * X + \
                pow(low,1.+eta),1./(1.+eta))
        return Y