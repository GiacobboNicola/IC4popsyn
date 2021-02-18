import argparse as arg
import numpy as np
import pandas as pd
import progressbar as pb

class Population:
# add model to choose between Moe di Stefeno and Sana
    def __init__(self, number_of_binaries, model, **kwards):

        self.Nbin = number_of_binaries
        self.Nstars = int(2.*self.Nbin)
        self.model = model
        self.setup_params(self, kwards)
        
        self.binaries = pd.DataFrame(columns=['primary','secondary','period','ecc','separation'])

        if self.model == 'sana':
            self.Sana_etal12()
        elif self.model == 'moedistefano':
            print('We will implent this model soon')
        elif self.model == None:
            print('Binaries are ready')

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
    
    def Sana_etal12(self):
        # The beginning of the main function
        pbar = pb.ProgressBar().start()
        self.binaries['primary'] = self.IMF(self.Nbin, self.mass_range, self.alphas)
        pbar.update((1/5)*100)
        q = self.mass_ratio(self.Nbin, self.q_min, self.q_max, self.q_slope)
        pbar.update((2/5)*100)
        secondary = self.binaries['primary'] * q
        self.binaries['secondary'] = np.where(secondary < self.mass_min, self.mass_min, secondary)
        pbar.update((3/5)*100)
        self.binaries['period'] = self.period(self.Nbin, self.P_min, self.P_max, self.P_slope)
        pbar.update((4/5)*100)
        self.binaries['ecc'] = self.eccentricity(self.Nbin, self.e_min, self.e_max, self.e_slope)
        pbar.finish()   
        # The end

    @staticmethod
    def IMF(number_of_stars, mass_bouders=[0.1,0.5,150], alphas=[-1.3,-2.3]):
        """
        default: f(m) = m^alpha (Kroupa01)
                0.1 <= m <= 0.5   alpha = -1.3
                0.5 <= m <= 150   alpha = -2.3
        """
        random = np.random.random

        breaks = np.array(mass_bouders)
        slopes = np.array(alphas)

        mass_max = breaks[-1]

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
    def period(Nbin, lower_lim=0.15, upper_lim=3.5, pi=-0.55):
        """
        default: f(logP) = logP^pi (Sana+12)
                lower_lim = 0.15
                upper_lim = 5.5
                pi = -0.55
        """
        X = np.random.random(Nbin)
        logP = pow((pow(upper_lim,1.+pi) - pow(lower_lim,1.+pi)) * X + \
                pow(lower_lim,1.+pi),1./(1.+pi))
        period = pow(10.,logP)
        return period

    @staticmethod
    def eccentricity(Nbin, lower_lim=0.0, upper_lim=0.9999, eta=-0.45):
        """
        default: f(e) = e^eta (Sana+12)
                lower_lim = 0.0
                upper_lim = 0.9999
                eta = -0,45 
        """
        X = np.random.random(Nbin)
        ecc = pow((pow(upper_lim,1.+eta) - pow(lower_lim,1.+eta)) * X + \
                pow(lower_lim,1.+eta),1./(1.+eta))
        return ecc

    @staticmethod
    def mass_ratio(Nbin, lower_lim=0.1, upper_lim=1., kappa=-0.1):
        """
        default: f(q) = q^kappa (Sana+12)
                lower_lim = 0.1
                upper_lim = 1.
                kappa = -0.1 
        """
        X = np.random.random(Nbin)
        q = pow((pow(upper_lim,1.+kappa) - pow(lower_lim,1.+kappa)) * X + \
                pow(lower_lim,1.+kappa),1./(1.+kappa))
        return q

#class Constants:
