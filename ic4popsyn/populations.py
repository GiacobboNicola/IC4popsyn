"""
Classes to build populations of stars (singles/bineries/triples)
"""
import argparse as arg
import numpy as np
from numpy.lib.shape_base import hsplit
import pandas as pd
import progressbar as pb
# this package and necessary for Binaries
from ic4popsyn import tools

class Binaries:
    def __init__(self, number_of_binaries, single_pop=False, **kwards):

        self.Nbin = number_of_binaries
        self.Nstars = int(2.*self.Nbin)
        self.single_pop = single_pop
        # Update paramenters needed for draw distributions
        self.setup_params(self, kwards)
        # Upload the method used to save inputs for our pop-syn
        self.save_mobse_input = self._save_mobse_input
        self.save_sevn_input = self._save_sevn_input
        # Check for IC model
        if 'model' in kwards:
            self.model = kwards['model'].lower()
        else:
            self.model = 'None'
        
        # Build dataframe with minimal info
        if single_pop:
            self.population = pd.DataFrame(columns=['m1'])
        else:
            self.population = pd.DataFrame(columns=['m1','m2','p','ecc','a'])

        # Select model (if specified)
        if self.model in ['sana12','sana_eccm&ds']:
            print('Building a population of binaries based on Sana+2012 and Kroupa2001')
            if self.model == 'sana_eccm&ds':
                print('Building a population of binaries based on Sana+2012, Kroupa2001 and Moe&DiStefano 2017')
            self.Sana_etal12(self)
        elif self.model == 'M&DS17':
            print('Be patient, we will implement this model soon')
        else:
            print('No model selected (binaries is empty)! Build your own model.')

    @classmethod
    def setup_params(cls, self, kwards):
        params = {
            'alphas':[-1.3,-2.3],
            'mass_range':[0.1,0.5,150],
            'logP_min': 0.15,
            'logP_max': 5.5,
            'logP_slope': -0.55,
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
        """
        It draws a population of binaries based on Sana+12.
        """
        # The beginning of the main function
        pbar = pb.ProgressBar().start()
        self.population['m1'] = tools.IMF(self.Nbin, self.mass_range, self.alphas)
        
        if self.single_pop: #if True returns only population['m1']
            pbar.finish()
            return

        pbar.update((1/6)*100)
        q = tools.power_law(self.Nbin, self.q_min, self.q_max, self.q_slope)
        pbar.update((2/6)*100)
        secondary = self.population['m1'] * q
        self.population['m2'] = np.where(secondary < self.mass_min, self.mass_min, secondary)
        pbar.update((3/6)*100)
        self.population['p'] = pow(10.,tools.power_law(self.Nbin, self.logP_min, self.logP_max, self.logP_slope))
        pbar.update((4/6)*100)
        if self.model == 'sana_eccm&ds':
            Pecc = np.where(self.population['p'] < 2., 2.1, self.population['p'])
            eccMax = tools.eccvsP(Pecc)
            self.population['ecc'] = tools.vec_power_law(self.e_min, eccMax, self.e_slope)
        else:
            self.population['ecc'] = tools.power_law(self.Nbin, self.e_min, self.e_max, self.e_slope)
        pbar.update((5/6)*100)
        self.population['a'] = tools.p2a(self.population['p'], self.population['m1'], self.population['m2'])        
        pbar.finish()   
        # The end
    
    # To use pandas dataframe use to define the object
    def _save_mobse_input(self, name, met, tmax, backup=10):
        """
        It saves a population in a file (MOBSE format).
        Input: 
            fname = name of the output file
            met = metallicity (float)
            tmax = max time (float)
        Output:
            __________________________________
            | nd | m1 | m2 | p | ecc | Z | t |
        """

        if self.single_pop:
            print("\nSingle population print to file only supported for SEVN input style.")
            print("Please consider to turn single_pop argument to True.\n")
            return

        Z = np.full(self.Nbin, met)
        t = np.full(self.Nbin, tmax)
        ind = np.arange(1,self.Nbin+1,1)            
        m1, m2, p, ecc = np.hsplit(self.population[['m1','m2','p','ecc']].to_numpy(),4)
        np.savetxt(name+"_"+str(met)+".in", 
                np.c_[ind, m1, m2, p, ecc, Z, t], 
                fmt=('%i %4.4f %4.4f %10.4f %1.4f %1.4f %6.1f'),
                header=str(self.Nbin-backup), delimiter=' ', comments='')

    def _save_sevn_input(self, name, z1=None, z2=None, o1=0.0, o2=0.0, tend=None, tstart1=None, \
                         sn1=None, sn2=None, dtout=None, tstart2=None, dt=None, Sevn_v=2):
        """
        It saves a population in a file (SEVN format).
        Input: 
            z1, z2 = metallicity of the stars 1 and 2 (arrays)
            o1, o2 = spin of the stars 1 and 2 [1/yr] (arrays)
            tend = ending time of the simulation [Myr] (float)
            tstart1 = starting time of the simulation [Myr] (float)
            dt = time step of the integration [Myr] (not used)
            sn1, sn2 = supernova explosion mechanism of the stars (string: _delayed_ / _rapid_ / _startrack_)
            dtout = time step for printing outputs [Myr] (float)
        Output:
            ...
        """

        #Set placeholder
        placeholder="xxx"
        if z1 is None:
            z1 = placeholder
        if z2 is None:
            z2 = placeholder
        if tend is None:
            tend = placeholder
        if tstart1 is None:
            tstart1 = placeholder
        if sn1 is None:
            sn1 = placeholder
        if sn2 is None:
            sn2 = placeholder
        if dtout is None:
            dtout = placeholder

        # transform values to string since in SEVN they can be either strings or floats
        tostring = [z1,z2,o1,o2,tend,tstart1,sn1,sn2,dtout,tstart2]

        for ith, i in enumerate(tostring):
            if type(i) == float:
                tostring[ith] = str(round(i, 3))

        z1,z2,o1,o2,tend,tstart1,sn1,sn2,dtout,tstart2 = tostring

        M1 = self.population['m1'].values
        if not self.single_pop:
            M2 = self.population['m2'].values
            A = self.population['a'].values
            E = self.population['ecc'].values

            # to remove eventual 1 due to formatting round process
            E[E > 0.999] = 0.999

        
        with open(name+"_"+z1+".in", mode='w') as f:
            
            if self.single_pop:
                for i in range(self.Nbin):
                    line = f'{M1[i]:10.3f} {z1:>10}' \
                           f'{o1:>10} {sn1:>10}' \
                           f'{tstart1:>10}' \
                           f'{tend:>10} {dtout:>10}'

                    print(line, file=f)
            
            elif Sevn_v == 1:
                
                if not dt:
                    dt = 0.1
                
                for i in range(self.Nbin):
                    line = f'{M1[i]:10.3f} {M2[i]:10.3f}' \
                           f'{z1:>10} {z2:>10}' \
                           f'{o1:>10} {o2:>10}' \
                           f'{A[i]:10.3f} {E[i]:10.3f}' \
                           f'{tend:>10} {tstart1:>10}' \
                           f'{dt:10.3f} {sn1:>10}' \
                           f'{sn2:>10} {dtout:>10}'

                    print(line, file=f)
            
            elif Sevn_v == 2:
                
                if not tstart2:
                    tstart2 = tstart1

                for i in range(self.Nbin):
                    line = f'{M1[i]:10.3f} {z1:>10}' \
                           f'{o1:>10} {sn1:>10}' \
                           f'{tstart1:>10} {M2[i]:10.3f}' \
                           f'{z2:>10} {o2:>10}' \
                           f'{sn2:>10} {tstart2:>10}' \
                           f'{A[i]:12.3g} {E[i]:12.3g}' \
                           f'{tend:>10} {dtout:>10}'

                    print(line, file=f)

        # SN1 = np.array(SN1, dtype=str)
        # tmp_arr = np.array(np.vstack([M1, M2, A, E, SN1]).T, dtype=str)
        # print(tmp_arr)
        # np.savetxt(name+".in", 
        #     tmp_arr, #, SN2, DTOUT], 
        #     fmt='%4.4f %4.4f %4.4f %4.4f %10s',
        #     delimiter=' ', comments='')
        return

    @staticmethod
    def save_mobse_input(name, m1, m2, p, ecc, met, tmax, backup=10):
        """
        It saves a population in a file (MOBSE format).
        Input: 
            fname = name of the output file
            m1 = primary masses (array)
            m2 = secondary masses (array)
            p = periods (array)
            ecc = eccentricity (array)
            met = metallicity (float)
            tmax = max time (float)
        Output:
            __________________________________
            | nd | m1 | m2 | p | ecc | Z | t |
        """
        Nbin = len(m1)
        Z = np.full(Nbin,float(met))
        t = np.full(Nbin,tmax)
        ind = np.arange(1,Nbin+1,1)
        np.savetxt(name+"_"+str(met)+".in", 
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
