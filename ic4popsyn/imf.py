"""
Classes for IMF
"""
import numpy as np
from abc import ABC, abstractmethod

class IMF(ABC):
    """
    Base template class
    """

    def __init__(self,mass_range):
        self._pdf_norm=None

    def generate(self,number_of_stars):
        """
        Function to generate stars from the IMF
        """
        number_of_stars=int(number_of_stars)
        return self._generate(number_of_stars)

    def pdf(self,mass):
        """
        Return the IMF pdf at mass
        """
        return self._pdf(mass)

    def cdf(self,mass):
        """
        Return the IMF cdf for M<mass
        """
        return self._cdf(mass)


    @abstractmethod
    def _pdf(self,mass):
        """
        Specialised Function to estimate pdf
        """
        pass

    @abstractmethod
    def _cdf(self,mass):
        """
        Specialised Function to estimate cdf
        """
        pass

    @abstractmethod
    def _generate(self,number_of_stars):
        """
        Specialised Function to generate stars from the IMF
        """
        pass

class PowerLaw(IMF):
    """
    It samples stellar masses from a power-law M^alpha
    Input:
        mass_ranges= tuples containing the minimum and maximum mass
        alpha= power-law slope
    """
    def __init__(self,mass_range,  alpha):
        """
        It samples stellar masses from a power-law M^alpha
        Input:
            mass_ranges= tuples containing the minimum and maximum mass
            alpha= power-law slope
        """
        self._mass_range=mass_range
        self._alpha=alpha
        self._mmin=min(mass_range)
        self._mmax=max(mass_range)


    def _generate(self,number_of_stars):
        """
        Specialised Function to generate stars from the IMF
        """
        #Inverse sampling
        u=np.random.uniform(0,1,number_of_stars)
        if self._alpha==-1:
            mmin_log=np.log(self._mmin)
            masses = np.exp( (np.log(self._mmax) - mmin_log)*u + mmin_log )
        else:
            slope = 1+self._alpha
            mmin_slope = self._mmin**slope
            masses = ( (self._mmax**slope - mmin_slope)*u + mmin_slope )**(1/slope)

        return masses

    def _pdf(self,mass):
        """
        Specialised Function to estimate pdf
        """

        if self._alpha==-1:
            norm=1/(np.log(self._mmax)-np.log(self._mmin))
            return norm/mass
        else:
            slope = 1+self._alpha
            norm = slope/(self._mmax**slope - self._mmin**slope)
            return norm*mass**self._alpha

    def _cdf(self,mass):
        """
        Specialised Function to estimate cdf
        """
        if self._alpha==-1:
            norm = 1 / (np.log(self._mmax) - np.log(self._mmin))
            return norm*(np.log(mass)-np.log(self._mmin))
        else:
            slope=1+self._alpha
            norm = 1/(self._mmax**slope - self._mmin**slope)
            return norm*(mass**slope - self._mmin**slope)



class BrokenPowerLaw(IMF):
    """
    It samples stellar masses from a broken-power-law with N pieces.
    Input:
        mass_ranges= tuple with N+1 edges of the broken power law pieces
        alphas= tuple with N values of power low index in each of the N pieces
    """
    def __init__(self,mass_ranges, breakpoints, alphas):
        """
        It samples stellar masses from a broken-power-law with N pieces.
        Input:
            mass_ranges= tuples containing the minimum and maximum mass
            breakpoints= points at which the power law change
            alphas= tuple with N values of power low index in each of the N pieces
        """
        self.mass_ranges=mass_ranges
        self.alphas=alphas



    def _generate(self,number_of_stars):
        """
        Specialed Function to generate stars from the IMF
        """

        random = np.random.random

        breaks = np.array(self.mass_ranges)
        slopes = np.array(self.alphas)

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

    def _continuity_constants(self,):
        """
        Given a broken power law pdf of the kind
        K0 x^-a0  for x<=b_0
        K1 x^-a1  for b_0<=x<=b_1

        Kn x^-an for x>=b_n-1
        with breaks points b_0..b_n-1 and slopes a0...an
        this generator returns the constant Kn to satisfy the continuity at the breakpoints.
        The last K is always set to 1, therefore the pdf will not integrate to 1!
        To be sure that the pdf integrates to 1 estimate the integral I and rescale of all the keys
        by I, i.e. Kn=Kn/I.
        """
        current_K=1
        iter = 0
        while (iter<len(slopes)):
            yield current_K
            if iter<len(breakpoints): current_K = current_K*(breakpoints[iter])**(slopes[iter+1]-slopes[iter])
            iter+=1


class Salpeter(BrokenPowerLaw):
    """
    Salpeter IMF dN/dM propto M^-2.35
    """
    def __init__(self,mass_ranges=(0.08,150)):

        if len(mass_ranges)!=2:
            raise ValueError(f"mass ranges in Salpeter IMF must contains exactly two values (lowest and highest mass)")

        super().__init__(mass_ranges=mass_ranges, alphas=(-2.35,))

class Kroupa(BrokenPowerLaw):
    """

    """
    def __init__(self, mass_ranges=(0.08,150)):
        if len(mass_ranges) != 2:
            raise ValueError(f"mass ranges in Kroupa IMF must contains exactly two values (lowest and highest mass)")

        mmin=min(mass_ranges)
        mmax=max(mass_ranges)

        kroupa_mass_ranges = (0.01,0.08,0.5,150)
        kroupa_alphas=(-0.3,-1.3,-2.3)

        _mass_ranges=None
        _alphas=None

        # This part below is need to get the right mass_ranges and alphas to incldue the input mass ranges
        # Limit cases:
        # Maximum range below minimum Kroupa edges
        if mmax<=kroupa_mass_ranges[0]:
            _mass_ranges=(mmin,mmax)
            _alphas=(kroupa_alphas[0],)
        # Maximum range above maximum Kroupa edges
        elif mmin>=kroupa_mass_ranges[-1]:
            _mass_ranges = (mmin, mmax)
            _alphas = (kroupa_alphas[-1],)
        # All the other siutations
        else:
            _mass_ranges = None
            _alphas = None
            #Setting the starting point
            imin=0
            imax=len(kroupa_mass_ranges)-1
            # Not consider the last item since it will be anyway overriden if mmax is larger than the -2 item
            for i,km in enumerate(kroupa_mass_ranges[:-1]):
                if mmin>=km: imin=i
                if mmax<=km:
                    imax=i
                    break # Maximum found we can stop

            _mass_ranges = [m for m in kroupa_mass_ranges[imin:imax+1]]
            _mass_ranges[0]  = mmin
            _mass_ranges[-1] = mmax
            _alphas = [m for m in kroupa_alphas[imin:imax]]


        super().__init__(mass_ranges=tuple(_mass_ranges), alphas=tuple(_alphas))



