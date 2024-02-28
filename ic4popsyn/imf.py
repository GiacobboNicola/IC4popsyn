"""
Classes for IMF
"""
import numpy as np
from scipy import integrate
from scipy import special
from scipy import stats

class IMF():
    """
    Base template class
    Input:
        mass_range= tuples containing the minimum and maximum mass
    """

    def __init__(self,mass_range):

        if len(mass_range)!=2:
            raise ValueError(f"mass ranges in Salpeter IMF must contains exactly two values (lowest and highest mass)")

        self._mass_range=mass_range
        self._mmin=min(mass_range)
        self._mmax=max(mass_range)
        if self._mmin<=0: raise ValueError("minimum mass cannot be lower or equal 0")

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

        if isinstance(mass,float) or isinstance(mass,int):
            if mass<self._mmin or mass>self._mmax: return 0
            else: return self._pdf(mass)
        else:
            pdf_array=np.where( (mass>=self._mmin) & (mass<=self._mmax), self._pdf(mass), 0)
            return pdf_array


    def cdf(self,mass):
        """
        Return the IMF cdf for M<mass
        """

        if isinstance(mass,float) or isinstance(mass,int):
            if mass<self._mmin: return 0
            elif  mass>self._mmax: return 1
            else: return self._cdf(mass)
        else:
            cdf_array=np.zeros_like(mass)
            cdf_array[mass>self._mmax]=1
            idx=(mass>=self._mmin) & (mass<=self._mmax)
            cdf_array[idx] = self._cdf(mass[idx])
            return cdf_array


    def mcdf(self,mass,mtot=1):
        """
        Return the integral of m times the IMF up to mass in input.
        Notice it is normalised so that the total mass is equal to mtot
        Input:
            mass= number or array of numbers with the mass where to calculate the enclosed mass
            mtot= normalised total mass
        """
        mcdf_max = self._mcdf(self._mmax)
        norm = mtot/mcdf_max

        return norm*self._mcdf(mass)

    def mpdf(self,mass):
        """
        Pdf of the integral m*IMF
        """
        norm = self.mcdf(self._mmax)
        return norm*mass*self.pdf(mass)


    def _pdf(self,mass):
        """
        Specialised Function to estimate pdf
        """
        raise NotImplementedError("Specialised method _pdf not implemented")

    def _cdf(self,mass):
        """
        Specialised Function to estimate cdf.
        If the specialised function is not provided, use a numerical integration
        """
        int_function = lambda m:  self.pdf(m)

        if isinstance(mass,int) or isinstance(mass,float):
            cdf = integrate.quad(int_function,self._mmin,mass)[0]
        else:
            cdf = np.zeros_like(mass)
            for i,m in enumerate(mass):
                cdf[i] = integrate.quad(int_function,self._mmin,m)[0]

        return cdf


    def _mcdf(self,mass):
        """
        Specialised Function to estimate m times the IMF up to mass in input
        If the specialised function is not provided, use a numerical integration
        """
        int_function = lambda  m: m*self.pdf(m)

        if isinstance(mass, int) or isinstance(mass, float):
            mcdf = integrate.quad(int_function, self._mmin, mass)[0]
        else:
            mcdf = np.zeros_like(mass)
            for i, m in enumerate(mass):
                mcdf[i] = integrate.quad(int_function, self._mmin, m)[0]

        return mcdf


    def _generate(self,number_of_stars):
        """
        Specialised Function to generate stars from the IMF
        """
        raise NotImplementedError("Specialised method _generate not implemented")

class PowerLaw(IMF):
    """
    It samples stellar masses from a power-law M^alpha
    Input:
        mass_range= tuples containing the minimum and maximum mass
        alpha= power-law slope
    """
    def __init__(self,mass_range,  alpha):
        """
        It samples stellar masses from a power-law M^alpha
        Input:
            mass_ranges= tuples containing the minimum and maximum mass
            alpha= power-law slope
        """
        self._alpha=alpha
        super(PowerLaw, self).__init__(mass_range=mass_range) #Initialise mass_range


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

    def _xcdf(self,mass):
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
        mass_range= tuples containing the minimum and maximum mass
        break_points= Points
        alphas= tuple with N values of power low index in each of the N pieces
    """

    # To be built


class Salpeter(PowerLaw):
    """
    Salpeter IMF dN/dM propto M^-2.35
    """
    def __init__(self, mass_range=(0.08,150)):

        super(Salpeter, self).__init__(mass_range=mass_range, alpha=-2.35)

class LogGau(IMF):
    """
    An IMF in which
    dN/dM prop to a log normal in m
    """

    def __init__(self,mass_range,m_mean,logm_std):
        """
        Input:
            mass_range= tuples containing the minimum and maximum mass
            m_mean = mass of the peak of the lognormal
            logm_std= std if tge log normal
        """

        # Since we are depends on scipy
        # Here we use the stat

        # We assume that the pdf is A*exp-(lm - lm_mean)^2/(2*std^2)
        # where A is the normalisation constant (to be found)
        # this is the pdf in the logspace log10(m), in the linear space
        # we have to add the Jacobian d log10(m) / dm  = 1/(m*ln(10))

        self._m_mean   = m_mean
        if self._m_mean<=0: raise ValueError("m_mean mass must be strictly positive")
        self._logm_std = logm_std
        if self._m_mean<=0: raise ValueError("logm_std  must be strictly positive")
        super(LogGau, self).__init__(mass_range)

        # We set not the normalisation constant  A
        # considering our pdf as A * _pdf where
        # _pdf is exp-(lm - lm_mean)^2/(2*std^2)
        # and A is found integrating such pdf in the mass range
        self._pdf_norm = 1/self._cdf_not_norm(self._mmin,self._mmax)

    def _pdf(self,mass):
        """
        This is simply a truncated Gaussian in the log10 space
        """
        logmean = np.log10(self._m_mean)
        logmass = np.log10(mass)
        std = self._logm_std
        t = (logmass-logmean)/std
        Jacobian = 1/(mass*np.log(10))
        # The Jacobian takes into account that the Gaussian is defined in the log10(m) space
        # but the pdf is in the linear m space (d log10(m) = d m /(m*ln(10))

        return self._pdf_norm*Jacobian*np.exp(-0.5*t*t)

    def _cdf(self,mass):
        """
        This is simply a truncated Gaussian in the log10 space
        see _cdf_not_norm
        """

        return self._pdf_norm*self._cdf_not_norm(self._mmin,mass)

    def _cdf_not_norm(self,mmin,mmax):
        """
        This the method to find the total integral of the pdf given the mass range assuming norm A=1
        in  A*exp-(lm - lm_mean)^2/(2*std^2). If we introduce the variable t=(lm - lm_mean)/(np.sqrt(2) std)
        we have pdf(t) = exp(-t^2) with Jacobian d lm / dt = sqrt(2) std,
        Considering that the error function is defined as Errf(t) = 2/sqrt(pi) * int^z_0 e^-t^2 dt
        the solution of our integral between tmin=(log(xmin) -lm_mean)/(np.sqrt(2) std) and xmax=(log(xmax) -lm_mean)/(np.sqrt(2) std) is
        sqrt(pi/2) * std (Errf(tmax)-Errf(tmin)), Errf is taken fro scipy
        """

        ft = lambda m: (np.log10(m) - np.log10(self._m_mean))/(self._logm_std*np.sqrt(2))

        tmin = ft(mmin)
        tmax = ft(mmax)

        return np.sqrt(np.pi/2)*self._logm_std*(special.erf(tmax) - special.erf(tmin))

    def _generate(self,number_of_stars):
        """
        This is a truncated Gaussian, we can in principle use rejection sampling
        or create an interpolation and then invert it.
        But since we are already depending on scipy, we can use the truncated normal distribution
        in scipy stats so that we can exploit vectorisation
        @TODO: for consistenscy we should use the same scipy utilities also for the pdf and cdf
        """

        # In stats the limit need to be given with y=(lm - lm_mean)/(std) (similar to what we use in th cdf)

        logmean = np.log10(self._m_mean)
        ft = lambda m: (np.log10(m) - logmean)/(self._logm_std)

        tmin = ft(self._mmin)
        tmax = ft(self._mmax)

        rv = stats.truncnorm(tmin, tmax, loc=logmean, scale=self._logm_std) #@TODO when using scipy everywhere we can initialise this obj in the init
        logm=rv.rvs(number_of_stars)
        return 10**logm


class LogGauPowerLaw(IMF):
    """
    An IMF in which
    dN/dM prop to a log m Gaussian m<mbreak
    dN/dM prop to  m**alpha for m>mbreak
    """




class _BrokenPowerLaw(IMF):
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


class _Kroupa(BrokenPowerLaw):
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




