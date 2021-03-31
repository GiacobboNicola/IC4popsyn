"""
Test the sample of the Sana+12 distributions with a Kroupa IMF:
f(m) ~ m^{-2.3} 
f(q) ~ q^{-0.1}
f(logP) ~ logP^{-0.55}
f(ecc) ~ ecc^{-0.45}
"""
import numpy as np
from scipy import integrate
from scipy.stats import chisquare
from scipy.optimize import curve_fit
from ic4popsyn_pkg import ic
from ic4popsyn_pkg import populations as pop

def Sana_population(N):
    # create population of binaries
    binaries = pop.Binaries(N, model='SaNa12',mass_range=[5.,150],alphas=[-2.3])
    binaries.population['logP'] = np.log10(binaries.population['p'])
    binaries.population['q'] = binaries.population['m2']/binaries.population['m1']

    keys = ['m1','q','logP','ecc']
    alphas = [-2.3,-0.1,-0.55,-0.45]
    binss = [np.logspace(np.log10(5.),np.log10(150),100),100,100,100]

    chis = []

    # theoretical function
    powerlaw = lambda x, alpha, b : b*x**(alpha)

    for key,a,bins in zip(keys,alphas,binss):

        # compute fits of the distributions
        hist, edges = np.histogram(binaries.population[key], density=True, bins=bins)
        centers = 0.5*(edges[1:]+ edges[:-1])
        pars, cov = curve_fit(powerlaw, centers, hist)

        # theoratical distributions
        k = integrate.quad(powerlaw, centers[0], centers[-1], args=(a,1))
        y = powerlaw(centers,a,1/k[0])

        chi, _ = chisquare(hist, y)
        chis.append(chi)
    return chis

def test_answer():
    assert Sana_population(10000) < [0.1, 0.1, 0.1, 0.1]