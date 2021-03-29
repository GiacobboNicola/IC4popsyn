import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.stats import chisquare
from scipy.optimize import curve_fit
from classes.ic import Population as pop

def Sana_population(N):
    # create population of binaries
    sys1 = pop(N, model='SaNa12',mass_range=[5.,150],alphas=[-2.3])
    sys1.binaries['logP'] = np.log10(sys1.binaries['p'])
    sys1.binaries['q'] = sys1.binaries['m2']/sys1.binaries['m1']

    # plot distributions and compare them with theoretical distributions
    _, axs = plt.subplots(nrows=2, ncols=2, figsize=[8,8])
    axs = axs.ravel()

    keys = ['m1','q','logP','ecc']
    colors = ['C0','C1','C2','C3']
    alphas = [-2.3,-0.1,-0.55,-0.45]
    bins = [np.logspace(np.log10(5.),np.log10(150),100),100,100,100]
    scales = [['log','log'],['linear','linear'],['linear','linear'],['linear','linear']]
    texts = [r'$f(m) \propto m^{-2.3}$ (Kroupa01)', r'$f(q) \propto q^{-0.1}$ (Sana+12)', 
            r'$f(logP) \propto logP^{-0.55}$ (Sana+12)', r'$f(ecc) \propto ecc^{-0.45}$ (Sana+12)']

    chis = []

    # theoretical function
    powerlaw = lambda x, alpha, b : b*x**(alpha)

    for ax,key,bin,scale,a,c,t in zip(axs,keys,bins,scales,alphas,colors,texts):

        # compute fits of the distributions
        hist, edges = np.histogram(sys1.binaries[key], density=True, bins=bin)
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