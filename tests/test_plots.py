import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ic4popsyn import populations as pop
from ic4popsyn import tools
from scipy import integrate

def Sana_plots(Nbin):
    # create population of binaries
    binaries = pop.Binaries(Nbin, model='sana12',mass_ranges=[5.,150],alphas=[-2.3])
    binaries.population['logP'] = np.log10(binaries.population['p'])
    binaries.population['q'] = binaries.population['m2']/binaries.population['m1']

    # plot distributions and compare them with theoretical distributions
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[8,8])
    axs = axs.ravel()

    keys = ['m1','q','logP','ecc']
    colors = ['C0','C1','C2','C3']
    alphas = [-2.3,-0.1,-0.55,-0.45]
    bins = [np.logspace(np.log10(5.),np.log10(150),30),30,30,30]
    scales = [['log','log'],['linear','linear'],['linear','linear'],['linear','linear']]
    texts = [r'$f(m) \propto m^{-2.3}$ (Kroupa01)', r'$f(q) \propto q^{-0.1}$ (Sana+12)', 
            r'$f(logP) \propto logP^{-0.55}$ (Sana+12)', r'$f(ecc) \propto ecc^{-0.45}$ (Sana+12)']
    txs = [0.06,0.3,0.2,0.2]
    tys = [0.75,0.95,0.9,0.9]
    # theoretical function
    powerlaw = lambda x, alpha, b : b*x**(alpha)

    for ax,key,bi,scale,a,c,t,tx,ty in zip(axs,keys,bins,scales,alphas,colors,texts,txs,tys):
        low = binaries.population[key].min()
        if key == 'ecc':
            low = 0.001

        # theoratical distributions
        up = binaries.population[key].max()
        x = np.linspace(low,up)
        k = integrate.quad(powerlaw, low, up, args=(a,1))
        y = powerlaw(x,a,1/k[0]) # normalization
        ax.plot(x,y,c='k', label='Theo. distr.')
        if key == 'm1':
            ax.text(tx*up, ty*np.max(y), t)
        else:
            ax.text(tx*up, ty*np.max(y), t)
        
        # samples
        ax.hist(binaries.population[key],histtype='step',density=True,bins=bi,color=c, label='simulation')
        ax.set_xscale(scale[0])
        ax.set_yscale(scale[1])
        ax.set_ylabel('PDF')
        ax.set_xlabel(key)
        ax.legend(loc='center right')
    plt.show(block=True)

def SanaMDS_plots(Nbin):
    # create population of binaries
    binaries = pop.Binaries(Nbin, model='sana_eccM&DS',mass_range=[5.,150],alphas=[-2.3])

    # plot distributions and compare them with theoretical distributions
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[10,5])

    # theoretical function
    x = np.logspace(np.log10(2),5.5,1000)
    y = tools.eccvsP(x)
    ax.plot(x,y,c='r', lw=3, label='Theo. distr.')
    # samples
    ax.scatter(binaries.population['p'], binaries.population['ecc'], s=3, c='k', alpha=0.15, label='simulation')
    ax.set_xscale('log')
    ax.set_ylabel('ecc')
    ax.set_xlabel('$P$ [days]')
    ax.legend(loc='upper left')
    plt.show(block=True)

def test_plots():
    Sana_plots(100000)
    SanaMDS_plots(100000)
    assert True