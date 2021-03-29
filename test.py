import numpy as np
import pandas as pd
from ic import Population as pop
from scipy.optimize import curve_fit

def func(x, a, b):
    return b*x**a

# create population as object
Nbin = 100000
sys1 = pop(Nbin, model='SaNa12',mass_range=[5.,150],alphas=[-2.3])
# save data
sys1.save_mobse_input('mobse_object','0.02',14000)
# compute log(P) and q
sys1.binaries['logP'] = np.log10(sys1.binaries['p'])
sys1.binaries['q'] = sys1.binaries['m2']/sys1.binaries['m1']

# create population from single functions
primary = pop.IMF(Nbin)
q = pop.power_law(Nbin, 0.1, 1, -0.1)
secondary = primary * q
secondary = np.where(secondary < 0.1, 0.1, secondary)
period = pow(10.,pop.power_law(Nbin, 0.15, 5.5, -0.55))
ecc = pop.power_law(Nbin, 0., 1, -0.45)
# save data
pop.save_mobse_input('mobse_functions',primary,secondary,period,ecc,'0.02',14000)
sys2 = pd.DataFrame({'m1': primary, 'q':q, 'logP':np.log10(period), 'ecc':ecc})

# plot distributions and compare them with theoretical distributions
import matplotlib.pyplot as plt
from scipy import integrate

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[8,8])
axs = axs.ravel()

keys = ['m1','q','logP','ecc']
colors = ['C0','C1','C2','C3']
consts = [8., 1, 0.22, 1.]
alphas = [-2.3,-0.1,-0.55,-0.45]
bins = [np.logspace(np.log10(5.),np.log10(150),30),30,30,30]
scales = [['log','log'],['linear','linear'],['linear','linear'],['linear','linear']]
texts = [r'$f(m) \propto m^{-2.3}$ (Kroupa01)', r'$f(q) \propto q^{-0.1}$ (Sana+12)', 
         r'$f(logP) \propto logP^{-0.55}$ (Sana+12)', r'$f(ecc) \propto ecc^{-0.45}$ (Sana+12)']

# theoretical function
powerlaw = lambda x, alpha : x**(alpha)

for ax,key,bin,scale,a,c,t in zip(axs,keys,bins,scales,alphas,colors,texts):
    low = sys1.binaries[key].min()
    if key == 'ecc':
        low = 0.001

    # theoratical distributions
    up = sys1.binaries[key].max()
    x = np.linspace(low,up)
    k = integrate.quad(powerlaw, low, up, args=(a))
    y = powerlaw(x,a)
    ax.plot(x,1/k[0]*y,c='k')
    if key == 'm1':
       ax.text(0.1*up, 0.8*np.max(1/k[0]*y), t)
    else:
       ax.text(0.2*up, 0.95*np.max(1/k[0]*y), t)

    # compute fits of the distributions
    hist, edges = np.histogram(sys1.binaries[key],density=True,bins=bin)
    centers = 0.5*(edges[1:]+ edges[:-1])
    pars, cov = curve_fit(func, centers, hist)
    ax.plot(centers, func(centers, *pars))

    # samples
    ax.hist(sys1.binaries[key],histtype='step',density=True,bins=bin,color=c, label='object')
    ax.hist(sys2[key],histtype='stepfilled',alpha=0.2,density=True,bins=bin,color=c, label='functions')
    ax.set_xscale(scale[0])
    ax.set_yscale(scale[1])
    ax.set_ylabel('PDF')
    ax.set_xlabel(key)
    ax.legend(loc='center right')

#todo: add check fitting the distributions!
# assert alpha == fitted exponent
plt.show()