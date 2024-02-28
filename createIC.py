import numpy
from ic4popsyn import populations as pop

# Number of systems (1 systems are for backup)
Nbin = 101  
backup = 1
# create a population of binaries 
for z in [0.02, 0.008, 0.002, 0.0008, 0.0002]:

    binSana = pop.Binaries(Nbin, model='sana12', mass_range=[5.,150.], alphas=[-2.3]) # available: sana12 / sana_eccm&ds
    binSana.generate()
    # save the population as input for MOBSE
    binSana.save_mobse_input('mobse', z, 13600, backup)

    type1 = [1]*Nbin
    type2 = [1]*Nbin
    tini = [0.0]*Nbin
    m1 = binSana.population.m1
    m2 = binSana.population.m2
    p = binSana.population.p
    e = binSana.population.ecc

    numpy.savetxt("petar_"+str(z)+".in", 
        numpy.c_[m1, m2, type1, type2, p/(365.24*10.**6), e, tini], 
        fmt="%4.4f %4.4f %i %i %15.9f %1.4f %1.2f",
        header=str(Nbin-backup),comments='')