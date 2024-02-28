import numpy
from ic4popsyn import populations as pop

# Number of systems (1 systems are for backup)
Nbin = 10000000
fimf = []
fimf2 = []
# create a population of binaries 
for z in range(10):

    print(z)
    binSana = pop.Binaries(Nbin, model='sana12') # available: sana12 / sana_eccm&ds
    binSana.generate()
    # save the population as input for MOBSE 
    nHigh = binSana.population[binSana.population['m1'] >= 5]
    highMass = nHigh.m1.sum() #+ nHigh.m2.sum()
    totMass = binSana.population.m1.sum() #+ binSana.population.m2.sum()
    fimf += [highMass / totMass * 100.]
    highMass = nHigh.m1.sum() + nHigh.m2.sum()
    totMass = binSana.population.m1.sum() + binSana.population.m2.sum()
    fimf2 += [highMass / totMass * 100.]

print('fimf = {}, {}'.format(numpy.mean(fimf),numpy.mean(fimf2)))
