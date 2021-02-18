import Classes

# test object
systems = Classes.Population(10000000, 'sana', alphas=[-1.0,-2.7])
print(systems.alphas, systems.Nstars)
print(systems.binaries)

from Classes import Population as pop
import numpy as np

# test single functions
primary = pop.IMF(1000)
q = pop.mass_ratio(1000)
secondary = primary * q
secondary = np.where(secondary < 0.1, 0.1, secondary)
period = pop.period(1000)
ecc = pop.eccentricity(1000)
print(primary, q, secondary, period, ecc)

