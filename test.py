import Classes
import numpy 

# test object
sys = Classes.Population(110, model='SaNa12', alphas=[-1.3,-2.3])
print(sys.alphas, sys.model, sys.Nstars)
sys.save_mobse_input('sanaCLASS','0.02',14000)

from Classes import Population as pop

# test single functions
primary = pop.IMF(20)
q = pop.mass_ratio(20)
secondary = primary * q
secondary = numpy.where(secondary < 0.1, 0.1, secondary)
period = pop.period(20)
ecc = pop.eccentricity(20)
pop.save_mobse_input('sanaOWN',primary,secondary,period,ecc,'0.02',14000)
