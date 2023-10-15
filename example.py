from ic4popsyn import populations as pop
from ic4popsyn import tools

# Number of systems (10 systems are for backup)
Nbin = 100010  
# create a population of binaries 
binSana = pop.Binaries(Nbin, model='sana12', mass_ranges=[5.,150.], alphas=[-2.3]) # available: sana12 / sana_eccm&ds
binSanaMDS = pop.Binaries(Nbin, model='sana_eccM&DS', mass_ranges=[5.,150.], alphas=[-2.3]) # available: sana12 / sana_eccm&ds
SinglePop = pop.Binaries(Nbin, model='sana_eccM&DS', single_pop=True, mass_ranges=[5.,150.], alphas=[-2.3]) # available: sana12 / sana_eccm&ds

# save the population as input for MOBSE 
binSana.save_mobse_input('mobseInputS_Z', 0.02, 13600, 10) # with 10 system used for backup
binSanaMDS.save_mobse_input('mobseInputMDS_Z', 0.02, 13600, 10) # with 10 system used for backup

# save the population as input for SEVN
binSanaMDS.save_sevn_input('SEVNInputMDS_Z', 0.02, 0.02, 0.0, 0.0, 'end', 0.6, 'delayed', 'delayed', 'end')
binSanaMDS.save_sevn_input('SEVNInputSingleMDS_Z', tend='end', dtout='end')

import os 
# download MOBSE
codeName = 'mobse'
getMOBSE = tools.GetFromUrl(codeName)
getMOBSE.main()
print(f'MOBSE is ready: {os.path.exists(codeName)}')
