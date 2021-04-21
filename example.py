from ic4popsyn import populations as pop

# Number of systems (10 systems are for backup)
Nbin = 1000010  
# create a population of binaries 
binSana = pop.Binaries(Nbin, model='sana12', mass_range=[5.,150.], alphas=[-2.3]) # available: sana12 / sana_eccm&ds
binSanaMDS = pop.Binaries(Nbin, model='sana_eccM&DS', mass_range=[5.,150.], alphas=[-2.3]) # available: sana12 / sana_eccm&ds
# save the population as input for MOBSE 
binSana.save_mobse_input('mobseInputS_Z', 0.02, 13600, 10) # with 10 system used for backup
binSanaMDS.save_mobse_input('mobseInputMDS_Z', 0.02, 13600, 10) # with 10 system used for backup

binSanaMDS.save_sevn_input('SEVNInputMDS_Z', 0.02, 0.02, 0.0, 0.0, 'end', 0.6, 'delayed', 'delayed', 'end')
