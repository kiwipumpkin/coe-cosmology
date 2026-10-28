# Derived: I_STOKES of SMICA 2048 downgraded to Nside 1024 with healpy.ud_grade (not an official PLA product)
import healpy as hp
m=hp.read_map('COM_CMB_IQU-smica_2048_R3.00_full.fits',field=0)
hp.write_map('smica_I_nside1024_udgrade_from2048.fits',hp.ud_grade(m,1024),overwrite=True,column_units='K_CMB')
