#
#
from model import *

# Load isochrone.
age  = 4.0
Fe_H = 0.0
a_Fe = 0.0

iso = Isochrone(age, Fe_H)
iso.load()


# Initialize log file.
bc.utils.log_init('example.log')

# Initialize bolometric correction table at fixed [Fe/H] and [a/Fe].
brand = 'marcs'
filters = ['r', 'J', 'H', 'K']
bc.bolcorrection.bc_init(Fe_H, a_Fe, brand, filters)

# Compute magnitudes.
magnitudes = []
for i in np.arange(len(iso.Teff)):
	if i == 0:
		magnitudes = iso.colorize(iso.Teff[i], iso.logg[i], iso.logL[i], filters=filters)
	else:
		magnitudes = np.column_stack((magnitudes,
				iso.colorize(iso.Teff[i], iso.logg[i], iso.logL[i], filters=filters)))

# Save magnitudes in a new file with previously loaded quantities.
iso.save_unspotted(magnitudes)


# Free parameters.
zeta    = np.array([1]) # Luminosity ratio
epsilon = np.array([1])	# Surface ratio
rho     = np.array([0])	# Spot coverage
pi      = np.array([1]) # Tspot/Tphot

for m in np.arange(len(zeta)):
	for n in np.arange(len(epsilon)):
		for p in np.arange(len(rho)):
			for q in np.arange(len(pi)):
				# Compute the effect of spots on stars.
				spots_params = [zeta[m], epsilon[n], rho[p], pi[q]]
				isodata_spots, nlog_g, Tphot, log_Lphot, Tspot, log_Lspot = \
						iso.add_spots(spots_params)
				# Compute magnitudes of spotted stars.
				magnitudes_spots = []
				for i in np.arange(len(Tphot)):
					if i != 0 and rho[p] != 0 and pi[q] != 0:
						mag_spot = iso.colorize(Tspot[i], nlog_g[i],
								log_Lspot[i], filters=filters)
						mag_phot = iso.colorize(Tphot[i], nlog_g[i], 
								log_Lphot[i], filters=filters)
						mag_star = mag_tot(mag_spot, mag_phot)

						magnitudes_spots = np.column_stack((magnitudes_spots,
								mag_star))
					elif i == 0 and rho[p] != 0 and pi[q] != 0:
						mag_spot = iso.colorize(Tspot[i], nlog_g[i],
								log_Lspot[i], filters=filters)
						mag_phot = iso.colorize(Tphot[i], nlog_g[i], 
								log_Lphot[i], filters=filters)

						magnitudes_spots = mag_tot(mag_spot, mag_phot)
					elif i != 0:
						magnitudes_spots = np.column_stack((magnitudes_spots,
								iso.colorize(Tphot[i], nlog_g[i],
								log_Lphot[i], filters=filters)))
					else:
						magnitudes_spots = iso.colorize(Tphot[i],
							nlog_g[i], log_Lphot[i], filters=filters)
											
				# Save spotted isochrone
				iso.save_spotted(spots_params, isodata_spots,
						magnitudes_spots)

# Release allocated memory and close log file
bc.bolcorrection.bc_clean()
bc.utils.log_close()
