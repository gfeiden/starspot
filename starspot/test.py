#
#
from model import *

# Load isochrone.
age = 120.0
Fe_H = 0.0
a_Fe = 0.0
isochrone=Isochrone(age, Fe_H)
isodata, Teff, log_g, log_L = isochrone.load()


# Initialize log file.
bc.utils.log_init('example.log')

# Initialize bolometric correction table at fixed [Fe/H] and [a/Fe].
brand = 'marcs'
filters = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
bc.bolcorrection.bc_init(Fe_H, a_Fe, brand, filters)

# Compute magnitudes.
magnitudes = []
for i in np.arange(len(Teff)):
	if i == 0:
		magnitudes = isochrone.colorize(Teff[i], log_g[i], log_L[i])
	else:
		magnitudes = np.column_stack((magnitudes,
				isochrone.colorize(Teff[i], log_g[i], log_L[i])))

# Save magnitudes in a new file with previously loaded quantities.
isochrone.save_unspotted(isodata, magnitudes)


# Free parameters.
zeta = np.arange(1,2)	# Luminosity ratio
epsilon = np.arange(1,2)	# Surface ratio
rho = np.arange(0,1)	# Spot coverage
pi = np.arange(1,2) # Tspot/Tphot

for m in np.arange(len(zeta)):
	for n in np.arange(len(epsilon)):
		for p in np.arange(len(rho)):
			for q in np.arange(len(pi)):
				# Compute the effect of spots on stars.
				spots_params = [zeta[m], epsilon[n], rho[p], pi[q]]
				isodata_spots, nlog_g, Tphot, log_Lphot, Tspot, log_Lspot = \
						isochrone.add_spots(isodata, spots_params)
				# Compute magnitudes of spotted stars.
				magnitudes_spots = []
				for i in np.arange(len(Tphot)):
					if i != 0 and rho[p] != 0:
						magnitudes_spots = np.column_stack((magnitudes_spots,
								isochrone.colorize(Tspot[i], nlog_g[i],
								log_Lspot[i]) +	isochrone.colorize(Tphot[i],
								nlog_g[i], log_Lphot[i])))
					elif i == 0 and rho[p] != 0:
						magnitudes_spots = isochrone.colorize(Tspot[i],
								nlog_g[i], log_Lspot[i]) + \
								isochrone.colorize(Tphot[i], nlog_g[i],
								log_Lphot[i])
					elif i == 0:
						magnitudes_spots = magnitudes
					else:
						magnitudes_spots = magnitudes
				# Save spotted isochrone
				isochrone.save_spotted(spots_params, isodata_spots,
						magnitudes_spots)

# Release allocated memory and close log file
bc.bolcorrection.bc_clean()
bc.utils.log_close()
