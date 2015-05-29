#
#
from model import *

#Loading isochrone
age=120	
Fe_H=0

isochrone=Isochrone(age,Fe_H)
isodata, Teff, log_g, log_L = isochrone.load()


#Compute magnitudes
magnitudes=[]
for i in np.arange(len(Teff)):
	if i==0:	
		magnitudes = isochrone.colorize(Teff[i],log_g[i],log_L[i])
	else:
		magnitudes = np.column_stack((magnitudes,
				isochrone.colorize(Teff[i],log_g[i],log_L[i])))


#Save magnitudes in a new file with previous quantities
isochrone.save_unspotted(isodata,magnitudes)


#Compute the effect of spots on stars

zeta=1	#luminosity ratio
epsilon=1	#surface ratio
rho=0.001	#spot coverage
pi=	1	#Tspot/Tphot

spots_params = [zeta,epsilon,rho,pi]
isodata_spots, nlog_g, Tphot, log_Lphot, Tspot, log_Lspot = \
		isochrone.add_spots(isodata,spots_params)

#Compute magnitudes of the photosphere of spotted stars
mmagnitudes_spots=[]
for i in np.arange(len(Tspot)):
	if i == 0:
		magnitudes_spots =isochrone.colorize(Tspot[i],nlog_g[i],
				log_Lspot[i]) + \
				isochrone.colorize(Tphot[i],nlog_g[i],log_Lphot[i])
	else:
		magnitudes_spots = np.column_stack((magnitudes,
				isochrone.colorize(Tspot[i],nlog_g[i],log_Lspot[i]) +
				isochrone.colorize(Tphot[i],nlog_g[i],log_Lphot[i])))

#Save spotted isochrone
isochrone.save_spotted(spots_params,isodata_spots,magnitudes_spots)
