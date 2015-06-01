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

#values of free parameters

zeta=np.arange(1,2)	#luminosity ratio
epsilon=np.arange(1,2)	#surface ratio
rho=np.arange(0,1)	#spot coverage
pi=np.arange(1,2) #Tspot/Tphot

for m in np.arange(len(zeta)) :
	for n in np.arange(len(epsilon)) :
		for p in np.arange(len(rho)) :
			for q in np.arange(len(pi)) :
				#Compute the effect of spots on stars
				spots_params = [zeta[m],epsilon[n],rho[p],pi[q]]
				isodata_spots, nlog_g, Tphot, log_Lphot, Tspot, log_Lspot = \
						isochrone.add_spots(isodata,spots_params)

				#Compute magnitudes of spotted stars
				magnitudes_spots=[]
				for i in np.arange(len(Tphot)):
					if i == 0 and rho != 0 :
						magnitudes_spots = isochrone.colorize(Tspot[i],
								nlog_g[i],log_Lspot[i]) + \
								isochrone.colorize(Tphot[i],nlog_g[i],
								log_Lphot[i])
					elif i != 0 and rho != 0 :
						magnitudes_spots = np.column_stack((magnitudes_spots,
								isochrone.colorize(Tspot[i],nlog_g[i],
								log_Lspot[i]) +	isochrone.colorize(Tphot[i],
								nlog_g[i],log_Lphot[i])))
					elif i == 0 :
						magnitudes_spots = isochrone.colorize(Tphot[i],
								nlog_g[i],log_Lphot[i])
					else :
						magnitudes_spots = np.column_stack((magnitudes_spots,
								isochrone.colorize(Tphot[i],nlog_g[i],
								log_Lphot[i])))

				#Save spotted isochrone
				isochrone.save_spotted(spots_params,isodata_spots,
						magnitudes_spots)
