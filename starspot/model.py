#
#
import numpy as np
from astropy import constants as const

from os import mkdir
from os.path import exists

import color.bolcor as bc

#Constants
G=const.G.cgs.value
L_sun=const.L_sun.cgs.value
M_sun=const.M_sun.cgs.value
R_sun=const.R_sun.cgs.value
sigma=const.sigma_sb.cgs.value

class Isochrone(object):
	""" provides instance of Dartmouth stellar evolution isochrone """
    
	def __init__(self, age, Fe_H, a_Fe = 0.0, mix = 'gas07', a_mlt = 'solar'):
		""" """
		self.age = age
		self.Fe_H = Fe_H
		self.a_Fe = a_Fe
		self.mix = mix
		self.a_mlt = a_mlt

	def load(self):
		""" loads data from isochrone models """
		isofile = 'isochrones/dmestar_00{}myr_z+{}_a+{}_marcs.iso'.format(
				"%.1f" % self.age, "%.2f" % self.Fe_H, "%.2f" % self.a_Fe)
		isodata = np.genfromtxt(isofile,unpack=True)
		m = isodata[0]
		Teff = 10**isodata[1]
		log_g = isodata[2]
		log_L = isodata[3]
		log_R = isodata[4]
		a_Li = isodata[5]
		return isodata, Teff, log_g, log_L

	def unload(self):
		""" """

	def add_spots(self,isodata,spots_params):
		""" adds spots to isochrone models assuming a two-temperature model """
		m = isodata[0]
		Teff = 10**isodata[1]
		log_g = isodata[2]
		log_L = isodata[3]
		log_R = isodata[4]
		a_Li = isodata[5]

		zeta=spots_params[0]
		epsilon=spots_params[1]
		rho=spots_params[2]
		pi=spots_params[3]
		
		nlog_L=np.log10(zeta)+log_L
		nlog_R=0.5*np.log10(epsilon)+log_R
		nR=10**nlog_R
		nlog_g=np.log10(G*m*M_sun/(R_sun*nR)**2)

		Tphot=Teff*(zeta/(epsilon*(1-rho*(1-pi**4))))**0.25
		Sphot=4*np.pi*(1-rho)*(R_sun*nR)**2
		log_Lphot=np.log10((sigma*Sphot*Tphot**4)/L_sun)

		Tspot=pi*Tphot
		if rho==0 :
			log_Lspot=np.zeros(len(m))
			log_Lspot.fill(np.nan)
		else :
			Sspot=4*np.pi*rho*(R_sun*nR)**2
			log_Lspot=np.log10((sigma*Sspot*Tspot**4)/L_sun)

		Tavg=(L_sun*10**nlog_L/(4*np.pi*sigma*(R_sun*nR)**2))**0.25

		isodata_spots=np.vstack((m, Tavg, nlog_g, nlog_L, nlog_R, a_Li, Tphot,
				Tspot, log_Lphot, log_Lspot))
		return isodata_spots, nlog_g, Tphot, log_Lphot, Tspot, log_Lspot

	def colorize(self,Teff,log_g,log_L):
		""" calculates magnitudes from isochrone data """
	   	# initialize log file
		bc.utils.log_init('example.log')
		
		# initialize bolometric correction table at fixed [Fe/H] and [a/Fe]
		brand='marcs'
		filters=['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
		bc.bolcorrection.bc_init(self.Fe_H, self.a_Fe, brand, filters)
		
		# transform stellar parameters into UBVRIJHK magnitudes
		magnitudes = bc.bolcorrection.bc_eval(Teff, log_g, log_L, 
				len(filters))
			
		# release allocated memory and close log file
		bc.bolcorrection.bc_clean()
		bc.utils.log_close()

		return magnitudes

	def save_unspotted(self,isodata,magnitudes):
		""" """
		m = isodata[0]
		log_Teff = isodata[1]
		log_g = isodata[2]
		log_L = isodata[3]
		log_R = isodata[4]
		a_Li = isodata[5]
		U=magnitudes[0]
		B=magnitudes[1]
		V=magnitudes[2]
		R=magnitudes[3]
		I=magnitudes[4]
		J=magnitudes[5]
		H=magnitudes[6]
		K=magnitudes[7]

		if exists('age_{}+z_{}'.format("%.1f" % self.age, 
				"%.2f" % self.Fe_H))==False:
			mkdir('age_{}+z_{}'.format("%.1f" % self.age, 
				"%.2f" % self.Fe_H))

		isofile=open('age_{}+z_{}/isochrone_{}myr_z+{}_a+{}_marcs.iso'.format(
				"%.1f" % self.age, "%.2f" % self.Fe_H,
				"%.1f" % self.age, "%.2f" % self.Fe_H, 
				"%.2f" % self.a_Fe),'w')
		isofile.write('#'+'\n'+'#'+'    '+
				'age = '+str("%.1f" % self.age)+'    '+
				'[Fe/H] = '+str("%.2f" % self.Fe_H)+'    '+
				'a(Fe) = '+str("%.2f" % self.a_Fe)+'    '+
				'\n'+'#'+'\n'+'#'+'     '+
				'Mass'+'          '+
				'log(Teff)'+'     '+
				'log(g)'+'       '+
				'log(L)'+'        '+
				'log(R)'+'        '+
				'A(Li)'+'        '+
				'U'+'             '+
				'B'+'             '+
				'V'+'             '+
				'Rc'+'            '+
				'Ic'+'             '+
				'J'+'             '+
				'H'+'             '+
				'K'+'\n')
		for i in np.arange(len(m)):
			isofile.write('    '+str("%10.6f" % m[i])+'    '+
					str("%10.6f" % log_Teff[i])+'    '+
					str("%10.6f" % log_g[i])+'    '+
					str("%10.6f" % log_L[i])+'    '+
					str("%10.6f" % log_R[i])+'    '+
					str("%10.6f" % a_Li[i])+'   '+
					str("%10.6f" % U[i])+'    '+
					str("%10.6f" % B[i])+'    '+
					str("%10.6f" % V[i])+'    '+
					str("%10.6f" % R[i])+'    '+
					str("%10.6f" % I[i])+'    '+
					str("%10.6f" % J[i])+'    '+
					str("%10.6f" % H[i])+'    '+
					str("%10.6f" % K[i])+'\n')
		isofile.close()

	def save_spotted(self,spots_params,isodata_spots,magnitudes_spots):
		""" """
		zeta=spots_params[0]
		epsilon=spots_params[1]
		rho=spots_params[2]
		pi=spots_params[3]
		#Creates a file with data comparable to observations. Included :
		# m, Tavg, log_g (new), log_L (new), log_R (new), A(Li)
		# U B V Rc Ic J H K (2MASS)
		m=isodata_spots[0]
		log_Tavg=np.log10(isodata_spots[1])
		nlog_g=isodata_spots[2]
		nlog_L=isodata_spots[3]
		nlog_R=isodata_spots[4]
		a_Li=isodata_spots[5]
		U=magnitudes_spots[0]
		B=magnitudes_spots[1]
		V=magnitudes_spots[2]
		R=magnitudes_spots[3]
		I=magnitudes_spots[4]
		J=magnitudes_spots[5]
		H=magnitudes_spots[6]
		K=magnitudes_spots[7]

		if exists('age_{}+z_{}'.format("%.1f" % self.age, 
				"%.2f" % self.Fe_H))==False:
			mkdir('age_{}+z_{}'.format("%.1f" % self.age, 
				"%.2f" % self.Fe_H))

		magfile = open('age_{}+z_{}/mag_zet+{}_eps+{}_rho+{}_pi+{}.dat'\
				.format("%.1f" % self.age, "%.2f" % self.Fe_H,zeta, epsilon, 
				rho, pi), 'w')
		magfile.write('#'+'\n'+'#'+'    '+
				'zeta = '+str(zeta)+'    '+
				'epsilon = '+str(epsilon)+'    '+
				'rho = '+str(rho)+'    '+
				'pi = '+str(pi)+
				'\n'+'#'+'\n'+'#'+'     '+
				'Mass'+'          '+
				'log(Tavg)'+'     '+
				'log(g)'+'       '+
				'log(L)'+'        '+
				'log(R)'+'        '+
				'A(Li)'+'        '+
				'U'+'             '+
				'B'+'             '+
				'V'+'             '+
				'Rc'+'            '+
				'Ic'+'             '+
				'J'+'             '+
				'H'+'             '+
				'K'+'\n')
		for i in np.arange(len(m)):
			magfile.write('    '+str("%10.6f" % m[i])+'    '+
					str("%10.6f" % log_Tavg[i])+'    '+
					str("%10.6f" % nlog_g[i])+'    '+
					str("%10.6f" % nlog_L[i])+'    '+
					str("%10.6f" % nlog_R[i])+'    '+
					str("%10.6f" % a_Li[i])+'   '+
					str("%10.6f" % U[i])+'    '+
					str("%10.6f" % B[i])+'    '+
					str("%10.6f" % V[i])+'    '+
					str("%10.6f" % R[i])+'    '+
					str("%10.6f" % I[i])+'    '+
					str("%10.6f" % J[i])+'    '+
					str("%10.6f" % H[i])+'    '+
					str("%10.6f" % K[i])+'\n')
		magfile.close()

		#Creates a file with data related to spot modelisation. Included :
		# m, Tphot, Tspot, log_Lphot, log_Lspot
		m=isodata_spots[0]
		log_Tphot=np.log10(isodata_spots[6])
		log_Tspot=np.log10(isodata_spots[7])
		log_Lphot=isodata_spots[8]
		log_Lspot=isodata_spots[9]
		modfile = open('age_{}+z_{}/spots_zet+{}_eps+{}_rho+{}_pi+{}.dat'\
				.format("%.1f" % self.age, "%.2f" % self.Fe_H, zeta, epsilon, 
				rho, pi), 'w')
		modfile.write('#'+'\n'+'#'+'    '+
				'zeta = '+str(zeta)+'    '+
				'epsilon = '+str(epsilon)+'    '+
				'rho = '+str(rho)+'    '+
				'pi = '+str(pi)+
				'\n'+'#'+'\n'+'#'+'     '+
				'Mass'+'          '+
				'log(T_phot)'+'   '+
				'log(T_spot)'+'  '+
				'log(L_phot)'+'   '+
				'log(L_spot)'+'\n')
		for i in np.arange(len(m)):
			modfile.write('    '+str("%10.6f" % m[i])+'    '+
					str("%10.6f" % log_Tphot[i])+'    '+
					str("%10.6f" % log_Tspot[i])+'    '+
					str("%10.6f" % log_Lphot[i])+'    '+
					str("%10.6f" % log_Lspot[i])+'\n')
		modfile.close()
        

class MassTrack(object):
	""" provides instance of Dartmouth stellar evolution mass track """

	def __init__(self, mass, Fe_H, a_Fe = 0.0, mix = 'gas07', a_mlt = 'solar'):
		""" """
		self.mass = mass
		self.Fe_H = Fe_H
		self.a_Fe = a_Fe
		self.mix = mix
		self.a_mlt = a_mlt

	def load(self):
		""" """

	def unload(self):
		""" """

