#
#
import numpy as np
import bolcor as bc

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
		isofile = 'dmestar_00{}myr_z+{}_a+{}_marcs.iso'.format(
				"%.1f" % self.age, "%.2f" % self.Fe_H, "%.2f" % self.a_Fe)
		isodata = np.genfromtxt(isofile,unpack=True)
		m = isodata[0]
		Teff = isodata[1]
		log_g = isodata[2]
		log_L = isodata[3]
		log_R = isodata[4]
		a_Li = isodata[5]
		return isodata, m, Teff, log_g, log_L, log_R, a_Li

	def unload(self):
		""" """

	def spot(self):
		""" """
		return isodata_spots

	def colorize(self,Teff,log_g,log_L):
		""" calculates magnitudes from isochrone data """
	   	# initialize log file
		bc.utils.log_init('example.log')
		
		# initialize bolometric correction table at fixed [Fe/H] and [a/Fe]
		brand='marcs'
		filters=['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
		bc.bolcorrection.bc_init(self.Fe_H, self.a_Fe, brand, filters)
		
		# transform stellar parameters into UBVRIJHK magnitudes
		magnitudes = bc.bolcorrection.bc_eval(Teff, log_g, log_L, len(filters))
			
		# release allocated memory and close log file
		bc.bolcorrection.bc_clean()
		bc.utils.log_close()

		return magnitudes

	def save_unspotted(self,isodata,magnitudes):
		""" """
		m = isodata[0]
		Teff = isodata[1]
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
		isofile=open('isochrone_{}myr_z+{}_a+{}_marcs.iso'.format(
				"%.1f" % self.age, "%.2f" % self.Fe_H, "%.2f" % self.a_Fe),'w')
		isofile.write('#'+'\n'+'#'+'    '+
				'age = '+str("%.1f" % self.age)+'    '+
				'[Fe/H] = '+str("%.2f" % self.Fe_H)+'    '+
				'a(Fe) = '+str("%.2f" % self.a_Fe)+'    '+
				'\n'+'#'+'    '+
				'Mass'+'    '+'Teff'+'    '+
				'log(g)'+'    '+'log_L'+'    '+
				'log(R)'+'    '+'A(Li)'+'    '+
				'U'+'    '+'B'+'    '+'V'+'    '+
				'Rc'+'    '+'Ic'+'    '+
				'J'+'    '+'H'+'    '+
				'K'+'    '+'\n')
		for i in np.arange(len(m)):
			isofile.write('    '+str(m[i])+'    '+str(Teff[i])+'    '+
					str(log_g[i])+'    '+str(log_L[i])+'    '+
					str(log_R[i])+'    '+str(a_Li[i])+'   '+
					str(U[i])+'    '+str(B[i])+'    '+str(V[i])+'    '+
					str(R[i])+'    '+str(I[i])+'    '+str(J[i])+'    '+
					str(H[i])+'    '+str(K[i])+'\n')
		magfile.close()

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
		Tavg=isodata_spots[1]
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
		magfile = open('mag_zet+{}_eps+{}_rho+{}_pi+{}.dat'.format(zeta,
                       epsilon, rho, pi), 'w')
		magfile.write('#'+'\n'+'#'+'    '+
				'zeta = '+str(zeta)+'    '+
				'epsilon = '+str(epsilon)+'    '+
				'rho = '+str(rho)+'    '+
				'pi = '+str(pi)+
				'\n'+'#'+'\n'+'#'+'   '+
				'Mass'+'    '+'T_average'+'    '+
				'log(g)'+'    '+'log(L)'+'    '+
				'log(R)'+'    '+'A(Li)'+'    '+
				'U'+'    '+'B'+'    '+'V'+'    '+
				'Rc'+'    '+'Ic'+'    '+
				'J'+'    '+'H'+'    '+
				'K'+'    '+'\n')
		for i in np.arange(len(m)):
			magfile.write('    '+str(m[i])+'    '+str(Tavg[i])+'    '+
					str(nlog_g[i])+'    '+str(nlog_L[i])+'    '+
					str(nlog_R[i])+'    '+str(a_Li[i])+'   '+
					str(U[i])+'    '+str(B[i])+'    '+str(V[i])+'    '+
					str(R[i])+'    '+str(I[i])+'    '+str(J[i])+'    '+
					str(H[i])+'    '+str(K[i])+'\n')
		magfile.close()

		#Creates a file with data related to spot modelisation. Included :
		# m, Tphot, Tspot, log_Lphot, log_Lspot
		m=isodata_spots[0]
		Tphot=isodata_spots[6]
		Tspot=isodata_spots[7]
		log_Lphot=isodata_spots[8]
		log_Lspot=isodata_spots[9]
		modfile = open('spots_zet+{}_eps+{}_rho+{}_pi+{}.dat'.format(zeta,
                       epsilon, rho, pi), 'w')
		modfile.write('#'+'\n'+'#'+'    '+
				'zeta = '+str(zeta)+'    '+
				'epsilon = '+str(epsilon)+'    '+
				'rho = '+str(rho)+'    '+
				'pi = '+str(pi)+
				'\n'+'#'+'\n'+'#'+'    '+
				'Mass'+'    '+'T_phot'+'    '+
				'T_spot'+'    '+'log(L_phot)'+'    '+
				'log(L_spot)'+'   '+'\n')
		for i in np.arange(len(m)):
			modfile.write('    '+str(m[i])+'    '+str(Tphot[i])+'    '+
					str(Tspot[i])+'    '+str(log_Lphot[i])+'    '+
					str(log_Lspot[i])+'\n')
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

