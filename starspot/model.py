#
#
import numpy as np
import t2color

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
        A_Li = isodata[5]
        return isodata, m, Teff, log_g, log_L, log_R, A_Li

    def unload(self):
        """ """

    def colorize(self,*kind):
        """ calculates magnitudes from isochrone data """
        #Fortran subroutine
        magnitude = t2color()
        #Creates a file with data comparable to observations. Included :
        # m, Tavg, log_g (new), log_L (new), log_R (new), A(Li)
        # B V Rc Ic J H K (2MASS)
        magfile = open('mag_zet+{}_eps+{}_rho+{}_pi+{}.dat'.format(kind[1],
                       kind[2], kind[3], kind[4]), 'a')
        magfile.write('#'+'\n'+'#'+'    '+
                'kind = '+kind[0]+'    '+
                'zeta = '+str(kind[1])+'    '+
                'epsilon = '+str(kind[2])+'    '+
                'rho = '+str(kind[3])+'    '+
                'pi = '+str(kind[4])+
                '\n'+'#'+'\n'+'#'+'   '+
                'Mass'+'    '+'T_average'+'    '+
                'log(g)'+'    '+'log_L'+'    '+
                'log(R)'+'    '+'A(Li)'+'    '+
                'B'+'    '+'V'+'    '+
                'Rc'+'    '+'Ic'+'    '+
                'J'+'    '+'H'+'    '+
                'K'+'    ')
        magfile.close()
        #Creates a file with data related to spot modelisation. Included :
        # m, Tphot, Tspot, log_Lphot, log_Lspot
        modfile = open('spots_zet+{}_eps+{}_rho+{}_pi+{}.dat'.format(kind[1],
                       kind[2], kind[3], kind[4]), 'a')
        modfile.write('#'+'\n'+'#'+'    '+
                'kind = '+kind[0]+'    '+
                'zeta = '+str(kind[1])+'    '+
                'epsilon = '+str(kind[2])+'    '+
                'rho = '+str(kind[3])+'    '+
                'pi = '+str(kind[4])+
                '\n'+'#'+'\n'+'#'+'    '+
                'Mass'+'    '+'T_phot'+'    '+
                'T_spot'+'    '+'log(L_phot)'+'    '+
                'log(L_spot)'+'   ')
        modfile.close()
        return magnitude
        

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