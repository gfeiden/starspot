#
#
import numpy as np
from astropy import constants as const

import os
from os.path import exists

from color import bolcor as bc

# Constants
G = const.G.cgs.value
L_SUN = const.L_sun.cgs.value
M_SUN = const.M_sun.cgs.value
R_SUN = const.R_sun.cgs.value
SIGMA = const.sigma_sb.cgs.value

class Isochrone(object):
    """ provides instance of Dartmouth stellar evolution isochrone """
    
    def __init__(self, age, Fe_H, a_Fe=0.0, mix='gas07', a_mlt='solar', 
                 magnetic=False, mag_field=0.0e0):
        """ """
        self.age = age
        self.Fe_H = Fe_H
        self.a_Fe = a_Fe
        self.mix = mix
        self.a_mlt = a_mlt
        
        self.magnetic  = magnetic
        self.mag_field = mag_field

    def load(self):
        """ loads data from isochrone models """
        
        # establish the proper model directory accounting for B field and Z
        base_directory = os.environ["DARTMOUTH_GRID"]
        if self.magnetic:
            base_directory += "/mag/iso"
        else:
            base_directory += "/std/iso"
        base_directory += "/{:+04.0f}".format(self.Fe_H*100.).replace('-', 'm').replace('+', 'p')
        
        #
        # place holder for [alpha/Fe]
        #
        
        # format isochrone file name
        isofile = "dmestar_{:07.1f}myr_z{:+5.2f}_a{:+5.2f}_gs98_phx".format(self.age, 
                   self.Fe_H, self.a_Fe)
        if self.magnetic:
            if type(self.magnetic) == 'float':
                isofile += "_mag{:3.1f}.iso".format(self.mag_field)
            elif type(self.magnetic) == 'str':
                isofile += "_mag{:s}.iso".format(self.mag_field)
            else:
                raise TypeError("Invalid Magnetic Field Strength Data Type")
        else:
            isofile += ".iso"
            
        isodata = np.genfromtxt(base_directory + '/' + isofile, unpack=True)
        m     = isodata[0]
        Teff  = 10**isodata[1]
        log_g = isodata[2]
        log_L = isodata[3]
        log_R = isodata[4]
        A_Li  = isodata[5]
        
        return isodata, Teff, log_g, log_L

    def unload(self):
        """ """

    def add_spots(self, isodata, spots_params):
        """ adds spots to isochrone model assuming a two-temperature model """
        m = isodata[0]
        Teff = 10**isodata[1]
        log_g = isodata[2]
        log_L = isodata[3]
        log_R = isodata[4]
        a_Li = isodata[5]

        zeta = spots_params[0]
        epsilon = spots_params[1]
        rho = spots_params[2]
        pi = spots_params[3]

        # Properties
        nlog_L = np.log10(zeta)+log_L
        nlog_R = 0.5*np.log10(epsilon)+log_R
        nR = 10**nlog_R
        nlog_g = np.log10(G*m*M_SUN/(R_SUN*nR)**2)
        Tavg = (L_SUN*10**nlog_L/(4*np.pi*SIGMA*(R_SUN*nR)**2))**0.25
        # Photosphere
        Tphot = Teff*(zeta/(epsilon*(1-rho*(1-pi**4))))**0.25
        Sphot = 4*np.pi*(1-rho)*(R_SUN*nR)**2
        log_Lphot = np.log10((SIGMA*Sphot*Tphot**4)/L_SUN)
        # Spots
        Tspot = pi*Tphot
        if rho != 0 and pi != 0:
            Sspot = 4*np.pi*rho*(R_SUN*nR)**2
            log_Lspot = np.log10((SIGMA*Sspot*Tspot**4)/L_SUN)
        else:
            log_Lspot = np.zeros(len(m))
            log_Lspot.fill(np.nan)
        # Make an array of all computed data.
        isodata_spots = np.vstack((m, Tavg, nlog_g, nlog_L, nlog_R, a_Li,
                Tphot, Tspot, log_Lphot, log_Lspot))
        return isodata_spots, nlog_g, Tphot, log_Lphot, Tspot, log_Lspot

    def colorize(self, Teff, log_g, log_L, filters=None):
        """ transforms stellar parameters into UBVRIJHK magnitudes """
        if filters == None:
            filters = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
        else:
            pass

        self.filters = filters

        magnitudes = bc.bolcorrection.bc_eval(Teff, log_g, log_L, 
                len(filters))

        return magnitudes

    def save_unspotted(self, isodata, magnitudes):
        """ """
        m = isodata[0]
        log_Teff = isodata[1]
        log_g = isodata[2]
        log_L = isodata[3]
        log_R = isodata[4]
        a_Li = isodata[5]

        # Create the directory where data will be saved if needed
        color_iso_directory = "./iso/age_{:07.1f}myr_z{:+05.2f}".format(self.age, self.Fe_H)
        #if exists(color_iso_directory) == False:
        #    os.mkdir(color_iso_directory)
        #else:
        #    pass
        color_iso_directory = "./iso"

        # format isochrone file name
        isofile = "dmestar_{:07.1f}myr_z{:+5.2f}_a{:+5.2f}_gs98_phx".format(self.age,
                   self.Fe_H, self.a_Fe)
        if self.magnetic:
            if type(self.magnetic) == 'float':
                isofile += "_mag{:3.1f}.phot".format(self.mag_field)
            elif type(self.magnetic) == 'str':
                isofile += "_mag{:s}.phot".format(self.mag_field)
            else:
                raise TypeError("Invalid Magnetic Field Strength Data Type")
        else:
            isofile += ".phot"

	iso_io = open(color_iso_directory + "/" + isofile, "w")

        # format header 
        iso_io.write('# Dartmouth Stellar Evolution Isochrone with Synthetic Photometry\n' + 
                '#'+'\n'+'#'+'    '+
                'Age = '+str("%7.1f" % self.age)+' Myr    '+
                '[Fe/H] = '+str("%+5.2f" % self.Fe_H)+' dex    '+
                '[a/Fe] = '+str("%+5.2f" % self.a_Fe)+' dex    '+
                '\n'+'#'+'\n')

        # column labels
        iso_io.write("# {:8s}{:10s}{:10s}{:10s}{:10s}{:10s}".format(
                     "Mass", "log(Teff)", "log(g)", "log(L)", "log(R)", "A(Li)"))
        for f in self.filters:
            iso_io.write("{:10s}".format(f))
        iso_io.write("\n")

        # units
        iso_io.write("# {:8s}{:10s}{:10s}{:10s}{:10s}{:10s}".format(
                     "[Msun]", "[K]", "[cm s^-2]", "[Lsun]", "[Rsun]", "[dex]"))
        for f in self.filters:
            iso_io.write("{:10s}".format("[mag]"))
        iso_io.write("\n")


        # Data
        for i in np.arange(len(m)):
            iso_io.write("{:8.2f}{:10.6f}{:10.6f}{:10.6f}{:10.6f}{:10.6f}".format(
                    m[i], log_Teff[i], log_g[i], log_L[i], log_R[i], a_Li[i]))
            for mag in magnitudes:
                iso_io.write("{:10.6f}".format(mag[i]))
            iso_io.write("\n")
            
        iso_io.close()

    def save_spotted(self, spots_params, isodata_spots, magnitudes_spots):
        """ """
        zeta = spots_params[0]
        epsilon = spots_params[1]
        rho = spots_params[2]
        pi = spots_params[3]

        # Create the directory where data will be saved if needed.
        if exists('age_{}+z_{}'.format("%.1f" % self.age, 
                "%.2f" % self.Fe_H)) == False:
            os.mkdir('age_{}+z_{}'.format("%.1f" % self.age, 
                "%.2f" % self.Fe_H))

        # Create a file with data comparable to observations. Included :
        # m, Tavg, log_g (new), log_L (new), log_R (new), A(Li)
        # U B V Rc Ic J H K (2MASS)
        m = isodata_spots[0]
        log_Tavg = np.log10(isodata_spots[1])
        nlog_g = isodata_spots[2]
        nlog_L = isodata_spots[3]
        nlog_R = isodata_spots[4]
        a_Li = isodata_spots[5]
        U = magnitudes_spots[0]
        B = magnitudes_spots[1]
        V = magnitudes_spots[2]
        R = magnitudes_spots[3]
        I = magnitudes_spots[4]
        J = magnitudes_spots[5]
        H = magnitudes_spots[6]
        K = magnitudes_spots[7]
        # Header
        magfile = open('age_{}+z_{}/mag_zet+{}_eps+{}_rho+{}_pi+{}.dat'\
                .format("%.1f" % self.age, "%.2f" % self.Fe_H,
                "%.2f" % zeta, "%.2f" % epsilon, "%.2f" % rho,
                "%.2f" % pi), 'w')
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
        # Data
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

        # Create a file with data related to spot modelisation. Included :
        # m, Tphot, Tspot, log_Lphot, log_Lspot
        m = isodata_spots[0]
        log_Tphot = np.log10(isodata_spots[6])
        log_Tspot = np.log10(isodata_spots[7])
        log_Lphot = isodata_spots[8]
        log_Lspot = isodata_spots[9]
        # Header
        modfile = open('age_{}+z_{}/spots_zet+{}_eps+{}_rho+{}_pi+{}.dat'\
                .format("%.1f" % self.age, "%.2f" % self.Fe_H, 
                "%.2f" % zeta, "%.2f" % epsilon, "%.2f" % rho,
                "%.2f" % pi), 'w')
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
        # Data
        for i in np.arange(len(m)):
            modfile.write('    '+str("%10.6f" % m[i])+'    '+
                    str("%10.6f" % log_Tphot[i])+'    '+
                    str("%10.6f" % log_Tspot[i])+'    '+
                    str("%10.6f" % log_Lphot[i])+'    '+
                    str("%10.6f" % log_Lspot[i])+'\n')
        modfile.close()
        

class MassTrack(object):
    """ provides instance of Dartmouth stellar evolution mass track """

    def __init__(self, mass, Fe_H, a_Fe=0.0, mix='gas07', a_mlt='solar'):
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

def mag_tot(mag1, mag2):
    """ computes total magnitude from two contributions """
    # Total flux
    flux1 = 10**(-mag1/2.5)
    flux2 = 10**(-mag2/2.5)
    flux_tot = flux1+flux2
    # Corresponding magnitude
    m_tot = -2.5*np.log10(flux_tot)
    return m_tot

