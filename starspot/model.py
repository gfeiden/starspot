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
        """ Initialize an instance of Isochrone  """
        self.age   = age
        self.Fe_H  = Fe_H
        self.a_Fe  = a_Fe
        self.mix   = mix
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

        self.mass = isodata[0]
        self.logT = isodata[1]
        self.logg = isodata[2]
        self.logL = isodata[3]
        self.logR = isodata[4]
        self.A_Li = isodata[5]
        
        return isodata, self.logT, self.logg, self.logL

    def unload(self):
        """ """

    def add_spots(self, isodata, spots_params):
        """ adds spots to isochrone model assuming a two-temperature model """
        self.Teff = 10**self.logT

        self.zeta    = spots_params[0] # fractional luminosity brightening/dimming
        self.epsilon = spots_params[1] # fractional radius de/inflation
        self.rho     = spots_params[2] # spot surface coverage
        self.pi      = spots_params[3] # spot temperature constrast

        # Properties
        nlog_L = np.log10(self.zeta) + self.logL         # spotted star log(L/Lsun)
        nlog_R = 0.5*np.log10(self.epsilon) + self.logR  # spotted star log(R/Rsun)
        nR     = 10**nlog_R                         # new radius [Rsun]
        nlog_g = np.log10(G*self.mass*M_SUN/(R_SUN*nR)**2)  # updated log(g)
        Tavg   = (L_SUN*10**nlog_L/(4*np.pi*SIGMA*(R_SUN*nR)**2))**0.25  # new Teff

        # Photosphere
        Tphot     = self.Teff*(self.zeta/(self.epsilon*(1 - self.rho*(1 - self.pi**4))))**0.25 # background phot temp
        Sphot     = 4*np.pi*(1 - self.rho)*(R_SUN*nR)**2                   # surface w/o spots
        log_Lphot = np.log10((SIGMA*Sphot*Tphot**4)/L_SUN)            # luminosity from unspotted regions

        # Spots
        Tspot = self.pi*Tphot  # spot temperature (umbral) [K]
        if self.rho != 0 and self.pi != 0:
            Sspot     = 4*np.pi*self.rho*(R_SUN*nR)**2         # fractional surface coverage of spots
            log_Lspot = np.log10((SIGMA*Sspot*Tspot**4)/L_SUN) # luminosity from spotted regions
        else:
            log_Lspot = np.zeros(len(self.mass))
            log_Lspot.fill(np.nan)

        # Make an array of all computed data.
        isodata_spots = np.vstack((self.mass, Tavg, nlog_g, nlog_L, nlog_R, self.A_Li,
                Tphot, Tspot, log_Lphot, log_Lspot))

        return isodata_spots, nlog_g, Tphot, log_Lphot, Tspot, log_Lspot

    def colorize(self, Teff, log_g, log_L, filters=None):
        """ transforms stellar parameters into UBVRIJHK magnitudes """
        if filters == None:
            filters = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
        else:
            pass

        self.filters = filters

        return bc.bolcorrection.bc_eval(Teff, log_g, log_L, len(filters))
        

    def save_unspotted(self, isodata, magnitudes):
        """ Write out a copy of the unspotted isochrone with synthetic photometry """

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


        # data
        for i in np.arange(len(self.mass)):
            iso_io.write("{:8.2f}{:10.6f}{:10.6f}{:10.6f}{:10.6f}{:10.6f}".format(
                           self.mass[i], self.logT[i], self.logg[i], self.logL[i], 
                           self.logR[i], self.A_Li[i]))
            for mag in magnitudes:
                iso_io.write("{:10.6f}".format(mag[i]))
            iso_io.write("\n")
            
        iso_io.close()

    def save_spotted(self, spots_params, isodata_spots, magnitudes_spots):
        """ Write out properties of a spotted isochrone with synthetic photometry """
        #zeta    = spots_params[0]
        #epsilon = spots_params[1]
        #rho     = spots_params[2]
        #pi      = spots_params[3]

        # Create the directory where data will be saved if needed
        color_iso_directory = "./iso/age_{:07.1f}myr_z{:+05.2f}".format(self.age, self.Fe_H)
        #if exists(color_iso_directory) == False:
        #    os.mkdir(color_iso_directory)
        #else:
        #    pass
        color_iso_directory = "./iso"

        # Create a file with data comparable to observations. Included :
        # m, Tavg, log_g (new), log_L (new), log_R (new), A(Li)
        # U B V Rc Ic J H K (2MASS)
        logTavg = np.log10(isodata_spots[1])
        nlogg = isodata_spots[2]
        nlogL = isodata_spots[3]
        nlogR = isodata_spots[4]
        
        # format isochrone file name
        magfile = "dmestar_{:07.1f}myr_z{:+5.2f}_a{:+5.2f}_gs98_phx".format(self.age,
                   self.Fe_H, self.a_Fe)
        if self.magnetic:
            if type(self.magnetic) == 'float':
                magfile += "_mag{:3.1f}".format(self.mag_field)
            elif type(self.magnetic) == 'str':
                magfile += "_mag{:s}".format(self.mag_field)
            else:
                raise TypeError("Invalid Magnetic Field Strength Data Type")
            magfile += "_zet{:+5.2f}_eps{:+5.2f}_rho{:+5.2f}_pi{:+5.2f}".format(
                        self.zeta, self.epsilon, self.rho, self.pi)
        else:
            magfile += "_zet{:+5.2f}_eps{:+5.2f}_rho{:+5.2f}_pi{:+5.2f}".format(
                        self.zeta, self.epsilon, self.rho, self.pi)
        
        iso_io = open(color_iso_directory + "/" + magfile + ".phot", "w")

        # format header 
        iso_io.write('# Dartmouth Stellar Evolution Isochrone with Synthetic Photometry & Starspots\n' + 
                '#'+'\n'+'#'+'    '+
                'Age = '+str("%7.1f" % self.age)+' Myr    '+
                '[Fe/H] = '+str("%+5.2f" % self.Fe_H)+' dex    '+
                '[a/Fe] = '+str("%+5.2f" % self.a_Fe)+' dex    '+
                '\n'+'#'+'\n')
        iso_io.write("#    zeta = {:+5.2f}, epsilon = {:+5.2f},".format(self.zeta, self.epsilon))
        iso_io.write(" rho = {:+5.2f}, pi = {:+5.2f} \n".format(self.rho, self.pi))
        iso_io.write("# \n")
        
        # column labels
        iso_io.write("# {:8s}{:10s}{:10s}{:10s}{:10s}{:10s}".format(
                     "Mass", "log(Tavg)", "log(g)", "log(L)", "log(R)", "A(Li)"))
        for f in self.filters:
            iso_io.write("{:10s}".format(f))
        iso_io.write("\n")


        # units
        iso_io.write("# {:8s}{:10s}{:10s}{:10s}{:10s}{:10s}".format(
                     "[Msun]", "[K]", "[cm s^-2]", "[Lsun]", "[Rsun]", "[dex]"))
        for f in self.filters:
            iso_io.write("{:10s}".format("[mag]"))
        iso_io.write("\n")


        # data
        for i in np.arange(len(self.mass)):
            iso_io.write("{:8.2f}{:10.6f}{:10.6f}{:10.6f}{:10.6f}{:10.6f}".format(
                           self.mass[i], logTavg[i], nlogg[i], nlogL[i], nlogR[i], 
                           self.A_Li[i]))
            for mag in magnitudes_spots:
                iso_io.write("{:10.6f}".format(mag[i]))
            iso_io.write("\n")

        iso_io.close()

        # Create a file with data related to spot modelisation. Included :
        # m, Tphot, Tspot, log_Lphot, log_Lspot
        log_Tphot = np.log10(isodata_spots[6])
        log_Tspot = np.log10(isodata_spots[7])
        log_Lphot = isodata_spots[8]
        log_Lspot = isodata_spots[9]
        
        # Header
        spotfile = open(color_iso_directory + "/" + magfile + ".spot", "w")
        
        # column labels
        spotfile.write("# {:8s}{:10s}{:10s}{:10s}{:10s}{:10s}".format(
                     "Mass", "log(Tavg)", "log(g)", "log(L)", "log(R)", "A(Li)"))
        spotfile.write("{:12s}{:12s}{:12s}{:12s}".format("log(Tphot)", 
                     "log(Tspot)", "log(Lphot)", "log(Lspot)"))
        spotfile.write("\n")

        # units
        spotfile.write("# {:8s}{:10s}{:10s}{:10s}{:10s}{:10s}".format(
                     "[Msun]", "[K]", "[cm s^-2]", "[Lsun]", "[Rsun]", "[dex]"))
        for i in range(4):
            spotfile.write("{:12s}".format("[K]"))
        spotfile.write("\n")

        # data
        for i in np.arange(len(self.mass)):
            spotfile.write("{:8.2f}{:10.6f}{:10.6f}{:10.6f}{:10.6f}{:10.6f}".format(
                           self.mass[i], logTavg[i], nlogg[i], nlogL[i], nlogR[i], 
                           self.A_Li[i]))
            spotfile.write("{:12.6f}{:12.6f}{:12.6f}{:12.6f}".format(log_Tphot[i], 
                           log_Tspot[i], log_Lphot[i], log_Lspot[i]))
            
            spotfile.write("\n")
            
        spotfile.close()
        

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

    # Find total flux by summing contributions
    flux_tot = flux1 + flux2

    # Corresponding magnitude
    return -2.5*np.log10(flux_tot)
