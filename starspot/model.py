#
#
import numpy as np

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
        isofile = 'dmestar_00{}myr_z+{}0_a+{}0_marcs.iso'.format(self.age,
                self.Fe_H,self.a_Fe)
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