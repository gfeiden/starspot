#
#

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
        """ """

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
