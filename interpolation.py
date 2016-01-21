import numpy as np


def interpolate(self, field, particles=[]):

    if particles == []:
        particles = self.particles
    
    if 'triinterp' in self.opt.interpolation:
        vel = _triinterp(self, field, particles)    
    
    return vel



def _triinterp(self, field, particles):
    
    return np.zeros((len(particles.lon1),))
    
    
    
def scattered():
    pass
