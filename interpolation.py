from __future__ import division, print_function
import numpy as np
import sys


def interpolate(self, field, particles=[]):

    if particles == []:
        particles = self.particles
    
    if 'triinterp' in self.opt.interpolation:
        vel = _triinterp(self, field, particles)    
    
    return vel

def _triinterp(self, field, particles):
    
    
    grid = self.grid
    
    # Finder element particle is in
    hosts = grid.finder.__call__(particles.xpt, particles.ypt)
    
    # Find layer
    if '2D' in self.opt.gridDim:
        layer = None
    else:
        # code to find layers above and below particle
        print('Still need to added 3D interpolation.')
        pass
    
    
    
    # Get distance from element center
    x0c = particles.xpt - grid.xc[hosts]
    y0c = particles.ypt - grid.yc[hosts] 
    
    # Get neighbours    
    e0=grid.nbe[hosts, 0]
    e1=grid.nbe[hosts, 1]
    e2=grid.nbe[hosts, 2]
    
    var_e = (field[layer, hosts]).flatten()   
    var_0 = (field[layer, e0]).flatten()
    var_1 = (field[layer, e1]).flatten()
    var_2 = (field[layer, e2]).flatten()
    var_0[e0==-1] = 0
    var_1[e1==-1] = 0
    var_2[e2==-1] = 0        
    
    dvardx = grid.a1u[0, hosts] * var_e + grid.a1u[1, hosts] * var_0 +\
             grid.a1u[2, hosts] * var_1 + grid.a1u[3, hosts] * var_2
    dvardy = grid.a2u[0, hosts] * var_e + grid.a2u[1, hosts] * var_0 +\
             grid.a2u[2, hosts] * var_1 + grid.a2u[3, hosts] * var_2
    
    var = var_e + dvardx * x0c + dvardy * y0c
    
    return var
    
    
    
def scattered():
    pass
