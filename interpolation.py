from __future__ import division, print_function
import numpy as np
import sys


def intriangle(xt, yt, x0, y0):
    """
    ** Function checks if points are in a specified triangle**
    
    Inputs:
      - xt, yt - n X 3 array defining the x,y locations of each triangles node
      - x0, y0 - The n points to check     
           
    """
    
    f1 = (y0 - yt[:, 0]) * (xt[:, 1] - xt[:, 0]) -\
         (x0 - xt[:, 0]) * (yt[:, 1] - yt[:, 0])
    f2 = (y0 - yt[:, 2]) * (xt[:, 0] - xt[:, 2]) -\
         (x0 - xt[:, 2]) * (yt[:, 0] - yt[:, 2])
    f3 = (y0 - yt[:, 1]) * (xt[:, 2] - xt[:, 1]) -\
         (x0 - xt[:, 1]) * (yt[:, 2] - yt[:, 1])    
    
    return ((f1 * f3 >= 0) & (f3 * f2 >= 0))
    
    
def interpolate(self, field, particles=[]):
    """ 
    ** Interpolates particle velocities for the given field. 
       Using the interpolation method specified in options.**
    
    Inputs:
      - self - pyticleClass
      - field - the field to use for interpolation
      - particles - the particles to get the velocities for
    """

    if particles == []:
        particles = self.particles
    
    if 'triinterp' in self.opt.interpolation:
        vel = _triinterpE(self, field, particles)    
    
    return vel


def _triinterpE(self, field, particles):
    """ 
    ** FVCOM interpolation method for element data using grid parameters.**
    
    Inputs:
      - self - pyticleClass
      - field - the field to use for interpolation
      - particles - the particles to get the velocities for
    """
    
    
    grid = self.grid
       
    # Find particles indomain
    indom = particles.indomain!=-1
    # Get node numbers for each particles element
    ns = grid.nv[particles.indomain[indom], :]
    # Check if each particle is in each element if not find the particles element
    intri = intriangle(grid.x[ns], grid.y[ns], particles.xpt[indom], particles.ypt[indom])
    if intri.all():
        hosts = particles.indomain
    else:
        # Figure out which particles aren't in there element and find there element
        check = indom*False
        check[indom]=~intri
        particles.indomain[check] = grid.finder.__call__(particles.xpt[check], particles.ypt[check])
        hosts = particles.indomain

    
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
