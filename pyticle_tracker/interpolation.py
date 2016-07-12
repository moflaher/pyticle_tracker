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
    
def find_layer(grid, particles, fieldshape):
    """
    ** Function to find which layer particles are in
       and calculate the layer weights**

    Inputs:
      - grid - grid object
      - particles - particle object
      - fieldshape - determine if data is siglev or siglay

    """
    # Find the upper layer
    layer = grid.siglaylen -1 - np.sum(particles.sigpt > grid.siglayrep, axis=0)
    # Find lower layer
    layer2 = layer + 1
    # If lower layer is below bottom use upper layer
    # Which layer is used won't matter because the weights will handle it
    layer2[layer2==grid.siglaylen] = -1

    # Initialize layers
    f1 = layer*0 + 0.0
    f2 = layer*0 + 0.0

    # Set weights between layers
    f1t = (particles.sigpt - grid.siglay[layer]) /\
          (grid.siglay[layer] - grid.siglay[layer2])
    f1[layer!=-1] = f1t[layer!=-1]
    f2[layer2!=-1] = 1 - f1t[layer2!=-1]

    # Set weights for above top layer (completely use lower layer
    # as it doesn't have an upper layer)
    f2[layer==-1] = 1
    # Set weights for below botton layer (same calculation as above
    # however the "bottom layer" siglay is -1)
    f1t = (particles.sigpt - grid.siglay[layer]) /\
          (grid.siglay[layer] - -1)
    f1[layer2==-1] = f1t[layer2==-1]
    
    print(grid.siglayrep.shape,grid.siglevrep.shape)
    print(particles.sigpt,f1,f2,layer,layer2)
    
    if fieldshape[0]==grid.siglevlen:
        #ttest siglev
        # Find the upper layer
        layer = grid.siglevlen -1 - np.sum(particles.sigpt > grid.siglevrep, axis=0)
        # Find lower layer
        layer2 = layer + 1
        # If lower layer is below bottom use upper layer
        # Which layer is used won't matter because the weights will handle it
        layer2[layer2==grid.siglevlen] = -1

        # Initialize layers
        f1 = layer*0 + 0.0
        f2 = layer*0 + 0.0

        # Set weights between layers
        f1t = (particles.sigpt - grid.siglev[layer]) /\
              (grid.siglev[layer] - grid.siglev[layer2])
        f1[layer!=-1] = f1t[layer!=-1]
        f2[layer2!=-1] = 1 - f1t[layer2!=-1]

        # Set weights for above top layer (completely use lower layer
        # as it doesn't have an upper layer)
        f2[layer==-1] = 1
        # Set weights for below botton layer (same calculation as above
        # however the "bottom layer" siglay is -1)
        f1t = (particles.sigpt - grid.siglay[layer]) /\
              (grid.siglay[layer] - -1)
        f1[layer2==-1] = f1t[layer2==-1]
    
        print(particles.sigpt,f1,f2,layer,layer2)
    
    return f1, f2, layer, layer2


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

    # Initialize varx and vary so there is always something to return.
    varx, vary = None, None

    if 'triinterp' in self.opt.interpolation:
        if self.grid.nele in field.shape:
            var = _triinterpE(self, field, particles)
        elif self.grid.node in field.shape:
            var, varx, vary = _triinterpN(self, field, particles)
        else:
            print("Given field doesn'''t match compatible dimensions.")
            sys.exit()

    return var, varx, vary

def _triinterpE(self, field, particles):
    """
    ** FVCOM interpolation method for element data using grid parameters.**

    Inputs:
      - self - pyticleClass
      - field - the field to use for interpolation
      - particles - the particles to get the velocities for
    """


    grid = self.grid

    hosts = __find_hosts(grid, particles)

    # Find layer
    if '2D' in self.opt.gridDim:
        layer = None
    else:
        # Find the upper layer
        f1, f2, layer, layer2 = find_layer(grid, particles, field.shape)

    # Get distance from element center
    x0c = particles.xpt - grid.xc[hosts]
    y0c = particles.ypt - grid.yc[hosts]

    # Get neighbouring elements
    e0 = grid.nbe[hosts, 0]
    e1 = grid.nbe[hosts, 1]
    e2 = grid.nbe[hosts, 2]

    def layer_vel(layer):
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

        vel = var_e + dvardx * x0c + dvardy * y0c

        return vel

    # Calculate velocities for 2D or upper layer in 3D
    vel = layer_vel(layer)

    if '3D' in self.opt.gridDim:
        # Interpolation lower layer
        vel2 = layer_vel(layer2)
        # Multiply by layer weights
        vel = vel * f1 + vel2 * f2

    # Zero velocity for any particles on land
    vel[hosts==-1] = 0
    # Zero velocity for any particles on the bottom
    vel[particles.inwater==0] = 0

    return vel

def _triinterpN(self, field, particles):
    """
    ** FVCOM interpolation method for node data using grid parameters.**

    Inputs:
      - self - pyticleClass
      - field - the field to use for interpolation
      - particles - the particles to get the depth/elevation for
    """

    grid = self.grid

    hosts = __find_hosts(grid, particles)

    # Find layer
    if len(field.shape)<2:
        layer = None
    else:
        # Find the upper layer
        f1, f2, layer, layer2 = find_layer(grid, particles, field.shape)
    

    # Get distance from element center
    x0c = particles.xpt - grid.xc[hosts]
    y0c = particles.ypt - grid.yc[hosts]

    # Get neighbouring elements
    n0=grid.nv[hosts, 0]
    n1=grid.nv[hosts, 1]
    n2=grid.nv[hosts, 2]

    def layer_var(layer):
        var_0 = (field[layer, n0]).flatten()
        var_1 = (field[layer, n1]).flatten()
        var_2 = (field[layer, n2]).flatten()

        var0 = grid.aw0[0, hosts] * var_0 + grid.aw0[1, hosts] * \
                var_1 + grid.aw0[2, hosts] * var_2
        varx = grid.awx[0, hosts] * var_0 + grid.awx[1, hosts] * \
                var_1 + grid.awx[2, hosts] * var_2
        vary = grid.awy[0, hosts] * var_0 + grid.awy[1, hosts] * \
                var_1 + grid.awy[2, hosts] * var_2
        var = var0  +  varx * x0c  +  vary * y0c

        return var, varx, vary

    # Calculate velocities for 2D or upper layer in 3D
    var, varx, vary = layer_var(layer)
    
    if len(field.shape)>1:
        # Interpolation lower layer
        var2, varx2, vary2 = layer_var(layer2)
        # Multiply by layer weights
        var = var * f1 + var2 * f2
        varx = varx * f1 + varx2 * f2
        vary = vary2 * f1 + vary2 * f2

    # Add code to not update land particles?

    return var, varx, vary
    
    
def __triinterpN(self, field, particles):
    """
    ** FVCOM interpolation method for node data using grid parameters.**

    Inputs:
      - self - pyticleClass
      - field - the field to use for interpolation
      - particles - the particles to get the depth/elevation for
    """

    grid = self.grid

    hosts = __find_hosts(grid, particles)

    # Get distance from element center
    x0c = particles.xpt - grid.xc[hosts]
    y0c = particles.ypt - grid.yc[hosts]

    # Get neighbouring elements
    n0=grid.nv[hosts, 0]
    n1=grid.nv[hosts, 1]
    n2=grid.nv[hosts, 2]

    var_0 = (field[n0]).flatten()
    var_1 = (field[n1]).flatten()
    var_2 = (field[n2]).flatten()

    var0 = grid.aw0[0, hosts] * var_0 + grid.aw0[1, hosts] * \
            var_1 + grid.aw0[2, hosts] * var_2
    varx = grid.awx[0, hosts] * var_0 + grid.awx[1, hosts] * \
            var_1 + grid.awx[2, hosts] * var_2
    vary = grid.awy[0, hosts] * var_0 + grid.awy[1, hosts] * \
            var_1 + grid.awy[2, hosts] * var_2
    var = var0  +  varx * x0c  +  vary * y0c

    # Add code to not update land particles?

    return var


def __find_hosts(grid, particles):
    """
    *Given an FVCOM grid and particles find the host elements of those particles*

    Inputs:
      - grid - an FVCOM grid
      - particles - the particles to get the depth/elevation for
    """

    # Find particles indomain
    indom = particles.indomain!=-1
    # Get node numbers for each particles element
    ns = grid.nv[particles.indomain[indom], :]
    # Check if each particle is in each element if not find the particles element
    intri = intriangle(grid.x[ns], grid.y[ns], particles.xpt[indom], \
                        particles.ypt[indom])

    if intri.all():
        hosts = particles.indomain
    else:
        # Figure out which particles aren't in there element and find their element
        check = indom*False
        check[indom]=~intri
        particles.indomain[check] = grid.finder.__call__(particles.xpt[check], \
                particles.ypt[check])
        hosts = particles.indomain

    return hosts


def scattered():
    pass
