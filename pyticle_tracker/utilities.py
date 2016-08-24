from __future__ import division, print_function
import pyproj as pyp
from scipy.io import netcdf
import six
import matplotlib.tri as mplt
import numpy as np
import sys
from interpolation import interpolate

# Make a dummy class to instantiate an object that is extendable
class container(object):
    pass


def _set_grid(self, data, locations):
    """
    ** Initializes grid by loading required fields from
       data for the particular model being used.**

    Inputs:
      - data can be a dict with required fields or a path to an ncfile with the required fields
    """

    if 'FVCOM' in self.opt.model:
        grid=__load_fvcom(data, self.opt, locations, self._debug)

    return grid


def _set_time(self):
    """
    ** Initializes all information related to time variables.**

    Inputs:
      - pyticleClass
    """

    time = container()

    time.starttime = self.opt.starttime
    # Special handling of default case.
    if self.opt.endtime == -2:
        time.endtime = len(self.grid.time)-2
    else:
        time.endtime = self.opt.endtime

    time.interp = self.opt.interpolationratio
    time.out = self.opt.outputratio

    time.time = self.grid.time[time.starttime]

    # Assuming days
    # Need to improve time handling in general
    time.dt = 24*60*60 * (self.grid.time[1]-self.grid.time[0]) / time.interp

    time.timesteps = 1 + ((time.endtime - time.starttime) * time.interp/time.out)
    time.totalsteps = 1 + ((time.endtime - time.starttime) * time.interp)

    return time

def _randm_locs(self):
    """
    *** Generates n random locations in a circle of radius
    r centered at a coordinate ***

    Inputs:
      - pyticleClass (tuple coordinate, radius (float), and number of particles)
    """

    # note: 1 deg lon ~ 111117 m
    #       1 def lat ~ 79868 m

    if self.opt.useLL:
        lon = np.random.uniform(self.opt.centre[0] - self.opt.radius/79868.0, \
                      self.opt.centre[0] + self.opt.radius/78868.0, self.opt.np)
        lat = np.random.uniform(self.opt.centre[1] - self.opt.radius/111117.0, \
                      self.opt.centre[1] + self.opt.radius/111117.0, self.opt.np)

        return lon, lat

    else:
        pass


def _set_particles(self, locations):
    """
    ** Initializes particles location, velocity, and additional fields.**

    Inputs:
      - pyticleClass and an array of particle locations
    """

    particles = container()

    if not self.opt.randomStart:
        locations = np.atleast_2d(locations)

    if self.opt.useLL:
        if self.opt.randomStart:
            particles.lon, particles.lat = _randm_locs(self)
        else:
            particles.lon = locations[:, 0]
            particles.lat = locations[:, 1]
        x, y = self.grid.proj(particles.lon, particles.lat)
        particles.x = x
        particles.y = y
    else:
        # if self.opt.randomStart:
        #     particles.x, particles.y = _randm_locs()
        # else:
        particles.x = locations[:, 0]
        particles.y = locations[:, 1]

    particles.xpt = particles.x
    particles.ypt = particles.y

    particles.indomain = self.grid.finder.__call__(particles.x, particles.y)
    if np.sum(particles.indomain != -1) == 0:
        print('No particles are initially in the domain.')
        sys.exit()

    if '3D' in self.opt.gridDim:
        if self.opt.randomStart:
            particles.z = np.tile(-5, self.opt.np)
        else:
            particles.z = locations[:,2]
        particles.zpt = particles.z

        #Find particle height
        particles.hpt = interpolate(self, self.grid.h, particles)[0]
        particles.ept = interpolate(self, self.grid.zeta[self.time.starttime,], \
                                    particles)[0]

        # If particles are above the water place them in the water
        particles.zpt = np.min([particles.zpt, particles.ept], axis=0)

        # If a particles is within cutoff (default - 1cm) of bottom stop movement
        particles.inwater = (particles.zpt + particles.hpt) > self.opt.cutoff
        # And if they are at the bottom put them at the bottom not below
        particles.zpt = np.max([particles.zpt, -particles.hpt], axis=0)

        # Finally update the sigma position of the particle
        # for layer interpolation of the velocity
        particles.sigpt = np.divide(particles.zpt, \
                                    -1*(particles.hpt + particles.ept))
    else:
        # If 2D then particles are always *vertically* in the water column
        particles.inwater = (particles.x * 0 + 1).astype(bool)

    particles.time = self.time.time
    particles.npts = len(particles.x)

    # Run interp code here to get particle velocities here
    particles.u = interpolate(self, self.grid.u[self.time.starttime,], particles)[0]
    particles.v = interpolate(self, self.grid.v[self.time.starttime,], particles)[0]

    if '3D' in self.opt.gridDim:
        particles.w = interpolate(self, self.grid.ww[self.time.starttime,], \
                                    particles)[0]

    if self.opt.diffusion:
        particles.viscofhp, particles.viscofhx, particles.viscofhy = \
            interpolate(self, self.grid.viscofh[self.time.starttime,], particles)[:3]
        # Some times this can be negative because of the interpolation method, force positive. Not sure how to properly handle this.
        particles.viscofhp = np.fabs(particles.viscofhp)
        particles.khp, _, _, particles.khz = \
            interpolate(self, self.grid.kh[self.time.starttime,], particles)
        particles.randomstate = self.opt.seed
        np.random.seed(self.opt.seed)
        particles.wiener = np.sqrt(self.time.dt)*np.random.randn(4 * \
                self.time.totalsteps, 1)
        particles.fudgefactor = self.opt.diffusionfudgefactor

    # progress counter
    particles.loop = 0
    particles.count = 1

    return particles


def __load_fvcom(data, options, locations, debug):
    """
    ** Loaded required grid data for FVCOM**

    Inputs:
      - data a dict or path to ncfile
      - options pyticleClass.opt
      - debug pyticleClass._debug
    """

    if debug: print(' Loading FVCOM data')

    if isinstance(data, six.string_types):
        filepath = data
        data={}
        ncid = netcdf.netcdf_file(filepath, 'r', mmap=True)
        for key in ncid.variables.keys():
            data[key] = ncid.variables[key].data
        if ('nv' in data):
            data['nv']=data['nv'].astype(int).T-1
        if ('nbe' in data):
            data['nbe']=data['nbe'].astype(int).T-1
        if ('nele' in ncid.dimensions):
            data['nele'] = ncid.dimensions['nele']
        if ('node' in ncid.dimensions):
            data['node'] = ncid.dimensions['node']

    # Initialize empty object
    grid = container()

    # Add the raw data as hidden
    grid._data = data

    # Load the proper data
    for key in options.reqvar:
        setattr(grid, key, data[key])

    # Initialize Traiangulation and finder
    grid.trigridxy = mplt.Triangulation(grid.x, grid.y, grid.nv)
    grid.finder = grid.trigridxy.get_trifinder()
    if options.useLL:
        grid.trigrid = mplt.Triangulation(grid.lon, grid.lat, grid.nv)

    # Handle 2D cases
    if ('2D' in options.gridDim) and ('da' in str(options.layer)):
        grid.u = grid.ua
        grid.v = grid.va
    elif '2D' in options.gridDim:
        grid.u = grid.u[:,options.layer,:]
        grid.v = grid.v[:,options.layer,:]
    else:
        # Create an array of siglay the size of the number of particles
        # This is so the interpolation code can find the particles layer
        grid.siglay = -1*grid.siglay[:,0]
        if options.randomStart:
            npts = options.np
        else:
            npts = len(locations[:,0])
        grid.siglaylen = len(grid.siglay)
        grid.siglayrep = grid.siglay.repeat(npts).reshape(grid.siglaylen, npts)
        # Same for siglev
        grid.siglev = -1*grid.siglev[:,0]
        grid.siglevlen = len(grid.siglev)
        grid.siglevrep = grid.siglev.repeat(npts).reshape(grid.siglevlen, npts)

    if (options.useLL) and (options.projstr==[]):
        # Define the lcc projection
        xmax = np.nanmax(grid.lon)
        xmin = np.nanmin(grid.lon)
        ymax = np.nanmax(grid.lat)
        ymin = np.nanmin(grid.lat)
        xavg = ( xmax + xmin ) * 0.5;
        yavg = ( ymax + ymin ) * 0.5;
        ylower = ( ymax - ymin ) * 0.25 + ymin;
        yupper = ( ymax - ymin ) * 0.75 + ymin;

        grid.projstr = 'lcc +lon_0='+str(xavg)+' +lat_0='+str(yavg) \
                + ' +lat_1='+str(ylower)+' +lat_2='+str(yupper)

    if options.useLL:
        grid.proj = pyp.Proj(proj=options.projstr)

    if options.diffusion:
        grid.hele = np.sum(grid.h[grid.nv],axis=1)

    return grid
