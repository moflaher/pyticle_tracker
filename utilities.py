from __future__ import division, print_function
import pyproj as pyp
from scipy.io import netcdf
import six
import matplotlib.tri as mplt
import numpy as np
from interpolation import interpolate

# Make a dummy class to instantiate an object that is extendable
class container(object):
    pass
    
def set_particles(self, locations):
    
    particles = container()
    
    locations = np.atleast_2d(locations)
    
    if self.opt.useLL:
        particles.lon = locations[:, 0]
        particles.lat = locations[:, 1]
        x, y = self.grid.proj(particles.lon, particles.lat)
        particles.x = x
        particles.y = y
    else:
        particles.x = x
        particles.y = y
        
    particles.xpt = particles.x
    particles.ypt = particles.y           
        
    if '3D' in self.opt.gridDim:
        particles.z = locations[:,0]
        particles.zpt = particles.z
    
    particles.indomain = self.grid.finder.__call__(particles.x, particles.y)
    particles.time = self.time.time
        
    # Run interp code here to get particle velocities here
    particles.u = interpolate(self, self.grid.u[particles.time,], particles)
    particles.v = interpolate(self, self.grid.v[particles.time,], particles)   
    if '3D' in self.opt.gridDim:
        particles.w = interpolate(self, self.grid.ww[particles.time,], particles)    
        
    particles.npts = len(particles.x)
    particles.loop = 0
    
    return particles           
    

def set_grid(self, data):
    """ 
        Function to pass off loading of the grid to correct function
        depending on model type.
    """
    
    if 'FVCOM' in self.opt.model:
        grid=_load_fvcom(data, self.opt, self._debug)
        
    return grid
    

def _load_fvcom(data, options, debug):
    """ 
        Function to load fvcom input data.
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
    
    # Define the lcc projection
    xmax = np.nanmax(grid.lon)
    xmin = np.nanmin(grid.lon)
    ymax = np.nanmax(grid.lat)
    ymin = np.nanmin(grid.lat)
    xavg = ( xmax + xmin ) * 0.5;
    yavg = ( ymax + ymin ) * 0.5;
    ylower = ( ymax - ymin ) * 0.25 + ymin;
    yupper = ( ymax - ymin ) * 0.75 + ymin;
    
    grid.projstr = 'lcc +lon_0='+str(xavg)+' +lat_0='+str(yavg)+' +lat_1='+str(ylower)+' +lat_2='+str(yupper)
    grid.proj = pyp.Proj(proj=grid.projstr)
    
    return grid

        
        
def set_time(self):
    time = container()
    
    time.starttime = self.opt.starttime
    # Special handling of default case.
    if self.opt.endtime == -2:
        time.endtime = len(self.grid.time)-2
    else:
        time.endtime = self.opt.endtime

    time.interp = self.opt.interpolationratio
    time.out = self.opt.outputratio
    
    time.timestepin1 = time.starttime
    time.timestepin2 = time.starttime+1    
    time.timein1 = self.grid.time[time.starttime]
    time.timein2 = self.grid.time[time.starttime+1]  
    
    time.timestep1 = time.timestepin1
    time.timestep2 = time.timestepin1*((1-time.interp)/-time.interp) +\
                     time.timestepin2*(1/time.interp)    
    time.time = time.timein1
    time.time2 = time.timein1*((1-time.interp)/-time.interp) +\
                 time.timein2*(1/time.interp)   
    # Assuming days 
    # Need to improve time handling in general
    time.dt = (time.time2 - time.time) * 24*60*60
    
    time.timesteps = 1 + ((time.endtime - time.starttime) * time.interp/time.out)
    time.totalsteps = 1 + ((time.endtime - time.starttime) * time.interp)

    return time
    
        
        
