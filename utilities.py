from __future__ import print_function
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
        particles.lon1 = locations[:, 0]
        particles.lat1 = locations[:, 1]
        x, y = self.grid.proj(particles.lon1, particles.lat1)
        particles.x1 = x
        particles.y1 = y
    else:
        particles.x1 = x
        particles.y1 = y
        
    if '3D' in self.opt.gridDim:
        particles.z1 = locations[:,0]
    
    #particles.time1 = self.time.
    
    # Run interp code here to get particle velocities here
    particles.u1 = interpolate(self, 'u', particles)
    particles.v1 = interpolate(self, 'v', particles)   
    if '3D' in self.opt.gridDim:
        particles.w1 = interpolate(self, 'w', particles)    
    
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
           
    grid = container()        
    for key in options.reqvar:
        setattr(grid, key, data[key])   
    grid.finder = mplt.Triangulation(grid.x, grid.y, grid.nv).get_trifinder()
    
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
    
    return time
    
        
        
