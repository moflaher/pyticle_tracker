from __future__ import print_function
import pyproj as pyp
from scipy.io import netcdf
import six
import matplotlib.tri as mplt
import numpy as np

# Make a dummy class to instantiate an object that is extendable
class container(object):
    pass
    
def set_particles(grid, locations, useLL, debug):
    
    particles = container()
    
    locations = np.atleast_2d(locations)
    
    if useLL:
        particles.lon1 = locations[:, 0]
        particles.lat1 = locations[:, 1]
        x, y = grid.proj(particles.lon1, particles.lat1)
        particles.x = x
        particles.y = y
    else:
        particles.x = x
        particles.y = y
    
    return particles           
    

def load_grid(data, options, debug):
    """ 
        Function to pass off loading of the grid to correct function
        depending on model type.
    """
    
    if 'FVCOM' in options.model:
        grid=_load_fvcom(data, options, debug)
        
    return grid
    

def _load_fvcom(data, options, debug):
    """ 
        Function to load fvcom input data.
    """
    
    if debug: print(' Loading FVCOM data')
    
    if isinstance(data, six.string_types):
        filepath = data
        data={}
        ncid = netcdf.netcdf_file(filepath, 'r',mmap=True)
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

        
        
def time_setup(grid, options, debug):
    time = container()
    
    return time
    
        
        
