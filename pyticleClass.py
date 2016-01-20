from __future__ import print_function
from defaults import *
from fileIO import *
from utilities import *

class pyticle:
    """
    ** A class to hold the particles that will be tracked**
    
    Inputs:
      - A dict or netcdf filename, either of which contain the required data.
      - A list of locations to place particles
      - A filepath to save the output
      - Additional dict specifying options
    """
    
    def __init__(self, data, locations, outfile='test.nc', options={}, debug=False):
        """ Initialize pyticle class"""

        if debug: print('-Debug mode on-')
        self._debug=debug         
    
        # Load and set the options for the specified model type
        if debug: print(' Setting options')
        self.opt = model_options(options)  
            
        # Load the required variables for specified model type
        if debug: print(' Loading model grid')
        self.grid = load_grid(data, self.opt, debug)
           
        # Set initial particle data
        if debug: print(' Setting particle locations')
        self.particles = set_particles(self.grid, locations, self.opt.useLL, debug)
        
        # Deal with time setup
        if debug: print(' Setup time')
        self.time = time_setup(self.grid, self.opt, debug)
            
        # Initialize output file
        if self.opt.saveOutput:
            if debug: print(' Initializing netcdf output')
            init_netcdf(self, outfile, debug)
            self.opt.outfile = outfile       
            
        if debug: print(' pyticle initialized!')
        
        return
        
            
        
            
            
        
        
