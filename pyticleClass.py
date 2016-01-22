from __future__ import division, print_function
from defaults import *
from fileIO import *
from utilities import *
from solvers import *

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
        self.grid = set_grid(self, data)
        
        # Deal with time setup
        if debug: print(' Setup time')
        self.time = set_time(self)
           
        # Set initial particle data
        if debug: print(' Setting particle locations')
        self.particles = set_particles(self, locations)
            
        # Initialize output file
        if self.opt.saveOutput:
            if debug: print(' Initializing netcdf output')
            init_netcdf(self, outfile)
            self.opt.outfile = outfile       
            
        if debug: print(' pyticle initialized!')
        
        return
        
     
    def run(self):
        self.grid.u1 = self.grid.u[self.time.starttime]
        self.grid.v1 = self.grid.v[self.time.starttime]
        if '3D' in self.opt.gridDim:
            self.grid.w1 = self.grid.ww[self.time.starttime]
        
        # Progress counter
        cnt = 0
        
        for ncstep in range(self.time.starttime, self.time.endtime):
            for interpstep in range(self.time.interp):
                
                # Linearly interpolate fields in time
                f1 = (interpstep - self.time.interp + 1) / -self.time.interp
                f2 = (interpstep + 1) / self.time.interp
                
                self.grid.u2 = self.grid.u[ncstep,:]*f1 + self.grid.u[ncstep + 1,:]
                self.grid.v2 = self.grid.v[ncstep,:]*f1 + self.grid.v[ncstep + 1,:]
                if '3D' in self.opt.gridDim:
                    self.grid.w2 = self.grid.ww[ncstep,:]*f1 + self.grid.ww[ncstep + 1,:]
                
                # Move particles
                self.particles = rungekutta(self)
            
                # Overwrite old velocity field with current velocity field
                self.grid.u1 = self.grid.u2
                self.grid.v1 = self.grid.v2
                if '3D' in self.opt.gridDim:
                    self.grid.w1 = self.grid.w2
            
                self.particles.time = self.grid.time[ncstep] * f1 +\
                                      self.grid.time[ncstep + 1] * f2
                
                # Only save output when specified based on outputratio
                if np.mod(interpstep, self.time.out) == 0:
                    self.particles.loop += 1
                    save_netcdf(self)
                    
                cnt += 1
                                
                # The code starts at "step 2" as step one happens during initialization
                print('Completed step {}/{}'.format(cnt +1 , self.time.totalsteps))
            
            
        
        
