from __future__ import division, print_function
import os
import sys
import netCDF4 as n4
import time

def init_netcdf(self,outfile):
    """ 
    ** Initialize an ncfile in which to save output and save first step.**
    
    Inputs:
      - self - pyticleClass
      - outfile -  filepath and name
    """
   
    if os.path.isfile(outfile):
        print('File ' + outfile + ' already exists.\nPlease remove file and rerun.\n')
        sys.exit()
    else:
        ncid=n4.Dataset(outfile, 'w', format=self.opt.ncformat)
        ncid.createDimension('time', self.time.timesteps)        
        ncid.createDimension('npts', self.particles.npts)

        ncid.createVariable('x', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)
        ncid.createVariable('y', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)
        ncid.createVariable('u', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)
        ncid.createVariable('v', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)
        ncid.createVariable('indomain', 'i4', ('time','npts'), zlib=self.opt.zlib)
        ncid.createVariable('time', 'd', ('time'), zlib=self.opt.zlib)
        
        if self.opt.useLL:        
            # Add 3 to least significant digit for long lat
            ncid.createVariable('lon', 'd', ('time','npts'), zlib=self.opt.zlib)
            ncid.createVariable('lat', 'd', ('time','npts'), zlib=self.opt.zlib)
            
        if '3D' in self.opt.gridDim:
            ncid.createVariable('z', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)
            ncid.createVariable('w', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)

        ncid.__setattr__('coordinateprojection', self.grid.projstr)
        ncid.__setattr__('history', 'Created on ' +time.ctime(time.time()) + ' by pyticle_tracker.' )
        ncid.__setattr__('options', str(self.opt.__dict__).replace("'", '+'))       


    #save initial particle positions and velocity
    ncid.variables['x'][0,:] = self.particles.x
    ncid.variables['y'][0,:] = self.particles.y
    ncid.variables['u'][0,:] = self.particles.u
    ncid.variables['v'][0,:] = self.particles.v

    if self.opt.useLL:       
        ncid.variables['lon'][0,:] = self.particles.lon
        ncid.variables['lat'][0,:] = self.particles.lat
        
    if '3D' in self.opt.gridDim:
        ncid.variables['z'][0,:] = self.particles.z
        ncid.variables['w'][0,:] = self.particles.w
        
    ncid.variables['indomain'][0,:] = self.particles.indomain
    ncid.variables['time'][0] = self.particles.time       



    ncid.close()

    return


def save_netcdf(self):
    """ 
    ** Save the current timestep particle data to an ncfile.**
    
    Inputs:
      - self - pyticleClass
    """
    
    loop = self.particles.loop

    ncid=n4.Dataset(self.opt.outfile,'r+', format=self.opt.ncformat)

    #save particle positions and velocities
    ncid.variables['x'][loop, :]=self.particles.x
    ncid.variables['y'][loop, :]=self.particles.y
    if self.opt.useLL:
        lon, lat = self.grid.proj(self.particles.x, self.particles.y, inverse=True)
        ncid.variables['lon'][loop, :] = lon
        ncid.variables['lat'][loop, :] = lat
    ncid.variables['u'][loop, :] = self.particles.u
    ncid.variables['v'][loop, :] = self.particles.v
    
    if '3D' in self.opt.gridDim:
        ncid.variables['z'][loop,:] = self.particles.z
        ncid.variables['w'][loop,:] = self.particles.w
    
    ncid.variables['indomain'][loop,:] = self.particles.indomain
    ncid.variables['time'][loop] = self.particles.time

    ncid.close()

    return
