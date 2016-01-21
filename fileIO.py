import os
import sys
import netCDF4 as n4
import time

def init_netcdf(self,outfile):
    
            
    npts = len(self.particles.lon1)
    
    if os.path.isfile(outfile):
        print('File ' + outfile + ' already exists.\nPlease remove file and rerun.\n')
        sys.exit()
    else:
        ncid=n4.Dataset(outfile,'w',format='NETCDF4_CLASSIC')
        ncid.createDimension('time',None)        
        ncid.createDimension('npts',npts)

        ncid.createVariable('x', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)
        ncid.createVariable('y', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)
        ncid.createVariable('u', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)
        ncid.createVariable('v', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)
        ncid.createVariable('time', 'd', ('time'))
        
        if self.opt.useLL:        
            ncid.createVariable('lon', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)
            ncid.createVariable('lat', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)
            
        if '3D' in self.opt.gridDim:
            ncid.createVariable('z', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)
            ncid.createVariable('w', 'd', ('time','npts'), zlib=self.opt.zlib, least_significant_digit=self.opt.lsd)

        ncid.__setattr__('coordinateprojection', self.grid.projstr)
        ncid.__setattr__('history', 'Created on ' +time.ctime(time.time()) + ' by pyticle_tracker.' )
        ncid.__setattr__('options', str(self.opt.__dict__))

    #save initial particle positions and velocity
    ncid.variables['x'][0,:] = self.particles.x1
    ncid.variables['y'][0,:] = self.particles.y1
    ncid.variables['u'][0,:] = self.particles.u1
    ncid.variables['v'][0,:] = self.particles.v1

    if self.opt.useLL:       
        ncid.variables['lon'][0,:] = self.particles.lon1
        ncid.variables['lat'][0,:] = self.particles.lat1
        
    if '3D' in self.opt.gridDim:
        ncid.variables['z'][0,:] = self.particles.z1
        ncid.variables['w'][0,:] = self.particles.w1
        
    ncid.variables['time'][0] = self.particles.time1       

    ncid.close()

    return lag
