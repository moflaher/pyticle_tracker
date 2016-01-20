

def init_netcdf(self,outfile):
    npts=len(locs)
    
    if os.path.isfile('output/'+filename):
        print('File '+ 'output/'+filename+ ' already exists.\nPlease remove file and rerun.\n')
        sys.exit()
    else:
        ncid=n4.Dataset('output/'+filename,'w',format='NETCDF3_CLASSIC')
        ncid.createDimension('time',None)        
        ncid.createDimension('npts',npts)

        ncid.createVariable('x','d',('time','npts'))
        ncid.createVariable('y','d',('time','npts'))
        ncid.createVariable('lon','d',('time','npts'))
        ncid.createVariable('lat','d',('time','npts'))
        ncid.createVariable('u','d',('time','npts'))
        ncid.createVariable('v','d',('time','npts'))

        ncid.createVariable('time','d',('time'))

        ncid.__setattr__('coordinateprojection',grid['projstr'])
        ncid.__setattr__('history','Created on ' +time.ctime(time.time()) + ' by scatter tracker.' )
        
    

    lag={}
    #last step locations
    lag['xl'],lag['yl']=grid['proj'](locs[:,0],locs[:,1])
    lag['lonl']=locs[:,0]
    lag['latl']=locs[:,1]

    #current step locations
    lag['xn']=np.empty((npts,))
    lag['yn']=np.empty((npts,))
    lag['lonn']=np.empty((npts,))
    lag['latn']=np.empty((npts,))

    #time stuff
    lag['timel']=np.empty((1,))
    lag['timen']=np.empty((1,))
    lag['timestep']=24*3600*(grid['time'][1]-grid['time'][0])/settings['irange'] 
    lag['loop']=0

    #other
    lag['u']=itp.gridded(grid,grid['ustep1'],lag['xl'],lag['yl'])
    lag['v']=itp.gridded(grid,grid['vstep1'],lag['xl'],lag['yl'])
    lag['npts']=npts
    lag['filename']='output/'+filename
    lag['proj']=grid['proj']


    #save initial particle positions and velocity
    ncid.variables['lon'][0,:]=locs[:,0]
    ncid.variables['lat'][0,:]=locs[:,1]
    ncid.variables['x'][0,:]=lag['xl']
    ncid.variables['y'][0,:]=lag['yl']
    ncid.variables['u'][0,:]=lag['u']
    ncid.variables['v'][0,:]=lag['u']

    ncid.close()

    return lag
