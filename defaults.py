from utilities import *


def model_options(options):
    
    if 'model' in options:
        model = options['model']
    else:
        model = 'FVCOM'
    
    if 'FVCOM' in model:
        defaults =  _fvcom_options(options)
    
    defaults.model = model
    
    return defaults
    
    
def _fvcom_options(options):
    defaults=container()
     
    defaults.useLL = True
    defaults.gridDim = '2D'
    defaults.layer='da'
    defaults.interpolation = 'triinterp'
    defaults.saveOutput = False
     
    for key in options:
        setattr(defaults, key, options[key])
     
    # Define required vars based on defaults and/or options
    defaults.reqvar = ['x', 'y']
     
    if defaults.useLL:
        defaults.reqvar += ['lon', 'lat']
    
    if ('2D' in defaults.gridDim) and ('da' in str(defaults.layer)):
        defaults.reqvar += ['ua', 'va']
    elif '2D' in defaults.gridDim:
        defaults.reqvar += ['u', 'v']
    elif '3D' in defaults.gridDim:
        defaults.reqvar += ['u', 'v', 'ww', 'siglay', 'siglev']
    else:
        print('gridDim must be 2D or 3D!')
    
    if 'triinterp' in defaults.interpolation:
        defaults.reqvar += ['nv', 'nbe', 'a1u', 'a2u', 'aw0', 'awx', 'awy']
        
    return defaults
        
