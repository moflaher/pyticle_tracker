from __future__ import division, print_function
from .utilities import *


def model_options(options):
    """
    ** Select which default options to run with given specified model type.
       Also loads options that user specified.**

    Inputs:
      - options - dict given when initializing pyticleClass
    """

    if 'model' in options:
        model = options['model']
    else:
        model = 'FVCOM'

    if 'FVCOM' in model:
        defaults =  _fvcom_options(options)

    defaults.model = model

    return defaults


def _fvcom_options(options):
    """
    ** Setup default and user specified options for FVCOM **

    Inputs:
      - options - dict given when initializing pyticleClass
    """
    defaults=container()

    defaults.useLL = True
    defaults.gridDim = '2D'
    defaults.layer = 0
    defaults.interpolation = 'triinterp'
    defaults.saveOutput = True
    defaults.ncformat = 'NETCDF3_64BIT'
    defaults.awgn = False

    # Specifiy the proj projstr so that the code does not default to an LCC project.
    defaults.projstr=[]

    # Using zlib has a huge performance impact as it has to recompress every save
    defaults.zlib = False
    defaults.lsd = None

    # Default time runs for whole run.
    defaults.starttime = 0
    defaults.endtime = -2

    # Input timestep is divided by interpolationratio
    defaults.interpolationratio = 12

    # Output timestep every x many interpolationratio
    # EX: input time is 1 hour, interpolation ratio is 20 (3 minutes)
    #     outputratio is 2 then data is saved every 6 minutes.
    defaults.outputratio = 1

    # Stop particle motion when within X of bottom
    defaults.cutoff = 0.01

    # Stop particle vertical motion if water depth is less them X
    defaults.limitv = 0.01
    
    # Default particle sinking speed (postive up, m/s)
    defaults.sinkspeed = 0.0

    # Load lon/lat data from txt files if true    
    defaults.useLLtxt = False
    
    # Specify files to use for lon and lat
    defaults.lonTxt = ''
    defaults.latTxt = ''


    for key in options:
        setattr(defaults, key, options[key])

    # Define required vars based on defaults and/or options
    defaults.reqvar = ['x', 'y', 'xc', 'yc', 'time', 'nele', 'node']

    if defaults.useLL:
        defaults.reqvar += ['lon', 'lat']

    if ('2D' in defaults.gridDim) and ('da' in str(defaults.layer)):
        defaults.reqvar += ['ua', 'va']
    elif '2D' in defaults.gridDim:
        defaults.reqvar += ['u', 'v']
    elif '3D' in defaults.gridDim:
        defaults.reqvar += ['u', 'v', 'ww', 'siglay', 'h', 'zeta']
    else:
        print('gridDim must be 2D or 3D!')

    if 'triinterp' in defaults.interpolation:
        defaults.reqvar += ['nv', 'nbe', 'a1u', 'a2u', 'aw0', 'awx', 'awy']

    return defaults

