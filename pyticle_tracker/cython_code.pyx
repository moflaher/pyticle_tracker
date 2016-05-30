from __future__ import division,print_function
import numpy as np
cimport numpy as np


def intriangle(np.ndarray[double, ndim=2] xt, np.ndarray[double, ndim=2] yt, np.ndarray[double, ndim=1] x0, np.ndarray[double, ndim=1] y0):
    """
    ** Function checks if points are in a specified triangle**
    
    Inputs:
      - xt, yt - n X 3 array defining the x,y locations of each triangles node
      - x0, y0 - The n points to check     
           
    """
    
    cdef np.ndarray f1 = (y0 - yt[:, 0]) * (xt[:, 1] - xt[:, 0]) -\
         (x0 - xt[:, 0]) * (yt[:, 1] - yt[:, 0])
    cdef np.ndarray f2 = (y0 - yt[:, 2]) * (xt[:, 0] - xt[:, 2]) -\
         (x0 - xt[:, 2]) * (yt[:, 0] - yt[:, 2])
    cdef np.ndarray f3 = (y0 - yt[:, 1]) * (xt[:, 2] - xt[:, 1]) -\
         (x0 - xt[:, 1]) * (yt[:, 2] - yt[:, 1])  
    
    cdef np.ndarray rt = ((f1 * f3 >= 0) & (f3 * f2 >= 0))
    
    return rt
