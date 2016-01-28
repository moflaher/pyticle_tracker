from __future__ import division, print_function
from interpolation import *


def rungekutta(self):
    """
    ** Rungekutta4 solver. Moves particles one timestep.**
    
    Inputs:
      - pyticleClass
    """
    
    grid = self.grid
    particles = self.particles

    mstage = 4
    a_rk = [0, 0.5, 0.5, 1]
    b_rk = [1/6, 1/3, 1/3, 1/6]
    c_rk = [0, 0.5, 0.5, 1]

    chix = np.zeros((particles.npts, 4))
    chiy = np.zeros((particles.npts, 4))
    if '3D' in self.opt.gridDim:
        chiz = np.zeros((particles.npts, 4))

    for ns in range(1, mstage):
        particles.xpt  = particles.x  + (a_rk[ns] * self.time.dt) * chix[:, ns-1]        
        particles.ypt  = particles.y  + (a_rk[ns] * self.time.dt) * chiy[:, ns-1]
        if '3D' in self.opt.gridDim:
            particles.zpt  = particles.z  + (a_rk[ns] * self.time.dt) * chiz[:, ns-1]

        uin  = ((1-c_rk[ns]) * grid.u1 + c_rk[ns] * grid.u2)
        vin  = ((1-c_rk[ns]) * grid.v1 + c_rk[ns] * grid.v2)
        if '3D' in self.opt.gridDim:
            win  = ((1-c_rk[ns]) * grid.w1 + c_rk[ns] * grid.w2)
            
        usam = interpolate(self, uin, particles)
        vsam = interpolate(self, vin, particles)
        if '3D' in self.opt.gridDim:
            wsam = interpolate(self, win, particles)

        chix[:, ns] = usam
        chiy[:, ns] = vsam
        if '3D' in self.opt.gridDim:
            chiz[:, ns] = wsam
 
    particles.xpt=particles.x
    particles.ypt=particles.y
    if '3D' in self.opt.gridDim:
        particles.zpt  = particles.z
        
    for ns in range(0, mstage):
        particles.xpt = particles.xpt + self.time.dt * b_rk[ns] * chix[:,ns]
        particles.ypt = particles.ypt + self.time.dt * b_rk[ns] * chiy[:,ns]
        if '3D' in self.opt.gridDim:
            particles.zpt = particles.zpt + self.time.dt * b_rk[ns] * chiz[:,ns]

    # Set particles positions and velocities for this timestep
    particles.x = particles.xpt
    particles.u = interpolate(self, uin, particles)
    particles.y = particles.ypt
    particles.v = interpolate(self, vin, particles)
    if '3D' in self.opt.gridDim:
        particles.z = particles.zpt
        particles.w = interpolate(self, win, particles)
    particles.indomain = grid.finder.__call__(particles.x, particles.y)

    return particles

