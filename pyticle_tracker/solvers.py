from __future__ import division, print_function
from interpolation import *
from interpolation import __find_hosts

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

    chix[:, 0] = particles.u
    chiy[:, 0] = particles.v
    if '3D' in self.opt.gridDim:
        chiz[:, 0] = particles.w

    if self.opt.diffusion:
        ff = particles.fudgefactor
        diffh = np.zeros((particles.npts))
        diffx = np.zeros((particles.npts))
        diffy = np.zeros((particles.npts))
        if '3D' in self.opt.gridDim:
            diffv = np.zeros((particles.npts))
            diffz = np.zeros((particles.npts))

        diffh = particles.viscofhp / ff
        diffx = particles.viscofhx / ff
        diffy = particles.viscofhy / ff
        if '3D' in self.opt.gridDim:
            diffv = particles.khp / ff
            diffz = particles.khz / ff

    # Loop over RK stages
    for ns in range(1, mstage):
        # Update particle positions at stage n
        if self.opt.diffusion:
            particles.xpt = particles.x + (a_rk[ns]*self.time.dti) * chix[:,ns-1] \
                    + diffx[:,ns-1] + np.sqrt(2*diffh[:,ns-1]) * \
                    particles.wiener[(4*(particles.count - 1)) + ns - 1]
            particles.ypt = particles.y + (a_rk[ns]*self.time.dti) * chiy[:,ns-1] \
                    + diffy[:,ns-1] + np.sqrt(2*diffh[:,ns-1]) * \
                    particles.wiener[(4*(particles.count - 1)) + ns - 1]
            if '3D' in self.opt.gridDim:
                particles.zpt = particles.z + (a_rk[ns]*self.time.dti) \
                    * chiz[:,ns-1] + diffz[:,ns-1] + np.sqrt(2*diffv[:,ns-1]) * \
                    particles.wiener[(4*(particles.count - 1)) + ns - 1]
        else:
            particles.xpt  = particles.x + (a_rk[ns]*self.time.dt) * chix[:, ns-1]
            particles.ypt  = particles.y + (a_rk[ns]*self.time.dt) * chiy[:, ns-1]
            if '3D' in self.opt.gridDim:
                particles.zpt  = particles.z + (a_rk[ns]*self.time.dt) * \
                                    chiz[:, ns-1]

        # Update velocity and elevation fields
        uin  = ((1-c_rk[ns]) * grid.u1 + c_rk[ns] * grid.u2)
        vin  = ((1-c_rk[ns]) * grid.v1 + c_rk[ns] * grid.v2)

        if self.opt.diffusion:
            viscofhin = ((1-c_rk[ns])*grid.viscofh1 + c_rk[ns]*grid.viscofh2)
            khin = ((1-c_rk[ns])*grid.kh1 + c_rk[ns]*grid.kh2)

        if '3D' in self.opt.gridDim:
            win = ((1-c_rk[ns]) * grid.w1 + c_rk[ns] * grid.w2)
            zin = ((1-c_rk[ns]) * grid.z1 + c_rk[ns] * grid.z2)

            #Find particle height
            particles.hpt = interpolate(self, grid.h, particles)[0]
            particles.ept = interpolate(self, zin, particles)[0]

            # If particles are above the water place them in the water
            particles.zpt = np.min([particles.zpt, particles.ept], axis=0)

            # If a particles is within cutoff (default-1cm) of bottom stop movement
            particles.inwater = (particles.zpt + particles.hpt) > self.opt.cutoff
            # And if they are at the bottom put them at the bottom not below
            particles.zpt = np.max([particles.zpt, -particles.hpt], axis=0)

            # Finally update the sigma position of the particle
            # for layer interpolation of the velocity
            particles.sigpt = np.divide(particles.zpt, \
                                        - 1 * (particles.hpt + particles.ept))

        usam = interpolate(self, uin, particles)[0]
        vsam = interpolate(self, vin, particles)[0]
        if '3D' in self.opt.gridDim:
            wsam = interpolate(self, win, particles)[0]
            t,x,y=interpolate(self,grid.kh[10,:,:],particles)
        #if self.opt.diffusion:
            # interpolation goes here

        chix[:, ns] = usam
        chiy[:, ns] = vsam
        if '3D' in self.opt.gridDim:
            chiz[:, ns] = wsam

            # If the particle is in shallow water then limit the vertical motion
            # Default value is 1 cm
            chiz[(particles.hpt + particles.ept) < self.opt.limitv, ns] = 0

    particles.xpt=particles.x
    particles.ypt=particles.y
    if '3D' in self.opt.gridDim:
        particles.zpt  = particles.z

    for ns in range(0, mstage):
        if self.opt.diffusion:
            particles.xpt = particles.xpt + (self.time.dt*b_rk[ns] * (chix[:,ns] \
                    + diffx[:,ns])) + particles.indomain[:] * (np.sqrt(2 * \
                    diffh[:,ns]) * particles.wiener[4*(particles.count - 1) + ns])
            particles.ypt = particles.ypt + (self.time.dt*b_rk[ns] * (chiy[:,ns] \
                    + diffy[:,ns])) + particles.indomain[:] * (np.sqrt(2 * \
                    diffh[:,ns]) * particles.wiener[4*(particles.count - 1) + ns])
            if '3D' in self.opt.gridDim:
                particles.zpt = particles.zpt + (self.time.dt*b_rk[ns] * \
                    (chiz[:,ns] + diffz[:,ns])) + particles.indomain[:] * \
                    (np.sqrt(2 * diffv[:,ns]) * \
                    particles.wiener[4*(particles.count - 1) + ns])
        else:
            particles.xpt = particles.xpt + self.time.dt * b_rk[ns] * chix[:,ns]
            particles.ypt = particles.ypt + self.time.dt * b_rk[ns] * chiy[:,ns]
            if '3D' in self.opt.gridDim:
                particles.zpt = particles.zpt + self.time.dt * b_rk[ns]*chiz[:,ns]

    # Set particles positions and velocities for this timestep
    # Unless the particle is on the bottom
    particles.x[particles.inwater] = particles.xpt[particles.inwater]
    particles.y[particles.inwater] = particles.ypt[particles.inwater]
    particles.u = interpolate(self, uin, particles)[0]
    particles.v = interpolate(self, vin, particles)[0]
    if '3D' in self.opt.gridDim:
        particles.z[particles.inwater] = particles.zpt[particles.inwater]
        particles.w = interpolate(self, win, particles)[0]

    #if self.opt.diffusion:
        # interpolation goes here!

    #particles.indomain = grid.finder.__call__(particles.x, particles.y)
    particles.indomain = __find_hosts(grid, particles)

    return particles

