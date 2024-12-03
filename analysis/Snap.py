""" 
    LOAD SNAPSHOTS AND PERFORM VARIOUS MANIPULATIONS,
    E.G., INTERPOLATING AND TIME-AVERAGING.
    ALSO PLOT MAPS OF QUANTITIES SUCH AS;
    DENSITY, TORQUE, EBDOT, ABDOT, ...

    created by Magdalena Siwek, last modified 10/2024
"""

import numpy as np  
import misc
import gadget
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import plot_utils
from scipy.interpolate import LinearNDInterpolator

figwidth = 7
figheight = 6
fs = 20
ticksize = 8
tickwidth = 1.5
linewidth = 2

class Snap:
    def __init__(self, ext, i, BoxSize=300):
        Snap = h5.File(ext +  '/snap_%03d.hdf5' %i, 'r')
        self.Snap = Snap
        SnapSim = gadget.Simulation(ext + '/snap_%03d.hdf5' %i)
        self.SnapSim = SnapSim
        self.BoxSize = BoxSize

    def excise(self, x_range, y_range):
        ind_x_l = self.Snap['PartType0']['Coordinates'][:,0] >= x_range[0]
        ind_x_u = self.Snap['PartType0']['Coordinates'][:,0] <= x_range[1]
        ind_x = ind_x_l & ind_x_u

        ind_y_l = self.Snap['PartType0']['Coordinates'][:,1] >= y_range[0]
        ind_y_u = self.Snap['PartType0']['Coordinates'][:,1] <= y_range[1]
        ind_y = ind_y_l & ind_y_u

        ind_range = ind_x & ind_y

        """ NOW APPLY MASK TO ALL HYDRO QUANTITIES WE NEED """
        for key in self.Snap['PartType0'].keys():
            print("key = ", key)




    def interp(self, x, y, z):
        """ MAP THE VORONOI MESH TO THE SPECIFIED GRID (LIKELY CARTESIAN).
            x, y: COORDINATES OF NEW MESH 
            z: QUANTITY WE WANT TO MAP TO CARTESIAN MESH
            ALL QUANTITIES ARE INTERPOLATED USING RectBivariateSpline,
            NOTE THAT WE MUST USE kx=1, ky=1 TO ENSURE LINEAR INTERPOLATION """
        x_voronoi = self.Snap['PartType0']['Coordinates'][:,0] - self.BoxSize/2.
        y_voronoi = self.Snap['PartType0']['Coordinates'][:,1] - self.BoxSize/2.
        self.x_interp = x
        self.y_interp = y
        self.X_interp, self.Y_interp = np.meshgrid(self.x_interp, self.y_interp)
        interp = LinearNDInterpolator(list(zip(x_voronoi,y_voronoi)), z)
        self.z_interp = interp(self.X_interp, self.Y_interp)
        self.Z_interp = self.z_interp.reshape(len(self.x_interp), len(self.y_interp))
    
    def stack_grid(self, grid=None):
        if grid is None:
            #if no grid is supplied, we are starting with the stacking at this snapshot
            grid = self.Z_interp.reshape(len(self.Z_interp), len(self.Z_interp), 1)
        else:
            #if a grid is supplied, we are 
            Z_interp_reshape = self.Z_interp.reshape(len(self.Z_interp), len(self.Z_interp), 1)
            grid = np.dstack([grid, Z_interp_reshape])
        return(grid)
