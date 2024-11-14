""" 
    LOAD SNAPSHOTS AND PERFORM VARIOUS MANIPULATIONS,
    E.G., INTERPOLATING AND TIME-AVERAGING.
    ALSO PLOT MAPS OF QUANTITIES SUCH AS;
    DENSITY, TORQUE, EBDOT, ABDOT, ...

    created by Magdalena Siwek, last modified 10/2024
"""

import numpy as np  
import torques_eb_lb as tq
import misc
import sys
import time
import gadget
import h5py as h5
import scipy
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
        
    def plot_interp(self, **kwargs):
        if 'ax' in kwargs:
            ax = kwargs['ax']
        else:
            fig = plt.figure(figsize = (figwidth, figheight))
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            ax = fig.add_subplot(1,1,1)
        
        v_range = 1.e-8
        cbar_norm_log = colors.SymLogNorm(linthresh=1.e-10, linscale=0.03, \
                                    vmin=-v_range, vmax=v_range)
        
        #norm=cbar_norm_log, \
        im = ax.pcolormesh(self.X_interp, self.Y_interp, self.Z_interp, norm=cbar_norm_log, \
                            cmap='RdBu_r')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        try:
            cbar.set_label(kwargs['cbar_label'], rotation=90, fontsize=fs)
        except:
            print("Warning: Did not supply cbar_label in kwargs.")

        ''' Add a marker for the sink particles '''
        ax = plot_utils.plot_sinks(ax, sn = self.Snap, size=0.1)

        ax.set_xlabel(r'$x/a_{\rm b}$', fontsize=fs)
        ax.set_ylabel(r'$y/a_{\rm b}$', fontsize=fs)

        plt.tight_layout()

        if 'ax' in kwargs:
            return(ax)
        else:
            misc.adjust_axes(fig, show_grid=False, ticksize=ticksize, fs=fs, tickwidth=tickwidth)
            plt.tight_layout()
            try:
                plt.savefig(kwargs['figpath'] + kwargs['fname'] + '.png')
            except:
                print("No figure path or figure name given! Saving in current path as mesh.png")
                plt.savefig('mesh.png')
                plt.close()


    def plot(self, z, box=10, plot_cbar=True, **kwargs):
        if 'ax' in kwargs:
            ax = kwargs['ax']
        else:
            fig = plt.figure(figsize = (figwidth, figheight))
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            ax = fig.add_subplot(1,1,1)
        
        v_range = 1.e-8
        cbar_norm_log = colors.SymLogNorm(linthresh=1.e-10, linscale=0.03, \
                                    vmin=-v_range, vmax=v_range)
        
        im = self.SnapSim.plot_Aslice(z, center = [150,150], \
                                contour=False, colorbar = False, \
                                axes = ax, box = box, \
                                norm=cbar_norm_log, \
                                vmin=-v_range, vmax=v_range, \
                                cmap='RdBu_r')
    
        ax.set_xlabel(r'$x/a_{\rm b}$')
        ax.set_ylabel(r'$y/a_{\rm b}$')

        if plot_cbar:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = plt.colorbar(im, cax=cax)
            try:
                cbar.set_label(kwargs['cbar_label'], rotation=90, fontsize=fs)
            except:
                print("Warning: Did not supply cbar_label in kwargs.")
        
        
        ''' Add a marker for the sink particles '''
        ax = plot_utils.plot_sinks(ax, sn = self.Snap, size=0.1)

        plt.tight_layout()

        if 'ax' in kwargs:
            return(ax)
        else:
            misc.adjust_axes(fig, show_grid=False, ticksize=ticksize, fs=fs, tickwidth=tickwidth)
            plt.tight_layout()
            try:
                plt.savefig(kwargs['figpath'] + kwargs['fname'] + '.png')
            except:
                print("No figure path or figure name given! Saving in current path as snap.png")
                plt.savefig('snap.png')
                plt.close()
            


        






    



