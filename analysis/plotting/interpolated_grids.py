
import numpy as np
import matplotlib.pyplot as plt
import misc
import Snap
import BinEvol as BE
import plot_utils
import h5py as h5

figwidth = 7
figheight = 6
fs = 20
ticksize = 8
tickwidth = 1.5
linewidth = 2

Pb = 2.*np.pi

def plot(all_fps, figpath, titles, param, box=10, ncells = 500, BoxSize=300, \
    dim='2D', SemiMajorAxis=1, fname='ebdot', **kwargs):

    i_initial, i_final = misc.snap_inds(all_fps, kwargs)

    sn = h5.File(all_fps[0] +  '/snap_%03d.hdf5' %i_initial, 'r')
    t_initial = sn['/Header'].attrs['Time']
    sn = h5.File(all_fps[0] +  '/snap_%03d.hdf5' %i_final, 'r')
    t_final = sn['/Header'].attrs['Time']

    x_interp = np.linspace(-box/2., box/2., ncells)
    y_interp = x_interp
    X_interp, Y_interp = np.meshgrid(x_interp, y_interp)

    for fp in all_fps:
        grid = None
        for i in range(i_initial, i_final+1):
            snap = Snap.Snap(fp, i, BoxSize=BoxSize)
            siminit = BE.SimInit(fp)
            binevol = BE.BinEvol(siminit, from_txt=False, from_snap=True, sn=snap.Snap)

            if param == 'ebdot':
                z = binevol.ebdot()
            if param == 'abdot':
                z = binevol.abdot()

            #interpolate grid
            snap.interp(x_interp, y_interp, z)
            #now stack the interpolated grid onto the last one
            grid = snap.stack_grid(grid=grid)

            # snap.plot(z, box=box, plot_cbar=True)
        
        snap.plot_interp(figpath = figpath, fname = fname + '_interp_%.3d' %i)
        title = r'$%.2f < t < %.2f$' %(t_initial/Pb, t_final/Pb)
        fname += '_stacked_grid' + '_%.2f<t<%.2f' %(t_initial/Pb, t_final/Pb)
        plot_utils.plot_mesh(X_interp, Y_interp, np.mean(grid, axis=-1), \
                             figpath = figpath, fname = fname,\
                             title=title)


            


            



