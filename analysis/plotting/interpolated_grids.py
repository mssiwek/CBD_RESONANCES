
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

def plot(all_fps, figpath, titles, param, \
         vlims, \
         cbar_symlog=True, \
         box=10, ncells = 500, BoxSize=300, \
         dim='2D', SemiMajorAxis=1, fname='ebdot', \
         **kwargs):

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
            if param == 'rho':
                z = snap.Snap['PartType0']['Density']

            #interpolate grid
            snap.interp(x_interp, y_interp, z)
            #now stack the interpolated grid onto the last one
            grid = snap.stack_grid(grid=grid)
        
        snap.plot_interp(figpath = figpath, vlims=vlims, \
                         cbar_symlog=cbar_symlog, \
                         fname = fname + '_interp_%.3d' %i)
        title = r'$%.2f < t < %.2f$' %(t_initial/Pb, t_final/Pb)
        fname += '_stacked_grid' + '_%.2f<t<%.2f' %(t_initial/Pb, t_final/Pb)
        plot_utils.plot_mesh(X_interp, Y_interp, np.mean(grid, axis=-1), \
                             figpath = figpath, fname = fname,\
                             title=title, vlims=vlims)

            
if __name__ == "__main__":
    #cbd_ebqb_study_rs_videos
    fp_root = "/n/holylfs05/LABS/hernquist_lab/msiwek/arepo/cbd_ebqb_study_rs_videos/res=5.000e-10/alpha=0.100/dim=2/"
    figpath = "../figures/"

    param_dict = \
        {\
        'e': [0.8],
        'q': [0.1],
        'aspect_ratio': [0.10],
        'phi': [0.0],
        'AccretionRadius': [0.005],
        'SinkRate': [0.50],\
        'keys': ['e', 'q', 'aspect_ratio', 'phi', 'AccretionRadius', 'SinkRate']
        }
    all_fps = misc.get_fpaths_params(fp_root, param_dict, test_param = 'e')
    titles = []
    for eb in param_dict['e']:
        titles.append(r'$e_{\rm b} = $' + '%.2f' %eb)

    param = 'rho'
    if param=='rho':
        vlims=[1.e-7, 1.e-4]
        cbar_symlog=False
    else:
        vlims=[-1.e-8,1.e-8]
        cbar_symlog=True
    fname = param + '_eb=%.2f' %param_dict['e'][0] + '_qb=%.2f' %param_dict['q'][0]
    plot(all_fps, figpath + '/stacked_grids/', titles, param, vlims, \
         cbar_symlog = cbar_symlog, box=10, ncells=500, BoxSize=300, \
         SemiMajorAxis=1, fname=fname, i_initial=300, i_final=601)


            



