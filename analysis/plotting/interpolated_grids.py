
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

def plot(all_fps, figpath, param, \
         vlims, \
         cbar_symlog=True, \
         make_video = True, \
         box=10, ncells = 500, BoxSize=300, \
         fname='ebdot', \
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
                cbar_label = r'$\dot{e}_{\rm b}$'
            if param == 'abdot':
                z = binevol.abdot()
                cbar_label = r'$\dot{a}_{\rm b}$'
            if param == 'rho':
                z = snap.Snap['PartType0']['Density']
                cbar_label = r'$\Sigma$'

            #interpolate grid
            snap.interp(x_interp, y_interp, z)
            #now stack the interpolated grid onto the last one
            grid = snap.stack_grid(grid=grid)
        
            if make_video:
                title = r'$t = %.2f\,P_{\rm b}$' %(snap.Snap['/Header'].attrs['Time']/Pb)
                plot_utils.plot_mesh(snap.X_interp, snap.Y_interp, snap.Z_interp, \
                             cbar_symlog=cbar_symlog, \
                             figpath = figpath, \
                             fname = fname + '_interp_%.3d' %i,\
                             title=title, \
                             vlims=vlims, \
                             cbar_label=cbar_label)
        
        title = r'$%.2f\,P_{\rm b} < t < %.2f\,P_{\rm b}$' %(t_initial/Pb, t_final/Pb)
        fname_stacked = fname + '_stacked_grid' + '_%.2f<t<%.2f' %(t_initial/Pb, t_final/Pb)
        plot_utils.plot_mesh(X_interp, Y_interp, np.mean(grid, axis=-1), \
                             cbar_symlog=cbar_symlog, \
                             figpath = figpath, \
                             fname = fname_stacked,\
                             title=title, \
                             vlims=vlims, \
                             cbar_label=cbar_label)
        
        if make_video:
            misc.make_video(figpath, fname + '_interp')

            
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

    param = 'ebdot'
    if param=='rho':
        vlims=[1.e-7, 1.e-4]
        cbar_symlog=False
    else:
        vlims=[-1.e-8,1.e-8]
        cbar_symlog=True
    fname = param + '_eb=%.2f' %param_dict['e'][0] + '_qb=%.2f' %param_dict['q'][0]
    plot(all_fps, figpath + '/stacked_grids/', param, vlims, \
         cbar_symlog = cbar_symlog, make_video=True, box=10, ncells=500, BoxSize=300, \
         SemiMajorAxis=1, fname=fname, i_initial=300, i_final=310)


            



