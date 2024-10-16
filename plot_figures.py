import misc
import orbital_evolution_maps as oem

fp_root = "/n/holylfs05/LABS/hernquist_lab/msiwek/arepo/cbd_ebqb_study_rs_videos//res=5.000e-10/alpha=0.100/dim=2/"
figpath = "./figures/"

param_dict = \
    {\
    'e': [0.8],
    'q': [1.0],
    'aspect_ratio': [0.10],
    'phi': [0.0],
    'AccretionRadius': [0.03],
    'SinkRate': [0.50],\
    'keys': ['e', 'q', 'aspect_ratio', 'phi', 'AccretionRadius', 'SinkRate']
    }

all_fps = misc.get_fpaths_params(fp_root, param_dict, test_param = 'e')
titles = []
for eb in param_dict['e']:
    titles.append(r'$e_{\rm b} = $' + '%.2f' %eb)

param = 'abdot'
for param in ['ebdot']:
    fname = '%s_slow_motion' %param

    oem.plot(all_fps, figpath, titles, param, box=20, BoxSize=300, \
        dim='2D', SemiMajorAxis=1, fname=fname, contour=False, ratio = False,\
        vmin = -1, vmax=1, i_initial=1000, i_final=1500)

