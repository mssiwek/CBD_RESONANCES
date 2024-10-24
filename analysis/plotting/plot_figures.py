import numpy as np
import misc
import ebdot_abdot_timeseries as eat
import ebdot_abdot_powerspec as eap
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

# for param in ['ebdot']:
#     fname = '%s_slow_motion_log' %param

#     oem.plot(all_fps, figpath, titles, param, box=20, BoxSize=300, \
#         dim='2D', SemiMajorAxis=1, fname=fname, contour=False, ratio = False,\
#         vmin = -1, vmax=1, i_initial=1000, i_final=1500)

Pb = 2.*np.pi
#compare 2000, 4000, 5000, 6000
tmin = 3000*Pb
tmax = 3100*Pb

#note: give xlim in units of Pb
ylim = [-5.e-6, 5.e-6]#[-5.e-6, 5.e-6]

#plot ebdot or abdot
param = 'ebdot'

# mass accretion included or not
accs = [True]
# linear momentum accretion included or not
f_accs = [True]

fname = 'timeseries_%s_qb%.2f' %(param,param_dict['q'][0])
for eb in param_dict['e']:
    fname += '_%.2feb' %eb 

fname = fname + "_accrad=%.3f" %(param_dict['AccretionRadius'][0]) \
        + '_%d<t<%d' %(tmin/Pb, tmax/Pb)

eat.plot(all_fps, figpath, titles, param, \
         accs, f_accs, \
         fname=fname, ylim=ylim, \
         tmin=tmin, tmax=tmax)

fname = 'powerspec_%s_qb%.2f' %(param,param_dict['q'][0])
for eb in param_dict['e']:
    fname += '_%.2feb' %eb 

fname = fname + "_accrad=%.3f" %(param_dict['AccretionRadius'][0]) \
        + '_%d<t<%d' %(np.ceil(tmin/Pb), np.ceil(tmax/Pb))


# mass accretion included or not
accs = [True]
# linear momentum accretion included or not
f_accs = [True]

xlim = [1.e-1, 3]
ylim = [1.e-17, 1.e-7]

eap.plot(all_fps, figpath, titles, param, \
         accs, f_accs, \
         fname=fname,  \
         tmin=tmin, tmax=tmax, \
         xlim=xlim, ylim=ylim)

