import numpy as np
import misc
import ebdot_abdot_timeseries as eat
import ebdot_abdot_powerspec as eap
import interpolated_grids as ig
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
