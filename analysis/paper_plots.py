import numpy as np
import os
import timeseries as ts
import powerspec as ps
import interpolated_grids as ig
import utils
Pb = 2*np.pi
#cbd_ebqb_study_rs_videos
fp_root = "/n/holylfs05/LABS/hernquist_lab/msiwek/arepo/cbd_ebqb_study/res=5.000e-10/alpha=0.100/dim=2/"
figpath = "../figures/fiducial/"

""" -------------------------------------------------------------- """
""" Plot timeseries of ebdot, showing that mbdot is insignificant,
    but f_acc can be important (in higher eccentricity binaries). """
""" -------------------------------------------------------------- """
fig1 = False
""" -------------------------------------------------------------- """
""" Plot timeseries of ebdot, showing that mbdot is insignificant,
    but f_acc can be important (in higher eccentricity binaries). """
""" -------------------------------------------------------------- """



""" -------------------------------------------------------------------------- """
""" Plot timeseries of ebdot, showing difference between cavity and outer CBD. """
""" -------------------------------------------------------------------------- """
fig2 = False
""" -------------------------------------------------------------------------- """
""" Plot timeseries of ebdot, showing difference between cavity and outer CBD. """
""" -------------------------------------------------------------------------- """



""" ------------------------------------------------------------------------------- """
"""  Plot power spectrum of ebdot, comparing simulations where ebdot>0 and ebdot<0  """
""" ------------------------------------------------------------------------------- """
fig3 = False
""" ------------------------------------------------------------------------------- """
"""  Plot power spectrum of ebdot, comparing simulations where ebdot>0 and ebdot<0  """
""" ------------------------------------------------------------------------------- """



""" ------------------------------------------------------------------------------- """
"""   Plot power spectrum of ebdot, comparing prograde and retrograde simulations   """
""" ------------------------------------------------------------------------------- """
fig4 = False
""" ------------------------------------------------------------------------------- """
"""   Plot power spectrum of ebdot, comparing prograde and retrograde simulations   """
""" ------------------------------------------------------------------------------- """



""" ------------------------------------------------------------------------------- """
"""     Plot time series of ebdot, comparing prograde and retrograde simulations    """
""" ------------------------------------------------------------------------------- """
fig5 = True
""" ------------------------------------------------------------------------------- """
"""     Plot time series of ebdot, comparing prograde and retrograde simulations    """
""" ------------------------------------------------------------------------------- """




""" -------------------------------------------------------------- """
""" Plot timeseries of ebdot, showing that mbdot is insignificant,
    but f_acc can be important (in higher eccentricity binaries). """
""" -------------------------------------------------------------- """
if fig1:
    fname = 'fig1'
    test_param = 'e'
    ebs = [0.6,0.8]
    param_dict = \
        {\
        'e': ebs,
        'q': [0.1],
        'aspect_ratio': [0.10],
        'phi': [0.0],
        'AccretionRadius': [0.03],
        'SinkRate': [0.50],\
        'keys': ['e', 'q', 'aspect_ratio', 'phi', 'AccretionRadius', 'SinkRate']
        }
    all_fps = utils.get_fpaths_params(fp_root, param_dict, test_param = test_param)

    titles = []
    for eb in ebs:
        titles.append(r'$e_{\rm b} = %.2f,\ q_{\rm b} = %.2f$' %(eb, param_dict['q'][0]))

    tmin = 5000*Pb
    tmax = 5003*Pb
    ylim = [-5.e-6, 5.e-6]
    param = 'ebdot'
    # mass accretion included or not
    accs = [False, True]
    # linear momentum accretion included or not
    f_accs = [False, True]
    # Include or exclude cavity?
    fg_cavs = [True]

    labels = []
    for i,acc in enumerate(accs):
        for j,f_acc in enumerate(f_accs):
            for l,fg_cav in enumerate(fg_cavs):
                label = ''                    
                if fg_cav:
                    if acc:
                        label += r'$\dot{M}_{\rm b} \geq 0 $'
                    else:
                        label += r'$\dot{M}_{\rm b} = 0 $'
                    if f_acc: 
                        label += r'$, \ f_{\rm acc} \ \rm{incl}$'
                    else:
                        label += r'$, \ f_{\rm acc} \ \rm{excl}$'
                else:
                    label += r'$ f_{\rm g, \, r>a}$'
                labels.append(label)
    
    ts.plot(all_fps, labels, param, \
        accs, f_accs, fg_cavs, \
        all_in_one_ax = False, \
        raw = True, mean=True, smoothed=True, \
        figpath = figpath, fname=fname, \
        tmin=tmin, tmax=tmax, ylim=ylim,\
        titles=titles)

""" -------------------------------------------------------------- """
""" Plot timeseries of ebdot, showing that mbdot is insignificant. """
""" -------------------------------------------------------------- """



""" -------------------------------------------------------------------------- """
""" Plot timeseries of ebdot, showing difference between cavity and outer CBD. """
""" -------------------------------------------------------------------------- """
if fig2:
    fname = 'fig2'
    test_param = 'phi'
    phis = [0]
    param_dict = \
        {\
        'e': [0.3],
        'q': [1.0],
        'aspect_ratio': [0.10],
        'phi': phis,
        'AccretionRadius': [0.03],
        'SinkRate': [0.50],\
        'keys': ['e', 'q', 'aspect_ratio', 'phi', 'AccretionRadius', 'SinkRate']
        }
    all_fps = utils.get_fpaths_params(fp_root, param_dict, test_param = test_param)

    titles = []
    for phi in phis:
        titles.append(r'$\Phi = %d, \ e_{\rm b} = %.2f,\ q_{\rm b} = %.2f$' %(phi, param_dict['e'][0], param_dict['q'][0]))

    tmin = 5000*Pb
    tmax = 5010*Pb
    ylim = [-300,300]#[-5.e-6, 5.e-6]
    param = 'ebdot'
    # mass accretion included or not
    accs = [False]
    # linear momentum accretion included or not
    f_accs = [False]
    # Include or exclude cavity?
    fg_cavs = [True, False]

    labels = []
    for i,acc in enumerate(accs):
        for j,f_acc in enumerate(f_accs):
            for l,fg_cav in enumerate(fg_cavs):
                label = ''                    
                if fg_cav:
                    if acc:
                        label += r'$\dot{M}_{\rm b} \geq 0 $'
                    else:
                        label += r'$\dot{M}_{\rm b} = 0 $'
                    if f_acc: 
                        label += r'$, \ f_{\rm acc} \ \rm{incl}$'
                    else:
                        label += r'$, \ f_{\rm acc} \ \rm{excl}$'
                else:
                    label += r'$ f_{\rm g, \, r>a}$'
                labels.append(label)
    
    ts.plot(all_fps, labels, param, \
        accs, f_accs, fg_cavs, \
        all_in_one_ax = False, \
        raw = True, mean=True, smoothed=True, \
        figpath = figpath, fname=fname, \
        tmin=tmin, tmax=tmax, ylim=ylim,\
        titles=titles)
""" -------------------------------------------------------------------------- """
""" Plot timeseries of ebdot, showing difference between cavity and outer CBD. """
""" -------------------------------------------------------------------------- """



""" ------------------------------------------------------------------------------- """
"""  Plot power spectrum of ebdot, comparing simulations where ebdot>0 and ebdot<0  """
""" ------------------------------------------------------------------------------- """
if fig3:
    fname = 'fig3' 

    ylim = [1.e-15, 1.e-3]
    xlim = [0.1, 10]
    tmin = 5000*Pb
    tmax = 10000*Pb
    param = 'ebdot'

    test_param = 'q'
    param_dict = \
        {\
        'e': [0.3],
        'q': [0.1,1.0],
        'aspect_ratio': [0.10],
        'phi': [0.0],
        'AccretionRadius': [0.03],
        'SinkRate': [0.50],\
        'keys': ['e', 'q', 'aspect_ratio', 'phi', 'AccretionRadius', 'SinkRate']
        }
    all_fps = utils.get_fpaths_params(fp_root, param_dict, test_param = test_param)

    # mass accretion included or not
    accs = [False]
    # linear momentum accretion included or not
    f_accs = [False]
    # Include or exclude cavity? True: include cavity, False: exclude cavity
    fg_cavs = [True]

    labels = []
    for tp in param_dict[test_param]:
        if test_param == 'phi':
            labels.append(r'$\Phi = $' + '%d' %(tp))
        if test_param == 'e':
            labels.append(r'$e_{\rm b} = $' + '%.2f' %(tp))
        if test_param == 'q':
            labels.append(r'$q_{\rm b} = $' + '%.2f' %(tp))
    
    title = r'$e_{\rm b} = %.2f$' %(param_dict['e'][0])
    
    ps.plot(all_fps, labels, param, \
        accs, f_accs, fg_cavs, \
        all_in_one_ax = True, \
        figpath = figpath, fname=fname, \
        tmin=tmin, tmax=tmax, ylim=ylim, xlim=xlim, \
        title=title)
""" ------------------------------------------------------------------------------- """
"""  Plot power spectrum of ebdot, comparing simulations where ebdot>0 and ebdot<0  """
""" ------------------------------------------------------------------------------- """



""" ------------------------------------------------------------------------------- """
"""   Plot power spectrum of ebdot, comparing prograde and retrograde simulations   """
""" ------------------------------------------------------------------------------- """
if fig4:
    fname = 'fig4'

    ylim = [1.e-15, 1.e-3]
    xlim = [0.4, 10]
    tmin = 1500*Pb
    tmax = 2000*Pb
    param = 'ebdot'
    
    test_param = 'phi'
    param_dict = \
        {\
        'e': [0.5],
        'q': [1.0],
        'aspect_ratio': [0.10],
        'phi': [0.0,180],
        'AccretionRadius': [0.03],
        'SinkRate': [0.50],\
        'keys': ['e', 'q', 'aspect_ratio', 'phi', 'AccretionRadius', 'SinkRate']
        }
    all_fps = utils.get_fpaths_params(fp_root, param_dict, test_param = test_param)

    # mass accretion included or not
    accs = [False]
    # linear momentum accretion included or not
    f_accs = [False]
    # Include or exclude cavity? True: include cavity, False: exclude cavity
    fg_cavs = [False]

    labels = []
    for tp in param_dict[test_param]:
        if test_param == 'phi':
            labels.append(r'$\Phi = $' + '%d' %(tp))
        if test_param == 'e':
            labels.append(r'$e_{\rm b} = $' + '%.2f' %(tp))
        if test_param == 'q':
            labels.append(r'$q_{\rm b} = $' + '%.2f' %(tp))
    
    title = r'$e_{\rm b} = %.2f, q_{\rm b} = %.2f $' %(param_dict['e'][0],param_dict['q'][0])
    
    ps.plot(all_fps, labels, param, \
        accs, f_accs, fg_cavs, \
        all_in_one_ax = True, \
        figpath = figpath, fname=fname, \
        tmin=tmin, tmax=tmax, ylim=ylim, xlim=xlim, \
        title=title)
""" ------------------------------------------------------------------------------- """
"""   Plot power spectrum of ebdot, comparing prograde and retrograde simulations   """
""" ------------------------------------------------------------------------------- """



""" ------------------------------------------------------------------------------- """
"""     Plot time series of ebdot, comparing prograde and retrograde simulations    """
""" ------------------------------------------------------------------------------- """
if fig5:
    fname = 'fig5'
    test_param = 'phi'
    phis = [0,180]
    param_dict = \
        {\
        'e': [0.5],
        'q': [1.0],
        'aspect_ratio': [0.10],
        'phi': phis,
        'AccretionRadius': [0.03],
        'SinkRate': [0.50],\
        'keys': ['e', 'q', 'aspect_ratio', 'phi', 'AccretionRadius', 'SinkRate']
        }
    all_fps = utils.get_fpaths_params(fp_root, param_dict, test_param = test_param)

    tmin = 2000*Pb
    tmax = 2010*Pb
    ylim = [-300,300]#[-5.e-6, 5.e-6]
    param = 'ebdot'
    # mass accretion included or not
    accs = [True]
    # linear momentum accretion included or not
    f_accs = [True]
    # Include or exclude cavity?
    fg_cavs = [True]

    labels = []
    for tp in param_dict[test_param]:
        if test_param == 'phi':
            labels.append(r'$\Phi = $' + '%d' %(tp))
        if test_param == 'e':
            labels.append(r'$e_{\rm b} = $' + '%.2f' %(tp))
        if test_param == 'q':
            labels.append(r'$q_{\rm b} = $' + '%.2f' %(tp))
    
    title = r'$ e_{\rm b} = %.2f,\ q_{\rm b} = %.2f$' %(param_dict['e'][0], param_dict['q'][0])

    ts.plot(all_fps, labels, param, \
        accs, f_accs, fg_cavs, \
        all_in_one_ax = True, \
        raw = True, mean=True, smoothed=True, \
        figpath = figpath, fname=fname, \
        tmin=tmin, tmax=tmax, ylim=ylim,\
        title=title)
""" ------------------------------------------------------------------------------- """
"""     Plot time series of ebdot, comparing prograde and retrograde simulations    """
""" ------------------------------------------------------------------------------- """