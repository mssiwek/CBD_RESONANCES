import os
import pickle as pkl
import BinEvol as BE
import pickle

def load_binevol(fp, acc, f_acc, fg_cav):
    binevol_name = 'BinEvoltxt'
    if fg_cav:
        if acc:
            binevol_name += '_acc'
        else:
            binevol_name += '_no_acc'
        if f_acc:
            binevol_name += '_f_acc'
        else:
            binevol_name += '_no_f_acc'
    else:
        binevol_name += '_cav_excl'
    
    filename = binevol_name + '.pkl'
    
    if not os.path.exists(fp):
        os.makedirs(fp)
    
    if os.path.isfile(fp + filename):
        filehandler = open(fp + filename, 'rb')
        binevol = pkl.load(filehandler)
    else:
        siminit = BE.SimInit(fp + '/../', acc=acc, f_acc=f_acc, fg_cav=fg_cav)
        binevol = BE.BinEvol(siminit, from_txt=True, from_snap=False)
        filehandler = open(fp + filename, 'wb') 
        pkl.dump(binevol, filehandler)
        
    return(binevol)


def load_param_txt(fp):
    fp_list = fp.split('/')
    if fp_list[-1] == 'output':
        fp += '..'
    try:
        all_param_fp = open(fp + '/all_param.pkl', 'rb')
        all_param = pickle.load(all_param_fp)
        return(all_param)
    except FileNotFoundError:
        print("no param.pkl file found in %s" %fp)
        try:
            all_param_fp = open(fp + '/../' + '/all_param.pkl', 'rb')
            all_param = pickle.load(all_param_fp)
            return(all_param)
        except FileNotFoundError:
            from pathlib import Path
            path = Path(fp)
            fp2 = str(path.parent.absolute())
            print("fp2 = %s" %fp2)
            fo = open(fp2 + '/param.txt', "r")
        except FileNotFoundError:
            print("No param.txt in %s" %(fp2))
            try:
                fo = open(fp + '/param.txt', "r")
            except FileNotFoundError:
                print("No param.txt in %s or %s" %(fp2, fp))
    keys = []
    all_param = {}
    for line in fo.readlines():
        a = line.split()
        if len(a) > 1:
            key = a[0]
            key_val = a[1]
            keys.append(key)
            try:
                new_key_val = float(key_val)
            except ValueError:
                new_key_val = key_val
            all_param[key] = new_key_val
    return(all_param)