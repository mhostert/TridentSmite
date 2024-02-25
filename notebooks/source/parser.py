import numpy as np
import scipy
import re
import os

from . import Cfv
from . import const


def get_material_mass(exp, myZ):
    m=exp.mass[(abs(exp.Z-myZ)<1)]
    if np.size(m)>1:
        print("More than one entry found for Z=", myZ)
        return 0
    elif np.size(m)==0:
        print("No entry found for Z=", myZ)
        return 0        
    else:
        return m


def load_files_to_dic(exp, identifier=False):
    dic={}
    dic['A'] = []
    dic['Z'] = []
    dic['PB'] = []#,dtype='bool'
    dic['cp'] = []
    dic['data'] = []
        
    # identifier of the events to load -- currently doing it by BSM params
    if not identifier and exp.new_physics:  
        identifier = f"gv_{exp.gprimeV:.6f}_ga_{exp.gprimeA:.6f}_mz_{exp.mzprime:.6f}.npy"
    else:  
        identifier = f"gv_{0.0:.6f}_ga_{0.0:.6f}_mz_{0.0:.6f}.npy"

    for fname in os.listdir("../"+exp.path_to_eventfiles):    # change directory as needed
        
        if identifier in fname:
            m = re.search('-_(.+?)_gv', fname)
            material = m.group(1)
           
            if material[0] == "p":
                dic['A'].append(1)
                dic['Z'].append(1)
            elif material[0] == "n":
                dic['A'].append(1)
                dic['Z'].append(0)
            else:
                dic['Z'].append(float(re.search('coh_(.+?)_', material).group(1)))
                dic['A'].append(float(re.search('_(.+?)', material).group(1)))
            

            if "bar" in fname:
                dic['cp'].append(-1)
            else:
                dic['cp'].append(1)


            if "noPB" in material:
                dic['PB'].append(False)
            elif "_PB" in material:
                dic['PB'].append(True)
            else:
                dic['PB'].append(False)

            d = np.load("../"+exp.path_to_eventfiles+fname)
            dic['data'].append(d.T)

    return dic

def combine_MC_outputs(dic,exp):
    new_dic={}
    
    new_dic['Enu'] = np.empty(0)
    new_dic['Q2'] = np.empty(0)
    new_dic['Pnu'] = np.empty((4,0))
    new_dic['Pplus'] = np.empty((4,0))
    new_dic['Pminus'] = np.empty((4,0))
    new_dic['PB'] = np.empty(0)
    new_dic['cp'] = np.empty(0)
    new_dic['A'] = np.empty(0)
    new_dic['Z'] = np.empty(0)
    new_dic['weights'] = np.empty(0)

    for i in range(np.size(dic['Z'])):
        # get from exp definition
        if dic['A'][i]>1:
            mass=get_material_mass(exp, dic['Z'][i])
        elif (dic['A'][i]==1 and not dic['PB'][i]):
            mass=get_material_mass(exp, dic['Z'][i])
        elif dic['A'][i]==1 and dic['PB'][i]:
            if dic['Z'][i]==1:
                mass=np.sum(exp.Z*exp.mass/exp.A)
            elif dic['Z'][i]==0:
                mass=np.sum( (exp.A-exp.Z) * exp.mass/exp.A)

        else:
            print(f"Error! Could not find normalization for case A={dic['A']}, Z={dic['Z']}.")

        nevents = np.size(dic['data'][i][0])

        new_dic['Z'] = np.concatenate([new_dic['Z'], dic['Z'][i]*np.ones(nevents)], axis=0)
        new_dic['A'] = np.concatenate([new_dic['A'], dic['A'][i]*np.ones(nevents)], axis=0)
        new_dic['PB'] = np.concatenate([new_dic['PB'], dic['PB'][i]*np.ones(nevents)], axis=0)
        new_dic['cp'] = np.concatenate([new_dic['cp'], dic['cp'][i]*np.ones(nevents)], axis=0)

        new_dic['Enu'] = np.concatenate([new_dic['Enu'], dic['data'][i][0]],axis=0)
        new_dic['Q2'] = np.concatenate([new_dic['Q2'], dic['data'][i][1]],axis=0)
        new_dic['weights'] = np.concatenate([new_dic['weights'], dic['data'][i][-1]/dic['A'][i]*mass],axis=0)

        new_dic['Pnu'] = np.concatenate([new_dic['Pnu'],(dic['data'][i][2:6])],axis=1)
        new_dic['Pplus'] = np.concatenate([new_dic['Pplus'],(dic['data'][i][6:10])],axis=1)
        new_dic['Pminus'] = np.concatenate([new_dic['Pminus'],(dic['data'][i][10:14])],axis=1)

    return new_dic

