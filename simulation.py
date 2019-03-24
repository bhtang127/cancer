import random
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import pandas as pd
from random_methods import *


def data_init(molecule_per_well, mutated_per_well):
    data = []
    for i, nm in enumerate(molecule_per_well):
        for j in range(nm):
            data.append({"well_id": i,
                         "molecule_id": j,
                         "mutated_before": j < mutated_per_well[i]})
    data = pd.DataFrame(data)
    print("        Data Initialized")
    return data

def first_dilution(data, molecule_per_well, well, c_max, p, q, dilution_rate, bases_per_amplicon, error_rate):

    # dilution_to = [int(molecule_per_well[i] * 2**cycles * dilution_rate) for i in range(well)]
    nsamples_per_molecule, cycles = nsamples_generator(molecule_per_well, d_rate = dilution_rate, 
                                                       c_max=c_max, p=p, q=q)
    
    indexes = expand(nsamples_per_molecule)
    data_expand = data.iloc[indexes].reset_index(drop=True)

    nsamples_per_molecule = np.array(nsamples_per_molecule, dtype=np.int16)
    nsamples_per_molecule = nsamples_per_molecule[nsamples_per_molecule > 0]
    cycles = cycles[nsamples_per_molecule > 0]
    
    mutations = poisson_mutation(nsamples_per_molecule, cycles, bases_per_amplicon, error_rate, tag="UID")

    if data_expand.shape[0] != mutations.shape[0]:
        raise(ValueError)
    data = pd.concat([data_expand, mutations], axis=1)

    print("        First dilution done")
    return data

def second_dilution(data, well, c_max, p, q, sequencer_reads, bases_per_amplicon, error_rate):

    molecule_per_well = []
    for i in range[well]:
        n = data[data["well_id"] == i].shape[0]
        molecule_per_well.append(n)
    nsamples, cycles = nsamples_generator(molecule_per_well, d_to = int(sequencer_reads/2),
                                          c_max=c_max, p=p, q=q)
    
    indexes = expand(nsamples)
    data_expand = data.iloc[indexes].reset_index(drop=True)
    
    nsamples = nsamples[nsamples > 0]
    cycles = cycles[nsamples > 0]

    mutations = poisson_mutation(nsamples, cycles, bases_per_amplicon, error_rate, tag="WBC")

    if data_expand.shape[0] != mutations.shape[0]:
        raise(ValueError)
    
    data = pd.concat([data_expand, mutations], axis=1)

    print("        Second dilution done")
    return data

def one_run(UID_cycles, WBC_cycles, total_molecule, mutated_count, wells, first_dilution_rate, sequencer_reads, bases_per_amplicon, error_rate, cut_ratio=0):

    start = time.time()

    molecule_per_well = [int(total_molecule / wells) for i in range(wells)]
    for i in range(wells):
        if sum(molecule_per_well) == total_molecule:
            break
        else:
            molecule_per_well[i] += 1

    mutated_per_well = [0 for i in range(wells)]
    for i in range(mutated_count):
        loc = np.random.randint(wells)
        if mutated_per_well[loc] == molecule_per_well[loc]:
            continue
        else:
            mutated_per_well[loc] += 1

    init = data_init(molecule_per_well, mutated_per_well)

    fd_data = first_dilution(init, molecule_per_well, wells, UID_cycles, 
                                            first_dilution_rate, bases_per_amplicon, error_rate)

    data = second_dilution(fd_data, WBC_cycles, sequencer_reads, bases_per_amplicon, error_rate)
    # print( np.sum(data["leftUID"] == data["rightUID"]) / data.shape[0] )

    data["UID"] = data["left_side"] * data["leftUID"] + (1 - data["left_side"]) * data["rightUID"]

    data["left_mutation"] = np.array(data["mutated_before"] | (data["WBC_left_mutation"] > 0) | ((data["left_side"] * data["UID_left_mutation"] + (1 - data["left_side"]) * data["UID_right_mutation"]) > 0), dtype=np.int8)
    data["right_mutation"] = np.array(data["mutated_before"] | (data["WBC_right_mutation"] > 0) | ((data["left_side"] * data["UID_left_mutation"] + (1 - data["left_side"]) * data["UID_right_mutation"]) > 0), dtype=np.int8)

    data = data.drop(columns = ["mutated_before","leftUID","rightUID","UID_left_mutation","UID_right_mutation","WBC_left_mutation","WBC_right_mutation","left_side"])
    
    data = pd.concat([data[data["leftWBC"]>0][["well_id","molecule_id","UID","left_mutation"]].rename(columns={"left_mutation":"mutation"}),
                      data[data["rightWBC"]>0][["well_id","molecule_id","UID","right_mutation"]].rename(columns={"right_mutation":"mutation"})], axis=0).sort_values(["well_id","molecule_id","UID"]).reset_index(drop=True)

    if(cut_ratio == 0):
        data = data.groupby("well_id")['mutation'].agg([np.mean,"count"])
    else:
        data = data.groupby(["well_id","molecule_id","UID"])['mutation'].agg([np.mean, "count"])
        data["voted_mutate"] = np.array(data["mean"] > cut_ratio, dtype=np.int8)
        data = data.groupby("well_id")["voted_mutate"].agg([np.mean, "count"])
        
    end = time.time()
    print("    one run time:", end-start)
    return data

if __name__ == '__main__':
    wells = int(sys.argv[1])
    mutated_count = int(sys.argv[2])
    cut_off = float(sys.argv[3])

    bases_per_amplicon = 33
    error_rate = 1e-6
    total_molecule = 30 * 308 
    sequencer_reads = 50000000

    if wells == 94:
        UID_cycles = 30 
        first_dilution_rate = 0.01*0.000185 
        WBC_cycles = 4 
    elif wells == 6:
        UID_cycles = 15 
        first_dilution_rate = 0.01 
        WBC_cycles = 18
    else:
        raise(ValueError)

    data = one_run(UID_cycles,WBC_cycles,total_molecule, mutated_count, wells,first_dilution_rate,sequencer_reads,bases_per_amplicon,error_rate,cut_off)
    fname = "m"+str(mutated_count)+"w"+str(wells)+"c"+str(cut_off)+".csv"

    if not os.path.isfile(fname):
        data.to_csv(fname, mode='a')
    else:
        data.to_csv(fname, mode='a', header=False)



    # data = one_run(UID_cycles,WBC_cycles,total_molecule, 0, wells,first_dilution_rate,sequencer_reads,bases_per_amplicon,error_rate,0.9)
    # data.to_csv("m0w94c09.csv")

    # data = []
    # for i in range(10):
    #     data.append( one_run(15,18,total_molecule, 0, 6, 0.01,sequencer_reads,bases_per_amplicon,error_rate) )
    # data = pd.concat(data, axis=0)
    # data.to_csv("m0w6c0.csv")

