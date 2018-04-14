import random
import numpy as np
import matplotlib.pyplot as plt
import sys
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

def first_dilution(data, molecule_per_well, well, cycles, dilution_rate, bases_per_amplicon, error_rate):

    dilution_to = [int(molecule_per_well[i] * 2**cycles * dilution_rate) for i in range(well)]
    nsamples_per_molecule = []
    for i in range(well):
        nsamples_per_molecule += list(nsamples_generator(molecule_per_well[i], dilution_to[i]))
    
    indexes = expand(nsamples_per_molecule)
    data_expand = data.iloc[indexes].reset_index(drop=True)

    nsamples_per_molecule = np.array(nsamples_per_molecule, dtype=np.int16)
    nsamples_per_molecule = nsamples_per_molecule[nsamples_per_molecule > 0]
    
    mutations = poisson_mutation(nsamples_per_molecule, cycles, bases_per_amplicon, error_rate, tag="UID")

    if data_expand.shape[0] != mutations.shape[0]:
        raise(ValueError)
    data = pd.concat([data_expand, mutations], axis=1)
    data[:10000].to_csv("data.csv")

    print("        First dilution done")
    return data

def second_dilution(data, cycles, sequencer_reads, bases_per_amplicon, error_rate):

    n_amplicon = data.shape[0]
    nsamples = nsamples_generator(n_amplicon, int(sequencer_reads/2), restriction=2**cycles)
    
    indexes = expand(nsamples)
    data_expand = data.iloc[indexes].reset_index(drop=True)
    
    nsamples = nsamples[nsamples > 0]

    mutations = poisson_mutation(nsamples, cycles, bases_per_amplicon, error_rate, tag="WBC")

    if data_expand.shape[0] != mutations.shape[0]:
        raise(ValueError)
    
    data = pd.concat([data_expand, mutations], axis=1)

    print("        Second dilution done")
    return data

def one_run(UID_cycles, WBC_cycles, total_molecule, mutated_count, wells, first_dilution_rate, sequencer_reads, bases_per_amplicon, error_rate):

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
    data[:10000].to_csv("data0.csv")    

    data["UID"] = data["left_side"] * data["leftUID"] + (1 - data["left_side"]) * data["rightUID"]

    data["left_mutation"] = np.array(data["mutated_before"] | (data["WBC_left_mutation"] > 0) | ((data["left_side"] * data["UID_left_mutation"] + (1 - data["left_side"]) * data["UID_right_mutation"]) > 0), dtype=np.int8)
    data["right_mutation"] = np.array(data["mutated_before"] | (data["WBC_right_mutation"] > 0) | ((data["left_side"] * data["UID_left_mutation"] + (1 - data["left_side"]) * data["UID_right_mutation"]) > 0), dtype=np.int8)

    data = data.drop(columns = ["mutated_before","leftUID","rightUID","UID_left_mutation","UID_right_mutation","WBC_left_mutation","WBC_right_mutation","left_side"])
    
    data = pd.concat([data[data["leftWBC"]>0][["well_id","molecule_id","UID","left_mutation"]].rename(columns={"left_mutation":"mutation"}),
                      data[data["rightWBC"]>0][["well_id","molecule_id","UID","right_mutation"]].rename(columns={"right_mutation":"mutation"})], axis=0).sort_values(["well_id","molecule_id","UID"]).reset_index(drop=True)
    data[:10000].to_csv("data1.csv")

    # summary = [wells, total_molecule, mutated_count]
    # for w in range(wells):
    #     aux = data[data.well_id == wells]
        
    end = time.time()
    print("    one run time:", end-start)

if __name__ == '__main__':
    UID_cycles = 15 #15 #30
    bases_per_amplicon = 33
    error_rate = 1e-6

    total_molecule = 30 * 308 
    mutated_count = 1
    wells = 6 #6 #94
    first_dilution_rate = 0.01 #0.01 #0.000185*0.01

    WBC_cycles = 18 #18 #4
    sequencer_reads = 50000000

    # ["original_mutation", "n_molecule", "well_i_FC", "well_i_UID", "well_i_mutation", "well_i_UID", "well_i_mutation"]

    # print("Begin simulation for:", runs, "runs with cycles:", UID_cycles, WBC_cycles)

    # names = ["n_well", "n_molecule", "original_mutation"]
    # for i in range(wells):
    #     names.append( "well_" + str(i) + "_FC" )
    #     names.append( "well_" + str(i) + "_UID" )
    #     names.append( "well_" + str(i) + "_mutation" )
    #     names.append( "well_" + str(i) + "_UID_fd" )
    #     names.append( "well_" + str(i) + "_mutation_fd" )

    one_run(UID_cycles,WBC_cycles,total_molecule,mutated_count,wells,first_dilution_rate,sequencer_reads,bases_per_amplicon,error_rate)

    # one_run(15,15,900, 1, 6, 0.01, 5000000, bases_per_amplicon, error_rate)


