import random
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import pandas as pd
from random_methods import *


def data_init(molecule_per_well, mutated_per_well):
    start = time.time()

    data = []
    for i, nm in enumerate(molecule_per_well):
        for j in range(nm):
            data.append({"well_id": i,
                         "molecule_id": j,
                         "mutated_before": j < mutated_per_well[i]})
    data = pd.DataFrame(data)

    end = time.time()
    print("Data Prepare time:", end-start)
    return data

def first_dilution(data, molecule_per_well, well, cycles, dilution_rate, bases_per_amplicon, error_rate):
    start = time.time()

    dilution_to = [int(molecule_per_well[i] * 2**cycles * dilution_rate) for i in range(well)]
    nsamples_per_molecule = []
    for i in range(well):
        nsamples_per_molecule += list(nsamples_generator(molecule_per_well[i], dilution_to[i]))
    
    indexes = expand(nsamples_per_molecule)
    data_expand = data.iloc[indexes].reset_index(drop=True)

    nsamples_per_molecule = np.array(nsamples_per_molecule, dtype=np.int16)
    nsamples_per_molecule = nsamples_per_molecule[nsamples_per_molecule > 0]
    
    end = time.time()
    print("First dilution Prepare time:", end-start)
    
    mutations = poisson_mutation(nsamples_per_molecule, cycles, bases_per_amplicon, error_rate, tag="UID")

    if data_expand.shape[0] != mutations.shape[0]:
        raise(ValueError)
    data = pd.concat([data_expand, mutations], axis=1)

    end = time.time()
    print("First dilution time:", end-start)
    return data

def second_dilution(data, cycles, sequencer_reads, bases_per_amplicon, error_rate):
    start = time.time()

    n_amplicon = data.shape[0]
    nsamples = nsamples_generator(n_amplicon, int(sequencer_reads/2), restriction=2**cycles)
    
    well_tag = np.array(data.iloc[np.where(nsamples > 0)[0]]["well_id"])
    indexes = expand(nsamples)
    data_expand = data.iloc[indexes].reset_index(drop=True)
    
    nsamples = nsamples[nsamples > 0]
    end = time.time()
    print("Second dilution Prepare time:", end-start)
    mutations = poisson_mutation(nsamples, cycles, bases_per_amplicon, error_rate, tag="WBC")

    if data_expand.shape[0] != mutations.shape[0]:
        raise(ValueError)
    
    data = pd.concat([data_expand, mutations], axis=1)
    end = time.time()
    print("Second dilution time:", end-start)
    return data

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
    
    end = time.time()
    print("Allocation time:", end-start)

    init = data_init(molecule_per_well, mutated_per_well)

    summary_first_dilution = first_dilution(init, molecule_per_well, wells, UID_cycles, 
                                            first_dilution_rate, bases_per_amplicon, error_rate)

    summary_first_dilution.to_msgpack('first_dilution2.msg')
    # summary_first_dilution = pd.read_msgpack("first_dilution.msg")

    count = 0
    count += 2 * summary_first_dilution[summary_first_dilution["mutated_before"]==1].shape[0]
    nmdata = np.asarray(summary_first_dilution[summary_first_dilution["mutated_before"]==0][["UID_left_mutation","UID_right_mutation"]])
    count += np.sum(nmdata > 0)

    first_dilution_to = int( total_molecule * 2**UID_cycles * first_dilution_rate )
    print("mutated ratio:", count / (2.0*first_dilution_to))

    summary_second_dilution = second_dilution(summary_first_dilution, WBC_cycles, sequencer_reads, bases_per_amplicon, error_rate)
    summary_second_dilution.to_msgpack("second_dilution2.msg")
    
    count = 0
    count += 2 * summary_second_dilution[summary_second_dilution["mutated_before"]==1].shape[0]
    sliced = summary_second_dilution[summary_second_dilution["mutated_before"]==0]
    sliced_left = sliced[sliced["left_side"]==1]
    sliced_right = sliced[sliced["left_side"]==0]
    x = np.array(sliced_left[["WBC_left_mutation","WBC_right_mutation"]]) 
    x[:,0] += np.array(sliced_left["UID_left_mutation"])
    x[:,1] += np.array(sliced_left["UID_left_mutation"])
    y = np.array(sliced_right[["WBC_left_mutation","WBC_right_mutation"]]) 
    y[:,0] += np.array(sliced_right["UID_right_mutation"])
    y[:,1] += np.array(sliced_right["UID_right_mutation"])    
    count += np.sum(x > 0)
    count += np.sum(y > 0)

    print("mutated ratio:", count / sequencer_reads)

    # print(ancestor_to_kid(np.array([0,1,2,3,4,5,6,7]),0))
    #  x = pd.DataFrame([{"a":1,"b":2},{"a":2,"b":3}])
    #  y = pd.DataFrame([{"c":1,"d":2},{"c":2,"d":3}])