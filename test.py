import random
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import pandas as pd

def sampleUID(cycles, n_samples):
    ## cycles < 32
    if n_samples == 0:
        return None

    while True:
        ids = np.array(np.random.choice(2**cycles, n_samples), dtype=np.int32)
        if len(set(ids)) == n_samples:
            return np.sort(ids)

def nsamples_generator(counts, total_samples):
    pval = np.array(counts) / np.sum(counts)
    return np.random.multinomial(total_samples, pval)

def translate_aux(id, cycles):
    if cycles == 0:
        return [0, 0]
    if id % 2 == 1:
        uid = translate_aux( int(id / 2), cycles - 1 )[1]
        return [uid, uid]
    else:
        uid1 = translate_aux( int(id / 2), cycles - 1 )[0]
        uid2 = id / 2 + 2 ** (cycles - 1)
        return [uid1, uid2]

def translate(samples_id, cycles):
    if samples_id is None:
        return None
    n = len(samples_id)
    UID = np.array([ translate_aux(id, cycles) for id in samples_id ], dtype = np.int32)
    return UID

def translate_well_aux(id, cycles, well):
    if id == 0:
        return [0, well]
    elif 0 < id and id < 2**(cycles-1):
        return [well, well]
    else:
        return translate_well_aux( id - 2**(cycles-1), cycles-1, well )

def translate_well(samples_id, cycles, well):
    if samples_id is None:
        return None
    n = len(samples_id)
    translation = np.array([ translate_well_aux(id, cycles, well) for id in samples_id ], dtype = np.int16)
    return translation

def poisson_mutation(bases_per_amplicon, error_rate, influence_set, mutation_count):
    if isinstance( influence_set, list ):
        n = len(influence_set)
        poisson_mean = n * bases_per_amplicon * error_rate
        total_count = np.random.poisson(poisson_mean)
        if total_count == 0: return
        mutated_order = np.random.randint(n, size = (total_count,))
        for t in mutated_order:
            if influence_set[t] == 0:
                mutation_count[t, 1] += 1
            else:
                mutation_count[t, 0] += 1

    elif isinstance(influence_set, tuple):
        n = influence_set[0]
        ancestors = influence_set[1]
        is_influencial = influence_set[2]
        poisson_mean = n * bases_per_amplicon * error_rate
        total_count = np.random.poisson(poisson_mean)
        if total_count == 0: return
        mutated_order = np.random.randint(n, size=(total_count,))
        influence_list = np.sort(list(set(is_influencial * (ancestors + 1)) - set({0})))
        for loc in mutated_order:
            ancestor_plus_one = influence_list[loc]
            indexes = np.where( (ancestors+1) == ancestor_plus_one )[0]
            for i in indexes:
                mutation_count[i, 0] += 1 
                mutation_count[i, 1] += 1
    
    else:
        raise(TypeError)
    return

def influence(samples_id, depth):
    if depth == 0:
        return [id % 2 for id in samples_id]
    else:
        denominator = 2 ** (depth-1)
        prev_ancestors = np.asarray( np.floor(samples_id / denominator), dtype=np.int32 )
        ancestors = np.asarray( np.floor(prev_ancestors / 2), dtype=np.int32 )
        is_influencial = np.mod(ancestors + prev_ancestors, 2)

        influence_set = set(is_influencial * (ancestors + 1)) - set({0})
        count = len(influence_set)      
        return (count, ancestors, is_influencial)


def get_mutation_per_molecule(samples_id, cycles, bases_per_amplicon, error_rate):
    if samples_id is None:
        return None 

    n = len(samples_id)
    mutation_count = np.zeros([n, 2], dtype = np.int16)

    [ poisson_mutation(bases_per_amplicon, error_rate, influence(samples_id, i), mutation_count) for i in range(cycles) ]
    
    return mutation_count

def mylen(l):
    if l is None:
        return None
    else:
        return len(l)

def simulation_per_well(well, n_molecules, n_mutated, cycles, dilution_rate, bases_per_amplicon, error_rate):
    start = time.time()
    print("Processing Well:", well)
    dilution_to = int( n_molecules * 2**UID_cycles * first_dilution_rate )
    nsamples_per_molecule = nsamples_generator(np.ones(n_molecules), dilution_to)
    samples_id = [sampleUID(cycles, num) for num in nsamples_per_molecule]
    summary = []
    for i, ids in enumerate(samples_id):
        UIDs = translate(ids, cycles)
        mutations_count = get_mutation_per_molecule(ids, cycles, bases_per_amplicon, error_rate)
        if UIDs is not None:
            for uid, mutation in zip(UIDs, mutations_count):
                summary.append( {"leftUID": uid[0],
                                 "rightUID": uid[1],
                                 "mutated_before": i < n_mutated,
                                 "molecule_id": i,
                                 "well_id": well,
                                 "left_mutation_count": mutation[0],
                                 "right_mutation_count": mutation[1]} )
    summary = pd.DataFrame(summary)
    end = time.time()
    print("Well Processing time:", end-start)
    return summary

def second_dilution(summary, cycles, sequencer_reads, bases_per_amplicon, error_rate):
    start = time.time()

    n_wells = len(summary)
    molecule_per_well = nsamples_generator(np.ones(n_wells), sequencer_reads / 2)
    for i, n_molecules in enumerate(molecule_per_well):
        summary_this_well = summary[i]
        n_samples = [ s["n_samples"] for s in summary_this_well if s["n_samples"] is not None ]


    end = time.time()
    print("Second dilution time:", end-start)

if __name__ == '__main__':
    UID_cycles = 30
    bases_per_amplicon = 33
    error_rate = 1e-6

    total_molecule = 30 * 308
    mutated_count = 9
    wells = 94
    first_dilution_rate = 0.000185 * 0.01

    WBC_cycles = 4
    sequencer_reads = 50000000

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

    summary_first_dilution = \
    [simulation_per_well(w,n,m,UID_cycles,first_dilution_rate,bases_per_amplicon,error_rate)\
                        for w, n, m in zip(range(wells), molecule_per_well, mutated_per_well)]
    summary_first_dilution = pd.concat(summary_first_dilution)

    summary_first_dilution.to_csv("first_dilution_data.csv")

    count = 0
    for i, row in summary_first_dilution.iterrows():
        if row["mutated_before"] == 1:
            count += 2
        else:
            count += np.sum( row["mutation_count"] > 0 ) 
    
    print(count / (2 * total_molecule * 2**30 * first_dilution_rate))

    # print(ancestor_to_kid(np.array([0,1,2,3,4,5,6,7]),0))



    
    




