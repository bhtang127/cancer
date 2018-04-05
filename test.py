import random
import numpy as np
import matplotlib.pyplot as plt
import sys

def sampleUID(cycles, n_samples):
    ## cycles < 32
    if n_samples == 0:
        return None

    init = np.random.randint(2**cycles)
    maintain = [init]

    while len(maintain) < n_samples:
        sample = np.random.randint(2**cycles)
        if sample in maintain:
            continue
        else:
            maintain.append(sample)
    
    return np.sort(np.array(maintain, dtype=np.int32))

def nsamples_generator(counts, total_samples):
    events = np.sum(counts)
    l = len(counts)
    nsamples = [0 for i in range(l)]
    for i in range(total_samples):
        sample = np.random.randint(events)
        loc = 0
        s = 0
        while s < sample:
            s += counts[loc]
            loc += 1
        nsamples[loc] += 1
    
    return nsamples

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
    UID = np.array([ translate_aux(id, cycles) for id in samples_id], dtype = np.int32)
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

def poisson_mutation(bases_per_amplicon, error_rate, ancestor_to_kid, mutation_count):
    if isinstance( ancestor_to_kid, list ):
        n = len(ancestor_to_kid)
        poisson_mean = n * bases_per_amplicon * error_rate
        total_count = np.random.poisson(poisson_mean)
        if total_count == 0: return
        mutated_order = np.random.randint(n, size = (total_count,))
        for t in mutated_order:
            if ancestor_to_kid[t] == 0:
                mutation_count[t, 1] += 1
            else:
                mutation_count[t, 0] += 1
        return

    elif isinstance(ancestor_to_kid, dict):
        n = len(ancestor_to_kid)
        poisson_mean = n * bases_per_amplicon * error_rate
        total_count = np.random.poisson(poisson_mean)
        if total_count == 0: return
        mutated_order = np.random.randint(n, size=(total_count,))
        for loc in mutated_order:
            for id in ancestor_to_kid[loc]:
                mutation_count[id, 0] += 1 
                mutation_count[id, 1] += 1
    
    else:
        raise(TypeError)

def ancestor_to_kid(samples_id, depth):
    if depth == 0:
        return [id % 2 for id in samples_id]
    else:
        denominator = 2 ** (depth-1)
        ancestors = []
        ancestors_to_kids = {}
        for i in range(len(samples_id)):
            id = samples_id[i]
            prev_ancestor = int(id / denominator)
            ancestor = int(prev_ancestor / 2)
            if (ancestor + prev_ancestor) % 2 == 0:
                continue
            else:
                if ancestor in ancestors:
                    ancestors_to_kids[ancestors.index(ancestor)].append(i)
                else:
                    ancestors.append(ancestor)
                    ancestors_to_kids[ancestors.index(ancestor)] = [i]           
        return ancestors_to_kids


def get_mutation_per_molecule(samples_id, cycles, bases_per_amplicon, error_rate):
    if samples_id is None:
        return None 

    n = len(samples_id)
    mutation_count = np.zeros([n, 2], dtype = np.int16)

    [ poisson_mutation(bases_per_amplicon, error_rate, ancestor_to_kid(samples_id, i), mutation_count) for i in range(cycles) ]
        
    return mutation_count

def mylen(l):
    if l is None:
        return None
    else:
        return len(l)

def simulation_per_well(well, n_molecules, n_mutated, cycles, dilution_rate, bases_per_amplicon, error_rate):
    print("Processing Well:", well)
    dilution_to = int( n_molecules * 2**UID_cycles * first_dilution_rate )
    nsamples_per_molecule = nsamples_generator(np.ones(n_molecules), dilution_to)
    samples_id = [sampleUID(cycles, num) for num in nsamples_per_molecule]
    summary = [{"UIDs": translate(ids, cycles),
                "n_samples": mylen(ids),
                "mutated_before": i < n_mutated,
                "molecule_id": i,
                "well_id": well,
                "mutations_count": get_mutation_per_molecule(ids, cycles, bases_per_amplicon, error_rate)} for i, ids in enumerate(samples_id)]
    return summary

if __name__ == '__main__':
    UID_cycles = 30
    bases_per_amplicon = 33
    error_rate = 1e-6

    total_molecule = 30 * 308
    mutated_count = 9
    wells = 94
    first_dilution_rate = 0.000185 * 0.01
    sequencer_reads = 50000000

    print("Begin0")
    molecule_per_well = [int(total_molecule / wells) for i in range(wells)]
    for i in range(wells):
        if sum(molecule_per_well) == total_molecule:
            break
        else:
            molecule_per_well[i] += 1

    print("Begin1")
    mutated_per_well = [0 for i in range(wells)]
    for i in range(mutated_count):
        loc = np.random.randint(wells)
        if mutated_per_well[loc] == molecule_per_well[loc]:
            continue
        else:
            mutated_per_well[loc] += 1

    print("Begin2")
    summary_first_dilution = \
    [simulation_per_well(w,n,m,UID_cycles,first_dilution_rate,bases_per_amplicon,error_rate)\
                        for w, n, m in zip(range(wells), molecule_per_well, mutated_per_well)]

    count = 0
    for i in range(wells):
        summary_this_well = summary_first_dilution[i]
        for d in summary_this_well:
            if d["UIDs"] != None and d["mutated_before"] == 1:
                count += 2 * d["n_samples"]
            else:
                count += np.sum( d["mutations_count"] > 0 ) 
    
    print(count / (2 * total_molecule * 2**30 * first_dilution_rate))

    # print(ancestor_to_kid(np.array([0,1,2,3,4,5,6,7]),0))



    
    




