import random
import numpy as np
import matplotlib.pyplot as plt
import sys

def sampleUID(cycles, n_samples):
    ## cycles < 32
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
        loc = 0, s = 0
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
    n = len(samples_id)
    UID = np.array([ translate_aux(id, cycles) for id in samples_id], dtype = np.int32)
    return UID

def translate_well_aux(id, cycles):
    if id == 0:
        return [0, 1]
    elif 0 < id and id < 2**(cycles-1):
        return [1, 1]
    else:
        return translate_well_aux( id - 2**(cycles-1), cycles-1 )

def translate_well(samples_id, cycles):
    n = len(samples_id)
    well = np.array([ translate_well_aux(id, cycles) for id in samples_id], dtype = np.bool)
    return well

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
        order = 0
        for id in samples_id:
            kid = np.where(samples_id == id)[0][0]
            prev_ancestor = int(id / denominator)
            ancestor = int(prev_ancestor / 2)
            if (ancestor + prev_ancestor) % 2 == 0:
                continue
            else:
                if ancestor in ancestors:
                    ancestors_to_kids[order].append(kid)
                else:
                    ancestors_to_kids[order] = [kid] 
                    order += 1
                    ancestors.append(ancestor)


def get_mutation_per_molecule(samples_id, cycles, bases_per_amplicon, error_rate):
    mutation_count = np.zeros([n, 2], dtype = np.int16)

    [ poisson_mutation(bases_per_amplicon, error_rate, ancestor_to_kid(samples_id, i), mutation_count) for i in range(cycles) ]
        
    return mutation_count

if __name__ == '__main__':
    UID_cycles = 30
    bases_per_amplicon = 33
    error_rate = 1e-6

    total_molecule = 30 * 308
    mutated_count = 9
    wells = 94
    first_dilution_rate = 0.000185 * 0.01

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

    all_well_samples_per_molecule = []
    all_well_sampleUID = []
    all_well_mutations = []
    for i in range(wells):
        print("Processing Well:", i)
    
        molecule_this_well = molecule_per_well[i]
        mutated_this_well = mutated_per_well[i]

        first_dilution_to = int( molecule_this_well * 2**UID_cycles * first_dilution_rate )

        sampleUID_per_molecule = [0 for i in range(molecule_this_well)]
        for i in range(first_dilution_to):
            loc = np.random.randint(molecule_this_well)
            sampleUID_per_molecule[loc] += 1
        all_well_samples_per_molecule.append( sampleUID_per_molecule )
        
        samples = [ sampleUID(UID_cycles, s) for s in sampleUID_per_molecule ]
        all_well_sampleUID.append( samples )

        mutations = [ get_mutation(sample, UID_cycles, bases_per_amplicon, error_rate) for sample in samples ]
        all_well_mutations.append( mutations )

    count = 0
    for i in range(wells):
        mutation_count_this_well = mutated_per_well[i]
        molecule_this_well = molecule_per_well[i]
        mutations_this_well = all_well_mutations[i]
        nsamples_this_well = all_well_samples_per_molecule[i]
        for j in range(molecule_this_well):
            if j < mutation_count_this_well:
                count += 2 * nsamples_this_well[j]
            else:
                count += np.sum( mutations_this_well[j] > 0 ) 
    
    print(count / (2 * total_molecule * 2**30 * first_dilution_rate))



    
    




