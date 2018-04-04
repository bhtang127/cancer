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
    
    return np.array(maintain, dtype=np.int32) 

def get_mutation(samples_id, cycles, bases_per_amplicon, error_rate):
    n = len(samples_id)
    mutation_count = np.zeros([n, 2], dtype = np.int16)

    for i in range(cycles):
        if i == 0:
            poisson_mean = n * bases_per_amplicon * error_rate
            total_count = np.random.poisson(poisson_mean)
            if total_count == 0: continue
            mutated_order = np.random.randint(n, size = (total_count,))
            for t in mutated_order:
                if samples_id[t] % 2 == 0:
                    mutation_count[t, 1] += 1
                else:
                    mutation_count[t, 0] += 1
        
        else:
            denominator = 2 ** (i-1)
            ancestors_to_kids = {}
            
            for id in samples_id:
                prev_ancestor = int(id / denominator)
                ancestor = int(prev_ancestor / 2)
                if (ancestor + prev_ancestor) % 2 == 0:
                    continue
                else:
                    if ancestor in ancestors_to_kids.keys():
                        ancestors_to_kids[ancestor].append(id)
                    else:
                        ancestors_to_kids[ancestor] = [id]

            keys = list(ancestors_to_kids)
            n_ancestors = len(keys) 
            poisson_mean = n_ancestors * bases_per_amplicon * error_rate
            total_count = np.random.poisson(poisson_mean)
            if total_count == 0: 
                continue
            mutated_order = np.random.randint(n_ancestors, size=(total_count,))
            for o in mutated_order:
                for mutated_id in ancestors_to_kids[ keys[o] ]:
                    order = np.where(samples_id == mutated_id)[0][0]
                    mutation_count[order, 0] += 1 
                    mutation_count[order, 1] += 1

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



    
    




