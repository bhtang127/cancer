import random
import numpy as np
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
            print(poisson_mean, total_count)
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
                prev_ancestor = id / denominator
                ancestor = prev_ancestor / 2
                if (ancestor + prev_ancestor) % 2 == 0:
                    continue
                else:
                    if ancestor in ancestors_to_kids.keys():
                        ancestors_to_kids[ancestor].append(id)
                    else:
                        ancestors_to_kids[ancestor] = [id]

            keys = ancestors_to_kids.keys()
            n_ancestors = len(keys) 
            poisson_mean = n_ancestors * bases_per_amplicon * error_rate
            total_count = np.random.poisson(poisson_mean)
            print(poisson_mean, total_count)
            if total_count == 0: continue
            mutated_order = np.random.randint(n_ancestors, size=(total_count,))
            for o in mutated_order:
                mutated_id = ancestors_to_kids[ keys[o] ]
                order = np.where(samples_id == mutated_id)[0][0]
                mutation_count[order, 0] += 1 
                mutation_count[order, 1] += 1

    return mutation_count

if __name__ == '__main__':
    UID_cycles = 30
    bases_per_amplicon = 33
    error_rate = 1e-6

    first_dilution_to = 3973
    second_dilution_to = 2706

    is_mutated = 0

    SampleIDs = np.sort(sampleUID(UID_cycles, first_dilution_to))
    mutation_count = get_mutation(SampleIDs, UID_cycles, bases_per_amplicon, error_rate)
    
    mutated_uid_rate = np.sum(mutation_count > 0) / (2.0 * first_dilution_to)
    print(mutation_count, np.sum(mutation_count), mutated_uid_rate)

    
    




