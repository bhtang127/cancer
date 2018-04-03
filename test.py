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

def sample_mutation_list(cycles, bases_per_amplicon, error_rate):
    # list of elements like [n, m, loc] 
    # means mutation happens in cycle n, mth amplicon and in loc base  

    number_of_duplications = np.array([2**(c+1) for c in range(cycles)])
    possion_means = error_rate * number_of_duplications * bases_per_amplicon
    n_mutaions = np.random.poisson(possion_means)

    mutations = [[c+1, np.random.randint(2**(c+1)), np.random.randint(33)] for c in range(cycles) for mut in range(n_mutaions[c])]

    return mutations

def influence_bond(location, current_cycles, total_cycles):
    # current_cycles < total_cycles
    next_influence = 2 * location + (location % 2 == 0)
    lower_bound = next_influence
    higher_bound = next_influence
    for i in range(total_cycles - current_cycles - 1):
        lower_bound *= 2
        higher_bound = higher_bound * 2 + 1
    return [lower_bound, higher_bound]

def within(value, bond):
    if value <= bond[1] and value >= bond[0]:
        return 1
    else:
        return 0

if __name__ == '__main__':
    UID_cycles = 30
    bases_per_amplicon = 33
    error_rate = 1e-6

    first_dilution_to = 3973
    second_dilution_to = 2706

    is_mutated = 0

    SampleIDs = sampleUID(UID_cycles, first_dilution_to)
    mutation_list = sample_mutation_list(UID_cycles, bases_per_amplicon, error_rate)

    count = 0
    for id in SampleIDs:
        for mut in mutation_list:
            if mut[0] < UID_cycles and within(id, influence_bond(mut[1],mut[0],UID_cycles)):
                count += 2
                break
            else:
                if mut[0] == UID_cycles and mut[1] == id:
                    count += 1
                    break
    
    print(count / (2.0 * first_dilution_to))
    
    




