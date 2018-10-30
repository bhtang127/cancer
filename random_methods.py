import random
import numpy as np
import time
import pandas as pd

import fast_mutation as fm

def expand(n_samples):
    s = np.sum(n_samples)
    n_samples = np.asarray(n_samples, dtype=np.int32)
    indexes = -np.ones(s, dtype=np.int32)
    fm.inplace_expand( indexes, n_samples )
    if np.min(indexes) < 0:
        raise(ValueError)

    return indexes

def nsamples_generator(molecules, d_rate = -1, d_to = -1, c_max, p, q):
    cycles_group = []
    samples = []
    cycles = []

    for n in molecules:
        cycles_group += np.random.binomial(c_max, p, n) * np.random.binomial(1, q, 1)
    
    if d_rate >= 0:
        nsample_per_well = [ int(np.sum( 2**x ) * d_rate) for x in cycles_group ]
        for i in range(len(molecules)):
            samples_i = np.zeros(molecules[i], dtype=np.int32)
            fm.inplace_sampling(samples_i, cycles_group[i], nsample_per_well[i])
            samples += list(samples_i)
            cycles += list(cycles_group[i])

    else:
        samples = np.zeros( np.sum(molecules), dtype=np.int32 )
        for i in range(len(molecules)):
            cycles += list(cycles_group[i])
        max_mol = 2 ** np.asarray(cycles)
        if np.sum(max_mol) <= d_to:
            samples = max_mol
        else:
            fm.inplace_sampling(samples, cycles, d_to)

    return samples, cycles

def poisson_mutation(nsamples_per_molecule, cycles, bases_per_amplicon, error_rate, tag):
    rows = np.sum(nsamples_per_molecule)
    columns = 5 + (tag != "UID")
    nsamples_per_molecule = np.asarray(nsamples_per_molecule, dtype=np.uint32)
    data = np.zeros( [rows, columns] , dtype=np.uint32)

    fm.inplace_poisson_mutation( data, nsamples_per_molecule, cycles, bases_per_amplicon, error_rate )

    if tag == "UID":
        data = pd.DataFrame( data, columns=["leftUID", "rightUID", "id0",
                                            "UID_left_mutation", "UID_right_mutation"] )
    else:
        data = pd.DataFrame( data, columns=["leftWBC", "rightWBC", "id1", "left_side", 
                                            "WBC_left_mutation", "WBC_right_mutation"] )
    
    return data
