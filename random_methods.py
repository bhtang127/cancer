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

def nsamples_generator(counts, total_samples, restriction=0):
    if restriction is 0:
        pval = np.array(counts) / np.sum(counts)
        return np.random.multinomial(total_samples, pval)
    else:
        pval = np.array(counts) / np.sum(counts)
        l = len(counts)
        for a_try in range(20):
            samples = np.random.multinomial(total_samples, pval)
            if np.sum(samples > restriction) == 0:
                return samples
        n_samples = np.zeros(len(counts), dtype=np.int16)
        for i in range(total_samples):
            loc = np.random.choice(l, p=pval)
            while n_samples[loc] >= restriction[loc]:
                loc = np.random.choice(l, p=pval)
            n_samples[loc] += 1
        return n_samples

def poisson_mutation(nsamples_per_molecule, cycles, bases_per_amplicon, error_rate, tag):
    rows = np.sum(nsamples_per_molecule)
    columns = 4 + (tag != "UID")
    nsamples_per_molecule = np.asarray(nsamples_per_molecule, dtype=np.uint32)
    data = np.zeros( [rows, columns] , dtype=np.uint32)

    fm.inplace_poisson_mutation( data, nsamples_per_molecule, cycles, bases_per_amplicon, error_rate )

    if tag == "UID":
        data = pd.DataFrame( data, columns=["leftUID", "rightUID", 
                                            "UID_left_mutation", "UID_right_mutation"] )
    else:
        data = pd.DataFrame( data, columns=["leftWBC", "rightWBC", "left_side", 
                                            "WBC_left_mutation", "WBC_right_mutation"] )
    
    return data
