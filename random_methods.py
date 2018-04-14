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
    samples = np.zeros(counts, dtype=np.int32)
    fm.inplace_sampling(samples, total_samples, restriction)
    return samples

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
