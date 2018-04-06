import random
import numpy as np
import time
import pandas as pd
from scipy.stats import poisson

def expand(n_samples):
    indexes = []
    for i, ns in enumerate(n_samples):
        if ns > 0:
            indexes += list( i * np.ones(ns,dtype=np.int16) )
    return indexes

def sampleUID(cycles, n_samples):
    ## cycles < 32
    if n_samples == 0:
        return None
    
    if cycles < 20:
        ids = np.array(np.random.choice(2**cycles, n_samples, replace=False), dtype=np.int32)
        return np.sort(ids)

    while True:
        ids = np.array(np.random.choice(2**cycles, n_samples), dtype=np.int32)
        if len(set(ids)) == n_samples:
            return np.sort(ids)

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

def poisson_mutation(nsamples_per_molecule, cycles, bases_per_amplicon, error_rate, tag):
    [sampleUID(cycles, ns) for ns in nsamples_per_molecule]