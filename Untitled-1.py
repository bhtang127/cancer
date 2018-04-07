from libcpp.vector cimport vector
import numpy as np
cimport numpy as np
import pandas as pd
cimport cython

cdef extern from "mutation.cpp":
    vector[vector[unsigned int]] cpp_poisson_mutation(vector[int] nsamples_per_molecule, int cycles, int bases_per_amplicon, double error_rate, const char tag)
    vector[int] expand(vector[int] n_samples)

def expand(np.ndarray n_samples):
    cdef vector[int] result = expand()

