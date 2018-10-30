#include <cstdlib>
#include <set>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>

void inplace_poisson_mutation(unsigned int* retval, int rows, int columns,
                              unsigned int* nsamples_per_molecule, int n_molecule, int* cycles, int l, 
                              int bases_per_amplicon, double error_rate);

void inplace_expand(int* indexes, int olen, int* n_samples, int ilen);

void inplace_sampling(int* samples, int len, int* cycles, int l, int total_samples);