%module fast_mutation

%{
  #define SWIG_FILE_WITH_INIT
  #include "fast_mutation.h"
%}

/* Include the NumPy typemaps library */
%include "numpy.i"

%init %{
  import_array();
%}

/* Typemaps */
%apply (int* IN_ARRAY1, int DIM1) {(int* n_samples, int ilen)};
%apply (int* INPLACE_ARRAY1, int DIM1) {(int* indexes, int olen)};
%apply (unsigned int* IN_ARRAY1, int DIM1) {(unsigned int* nsamples_per_molecule, int n_molecule)};
%apply (unsigned int* INPLACE_ARRAY2, int DIM1, int DIM2) {(unsigned int* retval, int rows, int columns)};

%include "fast_mutation.h"


