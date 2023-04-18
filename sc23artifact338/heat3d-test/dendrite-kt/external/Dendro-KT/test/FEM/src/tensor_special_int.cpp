
#include "tensor_impl.h"


//
// template instantiations of KroneckerProduct().
//
template
void KroneckerProduct<2, int, true>(unsigned M, const int **A, const int **in, int **out, unsigned int ndofs);
template
void KroneckerProduct<3, int, true>(unsigned M, const int **A, const int **in, int **out, unsigned int ndofs);
template
void KroneckerProduct<4, int, true>(unsigned M, const int **A, const int **in, int **out, unsigned int ndofs);

/// template
/// void KroneckerProduct<2, int, false>(unsigned M, const int **A, const int **in, int **out, unsigned int ndofs);
/// template
/// void KroneckerProduct<3, int, false>(unsigned M, const int **A, const int **in, int **out, unsigned int ndofs);
/// template
/// void KroneckerProduct<4, int, false>(unsigned M, const int **A, const int **in, int **out, unsigned int ndofs);

