
#include "tensor_impl.h"


//
// template instantiations of KroneckerProduct().
//
template
void KroneckerProduct<2, float, true>(unsigned M, const float **A, const float **in, float **out, unsigned int ndofs);
template
void KroneckerProduct<3, float, true>(unsigned M, const float **A, const float **in, float **out, unsigned int ndofs);
template
void KroneckerProduct<4, float, true>(unsigned M, const float **A, const float **in, float **out, unsigned int ndofs);

/// template
/// void KroneckerProduct<2, float, false>(unsigned M, const float **A, const float **in, float **out, unsigned int ndofs);
/// template
/// void KroneckerProduct<3, float, false>(unsigned M, const float **A, const float **in, float **out, unsigned int ndofs);
/// template
/// void KroneckerProduct<4, float, false>(unsigned M, const float **A, const float **in, float **out, unsigned int ndofs);

