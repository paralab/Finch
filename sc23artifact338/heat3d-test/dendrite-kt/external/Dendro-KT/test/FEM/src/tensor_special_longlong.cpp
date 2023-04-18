
#include "tensor_impl.h"


//
// template instantiations of KroneckerProduct().
//
template
void KroneckerProduct<2, long long, true>(unsigned M, const long long **A, const long long **in, long long **out, unsigned int ndofs);
template
void KroneckerProduct<3, long long, true>(unsigned M, const long long **A, const long long **in, long long **out, unsigned int ndofs);
template
void KroneckerProduct<4, long long, true>(unsigned M, const long long **A, const long long **in, long long **out, unsigned int ndofs);

/// template
/// void KroneckerProduct<2, long long, false>(unsigned M, const long long **A, const long long **in, long long **out, unsigned int ndofs);
/// template
/// void KroneckerProduct<3, long long, false>(unsigned M, const long long **A, const long long **in, long long **out, unsigned int ndofs);
/// template
/// void KroneckerProduct<4, long long, false>(unsigned M, const long long **A, const long long **in, long long **out, unsigned int ndofs);


