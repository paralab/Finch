
#include "tensor_impl.h"


//
// template instantiations of KroneckerProduct().
//
template
void KroneckerProduct<2, unsigned, true>(unsigned M, const unsigned **A, const unsigned **in, unsigned **out, unsigned int ndofs);
template
void KroneckerProduct<3, unsigned, true>(unsigned M, const unsigned **A, const unsigned **in, unsigned **out, unsigned int ndofs);
template
void KroneckerProduct<4, unsigned, true>(unsigned M, const unsigned **A, const unsigned **in, unsigned **out, unsigned int ndofs);

/// template
/// void KroneckerProduct<2, unsigned, false>(unsigned M, const unsigned **A, const unsigned **in, unsigned **out, unsigned int ndofs);
/// template
/// void KroneckerProduct<3, unsigned, false>(unsigned M, const unsigned **A, const unsigned **in, unsigned **out, unsigned int ndofs);
/// template
/// void KroneckerProduct<4, unsigned, false>(unsigned M, const unsigned **A, const unsigned **in, unsigned **out, unsigned int ndofs);

