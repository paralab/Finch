
#include "tensor_impl.h"


//
// template instantiations of KroneckerProduct().
//
template
void KroneckerProduct<2, char, true>(unsigned M, const char **A, const char **in, char **out, unsigned int ndofs);
template
void KroneckerProduct<3, char, true>(unsigned M, const char **A, const char **in, char **out, unsigned int ndofs);
template
void KroneckerProduct<4, char, true>(unsigned M, const char **A, const char **in, char **out, unsigned int ndofs);

/// template
/// void KroneckerProduct<2, char, false>(unsigned M, const char **A, const char **in, char **out, unsigned int ndofs);
/// template
/// void KroneckerProduct<3, char, false>(unsigned M, const char **A, const char **in, char **out, unsigned int ndofs);
/// template
/// void KroneckerProduct<4, char, false>(unsigned M, const char **A, const char **in, char **out, unsigned int ndofs);
