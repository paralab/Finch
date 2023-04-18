
#include "tensor_impl.h"


//
// template instantiations of KroneckerProduct().
//
template
void KroneckerProduct<2, double, true>(unsigned M, const double **A, const double **in, double **out, unsigned int ndofs);
template
void KroneckerProduct<3, double, true>(unsigned M, const double **A, const double **in, double **out, unsigned int ndofs);
template
void KroneckerProduct<4, double, true>(unsigned M, const double **A, const double **in, double **out, unsigned int ndofs);

/// template
/// void KroneckerProduct<2, double, false>(unsigned M, const double **A, const double **in, double **out, unsigned int ndofs);
/// template
/// void KroneckerProduct<3, double, false>(unsigned M, const double **A, const double **in, double **out, unsigned int ndofs);
/// template
/// void KroneckerProduct<4, double, false>(unsigned M, const double **A, const double **in, double **out, unsigned int ndofs);
