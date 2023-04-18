//
// Created by milinda on 2/8/16.
//
/**
 * @brief:
 * contains main dendro data types and some basic definitions
 * 
*/
#ifndef DENDRO_KT_DENDRO_H
#define DENDRO_KT_DENDRO_H

#include <climits>
#include <complex>

#define RED "\e[1;31m"
#define BLU "\e[2;34m"
#define GRN "\e[0;32m"
#define YLW "\e[0;33m"
#define MAG "\e[0;35m"
#define CYN "\e[0;36m"
#define NRM "\e[0m"



#ifdef USE_64BIT_INDICES
#define DendroIntL long long
#define DendroIntLSpecifier %lld
#define DendroUIntLSpecifier %llu
#else
#define DendroIntL unsigned int
#define DendroIntLSpecifier %d
#define DendroUIntLSpecifier %u
#endif


#define DendroScalar double
#define DendroComplex std::complex<double>

// Use DendroScopeBegin() and DendroScopeEnd() in main().
// - DendroScopeBegin() after init'ing MPI, before declaring Dendro objects.
// - DendroScopeEnd() before finalizing MPI to de-scope automatic variables.
#define DendroScopeBegin() {
#define DendroScopeEnd() }


//#define DendroIntL unsigned int
typedef unsigned __int128 DendroUInt_128;

#endif //DENDRO_KT_DENDRO_H
