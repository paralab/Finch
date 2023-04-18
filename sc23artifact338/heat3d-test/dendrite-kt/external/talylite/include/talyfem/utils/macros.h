/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#pragma once

// Used in MPI communications that transfer PetscInts.
#ifdef PETSC_USE_64BIT_INDICES
#define MPI_TALYFEM_INT MPI_INT64_T
#else
#define MPI_TALYFEM_INT MPI_INT
#endif

// PETSCINT_F macro - use this when printfing PetscInts. Example:
// printf("%" PETSCINT_F " rows", n_rows);

// WHY?
// Depending on the compiler, PetscInt is typedefed to either
// 'int64_t' (guaranteed by the C standard to be exactly 64 bits),
// or 'long long int' (guaranteed to be at least 64 bits).
// (See petsc/include/petscsys.h for the exact typdef, around line 220.)

// HOWEVER, 'int64_t' is typedefed to 'long int' on some compilers.
// 'long int' is only guaranteed to be at least 32 bits, but on some compilers
// it happens to be 64 bits.
// The problem is that 'long int' and 'long long int' have different
// printf format specifiers (%ld vs %lld), even though they might be
// the same size.
// So, this macro exists so we can easily change the format specifier
// for PetscInts between 'long long int' and 'long int' as necessary.
#ifdef PETSC_USE_64BIT_INDICES
#define PETSCINT_F "ld"
#else
#define PETSCINT_F "d"
#endif

// Used to detect of a type has a member function for template metaprogramming.
// http://stackoverflow.com/a/23133787
#define DEFINE_HAS_SIGNATURE(traitsName, funcName, signature)               \
    template <typename U>                                                   \
    class traitsName                                                        \
    {                                                                       \
     private:                                                               \
        template<typename T, T> struct helper;                              \
        template<typename T>                                                \
        static std::uint8_t check(helper<signature, &funcName>*);           \
        template<typename T> static std::uint16_t check(...);               \
     public:                                                                \
        static                                                              \
        constexpr bool value = sizeof(check<U>(0)) == sizeof(std::uint8_t); \
    }

// https://en.wikibooks.org/wiki/More_C++_Idioms/Member_Detector
#define DEFINE_HAS_MEMBER(traitsName, X)                                            \
template<typename T>                                                                \
class traitsName {                                                                  \
    struct Fallback { int X; };                                                     \
    struct Derived : T, Fallback { };                                               \
                                                                                    \
    template<typename U, U> struct Check;                                           \
                                                                                    \
    typedef char ArrayOfOne[1];                                                     \
    typedef char ArrayOfTwo[2];                                                     \
                                                                                    \
    template<typename U> static ArrayOfOne & func(Check<int Fallback::*, &U::X> *); \
    template<typename U> static ArrayOfTwo & func(...);                             \
  public:                                                                           \
    typedef traitsName type;                                                        \
    enum { value = sizeof(func<Derived>(0)) == 2 };                                 \
};
