/**
 * @file arraySlice.h
 * @author Masado Ishii
 * @created 2019-03-31
 * @description Axis-aligned slices from multi-d arrays.
 */

#ifndef DENDRO_KT_ARRAY_SLICE_H
#define DENDRO_KT_ARRAY_SLICE_H

#include "mathUtils.h"

// Number of consecutive most-significant set bits in sub-vector. (Compile-time expression)
template <typename T>
constexpr unsigned char c_numConsecutiveOnesTop(T arg, unsigned char strLength = 8*sizeof(T), unsigned char prefix = 0)
{
  return (!strLength || !(arg & (1u << (strLength-1)))) ? prefix : c_numConsecutiveOnesTop(arg, strLength-1, prefix + 1);
}


//
// PowerLoop: Iterate over the volume of a hypercube, calling kernel on lexicographic index.
//
template <unsigned int dim>
struct PowerLoop
{
  template <typename Kernel>
  static void loop(const unsigned int M, Kernel &kernel)
  {
    const unsigned int numIterations = intPow(M, dim);
    for (int ii = 0; ii < numIterations; ii++)
      kernel(ii);
  }
};
template <>
struct PowerLoop<0u>
{
  template <typename Kernel>
  static void loop(const unsigned int M, Kernel &kernel) { kernel(0); }
};


//
// OffsetKernel: Defining '+' on index passed to kernel.
//
template <typename Kernel>
struct OffsetKernel
{
  unsigned int m_offset;
  Kernel &m_kernel;

  OffsetKernel(unsigned int offset, Kernel &kernel) : m_offset(offset), m_kernel(kernel) {}
  OffsetKernel() = delete;

  void operator() (unsigned int ii) { m_kernel(m_offset + ii); }
};


//
// PowerLoopSlice: Create nested loops when some iteration-space axes are skipped.
//
template <unsigned int dim, unsigned int orientFlag>
struct PowerLoopSlice
{
  static constexpr unsigned char numOnes = c_numConsecutiveOnesTop<unsigned int>(orientFlag, dim);
  static constexpr unsigned char numZeros = c_numConsecutiveOnesTop<unsigned int>(~orientFlag, dim-numOnes);

  static constexpr unsigned int outerStridePow = dim - numOnes;

  static constexpr unsigned char innerDim = dim - numOnes - numZeros;
  static constexpr unsigned int innerOrientFlag = orientFlag & ((1u << innerDim) - 1);

  template <typename Kernel>
  static void loop(const unsigned int M, Kernel &kernel)
  {
    using InnerKernel = OffsetKernel<Kernel>;

    struct PowerKernel
    {
      unsigned int m_M;
      InnerKernel m_innerKernel;

      PowerKernel(unsigned int M, Kernel &kernel) : m_M(M), m_innerKernel(0, kernel) {}
      PowerKernel() = delete;

      void operator() (unsigned int ii)
      {
        m_innerKernel.m_offset = ii * intPow(m_M, outerStridePow);
        PowerLoopSlice<innerDim, innerOrientFlag>::template loop(m_M, m_innerKernel);
      }
    };

    PowerKernel powerKernel(M, kernel);
    PowerLoop<numOnes>::template loop(M, powerKernel);
  }
};
template <unsigned int orientFlag>
struct PowerLoopSlice<0u, orientFlag>
{
  template <typename Kernel>
  static void loop(const unsigned int M, Kernel &kernel) { kernel(0); }
};


#endif//DENDRO_KT_ARRAY_SLICE_H
