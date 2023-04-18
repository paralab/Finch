
#ifndef DENDRO_KT_TENSOR_IMPL_H
#define DENDRO_KT_TENSOR_IMPL_H

#include "tensor.h"
#include <arraySlice.h>


/**
 * @tparam dim Dimension of element, i.e. order of tensor.
 * @tparam da  Datatype of vectors.
 * @tparam forward If true, axes are evaluated in increasing order.
 * @tparam ndofs Number of degrees of freedom in the vector, e.g. 3 for xyzxyz.
 * @param [in] M Size of tensor in 1D.
 * @param [in] A Array of pointers to interpolation matrices, ordered by axis.
 * @param [in] in Array of pointers to input buffers, ordered with source in position 0.
 * @param [in] out Array of pointers to output buffers, ordered with destination in position (dim-1).
 */
template <unsigned int dim, typename da, bool forward, unsigned int ndofs>
void KroneckerProductFixedDof(unsigned M, const da **A, const da **in, da **out);



template <typename da, unsigned int ndofs>
struct MatKernelAssign
{
  const da *X;
  da *Y;
  void operator()(da d, unsigned int x_ii, unsigned int y_ii) const
  {
    for (int dof = 0; dof < ndofs; dof++)
      Y[ndofs*y_ii + dof] = d * X[ndofs*x_ii + dof];
  }
};

template <typename da, unsigned int ndofs>
struct MatKernelAccum
{
  const da *X;
  da *Y;
  void operator()(da d, unsigned int x_ii, unsigned int y_ii) const
  {
    for (int dof = 0; dof < ndofs; dof++)
      Y[ndofs*y_ii + dof] += d * X[ndofs*x_ii + dof];
  }
};

//TODO there is a way to incorporate ndofs as a new axis of the tensor
//  instead of modifying MatKernelAssign and MatKernelAccum.
//  I have chosen the latter for now for simplicity.


// Local matrix multiplication over 1D of a general slice of tensor.
template <unsigned int dim, unsigned int sliceFlag, unsigned int tangent, unsigned int ndofs>
struct IterateBindMatrix;

// Specialization of IterateBindMatrix to the full tensor.
template <unsigned int dim, unsigned int tangent, unsigned int ndofs>
using IterateTensorBindMatrix = IterateBindMatrix<dim, (1u<<dim)-1, tangent, ndofs>;

// Specialization of IterateBindMatrix to a (K-1)-face of the tensor.
template <unsigned int dim, unsigned int face, unsigned int tangent, unsigned int ndofs>
using IterateFacetBindMatrix = IterateBindMatrix<dim, ((1u<<dim)-1) - (1u<<face), tangent, ndofs>;


//
// IterateBindMatrix
//
// Usage: IterateBindMatrix<dim, sliceFlag, tangent, ndofs>::template iterate_bind_matrix<da>(M, A, X, Y);
//
template <unsigned int dim, unsigned int sliceFlag, unsigned int tangent, unsigned int ndofs>
struct IterateBindMatrix
{
  static constexpr unsigned int upperOrientFlag = sliceFlag & (- (1u<<(tangent+1)));
  static constexpr unsigned int lowerOrientFlag = sliceFlag & ((1u<<tangent) - 1);

  using OuterLoop = PowerLoopSlice<dim, upperOrientFlag>;
  using InnerLoop = PowerLoopSlice<dim, lowerOrientFlag>;

  template <typename da>
  struct InnerKernel
  {
    MatKernelAssign<da, ndofs> &m_kernel1;
    MatKernelAccum<da, ndofs> &m_kernel2;
    struct SAssign { InnerKernel &p; void operator()(unsigned int innerIdx) { p.actuate(p.m_kernel1, innerIdx); } } Assign;
    struct SAccum { InnerKernel &p; void operator()(unsigned int innerIdx) { p.actuate(p.m_kernel2, innerIdx); } } Accum;
    unsigned int m_M;
    InnerKernel(MatKernelAssign<da, ndofs> &kernel1, MatKernelAccum<da, ndofs> &kernel2, unsigned int M)
        : m_kernel1(kernel1), m_kernel2(kernel2), Assign{*this}, Accum{*this}, m_M(M) {}
    InnerKernel() = delete;

    da m_d;
    unsigned int m_row;
    unsigned int m_col;
    unsigned int m_outerIdx;
    template <typename Kernel>
    inline void actuate(Kernel &kernel, unsigned int innerIdx);
  };

  template <typename da>
  struct OuterKernel
  {
    InnerKernel<da> &m_kernel;
    unsigned int m_M;
    const da *m_APtr;
    OuterKernel(InnerKernel<da> &kernel, unsigned int M, const da *A) : m_kernel(kernel), m_M(M), m_APtr(A) {}
    OuterKernel() = delete;

    unsigned int m_row;
    inline void operator()(unsigned int outerIdx);
  };

  template <typename da>
  inline static void iterate_bind_matrix(const unsigned int M, const da *A, const da *X, da *Y);
};


  //
  // iterate_bind_matrix()
  //
  template <unsigned int dim, unsigned int sliceFlag, unsigned int tangent, unsigned int ndofs>
  template <typename da>
  inline void IterateBindMatrix<dim,sliceFlag,tangent,ndofs>::iterate_bind_matrix(const unsigned int M, const da *A, const da *X, da *Y)
  {
    // For each row of the matrix, iterate through the hyperplane
    // of the index space that is normal to the axis `face'.
    // `dim' specifies the dimension of the index space.
    // `tangent' specifies which axis should be bound to the matrix row.

    MatKernelAssign<da, ndofs> matKernelAssign{X, Y};
    MatKernelAccum<da, ndofs> matKernelAccum{X, Y};

    InnerKernel<da> inner{matKernelAssign, matKernelAccum, M};
    OuterKernel<da> outer{inner, M, A};

    unsigned int &row = outer.m_row;
    for (row = 0; row < M; row++)
      OuterLoop::template loop(M, outer);
  }
  //
  // OuterKernel::operator()()
  //
  template <unsigned int dim, unsigned int sliceFlag, unsigned int tangent, unsigned int ndofs>
  template <typename da>
  inline void IterateBindMatrix<dim,sliceFlag,tangent,ndofs>::
       OuterKernel<da>::
       operator()(unsigned int outerIdx)
  {
    // Iterate down tangent axis (input and matrix).
    m_kernel.m_outerIdx = outerIdx;
    m_kernel.m_row = m_row;
    unsigned int &col = m_kernel.m_col;
    col = 0;
    m_kernel.m_d = m_APtr[m_row];
    InnerLoop::template loop(m_M, m_kernel.Assign);
    for (col = 1; col < m_M; col++)
    {
      m_kernel.m_d = m_APtr[col * m_M + m_row];
      InnerLoop::template loop(m_M, m_kernel.Accum);
    }
  }
  //
  // InnerKernel::actuate()
  //
  template <unsigned int dim, unsigned int sliceFlag, unsigned int tangent, unsigned int ndofs>
  template <typename da>
  template <typename Kernel>
  inline void IterateBindMatrix<dim,sliceFlag,tangent,ndofs>::
       InnerKernel<da>::
       actuate(Kernel &kernel, unsigned int innerIdx)
  {
    // Combine indices into x_ii and y_ii;
    const unsigned int tangentStride = intPow(m_M, tangent);
    unsigned int x_ii = m_outerIdx + m_col * tangentStride + innerIdx;
    unsigned int y_ii = m_outerIdx + m_row * tangentStride + innerIdx;
    kernel(m_d, x_ii, y_ii);
  }



template <unsigned int dim, unsigned int d, bool forward, unsigned int ndofs>
struct KroneckerProduct_loop { template <typename da> static void body(unsigned M, const da **A, const da **in, da **out) {
  constexpr unsigned int ii = (forward ? dim-1 - d : d);

  if (forward)
    IterateTensorBindMatrix<dim, d, ndofs>::template iterate_bind_matrix<da>(M, A[d], in[ii], out[ii]);

  KroneckerProduct_loop<dim, d-1, forward, ndofs>::template body<da>(M,A,in,out);

  if (!forward)
    IterateTensorBindMatrix<dim, d, ndofs>::template iterate_bind_matrix<da>(M, A[d], in[ii], out[ii]);
}};
template <unsigned int dim, bool forward, unsigned int ndofs>
struct KroneckerProduct_loop<dim, 0, forward, ndofs> { template <typename da> static void body(unsigned M, const da **A, const da **in, da **out) {
  constexpr unsigned int ii = (forward ? dim-1 : 0);
  IterateTensorBindMatrix<dim, 0, ndofs>::template iterate_bind_matrix<da>(M, A[0], in[ii], out[ii]);
}};

template <unsigned int dim, typename da, bool forward, unsigned int ndofs>
void KroneckerProductFixedDof(unsigned M, const da **A, const da **in, da **out)
{
  KroneckerProduct_loop<dim, dim-1, forward, ndofs>::template body<da>(M,A,in,out);
}

template <unsigned int dim, typename da, bool forward>
void KroneckerProduct(unsigned M, const da **A, const da **in, da **out, unsigned int ndofs)
{
  // Convert runtime argument to template argument.
  switch (ndofs)
  {
    // TODO add CMake options and #if macros for greater number of dofs.
    case 1:
      KroneckerProductFixedDof<dim, da, forward, 1>(M, A, in, out);
      break;
    case 2:
      KroneckerProductFixedDof<dim, da, forward, 2>(M, A, in, out);
      break;
    case 3:
      KroneckerProductFixedDof<dim, da, forward, 3>(M, A, in, out);
      break;
    case 4:
      KroneckerProductFixedDof<dim, da, forward, 4>(M, A, in, out);
      break;
    case 5:
      KroneckerProductFixedDof<dim, da, forward, 5>(M, A, in, out);
      break;
    case 6:
      KroneckerProductFixedDof<dim, da, forward, 6>(M, A, in, out);
      break;
    case 7:
      KroneckerProductFixedDof<dim, da, forward, 7>(M, A, in, out);
      break;
    case 8:
      KroneckerProductFixedDof<dim, da, forward, 8>(M, A, in, out);
      break;
    case 9:
      KroneckerProductFixedDof<dim, da, forward, 9>(M, A, in, out);
      break;
    case 10:
      KroneckerProductFixedDof<dim, da, forward, 10>(M, A, in, out);
      break;
    case 11:
      KroneckerProductFixedDof<dim, da, forward, 11>(M, A, in, out);
      break;
    case 12:
      KroneckerProductFixedDof<dim, da, forward, 12>(M, A, in, out);
      break;
    case 13:
      KroneckerProductFixedDof<dim, da, forward, 13>(M, A, in, out);
      break;
    case 14:
      KroneckerProductFixedDof<dim, da, forward, 14>(M, A, in, out);
      break;
    case 15:
      KroneckerProductFixedDof<dim, da, forward, 15>(M, A, in, out);
      break;
    case 16:
      KroneckerProductFixedDof<dim, da, forward, 16>(M, A, in, out);
      break;
    case 17:
      KroneckerProductFixedDof<dim, da, forward, 17>(M, A, in, out);
      break;
    case 18:
      KroneckerProductFixedDof<dim, da, forward, 18>(M, A, in, out);
      break;
    case 19:
      KroneckerProductFixedDof<dim, da, forward, 19>(M, A, in, out);
      break;
    case 20:
      KroneckerProductFixedDof<dim, da, forward, 20>(M, A, in, out);
      break;
    case 21:
      KroneckerProductFixedDof<dim, da, forward, 21>(M, A, in, out);
      break;
    case 22:
      KroneckerProductFixedDof<dim, da, forward, 22>(M, A, in, out);
      break;
    case 23:
      KroneckerProductFixedDof<dim, da, forward, 23>(M, A, in, out);
      break;
    case 24:
      KroneckerProductFixedDof<dim, da, forward, 24>(M, A, in, out);
      break;
    case 25:
      KroneckerProductFixedDof<dim, da, forward, 25>(M, A, in, out);
      break;
    case 26:
      KroneckerProductFixedDof<dim, da, forward, 26>(M, A, in, out);
      break;
    case 27:
      KroneckerProductFixedDof<dim, da, forward, 27>(M, A, in, out);
      break;
    case 28:
      KroneckerProductFixedDof<dim, da, forward, 28>(M, A, in, out);
      break;
    case 29:
      KroneckerProductFixedDof<dim, da, forward, 29>(M, A, in, out);
      break;
    case 30:
      KroneckerProductFixedDof<dim, da, forward, 30>(M, A, in, out);
      break;
    case 31:
      KroneckerProductFixedDof<dim, da, forward, 31>(M, A, in, out);
      break;
    case 32:
      KroneckerProductFixedDof<dim, da, forward, 32>(M, A, in, out);
      break;
    case 33:
      KroneckerProductFixedDof<dim, da, forward, 33>(M, A, in, out);
      break;
    case 34:
      KroneckerProductFixedDof<dim, da, forward, 34>(M, A, in, out);
      break;
    case 35:
      KroneckerProductFixedDof<dim, da, forward, 35>(M, A, in, out);
      break;
    case 36:
      KroneckerProductFixedDof<dim, da, forward, 36>(M, A, in, out);
      break;
    case 37:
      KroneckerProductFixedDof<dim, da, forward, 37>(M, A, in, out);
      break;
    case 38:
      KroneckerProductFixedDof<dim, da, forward, 38>(M, A, in, out);
      break;
    case 39:
      KroneckerProductFixedDof<dim, da, forward, 39>(M, A, in, out);
      break;
    case 40:
      KroneckerProductFixedDof<dim, da, forward, 40>(M, A, in, out);
      break;
    case 41:
      KroneckerProductFixedDof<dim, da, forward, 41>(M, A, in, out);
      break;
    case 42:
      KroneckerProductFixedDof<dim, da, forward, 42>(M, A, in, out);
      break;
    case 43:
      KroneckerProductFixedDof<dim, da, forward, 43>(M, A, in, out);
      break;
    case 44:
      KroneckerProductFixedDof<dim, da, forward, 44>(M, A, in, out);
      break;
    case 45:
      KroneckerProductFixedDof<dim, da, forward, 45>(M, A, in, out);
      break;
    case 46:
      KroneckerProductFixedDof<dim, da, forward, 46>(M, A, in, out);
      break;
    case 47:
      KroneckerProductFixedDof<dim, da, forward, 47>(M, A, in, out);
      break;
    case 48:
      KroneckerProductFixedDof<dim, da, forward, 48>(M, A, in, out);
      break;
    case 49:
      KroneckerProductFixedDof<dim, da, forward, 49>(M, A, in, out);
      break;
    case 50:
      KroneckerProductFixedDof<dim, da, forward, 50>(M, A, in, out);
      break;
    default:
      const bool isNumberOfDegreesOfFreedomSupported = false;
      assert(isNumberOfDegreesOfFreedomSupported);  // Need to add more cases.
      break;
  }
}













#endif//DENDRO_KT_TENSOR_IMPL_H
