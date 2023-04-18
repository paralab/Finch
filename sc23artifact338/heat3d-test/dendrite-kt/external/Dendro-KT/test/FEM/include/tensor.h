//
// Created by milinda on 1/15/17.
//

/**
 *
 * @author Milinda Fernando
 * @author Masado Ishii
 * @breif contains the utilities for tensor kronecker products for interpolations.
 *
 * */

#ifndef SFCSORTBENCH_DENDROTENSOR_H
#define SFCSORTBENCH_DENDROTENSOR_H


#include "mathUtils.h"


// TODO make a new namespace



/**
 * @tparam dim Dimension of element, i.e. order of tensor.
 * @tparam da  Datatype of vectors.
 * @tparam forward If true, axes are evaluated in increasing order.
 * @param [in] M Size of tensor in 1D.
 * @param [in] A Array of pointers to interpolation matrices, ordered by axis.
 * @param [in] in Array of pointers to input buffers, ordered with source in position 0.
 * @param [in] out Array of pointers to output buffers, ordered with destination in position (dim-1).
 * @param ndofs Number of degrees of freedom in the vector, e.g. 3 for xyzxyz.
 */
template <unsigned int dim, typename da, bool forward>
void KroneckerProduct(unsigned M, const da **A, const da **in, da **out, unsigned int ndofs);


/**
 * @brief Forall (i,j,k), Q_ijk *= A * P_ijk, where P_ijk == W_i * W_j * W_k;
 * @tparam dim Order of the tensor.
 * @tparam T Component type.
 *
 * Generalizes the following nested for-loop to any nesting level:
 *   int idx = 0;
 *   for (int k = 0; k < l; k++)
 *   {
 *     const double a2 = A*W[k];
 *     for (int j = 0; j < l; j++)
 *     {
 *       const double a1 = a2*W[j];
 *       for (int i = 0; i < l; i++)
 *       {
 *         const double a0 = a1*W[i];
 *         Q[idx++] *= a0;
 *       }
 *     }
 *   }
 */
template <typename T, unsigned int dim>
struct SymmetricOuterProduct
{
  inline static void applyHadamardProduct(unsigned int length1d, T *Q, const T *W1d, const T premult = 1);
};

template <typename T, unsigned int dim>
inline void SymmetricOuterProduct<T,dim>::applyHadamardProduct(
    unsigned int length1d, T *Q, const T *W1d, const T premult)
{
  const unsigned int stride = intPow(length1d, dim-1);
  for (unsigned int ii = 0; ii < length1d; ii++)
    SymmetricOuterProduct<T, dim-1>::applyHadamardProduct(length1d, &Q[stride * ii], W1d, premult * W1d[ii]);
}

template <typename T>
struct SymmetricOuterProduct<T, 1>
{
  inline static void applyHadamardProduct(unsigned int length1d, T *Q, const T *W1d, const T premult = 1)
  {
    for (unsigned int ii = 0; ii < length1d; ii++)
      Q[ii] *= premult * W1d[ii];
  }
};




/** Apply the 1D interpolation for the input vector x and output the interpolated values in the vector Y.
 *
 *
 *
 * \param [in]  M  size of the vector
 * \param [in]  A  interpolation matrix
 * \param [in]  X  input data for the interpolation
 * \param [out] Y  interpolated values.
 *

 */
void DENDRO_TENSOR_AIIX_APPLY_ELEM (const int M, const double*  A, const double*  X, double*  Y);



/** Apply the 1D interpolation for the input vector x and output the interpolated values in the vector Y.
 *
 *
 *
 * \param [in]  M  size of the vector
 * \param [in]  A  interpolation matrix
 * \param [in]  X  input data for the interpolation
 * \param [out] Y  interpolated values.
 *

 */
void DENDRO_TENSOR_IIAX_APPLY_ELEM(const int M, const double*  A, const double*  X, double*  Y);


/** Apply the 1D interpolation for the input vector x and output the interpolated values in the vector Y.
 *
 *
 *
 * \param [in]  M  size of the vector
 * \param [in]  A  interpolation matrix
 * \param [in]  X  input data for the interpolation
 * \param [out] Y  interpolated values.
 *

 */
void DENDRO_TENSOR_IAIX_APPLY_ELEM (const int M, const double*  A, const double*  X, double*  Y);




/** Apply the 1D interpolation for the input vector x and output the interpolated values in the vector Y.
 *
 *
 *
 * \param [in]  M  size of the vector
 * \param [in]  A  interpolation matrix
 * \param [in]  X  input data for the interpolation
 * \param [out] Y  interpolated values.
 *

 */
void DENDRO_TENSOR_IAX_APPLY_ELEM_2D(const int M, const double*  A, const double*  X, double*  Y);



/** Apply the 1D interpolation for the input vector x and output the interpolated values in the vector Y.
 *
 *
 *
 * \param [in]  M  size of the vector
 * \param [in]  A  interpolation matrix
 * \param [in]  X  input data for the interpolation
 * \param [out] Y  interpolated values.
 *

 */
void DENDRO_TENSOR_AIX_APPLY_ELEM_2D (const int M, const double*  A, const double*  X, double*  Y);


#endif //SFCSORTBENCH_DENDROTENSOR_H
