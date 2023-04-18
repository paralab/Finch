//
// Created by milinda on 3/31/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains useful math helper routines
*/
//

#ifndef SFCSORTBENCH_MATHUTILS_H
#define SFCSORTBENCH_MATHUTILS_H

#include <iostream>
#include <cmath>
#include <math.h>
#include <array>
#include "parUtils.h"

/**
 * @brief computes the l2 norm between two vectors.
 * @param[in] vec1 input vector 1
 * @param[in] vec2 input vector 2
 * @param[in] dimension of the vector.
 * @param[in] comm communicator
 * @return l2 norm between two vectors.
 * */
template <typename T>
T normL2(T * vec1,T* vec2, unsigned int n,MPI_Comm comm);

/**
 * @brief computes the l2 norm between two vectors.
 * @param[in] vec1 input vector 1
 * @param[in] vec2 input vector 2
 * @param[in] dimension of the vector.
 * @param[in] comm communicator
 * @return l1 norm between two vectors.
 * */
template <typename T>
T normLInfty(T *vec1, T *vec2, unsigned int n, MPI_Comm comm);


/**
 * @brief computes the l2 norm between two vectors.
 * @param[in] vec1 input vector 1
 * @param[in] vec2 input vector 2
 * @param[in] dimension of the vector.
 * @return l1 norm between two vectors.
 * */
template <typename T>
T normLInfty(T *vec1, T *vec2, unsigned int n);



/**
 * @brief computes the l2 norm of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @param[in] comm communicator
 * @return l2 norm of vec.
 * */
template <typename T>
T normL2(T * vec,unsigned int n, MPI_Comm comm);


/**
 * @brief computes the l2 norm between two vectors.
 * @param[in] vec1 input vector 1
 * @param[in] vec2 input vector 2
 * @param[in] dimension of the vector.
 * @return l2 norm between two vectors.
 * */
template <typename T>
T normL2(T * vec1,T* vec2, unsigned int n);

/**
 * @brief computes the l2 norm of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @return l2 norm of vec.
 * */
template <typename T>
T normL2(T * vec,unsigned int n);


/**
 * @brief computes the l_inf norm of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @return l2 norm of vec.
 * */
template <typename T>
T normLInfty(T * vec,unsigned int n);


/**
 * @brief computes the l_inf norm of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @param[in] comm communicator
 * @return l2 norm of vec.
 * */
template <typename T>
T normLInfty(T * vec,unsigned int n,MPI_Comm comm);

/**
 * @brief computes the min of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @return min of vec.
 * */
template <typename T>
T vecMin(T * vec,unsigned int n);


/**
 * @brief computes the min of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @return max of vec.
 * */
template <typename T>
T vecMax(T * vec,unsigned int n);


/**
 * @brief computes the min of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @return min of vec.
 * */
template <typename T>
T vecMin(T * vec,unsigned int n,MPI_Comm comm);


/**
 * @brief computes the min of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @return max of vec.
 * */
template <typename T>
T vecMax(T * vec,unsigned int n,MPI_Comm comm);

/**@brief : performs the dot product of any two given vectors
 * @param[in] v1: input vector 1
 * @param[in] v2: input vector 2
 * @param[in] n: size of the vector ( dof of vector)
 * @return  v1^tv2
 * */
template <typename T>
T dot(const T* v1, const T*v2,const unsigned int n);

/**@brief : performs the dot product of any two given distributed vectors
 * Assumption : v1 and v2 are partitioned in the sameway.
 * @param[in] v1: input vector 1
 * @param[in] v2: input vector 2
 * @param[in] n: size of the vector ( dof of vector)
 * @return  v1^tv2
 * */
template <typename T>
T dot(const T* v1, const T*v2,const unsigned int n,MPI_Comm comm);


/**
 * @brief: Scaler multiplication of a vector.
 * @param[in] alpha: scalar value
 * @param[in] v: input vector.
 * @param[in] n: size of the vector ( dof of vector)
 * @param[in] comm : MPI communicator.
 *
 * */
template<typename T>
void mul(const T alpha, const T* v, const unsigned int n, T* out);


/**
 * @brief : add two vectors.
 * @param[in] v1: input vector 1
 * @param[in] v2: input vector 2
 * @param[in] n: size of the vector ( dof of vector)
 * @param[out] out: out=v1+v2
 * */
template <typename T>
T add(const T* v1, const T*v2,const unsigned int n, T* out);

/**
* @brief : substract two vectors.
* @param[in] v1: input vector 1
* @param[in] v2: input vector 2
* @param[in] n: size of the vector ( dof of vector)
* @param[out] out: out=v1+v2
                          * */
template <typename T>
T subt(const T* v1, const T*v2,const unsigned int n, T* out);

/**
 @brief Compile-time computation of integer (A*pow(b,p)).
 @author Masado Ishii
 */
template <typename T>
constexpr T intPow(T b, unsigned p, T A = 1);

/**
 @brief Compile-time computation of factorial "f!".
 @author Masado Ishii
 */
template <typename T>
constexpr T intFactorial(T f, T A = 1);

/**
 @brief Compile-time computation of combination "n choose k".
 @note Computes numerator and denominator separately; these could get large, so n and k should be small or T should be large.
 @author Masado Ishii
 */
template <typename T>
constexpr T intCombination(T n, T k, T p = 1, T q = 1);

/**
 @brief Prefix sum of binomial coefficients "n choose k'" up to k'==k-1. (Exclusive prefix sum).
 @author Masado Ishii
 @note Assumes that 0 <= k <= n.
 */
template <typename T>
T intCombinationSum(T n, T k, T &combo);


/**
 @brief Increments an array of numbers as if they are digits in a base B number.
 @description The reason for writing this is to traverse multidimensional lattices
     without using nested for-loops; for example, when the dimension is a template parameter.
 @author Masado Ishii
 */
template <typename T, unsigned int L>
void incrementBaseB(std::array<T,L> &digits, T B, unsigned int lstart = 0);

/**
 * @brief Increment a multi-index in a nested for loop.
 * @note This version allows different limits in each axis.
 */
template <typename T, unsigned int L>
void incrementFor(std::array<T,L> &digits,
                  const std::array<T,L> &limits,
                  unsigned int d = 0);



#include "mathUtils.tcc"


#endif //SFCSORTBENCH_MATHUTILS_H
