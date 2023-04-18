#include "IMGA/MathOperation.h"

namespace MathOp {

/**
 * Return the transpose of a 3 * 3 matrix (represented by 3 ZEROPTVs)
 * @param m_in
 * @return
 */
std::vector<ZEROPTV> matTranspose(const std::vector<ZEROPTV> &m_in) {
  assert(m_in.size() == 3);
  std::vector<ZEROPTV> m_out;
  m_out.resize(3);
  for (int dim1 = 0; dim1 < 3; dim1++) {
    for (int dim2 = 0; dim2 < 3; dim2++) {
      m_out[dim1](dim2) = m_in[dim2](dim1);
    }
  }
  return m_out;
}

/**
 * return inverse of a 3 * 3 matrix
 * @param m_in
 * @return m_out
 */
std::vector<ZEROPTV> inverse_matrix(const std::vector<ZEROPTV> &m_in) {
  assert(m_in.size() == 3);
  std::vector<ZEROPTV> m_out;
  m_out.resize(3);

  double det = m_in[0](0) * (m_in[1](1) * m_in[2](2) - m_in[2](1) * m_in[1](2))
      - m_in[0](1) * (m_in[1](0) * m_in[2](2) - m_in[1](2) * m_in[2](0))
      + m_in[0](2) * (m_in[1](0) * m_in[2](1) - m_in[1](1) * m_in[2](0));

  double invdet = 1.0 / det;

  m_out[0](0) = (m_in[1](1) * m_in[2](2) - m_in[2](1) * m_in[1](2)) * invdet;
  m_out[0](1) = (m_in[0](2) * m_in[2](1) - m_in[0](1) * m_in[2](2)) * invdet;
  m_out[0](2) = (m_in[0](1) * m_in[1](2) - m_in[0](2) * m_in[1](1)) * invdet;
  m_out[1](0) = (m_in[1](2) * m_in[2](0) - m_in[1](0) * m_in[2](2)) * invdet;
  m_out[1](1) = (m_in[0](0) * m_in[2](2) - m_in[0](2) * m_in[2](0)) * invdet;
  m_out[1](2) = (m_in[1](0) * m_in[0](2) - m_in[0](0) * m_in[1](2)) * invdet;
  m_out[2](0) = (m_in[1](0) * m_in[2](1) - m_in[2](0) * m_in[1](1)) * invdet;
  m_out[2](1) = (m_in[2](0) * m_in[0](1) - m_in[0](0) * m_in[2](1)) * invdet;
  m_out[2](2) = (m_in[0](0) * m_in[1](1) - m_in[1](0) * m_in[0](1)) * invdet;
  return m_out;
}

/**
 * matrix matrix multiplication
 * @param m_in_1
 * @param m_in_2
 * @return
 */
std::vector<ZEROPTV> mat_mat_multi(const std::vector<ZEROPTV> &m_in_1, const std::vector<ZEROPTV> &m_in_2) {
  assert(m_in_1.size() == 3 && m_in_2.size() == 3);

  std::vector<ZEROPTV> m_out;
  m_out.resize(3);

  for (int dim1 = 0; dim1 < 3; dim1++) {
    for (int dim2 = 0; dim2 < 3; dim2++) {
      for (int k = 0; k < 3; k++) {
        m_out[dim1](dim2) += m_in_1[dim1](k) * m_in_2[k](dim2);
      }
    }
  }
  return m_out;
}

/**
* return matrix vector multiplication for 3 by 3 only
* @param m_in
* @param v_in
* @return v_out
*/
ZEROPTV mat_vec_multi(const std::vector<ZEROPTV> &m_in, const ZEROPTV &v_in) {
  assert(m_in.size() == 3);
  ZEROPTV vec_out;
  for (int dim = 0; dim < 3; dim++) {
    for (int k = 0; k < 3; k++) {
      vec_out(dim) += m_in[dim](k) * v_in(k);
    }
  }
  return vec_out;
}

/**
 * return the star matrix of a vector
 * @param vec_in
 * @return
 */
std::vector<ZEROPTV> Star(const ZEROPTV &vec_in) {
  std::vector<ZEROPTV> m_out;
  m_out.resize(3);
  m_out[0](0) = 0.0;
  m_out[1](1) = 0.0;
  m_out[2](2) = 0.0;
  m_out[0](1) = -vec_in(2);
  m_out[1](0) = vec_in(2);
  m_out[0](2) = vec_in(1);
  m_out[2](0) = -vec_in(1);
  m_out[1](2) = -vec_in(0);
  m_out[2](1) = vec_in(0);
  return m_out;
}

}

