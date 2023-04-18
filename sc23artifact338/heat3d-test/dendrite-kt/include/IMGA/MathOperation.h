#pragma once

#include <algorithm>
#include <vector>
#include <talyfem/grid/zeroptv.h>
#include <point.h>
#include <talyfem/common/exceptions.h>

namespace MathOp {
using TALYFEMLIB::ZEROPTV;
using namespace std;
using namespace TALYFEMLIB;
/**
 * transpose of a matrix (3*3 only)
 * @param m_in
 * @return
 */
std::vector<ZEROPTV> matTranspose(const std::vector<ZEROPTV> &m_in);

/**
 * inverse of a matrix (3*3 only)
 * @param m_in
 * @return m_out
 */
std::vector<ZEROPTV> inverse_matrix(const std::vector<ZEROPTV> &m_in);

/**
 * matrix matrix multiplication (3*3 only)
 * @param m_in_1
 * @param m_in_2
 * @return
 */
std::vector<ZEROPTV> mat_mat_multi(const std::vector<ZEROPTV> &m_in_1, const std::vector<ZEROPTV> &m_in_2);

/**
* matrix vector multiplication (3*3 only)
* @param m_in
* @param v_in
* @return v_out
*/
ZEROPTV mat_vec_multi(const std::vector<ZEROPTV> &m_in, const ZEROPTV &v_in);

/**
 * the star matrix of a vector (3*3 only)
 * @param vec_in
 * @return
 */
std::vector<ZEROPTV> Star(const ZEROPTV &vec_in);


}
