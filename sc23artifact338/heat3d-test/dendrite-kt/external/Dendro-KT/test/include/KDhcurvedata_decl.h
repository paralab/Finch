/**
 * @author Masado Ishii
 * @date   2019-01-18
 */

#ifndef DENDRO_KT_HCURVEDATA_DECL_H
#define DENDRO_KT_HCURVEDATA_DECL_H

#include "hcurvedata.h"

#include <cstring>

template <int pDim>
struct HilbertData
{
  static void copyData(char * &rotations, int * &HILBERT_TABLE);
  static const char m_rotations[_KD_ROTATIONS_SIZE(pDim)];
  static const int m_HILBERT_TABLE[_KD_HILBERT_TABLE_SIZE(pDim)];
};

template <int pDim>
void HilbertData<pDim>::copyData(char * &rotations, int * &HILBERT_TABLE)
{
  rotations = new char[sizeof(m_rotations)];
  memcpy(rotations, m_rotations, sizeof(m_rotations));

  HILBERT_TABLE = new int[sizeof(m_HILBERT_TABLE) / sizeof(int)];
  memcpy(HILBERT_TABLE, m_HILBERT_TABLE, sizeof(m_HILBERT_TABLE));
}

#endif // DENDRO_KT_HCURVEDATA_DECL_H
