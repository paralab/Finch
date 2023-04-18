/**
 * @file:nsort.cpp
 * @author: Masado Ishii  --  UofU SoC,
 * @date: 2019-02-20
 */


#include "nsort.h"

namespace ot
{

  std::ostream & operator<<(std::ostream &out, const ScatterMap &sm)
  {
    out << "{\n";

    out << "m_map: [";
    for (RankI locId : sm.m_map)
      out << locId << ", ";
    out << "]\n";

    out << "m_sendCounts: [";
    for (RankI c : sm.m_sendCounts)
      out << c << ", ";
    out << "]\n";

    out << "m_sendOffsets: [";
    for (RankI off : sm.m_sendOffsets)
      out << off << ", ";
    out << "]\n";

    out << "m_sendProc: [";
    for (int rank : sm.m_sendProc)
      out << rank << ", ";
    out << "]\n";

    out << "}";

    return out;
  }


  std::ostream & operator<<(std::ostream &out, const GatherMap &gm)
  {
    out << "{\n";

    out << "m_recvProc: [";
    for (int rank : gm.m_recvProc)
      out << rank << ", ";
    out << "]\n";

    out << "m_recvCounts: [";
    for (RankI c : gm.m_recvCounts)
      out << c << ", ";
    out << "]\n";

    out << "m_recvOffsets: [";
    for (RankI off : gm.m_recvOffsets)
      out << off << ", ";
    out << "]\n";

    out << "m_totalCount: " << gm.m_totalCount << "; ";
    out << "m_locCount: " << gm.m_locCount << "; ";
    out << "m_locOffset: " << gm.m_locOffset << "; ";

    out << "}";

    return out;
  }

}

