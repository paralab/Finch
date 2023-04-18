/**
 * @file:tsort.tcc
 * @author: Masado Ishii  --  UofU SoC,
 * @date: 2019-01-11
 */

namespace ot
{


//
// locTreeSort()
//
template<typename T, unsigned int D>
template <class PointType>
void
SFC_Tree<T,D>:: locTreeSort(PointType *points,
                          RankI begin, RankI end,
                          LevI sLev,
                          LevI eLev,
                          RotI pRot)
{
  //// Recursive Depth-first, similar to Most Significant Digit First. ////

  if (end <= begin) { return; }

  constexpr char numChildren = TreeNode<T,D>::numChildren;
  constexpr unsigned int rotOffset = 2*numChildren;  // num columns in rotations[].

  // Reorder the buckets on sLev (current level).
  std::array<RankI, numChildren+1> tempSplitters;
  RankI ancStart, ancEnd;
  /// SFC_bucketing(points, begin, end, sLev, pRot, tempSplitters, ancStart, ancEnd);
  SFC_bucketing_impl<KeyFunIdentity_Pt<PointType>, PointType, PointType>(
      points, begin, end, sLev, pRot,
      KeyFunIdentity_Pt<PointType>(), true, true,
      tempSplitters,
      ancStart, ancEnd);

  // The array `tempSplitters' has numChildren+1 slots, which includes the
  // beginning, middles, and end of the range of children.
  // Ancestor splitters are is ancStart and ancEnd, not tempSplitters.

  // Lookup tables to apply rotations.
  const ChildI * const rot_perm = &rotations[pRot*rotOffset + 0*numChildren];
  const RotI * const orientLookup = &HILBERT_TABLE[pRot*numChildren];

  if (sLev < eLev)  // This means eLev is further from the root level than sLev.
  {
    // Recurse.
    // Use the splitters to specify ranges for the next level of recursion.
    for (char child_sfc = 0; child_sfc < numChildren; child_sfc++)
    {
      // Columns of HILBERT_TABLE are indexed by the Morton rank.
      // According to Dendro4 TreeNode.tcc:199 they are.
      // (There are possibly inconsistencies in the old code...?
      // Don't worry, we can regenerate the table later.)
      ChildI child = rot_perm[child_sfc];
      RotI cRot = orientLookup[child];

      if (tempSplitters[child_sfc+1] - tempSplitters[child_sfc+0] <= 1)
        continue;

      if (sLev > 0)
      {
        locTreeSort(points,
            tempSplitters[child_sfc+0], tempSplitters[child_sfc+1],
            sLev+1, eLev,
            cRot);
      }
      else   // Special handling if we have to consider the domain boundary.
      {
        locTreeSort(points,
            tempSplitters[child_sfc+0], tempSplitters[child_sfc+1],
            sLev+1, eLev,
            pRot);
      }
    }
  }
}// end function()


//
// locTreeSort() (with parallel companion array)
//
template<typename T, unsigned int D>
template <class KeyFun, typename PointType, typename KeyType, typename Companion, bool useCompanions>
void
SFC_Tree<T,D>:: locTreeSort(PointType *points, Companion *companions,
                          RankI begin, RankI end,
                          LevI sLev,
                          LevI eLev,
                          RotI pRot,
                          KeyFun keyfun)
{
  //// Recursive Depth-first, similar to Most Significant Digit First. ////

  if (end <= begin) { return; }

  constexpr char numChildren = TreeNode<T,D>::numChildren;
  constexpr unsigned int rotOffset = 2*numChildren;  // num columns in rotations[].

  // Reorder the buckets on sLev (current level).
  std::array<RankI, numChildren+1> tempSplitters;
  RankI ancStart, ancEnd;
  SFC_bucketing_general<KeyFun, PointType, KeyType, Companion, useCompanions>(
      points, companions, begin, end, sLev, pRot,
      keyfun, true, true,
      tempSplitters,
      ancStart, ancEnd);

  // The array `tempSplitters' has numChildren+1 slots, which includes the
  // beginning, middles, and end of the range of children.
  // Ancestor splitters are is ancStart and ancEnd, not tempSplitters.

  // Lookup tables to apply rotations.
  const ChildI * const rot_perm = &rotations[pRot*rotOffset + 0*numChildren];
  const RotI * const orientLookup = &HILBERT_TABLE[pRot*numChildren];

  if (sLev < eLev)  // This means eLev is further from the root level than sLev.
  {
    // Recurse.
    // Use the splitters to specify ranges for the next level of recursion.
    for (char child_sfc = 0; child_sfc < numChildren; child_sfc++)
    {
      // Columns of HILBERT_TABLE are indexed by the Morton rank.
      // According to Dendro4 TreeNode.tcc:199 they are.
      // (There are possibly inconsistencies in the old code...?
      // Don't worry, we can regenerate the table later.)
      ChildI child = rot_perm[child_sfc];
      RotI cRot = orientLookup[child];

      if (tempSplitters[child_sfc+1] - tempSplitters[child_sfc+0] <= 1)
        continue;

      if (sLev > 0)
      {
        locTreeSort<KeyFun, PointType, KeyType, Companion, useCompanions>
            (points, companions,
            tempSplitters[child_sfc+0], tempSplitters[child_sfc+1],
            sLev+1, eLev,
            cRot,                         // This branch uses cRot.
            keyfun);
      }
      else   // Special handling if we have to consider the domain boundary.
      {
        locTreeSort<KeyFun, PointType, KeyType, Companion, useCompanions>
            (points, companions,
            tempSplitters[child_sfc+0], tempSplitters[child_sfc+1],
            sLev+1, eLev,
            pRot,                         // This branch uses pRot.
            keyfun);
      }
    }
  }
}// end function()


//
// SFC_bucketing_impl()
//
template <typename T, unsigned int D>
template <class KeyFun, typename PointType, typename KeyType>
void
SFC_Tree<T,D>:: SFC_bucketing_impl(PointType *points,
                          RankI begin, RankI end,
                          LevI lev,
                          RotI pRot,
                          KeyFun keyfun,
                          bool separateAncestors,
                          bool ancestorsFirst,
                          std::array<RankI, 1+TreeNode<T,D>::numChildren> &outSplitters,
                          RankI &outAncStart,
                          RankI &outAncEnd)
{
  // Call the "companion" implementation without giving or using companions.
  SFC_bucketing_general<KeyFun, PointType, KeyType, int, false>(
      points, nullptr,
      begin, end, lev, pRot,
      keyfun, separateAncestors, ancestorsFirst,
      outSplitters, outAncStart, outAncEnd);
}


//
// SFC_bucketing_general()
//
template <typename T, unsigned int D>
template <class KeyFun, typename PointType, typename KeyType, typename Companion, bool useCompanions>
void
SFC_Tree<T,D>:: SFC_bucketing_general(PointType *points, Companion* companions,
                          RankI begin, RankI end,
                          LevI lev,
                          RotI pRot,
                          KeyFun keyfun,
                          bool separateAncestors,
                          bool ancestorsFirst,
                          std::array<RankI, 1+TreeNode<T,D>::numChildren> &outSplitters,
                          RankI &outAncStart,
                          RankI &outAncEnd)
{
  using TreeNode = TreeNode<T,D>;
  constexpr char numChildren = TreeNode::numChildren;
  constexpr char rotOffset = 2*numChildren;  // num columns in rotations[].

  // -- Reconstruct offsets and bucketEnds from returned splitters. -- //
  SFC_locateBuckets_impl<KeyFun,PointType,KeyType>(
      points, begin, end, lev, pRot,
      keyfun, separateAncestors, ancestorsFirst,
      outSplitters, outAncStart, outAncEnd);

  std::array<RankI, numChildren+1> offsets, bucketEnds;  // Last idx represents ancestors.

  const ChildI *rot_perm = &rotations[pRot*rotOffset + 0*numChildren];
  for (ChildI child_sfc = 0; child_sfc < numChildren; child_sfc++)
  {
    ChildI child = rot_perm[child_sfc];
    offsets[child] = outSplitters[child_sfc];
    bucketEnds[child] = outSplitters[child_sfc+1];
  }
  offsets[numChildren] = outAncStart;
  bucketEnds[numChildren] = outAncEnd;

  // -- Movement phase. -- //
  std::array<PointType, numChildren+1> unsortedBuffer;
  std::array<Companion, numChildren+1> unsortedBufferComp;
  int bufferSize = 0;

  for (char bucketId = 0; bucketId <= numChildren; bucketId++)
  {
    if (offsets[bucketId] < bucketEnds[bucketId])
    {
      unsortedBuffer[bufferSize]     = points[offsets[bucketId]];  // Copy TreeNode.
      if (useCompanions)
        unsortedBufferComp[bufferSize] = companions[offsets[bucketId]];  // Copy companion.
      bufferSize++;
    }
  }

  while (bufferSize > 0)
  {
    PointType *bufferTop     = &unsortedBuffer[bufferSize-1];
    Companion *bufferTopComp;
    if (useCompanions)
      bufferTopComp = &unsortedBufferComp[bufferSize-1];
    unsigned char destBucket
      = (separateAncestors && keyfun(*bufferTop).getLevel() < lev) ? numChildren : keyfun(*bufferTop).getMortonIndex(lev);
    // destBucket is used to index into offsets[] and bucketEnds[], for which
    // ancestors are represented in [numChildren] regardless of `ancestorsFirst'.

    points[offsets[destBucket]]     = *bufferTop;  // Set down the TreeNode.
    if (useCompanions)
      companions[offsets[destBucket]] = *bufferTopComp;  // Set down the companion.
    offsets[destBucket]++;

    if (offsets[destBucket] < bucketEnds[destBucket])
    {
      *bufferTop     = points[offsets[destBucket]];    // Copy TreeNode.
      if (useCompanions)
        *bufferTopComp = companions[offsets[destBucket]];    // Copy companion.
    }
    else
      bufferSize--;
  }
}


//
// SFC_locateBuckets_impl()
//
template <typename T, unsigned int D>
template <class KeyFun, typename PointType, typename KeyType>
void
SFC_Tree<T,D>:: SFC_locateBuckets_impl(const PointType *points,
                          RankI begin, RankI end,
                          LevI lev,
                          RotI pRot,
                          KeyFun keyfun,
                          bool separateAncestors,
                          bool ancestorsFirst,
                          std::array<RankI, 1+TreeNode<T,D>::numChildren> &outSplitters,
                          RankI &outAncStart,
                          RankI &outAncEnd)
{
  using TreeNode = TreeNode<T,D>;
  constexpr char numChildren = TreeNode::numChildren;
  constexpr char rotOffset = 2*numChildren;  // num columns in rotations[].

  std::array<int, numChildren> counts;
  counts.fill(0);
  int countAncestors = 0;   // Special bucket to ensure ancestors bucketed properly.
  for (const PointType *pt = points + begin; pt < points + end; pt++)
  {
    const KeyType &tn = keyfun(*pt);
    if (separateAncestors && tn.getLevel() < lev)
      countAncestors++;
    else
      counts[tn.getMortonIndex(lev)]++;
  }

  /// std::array<RankI, numChildren+1> offsets, bucketEnds;  // Last idx represents ancestors.
  RankI accum;
  if (ancestorsFirst)
    accum = begin + countAncestors;                  // Ancestors belong in front.
  else
    accum = begin;                                   // Else first child is front.

  const ChildI *rot_perm = &rotations[pRot*rotOffset + 0*numChildren];
  ChildI child_sfc = 0;
  for ( ; child_sfc < numChildren; child_sfc++)
  {
    ChildI child = rot_perm[child_sfc];
    outSplitters[child_sfc] = accum;
    /// offsets[child] = accum;           // Start of bucket. Moving marker.
    accum += counts[child];
    /// bucketEnds[child] = accum;        // End of bucket. Fixed marker.
  }
  outSplitters[child_sfc] = accum;  // Should be the end of siblings..

  if (ancestorsFirst)
  {
    /// offsets[numChildren] = begin;
    /// bucketEnds[numChildren] = begin + countAncestors;
    outAncStart = begin;
    outAncEnd = begin + countAncestors;
  }
  else
  {
    /// offsets[numChildren] = accum;
    /// bucketEnds[numChildren] = accum + countAncestors;
    outAncStart = accum;
    outAncEnd = accum + countAncestors;
  }
}



//
// SFC_Tree::dist_bcastSplitters()
//
template <typename T, unsigned int dim>
std::vector<TreeNode<T,dim>> SFC_Tree<T,dim>::dist_bcastSplitters(const TreeNode<T,dim> *start, MPI_Comm comm)
{
  int nProc, rProc;
  MPI_Comm_rank(comm, &rProc);
  MPI_Comm_size(comm, &nProc);

  using TreeNode = TreeNode<T,dim>;
  std::vector<TreeNode> splitters(nProc);
  splitters[rProc] = *start;

  for (int turn = 0; turn < nProc; turn++)
    par::Mpi_Bcast<TreeNode>(&splitters[turn], 1, turn, comm);

  return splitters;
}


//
// SFC_Tree::dist_bcastSplitters() (when activeComm != globalComm)
//
template <typename T, unsigned int dim>
std::vector<TreeNode<T, dim>> SFC_Tree<T, dim>::dist_bcastSplitters(
    const TreeNode<T, dim> *start,
    MPI_Comm globalComm,
    MPI_Comm activeComm,
    bool isActive,
    std::vector<int> &activeList)
{
  int rProc, nProc;
  MPI_Comm_size(globalComm, &nProc);
  MPI_Comm_rank(globalComm, &rProc);

  int activeSize, activeRank;;
  if (isActive)
  {
    MPI_Comm_size(activeComm, &activeSize);
    MPI_Comm_rank(activeComm, &activeRank);
  }

  // Decide on a global root, who must also be an active rank.
  const int voteRoot = (isActive ? rProc : nProc);
  int globalRoot = -1, activeRoot = -1;
  par::Mpi_Allreduce(&voteRoot, &globalRoot, 1, MPI_MIN, globalComm);
  if (rProc == globalRoot)
    activeRoot = activeRank;
  par::Mpi_Bcast(&activeRoot, 1, globalRoot, globalComm); // For active ranks.
  par::Mpi_Bcast(&activeSize, 1, globalRoot, globalComm); // For inactive ranks.

  activeList.clear();
  activeList.resize(activeSize);
  std::vector<TreeNode<T, dim>> activeSplitters(activeSize);

  if (isActive)
  {
    // Collect from active ranks to active root.
    par::Mpi_Gather(&rProc, &activeList[0], 1, activeRoot, activeComm);
    par::Mpi_Gather(start, &activeSplitters[0], 1, activeRoot, activeComm);
  }
  // Share with everyone.
  par::Mpi_Bcast(&activeList[0], activeSize, globalRoot, globalComm);
  par::Mpi_Bcast(&activeSplitters[0], activeSize, globalRoot, globalComm);

  return activeSplitters;
}


} // end namespace ot
