/**
 * @file distTree.h
 * @author Masado Ishii, UofU SoC
 * @date 2019-10-04
 * @brief Struct to hold part of a distributed tree.
 */

#ifndef DENDRO_KT_DIST_TREE_H
#define DENDRO_KT_DIST_TREE_H


#include "treeNode.h"
#include "octUtils.h"
#include "mpi.h"

#include "filterFunction.h"

namespace ot
{
  /**
   * @brief Intermediate container for filtering trees before creating the DA.
   *
   * @note  DistTree takes ownership of the provided TreeNodes and empties the
   *        provided std::vector.
   *
   * @note  It is intended that, during construction of the DA, the tree
   *        vector held by DistTree will be destroyed.
   *
   * @note  Create a DistTree from a partitioned complete tree, i.e. taking
   *        union of the TreeNodes across all processors should be the entire
   *        unit hypercube. If you want to filter the domain to a subset of the
   *        unit hypercube, use DistTree to accomplish that.
   *
   * @note  A default filter function is provided and applied upon construction.
   *        To provide a custom filter function, call filterTree().
   *
   * @note  DistTree remembers the front and back TreeNode from the original partition.
   * @note  The partition cannot be changed without creating a new DistTree.
   */
  template <typename T, unsigned int dim>
  class DistTree
  {
    public:
      // Member functions.
      //
      DistTree();
      DistTree(std::vector<TreeNode<T, dim>> &treePart, MPI_Comm comm);
      DistTree(const DistTree &other) { this->operator=(other); }
      DistTree & operator=(const DistTree &other);

      static DistTree constructSubdomainDistTree(
          unsigned int finestLevel,
          MPI_Comm comm,
          double sfc_tol = 0.3);

      static DistTree constructSubdomainDistTree(
          unsigned int finestLevel,
          const ::ibm::DomainDecider &domainDecider_phys,
          MPI_Comm comm,
          double sfc_tol = 0.3);

      template <typename D>
      static DistTree constructDistTreeByFunc(
          std::function<void(const D *, D *)> func,
          unsigned int dofSz,
          MPI_Comm comm,
          unsigned int order,
          double interp_tol,
          double sfc_tol);


      /** distRemeshSubdomain
       *  @brief: Uses DistTree filter function to carve out subdomain from remeshed tree.
       */
      static void distRemeshSubdomain( const DistTree &inTree,
                                       const std::vector<OCT_FLAGS::Refine> &refnFlags,
                                       DistTree &outTree,
                                       DistTree &surrogateTree,
                                       double loadFlexibility );

      // generateGridHierarchyUp()
      //   TODO should return surrogate DistTree
      void generateGridHierarchyUp(bool isFixedNumStrata,
                                 unsigned int lev,
                                 double loadFlexibility);

      // generateGridHierarchyDown()
      //
      //   Replaces internal single grid with internal list of grids.
      //   Returns a DistTree with the surrogate grid at same level as each
      //     coarse grid (the finest grid has no surrogate).
      DistTree generateGridHierarchyDown(unsigned int numStrata, double loadFlexibility);

      // filterTree() requires the input decider to be in physical coordinate form.
      // Discards 'IN' elements, keeps 'OUT' and 'INTERCEPTED' elements,
      //   and calls setIsOnTreeBdry(true) on 'INTERCEPTED' elements.
      void filterTree(const ::ibm::DomainDecider &domainDecider);

      void destroyTree();

      const ::ibm::DomainDecider & getDomainDecider() const;

      const std::function<::ibm::Partition(const TreeNode<T, dim> &treeNodeElem)>
        & getDomainDeciderTN_asCell() const;

      const std::function<::ibm::Partition(const TreeNode<T, dim> &treeNodeElem)>
        & getDomainDeciderTN_asPoint() const;

      // Tree accessors.
      const std::vector<TreeNode<T, dim>> & getTreePartFiltered(int stratum = 0) const;
      size_t getOriginalTreePartSz(int stratum = 0) const;
      size_t getFilteredTreePartSz(int stratum = 0) const;
      TreeNode<T, dim> getTreePartFront(int stratum = 0) const;
      TreeNode<T, dim> getTreePartBack(int stratum = 0) const;


      int getNumStrata() { return m_numStrata; }


      // These deciders can be called directly.

      // Default domainDecider (physical)
      static ::ibm::Partition defaultDomainDecider(const double * physCoords, double physSize)
      {
        // Default: keep the unit cube.

        // In = discard, Out = retain

        // In_1d:    (coord + sz) <= 0.0  or  coord >= 1.0  // Closed "In" contains closure of octant.
        // Out_1d:   coord  > 0.0  and (coord + sz) < 1.0   // Open "Out" contains closure of octant.

        // In this case of the box,
        //   In_Nd:  OR_{over dims} { In_1d() }
        //   Out_Nd: AND_{over dims} { Out_1d() }

        bool isIn = false;
        bool isOut = true;

        for (int d = 0; d < dim; ++d)
        {
          isIn  |= (physCoords[d] + physSize <= 0.0 or physCoords[d] >= 1.0);
          isOut &= (physCoords[d] > 0.0 and physCoords[d] + physSize < 1.0);
        }

        // Can be isIn, or isOut, or neither, but NEVER both.

        if (isIn && !isOut)
          return ::ibm::IN;
        else if (isOut && !isIn)
          return ::ibm::OUT;
        else if (!isIn && !isOut)
          return ::ibm::INTERCEPTED;
        else
          throw std::logic_error("Filter function returned both isIn and isOut.");
      }


      // BoxDecider
      struct BoxDecider
      {
        BoxDecider(const std::array<double, dim> &physDims)
          : m_maxs(physDims)
        {
          m_mins.fill(0.0);
        }

        BoxDecider(const std::array<double, dim> &physMins,
                   const std::array<double, dim> &physMaxs)
          : m_mins(physMins), m_maxs(physMaxs)
        {}

        std::array<double, dim> m_maxs;
        std::array<double, dim> m_mins;


        ::ibm::Partition operator()(const double * physCoords, double physSize) const
        {
          bool isIn = false;
          bool isOut = true;
  
          for (int d = 0; d < dim; ++d)
          {
            isIn  |= (physCoords[d] + physSize <= m_mins[d] or physCoords[d] >= m_maxs[d]);
            isOut &= (physCoords[d] > m_mins[d] and physCoords[d] + physSize < m_maxs[d]);
          }
  
          // Can be isIn, or isOut, or neither, but NEVER both.
  
          if (isIn && !isOut)
            return ::ibm::IN;
          else if (isOut && !isIn)
            return ::ibm::OUT;
          else if (!isIn && !isOut)
            return ::ibm::INTERCEPTED;
          else
            throw std::logic_error("Filter function returned both isIn and isOut.");
        }
      };


    protected:
      // Member variables.
      //
      MPI_Comm m_comm;
      ::ibm::DomainDecider m_domainDecider;
      std::function<::ibm::Partition(const TreeNode<T, dim> &treeNodeElem)> m_domainDeciderTN_asCell;
      std::function<::ibm::Partition(const TreeNode<T, dim> &treeNodeElem)> m_domainDeciderTN_asPoint;

      bool m_hasBeenFiltered;

      // Multilevel grids, finest first. Must initialize with at least one level.
      std::vector<std::vector<TreeNode<T, dim>>> m_gridStrata;
      std::vector<TreeNode<T, dim>> m_tpFrontStrata;
      std::vector<TreeNode<T, dim>> m_tpBackStrata;
      std::vector<size_t> m_originalTreePartSz;
      std::vector<size_t> m_filteredTreePartSz;


      int m_numStrata;

      // Protected accessors return a reference to the 0th stratum.
      std::vector<TreeNode<T, dim>> &get_m_treePartFiltered() { return m_gridStrata[0]; }
      const std::vector<TreeNode<T, dim>> &get_m_treePartFiltered() const { return m_gridStrata[0]; }
      TreeNode<T, dim> &get_m_treePartFront() { return m_tpFrontStrata[0]; }
      const TreeNode<T, dim> &get_m_treePartFront() const { return m_tpFrontStrata[0]; }
      TreeNode<T, dim> &get_m_treePartBack() { return m_tpBackStrata[0]; }
      const TreeNode<T, dim> &get_m_treePartBack() const { return m_tpBackStrata[0]; }


      void assignDomainDecider(const ::ibm::DomainDecider &domainDecider);


      //
      // Intrinsic Deciders (not callable directly).
      //

      // If given a decider on phys coords, can still test treeNodes.
      ::ibm::Partition conversionDomainDeciderTN_asCell(const TreeNode<T, dim> &tn)
      {
        double physCoords[dim];
        double physSize;
        treeNode2Physical(tn, physCoords, physSize);

        return m_domainDecider(physCoords, physSize);
      }

      ::ibm::Partition conversionDomainDeciderTN_asPoint(const TreeNode<T, dim> &tn)
      {
        double physCoords[dim];
        double physSize;
        treeNode2Physical(tn, physCoords, physSize);

        return m_domainDecider(physCoords, 0.0);
      }

  };


  //
  // DistTree() - default constructor
  //
  template <typename T, unsigned int dim>
  DistTree<T, dim>::DistTree()
  : m_comm(MPI_COMM_NULL),
    m_gridStrata(m_uiMaxDepth+1),
    m_tpFrontStrata(m_uiMaxDepth+1),
    m_tpBackStrata(m_uiMaxDepth+1),
    m_originalTreePartSz(m_uiMaxDepth+1, 0),
    m_filteredTreePartSz(m_uiMaxDepth+1, 0),
    m_numStrata(0)
  {
    this->assignDomainDecider(DistTree::defaultDomainDecider);

    m_hasBeenFiltered = false;
  }


  //
  // DistTree() - constructor
  //
  template <typename T, unsigned int dim>
  DistTree<T, dim>::DistTree(std::vector<TreeNode<T, dim>> &treePart, MPI_Comm comm)
  : m_comm(comm),
    m_gridStrata(m_uiMaxDepth+1),
    m_tpFrontStrata(m_uiMaxDepth+1),
    m_tpBackStrata(m_uiMaxDepth+1),
    m_originalTreePartSz(m_uiMaxDepth+1, 0),
    m_filteredTreePartSz(m_uiMaxDepth+1, 0),
    m_numStrata(1)
  {
    if (comm != MPI_COMM_NULL)
      SFC_Tree<T, dim>::distCoalesceSiblings(treePart, comm);

    m_originalTreePartSz[0] = treePart.size();
    m_filteredTreePartSz[0] = treePart.size();

    if (treePart.size())
    {
      get_m_treePartFront() = treePart.front();
      get_m_treePartBack() = treePart.back();
    }

    get_m_treePartFiltered().clear();
    std::swap(get_m_treePartFiltered(), treePart);  // Steal the tree vector.

    this->filterTree(DistTree::defaultDomainDecider);
  }



  template <typename C, unsigned int dim>
  template <typename T>
  DistTree<C, dim> DistTree<C, dim>::constructDistTreeByFunc(
      std::function<void(const T *, T *)> func,
      unsigned int dofSz,
      MPI_Comm comm,
      unsigned int order,
      double interp_tol,
      double sfc_tol)
  {
    std::vector<unsigned int> varIndex(dofSz);
    for (unsigned int ii = 0; ii < dofSz; ii++)
      varIndex[ii] = ii;

    // Get a complete tree sufficiently granular to represent func with accuracy interp_tol.
    std::vector<ot::TreeNode<C,dim>> completeTree;
    function2Octree<C,dim>(func, dofSz, &(*varIndex.cbegin()), dofSz, completeTree, m_uiMaxDepth, interp_tol, sfc_tol, order, comm);

    // Make the tree balanced, using completeTree as a minimal set of TreeNodes.
    // Calling distTreeBalancing() on a complete tree with ptsPerElement==1
    // should do exactly what we want.
    std::vector<ot::TreeNode<C,dim>> balancedTree;
    ot::SFC_Tree<C,dim>::distTreeBalancing(completeTree, balancedTree, 1, sfc_tol, comm);

    ot::DistTree<C,dim> distTree(balancedTree, comm);   // Uses default domain decider.

    return distTree;
  }





  //
  // operator=()
  //
  template <typename T, unsigned int dim>
  DistTree<T, dim> &  DistTree<T, dim>::operator=(const DistTree &other)
  {
    m_comm =                  other.m_comm;
    m_hasBeenFiltered =       other.m_hasBeenFiltered;
    m_gridStrata =            other.m_gridStrata;
    m_tpFrontStrata =         other.m_tpFrontStrata;
    m_tpBackStrata =          other.m_tpBackStrata;
    m_originalTreePartSz =    other.m_originalTreePartSz;
    m_filteredTreePartSz =    other.m_filteredTreePartSz;
    m_numStrata =             other.m_numStrata;

    this->assignDomainDecider(other.m_domainDecider);

    return *this;
  }


  //
  // destroyTree()
  //
  template <typename T, unsigned int dim>
  void DistTree<T, dim>::destroyTree()
  {
    for (std::vector<TreeNode<T, dim>> &gridStratum : m_gridStrata)
    {
      gridStratum.clear();
      gridStratum.shrink_to_fit();
    }
  }



  //
  // assignDomainDecider()
  //
  template <typename T, unsigned int dim>
  void DistTree<T, dim>::assignDomainDecider(const ibm::DomainDecider &domainDecider)
  {
    m_domainDecider = domainDecider;
    {
      using namespace std::placeholders;
      m_domainDeciderTN_asCell  = std::bind(&DistTree<T,dim>::conversionDomainDeciderTN_asCell, this, _1);
      m_domainDeciderTN_asPoint = std::bind(&DistTree<T,dim>::conversionDomainDeciderTN_asPoint, this, _1);
    }
  }


  //
  // filterTree()
  //
  template <typename T, unsigned int dim>
  void DistTree<T, dim>::filterTree( const ::ibm::DomainDecider &domainDecider)
  {
    this->assignDomainDecider(domainDecider);

    for (int l = 0; l < m_numStrata; ++l)
    {
      const size_t oldSz = m_gridStrata[l].size();
      size_t ii = 0;

      m_filteredTreePartSz[l] = ii;

      ::ibm::Partition subdomain;

      // Keep finding and deleting elements.
      for ( ; ii < oldSz ; ii++)
        if ( (subdomain = this->m_domainDeciderTN_asCell(m_gridStrata[l][ii])) != ::ibm::IN )
        {
          m_gridStrata[l][m_filteredTreePartSz[l]] = std::move(m_gridStrata[l][ii]);
          m_gridStrata[l][m_filteredTreePartSz[l]].setIsOnTreeBdry(subdomain == ::ibm::INTERCEPTED);
          m_filteredTreePartSz[l]++;
        }

      m_gridStrata[l].resize(m_filteredTreePartSz[l]);
    }

    m_hasBeenFiltered = true;
  }



  //
  // getDomainDecider()
  //
  template <typename T, unsigned int dim>
  const ::ibm::DomainDecider & DistTree<T, dim>::getDomainDecider() const
  {
    return m_domainDecider;
  }

  //
  // getDomainDeciderTN_asCell()
  //
  template <typename T, unsigned int dim>
  const std::function<::ibm::Partition(const TreeNode<T, dim> &treeNodeElem)> &
      DistTree<T, dim>::getDomainDeciderTN_asCell() const
  {
    return m_domainDeciderTN_asCell;
  }

  //
  // getDomainDeciderTN_asPoint()
  //
  template <typename T, unsigned int dim>
  const std::function<::ibm::Partition(const TreeNode<T, dim> &treeNodePoint)> &
      DistTree<T, dim>::getDomainDeciderTN_asPoint() const
  {
    return m_domainDeciderTN_asPoint;
  }


  //
  // getTreePartFiltered()
  //
  template <typename T, unsigned int dim>
  const std::vector<TreeNode<T, dim>> &
      DistTree<T, dim>::getTreePartFiltered(int stratum) const
  {
    return m_gridStrata[stratum];
  }


  //
  // getOriginalTreePartSz()
  //
  template <typename T, unsigned int dim>
  size_t DistTree<T, dim>::getOriginalTreePartSz(int stratum) const
  {
    return m_originalTreePartSz[stratum];
  }


  //
  // getFilteredTreePartSz()
  //
  template <typename T, unsigned int dim>
  size_t DistTree<T, dim>::getFilteredTreePartSz(int stratum) const
  {
    return m_filteredTreePartSz[stratum];
  }


  //
  // getTreePartFront()
  //
  template <typename T, unsigned int dim>
  TreeNode<T, dim> DistTree<T, dim>::getTreePartFront(int stratum) const
  {
    return m_tpFrontStrata[stratum];
  }


  //
  // getTreePartBack()
  //
  template <typename T, unsigned int dim>
  TreeNode<T, dim> DistTree<T, dim>::getTreePartBack(int stratum) const
  {
    return m_tpBackStrata[stratum];
  }



}//namespace ot



#endif//DENDRO_KT_DIST_TREE_H
