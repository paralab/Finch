////
//// Created by maksbh on 7/12/20.
////
//
//#ifndef DENDRITEKT_IBMSOLVER_H
//#define DENDRITEKT_IBMSOLVER_H
//#include <IMGA/Geometry/Geometry.h>
//#include <NodeAndValues.h>
//#include <point.h>
//#include <OctToPhysical.h>
//class IBMSolver{
//  /**
//   * @brief fill the points vector with the Gauss points. This function does not require any
//   * communication but does repetitive work
//   * @param points points
//   */
//  void fillGaussPointsWithoutCommunication(std::vector<Point<DIM>> & points);
//  /**
//   * @brief fill the points vector with the Gauss points. This function require communication
//   * @param points
//   */
//  void fillGaussPointsWithCommunication(std::vector<Point<DIM>> & points);
//
// protected:
//    const STL::Geometry * m_geometry; /// geometry
//    std::vector<NodeAndValues<DENDRITE_REAL,DIM>> gaussPoints_; /// vector with Gauss points
//    const DA * m_octDA; /// DA
//    const DomainInfo m_domain; /// domain
//
//    /**
//     * @brief Find the correspondong tree nodes
//     * @param domainInfo the domain boundary information
//     * @param points Gauss points (basically any points, for us its the Gauss points)
//     * @param [out] treeNodes return the tree nodes
//     */
//    static void computeTreeNodesAndCordinates(const DomainInfo & domainInfo, const std::vector<Point<DIM>> & points,std::vector<TREENODE>& treeNodes) ;
//
// public:
//  /**
//   * @brief Constructor: Once called this function automatically fills the Gauss points and partition them
//   * @param octDA DA
//   * @param geometry geometry
//   * @param domainInfo domainInfo
//   */
//    IBMSolver(const DA * octDA, const STL::Geometry * geometry, const DomainInfo & domainInfo);
//
//    /**
//     * @brief fill the Gauss point vector and partition them
//     */
//    void fillAndPartitionGaussPoints();
//
//    /**
//     * @brief Finds the background elements
//     */
//    void getBackGroundElement();
//
//    /**
//     *
//     * @return return the vector with Gauss points
//     */
//    const std::vector<NodeAndValues<DENDRITE_REAL,DIM>> & getGaussPoints() const;
//
//
//
//
//};
//
//#endif //DENDRITEKT_IBMSOLVER_H
