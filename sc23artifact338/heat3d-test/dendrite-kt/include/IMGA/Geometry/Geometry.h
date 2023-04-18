////
//// Created by maksbh on 7/12/20.
////
//
//#ifndef DENDRITEKT_GEOMETRY_H
//#define DENDRITEKT_GEOMETRY_H
//
//#include <DataTypes.h>
//#include <IMGA/Geometry/GeomType.h>
//
//#if (DIM == 3)
//namespace STL{
//
//  /**
//   * Every triangle is identified uniquely by three coordinates and one normal
//   * In STL files the connectivity is fixed. So, no additional information is required.
//   */
//  struct Triangles{
//    DENDRITE_REAL triangleCoord[3][3];
//    DENDRITE_REAL normal[3];
//  };
//
//  /**
//   * @brief class for Geometry class
//   */
//  class Geometry{
//   protected:
//    std::vector<Triangles> m_triangles; /// List of triangles and associated normals
//    const DENDRITE_UINT m_id; /// id of the geometry, if multiple exists
//    const GeomType m_geomType;
//   public:
//    /**
//     * @brief Constructor. The function loads the STL
//     * @param filePath filepath
//     * @param GeomType type of geometry : useful for analytical surface
//     * @param id id of the STL
//     */
//    Geometry(const std::string filePath,const GeomType & geomType, const DENDRITE_UINT id = 0);
//
//    /**
//     * Load STL
//     * @param filePath filepath
//     */
//    void loadSTL(const std::string filePath);
//
//    /**
//     * @brief clean up the geomwtry if not needed.
//     */
//    void cleanupGeometry();
//  };
//}
//#endif
//#endif //DENDRITEKT_GEOMETRY_H
