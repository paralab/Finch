#include <talyfem/grid/elem_types/elem4dtesseract.h>

#include <talyfem/grid/grid_types/grid.h>
#include <talyfem/grid/node.h>


namespace TALYFEMLIB {

const int* ELEM4dTesseract::GetSurfaceCheckArray() const {
  /// NODE ORDERING OF A 3-D LINEAR SURFACE ELEMENT (Hilbert ordering)
  /// when the node order of the volume element is in HILBERT-ORDER
  /*static const int B16n4DCheckArray[] = {
      // {surface_id, node_id1, node_id2, node_id3, node_id4 }
      -1, 0, 3, 7, 4, 8, 11, 15, 12,
      +1, 1, 2, 6, 5, 9, 10, 14, 13,
      -2, 0, 4, 12, 8, 1, 5, 13, 9,
      +2, 3, 7, 15, 11, 2, 6, 14, 10,
      -3, 0, 8, 9, 1, 3, 11, 10, 2,
      +3, 4, 12, 13, 5, 7, 15, 14, 6,
      -4, 0, 1, 2, 3, 4, 5, 6, 7,
      +4, 8, 9, 10, 11, 12, 13, 14, 15,
  };*/

  /// NODE ORDERING OF A 3-D LINEAR SURFACE ELEMENT (Hilbert ordering)
  /// when the node order of the volume element is in Z-ORDER
  static const int B16n4DCheckArray[] = {
      // {surface_id, node_id1, node_id2, node_id3, node_id4 }
      -1, 0, 2, 6, 4, 8, 10, 14, 12,
      +1, 1, 3, 7, 5, 9, 11, 15, 13,
      -2, 0, 4, 12, 8, 1, 5, 13, 9,
      +2, 2, 6, 14, 10, 3, 7, 15, 11,
      -3, 0, 8, 9, 1, 2, 10, 11, 3,
      +3, 4, 12, 13, 5, 6, 14, 15, 7,
      -4, 0, 1, 3, 2, 4, 5, 7, 6,
      +4, 8, 9, 11, 10, 12, 13, 15, 14,
  };

  static const int B81n4DCheckArray[] = {
      // {surface_id, node_id1, node_id2, node_id3, node_id4 }
      -1, 0, 6, 24, 18, 54, 60, 78, 72, 27, 33, 51, 45, 3, 15, 21, 9, 12, 57, 69, 75, 63, 66, 30, 42, 48, 36, 39,
      +1, 2, 8, 26, 20, 56, 62, 80, 74, 29, 35, 53, 47, 5, 17, 23, 11, 14, 59, 71, 77, 65, 68, 32, 44, 50, 38, 41,
      -2, 0, 18, 72, 54, 2, 20, 74, 56, 1, 19, 73, 55, 9, 45, 63, 27, 36, 11, 47, 65, 29, 38, 10, 46, 64, 28, 37,
      +2, 6, 24, 78, 60, 8, 26, 80, 62, 7, 25, 79, 61, 15, 51, 69, 33, 42, 17, 53, 71, 35, 44, 16, 52, 70, 34, 43,
      -3, 0, 54, 56, 2, 6, 60, 62, 8, 3, 57, 59, 5, 27, 55, 29, 1, 28, 33, 61, 35, 7, 34, 30, 58, 32, 4, 31,
      +3, 18, 72, 74, 20, 24, 78, 80, 26, 21, 75, 77, 23, 45, 73, 47, 19, 46, 51, 79, 53, 25, 52, 48, 76, 50, 22, 49,
      -4, 0, 2, 8, 6, 18, 20, 26, 24, 9, 11, 17, 15, 1, 5, 7, 3, 4, 19, 23, 25, 21, 22, 10, 14, 16, 12, 13,
      +4, 54, 56, 62, 60, 72, 74, 80, 78, 63, 65, 71, 69, 55, 59, 61, 57, 58, 73, 77, 79, 75, 76, 64, 68, 70, 66, 67,
  };

  static const int B256n4DCheckArray[] = {
      // {surface_id, node_id1, node_id2, node_id3, node_id4 }
      -1, 0, 12, 60, 48, 192, 204, 252, 240, 64, 76, 124, 112, 128, 140, 188, 176, 4, 8, 28, 44, 56, 52, 32, 16, 20, 24,
      40, 36, 196, 200, 220, 236, 248, 244, 224, 208, 212, 216, 232, 228, 68, 72, 92, 108, 120, 116, 96, 80, 84, 88,
      104, 100, 132, 136, 156, 172, 184, 180, 160, 144, 148, 152, 168, 164,

      +1, 3, 15, 63, 51, 195, 207, 255, 243, 67, 79, 127, 115, 131, 143, 191, 179, 7, 11, 31, 47, 59, 55, 35, 19, 23,
      27, 43, 39, 199, 203, 223, 239, 251, 247, 227, 211, 215, 219, 235, 231, 71, 75, 95, 111, 123, 119, 99, 83, 87, 91,
      107, 103, 135, 139, 159, 175, 187, 183, 163, 147, 151, 155, 171, 167,

      -2, 0, 48, 240, 192, 3, 51, 243, 195, 1, 49, 241, 193, 2, 50, 242, 194, 16, 32, 112, 176, 224, 208, 128, 64, 80,
      96, 160, 144, 19, 35, 115, 179, 227, 211, 131, 67, 83, 99, 163, 147, 17, 33, 113, 177, 225, 209, 129, 65, 81, 97,
      161, 145, 18, 34, 114, 178, 226, 210, 130, 66, 82, 98, 162, 146,

      +2, 12, 60, 252, 204, 15, 63, 255, 207, 13, 61, 253, 205, 14, 62, 254, 206, 28, 44, 124, 188, 236, 220, 140, 76,
      92, 108, 172, 156, 31, 47, 127, 191, 239, 223, 143, 79, 95, 111, 175, 159, 29, 45, 125, 189, 237, 221, 141, 77,
      93, 109, 173, 157, 30, 46, 126, 190, 238, 222, 142, 78, 94, 110, 174, 158,

      -3, 0, 192, 195, 3, 12, 204, 207, 15, 4, 196, 199, 7, 8, 200, 203, 11, 64, 128, 193, 194, 131, 67, 2, 1, 65, 129,
      130, 66, 76, 140, 205, 206, 143, 79, 14, 13, 77, 141, 142, 78, 68, 132, 197, 198, 135, 71, 6, 5, 69, 133, 134, 70,
      72, 136, 201, 202, 139, 75, 10, 9, 73, 137, 138, 74,

      +3, 48, 240, 243, 51, 60, 252, 255, 63, 52, 244, 247, 55, 56, 248, 251, 59, 112, 176, 241, 242, 179, 115, 50, 49,
      113, 177, 178, 114, 124, 188, 253, 254, 191, 127, 62, 61, 125, 189, 190, 126, 116, 180, 245, 246, 183, 119, 54,
      53, 117, 181, 182, 118, 120, 184, 249, 250, 187, 123, 58, 57, 121, 185, 186, 122,

      -4, 0, 3, 15, 12, 48, 51, 63, 60, 16, 19, 31, 28, 32, 35, 47, 44, 1, 2, 7, 11, 14, 13, 8, 4, 5, 6, 10, 9, 49, 50,
      55, 59, 62, 61, 56, 52, 53, 54, 58, 57, 17, 18, 23, 27, 30, 29, 24, 20, 21, 22, 26, 25, 33, 34, 39, 43, 46, 45,
      40, 36, 37, 38, 42, 41,

      +4, 192, 195, 207, 204, 240, 243, 255, 252, 208, 211, 223, 220, 224, 227, 239, 236, 193, 194, 199, 203, 206, 205,
      200, 196, 197, 198, 202, 201, 241, 242, 247, 251, 254, 253, 248, 244, 245, 246, 250, 249, 209, 210, 215, 219, 222,
      221, 216, 212, 213, 214, 218, 217, 225, 226, 231, 235, 238, 237, 232, 228, 229, 230, 234, 233,
  };


  switch (n_nodes()) {
    case 16:
      return B16n4DCheckArray;
    case 81:
      return B81n4DCheckArray;
    case 256:
      return B256n4DCheckArray;
    default:
      throw NotImplementedException();
  }
}

int ELEM4dTesseract::GetNodesPerSurface() const {
  switch (n_nodes()) {
    case 16:
      return 8;
    case 81:
      return 27;
    case 256:
      return 64;
    default:
      throw NotImplementedException();
  }
}

ZEROPTV ELEM4dTesseract::CalculateNormal(const GRID* grid, int surf_id) const {
  const int* nodes = GetNodesInSurface(surf_id);
  const ZEROPTV &p1 = grid->GetNode(ElemToLocalNodeID(nodes[0]))->location();
  const ZEROPTV &p2 = grid->GetNode(ElemToLocalNodeID(nodes[1]))->location();
  const ZEROPTV &p3 = grid->GetNode(ElemToLocalNodeID(nodes[2]))->location();
  const ZEROPTV &p4 = grid->GetNode(ElemToLocalNodeID(nodes[3]))->location();
  ZEROPTV a(p2 - p1);
  ZEROPTV b(p3 - p1);
  ZEROPTV c(p4 - p1);

  ZEROPTV normal;
  normal.crossProduct(p2 - p1, p3 - p1, p4 - p1);

  normal.Normalize();
  return normal;
  // throw NotImplementedException();
}

kBasisFunction ELEM4dTesseract::basis_function() const {
  switch (n_nodes()) {
    case 16:
      return BASIS_LINEAR;
    case 81:
      return BASIS_QUADRATIC;
    case 256:
      return BASIS_CUBIC;
    default:
      throw NotImplementedException();
  }
}

}  // namespace TALYFEMLIB
