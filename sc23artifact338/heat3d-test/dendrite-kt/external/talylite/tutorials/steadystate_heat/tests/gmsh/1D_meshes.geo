Point(1) = {0, 0, 0, 0.002};
Point(2) = {1, 0, 0, 0.020};
Line(1) = {1, 2};
Physical Line(3) = {1};
Physical Point(1) = {1};
Physical Point(2) = {2};

Mesh.ElementOrder = 1;
// Rotate {{0, 1, 1}, {0, 0, 0}, Pi/4} {
//   Line{1};
// }
Rotate {{0, 1, 0}, {0, 0, 0}, -Pi/2} {
  Line{1};
}
