Point(1) = {0, 0, 0, 0.008};
Point(2) = {1, 0, 0, 0.02};
Point(3) = {1, 1, 0, 0.08};
Point(4) = {0, 1, 0, 0.02};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};

Physical Surface(7) = {6};

rotAngle = 0.204832765*Pi;
// rotAngle = 0;
Mesh.ElementOrder = 1;
 
// Recombine Surface {6};
// 	Physical Line(1) = {4};
// 	Physical Line(2) = {2};
// 	Physical Line(4) = {3};
// 	Physical Line(3) = {1};
// Rotate {{0, 1, 0}, {0, 0, 0}, 0} {
//   Surface{6};
// }


// Physical Line(5) = {4};
// Physical Line(6) = {2};
// Physical Line(4) = {3};
// Physical Line(3) = {1};
// Rotate {{0, 1, 0}, {0, 0, 0}, -Pi/2} {
//   Surface{6};
// }

Physical Line(1) = {4};
Physical Line(2) = {2};
Physical Line(6) = {3};
Physical Line(5) = {1};
Rotate {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Surface{6};
}



