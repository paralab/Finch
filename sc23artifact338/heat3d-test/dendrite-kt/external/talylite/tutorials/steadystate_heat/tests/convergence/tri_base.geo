Point(1) = {0, 0, 0, cl};
Point(2) = {1, 0, 0, cl};
Point(3) = {1, 1, 0, cl};
Point(4) = {0, 1, 0, cl};
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};

//Transfinite Line "*" = 20;
//Transfinite Surface "*";
//Transfinite Volume "*";

Physical Line(1) = {1};
Physical Line(3) = {2};
Physical Line(2) = {3};
Physical Line(4) = {4};
Physical Surface(7) = {6};
