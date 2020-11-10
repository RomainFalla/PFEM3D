//+
d = 0.25;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {1, 0, 0, d};
//+
Point(3) = {0, 1, 0, d};
//+
Point(4) = {0, 0, 1, d};
//+
Point(5) = {1, 1, 0, d};
//+
Point(6) = {0, 1, 1, d};
//+
Point(7) = {1, 0, 1, d};
//+
Point(8) = {1, 1, 1, d};
//+
Line(1) = {4, 6};
//+
Line(2) = {6, 3};
//+
Line(3) = {3, 1};
//+
Line(4) = {1, 4};
//+
Line(5) = {4, 7};
//+
Line(6) = {7, 2};
//+
Line(7) = {2, 1};
//+
Line(8) = {3, 5};
//+
Line(9) = {5, 2};
//+
Line(10) = {6, 8};
//+
Line(11) = {8, 7};
//+
Line(12) = {8, 5};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, -11, -10, -1};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 6, -9, -12};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {8, -12, -10, 2};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {3, -7, -9, -8};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {4, 5, 6, 7};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {1, 2, 6, 3, 5, 4};
//+
Volume(1) = {1};
//+
Physical Surface("FreeSurface") = {1, 2, 5, 3, 4, 6};
//+
Physical Volume("Fluid") = {1};
//+
Transfinite Surface{1, 2, 3, 4, 5, 6};
//+
Transfinite Volume{1};
