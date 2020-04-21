//+
d = 0.025;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {1, 0, 0, d};
//+
Point(3) = {1, 1, 0, d};
//+
Point(4) = {0, 1, 0, d};
//+
Point(5) = {0, 1.2, 0, d};
//+
Point(6) = {1, 1.2, 0, d};
//+
Point(7) = {0.5, 0.15, 0, d};
//+
Point(9) = {0.55, 0.2, 0, d};
//+
Point(10) = {0.5, 0.25, 0, d};
//+
Point(11) = {0.45, 0.2, 0, d};
//+
Point(12) = {0.5, 0.2, 0, d};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 5};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 2};
//+
Line(5) = {2, 1};
//+
Line(6) = {3, 6};
//+
Circle(7) = {11, 12, 10};
//+
Circle(8) = {10, 12, 9};
//+
Circle(9) = {9, 12, 7};
//+
Circle(10) = {7, 12, 11};
//+
Curve Loop(1) = {1, 3, 4, 5};
//+
Curve Loop(2) = {7, 8, 9, 10};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("Boundary") = {2, 1, 5, 4, 6};
//+
Physical Curve("FreeSurface") = {3};
//+
Physical Curve("DiskBoundary") = {7, 10, 9, 8};
//+
Physical Surface("Fluid") = {1};
