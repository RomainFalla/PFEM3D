//+
d = 0.05;
//+
Point(1) = {0, 0, 1, d};
//+
Point(2) = {0, 0.5, 0.5, d};
//+
Point(3) = {0, -0.5, 0.5, d};
//+
Point(4) = {0, 0, 0, d};
//+
Point(5) = {0, 0, 0.5, d};
//+
Point(6) = {4, 0, 1, d};
//+
Point(7) = {4, 0.5, 0.5, d};
//+
Point(8) = {4, -0.5, 0.5, d};
//+
Point(9) = {4, 0, 0, d};
//+
Point(10) = {4, 0, 0.5, d};
//+
Point(11) = {0.2, 0, 1, d};
//+
Point(12) = {0.2, 0.5, 0.5, d};
//+
Point(13) = {0.2, -0.5, 0.5, d};
//+
Point(14) = {0.2, 0, 0, d};
//+
Point(15) = {0.2, 0, 0.5, d};
//+
Circle(1) = {3, 5, 1};
//+
Circle(2) = {1, 5, 2};
//+
Circle(3) = {2, 5, 4};
//+
Circle(4) = {4, 5, 3};
//+
Circle(5) = {13, 15, 11};
//+
Circle(6) = {11, 15, 12};
//+
Circle(7) = {12, 15, 14};
//+
Circle(8) = {14, 15, 13};
//+
Circle(9) = {8, 10, 6};
//+
Circle(10) = {6, 10, 7};
//+
Circle(11) = {7, 10, 9};
//+
Circle(12) = {9, 10, 8};
//+
Line(13) = {1, 11};
//+
Line(14) = {11, 6};
//+
Line(15) = {3, 13};
//+
Line(16) = {4, 14};
//+
Line(17) = {2, 12};
//+
Line(18) = {13, 8};
//+
Line(19) = {12, 7};
//+
Line(20) = {14, 9};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(4) = {2, 17, -6, -13};
//+
Surface(4) = {4};
//+
Curve Loop(5) = {1, 13, -5, -15};
//+
Surface(5) = {5};
//+
Curve Loop(6) = {15, -8, -16, 4};
//+
Surface(6) = {6};
//+
Curve Loop(7) = {16, -7, -17, 3};
//+
Surface(7) = {7};
//+
Curve Loop(8) = {5, 14, -9, -18};
//+
Surface(8) = {8};
//+
Curve Loop(9) = {8, 18, -12, -20};
//+
Surface(9) = {9};
//+
Curve Loop(10) = {19, -10, -14, 6};
//+
Surface(10) = {10};
//+
Curve Loop(11) = {11, -20, -7, 19};
//+
Surface(11) = {11};
//+
Surface Loop(1) = {5, 1, 4, 7, 6, 2};
//+
Volume(1) = {1};
//+
Physical Surface("FluidInput") = {1};
//+
Physical Surface("Boundary") = {4, 10, 6, 9, 7, 5, 8, 11};
//+
Physical Surface("FreeSurface") = {2};
//+
Physical Volume("Fluid") = {1};
//+
Transfinite Surface{4, 5, 6, 7, 8, 9, 10, 11};