//+
d = 0.1;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {0.25, 0, 0, d};
//+
Point(3) = {-0.25, 0, 0, d};
//+
Point(4) = {0, 0.25, 0, d};
//+
Point(5) = {0, -0.25, 0, d};
//+
Point(6) = {0, -1, 0, d};
//+
Point(7) = {0, 1, 0, d};
//+
Point(8) = {1, 0, 0, d};
//+
Point(9) = {-1, 0, 0, d};
//+
Circle(1) = {3, 1, 4};
//+
Circle(2) = {4, 1, 2};
//+
Circle(3) = {2, 1, 5};
//+
Circle(4) = {5, 1, 3};
//+
Circle(5) = {6, 1, 9};
//+
Circle(6) = {9, 1, 7};
//+
Circle(7) = {7, 1, 8};
//+
Circle(8) = {8, 1, 6};
//+
Line(9) = {9, 3};
//+
Line(10) = {4, 7};
//+
Line(11) = {2, 8};
//+
Line(12) = {5, 6};
//+
Curve Loop(1) = {6, -10, -1, -9};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {2, 11, -7, -10};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {3, 12, -8, -11};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {9, -4, 12, 5};
//+
Surface(4) = {4};
//+
Physical Curve("MovingBoundary") = {6, 5, 8, 7};
//+
Physical Curve("StaticBoundary") = {1, 4, 3, 2};
//+
Physical Surface("Fluid") = {4, 1, 2, 3};
