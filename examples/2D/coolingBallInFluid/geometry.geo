//+
d = 0.02;
//+
ballH = 1.25;
//+
ballR = 0.125;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {0, 1, 0, d};
//+
Point(3) = {4, 0, 0, d};
//+
Point(4) = {4, 1, 0, d};
//+
Point(5) = {0, 2.5, 0, d};
//+
Point(6) = {4, 2.5, 0, d};
//+
Point(7) = {1, ballH, 0, d};
//+
Point(8) = {2, ballH, 0, d};
//+
Point(9) = {3, ballH, 0, d};
//+
Point(10) = {1 - ballR, ballH, 0, d};
//+
Point(11) = {1 + ballR, ballH, 0, d};
//+
Point(12) = {2 - ballR, ballH, 0, d};
//+
Point(13) = {2 + ballR, ballH, 0, d};
//+
Point(14) = {3 - ballR, ballH, 0, d};
//+
Point(15) = {3 + ballR, ballH, 0, d};
//+
Line(1) = {5, 2};
//+
Line(2) = {2, 1};
//+
Line(3) = {1, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 6};
//+
Line(6) = {2, 4};
//+
Circle(7) = {10, 7, 11};
//+
Circle(8) = {11, 7, 10};
//+
Circle(9) = {12, 8, 13};
//+
Circle(10) = {13, 8, 12};
//+
Circle(11) = {14, 9, 15};
//+
Circle(12) = {15, 9, 14};
//+
Curve Loop(1) = {2, 3, 4, -6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {8, 7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, 9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, 11};
//+
Plane Surface(4) = {4};
//+
Physical Curve("Boundary") = {1, 2, 3, 4, 5};
//+
Physical Curve("FreeSurface") = {6, 7, 8, 10, 9, 11, 12};
//+
Physical Surface("Fluid") = {1, 2, 3, 4};
//+
Transfinite Surface{1};
