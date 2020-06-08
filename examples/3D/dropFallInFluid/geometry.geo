//+
d = 0.025;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {0.5, 0, 0, d};
//+
Point(3) = {-0.5, 0, 0, d};
//+
Point(4) = {0, 0.5, 0, d};
//+
Point(5) = {0, -0.5, 0, d};
//+
Point(6) = {0, 0, 1, d};
//+
Point(7) = {0.5, 0, 1, d};
//+
Point(8) = {-0.5, 0, 1, d};
//+
Point(9) = {0, 0.5, 1, d};
//+
Point(10) = {0, -0.5, 1, d};
//+
Point(11) = {0, 0, 0.3, d};
//+
Point(12) = {0.5, 0, 0.3, d};
//+
Point(13) = {-0.5, 0, 0.3, d};
//+
Point(14) = {0, 0.5, 0.3, d};
//+
Point(15) = {0, -0.5, 0.3, d};
//+
Point(16) = {0, 0, 0.7, d};
//+
Point(17) = {-0.1, 0, 0.7, d};
//+
Point(18) = {0.1, 0, 0.7, d};
//+
Point(19) = {0, 0.1, 0.7, d};
//+
Point(20) = {0, -0.1, 0.7, d};
//+
Point(21) = {0, 0, 0.8, d};
//+
Point(22) = {0, 0, 0.6, d};
//+
Circle(1) = {5, 1, 3};
//+
Circle(2) = {3, 1, 4};
//+
Circle(3) = {4, 1, 2};
//+
Circle(4) = {2, 1, 5};
//+
Circle(5) = {15, 11, 13};
//+
Circle(6) = {13, 11, 14};
//+
Circle(7) = {14, 11, 12};
//+
Circle(8) = {12, 11, 15};
//+
Circle(9) = {10, 6, 8};
//+
Circle(10) = {8, 6, 9};
//+
Circle(11) = {9, 6, 7};
//+
Circle(12) = {7, 6, 10};
//+
Circle(13) = {20, 16, 17};
//+
Circle(14) = {17, 16, 19};
//+
Circle(15) = {19, 16, 18};
//+
Circle(16) = {18, 16, 20};
//+
Circle(17) = {20, 16, 21};
//+
Circle(18) = {21, 16, 19};
//+
Circle(19) = {19, 16, 22};
//+
Circle(20) = {22, 16, 20};
//+
Circle(21) = {18, 16, 21};
//+
Circle(22) = {21, 16, 17};
//+
Circle(23) = {17, 16, 22};
//+
Circle(24) = {22, 16, 18};
//+
Line(25) = {5, 15};
//+
Line(26) = {15, 10};
//+
Line(27) = {2, 12};
//+
Line(28) = {12, 7};
//+
Line(29) = {4, 14};
//+
Line(30) = {14, 9};
//+
Line(31) = {3, 13};
//+
Line(32) = {8, 13};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {29, -6, -31, 2};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {1, 31, -5, -25};
//+
Surface(4) = {4};
//+
Curve Loop(5) = {25, -8, -27, 4};
//+
Surface(5) = {5};
//+
Curve Loop(6) = {27, -7, -29, 3};
//+
Surface(6) = {6};
//+
Curve Loop(7) = {28, -11, -30, 7};
//+
Surface(7) = {7};
//+
Curve Loop(8) = {30, -10, 32, 6};
//+
Surface(8) = {8};
//+
Curve Loop(9) = {32, -5, 26, 9};
//+
Surface(9) = {9};
//+
Curve Loop(10) = {8, 26, -12, -28};
//+
Surface(10) = {10};
//+
Curve Loop(11) = {22, -13, 17};
//+
Surface(11) = {11};
//+
Curve Loop(12) = {23, 20, 13};
//+
Surface(12) = {12};
//+
Curve Loop(13) = {17, -21, 16};
//+
Surface(13) = {13};
//+
Curve Loop(14) = {20, -16, -24};
//+
Surface(14) = {14};
//+
Curve Loop(15) = {21, 18, 15};
//+
Surface(15) = {15};
//+
Curve Loop(16) = {24, -15, 19};
//+
Surface(16) = {16};
//+
Curve Loop(17) = {18, -14, -22};
//+
Surface(17) = {17};
//+
Curve Loop(18) = {19, -23, 14};
//+
Surface(18) = {18};
//+
Surface Loop(1) = {1, 4, 3, 6, 5, 2};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {16, 14, 12, 18, 17, 15, 13, 11};
//+
Volume(2) = {2};
//+
Physical Surface("Boundary") = {4, 9, 10, 5, 7, 6, 8, 3, 1};
//+
Physical Surface("FreeSurface") = {2, 14, 17, 15, 13, 11, 18, 16, 12};
//+
Physical Volume("Fluid") = {2, 1};
//+
Transfinite Surface{3, 4, 5, 6, 7, 8, 9, 10};
