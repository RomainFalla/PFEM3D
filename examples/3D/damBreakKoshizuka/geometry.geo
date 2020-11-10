//+
L = 0.146;
//+
b = 0.175;
//+
d = 0.146/10; 
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {L, 0, 0, d};
//+
Point(3) = {4*L, 0, 0, d};
//+
Point(4) = {0, 0, 2*L, d};
//+
Point(5) = {L, 0, 2*L, d};
//+
Point(6) = {L, 0, 4*L, d};
//+
Point(7) = {4*L, 0, 4*L, d};
//+
Point(8) = {0, 0, 4*L, d};
//+
Point(9) = {0, b, 0, d};
//+
Point(10) = {L, b, 0, d};
//+
Point(11) = {4*L, b, 0, d};
//+
Point(12) = {0, b, 2*L, d};
//+
Point(13) = {L, b, 2*L, d};
//+
Point(14) = {L, b, 4*L, d};
//+
Point(15) = {4*L, b, 4*L, d};
//+
Point(16) = {0, b, 4*L, d};
//+
Point(17) = {4*L, 0, 2*L, d};
//+
Point(18) = {4*L, b, 2*L, d};
//+
Line(1) = {9, 12};
//+
Line(2) = {12, 16};
//+
Line(3) = {16, 14};
//+
Line(4) = {14, 13};
//+
Line(5) = {13, 10};
//+
Line(6) = {10, 9};
//+
Line(7) = {12, 13};
//+
Line(8) = {13, 18};
//+
Line(9) = {18, 15};
//+
Line(10) = {15, 14};
//+
Line(11) = {18, 11};
//+
Line(12) = {11, 10};
//+
Line(13) = {1, 4};
//+
Line(14) = {4, 8};
//+
Line(15) = {8, 6};
//+
Line(16) = {6, 7};
//+
Line(17) = {7, 17};
//+
Line(18) = {17, 3};
//+
Line(19) = {8, 16};
//+
Line(20) = {10, 2};
//+
Line(21) = {4, 12};
//+
Line(22) = {4, 5};
//+
Line(23) = {5, 13};
//+
Line(24) = {5, 17};
//+
Line(25) = {1, 9};
//+
Line(26) = {1, 2};
//+
Line(27) = {2, 3};
//+
Line(28) = {3, 11};
//+
Line(29) = {17, 18};
//+
Line(30) = {7, 15};
//+
Line(31) = {5, 2};
//+
Line(32) = {6, 5};

//+
Curve Loop(1) = {13, 21, -1, -25};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {13, 22, 31, -26};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {24, 18, -27, -31};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {5, 20, -31, 23};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {6, -25, 26, -20};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {7, -23, -22, 21};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {14, 15, 32, -22};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {2, -19, -14, 21};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {24, -17, -16, 32};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {7, -4, -3, -2};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {10, 4, 8, 9};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {20, 27, 28, 12};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {29, 11, -28, -18};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {17, 29, 9, -30};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {11, 12, -5, 8};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {5, 6, 1, 7};
//+
Plane Surface(16) = {16};
//+
Surface Loop(1) = {6, 4, 1, 2, 5, 16};
//+
Volume(1) = {1};
//+
Physical Surface("Boundary") = {2, 3, 9, 7, 10, 11, 15, 16, 5, 12, 1, 8, 14, 13};
//+
Physical Surface("FreeSurface") = {6, 4};
//+
Physical Volume("Fluid") = {1};
//+
Transfinite Surface{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
//+
Transfinite Volume{1};
