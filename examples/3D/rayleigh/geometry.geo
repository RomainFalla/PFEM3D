//+
d = 0.1;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {0, 4, 0, d};
//+
Point(3) = {4, 0, 0, d};
//+
Point(4) = {4, 4, 0, d};
//+
Point(5) = {0, 0, 1, d};
//+
Point(6) = {0, 4, 1, d};
//+
Point(7) = {4, 0, 1, d};
//+
Point(8) = {4, 4, 1, d};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 2};
//+
Line(5) = {6, 5};
//+
Line(6) = {5, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 6};
//+
Line(9) = {6, 2};
//+
Line(10) = {5, 1};
//+
Line(11) = {7, 3};
//+
Line(12) = {8, 4};
//+
Curve Loop(1) = {5, 10, -1, -9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, -11, -6, 10};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 3, -12, -7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {8, 9, -4, -12};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {5, 6, 7, 8};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {1, 2, 3, 4};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {1, 5, 2, 6, 3, 4};
//+
Volume(1) = {1};
//+
Physical Surface("Top") = {5};
//+
Physical Surface("Bottom") = {6};
//+
Physical Surface("Lateral") = {1, 4, 3, 2};
//+
Physical Volume("Fluid") = {1};
