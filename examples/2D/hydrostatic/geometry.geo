//+
d = 0.25;
//+
Point(1) = {10, 0, 0, d};
//+
Point(2) = {0, 0, 0, d};
//+
Point(3) = {0, 10, 0, d};
//+
Point(4) = {10, 10, 0, d};
//+
Point(5) = {0, 5, 0, d};
//+
Point(6) = {10, 5, 0, d};
//+
Line(1) = {3, 5};
//+
Line(2) = {5, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 6};
//+
Line(5) = {6, 4};
//+
Line(6) = {5, 6};
//+
Physical Curve("Boundary") = {1, 2, 3, 4, 5};
//+
Physical Curve("FreeSurface") = {6};
//+
Curve Loop(1) = {6, -4, -3, -2};
//+
Plane Surface(1) = {1};
//+
Physical Surface("Fluid") = {1};
//+
Transfinite Surface{1};
