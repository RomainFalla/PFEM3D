//+
d = 0.1;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {0, 1, 0, d};
//+
Point(3) = {4, 0, 0, d};
//+
Point(4) = {4, 1, 0, d};
//+
Point(5) = {0.2, 0, 0, d};
//+
Point(6) = {0.2, 1, 0, d};
//+
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 5};
//+
Line(3) = {5, 6};
//+
Line(4) = {6, 2};
//+
Line(5) = {6, 4};
//+
Line(6) = {3, 5};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Surface(1) = {1};
//+
Physical Curve("Boundary") = {5, 6, 2, 4};
//+
Physical Curve("FluidInput") = {1};
//+
Physical Curve("FreeSurface") = {3};
//+
Physical Surface("Fluid") = {1};
//+
Transfinite Surface{1};
