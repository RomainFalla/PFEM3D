//+
d = 0.1;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {0, 1, 0, d};
//+
Point(3) = {1, 0, 0, d};
//+
Point(4) = {1, 1, 0, d};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Surface(1) = {1};
//+
Physical Curve("FreeSurface") = {1};
//+
Physical Surface("Fluid") = {1};
//+
Transfinite Surface{1};
