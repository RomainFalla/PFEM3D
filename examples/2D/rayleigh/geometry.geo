//+
d = 0.02;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {0, 1, 0, d};
//+
Point(3) = {4, 0, 0, d};
//+
Point(4) = {4, 1, 0, d};
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
Plane Surface(1) = {1};
//+
Physical Curve("Top") = {2};
//+
Physical Curve("Bottom") = {4};
//+
Physical Curve("Left") = {1};
//+
Physical Curve("Right") = {3};
//+
Physical Surface("Fluid") = {1};
//+
Transfinite Surface{1};
