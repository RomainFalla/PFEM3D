//+
d = 0.5;
//+
L = 1;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {0, L, 0, d};
//+
Point(3) = {L, 0, 0, d};
//+
Point(4) = {L, L, 0, d};
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
Physical Curve("Down") = {1};
//+
Physical Curve("Right") = {2};
//+
Physical Curve("Up") = {3};
//+
Physical Curve("Left") = {4};
//+
Physical Surface("Fluid") = {1};
//+
Transfinite Surface{1};
