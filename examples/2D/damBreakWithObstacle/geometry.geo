//+
L = 0.146;
//+
s = 0.024;
//+
d=s/4;
//+
Point(1) = {0, 0, 0, d};
//+
Point(2) = {4*L+s, 0, 0, d};
//+
Point(3) = {4*L+s, 4*L, 0, d};
//+
Point(4) = {0, 2.175*L, 0, d};
//+
Point(5) = {L, 0, 0, d};
//+
Point(6) = {L, 2*L, 0, d};
//+
Point(7) = {0, 2*L, 0, d};
//+
Point(8) = {2*L, 0, 0, d};
//+
Point(9) = {2*L, 2*s, 0, d};
//+
Point(10) = {2*L+s, 2*s, 0, d};
//+
Point(11) = {2*L+s, 0, 0, d};
//+
Line(1) = {4, 7};
//+
Line(2) = {7, 1};
//+
Line(3) = {1, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 7};
//+
Line(6) = {5, 8};
//+
Line(7) = {8, 9};
//+
Line(8) = {9, 10};
//+
Line(9) = {10, 11};
//+
Line(10) = {11, 2};
//+
Line(11) = {2, 3};
//+
Curve Loop(1) = {2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface{1};
//+
Physical Curve("Boundary") = {1, 2, 3, 6, 7, 8, 9, 10, 11};
//+
Physical Curve("FreeSurface") = {5, 4};
//+
Physical Surface("Fluid") = {1};
