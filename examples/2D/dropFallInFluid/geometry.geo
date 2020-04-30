//+
L = 0.584;
//+
H = L/4;
//+
d = H/20;
//+
Point(1) = {L, 0, 0, d};
//+
Point(2) = {0, 0, 0, d};
//+
Point(3) = {0, L, 0, d};
//+
Point(4) = {L, L, 0, d};
//+
Point(5) = {L/2, 2.5*H, 0, d};
//+
Point(6) = {L/2+H/2, 2.5*H, 0, d};
//+
Point(7) = {L/2-H/2, 2.5*H, 0, d};
//+
Point(8) = {0, H, 0, d};
//+
Point(9) = {L, H, 0, d};
//+
Line(1) = {3, 8};
//+
Line(2) = {8, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 9};
//+
Line(5) = {9, 4};
//+
Line(6) = {8, 9};
//+
Circle(7) = {7, 5, 6};
//+
Circle(8) = {6, 5, 7};
//+
Physical Curve("Boundary") = {1, 2, 3, 4, 5};
//+
Physical Curve("FreeSurface") = {6, 8, 7};
//+
Curve Loop(1) = {2, 3, 4, -6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, 8};
//+
Plane Surface(2) = {2};
//+
Physical Surface("Fluid") = {1, 2};
//+
Transfinite Surface(1);
