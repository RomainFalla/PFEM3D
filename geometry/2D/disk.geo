//+
d = 0.05;
//+
Point(1) = {10, 0, 0, d};
//+
Point(2) = {0, 0, 0, d};
//+
Point(3) = {0, 10, 0, d};
//+
Point(4) = {10, 10, 0, d};

//+
Line(1) = {3, 2};
//+
Line(2) = {2, 1};
//+
Line(3) = {1, 4};
//+
Physical Curve("Boundary") = {1, 2, 3};

//+
Point(5) = {5, 7.5, 0, d};
//+
Point(6) = {6, 7.5, 0, d};
//+
Point(7) = {4, 7.5, 0, d};
//+
Circle(4) = {7, 5, 6};
//+
Circle(5) = {6, 5, 7};
//+
Curve Loop(1) = {5, 4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("FreeSurface") = {5, 4};
//+
Physical Surface("Fluid") = {1};
