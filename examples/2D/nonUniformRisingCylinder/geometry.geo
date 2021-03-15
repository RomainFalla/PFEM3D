
D=0.05;
N=1.;
R_inf=10.;
prof_cyl=5.;

Point(1) = {0, 0, 0, D/(20*N)};
Point(2) = {-R_inf*D, 0, 0, D/(N)};
Point(3) = {-R_inf*D, -(R_inf+prof_cyl)*D, 0, D/(N)};
Point(4) = {0, -(R_inf+prof_cyl)*D, 0, D/(N)};
Point(5) = {R_inf*D, -(R_inf+prof_cyl)*D, 0, D/(N)};
Point(6) = {R_inf*D, 0, 0, D/(N)};

Point(7) = {0, -prof_cyl*D, 0, D/(20*N)};
Point(8) = {0, (-prof_cyl-0.5)*D, 0, D/(20*N)};
Point(9) = {0, (-prof_cyl+0.5)*D, 0, D/(20*N)};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Circle(7) = {8, 7, 9};
Circle(8) = {9, 7, 8};
Line Loop(9) = {1,2,3,4,5,6};
Line Loop(10) = {7,8};
Plane Surface(11)={9,10};
Physical Curve("Boundary")={2,3,4,5};
Physical Curve("DiskBoundary")={7,8};
Physical Curve("FreeSurface")={6,1};
Physical Surface("Fluid")={11};
//+


