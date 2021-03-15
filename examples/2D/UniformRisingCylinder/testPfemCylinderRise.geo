
D=0.05;
N=10.;
R_inf=5.5;
prof_cyl=2.5;

Point(1) = {0, 0, 0, D/(N)};
Point(2) = {-R_inf*D, 0, 0, D/(N)};
Point(3) = {-R_inf*D, -(R_inf+prof_cyl)*D, 0, D/(N)};
Point(4) = {0, -(R_inf+prof_cyl)*D, 0, D/(N)};
Point(5) = {R_inf*D, -(R_inf+prof_cyl)*D, 0, D/(N)};
Point(6) = {R_inf*D, 0, 0, D/(N)};

Point(7) = {0, -prof_cyl*D, 0, D/(N)};
Point(8) = {0, (-prof_cyl-0.5)*D, 0, D/(N)};
Point(9) = {0, (-prof_cyl+0.5)*D, 0, D/(N)};

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
Physical Line(12)={2,3,4,5};
Physical Line(13)={7,8};
Physical Surface(14)={11};
//+



