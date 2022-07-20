/*********************************************************************
 *
 *  hollow_cylinder
 *
 *********************************************************************/

a = 1.0;
b = 2.0;

lc = 1;
ndiv = 10;

Point(5) = {0 , 0., 0.};
Point(1) = {a , 0., 0.};
Point(2) = {b , 0., 0.};
Point(3) = {0., b,  0.};
Point(4) = {0., a,  0.};

Line(1)  = {1,2};
Circle(2)  = {2,5,3};
Line(3)  = {3,4};
Circle(4)  = {4,5,1};

Line Loop(5) = { 1, 2, 3, 4};

Plane Surface(1) = {5};

Physical Curve("Γᵍ",2) = {1,3};
Physical Curve("Γᴹ",3) = {2,4};
Physical Surface("Ω",1) = {1};
Transfinite Curve{1} = ndiv+1;
Transfinite Curve{2} = 2*ndiv+1;
Transfinite Curve{3} = ndiv+1;
Transfinite Curve{4} = 2*ndiv+1;
Transfinite Surface{1} = {1,2,3,4};
Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;