/*********************************************************************
 *
 *  hollow_cylinder
 *
 *********************************************************************/

a = 1.0;
b = 2.0;

ndiv = 32;

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

Physical Curve("Γᶿ") = {1,3};
Physical Curve("Γʳ") = {2,4};
Physical Surface("Ω") = {1};

// 10 -> 1.07177, 20 -> 1.035, 40 -> 1.01747, 80 -> 1.0087
// 8 -> 1.0905, 16 -> 1.04427378, 32 -> 1.021897149, 64 -> 1.010889286
Transfinite Curve{1,-3} = ndiv+1 Using Progression 1.021897149;
// Transfinite Curve{1,-3} = ndiv+1;
Transfinite Curve{2,4} = 2*ndiv+1;
Transfinite Surface{1} = {1,2,3,4};
Mesh.Algorithm = 8;
Mesh.MshFileVersion = 2;
Mesh 2;
