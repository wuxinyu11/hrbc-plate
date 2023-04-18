
a = 10.0;
b = 10.0;
n = 80;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {  a, 0.0, 0.0};
Point(3) = {  a,   b, 0.0};
Point(4) = {0.0,   b, 0.0};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Transfinite Curve{1,2,3,4} = n+1;
Transfinite Surface{1};
// Physical Curve("Γᵗ₁") = {2};
// Physical Curve("Γᵗ₂") = {3};
Physical Point("Γₚ₁") = {1};
Physical Point("Γₚ₂") = {2};
Physical Point("Γₚ₃") = {3};
Physical Point("Γₚ₄") = {4};
Physical Curve("Γ₁") = {1};
Physical Curve("Γ₂") = {2};
Physical Curve("Γ₃") = {3};
Physical Curve("Γ₄") = {4};
Physical Surface("Ω") = {1};

Mesh.Algorithm = 9;
Mesh.MshFileVersion = 2;
Mesh 2;
// RecombineMesh;
