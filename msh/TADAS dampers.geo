a = 12.0;
b = 3.0;
c = 5.0;
d = 28.0;
e = 14.2;
n =3;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {   120,0.0,0.0};
Point(3) = {   120,  30.0,0.0};
Point(4) = {90,30,0.0};
Point(5) = {75,80,0.0};
Point(6) = {131,360,0.0};
Point(7) = {-11,360,0.0};
Point(8) = {45,80,0.0};
Point(9) = {30,30,0.0};
Point(10) = {0.0,30,0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,1};

Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10};
Plane Surface(1) = {1};

Transfinite Curve{1} = 4*n;
Transfinite Curve{2,10} = n;
Transfinite Curve{3,9} = n;
Transfinite Curve{4,8} = 2*n;
Transfinite Curve{5,7} = 9*n;
Transfinite Curve{6} = 4*n;

Physical Curve("Γᵍ") = {6};
Physical Curve("Γᵗ") = {1};


Physical Curve("Γ2") = {2};
Physical Curve("Γ3") = {3};
Physical Curve("Γ4") = {4};
Physical Curve("Γ5") = {5};
Physical Curve("Γ7") = {7};
Physical Curve("Γ8") = {8};
Physical Curve("Γ9") = {9};
Physical Curve("Γ10") = {10};


Physical Point("Γₚ₁") = {1};
Physical Point("Γₚ₂") = {2};
Physical Point("Γₚ₃") = {3};
Physical Point("Γₚ₄") = {4};
Physical Point("Γₚ₅") = {5};
Physical Point("Γₚ₆") = {6};
Physical Point("Γₚ₇") = {7};
Physical Point("Γₚ₈") = {8};
Physical Point("Γₚ₉") = {9};
Physical Point("Γₚ₁₀") = {10};
Physical Surface("Ω") = {1};


Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;