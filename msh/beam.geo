
l = 10.0;
n = 80;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {  l, 0.0, 0.0};

Line(1) = {1,2};

Physical Curve("Ω") = {1};
Physical Point("Γ") = {1,2};

Transfinite Curve{1} = n+1;

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 1;
