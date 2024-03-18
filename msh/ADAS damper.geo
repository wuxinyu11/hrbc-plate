Point(1) = {0,0,0};
Point(2)= {240,0,0};
Point(3)= {240,40,0};
Point(4) =  {180,40,0};
Point(5)= {140,180,0};
Point(6) ={180,320,0};
Point(7) = {240,320,0};
Point(8)= {240,360,0};
Point(9) = {0,360,0};
Point(10) ={0,320,0};
Point(11) ={60,320,0};
Point(12) ={100,180,0};
Point(13) ={60,40,0};
Point(14) ={0,40,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13, 14};
Line(14) = {14, 1};


Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
Plane Surface(1) = {1};

Physical Curve("Γᵍ") = {8};
Physical Curve("Γᵗ") = {1};
Physical Surface("Ω") = {1};


Mesh.Algorithm = 8;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;




