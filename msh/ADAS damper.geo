Point(1) = {0,0,0};
Point(2)= {320,0,0};
Point(3)= {320,40,0};
Point(4) =  {280,40,0};
Point(5)= {250,70,0};
Point(6)= {180,160,0};
Point(7) ={250,250,0};
Point(8) = {280,280,0};
Point(9)= {320,280,0};
Point(10) = {320,320,0};
Point(11) ={0,320,0};
Point(12) ={0,280,0};
Point(13) ={40,280,0};
Point(14) ={70,250,0};
Point(15) ={140,160,0};
Point(16) = {70,70,0};
Point(17) =  {40,40,0};
Point(18) = {0,40,0};



Point(23) = {40,70,0};
Point(22) = {40,250,0};
Point(21) = {280,250,0};
Point(20) = {280,70,0};



Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Circle(4) = {4,20,5};
Line(5) = {5,6};
Line(6) = {6,7};
Circle(7) = {7,21,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,13};
Circle(13) = {13,22,14};
Line(14) = {14,15};
Line(15) = {15,16};
Circle(16) = {16, 23, 17};
Line(17) = {17, 18};
Line(18) = {18, 1};


Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
Plane Surface(1) = {1};

Physical Curve("Γᵍ") = {10};
Physical Curve("Γᵗ") = {1};
Physical Surface("Ω") = {1};


Mesh.Algorithm = 8;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;



