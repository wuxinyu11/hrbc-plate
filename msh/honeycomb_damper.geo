Point(1) = {0,0,0};
Point(2)= {560,0,0};
Point(3)= {560,40,0};
Point(4) =  {500,40,0};
Point(5)= {460,180,0};
Point(6) ={500,320,0};
Point(7) = {560,320,0};
Point(8)= {560,360,0};
Point(9) = {0,360,0};
Point(10) ={0,320,0};
Point(11) ={60,320,0};
Point(12) ={100,180,0};
Point(13) ={60,40,0};
Point(14) ={0,40,0};

Point(15) = {180,40,0};
Point(16) = {140,180,0};
Point(17) = {180,320,0};
Point(18) = {220,320,0};
Point(19) = {260,180,0};
Point(20) = {220,40,0};

Point(21) = {340,40,0};
Point(22) = {300,180,0};
Point(23) = {340,320,0};
Point(24) = {380,320,0};
Point(25) = {420,180,0};
Point(26) = {380,40,0};

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
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};
Line(20) = {20, 15};
Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 24};
Line(24) = {24, 25};
Line(25) = {25, 26};
Line(26) = {26, 21};


Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
Curve Loop(2) = {15,16,17,18,19,20};
Curve Loop(3) = {21,22,23,24,25,26};
Plane Surface(1) = {1,2,3};

Transfinite Curve{1,8} = 100;
Transfinite Curve{13,26,20,3,2,14,10,17,23,6,7,9} =15;
Transfinite Curve{12,15,19,21,25,4,5,24,22,18,16,11} = 75;


Physical Curve("Γᵍ") = {8};
Physical Curve("Γᵗ") = {1};
Physical Surface("Ω") = {1};


Mesh.Algorithm =2;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;




