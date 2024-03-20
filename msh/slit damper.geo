//边框
Point(1) = {0.0,0.0,0.0};
Point(2) = { 420-180,0,0};
Point(3) = { 420-180,20,0};
Point(4) = {400-180,40,0};
Point(5) = {400-180,140,0};
Point(6) = {420-180,160,0};
Point(7) = {420-180,180,0};
Point(8) = {0,180,0};
Point(9) = {0,160,0};
Point(10)  = (20,140,0);
Point(11)  = (20,40,0);
Point(12)  = (0,20,0);

Point(100)  = (420-180,40,0);
Point(200)  = (420-180,140,0);
Point(300)  = (0,140,0);
Point(400)  = (0,40,0);

Line(1) = {1,2};
Line(2) = {2,3};
Circle(3)  = {3,100,4};
Line(4) = {4,5};
Circle(5)  = {5,200,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Circle(9)  = {9,300,10};
Line(10) = {10,11};
Circle(11)  = {11,400,12};
Line(12) = {12,1};

Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12};


//内部1

Point(21) = {100,40,0};
Point(22)  = (60,40,0);
Point(23)  = (60,140,0);
Point(24)  = (100,140,0);
Point(500)  = (80,40,0);
Point(600)  = (80,140,0);

Circle(21)  = {22,500,21};
Line(22) = {22,23};
Circle(23)  = {24,600,23};
Line(24) = {24,21};

Curve Loop(2) = {-21,22,-23,24};
//Plane Surface(1) = {1,2};


//内部2
Point(31) = {180,40,0};
Point(32)  = (140,40,0);
Point(33)  = (140,140,0);
Point(34)  = (180,140,0);
Point(700)  = (160,40,0);
Point(800)  = (160,140,0);

Circle(31)  = {32,700,31};
Line(32) = {32,33};
Circle(33)  = {34,800,33};
Line(34) = {34,31};

Curve Loop(3) = {-31,32,-33,34};
Plane Surface(1) = {1,2,3};

Transfinite Curve{1,7} = 80;
Transfinite Curve{10,22,24,32,34,4} = 40;
Transfinite Curve{9,5,11,3} = 8;
Transfinite Curve{23,33,21,31} = 32;


Physical Curve("Γᵍ") = {7};
Physical Curve("Γᵗ") = {1};
Physical Surface("Ω") = {1};

Mesh.Algorithm = 2;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;