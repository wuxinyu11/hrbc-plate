
Point(1) = {0.0,0.0,0.0};
Point(2) = { 20,20,0};
Point(3) = {  20,120,0};
Point(4) = {0,140,0};
Point(5) = {80,140,0};
Point(6) = {60,120,0};
Point(7) = {60,20,0};
Point(8) = {80,0,0};


Point(10)  = (0,20,0);
Point(11)  = (0,120,0);
Point(12)  = (80,120,0);
Point(13)  = (80,20,0);

Circle(1)  = {2,10,1};
Line(2) = {2,3};
Circle(3)  = {3,11,4};
Line(4) = {4,5};
Circle(5)  = {6,12,5};
Line(6) = {6,7};
Circle(7)  = {7,13,8};
Line(8) = {8,1};



Curve Loop(1) = {1,2,3,4,5,6,7,8};
Plane Surface(1) = {1};
Transfinite Surface {1};


Mesh.Algorithm = 8;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;