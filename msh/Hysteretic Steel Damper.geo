Point(1) = {0,0,0};
Point(2) = {0,40,0};
Point(3) = {40,40,0};
Point(4) = {70,70,0};
Point(5) = {140,160,0};
Point(6) = {70,250,0};
Point(7) = {40,280,0};
Point(8) = {0,280,0};
Point(9) = {0,320,0};
Point(10) = {320,320,0};
Point(11) = {320,280,0};
Point(12) = {280,280,0};
Point(13) = {250,250,0};
Point(14) = {180,160,0};
Point(15) = {250,70,0};
Point(16) = {280,40,0};
Point(17) = {320,40,0};
Point(18) = {320,0,0};

Point(20) = {40,70,0};
Point(21) = {40,250,0};
Point(22) = {280,250,0};
Point(23) = {280,70,0};





//连接点成直线或者曲线，使用Line和Circle arc功能。
Line(1) = {1, 2};
Line(2) = {2, 3};
Circle(3) = {4, 20, 3};
Line(4) = {4,5};
Line(5) = {5,6};
Circle(6) = {6,21,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Circle(12) = {12,22,13};
Line(13) = {13,14};
Line(14) = {14,15};
Circle(15) = {15,23,16};
Line(16) = {16,17};
Line(17) = {17,18};
Line(18) = {18,1};


//给前面定义的Line划分间距。这里的意思是将线等间距的添加10个点，即9个网格。命令在Mesh/Define/Transfinite/Curve。需要注意的是.geo文件中的Curve在openFoam转换的时候识别不是出来，只需把文件粘贴在算例文件夹之后，编辑文件时使用替换功能将所有的Curve改成Line。
Transfinite Curve {1} = 10 Using Progression 1;
Transfinite Curve {2} = 10 Using Progression 1;
Transfinite Curve {4} = 10 Using Progression 1;
Transfinite Curve {5} = 10 Using Progression 1;
Transfinite Curve {6} = 10 Using Progression 1;
Transfinite Curve {7} = 10 Using Progression 1;
Transfinite Curve {3} = 10 Using Progression 1;
Transfinite Curve {10} = 10 Using Progression 1;
Transfinite Curve {9} = 10 Using Progression 1;
Transfinite Curve {8} = 10 Using Progression 1;
Transfinite Curve {11} = 10 Using Progression 1;
Transfinite Curve {15} = 10 Using Progression 1;
Transfinite Curve {12} = 10 Using Progression 1;
Transfinite Curve {13} = 10 Using Progression 1;
Transfinite Curve {14} = 10 Using Progression 1;
Transfinite Curve {16} = 10 Using Progression 1;
Transfinite Curve {17} = 10 Using Progression 1;
Transfinite Curve {18} = 10 Using Progression 1;

//定义面，因为一开始画几何体的时候只画了二维的几何图形，其他部分拉通过拉伸得到。
Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
Plane Surface(1) = {1};

//为前面定义的面划分网格。
Transfinite Surface {1};


//生成结构性网格使用Recombine。
Recombine Surface {1};





