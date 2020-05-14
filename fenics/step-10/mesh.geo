n = 20;
h = 1.0/(n-1);

Point(1) = {-1,-1,0,h};
Point(2) = { 1,-1,0,h};
Point(3) = { 1, 1,0,h};
Point(4) = {-1, 1,0,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Point(5) = { 0.0, 0.0, 0.0, h};
Point(6) = {-0.5, 0.0, 0.0, h};
Point(7) = { 0.0, -0.5, 0.0, h};
Point(8) = { 0.5, 0.0, 0.0, h};
Point(9) = { 0.0, 0.5, 0.0, h};

Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Line{5,6,7,8} In Surface{1};

Physical Surface(100) = {1};
Physical Line(200) = {1,2,3,4};
