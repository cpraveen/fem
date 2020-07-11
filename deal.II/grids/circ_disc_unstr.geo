Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 2; // 2 or 3
Mesh.Algorithm = 8;

n = 10;
h = 1.0/(n-1);

Point(1) = {-1,0,0,h};
Point(2) = { 0,-1,0,h};
Point(3) = { 1, 0,0,h};
Point(4) = { 0, 1,0,h};

Point(5) = { 0.0, 0.0, 0.0, h};

Circle(1) = {1,5,2};
Circle(2) = {2,5,3};
Circle(3) = {3,5,4};
Circle(4) = {4,5,1};

Point(6) = {-0.5, 0.0, 0.0, h};
Point(7) = { 0.0, -0.5, 0.0, h};
Point(8) = { 0.5, 0.0, 0.0, h};
Point(9) = { 0.0, 0.5, 0.0, h};

Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};

Line Loop(1) = {1,2,3,4}; // inner circle
Line Loop(2) = {5,6,7,8}; // outer circle

Plane Surface(1) = {1,2}; // annulus region
Plane Surface(2) = {2};   // inner disc region

Physical Surface(100) = {1,2};
Physical Line(200) = {1,2,3,4};
