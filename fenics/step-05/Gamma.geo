n = 11;
h = 1/(n-1);

Point(1) = { 0, 0, 0, h};
Point(2) = { 1, 0, 0, h};
Point(3) = { 1, 1, 0, h};
Point(4) = {-1, 1, 0, h};
Point(5) = {-1,-1, 0, h};
Point(6) = { 0,-1, 0, h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

Line Loop(1) = {1,2,3,4,5,6};
Plane Surface(1) = {1};

Physical Surface(100000) = {1};
Physical Line   (100001) = {1,2,3,4,5,6};
