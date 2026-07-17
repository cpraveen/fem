// Generate mesh for a rectangulat domain (xmin,xmax) x (ymin,ymax)

// Domain
xmin = -10.0;
xmax =  10.0;
ymin = -10.0;
ymax =  10.0;

// number of points on each side
n = 101;

Point(1) = {xmin, ymin, 0};
Point(2) = {xmax, ymin, 0};
Point(3) = {xmax, ymax, 0};
Point(4) = {xmin, ymax, 0};

Line(1) = {1,2}; // bottom
Line(2) = {2,3}; // right
Line(3) = {3,4}; // top
Line(4) = {4,1}; // left

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Transfinite Line{1,2,3,4} = n;

// Uncomment next line to get cartesian grid
//Transfinite Surface{1};

// Comment next line to get triangles only
Recombine Surface{1};

Periodic Curve {4} = {-2};
Periodic Curve {1} = {-3};

Physical Surface(100) = {1};

Physical Line(0) = {4}; // left
Physical Line(1) = {2}; // right
Physical Line(2) = {1}; // bottom
Physical Line(3) = {3}; // top
