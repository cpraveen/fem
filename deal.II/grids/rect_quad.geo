// Generate mesh for a rectangulat domain (xmin,xmax) x (ymin,ymax)

// Domain
xmin = 0.0;
xmax = 1.0;
ymin = 0.0;
ymax = 1.0;

// number of points on each side
n = 51;

Point(1) = {xmin, ymin, 0};
Point(2) = {xmax, ymin, 0};
Point(3) = {xmax, ymax, 0};
Point(4) = {xmin, ymax, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Transfinite Line{1,2,3,4} = n;

// Comment next line to get triangles only
Recombine Surface{1};

Physical Line(0) = {4}; // left
Physical Line(1) = {2}; // right
Physical Line(2) = {1}; // bottom
Physical Line(3) = {3}; // top
