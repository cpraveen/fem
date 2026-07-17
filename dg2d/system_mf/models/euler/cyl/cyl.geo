// Generate mesh for annular region

// Comment following line to generate triangles
Mesh.RecombineAll=1;

ri = 1.0;      // inner radius
ro = 100.0*ri;  // outer radius

nin = 200;
hin = 2*Pi*ri/nin;

nout = 50;
hout = 2*Pi*ro/nout;

Point(1) = {0, 0, 0, hin};

Point(2) = {ri, 0, 0, hin};
Point(3) = {0, ri, 0, hin};
Point(4) = {-ri, 0, 0, hin};
Point(5) = {0, -ri, 0, hin};

Point(6) = {ro, 0, 0, hout};
Point(7) = {0, ro, 0, hout};
Point(8) = {-ro, 0, 0, hout};
Point(9) = {0, -ro, 0, hout};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};

Line(9) = {2,6};
Line(10) = {3,7};
Line(11) = {4,8};
Line(12) = {5,9};

Line Loop(1) = {9, 5, -10, -1};
Plane Surface(1) = {1};

Line Loop(2) = {10, 6, -11, -2};
Plane Surface(2) = {2};

Line Loop(3) = {11, 7, -12, -3};
Plane Surface(3) = {3};

Line Loop(4) = {12, 8, -9, -4};
Plane Surface(4) = {4};

//Transfinite Line{1,2,3,4} = nin/4;     // inner boundary
//Transfinite Line{5,6,7,8} = nout/4;    // outer boundary

Physical Surface(100) = {1,2,3,4};
Physical Line(101) = {1,2,3,4}; // inner circlular boundary
Physical Line(102) = {5,6,7,8}; // outer circlular boundary

//Geometry.Normals = 100;
