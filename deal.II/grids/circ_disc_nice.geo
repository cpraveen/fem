l = 0.1;
r1 = 0.5;  // circle of discontinuity
r2 = 1.0;

Point(newp) = {0, 0, 0};

// inner square
Point(newp) = { l,-l, 0};
Point(newp) = { l, l, 0};
Point(newp) = {-l, l, 0};
Point(newp) = {-l,-l, 0};

Line(newl) = {2,3};
Line(newl) = {3,4};
Line(newl) = {4,5};
Line(newl) = {5,2};

Point(newp) = {r1*Cos(-Pi/4), r1*Sin(-Pi/4), 0};
Point(newp) = {r1*Cos(+Pi/4), r1*Sin(+Pi/4), 0};
Line(newl)   = {2,6};
Line(newl)   = {3,7};
Circle(newl) = {6,1,7};

// right inner annulus
Line Loop(1) = {5,7,-6,-1};
Plane Surface(1) = {1};
Transfinite Surface(1) = {2,6,7,3};

// left inner annulus
Symmetry{1,0,0,0}{ Duplicata{Surface{1};} }
Transfinite Surface(8) = {4,14,9,5};
Reverse Surface{8};

// top inner annulus
Circle(newl) = {7,1,14};
Line Loop(2) = {6,12,11,-2};
Plane Surface(2) = {2};
Transfinite Surface(2) = {3,7,14,4};

// bottom inner annulus
Symmetry{0,1,0,0}{ Duplicata{Surface{2};} }
Transfinite Surface(13) = {5,9,6,2};
Reverse Surface{13};

// central square
Line Loop(3) = {1,2,3,4};
Plane Surface(3) = {3};
Transfinite Surface(3) = {2,3,4,5};

// right outer annulus
Point(newp) = {r2*Cos(-Pi/4), r2*Sin(-Pi/4), 0};
Point(newp) = {r2*Cos(+Pi/4), r2*Sin(+Pi/4), 0};
Line(newl) = {6,15};
Line(newl) = {7,16};
Circle(newl) = {15,1,16};

Line Loop(4) = {16,18,-17,-7};
Plane Surface(4) = {4};
Transfinite Surface(4) = {6,15,16,7};

// left outer annulus
Symmetry{1,0,0,0}{ Duplicata{Surface{4};} }
Transfinite Surface(19) = {14,23,18,9};
Reverse Surface{19};

// top outer annulus
Circle(newl) = {16,1,23};
Line Loop(5) = {17,23,22,-12};
Plane Surface(5) = {5};
Transfinite Surface(5) = {7,16,23,14};

// bottom outer annulus
Symmetry{0,1,0,0}{ Duplicata{Surface{5};} }
Transfinite Surface(24) = {18,15,6,9};
Reverse Surface{24};

Transfinite Line{1,2,3,4,7,12,10,15,18,23,21,26} = 4;
Transfinite Line{5,6,11,9} = 4;
Transfinite Line{16,17,22,20} = 4;

Physical Line(0) = {18,23,21,26}; // outer boundary
Physical Surface(1) = {1,2,3,8,13,4,5,19,24}; // domain

Mesh.RecombineAll=1;
Geometry.Normals = 100;
