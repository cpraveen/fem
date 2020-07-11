r = 1.0;  // cylinder radius
R = 2.0; // outer boundary radius

n1 = 20; // number of points along radius
n2 = 20;  // number of points along angular direction in (0,pi/2)
p  = 1.0; // progression ration along radius

Point(1) = { 0, 0, 0};
Point(2) = {-r, 0, 0};
Point(3) = { 0, r, 0};

Point(4) = {-R, 0, 0};
Point(5) = { 0, R, 0};

Circle(1) = {2,1,3};
Circle(2) = {4,1,5};
Line(3)   = {3,5};
Line(4)   = {4,2};

// Make one quadrant
Line Loop(1) = {1,3,-2,4};
Plane Surface(1) = {1};
Transfinite Surface(1) = {2,3,5,4};
Recombine Surface(1);

Transfinite Line {3,-4} = n1 Using Progression p;
Transfinite Line {1,2} = n2;

// Transform above quadrant to three others
Symmetry{1,0,0,0}{ Duplicata{Surface{1};} }
Transfinite Surface(5) = {6,17,5,3};
Recombine Surface(5);

Symmetry{0,1,0,0}{ Duplicata{Surface{1};} }
Transfinite Surface(10) = {4,24,20,2};
Recombine Surface(10);

Symmetry{0,1,0,0}{ Duplicata{Surface{5};} }
Transfinite Surface(14) = {6,20,24,17};
Recombine Surface(14);

Transfinite Line {-9,12} = n1 Using Progression p;
Transfinite Line {6,11,15,17,13,8} = n2;

// Inner circular disc
Line Loop(2) = {1,-6,15,-11};
Plane Surface(2) = {2};
Recombine Surface(2);

Physical Line(1) = {1,6,11,15}; // cylinder surface
Physical Line(2) = {2,8,17,13}; // outer boundary
Physical Surface(3) = {1,5,14,10}; 

//Geometry.Normals = 100;
//Mesh.MshFileVersion = 2.1;
