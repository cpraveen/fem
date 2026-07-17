// Generate quads
Mesh.RecombineAll = 1;

// Comment both the above lines to get triangular prisms.
// NOTE: Make sure you dont get any tetrahedra.

h = -0.1;
nlayers = 1;

n = 100;     // points on upper surface of airfoil
m = 2*n - 2; // total number of points on airfoil without repetition
             // LE and TE points are common to upper/lower surface
nc = 25;     // points on each quarter of intermediate circle
R = 50;      // Radius of outer boundary
r = 2.0;

cl1 = 1.0/400;
cl2 = 1.0/50;
cl3 = 5.0; // length scale on outer boundary

naf = 100; // no. of points on airfoil (each side)
raf = 0.1; // clustering ratio on airfoil

// NACA0012 profile
// formula taken from http://turbmodels.larc.nasa.gov/naca0012_val.html
Macro NACA0012
   x2 = x * x;
   x3 = x * x2;
   x4 = x * x3;
   y = 0.594689181*(0.298222773*Sqrt(x) 
       - 0.127125232*x - 0.357907906*x2 + 0.291984971*x3 - 0.105174606*x4);
Return

// put points on upper surface of airfoil
For i In {1:n}
   theta = Pi * (i-1) / (n-1);
   x = 0.5 * (Cos(theta) + 1.0);
   Call NACA0012;
   cl[i] = cl1 + 4*(cl2-cl1)*x*(1-x);
   Point(i) = {x, y, 0.0};
   xx[i] = x;
   yy[i] = y;
EndFor

// put points on lower surface of airfoil, use upper surface points and reflect
For i In {n+1:m}
   Point(i) = {xx[2*n-i], -yy[2*n-i], 0.0};
EndFor

Spline(1) = {1:n};
Spline(2) = {n:m, 1};

// Points on outer boundary
Point(1000) = {R+0.5, 0, 0, cl3};
Point(1001) = {0.5, R, 0, cl3};
Point(1002) = {-R+0.5, 0, 0, cl3};
Point(1003) = {0.5, -R, 0, cl3};

Point(1004) = {0.5, 0, 0};

Point(1005) = {0.5+r, 0, 0};
Point(1006) = {0.5, r, 0};
Point(1007) = {0.5-r, 0, 0};
Point(1008) = {0.5, -r, 0};

Circle(3) = {1000, 1004, 1001};
Circle(4) = {1001, 1004, 1002};
Circle(5) = {1002, 1004, 1003};
Circle(6) = {1003, 1004, 1000};

Circle(7) = {1005, 1004, 1006};
Circle(8) = {1006, 1004, 1007};
Circle(9) = {1007, 1004, 1008};
Circle(10) = {1008, 1004, 1005};

Line Loop(1) = {1,2};
Line Loop(2) = {3,4,5,6};
Line Loop(3) = {7,8,9,10};

Transfinite Line{1,2} = naf Using Bump raf; // airfoil
Transfinite Line{7,8,9,10} = nc;            // first circular boundary

Plane Surface(1) = {1,3};
Plane Surface(2) = {2,3};

Physical Surface(100) = {1,2};
Physical Line(0) = {1,2};     // airfoil
Physical Line(1) = {3,4,5,6}; // outer boundary
