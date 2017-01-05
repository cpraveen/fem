/*
 * Solve 
 *      -u'' = (2*pi)^2 sin(2*pi*x) in (0,1)
 *      u(0) = u0, u(1) = u1
 * Exact solution
 *       u(x) = sin(2*pi*x) + (1-x)*u0 + x*u1
 */
#include <iostream>
#include <fstream>
#include <cmath>

#include "Vector.h"
#include "sparse_matrix.h"
#include "cg_solver.h"
#include "jacobi_solver.h"
#include "sor_solver.h"
#include "ssor_solver.h"

using namespace std;

//------------------------------------------------------------------------------
// Problem definition and exact solution
//------------------------------------------------------------------------------
double u0 = 0, u1 = 1;
double xmin = 0, xmax = 1;
double exact_solution(const double x)
{
   return sin(2*M_PI*x) + (1-x)*u0 + x*u1;
}

//------------------------------------------------------------------------------
// Main program
//------------------------------------------------------------------------------
int main ()
{
   const unsigned int n = 10;
   const double h = (xmax - xmin) / (n - 1);

   // Create 1-d mesh
   Vector<double> x(n);
   for(unsigned int i=0; i<n; ++i)
      x(i) = xmin + i*h;

   const double a =  2.0/(h*h);
   const double b = -1.0/(h*h);

   // Construct matrix
   SparseMatrix<double> A(n);
   A.set(0, 0, a); 
   for(unsigned int i=1; i<n-1; ++i)
   {
      A.set(i, i  , a); // diagonal element should be first
      A.set(i, i-1, b);
      A.set(i, i+1, b);
   }
   A.set(n-1, n-1, a);
   A.close();

   // Construct right hand side vector
   Vector<double> f(n);
   for(unsigned int i=0; i<n; ++i)
      f(i) = pow(2*M_PI,2) * sin(2*M_PI*x(i));

   // Apply boundary conditions
   f(0)   = A(0,0) * u0;
   f(1)  -= A(1,0) * u0;
   f(n-2)-= A(n-2,n-1) * u1;
   f(n-1) = A(n-1,n-1) * u1;

   A(1,0)     = 0;
   A(n-2,n-1) = 0;

   cout << A << endl;
   cout << f << endl;

   // Solution vector
   Vector<double> u(n);
   u      = 0;
   u(0)   = u0;
   u(n-1) = u1;

   unsigned int max_iter = 1000;
   double tol = 1.0e-6;
   //JacobiSolver<double> solver (max_iter, tol);
   //SORSolver<double> solver (max_iter, tol);
   SSORSolver<double> solver (max_iter, tol);
   //CGSolver<double> solver (max_iter, tol);
   unsigned int iter = solver.solve (A, u, f, 1.5);
   //unsigned int iter = solver.solve (A, u, f);

   cout << "Number of iterations = " << iter << endl;

   // Save solution to file
   ofstream fsol("sol.dat");
   for(unsigned int i=0; i<n; ++i)
      fsol << x(i) << "  " 
           << u(i) << "  " 
           << exact_solution(x(i)) << endl;
   fsol.close ();
}
