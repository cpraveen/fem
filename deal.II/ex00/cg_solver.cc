#include <cassert>
#include <cmath>

#include "cg_solver.h"
#include "math_functions.h"

//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------
template <class T>
CGSolver<T>::CGSolver (unsigned int max_iter,
                       T            tol)
   :
      max_iter (max_iter),
      tol (tol)
{
   assert (max_iter > 0);
   assert (tol > 0);
}

//-----------------------------------------------------------------------------
// Solves A*x = f for x
// We assume that x has already been initialized.
//-----------------------------------------------------------------------------
template <class T>
unsigned int CGSolver<T>::solve (const SparseMatrix<T>& A,
                                             Vector<T>& x,
                                 const       Vector<T>& f) const
{
   const unsigned int n = x.size();
   assert (n == A.size());
   assert (n == f.size());

   Vector<T> r(n), v(n);

   // initial residual: r = f - A*x
   A.multiply(x, r, -1); // r = -A*x
   r += f;               // r = r + f

   // initial direction
   Vector<T> d (r);

   std::vector<T> r2 (max_iter);
   r2[0] = dot<T>(r, r);

   unsigned int iter = 0;

   while ( sqrt(r2[iter]/r2[0]) > tol && iter < max_iter)
   {
      if(iter >= 1) // update descent direction
      {             // d = r + beta * d
         const T beta = r2[iter] / r2[iter-1];
         d   *= beta;
         d   += r;
      }

      // v = A*d
      A.multiply (d, v);
      const T omega = r2[iter] / dot<T>(d, v);

      // update x: x = x + omega * d
      x.add (omega, d);

      // update residual: r = r - omega * v
      r.add (-omega, v);

      ++iter;

      r2[iter] = dot<T> (r, r);
   }

   if ( sqrt(r2[iter]/r2[0]) > tol && iter == max_iter)
   {
      T r0 = sqrt(r2[0]);
      T r1 = sqrt(r2[iter]);
      std::cout << "CGSolver did not converge !!!\n";
      std::cout << "   No. of iterations= " << iter << std::endl;
      std::cout << "   Initial residual = " << r0 << std::endl;
      std::cout << "   Final   residual = " << r1 << std::endl;
      std::cout << "   Final/Initial    = " << r1/r0 << std::endl;
   }

   return iter;
}

//-----------------------------------------------------------------------------
// Instantiation
//-----------------------------------------------------------------------------
template class CGSolver<float>;
template class CGSolver<double>;
