#include <cassert>
#include <cmath>

#include "ssor_solver.h"
#include "math_functions.h"

//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------
template <class T>
SSORSolver<T>::SSORSolver (unsigned int max_iter,
                           T            tol,
                           T            omg)
   :
      max_iter (max_iter),
      tol (tol),
      omg (omg)
{
   assert (max_iter > 0);
   assert (tol > 0);
   assert (omg >= 1.0 && omg <= 2.0);
   std::cout << "Using SSOR solver\n";
}

//-----------------------------------------------------------------------------
// Solves A*x = f for x
// We assume that x has already been initialized.
//-----------------------------------------------------------------------------
template <class T>
unsigned int SSORSolver<T>::solve (const SparseMatrix<T>& A,
                                               Vector<T>& x,
                                   const       Vector<T>& f) const
{
   const unsigned int n = x.size();
   assert (n == A.size());
   assert (n == f.size());

   T res, res0;
   res = res0 = 1;
   unsigned int iter = 0;

   while (res/res0 > tol && iter < max_iter)
   {
      // forward loop
      for(unsigned int i=0; i<n; ++i)
      {
         T r = f(i) - A.multiply(x, i);
         x(i) += omg * r / A.diag(i);
      }
      
      // backward loop
      res = 0;
      for(int i=n-1; i>=0; --i)
      {
         T r = f(i) - A.multiply(x, i);
         x(i) += omg * r / A.diag(i);
         res  += r * r;
      }
      
      res = std::sqrt( res );
      if(iter==0) res0 = res;
      ++iter;
   }

   if (res/res0 > tol)
   {
      std::cout << "SSORSolver did not converge !!!\n";
      std::cout << "   No. of iterations= " << iter << std::endl;
      std::cout << "   Initial residual = " << res0 << std::endl;
      std::cout << "   Final   residual = " << res  << std::endl;
      std::cout << "   Final/Initial    = " << res/res0 << std::endl;
   }

   return iter;
}

//-----------------------------------------------------------------------------
// Instantiation
//-----------------------------------------------------------------------------
template class SSORSolver<float>;
template class SSORSolver<double>;
