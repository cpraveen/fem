#include <cassert>
#include <cmath>

#include "sor_solver.h"
#include "math_functions.h"

//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------
template <class T>
SORSolver<T>::SORSolver (unsigned int max_iter,
                         T            tol,
                         T            omg)
   :
      max_iter (max_iter),
      tol (tol),
      omg (omg)
{
   assert (max_iter > 0);
   assert (tol > 0);
   assert (omg >= 1 && omg <= 2);
   std::cout << "Using SOR solver\n";
}

//-----------------------------------------------------------------------------
// Solves A*x = f for x
// We assume that x has already been initialized.
//-----------------------------------------------------------------------------
template <class T>
unsigned int SORSolver<T>::solve (const SparseMatrix<T>& A,
                                              Vector<T>& x,
                                  const       Vector<T>& f) const
{
   const unsigned int n = x.size();
   assert (n == A.size());
   assert (n == f.size());
   assert (omg > 0 && omg < 2);

   T res, res0;
   res = res0 = 1;
   unsigned int iter = 0;

   while (res/res0 > tol && iter < max_iter)
   {      
      res = 0;
      for(unsigned int i=0; i<n; ++i)
      {
         T r = f(i) - A.multiply(x,i);
         x(i) += omg * r / A.diag(i);
         res  += r * r;
      }
      
      res = std::sqrt( res );
      if(iter==0) res0 = res;
      ++iter;
   }

   if (res/res0 > tol)
   {
      std::cout << "SORSolver did not converge !!!\n";
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
template class SORSolver<float>;
template class SORSolver<double>;
