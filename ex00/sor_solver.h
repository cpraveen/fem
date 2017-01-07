#ifndef __SOR_SOLVER_H__
#define __SOR_SOLVER_H__

#include "sparse_matrix.h"
#include "Vector.h"

template <class T>
class SORSolver
{
   public:
      SORSolver (unsigned int max_iter,
                 T            tol,
                 T            omg);
      ~SORSolver () {};
      unsigned int solve (const SparseMatrix<T>& A,
                                      Vector<T>& x, 
                          const       Vector<T>& f) const;

   private:
      unsigned int max_iter;
      T            tol;
      T            omg;
};

#endif
