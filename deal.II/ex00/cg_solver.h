#ifndef __CG_SOLVER_H__
#define __CG_SOLVER_H__

#include "sparse_matrix.h"
#include "Vector.h"

template <class T>
class CGSolver
{
   public:
      CGSolver (unsigned int max_iter,
                T            tol);
      ~CGSolver () {};
      unsigned int solve (const SparseMatrix<T>& A,
                                      Vector<T>& x, 
                          const       Vector<T>& f) const;

   private:
      unsigned int max_iter;
      T            tol;
};

#endif
