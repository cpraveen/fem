#ifndef __SSOR_SOLVER_H__
#define __SSOR_SOLVER_H__

#include "sparse_matrix.h"
#include "Vector.h"

template <class T>
class SSORSolver
{
   public:
      SSORSolver (unsigned int max_iter,
                  T            tol);
      ~SSORSolver () {};
      unsigned int solve (const SparseMatrix<T>& A,
                                      Vector<T>& x, 
                          const       Vector<T>& f,
                          const               T  omg = 1) const;

   private:
      unsigned int max_iter;
      T            tol;
};

#endif
