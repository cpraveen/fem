#ifndef __JACOBI_SOLVER_H__
#define __JACOBI_SOLVER_H__

#include "sparse_matrix.h"
#include "Vector.h"

template <class T>
class JacobiSolver
{
   public:
      JacobiSolver (unsigned int max_iter,
                    T            tol);
      ~JacobiSolver () {};
      unsigned int solve (const SparseMatrix<T>& A,
                                      Vector<T>& x, 
                          const       Vector<T>& f) const;

   private:
      unsigned int max_iter;
      T            tol;
};

#endif
