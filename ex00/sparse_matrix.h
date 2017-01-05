#ifndef __SPARSE_MATRIX_H__
#define __SPARSE_MATRIX_H__

#include <iomanip>
#include <vector>

#include "Vector.h"

enum MatrixState { OPEN, CLOSED };

template <class T>
class SparseMatrix
{
   public:
      SparseMatrix (std::vector<unsigned int>& row_ptr, 
                    std::vector<unsigned int>& col_ind, 
                               std::vector<T>& val);
      SparseMatrix (unsigned int nrow);
      ~SparseMatrix() {};
      unsigned int size () const 
      {
         return nrow;
      }
      void set (const unsigned int i, 
                const unsigned int j, 
                const T            value);
      void close ();
      void multiply(const Vector<T>& x, 
                          Vector<T>& y,
                    const T          scalar = 1) const;
      T operator()(unsigned int i, 
                   unsigned int j) const;
      T& operator()(unsigned int i, 
                    unsigned int j);
      T diag (unsigned int i) const;
      T Jacobi_step(Vector<T>&       x, 
                    const Vector<T>& rhs) const;
      T SOR_step(Vector<T>&       x, 
                 const Vector<T>& rhs,
                 const T          omg) const;
      T SSOR_step(Vector<T>&       x, 
                  const Vector<T>& rhs,
                  const T          omg) const;
      friend std::ostream& operator<< (std::ostream&          os, 
                                       const SparseMatrix<T>& A)
      {
         for(unsigned int i=0; i<A.size(); ++i)
         {
            for(unsigned int j=0; j<A.size(); ++j)
               os << std::setw(10) << A(i,j);
            os << std::endl;
         }
         return os;
      }
      
   private:
      unsigned int nrow;
      std::vector<unsigned int> row_ptr, col_ind;
      std::vector<T> val;
      MatrixState state;
};

#endif
