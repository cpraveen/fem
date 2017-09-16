#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>

#include "sparse_matrix.h"
#include "Vector.h"
#include "math_functions.h"

//-----------------------------------------------------------------------------
// Constructor, fully initializes matrix.
//-----------------------------------------------------------------------------
template <class T>
SparseMatrix<T>::SparseMatrix (std::vector<unsigned int>& row_ptr, 
                               std::vector<unsigned int>& col_ind, 
                                          std::vector<T>& val)
   :
   nrow (row_ptr.size()-1),
   row_ptr (row_ptr),
   col_ind (col_ind),
   val (val)
{
   assert (row_ptr.size() >= 2);
   assert (col_ind.size() > 0);
   assert (col_ind.size() == val.size());
   assert (row_ptr[nrow] == val.size());

   for(unsigned int i=0; i<col_ind.size(); ++i)
      assert (col_ind[i] >= 0 && col_ind[i] < nrow);
   
   state = CLOSED;
}

//------------------------------------------------------------------------------
// Constructor, only size if set, values are not set.
//------------------------------------------------------------------------------
template <class T>
SparseMatrix<T>::SparseMatrix (unsigned int nrow)
   :
      nrow (nrow)
{
   assert (nrow > 0);
   row_ptr.resize (nrow+1, 0);
   cols.resize(nrow);
   vals.resize(nrow);
   state = OPEN;
}

//------------------------------------------------------------------------------
// Create sparsity pattern row by row
// Insert "value" into location (i,j)
//------------------------------------------------------------------------------
template <class T>
void SparseMatrix<T>::set (const unsigned int i,
                           const unsigned int j,
                           const T            value)
{
   assert (state == OPEN);
   assert (i >= 0 && i < nrow);
   assert (j >= 0 && j < nrow);

   cols[i].push_back(j);
   vals[i].push_back(value);
}

//------------------------------------------------------------------------------
// Close sparsity pattern. After this point, sparsity pattern cannot be changed.
//------------------------------------------------------------------------------
template <class T>
void SparseMatrix<T>::close ()
{
   unsigned int nval = 0;
   for(unsigned int i=0; i<nrow; ++i)
   {
      nval += cols[i].size();
      assert(cols[i].size() > 0); // i'th row is empty
      assert(i == cols[i][0]); // set diagonal first
   }
   col_ind.resize(nval);
   val.resize(nval);

   unsigned int c = 0;
   row_ptr[0] = 0;
   for(unsigned int i=0; i<nrow; ++i)
   {
      for(unsigned int j=0; j<cols[i].size(); ++j)
      {
         col_ind[c] = cols[i][j];
         val[c] = vals[i][j];
         ++c;
      }
      row_ptr[i+1] = c;
   }

   // Free up memory
   for(unsigned int i=0; i<nrow; ++i)
   {
      cols[i].resize(0);
      vals[i].resize(0);
   }
   cols.resize(0);
   vals.resize(0);

   state = CLOSED;
}

//-----------------------------------------------------------------------------
// Zero all elements in i'th row and i'th column except diagonal
//-----------------------------------------------------------------------------
template <class T>
void SparseMatrix<T>::zero_off_diag(const unsigned int i)
{
   unsigned int row_beg = row_ptr[i];
   unsigned int row_end = row_ptr[i+1];
   for(unsigned int d=row_beg+1; d<row_end; ++d)
   {
      val[d] = 0; // a(i,j)
      (*this)(col_ind[d],i) = 0; // a(j,i)
   }
}

//-----------------------------------------------------------------------------
// Get element value of A(i,j)
//-----------------------------------------------------------------------------
template <class T>
T SparseMatrix<T>::operator()(unsigned int i, 
                              unsigned int j) const
{
   unsigned int row_beg = row_ptr[i];
   unsigned int row_end = row_ptr[i+1];
   for(unsigned int d=row_beg; d<row_end; ++d)
      if(col_ind[d] == j) 
         return val[d];
   return 0;
}

//-----------------------------------------------------------------------------
// Get reference to element A(i,j)
//-----------------------------------------------------------------------------
template <class T>
T& SparseMatrix<T>::operator()(unsigned int i, 
                               unsigned int j)
{
   unsigned int row_beg = row_ptr[i];
   unsigned int row_end = row_ptr[i+1];
   for(unsigned int d=row_beg; d<row_end; ++d)
      if(col_ind[d] == j) 
         return val[d];
   std::cout << "Element " << i << ", " << j << " does not exist\n";
   abort ();
}

//-----------------------------------------------------------------------------
// y = scalar * A * x
//-----------------------------------------------------------------------------
template <class T>
void SparseMatrix<T>::multiply(const Vector<T>& x, 
                                     Vector<T>& y,
                               const T          scalar) const
{
   assert (x.size() == nrow);
   assert (x.size() == y.size());
   for(unsigned int i=0; i<nrow; ++i)
   {
      y(i) = 0;
      unsigned int row_beg = row_ptr[i];
      unsigned int row_end = row_ptr[i+1];
      for(unsigned int j=row_beg; j<row_end; ++j)
         y(i) += val[j] * x(col_ind[j]);
      y(i) *= scalar;
   }
}

//-----------------------------------------------------------------------------
// Computes only i'th component of A*x, i.e., dot product between
// i'th row of A and vector x.
//-----------------------------------------------------------------------------
template <class T>
T SparseMatrix<T>::multiply(const Vector<T>&   x,
                            const unsigned int i) const
{
   assert (x.size() == nrow);
   
   unsigned int row_beg = row_ptr[i];
   unsigned int row_end = row_ptr[i+1];
   T r = 0;
   for(unsigned int j=row_beg; j<row_end; ++j)
      r += val[j] * x(col_ind[j]);
   return r;
}

//-----------------------------------------------------------------------------
// Return i'th diagonal 
//-----------------------------------------------------------------------------
template <class T>
T SparseMatrix<T>::diag (const unsigned int i) const
{
   return val[ row_ptr[i] ];
}

//-----------------------------------------------------------------------------
// Instantiation
//-----------------------------------------------------------------------------
template class SparseMatrix<int>;
template class SparseMatrix<float>;
template class SparseMatrix<double>;
