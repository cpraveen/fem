#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>

#include "sparse_matrix.h"
#include "Vector.h"
#include "math_functions.h"

//-----------------------------------------------------------------------------
// Constructor
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
}

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template <class T>
SparseMatrix<T>::SparseMatrix (unsigned int nrow)
   :
      nrow (nrow)
{
   assert (nrow > 0);
   row_ptr.resize (nrow+1, 0);
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

   static unsigned int count = 0;
   
   if(row_ptr[i+1] == 0) 
      row_ptr[i+1] = count;
   
   ++count;
   ++row_ptr[i+1];
   col_ind.push_back (j);
   val.push_back (value);
}

//------------------------------------------------------------------------------
// Close sparsity pattern. After this point, sparsity pattern cannot be changed.
//------------------------------------------------------------------------------
template <class T>
void SparseMatrix<T>::close ()
{
   // Take care of empty row
   for(unsigned int i=1; i<nrow+1; ++i)
   {
      if(row_ptr[i] == 0) 
      {
         row_ptr[i] = row_ptr[i-1];
         std::cout << "Warning: Row " << i << " is empty\n";
      }
   }
   
   state = CLOSED;
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
// Return i'th diagonal 
//-----------------------------------------------------------------------------
template <class T>
T SparseMatrix<T>::diag (const unsigned int i) const
{
   return val[ row_ptr[i] ];
}

//-----------------------------------------------------------------------------
// Perform one step of Jacobi
//-----------------------------------------------------------------------------
template <class T>
T SparseMatrix<T>::Jacobi_step(Vector<T>&       x, 
                               const Vector<T>& rhs) const
{
   Vector<T> r(nrow);
   multiply(x, r, -1); // r = -A*x
   r += rhs;           // r = r + rhs

   for(unsigned int i=0; i<nrow; ++i)
      x(i) += r(i) / diag(i);

   return std::sqrt( dot<T>(r,r) );
}

//-----------------------------------------------------------------------------
// Perform one step of SOR
//-----------------------------------------------------------------------------
template <class T>
T SparseMatrix<T>::SOR_step(Vector<T>&       x, 
                            const Vector<T>& rhs,
                            const T          omg) const
{
   T res = 0;
   for(unsigned int i=0; i<nrow; ++i)
   {
      T r = rhs(i);
      unsigned int row_beg = row_ptr[i];
      unsigned int row_end = row_ptr[i+1];
      for(unsigned int j=row_beg; j<row_end; ++j)
         r -= val[j] * x(col_ind[j]);
      x(i) += omg * r / diag(i);
      res  += r * r;
   }

   return std::sqrt( res );
}

//-----------------------------------------------------------------------------
// Perform one step of SSOR
//-----------------------------------------------------------------------------
template <class T>
T SparseMatrix<T>::SSOR_step(Vector<T>&       x, 
                             const Vector<T>& rhs,
                             const T          omg) const
{
   // forward loop
   for(unsigned int i=0; i<nrow; ++i)
   {
      T r = rhs(i);
      unsigned int row_beg = row_ptr[i];
      unsigned int row_end = row_ptr[i+1];
      for(unsigned int j=row_beg; j<row_end; ++j)
         r -= val[j] * x(col_ind[j]);
      x(i) += omg * r / diag(i);
   }

   // backward loop
   T res = 0;
   for(int i=nrow-1; i>=0; --i)
   {
      T r = rhs(i);
      unsigned int row_beg = row_ptr[i];
      unsigned int row_end = row_ptr[i+1];
      for(unsigned int j=row_beg; j<row_end; ++j)
         r -= val[j] * x(col_ind[j]);
      x(i) += omg * r / diag(i);
      res  += r * r;
   }

   return std::sqrt( res );
}

//-----------------------------------------------------------------------------
// Instantiation
//-----------------------------------------------------------------------------
template class SparseMatrix<int>;
template class SparseMatrix<float>;
template class SparseMatrix<double>;
