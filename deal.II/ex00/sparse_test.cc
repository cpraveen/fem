#include <iostream>

#include "sparse_matrix.h"
#include "Vector.h"

using namespace std;

int main ()
{
   unsigned int nrow=4, nval=8;
   vector<unsigned int> row_ptr(nrow+1), col_ind(nval);
   vector<double> val(nval);

   row_ptr[0] = 0; 
   row_ptr[1] = 2;
   row_ptr[2] = 4;
   row_ptr[3] = 7;
   row_ptr[4] = 8;

   col_ind[0] = 0;
   col_ind[1] = 3;
   col_ind[2] = 1;
   col_ind[3] = 2;
   col_ind[4] = 0;
   col_ind[5] = 2;
   col_ind[6] = 3;
   col_ind[7] = 3;

   val[0] = 10;
   val[1] = 7;
   val[2] = 1;
   val[3] = 2;
   val[4] = 3;
   val[5] = 5;
   val[6] = 9;
   val[7] = 1;

   SparseMatrix<double> A(row_ptr, col_ind, val);
   cout << A << endl;
   cout << "A(2,3) = " << A(2,3) << endl << endl;
   
   Vector<double> x(nrow), y(nrow);

   x = 1.0;
   cout << "x = \n" << x << endl;

   A.multiply(x, y);
   cout << "y = \n" << y;
}
