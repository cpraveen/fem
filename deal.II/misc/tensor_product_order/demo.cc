/*
 Check the order of shape functions in FE_DGQLegendre and FE_DGQArbitraryNodes
 */
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/point.h>
#include <deal.II/fe/fe_dgq.h>

#include <iostream>
#include <cmath>
#include <vector>

using namespace dealii;

void order_fe_dgqarbitrary(const int degree)
{
   const int n = degree + 1;
   const QGauss<1> quad1d(n);
   const FE_DGQArbitraryNodes<1> fe1(quad1d);
   const FE_DGQArbitraryNodes<3> fe3(quad1d);
   const QGauss<3> quad3d(n);
   const std::vector<Point<3>>& points3d = quad3d.get_points();

   int index3d = 0;
   for(int k=0; k<n; k++)
      for(int j=0; j<n; j++)
         for(int i=0; i<n; i++)
         {
            //int index3d = i + j*n + k*n*n;
            // compute error between 'index3d'-th polynomial of 'fe3' and the
            // product of i-th, j-th and k-th polynomials of 'fe1'
            double error = 0.0;
            for(unsigned int q=0; q<quad3d.size(); q++)
            {
               double shape3 = fe3.shape_value(index3d, points3d[q]);
               double shape1 = (fe1.shape_value(i, Point<1>(points3d[q][0]))*
                                fe1.shape_value(j, Point<1>(points3d[q][1]))*
                                fe1.shape_value(k, Point<1>(points3d[q][2])));
               error = fmax(error, fabs(shape1-shape3));
            }
            if(error < 1.0e-14)
            std::cout << "Tensor indices: " << i << " " << j << " "
                      << k << ", " << "3d index: " << index3d << ", error: "
                      << error << "\n";
            ++index3d;
         }
}

void order_fe_dgqlegendre(const int degree)
{
   const int n = degree + 1;
   const FE_DGQLegendre<1> fe1(degree);
   const FE_DGQLegendre<3> fe3(degree);
   const QGauss<3> quad3d(n);
   const std::vector<Point<3>>& points3d = quad3d.get_points();

   int index3d = 0;
   for(int k=0; k<n; k++)
      for(int j=0; j<n; j++)
         for(int i=0; i<n; i++)
         {
            //int index3d = i + j*n + k*n*n;
            // compute error between 'index3d'-th polynomial of 'fe3' and the
            // product of i-th, j-th and k-th polynomials of 'fe1'
            double error = 0.0;
            for(unsigned int q=0; q<quad3d.size(); q++)
            {
               double shape3 = fe3.shape_value(index3d, points3d[q]);
               double shape1 = (fe1.shape_value(i, Point<1>(points3d[q][0]))*
                                fe1.shape_value(j, Point<1>(points3d[q][1]))*
                                fe1.shape_value(k, Point<1>(points3d[q][2])));
               error = fmax(error, fabs(shape1-shape3));
            }
            if(error < 1.0e-14)
            std::cout << "Tensor indices: " << i << " " << j << " "
                      << k << ", " << "3d index: " << index3d << ", error: "
                      << error << "\n";
            ++index3d;
         }
}

int main()
{
   for(int i=1; i<5; ++i)
   {
      std::cout << "----- DGQLegendre: degree = " << i << " -----" << std::endl;
      order_fe_dgqlegendre(i);
   }

   for(int i=1; i<5; ++i)
   {
      std::cout << "----- DGQArbitraryNodes: degree = " << i << " -----" << std::endl;
      order_fe_dgqarbitrary(i);
   }
   return 0;
}
