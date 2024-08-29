//------------------------------------------------------------------------------
// Linear advection with velocity = (1,1) in [-1,1] x [-1,1] and periodic bc.
//------------------------------------------------------------------------------

#ifndef __TEST_DATA_H__
#define __TEST_DATA_H__

using namespace dealii;

//------------------------------------------------------------------------------
// Bring x into [xmin,xmax] by periodicity
//------------------------------------------------------------------------------
void periodic(const double xmin, const double xmax, double& x)
{
   const double L = xmax-xmin;
   while(x < xmin) x += L;
   while(x > xmax) x -= L;
}

//------------------------------------------------------------------------------
// Exact solution
//------------------------------------------------------------------------------
template <int dim>
class Solution : public Function<dim>
{
public:
   Solution(double time)
      :
      Function<dim>(),
      time(time)
   {}

   double value(const Point<dim>&    p,
                const unsigned int  component = 0) const override;
   Tensor<1, dim> gradient(const Point<dim>&    p,
                           const unsigned int  component = 0) const override;
private:
   const double xmin = -1.0;
   const double xmax =  1.0;
   const double ymin = -1.0;
   const double ymax =  1.0;
   const double alpha = 20.0;
   const double time;
};

//------------------------------------------------------------------------------
// Exact solution works correctly only for periodic case
//------------------------------------------------------------------------------
template<int dim>
double
Solution<dim>::value(const Point<dim>&    p,
                     const unsigned int) const
{
   double x = p[0] - time;
   double y = p[1] - time;

   // Apply periodicity
   periodic(xmin, xmax, x);
   periodic(ymin, ymax, y);

   double f = exp(-alpha*(x*x + y*y));
   return f;
}

//------------------------------------------------------------------------------
// Exact solution works correctly only for periodic case
//------------------------------------------------------------------------------
template<int dim>
Tensor<1, dim>
Solution<dim>::gradient(const Point<dim>&    p,
                        const unsigned int) const
{
   double f = value(p);

   double x = p[0] - time;
   double y = p[1] - time;

   // Apply periodicity
   periodic(xmin, xmax, x);
   periodic(ymin, ymax, y);

   Tensor<1,dim> g;
   g[0] = f * (-2 * alpha * x);
   g[1] = f * (-2 * alpha * y);
   return g;
}

#endif
