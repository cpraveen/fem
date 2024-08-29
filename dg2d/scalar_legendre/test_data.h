//------------------------------------------------------------------------------
// Linear advection with velocity = (1,1) in [-1,1] x [-1,1] and periodic bc.
//------------------------------------------------------------------------------

#ifndef __TEST_DATA_H__
#define __TEST_DATA_H__

using namespace dealii;

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
   x = std::fmod(x - xmin, xmax - xmin) + xmin;
   y = std::fmod(y - ymin, ymax - ymin) + ymin;

   double value = exp(-alpha*(x*x + y*y));
   return value;
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
   x = std::fmod(x - xmin, xmax - xmin) + xmin;
   y = std::fmod(y - ymin, ymax - ymin) + ymin;

   Tensor<1,dim> value;
   value[0] = f * (-2 * alpha * x);
   value[1] = f * (-2 * alpha * y);
   return value;
}

#endif
