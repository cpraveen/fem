//------------------------------------------------------------------------------
// Linear advection with velocity = (1,1) in [-1,1] x [-1,1] and periodic bc.
// Time period = 2
//------------------------------------------------------------------------------

#ifndef __TEST_DATA_H__
#define __TEST_DATA_H__

#define XMIN       (-1.0)
#define XMAX       (+1.0)
#define YMIN       (-1.0)
#define YMAX       (+1.0)
#define FINAL_TIME (2.0)

using namespace dealii;

//------------------------------------------------------------------------------
// 2d velocity field
//------------------------------------------------------------------------------
void velocity(const Point<2> & /*p*/, Tensor<1, 2> &v)
{
   v[0] = 1.0;
   v[1] = 1.0;
}

//------------------------------------------------------------------------------
// Bring x into [a,b] by periodicity
//------------------------------------------------------------------------------
void periodic(const double a, const double b, double& x)
{
   const double L = b - a;
   while(x < a) x += L;
   while(x > b) x -= L;
}

//------------------------------------------------------------------------------
// Exact solution
//------------------------------------------------------------------------------
template <int dim>
class Solution : public Function<dim>
{
public:
   Solution() = default;
   double value(const Point<dim>&    p,
                const unsigned int  component = 0) const override;
   Tensor<1, dim> gradient(const Point<dim>&    p,
                           const unsigned int  component = 0) const override;
private:
   const double xmin = XMIN;
   const double xmax = XMAX;
   const double ymin = YMIN;
   const double ymax = YMAX;
   const double alpha = 20.0;
};

//------------------------------------------------------------------------------
// Exact solution works correctly only for periodic case
//------------------------------------------------------------------------------
template<int dim>
double
Solution<dim>::value(const Point<dim>&    p,
                     const unsigned int) const
{
   double time = this->get_time();
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
   double time = this->get_time();
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
