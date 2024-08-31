//------------------------------------------------------------------------------
// Linear advection with velocity = (-y,x) in [-1,1] x [-1,1].
// Time period = 2*pi
//------------------------------------------------------------------------------

#ifndef __TEST_DATA_H__
#define __TEST_DATA_H__

#define XMIN       (-1.0)
#define XMAX       (+1.0)
#define YMIN       (-1.0)
#define YMAX       (+1.0)
#define FINAL_TIME (2.0*M_PI)

using namespace dealii;

//------------------------------------------------------------------------------
// 2d velocity field
//------------------------------------------------------------------------------
void velocity(const Point<2> & p, Tensor<1, 2> &v)
{
   v[0] = -p[1];
   v[1] =  p[0];
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
   const double alpha = 50.0;
   const double x0 = 0.5;
   const double y0 = 0.0;
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
   double x = p[0];
   double y = p[1];
   double r = p.norm();
   double theta = atan2(y, x) - time;
   x = r * cos(theta);
   y = r * sin(theta);

   x -= x0;
   y -= y0;

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
   double x = p[0];
   double y = p[1];
   double r = p.norm();
   double theta = atan2(y, x) - time;
   x = r * cos(theta);
   y = r * sin(theta);

   x -= x0;
   y -= y0;

   double f = exp(-alpha*(x*x + y*y));

   Tensor<1,dim> gp, g;
   gp[0] = f * (-2 * alpha * x);
   gp[1] = f * (-2 * alpha * y);

   g[0] = cos(time) * gp[0] - sin(time) * gp[1];
   g[1] = sin(time) * gp[0] + cos(time) * gp[1];
   return g;
}

#endif
