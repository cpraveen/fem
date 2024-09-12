using namespace dealii;

//------------------------------------------------------------------------------
// Initial condition
//------------------------------------------------------------------------------
template <int dim>
class InitialCondition : public Function<dim>
{
public:
   InitialCondition() : Function<dim>(nvar) 
   {
      xmin = 0.0;
      xmax = 1.0;

      x0 = 0.75;
      beta = 100.0;
   }
   void vector_value(const Point<dim>& p,
                     Vector<double>&   values) const override;
   double xmin, xmax;
   double x0, beta;
};

//------------------------------------------------------------------------------
// Initial condition
//------------------------------------------------------------------------------
template<int dim>
void
InitialCondition<dim>::vector_value(const Point<dim>& p,
                                    Vector<double>&   values) const
{
   double x = p[0];
   values[0] = exp(-beta*pow(x-x0,2));
   values[1] = 0.0;
}
