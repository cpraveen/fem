//------------------------------------------------------------------------------
namespace ProblemData
{
   const std::string name = "2D Riemann problem";
   const double xmin = 0.0;
   const double xmax = 1.0;
   const double ymin = 0.0;
   const double ymax = 1.0;
   const double final_time = 0.25;
   const bool periodic_x = false;
   const bool periodic_y = false;

   const double gamma = 1.4;
}

//------------------------------------------------------------------------------
template <int dim>
struct Problem : ProblemBase<dim>
{
   const double gamma = ProblemData::gamma;

   //---------------------------------------------------------------------------
   void initial_value(const Point<dim>& /*p*/,
                      Vector<double>&   u) const override
   {
      // TODO: Implement 4 states
      const double rho = 0.5313;
      const double pre = 0.4;
      Tensor<1,dim> vel;
      vel[0] = 0.0;
      vel[1] = 0.0;

      PDE::prim2con(rho, vel, pre, u);
   }

   //---------------------------------------------------------------------------
   void boundary_value(const int             /*boundary_id*/,
                       const Point<dim>&     /*p*/,
                       const double          /*t*/,
                       const Tensor<1,dim>&  /*normal*/,
                       const Vector<double>& Uint,
                       Vector<double>&       Uout) const override
   {
      // All boundaries are Neumann
      Uout = Uint;
   }

};
