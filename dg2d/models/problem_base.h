//------------------------------------------------------------------------------
// This should be set by user in a problem.h file
//------------------------------------------------------------------------------
namespace ProblemData
{
   extern const std::string name;
   extern const double xmin;
   extern const double xmax;
   extern const double ymin;
   extern const double ymax;
   extern const double final_time;
   extern const bool periodic_x;
   extern const bool periodic_y;
}

//------------------------------------------------------------------------------
template <int dim>
struct ProblemBase
{
   ProblemBase()
   {}

   virtual void make_grid(Triangulation<dim>& /*triangulation*/) const
   {
      AssertThrow(false, ExcNotImplemented());
   }

   // Transform grid but cells must remain rectangular if you are using FE_DGP. 
   // E.g., you can do some local grid refinement.
   virtual void transform_grid(Triangulation<dim>& /*triangulation*/) const
   {
      // can be empty
   }

   // Attach manifolds to define curved boundaries
   virtual void set_manifolds(Triangulation<dim>& /*triangulation*/) const
   {
      // can be empty
   }

   virtual void initial_value(const Point<dim>& p,
                              Vector<double>&   u) const = 0;

   virtual void boundary_value(const int             /*boundary_id*/,
                               const Point<dim>&     /*p*/,
                               const double          /*t*/,
                               const Tensor<1,dim>&  /*normal*/,
                               const Vector<double>& /*Uint*/,
                               Vector<double>&       /*Uout*/) const
   {
      AssertThrow(false, ExcNotImplemented());
   }

   virtual std::string get_name()
   {
      return ProblemData::name;
   }

   virtual double get_xmin()
   {
      return ProblemData::xmin;
   }

   virtual double get_xmax()
   {
      return ProblemData::xmax;
   }

   virtual double get_ymin()
   {
      return ProblemData::ymin;
   }

   virtual double get_ymax()
   {
      return ProblemData::ymax;
   }

   virtual bool get_periodic_x()
   {
      return ProblemData::periodic_x;
   }

   virtual bool get_periodic_y()
   {
      return ProblemData::periodic_y;
   }

   virtual bool get_periodic()
   {
      return (ProblemData::periodic_x && ProblemData::periodic_y);
   }

   virtual double get_final_time()
   {
      return ProblemData::final_time;
   }
};
