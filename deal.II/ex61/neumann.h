namespace PrescribedSolution
{
   //---------------------------------------------------------------------------
   template <int dim>
   class ExactSolution : public Function<dim>
   {
   public:
      ExactSolution() : Function<dim>(dim+1) {}

      double value(const Point<dim> &p,
                   const unsigned int component) const override;
   };

   template <>
   double ExactSolution<2>::value(const Point<2> &p,
                                  const unsigned int component) const
   {
      if(component == 0)      // p_x
         return -2 * M_PI * sin(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]);
      else if(component == 1) // p_y
         return -2 * M_PI * cos(2 * M_PI * p[0]) * sin(2 * M_PI * p[1]);
      else                    // p
         return cos(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]);
   }

   //---------------------------------------------------------------------------
   template <int dim>
   class RHSFunction : public Function<dim>
   {
   public:
      RHSFunction() : Function<dim>() {}

      double value(const Point<dim> &p,
                   const unsigned int component = 0) const override;
   };

   template <>
   double RHSFunction<2>::value(const Point<2> &p,
                                const unsigned int /*component*/) const
   {
      return 8 * M_PI * M_PI * cos(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]);
   }
} // namespace PrescribedSolution
