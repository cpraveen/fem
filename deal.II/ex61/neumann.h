namespace PrescribedSolution
{
   //---------------------------------------------------------------------------
   template <int dim>
   class ExactSolution : public Function<dim>
   {
   public:
      ExactSolution() : Function<dim>() {}

      double value(const Point<dim> &p,
                   const unsigned int component = 0) const override;
      Tensor<1, dim> gradient(const Point<dim> &p,
                              const unsigned int component = 0) const override;
   };

   template <>
   double ExactSolution<2>::value(const Point<2> &p,
                                  const unsigned int /*component*/) const
   {
      return cos(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]);
   }

   template <>
   Tensor<1, 2> ExactSolution<2>::gradient(const Point<2> &p,
                                           const unsigned int) const
   {
      Tensor<1, 2> values;
      values[0] = -2 * M_PI * sin(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]);
      values[1] = -2 * M_PI * cos(2 * M_PI * p[0]) * sin(2 * M_PI * p[1]);
      return values;
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
