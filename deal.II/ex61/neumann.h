namespace PrescribedSolution
{
   const double pi = M_PI;

   //---------------------------------------------------------------------------
   template <int dim>
   class ExactSolution : public Function<dim>
   {
   public:
      ExactSolution() : Function<dim>(dim+1) {}

      void vector_value(const Point<dim>& p,
                        Vector<double>&   value) const override;
      void vector_gradient(const Point<dim>&           p,
                           std::vector<Tensor<1,dim>>& value) const override;
   };

   template <>
   void ExactSolution<2>::vector_value(const Point<2>& p,
                                       Vector<double>& value) const
   {
      const double x = p[0];
      const double y = p[1];
      const double cx = cos(2 * pi * x);
      const double sx = sin(2 * pi * x);
      const double cy = cos(2 * pi * y);
      const double sy = sin(2 * pi * y);

      value[0] = -2 * pi * sx * cy;
      value[1] = -2 * pi * cx * sy;
      value[2] = cx * cy;
   }

   template <>
   void ExactSolution<2>::vector_gradient(const Point<2>&           p,
                                          std::vector<Tensor<1,2>>& value) const
   {
      const double x = p[0];
      const double y = p[1];
      const double cx = cos(2 * pi * x);
      const double sx = sin(2 * pi * x);
      const double cy = cos(2 * pi * y);
      const double sy = sin(2 * pi * y);
      const double a  = pow(2 * pi, 2);
      const double b  = 2 * pi;

      value[0][0] = -a * cx * cy;
      value[0][1] =  a * sx * sy;

      value[1][0] =  a * sx * sy;
      value[1][1] = -a * cx * cy;

      value[2][0] = -b * sx * cy;
      value[2][1] = -b * cx * sy;
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
      return 8 * pi * pi * cos(2 * pi * p[0]) * cos(2 * pi * p[1]);
   }
} // namespace PrescribedSolution
