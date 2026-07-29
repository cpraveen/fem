//------------------------------------------------------------------------------
// Euler equations for compressible flows
//------------------------------------------------------------------------------
#ifndef __PDE_H__
#define __PDE_H__

#include <deal.II/base/vectorization.h>
#include <deal.II/numerics/data_postprocessor.h>

using namespace dealii;

constexpr unsigned int nvar = 1;

// Numerical flux functions
enum class FluxType {upwind, none};

std::map<std::string, FluxType> FluxTypeList{{"upwind", FluxType::upwind},
                                             {"none",   FluxType::none}};

//------------------------------------------------------------------------------
template <int dim, typename Number>
struct FluxData
{
   Point<dim, Number> p;                                 // quadrature point
   double t;                                             // time
   Tensor<1, nvar, Number> ul;                           // left  cell average
   Tensor<1, nvar, Number> ur;                           // right cell average
};

//------------------------------------------------------------------------------
// This should be set by user in a problem.h file
//------------------------------------------------------------------------------
namespace ProblemData
{
   template <int dim, typename Number>
   void velocity(const Point<dim, Number>& p, Tensor<1, dim, Number>& v);
}

//------------------------------------------------------------------------------
namespace PDE
{

   const std::string name = "2D linear advection equation";
   using ProblemData::velocity;

   template <int dim, typename Number>
   using FluxMatrix = Tensor<1, dim, Number>;

   template <int dim, typename Number>
   using NormalFlux = Number;

   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE void
   upwind_flux(const Number&                            ul,
               const Number&                            ur,
               const Tensor<1, dim, Number>&            normal,
               const FluxData<dim, Number>&             data,
               NormalFlux<dim, Number>&                 flux)
   {
      Tensor<1, dim, Number> vel;
      velocity(data.p, vel);
      const Number vn = vel * normal;
      const Number u = compare_and_apply_mask<SIMDComparison::greater_than>(
         vn, Number(0.0), ul, ur);
      flux = vn * u;
   }

   //---------------------------------------------------------------------------
   // Following functions are directly called from DG solver
   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE void
   max_speed(const Tensor<1, nvar, Number>& /*u*/,
             const Point<dim, Number>&       p,
             Tensor<1, dim, Number>&         speed)
   {
      velocity(p, speed);
   }

   //---------------------------------------------------------------------------
   // Flux of the PDE model: f(u,x)
   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE void
   physical_flux(const Number&                         u,
                 const FluxData<dim, Number>&          data,
                 FluxMatrix<dim, Number>&              flux)
   {
      Tensor<1, dim, Number> vel;
      velocity(data.p, vel);
      for(unsigned int d = 0; d < dim; ++d)
         flux[d] = vel[d] * u;
   }

   //---------------------------------------------------------------------------
   // Compute flux across cell faces
   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE void
   numerical_flux(const FluxType                            flux_type,
                  const Number&                             ul,
                  const Number&                             ur,
                  const Tensor<1, dim, Number>&             normal,
                  const FluxData<dim, Number>&              data,
                  NormalFlux<dim, Number>&                  flux)
   {
      switch(flux_type)
      {
         case FluxType::upwind:
            upwind_flux(ul, ur, normal, data, flux);
            break;

         default:
            AssertThrow(false, ExcMessage("Unknown numerical flux"));
      }
   }

   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE void
   boundary_flux(const Number&                            ul,
                 const Number&                            ur,
                 const Tensor<1, dim, Number>&            normal,
                 const FluxData<dim, Number>&             data,
                 NormalFlux<dim, Number>&                 flux)
   {
      upwind_flux(ul, ur, normal, data, flux);
   }

   //---------------------------------------------------------------------------
   void print_info()
   {
   }

   //---------------------------------------------------------------------------
   template <int dim>
   class Postprocessor : public DataPostprocessor<dim>
   {
   public:
      void
      evaluate_scalar_field(const DataPostprocessorInputs::Scalar<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const override
      {
         const std::vector<double> &uh = input_data.solution_values;
         for (unsigned int q = 0; q < uh.size(); ++q)
         {
            computed_quantities[q][0] = uh[q];
         }
      }
      std::vector<std::string> get_names() const override
      {
         return {"Solution"};
      }
      UpdateFlags get_needed_update_flags() const override
      {
         return update_values;
      }
   };

}
#endif
