//------------------------------------------------------------------------------
// Euler equations for compressible flows
//------------------------------------------------------------------------------
#ifndef __PDE_H__
#define __PDE_H__

#include <deal.II/numerics/data_postprocessor.h>

using namespace dealii;

constexpr unsigned int nvar = 1;

// Numerical flux functions
enum class FluxType {upwind, none};

std::map<std::string, FluxType> FluxTypeList{{"upwind", FluxType::upwind},
                                             {"none",   FluxType::none}};

//------------------------------------------------------------------------------
// This should be set by user in a problem.h file
//------------------------------------------------------------------------------
namespace ProblemData
{
   template <int dim>
   extern void velocity(const Point<dim>& p, Tensor<1,dim>& v);
}

//------------------------------------------------------------------------------
namespace PDE
{

   const std::string name = "2D linear advection equation";
   using ProblemData::velocity;

   //---------------------------------------------------------------------------
   template <int dim>
   void
   upwind_flux(const Vector<double>&  ul,
               const Vector<double>&  ur,
               const Point<dim>&      p,
               const Tensor<1, dim>&  normal,
               Vector<double>&        flux)
   {
      Tensor<1,dim> vel;
      velocity(p, vel);
      const auto vn = vel * normal;
      flux[0] = vn * ((vn > 0.0) ? ul[0] : ur[0]);
   }

   //---------------------------------------------------------------------------
   // Following functions are directly called from DG solver
   //---------------------------------------------------------------------------
   template <int dim>
   void
   max_speed(const Vector<double>& /*u*/,
             const Point<dim>&     p,
             Tensor<1, dim>&       speed)
   {
      velocity(p, speed);
   }

   //---------------------------------------------------------------------------
   // Flux of the PDE model: f(u,x)
   //---------------------------------------------------------------------------
   template <int dim>
   void
   physical_flux(const Vector<double>&       u,
                 const Point<dim>&           p,
                 ndarray<double, nvar, dim>& flux)
   {
      Tensor<1,dim> vel;
      velocity(p, vel);
      flux[0][0] = vel[0] * u[0];
      flux[0][1] = vel[1] * u[0];
   }

   //---------------------------------------------------------------------------
   // Compute flux across cell faces
   //---------------------------------------------------------------------------
   template <int dim>
   void
   numerical_flux(const FluxType        flux_type,
                  const Vector<double>& ul,
                  const Vector<double>& ur,
                  const Point<dim>&     p,
                  const Tensor<1, dim>& normal,
                  Vector<double>&       flux)
   {
      switch(flux_type)
      {
         case FluxType::upwind:
            upwind_flux(ul, ur, p, normal, flux);
            break;

         default:
            AssertThrow(false, ExcMessage("Unknown numerical flux"));
      }
   }

   //---------------------------------------------------------------------------
   template <int dim>
   void
   boundary_flux(const Vector<double>& ul,
                 const Vector<double>& ur,
                 const Point<dim>&     p,
                 const Tensor<1, dim>& normal,
                 Vector<double>&       flux)
   {
      upwind_flux(ul, ur, p, normal, flux);
   }

   //---------------------------------------------------------------------------
   void
   char_mat(const Vector<double>& /*sol*/,
            const Point<2>&       /*p*/,
            const Tensor<1, 2>&   /*ex*/,
            const Tensor<1, 2>&   /*ey*/,
            FullMatrix<double>&   Rx,
            FullMatrix<double>&   Lx,
            FullMatrix<double>&   Ry,
            FullMatrix<double>&   Ly)
   {
      Rx[0][0] = 1.0;
      Ry[0][0] = 1.0;
      Lx[0][0] = 1.0;
      Ly[0][0] = 1.0;
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
