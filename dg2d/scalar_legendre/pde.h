#ifndef __PDE_H__
#define __PDE_H__

using namespace dealii;

// Numerical flux functions
enum class FluxType {central, upwind};

std::map<std::string, FluxType> FluxTypeList{{"central", FluxType::central}, 
                                             {"upwind",  FluxType::upwind}};

//------------------------------------------------------------------------------
// 2d velocity field
//------------------------------------------------------------------------------
void velocity(const Point<2> &p, Tensor<1, 2> &v)
{
   v[0] = 1.0;
   v[1] = 1.0;
}

//------------------------------------------------------------------------------
// Flux of the PDE model: f(u)
//------------------------------------------------------------------------------
template <int dim>
void
physical_flux(const double u, const Point<dim>& p, Tensor<1,dim>& flux)
{
   Tensor<1,dim> v;
   velocity(p, v);
   for (int d=0; d<dim; ++d)
      flux[d] = v[d] * u;
}

//------------------------------------------------------------------------------
template <int dim>
void
flux_jacobian(const double /*u*/, const Point<dim>& p, Tensor<1,dim>& jac)
{
   velocity(p, jac);
}

//------------------------------------------------------------------------------
// Central flux
//------------------------------------------------------------------------------
template <int dim>
void
CentralFlux(const double         ul,
            const double         ur,
            const Point<dim>&    p,
            const Tensor<1,dim>& normal,
            double&              flux)
{
   Tensor<1,dim> fl, fr;
   physical_flux(ul, p, fl);
   physical_flux(ur, p, fr);
   flux = 0.5 * (fl + fr) * normal;
}

//------------------------------------------------------------------------------
// Upwind flux
//------------------------------------------------------------------------------
template <int dim>
void
UpwindFlux(const double         ul,
           const double         ur,
           const Point<dim>&    p,
           const Tensor<1,dim>& normal,
           double&              flux)
{
   Tensor<1,dim> v;
   velocity(p, v);
   double vn = v * normal;
   flux = vn * ((vn > 0.0) ? ul : ur);
}

//------------------------------------------------------------------------------
// Compute flux across cell faces
//------------------------------------------------------------------------------
template <int dim>
void
numerical_flux(const FluxType       flux_type,
               const double         ul,
               const double         ur,
               const Point<dim>&    p,
               const Tensor<1,dim>& normal,
               double& flux)
{
   switch(flux_type)
   {
      case FluxType::central:
         CentralFlux(ul, ur, p, normal, flux);
         break;

      case FluxType::upwind:
         UpwindFlux(ul, ur, p, normal, flux);
         break;

      default:
         AssertThrow(false, ExcMessage("Unknown numerical flux"));
   }
}

#endif
