#define sign(a)   (((a) > 0.0) ? 1 : -1)

using namespace dealii;

// Advection speed
const double speed = 1.0;

// Numerical flux functions
enum class FluxType {central, upwind};

//------------------------------------------------------------------------------
// Flux of the PDE model: f(u)
//------------------------------------------------------------------------------
double
physical_flux(const double u)
{
   return speed * u;
}

//------------------------------------------------------------------------------
// Maximum wave speed: |df/du(u)|
//------------------------------------------------------------------------------
double
max_speed(const double /*u*/)
{
   return fabs(speed);
}

//------------------------------------------------------------------------------
// Central flux
//------------------------------------------------------------------------------
void
CentralFlux(const double left_state,
            const double right_state,
            double& flux)
{
   flux = 0.5 * speed * (left_state + right_state);
}

//------------------------------------------------------------------------------
// Upwind flux
//------------------------------------------------------------------------------
void
UpwindFlux(const double left_state,
           const double right_state,
           double& flux)
{
   flux  = speed * ((speed > 0) ? left_state : right_state);
}

//------------------------------------------------------------------------------
// Compute flux across cell faces
//------------------------------------------------------------------------------
void
numerical_flux(const FluxType flux_type,
               const double left_state,
               const double right_state,
               double& flux)
{
   switch(flux_type)
   {
      case FluxType::central:
         CentralFlux(left_state, right_state, flux);
         break;

      case FluxType::upwind:
         UpwindFlux(left_state, right_state, flux);
         break;

      default:
         std::cout << "Unknown flux_type !!!\n";
         abort();
   }
}
