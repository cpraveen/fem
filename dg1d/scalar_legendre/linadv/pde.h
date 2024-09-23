//------------------------------------------------------------------------------
// PDE is
//    u_t + (a u)_x = 0,  a = constant
// "a" is defined as "speed" below.
//------------------------------------------------------------------------------

using namespace dealii;

// Advection speed
const double speed = 1.0;

// Numerical flux functions
enum class FluxType {central, upwind};

std::map<std::string, FluxType> FluxTypeList{{"central", FluxType::central}, 
                                             {"upwind",  FluxType::upwind}};

//------------------------------------------------------------------------------
// Flux of the PDE model: f(u)
//------------------------------------------------------------------------------
double
physical_flux(const double u, const Point<1>& /*p*/)
{
   return speed * u;
}

//------------------------------------------------------------------------------
// Maximum wave speed: |df/du(u)|
//------------------------------------------------------------------------------
double
max_speed(const double /*u*/, const Point<1>& /*p*/)
{
   return fabs(speed);
}

//------------------------------------------------------------------------------
// Central flux
//------------------------------------------------------------------------------
void
CentralFlux(const double    ul,
            const double    ur,
            const Point<1>& /*p*/,
            double&         flux)
{
   flux = 0.5 * speed * (ul + ur);
}

//------------------------------------------------------------------------------
// Upwind flux
//------------------------------------------------------------------------------
void
UpwindFlux(const double    ul,
           const double    ur,
           const Point<1>& /*p*/,
           double&         flux)
{
   flux  = speed * ((speed > 0) ? ul : ur);
}

//------------------------------------------------------------------------------
// Compute flux across cell faces
//------------------------------------------------------------------------------
void
numerical_flux(const FluxType  flux_type,
               const double    ul,
               const double    ur,
               const Point<1>& p,
               double& flux)
{
   switch(flux_type)
   {
      case FluxType::central:
         CentralFlux(ul, ur, p, flux);
         break;

      case FluxType::upwind:
         UpwindFlux(ul, ur, p, flux);
         break;

      default:
         AssertThrow(false, ExcMessage("Unknown flux type"));
   }
}
