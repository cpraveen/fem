//------------------------------------------------------------------------------
// PDE is
//    u_t + a(x) u_x = 0
// and is written as conservation law with source term
//    u_t + (a(x) u)_x = a'(x) u
//------------------------------------------------------------------------------

using namespace dealii;

// Numerical flux functions
enum class FluxType {central, upwind};

std::map<std::string, FluxType> FluxTypeList{{"central", FluxType::central}, 
                                             {"upwind",  FluxType::upwind}};

//------------------------------------------------------------------------------
// Trefethen: Spectral Methods
// Period = 
//------------------------------------------------------------------------------
double speed(const Point<1>& p)
{
   return 0.2 + pow(sin(p[0]-1.0),2);
}

//------------------------------------------------------------------------------
// Flux of the PDE model: f(u)
//------------------------------------------------------------------------------
double
physical_flux(const double u, const Point<1>& p)
{
   return speed(p) * u;
}

//------------------------------------------------------------------------------
// Maximum wave speed: |df/du(u)|
//------------------------------------------------------------------------------
double
max_speed(const double /*u*/, const Point<1>& p)
{
   return fabs(speed(p));
}

//------------------------------------------------------------------------------
// Central flux
//------------------------------------------------------------------------------
void
CentralFlux(const double    ul,
            const double    ur,
            const Point<1>& p,
            double&         flux)
{
   flux = 0.5 * speed(p) * (ul + ur);
}

//------------------------------------------------------------------------------
// Upwind flux
//------------------------------------------------------------------------------
void
UpwindFlux(const double    ul,
           const double    ur,
           const Point<1>& p,
           double&         flux)
{
   double a = speed(p);
   flux  = a * ((a > 0) ? ul : ur);
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

//------------------------------------------------------------------------------
double source(const double u, const Point<1>& p)
{
   double da = 2.0 * sin(p[0]-1.0) * cos(p[0]-1.0);
   return da * u;
}
