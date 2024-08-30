//------------------------------------------------------------------------------
// PDE is
//    u_t + (a(x) u)_x = 0
//------------------------------------------------------------------------------

using namespace dealii;

// Numerical flux functions
enum class FluxType {central, upwind};

std::map<std::string, FluxType> FluxTypeList{{"central", FluxType::central}, 
                                             {"upwind",  FluxType::upwind}};

//------------------------------------------------------------------------------
// Hesthaven DG book: Example 5.3
//------------------------------------------------------------------------------
double speed(const Point<1>& p)
{
   double x2 = pow(p[0],2);
   return 1.0 + pow(1-x2,5);
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
