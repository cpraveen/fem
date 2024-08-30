//------------------------------------------------------------------------------
// PDE is
//    u_t + f(u)_x = 0,   f(u) = u^2/2
//------------------------------------------------------------------------------
using namespace dealii;

// Numerical flux functions
enum class FluxType {central, roe, godunov, rusanov};

std::map<std::string, FluxType> FluxTypeList{{"central", FluxType::central},
                                             {"roe",     FluxType::roe},
                                             {"godunov", FluxType::godunov},
                                             {"rusanov", FluxType::rusanov}};
//------------------------------------------------------------------------------
// Flux of the PDE model: f(u)
//------------------------------------------------------------------------------
double
physical_flux(const double u, const Point<1>& p = Point<1>(0.0))
{
   return 0.5 * std::pow(u, 2);
}

//------------------------------------------------------------------------------
// Maximum wave speed: |df/du(u)|
//------------------------------------------------------------------------------
double
max_speed(const double u, const Point<1>& /*p*/)
{
   return std::fabs(u);
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
   flux = 0.5 * (physical_flux(ul) + physical_flux(ur));
}

//------------------------------------------------------------------------------
// Roe flux
//------------------------------------------------------------------------------
void
RoeFlux(const double    ul,
        const double    ur,
        const Point<1>& /*p*/,
        double&         flux)
{
   double speed = 0.5 * (ul + ur);
   flux = 0.5 * (physical_flux(ul) + physical_flux(ur))
          - 0.5 * std::fabs(speed) * (ur - ul);
}

//------------------------------------------------------------------------------
// Godunov flux
//------------------------------------------------------------------------------
void
GodunovFlux(const double    ul,
            const double    ur,
            const Point<1>& /*p*/,
            double&         flux)
{
   double u1 = std::max(0.0, ul);
   double u2 = std::min(0.0, ur);
   flux = std::max(physical_flux(u1), physical_flux(u2));
}

//------------------------------------------------------------------------------
// Rusanov flux
//------------------------------------------------------------------------------
void
RusanovFlux(const double    ul,
            const double    ur,
            const Point<1>& /*p*/,
            double&         flux)
{
   double lambda = std::max(std::fabs(ul), std::fabs(ur));
   flux = 0.5 * (physical_flux(ul) + physical_flux(ur))
          - 0.5 * lambda * (ur - ul);
}

//------------------------------------------------------------------------------
// Compute flux across cell faces
//------------------------------------------------------------------------------
void
numerical_flux(const FluxType  flux_type,
               const double    ul,
               const double    ur,
               const Point<1>& p,
               double&         flux)
{
   switch(flux_type)
   {
      case FluxType::central:
         CentralFlux(ul, ur, p, flux);
         break;

      case FluxType::roe:
         RoeFlux(ul, ur, p, flux);
         break;

      case FluxType::godunov:
         GodunovFlux(ul, ur, p, flux);
         break;

      case FluxType::rusanov:
         RusanovFlux(ul, ur, p, flux);
         break;

      default:
         AssertThrow(false, ExcMessage("Unknown flux type"));
   }
}
