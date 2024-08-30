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
physical_flux(const double u)
{
   return 0.5 * std::pow(u, 2);
}

//------------------------------------------------------------------------------
// Maximum wave speed: |df/du(u)|
//------------------------------------------------------------------------------
double
max_speed(const double u)
{
   return std::fabs(u);
}

//------------------------------------------------------------------------------
// Central flux
//------------------------------------------------------------------------------
void
CentralFlux(const double left_state,
            const double right_state,
            double&      flux)
{
   flux = 0.5 * (physical_flux(left_state) + physical_flux(right_state));
}

//------------------------------------------------------------------------------
// Roe flux
//------------------------------------------------------------------------------
void
RoeFlux(const double left_state,
        const double right_state,
        double&      flux)
{
   double speed = 0.5 * (left_state + right_state);
   flux = 0.5 * (physical_flux(left_state) + physical_flux(right_state))
          - 0.5 * std::fabs(speed) * (right_state - left_state);
}

//------------------------------------------------------------------------------
// Godunov flux
//------------------------------------------------------------------------------
void
GodunovFlux(const double left_state,
            const double right_state,
            double&      flux)
{
   double u1 = std::max(0.0, left_state);
   double u2 = std::min(0.0, right_state);
   flux = std::max(physical_flux(u1), physical_flux(u2));
}

//------------------------------------------------------------------------------
// Rusanov flux
//------------------------------------------------------------------------------
void
RusanovFlux(const double left_state,
            const double right_state,
            double&      flux)
{
   double lambda = std::max(std::fabs(left_state), std::fabs(right_state));
   flux = 0.5 * (physical_flux(left_state) + physical_flux(right_state))
          - 0.5 * lambda * (right_state - left_state);
}

//------------------------------------------------------------------------------
// Compute flux across cell faces
//------------------------------------------------------------------------------
void
numerical_flux(const FluxType flux_type,
               const double   left_state,
               const double   right_state,
               double&        flux)
{
   switch(flux_type)
   {
      case FluxType::central:
         CentralFlux(left_state, right_state, flux);
         break;

      case FluxType::roe:
         RoeFlux(left_state, right_state, flux);
         break;

      case FluxType::godunov:
         GodunovFlux(left_state, right_state, flux);
         break;

      case FluxType::rusanov:
         RusanovFlux(left_state, right_state, flux);
         break;

      default:
         AssertThrow(false, ExcMessage("Unknown flux type"));
   }
}
