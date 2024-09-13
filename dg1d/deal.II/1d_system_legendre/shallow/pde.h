//------------------------------------------------------------------------------
// PDE is linear acoustics
//    p_t + K * v_x = 0
//    v_t + (1/rho) * p_x = 0
//------------------------------------------------------------------------------

using namespace dealii;

// Number of PDE in the system
const unsigned int nvar = 2;

//------------------------------------------------------------------------------
// Extend namespace Problem with data specific to this pde
// This data must be set in problem_data.h file.
//------------------------------------------------------------------------------
namespace Problem 
{
   extern double g;
}

//------------------------------------------------------------------------------
// Linear acoustics
//------------------------------------------------------------------------------
namespace PDE
{

const double g = Problem::g;

// Numerical flux functions
enum class FluxType {roe, rusanov};

std::map<std::string, FluxType> FluxTypeList{{"roe",     FluxType::roe}, 
                                             {"rusanov", FluxType::rusanov}};

//------------------------------------------------------------------------------
// Flux of the PDE model: f(u)
//------------------------------------------------------------------------------
void
physical_flux(const Vector<double>& u, 
              const Point<1>&       /*p*/, 
              Vector<double>&       flux)
{
   const double h = u[0];
   const double v = u[1] / h;
   flux[0] = h * v;
   flux[1] = h * pow(v,2) + 0.5 * g * pow(h,2);
}

//------------------------------------------------------------------------------
// Maximum wave speed: |df/du(u)|
//------------------------------------------------------------------------------
double
max_speed(const Vector<double>& u, 
          const Point<1>&       /*p*/)
{
   const double h = u[0];
   const double v = u[1] / h;
   return abs(v) + sqrt(g * h);
}

//------------------------------------------------------------------------------
// R = matrix of right eigenvectors, columns are right eigenvectors
// L = matrix of left eigenvectors = R^(-1), rows are left eigenvectors
//------------------------------------------------------------------------------
void 
char_mat(const Vector<double>& /*u*/,
         const Point<1>&       /*p*/,
         FullMatrix<double>&   R,
         FullMatrix<double>&   L)
{
   R(0,0) = 1.0; R(0,1) = 0.0;
   R(1,0) = 0.0; R(1,1) = 1.0;

   L(0,0) = 1.0; L(0,1) = 0.0;
   L(1,0) = 0.0; L(1,1) = 1.0;
}

//------------------------------------------------------------------------------
// Compute flux across cell faces
//------------------------------------------------------------------------------
void roe_flux(const Vector<double> &/*ul*/,
              const Vector<double> &/*ur*/,
              const Point<1>&       /*p*/,
              Vector<double> &/*flux*/)
{
   AssertThrow(false, ExcNotImplemented());
}

//------------------------------------------------------------------------------
// Compute flux across cell faces
//------------------------------------------------------------------------------
void rusanov_flux(const Vector<double> &ul,
                  const Vector<double> &ur,
                  const Point<1>&       p,
                  Vector<double> &flux)
{
   Vector<double> fl(nvar);
   physical_flux(ul, p, fl);

   Vector<double> fr(nvar);
   physical_flux(ur, p, fr);

   double cl = max_speed(ul, p);
   double cr = max_speed(ur, p);
   double lam = std::max(cl, cr);

   flux[0] = 0.5*(fl[0] + fr[0]) - 0.5 * lam * (ur[0] - ul[0]);
   flux[1] = 0.5*(fl[1] + fr[1]) - 0.5 * lam * (ur[1] - ul[1]);
}

//------------------------------------------------------------------------------
// Compute flux across cell faces
//------------------------------------------------------------------------------
void numerical_flux(const FluxType        flux_type,
                    const Vector<double>& ul,
                    const Vector<double>& ur,
                    const Point<1>&       p,
                    Vector<double>&       flux)
{
   switch (flux_type)
   {
   case FluxType::roe:
      roe_flux(ul, ur, p, flux);
      break;

   case FluxType::rusanov:
      rusanov_flux(ul, ur, p, flux);
      break;

   default:
      AssertThrow(false, ExcMessage("Unknown flux type"));
   }
}

}
