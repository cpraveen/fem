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
   extern double rho;
   extern double bulk;
}

//------------------------------------------------------------------------------
// Linear acoustics
//------------------------------------------------------------------------------
namespace PDE
{

const double rho = Problem::rho;
const double bulk = Problem::bulk;
const double zz = sqrt(rho*bulk);
const double cc = sqrt(bulk/rho);

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
   const double p = u[0];
   const double v = u[1];
   flux[0] = bulk * v;
   flux[1] = p / rho;
}

//------------------------------------------------------------------------------
// Maximum wave speed: |df/du(u)|
//------------------------------------------------------------------------------
double
max_speed(const Vector<double>& /*u*/, 
          const Point<1>&       /*p*/)
{
   return cc;
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
   R(0,0) = -zz; R(0,1) = zz;
   R(1,0) = 1.0; R(1,1) = 1.0;

   L(0,0) = -0.5/zz; L(0,1) = 0.5;
   L(1,0) =  0.5/zz; L(1,1) = 0.5;
}

//------------------------------------------------------------------------------
// Compute flux across cell faces
//------------------------------------------------------------------------------
void roe_flux(const Vector<double> &ul,
              const Vector<double> &ur,
              const Point<1>&       p,
              Vector<double> &flux)
{
   // Intermediate state
   double pm = 0.5 * (ul[0] + ur[0]) - 0.5 * zz * (ur[1] - ul[1]);
   double vm = 0.5 * (ul[1] + ur[1]) - 0.5 / zz * (ur[0] - ul[0]);

   Vector<double> fl(nvar);
   physical_flux(ul, p, fl);
   flux[0] = fl[0] - cc * (pm - ul[0]);
   flux[1] = fl[1] - cc * (vm - ul[1]);
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

   flux[0] = 0.5*(fl[0] + fr[0]) - 0.5 * cc * (ur[0] - ul[0]);
   flux[1] = 0.5*(fl[1] + fr[1]) - 0.5 * cc * (ur[1] - ul[1]);
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
      DEAL_II_NOT_IMPLEMENTED();
   }
}

}
