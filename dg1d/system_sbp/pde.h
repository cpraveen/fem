//------------------------------------------------------------------------------
// 1d compressible Euler equations
//------------------------------------------------------------------------------

using namespace dealii;

// Number of PDE in the system
const unsigned int nvar = 3;

//------------------------------------------------------------------------------
// Extend namespace Problem with data specific to this pde
// This data must be set in problem_data.h file.
//------------------------------------------------------------------------------
namespace Problem
{
   extern double gamma;
}

//------------------------------------------------------------------------------
// Linear acoustics
//------------------------------------------------------------------------------
namespace PDE
{

   const double gamma = Problem::gamma;

   // Numerical flux functions
   enum class FluxType {roe, rusanov};
   std::map<std::string, FluxType> FluxTypeList{{"roe",     FluxType::roe},
                                                {"rusanov", FluxType::rusanov}};

   // Two point volume flux
   enum class VFluxType {central, keeppe};
   std::map<std::string, VFluxType> VFluxTypeList{{"central", VFluxType::central},
                                                  {"keeppe",  VFluxType::keeppe}};
   //---------------------------------------------------------------------------
   // Flux of the PDE model: f(u)
   //---------------------------------------------------------------------------
   void
   physical_flux(const Vector<double>& u,
                 const Point<1>&       /*p*/,
                 Vector<double>&       flux)
   {
      const double rho = u[0];
      const double vel = u[1] / rho;
      const double pre = (gamma - 1.0) * (u[2] - 0.5 * rho * pow(vel, 2));
      flux[0] = rho * vel;
      flux[1] = pre + rho * pow(vel, 2);
      flux[2] = (u[2] + pre) * vel;
   }

   //---------------------------------------------------------------------------
   // Maximum wave speed: |df/du(u)|
   //---------------------------------------------------------------------------
   double
   max_speed(const Vector<double>& u,
             const Point<1>&       /*p*/)
   {
      const double rho = u[0];
      const double vel = u[1] / rho;
      const double pre = (gamma - 1.0) * (u[2] - 0.5 * rho * pow(vel, 2));
      return abs(vel) + sqrt(gamma * pre / rho);
   }

   //---------------------------------------------------------------------------
   // R = matrix of right eigenvectors, columns are right eigenvectors
   // L = matrix of left eigenvectors = R^(-1), rows are left eigenvectors
   //---------------------------------------------------------------------------
   void
   char_mat(const Vector<double>& u,
            const Point<1>&       /*p*/,
            FullMatrix<double>&   R,
            FullMatrix<double>&   L)
   {
      const double rho = u[0];
      const double vel = u[1] / rho;
      const double pre = (gamma - 1.0) * (u[2] - 0.5 * rho * pow(vel, 2));
      const double H = (u[2] + pre) / rho;
      const double a = sqrt(gamma * pre / rho);

      R(0, 0) = 1.0;     R(0, 1) = 1.0;            R(0, 2) = 1.0;
      R(1, 0) = vel-a;   R(1, 1) = vel;            R(1, 2) = vel+a;
      R(2, 0) = H-vel*a; R(2, 1) = 0.5*pow(vel,2); R(2, 2) = H+vel*a;

      const double M = vel/a;
      const double M2 = pow(M,2);
      const double g1 = gamma - 1.0;

      L(0, 0) = 0.25*g1*M2 + 0.5*M;
      L(1, 0) = 1.0 - 0.5*g1*M2;
      L(2, 0) = 0.25*g1*M2 - 0.5*M;

      L(0,1) = -0.5*g1*M/a - 0.5/a;
      L(1,1) = g1*M/a;
      L(2,1) = -0.5*g1*M/a + 0.5/a;

      L(0,2) = 0.5*g1/pow(a,2);
      L(1,2) = -g1/pow(a,2);
      L(2,2) = 0.5*g1/pow(a,2);
   }


   //---------------------------------------------------------------------------
   // Compute flux across cell faces
   //---------------------------------------------------------------------------
   void
   roe_flux(const Vector<double>& /*ul*/,
            const Vector<double>& /*ur*/,
            const Point<1>&       /*p*/,
            Vector<double>&       /*flux*/)
   {
      DEAL_II_NOT_IMPLEMENTED();
   }

   //---------------------------------------------------------------------------
   // Compute flux across cell faces
   //---------------------------------------------------------------------------
   void
   rusanov_flux(const Vector<double>& ul,
                const Vector<double>& ur,
                const Point<1>&       p,
                Vector<double>&       flux)
   {
      Vector<double> fl(nvar);
      physical_flux(ul, p, fl);

      Vector<double> fr(nvar);
      physical_flux(ur, p, fr);

      double cl = max_speed(ul, p);
      double cr = max_speed(ur, p);
      double lam = std::max(cl, cr);

      for(unsigned int i=0; i<nvar; ++i)
         flux[i] = 0.5 * (fl[i] + fr[i]) - 0.5 * lam * (ur[i] - ul[i]);
   }

   //---------------------------------------------------------------------------
   // Compute flux across cell faces
   //---------------------------------------------------------------------------
   void
   numerical_flux(const FluxType        flux_type,
                  const Vector<double>& ul,
                  const Vector<double>& ur,
                  const Point<1>&       p,
                  Vector<double>&       flux)
   {
      switch(flux_type)
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

   //---------------------------------------------------------------------------
   // Following are needed only for SBP schemes
   //---------------------------------------------------------------------------
   // Average flux
   //---------------------------------------------------------------------------
   void
   central_flux(const Vector<double> &ul,
                const Vector<double> &ur,
                Vector<double> &flux)
   {
      Point<1> p(0); // dummy

      Vector<double> fl(nvar);
      physical_flux(ul, p, fl);

      Vector<double> fr(nvar);
      physical_flux(ur, p, fr);

      for (unsigned int i = 0; i < nvar; ++i)
         flux[i] = 0.5 * (fl[i] + fr[i]);
   }

   //---------------------------------------------------------------------------
   // KEEP-PE flux of Shima, et al.
   // https://doi.org/10.1016/j.jcp.2020.110060
   //---------------------------------------------------------------------------
   void
   keeppe_flux(const Vector<double>& ul,
               const Vector<double>& ur,
               Vector<double>&       flux)
   {
      const double dl = ul[0];
      const double vl = ul[1] / dl;
      const double pl = (gamma - 1.0) * (ul[2] - 0.5 * dl * pow(vl, 2));

      const double dr = ur[0];
      const double vr = ur[1] / dr;
      const double pr = (gamma - 1.0) * (ur[2] - 0.5 * dr * pow(vr, 2));

      const double rho = 0.5 * (dl + dr);
      const double vel = 0.5 * (vl + vr);
      const double pre = 0.5 * (pl + pr);
      const double v2  = vl * vr;
      const double pv  = 0.5 * (pl * vr + pr * vl);

      flux[0] = rho * vel;
      flux[1] = pre + rho * pow(vel, 2);
      flux[2] = pre * vel / (gamma - 1.0) + 0.5 * rho * vel * v2 + pv;
   }

   //---------------------------------------------------------------------------
   // Compute flux across cell faces
   //---------------------------------------------------------------------------
   void
   volume_flux(const VFluxType       flux_type,
               const Vector<double>& ul,
               const Vector<double>& ur,
               Vector<double>&       flux)
   {
      switch (flux_type)
      {
      case VFluxType::central:
         central_flux(ul, ur, flux);
         break;

      case VFluxType::keeppe:
         keeppe_flux(ul, ur, flux);
         break;

      default:
         DEAL_II_NOT_IMPLEMENTED();
      }
   }
}
