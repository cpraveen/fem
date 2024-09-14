//------------------------------------------------------------------------------
// Euler equations for compressible flows
//------------------------------------------------------------------------------
#ifndef __PDE_H__
#define __PDE_H__

using namespace dealii;

const unsigned int nvar = 4;

// Numerical flux functions
enum class FluxType {rusanov, none};

std::map<std::string, FluxType> FluxTypeList{{"rusanov", FluxType::rusanov},
                                             {"none",    FluxType::none}};

//------------------------------------------------------------------------------
// This should be set by user in a problem.h file
//------------------------------------------------------------------------------
namespace ProblemData
{
   extern const double gamma;
}

//------------------------------------------------------------------------------
//                 | rho       |  0
// u = conserved = | rho * vel |  1,...,dim
//                 |    E      |  dim+1
//
//                 | rho |  0
// q = primitive = | vel |  1,...,dim
//                 | pre |  dim+1
//------------------------------------------------------------------------------
namespace PDE
{

   const std::string name = "2D Euler equations";
   const double gamma = ProblemData::gamma;

   //---------------------------------------------------------------------------
   template <int dim>
   inline void
   con2prim(const Vector<double>& u,
            double&               rho,
            Tensor<1,dim>&        vel,
            double&               pre)
   {
      rho = u[0];

      double v2 = 0.0;
      for(unsigned int d = 0; d < dim; ++d)
      {
         vel[d] = u[d + 1] / rho;
         v2 += pow(vel[d], 2);
      }

      const double E = u[dim + 1];
      pre = (gamma - 1.0) * (E - 0.5 * rho * v2);
   }

   //---------------------------------------------------------------------------
   template <int dim>
   inline void
   con2prim(const Vector<double>& u, Vector<double>& q)
   {
      q[0] = u[0];

      double v2 = 0.0;
      for(unsigned int d = 1; d <= dim; ++d)
      {
         q[d] = u[d] / u[0];
         v2 += pow(q[d], 2);
      }

      q[dim+1] = (gamma - 1.0) * (u[dim+1] - 0.5 * u[0] * v2);
   }

   //---------------------------------------------------------------------------
   template <int dim>
   inline void
   prim2prim(const Vector<double>& q,
             double&               rho,
             Tensor<1,dim>&        vel,
             double&               pre)
   {
      rho = q[0];
      pre = q[dim+1];

      for(unsigned int d=0; d<dim; ++d)
         vel[d] = q[d+1];
   }

   //---------------------------------------------------------------------------
   template <int dim>
   void
   physical_flux(const Vector<double>& q,
                 const Tensor<1, dim>& normal,
                 Vector<double>&       flux)
   {
      double vn = 0.0, v2 = 0.0;
      for(unsigned int d = 0; d < dim; ++d)
      {
         vn += q[d+1] * normal[d];
         v2 += pow(q[d+1], 2);
      }

      flux[0] = q[0] * vn;
      for(unsigned int d = 0; d < dim; ++d)
         flux[d+1] = q[dim+1] * normal[d] + q[0] * q[d+1] * vn;

      const double E = q[dim+1] / (gamma - 1.0) + 0.5 * q[0] * v2;
      flux[dim + 1] = (E + q[dim+1]) * vn;
   }

   //---------------------------------------------------------------------------
   template <int dim>
   inline double
   max_speed(const Vector<double>&  q,
             const Tensor<1, dim>&  normal)
   {
      double vn = 0.0;
      for(unsigned int d = 0; d < dim; ++d)
         vn += q[d + 1] * normal[d];

      return abs(vn) + sqrt(gamma * q[dim + 1] / q[0]);
   }

   //---------------------------------------------------------------------------
   template <int dim>
   void
   rusanov_flux(const Vector<double>&  ul,
                const Vector<double>&  ur,
                const Tensor<1, dim>&  normal,
                Vector<double>&        flux)
   {
      Vector<double> ql(nvar), qr(nvar);
      con2prim<dim>(ul, ql);
      con2prim<dim>(ur, qr);

      Vector<double> fl(nvar), fr(nvar);
      physical_flux(ql, normal, fl);
      physical_flux(qr, normal, fr);

      const double al = max_speed(ql, normal);
      const double ar = max_speed(qr, normal);
      const double lam = std::max(al, ar);

      for(unsigned int i = 0; i < nvar; ++i)
         flux[i] = 0.5 * (fl[i] + fr[i] - lam * (ur[i] - ul[i]));
   }

   //---------------------------------------------------------------------------
   // Following functions are directly called from DG solver
   //---------------------------------------------------------------------------
   template <int dim>
   void
   max_speed(const Vector<double>& u,
             const Point<dim>&     /*p*/,
             Tensor<1, dim>&       speed)
   {
      double rho, pre;
      Tensor<1,dim> vel;
      con2prim<dim>(u, rho, vel, pre);
      const double c = sqrt(gamma * pre / rho);

      for(unsigned int d = 0; d < dim; ++d)
         speed[d] = abs(vel[d]) + c;
   }

   //---------------------------------------------------------------------------
   // Flux of the PDE model: f(u,x)
   //---------------------------------------------------------------------------
   template <int dim>
   void
   physical_flux(const Vector<double>&       u,
                 const Point<dim>&           /*p*/,
                 ndarray<double, nvar, dim>& flux)
   {
      double rho, pre;
      Tensor<1,dim> vel;
      con2prim<dim>(u, rho, vel, pre);

      const double E = u[dim + 1];

      for(unsigned int d = 0; d < dim; ++d)
      {
         flux[0][d] = u[d+1]; // mass flux

         // momentum flux
         for(unsigned int e=1; e<=dim; ++e)
            flux[e][d] = u[e] * vel[d];

         flux[d+1][d] += pre;

         // energy flux
         flux[dim + 1][d] = (E + pre) * vel[d];
      }
   }

   //---------------------------------------------------------------------------
   // Compute flux across cell faces
   //---------------------------------------------------------------------------
   template <int dim>
   void
   numerical_flux(const FluxType        flux_type,
                  const Vector<double>& ul,
                  const Vector<double>& ur,
                  const Point<dim>&     /*p*/,
                  const Tensor<1, dim>& normal,
                  Vector<double>&       flux)
   {
      switch(flux_type)
      {
         case FluxType::rusanov:
            rusanov_flux(ul, ur, normal, flux);
            break;

         default:
            AssertThrow(false, ExcMessage("Unknown numerical flux"));
      }
   }

   //---------------------------------------------------------------------------
   template <int dim>
   void
   boundary_flux(const Vector<double>& /*ul*/,
                 const Vector<double>& /*ur*/,
                 const Point<dim>&     /*p*/,
                 const Tensor<1, dim>& /*normal*/,
                 Vector<double>&       /*flux*/)
   {
      AssertThrow(false, ExcNotImplemented());
   }

   //---------------------------------------------------------------------------
   void print_info()
   {
      std::cout << "Ratio of specific heats, gamma = " << gamma << std::endl;
   }
}

#endif
