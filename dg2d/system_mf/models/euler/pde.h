
//------------------------------------------------------------------------------
// Euler equations for compressible flows
//------------------------------------------------------------------------------
#ifndef __PDE_H__
#define __PDE_H__

#include <deal.II/base/vectorization.h>
#include <deal.II/numerics/data_postprocessor.h>

using namespace dealii;

constexpr unsigned int nvar = 4;

// Numerical flux functions
enum class FluxType {rusanov, roe, steger_warming, none};

std::map<std::string, FluxType> 
FluxTypeList{{"rusanov",        FluxType::rusanov},
             {"roe",            FluxType::roe},
             {"steger_warming", FluxType::steger_warming},
             {"none",           FluxType::none}};

//------------------------------------------------------------------------------
template <int dim, typename Number>
struct FluxData
{
   Point<dim, Number> p;                                 // quadrature point
   double t;                                             // time
   Tensor<1, nvar, Number> ul;                           // left  cell average
   Tensor<1, nvar, Number> ur;                           // right cell average
};

//------------------------------------------------------------------------------
// This should be set by user in a problem.h file
//------------------------------------------------------------------------------
namespace ProblemData
{
   extern const double gamma;
}

//------------------------------------------------------------------------------
//                 |    rho    |  0
// u = conserved = | rho * vel |  1,...,dim
//                 |     E     |  dim+1
//
//                 | rho |  0
// q = primitive = | vel |  1,...,dim
//                 | pre |  dim+1
//------------------------------------------------------------------------------
namespace PDE
{
   const std::string name = "2D Euler equations";
   const double gamma = ProblemData::gamma;

   constexpr unsigned int irho = 0;
   constexpr unsigned int iE = nvar-1;
   constexpr unsigned int ipre = nvar-1;

   template <int dim, typename Number>
   using FluxMatrix = Tensor<1, nvar, Tensor<1, dim, Number>>;

   template <int dim, typename Number>
   using NormalFlux = Tensor<1, nvar, Number>;

   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE void
   con2prim(const Tensor<1, nvar, Number>& u,
            Number&                        rho,
            Tensor<1, dim, Number>&        vel,
            Number&                        pre)
   {
      rho = u[0];

      Number v2 = 0.0;
      for(unsigned int d = 0; d < dim; ++d)
      {
         vel[d] = u[d + 1] / rho;
         v2 += vel[d] * vel[d];
      }

      const Number E = u[dim + 1];
      pre = (gamma - 1.0) * (E - 0.5 * rho * v2);
   }

   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   inline void
   prim2con(const Number                        rho,
            const Tensor<1, dim, Number>&       vel,
            const Number                        pre,
            Tensor<1, nvar, Number>&            u)
   {
      u[0] = rho;
      u[dim+1] = pre/(gamma - 1.0) + 0.5 * rho * vel.norm_square();

      for(unsigned int d = 0; d < dim; ++d)
         u[d+1] = rho * vel[d];
   }

   //---------------------------------------------------------------------------
   template <int dim>
   inline void
   con2prim(const Vector<double>& u, Vector<double>& q)
   {
      // density
      q[0] = u[0];

      // velocity
      double v2 = 0.0;
      for(unsigned int d = 1; d <= dim; ++d)
      {
         q[d] = u[d] / u[0];
         v2 += q[d] * q[d];
      }
      q[dim+1] = (gamma - 1.0) * (u[dim+1] - 0.5 * u[0] * v2);
   }

   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE void
   con2prim(const Tensor<1, nvar, Number>& u,
            Tensor<1, nvar, Number>&       q)
   {
      q[0] = u[0];

      Number v2 = 0.0;
      for(unsigned int d = 1; d <= dim; ++d)
      {
         q[d] = u[d] / u[0];
         v2 += q[d] * q[d];
      }


      // pressure
      q[dim+1] = (gamma - 1.0) * (u[dim+1] - 0.5 * u[0] * v2);
   }

   //---------------------------------------------------------------------------
   // q = primitive
   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   inline void
   prim2prim(const Tensor<1, nvar, Number>& q,
             Number&                        rho,
             Tensor<1, dim, Number>&        vel,
             Number&                        pre)
   {
      rho = q[0];
      pre = q[dim+1];

      for(unsigned int d=0; d<dim; ++d)
         vel[d] = q[d+1];
   }

   //---------------------------------------------------------------------------
   // q = primitive
   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE void
   physical_flux(const Tensor<1, nvar, Number>& q,
                 const Tensor<1, dim, Number>&  normal,
                 NormalFlux<dim, Number>&        flux)
   {
      Number vn = 0.0, v2 = 0.0;
      for(unsigned int d = 0; d < dim; ++d)
      {
         vn += q[d+1] * normal[d];
         v2 += q[d+1] * q[d+1];
      }

      flux[0] = q[0] * vn;
      for(unsigned int d = 0; d < dim; ++d)
         flux[d+1] = q[dim+1] * normal[d] + q[0] * q[d+1] * vn;

      const Number E = q[dim+1] / (gamma - 1.0) + 0.5 * q[0] * v2;
      flux[dim + 1] = (E + q[dim+1]) * vn;
   }

   //---------------------------------------------------------------------------
   // q = primitive
   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE Number
   max_speed(const Tensor<1, nvar, Number>& q,
             const Tensor<1, dim, Number>&  normal)
   {
      Number vn = 0.0;
      for(unsigned int d = 0; d < dim; ++d)
         vn += q[d + 1] * normal[d];

      if constexpr (std::is_same_v<Number, double>)
      {
         if(q[0] <= 0.0 || q[dim+1] <= 0.0)
         {
            std::cout << "Non-physical trace: rho, pre = " << q[0] << " " 
                      << q[dim+1] << std::endl;
         }
      }

      return std::abs(vn) + std::sqrt(gamma * q[dim + 1] / q[0]);
   }

   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE void
   rusanov_flux(const Tensor<1, nvar, Number>& ul,
                const Tensor<1, nvar, Number>& ur,
                const Tensor<1, dim, Number>&  normal,
                const FluxData<dim, Number>& data,
                NormalFlux<dim, Number>&       flux)
   {
      Tensor<1, nvar, Number> ql, qr;
      con2prim<dim>(ul, ql);
      con2prim<dim>(ur, qr);

      NormalFlux<dim, Number> fl, fr;
      physical_flux(ql, normal, fl);
      physical_flux(qr, normal, fr);

      // Wave speed from cell averages
      Tensor<1, nvar, Number> ql_avg, qr_avg;
      con2prim<dim>(data.ul, ql_avg);
      con2prim<dim>(data.ur, qr_avg);
      const Number al = max_speed(ql_avg, normal);
      const Number ar = max_speed(qr_avg, normal);
      const Number lam = std::max(al, ar);

      for(unsigned int i = 0; i < nvar; ++i)
         flux[i] = 0.5 * (fl[i] + fr[i] - lam * (ur[i] - ul[i]));
   }

   //---------------------------------------------------------------------------
   // Roe flux, taken from dflo code
   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   void
   roe_flux(const Tensor<1, nvar, Number>& ul,
            const Tensor<1, nvar, Number>& ur,
            const Tensor<1, dim, Number>&  normal,
            NormalFlux<dim, Number>&        flux)
   {
      Number rho_l, rho_r, p_l, p_r;
      Tensor<1, dim, Number> v_l, v_r;
      con2prim<dim>(ul, rho_l, v_l, p_l);
      con2prim<dim>(ur, rho_r, v_r, p_r);

      const Number rho_l_sqrt = std::sqrt(rho_l);
      const Number rho_r_sqrt = std::sqrt(rho_r);
      const Number fact_l = rho_l_sqrt / (rho_l_sqrt + rho_r_sqrt);
      const Number fact_r = 1.0 - fact_l;

      const Tensor<1, dim, Number> velocity = v_l * fact_l + v_r * fact_r;
      const Number v2 = velocity * velocity;
      const Number vel_normal = velocity * normal;

      const Number v_l_normal = v_l * normal;
      const Number v_r_normal = v_r * normal;

      const Tensor<1, dim, Number> dv = v_r - v_l;
      const Number v_dot_dv = velocity * dv;

      const Number v2_l = v_l * v_l;
      const Number v2_r = v_r * v_r;

      const Number h_l = gamma * p_l / (rho_l * (gamma-1.0)) + 0.5 * v2_l;
      const Number h_r = gamma * p_r / (rho_r * (gamma-1.0)) + 0.5 * v2_r;

      const Number density = rho_l_sqrt * rho_r_sqrt;
      const Number h = h_l * fact_l + h_r * fact_r;
      const Number c = std::sqrt((gamma-1.0) * (h - 0.5*v2));
      const Number drho = rho_r - rho_l;
      const Number dp = p_r - p_l;
      const Number dvn = v_r_normal - v_l_normal;

      const Number a1 = (dp - density * c * dvn) / (2.0*c*c);
      const Number a2 = drho - dp / (c*c);
      const Number a3 = (dp + density * c * dvn) / (2.0*c*c);

      Number l1 = std::abs(vel_normal - c);
      Number l2 = std::abs(vel_normal);
      Number l3 = std::abs(vel_normal + c);

      // entropy fix
      const Number delta = 0.1 * c;
      if constexpr (std::is_same_v<Number, double>)
      {
         if(l1 < delta) l1 = 0.5 * (l1*l1/delta + delta);
         if(l3 < delta) l3 = 0.5 * (l3*l3/delta + delta);
      }
      else
      {
         const Number l1_smooth = 0.5 * (l1 * l1 / delta + delta);
         l1 = compare_and_apply_mask<SIMDComparison::less_than>(l1, delta, l1_smooth, l1);
         const Number l3_smooth = 0.5 * (l3 * l3 / delta + delta);
         l3 = compare_and_apply_mask<SIMDComparison::less_than>(l3, delta, l3_smooth, l3);
      }

      NormalFlux<dim, Number> Dflux;
      Dflux[irho] = l1 * a1 + l2 * a2 + l3 * a3;
      Dflux[iE] =  l1 * a1 * (h - c * vel_normal)
                 + l2 * a2 * 0.5 * v2
                 + l2 * density * (v_dot_dv - vel_normal * dvn)
                 + l3 * a3 * (h + c * vel_normal);
      flux[irho] = 0.5 * (rho_l * v_l_normal + rho_r * v_r_normal
                          - Dflux[irho]);
      flux[iE] = 0.5 * (rho_l * h_l * v_l_normal + rho_r * h_r * v_r_normal
                        - Dflux[iE]);
      const Number p_avg = 0.5 * (p_l + p_r);
      for(unsigned int d=0; d<dim; ++d)
      {
         Dflux[d+1] = (velocity[d] - normal[d] * c) * l1 * a1
                  + velocity[d] * l2 * a2
                  + (dv[d] - normal[d] * dvn) * l2 * density
                  + (velocity[d] + normal[d] * c) * l3 * a3;
         flux[d+1] = normal[d] * p_avg
                   + 0.5 * (ul[d+1] * v_l_normal + ur[d+1] * v_r_normal)
                   - 0.5 * Dflux[d+1];
      }
   }

   //---------------------------------------------------------------------------
   // steger-warming flux
   // TODO: Add reference
   // See 
   //   Toro, Section 8.4.2
   //   Steger & Warming, JCP, 1981, Eq. (B9)
   //---------------------------------------------------------------------------
   template <int dim,typename Number>
   void
   steger_warming_flux(const Tensor<1, nvar, Number>& ul,
                       const Tensor<1, nvar, Number>& ur,
                       const Tensor<1, dim, Number>&  normal,
                       NormalFlux<dim, Number>&        flux)
   {
      Number rho_l, rho_r, pre_l, pre_r;
      Tensor<1, dim, Number> vel_l, vel_r;
      con2prim<dim>(ul, rho_l, vel_l, pre_l);
      con2prim<dim>(ur, rho_r, vel_r, pre_r);

      const Number c_l = std::sqrt(gamma * pre_l / rho_l);
      const Number c_r = std::sqrt(gamma * pre_r / rho_r);
      const Number vn_l = vel_l * normal;
      const Number vn_r = vel_r * normal;

      const Number zero = Number(0.0);
      // positive flux
      const Number l1p = std::max(vn_l,       zero);
      const Number l2p = std::max(vn_l + c_l, zero);
      const Number l3p = std::max(vn_l - c_l, zero);
      const Number ap  = 2.0 * (gamma - 1.0) * l1p + l2p + l3p;
      const Number fp  = 0.5 * rho_l / gamma;

      NormalFlux<dim, Number> pflux;
      pflux[0] = ap;
      for(unsigned int d=0; d<dim; ++d)
         pflux[d+1] = ap * vel_l[d] + c_l * (l2p - l3p) * normal[d];
      pflux[dim+1] = 0.5 * ap * vel_l.norm_square() +
                     c_l * vn_l * (l2p - l3p) +
                     c_l * c_l * (l2p + l3p) / (gamma - 1.0);

      // negative flux
      const Number l1m = std::min(vn_r,       zero);
      const Number l2m = std::min(vn_r + c_r, zero);
      const Number l3m = std::min(vn_r - c_r, zero);
      const Number am  = 2.0 * (gamma - 1.0) * l1m + l2m + l3m;
      const Number fm  = 0.5 * rho_r / gamma;

      NormalFlux<dim, Number> mflux;
      mflux[0] = am;
      for(unsigned int d=0; d<dim; ++d)
         mflux[d+1] = am * vel_r[d] + c_r * (l2m - l3m) * normal[d];
      mflux[dim+1] = 0.5 * am * vel_r.norm_square() +
                     c_r * vn_r * (l2m - l3m) +
                     c_r * c_r * (l2m + l3m) / (gamma - 1.0);

      // Total flux
      for(unsigned int i=0; i<nvar; ++i)
         flux[i] = fp * pflux[i] + fm * mflux[i];
   }

   //---------------------------------------------------------------------------
   // Following functions are directly called from DG solver
   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   void
   max_speed(const Tensor<1, nvar, Number>& u,
             const Point<dim, Number>&       /*p*/,
             Tensor<1, dim, Number>&         speed)
   {
      Number rho, pre;
      Tensor<1, dim, Number> vel;
      con2prim<dim>(u, rho, vel, pre);

      if constexpr (std::is_same_v<Number, double>)
      {
         if(rho <= 0.0 || pre <= 0.0)
         {
            std::cout << "Non-physical avg: rho, pre = " << rho << " " 
                      << pre << std::endl;
         }
      }

      const Number c = std::sqrt(gamma * pre / rho);

      for(unsigned int d = 0; d < dim; ++d)
         speed[d] = std::abs(vel[d]) + c;
   }

   //---------------------------------------------------------------------------
   // Flux of the PDE model: f(u,x)
   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE void
   physical_flux(const Tensor<1, nvar, Number>& u,
                 const FluxData<dim, Number>&           /*data*/,
                 FluxMatrix<dim, Number>&               flux)
   {
      Number rho, pre;
      Tensor<1, dim, Number> vel;
      con2prim<dim>(u, rho, vel, pre);

      const Number E = u[dim + 1];

      for(unsigned int d = 0; d < dim; ++d)
      {
         flux[0][d] = u[d+1];

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
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE void
   numerical_flux(const FluxType                        flux_type,
                  const Tensor<1, nvar, Number>&        ul,
                  const Tensor<1, nvar, Number>&        ur,
                  const Tensor<1, dim, Number>&         normal,
                  const FluxData<dim, Number>& data,
                  NormalFlux<dim, Number>&               flux)
   {
      switch(flux_type)
      {
         case FluxType::rusanov:
            rusanov_flux(ul, ur, normal, data, flux);
            break;

         case FluxType::roe:
            roe_flux(ul, ur, normal, flux);
            break;

         case FluxType::steger_warming:
            steger_warming_flux(ul, ur, normal, flux);
            break;

         default:
            AssertThrow(false, ExcMessage("Unknown numerical flux"));
      }
   }

   //---------------------------------------------------------------------------
   template <int dim, typename Number>
   void
   boundary_flux(const Tensor<1, nvar, Number>&        ul,
                 const Tensor<1, nvar, Number>&        ur,
                 const Tensor<1, dim, Number>&         normal,
                 const FluxData<dim, Number>&                  /*data*/,
                 NormalFlux<dim, Number>&               flux)
   {
      steger_warming_flux(ul, ur, normal, flux);
   }

   //---------------------------------------------------------------------------
   void print_info()
   {
      std::cout << "Ratio of specific heats, gamma = " << gamma << std::endl;
   }

   //---------------------------------------------------------------------------
   template <int dim>
   class Postprocessor : public DataPostprocessor<dim>
   {
   public:
      void
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const override
      {
         const std::vector<Vector<double>> &uh = input_data.solution_values;
         Assert(uh[0].size() == nvar && computed_quantities[0].size() == nvar,
                ExcInternalError());
         for (unsigned int q = 0; q < uh.size(); ++q)
         {
            con2prim<dim>(uh[q], computed_quantities[q]);
         }
      }
      std::vector<std::string> get_names() const override
      {
         if (dim == 2)
            return {"Density", "XVelocity", "YVelocity", "Pressure"};
         else
            return {"Density", "XVelocity", "YVelocity", "ZVelocity", "Pressure"};
      }
      UpdateFlags get_needed_update_flags() const override
      {
         return update_values;
      }
   };

}
#endif

