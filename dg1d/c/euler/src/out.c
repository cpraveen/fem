#include "dg.h"
#include "dg1d.h"
void ausmflux(REAL * ql, REAL * qr, REAL * fl)
{

   UINT i;
   REAL rho_left, rho_right, u_left, u_right,
      p_left, p_right, c_left, c_right, c_mean,
      h_left, h_right, amach_l, amach_r,
      alpha_l, alpha_r, temp1, temp2, u_lp, u_rm,
      p_lp, p_rm, ak_half, p_half, s_half, flux_v1,
      flux_d1, flux_p, flux_n, p_min, a_min, gam, gam1;

   gam = 1.4;
   gam1 = gam - 1.0;

   rho_left = ql[0];
   u_left = ql[1] / ql[0];
   p_left = gam1 * (ql[2] - 0.5 * ql[1] * ql[1] / ql[0]);

   rho_right = qr[0];
   u_right = qr[1] / qr[0];
   p_right = gam1 * (qr[2] - 0.5 * qr[1] * qr[1] / qr[0]);

   c_left = sqrt(gam * p_left / rho_left);
   h_left = c_left * c_left / gam1 + 0.5 * u_left * u_left;

//......right state 
   c_right = dsqrt(gam * p_right / rho_right);
   h_right = c_right * c_right / gam1 + 0.5 * u_right * u_right;

//......mean speed of sound
   c_mean = (c_left > c_right) ? c_left : c_right;

//......left and right mach no.
   amach_l = fabs(u_left) / c_mean;
   amach_r = fabs(u_right) / c_mean;

//......alpha
   temp1 = p_left / rho_left;
   temp2 = p_right / rho_right;

   alpha_l = 2 * temp1 / (temp1 + temp2);
   alpha_r = 2 * temp2 / (temp1 + temp2);

//......u_left+
   temp1 = u_left + c_mean;
   temp2 = u_left + dabs(u_left);
   if(amach_l <= 1)
      u_lp = alpha_l * (0.25 * temp1 * temp1 / c_mean - 0.5 * temp2)
         + 0.5 * temp2;
   else
      u_lp = 0.5 * temp2;

//......u_right-
   temp1 = u_right - c_mean;
   temp2 = u_right - dabs(u_right);
   if(amach_r <= 1)
      u_rm = alpha_r * (-0.25 * temp1 * temp1 / c_mean - 0.5 * temp2)
         + 0.5 * temp2;
   else
      u_rm = 0.5 * temp2;

//......p_left+
   temp1 = u_left / c_mean + 1;
   temp2 = 2 - u_left / c_mean;
   if(amach_l <= 1)
      p_lp = 0.25 * p_left * temp1 * temp1 * temp2;
   else
      p_lp = 0.5 * p_left * (u_left + fabs(u_left)) / (u_left);

//......p_right-
   temp1 = u_right / c_mean - 1;
   temp2 = 2 + u_right / c_mean;
   if(amach_r <= 1)
      p_rm = 0.25 * p_right * temp1 * temp1 * temp2;
   else
      p_rm = 0.5 * p_right * (u_right - fabs(u_right)) / (u_right);

   ak_half = 10.0;
   p_half = p_lp + p_rm;
   p_min = (p_left < p_right) ? p_left : p_right;
   a_min = ak_half * s_half;
   a_min = (a_min < 1.0) ? a_min : 1.0;
   s_half = fabs(p_right - p_left) / p_min;
   s_half = 0.5 * a_min;

//......calculating the inviscid interface fluxes

//......mass flux
   fl[0] = u_lp * rho_left + u_rm * rho_right;

//......enthalpy flux
   fl[2] = 0.5 * (fl[0] * (h_left + h_right)
                  - fabs(fl[0]) * (h_right - h_left));

//......momentum flux

   flux_v1 = rho_left * u_left * u_lp + rho_right * u_right * u_rm;

   flux_d1 = 0.5 * (fl[0] * (u_left + u_right) -
                    fabs(fl[0]) * (u_right - u_left));

   flux_p = p_half;

   flux_n = (0.5 + s_half) * flux_v1 + (0.5 - s_half) * flux_d1;

   fl[1] = flux_n + flux_p;
}
