#include <math.h>
#include "dg.h"
#include "dg1d.h"

void Project(CELL * cell)
{
   REAL minmod(REAL, REAL, REAL);
   REAL minmod2(REAL, REAL, REAL, REAL, REAL);
   UINT i, j, k;
   REAL u[NVAR], ux[NVAR], uxb[NVAR], dul[NVAR], dur[NVAR], R[NVAR][NVAR],
      Ri[NVAR][NVAR], fact;

   fact = sqrt(3.0);

   for(i = 1; i < NC - 1; i++) {
      for(j = 0; j < NVAR; j++) {
         dul[j] = cell[i].Un[j][0] - cell[i - 1].Un[j][0];
         dur[j] = cell[i + 1].Un[j][0] - cell[i].Un[j][0];
         u[j] = cell[i].Un[j][0];
         ux[j] = fact * cell[i].Un[j][1];
      }

      EigMat(u, R, Ri);
      Multi(Ri, ux);
      Multi(Ri, dul);
      Multi(Ri, dur);
      for(j = 0; j < NVAR; j++) {
         if(fabs(ux[j]) <= Mfact * cell[i].h * cell[i].h)
            uxb[j] = ux[j];
         else
            uxb[j] = minmod(ux[j], dul[j], dur[j]);
      }
      Multi(R, uxb);

      for(j = 0; j < NVAR; j++) {
         uxb[j] = uxb[j] / fact;
         if(cell[i].Un[j][1] != uxb[j]) {
            cell[i].Un[j][1] = uxb[j];
            for(k = 2; k < cell[i].p; k++)
               cell[i].Un[j][k] = 0.0;
         }

      }
   }
}

/* minmod limiter function */
REAL minmod(REAL a, REAL b, REAL c)
{
   REAL sgn, m;

   if(a * b <= 0.0 || b * c <= 0.0)
      return 0.0;

   sgn = (a > 0.0) ? 1.0 : -1.0;
   a = fabs(a);
   b = fabs(b);
   c = fabs(c);
   m = (a < b) ? a : b;
   m = (c < m) ? c : m;
   return sgn * m;

}

/* Eigenvector matrix */
void EigMat(REAL * U, REAL R[][3], REAL Ri[][3])
{
   REAL d, v, p, c, h, f, g1, g2;

   g1 = GAMMA - 1.0;
   g2 = g1 / 2.0;

   d = U[0];
   v = U[1] / d;
   p = (GAMMA - 1.0) * (U[2] - 0.5 * d * v * v);
   c = sqrt(GAMMA * p / d);
   h = c * c / g1 + 0.5 * v * v;
   f = d / c / 2.0;

   /* Inverse eigenvector-matrix */
   Ri[0][0] = 1.0 - g2 * v * v / c / c;
   Ri[1][0] = (g2 * v * v - v * c) / d / c;
   Ri[2][0] = -(g2 * v * v + v * c) / d / c;

   Ri[0][1] = g1 * v / c / c;
   Ri[1][1] = (c - g1 * v) / d / c;
   Ri[2][1] = (c + g1 * v) / d / c;

   Ri[0][2] = -g1 / c / c;
   Ri[1][2] = g1 / d / c;
   Ri[2][2] = -g1 / d / c;

   /* Eigenvector matrix */
   R[0][0] = 1.0;
   R[1][0] = v;
   R[2][0] = v * v / 2.0;

   R[0][1] = f;
   R[1][1] = (v + c) * f;
   R[2][1] = (h + v * c) * f;

   R[0][2] = -f;
   R[1][2] = -(v - c) * f;
   R[2][2] = -(h - v * c) * f;

}

/* Multiply matrix R and vector U */
void Multi(REAL R[][3], REAL * U)
{
   UINT i, j;
   REAL Ut[NVAR];

   for(i = 0; i < NVAR; i++)
      Ut[i] = U[i];

   for(i = 0; i < NVAR; i++) {
      U[i] = 0.0;
      for(j = 0; j < NVAR; j++)
         U[i] += R[i][j] * Ut[j];
   }
}
