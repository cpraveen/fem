#include<stdio.h>
#include<math.h>
#include "dg.h"
#include "dg1d.h"

/* Gauss integration points and weights */
void GaussInit()
{
   UINT n = 10, i, j;

   printf("Calculating Gauss integration points and weights ...\n");

   for(i = 0; i < n; i++)
      for(j = 0; j < n; j++) {
         xg[i][j] = 0.0;
         wg[i][j] = 0.0;
      }

   /* Gauss integration points */
   xg[0][0] = 0.0;

   xg[1][0] = -1.0 / sqrt(3.0);
   xg[1][1] = 1.0 / sqrt(3.0);

   xg[2][0] = -sqrt(15.0) / 5.0;
   xg[2][1] = 0.0;
   xg[2][2] = sqrt(15.0) / 5.0;

   xg[3][0] = -sqrt(525.0 + 70.0 * sqrt(30.0)) / 35.0;
   xg[3][1] = -sqrt(525.0 - 70.0 * sqrt(30.0)) / 35.0;
   xg[3][2] = sqrt(525.0 - 70.0 * sqrt(30.0)) / 35.0;
   xg[3][3] = sqrt(525.0 + 70.0 * sqrt(30.0)) / 35.0;

   xg[4][0] = -sqrt(245.0 + 14.0 * sqrt(70.0)) / 21.0;
   xg[4][1] = -sqrt(245.0 - 14.0 * sqrt(70.0)) / 21.0;
   xg[4][2] = 0.0;
   xg[4][3] = sqrt(245.0 - 14.0 * sqrt(70.0)) / 21.0;
   xg[4][4] = sqrt(245.0 + 14.0 * sqrt(70.0)) / 21.0;

   /* Gaussian weights */
   wg[0][0] = 2.0;

   wg[1][0] = 1.0;
   wg[1][1] = 1.0;

   wg[2][0] = 5.0 / 9.0;
   wg[2][1] = 8.0 / 9.0;
   wg[2][2] = 5.0 / 9.0;

   wg[3][0] = (18.0 - sqrt(30.0)) / 36.0;
   wg[3][1] = (18.0 + sqrt(30.0)) / 36.0;
   wg[3][2] = (18.0 + sqrt(30.0)) / 36.0;
   wg[3][3] = (18.0 - sqrt(30.0)) / 36.0;

   wg[4][0] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
   wg[4][1] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
   wg[4][2] = 128.0 / 225.0;
   wg[4][3] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
   wg[4][4] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;

}

/* Find physical coordinates of Gauss points in each cell */
void GaussPoints(CELL * cell)
{
   REAL xl, xr;
   UINT i;

   xl = cell->xl;
   xr = cell->xr;

   for(i = 0; i < cell->ng; i++)
      cell->xg[i] = 0.5 * (xl * (1.0 - xg[cell->ng - 1][i]) +
                           xr * (1.0 + xg[cell->ng - 1][i]));

}

/* Perform Gauss quadrature */
REAL GaussQuadrature(REAL * f, UINT ng)
{
   REAL integral = 0.0;
   UINT i;

   for(i = 0; i < ng; i++)
      integral += f[i] * wg[ng - 1][i];

   return integral;
}
