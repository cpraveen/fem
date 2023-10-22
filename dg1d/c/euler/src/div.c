#include <stdio.h>
#include <stdlib.h>
#include "dg.h"
#include "dg1d.h"

void Flux(CELL * cell)
{
   UINT i, j, k, l, cl, cr;
   REAL UL[NVAR], UR[NVAR], fl[NC + 1][NVAR], UG[NVAR], flg[NVAR], vl, vr, vx;

   for(i = 0; i < NC; i++)
      for(j = 0; j < NVAR; j++)
         for(k = 0; k < cell[i].p; k++)
            cell[i].Re[j][k] = 0.0;

   /* Loop over cell faces and find flux, periodic bc */
   for(i = 0; i <= NC; i++) {

      /* Use this for periodic bc */
      /*
      cl = (i == 0) ? (NC-1) : (i - 1);
      cr = (i == NC) ? 0 : i;
      */

      cl = (i == 0) ? i : (i - 1);
      cr = (i == NC) ? (i - 1) : i;

      Uvect(&cell[cl], cell[cl].xr, UL);
      Uvect(&cell[cr], cell[cr].xl, UR);

      switch (FLUX) {
         case LF:
            LFFlux(UL, UR, fl[i]);
            break;
         case ECUSP:
            ECUSPFlux(UL, UR, fl[i]);
            break;
         case HLLC:
            HLLCFlux(UL, UR, fl[i]);
            break;
         case AUSMDV:
            AUSMDVFlux(UL, UR, fl[i]);
            break;
         case LFC:
            LFCFlux(UL, UR, fl[i]);
            break;
         default:
            printf("Error: Flux number %d not defined\n", FLUX);
            exit(0);
      }
      //printf("%d %f %f %f\n", i, fl[i][0], fl[i][1], fl[i][2]);
   }
   //exit(0);

   /* Add interface flux to the cells */
   for(i = 0; i < NC; i++)
      for(j = 0; j < NVAR; j++)
         for(k = 0; k < cell[i].p; k++) {
            vl = ShapeFun(cell[i].xl, &cell[i], k);
            vr = ShapeFun(cell[i].xr, &cell[i], k);
            cell[i].Re[j][k] += fl[i + 1][j] * vr - fl[i][j] * vl;
         }

   /* Flux quadrature */
   for(i = 0; i < NC; i++)
      for(j = 0; j < cell[i].ng; j++) {
         Uvect(&cell[i], cell[i].xg[j], UG);
         EulerFlux(UG, flg);
         for(k = 0; k < cell[i].p; k++) {
            vx = ShapeFunDeriv(cell[i].xg[j], &cell[i], k);
            for(l = 0; l < NVAR; l++)
               cell[i].Re[l][k] -=
                  0.5 * cell[i].h * flg[l] * vx * wg[cell[i].ng - 1][j];
         }
      }

}
