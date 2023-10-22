#include <stdio.h>
#include <math.h>
#include "dg.h"
#include "dg1d.h"

void Result(CELL * cell)
{
   FILE *fp1, *fp2;
   UINT i, j;
   REAL dx, x, U[NVAR], d, u, p, c, m;

   fp1 = fopen("out", "w");
   fp2 = fopen("sol", "w");
   for(i = 0; i < NC; i++) {
      if(NPLT == 1)
         dx = 0.5 * cell[i].h;
      else
         dx = cell[i].h / (NPLT - 1);
      for(j = 0; j < NPLT; j++) {
         if(NPLT == 1)
            x = cell[i].xl + dx;
         else
            x = cell[i].xl + dx * j;
         Uvect(&cell[i], x, U);
         d = U[0];
         u = U[1] / d;
         p = (GAMMA - 1.0) * (U[2] - 0.5 * d * u * u);
         c = sqrt(GAMMA * p / d);
         m = u / c;
         fprintf(fp1, "%f %f %f %f\n", x, U[0], U[1], U[2]);
         fprintf(fp2, "%f %f %f %f %f\n", x, d, u, p, m);
      }
      fprintf(fp1, "\n");
      fprintf(fp2, "\n");
   }
   fclose(fp1);
   fclose(fp2);
}
