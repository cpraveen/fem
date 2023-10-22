#include "dg.h"
#include "dg1d.h"

void Update(UINT rk, CELL * cell)
{
   UINT i, j, k;
   REAL f;

   for(i = 0; i < NC; i++) {
      f = dt / cell[i].h;
      for(j = 0; j < NVAR; j++)
         for(k = 0; k < cell[i].p; k++)
            cell[i].Un[j][k] =
               ark[rk] * cell[i].Uo[j][k] + brk[rk] * (cell[i].Un[j][k] -
                                                       f * cell[i].Re[j][k]);
   }
}
