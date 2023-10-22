#include "dg.h"
#include "dg1d.h"

void UatGauss(CELL * cell, REAL ** U)
{
   UINT iv, ig, ip;

   for(ig = 0; ig < cell->ng; ig++)
      for(iv = 0; iv < NVAR; iv++) {
         U[ig][iv] = 0.0;
         for(ip = 0; ip < cell->p; ip++)
            U[ig][iv] += cell->Uo[iv][ip] * ShapeFun(cell->xg[ig], cell, ip);
      }
}

void Uvect(CELL * cell, REAL x, REAL * U)
{
   UINT iv, ip;

   for(iv = 0; iv < NVAR; iv++) {
      U[iv] = 0.0;
      for(ip = 0; ip < cell->p; ip++)
         U[iv] += cell->Un[iv][ip] * ShapeFun(x, cell, ip);
   }
}
