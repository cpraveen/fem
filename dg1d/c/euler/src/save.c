#include "dg.h"
#include "dg1d.h"

void SaveSol(CELL * cell)
{
   UINT i, j, k;

   for(i = 0; i < NC; i++)
      for(j = 0; j < NVAR; j++)
         for(k = 0; k < cell[i].p; k++)
            cell[i].Uo[j][k] = cell[i].Un[j][k];
}
