#include <math.h>
#include "dg.h"
#include "dg1d.h"

void TimeStep(CELL * cell)
{
   UINT i;
   REAL d, u, p, c, t;

   dt = 1.0e20;

   for(i = 0; i < NC; i++) {
      d = cell[i].Un[0][0];
      u = cell[i].Un[1][0] / d;
      p = (GAMMA - 1.0) * (cell[i].Un[2][0] - 0.5 * d * u * u);
      c = sqrt(GAMMA * p / d);
      t = cell[i].h / (fabs(u) + c);
      dt = (t < dt) ? t : dt;
   }

   dt = cfl * dt;


}
