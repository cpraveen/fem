#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "dg.h"
#include "dg1d.h"

/* Shape function */
REAL ShapeFun(REAL x, CELL * cell, UINT nshape)
{
   REAL Legendre(REAL, UINT);
   REAL f;

   x = 2.0 * (x - cell->x) / cell->h;
   f = sqrt(2.0 * nshape + 1);
   return f * Legendre(x, nshape);
}

/* Derivative of Shape function */
REAL ShapeFunDeriv(REAL x, CELL * cell, UINT nshape)
{
   REAL LegendreDeriv(REAL, UINT);
   REAL f;

   x = (x - cell->x) / cell->h;
   f = 2.0 * sqrt(2.0 * nshape + 1) / cell->h;
   return f * LegendreDeriv(x, nshape);
}

/* Legendre polynomials */
REAL Legendre(REAL x, UINT n)
{

   switch (n) {
      case 0:
         return 1.0;
      case 1:
         return x;
      case 2:
         return 0.5 * (3.0 * x * x - 1.0);
      case 3:
         return 0.5 * (5.0 * x * x * x - 3.0 * x);
      case 4:
         return (35.0 * x * x * x * x - 30.0 * x * x + 3.0) / 8.0;
      default:
         printf("Legendre: Unknown Legendre function no = %d\n", n);
         exit(0);
   }
}

/* Derivative of Legendre polynomial */
REAL LegendreDeriv(REAL x, UINT n)
{

   switch (n) {
      case 0:
         return 0.0;
      case 1:
         return 1.0;
      case 2:
         return 3.0 * x;
      case 3:
         return 0.5 * (15.0 * x * x - 3.0);
      case 4:
         return 0.5 * (35.0 * x * x * x - 15.0 * x);
      default:
         printf("LegendreDeriv: Unknown Legendre function no = %d\n", n);
         exit(0);
   }
}
