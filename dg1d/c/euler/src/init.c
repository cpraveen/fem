#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "dg.h"
#include "dg1d.h"

/* Allocate memory for main structure and set initial conditions */
CELL *Init()
{
   void ReadInput();
   void InitCondEuler(REAL, REAL *);
   UINT i, j, k, l;
   REAL dx, U[NVAR], v;
   CELL *cell;

   ReadInput();

   /* Coefficients for RK3 */
   ark[0] = 0.0;
   ark[1] = 3.0 / 4.0;
   ark[2] = 1.0 / 3.0;

   brk[0] = 1.0;
   brk[1] = 1.0 / 4.0;
   brk[2] = 2.0 / 3.0;

   NG = PORD + 2;

   printf("Allocating memory and setting initial condition ...\n");

   cell = (CELL *) calloc(NC, sizeof(CELL));
   if(cell == NULL) {
      printf("Init: Could not allocate cell\n");
      exit(0);
   }

   dx = (xmax - xmin) / NC;
   printf("No of cells = %d\n", NC);
   printf("         dx = %f\n", dx);

   for(i = 0; i < NC; i++) {
      cell[i].xl = xmin + i * dx;
      cell[i].xr = cell[i].xl + dx;
      cell[i].x = 0.5 * (cell[i].xl + cell[i].xr);
      cell[i].h = cell[i].xr - cell[i].xl;


      cell[i].p = PORD;

      cell[i].ng = NG;
      cell[i].xg = (REAL *) calloc(cell[i].ng, sizeof(REAL));
      GaussPoints(&cell[i]);

      cell[i].Un = (REAL **) calloc(NVAR, sizeof(REAL *));
      cell[i].Uo = (REAL **) calloc(NVAR, sizeof(REAL *));
      cell[i].Re = (REAL **) calloc(NVAR, sizeof(REAL *));
      for(j = 0; j < NVAR; j++) {
         cell[i].Un[j] = (REAL *) calloc(cell[i].p, sizeof(REAL));
         cell[i].Uo[j] = (REAL *) calloc(cell[i].p, sizeof(REAL));
         cell[i].Re[j] = (REAL *) calloc(cell[i].p, sizeof(REAL));
      }
   }

   /* Set initial condition by L2 projection */
   for(i = 0; i < NC; i++) {

      for(j = 0; j < NVAR; j++)
         for(k = 0; k < cell[i].p; k++)
            cell[i].Un[j][k] = 0.0;

      for(j = 0; j < cell[i].p; j++)
         for(k = 0; k < cell[i].ng; k++) {
            InitCondEuler(cell[i].xg[k], U);
            v = ShapeFun(cell[i].xg[k], &cell[i], j);
            for(l = 0; l < NVAR; l++)
               cell[i].Un[l][j] += 0.5 * U[l] * v * wg[cell[i].ng - 1][k];
         }
   }

   return cell;
}

void ReadInput()
{
   FILE *fp;
   char dummy[100];
   fp = fopen("inp.dat", "r");
   if(fp == NULL) {
      printf("Error: Could not open inp.dat\n");
      exit(0);
   }

   fscanf(fp, "%s%lf", dummy, &cfl);
   fscanf(fp, "%s%lf", dummy, &finaltime);
   fscanf(fp, "%s%d", dummy, &NC);
   fscanf(fp, "%s%d", dummy, &PORD);
   fscanf(fp, "%s%d", dummy, &NPLT);
   fscanf(fp, "%s%d", dummy, &FLUX);
   fscanf(fp, "%s%lf", dummy, &Mfact);
   fscanf(fp, "%s%lf%lf", dummy, &xmin, &xmax);
   fscanf(fp, "%s%lf", dummy, &XS);
   fscanf(fp, "%s%lf%lf%lf", dummy, &d_left, &u_left, &p_left);
   fscanf(fp, "%s%lf%lf%lf", dummy, &d_right, &u_right, &p_right);
   fclose(fp);
}

/* Initial condition for Burgers equation */
REAL InitCondBurger(REAL x)
{
   if(x < 0.5)
      return 1.0;
   else
      return 0.0;
}

/* Initial condition for Euler equation */
void InitCondEuler(REAL x, REAL * U)
{
   REAL d, u, p;

   if(x < XS) {
      d = d_left;
      u = u_left;
      p = p_left;
   }
   else {
      d = d_right;
      u = u_right;
      p = p_right;
   }

   U[0] = d;
   //U[0] = sin(M_PI*x);
   U[1] = d * u;
   U[2] = p / (GAMMA - 1.0) + 0.5 * d * u * u;
}
