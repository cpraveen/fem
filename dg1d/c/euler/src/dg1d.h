#ifndef _DG1D

#define _DG1D  1

#define LINCON 1
#define BURGER 2
#define EULER  3

REAL mass0, mass1[2][2], mass2[3][3], mass3[4][4], mass4[5][5], mass5[6][6];
REAL cfl, dt, finaltime;
REAL XS;                        /* Shock position */
REAL xmin, xmax;
REAL d_left, u_left, p_left;
REAL d_right, u_right, p_right;
REAL Mfact;

typedef struct
{
   REAL x, xl, xr, h, *xg;
   UINT p, ng;
   REAL **Un, **Uo, **Re;
} CELL;

void GaussInit ();
void GaussPoints (CELL *);
REAL ShapeFun (REAL, CELL *, UINT);
REAL ShapeFunDeriv (REAL, CELL *, UINT);

void Uvect (CELL * cell, REAL x, REAL * U);
void EulerFlux (REAL * U, REAL * flux);
void RoeFlux (REAL * Ul, REAL * Ur, REAL * flux);
void LFFlux (REAL * Ul, REAL * Ur, REAL * flux);
void ECUSPFlux (REAL * Ul, REAL * Ur, REAL * flux);
void HLLCFlux (REAL * Ul, REAL * Ur, REAL * flux);
void AUSMDVFlux (REAL * Ul, REAL * Ur, REAL * flux);
void LFCFlux (REAL * Ul, REAL * Ur, REAL * flux);

void EigMat (REAL *, REAL[][3], REAL[][3]);
void Multi (REAL[][3], REAL *);

#endif
