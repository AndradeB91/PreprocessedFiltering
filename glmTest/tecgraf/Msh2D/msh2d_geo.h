/*
** ----------------------------------------------------------------------
**
** msh2d_geo.h - Header file for auxiliar 2D meshing routine. 
**
** ----------------------------------------------------------------------
**
** Created:      June-2010      Antonio C.O. Miranda
**
** ----------------------------------------------------------------------
*/


#ifndef _MSH2D_GEO_H_
#define _MSH2D_GEO_H_

// basic operations
double Msh2DCrossProd     (double *u, double *v, double *w);
void   Msh2DDiffVector    (double *a, double *b, double *c);
double Msh2DCompr         (double *u);
double Msh2DSquareCompr   (double *u);


// compute 2d area
double Msh2DArea2D        (int nNodes, double *boundpts2d);

// find border stones in a loop composed by pts
void   Msh2DFindBorders   (int nNodes, double *pts, int *nbordes, int **borders );

// angle between "a" point
double Msh2DAngleAboutPt  (double a[2], double b[2], double c[2]);

// get bilinear Cornes if possible
int Msh2DAutomaticBilinearCornes  (int np, double *pts, int *m, int *n);
int Msh2DAutomaticTrilinearCornes (int np, double *pts, int *m);


#endif

