/*
** ----------------------------------------------------------------------
**
** mshsurf_geo.h - Header file for auxiliar 3D meshing routine. 
**
** ----------------------------------------------------------------------
**
** Created:      April-2010      Antonio C.O. Miranda
**
** ----------------------------------------------------------------------
*/


#ifndef _MSHSURF_GEO_H_
#define _MSHSURF_GEO_H_

// basic operations
void   MshSurfCrossProd     (double *u, double *v, double *w);
void   MshSurfDiffVector    (double *a, double *b, double *c);
void   MshSurfCrossProdNorm (double *u, double *v, double *w);
double MshSurfCompr         (double *u);
double MshSurfSquareCompr   (double *u);

// get plane
int    MshSurfLstSqrPlanEqn (int num_pts, double *pts, double plan_eqn[4]);
int    MshSurfLstSqrPlanMtx (int num_pts, double *pts, double BaseMatrix[4][4]);
void   MshSurfGenPlanEqn    (double pt0[3], double pt1[3], double pt2[3], 
                             double plan_eqn[4]);
void   MshSurfCoord_xy      (double x, double y, double z,
                             double BaseMatrix[4][4], double *p, double *h);
void   MshSurfCoord_xyz     (double x, double y, double BaseMatrix[4][4],
                             double *p);

void   MshSurfTransf_xyx    (double *p_in, double BaseMatrix[4][4], double *p_out);
void   MshSurfInvTransf_xyz (double *p_in, double BaseMatrix[4][4], double *p_out);

// compute 2d area
double MshSurfArea2D        (int nNodes, double *boundpts2d);
double MshSurfArea2DEdges   (int nNodes, double *boundpts2d, int nedge, int *edges);

// find border stones in a loop composed by pts
void   MshSurfFindBorders   (int nNodes, double *pts, int *nbordes, int **borders );

// angle between "a" point
double MshSurfAngleAboutPt  (double a[3], double b[3], double c[3]);


#endif

