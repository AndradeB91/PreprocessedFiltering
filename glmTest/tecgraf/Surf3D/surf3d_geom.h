/*
** ----------------------------------------------------------------------
**
** surf3d_geom.h - Functions to geometric tests and transformations. 
**
** ----------------------------------------------------------------------
**
** Created:       21-Jan-2007    Antonio C.O. Miranda (surf3D)
**
** ----------------------------------------------------------------------
**
*/

#ifndef _SURF3DGEOM_H_
#define _SURF3DGEOM_H_

#include "surf3d_def.h"

double Surf3DVectLen2         (double *u);
double Surf3DVectLen          (double *u);
void   Surf3DCrossprod3d      (double *u, double *v, double *w);
void   Surf3DCrossprodNorm    (double *u, double *v, double *w);
void   Surf3DGetMatrixTransfZ (double *n, double matrix[3][3]);
void   Surf3DGetMatrixTransfX (double *n, double matrix[3][3]);
void   Surf3DTransf_xy        (double orig[3], double BaseMatrix[3][3], 
                               double p3d[3], double p2d[2]);
void   Surf3DTransf_xyz       (double orig[3], double BaseMatrix[3][3], 
                               double p3d[3], double p3d_l[3]);
void   Surf3DTransf_XYZ       (double orig[3], double BaseMatrix[3][3], 
                               double p2d[2],double p3d[3]);

int Surf3DPlaneLineIntersec_xy (double orig[3], double BaseMatrix[3][3],
                                double tri[3][3], int *n, double p2d[3][2]);

double Surf3DVectArea         (double *u, double *v);
int    Surf3DCrossEdge2d      (int I1, int I2, int J1, int J2, 
                               Surf3DBdryNode *nodes, double BaseMatrix[3][3]);
double Surf3DCrossprod2d      (double I[2], double J[2], double K[2]);
double Surf3DAreaEdge_Node    (Surf3DBdryEdge *edge, int node_indx, Surf3DBdryNode *nodes);
double Surf3DVectAngle        (double *u, double *v);
double Surf3DMinDist2Lines    (double L1i[3], double L1j[3], double L2i[3], double L2j[3]);      
int    Surf3DPointInsideElem  (double pts_tri[3][3], double pts_in[3], double base[3][3]);


#endif

