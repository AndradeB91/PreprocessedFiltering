/*
**  surf3d_auxfunc.h - This file contains auxiliar functions (AuxFuncGetElemSize,
    AuxFuncIdealPoints, AuxFuncSnapPoint) used in advance front algorithm.
*/



#ifndef _AUXFUNCTOSURF3D_H_
#define _AUXFUNCTOSURF3D_H_

#ifdef __cplusplus
extern "C" {
#endif


/* AuxFuncSetCurrentSurface - define the current surface using feedback functions */
int  AuxFuncSetCurrentSurface
(
/* this function obtains uv coord. from xyz coordinate    */
  int (*getUVfromXYZ)   (double xyz[3], double tol, double uv[2]),
/* this function obtains xyz coord. from xy coordinate    */
  int (*getXYZfromUV)   (double uv[2], double xyz[3]),
/* this function obtains the gradiente vector in uv point */
  int (*getGradVectors) (double uv[2], double Gu[3], double Gv[3])
);

/* AuxFuncSetCurrentMesh builds internal structures and AuxFuncRelease releases them */
int  AuxFuncSetCurrentBoundary (double max_size, int np, double *pts, 
                                int ne, int *edges);

int  AuxFuncGetNormal2      (double *pts, double hint_size, double *normal); 

void AuxFuncRelease2        (void);


/* Auxiliar functions:
**  AuxFuncGetElemSize - to obtain the size of new elements
**  AuxFuncIdealPoints - to obtain the position of ideal point
**  AuxFuncSnapPoint   - to obtain the position of a point onto old triangular surface
*/
void AuxFuncGetElemSize2 (double x, double y, double z, double *size);
int  AuxFuncIdealPoints2 (double pts[3], double n[3], double direction[3], double length, 
                         double ipts[3]);
void AuxFuncSnapPoint2  (double pts[3], double normal[3], double box, double snap_pts[3]);


void AuxFuncDrawOctree2 (void);

#ifdef __cplusplus
}
#endif

#endif

