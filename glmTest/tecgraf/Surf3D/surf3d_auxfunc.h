/*
**  surf3d_auxfunc.h - This file contains auxiliar functions (AuxFuncGetElemSize,
    AuxFuncIdealPoints, AuxFuncSnapPoint) used in advance front algorithm.
*/



#ifndef _AUXFUNCTOMESHSURF3D_H_
#define _AUXFUNCTOMESHSURF3D_H_

#ifdef __cplusplus
extern "C" {
#endif

/* AuxFuncSetCurrentMesh builds internal structures and AuxFuncRelease releases them */
int  AuxFuncSetCurrentMesh (void *surf,  double max_size, int np, double *pts, 
                            int ne, int *edges);
int  AuxFuncGetNormal      (double *pts, double hint_size, double *normal); 
void AuxFuncRelease        (void);


/* Auxiliar functions:
**  AuxFuncGetElemSize - to obtain the size of new elements
**  AuxFuncIdealPoints - to obtain the position of ideal point
**  AuxFuncSnapPoint   - to obtain the position of a point onto old triangular surface
*/
void AuxFuncGetElemSize (double x, double y, double z, double *size);
int  AuxFuncIdealPoints (double pts[3], double n[3], double direction[3], double length, 
                         double ipts[3]);
void AuxFuncSnapPoint   (double pts[3], double normal[3], double box, double snap_pts[3]);


void AuxFuncDrawOctree (void);

#ifdef __cplusplus
}
#endif

#endif

