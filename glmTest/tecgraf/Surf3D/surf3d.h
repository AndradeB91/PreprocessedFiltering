/*
** ----------------------------------------------------------------------
**
** surf3d.h - Header file for the Surf 3D meshing routine. 
**
** ----------------------------------------------------------------------
**
** Created:      14-Jan-2007      Antonio C.O. Miranda
**
** ----------------------------------------------------------------------
**
*/


#ifndef _SURFDIR3DMESH_H_
#define _SURFDIR3DMESH_H_

/* #include "../mshsurf3d.h" */

/* register auxiliar functions */
typedef void Surf3DSizeElement (double x, double y, double z, double *size);
typedef int  Surf3DIdealPoints (double pts[3], double n[3], double direction[3], 
                                double length, double ipts[3]);
typedef void Surf3DSnapPoint   (double pts[3], double normal[3], double box, double surfpts[3]);

void Surf3DRegFunc (
Surf3DSizeElement *Surf3D_size,      /* This function obtains the element size in a region. */
Surf3DIdealPoints *Surf3D_IPoints,   /* This function obtains one points on surface from a 
                                        initial direction and position  */
Surf3DSnapPoint   *Surf3D_SnapPts    /* This function obtains the real position of the point 
                                        on surface */ 
);


int Surf3DAdvFront ( 
int     n_pts,            /* # of points                          (in) */
double  *bdry_pts,        /* coordinates of all points            (in) */
double  *normal,          /* normal of points                     (in) */
int     bound_edge,       /* # of boundary edges                  (in) */
int     *edges,           /* vector with the edges                (in) */
void    (*mes_func) (char *), /* message function                     (in) */
int     *n_node,          /* counts # of pts for meshing         (out) */
double  **coords,         /* coordinate array used for meshing   (out) */
int     *n_elem,          /* number of elements generated        (out) */
int     **Conn            /* elem.connectivity list from meshing (out) */
);


void Surf3DDraw (void *dummy); 

#endif

