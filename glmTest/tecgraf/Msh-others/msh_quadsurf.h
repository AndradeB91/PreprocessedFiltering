/*
** ----------------------------------------------------------------------
**
** msh_quadsurf.h - Header file for quadtree meshing routines. 
**
** ----------------------------------------------------------------------
**
** Created:      Fev-2008      Antonio C.O. Miranda
**
** ----------------------------------------------------------------------
**
*/

#ifndef _MSH_QUAD_SURF_H_
#define _MSH_QUAD_SURF_H_


/*
** -----------------------------------------------------------------
** Public Functions:
*/

/* Generate a equivalente quadtree for surface mesh generation */
int MshSurfGenQuadTree  ( int num_org_nodes, 
                          int num_org_edges, 
                          double original_nodes[][2], 
                          int original_edges[][2],
                          void (*f_surf)(double, double, double *)) ;

/* set a maximum size of cells */
void MshSurfSetMaxQuadSize (double max_size);

/* release quadtree created with previous function */
void MshSurfQuadFreeAll();

/* Obtain an optimal node from quadtree */
double MshSurfOptimalNodes ( double edge_center[2], 
                             double tree_center[2], 
                             int  *tree_level ) ;

/* Obtain deviations held in quadtree */
void   MshSurfQuadTreeDeriv  ( double u, double v, double *dev, int interpol );


#endif
