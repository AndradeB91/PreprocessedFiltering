/** mshq2d.h - header files
* ------------------------------------------------------------------
*/

#ifndef MSHQ2DSEAM_H
#define MSHQ2DSEAM_H

int Msh2DSeamGeneration (
int     n_pts,         /* # of points                          (in) */
double  *bdry_pts,     /* coordinates of all points            (in) */
int     bound_edge,    /* # of boundary edges                  (in) */
int     inter_edge,    /* # of internal free edges             (in) */
int     *edges_bound,  /* vector with the edges                (in) */
int     type_mesh,     /*********************************************/ 
int     *n_node,       /* counts # of pts for meshing         (out) */
double  **coords,      /* coordinate array used for meshing   (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn         /* elem.connectivity list from meshing (out) */
);

int Msh2DSeamTest (
int     n_pts,         /* # of points                          (in) */
double  *bdry_pts,     /* coordinates of all points            (in) */
int     bound_edge,    /* # of boundary edges                  (in) */
int     inter_edge,    /* # of internal free edges             (in) */
int     *edges_bound,  /* vector with the edges                (in) */
int     type_mesh,     /*********************************************/ 
int     *n_node,       /* counts # of pts for meshing         (out) */
double  **coords,      /* coordinate array used for meshing   (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn         /* elem.connectivity list from meshing (out) */
);

#endif
