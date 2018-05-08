/*
** ----------------------------------------------------------------------
**
** mshsurf3d.h - Header file for the Surf 3D meshing routine. 
**
** ----------------------------------------------------------------------
**
** Created:      18-Jun-2007      Antonio C.O. Miranda
** - This file incorporates the old mesh generator using parametric space,
**   the topological copy of the mesh, and the new surface mesh generator
**   directly in 3D (developed in the period of Cornell) 
**
** ----------------------------------------------------------------------
**
*/


#ifndef _SURF3DMESH_H_
#define _SURF3DMESH_H_

#include "deflib.h"

#ifdef __cplusplus
extern "C" {
#endif


/* auxiliar functions types*/
typedef MSHSURF_API void MshSurfSizeElement (double u, double v, double *size);
typedef MSHSURF_API void Surf3DMessFunc      (char *message);


/*
** ----------------------------------------------------------------------
** Begin surface mesh generator using parametric 2D space
** ----------------------------------------------------------------------
*/

/*
** The MshSurfEdge function is based on paper bellow:
** ----------------------------------------------------------------------
** MIRANDA, A. C. O., MARTHA, L. F., 
** Mesh Generation on High-Curvature Surfaces Based on a Background Quadtree Structure 
** In: 11th International Meshing Roundtable, 2002, Ithaca, New York. 11th 
** International Meshing Roundtable. Albuquerque, New Mexico: Sandia National 
** Laboratories, 2002. v.11. p.333 - 341.
** ----------------------------------------------------------------------
*/


MSHSURF_API int MshSurfEdge ( 
int     n_pts,         /* # of points                          (in) */
double  *bdry_pts,     /* coordinates of all points            (in) */
int     bound_edge,    /* # of boundary free edges             (in) */
int     inter_edge,    /* # of internal free edges             (in) */
int     *edges,        /* edge vector (i0,j0; i1, j1; ...)     (in) */
int     type_mesh,     /* 3 -> T3;  6 -> T6                    (in) */ 
int     *n_node,       /* # of pts in the mesh                (out) */
double  **coords,      /* coordinate array of the mesh        (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn,        /* elem.connectivity list of the mesh  (out) */
void    (*f_surf)      (double, double, double *) 
                       /* function to compute the deviations   (in) */
);

/* function to obtain the element size in a region. If NULL,
it uses an Quadtree based on boundary edges              */
MSHSURF_API void MshSurfRegFunc (MshSurfSizeElement *mshsurf_size);


/* This function set options to surface mesh generator
**
** nf               - number of flags 
** flags[0] - (0,1) - refine mesh considering the max size of edge.
**                    Default is 1.
** flags[1] - (0,1) - refine mesh considering the max size of element
**                    given in param variable. Default is 0.
** flags[2] - (0,1) - refine locations of high curvature given a max
**                    angle in param variable. Default is 1.
** flags[3] - (0,1) - refine again locations of high curvature. Default is 0.
**
** np               - number of parameters
** param[0]         - max size of element. Default is the max size of edge.
** param[1]         - max angle between elements (radius). Default is PI/12.
**
*/
MSHSURF_API void MshSurfOptions (int nf, int *flags, int np, double *param);



/*
** ----------------------------------------------------------------------
** Begin copy topological mesh - source to a target mesh
** ----------------------------------------------------------------------
*/


/* Create copy mesh */
/* These function create a copy mesh under another mesh.
** Both meshes must have the same number of edges on the boundary,
** and holes and internal restricions are not considered. A base
** point in both mesh have to be informed. Using these points 
** the copy mesh is created.
** Obs: The input mesh must contain only triangular e quadrilateral
**      elements
*/

/* First Function - Set the source mesh */
MSHSURF_API void MshSurfSetSourceMesh (
double  base_pts[3],  /* point on boundary where inits surf copy (in) */
void    *surfmesh,    /* surface mesh using topology lib or the 
                      inserting a mesh as bellow              (in) */
int     n_node,       /* # of pts in the mesh                    (in) */
double  *coords,      /* coordinate array of the mesh            (in) */
int     n_elem,       /* number of elements generated            (in) */
int     *Conn         /* elem.connectivity list of the mesh      (in) */
); 

/* Second Function - Set the target mesh */
MSHSURF_API void MshSurfSetTargetMesh (
double  base_pts[3],  /* point on boundary where inits surf copy (in) */
void    *surfmesh,    /* surface mesh using topology lib or the 
                      inserting a mesh as bellow              (in) */
int     n_node,       /* # of pts in the mesh                    (in) */
double  *coords,      /* coordinate array of the mesh            (in) */
int     n_elem,       /* number of elements generated            (in) */
int     *Conn         /* elem.connectivity list of the mesh      (in) */
); 

/* Third Function - Get the copy mesh
** This function consider two cases: a boundary nodes as input or
**                                   the boundary of target surface
*/
MSHSURF_API int MshSurfGetCopyMesh(
int     n_bound_nodes, /* # of boundary nodes (optional or 0)  (in) */
double  *bound_nodes,  /* boundary nodes (optional or NULL)    (in) */
int     *n_node,       /* # of pts in the mesh                (out) */
double  **coords,      /* coordinate array of the mesh        (out) */
double  **normal,      /* normal of the coordinate            (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn         /* elem.connectivity list of the mesh  (out) */
);

/*
** ----------------------------------------------------------------------
** Begin 3D Mesh to Surface Mesh
** ----------------------------------------------------------------------
*/

/*
** This algorithm needs to call two functions:
** 1 - MshSurfSetTargetMesh - target support mesh (the same function above)
** 2 - MshSurfCreate3dBil  - create 3d bilinear mesh and move to support mesh
*/
MSHSURF_API int MshSurfCreate3dBil 
(
 int     np,            /* # of points                          (in) */
 double  *bry,          /* coordinates of all points            (in) */
 int     elem_type,     /* element type:  T3 (3), T6 (6), Q4 (4) ou  Q8 (8) (in) */
 int     diagtype,      /* Triangle option (option for cell diagonal) (in)
                        = 1 --> diagonal oriented to right direction
                        = 2 --> diagonal oriented to left direction
                        = 3 --> union jack alternation
                        = 4 --> optimum diagonal (smallest of two possible) */
int     *n_node,       /* # of pts in the mesh                (out) */
double  **coords,      /* coordinate array of the mesh        (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn         /* elem.connectivity list of the mesh  (out) */
);


/*
** ----------------------------------------------------------------------
** Begin surface mesh generator directly in 3D space using a support surface
** ----------------------------------------------------------------------
*/

/* First Function - Set the support surface. */

/* The support surface is a mesh */
MSHSURF_API void MshSurfSetSupportMesh (
void    *surfmesh,    /* surface mesh using topology lib or the 
                         inserting a mesh as bellow              (in) */
int     n_node,       /* # of pts in the mesh                    (in) */
double  *coords,      /* coordinate array of the mesh            (in) */
int     n_elem,       /* number of elements generated            (in) */
int     *Conn         /* elem.connectivity list of the mesh      (in) */
); 



/* The support surface is an analytical surface */
MSHSURF_API void MshSurfSetSupportAnalytical (
/* this function obtains uv coord. from xyz coordinate    */
  int (*getUVfromXYZ)   (double xyz[3], double tol, double uv[2]),
/* this function obtains xyz coord. from xy coordinate    */
  int (*getXYZfromUV)   (double uv[2], double xyz[3]),
/* this function obtains the gradiente vector in uv point */
  int (*getGradVectors) (double uv[2], double Gu[3], double Gv[3])
); 



/* Second Function - Create a new mesh using the currente support mesh */
MSHSURF_API int MshSurf3D ( 
int     n_pts,            /* # of points                          (in) */
double  *bdry_pts,        /* coordinates of all points            (in) */
int     bound_edge,       /* # of boundary edges                  (in) */
int     inter_edge,       /* # of internal free edges             (in) */
int     *edges,           /* vector with the edges                (in) */
double  max_elm_size,     /* maximum size of elements             (in) */
int     curvature,        /* compute curvature ?                  (in) */
Surf3DMessFunc *mes_func, /* message function                     (in) */
int     *n_node,          /* counts # of pts for meshing         (out) */
double  **coords,         /* coordinate array used for meshing   (out) */
int     *n_elem,          /* number of elements generated        (out) */
int     **Conn            /* elem.connectivity list from meshing (out) */
);

/* Third Function - Release current support mesh */
MSHSURF_API void MshSurfReleaseSupportMesh (void);



/*
** ----------------------------------------------------------------------
** Begin 2D to Surface Mesh
** ----------------------------------------------------------------------
*/

/*
** This algorithm needs to call two functions:
** 1 - MshSurfSetSupport2d - target support mesh
** 2 - MshSurfCreate2d     - create 2d mesh and move to support mesh
*/

/* First Function - Set the support surface. */

/*  target support mesh */
MSHSURF_API void MshSurfSetSupport2d (
int     n_node,       /* # of pts in the mesh                    (in) */
double  *coords,      /* coordinate array of the mesh            (in) */
int     n_elem,       /* number of elements generated            (in) */
int     *Conn         /* elem.connectivity list of the mesh      (in) */
); 


/* Second Function - Create a new mesh using the current support mesh */
/*                   This function use boundary from support mesh     */
MSHSURF_API int MshSurfCreate2d (
int     alg_2d,        /* 2d algorithm                         (in)
                          0 - quadrilateral mapping
                          1 - boundary contraction
                          2 - quadrilateral with quadtree
                          3 - quadrilateral seam                    */
int     *n_node,       /* # of pts in the mesh                (out) */
double  **coords,      /* coordinate array of the mesh        (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn         /* elem.connectivity list of the mesh  (out) */
);

/* Second Function - Create a new mesh using the current support mesh */
/*                   This function use boundary from edges            */
MSHSURF_API int MshSurfCreate2dEdge (
int     alg_2d,        /* 2d algorithm                         (in)
                          0 - boundary contraction
                          1 - quadrilateral seam                    */
int     n_pts,         /* # of points                          (in) */
double  *bdry_pts,     /* coordinates of all points            (in) */
int     bound_edge,    /* # of boundary free edges             (in) */
int     *edges,        /* edge vector (i0,j0; i1, j1; ...)     (in) */
int     *n_node,       /* # of pts in the mesh                (out) */
double  **coords,      /* coordinate array of the mesh        (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn         /* elem.connectivity list of the mesh  (out) */
);


/* Mesh Algorithm
   - Create surface mesh using Msh2DEdge, but the boundary points is used
     to create a plane. If the option to create internal points, these
     internal points is smoothed to 3D */
MSHSURF_API int MshSurfEdge2D ( 
int     n_pts,            /* # of points                          (in) */
double  *bdry_pts,        /* coordinates of all points            (in) */
int     bound_edge,       /* # of boundary edges                  (in) */
int     inter_edge,       /* # of internal free edges             (in) */
int     *edges,           /* vector with the edges                (in) */
int     internal_pts,     /* create internal points?              (in) */
int     *n_node,          /* counts # of pts for meshing         (out) */
double  **coords,         /* coordinate array used for meshing   (out) */
int     *n_elem,          /* number of elements generated        (out) */
int     **Conn            /* elem.connectivity list from meshing (out) */
);


/*
** ----------------------------------------------------------------------
** Begin Surface Mesh using TEMPLATES
** ----------------------------------------------------------------------
*/

/* Optional Function: MshSurfSetTargetMesh - target support mesh (the same 
 *                    function above). Set the support surface. The surface 
 *                    is used to project 3D points to support surface
*/

/* Mesh Algorithm
   - Create quadrilateral mesh in 3D surface. If the support mesh is given, all
     points is projected to this surface. 
*/
MSHSURF_API int MshSurfTemplate (
int     n_sides,          /* # of sides (2, 3, 4)                 (in) */
int     subdvision[4],    /* subdvision of each side              (in) */
double  *bdry_pts,        /* coordinates of all points            (in) */
int     *n_node,          /* counts # of pts for meshing         (out) */
double  **coords,         /* coordinate array used for meshing   (out) */
int     *n_elem,          /* number of elements generated        (out) */
int     **Conn            /* elem.connectivity list from meshing (out) */
);


/*
** ----------------------------------------------------------------------
** Smoothing Algorithm
** ----------------------------------------------------------------------
*/

MSHSURF_API int MshSurfSmooth (
int     n_node,       /* # of pts in the mesh                    (in) */
double  *coords,      /* coordinate array of the mesh        (in/out) */
int     n_elem,       /* number of elements generated            (in) */
int     *Conn,        /* elem.connectivity list of the mesh      (in) */
int     algorithm,    /* algorithm: 1-laplacian, 2 - mod. laplacian, 3 - taubin (in) */
int     num_steps,    /* number of steps                         (in) */
int     n_param,      /* number of parameters                    (in) */
double  *param        /* parameters                              (in) */
);

/** Adiciona/Obtêm parametros para os geradores de malha
    @param index   (in) indice do parametro 
                        0 - faz a suavização dos pontos quando a malha é
                            projeção 2D para uma superfície
    @param value   (in) valor do parametro de parametros
                        value = 0.0 (false) e value = 1.0 (true) 
*/
MSHSURF_API void   MshSurfSetParams (int index, double value);
MSHSURF_API double MshSurfGetParams (int index);


#ifdef __cplusplus
}
#endif

#endif
