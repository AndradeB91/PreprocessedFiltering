/*
** ---------------------------------------------------------------
**
** surf3d_main.c: Main driver to generate surface 3D mesh.
**
** ---------------------------------------------------------------
**
** Created:      Jan-2007      Antonio Miranda
**
** ---------------------------------------------------------------
**
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "surf3d_def.h"


/* --------------------------------------------------------------
** Public function: 
*/

/************************ Surf3DAdvFront *************************/
int Surf3DAdvFront ( 
int     n_pts,             /* # of points                          (in) */
double  *bdry_pts,         /* coordinates of all points            (in) */
double  *normal,           /* normal of points                     (in) */
int     bound_edge,        /* # of boundary edges                  (in) */
int     *edges,            /* vector with the edges                (in) */
void   (*mes_func) (char *),
int     *n_node,           /* counts # of pts for meshing         (out) */
double  **coords,          /* coordinate array used for meshing   (out) */
int     *n_elem,           /* number of elements generated        (out) */
int     **Conn             /* elem.connectivity list from meshing (out) */
)
{
 int status, i, num_gen_nodes, num_elems, *elems;
 double  *generated_nodes ;

 Surf3DMessFunction = mes_func;

 num_gen_nodes = 0;

 if (Surf3DSizeFunction    == NULL ||
     Surf3DIPointsFunction == NULL ||
     Surf3DSnapPtFunction  == NULL )
     return 0;

 /* 0.3 Generate Surface 3D mesh by Boundary Contraction operations */
 status = Surf3DBdryContraction (n_pts, bound_edge, (double (*)[3]) bdry_pts,           /* IN  */
                                (double (*)[3]) normal, (int (*)[2]) edges,             /* IN  */
                                &num_gen_nodes, &generated_nodes, &num_elems, &elems ); /* OUT */

 if (num_elems == 0)
   return 0;

 *n_node = num_gen_nodes;
 *coords = (double *)generated_nodes;

 *n_elem = num_elems;

 *Conn = (int *) calloc( 4*num_elems, sizeof(int) );
 for( i=0; i<num_elems; i++)
 {
   (*Conn)[i*4+0] = 3;
   (*Conn)[i*4+1] = elems[i*3+0];
   (*Conn)[i*4+2] = elems[i*3+1];
   (*Conn)[i*4+3] = elems[i*3+2];
 }

 free(elems);

 /* set the functions to NULL */
 Surf3DSizeFunction    = NULL;
 Surf3DIPointsFunction = NULL;
 Surf3DSnapPtFunction  = NULL;

 return status ;
}

/************************ Surf3DRegFunc *************************/
void Surf3DRegFunc (Surf3DSizeElement *Surf3D_size, Surf3DIdealPoints *Surf3D_IPoints,
                    Surf3DSnapPoint   *Surf3D_SnapPts)
{
  Surf3DSizeFunction    = Surf3D_size;
  Surf3DIPointsFunction = Surf3D_IPoints;
  Surf3DSnapPtFunction  = Surf3D_SnapPts;
}
