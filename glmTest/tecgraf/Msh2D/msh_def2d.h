
#ifndef _MSH_DEF2D_H_
#define _MSH_DEF2D_H_

/*
** ----------------------------------------------------------------------
**
** msh_def2d.h - Header file for the 2D meshing routines. 
**
** ----------------------------------------------------------------------
**
** Created:      01-Oct-97      Joaquim B.C. Neto (3D)
** Modify:       01-Mar-98      Antonio C.O. Miranda (2D)
**
** ----------------------------------------------------------------------
**
*/

#include "msh2d.h"

#define Msh2DMalloc  malloc
#define Msh2DRealloc realloc
#define Msh2DFree    free

/*
** -----------------------------------------------------------------
** Public Types:
*/
typedef struct shape_pt
{
   double x, y;
} Sh_pt;

typedef struct _Msh2DBdryEdgeRec {
  struct  _Msh2DBdryEdgeRec  *next ; /* next edge pointer       */
  struct  _Msh2DBdryEdgeRec  *prev ; /* previous edge pointer   */
  int      verts[2] ;                /* corner vertices         */
  int      use ;                     /* new node use status     */
  double   center[2] ;               /* center point coordinate */
  double   max[2] ;                  /* max x and y             */  
  double   min[2] ;                  /* min x and y             */
  double   r[2] ;                    /* vector from v1 to v2    */
  double   nrm[2];
  double   length ;                  /* length of the triangle  */
  double   key ;                     /* key value for sort      */
}
Msh2DBdryEdgeRec, Msh2DBdryEdge, Msh2DBdryStack ;

typedef struct _Msh2DPolyEdgeRec {
    struct _Msh2DPolyEdgeRec *next ; /* next pointer (or null)  */
    struct _Msh2DPolyEdgeRec *prev ; /* prev pointer (or null)  */
    Msh2DBdryEdge            *edge ; /* address of the edge     */
}
Msh2DPolyEdgeRec, Msh2DPolyEdge, Msh2DPolyEdgeList ;

typedef struct _Msh2DAdjEdgeRec {
    struct _Msh2DAdjEdgeRec *next ; /* next pointer (or null)  */
    Msh2DBdryEdge           *edge ; /* address of the edge     */
}
Msh2DAdjEdgeRec, Msh2DAdjEdge, Msh2DAdjEdgeList ;

typedef struct _Msh2DAdjIniEdgeRec {
    struct _Msh2DAdjIniEdgeRec *next ;     /* next pointer (or null)         */
    int                         verts[2] ; /* verts of initial adjacent edge */
    double                      r[2];
}
Msh2DAdjIniEdgeRec, Msh2DAdjIniEdge, Msh2DAdjIniEdgeList ;

typedef struct _Msh2DAdjElemRec {
    struct _Msh2DAdjElemRec *next ;    /* next pointer (or null) */
    int                      elem ;    /* address of the element */
} Msh2DAdjElemRec, Msh2DAdjElem, Msh2DAdjElemList ;

typedef struct _Msh2DBdryNodeRec {
    double               coord[2] ;     /* nodal coordinates          */
    int                  active_flag ;  /* active node flag           */
    Msh2DAdjElemList     *elems ;       /* adjacent elem list         */
    Msh2DAdjEdgeList     *edges ;       /* adjacent edge list         */
    Msh2DAdjIniEdgeList  *iedges ;      /* adjacent initial edge list */
}
Msh2DBdryNodeRec, Msh2DBdryNode, Msh2DBdryNodeList ;

typedef enum {
  GDB_VERTEX,
  GDB_EDGE,
  GDB_NEWEDGE,
  GDB_DELEDGE,
  GDB_TREE
} Top;

/*
** -----------------------------------------------------------------
** Public Function:
*/

#ifdef NO_PROTOTYPE
int Msh2DInternalNodes();
double Msh2DOptimalNodes();
#else

Msh2DSizeElement *Msh2DSizeFunction; 

int Msh2DGenQuadTree( int num_org_nodes, int num_org_edges, 
                      double original_nodes[][2], int original_edges[][2],
                      int *num_gen_nodes, double **generated_nodes ) ;

int Msh2DBdryContraction ( int num_org_nodes, int  num_org_edges,
                           double original_nodes[][2], int original_edges[][2],
                           int *num_int_nodes, double **internal_nodes, 
                           int *num_gen_elements, int **generated_elements);

double Msh2DOptimalNodes( double edge_center[2], 
                         double tree_center[2], 
                         int  *tree_level ) ;
#endif

#endif
