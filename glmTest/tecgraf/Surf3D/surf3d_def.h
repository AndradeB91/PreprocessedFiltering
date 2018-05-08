/*
** ----------------------------------------------------------------------
**
** surf3d_def.h - Header file for the surf3D meshing routines.
**
** ----------------------------------------------------------------------
**
** Created:      01-Oct-97      Joaquim B.C. Neto (3D)
** Modify:       01-Mar-98      Antonio C.O. Miranda (2D)
** Modify:       14-Jan-2007    Antonio C.O. Miranda (surf3D)
**
** ----------------------------------------------------------------------
**
*/

#ifndef _SURF3DDEFS_H_
#define _SURF3DDEFS_H_


#define Surf3DMalloc  malloc
#define Surf3DRealloc realloc
#define Surf3DFree    free

#include "surf3d.h"

/*
** -----------------------------------------------------------------
** Public Types:
*/


typedef struct _Surf3DBdryEdgeRec
{
  struct  _Surf3DBdryEdgeRec  *next ; /* next edge pointer       */
  struct  _Surf3DBdryEdgeRec  *prev ; /* previous edge pointer   */
  int      verts[2] ;                /* corner vertices         */
  int      use ;                     /* new node use status     */
  int      layer ;                   /* markos - layer in which the edge belongs (0 is the input front) */
  double   center[3] ;               /* center point coordinate */
  double   max[3] ;                  /* max x and y             */
  double   min[3] ;                  /* min x and y             */
  double   r[3] ;                    /* vector from v1 to v2    */
  double   nrm[3];
  double   length ;                  /* length of the triangle  */
  double   key ;                     /* key value for sort      */
}
Surf3DBdryEdgeRec, Surf3DBdryEdge, Surf3DBdryStack ;

typedef struct _Surf3DPolyEdgeRec
{
  struct _Surf3DPolyEdgeRec *next ; /* next pointer (or null)  */
  struct _Surf3DPolyEdgeRec *prev ; /* prev pointer (or null)  */
  Surf3DBdryEdge            *edge ; /* address of the edge     */
}
Surf3DPolyEdgeRec, Surf3DPolyEdge, Surf3DPolyEdgeList ;

typedef struct _Surf3DAdjEdgeRec
{
  struct _Surf3DAdjEdgeRec *next ; /* next pointer (or null)  */
  Surf3DBdryEdge           *edge ; /* address of the edge     */
}
Surf3DAdjEdgeRec, Surf3DAdjEdge, Surf3DAdjEdgeList ;

typedef struct _Surf3DAdjIniEdgeRec
{
  struct _Surf3DAdjIniEdgeRec *next ;     /* next pointer (or null)         */
  int                         verts[2] ; /* verts of initial adjacent edge */
  double                      nrm[3];
}
Surf3DAdjIniEdgeRec, Surf3DAdjIniEdge, Surf3DAdjIniEdgeList ;

typedef struct _Surf3DElemRec
{
  int       conn[3];               /* connectivity         */
  double    normal[3];             /* normal               */
  int       del;                   /* 1 - deleted element  */
} Surf3DElemRec, Surf3DElem, Surf3DElemList ;

typedef struct _Surf3DAdjElemRec
{
  struct _Surf3DAdjElemRec *next;    /* next pointer (or null) */
  int                       id;      /* address of the element */
} Surf3DAdjElemRec, Surf3DAdjElem, Surf3DAdjElemList ;

typedef struct _Surf3DBdryNodeRec
{
  double               coord[3] ;      /* nodal coordinates          */
  double               normal[3] ;     /* normal                     */
  int                  active_flag ;   /* active node flag           */
  Surf3DAdjElemList     *elems ;       /* adjacent elem list         */
  Surf3DAdjEdgeList     *edges ;       /* adjacent edge list         */
  Surf3DAdjIniEdgeList  *iedges ;      /* adjacent initial edge list */
}
Surf3DBdryNodeRec, Surf3DBdryNode, Surf3DBdryNodeList ;



/*
** -----------------------------------------------------------------
** Public Variables:
*/

extern  Surf3DSizeElement   *Surf3DSizeFunction;
extern  Surf3DIdealPoints   *Surf3DIPointsFunction;
extern  Surf3DSnapPoint     *Surf3DSnapPtFunction;
extern  void (*Surf3DMessFunction) (char *);


/*
** -----------------------------------------------------------------
** Public Function:
*/


int Surf3DBdryContraction ( int num_org_nodes, int  num_org_edges,
                           double original_nodes[][3], double normal[][3],
                           int original_edges[][2],
                           int *num_int_nodes, double **internal_nodes,
                           int *num_gen_elements, int **generated_elements);

#endif

