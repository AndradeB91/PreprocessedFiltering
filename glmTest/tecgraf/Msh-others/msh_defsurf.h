/*
** ----------------------------------------------------------------------
**
** msh_defsurf.h - Header file for the surface meshing routines. 
**
** ----------------------------------------------------------------------
**
** Created:      Fev-2000      Antonio C.O. Miranda
**
** ----------------------------------------------------------------------
**
*/

#ifndef _MSH_DEF_SURF_H_
#define _MSH_DEF_SURF_H_

#define Msh2DMalloc  malloc
#define Msh2DRealloc realloc
#define Msh2DFree    free


/* functions to get flags and parameters */
void   MshSurfSetDefaultFlagsParams    ( void );
int    MshSurfRefineByMaxEdgeSize_ask  ( void );
int    MshSurfRefineByMaxSetSize_ask   ( void );
int    MshSurfRefineCurvature_ask      ( void );
int    MshSurfRefineCurvature2_ask     ( void );
double MshSurfGetMaxSizeElement        ( void );
double MshSurfGetMaxAngle              ( void );

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
  double   length ;                  /* length of the edge      */
  double   key ;                     /* key value for sort      */
  double   real_length;              /* real length of the edge */
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
    double               E;
    double               F;
    double               G;
}
Msh2DBdryNodeRec, Msh2DBdryNode, Msh2DBdryNodeList ;


#endif
