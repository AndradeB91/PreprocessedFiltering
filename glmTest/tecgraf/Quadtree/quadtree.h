/*
** ---------------------------------------------------------------------------
**
** quadtree.h  -  This file contains defined types and prototypes  of the 
**                functions used by quatree generation modules, that is,
**                quadtree.c,  quadint.c,  quadtri.c  and quadqua.c  .
**
** ---------------------------------------------------------------------------
**
** Version     0-001
**
** Created:   27-Mai-93   -   Joaquim Bento Cavalcante Neto ( Joaquim )
**                            Stolen modules quadtree.h from Marcelo Tilio
**
** Modified:  28-Jul-94	  -   Joaquim Bento Cavalcante Neto
**			      Introduced a new field edepth in tree_node_t 
**			      structure that will be used in error refinement.
**
** Modified:  29-May-96	  -   Eduardo Setton Sampaio Silveira &&
**                            Joaquim Bento Cavalcante Neto
**                            Introduced function quadmesh_bound to allow 
**                            generation of quadrilateral mesh. Included
**                            a new field type in list structure as the 
**                            enumeration QUAD_BOUNDARY and INTERNAL for this
**                            type.
**
** Supervised by:            Luiz Fernando Martha
**
** ---------------------------------------------------------------------------
*/

#ifndef QUADTREE_H
#define QUADTREE_H

#include "deflib2d.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
** ---------------------------------------------------------------------------
** Coordinate typedefs:
*/

typedef MSH2D_API struct 
{ 
 double x, y;
} point;

typedef MSH2D_API point GeoPnt;

typedef MSH2D_API struct 
{ 
 double x, y;
} vector;

typedef MSH2D_API struct vertex_
{ 
 int	bdrypt;
 double  x;
 double  y;
} vertex_t;

typedef MSH2D_API struct 
{ 
 int ex, ey;
} new_coord;

/*
** ---------------------------------------------------------------------------
** Box typedefs:
*/

typedef MSH2D_API struct box2D_
{
 double xmin, ymin, xmax, ymax;
} box2D_t;                    /* bidimensional box */

typedef MSH2D_API struct win2D_
{
 int xmin, ymin, xmax, ymax;
} win2D_t;                    /* bidimensional windown( quadtree )*/

/*
** ---------------------------------------------------------------------------
** Geometric typedefs:
*/

typedef MSH2D_API enum   elistt_
{
 QUAD_BOUNDARY,
 QUAD_INTERNAL
} elist;

typedef MSH2D_API struct face_
{
 struct face_     *next;
 int              n;
 vertex_t         **vt;
} face_t;

typedef MSH2D_API struct listt_
{
 struct listt_   *next;
 int             vi;
 int             vj;
 int             type;
} list;

typedef MSH2D_API struct adjac_
{
 struct adjac_    *next;
 list             *p;    
} adjnode; 
  
typedef MSH2D_API struct adj_
{
 struct adj_  *next;
 int          elm;
} adj_t;

/*
** ---------------------------------------------------------------------------
** Fem typedefs:
*/

typedef MSH2D_API struct mesh_vertex_
{	       
  struct mesh_vertex_  *next;
  int                  x, y;
  int                  nno;
  int                  status;
  int                  bdrypt; 
  adj_t                *fadj;  
} mesh_vertex_t;

typedef MSH2D_API struct mesh_elem_
{
 struct mesh_elem_     *next;
 int                   n;
 int                   nel;
 mesh_vertex_t         **conect;
} mesh_elem_t;

/*
** ---------------------------------------------------------------------------
** Tree typedefs:
*/

typedef MSH2D_API enum
{
 NODE_UNDEF,
 NODE_VERTEX,
 NODE_BOUND,
 NODE_INTERIOR,
 NODE_OUT
} tree_code_t;                /* tree node code */

typedef MSH2D_API enum
{
 NO_SIDE = -1,
 BOTTOM  =  0,
 RIGHT   =  1,
 TOP     =  2,
 LEFT    =  3,
 CORNER0 =  4,
 CORNER1 =  5,
 CORNER2 =  6,
 CORNER3 =  7
} side_t;                      /* tree node side */

typedef MSH2D_API struct tree_vertex_
{		     
 struct tree_vertex_ *next;
 int                 id;
 new_coord           coord;
} tree_vertex_t;              
      
typedef MSH2D_API struct tree_node_
{
 tree_code_t        code;
 unsigned           depth;
 unsigned	        edepth;
 new_coord          coord;
 mesh_elem_t        *elmhead;
 tree_vertex_t      *vhead;
 struct tree_node_  **child;
} tree_node_t;                /* tree node */

typedef MSH2D_API struct neighbour_
{
 struct neighbour_ *next;
 tree_node_t       *node;
} neighbour_t;

typedef MSH2D_API struct nodeadj_
{ 
 struct nodeadj_  *next;
 side_t           side;
 tree_node_t      *node;
} nodeadj_t;
 
typedef MSH2D_API struct tree_
{
 box2D_t     box;
 tree_node_t *root;
} tree_t;                     /* tree node root */


/*
** ---------------------------------------------------------------------------
** Public Functions:
*/

#ifdef	QUADTREE_C

/*
MSH2D_API int                 msh_quadtree  ( int, int *, double *, int, int *,
												          double **, int *, int ** ); */
MSH2D_API void                param_form    ( double , double , new_coord * );
MSH2D_API void                invparam_form ( int, int, vertex_t * );
MSH2D_API neighbour_t         *locate_quad  ( new_coord *, tree_node_t * );
MSH2D_API neighbour_t         *find_adj_node( side_t, tree_node_t * );
MSH2D_API int                 pow_2         ( int  );
MSH2D_API double              cross_product ( vector *, vector * );
MSH2D_API double              pseudo_ang    ( vector *, vector * );
MSH2D_API double              dot_product   ( vector *, vector * );
MSH2D_API void                add_bdrypts_quad ( tree_node_t *, list ** );
#else

extern MSH2D_API int          msh_quadtree  ( int, int *, double *, int, int *,
														    double **, int *, int ** );
extern MSH2D_API void         param_form    ( double , double , new_coord * );
extern MSH2D_API void         invparam_form ( int, int, vertex_t * );
extern MSH2D_API neighbour_t  *locate_quad  ( new_coord *, tree_node_t * );
extern MSH2D_API neighbour_t  *find_adj_node( side_t, tree_node_t * );
extern MSH2D_API int          pow_2         ( int  );
extern MSH2D_API double       cross_product ( vector *, vector * );
extern MSH2D_API double       pseudo_ang    ( vector *, vector * );
extern MSH2D_API double       dot_product   ( vector *, vector * );
extern MSH2D_API void         add_bdrypts_quad ( tree_node_t *, list ** );

#endif

#ifdef	QUADINT_C

MSH2D_API void                trimesh_int   ( tree_node_t *, int, int *, int *,
												          int *, mesh_vertex_t ** );
MSH2D_API void                quadmesh_int  ( tree_node_t *, int, int *, int *,
															 int *, mesh_vertex_t ** );
MSH2D_API int                 get_nnoint    ( void );
MSH2D_API int                 get_nnogeomint( void );

#else

extern MSH2D_API void         trimesh_int   ( tree_node_t *, int, int *, int *,
												          int *, mesh_vertex_t ** );
extern MSH2D_API void         quadmesh_int  ( tree_node_t *, int, int *, int *,
													       int *, mesh_vertex_t ** );
extern MSH2D_API int          get_nnoint    ( void );
extern MSH2D_API int          get_nnogeomint( void );

#endif

#ifdef	QUADTRI_C

MSH2D_API int                 trimesh_bound ( tree_node_t *, face_t *, mesh_vertex_t  *, 
											             int, int, int *, mesh_vertex_t **,
															 mesh_elem_t ** );

#else

extern MSH2D_API int          trimesh_bound ( tree_node_t *, face_t *, mesh_vertex_t  *, 
						                            int, int, int *, mesh_vertex_t **,
								                      mesh_elem_t ** );

#endif

#ifdef  QUADQUA_C
 
MSH2D_API int                 quadmesh_bound( tree_node_t *, face_t *, mesh_vertex_t  *,
														    int *, int, int *, mesh_vertex_t **,
															 mesh_elem_t ** );
									    
#else
		     
extern MSH2D_API int          quadmesh_bound( tree_node_t *, face_t *, mesh_vertex_t  *,
						                            int *, int, int *, mesh_vertex_t **,
	 										             mesh_elem_t ** );

#endif

#ifdef __cplusplus
}
#endif

#endif

