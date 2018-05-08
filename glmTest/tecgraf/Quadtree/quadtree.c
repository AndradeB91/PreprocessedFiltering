/*
** ---------------------------------------------------------------------------
**
**  quadtree.c - This file contains routines for tree manipulations in the
**               quadtree finite element generation algorithm.
**
** ---------------------------------------------------------------------------
**
**  Public Mesh Generating Function :
**
** ---------------------------------------------------------------------------
**
** void Msh2DQuadBound( n_loops, loop_segs[], bdry_pts, flag_mesh, n_node,
**                     coords, n_elem, conn )
**
**   int     n_loops      - # of circuits (connect.bdry parts)         (in )
**   int     loop_segs[]  - # of segs (elem.side) on circuits        (in )
**   GeoPnt  *bdry_pts    - coordinates of points on boundary        (in )
**   int     flag_mesh    - flag for triangular or quadrilateral mesh  (in )
**   int     ref_quad     - 1 = refine quadtree, 0 == not refine       (in )
**   int     *n_node      - counts # of pts for meshing          (out)
**   GeoPnt  *coords      - coordinate array used for meshing          (out)
**   int     *n_elem      - number of elements generated         (out)
**   int     *conn        - elem.connectivity list from meshing        (out)
**
**   This routine is the main public mesh generation function, that is, it
**   receives a geometrical structure and get out a finite  element   mesh
**   structure.
**   The process is convert the geometrical input structure to a geometrical
**   local one, generate and manipulate the tree to prepare  interior  tree
**   nodes to make the templates by trimesh_int or quadmesh_int   functions
**   (depending of the input flag_mesh for triangular or quadrilateral mesh)
**   and then make the boundary contraction in the rest of the region  that
**   even has no finite elements using Delaunay tecnique, after this smooth
**   the mesh and  convert to the final output finite elements and vertexs.
**
** ---------------------------------------------------------------------------
**
** Private Tree Functions :
**
** ---------------------------------------------------------------------------
**
** void convert_geom( n_loops, loop_segs[], bdry_pts, f )
**
**   int n_loops          - # of circuits (connect.bdry parts)        (in )
**   int     loop_segs[]  - # of segs (elem.side) on circuits       (in )
**   GeoPnt  *bdry_pts    - coordinates of points on boundary       (in )
**   face_t  *f           - pointer of the local geometric structure  (out)
**
**   This routine converts the input geometrical structure to the local
**   geometric one, that will be used in this algorithm.
**
** ---------------------------------------------------------------------------
**
** Version: 0-001
**
** Created:  20-May-93    Joaquim Bento Cavalcante Neto
**    Stolen from Marcelo Tilio.
**
** Modified: 01-Jun-94          Joaquim Bento Cavalcante Neto
**   Introduced a new function div_error_tree, that divides more the tree, ba-
**   sed in the error of elements that belongs to that tree leaf.Local  varia-
**   bles elem_vector and node_vector were changed to el_vector and  no_vector
**   respectively.It was also inibed the function only_boundsize_tree, that is
**   not used anymore.
**
** Modified: 17-Jul-94    Joaquim Bento Cavalcante Neto
**   Introduced functions draw_interior_mesh and draw_mesh to draw the mesh
**   in interior of quadtree created by the templates, and final mesh after
**   boundary contraction without smooth.
**
** Modified: 28-Jul-94    Joaquim Bento Cavalcante Neto
**   Modified function div_error_tree to avoid use of some elements that the
**   norm error is greater the max_error * max_norm in the model. Introduced
**   also a new function new_error_tree with another criteria of refinement.
**   Tests will decide wich function is the best to be used.
**
** Modified: 12-Aug-94    Joaquim Bento Cavalcante Neto
**   Corrected function get_edge because when x or y coordinates are the same
**   the new coordinates should be corrected and the function should return,
**   and this was not done before.
**
** Modified: 01-Jun-96    Joaquim Bento Cavalcante Neto &&
**        Eduardo Setton Sampaio Silveira
**   Modified functions ini_sub_tree, create_tree, make_tree, ins_vertex_tree,
**   div_error_tree, classify_tree, traverse_tree and clip_face to include a
**   new parameter type_mesh to indicate if the mesh is a triangular mesh or a
**   quadrilateral one. The modification is when is a quadrilateral mesh, the
**   edges in the boundary should be taken considering two edges each time.
**
** Modified: 21-Jun-96    Joaquim Bento Cavalcante Neto &&
**   Modified function smooth_fem and new_position_fem to include a type_mesh
**   parameter to make a diferent smooth from triangular and quadrilateral
**   elements.
**
** Supervised:      Luiz Fernando Martha
**
** ---------------------------------------------------------------------------
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "msh2d.h"
#include "quadtree.h"
#include "quadsll.h"

#include "btree.h"			// mshaux library

#define QUADTREE_C

#define  MAX(a,b)  ( (a>b) ? a: b )
#define  MIN(a,b)  (((a)<(b))?(a):(b))
#define  NUMSTEP 4

/*
** ---------------------------------------------------------------------------
** Definitions for use of algorithms:
*/

#define QUAD_DEBUG 0
#define QUAD_ERROR 0
#define QUAD_BOUND 1

#define QUAD_REF_FACTOR  1.0
/*
** ---------------------------------------------------------------------------
** Includes for DEBUG:
*/

#if QUAD_DEBUG
#include "mquadraw.h"
#include "mquadgra.h"
static int flag_debug = 1;
FILE   *arq_quad;
#endif

/*
** ---------------------------------------------------------------------------
** Local functions prototypes:
*/

static void        convert_geom( int, int *, GeoPnt *, face_t ** , int);
static void        clear_geom( face_t **, mesh_vertex_t ** );
static void        ini_sub_tree( face_t *, int );
static box2D_t     find_box( face_t * );
static void        ini_tree( box2D_t * );
static void        create_tree ( face_t *, int );
static void        make_tree ( new_coord *, tree_node_t *, int, int );
static void        subdivide_tree( tree_node_t * );
static void        only_onelevel_tree( void );
static void        correct_level( tree_node_t * );
static void        ins_vertex_tree ( face_t *, int *, int );
static void        classify_tree ( face_t *, double, int );
static void        traverse_tree( face_t *, tree_node_t *, double, int );
static tree_code_t clip_face ( face_t *, tree_node_t *, double, int  );
static void        get_edge( int , int , int , new_coord *, new_coord * );
static void        codigo ( int , int , win2D_t *, int * );
static int         pick_face( int , int , face_t * );
static void        bound_size( tree_node_t * );
#if QUAD_BOUND
static void        only_boundsize_tree( void );
#endif
#if QUAD_ERROR
static void    div_error_tree( face_t *, int );
static void      new_error_tree( face_t * );
static void        error_size( tree_node_t * );
static int         error_adj( tree_node_t *, int * );
#endif
static void        smooth_fem( int, int, int, int, mesh_elem_t *, mesh_vertex_t * );
static void        adj_faces_fem( int, mesh_elem_t *, mesh_vertex_t * );
static void        new_position_fem( int, int, int, mesh_elem_t *, mesh_vertex_t * );
static void        convert_fem( mesh_vertex_t **, mesh_elem_t **, int, int,
                                GeoPnt *, GeoPnt **, int **, int *, int * , int );
static void        clear_tree( tree_node_t * );
static void        clear_struct_tree( tree_node_t * );
static neighbour_t *locate( new_coord *, tree_node_t * );
static neighbour_t *list_loc_nodes ( new_coord *, tree_node_t * );
static void        ins_loc_list ( tree_node_t * );
static neighbour_t *adj_nodes ( side_t, tree_node_t * );
static int         get_path( tree_node_t *, int * );
static void        path_node( unsigned, int, int, int *, int *, tree_node_t * );
static tree_node_t *get_node ( int, int * );
static neighbour_t *adjacent( side_t, int, int, int *, tree_node_t * );
static void        add_nodes ( side_t, tree_node_t * );
static int        fn_quad_sort (void *, void *);
static void       MshQuadratic ( GeoPnt **, int **, int *, int , Tbtree *, int, GeoPnt *);
static void       Quadratic ( GeoPnt **, int **, int *, int , int, int *, GeoPnt *, face_t *);
static void       fix_mesh_orientation (int nnode, double *coords, int n_elem, int *conn);


/*
** ---------------------------------------------------------------------------
** Local static variables:
*/

static tree_t       *tree;          /* head of tree structure */
static unsigned     max_depth;      /* maximum depth in the tree */
static int          adj_type;       /* type of adjacent in the tree */
static int          sizemaxbound;   /* maximum bound size tree node */
static neighbour_t  *head_neigh;    /* head of the adjacent tree node list */
static neighbour_t  *head_loc;      /* head of the locate tree node list */

/*
** ---------------------------------------------------------------------------
** Public Mesh Generation Function:
*/

/* ========================= msh_quadtree ================================ */

int  Msh2DQuadBound(
int     n_loops,       /* # of circuits (connect.bdry parts)          (in) */
int     loop_segs[],   /* # of segs (elem.side) on circuits           (in) */
double  *Bdry_pts,     /* coordinates of points on boundary           (in) */
int     type_mesh,     /* type for triangular or quadrilateral mesh   (in) */
int     ref_quad,      /* 1 = refine quadtree, 0 == not refine        (in) */
int     *n_node,       /* counts # of pts for meshing                (out) */
double  **Coords,      /* coordinate array used for meshing          (out) */
int     *n_elem,       /* number of elements generated               (out) */
int     **conn         /* elem.connectivity list from meshing        (out) */
)

{
 /* Local variables and symbols */

 face_t        *f = NULL;         /* pointer for local gometric structure */
 mesh_vertex_t *no_vector = NULL; /* finite element vertex structure      */
 mesh_elem_t   *el_vector = NULL; /* finite element structure             */
 mesh_vertex_t *vhead = NULL;     /* quadtree finite element vertexs number */
 int           nnogeom = 0;       /* geometry finite element vertexs number */
 int           nno = 0;           /* total number of finite element vertexs */
 int           nel = 0;           /* total number of finite elements */
 int           nelquadtree = 0;   /* number of quadtree finite elements */
 int           rid = 0;           /* returned indicated generated mesh */
 GeoPnt        *bdry_pts = (GeoPnt *)Bdry_pts;
 GeoPnt        *coords = NULL;

#if QUAD_DEBUG
arq_quad=fopen("debug_quad.gdb","r");
if(arq_quad!=NULL)
{
 fscanf(arq_quad,"%d",&flag_debug);
 fclose(arq_quad);
}
else
 flag_debug=0;
#endif

 /* Convert the input geometric structure to local geometric structure */

 convert_geom( n_loops, loop_segs, bdry_pts, &f , type_mesh);

 #if QUAD_DEBUG
 if( flag_debug )
 {
  draw_clear( 1 );
  draw_boundary( f );
  draw_contour( f );
  ItfUtlMsg( "Debug" );
 }
 #endif

 /* Inicializate the tree, based in the local geometric structure */

 ini_sub_tree( f, type_mesh );
 #if 0 /*QUAD_DEBUG */
 if( flag_debug )
 {
  draw_clear( 1 );
  draw_boundary( f );
  draw_contour( f );
  draw_tree( tree );
  ItfUtlMsg( "Debug" );
 }
 #endif

 /* Divide the tree, based in the error of the elements inside of each tree
    leaf */

 #if QUAD_ERROR
 div_error_tree( f, type_mesh );

 #if QUAD_DEBUG
 if( flag_debug )
 {
  draw_clear( 1 );
  draw_boundary( f );
  draw_contour( f );
  draw_tree( tree );
  ItfUtlMsg( "Debug" );
 }
 #endif
 #endif

 /* Insurance only one level difference between adjacent nodes into the tree */

 only_onelevel_tree( );

 #if 0 /*QUAD_DEBUG */
 if( flag_debug )
 {
  draw_clear( 1 );
  draw_boundary( f );
  draw_contour( f );
  draw_tree( tree );
  ItfUtlMsg( "Debug" );
 }
 #endif

 /* Update the tree nodes that have a vertex inside them */

 ins_vertex_tree( f, &nnogeom, type_mesh );

 /* Classify the tree nodes with a collapse percentage of 10% */

 classify_tree( f, 0.1, type_mesh );

 #if 0 /*QUAD_DEBUG */
 if( flag_debug )
 {
  draw_clear( 1 );
  draw_boundary( f );
  draw_contour( f );
  draw_tree( tree );
  draw_interior_tree( tree );
  ItfUtlMsg( "Debug" );
 }
 #endif

 /* Insurance that have no interior tree node bigger than the biggest boundary
   tree node */

 #if QUAD_BOUND
 if (ref_quad)
   only_boundsize_tree( );

 #if 0 /* QUAD_DEBUG */
 if( flag_debug )
 {
  draw_clear( 1 );
  draw_boundary( f );
  draw_contour( f );
  draw_tree( tree );
  draw_interior_tree( tree );
  ItfUtlMsg( "Debug" );
 }
 #endif
 #endif

 /* Make triangular or quadrilateral mesh depending of the client desire */

 switch( type_mesh )
 {
   case 3 :
   case 6 :
   /* Make triangular mesh in the interior tree nodes using templates */

   trimesh_int( tree->root, nnogeom, &nno, &nel, &nelquadtree, &vhead );

   #if QUAD_DEBUG
   if( 0 /*flag_debug*/ )
   {
    draw_clear( 1 );
    draw_boundary( f );
    draw_contour( f );
    draw_tree( tree );
    draw_interior_mesh( tree );
    ItfUtlMsg( "Debug" );
    draw_clear( 1 );
    draw_boundary( f );
    draw_contour( f );
    draw_interior_mesh( tree );
    ItfUtlMsg( "Debug" );
   }
   #endif

   /* Make triangular boundary contraction in the rest of the geometry */

   rid = trimesh_bound( tree->root, f, vhead, nno, nnogeom, &nel, &no_vector,
                        &el_vector );

   #if QUAD_DEBUG
   if( flag_debug && rid )
   {
    draw_clear( 1 );
    draw_boundary( f );
    draw_contour( f );
    draw_interior_mesh( tree );
    draw_mesh( nel, nelquadtree, flag_debug, no_vector, el_vector );
    ItfUtlMsg( "Debug" );
   }
   #endif

   /* Escape switch */

   break;

   case 4 :
   case 8 :

   /* Make quadrilateral mesh in the interior tree nodes using templates */

   quadmesh_int( tree->root, nnogeom, &nno, &nel, &nelquadtree, &vhead );

   #if QUAD_DEBUG
   if( 0 /* flag_debug */ )
   {
    draw_clear( 1 );
    draw_boundary( f );
    draw_contour( f );
    draw_tree( tree );
    draw_interior_mesh( tree );
    ItfUtlMsg( "Debug" );
    draw_clear( 1 );
    draw_boundary( f );
    draw_contour( f );
    draw_interior_mesh( tree );
    ItfUtlMsg( "Debug" );
   }
   #endif

   /* Make quadrilateral boundary contraction in the rest of the geometry */

   rid = quadmesh_bound( tree->root, f, vhead, &nno, nnogeom, &nel, &no_vector,
                         &el_vector );

   #if QUAD_DEBUG
   if( flag_debug && rid )
   {
    draw_clear( 1 );
    draw_boundary( f );
    draw_contour( f );
    draw_interior_mesh( tree );
    draw_mesh( nel, nelquadtree, flag_debug, no_vector, el_vector );
    ItfUtlMsg( "Debug" );
   }
   #endif

   /* Escape switch */

   break;
 }

 /* If it was possible to generate the mesh (rid == 1) then smooth it */


 if ( rid )
 {
  /* Smooth the mesh */

  smooth_fem( nel, nno, nnogeom, type_mesh, el_vector, no_vector );

  /* Convert local finite element structure to export finite element one */

  convert_fem( &no_vector, &el_vector, nno, nel, bdry_pts, &coords, conn,
               n_node, n_elem, type_mesh);

 }

 if (*n_node != 0 && n_elem != 0)
 {
   if(((type_mesh==6)||(type_mesh==8)) && rid)
     Quadratic ( &coords, conn, n_node, *n_elem, n_loops, loop_segs, bdry_pts, f);
   (*Coords) = (double *)coords;

   /* Fix orientation */
   fix_mesh_orientation (*n_node, *Coords, *n_elem, *conn);
 }
 else
   rid = 0;

 /* Release local finite element structure */
 free( no_vector );
 free( el_vector );

 /* Release local geometric structure and quadtree finite element vertexs */

 clear_geom( &f, &vhead );

 /* Release tree structure */

 clear_tree( tree->root );
 free( tree );

 /* Return indicated generated mesh */

 return rid;
}

/*
** ---------------------------------------------------------------------------
** Local quadtree functions:
*/

/* ========================= convert_geom ================================= */

static void convert_geom( int n_loops, int *loop_segs, GeoPnt *bdry_pts,
                          face_t **f ,int type_mesh)
{
 face_t *cnew;     /* new pointer for local geometric structure */
 int    ind = 0;   /* auxiliary vertex index */
 int    i;         /* auxiliary loop counter */
 int    j;         /* auxiliary loop counter */
 int    md = 1;

 /* Look trough the loops to convert to the local geometrical structure */

 if((type_mesh==6)||(type_mesh==8))
   md = 2;

 for( i = 0; i < n_loops; i++ )
  {
   /* Allocate a new element */

   SllAddEnd( (Sll *)f, sizeof( face_t ), (Sll *)&cnew );

   /* Update the new element vector */

   cnew->vt = ( vertex_t ** ) calloc( (loop_segs[i]/md), sizeof(vertex_t *) );

   /* Update the new element vector structure */

   for( j=0; j < (loop_segs[i]/md); j++ )
     {
      /* Allocate vertex and update its coordinates */

      cnew->vt[j] = ( vertex_t * ) calloc( 1, sizeof( vertex_t ) );
      cnew->vt[j]->x = bdry_pts[ind*md].x;
      cnew->vt[j]->y = bdry_pts[ind*md].y;
      cnew->vt[j]->bdrypt = ind;

      /* Update the auxiliary vertex index */

      ind++;
     }

   /* Update the number of nodes */

   cnew->n = loop_segs[i]/md;
  }
}

/* ==================== clear_geom ======================================== */

static void clear_geom( face_t **f, mesh_vertex_t **v )
{
 int i;

 /* Release local geometric structure */

 while( *f )
 {
  /* Release vector structure */

  for( i=0; i < (*f)->n; i++ )
    free( (*f)->vt[i] );

  free( (*f)->vt );

  /* Release element list */

  SllDelete( (Sll *)f, (Sll)*f );
 }

 /* Release quadtree finite element vertexs structure */

 SllDelAll( (Sll *)v );
}

/* ========================= ini_sub_tree ================================= */

static void ini_sub_tree( face_t *f, int type_mesh )
{
 box2D_t  box;             /* two-dimensional box to generate the tree */

 /* Find the tree bounding box */

 box = find_box( f );

 /* Inicializate the tree root */

 ini_tree( &box );

 /* Create the inicial tree */

 create_tree( f, type_mesh );
}

/* ========================== find_box ==================================== */

static box2D_t find_box( face_t *f )
{
 int     i;
 double   sizex, sizey, c;
 box2D_t  box;

 /* Inicializate the box structure */

 box.xmin = box.xmax = f->vt[0]->x;
 box.ymin = box.ymax = f->vt[0]->y;

 /* Update the box structure */

 for( i = 0; i < f->n; i++ )
 {
   if( f->vt[i]->x < box.xmin ) box.xmin = f->vt[i]->x;
   if( f->vt[i]->x > box.xmax ) box.xmax = f->vt[i]->x;
   if( f->vt[i]->y < box.ymin ) box.ymin = f->vt[i]->y;
   if( f->vt[i]->y > box.ymax ) box.ymax = f->vt[i]->y;
 }

 /* Return the final box structure */

 sizex = box.xmax - box.xmin;
 sizey = box.ymax - box.ymin;
 if( sizex > sizey )
 {
   c = ( box.ymax + box.ymin ) / 2.0;
   box.ymin = c - sizex / 2.0;
   box.ymax = c + sizex / 2.0;
 }
 else
 {
   c = ( box.xmax + box.xmin ) / 2.0;
   box.xmin = c - sizey / 2.0;
   box.xmax = c + sizey / 2.0;
 }
 return box;
}

/* ========================= ini_tree ===================================== */

static void ini_tree( box2D_t *lim )
{
 /* Inicializate the tree root */

 tree       = (tree_t *)      calloc( 1, sizeof(tree_t) );
 tree->root = (tree_node_t *) calloc ( 1, sizeof(tree_node_t) );

 tree->box              = *lim;
 tree->root->code       = NODE_UNDEF;
 tree->root->depth      = 0;
 tree->root->edepth     = 0;
 tree->root->coord.ex   = 0;
 tree->root->coord.ey   = 0;
 tree->root->vhead      = NULL;
 tree->root->elmhead    = NULL;
 tree->root->child      = NULL;

 /* Inicializate the maximum depth in the tree node */

 max_depth = 0;
}

/* ========================= create_tree ================================== */

static void  create_tree ( face_t *f, int type_mesh )
{
 int        i, ilength;
 double      xaux = 0, yaux = 0, length = 0;
 new_coord  p;
 double     reffactor = Msh2DGetRefFactor ( );

 /* Look trough all edges in the geometry to make the inicial tree */

 while( f )
 {
  i = 0;
  while( i < f->n ) /* for( i = 0; i < f->n; i++ ) */
  {
   /* Calculate the edge lenght */

   switch( type_mesh )
   {
    case 3 :
  case 6 :
     length = sqrt( ( f->vt[(i+1)%f->n]->x - f->vt[i]->x ) *
                    ( f->vt[(i+1)%f->n]->x - f->vt[i]->x ) +
                    ( f->vt[(i+1)%f->n]->y - f->vt[i]->y ) *
                    ( f->vt[(i+1)%f->n]->y - f->vt[i]->y ) );
    break;

    case 4 :
  case 8 :
     length = sqrt( ( f->vt[(i+2)%f->n]->x - f->vt[i]->x ) *
        ( f->vt[(i+2)%f->n]->x - f->vt[i]->x ) +
        ( f->vt[(i+2)%f->n]->y - f->vt[i]->y ) *
        ( f->vt[(i+2)%f->n]->y - f->vt[i]->y ) );
    break;
   }

   /* Translate the double value lenght to an integer value lenght */

   ilength = pow_2(14) * length / reffactor / ( tree->box.xmax - tree->box.xmin );

   /* Calculate the edge middle point */

   switch( type_mesh )
   {
    case 3 :
  case 6 :
     xaux = ( f->vt[i]->x + f->vt[((i+1)%f->n)]->x ) / 2. ;
     yaux = ( f->vt[i]->y + f->vt[((i+1)%f->n)]->y ) / 2. ;
    break;

    case 4 :
  case 8 :
     xaux = ( f->vt[i]->x + f->vt[((i+2)%f->n)]->x ) / 2. ;
     yaux = ( f->vt[i]->y + f->vt[((i+2)%f->n)]->y ) / 2. ;
    break;
   }

   /* Translate double edge middle point to an integer edge middle point */

   param_form( xaux, yaux, &p );

   /* Update the tree by the edge middle point */

   make_tree( &p, tree->root, ilength, type_mesh );

   /* Get next edge */

   switch( type_mesh )
   {
    case 3 :
  case 6 :
     i++;
    break;

    case 4 :
  case 8 :
     i += 2;
    break;
   }
  }

  /* Get the next hole geometry structure, if exists */

  f = f->next;
 }
}

/* ========================= make_tree ==================================== */

static void make_tree ( new_coord *p, tree_node_t *node,
      int length, int type_mesh )
{
 neighbour_t  *n;
 int          size;
 double        factor = .0f;

 /* Locate in the actual state of the tree, in what node the edge middle point
    is.  The function locate_quad return a list of tree nodes  where the point
    could be,that seldom is just one node but eventually could be more */

 n = locate_quad( p, node );

 /* Look trough the returned list.If the size node where the edge middle point
    is it's bigger than a edge lenght percentage( here 100% ),the node  should
    be divided and the process is recursively done until be less than */

 while ( n )
 {
  /* If the first tree node size in the n list is less than 100% edge lenght,
     update max depth, release n list to this call of make_tree and return.
     Note that is used a tolerance of 1 in the test of the size */

  size = pow_2( 14 - n->node->depth );

  switch( type_mesh )
  {
   case 3 :
   case 6 :
    factor = 1.2;
   break;

   case 4 :
   case 8 :
    factor = 1.2;
   break;
  }

  if( size <= (factor*length + 1) )
  {
   /* Update the maximum depth in the tree */

   max_depth = MAX(n->node->depth, max_depth);

   /* Release all n locate_quad list, even it has not be looked at the others
      tree nodes in the n list, because  like  a  tree node is less than 100%
      edge lenght,  it doesn't need to  look  in others tree nodes in  the  n
      list, because they are childs too and will be less than too, correspon-
      ding that they won't be subdivided */

   SllDelAll( (Sll *)&n );

   /* Return to the other call of make_tree */

   return;
  }

  /* Else divide the tree node and and the process is recursively done until
     the tree node size be less than edge length*/

  /* Subdivide tree node */

  subdivide_tree( n->node );

  /* Repeat the process */

  make_tree( p, n->node, length, type_mesh );

  /* Release the n locate_quad list used and get a new one */

  SllDelete( (Sll *)&n, (Sll)n );
 }
}

/* ========================= subdivide_tree =============================== */

static void subdivide_tree( tree_node_t *node )
{
 int i, j;
 int ex, ey;
 int n;
 int sizlev;

 /* Allocate memory for the four child tree node */

 node->child = ( tree_node_t ** ) calloc( 4, sizeof(tree_node_t *) );

 /* Calculate the tree node size */

 sizlev = pow_2( 13 - node->depth );

 /* Update the four child tree node */

 n = 0;
 for( i = 0, ey = node->coord.ey; i < 2; ey += sizlev, i++ )
   for( j = 0, ex = node->coord.ex; j < 2; ex += sizlev, j++ )
   {
     node->child[n] = ( tree_node_t * ) calloc( 1, sizeof(tree_node_t) );
     node->child[n]->code       = NODE_UNDEF;
     node->child[n]->depth      = node->depth  + 1;
     node->child[n]->edepth     = node->edepth + 1;
     node->child[n]->coord.ex   = ex;
     node->child[n]->coord.ey   = ey;
     node->child[n]->vhead      = NULL;
     node->child[n]->elmhead    = NULL;
     node->child[n]->child      = NULL;
     n++;
   }
}

/* ===================== div_error_tree =================================== */

#if QUAD_ERROR
static void div_error_tree( face_t *f, int type_mesh )
{
 int     i,j,l,k,bdry,np,narea;
 int     n;
 int       isize;
 int       esize;
 double    size;
 double           hmax;
 double    helm;
 double   *x, xp;
 double   *y, yp;
 double           ratio_error;
 new_coord   p;
 GeoPnt         *area;
 GeoPnt          pt[8];
 double    norm_error, norm_max, norm_disp, hold, a;

 /* Find the maximum norm error */

 #if 0
 norm_max = 0.0;
 for( i = 0; i < nelem; i++ )
 {
  ElmNormError( &elem_vector[i], &norm_error, &norm_disp );
  norm_max = MAX(norm_max,norm_error);
 }
 #endif

 /* Get the region area */

 narea = f->n;
 area = (GeoPnt *) calloc( narea, sizeof( GeoPnt ) );
 for( i = 0; i < narea; i++ )
 {
  area[i].x = f->vt[i]->x;
  area[i].y = f->vt[i]->y;
 }

 /* Update the tree based in the error of each element */

 for( i = 0; i < nelem; i++ )
 {
  /* Get the element middle point coordinate */

  np = 0;
  xp = yp = 0.0;
  for( l = 0; l < elem_vector[i].n; l++ )
  {
   xp += node_vector[elem_vector[i].inc[l]].coord.x;
   yp += node_vector[elem_vector[i].inc[l]].coord.y;
   np++;
  }
  xp /= np;
  yp /= np;

  /* Verify if the element is inside the area */

  if( PkArea( narea, area, xp, yp ) )
  {
   /* Verify if the norm is greater error */

   #if 0
   ElmNormError( &elem_vector[i], &norm_error, &norm_disp );
   #endif

   if( 1 )  /* option => ( norm_error >= max_error * norm_max ) */
   {
    /* Get the element error size */

    ElmSizeError( &elem_vector[i], &size );

    /* Translate the error size length to integer value length */

    isize = pow_2(14) * size / ( tree->box.xmax - tree->box.xmin );

    /* Get the element edges middles points */

    ElmEdgeError( &elem_vector[i], &n, &x, &y );

    /* Update the tree based in the element error size */

    for( j = 0; j < n; j++ )
    {
     /* Translate the edge middle point to integer middle point */

     param_form( x[j], y[j], &p );

     /* Update the tree by the edge middle point. The size should be
  duplicated in quadrilateral mesh as is done in the boundary
  when the initial mesh is being created */

     switch( type_mesh )
     {
      case 3 :
    case 6 :
       esize = isize;
      break;

      case 4 :
    case 8 :
       esize = 2 * isize;
      break;
     }

     make_tree( &p, tree->root, esize, type_mesh );
    }
   }
  }
 }
}
#endif

/* ====================== only_onelevel_tree ============================== */

static void only_onelevel_tree( void )
{
 /* Insurance only one level difference between adjacent nodes in the tree */

 correct_level( tree->root );
}

/* ====================== correct_level =================================== */

static void correct_level( tree_node_t *node )
{
 int i;
 int j;

 /* The routine traverse recursivilying the tree and if the tree node has  no
    child, that is, if the tree node is a leaf, insure that the tree adjacent
    nodes has only one difference level */

 /* Verify if the tree node is a leaf */

 if( node->child != NULL )
 {
  for( i = 0; i < 4; i++ )
   correct_level( node->child[i] );
 }

 /* Insurance only one level difference */

 else
 {
  int          nr = 1, nl = 1, nb = 1, nt = 1;
  int          type[4] = { 0, 0, 0, 0 };
  neighbour_t  *h[4], *n;

  /* Get number of adjacent tree nodes of the current tree node from each side,
     that is, nr - RIGHT , nl - LEFT, nb - BOTTOM, nt - TOP */

  h[0] = adj_nodes ( RIGHT, node );
  if ( adj_type ) type[0] = 1;
  for ( nr = 0, n = h[0]; n; nr++, n = n->next );
  h[1] = adj_nodes ( LEFT, node );
  if ( adj_type ) type[1] = 1;
  for ( nl = 0, n = h[1]; n; nl++, n = n->next );
  h[2] = adj_nodes ( BOTTOM, node );
  if ( adj_type ) type[2] = 1;
  for ( nb = 0, n = h[2]; n; nb++, n = n->next );
  h[3] = adj_nodes ( TOP, node );
  if ( adj_type ) type[3] = 1;
  for ( nt = 0, n = h[3]; n; nt++, n = n->next );

  /* Insurance only one level difference between adjacent tree nodes, that is,
     number of adjacent tree nodes could not be more than two tree nodes */

  if( nr > 2 || nl > 2 || nb > 2 || nt > 2 )
  {
    /* If any adjacent has more than two tree nodes, it means that no have only
       one level difference, and so the current tree node should be divided */

    subdivide_tree( node );
    node->code = NODE_UNDEF;

    /* Ensure only one level for the tree nodes on the adjacent list too */

    for( i= 0; i < 4; i++ )
    {
     if( type[i] == 1 ) correct_level( h[i]->node );
    }

    /* Ensure only one level for the tree nodes childs too */

    for( i = 0; i < 4; i++ )
    {
      correct_level( node->child[i] );
    }
  }

  /* Release tree nodes adjacent list */

  for( j=0; j < 4; j++ )
    SllDelAll( (Sll *)&h[j] );
 }
}

/* ====================== ins_vertex_tree ================================= */

static void ins_vertex_tree ( face_t *f , int *nnogeom, int type_mesh )
{
 int           i;
 int           nnoaux = 0;
 int           sizeb;
 new_coord     p;
 neighbour_t   *loc;
 neighbour_t   *loc_head;

 /* Inicializate bigger bound tree node size */

 sizemaxbound = 0;

 /* Look trough the geometry to classify the tree nodes with a vertex inside */

 while ( f )
 {
  i = 0;
  while( i < f->n )
  {
   /* Get an geometrical vertex from the geometric structure */

   param_form( f->vt[i]->x, f->vt[i]->y, &p );

   /* Locate the tree nodes list where the vertex is */

   loc = loc_head = locate_quad( &p, tree->root );

   /* Classify the tree nodes like NODE_VERTEX and update the tree node data,
      already geting the bigger tree node bound size that will be used in the
      only_boundsize_tree function */

   while (loc)
   {
    /* Insert new vertex into the tree node vertex head list */

    if( SllAddTop((Sll *)&loc->node->vhead, sizeof( tree_vertex_t )) )
     {
      /* Update the tree node data and tree node vertex head list */

      loc->node->code = NODE_VERTEX;
      loc->node->vhead->id = nnoaux + i;
      loc->node->vhead->coord = p;
     }

    /* Update the bigger bound tree node size */

    sizeb = pow_2( 14 - loc->node->depth );
    if( sizeb > sizemaxbound ) sizemaxbound = sizeb;

    /* Get a new tree node and release the tree node used */

    loc = loc->next;
   }

   /* Release loc list */

   SllDelAll( (Sll *)&loc_head );

   /* Get next edge */

   switch( type_mesh )
   {
    case 3 :
  case 6 :
     i++;
    break;

    case 4 :
  case 8 :
     i += 2;
    break;
   }
  }

  /* Update the number of finite element vertexs in the geometry */

  nnoaux += f->n;

  /* Get the hole geometry structure, if exists */

  f = f->next;
 }

 /*Update final nnogeom*/

 *nnogeom = nnoaux;
}

/* ==================== classify_tree ===================================== */

static void classify_tree ( face_t *f, double codtol, int type_mesh )
{
 /* Traverse the tree and classify the tree nodes that even was  not   classi-
    fied( NODE_UNDEF ) in the first call, or collapse the tree node classified
    in the interior of geometry( NODE_INTERIOR ) in the second call */

 traverse_tree ( f, tree->root, codtol, type_mesh );
}

/* ==================== traverse_tree ===================================== */

static void traverse_tree( face_t *f, tree_node_t *node, double codtol,
         int type_mesh )
{
 int i;

 /* This routine is recursively done until all tree nodes be classified */

 /* Verify if the tree node has no child,that is, if it is a leaf */

 if( node->child != NULL )
  for( i = 0; i < 4; i++ )
   traverse_tree( f, node->child[i], codtol, type_mesh );

 /* Then classify the tree nodes that even has no classification( NODE_UNDEF )
    or the tree nodes that are NODE_INTERIOR( this is because that in the  se-
    cond time when this function is used, it's to collapse interior tree nodes
    nearly  the geometry ) */

 else if( node->code == NODE_UNDEF )
       node->code = clip_face ( f, node, codtol, type_mesh );
}

/* ==================== clip_face ========================================= */

static tree_code_t clip_face ( face_t *f, tree_node_t *node, double codtol,
             int type_mesh )
{
 int        i, j, code;
 int        cod0[4], cod1[4];
 int        trivial_key;
 double      length = .0f;
 int        ilength;
 int        size;
 int        sizeb;
 int        tolx;
 int        toly;
 int        k;
 new_coord  pa,pb;
 new_coord  p1 = {0, 0}, p2 = {0, 0};
 win2D_t    wn;
 face_t     *faux;

 /* Calculate the tree node size and their corners */

 size = pow_2( 14 -  node->depth );
 wn.xmin = node->coord.ex;
 wn.xmax = node->coord.ex + size;
 wn.ymin = node->coord.ey;
 wn.ymax = node->coord.ey + size;

 /* Test all geometrical edges against the tree node.If the edge cross the no-
    de, it's NODE_BOUND. Else it will be NODE_INTERIOR or NODE_OUT   depending
    of their position */

 /* Get the inicial geometry pointer */

 faux = f;

 while ( f )
 {
  i = 0;
  while( i < f->n )
  {
   /* Calculate the edge lenght */

   switch( type_mesh )
   {
    case 3 :
  case 6 :
     length = sqrt( ( f->vt[(i+1)%f->n]->x - f->vt[i]->x ) *
                    ( f->vt[(i+1)%f->n]->x - f->vt[i]->x ) +
                    ( f->vt[(i+1)%f->n]->y - f->vt[i]->y ) *
                    ( f->vt[(i+1)%f->n]->y - f->vt[i]->y ) );
    break;

    case 4 :
  case 8 :
     length = sqrt( ( f->vt[(i+2)%f->n]->x - f->vt[i]->x ) *
        ( f->vt[(i+2)%f->n]->x - f->vt[i]->x ) +
        ( f->vt[(i+2)%f->n]->y - f->vt[i]->y ) *
        ( f->vt[(i+2)%f->n]->y - f->vt[i]->y ) );
    break;
   }

   /* Translate the double value lenght to an integer value lenght */

   ilength = pow_2(14) * length / ( tree->box.xmax - tree->box.xmin );

   /* Calculate the integer value for the two edge vertexs */

   switch( type_mesh )
   {
    case 3 :
  case 6 :
     param_form( f->vt[i]->x, f->vt[i]->y, &p1 );
     param_form( f->vt[((i+1)%f->n)]->x, f->vt[((i+1)%f->n)]->y, &p2 );
    break;

    case 4 :
  case 8 :
     param_form( f->vt[i]->x, f->vt[i]->y, &p1 );
     param_form( f->vt[((i+2)%f->n)]->x, f->vt[((i+2)%f->n)]->y, &p2 );
    break;
   }

   /* Calculate the two edge projections( in x and y ) */

   tolx  = codtol * ilength;
   toly  = codtol * ilength;

   /* Test the tree node against the edge.In first time that this function  is
      used, then it tests against the real edge,and the second time for colla-
      pse the tree nodes nearly the geometry, the test is done against two new
      edges  with  the width  enlarged  20%  for  the  two sides to avoid this
      proximity */

   for( k = -1; k < 2; k++ )
   {
    /* Get the new two edges.This  step  has  no  meaning  in  the  first  call
       to the function( when codtol = 0.0 ) because the tests will be done aga-
       inst the real edge and not against the new edges */

    pa.ex = p1.ex;
    pa.ey = p1.ey;
    pb.ex = p2.ex;
    pb.ey = p2.ey;
    get_edge( k, tolx, toly, &pa, &pb );

    /* Get the second edge vertex code */

    trivial_key = 0;
    codigo ( pb.ex, pb.ey, &wn, cod1 );

    /* Stay in this infinite loop until the edge is out of the tree   node  or
       the edge cross the tree node( NODE_BOUND ) */

    while ( 1 )
     {
      codigo ( pa.ex, pa.ey, &wn, cod0 );

      for ( j = 0; j < 4; j++ )
       if (cod0[j] && cod1[j]) /* test no-pick trivial */
        {
         trivial_key = 1;  /* edge tottaly left, right,top or down */
         break;
        }
       if ( trivial_key ) break;

      /* move point 0 for a window limit and stop if the edge cross the node */

      if ( cod0[0] )
       {
        pa.ey += ( wn.xmin - pa.ex ) * ( pb.ey - pa.ey ) / ( pb.ex - pa.ex );
        pa.ex = wn.xmin;
       }
      else if ( cod0[1] )
            {
             pa.ey += ( wn.xmax - pa.ex ) * ( pb.ey - pa.ey ) /
                      ( pb.ex - pa.ex );
             pa.ex = wn.xmax;
            }
           else if ( cod0[2] )
                 {
                  pa.ex += ( wn.ymin - pa.ey ) * ( pb.ex - pa.ex ) /
                           ( pb.ey - pa.ey );
                  pa.ey = wn.ymin;
                 }
                else if ( cod0[3] )
                      {
                       pa.ex += ( wn.ymax - pa.ey ) * ( pb.ex - pa.ex ) /
                                ( pb.ey - pa.ey );
                       pa.ey = wn.ymax;
                      }
                     else
                      {
                       /*Update maximum tree node size bound if the edge is
                         not an *alargada* edge*/

                       if( k == 0 )
                        {
                         sizeb = pow_2( 14 - node->depth );
                         if( sizeb > sizemaxbound ) sizemaxbound = sizeb;
                        }

                       /*Return tree node boundary*/

                       return NODE_BOUND;
                      }
     }
   }

   /* Get next edge */

   switch( type_mesh )
   {
    case 3 :
  case 6 :
     i++;
    break;

    case 4 :
  case 8 :
     i += 2;
    break;
   }
  }

  /* Get the hole geometry pointer to test against the hole edges */

  f = f->next;
 }

 /* Test if box contains face */

 param_form( faux->vt[0]->x, faux->vt[0]->y, &p1 );
 if ( (p1.ex > wn.xmin) && (p1.ex < wn.xmax) &&
      (p1.ey > wn.ymin) && (p1.ey < wn.ymax) )
  {
   code = NODE_BOUND;

   /*Update maximum tree node bound*/

   sizeb = pow_2( 14 - node->depth );
   if( sizeb > sizemaxbound ) sizemaxbound = sizeb;
  }

 /* Test if face contains box */

 if ( pick_face ( wn.xmin+size/2.0, wn.ymin+size/2.0,  faux ) )
  code = NODE_INTERIOR;
 else
  code = NODE_OUT;

 /* Return the tree node code */

 return code;
}

/* ==================== get_edge ========================================= */

static void get_edge( int k, int tolx, int toly, new_coord *pa, new_coord *pb )
{
 /* Get the two *alargadas* edges that will be tested against the tree node in
    the second call of classify_tree to avoid tree nodes nearly to the  geome-
    try */

 /* Get the large edges for a same x coordinate edge */

 if( pa->ex == pb->ex )
  {
   pa->ex += ( k * toly );
   pb->ex += ( k * toly );
   return;
  }

 /* Get the large edges for a same y coordinate edge */

 if( pa->ey == pb->ey )
  {
   pa->ey += ( k * tolx );
   pb->ey += ( k * tolx );
   return;
  }

 /* Get the large edges for an no-particular edge */

 if( ( (pa->ex > pb->ex) && (pa->ey < pb->ey) )  ||
     ( (pa->ex < pb->ex) && (pa->ey > pb->ey) )   )
  {
   pa->ex += ( k * tolx );
   pa->ey += ( k * toly );
   pb->ex += ( k * tolx );
   pb->ey += ( k * toly );
   return;
  }
 else
  {
   pa->ex += ( -1.0 * k * tolx );
   pa->ey += ( k * toly );
   pb->ex += ( -1.0 * k * tolx );
   pb->ey += ( k * toly );
   return;
  }
}

/* ==================== codigo =========================================== */

static void codigo ( int x, int y, win2D_t *wn, int cod[4] )
{
 /* Get the code box for the tree node */

 cod[0] = ( x < wn->xmin );
 cod[1] = ( x > wn->xmax );
 cod[2] = ( y < wn->ymin );
 cod[3] = ( y > wn->ymax );
}

/* ==================== pick_face ======================================== */

static int pick_face( int x, int y, face_t *f )
{
 int       i;
 int        ni = 0;     /* intersections number */
 double      xc;
 new_coord  p1, p2;

 /* Pick a face to verify intesection( NODE_INTERIOR verification ) */

 /* Get the inicial geometry pointer */

 while( f )
 {
  for( i = 0; i < f->n; i++ )
  {
   param_form( f->vt[i]->x, f->vt[i]->y, &p1 );
   param_form( f->vt[((i+1)%f->n)]->x, f->vt[((i+1)%f->n)]->y, &p2 );

   if( !(  p1.ey == p2.ey ) &&             /* discard horizontals   */
       !( (p1.ey > y) && (p2.ey > y) ) &   /* discard edges above   */
       !( (p1.ey < y) && (p2.ey < y) ) &&  /* discard edges bellow  */
       !( (p1.ex < x) && (p2.ex < x) ) )   /* discard edges left    */
    {
     if( p1.ey == y )        /* first point in the same cote  */
     {
      if( ( p1.ex > x ) && ( p2.ey > y ) )
       ni++;                 /* on right and above of point   */
     }
     else
     {
      if( p2.ey == y )               /* second point of the same cote */
       {
  if( ( p2.ex > x ) && ( p1.ey > y ) )
   ni++;                     /* on right and above of point   */
       }
      else
       {
        if( ( p1.ex > x ) && ( p2.ex > x ) )
         ni++;                      /* tottaly on right */
  else                            /* verify intersection point */
     {
    double dx = p1.ex - p2.ex;
          xc = p1.ex;
          if( dx != 0.0 )
          xc += ( y - p1.ey ) * dx / ( p1.ey - p2.ey );
    if( xc > x )
     ni++;
         }
       }
     }
    }
  }

  /* Get the hole geometry structure, if exists */

  f = f->next;
 }

 /* Return pick face code */

 return( ni % 2 );
}

/* ==================== only_boundsize_tree =============================== */

#if QUAD_BOUND
static void only_boundsize_tree( void )
{
 /* Insurance that have no interior tree node bigger than the biggest boundary
   tree node */

 bound_size( tree->root );
}
#endif

/* ==================== bound_size ======================================== */

static void bound_size( tree_node_t *node )
{
 int          i;
 int          n;
 int          node_size;
 tree_code_t  code;

 /* The routine traverse recursivilying the tree and if the tree node has  no
    child, that is, if the tree node is a leaf, insure that has  no  interior
    tree node bigger than the biggest boundary tree node */

 if( node->child != NULL )
  for( i = 0; i < 4; i++ )
   bound_size( node->child[i] );

 else if( node->code == NODE_INTERIOR )
       {
        node_size = pow_2( 14 - node->depth );
        if( node_size > sizemaxbound )
         {
          subdivide_tree( node );
          code = node->code;
          node->code = NODE_UNDEF;
          for( n = 0; n < 4; n++ )
            {
             node->child[n]->code = code;
             bound_size( node->child[n] );
            }
         }

       }
}

/* ========================== new_error_tree ============================== */

#if QUAD_ERROR
static void new_error_tree( face_t *f )
{
 int     i,j,l,k,bdry,np,narea;
 int     n;
 int       isize;
 int             node_size;
 double    size;
 double           hmax;
 double    helm;
 double   *x, xp;
 double   *y, yp;
 new_coord   p;
 GeoPnt         *area;
 GeoPnt          pt[8];
 double    norm_error, norm_max, norm_disp, hold, a, det;
 neighbour_t     *nodes;

 /* Get the region area */

 narea = f->n;
 area = (GeoPnt *) calloc( narea, sizeof( GeoPnt ) );
 for( i = 0; i < narea; i++ )
 {
  area[i].x = f->vt[i]->x;
  area[i].y = f->vt[i]->y;
 }

 /* Update the tree based in the error of each element */

 for( i = 0; i < nelem; i++ )
 {
  /* Get the element middle point coordinate */

  np = 0;
  xp = yp = 0.0;
  for( l = 0; l < elem_vector[i].n; l++ )
  {
   xp += node_vector[elem_vector[i].inc[l]].coord.x;
   yp += node_vector[elem_vector[i].inc[l]].coord.y;
   np++;
  }
  xp /= np;
  yp /= np;

  /* Verify if the element is inside the area */

  if( PkArea( narea, area, xp, yp ) )
  {
    /* Get the element error size */

    ElmSizeError( &elem_vector[i], &size );

    /* Translate the error size length to integer value length */

    isize = pow_2(14) * size / ( tree->box.xmax - tree->box.xmin );

    /* Get the element edges middles points */

    ElmEdgeError( &elem_vector[i], &n, &x, &y );

    for( j = 0; j < n; j++ )
    {
     /* Translate the edge middle point to integer middle point */

     param_form( x[j], y[j], &p );

     /* Update the tree by the edge middle point */

     nodes = locate_quad( &p, tree->root );

     while( nodes )
     {
      node_size = pow_2( 14 - nodes->node->depth );
      if( node_size <= isize )
       nodes->node->edepth = nodes->node->depth;
      else
      {
       det = (node_size/isize)/2.0 + 0.5;
       nodes->node->edepth = nodes->node->depth +  (int)det;
      }
      nodes = nodes->next;
     }
    }
   }
  }

 /* Traverse tree to update the level */

 error_size( tree->root );
}

/* ====================== error_size ====================================== */

static void error_size( tree_node_t *node )
{
 int i;
 int n;
 int k;
 int adepth;
 int edepth;
 int equal;
 tree_code_t code;

 /* The routine traverse recursivilying the tree and if the tree node has  no
    child, that is, if the tree node is a leaf, insure that has  no  interior
    tree node bigger than the biggest boundary tree node */

 if( node->child != NULL )
  for( i = 0; i < 4; i++ )
   error_size( node->child[i] );

 else if( node->edepth != 0 )
 {
  adepth = error_adj( node, &equal );
  edepth = MIN( node->edepth, (adepth + 1) ); /* rule 2 */
  if( node->depth >= edepth )
  {
   max_depth = MAX(node->depth, max_depth);
   return;
  }
  else
  {
   subdivide_tree( node );
   code = node->code;
   node->code = NODE_UNDEF;
   for( k = 0; k < 4; k++ )
    node->child[k]->edepth = edepth;
   for( n = 0; n < 4; n++ )
   {
    node->child[n]->code = code;
    error_size( node->child[n] );
   }
  }
 }
}

/* ========================== error_adj =================================== */

static int error_adj( tree_node_t *node, int *equal )
{
 int          depth;
 int        code;
 neighbour_t *adj;

 depth = 0;
 code  = 1;

 adj = find_adj_node( RIGHT, node );
 while( adj )
 {
  if( code )
   if( adj->node->edepth != node->edepth )
    code = 0;
  depth = MAX( adj->node->edepth, depth );
  SllDelete( (Sll *)&adj, (Sll)adj );
 }

 adj = find_adj_node( LEFT, node );
 while( adj )
 {
  if( code )
   if( adj->node->edepth != node->edepth )
    code = 0;
  depth = MAX( adj->node->edepth, depth );
  SllDelete( (Sll *)&adj, (Sll)adj );
 }

 adj = find_adj_node( TOP, node );
 while( adj )
 {
  if( code )
   if( adj->node->edepth != node->edepth )
    code = 0;
  depth = MAX( adj->node->edepth, depth );
  SllDelete( (Sll *)&adj, (Sll)adj );
 }

 adj = find_adj_node( BOTTOM, node );
 while( adj )
 {
  if( code )
   if( adj->node->edepth != node->edepth )
    code = 0;
  depth = MAX( adj->node->edepth, depth );
  SllDelete( (Sll *)&adj, (Sll)adj );
 }

 adj = find_adj_node( CORNER0, node );
 while( adj )
 {
  if( code )
   if( adj->node->edepth != node->edepth )
    code = 0;
  depth = MAX( adj->node->edepth, depth );
  SllDelete( (Sll *)&adj, (Sll)adj );
 }

 adj = find_adj_node( CORNER1, node );
 while( adj )
 {
  if( code )
   if( adj->node->edepth != node->edepth )
    code = 0;
  depth = MAX( adj->node->edepth, depth );
  SllDelete( (Sll *)&adj, (Sll)adj );
 }

 adj = find_adj_node( CORNER2, node );
 while( adj )
 {
  if( code )
   if( adj->node->edepth != node->edepth )
    code = 0;
  depth = MAX( adj->node->edepth, depth );
  SllDelete( (Sll *)&adj, (Sll)adj );
 }

 adj = find_adj_node( CORNER3, node );
 while( adj )
 {
  if( code )
   if( adj->node->edepth != node->edepth )
    code = 0;
  depth = MAX( adj->node->edepth, depth );
  SllDelete( (Sll *)&adj, (Sll)adj );
 }

 *equal = code;
 return depth;
}
#endif

/* ========================= smooth_fem =================================== */

static void smooth_fem(

int           nel,         /* total finite elements number             (in) */
int           nno,         /* total finite element vertexs number      (in) */
int           nnogeom,     /* geometric finite element vertexs number  (in) */
int           type_mesh,   /* type for triangular or quadrilateral mesh(in) */
mesh_elem_t   *el_vector,  /* finite elements structure                (in) */
mesh_vertex_t *no_vector   /* finite element vertexs structure     (in/out) */
)

{
  /* Make the adjacent finite element structure of the mesh */

  adj_faces_fem( nel, el_vector, no_vector );

  /* Suavize the mesh based in a lagrangean interpolation */

  new_position_fem( nno, nnogeom, type_mesh, el_vector, no_vector );
}

/* ========================= adj_faces_fem ================================ */

static void adj_faces_fem(

int            nel,        /* total finite elements number             (in) */
mesh_elem_t   *el_vector,  /* finite elements structure                (in) */
mesh_vertex_t *no_vector   /* finite element vertexs structure     (in/out) */
)

{
 int    i, id, j;
 adj_t  *fadj;

 /* Make a loop trough all finite elements to make the adjacent structure */
 for( i = 0; i < nel; i++ )
 {
  /* Get the finite element id */

  id = el_vector[i].nel;

  /* For each finite element, update the adjacent structure of each finite
     element vertex of it */

  for( j = 0; j < el_vector[i].n; j++ )
  {
    fadj = ( adj_t * ) calloc( 1, sizeof( adj_t ) );
    fadj->elm = id;
    fadj->next = no_vector[el_vector[i].conect[j]->nno].fadj;
    no_vector[el_vector[i].conect[j]->nno].fadj = fadj;
  }
 }
}

/* ======================= new_position_fem =============================== */

static void new_position_fem(

int           nno,         /* total finite element vertexs number      (in) */
int           nnogeom,     /* geometric finite element vertexs number  (in) */
int           type_mesh,   /* type for triangular or quadrilateral mesh(in) */
mesh_elem_t   *el_vector,  /* finite elements structure                (in) */
mesh_vertex_t *no_vector   /* finite element vertexs structure     (in/out) */
)

{
 int     i, j, k, n, id;
 int     pos = 0, before, after;
 point   cnew;
 adj_t   *fadj;

 /* Do the suavization about five times */

 for( i = 0; i < NUMSTEP; i++ )
 {
  /* Smooth the internal finite elements vertexs( quadtree vertexs ),
     because we can't change the geometrical finite element vertexs */

  for( j = 0; j < nno; j++ )
  {
    /* Smooth only the internal nodes */

    if( no_vector[j].bdrypt == -1 )
    {
     /* Inicializate the number of finite element adjacent vertexs of each
        vertex and the new coordinate of this vertex */

     n = 0;
     cnew.x = cnew.y = 0.0;

     /* Get the adjacent finite elements list of this vertex */

     fadj = no_vector[j].fadj;

     /* Look trough this list to make the smooth */

     while( fadj )
     {
      switch( type_mesh )
      {
       case 3 :
       case 6 :
       case 4 :
       case 8 :

  /* Get the contribution to vertex coordinates */

        for( k = 0; k < el_vector[fadj->elm].n; k++ )
        {
         /* Get the finite element vertex id */

         id = el_vector[fadj->elm].conect[k]->nno;

         /* Update the new position if the vertex it's not itself */

         if( id != no_vector[j].nno )
         {
           cnew.x += no_vector[id].x;
           cnew.y += no_vector[id].y;
           n++;
         }
        }

       break;

       case 5 :  /* case 4 : */

  /* Get the position of vertex and the adjacent vertexs of it */

        for( k = 0; k < el_vector[fadj->elm].n; k++ )
         if( el_vector[fadj->elm].conect[k]->nno == no_vector[j].nno )
            pos = k;
        before = pos - 1;
        after  = pos + 1;
        if( pos == 0 )
           before = el_vector[fadj->elm].n - 1;
        if( pos == (el_vector[fadj->elm].n - 1) )
          after  = 0;

  /* Get the contribution to vertex coordinates */

        for( k = 0; k < el_vector[fadj->elm].n; k++ )
        {
         /* Get the finite element vertex id */

         id = el_vector[fadj->elm].conect[k]->nno;

         /* Update the new position if the vertex it's not itself */

         if( id != no_vector[j].nno )
         {
          if( k == before || k == after )
          {
           cnew.x += no_vector[id].x;
           cnew.y += no_vector[id].y;
           n++;
          }
          else
          {
           cnew.x += 2 * no_vector[id].x;
           cnew.y += 2 * no_vector[id].y;
           n += 2;
          }
         }
        }

       break;
      }

      /* Get the next adjacent finite element of this vertex  in the list */

      fadj = fadj->next;
     }

     /* Update new vertex coordinates dividing by the number of adjacent ver-
        texs of this vertex */

     no_vector[j].x = cnew.x / n;
     no_vector[j].y = cnew.y / n;
    }
  }
 }
}

/* ========================= convert_fem ================================== */

static void convert_fem(

mesh_vertex_t **no_vector,  /* finite element vertexs structure        (in) */
mesh_elem_t   **el_vector,  /* finite elements structure               (in) */
int           nno,          /* total finite element vertexs number     (in) */
int           nel,          /* total finite elements number            (in) */
GeoPnt        *bdry_pts,    /* input boundary points vector            (in) */
GeoPnt        **coords,     /* coordinate array used for meshing      (out) */
int           **conn,     /* elem.connectivity list from meshing    (out) */
int           *n_node,      /* counts # of pts for meshing          (out) */
int           *n_elem,      /* number of elements generated         (out) */
int           type_mesh
)

{
 int       i;         /* loop counter */
 int       j;         /* loop counter */
 new_coord vert;      /* integer no_vector coordinates */
 vertex_t  verta;     /* double no_vector coordinates */
 int       size = 0;  /* elem.connectivity list from meshing size */
 int       index = 0; /* elem.connectivity list from meshing index */
 int       md=1;

 if((type_mesh==6)||(type_mesh==8))
   md=2;

 /* Update counts # of pts for meshing */

 *n_node = nno;

 /* Update number of elements generated */

 *n_elem = nel;

 /* Update coordinate array used for meshing */

 *coords = ( GeoPnt * ) calloc( nno, sizeof( GeoPnt ) );

 for( i=0; i < nno; i++ )
  {
   /* Verify if the vertex is a boundary one, to get it from original boundary
      vector */

   if( (*no_vector)[i].bdrypt != -1 )
   {
    /* Get the vertex from original boundary vector */

    (*coords)[i].x = bdry_pts[(*no_vector)[i].bdrypt*md].x;
    (*coords)[i].y = bdry_pts[(*no_vector)[i].bdrypt*md].y;
   }

   /* Else get the vertex from node vector */

   else
   {
    /* Convert integer no_vector coordinates to double coordinates */

    vert.ex = (*no_vector)[i].x;
    vert.ey = (*no_vector)[i].y;
    invparam_form( vert.ex, vert.ey, &verta );

    /* Update final output vertex coordinates */

    (*coords)[i].x = verta.x;
    (*coords)[i].y = verta.y;
   }

   /*Release node vector adjacent finite elements list*/

   SllDelAll( (Sll *)&(*no_vector)[i].fadj );
  }

 /* Update output elem.connectivity list from meshing. This list contains
    the number of vertexs followed by the connectivity for each element */

 for( i=0; i < nel; i++ )
   size += ((*el_vector)[i].n + 1);

 *conn = ( int * ) calloc( size, sizeof( int ) );

 for( i=0; i < nel ; i++ )
  {
   /* Update number of vertexs for the element */

   (*conn)[index++] = (*el_vector)[i].n;

   /* Update connectivity for the element */

   for( j=0; j < (*el_vector)[i].n; j++ )
     (*conn)[index++] = ((*el_vector)[i].conect[j])->nno;

   /* Release element vector connectivity memory */

   free( (*el_vector)[i].conect );
  }
}

/* ==================== clear_tree ======================================== */

static void clear_tree(

tree_node_t *node          /* quadtree tree root (in/out) */
)

{
 int i;  /* loop counter */

  /* This routine is recursively done until all tree be released */

  if( node->child != NULL )
     {
      for( i = 0; i < 4; i++ )
       clear_tree( node->child[i] );

      /* Release tree node child vector and the tree node */

      free( node->child );
      free( node );
     }

  else
   {
    /* Release structures inside the tree node and the tree node allocated */

    clear_struct_tree( node );
    free( node );
   }
}

/* ==================== clear_struct_tree ================================= */

static void clear_struct_tree(

tree_node_t *node    /* quadtree tree root (in/out) */
)

{
 mesh_elem_t  *elem;

    /* Release tree node structure */

    switch( node->code )
    {
     case NODE_VERTEX :
      {
       /* Release list of vertexs inside the tree node */

       SllDelAll( (Sll *)&(node->vhead) );

       /* Escape switch */

       break;
      }

     case NODE_INTERIOR :
      {
       /* Release list of finite elements inside the tree node */

       /* Get the first finite element */

       elem = node->elmhead;
       while( elem )
        {
         /* Release connectivity vector*/

         free( elem->conect );

         /* Release the finite element */

         SllDelete( (Sll *)&elem, (Sll)elem );
        }

       /* Escape switch */

       break;
      }

     default :
       break;

    }
}

/*
** ---------------------------------------------------------------------------
** Public Auxiliary quadtree functions:
*/

/* ==================== param_form ======================================== */

void param_form( double x, double y, new_coord *p )
{
 int     sizlev;
 double   dx, dy;

 /* This function converts double values to integer values of a tree point */

 dx        = tree->box.xmax - tree->box.xmin;
 dy        = tree->box.ymax - tree->box.ymin;
 sizlev    = pow_2( 14 );
 p->ex     = sizlev * ( x - tree->box.xmax ) / dx + sizlev;
 p->ey     = sizlev * ( y - tree->box.ymax ) / dy + sizlev;
}


/* ==================== invparam_form ===================================== */

void invparam_form( int ex, int ey, vertex_t *v )
{
 int     sizlev;
 double   dx, dy;

 /* This function converts integer values to double values of a tree point */

 dx        = tree->box.xmax - tree->box.xmin;
 dy        = tree->box.ymax - tree->box.ymin;
 sizlev    = pow_2( 14 );
 v->x      = dx * ( ex - sizlev ) / sizlev + tree->box.xmax;
 v->y      = dy * ( ey - sizlev ) / sizlev + tree->box.ymax;
}

/* ==================== locate_quad ======================================= */

neighbour_t *locate_quad( new_coord *p, tree_node_t *node )
{
 /* This function returns a list of the tree nodes where point p is */

 /* Inicializate the head of the list */

 head_loc = NULL;

 /* return the list of the tree nodes */

 return locate( p, node );
}

/* ==================== locate  =========================================== */

static neighbour_t *locate( new_coord *p, tree_node_t *node )
{
 int  sizlev;

 /* Verify if the tree node is a leaf. If it isn't, traverse the tree to get
    the first tree node where the point p is */

 if( node->child != NULL )
 {
   /* Calculate the tree node size */

   sizlev = pow_2( 13 -  node->depth );

   /* Traverse recursively the  tree  until  get  the  first tree node  where
      the point p is.It could be child[0], child[1], etc.. depending of the p
      position */

   if( p->ex < node->coord.ex + sizlev )
   {
    if( p->ey < node->coord.ey + sizlev )
     return locate( p, node->child[0] );
    else
     return locate( p, node->child[2] );
   }
   else
   {
    if( p->ey < node->coord.ey + sizlev )
     return locate( p, node->child[1] );
    else
     return locate( p, node->child[3] );
   }
 }

 /* Then with the first tree leaf node where the point p is, verify if there's
    another tree nodes where the point p is too ( in some special cases ) */

 else
  return list_loc_nodes (p, node);
}

/* =================== list_loc_nodes ===================================== */

static neighbour_t *list_loc_nodes ( new_coord *p, tree_node_t *node )
{
 neighbour_t *adj;
 int         size;

 /* Insert the first tree leaf node geted by locate function in locate_quad
    tree nodes list returned */

 ins_loc_list ( node );

 /* Verify if there's another tree nodes where the point p is too (some
    special cases) */

 /* First verify left possibilities, after bottom and corner0 */

 if ( p->ex == node->coord.ex && p->ey == node->coord.ey )
 {
  /* Get the adjacent tree left nodes */

  adj = find_adj_node( LEFT, node );

  /* Look trough the adjacent list to see what nodes should be inserted */

  while( adj != NULL )
   {
    /* Calculate the tree node size */

    size = pow_2( 14 - adj->node->depth );

    /* Verify if point p really is inside of this tree node and if so, insert
       into yhe list */

    if( (p->ey >= adj->node->coord.ey) &&
        (p->ey <= adj->node->coord.ey + size) )
     ins_loc_list ( adj->node );

    /* Get the next adjacent left tree nodes and release the actual one */

    SllDelete( (Sll *)&adj, (Sll)adj );
   }

  /* Get the adjacent tree bottom nodes */

  adj = find_adj_node( BOTTOM, node );

  /* Look trough the adjacent list to see what nodes should be inserted */

  while( adj != NULL )
   {
    /* Calculate the tree node size */

    size = pow_2( 14 - adj->node->depth );

    /* Verify if point p really is inside of this tree node and if so, insert
       into the list */

    if( (p->ex >= adj->node->coord.ex) &&
        (p->ex <= adj->node->coord.ex + size) )
     ins_loc_list ( adj->node );

    /* Get the next adjacent bottom tree nodes and release the actual one */

    SllDelete( (Sll *)&adj, (Sll)adj );
   }

  /* Get the adjacent tree corner0 nodes */

  adj = find_adj_node( CORNER0, node );

  /* Look trough the adjacent list to see what nodes should be inserted */

  while( adj != NULL )
   {
    /* Insert the tree node */

    ins_loc_list (adj->node);

    /* Get the next adjacent corner0 tree nodes and release the actual one */

    SllDelete( (Sll *)&adj, (Sll)adj );
   }
 }

 /* Else verify only left possibilities */

 else if ( p->ex == node->coord.ex )
 {
  /* Get the adjacent tree left nodes */

  adj = find_adj_node( LEFT, node );

  /* Look trough the adjacent list to see what nodes should be inserted */

  while( adj != NULL )
   {
    /* Calculate the tree node size */

    size = pow_2( 14 - adj->node->depth );

    /* Verify if point p really is inside of this tree node and if so, insert
       into the list */

    if( (p->ey >= adj->node->coord.ey) &&
        (p->ey <= adj->node->coord.ey + size) )
     ins_loc_list ( adj->node );

    /* Get the next adjacent left tree nodes and release the actual one */

    SllDelete( (Sll *)&adj, (Sll)adj );
   }
 }

 /* Else verify only bottom posssibilities */

 else if ( p->ey == node->coord.ey )
 {
  /* Get the adjacent tree bottom nodes */

  adj = find_adj_node( BOTTOM, node );

  /* Look trough the adjacent list to see what nodes should be inserted */

  while( adj != NULL )
   {
    /* Calculate the tree node size */

    size = pow_2( 14 - adj->node->depth );

    /* Verify if point p really is inside of this tree node and if so, insert
       into the list */

    if( (p->ex >= adj->node->coord.ex) &&
        (p->ex <= adj->node->coord.ex + size) )
     ins_loc_list ( adj->node );

    /* Get the next adjacent bottom tree nodes and release the actual one */

    SllDelete( (Sll *)&adj, (Sll)adj );
   }
 }

 /* Return the head of locate_quad tree nodes list */

 return head_loc;
}

/* =================== ins_loc_list ======================================= */

static void ins_loc_list ( tree_node_t *node )
{
 neighbour_t  *n;

 /* Update locate_quad tree node list that will be returned by the locate_quad
    function */

 n = (neighbour_t *) calloc ( 1, sizeof(neighbour_t) );
 n->node = node;
 n->next = head_loc;
 head_loc = n;
}

/* =================== find_adj_node ====================================== */

neighbour_t *find_adj_node( side_t side, tree_node_t *node )
{
 int         size;
 neighbour_t *adj = NULL;
 neighbour_t *adj_corner = NULL;
 neighbour_t *adj_out = NULL;

 /* This function gets adjacent tree nodes list of a input tree node relative
    an input adjacent desired side( if the input desired side is any  corner,
    then it's returned only one tree node, not a list of them ) */

 /* Calculate the input tree node size */

 size = pow_2( 14 - node->depth );

 /* Find adjacent tree node if the input adjacent desired side is corner0 */

 if( side == CORNER0 )
 {
  /* Get the adjacent left tree node list of the input tree node */

  adj = adj_nodes( LEFT, node );

  /* Look trough the list to get adj tree node that will be returned by this
     function, if the list was not empty */

  if( adj != NULL )
  {
   while ( adj )
   {
    /* Verify if the node in the adjacent left list has the same y cote of the
       input tree node */

    if ( adj->node->coord.ey == node->coord.ey )
    {
     /* Get adjacent bottom tree node list of the adjacent left tree node */

     adj_corner = adj_nodes( BOTTOM, adj->node );

     /* Look trough the list to get adj tree node that will be  returned  by
        this function, if the list was not empty */

     if( adj_corner != NULL )
     {
      while ( adj_corner )
      {
       /* Verify if the tree node in the list has the same x cote of the input
          tree node. If it was, this is the corner0 adjacent desired tree node
          of the input tree node, and then break this loop to update adj, that
          will be the tree node returned by this function */

       if ( node->coord.ex ==
            adj_corner->node->coord.ex + pow_2(14-adj_corner->node->depth ) )
        break;

       /* If the loop was not breaked, get the next tree node of the adjacent
          bottom list and release the actual one */

       SllDelete( (Sll *)&adj_corner, (Sll)adj_corner );
      }
     }

     /* Update the adj tree node, that will be the tree node returned  by this
        function */

     if( adj_corner == NULL )  adj_out = adj_corner;
     else
      {
       adj_out = ( neighbour_t * ) calloc( 1, sizeof( neighbour_t ) );
       adj_out->node = adj_corner->node;
       adj_out->next = NULL;
      }

     /* Release adj tree node list and adj_corner tree node list */

     SllDelAll( (Sll *)&adj );
     SllDelAll( (Sll *)&adj_corner );

     /* Break the loop because the tree node was already found */

     break;
    }

    /* If the loop was not breaked, get the next tree node in the adjacent
       left list and release the actual one*/

    SllDelete( (Sll *)&adj, (Sll)adj );
   }
  }
 }

 /* Find adjacent tree node if the input adjacent desired side is corner1 */

 else if( side == CORNER1 )
 {
  /* Get the adjacent right tree node list of the input tree node */

  adj = adj_nodes( RIGHT, node );

  /* Look trough the list to get adj tree node that will be returned by this
     function, if the list was not empty */

  if( adj != NULL )
  {
   while ( adj )
   {
    /* Verify if the node in the adjacent right list has the same y cote of
       the input tree node */

    if ( adj->node->coord.ey == node->coord.ey )
    {
     /* Get adjacent bottom tree node list of the adjacent right tree node */

     adj_corner = adj_nodes( BOTTOM, adj->node );

     /* Look trough the list to get adj tree node that will be  returned  by
        this function, if the list was not empty */

     if( adj_corner != NULL )
     {
      while ( adj_corner )
      {
       /* Verify if the tree node in the list has the same x cote of the input
          tree node. If it was, this is the corner1 adjacent desired tree node
          of the input tree node, and then break this loop to update adj, that
          will be the tree node returned by this function  */

       if ( adj_corner->node->coord.ex == node->coord.ex + size )
        break;

       /* If the loop was not breaked, get the next tree node of the adjacent
         bottom list and release the actual one */

       SllDelete( (Sll *)&adj_corner, (Sll)adj_corner );
      }
     }

     /* Update the adj tree node, that will be the tree node returned by this
        function */

     if( adj_corner == NULL )  adj_out = adj_corner;
     else
      {
       adj_out = ( neighbour_t * ) calloc( 1, sizeof( neighbour_t ) );
       adj_out->node = adj_corner->node;
       adj_out->next = NULL;
      }

     /* Release adj tree node list and adj_corner tree node list */

     SllDelAll( (Sll *)&adj );
     SllDelAll( (Sll *)&adj_corner );

     /* Break the loop because the tree node was already found */

     break;
    }

    /* If the loop was not breaked, get the next tree node in the adjacent left
       list and release the actual one */

    SllDelete( (Sll *)&adj, (Sll)adj );
   }
  }
 }

 /* Find adjacent tree node if the input adjacent desired side is corner2 */

 else if( side == CORNER2 )
 {
  /* Get the adjacent left tree node list of the input tree node */

  adj = adj_nodes( LEFT, node );

  /* Look trough the list to get adj tree node that will be returned by this
     function, if the list was not empty */

  if( adj != NULL )
  {
   while ( adj )
   {
    /* Verify if the node in the adjacent left list has the same y cote of the
       input tree node */

    if ( node->coord.ey + size ==
         adj->node->coord.ey + pow_2(14-adj->node->depth ) )
    {
     /* Get adjacent top tree node list of the adjacent left tree node */

     adj_corner = adj_nodes( TOP, adj->node );

     /* Look trough the list to get adj tree node that will be  returned  by
        this function, if the list was not empty */

     if( adj_corner != NULL )
     {
      while ( adj_corner )
      {
       /* Verify if tree node in that list has the same x cote of the input
          tree node. If it was, this is the corner2 adjacent desired tree node
          of the input tree node, and then break this loop to update adj, that
          will be the tree node returned by this function  */

       if ( node->coord.ex ==
            adj_corner->node->coord.ex + pow_2(14-adj_corner->node->depth ) )
        break;

       /*If the loop was not breaked, get the next tree node of the adjacent
         bottom list and release the actual one */

       SllDelete( (Sll *)&adj_corner, (Sll)adj_corner );
      }
     }

     /* Update adj tree node, that will be the tree node returned  by  this
        function */

     if( adj_corner == NULL )  adj_out = adj_corner;
     else
      {
       adj_out = ( neighbour_t * ) calloc( 1, sizeof( neighbour_t ) );
       adj_out->node = adj_corner->node;
       adj_out->next = NULL;
      }

     /* Release adj tree node list and adj_corner tree node list */

     SllDelAll( (Sll *)&adj );
     SllDelAll( (Sll *)&adj_corner );

     /* Break the loop because the tree node was already found */

     break;
    }

    /* If the loop was not breaked, get next tree node in the adjacent left
       list and release the actual one */

    SllDelete( (Sll *)&adj, (Sll)adj );
   }
  }
 }

 /* Find adjacent tree node if the input adjacent desired side is corner3 */

 else if( side == CORNER3 )
 {
  /* Get the adjacent right tree node list of the input tree node */

  adj = adj_nodes( RIGHT, node );

  /* Look trough the list to get adj tree node that will be returned by this
     function, if the list was not empty */

  if( adj != NULL )
  {
   while ( adj )
   {
    /* Verify if node in the adjacent right list has the same y cote of the
       input tree node */

    if ( node->coord.ey + size ==
         adj->node->coord.ey + pow_2(14-adj->node->depth ) )
    {
     /* Get the adjacent top tree node list of the adjacent left tree node */

     adj_corner = adj_nodes( TOP, adj->node );

     /* Look trough the list to get adj tree node that will be  returned  by
        this function, if the list was not empty */

     if( adj_corner != NULL )
     {
      while ( adj_corner )
      {
       /* Verify if tree node in that list has the same x cote of the input
          tree node. If it was, this is the corner3 adjacent desired tree node
          of the input tree node, and then break this loop to update adj, that
          will be the tree node returned by this function  */

       if ( adj_corner->node->coord.ex == node->coord.ex + size )
        break;

       /* If the loop was not breaked, get the next tree node of the adjacent
          top list and release the actual one */

       SllDelete( (Sll *)&adj_corner, (Sll)adj_corner );
      }
     }

     /* Update adj tree node, that will be the tree node returned  by  this
        function */

     if( adj_corner == NULL )  adj_out = adj_corner;
     else
      {
       adj_out = ( neighbour_t * ) calloc( 1, sizeof( neighbour_t ) );
       adj_out->node = adj_corner->node;
       adj_out->next = NULL;
      }

     /* Release adj tree node list and adj_corner tree node list */

     SllDelAll( (Sll *)&adj );
     SllDelAll( (Sll *)&adj_corner );

     /* Break the loop because the tree node was already found */

     break;
    }

    /* If the loop was not breaked, get next tree node in the adjacent right
       list and release the actual one */

    SllDelete( (Sll *)&adj, (Sll)adj );
   }
  }
 }


 /* Find adjacent tree node list if the input adjacent desired side  is  not
    any corner, that is, it's right,left,top or bottom */

 else
  adj_out = adj_nodes( side, node );

 /* Return adjacent tree node list of the input tree node relative of an input
    adjacent desired side( if it's any corner, it's only one tree node,  not
    a list ) */

 return adj_out;
}

/* =================== adj_nodes ========================================== */

static neighbour_t *adj_nodes ( side_t side, tree_node_t *node )
{
 int i;
 int depth;
 int *path;

 /* Inicializate the head of the tree adjacent node list that will be returned,
    allocate path and inicaializate adj_type */

 head_neigh = NULL;
 path = ( int * ) calloc ( node->depth, sizeof(int) );
 adj_type = 0;

 /* Get the tree node depth and path in the tree */

 depth = get_path ( node, path );

 /* Get the tree node adjacent list of the input desired side( for right, left,
    bottom,top ) */

 switch ( side )
 {
  /* Get the tree node adjacent list for right input desired side */

  case RIGHT:

   /* Go trough the tree depths to get the tree nodes */

   for ( i = depth-1; i >= 0; i-- )
   {
    if ( path[i] == 0 )
    {
     /* Get tree node related to its path */

     node = get_node ( i, path );

     /* Update tree adjacent node list */

     head_neigh = adjacent ( RIGHT, i, depth, path, node->child[1] );

     /* Release path vector memory */

     free( path );

     /* Return tree adjacent node list */

     return head_neigh;
    }
    else if ( path[i] == 2 )
    {
     /* Get tree node related to its path */

     node = get_node ( i, path );

     /* Update tree adjacent node list */

     head_neigh = adjacent ( RIGHT, i, depth, path, node->child[3] );

     /* Release path vector memory */

     free( path );

     /* Return tree adjacent node list */

     return head_neigh;
    }
   }
   break;

  /* Get the tree node adjacent list for left input desired side */

  case LEFT:

   /* Go trough the tree depths to get the tree nodes */

   for ( i = depth-1; i >= 0; i-- )
   {
    if ( path[i] == 1 )
    {
     /* Get tree node related to its path */

     node = get_node ( i, path );

     /* Update the tree adjacent node list */

     head_neigh = adjacent ( LEFT, i, depth, path, node->child[0] );

     /* Release path vector memory */

     free( path );

     /* Return tree adjacent node list */

     return head_neigh;
    }
    else if ( path[i] == 3 )
    {
     /* Get tree node related to its path */

     node = get_node ( i, path );

     /* Update the tree adjacent node list */

     head_neigh = adjacent ( LEFT, i, depth, path, node->child[2] );

     /* Release path vector memory */

     free( path );

     /* Return tree adjacent node list */

     return head_neigh;
    }
   }
   break;

  /* Get the tree node adjacent list for bottom input desired side */

  case BOTTOM:

   /* Go trough the tree depths to get the tree nodes */

   for ( i = depth-1; i >= 0; i-- )
   {
    if ( path[i] == 2 )
    {
     /* Get tree node related to its path */

     node = get_node ( i, path );

     /* Update the tree adjacent node list */

     head_neigh = adjacent ( BOTTOM, i, depth, path, node->child[0] );

     /* Release path vector memory */

     free( path );

     /* Return tree adjacent node list */

     return head_neigh;
    }
    else if ( path[i] == 3 )
    {
     /* Get tree node related to its path */

     node = get_node ( i, path );

     /* Update the tree adjacent node list */

     head_neigh = adjacent ( BOTTOM, i, depth, path, node->child[1] );

     /* Release path vector memory */

     free( path );

     /* Return tree adjacent node list */

     return head_neigh;
    }
   }
   break;

  /* Get the tree node adjacent list for top input desired side */

  case TOP:

   /* Go trough the tree depths to get the tree nodes */

   for ( i = depth-1; i >= 0; i-- )
   {
    if ( path[i] == 0 )
    {
     /* Get tree node related to its path */

     node = get_node ( i, path );

     /* Update the tree adjacent node list */

     head_neigh = adjacent ( TOP, i, depth, path, node->child[2] );

     /* Release path vector memory */

     free( path );

     /* Return tree adjacent node list */

     return head_neigh;
    }
    else if ( path[i] == 1 )
    {
     /* Get tree node related to its path */

     node = get_node ( i, path );

     /* Update the tree adjacent node list */

     head_neigh = adjacent ( TOP, i, depth, path, node->child[3] );

     /* Release path vector memory */

     free( path );

     /* Return tree adjacent node list */

     return head_neigh;
    }
   }
   break;

  /* Default scape */

  default:
   break;
  }

 /* Release path vector memory */

 free( path );

 /* Return NULL if no tree adjacent node was found */

 return NULL;
}

/* =================== get_path =========================================== */

static int get_path( tree_node_t *node, int *path )
{
 int d;

 /* Inicializate the depth */

 d = 0;

 /* Get the tree node depth and path */

 path_node( node->depth, node->coord.ex, node->coord.ey, &d, path,
            tree->root );

 /* Return tree node depth */

 return d;
}

/* =================== path_node ========================================== */

static void path_node( unsigned level, int ex, int ey, int *d,
                       int *path, tree_node_t *node)
{
 /* Verify if the node->depth is less than level */

 if( node->depth < level )
 {
   int size;

   /* Calculate the tree node child size */

   size = pow_2( 13 - node->depth );

   /* Get the tree node path acording with its coordinates */

   if( ex < node->coord.ex + size )
   {
     if( ey < node->coord.ey + size )
     {
       /* Update path */

       path[(*d)++] = 0;
       path_node ( level, ex, ey, d, path, node->child[0] );
     }
     else
     {
       /* Update path */

       path[(*d)++] = 2;
       path_node ( level, ex, ey, d, path, node->child[2] );
     }
   }
   else
   {
     if( ey < node->coord.ey + size )
     {
       /* Upadte path */

       path[(*d)++] = 1;
       path_node ( level, ex, ey, d, path, node->child[1] );
     }
     else
     {
       /* Upadte path */

       path[(*d)++] = 3;
       path_node ( level, ex, ey, d, path, node->child[3] );
     }
   }
 }
}

/* =================== get_node =========================================== */

static tree_node_t *get_node ( int n, int *path )
{
 int         i;
 tree_node_t *node;

 /* Get tree node relative a input depth( n ) and a input path traversing the
    depths into the tree  */

 node = tree->root;
 for ( i = 0; i < n; i++ )
  node = node->child[path[i]];

 /* Return tree node */

 return node;
}

/* =================== adjacent =========================================== */

static neighbour_t *adjacent( side_t side, int n, int depth, int *path,
                              tree_node_t *node )
{
 int i;

 /* Go trough the depths */

 for ( i = n + 1; i < depth; i++ )
 {
  /* If tree node is a leaf, update adj_type and break the loop */

  if ( node->child == NULL )
  {
   adj_type = 1;
   break;
  }

  /* Get tree node relative to the path and side */

  switch ( side )
  {
   case RIGHT:
    if ( path[i] == 1 )
     node = node->child[0];
    else
     node = node->child[2];
    break;
   case LEFT:
    if ( path[i] == 0 )
     node = node->child[1];
    else
     node = node->child[3];
    break;
   case BOTTOM:
    if ( path[i] == 0 )
     node = node->child[2];
    else
     node = node->child[3];
    break;
   case TOP:
    if ( path[i] == 2 )
     node = node->child[0];
    else
     node = node->child[1];
    break;
   default:
    break;
  }
 }

 /* Add the tree node into the tree node adjacent list that will be returned */

 add_nodes ( side, node );

 /* Return the tree node adjacent list( its head pointer ) */

 return head_neigh;
}

/* =================== add_nodes ========================================== */

static void add_nodes ( side_t side, tree_node_t *node )
{
 neighbour_t *adj;

 /* If tree node is a leaf, allocate memory and put into tree adjacent node
    list */

 if ( node->child == NULL )
 {
  adj = ( neighbour_t * ) calloc ( 1, sizeof(neighbour_t) );
  adj->node = node;
  adj->next = head_neigh;
  head_neigh = adj;
 }

 /* Else done recursively this function for tree node childs to update  the
    tree adjacent node list */

 else
 {
  switch ( side )
  {
   case RIGHT:
    add_nodes ( RIGHT, node->child[0] );
    add_nodes ( RIGHT, node->child[2] );
    break;
   case LEFT:
    add_nodes ( LEFT, node->child[1] );
    add_nodes ( LEFT, node->child[3] );
    break;
   case BOTTOM:
    add_nodes ( BOTTOM, node->child[2] );
    add_nodes ( BOTTOM, node->child[3] );
    break;
   case TOP:
    add_nodes ( TOP, node->child[0] );
    add_nodes ( TOP, node->child[1] );
    break;
   default:
    break;
  }
 }
}

/*
** ---------------------------------------------------------------------------
** Public Auxiliary quadtree math functions:
*/

/* =================== pow_2 =============================================== */

int pow_2( int n )
{
 int p = 1;

 while( n )
 {
   p *= 2;
   n--;
 }
 return p;
}

/* =================== cross_product ====================================== */

double cross_product( vector *v1, vector *v2 )
{
 return ( v1->x * v2->y - v2->x * v1->y );
}

/* =================== pseudo_ang ========================================= */

double pseudo_ang ( vector *v1, vector *v2 )
{
 double arg;

 arg = ( v1->x*v2->x + v1->y*v2->y ) /
       ( sqrt( v1->x*v1->x + v1->y*v1->y ) *
         sqrt( v2->x*v2->x + v2->y*v2->y ) );
 return  cross_product( v1, v2 ) >= 0. ? ( 1. - arg ) : ( 3. + arg );
}

/* =================== dot_product ======================================== */

double dot_product( vector *v1, vector *v2 )
{
  return ( v1->x * v2->x ) + ( v1->y * v2->y );
}



typedef struct Tinfo_
{
  int id;
  int i;
  int j;
  GeoPnt coord;
} Tinfo;

/*
** Function to compare node x-y coordinates to quadratic triangulation and
** quadtree generation.
*/
static int fn_quad_sort (void *c1, void *c2)
{
 Tinfo *info1 = c1, *info2 = c2;

 /* if (info1 == NULL || info2 == NULL)          return -1; */

 if      (info1->i < info2->i)      return -1;
 else if (info1->i > info2->i)      return  1;
 else if (info1->j < info2->j)      return -1;
 else if (info1->j > info2->j)      return  1;
 else                               return  0;
}

/*
** Quadratic mesh generation.
**
**  This function updates the mesh structure from triangulation and quadtree
**  generation, to assure quadratic mesh. It is included quadratic mid nodes
**  in each element side, and in a binary tree to an eficient search.
**  it returns 1 if all quadratic generation it was done.
*/
static void MshQuadratic ( GeoPnt **coords, int **conn, int *nn, int nel,
                           Tbtree *root, int newnn, GeoPnt  *midcoords)
{
 Tinfo   *f;
 Tinfo   *fe;
 Tinfo   *finfo;
 int     *midconn;
 int      prenel, pronel;
 int      new_, i, j;
 int      first, second;

 /* alloc auxiliary mid coords and conn structure */
 /* midcoords = (CooPnt *) calloc (2*nel*((*conn)[0]+1), sizeof(CooPnt));  */
 midconn   = (int *) calloc (2*nel*((*conn)[0]+1), sizeof (int));

 /* traverse all elements to update auxiliary mid nodes and coords structure */
 new_ = prenel = 0;
 pronel = prenel + (*conn)[prenel] + 1;
 for (i = 0; i < nel; i++)
 {
  midconn[new_++] = 2*(*conn)[prenel];
  for (j = prenel+1; j < pronel; j++)
  {
   /* find element mshtri_edge middle coordinates */
   first  = j;
   second = (j+1)==pronel ? prenel+1 : j+1;
   f  = (Tinfo *) calloc (1, sizeof(Tinfo));
   f->i = (*conn)[first]<(*conn)[second]?(*conn)[first]:(*conn)[second];
   f->j = (*conn)[first]>(*conn)[second]?(*conn)[first]:(*conn)[second];
   f->coord.x = ((*coords)[(*conn)[first]].x + (*coords)[(*conn)[second]].x) / 2.0;
   f->coord.y = ((*coords)[(*conn)[first]].y + (*coords)[(*conn)[second]].y) / 2.0;

   /* look in the tree to see if a new mshtri_edge middle point that will be inserted
      already exists */
   fe = BtreeFind (root, f);

   /* update btree and auxiliary coords and conn structure */
   if (fe == NULL)
   {
    finfo  = (Tinfo *) calloc (1, sizeof(Tinfo));
    finfo->i = f->i;
    finfo->j = f->j;
    finfo->coord.x = f->coord.x;
    finfo->coord.y = f->coord.y;
    finfo->id = newnn;
    midcoords[newnn-(*nn)].x  = finfo->coord.x;
    midcoords[newnn-(*nn)].y  = finfo->coord.y;
    midconn[new_++]   = (*conn)[j];
    midconn[new_++]   = newnn;
    root  = BtreeInsert (root, finfo);
    newnn++;
   }
   else
   {
    midconn[new_++]   = (*conn)[j];
    midconn[new_++]   = fe->id;
   }
   if (f)
   {
    free (f);
    f = 0L;
   }
  }

  /*update present nel and proxim nel */
  prenel = pronel;
  pronel = prenel + (*conn)[prenel] + 1;
 }

 /* update real coords and conn structure */
 *coords = (GeoPnt *) realloc ( *coords, newnn*sizeof(GeoPnt) );
 for (i = *nn; i < newnn; i++)
   (*coords)[i] = midcoords[i-(*nn)];
 *conn = realloc ( *conn, new_*sizeof(int) );
 for (i = 0; i < new_; i++)
   (*conn)[i] = midconn[i];

 /* update total number of nodes */
 *nn = newnn;

 /* free memory */
 if (midcoords)
 {
  free (midcoords);
  midcoords = 0L;
 }
 if (midconn)
 {
  free (midconn);
  midconn = 0L;
 }
}


/*=========================  Quadratic =============================*/
static void Quadratic(GeoPnt **coords, int **conn, int *n_node, int n_elem,
                int n_loops, int *loop_segs, GeoPnt *bdry_pts, face_t *cnew)
{
  Tbtree    *root = NULL;
  Tinfo     *finfo;
  GeoPnt    *midcoords;
  int        j,i;
  int        newnn=*n_node;

  /* inicializate btree comparative function */
  BtreeInit (fn_quad_sort);
  midcoords = (GeoPnt *) calloc (2*n_elem*((*conn)[0]+1), sizeof(GeoPnt));

  for( i=0; i<n_loops; i++)
  {
    for( j=1; j<(loop_segs[i]/2); j++)
    {
      finfo  = (Tinfo *) calloc (1, sizeof(Tinfo));
      finfo->i = cnew->vt[j-1]->bdrypt;
      finfo->j = cnew->vt[j]->bdrypt;
      finfo->coord.x = midcoords[newnn-(*n_node)].x =  bdry_pts[(cnew->vt[j]->bdrypt*2)-1].x ;
      finfo->coord.y = midcoords[newnn-(*n_node)].y =  bdry_pts[(cnew->vt[j]->bdrypt*2)-1].y ;
      finfo->id = newnn;
      root  = BtreeInsert (root, finfo);
      newnn++;
    }
  finfo  = (Tinfo *) calloc (1, sizeof(Tinfo));
    finfo->i = cnew->vt[0]->bdrypt;
    finfo->j = cnew->vt[cnew->n-1]->bdrypt;
    finfo->coord.x = midcoords[newnn-(*n_node)].x =  bdry_pts[(cnew->vt[j-1]->bdrypt*2)+1].x ;
    finfo->coord.y = midcoords[newnn-(*n_node)].y =  bdry_pts[(cnew->vt[j-1]->bdrypt*2)+1].y ;
    finfo->id = newnn;
    root  = BtreeInsert (root, finfo);
    newnn++;
  cnew=cnew->next;
  }
/*      free(Bdry_pts); */
  MshQuadratic ( coords, conn, n_node, n_elem, root, newnn, midcoords);
}





static Tbtree *root_quad;
static list *cedgequad;

/* =========================== NewEdge ========================== */
static void NewEdge(int ed1, int ed2)
{
 Tinfo *ed,copy;
 copy.i=(ed1<ed2)?ed1:ed2;
 copy.j=(ed1>ed2)?ed1:ed2;
 copy.id=0;
 if((ed=BtreeFind(root_quad,&copy))==NULL)
 {
  ed=(Tinfo *) calloc (1, sizeof(Tinfo));
  *ed=copy;
  root_quad=BtreeInsert(root_quad,ed);
 }
 ed->id=ed->id+1;
}

/* =========================== bdrypts_quad ========================== */
static void EdgeQuad(void *p)
{
 Tinfo *ed=p;
 if(ed->id==1)
 {
  list   *mnew = NULL;
  mnew=(list *)calloc(1,sizeof(list));
  mnew->vi=ed->i;
  mnew->vj=ed->j;
  mnew->next=cedgequad;
  cedgequad=mnew;
 }
}

/* =========================== bdrypts_quad ========================== */
static void bdrypts_quad( tree_node_t *node)
{

 int           i, j;
  mesh_elem_t  *f;

  if( node->child != NULL )
    for( i = 0; i < 4; i++ ) bdrypts_quad( node->child[i]);
  else
  {
    f = node->elmhead;
    while( f )
    {
      /* Add edges of the element */
      for( j = 0; j < f->n; j++ )
       NewEdge(f->conect[j]->nno,f->conect[(j+1)%f->n]->nno);

      f = f->next;
    }
  }
}

/* =========================== add_bdrypts_quad ========================== */
void add_bdrypts_quad(
tree_node_t *node,
list    **head         /* cedge boundary of quad structure           (out) */
)
{

 cedgequad=NULL;
 root_quad=NULL;
 /* Initialize Btree */
 BtreeInit (fn_quad_sort);
 /* Build Btree */
 bdrypts_quad(node);
 /* Build boundary of quad */
 if(root_quad==NULL)
 {
  (*head)=NULL;
 }
 else
 {
  BtreeTraverse(root_quad,EdgeQuad);
  (*head)=cedgequad;
 }

}

/* =========================== fix_mesh_orientation ========================== */
static void
fix_mesh_orientation (int nnode, double *coords, int n_elem, int *conn)
{
  int i, j, index;

  /* for all elements */
  index = 0;
  for (i = 0; i < n_elem; ++i)
  {
    double area = 0.0;

    /* obtain area */
    for (j = 0; j < conn[index]; ++j)
    {
      int id_i = conn[index+j+1];
      int id_j = conn[index+((j+1)%conn[index])+1];
      area = area + coords[id_i*2+0]*coords[id_j*2+1] -
                    coords[id_j*2+0]*coords[id_i*2+1];
    }

    /* invert orientation of area in case of wrong area */
    if (area > 0)
    {
      int pts_swap;
      int type_mesh = conn[index];
      for (j = 0; j < type_mesh/2; ++j)
      {
        pts_swap = conn[index+j+1];
        conn[index+j+1] = conn[index+type_mesh-j];
        conn[index+type_mesh-j] = pts_swap;
      }
    }

    /* update index */
    index = index + conn[index] + 1;
  }

}
