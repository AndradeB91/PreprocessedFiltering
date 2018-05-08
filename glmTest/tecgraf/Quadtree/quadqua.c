/*
** ---------------------------------------------------------------------------
**
** quadqua.c - This module contains routines to generate quadrangular elements
**             in the part between quadtree and boundary by a procedure of
**             boundary contraction.
**
** ---------------------------------------------------------------------------
**
** Public Function:
**
** ---------------------------------------------------------------------------
**
** Version: 0-001
**
** Created:  01-Jun-96          Eduardo Setton Sampaio Silveira &&
**                              Joaquim Bento Cavalcante Neto
**
** Supervised:                  Luiz Fernando Martha
**
** ---------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "quadtree.h"
#include "quadsll.h"

#define QUADQUA_C

#define API acos(-1.0)
#define ATOL 0.001
#define RTOL 1
#define BTOL 1e-05

#define TOPL		0
#define END		1

#define ABS(n)  (((n)>=0)?(n):(-(n)))
#define DIST(cx,cy,px,py) sqrt(((cx-px)*(cx-px))+((cy-py)*(cy-py)))

/*
** ---------------------------------------------------------------------------
** Definition for use of algorithms:
*/

#define QUAD_NORMAL  1
#define QUAD_DEBUG   0
#define QUAD_ALLQUAD 0
#define QUAD_INTERS  0

/*
** ---------------------------------------------------------------------------
** Includes for DEBUG:
*/

#if QUAD_DEBUG
#include "mquadraw.h"
#include "mquadgra.h"
static int  flag_debug = 1;
#endif

/*
** ---------------------------------------------------------------------------
** Local functions prototypes:
*/

static void   msh_qua_vertex( mesh_vertex_t *,face_t *, int, mesh_vertex_t ** );
static void   msh_add_bdrypts( face_t *, adjnode **, list ** );
static int    msh_gen_regintpts( list *, list *, adjnode **, tree_node_t *,
                                 mesh_vertex_t **, mesh_elem_t * );
static void   msh_gen_regintpts_elem( list *, list *, list *, adjnode **,tree_node_t *,
                                      mesh_vertex_t **, mesh_elem_t * );
static void   msh_quad_best( list *, mesh_vertex_t ** );
static void   msh_quad_mesh( list *, mesh_vertex_t **, int *, mesh_elem_t ** );
static list  *msh_cedge_end( list * );
static list  *msh_cedge_add( int, int, list * );
static void   msh_update_cedge( tree_node_t *, int, adjnode **, list ** );
static void   msh_bound_qua_mesh( face_t *, mesh_vertex_t *, mesh_elem_t **,
				  int *,int *,mesh_vertex_t **,mesh_elem_t ** );
static void   msh_quad_qua_mesh( int *,int *, tree_node_t *, mesh_vertex_t **,
                                 mesh_elem_t ** );
static void   insert_tree_list( tree_node_t *, nodeadj_t ** );
#if 0
static void   insert_list( neighbour_t *, nodeadj_t ** );
static void   coord_int( tree_node_t * , side_t, nodeadj_t ** );
#endif
static void   insert_prob_list( neighbour_t *, nodeadj_t ** );
static void   coord_prob_int( tree_node_t * , side_t, nodeadj_t ** );
static void   locate_prob( neighbour_t *, neighbour_t ** );
static void   locate_prob_list( tree_node_t *, side_t, neighbour_t ** );
static void   insert_prob_tree_list( tree_node_t *node, nodeadj_t **list_head);
static int    find_list( tree_node_t *, side_t, nodeadj_t ** );
static void   int_quadtree( list *, list *, nodeadj_t *, nodeadj_t *, adjnode **,
                            mesh_vertex_t **, mesh_elem_t * );
static void   find_vert_id( list *, list *, nodeadj_t *, nodeadj_t *, tree_node_t *,
                            adjnode **, mesh_vertex_t **, mesh_elem_t * );
static void   find_quad_id( list *, list *, nodeadj_t *, nodeadj_t *, tree_node_t *,
                            side_t, adjnode **, mesh_vertex_t **, mesh_elem_t * );
#if !QUAD_NORMAL
static int    filtro_quad( tree_node_t *, side_t );
static int    filtro_quad_test( tree_node_t *, side_t, new_coord,
                                int, int, int, int, int, int );
#endif
#if QUAD_INTERS
static int    intersec( adjnode **, mesh_vertex_t ** );
#else
static int    intersec( list *,list *, mesh_vertex_t ** );
#endif
#if defined(QUAD_INQUAD) || QUAD_ALLQUAD
static int    interqtree( nodeadj_t *,   mesh_vertex_t ** );
#endif
#if QUAD_INTERS
static int    intersec_list( nodeadj_t * ,adjnode **,mesh_vertex_t **,list ** );
static int    intersec_vert( int , adjnode **, mesh_vertex_t **, list ** );
static int    intersec_quad( int, int, mesh_vertex_t ** );
static int    find_intersec( int, int, list ** );
#endif
static int    intersec_cross( int ,int  ,int, int, mesh_vertex_t ** );
static float  intersec_cross_value( int ,int ,int, mesh_vertex_t ** );
#ifdef QUAD_INSIDE
static int    inside( mesh_vertex_t **, mesh_elem_t * );
#endif
static int    inside_polygon( int, vector *, float, float );
static int    boundr_polygon( int, vector *, float, float );
static double geoangle ( void );
static void   quad_edge( tree_node_t *, new_coord, new_coord, int, int,
                         adjnode **, list ** );
static void   update_edge( int, int, adjnode **, list ** ,int);
#if QUAD_INTERS
static void   add_adjnode( list  *, adjnode ** );
static void   del_adjnode( int, int, adjnode ** );
#endif
static int    test_if_isbound();
static int    find_edge_boundary(list *,int,int,int *);
static int    find_meshv( int, int, int, int, face_t *, mesh_vertex_t *,
			  mesh_vertex_t ** );
static double calcangle ( int, int, int, new_coord [4]);
static double getsvalue ( new_coord [4]);
static int    quadtree_edge(int,int);

/*
** ---------------------------------------------------------------------------
** Local static variables:
*/

static list      *listquad   = NULL; /*auxiliary listquad list*/
static nodeadj_t *listnode   = NULL; /*tree nodes list*/
static nodeadj_t *listvertex = NULL; /*base edge vertex adjacent tree nodes list*/
static double      maxang;           /*current maximum angle for each edge*/
static double      ang;              /*computed current angle*/
static new_coord   first;            /*integer coordinates of first vertex base edge*/
static new_coord   second;           /*integer coordinates of second vertex base edge*/
static new_coord   candidate;        /*integer coordinates of candidate vertex*/
static new_coord   iddate;           /*integer coordinates of choosen vertex*/
static int         firstid;          /*id number of the first vertex base edge*/
static int         secondid;         /*id number of the second vertex base edge*/
static int         beforeid;         /*id number of the second vertex base edge*/
static int         candid;           /*id number of the candidate vertex*/
static int         id;               /*id number of the choosen vertex*/
static int         nnos;
static int         nelems;
static int         nelemint;
static int         step;

static int         fid;
static int         sid;
static int         tid;
static int         fic;
static int         sic;
static int         tic;
static int         fquadid;
static new_coord   fquadiddate; /*integer coordinates of choosen vertex*/
static int         squadid;
static new_coord   squadiddate; /*integer coordinates of choosen vertex*/
static int         curfirstid;
static int         cursecondid;
static new_coord   curfirst;
static new_coord   cursecond;

int nnogeometric;

/*
** ---------------------------------------------------------------------------
** Public Function:
*/

/* =========================  quadmesh_bound  ============================= */

int quadmesh_bound(

tree_node_t    *node,      /* quadtree tree root                       (in) */
face_t         *face,      /* geometrical boundary structure           (in) */
mesh_vertex_t  *vhe,       /* quadtree finite element vertexs structure(in) */
int            *nno,       /* total finite element vertexs number      (in) */
int            nnogeom,    /* geometry finite element vertexs number   (in) */
int            *nel,       /* total finite elements number         (in/out) */
mesh_vertex_t  **no_vector,/* finite element vertexs structure        (out) */
mesh_elem_t    **el_vector /* finite elements structure               (out) */

)

{
 /*Local variables and symbols*/

 list           *cedge = NULL;       /*cedge boundary structure*/
 list           *qedge = NULL;       /* cedge boundary quad structure */
 adjnode        **adj_no_vector;     /*adjacent edges vertex structure*/
 mesh_elem_t    *list_vector = NULL; /*auxiliary finite elements structure*/
 int            i;                   /*loop counter*/

 nelemint = *nel;

 nnogeometric = nnogeom;
 nnos = nelems = 0;

 /* The boundary contraction is first done  for the  geometrical  contour
    and after for the holes that could exist  in the geometry.The routine
    finishes when all the triangulation of the geometry was done. */

  /* Update final output finite element vertexs structure */

  msh_qua_vertex( vhe, face, *nno, no_vector );

  /* Inicializate adjacent finite element vertexs structure */

  adj_no_vector = ( adjnode ** )calloc( *nno,sizeof( adjnode * ) );
  for( i = 0; i < *nno; i++ )   adj_no_vector[i] = NULL;

  /* This routine only first put the geometrical informations in the auxi-
     liary cedge(contract edge) structure, that will be used  to   do  the
     boundary contraction,where cedge is the head of cedge structure. */

  msh_add_bdrypts( face, adj_no_vector, &cedge );
  add_bdrypts_quad( node, &qedge);

  /* For each edge in the cedge structure, get the internal points    and
     makes the boundary contraction generating the finite elements, atua-
     lizating the cedge structure. This step finishes  when  the    cedge
     structure is empty */

  while( cedge != NULL )
   {

     if( cedge->next == NULL )
      i = 1;

    /* Generate the possible internal points,  choose  the best  triangle,
       atualizate the id for the choosen point that will be used  to  make
       the new finite element in the element_vector structure and   return
       the pointer for the edge auxiliary structure that contains the  two
       new edges,one or none to atualizate the cegde boundary structure */

    if( !msh_gen_regintpts( cedge, qedge, adj_no_vector, node, no_vector,
                            list_vector ))
    {
     SllDelAll( (Sll *)&cedge );
     free( adj_no_vector );
     return 0;
    }

    /* Choose the best quadrilateral if there's more then one */

    msh_quad_best( cedge, no_vector );

    /* With the choosen point and your identificator( id ) updated, atua-
       lizate the list for  the   finite  elements that in the final will
       make the output finite element structure */

    msh_quad_mesh( cedge, no_vector, nel, &list_vector );

    /* Atualizate the cedge boundary structure, that is, insert the   two
       new created edges, only one or none, and delete the base the  head
       cedge structure used as a base edge */

    msh_update_cedge( node, nnogeom, adj_no_vector, &cedge );

    /*Get the new edge in the cedge structure */

   }

  /* Now update the final output finite element vector, adding  the  finite
     elements done in the quadtree and finite elements done by the boundary
     contraction with Delaunay */

    /* Update finite elements done by boundary contraction */

    msh_bound_qua_mesh( face, vhe, &list_vector, nel, nno, no_vector, el_vector );

    /* Update finite elements done by quadtree */

    msh_quad_qua_mesh( nno, nel, node, no_vector, el_vector );

  /* Release static global auxiliary lists used */

  SllDelAll( (Sll *)&listquad );
  SllDelAll( (Sll *)&listnode );

  /* Release adjacent finite element vertexs structure */

  free( adj_no_vector );

  /* Return the indicated generated mesh */

  return 1;
}

/*
** ---------------------------------------------------------------------------
** Local functions:
*/

/*=========================  msh_qua_vertex ==================================*/

static void msh_qua_vertex(

mesh_vertex_t *v,
face_t        *f,
int           nno,
mesh_vertex_t **no_vector
)

{
 int           nnogeomaux = 0;
 int           j;
 new_coord     p;
#if 0
 vertex_t      va;   /* tirar depois */
 char          text[5];
#endif
 int           maxalloc = 100000;


    /* Allocate finite element vertexs structure */

    /* *no_vector = (mesh_vertex_t *) calloc( nno, sizeof( mesh_vertex_t ) ); */
    *no_vector = (mesh_vertex_t *) calloc( maxalloc, sizeof( mesh_vertex_t ) );

    /* Put the quadtree vertexs in the finite element vertexs structure */

    while( v )
    {
      /* Update vertex */
      if (v->nno >=  maxalloc)
      {
        maxalloc = v->nno * 2;
        *no_vector = (mesh_vertex_t *) realloc ( *no_vector, maxalloc * sizeof( mesh_vertex_t ) );
      }


      (*no_vector)[v->nno].x      = v->x;
      (*no_vector)[v->nno].y      = v->y;
      (*no_vector)[v->nno].status = 1;
      (*no_vector)[v->nno].bdrypt = -1;
      (*no_vector)[v->nno].fadj   = NULL;
      (*no_vector)[v->nno].next   = NULL;
      (*no_vector)[v->nno].nno    = v->nno;

      /* Draw node - take off later */

      #if 0
      if( v->nno < (get_nnoint() + get_nnogeomint()) )
      {
       p.ex = v->x;
       p.ey = v->y;
       invparam_form( p.ex, p.ey, &va );
       sprintf( text, "%4d", v->nno );
       for( i = 0; text[i] == ' '; i++ );
       GraText( va.x, va.y, &text[i] );
      }
      #endif

      /* Get next vertex */

      v = v->next;
    }

    /* Put the geometry nodes in the finite element vertexs structure */

    while( f != NULL )
    {
     for( j = 0; j < f->n; j++ )
     {
      /* Transform vertex float in integer float as in quadtree vertexs */

      param_form( f->vt[j]->x, f->vt[j]->y, &p );

      /* Update vertex */
      if (nnogeomaux >=  maxalloc)
      {
        maxalloc = nnogeomaux * 2;
        *no_vector = (mesh_vertex_t *) realloc ( *no_vector, maxalloc * sizeof( mesh_vertex_t ) );
      }

      (*no_vector)[nnogeomaux].x       = p.ex;
      (*no_vector)[nnogeomaux].y       = p.ey;
      (*no_vector)[nnogeomaux].status  = 1;
      (*no_vector)[nnogeomaux].bdrypt  = f->vt[j]->bdrypt;
      (*no_vector)[nnogeomaux].fadj    = NULL;
      (*no_vector)[nnogeomaux].next    = NULL;
      (*no_vector)[nnogeomaux].nno     = nnogeomaux;
      nnogeomaux++;

      /* Draw node - take off later */

      #if 0
      va.x = f->vt[j]->x;
      va.y = f->vt[j]->y;
      sprintf( text, "%4d", f->vt[j]->bdrypt );
      for( i = 0; text[i] == ' '; i++ );
      GraText( va.x, va.y, &text[i] );
      #endif
     }

     /* Get the next hole, if exists */

     f = f->next;
    }
}

/* =========================  msh_add_bdrypts  =============================*/

static void msh_add_bdrypts(

face_t  *f,            /* geometrical boundary structure      (in) */
adjnode **adj_vector,  /* adjacent edges vertex structure (in/out) */
list    **head         /* cedge boundary structure           (out) */
)

{
 /*Local variables and symbols*/

 list   *mnew = NULL; /* new cedge boundary structure element */
 int    nnoaux = 0;   /* auxiliary vertex number */
 int    j;            /* loop counter */


 /* This routine only first put the geometrical informations in the   auxi-
   liary cedge(contract edge) structure, that will be used to do the boun-
   dary contraction.The cedge is the pointer for the cegde structure  that
   is a single-linked list */

  /* Make the inicial cedge boundary structure and the adjacent edges vertex
     structure */

  while( f != NULL )
  {
   /* For each vertex of the geometrical structure, make the cedge struc-
      ture and its identificator( id ) */

   for( j = 0; j < f->n ; j++ )
    {
     /* Allocate the new cell */

     SllAddEnd( (Sll *)head, sizeof( list ), (Sll *)&mnew );

     /* Insert the edge being created in the cedge structure */

     mnew->vi = ( j % f->n ) + nnoaux;
     mnew->vj = ( (j+1) % f->n ) + nnoaux;
     mnew->type = QUAD_BOUNDARY;

     /* Update adjacent edges vertex structure */

     #if QUAD_INTERS
     add_adjnode( mnew, &adj_vector[mnew->vi] );
     add_adjnode( mnew, &adj_vector[mnew->vj] );
     #endif
    }

   /*Increment number of vertexs */

   nnoaux = nnoaux + f->n;

   /*Get the next hole if exist */

   f = f->next;
  }
}

/* ======================== msh_gen_regintpts ============================= */

static int msh_gen_regintpts (

list        *head ,          /* cedge boundary structure              (in) */
list        *head_quad ,     /* qedge boundary structure              (in) */
adjnode     **adj_vector,    /* adjacent edges vertex structure       (in) */
tree_node_t *node,           /* root quadtree tree node               (in) */
mesh_vertex_t **node_vector, /* finite element vertex structure       (in) */
mesh_elem_t    *list_vector  /* auxiliay finite element structure     (in) */
)

{
 int   mid, vj;
 list *elem;
 int   rid;

 /* Initialize the step for getting internal points */

 step = 0;

 /* Initialize the insertion of temporary edges in cedge to find the ids */

 fic = sic = tic = 0;

 /* Initialize the first, second and third ids to be found */

 fid = sid = tid = -1;

 /* Find the first, second and third ids if it's possible */

 while( step < 3 )
 {
  switch (step)
  {
   case 0:
    msh_gen_regintpts_elem( head, head_quad, head, adj_vector, node, node_vector,
                            list_vector );
    fid = id;
    step++;
   break;
   case 1:
    if ( !find_edge_boundary( head, head->vi, fid, &mid ) &&
         !quadtree_edge( head->vi,fid )  )
    {
     elem = msh_cedge_add( head->vi, fid, head );
     if( elem != NULL )
     {
      msh_gen_regintpts_elem( head, head_quad, elem, adj_vector, node, node_vector,
                              list_vector );
      sid = id;
      sic = 1;
     }
    }
    step++;
   break;
   case 2:
    if( head->type == QUAD_BOUNDARY )
     vj = head->next->vj;
    else
     vj = head->vj;
    if ( !find_edge_boundary( head, fid, vj, &mid ) &&
         !quadtree_edge( fid, vj ) )
    {
     elem = msh_cedge_add( fid, vj, head );
     if( elem != NULL )
     {
      msh_gen_regintpts_elem( head, head_quad, elem, adj_vector, node, node_vector,
                              list_vector );
      tid = id;
      tic = 1;
     }
    }
    step++;
   break;
  }
 }

 /* If the base edge didn't get the first id, it's not possible to generate
    the mesh, that is, return zero */

 if( fid == -1 )
  rid = 0;
 else
  rid = 1;

 /* Return the indicated generated polygon */

 return rid;
}

/* ==========================  msh_gen_regintpts_elem  ==================== */

static void msh_gen_regintpts_elem (

list        *head ,          /* cedge boundary structure              (in) */
list        *head_quad ,     /* cedge boundary quad structure         (in) */
list        *elem ,          /* cedge boundary structure element      (in) */
adjnode     **adj_vector,    /* adjacent edges vertex structure       (in) */
tree_node_t *node,           /* root quadtree tree node               (in) */
mesh_vertex_t **no_vector,   /* finite element vertex structure       (in) */
mesh_elem_t    *list_vector  /* auxiliay finite element structure     (in) */
)

{
 /*Local variables and symbols*/

 neighbour_t  *locfirst = NULL;
 neighbour_t  *locfirst_prob = NULL;
 neighbour_t  *locsecond = NULL;
 neighbour_t  *locsecond_prob = NULL;

  /* Inicializate tree nodes candidates vertex list */

  listvertex = NULL;

  /* Inicializate the maximum angle */

  maxang = 0.0;

  /* Inicializate the id */

  id = -1;

  /* Get the first and second coordinates and ids of the base edge */

  if (elem->type == QUAD_INTERNAL)
  {
   curfirst.ex = (*no_vector)[elem->vi].x;
   curfirst.ey = (*no_vector)[elem->vi].y;
   curfirstid = elem->vi;

   cursecond.ex = (*no_vector)[elem->vj].x;
   cursecond.ey = (*no_vector)[elem->vj].y;
   cursecondid = elem->vj;
  }
  else
  {
   curfirst.ex = (*no_vector)[elem->vi].x;
   curfirst.ey = (*no_vector)[elem->vi].y;
   curfirstid = elem->vi;

   beforeid =  elem->vj;

   cursecond.ex = (*no_vector)[elem->next->vj].x;
   cursecond.ey = (*no_vector)[elem->next->vj].y;
   cursecondid = elem->next->vj;
  }

  if (step == 0)
  {
   first.ex = curfirst.ex;
   first.ey = curfirst.ey;
   second.ex = cursecond.ex;
   second.ey = cursecond.ey;
   firstid  = curfirstid;
   secondid = cursecondid;
  }

  /* Get the tree nodes list where the first and second vertexs are */

  locfirst  = locate_quad( &curfirst, node );
  locsecond = locate_quad( &cursecond, node );

  /* Make a list for all tree nodes if wasn't done yet */

  if( listnode == NULL )
   insert_tree_list( node, &listnode );

  /* If the base edge didn't have any id yet, get the adjacent tree nodes list
     of the first and second vertexs based on the tree nodes list where the
     vertexs are */

  if( id == -1 )
  {
   #if 0
   /* Make the new list following description above */

   insert_list( locfirst, &listvertex );
   insert_list( locsecond, &listvertex );

   /* Look trough the list to get the id for the best triangle of the base
     edge, based on a Delaunay triangulation */

   int_quadtree( head, listvertex, listnode, adj_vector, no_vector, list_vector );
   #endif
  }

  /* If the  base  edge  didn't get an id, make the best triangle  considering
     a new list of tree nodes to get the best vertex to do it. This  list   is
     maked in that way: get the first and second base edge vertexs tree  nodes
     and go through both tree nodes all adjacents only stoping when get a tree
     node NODE_INTERIOR  or  NODE_VERTEX that is updated in the new list */

  if( id == -1 )
  {
   #if 0
   /* Release the list before used and that didn't get and id */

   SllDelAll( (Sll *)&listvertex );

   /* Make the new list following description above */

   insert_prob_list( locfirst, &listvertex );
   insert_prob_list( locsecond, &listvertex );

   /* Look trough the lists to get the id for the best triangle of the base
      edge, based on a Delaunay triangulation */

   int_quadtree( head, listvertex, listnode, adj_vector, no_vector, list_vector );
   #endif
  }

  /* If the  base  edge  didn't get an id, make the best triangle  considering
     a new list of tree nodes to get the best vertex to do it. This  list   is
     maked in that way: get the first and second base edge vertexs tree  nodes
     and their tree nodes adjacents right,left, top and botton.Then go through
     those  tree  nodes  all  adjacents  only stoping  when  get  a  tree node
     NODE_INTERIOR  or  NODE_VERTEX that is updated in the new list */

  if( id == -1 )
  {
   /* Release the list before used and that didn't get and id */

   SllDelAll( (Sll *)&listvertex );

   /* Get a new locfirst and locsecond list */

   locate_prob( locfirst, &locfirst_prob );
   locate_prob( locsecond, &locsecond_prob );

   /* Make the new list following description  above */

   insert_prob_list( locfirst_prob, &listvertex );
   insert_prob_list( locsecond_prob, &listvertex );

   /* Look trough the lists to get the id for the best triangle of the base
      edge, based on a Delaunay triangulation */

   int_quadtree( head, head_quad, listvertex, listnode, adj_vector, no_vector, list_vector );
  }

  /* If the base edge even didn't get and a id, you now should look trough a
     new list  that  have all quadtree tree nodes  that  are   NODE_INTERIOR
     or NODE_VERTEX */

  if ( id == -1 )
  {
   /* Release the list before used and that didn't get and id */

   SllDelAll( (Sll *)&listvertex );

   /* Make the new list following description above */

   insert_prob_tree_list ( node, &listvertex );

   /* Look trough the lists to get the id for the best triangle of the base
     edge, based on a Delaunay triangulation */

   int_quadtree( head, head_quad, listvertex, listnode, adj_vector, no_vector, list_vector );
  }

  /* Release memory used */

  SllDelAll( (Sll *)&locfirst );
  SllDelAll( (Sll *)&locfirst_prob );
  SllDelAll( (Sll *)&locsecond_prob );
  SllDelAll( (Sll *)&locsecond );
  SllDelAll( (Sll *)&listvertex );
}

/* ======================= msh_quad_best ================================== */

static void msh_quad_best(

list          *head,         /* cedge boundary structure             (in) */
mesh_vertex_t **no_vector    /* finite element vertexs structure     (in) */

)

{
 int   code = 0;
 list *elem;
 new_coord quad[4];
 double firstquad, secondquad;

 /* Find the code for the case of polygon */

 if( (tid == -1) && (sid == -1) )
  code = 0;                              /* one triangular   polygon  */
 else if( (tid == -1) )
  code = 1;                              /* one quadrangular polygon  */
 else if( (sid == -1) )
  code = 2;                              /* one quadrangular polygon  */
 else if( (tid != -1) && (sid != -1) )
  code = 3;                              /* two quadrangular polygons */

 /* Choose the polygon depending on the case */

 switch( code )
 {
  case 0:
   /* Choose the triangular polygon */

   fquadid = fid;
   fquadiddate.ex = (*no_vector)[fid].x;
   fquadiddate.ey = (*no_vector)[fid].y;
   squadid = -1;
  break;

  case 1:
   /* Look the form of quadrangular polygon */

   quad[0].ex = (*no_vector)[firstid].x;
   quad[0].ey = (*no_vector)[firstid].y;
   quad[1].ex = (*no_vector)[sid].x;
   quad[1].ey = (*no_vector)[sid].y;
   quad[2].ex = (*no_vector)[fid].x;
   quad[2].ey = (*no_vector)[fid].y;
   quad[3].ex = (*no_vector)[secondid].x;
   quad[3].ey = (*no_vector)[secondid].y;
   firstquad = getsvalue( quad );

   /* Choose the quadrangular or triangular polygon */

   if( firstquad <= 1.5 ) {
    fquadid = fid;
    fquadiddate.ex = (*no_vector)[fid].x;
    fquadiddate.ey = (*no_vector)[fid].y;
    squadid = sid;
    squadiddate.ex = (*no_vector)[sid].x;
    squadiddate.ey = (*no_vector)[sid].y;
   }
   else {
    fquadid = fid;
    fquadiddate.ex = (*no_vector)[fid].x;
    fquadiddate.ey = (*no_vector)[fid].y;
    squadid = -1;
   }
  break;

  case 2:
   /* Look the form of quadrangular polygon */

   quad[0].ex = (*no_vector)[firstid].x;
   quad[0].ey = (*no_vector)[firstid].y;
   quad[1].ex = (*no_vector)[fid].x;
   quad[1].ey = (*no_vector)[fid].y;
   quad[2].ex = (*no_vector)[tid].x;
   quad[2].ey = (*no_vector)[tid].y;
   quad[3].ex = (*no_vector)[secondid].x;
   quad[3].ey = (*no_vector)[secondid].y;
   firstquad = getsvalue( quad );

   /* Choose the quadrangular or triangular polygon */

   if( firstquad <= 1.5 ) {
    fquadid = tid;
    fquadiddate.ex = (*no_vector)[tid].x;
    fquadiddate.ey = (*no_vector)[tid].y;
    squadid = fid;
    squadiddate.ex = (*no_vector)[fid].x;
    squadiddate.ey = (*no_vector)[fid].y;
   }
   else {
    fquadid = fid;
    fquadiddate.ex = (*no_vector)[fid].x;
    fquadiddate.ey = (*no_vector)[fid].y;
    squadid = -1;
   }
  break;

  case 3:
   /* Look the form of first quadrangular polygon */

   quad[0].ex = (*no_vector)[firstid].x;
   quad[0].ey = (*no_vector)[firstid].y;
   quad[1].ex = (*no_vector)[sid].x;
   quad[1].ey = (*no_vector)[sid].y;
   quad[2].ex = (*no_vector)[fid].x;
   quad[2].ey = (*no_vector)[fid].y;
   quad[3].ex = (*no_vector)[secondid].x;
   quad[3].ey = (*no_vector)[secondid].y;
   firstquad = getsvalue( quad );

   /* Look the form of second quadrangular polygon */

   quad[0].ex = (*no_vector)[firstid].x;
   quad[0].ey = (*no_vector)[firstid].y;
   quad[1].ex = (*no_vector)[fid].x;
   quad[1].ey = (*no_vector)[fid].y;
   quad[2].ex = (*no_vector)[tid].x;
   quad[2].ey = (*no_vector)[tid].y;
   quad[3].ex = (*no_vector)[secondid].x;
   quad[3].ey = (*no_vector)[secondid].y;
   secondquad = getsvalue( quad );

   /* Choose the best quadrangular or the triangular polygon */

   if(  (firstquad <= secondquad) && (firstquad <= 1.5) ) {
    fquadid = fid;
    fquadiddate.ex = (*no_vector)[fid].x;
    fquadiddate.ey = (*no_vector)[fid].y;
    squadid = sid;
    squadiddate.ex = (*no_vector)[sid].x;
    squadiddate.ey = (*no_vector)[sid].y;
   }
   else if( (secondquad <= 1.5) ) {
    fquadid = tid;
    fquadiddate.ex = (*no_vector)[tid].x;
    fquadiddate.ey = (*no_vector)[tid].y;
    squadid = fid;
    squadiddate.ex = (*no_vector)[fid].x;
    squadiddate.ey = (*no_vector)[fid].y;
   }
   else {
    fquadid = fid;
    fquadiddate.ex = (*no_vector)[fid].x;
    fquadiddate.ey = (*no_vector)[fid].y;
    squadid = -1;
   }
  break;
 }

 /* Delete the two temporary edges if they were included in the end of
    cedge to find the ids */

 if( sic ) {
  elem = msh_cedge_end( head );
  SllDelete( (Sll *)&head, (Sll)elem );
 }
 if( tic ) {
  elem = msh_cedge_end( head );
  SllDelete( (Sll *)&head, (Sll)elem );
 }
}

/* ======================= msh_quad_mesh ================================== */

static void msh_quad_mesh(

list          *head,         /* cedge boundary structure               (in) */
mesh_vertex_t **no_vector,   /* finite element vertexs structure       (in) */
int           *nel,          /* total number of finite elements    (in/out) */
mesh_elem_t   **list_vector  /* auxiliary finite element structure (in/out) */

)

{
 mesh_elem_t       *mnew;
 int               i;
 new_coord         p;
 vertex_t          va;
#if QUAD_DEBUG
 point             *el_coord[4], mesh_coord[4];
#endif

 /* Create the new finite element polygon and put in the auxiliary
    finite element polygon list */

 /* Allocate the new finite element polygon */

 if( SllAddEnd( (Sll *)list_vector, sizeof(mesh_elem_t), (Sll *)&mnew) )
  {
   /* Allocate finite element polygon connectivity */

   if( squadid != -1 )
    mnew->conect = ( mesh_vertex_t ** ) calloc( 4, sizeof(mesh_vertex_t *) );
   else
    mnew->conect = ( mesh_vertex_t ** ) calloc( 3, sizeof(mesh_vertex_t *) );

   /* Update connectivity in the new finite element polygon */

   if( squadid != -1 )
   {
    mnew->conect[0] = &((*no_vector)[firstid]);
    mnew->conect[1] = &((*no_vector)[squadid]);
    mnew->conect[2] = &((*no_vector)[fquadid]);
    mnew->conect[3] = &((*no_vector)[secondid]);
   }
   else
   {
    mnew->conect[0] = &((*no_vector)[firstid]);
    mnew->conect[1] = &((*no_vector)[fquadid]);
    mnew->conect[2] = &((*no_vector)[secondid]);
   }

   /* Update number of vertexs in the new finite element polygon */

   if( squadid != -1 )
    mnew->n = 4;
   else
    mnew->n = 3;

   /* Update the new finite element number and total number of finite
      elements polygons */

   mnew->nel = (*nel)++;

   /* Update total number of elements polygons */

   nelems += mnew->n;

   /* Draw element test */

   for( i = 0; i < mnew->n; i++ )
    {
     p.ex = mnew->conect[i]->x;
     p.ey = mnew->conect[i]->y;
     invparam_form( p.ex, p.ey, &va );
#if QUAD_DEBUG
     mesh_coord[i].x = va.x;
     mesh_coord[i].y = va.y;

     el_coord[i] = &mesh_coord[i];
#endif
    }

    #if QUAD_DEBUG
    if( flag_debug )
     {
      Draw_Face( HOLLOW1, 5, mnew->n, (point **) el_coord );
      ItfUtlMsg( "Debug" );
     }
    #endif
  }
}

/* ======================== msh_update_cedge ============================== */

static void msh_update_cedge(

tree_node_t *node,       /* quadtree tree node root                    (in) */
int         nnogeom,     /* geometry finite element vertexs number     (in) */
adjnode     **adj_vector,/* adjacent edges vertex structure        (in/out) */
list        **head       /* cedge boundary structure               (in/out) */
)

{
 int pmidle;

 /* Makes the atualization( contraction )  of  the  cegde boundary structure.
    There are two possibilities for each edge: the first is  that  if the two
    ids are equal or bigger than nnogeom,it indicates that this is a internal
    edge of the quadtree and should  not be included in  the cedge  structure
    by the update_edge function,but when there are two nodes between the edge
    and they are not NODE_INTERIOR,then even both ids  being equal or  bigger
    then nnogeom, the  edge  should  be  included  in  the  cedge  structure;
    the second possibility is when one or  other or  both ids are lower  then
    nnogeom, then the edge should be  included directly in the cedge boundary
    structure by the update_edge function.In relation of the base edge, every
    time it should be deleted. */

   /* Delete the firstid and secondid edge ( base edge ) */

   if ((*head)->type == QUAD_BOUNDARY )
   {
    update_edge( firstid, beforeid , adj_vector, head, TOPL);
    update_edge( beforeid, secondid, adj_vector, head, TOPL);
   }
   else
   {
    update_edge( firstid, secondid , adj_vector, head, TOPL);
   }

   /* Update other edges depending if there's squadid or not */

   if( squadid != -1 )
   {
    /* Update the firstid and squadid edge */

    if( firstid >= nnogeom && squadid >= nnogeom )
      quad_edge( node, first, squadiddate, firstid, squadid, adj_vector, head );
    else if (find_edge_boundary(*head,firstid,squadid,&pmidle))
    {
     update_edge( firstid, pmidle , adj_vector, head, TOPL);
     update_edge( pmidle, squadid, adj_vector, head, TOPL);
    }
    else
     update_edge( firstid, squadid, adj_vector, head, TOPL);

    /* Update the fquadid e squadid edge */

    if( fquadid >= nnogeom && squadid >= nnogeom )
     quad_edge( node, squadiddate, fquadiddate, squadid, fquadid, adj_vector,
                head );
    else if (find_edge_boundary(*head,squadid,fquadid,&pmidle))
    {
     update_edge( squadid, pmidle , adj_vector, head, TOPL);
     update_edge( pmidle, fquadid, adj_vector, head, TOPL);
    }
    else
     update_edge( squadid, fquadid, adj_vector, head, TOPL);

    /* Update the secondid e squadid edge */

    if( secondid >= nnogeom && fquadid >= nnogeom )
     quad_edge( node, fquadiddate, second, fquadid, secondid, adj_vector,
                head );
    else if (find_edge_boundary(*head,fquadid,secondid,&pmidle))
    {
     update_edge( fquadid, pmidle , adj_vector, head, TOPL);
     update_edge( pmidle, secondid, adj_vector, head, TOPL);
    }
    else
     update_edge( fquadid, secondid, adj_vector, head, TOPL);
   }
   else
   {
    /* Update the firstid e fquadid edge */

    if( firstid >= nnogeom && fquadid >= nnogeom )
      quad_edge( node, first, fquadiddate, firstid, fquadid, adj_vector, head );
    else if (find_edge_boundary(*head,firstid,fquadid,&pmidle))
    {
     update_edge( firstid, pmidle , adj_vector, head, TOPL);
     update_edge( pmidle, fquadid, adj_vector, head, TOPL);
    }
    else
     update_edge( firstid, fquadid, adj_vector, head, TOPL);

    /* Update the secondid e fquadid edge */

    if( secondid >= nnogeom && fquadid >= nnogeom )
     quad_edge( node, fquadiddate, second, fquadid, secondid, adj_vector,
                head );
    else if (find_edge_boundary(*head,fquadid,secondid,&pmidle))
    {
     update_edge( fquadid, pmidle , adj_vector, head, TOPL);
     update_edge( pmidle, secondid, adj_vector, head, TOPL);
    }
    else
     update_edge( fquadid, secondid, adj_vector, head, TOPL);
   }
}

/* ======================= msh_bound_qua_mesh ============================= */

static void msh_bound_qua_mesh(

face_t         *f,          /* geometrical boundary structure          (in) */
mesh_vertex_t  *vhe,        /* quadtree finite element vertex structure(in) */
mesh_elem_t   **list_vector,/* auxiliary finite elements structure     (in) */
int            *nel,        /* total finite elements number            (in) */
int            *nno,        /* total finite element vertexs number (in/out) */
mesh_vertex_t **no_vector,  /* output finite element vertexs structure (in) */
mesh_elem_t   **el_vector   /* output finite elements structure    (in/out) */
)

{
 int  j,k,l,m,n,idnvertex,idelem,npts,idcenter,incenter,node[9];
 int  e1,e2,e3,e4 = 0,ef;
 int  vi, vj, v[4];
 int  sxb,syb,xb,yb,nx,ny;
 float dist1, dist2,s;
 vector pl[4];
 int   nnoaux;
 new_coord         p;
 vertex_t          va;
#if QUAD_DEBUG
#if 0
 point             *el_coord[5],mesh_coord[5];
#endif
#endif
 int               ntotalelems;

 /* Get number of nodes */

 nnoaux    = *nno;
 nnos      = nnoaux;

 /* Get number of elements */

 ntotalelems = nelems + nelemint;

 /* Allocate output finite element structure */

 *el_vector = (mesh_elem_t *) calloc( 4*ntotalelems, sizeof( mesh_elem_t ) );

 /* Look trough list vector to update elem_vector */

 idelem = nelemint;
 while( *list_vector )
 {
  npts   = 4;

  /* Update connectivity in the the new output finite element */

  /* Compute baricent and allocate midle edge point if doesn't exist */

  sxb = 0.0;
  syb = 0.0;
  for (j = 0; j < (*list_vector)->n; j++)
  {
   node[j] = (*list_vector)->conect[j]->nno;
   sxb += (*list_vector)->conect[j]->x;
   syb += (*list_vector)->conect[j]->y;
   nx = ((*list_vector)->conect[j]->x +
         (*list_vector)->conect[(j+1)%((*list_vector)->n)]->x) / 2;
   ny = ((*list_vector)->conect[j]->y +
         (*list_vector)->conect[(j+1)%((*list_vector)->n)]->y) / 2;
   vi = (*list_vector)->conect[j]->nno;
   vj = (*list_vector)->conect[(j+1)%((*list_vector)->n)]->nno;
   idnvertex = find_meshv( nx, ny, vi, vj, f, vhe, no_vector );
   if (idnvertex == -1)
   {
    /* Update vertex */

    (*no_vector)[nnoaux].x       = nx;
    (*no_vector)[nnoaux].y       = ny;
    (*no_vector)[nnoaux].status  = 1;
    (*no_vector)[nnoaux].bdrypt  = -1;
    (*no_vector)[nnoaux].fadj    = NULL;
    (*no_vector)[nnoaux].next    = NULL;
    (*no_vector)[nnoaux].nno     = nnoaux;
    nnoaux++;
    node[j+(*list_vector)->n] = nnoaux-1;
    sxb +=nx;
    syb +=ny;
   }
   else
   {
    node[j+(*list_vector)->n] = idnvertex;
    sxb +=(*no_vector)[idnvertex].x;
    syb +=(*no_vector)[idnvertex].y;
   }
  }

  xb = sxb/(*list_vector)->n/2;
  yb = syb/(*list_vector)->n/2;

  /* Allocate baricent of polygon */

  if( (*list_vector)->n == 4 )
  {
   for( m = 0; m < (*list_vector)->n; m++ )
   {
    pl[m].x = (*list_vector)->conect[m]->x;
    pl[m].y = (*list_vector)->conect[m]->y;
   }
   incenter = inside_polygon( (*list_vector)->n, pl, xb, yb );
   idcenter = boundr_polygon( (*list_vector)->n, pl, xb, yb );
   if( !incenter || idcenter )
   {
    /* Compute the new baricent from the lowest diagonal of polygon */

    dist1 = DIST((*list_vector)->conect[0]->x,(*list_vector)->conect[0]->y,
                 (*list_vector)->conect[2]->x,(*list_vector)->conect[0]->y);
    dist2 = DIST((*list_vector)->conect[1]->x,(*list_vector)->conect[1]->y,
                 (*list_vector)->conect[3]->x,(*list_vector)->conect[3]->y);
    if( dist1 < dist2 )
    {
     xb = ((*list_vector)->conect[0]->x + (*list_vector)->conect[2]->x) / 2.0;
     yb = ((*list_vector)->conect[0]->y + (*list_vector)->conect[2]->y) / 2.0;
    }
    else
    {
     xb = ((*list_vector)->conect[1]->x + (*list_vector)->conect[3]->x) / 2.0;
     yb = ((*list_vector)->conect[1]->y + (*list_vector)->conect[3]->y) / 2.0;
    }
   }
   /* Update vertex */

   (*no_vector)[nnoaux].x       = xb;
   (*no_vector)[nnoaux].y       = yb;
   (*no_vector)[nnoaux].status  = 1;
   (*no_vector)[nnoaux].bdrypt  = -1;
   (*no_vector)[nnoaux].fadj    = NULL;
   (*no_vector)[nnoaux].next    = NULL;
   (*no_vector)[nnoaux].nno     = nnoaux;
   nnoaux++;
   node[2*(*list_vector)->n] = nnoaux - 1;
  }
  else
  {
   for( m = 0; m < (*list_vector)->n; m++ )
   {
    pl[m].x = (*list_vector)->conect[m]->x;
    pl[m].y = (*list_vector)->conect[m]->y;
   }
   idcenter = boundr_polygon( (*list_vector)->n, pl, xb, yb );
   if( !idcenter )
   {
    /* Update vertex */

    (*no_vector)[nnoaux].x       = xb;
    (*no_vector)[nnoaux].y       = yb;
    (*no_vector)[nnoaux].status  = 1;
    (*no_vector)[nnoaux].bdrypt  = -1;
    (*no_vector)[nnoaux].fadj    = NULL;
    (*no_vector)[nnoaux].next    = NULL;
    (*no_vector)[nnoaux].nno     = nnoaux;
    nnoaux++;
    node[2*(*list_vector)->n] = nnoaux - 1;
   }
  }

  /* Make the new finite elements from finite elements polygons */

  if( (*list_vector)->n == 4 )
  {
   e1 = idelem;
   e2 = idelem+1;
   e3 = idelem+2;
   e4 = idelem+3;

   /* Allocate connectivity in the the new output finite element */

   (*el_vector)[e1].conect =
     ( mesh_vertex_t **) calloc(npts,sizeof( mesh_vertex_t *));

   /* Verify  orientation because this algorithm gets finite
     elements in  orientation */

   v[0] = 0;
   v[1] = 7;
   v[2] = 8;
   v[3] = 4;
   s = 0.0;
   for( n = 0; n < (*list_vector)->n; n++ )
   {
    vi = node[v[n]];
    vj = node[v[(n+1)%((*list_vector)->n)]];
    s += (*no_vector)[vi].x * (*no_vector)[vj].y -
         (*no_vector)[vj].x * (*no_vector)[vi].y ;
   }

   /* Update new output finite element in  orientation */

   if( s > 0.0 )
   {
    (*el_vector)[e1].conect[0]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e1].conect[1]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e1].conect[2]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e1].conect[3]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e1].nel         = e1;
    (*el_vector)[e1].n           = 4;
   }
   else
   {
    (*el_vector)[e1].conect[0]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e1].conect[1]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e1].conect[2]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e1].conect[3]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e1].nel         = e1;
    (*el_vector)[e1].n           = 4;
   }

   /* Allocate connectivity in the the new output finite element */

   (*el_vector)[e2].conect =
     ( mesh_vertex_t **) calloc(npts,sizeof( mesh_vertex_t *));

   /* Verify  orientation because this algorithm gets finite
     elements in  orientation */

   v[0] = 4;
   v[1] = 8;
   v[2] = 5;
   v[3] = 1;
   s = 0.0;
   for( n = 0; n < (*list_vector)->n; n++ )
   {
    vi = node[v[n]];
    vj = node[v[(n+1)%((*list_vector)->n)]];
    s += (*no_vector)[vi].x * (*no_vector)[vj].y -
         (*no_vector)[vj].x * (*no_vector)[vi].y ;
   }

   /* Update new output finite element in  orientation */

   if( s > 0.0 )
   {
    (*el_vector)[e2].conect[0]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e2].conect[1]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e2].conect[2]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e2].conect[3]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e2].nel         = e2;
    (*el_vector)[e2].n           = 4;
   }
   else
   {
    (*el_vector)[e2].conect[0]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e2].conect[1]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e2].conect[2]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e2].conect[3]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e2].nel         = e2;
    (*el_vector)[e2].n           = 4;
   }

   /* Allocate connectivity in the the new output finite element */

   (*el_vector)[e3].conect =
     ( mesh_vertex_t **) calloc(npts,sizeof( mesh_vertex_t *));

   /* Verify  orientation because this algorithm gets finite
     elements in  orientation */

   v[0] = 8;
   v[1] = 6;
   v[2] = 2;
   v[3] = 5;
   s = 0.0;
   for( n = 0; n < (*list_vector)->n; n++ )
   {
    vi = node[v[n]];
    vj = node[v[(n+1)%((*list_vector)->n)]];
    s += (*no_vector)[vi].x * (*no_vector)[vj].y -
         (*no_vector)[vj].x * (*no_vector)[vi].y ;
   }

   /* Update new output finite element in  orientation */

   if( s > 0.0 )
   {
    (*el_vector)[e3].conect[0]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e3].conect[1]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e3].conect[2]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e3].conect[3]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e3].nel         = e3;
    (*el_vector)[e3].n           = 4;
   }
   else
   {
    (*el_vector)[e3].conect[0]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e3].conect[1]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e3].conect[2]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e3].conect[3]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e3].nel         = e3;
    (*el_vector)[e3].n           = 4;
   }

   /* Allocate connectivity in the the new output finite element */

   (*el_vector)[e4].conect =
     ( mesh_vertex_t **) calloc(npts,sizeof( mesh_vertex_t *));

   /* Verify  orientation because this algorithm gets finite
     elements in  orientation */

   v[0] = 7;
   v[1] = 3;
   v[2] = 6;
   v[3] = 8;
   s = 0.0;
   for( n = 0; n < (*list_vector)->n; n++ )
   {
    vi = node[v[n]];
    vj = node[v[(n+1)%((*list_vector)->n)]];
    s += (*no_vector)[vi].x * (*no_vector)[vj].y -
         (*no_vector)[vj].x * (*no_vector)[vi].y ;
   }

   /* Update new output finite element in  orientation */

   if( s > 0.0 )
   {
    (*el_vector)[e4].conect[0]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e4].conect[1]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e4].conect[2]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e4].conect[3]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e4].nel         = e4;
    (*el_vector)[e4].n           = 4;
   }
   else
   {
    (*el_vector)[e4].conect[0]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e4].conect[1]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e4].conect[2]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e4].conect[3]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e4].nel         = e4;
    (*el_vector)[e4].n           = 4;
   }
  }
  else
  {
   e1 = idelem;
   e2 = idelem+1;
   e3 = idelem+2;

   /* Allocate connectivity in the the new output finite element */

   (*el_vector)[e1].conect =
     ( mesh_vertex_t **) calloc(npts,sizeof( mesh_vertex_t *));

   /* Verify  orientation because this algorithm gets finite
     elements in  orientation */

   v[0] = 0;
   v[1] = 5;
   v[2] = 6;
   v[3] = 3;
   s = 0.0;
   for( n = 0; n < (*list_vector)->n; n++ )
   {
    vi = node[v[n]];
    vj = node[v[(n+1)%((*list_vector)->n)]];
    s += (*no_vector)[vi].x * (*no_vector)[vj].y -
         (*no_vector)[vj].x * (*no_vector)[vi].y ;
   }

   /* Update new output finite element in  orientation */

   if( s > 0.0 )
   {
    (*el_vector)[e1].conect[0]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e1].conect[1]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e1].conect[2]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e1].conect[3]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e1].nel         = e1;
    (*el_vector)[e1].n           = 3;
   }
   else
   {
    (*el_vector)[e1].conect[0]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e1].conect[1]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e1].conect[2]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e1].conect[3]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e1].nel         = e1;
    (*el_vector)[e1].n           = 3;
   }

   /* Allocate connectivity in the the new output finite element */

   (*el_vector)[e2].conect =
     ( mesh_vertex_t **) calloc(npts,sizeof( mesh_vertex_t *));

   /* Verify  orientation because this algorithm gets finite
     elements in  orientation */

   v[0] = 3;
   v[1] = 6;
   v[2] = 4;
   v[3] = 1;
   s = 0.0;
   for( n = 0; n < (*list_vector)->n; n++ )
   {
    vi = node[v[n]];
    vj = node[v[(n+1)%((*list_vector)->n)]];
    s += (*no_vector)[vi].x * (*no_vector)[vj].y -
         (*no_vector)[vj].x * (*no_vector)[vi].y ;
   }

   /* Update new output finite element in  orientation */

   if( s > 0.0 )
   {
    (*el_vector)[e2].conect[0]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e2].conect[1]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e2].conect[2]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e2].conect[3]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e2].nel         = e2;
    (*el_vector)[e2].n           = 3;
   }
   else
   {
    (*el_vector)[e2].conect[0]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e2].conect[1]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e2].conect[2]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e2].conect[3]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e2].nel         = e2;
    (*el_vector)[e2].n           = 3;
   }

   /* Allocate connectivity in the the new output finite element */

   (*el_vector)[e3].conect =
     ( mesh_vertex_t **) calloc(npts,sizeof( mesh_vertex_t *));

   /* Verify  orientation because this algorithm gets finite
     elements in  orientation */

   v[0] = 6;
   v[1] = 5;
   v[2] = 2;
   v[3] = 4;
   s = 0.0;
   for( n = 0; n < (*list_vector)->n; n++ )
   {
    vi = node[v[n]];
    vj = node[v[(n+1)%((*list_vector)->n)]];
    s += (*no_vector)[vi].x * (*no_vector)[vj].y -
         (*no_vector)[vj].x * (*no_vector)[vi].y ;
   }

   /* Update new output finite element in  orientation */

   if( s > 0.0 )
   {
    (*el_vector)[e3].conect[0]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e3].conect[1]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e3].conect[2]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e3].conect[3]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e3].nel         = e3;
    (*el_vector)[e3].n           = 3;
   }
   else
   {
    (*el_vector)[e3].conect[0]   = &((*no_vector)[node[v[0]]]);
    (*el_vector)[e3].conect[1]   = &((*no_vector)[node[v[1]]]);
    (*el_vector)[e3].conect[2]   = &((*no_vector)[node[v[2]]]);
    (*el_vector)[e3].conect[3]   = &((*no_vector)[node[v[3]]]);
    (*el_vector)[e3].nel         = e3;
    (*el_vector)[e3].n           = 3;
   }
  }

  /* Draw element test */

  if( (*list_vector)->n == 4 ) ef = e4;
  else                         ef = e3;
  for( l = e1; l <= ef; l++ )
  {
   for( k = 0; k < 4; k++ )
    {
     p.ex =  (*el_vector)[l].conect[k]->x;
     p.ey =  (*el_vector)[l].conect[k]->y;
     invparam_form( p.ex, p.ey, &va );
#if QUAD_DEBUG
#if 0
     mesh_coord[k].x = va.x;
     mesh_coord[k].y = va.y;
     el_coord[k] = &mesh_coord[k];
#endif
#endif
    }

    #if QUAD_DEBUG
    #if 0
    if( flag_debug )
     {
      Draw_Face( HOLLOW1, 5, 4, (point **) el_coord );
      ItfUtlMsg( "Debug" );
     }
    #endif
    #endif

   /* Update number of vertexs in the new output finite element */

   (*el_vector)[l].n = npts;

   /* Update the new output finite element number */

   (*el_vector)[l].nel = l;
  }

  /* Update number of node and elements */

  nnos = nnoaux - 1;
  *nno = nnos+1;
  idelem += (*list_vector)->n;

  /* Release list_vector element */

  free( (*list_vector)->conect );
  SllDelete( (Sll *)list_vector, (Sll)*list_vector );
 }
 *nel = ntotalelems;
}

/* ======================= msh_quad_qua_mesh ============================== */

static void msh_quad_qua_mesh(

int         *nno,
int         *nel,
tree_node_t *node,         /*quadtree tree node root             (in)*/
mesh_vertex_t **no_vector,
mesh_elem_t **el_vector  /*output finite element structure (in/out)*/
)

{
 mesh_elem_t       *elme;
 int               i;
 int               j;
 int               vi;
 int               vj;
 float             s;

 /* Traverse function recursivily until get all quadtree finite elements */

 /* Verify if the tree node has no child, that is, if it's a leaf and then put
    its finite element in output finite elements structure */

 if( node->child != NULL )
   for( i = 0; i < 4; i++  )
    msh_quad_qua_mesh( nno, nel, node->child[i], no_vector, el_vector );
 else if( node->elmhead != NULL )
      {
        elme = node->elmhead;
        while( elme )
        {
         /* Allocate connectivity in the the new output finite element */

         (*el_vector)[elme->nel].conect =
           ( mesh_vertex_t **) calloc( elme->n, sizeof( mesh_vertex_t *) );

         /* Verify  orientation because this algorithm gets finite
	    elements in  orientation */

         s = 0.0;
	 for( j = 0; j < elme->n; j++ )
	 {
	  vi = elme->conect[j]->nno;
	  vj = elme->conect[(j+1)%(elme->n)]->nno;

	  s += (*no_vector)[vi].x * (*no_vector)[vj].y -
	       (*no_vector)[vj].x * (*no_vector)[vi].y ;
         }

	 /* Update new output finite element in  orientation */

	 if (s > 0.0)
	 {
	  for( j = 0; j < elme->n; j++ )
	   (*el_vector)[elme->nel].conect[j] = elme->conect[(elme->n-1) - j];
         }
	 else
	 {
	  for( j = 0; j < elme->n; j++ )
	   (*el_vector)[elme->nel].conect[j] = elme->conect[j];
	 }

         /* Update number of vertexs in the new output finite element */

         (*el_vector)[elme->nel].n = elme->n;

         /* Update the new output finite element number */

         (*el_vector)[elme->nel].nel = elme->nel;

         /* Get another finite element in the tree node */

         elme = elme->next;
        }
      }
}

/*
** ---------------------------------------------------------------------------
** Local auxiliary functions:
*/

/* ======================== insert_tree_list ============================== */

static void insert_tree_list(

tree_node_t *node,
nodeadj_t   **list_head
)

{
 int               i;
 nodeadj_t        *le;

 /* Traverse all quadtree to insert in a list all tree nodes that have finite
    elements inside */

 if( node->child != NULL )
  for( i = 0; i < 4; i++ )
   insert_tree_list( node->child[i], list_head );

 /* Insert tree node that has finite element already done by quadtree inside */

 else if( node->elmhead != NULL )
 {
  SllAddEnd( (Sll *)list_head, sizeof(nodeadj_t), (Sll *)&le );
  le->node = node;
 }
}

#if 0
/* ======================== insert_list =================================== */

static void insert_list(

neighbour_t   *locaux,
nodeadj_t     **list_head
 )

{
 /* Look trough tree nodes where the vertex is, to get adjacent tree nodes */

 while( locaux != NULL )
  {
   /* The adjacent is only done if the tree node isn't NODE_INTERIOR */

   if( locaux->node->code != NODE_INTERIOR )
    {
     /* Get all kinds of tree nodes adjacents */

     coord_int( locaux->node , RIGHT  , list_head );
     coord_int( locaux->node , LEFT   , list_head );
     coord_int( locaux->node , TOP    , list_head );
     coord_int( locaux->node , BOTTOM , list_head );
     coord_int( locaux->node , CORNER0, list_head );
     coord_int( locaux->node , CORNER1, list_head );
     coord_int( locaux->node , CORNER2, list_head );
     coord_int( locaux->node , CORNER3, list_head );
    }

   /* Get the next tree node where the vertex is, if there's more than one */

   locaux = locaux->next;
  }
}

/* ======================== coord_int ===================================== */

static void coord_int(

tree_node_t *nodeaux,
side_t      side,
nodeadj_t   **list_head
)

{
  neighbour_t      *adj;
  nodeadj_t        *le;

  /* Get tree nodes adjacent of the nodeaux from the side desired */

  adj = find_adj_node( side, nodeaux );

  /* Look on the first tree node on the list */

  while( adj != NULL )
   {
    /* The tree node only will be considered if it's interior or vertex */

    if( (adj->node->code==NODE_INTERIOR) || (adj->node->code==NODE_VERTEX) )
     {
      /* The tree node only will be included in list if it was not included
         before to avoid repetions*/

      if( !find_list( adj->node, side, list_head ) )
       {
        /* Update the tree node on the list */

        SllAddEnd( (Sll *)list_head, sizeof(nodeadj_t), (Sll *)&le );
        le->side = side;
        le->node = adj->node;
       }
     }

    /* Get the next adjacent tree node and release the actual one */

    SllDelete( (Sll *)&adj, (Sll)adj );
   }
}
#endif

/* ====================== insert_prob_list ================================ */

static void insert_prob_list(

neighbour_t   *locaux,
nodeadj_t     **list_head
 )

{
 /* Look trough tree nodes where the vertex is, to get adjacent tree nodes */

 while( locaux != NULL )
  {
   /* Insert the adjacents tree nodes only if tree node isn't NODE_INTERIOR */

   if( locaux->node->code != NODE_INTERIOR )
    {
     /* Get all kinds of tree nodes adjacents */

     coord_prob_int( locaux->node , RIGHT  , list_head );
     coord_prob_int( locaux->node , LEFT   , list_head );
     coord_prob_int( locaux->node , TOP    , list_head );
     coord_prob_int( locaux->node , BOTTOM , list_head );
     coord_prob_int( locaux->node , CORNER0, list_head );
     coord_prob_int( locaux->node , CORNER1, list_head );
     coord_prob_int( locaux->node , CORNER2, list_head );
     coord_prob_int( locaux->node , CORNER3, list_head );
    }

   /* Insert the own tree node it is NODE_INTERIOR */

   #if 0
   else
    {
     nodeadj_t        *le;
     if( !find_list( locaux->node, TOP, list_head ) )
     {
      /* Then update the tree node adjacent list */

      SllAddEnd( (Sll *)list_head, sizeof(nodeadj_t), (Sll *)&le );
      le->side = TOP;
      le->node = locaux->node;
     }
    }
   #endif

   /* Get the next tree node where the vertex is, if there's more than one */

   locaux = locaux->next;
  }
}

/* ====================== coord_prob_int ================================== */

static void coord_prob_int(

tree_node_t *nodeaux,
side_t      side,
nodeadj_t   **list_head
)

{
 nodeadj_t        *le;
 neighbour_t      *adj;

 /* Get tree nodes adjacent of the nodeaux from the side desired */

 adj = find_adj_node( side, nodeaux );

 /* Look on the first tree node on the list */

 while( adj != NULL )
  {
   /* The routine is recursively done until a tree node interior or vertex
      be founded */

   if( (adj->node->code!=NODE_VERTEX) && (adj->node->code!=NODE_INTERIOR) )
    coord_prob_int( adj->node, side, list_head );

   else
    {
     /* The tree node only will be included in the list if it was not included
        before to avoid repetions */

     if( !find_list( adj->node, side, list_head ) )
      {
       /* Then update the tree node adjacent list */

       SllAddEnd( (Sll *)list_head, sizeof(nodeadj_t), (Sll *)&le );
       le->side = side;
       le->node = adj->node;
      }
    }

   /* Get the next adjacent tree node and release the actual one */

   SllDelete( (Sll *)&adj, (Sll)adj );
  }
}

/* ====================== locate_prob ===================================== */

static void locate_prob(

neighbour_t     *locaux,
neighbour_t     **list_head
)

{
 /* Look trough tree nodes where the vertex is, to get adjacent tree nodes */

 while( locaux != NULL )
  {
   locate_prob_list( locaux->node , RIGHT  , list_head );
   locate_prob_list( locaux->node , LEFT   , list_head );
   locate_prob_list( locaux->node , TOP    , list_head );
   locate_prob_list( locaux->node , BOTTOM , list_head );

   /* Get the next tree node where the vertex is, if there's more than one */

   locaux = locaux->next;
  }
}

/* ====================== locate_prob_list ================================ */

static void locate_prob_list(

tree_node_t   *nodeaux,
side_t        side,
neighbour_t   **list_head
)

{
 neighbour_t      *le;
 neighbour_t      *adj;

 /* Get tree nodes adjacent of the locaux from the side desired */

 adj = find_adj_node( side, nodeaux );

 /* Look on the first tree node on the list */

 while( adj != NULL )
  {
   /* Update the tree node on the list */

   SllAddEnd( (Sll *)list_head, sizeof(neighbour_t), (Sll *)&le );
   le->node = adj->node;

   /* Get the next adjacent tree node and release the actual one */

   SllDelete( (Sll *)&adj, (Sll)adj );
  }
}

/* ==================== insert_prob_tree_list ============================= */

static void insert_prob_tree_list(

tree_node_t    *node,
nodeadj_t     **list_head
)

{
 int               i;
 nodeadj_t        *le;

 /* This routine is recursively done until all quadtree tree nodes
    NODE_INTERIOR or NODE_VERTEX be included in list_head */

 /* Verify if the tree node has no child,that is, if it is a leaf */

 if( node->child != NULL )
  for( i = 0; i < 4; i++ )
   insert_prob_tree_list( node->child[i], list_head );

 /* Then if the tree node is NODE_INTERIOR or NODE_VERTEX it should be
    include in the list_head to chose after the best triangle */

 else if( node->code == NODE_INTERIOR || node->code == NODE_VERTEX )
      {
       /* Update the list_head */

       SllAddEnd( (Sll *)list_head, sizeof(nodeadj_t), (Sll *)&le );
       le->side = NO_SIDE;
       le->node = node;
      }
}

/* ====================== find_list ======================================= */

static int find_list(

tree_node_t *node,
side_t      side,
nodeadj_t   **list_head
)

{
 nodeadj_t *l;

 /* Inicializate local list_head */

 l = *list_head;

 /* Look trough the list to see if the tree node already exist inside it */

 while( l )
  {
   /* If the tree node is found, return 1 */

   if( (l->node == node)  )  return 1;

   /* Else verify with the next tree node in the list */

   l = l->next;
  }

 /* If the list_head is NULL or tree node isn't inside the list, return 0 */

 return 0;
}

/* ======================== int_quadtree ================================== */

static void int_quadtree(

list          *head,
list          *head_quad,
nodeadj_t     *list_vertex,
nodeadj_t     *list_node,
adjnode       **adj_vector,
mesh_vertex_t **no_vector,
mesh_elem_t    *list_vector
)

{
  nodeadj_t *list  = list_vertex;
  nodeadj_t *ilist = list_vertex;
  nodeadj_t *alist = list_node;

  /* Look trough the list to get the best id based in a Delaunay tecnique */

  while( list != NULL )
  {
   switch( list->node->code )
    {
     /* It has two kinds of procedures depending on if the tree node is
        NODE_VERTEX or NODE_INTERIOR */

     case NODE_VERTEX :
      {
       find_vert_id( head, head_quad, ilist, alist, list->node, adj_vector,
                     no_vector, list_vector );
       break;
      }
     case NODE_INTERIOR :
      {
       find_quad_id( head, head_quad, ilist, alist, list->node, list->side,
                     adj_vector, no_vector, list_vector );
       break;
      }
     default :
       break;
    }

   /* Get another tree node in the list */

   list = list->next;
  }
}

/* ===================== find_vert_id ===================================== */

static void find_vert_id(

list	      *head,
list          *head_quad,
nodeadj_t     *ilist,
nodeadj_t     *alist,
tree_node_t   *node,
adjnode       **adj_vector,
mesh_vertex_t **no_vector,
mesh_elem_t    *list_vector
)

{
  tree_vertex_t  *elme;
  int            intersecv, interseci, intersecb, interseca;
  int            nnosint,nnoint,nnogeomint;

  nnoint     = get_nnoint();
  nnogeomint = get_nnogeomint();

  nnosint = nnogeomint+nnoint;

  /* Look on the vertexs list in the NODE_VERTEX to see if it's the best id */

  /* Get the first vertex in the tree node */

  elme = node->vhead;

  while( elme )
  {
   /* Get the candidate id */

   candid = elme->id;

   /* If the candid is the same of the id, this vertex was already looked */

   if(( candid != id ) && ( candid < nnosint ))
    {
     /* Get the candidate coordinates */

     candidate = elme->coord;

     /* Calculate the angle */

     ang = geoangle( );

     /* If the angle is bigger than the maxang, update the id like the candid */

     if ( ( ang > maxang ) && ( test_if_isbound() ) )
      {
       /* Check intersection of the two new edges that will be formed( candid-
          firstid and candid-secondid ) against edges that already exits */

       #if QUAD_INTERS
       intersecv = intersec( adj_vector, no_vector );
       #else
       intersecv = intersec( head, head_quad, no_vector );
       #endif

       /* Check intersection of the two new edges that will be formed( candid-
          firstid and candid-secondid ) against quadtree boundary */

       #ifdef QUAD_INQUAD
       intersecb = interqtree( ilist, no_vector );
       #else
       intersecb = 1;
       #endif

       /* Check intersection of the two new edges that will be formed( candid-
          firstid and candid-secondid ) against all quadtree boundary */

       #if QUAD_ALLQUAD
       interseca = interqtree( alist, no_vector );
       #else
       interseca = 1;
       #endif

       /* Check if id choosen is inside of any finite element polygon already
          maked */

       #ifdef QUAD_INSIDE
       interseci = inside( no_vector, list_vector );
       #else
       interseci = 1;
       #endif

       /* If there's no intersection, update the candid,iddate and maxang */

       if( intersecv && interseci && intersecb && interseca )
        {
         id = candid;
         iddate = candidate;
         maxang = ang;
        }
      }
    }

   /* Get a new vertex in the tree node */

   elme = elme->next;
  }
}

/* ======================== find_quad_id ================================== */

static void find_quad_id
(
 list        *head,
 list        *head_quad,
 nodeadj_t   *ilist,
 nodeadj_t   *alist,
 tree_node_t *node,
 side_t      side,
 adjnode     **adj_vector,
 mesh_vertex_t **no_vector,
 mesh_elem_t    *list_vector
 )

{
 mesh_elem_t       *elme;
 int               j;
 int               intersecq, interseci, intersecb, interseca;
 int               filtroquad;
 int               nnosint,nnoint,nnogeomint;

  nnoint     = get_nnoint();
  nnogeomint = get_nnogeomint();

  nnosint = nnogeomint+nnoint;

  /* Look if the interior tree node has any finite element */

  if( node->elmhead != NULL )
      {
        /* Get the first finite element */

        elme = node->elmhead;
        while( elme )
        {
          /* All vertexs in the tree node finite element should be tested */

          for( j = 0; j < elme->n; j++ )
            {
             /* Inicializate the candid */

             candid = elme->conect[j]->nno;

             /* If candid is the same of the id, means that this vertex was
                already tested */

             if( ( candid != id ) && (candid<nnosint) )
              {
               /* Update candidate coordinates */

               candidate.ex = elme->conect[j]->x;
               candidate.ey = elme->conect[j]->y;

               /* Test quadtree filter to avoid intersection tests */

               #if QUAD_NORMAL
               filtroquad = 1;
               #else
               filtroquad = filtro_quad( node, side );
               #endif

               if( filtroquad )
                {
                 /* Compute ang for the candidate */

                 ang = geoangle( );

                 /* If ang is bigger than the maxang, update it and the id */

                 if( ang > maxang )
                  {
                   /* Test with there's no intersection against edges that
                      already exists */

                   #if QUAD_INTERS
                   intersecq = intersec( adj_vector, no_vector );
		   #else
                   intersecq = intersec( head, head_quad, no_vector );
                   #endif

                   /* Test with there's no intersection against quadtree boundary */

                   #ifdef QUAD_INQUAD
                   intersecb = interqtree( ilist, no_vector );
                   #else
                   intersecb = 1;
                   #endif

                   /* Test with there's no intersection against all quadtree boundary */

                   #if QUAD_ALLQUAD
                   interseca = interqtree( alist, no_vector );
                   #else
                   interseca = 1;
                   #endif

                   /* Check if id choosen is inside of any finite element
                      polygon already maked */

                   #ifdef QUAD_INSIDE
                   interseci = inside( no_vector, list_vector );
		   #else
		   interseci = 1;
                   #endif

                   /* If there's no intersection, upadte id,iddate and maxang */

                   if( intersecq && interseci && intersecb && interseca )
                    {
                     id = candid;
                     iddate = candidate;
                     maxang = ang;
                    }
                  }
                }
              }
            }

          /* Get a new finite element in the tree node */

          elme = elme->next;
        }
      }
 }

#if !QUAD_NORMAL
/* ===================== filtro_quad ====================================== */

static int filtro_quad(

tree_node_t *node,
side_t      side
)

{
 new_coord   base;
 new_coord   data;
 int         size;
 int         co;
 int         c1;
 int         c2;

 /* Compute the tree node size.Tree node here is the tree node in the list,
    that is, adjacent of tree node where first or second vertex is ,  about
    any side adjacent */

 size = pow_2( 14 - node->depth );

 /* Inicializate tree node base coordinate and filter code.The tree node base
    is tree node coord, that is the left inferior corner of the tree node */

 base = node->coord;
 co = 0;

 /* Do the filter for witch case of adjacent size */

 switch( side )
  {
   case RIGHT :
    {
     /* If the candidate has the same x coordinate of the left tree node side,
        then could be choosen */

     if( candidate.ex == base.ex )  co = 1;
     break;
    }
   case LEFT :
    {
     /* If the candidate has same x coordinate of the right tree node side,
        then could be choosen */

     if( candidate.ex == (base.ex + size) )  co = 1;
     break;
    }
   case TOP :
    {
     /* If the candidate has same y coordinate of the bottom tree node side,
        then could be choosen */

     if( candidate.ey == base.ey )  co = 1;
     break;
    }
   case BOTTOM :
    {
     /* If the candidate has the same y coordinate of the top tree node side,
        then could be choosen */

     if( candidate.ey == (base.ey + size) )  co = 1;
     break;
    }
   case CORNER0 :
    {
     /* If the point is really in the corner, code is 1 without test */

     if( (candidate.ex == (base.ex+size)) && (candidate.ey == (base.ey+size)) )
       return 1;

     /* If the adjacent is any corner, the candidate only could be choosen if
        it has the same coordinates of the data coordinate shown below for
        two possibilities */

     data.ex = base.ex + size;
     data.ey = base.ey + ( size/2 );
     c1 = filtro_quad_test( node, RIGHT, data, candidate.ex, candidate.ey,
              base.ex + size, base.ey, base.ey + size/2, base.ey + size );

     data.ex = base.ex + ( size/2 );
     data.ey = base.ey + size;
     c2 = filtro_quad_test( node, TOP, data, candidate.ey, candidate.ex,
               base.ey + size, base.ex, base.ex + size/2, base.ex + size );

     if( c1 == 1 || c2 == 1 )  co = 1;
     break;
    }
   case CORNER1 :
    {
     /* If the point is really in the corner, code is 1 without test */

     if( (candidate.ex == base.ex) && (candidate.ey == (base.ey+size)) )
       return 1;

     /* If the adjacent is any corner, the candidate only could be choosen if
        it has the same coordinates of the data coordinate shown below for
        two possibilities */

     data.ex = base.ex - ( size/2 );
     data.ey = base.ey + ( size/2 );
     c1 = filtro_quad_test( node, LEFT, data, candidate.ex, candidate.ey,
               base.ex, base.ey, base.ey + size/2, base.ey + size );

     data.ex = base.ex;
     data.ey = base.ey + size;
     c2 = filtro_quad_test( node, TOP, data, candidate.ey, candidate.ex,
               base.ey + size, base.ex + size, base.ex + size/2, base.ex );

     if( c1 == 1 || c2 == 1 )  co = 1;
     break;
    }
   case CORNER2 :
    {
     /* If the point is really in the corner, code is 1 without test */

     if( (candidate.ex == (base.ex+size)) && (candidate.ey == base.ey) )
      return 1;

     /* If the adjacent is any corner, the candidate only could be choosen if
        it has the same coordinates of the data coordinate shown below for
        two possibilities */

     data.ex = base.ex + size;
     data.ey = base.ey;
     c1 = filtro_quad_test( node, RIGHT, data, candidate.ex, candidate.ey,
               base.ex + size, base.ey + size, base.ey + size/2, base.ey );

     data.ex = base.ex + ( size/2 );
     data.ey = base.ey - ( size/2 );
     c2 = filtro_quad_test( node, BOTTOM, data, candidate.ey, candidate.ex,
               base.ey, base.ex, base.ex + size/2, base.ex + size );

     if( c1 == 1 || c2 == 1 )  co = 1;
     break;
    }
   case CORNER3 :
    {
     /* If the point is really in the corner, code is 1 without test */

     if( (candidate.ex == base.ex) && (candidate.ey == base.ey) )
       return 1;

     /* If the adjacent is any corner, the candidate only could be choosen if
        it has the same coordinates of the data coordinate shown below for
        two possibilities */

     data.ex = base.ex - ( size/2 );
     data.ey = base.ey;
     c1 = filtro_quad_test( node, LEFT, data, candidate.ex, candidate.ey,
               base.ex, base.ey + size, base.ey + size/2, base.ey );

     data.ex = base.ex;
     data.ey = base.ey - ( size/2 );
     c2 = filtro_quad_test( node, BOTTOM, data, candidate.ey, candidate.ex,
               base.ey, base.ex + size, base.ex + size/2, base.ex  );

     if( c1 == 1 || c2 == 1 )  co = 1;
     break;
    }
   default :
     /* The tree node has NO_SIDE definition */

     co = 1;
     break;
  }

  /* Return filtro_quad test */

  return( co );
}

/* ================== filtro_quad_test ==================================== */

static int filtro_quad_test
(

tree_node_t *node,
side_t      side,
new_coord   d,
int         coorda,
int         coordb,
int         li,
int         di,
int         dm,
int         df
)

{
 neighbour_t *ad;
 new_coord   cn;

 /* Filter tests relative to corners are little complicated.You have tree
    vertexs possibilities: two vertexs in each corner, and a a third one if
    tree node from the corner is a depth more than tree node who get this
    tree node corner.More, the tree vertexs only could be in the same x or y
    coordinate of the tree node who get the tree node corner, depending of
    the corner adjacent*/

 /* Get adjacent tree nodes list of the tree nodes from the corner */

 ad = find_adj_node( side, node );

 /* If adjacent list has only one tree node, it means that candidate could be
    only the two vertexs in each corner from the tree node corner, that  will
    be the vertex with the same coordinate of the tree node that gave    this
    tree node corner depending wich corner is */

 if( ad != NULL )
 {

  if( ad->next == NULL )
   {
    if( (ad->node->code == NODE_BOUND) && (li == coorda) )
     {
      /* Release adjacent tree node list */

      SllDelAll( (Sll *)&ad );

      /* Return filter code */

      return  1;
     }
   }

  /* Else candidate could have the tree possibilities told before */

  else
   {
    /* The first is to be in the more far corner of the tree node corner and
       has the same coordinate like the two possibilities shown before */

    if( (di == coordb) && (li == coorda) )
     if( (ad->node->code == NODE_BOUND) &&
         (ad->next->node->code == NODE_BOUND) )
      {
       /* Release adjacent tree node structure */

       SllDelAll( (Sll *)&ad );

       /* Return filter code */

       return 1;
      }

    /* The second is to be in the middle of the tree node corner and also
       has the same coordinate like the two possibilities shown before */

    if( (dm == coordb) && (li == coorda) )
     {
      while( ad != NULL )
       {
        cn = ad->node->coord;
        if( (d.ex == cn.ex) && (d.ey == cn.ey) )
         if( ad->node->code == NODE_BOUND )
          {
           /* Release adjacent tree node structure */

           SllDelAll( (Sll *)&ad );

           /* Return filter code */

           return 1;
          }

        /* If not return, get the next tree node and release the actual one */

        SllDelete( (Sll *)&ad, (Sll)ad );
       }
     }

    /* And the third is to be in the other near corner and also has the same
       the same coordinate like the two possibilities shown before too */

    if( (df == coordb) && (li == coorda) )
     {
      /* Release adjacent tree node structure */

      SllDelAll( (Sll *)&ad );

      /* Return filter code */

      return 1;
     }
   }

 }

 /* If candidate isn't in neither possibilities,then it could'nt be choosen */

 /* Release adjacent tree node structure */

 SllDelAll( (Sll *)&ad );

 /* Return filter code */

 return 0;
}
#endif

/*======================= intersec  ==========================================*/

#if QUAD_INTERS
static int intersec(

adjnode       **adj_vector,
mesh_vertex_t **no_vector
)

{
 int  cl;
 int  cf;
 int  cs;
 int  ce;
 int  cquad;
 list *listv;
 list *li;

 /* This routine test intersection against the edges in a posssible   triangle
    formed by firstid,secondid and candid,and edges already formed by  vertexs
    in tree nodes adjacent to the base edge. Note that the test is   not  done
    against the base edge,because it already exist and  then  can't  intersect
    or not was choosen in a anterior step.It returns  a code 0  if  intersects
    and 1 if not */

   /* Inicializate edges intersect auxiliary list */

   listv = NULL;

   /* Test intersection of new edges against edges formed by vertexs in the
      tree node list adjacent to first and second vertex on the base edge */

   cl = intersec_list( listvertex, adj_vector, no_vector, &listv );
   if( cl == 0 )
    {
     SllDelAll( (Sll *)&listv );
     return cl;
    }

   /* Test intersection of new edges against the edges formed by firstid */

   cf = intersec_vert( curfirstid, adj_vector, no_vector, &listv );
   if( cf == 0 )
    {
     SllDelAll( (Sll *)&listv );
     return cf;
    }

   /* Test intersection of new edges against the edges formed by secondid */

   cs = intersec_vert( cursecondid, adj_vector, no_vector, &listv );
   if( cs == 0 )
    {
     SllDelAll( (Sll *)&listv );
     return cs;
    }

   /* Test intersection of the new edges against edges formed by vertexs in  a
      quadtree vertex list of edges, that don't belongs to the  same tree node
      interior */

   if( listquad != NULL )
    {
     li = listquad;
     while( li )
      {
       ce = li->vi;
       cquad = intersec_vert( ce, adj_vector, no_vector, &listv );
       if( cquad == 0 )
        {
         SllDelAll( (Sll *)&listv );
         return cquad;
        }

       ce = li->vj;
       cquad = intersec_vert( ce, adj_vector, no_vector, &listv );
       if( cquad == 0 )
        {
         SllDelAll( (Sll *)&listv );
         return cquad;
        }

       li = li->next;
      }
    }

 /* No intersection was found, then returns code 1 */

 SllDelAll( (Sll *)&listv );
 return 1;
}
#else
static int intersec(

list          *head,
list          *head_quad,
mesh_vertex_t **no_vector
)

{
 int cff, cfs;

 while( head )
 {
  cff = intersec_cross( curfirstid, candid, head->vi, head->vj, no_vector );
  if( cff  == 0 ) return cff;
  cfs = intersec_cross( candid, cursecondid, head->vi, head->vj, no_vector );
  if( cfs == 0 ) return cfs;
  head = head->next;
 }
 while( head_quad )
 {
  cff = intersec_cross( firstid, candid, head_quad->vi, head_quad->vj, no_vector );
  if( cff  == 0 ) return cff;
  cfs = intersec_cross( candid, secondid, head_quad->vi, head_quad->vj, no_vector );
  if( cfs == 0 ) return cfs;
  head_quad = head_quad->next;
 }

 return 1;
}
#endif

/* ===================== interqtree ======================================= */

#if defined(QUAD_INQUAD) || QUAD_ALLQUAD
static int interqtree(

nodeadj_t      *qlist,
mesh_vertex_t **no_vector
)

{
 int j, vi, vj, cff, cfs;
 mesh_elem_t *elme;

 /* Test intersection against finite elements in the interior of quadtree */

 while( qlist != NULL )
 {
  if( qlist->node->elmhead != NULL )
  {
   elme = qlist->node->elmhead;
   while( elme )
   {
    for( j = 0; j < elme->n; j++ ) /* intersection against edge */
    {
     vi = elme->conect[j]->nno;
     vj = elme->conect[(j+1)%elme->n]->nno;

     cff = intersec_cross( curfirstid, candid, vi, vj, no_vector );
     if( cff  == 0 ) return cff;
     cfs = intersec_cross( candid, cursecondid, vi, vj, no_vector );
     if( cfs == 0 ) return cfs;
    }
    for( j = 0; j < 2; j++ ) /* intersection against diagonals */
    {
     if( elme->n == 4 ) /* diagonal only of quadrilaterals */
     {
      vi = elme->conect[j]->nno;
      vj = elme->conect[(j+2)]->nno;

      cff = intersec_cross( curfirstid, candid, vi, vj, no_vector );
      if( cff  == 0 ) return cff;
      cfs = intersec_cross( candid, cursecondid, vi, vj, no_vector );
      if( cfs == 0 ) return cfs;
     }
    }
    elme = elme->next;
   }
  }
  qlist = qlist->next;
 }

 return 1;
}
#endif

#if QUAD_INTERS
/* ===================== intersect_list =================================== */

static int intersec_list(

nodeadj_t     *heq,
adjnode       **adj_vector,
mesh_vertex_t **no_vector,
list          **listv
)

{
 int ca;
 int ce;
 int cq;
 int j;
 int auxid;
 int size;
 new_coord base;
 new_coord aux;
 tree_vertex_t *elme;
 mesh_elem_t   *elma;
 int           p[4];

 /* Get a tree node from the list */

 while( heq )
  {
   switch( heq->node->code )
    {
     /* Test intersection against edges formed by vertexs only in a tree node
        that is NODE_VERTEX */

     case NODE_VERTEX :
      {
       /* Get the first vertex in the tree node */

       elme = heq->node->vhead;
       while( elme )
        {
         /*Get the vertex id*/

         ce = elme->id;

         /* If this vertex is the candid, it should not be tested intersection
            against edges formed by it,because like the two new edges contains
            candid like one or its vertex, then  it isn't   possible intersect
            edges formed by candid. Another vertexs should be tested and if  a
            intersection is found, return with code 0. */

         if( ce != candid )
          {
           ca = intersec_vert( ce, adj_vector, no_vector, listv );
           if( ca == 0 )  return ca;
          }

         /* Get the next vertex in the tree node */

         elme = elme->next;
        }
       break;
      }

     /* Test intersection against edges formed by vertexs only in a tree node
        that is NODE_INTERIOR */

     case NODE_INTERIOR :
      {
       /* Calculate the tree node size and the left corner tree node */

       size = pow_2( 14 - heq->node->depth );
       base = heq->node->coord;

       /* Get the first finite element in the tree node */

       if( heq->node->elmhead != NULL )
        {
         elma = heq->node->elmhead;
         while( elma )
         {
          /* Look trough all vertexs in the finite elements to find the four
             vertexs in the corners of the tree node */

          for( j = 0; j < elma->n; j++ )
           {
            /*Get the vertex id and coordinates*/

            auxid  = elma->conect[j]->nno;
            aux.ex = elma->conect[j]->x;
            aux.ey = elma->conect[j]->y;

            /* Get the four corners of the tree node */

            if( aux.ex==base.ex&&aux.ey==base.ey )               p[0] = auxid;
            if( aux.ex==(base.ex+size)&&aux.ey==base.ey )        p[1] = auxid;
            if( aux.ex==base.ex&&aux.ey==(base.ey+size) )        p[2] = auxid;
            if( aux.ex==(base.ex+size)&&aux.ey==(base.ey+size) ) p[3] = auxid;
           }

          /* Get a new finite element in the tree node */

          elma = elma->next;
         }
        }

       /* Verify intersection against the diagonals of the tree node */

       cq = intersec_quad( p[0], p[3], no_vector );
       if( cq == 0 ) return cq;

       cq = intersec_quad( p[1], p[2], no_vector );
       if( cq == 0 ) return cq;

       break;
      }

     /* Another types of tree nodes don't have geometric vertexs inside */

     default :
       break;
    }

   /* Get the next tree node in the list */

   heq = heq->next;
  }

 /* No intersection was found, then returns code 1 */

 return 1;
}

/* ===================== intersec_vert =================================== */

static int intersec_vert(

int           auxid,
adjnode       **adj_vector,
mesh_vertex_t **no_vector,
list          **listv
)

{
 adjnode       *t;
 list          *inew;
 int           cff;
 int           codef;
 int           cfs;
 int           codes;

   /* Look in the adjacent edges vertex structure to get the first
      edge formed by the input vertex in the tree node NODE_VERTEX */

   t = adj_vector[auxid];

   while( t )
    {
     /* Test if this edge intersects the edges formed by firstid and  candid
        or the edge formed by secondid and candid.Note that isn't     tested
        against the base edge because it already exists and then don't cross
        this edge.The test will be done using geometrical computacional tec-
        niques with cross products */

     /* Inicilaizate edges intersection auxiliary list code */

     codef = codes = 0;

     /* Verify if the intersection against this edge was not tested before.
        If it was not, test intersection against the two new edges */

     if( !find_intersec( t->p->vi, t->p->vj, listv ) )
      {
       /* First test the edge against edge formed by firstid and candid.If the
          edge has one of their vertexs equal of one of the vertexs of the edge
          formed by firstid and candid then it not should be tested    because
          they won't intersect */

       if( ( (t->p->vi)!=curfirstid && (t->p->vi)!=candid ) &&
           ( (t->p->vj)!=curfirstid && (t->p->vj)!=candid )  )
        {
         codef = 1;
         cff = intersec_cross( curfirstid, candid,t->p->vi,t->p->vj,no_vector );
         if( cff == 0 )  return cff;
        }

       /* Second test the edge against edge formed by firstid and candid.If the
          edge has one of their vertexs equal of one of the vertexs of the edge
          formed by firstid and candid then it  not should be tested    because
          they won't intersect */

       if( ( (t->p->vi)!=cursecondid && (t->p->vi)!=candid ) &&
           ( (t->p->vj)!=cursecondid && (t->p->vj)!=candid )  )
        {
         codes = 1;
         cfs = intersec_cross( candid,cursecondid,t->p->vi,t->p->vj,no_vector );
         if( cfs == 0 )  return cfs;
        }

       /* Update edges list intersection if the edge was used */

       if( (codef == 1) || (codes == 1) )
        {
         SllAddEnd( (Sll *)listv, sizeof(list), (Sll *)&inew );
         inew->vi = t->p->vi;
         inew->vj = t->p->vj;
        }

      }

     /* Get the next edge structure */

     t = t->next;
    }

   /* At this point no intersection was found, then returns code 1 */

   return 1;
}

/* ===================== intersec_quad ==================================== */

static int intersec_quad(

int           df,
int           ds,
mesh_vertex_t **no_vector
)

{
 int           cff;
 int           cfs;

     /* First test the edge against edge formed by firstid and candid.If the
        edge has one of their vertexs equal of one of the vertexs of the edge
        formed by firstid and candid then it not should be tested    because
        they won't intersect */

     if( ( df!=curfirstid && df!=candid ) &&
         ( ds!=curfirstid && ds!=candid )  )
      {
       cff = intersec_cross( curfirstid, candid, df, ds, no_vector );
       if( cff == 0 )  return cff;
      }

     /* Second test the edge against edge formed by firstid and candid.If the
        edge has one of their vertexs equal of one of the vertexs of the edge
        formed by firstid and candid then it  not should be tested    because
        they won't intersect */

     if( ( df!=cursecondid && df!=candid ) &&
         ( ds!=cursecondid && ds!=candid )  )
      {
       cfs = intersec_cross( candid, cursecondid, df, ds, no_vector );
       if( cfs == 0 )  return cfs;
      }


     /* At this point no intersection was found, then returns code 1 */

     return 1;
}

/* ====================== find_intersec =================================== */

static int find_intersec(

int         vi,
int         vj,
list        **list_head
)

{
 list *l;

 /* Inicializate local list_head */

 l = *list_head;

 /* Look trough the list to see if the  edge already exist inside it */

 while( l )
  {
   /* If the edge is found, return 1 */

   if( (l->vi==vi || l->vi==vj) && (l->vj==vi || l->vj==vj) )
    return 1;

   /* Else verify with the next tree node in the list */

   l = l->next;
  }

 /* If list_head is NULL or the tree node isn't inside the list, return 0 */

 return 0;
}
#endif

/* ===================== intersec_cross =================================== */

static int intersec_cross(

int           pa,
int           qa,
int           vi,
int           vj,
mesh_vertex_t **no_vector

)

{
 float c0,c1,c2,c3;

 /* Compute first two cross products: the edge qa-pa and edge vi-pa, and
    edge qa-pa and edge vj-pa */

 c0 = intersec_cross_value( pa, qa, vi, no_vector );
 c1 = intersec_cross_value( pa, qa, vj, no_vector );
 if( (c0*c1) >= 0.0 )  return 1;

 /* Compute next two cross products:  the edge  vj-vi and  edge  pa-vi,
    and edge vj-vi and qa-vi */

 c2 = intersec_cross_value( vi, vj, pa, no_vector );
 c3 = intersec_cross_value( vi, vj, qa, no_vector );
 if( (c2*c3) >= 0.0 )  return 1;
 else                  return 0;
}

/* ===================== intersec_cross_value ============================= */

static float intersec_cross_value(

int           a,
int           b,
int           c,
mesh_vertex_t **no_vector
)

{
 vector v1;
 vector v2;

 /* Make the two vectors to the cross product */

 v1.x = (*no_vector)[b].x - (*no_vector)[a].x;
 v1.y = (*no_vector)[b].y - (*no_vector)[a].y;
 v2.x = (*no_vector)[c].x - (*no_vector)[a].x;
 v2.y = (*no_vector)[c].y - (*no_vector)[a].y;

 /* Compute a cross product of the two vectors */

 return cross_product( &v1, &v2 );
}

#ifdef QUAD_INSIDE
/*======================= inside =============================================*/

static int inside(

mesh_vertex_t **no_vector,
mesh_elem_t    *list_vector
)

{
 int          i, j, n;
 float        x, y;
 vector       p[4];
 mesh_elem_t *elem;

 /*Get candidate coordinates*/

 x = (*no_vector)[candid].x;
 y = (*no_vector)[candid].y;

 /*Look trough all finite element polygon structure*/

 elem = list_vector;
 while( elem )
 {
  /*Get number of vertexs*/

  n = elem->n;

  /*Get vertexs vector*/

  for( i = 0, j = n-1; i < n; i++, j-- )
  {
   p[i].x = elem->conect[j]->x;
   p[i].y = elem->conect[j]->y;
  }

  /*Verify if candidate is inside of polygon*/

  if( inside_polygon( n, p, x, y ) )
   return 0;

  /*Get next element*/

  elem = elem->next;
 }

 /*If we get here, then candidate is outside of all finite element polygons*/

 return 1;
}
#endif

/* ======================= inside_polygon =================================== */

/*
** Pick polygon:
**  Verify if the point (x,y) is inside the polygon. It returns:
**   1 - inside
**   0 - outside
*/
int inside_polygon( int n, vector *p, float x, float y )
{
 int i;
 int ni = 0; 		/* number of intersections */
 int fst = n - 1;	/* begin with the last node */
 vector p1, p2; 	/* edges points */
 float xc;

 for (i = 0; i < n; i++)
 {
  p1 = p[i];
  p2 = p[fst];
  if (!(p1.y == p2.y) &&	  	/* descarta horizontais */
      !((p1.y > y) && (p2.y > y)) &&	/* descarta retas acima */
      !((p1.y < y) && (p2.y < y)) && 	/* descarta retas abaixo */
      !((p1.x < x) && (p2.x < x)))  	/* descarta retas esquerda */
  {
   if (p1.y == y)			/* primeiro ponto na mesma cota */
   {
    if ((p1.x > x) && (p2.y > y))
     ni++;				/* a direita e acima do ponto */
   }
   else
   {
    if (p2.y == y)		/* segundo ponto na mesma cota */
    {
     if ((p2.x > x) && (p1.y > y))
      ni++;			/* a direita e acima do ponto */
    }
    else
    {
     if ((p1.x > x) && (p2.x > x))
      ni++; 			/* inteiramente a direita */
     else
     {                   	/* verifica ponto de intersecao */
      float dx = p1.x - p2.x;
      xc = p1.x;
      if ( dx != 0.0 )
	xc += ( y - p1.y ) * dx / ( p1.y - p2.y );
      if (xc > x)
	ni++;
     }
    }
   }
  }
  fst = i;  			/* ultimo ponto para proxima aresta */
 }
 return ( ni % 2 );
}

/* ========================= boundr_polygon =============================== */

/*
** Boundary polygon:
**  Verify if the point (x,y) is in the boundary of the polygon. It returns:
**   1 - in  the boundary
**   0 - out the boundary
*/
static int boundr_polygon( int n, vector *p, float x, float y )
{
 int i;
 float nx, ny;

 for( i = 0; i < n; i++ )
 {
  nx = ( p[(i+1)%n].x + p[i].x ) / 2.0;
  ny = ( p[(i+1)%n].y + p[i].y ) / 2.0;
  if( ABS(nx-x) < BTOL && ABS(ny-y) < BTOL )
   return 1;
 }

 return 0;
}

/* ====================== geoangle ======================================== */

static double geoangle ( void )
{
 double ikx, iky, jkx, jky;
 double angle;
 double xprod;

 /* Look if the side is on the left of the contour; if it's the angle is 0,that
    is, this vertex could not be used.This test is done using cross product */

 xprod = - (( cursecond.ex - curfirst.ex ) * ( candidate.ey - curfirst.ey ))
         + (( cursecond.ey - curfirst.ey ) * ( candidate.ex - curfirst.ex )) ;

 /* If the vertex is inside the contour, compute the angle */

 if( xprod >= 0.0 )
  {

   ikx = curfirst.ex - candidate.ex;
   iky = curfirst.ey - candidate.ey;
   jkx = cursecond.ex - candidate.ex;
   jky = cursecond.ey - candidate.ey;

   angle = ((ikx * jkx) + (iky * jky)) /
		 (sqrt(ikx*ikx + iky*iky) * sqrt(jkx*jkx + jky*jky));

   if( angle > 1.0 )  angle = 1.0;
#if 0
   if( ( ((angle+ATOL) > 1.0) && ((angle-ATOL) < 1.0) ) ||
       ( ((angle+ATOL) > -1.0) && ((angle-ATOL) < -1.0) ) )
     return 0.0;
#endif
   angle = acos (angle*API/180.0);

   /* Return the angle calculated */

   return( angle );
  }

 else
   return 0.0;
}

/* ======================== getsvalue ===================================== */

static double getsvalue ( new_coord q_vertex[4] )
{
 int i, indc, indf, inds;
 double bm, difmax, svalue;
 double ang[4];

 /* Get the polygon aspect form */

 bm =0.0;

 for ( i = 0; i < 4 ; i++ ) {
  switch ( i )
  {
   case 0 :
    indc = 0;
    indf = 1;
    inds = 3;
    ang[i] = calcangle ( indc, indf, inds, q_vertex );
   break;

   case 1 :
    indc = 1;
    indf = 2;
    inds = 0;
    ang[i] = calcangle ( indc, indf, inds, q_vertex );
   break;
   case 2 :
    indc = 2;
    indf = 3;
    inds = 1;
    ang[i] = calcangle ( indc, indf, inds, q_vertex );
   break;
   case 3 :
    indc = 3;
    indf = 0;
    inds = 2;
    ang[i] = calcangle ( indc, indf, inds, q_vertex );
   break;
  }
  difmax = ABS(API/2 - ang[i]);
  if ( difmax > bm ) bm = difmax;
 }

 svalue = 1.0 + (( 2.0/API )*bm);
 return (svalue);
}

/* ======================== calcangle ===================================== */

static double calcangle ( int c, int f, int s, new_coord q_vertex[4] )
{
 double ikx, iky, jkx, jky;
 double angle;


 /* compute angle */

 ikx = q_vertex[f].ex - q_vertex[c].ex;
 iky = q_vertex[f].ey - q_vertex[c].ey;
 jkx = q_vertex[s].ex - q_vertex[c].ex;
 jky = q_vertex[s].ey - q_vertex[c].ey;

 angle = ((ikx * jkx) + (iky * jky)) /
 	  (sqrt(ikx*ikx + iky*iky) * sqrt(jkx*jkx + jky*jky));

 if( angle > 1.0 )  angle = 1.0;

 if( ( ((angle+ATOL) > 1.0) && ((angle-ATOL) < 1.0) ) ||
     ( ((angle+ATOL) > -1.0) && ((angle-ATOL) < -1.0) ) )
  return 0.0;

 #if 0
 angle = acos (angle*API/180.0);
 #endif
 angle = acos (angle);

 /* Return the angle calculated */

 return( angle );
}

/* ======================== quad_edge ===================================== */

static void quad_edge(

 tree_node_t *node,        /* quadtree tree node root                (in) */
 new_coord   pac,          /* first edge vertex integer coordinates  (in) */
 new_coord   qac,          /* second edge vertex integer coordinates (in) */
 int         pa,           /* first edge vertex id                   (in) */
 int         qa,           /* second edge vertex id                  (in) */
 adjnode     **adj_vector, /* adjacent edges vertex structure    (in/out) */
 list        **head        /* cedge boundary structure           (in/out) */
)

{
 neighbour_t *loc1;   /* first vertex edge tree node list */
 neighbour_t *loc2;   /* second vertex edge tree node list */
 neighbour_t *loca;   /* auxiliary vertex edge tree node list pointer */
 neighbour_t *locb;   /* auxiliary vertex edge tree node list pointer */
 list        *l;      /* new cedge boundary element */

 /* Get the tree node list where the first vertex is */

 loc1 = locate_quad( &pac, node );

 /* Get the tree node list where the second vertex is */

 loc2 = locate_quad( &qac, node );

 /* Look trough the two lists to see if they have a common tree node that is
    a NODE_INTERIOR. If it happens, the edge is a internal quadtree edge and
    should not be updated in the cedge structure by update_edge function,that
    is, return to msh_update_cedge function and not update_edge or SllAddEnd */

 locb = loc2;
 while( locb != NULL )
  {
   loca = loc1;
   while( loca != NULL )
    {
     if( (loca->node==locb->node) && (loca->node->code==NODE_INTERIOR) )
      {
       SllDelAll( (Sll *)&loc1 );
       SllDelAll( (Sll *)&loc2 );
       return;
      }
     loca = loca->next;
    }
   locb = locb->next;
  }

 /* If the edge isn't a internal quadtree edge, it should be updated in cedge
    boundary structure and adjacent edges vertexs structure */

 update_edge( pa, qa, adj_vector, head, END );

 /* It also should be updated an auxiliary list with vertexs that are  vertexs
    in the internal quadtree, but make edges that aren't internal quadtree ed-
    ges */

 SllAddEnd( (Sll *)&listquad, sizeof(list), (Sll *)&l );
 l->vi = pa;
 l->vj = qa;

 /* Release locate_quad lists */

 SllDelAll( (Sll *)&loc1 );
 SllDelAll( (Sll *)&loc2 );
}

/* ======================== update_edge ==================================== */

#if QUAD_INTERS
static void update_edge(

int   pa,              /* first edge vertex id                (in) */
int   qa,              /* second edge vertex id               (in) */
adjnode **adj_vector,  /* adjacent edges vertex structure (in/out) */
list  **head,          /* cedge boundary structure        (in/out) */
int   pos

)

{
 list  *le;         /* auxiliary cedge boundary structure pointer */
 list  *elem;       /* auxiliary cedge boundary structure element */
 int   ref;
 adjnode *ed;


 /* In this update_edge, if the edge is a new one, it should be introduced in
    the cedge structure and if already exists, should be deleted of the cedge
    list */

 /* Get the lower id of the edge. It will be used to look trough the  adjacent
    edges vertex structure to see if this vertex already make this edge,avoiding
    searchs in all cedge structure and otimizing the search to see if  the  edge
    already exists */

 if( pa < qa )  ref = pa;
 else           ref = qa;

 /* Get the head of edge list of this vertex */

 ed = adj_vector[ref];

 /* Look in this list to see if this vertex already make the edge to be updated.
    If the edge already exists, get the pointer elem to delete directly in the
    cedge structure without search. If not, introduce the edge in cedge structu-
    re and upadte adj_vector */

 if( ed == NULL )  elem = NULL;
 else {
    while( ed != NULL )
    {
     if( ( ed->p->vi == pa || ed->p->vi == qa ) &&
         ( ed->p->vj == pa || ed->p->vj == qa )  )
      {
       elem = ed->p;
       break;
      }
     else
       elem = NULL;

     ed = ed->next;
    }
   }


 /* If edge already exists, delete directly in the cedge structure, else intro-
    duce in the cedge structure. In both cases update adjacent vector */

 if( elem != NULL )
  {
   /* Update cedge structure */

   SllDelete( (Sll *)head, (Sll)elem );

   /* Update adjacent edges vertex structure for the two vertexs. Here   has  a
      point: the  search is done only for the lower vertex but I should  update
      the adjacent structure for the two vertexs becuse the two will be used in
      the intersection test */

   del_adjnode( pa, qa, &adj_vector[pa] );
   del_adjnode( pa, qa, &adj_vector[qa] );
 }
 else
  {
   /* Update cedge structure. The edge not exist yet and should be introduced */

   if (pos)
   {
    SllAddEnd( (Sll *)head, sizeof(list),(Sll *)&le);
    le->type = QUAD_INTERNAL;
    le->vi = pa;
    le->vj = qa;

    /* Update the adjacent vector for the two vertexs */

    add_adjnode( le, &adj_vector[pa] );
    add_adjnode( le, &adj_vector[qa] );
   }
   else
   {
    SllAddTop( (Sll *)head, sizeof(list));
    (*head)->type = QUAD_INTERNAL;
    (*head)->vi = pa;
    (*head)->vj = qa;

    /* Update the adjacent vector for the two vertexs */

    add_adjnode( head, &adj_vector[pa] );
    add_adjnode( head, &adj_vector[qa] );
   }
  }
}
#else
static void update_edge(

int   pa,             /* first edge vertex id                (in) */
int   qa,             /* second edge vertex id               (in) */
adjnode **adj_vector, /* adjacent edges vertex structure (in/out) */
list  **head,         /* cedge boundary structure        (in/out) */
int   pos
)

{
 list  *le = NULL;         /* auxiliary cedge boundary structure pointer */
 list  *elem = NULL;       /* auxiliary cedge boundary structure element */

 le = *head;

 while( le )
 {
  if( ( le->vi == pa || le->vi == qa ) &&
      ( le->vj == pa || le->vj == qa )  )
  {
   elem = le;
   break;
  }
  else
   elem = NULL;

  le = le->next;
 }

 /* If edge already exists, delete directly in cedge structure, else intro-
    duce in the cedge structure. In both cases update adjacent vector */

 if( elem != NULL )
  {
   /* Update cedge structure */

   SllDelete( (Sll *)head, (Sll)elem );

  }
 else
  {
   /* Update cedge structure.The edge not exist yet and should be introduced */

  if (pos)
  {
   SllAddEnd( (Sll *)head, sizeof(list),(Sll *)&le);
   le->type = QUAD_INTERNAL;
   le->vi = pa;
   le->vj = qa;
  } else
  {
   SllAddTop( (Sll *)head, sizeof(list) );
   le = *head;
   le->type = QUAD_INTERNAL;
   le = *head;
   le->vi = pa;
   le->vj = qa;
  }
 }
}
#endif

#if QUAD_INTERS
/* ======================= add_adjnode ==================================== */

static void add_adjnode(

list    *l,
adjnode **head

)

{
 adjnode *anew;

 /* Allocate new element */

 SllAddEnd( (Sll *)head, sizeof(adjnode), (Sll *)&anew );

 /* Update the new element */

 anew->p = l;
}

/* ======================= del_adjnode ==================================== */

static void del_adjnode(

int     pa,
int     qa,
adjnode **head

)

{
 adjnode *hea;

 /* Inicializate the element */

 hea = *head;

 /* If the adjacent vector is empty return */

 if( hea == NULL ) return;

 /* Find the edge and delete it from the adjacent vector */

 while( hea )
  {
   if((hea->p->vi==pa || hea->p->vi==qa) && (hea->p->vj==pa || hea->p->vj==qa))
    {
     SllDelete( (Sll *)head, (Sll)hea );
     break;
    }
   hea = hea->next;
  }
}
#endif

/* ======================= msh_cedge_end  ================================= */

static list *msh_cedge_end(

list *head
)
{
 list *elem;

 /* Returns the last element of the cedge */

 elem = head;
 while( elem->next != NULL )
 {
  elem = elem->next;
 }
 return elem;
}

/* ======================= msh_cedge_add  ================================= */

static list *msh_cedge_add(

int   pa,
int   qa,
list *head

)
{
 list *le=NULL;
 list *elem=NULL;

 le = head;
 while( le )
 {
  if( ( le->vi == pa || le->vi == qa ) &&
      ( le->vj == pa || le->vj == qa )  )
  {
   elem = le;
   break;
  }
  else
   elem = NULL;

  le = le->next;
 }

 /* If edge already exist then no edge should be introduced in cedge
    structure */

 if( elem != NULL )
 {
  return NULL;
 }
 else
 {
  SllAddEnd( (Sll *)&head, sizeof(list),(Sll *)&le);
  le->type = QUAD_INTERNAL;
  le->vi = pa;
  le->vj = qa;
  return le;
 }
}

/* ======================= quadtree_edge  ================================= */

static int quadtree_edge(
 int  vi,
 int  vj
)
{
 if ( (vi>=nnogeometric) && (vj>=nnogeometric) ) return 1;

 return 0;
}

/* ========================= find_edge_boundary =========================== */

static int find_edge_boundary(list *head,int vi,int vj,
                              int *pmidle)
{
 list *firstedge = NULL;

 firstedge = head;

 while ( (firstedge != NULL) )
 {
  if (firstedge->type == QUAD_BOUNDARY)
  {
   if (firstedge->next != NULL )
   {
    if ( ( (firstedge->vi == vi) && (firstedge->next->vj == vj) ) ||
         ( (firstedge->vi == vj) && (firstedge->next->vj == vi) )  )
    {
     *pmidle = firstedge->vj;
     return 1;
    }
    firstedge = firstedge->next->next;
   }
  }
  else
  {
   firstedge = firstedge->next;
  }
 }

 return 0;
}

#if 0
static int find_edge_boundary(list *head,int vi,int vj,
                              int *pmidle)
{
 list *firstedge = NULL;
 int aux;

 if ( ( (vi>vj) && (vj!=0) )  || (vi==0) )
 {
  aux = vi;
  vi  = vj;
  vj  = aux;
 }

 firstedge = head;

 while ( (firstedge != NULL) )
 {
  if (firstedge->type == QUAD_BOUNDARY)
  {
   if (firstedge->next != NULL )
   {
    if ( (firstedge->vi == vi) && (firstedge->next->vj == vj) )
    {
     *pmidle = firstedge->vj;
     return 1;
    }
    firstedge = firstedge->next->next;
   }
  }
  else
  {
   firstedge = firstedge->next;
  }
 }

 return 0;
}
#endif

/* ========================= find_meshv =================================== */

static int find_meshv( int x, int y, int vi, int vj, face_t *f,
		       mesh_vertex_t *vhe, mesh_vertex_t **no_vector )
{
 int i;
 int nnoaux = 0;
 mesh_vertex_t *vertex;

 /* Verify if vertex already exists in the boundary */

 if( (vi < nnogeometric) && (vj < nnogeometric) )
 {
  /* Find mesh vertex in boundary */

  while( f != NULL )
  {
   /* Verify if this edge belongs to this face */

   if( ( (vi < (nnoaux + f->n)) && (vi >= nnoaux) ) &&
       ( (vj < (nnoaux + f->n)) && (vj >= nnoaux) )  )
   {
    if( ABS(vj - vi) != 2 )
    {
     if( ((vi == nnoaux) && (vj == (f->n + nnoaux - 2))) ||
         ((vj == nnoaux) && (vi == (f->n + nnoaux - 2)))  )
      return (f->n + nnoaux - 1);     /* last edge of boundary */
    }
    else
    {
     return ((vj + vi) / 2);          /* some edge of boundary */
    }
   }

   /* Increment number of vertexs */

   nnoaux = nnoaux + f->n;

   /* Get the next hole if exists */

   f = f->next;
  }
 }

 /* Verify if vertex already exists in the quadtree */

 if( (vi >= nnogeometric) && (vj >= nnogeometric) )
 {
  vertex = vhe;
  while( vertex )
  {
   if( vertex->x == x && vertex->y == y )
    return  vertex->nno;
   vertex = vertex->next;
  }
 }

 /* Verify if vertex already exists in node vector */

 for( i = 0; i < nnos; i++ )
 {
  if( ABS(x-(*no_vector)[i].x) < RTOL && ABS(y-(*no_vector)[i].y) < RTOL )
   return i;
 }

 /* If we get here, the vertex doesn't exists */

 return -1;
}

#if 0
static int find_meshv( int x, int y, int vi, int vj, mesh_vertex_t *vhe,
		       mesh_vertex_t **no_vector )
{
 int i;
 mesh_vertex_t *vertex;

 #if 0
 /* Verify if vertex already exists in the boundary */

 if( (vi < nnogeometric) && (vj < nnogeometric) )
 {
  if( ABS(vj - vi) != 2 ) /* last edge of the boundary. Doesn't work in holes */
   return (nnogeometric - 1);
  else
   return ((vj + vi) / 2);
 }
 #else
 /* Verify if vertex already exists in the boundary */

 for( i = 0; i < nnogeometric; i++ )
 {
  if( ABS(x-(*no_vector)[i].x) < RTOL && ABS(y-(*no_vector)[i].y) < RTOL )
   return i;
 }
 #endif

 /* Verify if vertex already exists in the quadtree */

 if( (vi >= nnogeometric) && (vj >= nnogeometric) )
 {
  vertex = vhe;
  while( vertex )
  {
   if( vertex->x == x && vertex->y == y )
    return  vertex->nno;
   vertex = vertex->next;
  }
 }

 /* Verify if vertex already exists in node vector */

 for( i = 0; i < nnos; i++ )
 {
  if( ABS(x-(*no_vector)[i].x) < RTOL && ABS(y-(*no_vector)[i].y) < RTOL )
   return i;
 }

 /* If we get here, the vertex doesn't exists */

 return -1;
}
#endif

/* ========================== test_if_isbound ============================ */

static int test_if_isbound()
{
 int t1,t2;

 t1 = ABS((candid+1)%nnogeometric - (curfirstid+1)%nnogeometric);

 t2 = ABS((candid+1)%nnogeometric - (cursecondid+1)%nnogeometric);

 if ( (curfirstid<nnogeometric) && ( ((candid%2) != 0) || (t1 == 1)) )
  return 0;
 if ( (cursecondid<nnogeometric) && ( ((candid%2) != 0) || (t2 == 1)) )
  return 0;

 return 1;
}

