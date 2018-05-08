/*
** ---------------------------------------------------------------------------
**  
** quadtri.c - This module contains routines to generate triangular elements 
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
** Created:  20-May-93		Joaquim Bento Cavalcante Neto
**
** Supervised:			Luiz Fernando Martha
**
** ---------------------------------------------------------------------------
*/          

#include <stdio.h>	
#include <stdlib.h>	
#include <math.h>

#include "quadtree.h"
#include "quadsll.h"

#define QUADTRI_C

#define ATOL 0.01
#define API acos(-1.0)

/*
** ---------------------------------------------------------------------------
** Definition for use of algorithms:
*/

#define QUAD_NORMAL 1
#define QUAD_INTERS 0
#define QUAD_DEBUG  0

/*
** ---------------------------------------------------------------------------
** Includes for DEBUG:
*/

#if QUAD_DEBUG
#include "mquadraw.h"
#include "mquadgra.h"
#endif

/*
** ---------------------------------------------------------------------------
** Local functions prototypes:
*/

static void   msh_tri_vertex( mesh_vertex_t *,face_t *, int, mesh_vertex_t **);
static void   msh_add_bdrypts( face_t *, adjnode **, list ** );
static int    msh_gen_regintpts( list *, list *,adjnode **, tree_node_t *, 
                                 mesh_vertex_t ** );
static void   msh_tri_mesh( list *, mesh_vertex_t **, int *, mesh_elem_t ** );
static void   msh_update_cedge( tree_node_t *, int, adjnode **, list ** );
static void   msh_bound_tri_mesh( mesh_elem_t **, int *, mesh_vertex_t **,
                                  mesh_elem_t ** );
static void   msh_quad_tri_mesh( tree_node_t *, mesh_vertex_t **, 
                                 mesh_elem_t ** );
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
static void   int_quadtree( list *, list *,nodeadj_t *, adjnode **, mesh_vertex_t **);
static void   find_vert_id( list *, list *, tree_node_t *, adjnode **, 
			    mesh_vertex_t ** );
static void   find_quad_id( list *, list *, tree_node_t *, side_t, adjnode **,
                            mesh_vertex_t **  );
#if !QUAD_NORMAL
static int    filtro_quad( tree_node_t *, side_t );
static int    filtro_quad_test( tree_node_t *, side_t, new_coord,
                                int, int, int, int, int, int );
#endif
#if QUAD_INTERS
static int    intersec( adjnode **, mesh_vertex_t ** );
static int    intersec_list( nodeadj_t * ,adjnode **,mesh_vertex_t **,list **);
static int    intersec_vert( int , adjnode **, mesh_vertex_t **, list ** );
static int    intersec_quad( int, int, mesh_vertex_t ** );
static int    find_intersec( int, int, list ** );
#else
static int    intersec( list *, list *, mesh_vertex_t ** );
#endif
static int    intersec_cross( int ,int  , int, int, mesh_vertex_t ** );
static double  intersec_cross_value( int ,int ,int, mesh_vertex_t ** );
static double geoangle ( void );
static void   quad_edge( tree_node_t *, new_coord, new_coord, int, int,
                         adjnode **, list ** );   
static void   update_edge( int, int, adjnode **, list ** );
#if QUAD_INTERS
static void   add_adjnode( list  *, adjnode ** );
static void   del_adjnode( int, int, adjnode ** );
#endif

/*
** ---------------------------------------------------------------------------
** Local static variables:
*/

static list      *listquad;  /*auxiliary listquad list*/
static nodeadj_t *listvertex;/*base edge vertex adjacent tree nodes list*/
static double    maxang;     /*current maximum angle for each edge*/
static double    ang;        /*computed current angle*/
static new_coord first;      /*integer coordinates of first vertex base edge*/
static new_coord second;     /*integer coordinates of second vertex base edge*/
static new_coord candidate;  /*integer coordinates of candidate vertex*/
static new_coord iddate;     /*integer coordinates of choosen vertex*/
static int       firstid;    /*id number of the first vertex base edge*/
static int       secondid;   /*id number of the second vertex base edge*/
static int       candid;     /*id number of the candidate vertex*/
static int       id;         /*id number of the choosen vertex*/

/*
** ---------------------------------------------------------------------------
** Public Function:
*/

/* ======================= trimesh_bound ================================== */

int trimesh_bound(

tree_node_t    *node,      /* quadtree tree root                       (in) */
face_t         *face,      /* geometrical boundary structure           (in) */
mesh_vertex_t  *vhe,       /* quadtree finite element vertexs structure(in) */
int            nno,        /* total finite element vertexs number      (in) */
int            nnogeom,    /* geometry finite element vertexs number   (in) */
int            *nel,       /* total finite elements number         (in/out) */
mesh_vertex_t  **no_vector,/* finite element vertexs structure        (out) */
mesh_elem_t    **el_vector /* finite elements structure               (out) */

)

{
 /* Local variables and symbols */

 list           *cedge = NULL;       /* cedge boundary structure */     
 list           *qedge = NULL;       /* cedge boundary quad structure */
 adjnode        **adj_no_vector;     /* adjacent edges vertex structure */
 mesh_elem_t    *list_vector = NULL; /* auxiliary finite elements structure */ 
 int            i;                   /* loop counter */

 /* The boundary contraction is first done  for the  geometrical  contour
    and after for the holes that could exist  in the geometry.The routine
    finishes when all the triangulation of the geometry was done.  */
     
  /* Update final output finite element vertexs structure */

  msh_tri_vertex( vhe, face, nno, no_vector ); 

  /* Inicializate adjacent finite element vertexs structure */

  adj_no_vector = ( adjnode ** )calloc( nno,sizeof( adjnode * ) );
  for( i = 0; i < nno; i++ )    adj_no_vector[i] = NULL;
 
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
   /* Generate the possible internal points,  choose  the best  triangle,
      atualizate the id for the choosen point that will be used  to  make
      the new finite element in the element_vector structure and   return
      the pointer for the edge auxiliary structure that contains the  two
      new edges,one or none to atualizate the cegde boundary structure */      

   if (!msh_gen_regintpts( cedge, qedge, adj_no_vector, node, no_vector ))
   {
    SllDelAll ( (Sll *)&cedge );
    free( adj_no_vector );
    return 0;
   }

   /* With the choosen point and your identificator( id ) updated, atua-
      lizate the list for  the   finite  elements that in the final will  
      make the output finite element structure */
    
   msh_tri_mesh( cedge, no_vector, nel, &list_vector );

   /* Atualizate the cedge boundary structure, that is, insert the   two
      new created edges, only one or none, and delete the base the  head
      cedge structure used as a base edge */
       
   msh_update_cedge( node, nnogeom, adj_no_vector, &cedge );

   /* Get the new edge in the cedge structure */
  }

  /* Now update the final output finite element vector, adding  the  finite
     elements done in the quadtree and finite elements done by the boundary
     contraction with Delaunay */

    /* Update finite elements done by boundary contraction */

    msh_bound_tri_mesh( &list_vector, nel, no_vector, el_vector );

    /* Update finite elements done by quadtree */
  
    msh_quad_tri_mesh( node, no_vector, el_vector );

  /* Release static global auxiliary listquad used */

  SllDelAll( (Sll *)&listquad );
       
  /* Release adjacent finite element vertexs structure */

  free( adj_no_vector );

  /* Return the indicated generated mesh */
  
  return 1;
}

/*
** ---------------------------------------------------------------------------
** Local functions:
*/

/* ======================= msh_tri_vertex ================================= */

static void msh_tri_vertex( 

mesh_vertex_t *v,
face_t        *f,
int           nno,
mesh_vertex_t **no_vector 
)

{
 int           nnogeomaux = 0;
 int           j;
 new_coord     p;

    /* Allocate finite element vertexs structure */
 
    *no_vector = (mesh_vertex_t *) calloc( nno, sizeof( mesh_vertex_t ) );

    /* Put the quadtree vertexs in the finite element vertexs structure */

    while( v )
    {
      /* Update vertex */

      (*no_vector)[v->nno].x      = v->x;
      (*no_vector)[v->nno].y      = v->y;
      (*no_vector)[v->nno].status =  1;
      (*no_vector)[v->nno].bdrypt = -1; 
      (*no_vector)[v->nno].fadj   = NULL;
      (*no_vector)[v->nno].next   = NULL;
      (*no_vector)[v->nno].nno    = v->nno;

      /* Get next vertex */

      v = v->next;
    }  

    /* Put the geometry nodes in the finite element vertexs structure */
 
    while( f != NULL )
    {
     for( j = 0; j < f->n; j++ )
     {
      /* Transform vertex double in integer double as in quadtree vertexs */

      param_form( f->vt[j]->x, f->vt[j]->y, &p );

      /* Update vertex */

      (*no_vector)[nnogeomaux].x       = p.ex;
      (*no_vector)[nnogeomaux].y       = p.ey;
      (*no_vector)[nnogeomaux].status  = 1;
      (*no_vector)[nnogeomaux].bdrypt  = f->vt[j]->bdrypt;
      (*no_vector)[nnogeomaux].fadj    = NULL;
      (*no_vector)[nnogeomaux].next    = NULL;
      (*no_vector)[nnogeomaux].nno     = nnogeomaux;
      nnogeomaux++;
     }

     /* Get the next hole, if exists */

     f = f->next;
    }
}

/* ======================= msh_add_bdrypts =============================== */

static void msh_add_bdrypts(

face_t  *f,            /* geometrical boundary structure      (in) */
adjnode **adj_vector,  /* adjacent edges vertex structure (in/out) */
list    **head         /* cedge boundary structure           (out) */
)

{
 /* Local variables and symbols */

 list   *mnew = NULL; /* new cedge boundary structure element */
 int    nnoaux = 0;   /* auxiliary vertex number */
 int    j;            /* loop counter */
            
 /* This routine only first put the geometrical informations in the   auxi-
    liary cedge(contract edge) structure, that will be used to do the boun-
    dary contraction.The cedge is the pointer for the cegde structure  that 
    is a single-linked list */

  /* Make the inicial cedge boundary structure and the adjacent edges vertex 
     structure*/

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

     /* Update adjacent edges vertex structure */

     #if QUAD_INTERS
     add_adjnode( mnew, &adj_vector[mnew->vi] );
     add_adjnode( mnew, &adj_vector[mnew->vj] );  
     #endif
    }

   /* Increment number of vertexs */

   nnoaux = nnoaux + f->n;

   /* Get the next hole if exist */

   f = f->next;
  }
}

/* ======================== msh_gen_regintpts ============================= */

static int msh_gen_regintpts (

list        *head ,          /* cedge boundary structure              (in) */
list        *head_quad ,     /* qedge boundary structure              (in) */
adjnode     **adj_vector,    /* adjacent edges vertex structure       (in) */
tree_node_t *node,           /* root quadtree tree node               (in) */
mesh_vertex_t **no_vector    /* finite element vertex structure       (in) */
)

{
 /* Local variables and symbols */

 neighbour_t  *locfirst = NULL;
 neighbour_t  *locfirst_prob = NULL;
 neighbour_t  *locsecond = NULL;
 neighbour_t  *locsecond_prob = NULL;
 int           rid;


  /* Inicializate tree nodes candidates vertex list */

  listvertex = NULL;

  /* Inicializate the maximum angle */

  maxang = 0.0;  

  /* Inicializate the id */

  id = -1;

  /* Get the first and second coordinates and ids of the base edge */

  first.ex = (*no_vector)[head->vi].x;
  first.ey = (*no_vector)[head->vi].y;
  firstid = head->vi;

  second.ex = (*no_vector)[head->vj].x;
  second.ey = (*no_vector)[head->vj].y;
  secondid = head->vj;
  
  /* Get the tree nodes list where the first and second vertexs are */ 

  locfirst  = locate_quad( &first, node );
  locsecond = locate_quad( &second, node );

  #if 0
  /* Get the adjacent tree nodes list of the first and second vertexs based
    on the tree nodes list where the vertexs are */

  insert_list( locfirst, &listvertex );
  insert_list( locsecond, &listvertex ); 
    
  /* Look trough the list to get the id for the best triangle of the base 
    edge, based on a Delaunay triangulation */

  int_quadtree( head, head_quad,listvertex, adj_vector, no_vector );
  #endif 

  /* If the  base  edge  didn't get an id, make the best triangle  considering 
     a new list of tree nodes to get the best vertex to do it. This  list   is
     maked in that way: get the first and second base edge vertexs tree  nodes
     and go through both tree nodes all adjacents only stoping when get a tree 
     node NODE_INTERIOR  or  NODE_VERTEX that is updated in the new list */

#if 0
  if( id == -1 )
  {
   /* Release the list before used and that didn't get and id */

   SllDelAll( (Sll *)&listvertex );
   
   /* Make the new list following description above */

   insert_prob_list( locfirst, &listvertex );
   insert_prob_list( locsecond, &listvertex );    

   /* Look trough the lists to get the id for the best triangle of the base 
      edge, based on a Delaunay triangulation */

   int_quadtree( head, listvertex, adj_vector, no_vector );  
  }
#endif

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

   int_quadtree( head, head_quad, listvertex, adj_vector, no_vector );     
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

   int_quadtree( head, head_quad, listvertex, adj_vector, no_vector );     
  }
      
  /* If the base edge even didn't get and a id, it's not possible to generate
     the mesh, that is, return zero */

  if( id == -1 )
   rid = 0;
  else
   rid = 1;
 
  /* Release memory used */

  SllDelAll( (Sll *)&locfirst );
  SllDelAll( (Sll *)&locfirst_prob );
  SllDelAll( (Sll *)&locsecond_prob );
  SllDelAll( (Sll *)&locsecond );
  SllDelAll( (Sll *)&listvertex );

  /* Return the indicated generated triangle */

  return rid;
}

/* ======================== msh_tri_mesh ================================== */

static void msh_tri_mesh(

list          *head,         /* cedge boundary structure               (in) */
mesh_vertex_t **no_vector,   /* finite element vertexs structure       (in) */
int           *nel,          /* total number of finite elements    (in/out) */
mesh_elem_t   **list_vector  /* auxiliary finite element structure (in/out) */

)

{
 mesh_elem_t *mnew;
#if QUAD_DEBUG
 int          i;     
 new_coord    p;
 vertex_t     va;
 point       *el_coord[4],mesh_coord[4]; 
#endif

 /* Create new finite element and put in the auxiliary finite element list */

 /* Allocate the new finite element */

 if( SllAddEnd( (Sll *)list_vector, sizeof(mesh_elem_t), (Sll *)&mnew) )
  {
   /* Allocate finite element connectivity */

   mnew->conect = ( mesh_vertex_t ** ) calloc( 3, sizeof(mesh_vertex_t *) );

   /* Update finite element connectivity */

   mnew->conect[0] = &((*no_vector)[head->vi]);
   mnew->conect[1] = &((*no_vector)[id]);
   mnew->conect[2] = &((*no_vector)[head->vj]);
 
   /* Update number of vertexs in the new finite element */

   mnew->n = 3;
   
   /* Update new finite element number and total number of finite elements */

   mnew->nel = (*nel)++;

#if QUAD_DEBUG
   /* Draw element test */

   for( i = 0; i < 3; i++ )
    {
     p.ex = mnew->conect[i]->x;
     p.ey = mnew->conect[i]->y;
     invparam_form( p.ex, p.ey, &va );
     mesh_coord[i].x = va.x;
     mesh_coord[i].y = va.y;

     el_coord[i] = &mesh_coord[i];
    }

    if( 1 /*flag_debug*/ )
     {
      Draw_Face( HOLLOW1, 7, 3, (point **) el_coord );
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
    
  update_edge( firstid, secondid, adj_vector, head );

  /* Atualizate the secondid e id edge */
      
  if( secondid >= nnogeom && id >= nnogeom )
    quad_edge( node, iddate, second, id, secondid, adj_vector, head );
  else
    update_edge( id, secondid, adj_vector, head );

  /* Atualizate the firstid e id edge */

  if( firstid >= nnogeom && id >= nnogeom ) 
    quad_edge( node, first, iddate, firstid, id, adj_vector, head );
  else
    update_edge( firstid, id, adj_vector, head );
}

/* ======================= msh_bound_tri_mesh ============================= */

static void msh_bound_tri_mesh( 

mesh_elem_t   **list_vector,/* auxiliary finite elements structure     (in) */
int            *nel,        /* total finite elements number            (in) */
mesh_vertex_t **no_vector,  /* output finite element vertexs structure (in) */
mesh_elem_t   **el_vector   /* output finite elements structure    (in/out) */
)

{
 int   j;
 int   vi;
 int   vj;
 double s;

 /* Allocate output finite element structure */

 *el_vector = (mesh_elem_t *) calloc( *nel, sizeof( mesh_elem_t ) );

 /* Look trough list vector to update el_vector */

 while( *list_vector )
 {
  /* Allocate connectivity in the the new output finite element */

  (*el_vector)[(*list_vector)->nel].conect = 
  ( mesh_vertex_t **) calloc( (*list_vector)->n, sizeof( mesh_vertex_t *) );

  /* Verify clock wise orientation because this algorithm gets finite elements
    in clock wise orientation */

  s = 0.0;
  for( j = 0; j < (*list_vector)->n; j++ )
  {
   vi = (*list_vector)->conect[j]->nno;
   vj = (*list_vector)->conect[(j+1)%((*list_vector)->n)]->nno;

   s += (*no_vector)[vi].x * (*no_vector)[vj].y - 
        (*no_vector)[vj].x * (*no_vector)[vi].y ;
  }

  /* Update new output finite element in clock wise orientation */

  if (s > 0.0)
  {
   for( j = 0; j < (*list_vector)->n; j++ )
    (*el_vector)[(*list_vector)->nel].conect[j] = 
                 (*list_vector)->conect[((*list_vector)->n-1) - j];
  }
  else
  {
   for( j = 0; j < (*list_vector)->n; j++ )
    (*el_vector)[(*list_vector)->nel].conect[j] = (*list_vector)->conect[j];
  }
     
  /* Update number of vertexs in the new output finite element */

  (*el_vector)[(*list_vector)->nel].n = (*list_vector)->n;

  /* Update the new output finite element number */

  (*el_vector)[(*list_vector)->nel].nel = (*list_vector)->nel;
  
  /* Release list_vector element */

  free( (*list_vector)->conect );
  SllDelete( (Sll *)list_vector, (Sll)*list_vector );
 } 
}
  
/* ======================= msh_quad_tri_mesh ============================== */

static void msh_quad_tri_mesh( 

tree_node_t    *node,       /* quadtree tree node root                 (in) */
mesh_vertex_t **no_vector,  /* output finite element vertexs structure (in) */
mesh_elem_t   **el_vector   /* output finite element structure     (in/out) */
)

{
 mesh_elem_t       *elme;
 int               i;
 int               j;
 int               vi;
 int               vj;
 double             s;

 /* Traverse function recursivily until get all quadtree finite elements */
 
 /* Verify if the tree node has no child, that is, if it's a leaf and then put
    its finite element in output finite elements structure */

 if( node->child != NULL )
   for( i = 0; i < 4; i++  ) 
    msh_quad_tri_mesh( node->child[i], no_vector, el_vector );
 else if( node->elmhead != NULL )
      {
        elme = node->elmhead;
        while( elme )
        {
         /* Allocate connectivity in the the new output finite element */

         (*el_vector)[elme->nel].conect = 
           ( mesh_vertex_t **) calloc( elme->n, sizeof( mesh_vertex_t *) );

         /* Verify clock wise orientation because this algorithm gets finite 
            elements in clock wise orientation */

         s = 0.0;
         for( j = 0; j < elme->n; j++ )
         {
          vi = elme->conect[j]->nno;
          vj = elme->conect[(j+1)%(elme->n)]->nno;

          s += (*no_vector)[vi].x * (*no_vector)[vj].y - 
               (*no_vector)[vj].x * (*no_vector)[vi].y ;
         }

         /* Update new output finite element in clock wise orientation */

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
nodeadj_t     *list,
adjnode       **adj_vector,
mesh_vertex_t **no_vector 
)
 
{ 
  /* Look trough the list to get the best id based in a Delaunay tecnique */
  
  while( list != NULL ) 
  {
   switch( list->node->code )
    {
     /* It has two kinds of procedures depending on if the tree node is 
        NODE_VERTEX or NODE_INTERIOR */

     case NODE_VERTEX :
      {
       find_vert_id( head, head_quad, list->node, adj_vector, no_vector );
       break;
      }
     case NODE_INTERIOR :
      {
       find_quad_id( head, head_quad, list->node, list->side, adj_vector, no_vector );
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
tree_node_t   *node,
adjnode       **adj_vector,
mesh_vertex_t **no_vector
)

{
  tree_vertex_t  *elme;
  int            intersecv;

  /* Look on the vertexs list in the NODE_VERTEX to see if it's the best id */

  /* Get the first vertex in the tree node */

  elme = node->vhead;

  while( elme )
  {
   /* Get the candidate id */

   candid = elme->id;

   /* If the candid is the same of the id, this vertex was already looked */

   if( candid != id )
    {
     /* Get the candidate coordinates */

     candidate = elme->coord;

     /* Calculate the angle */

     ang = geoangle( );

     /* If angle is bigger than the maxang, update the id like the candid */

     if( ang > maxang )
      {
       /* Check intersection of the two new edges that will be formed( candid-
         firstid and candid-secondid ) against edges that already exits */

       #if QUAD_INTERS
       intersecv = intersec( adj_vector, no_vector );
       #else
       intersecv = intersec( head, head_quad, no_vector );
       #endif

       /* If there's no intersection, update the candid,iddate and maxang */

       if( intersecv )
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

 list	     *head,
 list        *head_quad,
 tree_node_t *node,
 side_t      side,
 adjnode     **adj_vector,
 mesh_vertex_t **no_vector
 )

{
 mesh_elem_t       *elme;
 int               j;
 int               intersecq;
 int               filtroquad;
 
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

             if( candid != id )
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

                   /* If there's no intersection,upadte id,iddate and maxang */

                   if( intersecq )
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
 
   cf = intersec_vert( firstid, adj_vector, no_vector, &listv );
   if( cf == 0 )  
    {
     SllDelAll( (Sll *)&listv );
     return cf;
    }

   /* Test intersection of new edges against the edges formed by secondid */

   cs = intersec_vert( secondid, adj_vector, no_vector, &listv );
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
  cff = intersec_cross( firstid, candid, head->vi, head->vj, no_vector );
  if( cff  == 0 ) return cff;
  cfs = intersec_cross( candid, secondid, head->vi, head->vj, no_vector );
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
       
       if( ( (t->p->vi)!=firstid && (t->p->vi)!=candid ) &&
           ( (t->p->vj)!=firstid && (t->p->vj)!=candid )  )
        {
         codef = 1;     
         cff = intersec_cross( firstid, candid,t->p->vi,t->p->vj, no_vector );
         if( cff == 0 )  return cff;
        }

       /* Second test the edge against edge formed by firstid and candid.If the
          edge has one of their vertexs equal of one of the vertexs of the edge
          formed by firstid and candid then it  not should be tested    because
          they won't intersect */      

       if( ( (t->p->vi)!=secondid && (t->p->vi)!=candid ) &&
           ( (t->p->vj)!=secondid && (t->p->vj)!=candid )  )
        {
         codes = 1;
         cfs = intersec_cross( candid,secondid,t->p->vi,t->p->vj, no_vector );
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
       
     if( ( df!=firstid && df!=candid ) &&
         ( ds!=firstid && ds!=candid )  )
      {
       cff = intersec_cross( firstid, candid, df, ds, no_vector );
       if( cff == 0 )  return cff;
      }
 
     /* Second test the edge against edge formed by firstid and candid.If the
        edge has one of their vertexs equal of one of the vertexs of the edge
        formed by firstid and candid then it  not should be tested    because
        they won't intersect */      

     if( ( df!=secondid && df!=candid ) &&
         ( ds!=secondid && ds!=candid )  )
      {
       cfs = intersec_cross( candid, secondid, df, ds, no_vector );
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
 double c0,c1,c2,c3;

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

static double intersec_cross_value(

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

/* ====================== geoangle ======================================== */

static double geoangle ( void )
{
 double ikx, iky, jkx, jky;
 double angle;
 double xprod;

 /* Look if the side is on the left of the contour;if it's the angle is 0,that
    is, this vertex could not be used.This test is done using cross product */

 xprod = - (( second.ex - first.ex ) * ( candidate.ey - first.ey ))
         + (( second.ey - first.ey ) * ( candidate.ex - first.ex )) ;

 /* If the vertex is inside the contour, compute the angle */

 if( xprod >= 0.0 )
  {

   ikx = first.ex - candidate.ex;
   iky = first.ey - candidate.ey;
   jkx = second.ex - candidate.ex;
   jky = second.ey - candidate.ey;

   angle = ((ikx * jkx) + (iky * jky)) /
		 (sqrt(ikx*ikx + iky*iky) * sqrt(jkx*jkx + jky*jky));

   if( angle < -1.0 )  return  0.0;
   if( angle > 1.0  )  angle = 1.0;

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

 update_edge( pa, qa, adj_vector, head );

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
       
/* ====================== update_edge ===================================== */

#if QUAD_INTERS
static void update_edge( 

int   pa,             /* first edge vertex id                (in) */
int   qa,             /* second edge vertex id               (in) */
adjnode **adj_vector, /* adjacent edges vertex structure (in/out) */
list  **head          /* cedge boundary structure        (in/out) */

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
    edges vertex structure to see if this vertex already make this edge, avoi-
    ding searchs in all cedge structure and otimizing the search to see if the
    edge already exists */
 
 if( pa < qa )  ref = pa;
 else           ref = qa;

 /* Get the head of edge list of this vertex */

 ed = adj_vector[ref];

 /* Look in list to see if this vertex already make the edge to be updated.
    If the edge already exists, get the pointer elem to delete directly in the
    cedge structure without search. If not, introduce the edge in cedge struc-
    ture and upadte adj_vector */

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

   /* Update adjacent edges vertex structure for the two vertexs. Here has a 
      point: the search is done only for the lower vertex but I should update 
      the adjacent structure for the two vertexs becuse the two will be used  
      in the intersection test */

   del_adjnode( pa, qa, &adj_vector[pa] );
   del_adjnode( pa, qa, &adj_vector[qa] );
 }
 else
  {
   /* Update cedge structure.The edge not exist yet and should be introduced */

   SllAddTop( (Sll *)head, sizeof(list) );
   le = *head;
   le->vi = pa;
   le->vj = qa;

   /* Update the adjacent vector for the two vertexs */

   add_adjnode( le, &adj_vector[pa] );
   add_adjnode( le, &adj_vector[qa] );
  }
}
#else
static void update_edge( 

int   pa,             /* first edge vertex id                (in) */
int   qa,             /* second edge vertex id               (in) */
adjnode **adj_vector, /* adjacent edges vertex structure (in/out) */
list  **head          /* cedge boundary structure        (in/out) */

)

{ 
 list  *le   = NULL;         /* auxiliary cedge boundary structure pointer */
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
   /*U pdate cedge structure */

   SllDelete( (Sll *)head, (Sll)elem );
  }
 else
  {
   /* Update cedge structure.The edge not exist yet and should be introduced */

   SllAddTop( (Sll *)head, sizeof(list) );
   le = *head;
   le->vi = pa;
   le->vj = qa;
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

