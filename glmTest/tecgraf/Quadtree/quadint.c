/*
** ---------------------------------------------------------------------------
**
** quadint.c  - This module contains routines to generate elements in the 
**              interior tree cells by using template in a quadtree algorithm.
**
** ---------------------------------------------------------------------------
**
** Version: 0-001
**
** Created:  20-May-93		Joaquim Bento Cavalcante Neto
**    Stolen from Marcelo Tilio.
**
** Modified: 01-Jun-96		Joaquim Bento Cavalcante Neto &&
**                              Eduardo Setton Sampaio Silveira
**    Created functions nnogetgoemitn and nnogemqua.
**
** Supervised:			Luiz Fernando Martha
**
** ---------------------------------------------------------------------------
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "quadtree.h"
#include "quadsll.h"

/* #include "cd.h" */

#define QUADINT_C

/*
** ---------------------------------------------------------------------------
** Local functions prototypes:
*/

static void create_trimesh( tree_node_t * );
static void internal_trimesh( tree_node_t * );
static void create_tri_elem( int, mesh_vertex_t *[][4], tree_node_t * ); 
static void create_polygons_quadmesh( tree_node_t * );
static void internal_quadmesh( tree_node_t * );
static void create_elem_quad( int, mesh_vertex_t *[][4], tree_node_t * );
static void create_elem_tri( int, mesh_vertex_t *[][4], tree_node_t * );
static void create_quadmesh( tree_node_t * );
static void quadrilaterals_quadmesh( mesh_elem_t * );
static mesh_vertex_t *find_meshv( int x, int y );
static void include_mid_node( side_t, int, tree_node_t *, int, int, 
                              mesh_vertex_t **, int *, int * );
 
/*
** ---------------------------------------------------------------------------
** Local Variables:
*/

static int nnoint;           /* finite element vertexs number in quadtree */
static int nnointquad;
static int nelint;           /* finite elements number in quadtree        */
static int nnogeomint;       /* finite element vertexs number in geometry */
static mesh_vertex_t *vhead; /* quadtree head finite element vertexs      */

/*
** ---------------------------------------------------------------------------
** Public Functions:
*/

/* ======================= trimesh_int ==================================== */

void trimesh_int( tree_node_t *node, int nnogeom,
                  int *nno, int *nel, int *nelquadtree, 
                  mesh_vertex_t **v )
{
 /* Inicializate static local variables*/

 nnoint = nnointquad = nelint = nnogeomint = 0;
 vhead = NULL;

 /* Inicializate geometry finite element vertexs number */

 nnogeomint = nnogeom; 
 
 /* Create quadtree interior tree nodes finite element structure */

 create_trimesh( node );

 /* Update total finite element vertexs number */

 *nno = nnoint + nnogeom;
 
 /* Update total finite elements number */

 *nel = *nelquadtree = nelint;

 /* Update quadtree finite element vertexs auxiliary list */

 *v = vhead;
}

/* ======================== quadmesh_int ================================== */

void quadmesh_int( tree_node_t *node, int nnogeom, 
                  int *nno, int *nel, int *nelquadtree, 
                  mesh_vertex_t **v )
{
 /* Inicializate static local variables*/

 nnoint = nnointquad = nelint = nnogeomint = 0;
 vhead = NULL;

 /* Inicializate geometry finite element vertexs number */

 nnogeomint = nnogeom;
 
 /* Create quadtree interior tree nodes finite element structure */

 create_polygons_quadmesh( node );
 create_quadmesh( node );

 /* Update total finite element vertexs number */

 *nno = nnoint + nnogeom;
 
 /* Update total finite elements number */

 *nel = *nelquadtree = nelint;

 /* Update quadtree finite element vertexs auxiliary list */

 *v = vhead;
}

/* =========================== get_nnoint ================================== */
 
int get_nnoint(void)
{
 return nnointquad;   
}
   
/* ========================== get_nnogeomint =============================== */

int get_nnogeomint(void)
{
 return nnogeomint;
}

/*
** ---------------------------------------------------------------------------
** Local functions:
*/

/* ======================= create_trimesh ================================= */

static void create_trimesh( tree_node_t *node )
{
 int     i;

 /* The routine traverse recursivilying the tree and if the tree node  has  no
    child, that is, if the tree node is a leaf and it is NODE_INTERIOR, create
    a finite element inside it */

 if( node->child != NULL )
   for( i = 0; i < 4; i++ ) create_trimesh( node->child[i] );
 else if( node->code == NODE_INTERIOR ) internal_trimesh( node );    
}

/* ======================= internal_trimesh =============================== */

static void internal_trimesh( tree_node_t *node )
{
 int            i, j, x, y, size, ind, aux;
 int            type = 0;
 int            mid_node[4];
 mesh_vertex_t  *v[9], *conect[6][4];

 /* Calculate the tree node size */

 size = pow_2( 14 - node->depth );

 /* Look trough the tree node childs to get an auxiliary v list with the 
    corner vertexs of the tree node */

 for( i = 0, y = node->coord.ey; i < 2; y += size, i++ )
   for( j = 0, x = node->coord.ex; j < 2; x += size, j++ )
   {
     ind = abs( 6*i - 2*j );  
     v[ind] = find_meshv( x, y );
   }

 /* Get the type of template that will be used to form a new finite element,
    depending on the tree node kind of adjacent, and update vlist with   the
    middle vertexs tree node adjacents, if exists */

 include_mid_node( BOTTOM, 1, node, node->coord.ex + size / 2,
                   node->coord.ey, v, &type, mid_node );
 include_mid_node( RIGHT, 3, node, node->coord.ex + size,
                   node->coord.ey + size / 2, v, &type, mid_node );
 include_mid_node( TOP, 5, node, node->coord.ex + size / 2,
                   node->coord.ey + size, v, &type, mid_node );
 include_mid_node( LEFT, 7, node, node->coord.ex,
                   node->coord.ey + size / 2 , v, &type, mid_node );

 /* With the template type, create a new finite element inside tree node */

 switch( type )
 {
   case 0:

     /* Get finite element connectivity */

     for( i = 0; i < 2; i++ )
       for( j = 0; j < 3; j++ )
       {
         ind = 4*i + 2*(j+1);
         if( ind > 7 ) ind -= 8;
         conect[i][j] = v[ind];
       }

     /* Create the finite element into the tree node */

     create_tri_elem( 2, conect, node );
     break;

   case 1:

     /* Get finite element connectivity */

     for( i = 0; i < 3; i++ )
     {
       conect[i][0] = v[mid_node[0]];
       for( j = 1; j < 3; j++ )
       {
         if( ( ind = mid_node[0] + 2*i + 2*j - 1 ) > 7 ) ind -= 8;
         conect[i][j] = v[ ind ];
       }
     }

     /* Create the finite element into the tree node */
    
     create_tri_elem( 3, conect, node );
     break;

   case 2:
     
     /* Get finite element connectivity */

     v[8] = find_meshv( node->coord.ex + size/2, node->coord.ey + size/2 );
     if( mid_node[1] - mid_node[0] != 4 )
     {
       if( mid_node[1] - mid_node[0] != 2 )
       {
         aux = mid_node[1];
         mid_node[1] = mid_node[0];
         mid_node[0] = aux;
       }
       for( i = 0; i < 3; i++ )
       {
         if( ( ind = mid_node[0] + i ) > 7 ) ind -= 8;
         conect[0][i] = v[ ind ];
       }
       for( i = 1; i < 6; i++ ) conect[i][0] = v[8];

       conect[1][1] = conect[5][2] = v[mid_node[1]];
       conect[5][1] = conect[4][2] = v[mid_node[0]];

       for( i = 1; i < 4; i++ )
       {
         if( ( ind = mid_node[0] + 2*i + 1 ) > 7 ) ind -= 8;
         conect[i][2] = conect[i+1][1] = v[ ind ];
       }
     }
     else
     {
       for( i = 0; i < 6; i++ ) conect[i][0] = v[8];

       conect[0][1] = conect[5][2] = v[mid_node[0]];
       conect[2][2] = conect[3][1] = v[mid_node[1]];

       for( i = 0; i < 2; i++ )
         for( j = 0; j < 2; j++ )
         {
           if( (ind = mid_node[i] + 2*j + 1) > 7 ) ind -= 8;
           conect[3*i+j][2] = conect[3*i+j+1][1] = v[ind];  
         }
     }

     /* Create the finite element into the tree node */

     create_tri_elem( 6, conect, node );
     break;

   case 3:

     /* Get finite element connectivity */

     if( mid_node[1] - mid_node[0] > 2 )
     {
       aux = mid_node[1];
       mid_node[1] = mid_node[0];
       mid_node[0] = aux;
       aux = mid_node[1];
       mid_node[1] = mid_node[2];
       mid_node[2] = aux;
     }
     else if( mid_node[2] - mid_node[1] > 2 )
          {
            aux = mid_node[1];
            mid_node[1] = mid_node[2];
            mid_node[2] = aux;
            aux = mid_node[1];
            mid_node[1] = mid_node[0];
            mid_node[0] = aux;
          }
     for( i = 0; i < 3; i++ )
     {
       conect[0][i] = v[mid_node[i]];
       for( j = 1; j < 3; j++ )
       {
         if( (ind = mid_node[j-1] + i) > 7 ) ind -= 8;
         conect[j][i] = v[ ind ];
       }
     }
     conect[3][0] = conect[4][0] = v[mid_node[2]];
     if( (ind = mid_node[2] + 3) > 7 ) ind -= 8;
     conect[3][2] = conect[4][1] = v[ind];
     if( (ind = mid_node[2] + 1) > 7 ) ind -= 8;
     conect[3][1] = v[ind];
     conect[4][2] = v[mid_node[0]];

     /* Create the finite element into the tree node */

     create_tri_elem( 5, conect, node );
     break;

   case 4:
    
     /* Get finite element connectivity */

     for( i = 0; i < 4; i++ )
       for( j = 0; j < 3; j++ )
       {
	 if( (ind = mid_node[i] + j) > 7 ) ind -= 8;
         conect[i][j] = v[ ind ];
       }
     conect[4][1] = conect[5][0] = v[mid_node[1]];
     conect[4][2] = conect[5][2] = v[mid_node[3]];
     conect[4][0] = v[mid_node[0]];
     conect[5][1] = v[mid_node[2]];

     /* Create the finite element into the tree node */

     create_tri_elem( 6, conect, node );
     break;
 }
}

/* ====================== create_tri_elem ================================= */

static void create_tri_elem( int n, mesh_vertex_t *conect[][4],
                             tree_node_t *node ) 
{
 int          i, j;
 mesh_elem_t  *f;

 /* With the finite element connectivity and the number of finite elements 
    that will be formed, make them and put into the tree node */

 /* Make n finite elements and put them into the tree node */

 for( i = 0; i < n; i++ )
 {
   /* Allocate finite element auxiliary structure */

   f = ( mesh_elem_t * ) calloc( 1, sizeof(mesh_elem_t) );
   f->conect = ( mesh_vertex_t ** ) calloc( 3, sizeof(mesh_vertex_t *) );

   /* Put the connectivity inside the finite element auxiliary structure */

   for( j = 0; j < 3; j++ )
     f->conect[j] = conect[i][j];

   /* Update number of vertexs and number of finite elements in the quadtree */

   f->n = 3;
   f->nel = (nelint++);
   
   /* Put the finite element auxiliary structure inside the tree node */
  
   f->next = node->elmhead;
   node->elmhead = f;   
 }
}

/* ================== create_polygons_quadmesh ============================ */

static void create_polygons_quadmesh( tree_node_t *node )
{
 int     i;

 /* The routine traverse recursivilying the tree and if the tree node  has  no
    child, that is, if the tree node is a leaf and it is NODE_INTERIOR, create
    polygons inside it */

 if( node->child != NULL )
   for( i = 0; i < 4; i++ ) create_polygons_quadmesh( node->child[i] );
 else if( node->code == NODE_INTERIOR ) internal_quadmesh( node );    
}

/* ==================== internal_quadmesh ================================= */

static void internal_quadmesh( tree_node_t *node )
{
 int            i, j, x, y, size, ind, aux;
 int            type = 0;
 int            mid_node[4];
 mesh_vertex_t  *v[9], *conect[4][4];


 /* Calculate the tree node size */

 size = pow_2( 14 - node->depth );

 /* Look trough the tree node childs to get an auxiliary v list with the 
    corner vertexs of the tree node */

 for( i = 0, y = node->coord.ey; i < 2; y += size, i++ )
   for( j = 0, x = node->coord.ex; j < 2; x += size, j++ )
   {
     ind = abs( 6*i - 2*j );
     v[ind] = find_meshv( x, y );
     nnointquad = nnoint;
   }    

 /* Get the type of template that will be used to form a new finite element,
    depending on the tree node kind of adjacent, and update v list with   the
    middle vertexs tree node adjacents, if exists */

 include_mid_node( BOTTOM, 1, node, node->coord.ex + size / 2,
                   node->coord.ey, v, &type, mid_node );
 include_mid_node( RIGHT, 3, node, node->coord.ex + size,
                   node->coord.ey + size / 2, v, &type, mid_node );
 include_mid_node( TOP, 5, node, node->coord.ex + size / 2,
                   node->coord.ey + size, v, &type, mid_node );
 include_mid_node( LEFT, 7, node, node->coord.ex,
                   node->coord.ey + size / 2 , v, &type, mid_node );

 /* With the template type, create a new finite element inside tree node */

 switch( type )
 {
   case 0:
  
     /* Get finite element connectivity */

     for( j = 0; j < 4; j++ )
     {
       conect[0][j] = v[2*j];
     }

     /* Create the finite element into the tree node */

     create_elem_quad( 1, conect, node );
     break;

   case 1:

     /* Get finite element connectivity */
 
     v[8] = find_meshv( node->coord.ex + size/2, node->coord.ey + size/2 );
     for( i = 0; i < 2; i++ )
     {
       conect[i][0] = v[8];
       for( j = 1; j < 4; j++ )
       { 
         if( j == 2 )
	 { 
           if( ( ind = mid_node[0] + 6* i + 1 ) > 7 ) ind -= 8;
	 }
	 else  
	 {
           if( ( ind = mid_node[0] + 5*i + 3*(j/2) ) > 7 ) ind -= 8;
	 }
         conect[i][j] = v[ ind ];
       }
     }

     /* Create the finite element into the tree node */

     create_elem_quad( 2, conect, node );
     
     /* Get finite element connectivity */

     conect[0][0] = v[8];
     for( j = 1; j < 3; j++ )
     {
       if( ( ind = mid_node[0] + 2*j + 1 ) > 7 ) ind -= 8;
       conect[0][j] = v[ind];
     }

     /* Create the finite element into the tree node */

     create_elem_tri( 1, conect, node );
     break;

   case 2:

     /* Get finite element connectivity */

     if( mid_node[1] - mid_node[0] != 4 )
     {
       v[8] = find_meshv( node->coord.ex + size/2, node->coord.ey + size/2);
       if( mid_node[1] - mid_node[0] != 2 )
       {
         aux = mid_node[1];
         mid_node[1] = mid_node[0];
         mid_node[0] = aux;
       }
       for( i = 0; i < 3; i++ ) conect[i][0] = v[8];
       for( i =  0; i < 2; i++ ) 
       {
         if( ( ind = mid_node[i] + 1 ) > 7 ) ind -= 8;
         conect[i][2] = v[ind];
       }
       conect[0][1] = conect[2][3] = v[mid_node[0]];
       conect[1][1] = conect[0][3] = v[mid_node[1]];
       if( ( ind = mid_node[1] + 3 ) > 7 ) ind -= 8;
       conect[1][3] = conect[2][1] = v[ind];
       if( ( ind = mid_node[1] + 5 ) > 7 ) ind -= 8;
       conect[2][2] = v[ind];

       /* Create the finite element into the tree node */

       create_elem_quad( 3, conect, node );
     }
     else
     {

       /* Get finite element connectivity */

       for( i = 0; i < 2; i++ )
       {
         conect[i][3] = v[mid_node[1-i]];  
         for( j = 0; j < 3; j++ )
         {
           if( (ind = mid_node[i] + j + (j/2) ) > 7 ) ind -= 8;
           conect[i][j] = v[ind];
         }
       }	
       create_elem_quad( 2, conect, node ); 
     }

     /* Create the finite element into the tree node */

     /* create_elem_quad( 2, conect, node ); */
     break;

   case 3:

      /* Get finite element connectivity */

     if( mid_node[1] - mid_node[0] > 2 )
     {
       aux = mid_node[1];
       mid_node[1] = mid_node[0];
       mid_node[0] = aux;
       aux = mid_node[1];
       mid_node[1] = mid_node[2];
       mid_node[2] = aux;
     }
     else if( mid_node[2] - mid_node[1] > 2 )
          {
            aux = mid_node[1];
            mid_node[1] = mid_node[2];
            mid_node[2] = aux;
            aux = mid_node[1];
            mid_node[1] = mid_node[0];
            mid_node[0] = aux;
          }  
	  
     conect[0][0] = v[mid_node[1]];	  
     if( ( ind = mid_node[1] + 1 ) > 7 ) ind -= 8;
     conect[0][1] = v[ind];	  
     conect[0][2] = v[mid_node[2]];

     /* Create the finite element into the tree node */

     create_elem_tri( 1, conect, node );
     
     /* Get finite element connectivity */

     for( j = 0; j < 3; j++ )
     {
       if( (ind = mid_node[0] + j) > 7 ) ind -= 8;
       conect[0][j] = v[ ind ];
     }	       
     for( j = 0; j < 3; j++ )
     {
       if( (ind = mid_node[2] + j + j/2 ) > 7 ) ind -= 8;
       conect[1][j+1] = v[ ind ];
     }	       
     conect[0][3] = v[mid_node[2]];
     conect[1][0] = v[mid_node[0]];
     
     /* Create the finite element into the tree node */

     create_elem_quad( 2, conect, node );
     break;

   case 4:

     /* Get finite element connectivity */

     v[8] = find_meshv( node->coord.ex + size/2, node->coord.ey + size/2 );
     for( i = 0; i < 4; i++ ) conect[i][3] = v[8];
     for( i = 0; i < 4; i++ )
       for( j = 0; j < 3; j++ )
       {
         if( (ind = mid_node[i] + j) > 7 ) ind -= 8;
         conect[i][j] = v[ ind ];
       }

     /* Create the finite element into the tree node */

     create_elem_quad( 4, conect, node );
     break;
 }



#if 0 /* para teste no elemento */
vis_interior_mesh( node, CD_MAGENTA);
#endif



}

/* ==================== create_elem_quad ================================== */

static void create_elem_quad( int n, mesh_vertex_t *conect[][4], 
                              tree_node_t *node ) 
{
 int          i, j;
 mesh_elem_t  *f;

 /* With the finite element connectivity and the number of finite elements 
    that will be formed, make them and put into the tree node */

 /* Make n finite elements and put them into the tree node */

 for( i = 0; i < n; i++ )
 {
   /* Allocate finite element auxiliary structure */

   f = ( mesh_elem_t * ) calloc( 1, sizeof(mesh_elem_t) );
   f->conect = ( mesh_vertex_t ** ) calloc( 4, sizeof(mesh_vertex_t *) );

   /* Put the connectivity inside the finite element auxiliary structure */

   for( j = 0; j < 4; j++ )
     f->conect[j] = conect[i][j];

   /* Update number of vertexs  */

   f->n = 4;
   
   /* Put the finite element auxiliary structure inside the tree node */

   f->next = node->elmhead;
   node->elmhead = f;
 }
}  

/* ==================== create_elem_tri =================================== */

static void create_elem_tri( int n, mesh_vertex_t *conect[][4], 
                             tree_node_t *node ) 
{
 int          i, j;
 mesh_elem_t  *f;

 /* With the finite element connectivity and the number of finite elements 
    that will be formed, make them and put into the tree node */

 /* Make n finite elements and put them into the tree node */

 for( i = 0; i < n; i++ )
 {
   /* Allocate finite element auxiliary structure */

   f = ( mesh_elem_t * ) calloc( 1, sizeof(mesh_elem_t) );
   f->conect = ( mesh_vertex_t ** ) calloc( 3, sizeof(mesh_vertex_t *) );

   /* Put the connectivity inside the finite element auxiliary structure */

   for( j = 0; j < 3; j++ )
     f->conect[j] = conect[i][j];

   /* Update number of vertexs  */

   f->n = 3;
   
   /* Put the finite element auxiliary structure inside the tree node */

   f->next = node->elmhead;
   node->elmhead = f;
 }
}

/* ================== create_quadmesh ===================================== */

static void create_quadmesh( tree_node_t *node )
{
 int     i;

 /* The routine traverse recursivilying the tree and if the tree node  has  no
    child, that is, if  the  tree  node  is  a  leaf  create  a finite element 
    inside it */

 if( node->child != NULL )
   for( i = 0; i < 4; i++ ) create_quadmesh( node->child[i] );
 else
   quadrilaterals_quadmesh( node->elmhead );
}

/* ==================== quadrilaterals_quadmesh =========================== */

static void quadrilaterals_quadmesh( mesh_elem_t *f )
{
 int            i, j, n;
 int            aux1, aux2;
 new_coord      cm;
 mesh_elem_t    *f1, *faux;
 mesh_vertex_t  *quad[4], *mp[5];
      
 /* At this point you have polygons inside the tree node and then they will
    be transformed in finite elements */

 /* Look trough all polygons inside the tree node */

 while( f )
 {  		       
  /* Make a finite element from the polygon, thais, work with the finite
     element connectivity */

  n = f->n;   
  faux = f->next;
  cm.ex = cm.ey = 0;
  for( i = 0; i < n; i++ )
  {
   cm.ex += f->conect[i]->x; 
   cm.ey += f->conect[i]->y;
   aux1 = ( f->conect[((i+1)%f->n)]->x + f->conect[i]->x ) / 2;
   aux2 = ( f->conect[((i+1)%f->n)]->y + f->conect[i]->y ) / 2;
   mp[i] = find_meshv( aux1 , aux2 );
  }
  cm.ex /= f->n;
  cm.ey /= f->n;
  quad[3] = find_meshv( cm.ex , cm.ey );

  /* Get the finite element connectivity and transform the polygon inside
     the tree node, in a finite element */

  for( i = 0; i < n; i++ )
  {					
   quad[0] = mp[i];
   quad[1] = f->conect[((i+1)%f->n)];
   quad[2] = mp[((i+1)%f->n)];

   if( i == n - 1 )
   {
    free( f->conect );
    f->conect = 
     ( mesh_vertex_t ** ) calloc( 4, sizeof(mesh_vertex_t *) );
    f->n = 4;
    f->nel = (nelint++);
    for( j = 0; j < 4; j++ ) f->conect[j] = quad[j];
   }
   else
   { 
    f1 = ( mesh_elem_t * ) calloc( 1, sizeof(mesh_elem_t) );
    f1->conect = 
     ( mesh_vertex_t ** ) calloc( 4, sizeof(mesh_vertex_t *) );
    f1->n = 4;	 
    f1->nel = (nelint++);
    for( j = 0; j < 4; j++ ) f1->conect[j] = quad[j];
    f1->next = f->next;
    f->next = f1;
   }
  } 
   
  /* Get the next polygon in the tree node */
  	       
  f = faux;
 }   	   
}  

/*
** ---------------------------------------------------------------------------
** Local auxiliary functions:
*/

/* ====================== find_meshv ====================================== */

static mesh_vertex_t *find_meshv( int x, int y )
{
 mesh_vertex_t *v, *vertex;

 /* Verify if the vertex already exists. If so, return */

 vertex = vhead;
 while( vertex )
 {
   if( vertex->x == x && vertex->y == y ) 
    return  vertex;
   vertex = vertex->next;
 }

 /* If the vertexs doesn't exist, allocate it and put it into the find_meshv
    quadtree finite element vertexs structure */

 v = ( mesh_vertex_t * ) calloc( 1, sizeof(mesh_vertex_t) );
 v->x       = x;
 v->y       = y;
 v->nno     = nnoint + nnogeomint;
 v->fadj    = NULL;
 v->status  = 1;

 /* Update number of finite elements vertexs in the quadtree */

 nnoint++;
 
/* Update static global quadtree finite element vertexs structure and 
   return find_meshv quadtree finite element vertexs structure*/

 v->next = vhead;
 vhead = v;
 return v;
}

/* ==================== include_mid_node ================================== */

static void include_mid_node( side_t side, int i, tree_node_t *node,
                              int x, int y, mesh_vertex_t **v,
                              int *type, int *mid_node )
{	    
 neighbour_t   *adj;

 /* Look at tree node adjacents*/

 adj = find_adj_node( side, node );

 /* If there's adjacent, update quadtree finite element vertexs structure,
    middle vertexs structure and type of template that will be used */

 if( adj != NULL && adj->next != NULL )
 {
   v[i] = find_meshv( x, y );
   nnointquad = nnoint;
   mid_node[*type] = i;
   (*type)++;
 }

 /* Release adjacent tree node structure */

 SllDelAll( (Sll *)&adj );
}

