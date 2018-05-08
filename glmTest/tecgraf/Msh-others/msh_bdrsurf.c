/*
** --------------------------------------------------------------
**
** msh_bdr2d.c - Routines to implement mesh2D generation based in
**               boundary contraction (advancing front) technique.
**
** ---------------------------------------------------------------
**
** Created:      01-Oct-97      Joaquim B.C. Neto (3D)
** Modify:       01-Mar-98      Antonio C.O. Miranda (2D)
**
** ---------------------------------------------------------------
**
*/


#define BDR_DRAW      0
#define ALG_OLD       0
#define QTREE_INTP    1
#define SMOOTH        0
#define DEBUG_MESSAGE 0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "msh_defsurf.h"
#include "msh_quadsurf.h"
#include "msh_bdrsurf.h"

#if  BDR_DRAW
#ifdef _WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#include "iup.h"
#include "iupgl.h"
//#include "vgliup.h"
//#include "vglint.h"
static int flag_debug = 0;
static double fs = 1.00;
#endif

#if SMOOTH
FILE *arq_smooth;
static int   flag_smooth = 1;
static float peso = 0.5;
static int   opcao = 1;
#endif


/* --------------------------------------------------------------
** Private definitons and data types:
*/

#define Msh2D_INSERT_NO_VALID_NODES     3
#define Msh2D_SMOOTH_STEP               5
#define MSH2D_NODE_QUANTUM              5000
#define MSH2D_POLY_FACE_QUANTUM         1000
#define MSH2D_BDRY_FACE_QUANTUM         1000
#define MSH2D_TEST_FACE_QUANTUM         1000
#define MSH2D_ADJ_ELEM_QUANTUM          1000
#define MSH2D_ADJ_FACE_QUANTUM          1000
#define MSH2D_ADJ_INI_FACE_QUANTUM      1000
#define MSH2D_ELEM_QUANTUM              1000
#define BDY_FACTOR                      0.30
#define SHA_IDEAL                       (4/sqrt(3))
#define MUL_MAX                         1.50
#define REL_FACTOR                      0.5
#define NEW_NODES                       2
#ifndef PI
#define PI                              3.141592653589793324
#endif
#define MAX(a,b)                        (((a)>(b))?(a):(b))
#define MIN(a,b)                        (((a)<(b))?(a):(b))
#define TRUE                            1
#define FALSE                           0

/* --------------------------------------------------------------
** Private variables:
*/

static int                num_elem_alloced = 0 ;
static int                num_node_alloced = 0 ;
static Msh2DBdryEdge      *test_block_ptr = 0 ;
static Msh2DBdryEdge      *test_tail = 0 ;
static Msh2DBdryEdge      *bdry_free = 0 ;
static Msh2DBdryEdge      *bdry_block_ptr = 0 ;
static Msh2DBdryEdge      *bdry_tail = 0 ;
static Msh2DBdryEdge      *bdry_cursor = 0 ;
static Msh2DBdryEdge      *bdry_try = 0 ;
static Msh2DAdjElem       *adj_elem_free = 0 ;
static Msh2DAdjElem       *adj_elem_block_ptr = 0 ;
static Msh2DAdjEdge       *adj_free = 0 ;
static Msh2DAdjEdge       *adj_block_ptr = 0 ;
static Msh2DAdjIniEdge    *adj_ini_free = 0 ;
static Msh2DAdjIniEdge    *adj_ini_block_ptr = 0 ;
static Msh2DBdryNodeList  *node_list  = NULL ;
static Msh2DBdryEdge      *bdry_stack = 0 ;
static int                nintnode  = 0 ;
static int                nbdrynode = 0 ;

static double             msh2d_toler;

MshSurfSizeElement        *mshsurf_reg_size;

/* --------------------------------------------------------------
** Private functions prototypes:
*/

static int             Msh2DSmooth( int, int, int, int*) ;
static int             Msh2DProxEdge( Msh2DBdryEdge *, Msh2DBdryEdge *,
                       Msh2DBdryNodeList *, int);

static int             Msh2DInsideElem (Msh2DBdryEdge *, Msh2DBdryEdge *,
                                        Msh2DBdryNodeList *, int);

static int             Msh2DCrossEdge( Msh2DBdryNodeList *, int, int, int, int);
static double          Msh2DCrossProd( Msh2DBdryNodeList *, int, int, int);
static double          Msh2DAreaTrian( Msh2DBdryEdge *, Msh2DBdryNodeList *, int) ;
static int             Msh2DCheckValid( Msh2DBdryEdge *, Msh2DBdryNodeList *, int,
                       int, int) ;
static int             Msh2DCheckArea( Msh2DBdryEdge *, Msh2DBdryNodeList *, int) ;
static int             Msh2DCheckOptim( Msh2DBdryEdge *, int, int, int,
                       double [][2], int [][2]) ;

static int             Msh2DDetNode( Msh2DBdryEdge *, Msh2DBdryNodeList *, int, int);
static void            Msh2DHeapInit( int ) ;
static void            Msh2DHeapDelete(void) ;
static void            Msh2DHeapInsert( double, int) ;
static int             Msh2DHeapExtract( double *) ;
static void            Msh2DTestFreeAll(void) ;
static Msh2DBdryEdge   *Msh2DEdgeAlloc(void) ;
static void            Msh2DBdryFree(Msh2DBdryEdge  *edge) ;
static void            Msh2DEdgeFreeAll(void) ;
static void            Msh2DBdryPush( Msh2DBdryEdge *) ;
static void            Msh2DBdryPushCorre( Msh2DBdryEdge *) ;
static void            Msh2DBdryPushSmall( Msh2DBdryEdge *) ;
static Msh2DBdryEdge   *Msh2DBdryPop(void) ;
static void            Msh2DBdryDelete(Msh2DBdryEdge  *edge) ;
static Msh2DAdjEdge    *Msh2DAdjEdgeAlloc(void) ;
static void            Msh2DAdjFree(Msh2DAdjEdge  *edge) ;
static void            Msh2DAdjFreeAll(void) ;
static Msh2DAdjIniEdge *Msh2DAdjIniEdgeAlloc(void) ;
static void            Msh2DAdjIniFreeAll(void) ;
static Msh2DAdjElem    *Msh2DAdjElemAlloc(void) ;
static void            Msh2DAdjElemFreeAll(void) ;
static void            Msh2DAdjNodeFreeAll(void) ;
static void            Msh2DAddElem(int *, int **, Msh2DBdryEdge *, int) ;
static void            Msh2DFindNode( Msh2DBdryEdge *, int, double *) ;
static int             Msh2DInsertNewNodesBar (Msh2DBdryEdge *, int *,
                                               Msh2DBdryNodeList *, Msh2DBdryNodeList  **) ;
static int             Msh2DChkScan( double [][2], int, int, int) ;

static void            Msh2DBdryReset(void) ;
static Msh2DBdryEdge   *Msh2DBdryNext(void) ;


/* function to surface mesh */
static double           SurfMetricShape(Msh2DBdryEdge *,
                        Msh2DBdryNodeList *,int ,
                        void (*f_surf)(double,double,double *));
static double           SurfDist3D(Msh2DBdryNodeList *,Msh2DBdryNodeList *);
static Msh2DBdryEdge   *SurfPushBdryEdge( Msh2DBdryNode *,int,int,int,int,
                        void (*f_surf)(double,double,double *));
static void             SurfAddEdges( Msh2DBdryEdge *, Msh2DBdryNodeList *, int,
                        void (*f_surf)(double,double,double *));
static void             SurfTreeNode( Msh2DBdryEdge *, int,
                        void (*f_surf)(double,double,double *),
                        Msh2DBdryNodeList *, double *) ;


/* ----------------------------------------------------------
** Msh2DError - return error information.
*/

typedef void RegSurfFunction (double u,double v, double *f_dev);

static        RegSurfFunction *LocalSurfFunc;
static double Msh2DBdryInit (int num_nodes, int num_edges, double _nodes[][2], int _edges[][2]);
static int    Msh2DBdryGetNode (Msh2DBdryEdge *, int, Msh2DBdryNodeList *, Msh2DBdryNodeList *);
static int    Msh2DBdyGetOptimalNode (Msh2DBdryEdge *, int *, Msh2DBdryNodeList *,
                                      Msh2DBdryNodeList *, int, int, double [][2], int [][2]);



#if  BDR_DRAW

extern void* vglGetCurrentCopy( void );

static void LocalFlush (void)
{
 //GLint buffer;

 //VglContext *ctx = (VglContext *) vglGetCurrentCopy( );

 //glGetIntegerv( GL_DRAW_BUFFER, &buffer );

 //if ( buffer == GL_BACK )
 //   IupGLSwapBuffers( vglIupGetCanvas( ctx ));
 //else
    glFlush();
}
#endif


/**************************************************************/
/* MshSurfRegFunc
*/
void MshSurfBdryRegFunc (MshSurfSizeElement *mshsurf_size)
{
  mshsurf_reg_size = mshsurf_size;
}


/* -------------------------------------------------------------
** Msh2DBdryContraction - main driver for the boundary contraction
**                        algorithm.
*/

int SurfBdryContraction
(
int                  num_org_nodes,
int                  num_org_edges,
double               original_nodes[][2],
int                  original_edges[][2],
int                  *num_int_nodes,
double               **internal_nodes,
int                  *num_gen_elements,
int                  **generated_elements,
void (*f_surf)(double u,double v, double *f_dev)
)
{
  Msh2DBdryEdge     *edge = NULL;
  int               num_nodes, node_indx, i;
  double            smallest_edge = 0.0;

  /* 2.1 build the node list, then push all the boundary edge
     descriptions onto a stack  */
  LocalSurfFunc = f_surf;
  *num_gen_elements = 0 ;
  nintnode  = (*num_int_nodes) ;
  nbdrynode = num_org_nodes ;
  num_nodes = num_org_nodes + (*num_int_nodes);
  smallest_edge  = Msh2DBdryInit (num_org_nodes, num_org_edges, original_nodes, original_edges);
  if (smallest_edge == 0.0)
  {
    printf ("There is a edge with zero length\n");
    Msh2DEdgeFreeAll() ;
    return (0);
  }
  msh2d_toler = smallest_edge / 100.0;



  /* 2.2 use the edge on the top of the stack.  If the stack is
          empty then we are done */
  while ( (edge = Msh2DBdryPop()) )
  {
    Msh2DBdryNodeList _new={{0,0},0,NULL,NULL,NULL,0,0,0};
#if DEBUG_MESSAGE
    printf ("Aresta a ser testada = %d %d - use %d\n", edge->verts[0]+1, edge->verts[1]+1, edge->use);
#endif
    /* Try to get a node from the current boundary */
    node_indx = Msh2DBdryGetNode (edge, num_nodes, node_list, &_new);

    /* Try to insert an optimal node */
    if (node_indx == -1 && edge->use == 0)
    {
      node_indx = Msh2DBdyGetOptimalNode (edge, &num_nodes, node_list, &_new,
                                          num_org_nodes, num_org_edges,
                                          original_nodes, original_edges);
      if (node_indx == -1)
      {
        num_nodes--;
        nintnode--;
      }
    }

    if (node_indx != -1) /* It found a valid node to create a new element */
    {
#if DEBUG_MESSAGE
      printf ("Elemnt = %5d, (%5d %5d %5d), use %d\n", *num_gen_elements+1, edge->verts[0]+1, edge->verts[1]+1, node_indx+1, edge->use);
#endif
      Msh2DAddElem (num_gen_elements, generated_elements, edge, node_indx);
      SurfAddEdges (edge, node_list, node_indx, f_surf) ;
#if BDR_DRAW
{
 LocalFlush ();
 IupMessage("Draw elem", "Pause" );
}
#endif
    }
    else  /* It did not find a valid node => element is not created */
    {
      if(edge->use != 2)
      {
        edge->use++;
        Msh2DBdryPush (edge) ;
        continue ;
      }
    }

  }   /* while ( (edge = Msh2DBdryPop()) ) */

#if 1
   /* 2.6 Smooth the internal nodes, because the boundary nodes can't be changed */
  Msh2DSmooth (num_nodes, nbdrynode, *num_gen_elements, *generated_elements);
#endif

  /* 2.7 Get the generated nodes */
  *num_int_nodes    = (nbdrynode+nintnode) ;
  (*internal_nodes) = (double *) calloc (2*(*num_int_nodes), sizeof(double));
  for (i = 0; i < (*num_int_nodes); i++)
  {
	  (*internal_nodes)[i*2+0] = node_list[i].coord[0];
	  (*internal_nodes)[i*2+1] = node_list[i].coord[1];
  }

  /* 2.8 release memory and return status indicating that a mesh was generated */
  Msh2DHeapDelete() ;
  Msh2DEdgeFreeAll() ;
  Msh2DTestFreeAll() ;
  Msh2DAdjFreeAll() ;
  Msh2DAdjIniFreeAll() ;
  Msh2DAdjElemFreeAll() ;
  Msh2DAdjNodeFreeAll() ;

  return (1) ;
}


/********************** Msh2DBdryInit ******************************/
double Msh2DBdryInit (int num_nodes, int num_edges, double _nodes[][2], int _edges[][2])
{
  Msh2DBdryEdge     *edge = NULL;
  Msh2DBdryNode     *node = NULL;
  Msh2DAdjEdge      *adj_edge = NULL;
  Msh2DAdjIniEdge   *adj_ini_edge = NULL;
  int               i, j, k;
  double            smaller_edge = -1.0;
  double            dev[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  num_elem_alloced  = 0 ;
  bdry_stack = NULL;

  node_list = (Msh2DBdryNodeList *) Msh2DMalloc(num_nodes * sizeof (Msh2DBdryNodeRec));
  num_node_alloced = num_nodes;
  for ( i=0 ; i < num_nodes ; i++ )
  {
    for ( j=0 ; j<2 ; j++ )
     node_list[i].coord[j] = _nodes[i][j] ;
#if QTREE_INTP
    MshSurfQuadTreeDeriv(node_list[i].coord[0],node_list[i].coord[1],dev,0);
#else
    LocalSurfFunc (node_list[i].coord[0],node_list[i].coord[1],dev);
#endif
    node_list[i].E = dev[0]*dev[0] + dev[1]*dev[1] + dev[2]*dev[2];
    node_list[i].F = dev[0]*dev[3] + dev[1]*dev[4] + dev[2]*dev[5];
    node_list[i].G = dev[3]*dev[3] + dev[4]*dev[4] + dev[5]*dev[5];
    node_list[i].active_flag = 1 ;
    node_list[i].edges  = NULL;
    node_list[i].iedges = NULL;
  }

  for (i = 0; i < num_edges; i++)
  {
    edge = SurfPushBdryEdge (node_list, _edges[i][0], _edges[i][1], 0, 0, LocalSurfFunc);
    if ((smaller_edge == -1.0) || (smaller_edge > edge->length))
      smaller_edge = edge->length;
  }

  /* build the list of adjacent edges for all the boundary nodes */
  Msh2DBdryReset();
  while ( (edge = Msh2DBdryNext()) )
  {
    for ( j=0; j<2; j++)
    {
      node = &node_list[edge->verts[j]];
      adj_edge = Msh2DAdjEdgeAlloc();
      adj_edge->next = node->edges;
      adj_edge->edge = edge;
      node->edges = adj_edge;
    }
  }

  /* build the list of initial adjacent edges for all the boundary nodes */
  Msh2DBdryReset ( );
  while ( (edge = Msh2DBdryNext ()) )
  {
    for ( j=0 ; j<2 ; j++ )
    {
      node = &node_list[edge->verts[j]];
      adj_ini_edge = Msh2DAdjIniEdgeAlloc ( );
      adj_ini_edge->next = node->iedges;
      for ( k=0; k<2; k++ )
      {
        adj_ini_edge->verts[k] = edge->verts[k] ;
        adj_ini_edge->r[k] = edge->r[k];
      }
      node->iedges = adj_ini_edge ;
    }
  }

  return smaller_edge;
}


/********************** Msh2DBdryGetNode ******************************/
static int Msh2DBdryGetNode (Msh2DBdryEdge *curr_edge, int n_nodes,
                             Msh2DBdryNodeList *nd_list, Msh2DBdryNodeList *_new)
{
  double h = 0.0, dist, metric, cand_vec[2], dot, next_metric;
  int    node_indx = -1, next_indx, check, i, j;

  /* examine all active nodes. Rank all nodes as to the goodness
     of the triangle they will form with the current edge */
  if (curr_edge->use == 0)
    SurfTreeNode (curr_edge, 1, LocalSurfFunc, _new, &h);
  dist = h * 0.85;

  Msh2DHeapInit (n_nodes) ;
  for (i = 0; i < n_nodes; i++)
  {
    if (nd_list[i].active_flag )
    {
      double d;

      /* find the vector from the center to the new point
          and make sure that this cross with the normal is
          positive */
      cand_vec[0] = nd_list[i].coord[0] - curr_edge->center[0] ;
      cand_vec[1] = nd_list[i].coord[1] - curr_edge->center[1] ;

      dot = cand_vec[0] * curr_edge->nrm[0] + cand_vec[1] * curr_edge->nrm[1];

      if ( dot <= 0.0 ) continue ;

      /* verify if node is out of sphere centered in optimal node
          for the edge and ratio equal of largest (smallest)
          distance from this node to any vertex of the edge.
          This should be doen only the fase of ideal elements. */
      if (curr_edge->use == 0)
      {
        d = SurfDist3D(&nd_list[i], _new);;
        if ( d >= dist ) continue ;
      }

      /* the metric we are currently using is the square of
          the distance from the center of the edge to the
          candidate node */

      if (curr_edge->use == 0)
      {
        metric = SurfMetricShape( curr_edge, nd_list, i, LocalSurfFunc ) ;
        if( metric > 0.0 )
          Msh2DHeapInsert( (1.0*metric), i ) ;
      }
      else
      {
        metric = 0.0 ;
        for ( j=0 ; j<2 ; j++ )
           metric += cand_vec[j] * cand_vec[j] ;

        metric = metric / dot ;
        // Msh2DHeapInsert( metric, i ) ;
        if( metric > 0.0 )
          Msh2DHeapInsert( (1.0*metric), i ) ;
      }
    }
  } /* end for - calcula metrica de escolha */

/* 2.4 start with the node with the best ranking.  Check to make
      sure it will form a valid element (i.e., does not intersect
      the current boundary).  If the element is invalid go on
      to the next.  If the element is valid, then look at the next
      candidate node.  If the two nodes have the same coordinates
      then use topology to determine the one we want. */

  while (1)
  {
    /* extract a node from the heap based in its metric */
    node_indx = Msh2DHeapExtract (&metric) ;

    if ((node_indx == curr_edge->verts[0]) || (node_indx == curr_edge->verts[1]))
      continue ;

    /* here a node was found to make the element (node_indx was
        extracted from the heap). Check its validity */
    if (node_indx != -1)
    {
      check = 1;

      /* check validity for node choosen */
      if (curr_edge->use == 0 && check)
        check = Msh2DCheckArea (curr_edge, nd_list, node_indx);
      if (check)
        check = Msh2DCheckValid (curr_edge, nd_list, node_indx, -1, 1);

      /* set validity for node choosen */
      if (check)
      {
        /* check if there is another node with equal metric - crack */
        next_indx = Msh2DHeapExtract (&next_metric);
        if (next_indx != -1)
        {
          if (fabs (next_metric - metric) < msh2d_toler &&
              fabs (nd_list[node_indx].coord[0] - nd_list[next_indx].coord[0]) < msh2d_toler &&
              fabs (nd_list[node_indx].coord[1] - nd_list[next_indx].coord[1]) < msh2d_toler )
          {
            node_indx = Msh2DDetNode (curr_edge, nd_list, node_indx, next_indx);
          }
        }
        return node_indx; /* it returns the valid node */
      }
    }
    else
      return -1;
  }

  return -1; /* it was not possible to find a valid node */
}

/********************** Msh2DBdyGetOptimalNode  ******************************/
static int Msh2DBdyGetOptimalNode (Msh2DBdryEdge *curr_edge, int *n_nodes,
                                   Msh2DBdryNodeList *nd_list, Msh2DBdryNodeList *_new,
                                   int num_org_nodes, int num_org_edges,
                                   double original_nodes[][2], int original_edges[][2])
{
  int node_indx = -1;
  int check;

  /* try insertion of new node */
  check = Msh2DInsertNewNodesBar (curr_edge, n_nodes, _new, &nd_list);

  /* check validity for node inserted */
  if (check)
    check = Msh2DCheckOptim (curr_edge, *n_nodes-1, num_org_nodes,
                             num_org_edges, original_nodes, original_edges) ;
  if (check)
    check = Msh2DCheckValid (curr_edge, nd_list, *n_nodes-1, -1, -1);

  /* set validity for node inserted */
  if (check)
  {
    node_indx = *n_nodes - 1 ;
    return node_indx;
  }

  return -1;
}

/* -------------------------------------------------------------------
** Msh2DSmooth - driver for nodal smoothing.
*/

static int Msh2DSmooth
(
int           num_nodes ,
int           bdr_nodes ,
int           num_elems ,
int           *elements
)
{
    Msh2DBdryNodeList *smoo_list ;
    Msh2DBdryNode     *node ;
    Msh2DAdjElem      *adj_elem ;
    Msh2DAdjElemList  *elems ;
    int               i, j, k, n, id ;
    double             x, y, z, w, wj ;
    double            dist_par, dist_3d;

/* 3.1 Allocate memory for smooth nodes vector */


    smoo_list = (Msh2DBdryNodeList *)Msh2DMalloc( num_nodes *
         sizeof(Msh2DBdryNodeRec) ) ;

/* 3.2 Initiate list of all elements adjacent to a given node and
       initiate list of smmoth nodes vector */

    for( i = 0; i < num_nodes; i++ )
    {
     node_list[i].elems = (Msh2DAdjElemList *)0 ;
     smoo_list[i].coord[0] = node_list[i].coord[0] ;
     smoo_list[i].coord[1] = node_list[i].coord[1] ;
    }

/* 3.3 Build a list of all elements adjacent to a given node */

    for( i = 0; i < num_elems; i++ ) {
     if( elements[i*3] == -1 ) continue ;
     for( j = 0; j < 3; j++ ) {
      node = &node_list[elements[(i*3)+j]] ;
      adj_elem = Msh2DAdjElemAlloc( ) ;
      adj_elem->next = node->elems ;
      adj_elem->elem = i ;
      node->elems = adj_elem ;
     }
    }

/* 3.4 Move each internal node to the centroid of all its adjacent
       nodes defined by its adjacent elements. Each internal node
       only will be moved if doesn't affect the consistency of its
       adjacent elements, so consistency tests are necessary after
       each move. */

    for( i = 0; i < Msh2D_SMOOTH_STEP; i++ )
    {
     for( j = 0; j < num_nodes; j++ )
     {
      if( j >= bdr_nodes )
      {
       n = 0 ;
       x = y = z = w = 0.0 ;
       elems = node_list[j].elems ;
       while( elems )
       {
         for( k = 0; k < 3; k++ )
         {
           id = elements[(elems->elem*3)+k] ;
           if( id != j )
           {
              dist_par = sqrt(((node_list[id].coord[0] - node_list[j].coord[0])  *
                               (node_list[id].coord[0] - node_list[j].coord[0])) +
                              ((node_list[id].coord[1] - node_list[j].coord[1])  *
                               (node_list[id].coord[1] - node_list[j].coord[1])) );
              dist_3d = SurfDist3D(&node_list[id],&node_list[j]);
              wj = (dist_3d / dist_par);;
              x += (wj * (node_list[id].coord[0] - node_list[j].coord[0])) ;
              y += (wj * (node_list[id].coord[1] - node_list[j].coord[1])) ;
              w += wj ;
              n++ ;
           }
         }
         elems = elems->next ;
       }
       elems = node_list[j].elems ;
       if( elems )
       {

        if( w > 0.0 )
         {
           smoo_list[j].coord[0] = node_list[j].coord[0] + (REL_FACTOR * (x / w)) ;
           smoo_list[j].coord[1] = node_list[j].coord[1] + (REL_FACTOR * (y / w)) ;
         }

       }
      }
     }
     for( j = 0; j < num_nodes; j++ )
     {
       node_list[j].coord[0] = smoo_list[j].coord[0] ;
       node_list[j].coord[1] = smoo_list[j].coord[1] ;
     }
    }

/* 3.5 Release memory for smooth nodes vector */

    Msh2DFree( smoo_list ) ;

/* 3.6 Return smooth status */

    return(1) ;
}

/* -------------------------------------------------------------------
** SurfPushBdryEdge - this routine pushes boundary edges onto the stack.
*/

static Msh2DBdryEdge *SurfPushBdryEdge
(
Msh2DBdryNode   *nodes ,
int             v1 ,
int             v2 ,
int             use,
int             init,
void (*f_surf)(double u,double v,double *f_dev)
)
{
    Msh2DBdryEdge   *edge ;
    register int    i ;

    edge = Msh2DEdgeAlloc() ;

    edge->verts[0] = v1 ; edge->verts[1] = v2 ;

    edge->use = use ;

    /* the center of the edge is the algebraic mean of the corners */

    edge->center[0] = (nodes[edge->verts[0]].coord[0] +
                       nodes[edge->verts[1]].coord[0])/2;
    edge->center[1] = (nodes[edge->verts[0]].coord[1] +
                       nodes[edge->verts[1]].coord[1])/2;

    /*  find the max and mins of the coordinates */

    for ( i=0 ; i<2 ; i++ )
    {
      edge->max[i] = nodes[v1].coord[i] ;
      edge->min[i] = nodes[v1].coord[i] ;
      if (edge->max[i] < nodes[v2].coord[i]) edge->max[i] = nodes[v2].coord[i] ;
      if (edge->min[i] > nodes[v2].coord[i]) edge->min[i] = nodes[v2].coord[i] ;
    }

    /* compute the r vector  */
    edge->r[0] = nodes[v2].coord[0] - nodes[v1].coord[0] ;
    edge->r[1] = nodes[v2].coord[1] - nodes[v1].coord[1] ;

    edge->length = sqrt((edge->r[0]*edge->r[0]) + (edge->r[1]*edge->r[1]));
    edge->r[0] = edge->r[0] / edge->length;
    edge->r[1] = edge->r[1] / edge->length;
    edge->nrm[0] =  edge->r[1];
    edge->nrm[1] = -edge->r[0];
    edge->real_length = SurfDist3D(&nodes[v1],&nodes[v2]);

    if (edge->length == 0.0)
    {
#if DEBUG_MESSAGE
      printf ("Edge (%d - %d) with zero length\n", edge->verts[0], edge->verts[1]);
#endif
    }

#if BDR_DRAW
{
 double nr[3], norm=edge->length/5.0;
 glColor3d (1.0, 0.0, 0.0);
 glBegin( GL_LINES );
   glVertex3d (nodes[v1].coord[0]*fs, nodes[v1].coord[1]*fs,0.0);
   glVertex3d (nodes[v2].coord[0]*fs, nodes[v2].coord[1]*fs,0.0);
 glEnd( );

 nr[0] = edge->center[0]+edge->nrm[0]*norm;
 nr[1] = edge->center[1]+edge->nrm[1]*norm;
 glColor3d (0.0, 1.0, 1.0);
 glBegin( GL_LINES );
   glVertex3d (edge->center[0]*fs, edge->center[1]*fs, 0.0);
   glVertex3d (nr[0]*fs,nr[1]*fs,0.0);
 glEnd( );
}
#endif

    if (init == 0)
      Msh2DBdryPushSmall( edge ) ;
    else
      Msh2DBdryPushCorre( edge ) ;

    return(edge) ;
}


/*************************** Msh2DInsertNewNodesBar ************************/
static int Msh2DInsertNewNodesBar
(
Msh2DBdryEdge       *edge ,
int                 *n ,
Msh2DBdryNodeList   *_new,
Msh2DBdryNodeList   **nd_list
)
{
 int               j ;
 int               num_node ;
 double            dot ;
 double            node[2] ;
/* double            _new[2] ; */
 double            cand_vec[2] ;

 /* initiate number of nodes and allocation */
 num_node = *n;

 /* add new coordinate for baricent of base edge */
 for (j = 0; j < 2; j++)
  node[j] = _new->coord[j];

 /* update node_list */
 num_node++;
 nintnode++;
 if (num_node > num_node_alloced)
 {
   num_node_alloced += MSH2D_NODE_QUANTUM;
   node_list = (Msh2DBdryNodeList *) Msh2DRealloc (node_list,
               num_node_alloced * sizeof (Msh2DBdryNodeRec));
 }
 for (j = 0; j < 2 ; j++)
  node_list[num_node-1].coord[j] = node[j];
 node_list[num_node-1].active_flag = 1;
 node_list[num_node-1].edges  = NULL;
 node_list[num_node-1].iedges = NULL;
 node_list[num_node-1].E = _new->E;
 node_list[num_node-1].F = _new->F;
 node_list[num_node-1].G = _new->G;

 *nd_list = node_list;

 /* update number of nodes and allocation */
 *n = num_node;

 /* the new point should only be considered if it's in the same semi-plane
    than the base edge, there is, in the direction of base edge's normal */

 cand_vec[0] = node[0] - edge->center[0] ;
 cand_vec[1] = node[1] - edge->center[1] ;

 dot = cand_vec[0]*edge->nrm[0] +
       cand_vec[1]*edge->nrm[1] ;

 if ( dot <= 0.0 ) return 0 ;
 else              return 1 ;
}


typedef struct _Msh2Dedge {
 int    flag;
 int    vnodes[2];
} Msh2Dedge ;



/*************************** Msh2DFindNode *********************************/
static void Msh2DFindNode
(
Msh2DBdryEdge       *edge ,
int                 debug ,
double              node[2]
)
{
 double             h ;

 /* get the polygon equilateral size of the same length of the edge */

 h = (edge->length * sqrt(3.0)) / 2.0;

 /* evaluate the new node coordinate */

 node[0] = edge->center[0] + (h/debug)*edge->nrm[0];
 node[1] = edge->center[1] + (h/debug)*edge->nrm[1];
}

/*************************** SurfTreeNode *********************************/
static void SurfTreeNode
(
Msh2DBdryEdge       *edge ,
int                 debug ,
void (*f_surf)(double u,double v,double *f_dev),
Msh2DBdryNodeList   *node,
double              *h0
)
{
 int                level ;
 double             h, h1 ;
 double             _new[2] ;
 double             dev[6];
 double             h3D, h2D, N2D[2], Nc2D[2], E, F, G, Vab[2];


 /* get the size of the cell where the center of the edge is */

 h = 0.0 ;
/* h += SurfOptimalNodes ( node_list[edge->verts[0]].coord, _new, &level ) ;
 h += SurfOptimalNodes ( node_list[edge->verts[1]].coord, _new, &level ) ;
 h /= 2.0 ; */

 if (mshsurf_reg_size != NULL)
 {
   (mshsurf_reg_size)(edge->center[0], edge->center[1], &h1);
   h = h1 * 0.7;
 }
 else
   h = MshSurfOptimalNodes ( edge->center, _new, &level ) ;

 if ( h > (edge->real_length*1.5) ) h = edge->real_length*1.5;

 /* evaluate the new coordinate using Riemannian metric */
 /* Advancing Front Surface Mesh Generation in Parametric Space
    Using a Riemanninan Surface Definition, 7th International
    Meshing Roundtable Proceeding, pp 429-445 */
#if 1
 h3D = h;
#else
 h3D = sqrt((edge->real_length*edge->real_length) -
            (edge->real_length*edge->real_length)/4.0);
#endif
/*
printf("H's -> Qtree = %f, Edge = %f, Real Len = %f \n",h,h3D,edge->real_length);
*/
 E = (node_list[edge->verts[0]].E + node_list[edge->verts[1]].E)/2.0;
 F = (node_list[edge->verts[0]].F + node_list[edge->verts[1]].F)/2.0;
 G = (node_list[edge->verts[0]].G + node_list[edge->verts[1]].G)/2.0;

 Vab[0] = edge->r[0] / edge->length;
 Vab[1] = edge->r[1] / edge->length;
 N2D[0] =  F*Vab[0] + G*Vab[1];
 N2D[1] = -E*Vab[0] - F*Vab[1];

 h2D = h3D/sqrt(E*N2D[0]*N2D[0] + 2*F*N2D[0]*N2D[1] + G*N2D[1]*N2D[1]);

 Nc2D[0] = edge->center[0] + N2D[0]*h2D;
 Nc2D[1] = edge->center[1] + N2D[1]*h2D;

 /* evaluate the new node coordinate */

 node->coord[0] = Nc2D[0];
 node->coord[1] = Nc2D[1];
#if QTREE_INTP
 MshSurfQuadTreeDeriv(Nc2D[0],Nc2D[1],dev,0);
/*
printf("Deriv 01 = %f %f %f and %f %f %f\n",dev[0],dev[1],dev[2],dev[3],dev[4],dev[5]);
(*f_surf)(Nc2D[0],Nc2D[1],dev);
 printf("Deriv 02 = %f %f %f and %f %f %f\n\n",dev[0],dev[1],dev[2],dev[3],dev[4],dev[5]);
*/
#else
 (*f_surf)(Nc2D[0],Nc2D[1],dev);
#endif

 node->E = dev[0]*dev[0] + dev[1]*dev[1] + dev[2]*dev[2];
 node->F = dev[0]*dev[3] + dev[1]*dev[4] + dev[2]*dev[5];
 node->G = dev[3]*dev[3] + dev[4]*dev[4] + dev[5]*dev[5];

 *h0 = h3D;
}


/* -------------------------------------------------------------------
** TemporaryEdge - this routine set temporary edge
*/
static void TemporaryEdge
(
double         x1, double y1,
double         x2, double y2,
Msh2DBdryEdge *edge
)
{

    /* the center of the edge is the algebraic mean of the corners */
    edge->center[0] = (x1 + x2)/2;
    edge->center[1] = (y1 + y2)/2;

    /* compute the r vector  */
    edge->r[0] = x2 - x1;
    edge->r[1] = y2 - y1;
    edge->length = sqrt((edge->r[0]*edge->r[0]) + (edge->r[1]*edge->r[1]));
    edge->r[0] = edge->r[0] / edge->length;
    edge->r[1] = edge->r[0] / edge->length;
    edge->nrm[0] = edge->r[1];
    edge->nrm[1] = -edge->r[0];
}


#if 0
/*************************** Msh2DProxTmpEdge *********************************/
static int Msh2DProxTmpEdge
(
Msh2DBdryEdge      *tmpedge1,
Msh2DBdryEdge      *tmpedge2,
double              x,
double              y,
double              div
)
{
 double   dot1, dot2, v[2];

 /* vetor partindo do centro */
 v[0] = x - tmpedge1->center[0];
 v[1] = y - tmpedge1->center[1];

 dot1 = v[0]*tmpedge1->nrm[0] + v[1]*tmpedge1->nrm[1]; /* direcao normal a aresta */
 dot2 = v[0]*tmpedge1->r[0]   + v[1]*tmpedge1->r[1];   /* direcao paralela a aresta */

 if ((dot1 < tmpedge1->length*BDY_FACTOR/div) &&  (dot1 > 0.0) &&
     (fabs(dot2) < tmpedge1->length*0.4))
   return (0);
/*
 if ((dot1 < tmpedge1->length*BDY_FACTOR/10.0) && (dot1 > 0.0) &&
     (fabs(dot2) < tmpedge1->length*0.51))
   return (0);
*/

 /* segunda areasta */
 v[0] = x - tmpedge2->center[0];
 v[1] = y - tmpedge2->center[1];

 dot1 = v[0]*tmpedge2->nrm[0] + v[1]*tmpedge2->nrm[1]; /* direcao normal a aresta */
 dot2 = v[0]*tmpedge2->r[0]   + v[1]*tmpedge2->r[1];   /* direcao paralela a aresta */

 if ((dot1 < tmpedge2->length*BDY_FACTOR/div) && (dot1 > 0.0) &&
     (fabs(dot2) < tmpedge2->length*0.4))
   return (0);
/*
 if ((dot1 < tmpedge2->length*BDY_FACTOR/10.0) && (dot1 > 0.0) &&
     (fabs(dot2) < tmpedge2->length*0.51))
   return (0);
*/
 return 1;
}
#endif


/* -------------------------------------------------------------------
** Msh2DAdjEdgsFormTriangle - check if the adjacent current front is a triangle.
*/
static int Msh2DAdjEdgsFormTriangle
(
Msh2DBdryEdge  *edge,
int             node_indx
)
{
  Msh2DAdjEdgeList *adj_edge;
  int               found_adj[2] = {0, 0};

  adj_edge = node_list[node_indx].edges;

  if (adj_edge == NULL)
    return 0;

  do
  {
    /* try to found edge->verts[0] */
    if (adj_edge->edge->verts[0] == edge->verts[0] ||
        adj_edge->edge->verts[1] == edge->verts[0] )
     found_adj[0] = 1;

    /* try to found edge->verts[1] */
    if (adj_edge->edge->verts[0] == edge->verts[1] ||
        adj_edge->edge->verts[1] == edge->verts[1] )
     found_adj[1] = 1;

    adj_edge = adj_edge->next;

  } while (adj_edge != NULL);

  return (found_adj[0] & found_adj[1]);
}



/* -------------------------------------------------------------------
** Msh2DCheckValid - this routine makes geometrical checks to make sure
**                   that no edges of a candidate element cross the
**                   the current boundary.
*/
static int  Msh2DCheckValid
(
Msh2DBdryEdge      *edge,
Msh2DBdryNodeList  *nd_list,
int                node_indx,
int                node_bdry,
int                found       /* it found one point on the boundary */
)
{
    int                valid, type, id_edge;
    Msh2DBdryEdge      *current=NULL ;
    double             x1,x2,x3,x4,y1,y2,y3,y4, area;
    double             zero_area = (edge->length*edge->length)/500.0/(edge->use+1);
    Msh2DBdryEdge      tmp_edge1, tmp_edge2;

  /* first of all, check if the adjacent current front is a triangle */
  if (found)
  {
    if (Msh2DAdjEdgsFormTriangle (edge, node_indx))
      return 1;
  }


   area = Msh2DAreaTrian (edge, nd_list, node_indx);
   if (area < zero_area)
     return 0;

    /* loop through all the edges in the current boundary and see
       if they intersect any of the edges of the proposed new element */
    valid = 1 ;
    Msh2DBdryReset() ;

    x1 = nd_list[edge->verts[0]].coord[0];
    y1 = nd_list[edge->verts[0]].coord[1];
    x2 = nd_list[edge->verts[1]].coord[0];
    y2 = nd_list[edge->verts[1]].coord[1];

    /* Set temporary edges */
    TemporaryEdge (x1, y1, nd_list[node_indx].coord[0],
                   nd_list[node_indx].coord[1], &tmp_edge1);
    TemporaryEdge (nd_list[node_indx].coord[0],
                   nd_list[node_indx].coord[1], x2, y2, &tmp_edge2);

    while ( valid && (current = Msh2DBdryNext()) )
    {
       /* test if the current edge is the same of given edge, because
          if so no intersection test is necessary. */
      if ( current == edge )
        continue ;

      x3 = nd_list[current->verts[0]].coord[0];
      y3 = nd_list[current->verts[0]].coord[1];
      x4 = nd_list[current->verts[1]].coord[0];
      y4 = nd_list[current->verts[1]].coord[1];
#if 0
      if (/*(current->verts[0] == edge->verts[0] && current->verts[1] == node_indx) || */
          (current->verts[0] == edge->verts[1] && current->verts[1] == node_indx) ||
           (current->verts[1] == edge->verts[0] && current->verts[0] == node_indx)/* ||
          (current->verts[1] == edge->verts[1] && current->verts[0] == node_indx)  */)
      {
        printf ("Nao precisa testar\n");
        continue;
      }
#endif

      /* test the intersection between the boundary edge we are
        currently checking with the two new edges. */

      /* check the tolerance between vertices, because there are
         some matters with cracks */
      /* teste the first vertice of base edge */
      if((fabs(x1-x3)<msh2d_toler) && (fabs(y1-y3)<msh2d_toler))
         id_edge = current->verts[0];
      else if ((fabs(x1-x4)<msh2d_toler) && (fabs(y1-y4)<msh2d_toler))
         id_edge = current->verts[1];
      else
         id_edge = edge->verts[0];
      type = !((Msh2DCrossEdge(nd_list, id_edge, node_indx,
                current->verts[0], current->verts[1])));
      if( type == 0 )
      {
#if DEBUG_MESSAGE
        printf ("Cross Edge %d %d\n", current->verts[0], current->verts[1]);
#endif
        return 0;
      }

      /* teste the second vertice of base edge */
      if((fabs(x2-x3)<msh2d_toler) && (fabs(y2-y3)<msh2d_toler))
         id_edge = current->verts[0];
      else if ((fabs(x2-x4)<msh2d_toler) && (fabs(y2-y4)<msh2d_toler))
         id_edge = current->verts[1];
      else
         id_edge = edge->verts[1];
      type = !((Msh2DCrossEdge(nd_list, id_edge, node_indx,
                current->verts[0], current->verts[1])));
      if( type == 0 )
      {
#if DEBUG_MESSAGE
        printf ("Cross Edge %d %d\n", current->verts[0], current->verts[1]);
#endif
        return 0;
      }

      if (found != 1)  /* the point is not on the boundary => new point */
      {
        /* test to see if new point is near from current edge. If it's
           then this point should not be considered */
        type = Msh2DProxEdge (current, edge, nd_list, node_indx);
        if( type == 0 )
        {
#if DEBUG_MESSAGE
          printf ("Msh2DProxEdge %d %d\n", current->verts[0], current->verts[1]);
#endif
          return 0;
        }
      }

#if 0
      if (edge->use != 2)
      {
        /* test to see if new edges are near from current points. If they are
           then the point should not be considered */
        type = Msh2DProxTmpEdge (&tmp_edge1, &tmp_edge2, x3, y3, 1.0);
        if ( type == 0 )
        {
#if DEBUG_MESSAGE
          printf ("Msh2DProxTmpEdge %d %d\n", current->verts[0], current->verts[1]);
#endif
          return 0;
        }
        type = Msh2DProxTmpEdge (&tmp_edge1, &tmp_edge2, x4, y4, 1.0);
        if ( type == 0 )
        {
#if DEBUG_MESSAGE
          printf ("Msh2DProxTmpEdge %d %d\n", current->verts[0], current->verts[1]);
#endif
          return 0;
        }
      }
#endif

      type = Msh2DInsideElem (current, edge, nd_list, node_indx);
      if ( type == 0 )
      {
#if DEBUG_MESSAGE
        printf ("Msh2DInsideElem %d %d\n", current->verts[0], current->verts[1]);
#endif
        return 0;
      }

    }

    /* if we get here, the new element is valid */
    return(valid) ;
}


/*************************** Msh2DProxEdge *********************************/
static int Msh2DProxEdge
(
Msh2DBdryEdge      *edge,
Msh2DBdryEdge      *bedge,
Msh2DBdryNodeList  *nd_list,
int                node_indx
)
{
 double   dot1, dot2, v[2];

 /* vetor partindo do centro */
 v[0] = nd_list[node_indx].coord[0] - edge->center[0];
 v[1] = nd_list[node_indx].coord[1] - edge->center[1];

 dot1 = v[0]*edge->nrm[0] + v[1]*edge->nrm[1]; /* direcao normal a aresta */
 dot2 = v[0]*edge->r[0]   + v[1]*edge->r[1];   /* direcao paralela a aresta */

/* if ((dot1 < edge->length*BDY_FACTOR) && (dot1 > 0.0) &&
     (fabs(dot2) < edge->length*0.4))
   return (0); */

 if ((dot1 < edge->length*BDY_FACTOR/5.0) &&  (dot1 >= 0.0) &&
     (fabs(dot2) < edge->length*0.5))
   return (0);

 if ((fabs(dot1) < edge->length*0.1/5.0) &&  (dot1 <= 0.0) &&
     (fabs(dot2) < edge->length*0.49))
   return (0);

 return 1;
}


/*/ CrossProd - Compute the Cross Prod.
///////////////////////////////////////////////////////////////////// */
static double CrossProd(double *Pi, double *Pj, double *Pk)
{
   double    kx, ky, jx, jy;

   /* Orientation:

   | J      / K
   |       /
   |      /
   |     /
   |    /
   |   /
   |  /
   | /
   |/ I
   *
 */

   kx = Pk[0] - Pi[0];
   ky = Pk[1] - Pi[1];
   jx = Pj[0] - Pi[0];
   jy = Pj[1] - Pi[1];

   return (kx*jy - ky*jx);
}

/*************************** Msh2DInisideElem ******************************/
static int Msh2DInsideElem
(
Msh2DBdryEdge      *bedge,
Msh2DBdryEdge      *edge,
Msh2DBdryNodeList  *nd_list,
int                node_indx
)
{
  double P[3][2] ={ {nd_list[edge->verts[1]].coord[0], nd_list[edge->verts[1]].coord[1]},
                    {nd_list[edge->verts[0]].coord[0], nd_list[edge->verts[0]].coord[1]},
                    {nd_list[node_indx].coord[0],      nd_list[node_indx].coord[1]} };
  double Pj[2] =    {nd_list[bedge->verts[0]].coord[0],nd_list[bedge->verts[0]].coord[1]};
  double *Pi, *Pk;
  int    i, count = 0;

   /* Orientation:
   *---------*
   | P[1]   / P[0]
   |       /
   |      /
   | *   /
   | PK /
   |   /
   |  /
   | /
   |/ P[2] - node indx
   *
 */

  if ((fabs (P[0][0] - Pj[0])) < msh2d_toler &&
      (fabs (P[0][1] - Pj[1])) < msh2d_toler )
    return 1;

  if ((fabs (P[1][0] - Pj[0])) < msh2d_toler &&
      (fabs (P[1][1] - Pj[1])) < msh2d_toler )
    return 1;

  if ((fabs (P[2][0] - Pj[0])) < msh2d_toler &&
      (fabs (P[2][1] - Pj[1])) < msh2d_toler )
    return 1;

  for (i = 0; i < 3; i++)
  {
    Pi = P[i];
    Pk = P[(i+1)%3];
    if (CrossProd (Pi, Pj, Pk) > -0.01*msh2d_toler)
      count++;
  }

  if (count == 3) return 0;

  return 1;
}




/* -------------------------------------------------------------------
** Msh2DDetNode  - this routine determines the node that should be
**                 taken in case of two nodes choosen has the same
**                 coordinates.
*/

static int Msh2DDetNode
(
Msh2DBdryEdge      *edge,
Msh2DBdryNodeList  *nd_list,
int                node_indx,
int                next_indx
)
{
  int             num_adj;
  double          vector[2];
  Msh2DAdjIniEdge *cur ;
  double          med_nrm[2] = {0.0, 0.0}, dot1, dot2;

  num_adj = 0;
  for( cur = nd_list[node_indx].iedges ; cur ; cur = cur->next )
  {
    med_nrm[0] += cur->r[1];
    med_nrm[1] += -cur->r[0];
    num_adj++;
  }
  med_nrm[0] /= num_adj;
  med_nrm[1] /= num_adj;
  vector[0] = edge->center[0] - nd_list[node_indx].coord[0];
  vector[1] = edge->center[1] - nd_list[node_indx].coord[1];
  dot1 = vector[0]*med_nrm[0] + vector[1]*med_nrm[1];


  num_adj = 0;
  for( cur = nd_list[next_indx].iedges ; cur ; cur = cur->next )
  {
    med_nrm[0] += cur->r[1];
    med_nrm[1] += -cur->r[0];
    num_adj++;
  }
  med_nrm[0] /= num_adj;
  med_nrm[1] /= num_adj;
  vector[0] = edge->center[0] - nd_list[next_indx].coord[0];
  vector[1] = edge->center[1] - nd_list[next_indx].coord[1];
  dot2 = vector[0]*med_nrm[0] + vector[1]*med_nrm[1];

  if (dot1 > dot2)
    return node_indx;
  else
    return next_indx;

  return -1 ;
}

/* -------------------------------------------------------------------
** SurfAddEdges - this routine updates the edge stack to include/
**                 delete the necessary boundary edges to include
**                 the new element.  It also updates the active
**                 flags in the node list.
*/

static void SurfAddEdges
(
Msh2DBdryEdge      *edge,
Msh2DBdryNodeList  *nd_list,
int                 node_indx,
void (*f_surf)(double u,double v,double *f_dev)
)
{
    int             v1i, v1, i ;
    int             found ;
    Msh2DAdjEdge   *cur, **temp, *save ;
    Msh2DBdryEdge  *edge_to_delete, *_new ;

    /* start with each of the base nodes.  See if they are adjacent
       to a edge that includes the cap node along with the previous
       vertex on the base.  If so, this edge is already part of the
       boundary, and it should be deleted. */

    for ( v1i=0 ; v1i<2 ; v1i++ )
    {
      v1 = edge->verts[v1i] ; found = 0;
      for ( cur=nd_list[v1].edges ; cur ; cur=cur->next )
      {
        if ( (( cur->edge->verts[0] == v1 ) &&
              ( cur->edge->verts[1] == node_indx ) && (v1i == 1)) ||
             (( cur->edge->verts[1] == v1 ) &&
              ( cur->edge->verts[0] == node_indx ) && (v1i == 0)) )
        {
          found = 1 ;
          break ;
        }
      }

       /* if we have found a edge, first we update all the
          adjacent edge lists so they no longer point to it.
          Then we delete the edge from the edge list */

      if ( found )
      {
        edge_to_delete = cur->edge ;
        for ( temp = &(nd_list[v1].edges) ; *temp ;
              temp = &((*temp)->next) )
        {
          if ( (*temp)->edge == edge_to_delete )
          {
            save = *temp ;
            *temp = (*temp)->next ;
            Msh2DAdjFree( save ) ;
            break ;
          }
        }
        for ( temp = &(nd_list[node_indx].edges) ; *temp ;
              temp = &((*temp)->next) )
        {
          if ( (*temp)->edge == edge_to_delete )
          {
            save = *temp ;
            *temp = (*temp)->next ;
            Msh2DAdjFree( save ) ;
            break ;
          }
        }

        Msh2DBdryDelete( cur->edge ) ;
      }

 /* if we did not find the edge, we must create it and update
    all the adjacent lists */

      else
      {
        if ( v1i == 0)
          _new = SurfPushBdryEdge( nd_list, v1, node_indx, edge->use, 1,f_surf) ;
        else
          _new = SurfPushBdryEdge( nd_list, node_indx, v1, edge->use, 1,f_surf) ;

        cur = Msh2DAdjEdgeAlloc() ;
        cur->edge = _new ;

        cur->next = nd_list[v1].edges ;
        nd_list[v1].edges = cur ;

        cur = Msh2DAdjEdgeAlloc() ;
        cur->edge = _new ;

        cur->next = nd_list[node_indx].edges ;
        nd_list[node_indx].edges = cur ;
      }
    }

    /* update all the adjacent edge lists so they no longer point
       to the base edge */

    for ( i=0 ; i<2 ; i++ )
    {
      for ( temp = &(nd_list[edge->verts[i]].edges) ; *temp ;
            temp = &((*temp)->next) )
      {
        if ( (*temp)->edge == edge )
        {
          save = *temp ;
          *temp = (*temp)->next ;
          Msh2DAdjFree( save ) ;
          break ;
        }
      }
    }

    /* update the active flags.  Any node on the new element that
       no longer adjacent to at least on edge on the boundary
       becomes inactive */

    for ( i=0 ; i<2 ; i++ )
    {
      if ( !nd_list[edge->verts[i]].edges )
            nd_list[edge->verts[i]].active_flag = 0 ;
    }
    if ( !nd_list[node_indx].edges )
     nd_list[node_indx].active_flag = 0 ;
    else
     nd_list[node_indx].active_flag = 1 ;

    /* give up the memory associated with this edge */

#if 0
if(flag_debug)
{
G3dForeground(G3D_RED);
G3dLine(node_list[edge->verts[0]].coord[0], node_list[edge->verts[0]].coord[1],
0.0,node_list[edge->verts[1]].coord[0], node_list[edge->verts[1]].coord[1],0.0);
}
#endif


    Msh2DBdryFree( edge ) ;
}


/* -------------------------------------------------------------------
** Msh2DHeap - these routines manage a priorty queue using a heap
**             data structure.
*/

typedef struct _Msh2DHeapEntry
{
    double    value ;
    int       indx ;
} Msh2DHeapEntry, *Msh2DHeap ;

static Msh2DHeap the_heap = 0 ;
static int       alloced_size ;
static int       heap_last ;

static void Msh2DHeapInit( int size )
{
    /* make sure we have enough room for the heap */

    if ( !the_heap ) {
 the_heap = (Msh2DHeap)Msh2DMalloc( size * sizeof(Msh2DHeapEntry) ) ;
 alloced_size = size ;
    }
    else if ( alloced_size < size ) {
 Msh2DFree( the_heap ) ;
 the_heap = (Msh2DHeap)Msh2DMalloc( size * sizeof(Msh2DHeapEntry) ) ;
 alloced_size = size ;
    }
    heap_last = -1 ;
}

static void Msh2DHeapDelete(void)
{
    Msh2DFree( the_heap ) ;
    the_heap = 0 ;
}

static void Msh2DHeapInsert
(
double    value,
int      indx
)
{
    int  cur, parent ;

    /* put the new values at the end of the heap */

    heap_last++ ;
    the_heap[heap_last].value = value ;
    the_heap[heap_last].indx = indx ;

    /* perform a shift up operation to restore the heap property */

    cur = heap_last ;
    while ( cur != 0 )
    {
      parent = cur >> 1 ;    /* fast integer division by 2 */
      if ( the_heap[parent].value <= the_heap[cur].value )
        break ;
      value = the_heap[parent].value ;  /* else swap */
      indx = the_heap[parent].indx ;
      the_heap[parent].value = the_heap[cur].value ;
      the_heap[parent].indx = the_heap[cur].indx ;
      the_heap[cur].value = value ;
      the_heap[cur].indx = indx ;
      cur = parent ;
    }
}

static int Msh2DHeapExtract( double  *value )
{
    int   rtn_indx, cur, child, indx ;
    double tmp_value ;

    if ( heap_last < 0 ) return( -1 ) ;

    /* extract the element at the top of the heap */

    rtn_indx = the_heap[0].indx ;
    *value = the_heap[0].value ;

    /* replace this element with the one at the bottom of the heap */

    the_heap[0].indx = the_heap[heap_last].indx ;
    the_heap[0].value = the_heap[heap_last].value ;
    heap_last-- ;

    /* shift down to restore the heap order */

    cur = 0 ;
    child = 2 * cur ;
    while ( child <= heap_last )
    {
      if ( child+1 <= heap_last )
      {
        if ( the_heap[child+1].value < the_heap[child].value )
          child = child + 1 ;
      }
      if ( the_heap[cur].value <= the_heap[child].value ) break ;
      tmp_value = the_heap[child].value ;  /* else swap */
      indx = the_heap[child].indx ;
      the_heap[child].value = the_heap[cur].value ;
      the_heap[child].indx = the_heap[cur].indx ;
      the_heap[cur].value = tmp_value ;
      the_heap[cur].indx = indx ;
      cur = child ;
      child = 2 * cur ;
    }

    return( rtn_indx ) ;
}


/* -------------------------------------------------------------------
** Msh2DEdgeStack - these routines manage a stack that is used to store
**                  the active boundary edge data.
*/

static Msh2DBdryEdge *Msh2DEdgeAlloc( )
{
    Msh2DBdryEdge  *new_block, *alloced ;
    int            i ;

    /* if the free pointer is null we need to allocate a new block
       of boundary nodes */

    if ( !bdry_free )
    {
     new_block = (Msh2DBdryEdge *)Msh2DMalloc(
      MSH2D_BDRY_FACE_QUANTUM * sizeof(Msh2DBdryEdgeRec) ) ;
     new_block[0].next = bdry_block_ptr ;
     bdry_block_ptr = new_block ;
     for ( i=1 ; i<(MSH2D_BDRY_FACE_QUANTUM-1) ; i++ )
     {
      new_block[i].next = &(new_block[i+1]) ;
     }
     new_block[MSH2D_BDRY_FACE_QUANTUM-1].next = 0 ;
     bdry_free = &(new_block[1]) ;
    }

    /* return the top thing on the free list */

    alloced = bdry_free ;
    bdry_free = bdry_free->next ;

    return( alloced ) ;
}

static void Msh2DBdryFree( Msh2DBdryEdge  *edge )
{
    /* put this edge back on the free list */

    edge->next = bdry_free ;
    bdry_free = edge ;
}

static void Msh2DEdgeFreeAll()
{
    Msh2DBdryEdge  *cur, *next ;

    /* free all blocks allocated to store edge information */

    if ( bdry_block_ptr ) cur = bdry_block_ptr ;
    else return ;

    while ( cur->next ) {
    next = cur->next ;
    Msh2DFree( cur ) ;
    cur = next ;
    }
    Msh2DFree( cur ) ;

    bdry_stack = 0 ;
    bdry_tail  = 0 ;
    bdry_free  = 0 ;
    bdry_try   = 0 ;
    bdry_cursor    = 0 ;
    bdry_block_ptr = 0 ;
}

static void Msh2DBdryPush( Msh2DBdryEdge  *edge )
{
    /* push this edge on the end of the stack that is implemented
       as a doubly linked list */
    edge->prev = bdry_tail ;
    edge->next = 0 ;
    if ( bdry_tail ) bdry_tail->next = edge ;
    bdry_tail = edge ;
    if ( !bdry_stack ) bdry_stack = edge ;
}

static void Msh2DBdryPushCorre( Msh2DBdryEdge  *edge )
{
    Msh2DBdryEdge  *bdry_cur = 0 ;

    /* find the bdry_try, that is, the first edge that went one time
       to the end of the stack */

    bdry_try = 0 ;
    bdry_cur = bdry_stack;
    do {
     if (bdry_cur->use > edge->use)
     {
      bdry_try = bdry_cur ;
      break ;
     }
     bdry_cur = bdry_cur->next ;
    } while ( bdry_cur != bdry_tail ) ;

    /* modify the list to push the edge in its position considering the
       list only from bdry_stack to bdry_try, because edges from bdry_try
       until bdry_tail are edges already tested and pushed to the end of
       the list */

    if (bdry_try == 0) {
     edge->prev = bdry_tail ;
     edge->next = 0 ;
     if ( bdry_tail ) bdry_tail->next = edge ;
     bdry_tail = edge ;
     if ( !bdry_stack ) bdry_stack = edge ;
    }
    else {
     edge->prev = bdry_try->prev ;
     edge->next = bdry_try ;
     if ( bdry_try->prev ) bdry_try->prev->next = edge ;
     bdry_try->prev = edge ;
     if ( bdry_try == bdry_stack ) bdry_stack = edge ;
    }

}

static void Msh2DBdryPushSmall( Msh2DBdryEdge  *edge )
{
    Msh2DBdryEdge  *bdry_cur = 0 ;
    Msh2DBdryEdge  *bdry_real_tail = 0 ;
    int            found_pos = 0 ;

    /* find the key value for the edge */
    edge->key = edge->real_length ;

    /* find the bdry_try, that is, the first edge that went one time
       to the end of the stack */

    bdry_try = 0 ;
    bdry_cur = bdry_stack;
    do {
     if( bdry_cur ) {
      if (bdry_cur->use > 0)
      {
       bdry_try = bdry_cur ;
       break ;
      }
      bdry_cur = bdry_cur->next ;
      if( !bdry_cur ) break ;
     }
    } while ( bdry_cur != bdry_tail ) ;

    /* modify the list to push the edge in its position considering the
       list only from bdry_stack to bdry_try, because edges from bdry_try
       until bdry_tail are edges already tested and pushed to the end of
       the list */

    if ( (bdry_try != 0) && (bdry_try->prev != 0) ) {
     if (bdry_real_tail == 0) {
      bdry_real_tail = bdry_tail ;
      bdry_tail = bdry_try->prev ;
      bdry_tail->next = 0 ;
     }
    }

    /* push this edge in its right position in the list acoording
       wiht its length. Here is only considered the list from the
       bdry_stack until bdry_try */

    if ( !found_pos ) {
     if ( (bdry_stack == 0) || (edge->key <= bdry_stack->key) ) {
      edge->prev = 0 ;
      edge->next = bdry_stack ;
      if ( bdry_stack ) bdry_stack->prev = edge ;
      bdry_stack = edge ;
      if ( !bdry_tail ) bdry_tail = edge ;
      found_pos = 1 ;
     }
    }

    if ( !found_pos ) {
     if ( (bdry_tail != 0) && (edge->key >= bdry_tail->key) ) {
      edge->prev = bdry_tail ;
      edge->next = 0 ;
      if ( bdry_tail ) bdry_tail->next = edge ;
      bdry_tail = edge ;
      found_pos = 1 ;
     }
    }

    if ( !found_pos ) {
     bdry_cur = bdry_stack;
     do {
      if ( edge->key == bdry_cur->key ) {
       edge->prev = bdry_cur->prev ;
       edge->next = bdry_cur ;
       if ( bdry_cur->prev ) bdry_cur->prev->next = edge ;
       bdry_cur->prev = edge ;
       found_pos = 1 ;
       break ;
      }
      else if ( (edge->key > bdry_cur->key) &&
  (edge->key < bdry_cur->next->key) ) {
       edge->prev = bdry_cur ;
       edge->next = bdry_cur->next ;
       if ( bdry_cur->next ) bdry_cur->next->prev = edge ;
       bdry_cur->next = edge ;
       found_pos = 1 ;
       break ;
      }
      bdry_cur = bdry_cur->next ;
     } while ( bdry_cur != bdry_tail ) ;
    }

    /* correct the list to push back edges from bdry_try until bdry_tail */

    if ( (bdry_try != 0) && (bdry_try->prev != 0) ) {
     if (bdry_real_tail != 0) {
      bdry_tail->next = bdry_try ;
      bdry_try->prev = bdry_tail ;
      bdry_tail = bdry_real_tail ;
     }
    }
}

static Msh2DBdryEdge *Msh2DBdryPop()
{
    Msh2DBdryEdge  *popped ;

    /* pop a edge from the front of the stack that is implemented
       as a doubly linked list */

    if (!bdry_stack) return(0) ;
    if ( bdry_stack == bdry_tail ) bdry_tail = 0 ;
    popped = bdry_stack ;
    bdry_stack = bdry_stack->next ;
    if ( bdry_stack ) bdry_stack->prev = 0 ;
    popped->next = popped->prev = 0 ;

    return( popped ) ;
}

static void Msh2DBdryDelete( Msh2DBdryEdge  *edge )
{
    /* delete this edge from the middle of the doubly linked list */

#if 0
if(flag_debug)
{
G3dForeground(G3D_RED);
G3dLine(node_list[edge->verts[0]].coord[0], node_list[edge->verts[0]].coord[1],
0.0,node_list[edge->verts[1]].coord[0], node_list[edge->verts[1]].coord[1],0.0);
}
#endif

    if ( bdry_stack == edge ) bdry_stack = edge->next ;
    if ( bdry_tail  == edge ) bdry_tail  = edge->prev ;
    if ( edge->next ) edge->next->prev = edge->prev ;
    if ( edge->prev ) edge->prev->next = edge->next ;
    Msh2DBdryFree( edge ) ;
}

static void Msh2DBdryReset(void)
{
    /* reset the cursor for scanning through the list */

    bdry_cursor = bdry_stack ;
}

static Msh2DBdryEdge *Msh2DBdryNext(void)
{
    Msh2DBdryEdge *current ;

    /* return the edge the cursor is pointing to, and increment the
       cursor */

    current = bdry_cursor ;
    if ( bdry_cursor ) bdry_cursor = bdry_cursor->next ;
    return( current ) ;
}

/* -------------------------------------------------------------------
** Msh2DAdjElemAlloc - these routines manage the allocation and freeing
**                     of adjacent elem list entries
*/

static Msh2DAdjElem *Msh2DAdjElemAlloc( )
{
    Msh2DAdjElem  *new_block, *alloced ;
    int           i ;

    /* if the free pointer is null we need to allocate a new block
       of elements */

    if ( !adj_elem_free )
    {
     new_block = (Msh2DAdjElem *)Msh2DMalloc(
                  MSH2D_ADJ_ELEM_QUANTUM * sizeof(Msh2DAdjElemRec) ) ;
     new_block[0].next = adj_elem_block_ptr ;
     adj_elem_block_ptr = new_block ;
     for ( i=1 ; i<(MSH2D_ADJ_ELEM_QUANTUM-1) ; i++ )
     {
       new_block[i].next = &(new_block[i+1]) ;
     }
      new_block[MSH2D_ADJ_ELEM_QUANTUM-1].next = 0 ;
      adj_elem_free = &(new_block[1]) ;
    }

    /* return the top thing on the free list */

    alloced = adj_elem_free ;
    adj_elem_free = adj_elem_free->next ;

    return( alloced ) ;
}


static void Msh2DAdjElemFreeAll()
{
    Msh2DAdjElem  *cur, *next ;

    /* free all blocks allocated to store elem information */

    if ( adj_elem_block_ptr ) cur = adj_elem_block_ptr ;
    else return ;

    while ( cur->next )
    {
     next = cur->next ;
     Msh2DFree( cur ) ;
     cur = next ;
    }
    Msh2DFree( cur ) ;

    adj_elem_free = 0 ;
    adj_elem_block_ptr = 0 ;
}

/* -------------------------------------------------------------------
** Msh2DAdjNodeAlloc - these routines manage the allocation and freeing
**                     of node list entries
*/

static void Msh2DAdjNodeFreeAll()
{
    Msh2DBdryNodeList *cur ;

    /* free all blocks allocated to store node information */

    cur = node_list ;

    Msh2DFree( cur ) ;

    node_list = NULL ;
}

/* -------------------------------------------------------------------
** Msh2DAdjEdgeAlloc - these routines manage the allocation and freeing
**                     of adjacent edge list entries
*/

static Msh2DAdjEdge *Msh2DAdjEdgeAlloc( )
{
    Msh2DAdjEdge  *new_block, *alloced ;
    int           i ;

    /* if the free pointer is null we need to allocate a new block
       of boundary nodes */

    if ( !adj_free )
    {
     new_block = (Msh2DAdjEdge *)Msh2DMalloc(
                 MSH2D_ADJ_FACE_QUANTUM * sizeof(Msh2DAdjEdgeRec) ) ;
     new_block[0].next = adj_block_ptr ;
     adj_block_ptr = new_block ;
     for ( i=1 ; i<(MSH2D_ADJ_FACE_QUANTUM-1) ; i++ )
     {
      new_block[i].next = &(new_block[i+1]) ;
     }
     new_block[MSH2D_ADJ_FACE_QUANTUM-1].next = 0 ;
     adj_free = &(new_block[1]) ;
    }

    /* return the top thing on the free list */

    alloced = adj_free ;
    adj_free = adj_free->next ;

    return( alloced ) ;
}

static void Msh2DAdjFree( Msh2DAdjEdge  *edge )
{
    /* put this edge back on the free list */

    edge->next = adj_free ;
    adj_free = edge ;
}

static void Msh2DAdjFreeAll()
{
    Msh2DAdjEdge  *cur, *next ;

    /* free all blocks allocated to store edge information */

    if ( adj_block_ptr ) cur = adj_block_ptr ;
    else return ;

    while ( cur->next )
    {
     next = cur->next ;
     Msh2DFree( cur ) ;
     cur = next ;
    }
    Msh2DFree( cur ) ;

    adj_free = 0 ;
    adj_block_ptr = 0 ;
}

/* -------------------------------------------------------------------
** Msh2DAdjIniEdgeAlloc - these routines manage the allocation and
**                        freeing of initial adjacent edge list entries
*/

static Msh2DAdjIniEdge *Msh2DAdjIniEdgeAlloc( )
{
    Msh2DAdjIniEdge  *new_block, *alloced ;
    int              i ;

    /* if the free pointer is null we need to allocate a new block
       of boundary nodes */

    if ( !adj_ini_free )
    {
     new_block = (Msh2DAdjIniEdge *)Msh2DMalloc(
       MSH2D_ADJ_INI_FACE_QUANTUM * sizeof(Msh2DAdjIniEdgeRec) ) ;
     new_block[0].next = adj_ini_block_ptr ;
     adj_ini_block_ptr = new_block ;
     for ( i=1 ; i<(MSH2D_ADJ_INI_FACE_QUANTUM-1) ; i++ )
     {
       new_block[i].next = &(new_block[i+1]) ;
     }
     new_block[MSH2D_ADJ_INI_FACE_QUANTUM-1].next = 0 ;
     adj_ini_free = &(new_block[1]) ;
    }

    /* return the top thing on the free list */

    alloced = adj_ini_free ;
    adj_ini_free = adj_ini_free->next ;

    return( alloced ) ;
}

static void Msh2DAdjIniFreeAll()
{
    Msh2DAdjIniEdge  *cur, *next ;

    /* free all blocks allocated to store edge information */

    if ( adj_ini_block_ptr ) cur = adj_ini_block_ptr ;
    else return ;

    while ( cur->next )
    {
     next = cur->next ;
     Msh2DFree( cur ) ;
     cur = next ;
    }
    Msh2DFree( cur ) ;

    adj_ini_free = 0 ;
    adj_ini_block_ptr = 0 ;
}

/* -------------------------------------------------------------------
** Msh2DTestStack - these routines manage a stack that is used to store
**                  the boundary edge data already done.
*/

static void Msh2DTestFreeAll()
{
    Msh2DBdryEdge  *cur, *next ;

    /* free all blocks allocated to store edge information */

    if ( test_block_ptr ) cur = test_block_ptr ;
    else return ;

    while ( cur->next )
    {
     next = cur->next ;
     Msh2DFree( cur ) ;
     cur = next ;
    }
    Msh2DFree( cur ) ;

    test_tail  = 0 ;
    test_block_ptr = 0 ;
}

/* -------------------------------------------------------------------
** Msh2DAddElem - add a new element to the element list.
*/

static void Msh2DAddElem
(
int           *num_elems,
int           **elements,
Msh2DBdryEdge *edge,
int           node_indx
)
{
    int    *slot ;

    /* allocate new space if we need it */

 if ( *num_elems >= num_elem_alloced )
 {
   if ( !num_elem_alloced )
   {
     num_elem_alloced = MSH2D_ELEM_QUANTUM ;
     *elements = (int *)Msh2DMalloc(
                 num_elem_alloced * 3 * sizeof(int)) ;
   }
   else
   {
     num_elem_alloced += MSH2D_ELEM_QUANTUM ;
     *elements = (int *)Msh2DRealloc( *elements,
                  num_elem_alloced * 3 * sizeof(int)) ;
   }
 }

    /* store away the vertex numbers */

    slot = &(*elements)[3*(*num_elems)] ;
    slot[0] = edge->verts[0] ;
    slot[1] = edge->verts[1] ;
    slot[2] = node_indx ;
    (*num_elems)++ ;

#if BDR_DRAW
/************* desenha o elemento  */
if(1)
{
 glColor3d (0.0, 0.0, 1.0);
 glBegin( GL_LINES );
   glVertex3d (node_list[edge->verts[0]].coord[0]*fs, node_list[edge->verts[0]].coord[1]*fs, 0.0);
   glVertex3d (node_list[edge->verts[1]].coord[0]*fs, node_list[edge->verts[1]].coord[1]*fs, 0.0);
 glEnd( );

 glColor3d (0.0, 0.0, 1.0);
 glBegin( GL_LINES );
   glVertex3d (node_list[edge->verts[1]].coord[0]*fs, node_list[edge->verts[1]].coord[1]*fs, 0.0);
   glVertex3d (node_list[node_indx].coord[0]*fs, node_list[node_indx].coord[1]*fs, 0.0);
 glEnd( );

 glColor3d (0.0, 0.0, 1.0);
 glBegin( GL_LINES );
   glVertex3d (node_list[edge->verts[0]].coord[0]*fs, node_list[edge->verts[0]].coord[1]*fs, 0.0);
   glVertex3d (node_list[node_indx].coord[0]*fs, node_list[node_indx].coord[1]*fs, 0.0);
 glEnd( );
}
#endif

}



/* -------------------------------------------------------------------
** Msh2DCrossProd - Compute the Cross Prod.
*/
static double Msh2DCrossProd
(
Msh2DBdryNodeList   *nd_list,
int                 I,
int                 J,
int                 K
)
{
   double    kx, ky, jx, jy;

   /* Orientation:

   | J      / K
   |       /
   |      /
   |     /
   |    /
   |   /
   |  /
   | /
   |/ I
   *
 */

   kx = nd_list[K].coord[0] - nd_list[I].coord[0];
   ky = nd_list[K].coord[1] - nd_list[I].coord[1];
   jx = nd_list[J].coord[0] - nd_list[I].coord[0];
   jy = nd_list[J].coord[1] - nd_list[I].coord[1];

   return (kx*jy - ky*jx);
}


/* -------------------------------------------------------------------
** Msh2DCrossEdge
*/
static int Msh2DCrossEdge
(
Msh2DBdryNodeList    *nd_list,
int                  I1,           /* Node pointer   */
int                  I2,           /* Node pointer   */
int                  J1,           /* Node pointer   */
int                  J2            /* Node pointer   */
)
{
  /* locals variables */
  int cross;
  /*  Line orientation:

\I1   /J2
 \   /
  \ /
   \
  / \
 /   \
/J1   \I2
  */
  cross = TRUE;

/*  First do some simple box checks to eliminate lines that are not
    close together */

  if ((nd_list[I1].coord[0] > nd_list[J1].coord[0]) &&
      (nd_list[I1].coord[0] > nd_list[J2].coord[0]) &&
      (nd_list[I2].coord[0] > nd_list[J1].coord[0]) &&
      (nd_list[I2].coord[0] > nd_list[J2].coord[0]))
  {
    cross = FALSE;
    return(cross);
  }

  if ((nd_list[I1].coord[1] > nd_list[J1].coord[1]) &&
      (nd_list[I1].coord[1] > nd_list[J2].coord[1]) &&
      (nd_list[I2].coord[1] > nd_list[J1].coord[1]) &&
      (nd_list[I2].coord[1] > nd_list[J2].coord[1]))
  {
    cross = FALSE;
    return(cross);
  }

/*  Now cross line I with J1 and line I with J2, if have the same sign they
   cannot cross */

  if (Msh2DCrossProd(nd_list,I1,I2,J1) *
      Msh2DCrossProd(nd_list,I1,I2,J2) >= 0.0)
  {
    cross = FALSE;
    return(cross);
  }

/*  Now cross line J with I1 and line J with I2, if have the same sign they
     cannot cross */

  if (Msh2DCrossProd(nd_list,J1,J2,I1) *
      Msh2DCrossProd(nd_list,J1,J2,I2) >= 0.0)
  {
    cross = FALSE;
    return(cross);
  }

  return(cross);
}



/* -------------------------------------------------------------------
** Msh2DCheckArea - Check the area of added element.
*/

static double Msh2DAreaTrian
(
Msh2DBdryEdge      *edge,
Msh2DBdryNodeList  *nd_list,
int                node_indx
)
{
 double area = 0.0;

 area = Msh2DCrossProd( nd_list, edge->verts[0], edge->verts[1], node_indx) / 2.0 ;

 return (area);
}

static int Msh2DCheckArea
(
Msh2DBdryEdge      *edge,
Msh2DBdryNodeList  *nd_list,
int                node_indx
)
{
    double             node[2], old[2] ;
    double             area = 0.0 ;
    double             min_area = 0.0 ;

    /* return if it's not the phase of ideal elements */

    if (edge->use != 0) return 1 ;

    /* compute the area of polygon */

    area = Msh2DAreaTrian( edge, nd_list, node_indx ) ;

    /* compute the minimal node */

    Msh2DFindNode (edge,10,node) ;

    /* save node_indx choosen */

    old[0] = nd_list[node_indx].coord[0] ;
    old[1] = nd_list[node_indx].coord[1] ;

    /* set new node_indx */

    nd_list[node_indx].coord[0] = node[0] ;
    nd_list[node_indx].coord[1] = node[1] ;

    /* compute the minimal area */

    min_area = Msh2DAreaTrian( edge, nd_list, node_indx ) ;

    /* restore node_indx */

    nd_list[node_indx].coord[0] = old[0] ;
    nd_list[node_indx].coord[1] = old[1] ;

    /* Compare the two areas */

    if( fabs(area) > fabs(min_area) )
     return 1 ;
    else
     return 0 ;

}


static int Msh2DChkScan
(
double              original_nodes[][2],
int                 v0,
int                 v1,
int                 indx
)
{
    double              x, y ;
    double              x0, y0, x1, y1;
    double              xa, ya, xb, yb;
    double              max_y, min_y;

    x0 = original_nodes[v0][0];
    y0 = original_nodes[v0][1];
    x1 = original_nodes[v1][0];
    y1 = original_nodes[v1][1];
    x = node_list[indx].coord[0];
    y = node_list[indx].coord[1];

    /*  find the max and min y coordinates   */

    max_y = MAX( y0, y1);
    min_y = MIN( y0, y1);

/*
    the horizontal line does not intersect the line to the
    left of location pt if:
*/
    /*  (1) if line is horizontal        */

    if( y0 == y1 ) return( 0 );

    /*  (2) if the y coordinate of point is outside the range
     max_y and min_y (but not including min_y)      */

    if( (y < min_y) || (y >= max_y) ) return( 0 );

    /*  (3) if given line is vertical and x coordinate of point is
     smaller than that of line        */

    if( (x0 == x1) && (x < x0) ) return( 0 );

    /*  (4) if inclined line is located to the right of given point */
    /*     (first reduce the problem to a case where yb >ya, always)*/

    if( x0 != x1 )
    {
      if( ((x1 > x0) && (y1 > y0)) || ((x1 < x0) && (y1 > y0)) )
      {
 xa = x0;    ya = y0;
 xb = x1;    yb = y1;
      }
      else
      {
 xa = x1;    ya = y1;
 xb = x0;    yb = y0;
      }

      if( ((y - ya) * (xb - xa)) > ((x - xa) * (yb - ya)) )
      {
 return( 0 );
      }
    }

/*
    if we get here that is because the horizontal line passing
    at the location pt intersects the given line to the left of the pt
*/

    return( 1 );

}



/* -------------------------------------------------------------------
** Msh2DCheckOptim - Check the valid for the optimal point inserted.
*/


int Msh2DCheckOptim
(
Msh2DBdryEdge     *edge,
int               indx,
int               num_org_nodes,
int               num_org_edges,
double            original_nodes[][2],
int               original_edges[][2]
)
{
    int           i, ncross ;

    /* test if the node is outside of the body */

    ncross= 0 ;

 for( i = 0; i < num_org_edges; i++ )
    {
     if ( Msh2DChkScan( original_nodes, original_edges[i][0],
      original_edges[i][1], indx )) ncross++ ;
    }
    if( (ncross % 2) == 0 )
     return 0 ;   /* node outside of the body */


 return 1 ;

}


/**************** SurfMetricShape **********************/
static double SurfMetricShape
(
Msh2DBdryEdge      *edge,
Msh2DBdryNodeList  *nd_list,
int                node_indx,
void (*f_surf)(double u,double v,double *f_dev)
)
{
 double la, lb, lc;
 double p;
 double ind_area, ind_rms, metric_sha;

 lc = edge->real_length;

 la = SurfDist3D(&nd_list[node_indx],&nd_list[edge->verts[0]]);

 lb = SurfDist3D(&nd_list[node_indx],&nd_list[edge->verts[1]]);

 p = (la + lb + lc) * 0.5;

 ind_area = sqrt(p*(p-la)*(p-lb)*(p-lc));

 ind_rms = 1/3.0 * (la*la + lb*lb + lc*lc);

 if(ind_area==0.0)
   metric_sha = 100000;
 else
   metric_sha = ind_rms / ind_area;

 return(metric_sha/SHA_IDEAL);
}


/**************** Dist3D *****************************/
static double SurfDist3D
(
Msh2DBdryNodeList  *A,
Msh2DBdryNodeList  *B
)
{
 double du = (A->coord[0] - B->coord[0]);
 double dv = (A->coord[1] - B->coord[1]);
 double distAB, distBA;

 distAB = sqrt(A->E*du*du + 2*A->F*du*dv + A->G*dv*dv);

 distBA = sqrt(B->E*du*du + 2*B->F*du*dv + B->G*dv*dv);

 return((distAB+distBA)/2.0);
}



