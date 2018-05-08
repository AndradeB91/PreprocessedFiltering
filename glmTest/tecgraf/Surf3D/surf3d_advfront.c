/*
** --------------------------------------------------------------
**
** surf3d_advfront.c - Routines to implement surf3D generation based in
**                     boundary contraction (advancing front) technique.
**
** ---------------------------------------------------------------
**
** Created:       14-Jan-2007      Antonio C.O. Miranda
**
** ---------------------------------------------------------------
**
*/
#define SURF3DDRAW      0
#define PRINT_DEBUG     0
#define PRINT_FEEDBACK  1


/* #include <time.h> */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "surf3d_def.h"
#include "surf3d_geom.h"
#include "amr3bind.h"

#if SURF3DDRAW
#include "surf3d_draw.h"
/* #include "vgl/vglc/vglc.h"  */
#endif

#define NUM_STEP_MESS  100

typedef struct point_3d
{
 double x,y,z;
} point_3d;

#undef vector_3d
#define vector_3d point_3d

/* --------------------------------------------------------------
** Private definitons and data types:
*/

#define Surf3D_INSERT_NO_VALID_NODES     3
#define Surf3D_SMOOTH_STEP               2
#define Surf3D_NODE_QUANTUM              5000
#define Surf3D_POLY_FACE_QUANTUM         1000
#define Surf3D_BDRY_FACE_QUANTUM         1000
#define Surf3D_TEST_FACE_QUANTUM         1000
#define Surf3D_ADJ_ELEM_QUANTUM          1000
#define Surf3D_ADJ_FACE_QUANTUM          1000
#define Surf3D_ADJ_INI_FACE_QUANTUM      1000
#define Surf3D_ELEM_QUANTUM              10000
#define BDY_FACTOR                      0.20
#define SHA_IDEAL                       (4/sqrt(3))
#define MUL_MAX                         1.50
#define REL_FACTOR                      1.0
#define NEW_NODES                       2
#define PI                              3.141592653589793324
#define MAX(a,b)                        (((a)>(b))?(a):(b))
#define MIN(a,b)                        (((a)<(b))?(a):(b))
#define TRUE                            1
#define FALSE                           0
#ifndef ABS
#define ABS(a)  (((a)>=0)?(a):(-(a)))
#endif
/* --------------------------------------------------------------
** Private variables:
*/

static int                num_elem_alloced = 0 ;
static int                num_node_alloced = 0 ;
static Surf3DBdryEdge      *test_block_ptr = 0 ;
static Surf3DBdryEdge      *test_tail = 0 ;
static Surf3DBdryEdge      *bdry_free = 0 ;
static Surf3DBdryEdge      *bdry_block_ptr = 0 ;
static Surf3DBdryEdge      *bdry_tail = 0 ;
static Surf3DBdryEdge      *bdry_cursor = 0 ;
static Surf3DBdryEdge      *bdry_try = 0 ;
static Surf3DAdjElem       *adj_elem_free = 0 ;
static Surf3DAdjElem       *adj_elem_block_ptr = 0 ;
static Surf3DAdjEdge       *adj_free = NULL ;
static Surf3DAdjEdge       *adj_block_ptr = 0 ;
static Surf3DAdjIniEdge    *adj_ini_free = 0 ;
static Surf3DAdjIniEdge    *adj_ini_block_ptr = 0 ;

static Surf3DBdryNodeList  *node_list  = NULL ;
static Surf3DBdryEdge      *bdry_stack = 0 ;
static int                  nintnode  = 0 ;
static int                  nbdrynode = 0 ;


#if SURF3DDRAW
static int             *ptr_nelem;
static int             *ptr_nnodes = NULL;
static Surf3DElemList **ptr_elems;
static Surf3DBdryEdge  *draw_edge;
static double          *opt_node;
static double          try_pts[3];
#endif

static double             Surf3D_toler;
Surf3DSizeElement         *Surf3DSizeFunction = NULL;
Surf3DIdealPoints         *Surf3DIPointsFunction = NULL;
Surf3DSnapPoint           *Surf3DSnapPtFunction = NULL;
void (*Surf3DMessFunction) (char *) = NULL;
static void               *Surf3DElemTree = NULL;
static void               *Surf3DEdgeTree = NULL;
static void               *Surf3DNodeTree = NULL;
static Surf3DElemList     *ptr_ElemList;   /* element list */


/* --------------------------------------------------------------
** Private functions prototypes:
*/

double                 Surf3DBdryInit          (int, int, double [][3], double [][3], int [][2]);
static int             Surf3DBdryGetNode       (Surf3DBdryEdge *, int, double *);
static int             Surf3DBdyGetOptimalNode (Surf3DBdryEdge *, int *, double *, int, int,
                                                double [][3], double [][3], int [][2]);
static int             Surf3DSmooth            (int, int, int, int*, Surf3DElemList *);
static Surf3DBdryEdge *Surf3DPushBdryEdge      (Surf3DBdryNode *, int, int, int, int, int) ; /* markos - added parameter for layer */
static int             Surf3DProxEdge          (int , int , int, double, double base[3][3]);
static int             Surf3DProxEdge3d        (int , int , int, double, double base[3][3]); /* markos - same as Surf3DProxEdge, but tests 3D bbox before it (currently unused) */
static int             Surf3DInsideElem        (double pts[3][3], double pts_idx[3], double base[3][3]);
static int             Surf3DInsideElem3d      (double pts[3][3], double pts_idx[3], double base[3][3]); /* markos - same as Surf3DInsideElem, but tests 3D bbox before it (currently unused) */
static int             Surf3DCheckValid        (Surf3DBdryEdge *, int, int) ;
static int             Surf3DCheckArea         (Surf3DBdryEdge *, int) ;
static int             Surf3DDetNode           (Surf3DBdryEdge *, int, int);
static void            Surf3DAddEdges          (Surf3DBdryEdge *, int) ;
static void            Surf3DAddElem           (int *, int **, Surf3DElemList **, Surf3DBdryEdge *, int) ;
static void            Surf3DUpdateAdjNodeElem (int node_indx, Surf3DElemList *elem_list, int new_elem);
static void            Surf3DFindNode          (Surf3DBdryEdge *, int, double *) ;
static int             Surf3DTreeNode          (Surf3DBdryEdge *, int, double *, double *) ;
static int             Surf3DInsertNewNodesBar (Surf3DBdryEdge *, int *, double *) ;
static int             Surf3DAdjEdgsFormTriangle (Surf3DBdryEdge  *,int );
static int             Surf3DCrossEdge3d       (int, int, int, int, Surf3DBdryNode *, double [3][3]); /* markos - same as Surf3DCrossEdge2d, but tests 3D bbox before it (currently unused) */
static int             Surf3DCrossTri3d        (int, int, int, int, int, int, Surf3DBdryNode *, double *, double *, double[3][3]); /* markos - triangle x triangle intersection tests (based on Msh3D) */


static void              Surf3DHeapInit        (int) ;
static void              Surf3DHeapDelete      (void) ;
static void              Surf3DHeapInsert      (double, int) ;
static int               Surf3DHeapExtract     (double *) ;
static double            Surf3DHeapAngle       (int, int, int) ;
static void              Surf3DTestFreeAll     (void) ;
static Surf3DBdryEdge   *Surf3DEdgeAlloc       (void) ;
static void              Surf3DBdryFree        (Surf3DBdryEdge  *edge) ;
static void              Surf3DEdgeFreeAll     (void) ;
static void              Surf3DBdryPush        (Surf3DBdryEdge *) ;
static void              Surf3DBdryPushCorre   (Surf3DBdryEdge *) ;
static void              Surf3DBdryPushSmall   (Surf3DBdryEdge *) ;
static void              Surf3DBdryPushLayer   (Surf3DBdryEdge *) ; /* markos - adding edge to the front based on use, layer and length, in that order */
static Surf3DBdryEdge *  Surf3DBdryPop         (void) ;
static void              Surf3DBdryDelete      (Surf3DBdryEdge  *edge) ;
static Surf3DAdjEdge  *  Surf3DAdjEdgeAlloc    (void) ;
static void              Surf3DAdjFree         (Surf3DAdjEdge  *edge) ;
static void              Surf3DAdjFreeAll      (void) ;
static Surf3DAdjIniEdge *Surf3DAdjIniEdgeAlloc (void) ;
static void              Surf3DAdjIniFreeAll   (void) ;
#if 0
static Surf3DAdjElem    *Surf3DAdjElemAlloc    (void) ;
#endif
static void              Surf3DAdjElemFreeAll  (void) ;
static void              Surf3DAdjNodeFreeAll  (void) ;
static void              Surf3DBdryReset       (void) ;
static Surf3DBdryEdge   *Surf3DBdryNext        (void) ;

#if 0
static double          Surf3DCrossProd( int, int, int);
static int             Surf3DCrossEdge( int, int, int, int);
#endif

/* markos - functions that test the triangle x triangle intersection, based on Msh3D */
#define ABOUT_ZERO(val,tol)             ((val<tol)&&(val>(-tol)))
static int     GeoTriIntersect(double a[3], double b[3], double c[3], int p[3], double u[3], double v[3], double w[3], int q[3], double tol_inters);
static int     GeoTriSegIntersect(double a[3], double b[3], double c[3], double u[3], double v[3], double tol_inters);
static int     GeoTriSegIntersectCheck(double p[3], double v[3], double tol_inters);
static int     GeoSegSegIntersect(double p1[3], double q1[3], double p2[3], double q2[3], double p[3], double tol_inters);
static double  GeoNorm(double v[3]);
static double *GeoCrossProd(double a[3], double b[3], double n[3]);
static double *GeoCross(double a[3], double b[3], double c[3], double n[3]);
static double *GeoCrossNorm(double a[3], double b[3], double c[3], double n[3]);
static double  GeoDot(double u[3], double v[3]);

#if SURF3DDRAW
/*
static void Delay( void )
{
  long double init_time = clock( );
  long double total_time = (1.0);
  // return;
  do { } while ((clock ()-init_time) / CLOCKS_PER_SEC  < total_time);
}*/
#endif


#if 1
void cross(vector_3d *u,vector_3d *v,vector_3d *w)
{
 w->x=u->y*v->z - u->z*v->y;
 w->y=u->z*v->x - u->x*v->z;
 w->z=u->x*v->y - u->y*v->x;
}

static int colinear(double aX, double aY, double aZ,
                    double bX, double bY, double bZ,
                    double cX, double cY, double cZ)
{
 vector_3d u,v,w;
 u.x = bX-aX; u.y=bY-aY; u.z=bZ-aZ;
 v.x = cX-aX; v.y=cY-aY; v.z=cZ-aZ;
 cross(&u,&v,&w);

 return ((fabs(w.x)+fabs(w.y)+fabs(w.z))< Surf3D_toler*0.01);
}
#endif

/* -------------------------------------------------------------
** Surf3DBdryContraction - main driver for the boundary contraction
**                        algorithm.
*/
int Surf3DBdryContraction
(
int                  num_org_nodes,
int                  num_org_edges,
double               original_nodes[][3],
double               original_normal[][3],
int                  original_edges[][2],
int                  *num_int_nodes,
double               **internal_nodes,
int                  *num_gen_elements,
int                  **generated_elements
)
{
  Surf3DBdryEdge    *edge = NULL;
  int               num_nodes = 0, node_indx, i, j;
  /*double            smallest_edge = 0.0;*/
  Surf3DElemList    *ElemList = NULL;   /* element list */
  char              message[256];
  int               curr_part_elem = 1;
/*
  long double       cpu_time;
*/
  int pos;
  double area;


  /* Init structures */
  /* 2.1 build the node list, then push all the boundary edge
          descriptions onto a stack  */
  *num_gen_elements = 0 ;
  nintnode  = (*num_int_nodes) ;
  nbdrynode = num_org_nodes ;
  num_nodes = num_org_nodes; /*  + (*num_int_nodes) ; */
  /*smallest_edge  = */Surf3DBdryInit (num_org_nodes, num_org_edges, original_nodes,
                                   original_normal, original_edges);

#if 0
  printf ("advance front\n");
  for (i = 0; i < num_org_nodes; ++i)
  {
   printf ("%f %f %f\n", original_normal[i][0], original_normal[i][1], original_normal[i][2]);
  }
  printf ("\n");

  printf ("node_list\n");
  for (i = 0 ; i < num_org_nodes ; i++)
  {
   printf ("%f %f %f\n", node_list[i].normal[0], node_list[i].normal[1], node_list[i].normal[2]);
  }
  printf ("\n");
#endif


#if SURF3DDRAW
  ptr_nelem  = num_gen_elements;
  ptr_elems  = &ElemList;
  ptr_nnodes = &num_nodes;
#endif


  /* 2.2 use the edge on the top of the stack.  If the stack is
        empty then we are done */
/*
  cpu_time = clock( );
*/
  pos=0;
  while ( (edge = Surf3DBdryPop()) )
  {
    /*printf("markos - testing edge %d %d, use %d\n", edge->verts[0]+1, edge->verts[1]+1, edge->use);*/

    double _new[3]; /* optimal node coord. */
    pos++;
#if PRINT_DEBUG
    fprintf (stdout, "pos %d\n",pos);
    fprintf (stdout, "\n\nTentando aresta %d-%d, use %d\n", edge->verts[0], edge->verts[1], edge->use);
#endif
    ptr_ElemList = ElemList;

    /* Try to get a node from the current boundary */
    node_indx = Surf3DBdryGetNode (edge, num_nodes, _new);

#if SURF3DDRAW
    if ((*num_gen_elements) > curr_part_elem*NUM_STEP_MESS)
    {
      draw_edge = edge;
      opt_node  = _new;
      Surf3DDraw (NULL);
    }
#endif

#if PRINT_DEBUG
    if (node_indx != -1 && node_indx != -2)
      fprintf (stdout, "pegou um ponto \n");
    else
      fprintf (stdout, "recusou pontos \n");
#endif

    /* Try to insert an optimal node */
    if (node_indx == -1 && edge->use == 0)
    {
      node_indx = Surf3DBdyGetOptimalNode (edge, &num_nodes, _new, num_org_nodes, num_org_edges,
                                           original_nodes, original_normal, original_edges);
#if PRINT_DEBUG
      if (node_indx != -1)
        fprintf (stdout, "gerou um ponto \n");
      else
        fprintf (stdout, "nao gerou ponto \n");
#endif

      if (node_indx == -1)
      {
        num_nodes--;
        nintnode--;
      }
    }

    if (node_indx != -1 && node_indx != -2) /* It found a valid node to create a new element */
    {
#if PRINT_FEEDBACK
      if (Surf3DMessFunction != NULL && (*num_gen_elements) > curr_part_elem*NUM_STEP_MESS)
      {
        sprintf (message, "%d elements generated...", *num_gen_elements);
        (*Surf3DMessFunction) (message);
        curr_part_elem++;
#if PRINT_DEBUG
   fflush (stdout);
#endif
      }
#endif
#if 0 /* markos - changed to zero, because collinearity test is now performed in Surf3DCheckValid() */
      area=colinear(node_list[edge->verts[0]].coord[0],node_list[edge->verts[0]].coord[1],node_list[edge->verts[0]].coord[2],
                    node_list[node_indx].coord[0],node_list[node_indx].coord[1],node_list[node_indx].coord[2],
                    node_list[edge->verts[1]].coord[0],node_list[edge->verts[1]].coord[1],node_list[edge->verts[1]].coord[2]);

      if (!area)
      {
        /*printf("markos - creating element %d: %d %d %d\n", (*num_gen_elements)+1, edge->verts[0]+1, edge->verts[1]+1, node_indx+1);*/
        Surf3DAddElem  (num_gen_elements, generated_elements, &ElemList, edge, node_indx);
        Surf3DAddEdges (edge, node_indx);
      }
      else
        continue;
#else
      /*printf("markos - creating element %d: %d %d %d\n", (*num_gen_elements)+1, edge->verts[0]+1, edge->verts[1]+1, node_indx+1);*/
      Surf3DAddElem  (num_gen_elements, generated_elements, &ElemList, edge, node_indx);
      Surf3DAddEdges (edge, node_indx);
#endif
    }
    else  /* It did not find a valid node => element is not created */
    {
      if(edge->use != 3) /* markos - added one more use for edge */
      {
        edge->use++;
        Surf3DBdryPushLayer (edge) ; /* markos - pushing based on layer */
        continue ;
      }
    }

  }

/*
  cpu_time = (clock( ) - cpu_time)/CLOCKS_PER_SEC;
  printf("\n");
  printf("\t\tCPU time boundary contraction............. %0.3f (s)\n", (double)cpu_time);
*/


#if SURF3DDRAW
  draw_edge = NULL;
  Surf3DDraw (NULL);
#endif

   /* 2.6 Smooth the internal nodes, because the boundary nodes can't be changed */
   if (Surf3DMessFunction != NULL)
  (*Surf3DMessFunction) ("Smoothing...");
/*
  cpu_time = clock( );
*/

  Surf3DSmooth (num_nodes, nbdrynode, *num_gen_elements, *generated_elements, ElemList);

/*
  cpu_time = (clock( ) - cpu_time)/CLOCKS_PER_SEC;
  printf("\n");
  printf("\t\tCPU time smoothing............. %0.3f (s)\n", (double)cpu_time);
*/

   /* 2.7 Get the generated nodes */
  *num_int_nodes    = (nbdrynode+nintnode) ;
  (*internal_nodes) = (double *) calloc (3*(*num_int_nodes), sizeof(double));
  for (i = 0; i < (*num_int_nodes); i++)
  {
    for (j = 0 ; j < 3 ; j++)
      (*internal_nodes)[i*3+j] = node_list[i].coord[j];
  }


  /* 2.8 release memory and return status indicating that a mesh was generated */
  Surf3DHeapDelete() ;
  Surf3DEdgeFreeAll() ;
  Surf3DTestFreeAll() ;
  Surf3DAdjFreeAll() ;
  Surf3DAdjIniFreeAll() ;
  Surf3DAdjElemFreeAll() ;
  Surf3DAdjNodeFreeAll() ;
  free (ElemList);
  RtreeDestroy (Surf3DElemTree);
  RtreeDestroy (Surf3DEdgeTree);
  RtreeDestroy (Surf3DNodeTree);
  Surf3DSizeFunction = NULL;
  Surf3DIPointsFunction = NULL;
  Surf3DSnapPtFunction  = NULL;

  return (1) ;
}



/********************** Surf3DAddNodeToRTree ******************************/
static void	Surf3DAddNodeToRTree (int i, double x, double y, double z, double size)
{
  double xmin, ymin, zmin, xmax, ymax, zmax;
  int    *index = (int *) malloc (sizeof (int));
  xmin = x - size*0.5;
  ymin = y - size*0.5;
  zmin = z - size*0.5;
  xmax = x + size*0.5;
  ymax = y + size*0.5;
  zmax = z + size*0.5;
  index[0] = i;

  RtreeInsert (Surf3DNodeTree, (int *) index, xmin, xmax, ymin, ymax, zmin, zmax);
}


/* --------------------------------------------------------------
** Private functions:
*/

/********************** Surf3DBdryInit ******************************/
double Surf3DBdryInit (int num_nodes, int num_edges, double _nodes[][3],
                       double _normal[][3], int _edges[][2])
{
  Surf3DBdryEdge     *edge = NULL;
  Surf3DBdryNode     *node = NULL;
  Surf3DAdjEdge      *adj_edge = NULL;
  Surf3DAdjIniEdge   *adj_ini_edge = NULL;
  int               i, j, k;
  double            smaller_edge=0;

  num_elem_alloced  = 0 ;
  bdry_stack = NULL;

  node_list = (Surf3DBdryNodeList *) Surf3DMalloc(num_nodes * sizeof (Surf3DBdryNodeRec));
  num_node_alloced = num_nodes;
  for (i = 0 ; i < num_nodes ; i++)
  {
    for ( j=0 ; j<3 ; j++ )
    {
      node_list[i].coord[j]  = _nodes [i][j];
      node_list[i].normal[j] = _normal[i][j];
    }
    node_list[i].active_flag = 1 ;
    node_list[i].edges  = NULL;
    node_list[i].iedges = NULL;
    node_list[i].elems  = NULL;
  }

  /* init edges and edgetree */
  Surf3DEdgeTree = RtreeCreate( );
  for (i = 0; i < num_edges; i++)
  {
    edge = Surf3DPushBdryEdge (node_list, _edges[i][0], _edges[i][1], 0, 0, 0); /* markos - added parameter for layer */
    if ((smaller_edge == 0.0) || (smaller_edge > edge->length))
      smaller_edge = edge->length;
  }

  Surf3D_toler = smaller_edge / 100.0;

  /* Add nodes to RTree */
  Surf3DNodeTree = RtreeCreate( );
  for (i = 0 ; i < num_nodes ; i++)
  {
    Surf3DAddNodeToRTree (i, node_list[i].coord[0],node_list[i].coord[1],
                          node_list[i].coord[2], Surf3D_toler);
  }


  /* build the list of adjacent edges for all the boundary nodes */
  Surf3DBdryReset();
  while ( (edge = Surf3DBdryNext()) )
  {
    for ( j=0; j<2; j++)
    {
      node = &node_list[edge->verts[j]];
      adj_edge = Surf3DAdjEdgeAlloc();
      adj_edge->next = node->edges;
      adj_edge->edge = edge;
      node->edges = adj_edge;
    }
  }

  /* build the list of initial adjacent edges for all the boundary nodes */
  Surf3DBdryReset ( );
  while ( (edge = Surf3DBdryNext ()) )
  {
    for ( j=0 ; j<2 ; j++ )
    {
      node = &node_list[edge->verts[j]];
      adj_ini_edge = Surf3DAdjIniEdgeAlloc ( );
      adj_ini_edge->next = node->iedges;
      for ( k=0; k<2; k++ )
        adj_ini_edge->verts[k] = edge->verts[k] ;
      for ( k=0; k<3; k++ )
        adj_ini_edge->nrm[k] = edge->nrm[k];
      node->iedges = adj_ini_edge ;
    }
  }

  /* init elem tree */
  Surf3DElemTree = RtreeCreate( );

  return smaller_edge;
}

/********************** Surf3DBdryGetOptimalPhase ***************************/
static int Surf3DCheckTriangle (Surf3DBdryEdge *edge)
{
    /* markos - changed the test if the edge has a missing triangle adjacent */

    /* this test checks if there are two front edges adjacent to the base edge
       that share a common node, i.e., if there is a missing triangle just
       waiting to be created.

       the previous test did not do it properly, because it used only the last
       adjacent front edges, instead of all adjacent front edges. */

    Surf3DAdjEdgeList *adj_edge1, *adj_edge2;
    int adj_vert_indx1, adj_vert_indx2;

    for (adj_edge1 = node_list[edge->verts[0]].edges; adj_edge1; adj_edge1 = adj_edge1->next)
    {
        if (adj_edge1 == edge)
        {
            continue;
        }

        adj_vert_indx1 = (adj_edge1->edge->verts[0] == edge->verts[0]) ? adj_edge1->edge->verts[1] : adj_edge1->edge->verts[0];

        for (adj_edge2 = node_list[edge->verts[1]].edges; adj_edge2; adj_edge2 = adj_edge2->next)
        {
            if (adj_edge2 == edge)
            {
                continue;
            }

            adj_vert_indx2 = (adj_edge2->edge->verts[0] == edge->verts[1]) ? adj_edge2->edge->verts[1] : adj_edge2->edge->verts[0];

            if (adj_vert_indx1 == adj_vert_indx2)
            {
                return adj_vert_indx1;
            }
        }
    }

    return -1;


  Surf3DAdjEdgeList *adj_edge;
  int               found_adj[2] = {-1, -1};

  adj_edge = node_list[edge->verts[0]].edges;
  do
  {
    /* try to found edge->verts[0] */
    if (edge != adj_edge->edge)
    {
      if (adj_edge->edge->verts[0] != edge->verts[0])
        found_adj[0] = adj_edge->edge->verts[0];
      else
        found_adj[0] = adj_edge->edge->verts[1];
    }
    adj_edge = adj_edge->next;
  } while (adj_edge != NULL);

  adj_edge = node_list[edge->verts[1]].edges;
  do
  {
    /* try to found edge->verts[0] */
    if (edge != adj_edge->edge)
    {
      if (adj_edge->edge->verts[0] != edge->verts[0])
        found_adj[1] = adj_edge->edge->verts[0];
      else
        found_adj[1] = adj_edge->edge->verts[1];
    }
    adj_edge = adj_edge->next;
  } while (adj_edge != NULL);

  if (found_adj[0] == found_adj[1] && edge->verts[0] != found_adj[0] &&
      edge->verts[1] != found_adj[0])
    return found_adj[0];
  return -1;

}


/********************** Surf3DBdryGetOptimalPhase ***************************/
static void Surf3DGeometryBasedNode (Surf3DBdryEdge *curr_edge, int n_nodes,
                                     double *_new, double dist)
{
  double metric, cand_vec[3], dot;
  int    /*node_indx = -1,*/ j;
  double xmin, ymin, zmin, xmax, ymax, zmax;
  double xmn, ymn, zmn, xmx, ymx, zmx;
  int    *index;

  xmin = _new[0] - dist;
  ymin = _new[1] - dist;
  zmin = _new[2] - dist;
  xmax = _new[0] + dist;
  ymax = _new[1] + dist;
  zmax = _new[2] + dist;

  RtreeInitSearchBox (Surf3DNodeTree, xmin, xmax, ymin, ymax, zmin, zmax);
  while ((index = (int *)
          RtreeSearchBox(Surf3DNodeTree, &xmn, &xmx, &ymn,
                                            &ymx, &zmn, &zmx)) != NULL)
  {
    if (*index == curr_edge->verts[0] || *index == curr_edge->verts[1])
      continue;

    if (node_list[*index].active_flag )
    {
      double d, dist_vec[3];

      /* find the vector from the center to the new point
          and make sure that this cross with the normal is
          positive */
      for (j = 0; j < 3; j++)
        cand_vec[j] = node_list[*index].coord[j] - curr_edge->center[j] ;

      dot = cand_vec[0] * curr_edge->nrm[0] +
            cand_vec[1] * curr_edge->nrm[1] +
            cand_vec[2] * curr_edge->nrm[2] ;

      if ( dot <= 0.0 ) continue ;

      /* verify if node is out of sphere centered in optimal node
          for the edge and ratio equal of largest (smallest)
          distance from this node to any vertex of the edge.
          This should be doen only the fase of ideal elements. */
      for (j = 0; j < 3; j++)
        dist_vec[j] = node_list[*index].coord[j] - _new[j] ;

      d = sqrt (dist_vec[0] * dist_vec[0] +
                dist_vec[1] * dist_vec[1] +
                dist_vec[2] * dist_vec[2]);

      if ( d >= dist ) continue ;

      /* the metric we are currently using is the square of
          the distance from the center of the edge to the
          candidate node */
      metric = Surf3DHeapAngle (curr_edge->verts[0], curr_edge->verts[1], *index);
      if( metric != 0.0 )    /* multiply by -1 because of heap */
        Surf3DHeapInsert( (-1.0*metric), *index) ;

    }
  }

}


/********************** Surf3DTopologyBasedNode ***************************/
static void Surf3DTopologyBasedNode (Surf3DBdryEdge *curr_edge, int n_nodes,
                                     double *_new, double dist)
{
  double metric, cand_vec[3], dot;
  int    i, j;

  for (i = 0; i < n_nodes; i++)
  {
    if (node_list[i].active_flag )
    {

      /* find the vector from the center to the new point
          and make sure that this cross with the normal is
          positive */
      for (j = 0; j < 3; j++)
        cand_vec[j] = node_list[i].coord[j] - curr_edge->center[j] ;

      dot = cand_vec[0] * curr_edge->nrm[0] +
            cand_vec[1] * curr_edge->nrm[1] +
            cand_vec[2] * curr_edge->nrm[2] ;

      if ( dot <= 0.0 ) continue ;

      /* the metric */
      /* markos - added one more use for edge, that also uses angle metric */
      if (curr_edge->use <= 2)
      {
        metric = Surf3DHeapAngle (curr_edge->verts[0], curr_edge->verts[1], i);
        if (metric != 0.0)
          Surf3DHeapInsert( -metric, i ) ;
      }
      else
      {
      /* markos - previous use for edge */
        metric = 0.0;
        for ( j=0 ; j<3 ; j++ )
            metric += cand_vec[j] * cand_vec[j] ;
        metric = metric / dot;
        Surf3DHeapInsert( metric, i ) ;
      }
    }
  }

}
/********************** Surf3DBdryGetNode ******************************/
static int Surf3DBdryGetNode (Surf3DBdryEdge *curr_edge, int n_nodes,
                             double *_new)
{
  double h = 0.0, dist, metric, next_metric;
  int    node_indx = -1, next_indx, check;
  int    topology_on_close_nodes = 0, i = 0; /* markos - new use for edge */

#if 1
  node_indx = Surf3DCheckTriangle (curr_edge);
  if (node_indx != -1)
    return node_indx;
#endif

  /* examine all active nodes. Rank all nodes as to the goodness
     of the triangle they will form with the current edge */
  if (curr_edge->use == 0)
  {
    if (!Surf3DTreeNode (curr_edge, 1, _new, &h))
      return -2;
    if (h <= 0.0)
      return -2;

  }
  dist = h * 0.85;

  /* markos - added one more use for edge
     use = 0 -> same as previous use 0, can add new node, uses rtree to search for near nodes, uses angle for metric
     use = 1 -> new use added, cannot add new node, uses rtree to search for near nodes, uses angle for metric
                (if no distance is found in Surf3DTreeNode(), does the same as use = 3)
     use = 2 -> modified from previous use 1, cannot add new node, does not use rtree to search for nodes, uses angle for metric
                (previous use 1 used dot product for metric)
     use = 3 -> same as previous use 2, cannot add new node, does not use rtree to search for nodes, uses dot product for metric
   */
  if (curr_edge->use == 1)
  {
    for (i = 1; i < 10; i++)
    {
      if (Surf3DTreeNode (curr_edge, i, _new, &h))
      {
        if (h > 0.0)
        {
          topology_on_close_nodes = 1;
          dist = h * i * 1.5;
          break;
        }
      }
    }

    if (topology_on_close_nodes == 0)
    {
      Surf3DIdealPoints *func = Surf3DIPointsFunction;
      Surf3DIPointsFunction = NULL;

      if (Surf3DTreeNode (curr_edge, 1, _new, &h))
      {
        if (h > 0.0)
        {
          topology_on_close_nodes = 1;
          dist = h * 1.5;
        }
      }

      Surf3DIPointsFunction = func;
    }
  }

  Surf3DHeapInit (n_nodes);

  /* get nodes to heap list */
  if (curr_edge->use == 0)
    Surf3DGeometryBasedNode (curr_edge, n_nodes, _new, dist);
  /* markos - new use for edge */
  else if ((curr_edge->use == 1) && (topology_on_close_nodes))
    Surf3DGeometryBasedNode (curr_edge, n_nodes, _new, dist);
  else
    Surf3DTopologyBasedNode (curr_edge, n_nodes, _new, dist);



/* 2.4 start with the node with the best ranking.  Check to make
      sure it will form a valid element (i.e., does not intersect
      the current boundary).  If the element is invalid go on
      to the next.  If the element is valid, then look at the next
      candidate node.  If the two nodes have the same coordinates
      then use topology to determine the one we want. */

  while (1)
  {
    /* extract a node from the heap based in its metric */
    node_indx = Surf3DHeapExtract (&metric) ;

    if ((node_indx == curr_edge->verts[0]) || (node_indx == curr_edge->verts[1]))
      continue ;

    /* here a node was found to make the element (node_indx was
        extracted from the heap). Check its validity */
    if (node_indx != -1)
    {
      check = 1;

      /* check validity for node choosen */
      if (curr_edge->use == 0 && check)
        check = Surf3DCheckArea (curr_edge, node_indx);

      if (check)
          check = Surf3DCheckValid (curr_edge, node_indx, 1);

      /* set validity for node choosen */
      if (check)
      {
        /* check if there is another node with equal metric - crack */
        next_indx = Surf3DHeapExtract (&next_metric);
        /* markos - both nodes must belong to the boundary (initial front) */
        if ((next_indx != -1) && (node_list[node_indx].iedges) && (node_list[next_indx].iedges))
        {
          if (fabs (next_metric - metric) < Surf3D_toler &&
              fabs (node_list[node_indx].coord[0] - node_list[next_indx].coord[0]) < Surf3D_toler &&
              fabs (node_list[node_indx].coord[1] - node_list[next_indx].coord[1]) < Surf3D_toler &&
              fabs (node_list[node_indx].coord[2] - node_list[next_indx].coord[2]) < Surf3D_toler)
          {
            node_indx = Surf3DDetNode (curr_edge, node_indx, next_indx);
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


/********************** Surf3DBdyGetOptimalNode  ******************************/
static int Surf3DBdyGetOptimalNode (Surf3DBdryEdge *curr_edge, int *n_nodes,
                                   double *_new, int num_org_nodes, int num_org_edges,
                                   double original_nodes[][3], double original_normal[][3],
                                   int original_edges[][2])
{
  int node_indx = -1;
  int check;

  /* try insertion of new node */
  check = Surf3DInsertNewNodesBar (curr_edge, n_nodes, _new);

#if PRINT_DEBUG
  fprintf (stdout, "Tentando inserir ponto (%f, %f, %f)\n", _new[0], _new[1], _new[2]);
  if (check)
    fprintf (stdout, "  Ponto inicialmente inserido no NewNodesBar\n");
  else
    fprintf (stdout, "  Ponto nao inserido - Nao passou no NewNodesBar\n");
#endif

  if (check)
  {
    check = Surf3DCheckValid (curr_edge, *n_nodes-1, -1);

#if PRINT_DEBUG
  if (check)
    fprintf (stdout, "  Passou no Surf3DCheckValid\n");
  else
    fprintf (stdout, "  Nao Passou no Surf3DCheckValid\n");
#endif
  }

  /* set validity for node inserted */
  if (check)
  {
    node_indx = *n_nodes - 1 ;
    Surf3DAddNodeToRTree (node_indx, _new[0], _new[1], _new[2], Surf3D_toler);
    #if PRINT_DEBUG
      fprintf (stdout, "Ponto inserido -%d- (%f, %f, %f)\n", node_indx, _new[0], _new[1], _new[2]);
    #endif
    return node_indx;
  }

  return -1;
}


/* -------------------------------------------------------------------
** Surf3DSmooth - driver for nodal smoothing.
*/

static int Surf3DSmooth
(
int            num_nodes,
int            bdr_nodes,
int            num_elems,
int            *elements,
Surf3DElemList *elemlist
)
{
  Surf3DBdryNodeList *smoo_list ;
  Surf3DAdjElemList  *elems ;
  int                i, j, k, n, m, id ;
  double             new_pts[3], xi[3], w, wj, box;
  double             dx, dy, dz, dist;

  /* 3.1 Allocate memory for smooth nodes vector */
  smoo_list = (Surf3DBdryNodeList *) Surf3DMalloc (num_nodes * sizeof (Surf3DBdryNodeRec));

  /* 3.2 Initiate list of all elements adjacent to a given node and
        initiate list of smmoth nodes vector */
  for( i = 0; i < num_nodes; i++ )
  {
    for (m = 0; m < 3; m++)
      smoo_list[i].coord[m] = node_list[i].coord[m] ;
  }


  /* 3.4 Move each internal node to the centroid of all its adjacent
         nodes defined by its adjacent elements. Each internal node
         only will be moved if doesn't affect the consistency of its
         adjacent elements, so consistency tests are necessary after
         each move. */
  wj = 1.0 ;
  for (i = 0; i < Surf3D_SMOOTH_STEP; i++)
  {
    for (j = bdr_nodes; j < num_nodes; j++)  /* for all internal nodes */
    {
      /* init variables */
      n = 0 ;
      xi[0] = xi[1] = xi[2] = w = 0.0 ;
      elems = node_list[j].elems ;

      /* for all adjacent elements */
      while (elems)
      {
        for (k = 0; k < 3; k++)
        {
          id = elemlist[elems->id].conn[k];
          if( id != j )
          {
            for (m = 0; m < 3; m++)
              xi[m] += (wj * (node_list[id].coord[m] - node_list[j].coord[m])) ;
            w += wj ;
            n++ ;
          }
        }
        elems = elems->next ;
      }

      /* compute smoothing coordinate */
      if (node_list[j].elems != NULL)
      {
        if (w > 0.0)
        {
          for (m = 0; m < 3; m++)
            smoo_list[j].coord[m] = REL_FACTOR * (xi[m] / w);
        }
      }
    }

    /* update node coordenates */
    for (j = bdr_nodes; j < num_nodes; j++)
    {
      if (Surf3DSnapPtFunction != NULL)
      {
        box = fabs(smoo_list[j].coord[0] + smoo_list[j].coord[1] +
                   smoo_list[j].coord[2]) / 3.0;
        xi[0] = node_list[j].coord[0] + smoo_list[j].coord[0];
        xi[1] = node_list[j].coord[1] + smoo_list[j].coord[1];
        xi[2] = node_list[j].coord[2] + smoo_list[j].coord[2];
        if (box == 0.0)
          box = 0.001;
        (*Surf3DSnapPtFunction) (xi, node_list[j].normal, box, new_pts);

        // test uvs difference at closest point(u,v coordinates maybe outside the surface)
#if 0
        if ((fabs(node_list[j].coord[0]-new_pts[0]) > 3.0) ||
            (fabs(node_list[j].coord[1]-new_pts[1]) > 3.0) ||
            (fabs(node_list[j].coord[2]-new_pts[2]) > 3.0))
        {
         for (m = 0; m < 3; m++)
           node_list[j].coord[m] = node_list[j].coord[m] + smoo_list[j].coord[m];
        }
        else
        {
        for (m = 0; m < 3; m++)
          node_list[j].coord[m] = new_pts[m];
      }
#else
        // compute distance
        dx = new_pts[0] - node_list[j].coord[0];
        dy = new_pts[1] - node_list[j].coord[1];
        dz = new_pts[2] - node_list[j].coord[2];
        dist = sqrt (dx*dx + dy *dy + dz*dz);

        // try avoiding to go to wrong positions
        if (dist < 2*box)
        {
          for (m = 0; m < 3; m++)
            node_list[j].coord[m] = new_pts[m];
        }
#endif
      }
      else
      {
        for (m = 0; m < 3; m++)
          node_list[j].coord[m] = node_list[j].coord[0] + smoo_list[j].coord[0];
      }
    }

    /* Update nodal normals */

  }

  /* 3.5 Release memory for smooth nodes vector */
  Surf3DFree (smoo_list);

  /* 3.6 Return smooth status */
  return(1) ;
}

/* -------------------------------------------------------------------
** Surf3DPushBdryEdge - this routine pushes boundary edges onto the stack.
*/

static Surf3DBdryEdge *Surf3DPushBdryEdge
(
Surf3DBdryNode   *nodes ,
int             v1 ,
int             v2 ,
int             use,
int             init,
int             layer     /* markos - added parameter for layer */
)
{
    Surf3DBdryEdge  *edge;
    double          tmp_normal[3], tmp_nrm[3], tmp_len;
    register  int   i ;

    edge = Surf3DEdgeAlloc() ;

    edge->verts[0] = v1 ; edge->verts[1] = v2 ;

    edge->use = use ;
    edge->layer = layer; /* markos - added parameter for layer */

    /* the center of the edge is the algebraic mean of the corners */
    for (i = 0 ; i < 3 ; i++)
    {
      edge->center[i] = (nodes[edge->verts[0]].coord[i] +
                         nodes[edge->verts[1]].coord[i]) * 0.5;
    }

    /*  find the max and mins of the coordinates */
    for (i = 0 ; i < 3 ; i++)
    {
      edge->max[i] = nodes[v1].coord[i] ;
      edge->min[i] = nodes[v1].coord[i] ;
      if (edge->max[i] < nodes[v2].coord[i]) edge->max[i] = nodes[v2].coord[i] ;
      if (edge->min[i] > nodes[v2].coord[i]) edge->min[i] = nodes[v2].coord[i] ;
    }

    /* compute the r vector  */
    for (i = 0 ; i < 3 ; i++)
      edge->r[i] = nodes[v2].coord[i] - nodes[v1].coord[i] ;

    /* edge length */
    edge->length = sqrt((edge->r[0]*edge->r[0]) + (edge->r[1]*edge->r[1]) + (edge->r[2]*edge->r[2]));

    for (i = 0 ; i < 3 ; i++)
    {
      tmp_normal[i] = (nodes[edge->verts[0]].normal[i] +
                       nodes[edge->verts[1]].normal[i]) * 0.5;
    }

    tmp_nrm[0] = tmp_normal[1]*edge->r[2] - edge->r[1]*tmp_normal[2];
    tmp_nrm[1] = tmp_normal[2]*edge->r[0] - edge->r[2]*tmp_normal[0];
    tmp_nrm[2] = tmp_normal[0]*edge->r[1] - edge->r[0]*tmp_normal[1];

    tmp_len = sqrt (tmp_nrm[0]*tmp_nrm[0] + tmp_nrm[1]*tmp_nrm[1] + tmp_nrm[2]*tmp_nrm[2]);

    // rgd
    if (tmp_len!=0)
    {
    for (i = 0 ; i < 3 ; i++)
      edge->nrm[i] = tmp_nrm[i] / tmp_len;
    }

    /* insert in rtree */
   RtreeInsert (Surf3DEdgeTree, (void *)edge, edge->min[0], edge->max[0],
                   edge->min[1], edge->max[1], edge->min[2], edge->max[2]);


#if PRINT_DEBUG
    fprintf (stdout, "\t*****Inserted Edge -> %d-%d = %f\n", edge->verts[0], edge->verts[1],
                     edge->length);
    fprintf (stdout, "\t*****Edge Normal -> %f %f %f\n", tmp_normal[0], tmp_normal[1],
                     tmp_normal[2]);
    fprintf (stdout, "\t*****Edge nrm -> %f %f %f\n", edge->nrm[0], edge->nrm[1],
                     edge->nrm[2]);
#endif

    /* markos - pushing based on layer, independently on parameter init, which now becomes unused */
    /*if (init == 0)
     Surf3DBdryPushSmall( edge ) ;
    else
     Surf3DBdryPushCorre( edge ) ;*/
    Surf3DBdryPushLayer( edge ) ;

    return(edge) ;
}


/*************************** Surf3DInsertNewNodesBar ************************/
static int Surf3DInsertNewNodesBar
(
Surf3DBdryEdge       *edge ,
int                 *n ,
double              *_new
)
{
  int               num_node ;
  double            dot, cand_vec[3];

  /* initiate number of nodes and allocation */
  num_node = *n;

  /* update node_list */
  num_node++;
  nintnode++;
  if (num_node > num_node_alloced)
  {
    num_node_alloced += Surf3D_NODE_QUANTUM;
    node_list = (Surf3DBdryNodeList *) Surf3DRealloc (node_list,
                num_node_alloced * sizeof (Surf3DBdryNodeRec));
  }
  node_list[num_node-1].coord[0] = _new[0];
  node_list[num_node-1].coord[1] = _new[1];
  node_list[num_node-1].coord[2] = _new[2];
  node_list[num_node-1].active_flag = 1;
  node_list[num_node-1].edges  = NULL;
  node_list[num_node-1].iedges = NULL;
  node_list[num_node-1].elems  = NULL;

  /* update number of nodes and allocation */
  *n = num_node;

  /* the new point should only be considered if it's in the same semi-plane
     than the base edge, there is, in the direction of base edge's normal */

  cand_vec[0] = _new[0] - edge->center[0] ;
  cand_vec[1] = _new[1] - edge->center[1] ;
  cand_vec[2] = _new[2] - edge->center[2] ;

  dot = cand_vec[0]*edge->nrm[0] +
        cand_vec[1]*edge->nrm[1] +
        cand_vec[2]*edge->nrm[2] ;

  if ( dot <= 0.0 ) return 0 ;
  else              return 1 ;
}

/*************************** Surf3DFindNode *********************************/
static void Surf3DFindNode
(
Surf3DBdryEdge       *edge ,
int                 debug ,
double              node[3]
)
{
  double             h ;

  /* get the polygon equilateral size of the same length of the edge */
  h = (edge->length * sqrt(3.0)) / 2.0;

  /* evaluate the new node coordinate */
  node[0] = edge->center[0] + (h/debug)*edge->nrm[0];
  node[1] = edge->center[1] + (h/debug)*edge->nrm[1];
  node[2] = edge->center[2] + (h/debug)*edge->nrm[2];
}

/*************************** Surf3DTreeNode *********************************/
static int Surf3DTreeNode
(
Surf3DBdryEdge      *edge ,
int                 debug ,
double              node[3] ,
double              *h0
)
{
 double  h, size, dt[3], dist;
 int     i, get;


 if (Surf3DSizeFunction == NULL)
 {
   /* get the size of the cell where the center of the edge is */
   h = edge->length;
 }
 else
 {
   (*Surf3DSizeFunction) (node_list[edge->verts[0]].coord[0],
                          node_list[edge->verts[0]].coord[1],
                          node_list[edge->verts[0]].coord[2], &size);
   h = size;
   (*Surf3DSizeFunction) (node_list[edge->verts[1]].coord[0],
                          node_list[edge->verts[1]].coord[1],
                          node_list[edge->verts[1]].coord[2],
                          &size);
   h += size;
   h /= 2.0 ;
 }

 if (h > (edge->length*1.5))
   h = edge->length*1.5;

 /* evaluate the new node coordinate */
 if (Surf3DIPointsFunction == NULL)
 {
   *h0=h;
    for (i = 0; i < 3; i++)
     node[i] = edge->center[i] + (h/debug)*edge->nrm[i];
 }
 else
 {
   double ipts[3], n[3];
   for (i = 0; i < 3; i++)
     n[i] = (node_list[edge->verts[0]].normal[i] + node_list[edge->verts[1]].normal[i]) * 0.5;

/*
   if (edge->verts[0] == 2633 && edge->verts[1] == 2634)
     printf ("Debug\n");
*/
   for (i = 0; i < 3; i++)
   {
     ipts[i] = edge->center[i] + (h/debug)*edge->nrm[i];
/*     try_pts[i] = edge->center[i] + (h/debug)*edge->nrm[i]; */
   }

   #if PRINT_DEBUG
     fprintf (stdout, "    Ponto tentativa (%f)= %f %f %f\n", h, ipts[0], ipts[1], ipts[2]);
     fprintf (stdout, "    Center Edge     (%f)= %f %f %f\n", h, edge->center[0], edge->center[1], edge->center[2]);
     fprintf (stdout, "    nrm Edge        (%f)= %f %f %f\n", h, edge->nrm[0], edge->nrm[1], edge->nrm[2]);
   #endif

   get = (*Surf3DIPointsFunction) (edge->center, n, edge->nrm, h, node);

   dt[0] = ipts[0] - node [0];
   dt[1] = ipts[1] - node [1];
   dt[2] = ipts[2] - node [2];
   dist = sqrt (dt[0]*dt[0] + dt[1]*dt[1] + dt[2]*dt[2]);

   if (!get || (dist > h*0.25))
   {
     h *= 0.5;
     get = (*Surf3DIPointsFunction) (edge->center, n, edge->nrm, h, node);
   }

   if (!get)
   {
     #if PRINT_DEBUG
        fprintf (stdout, "    Nao pegou ponto em Surf3DTreeNode -> %d %d\n", edge->verts[0], edge->verts[1]);
     #endif
     return 0;
   }

   #if PRINT_DEBUG
     fprintf (stdout, "    Ponto otimo (%f)= %f %f %f\n", h, node[0], node[1], node[2]);
   #endif

   *h0=h;

 }

 return 1;

}


/*************************** Surf3DGetMaxMin *********************************/
void Surf3DGetMaxMin (double p[3][3], double *xmax, double *xmin, double *ymax,
                      double *ymin, double *zmax, double *zmin, double scale)
{
  int    i;
  double center[3], delta[3], max;

  *xmin = *xmax = p[0][0];
  *ymin = *ymax = p[0][1];
  *zmin = *zmax = p[0][2];
  for (i = 1; i < 3; i++)
  {
    if (p[i][0] < *xmin)
      *xmin = p[i][0];
    if (p[i][0] > *xmax)
      *xmax = p[i][0];

    if (p[i][1] < *ymin)
      *ymin = p[i][1];
    if (p[i][1] > *ymax)
      *ymax = p[i][1];

    if (p[i][2] < *zmin)
      *zmin = p[i][2];
    if (p[i][2] > *zmax)
      *zmax = p[i][2];
  }
  /* center */
  center[0] = (*xmin + *xmax) * 0.5;
  center[1] = (*ymin + *ymax) * 0.5;
  center[2] = (*zmin + *zmax) * 0.5;
  delta[0] = (*xmax - *xmin) * 0.5;
  delta[1] = (*ymax - *ymin) * 0.5;
  delta[2] = (*zmax - *zmin) * 0.5;

  /* get the biggest size */
  if (delta[0] > delta[1])
  {
    if (delta[0] > delta[2])
      max = delta[0];
    else
      max = delta[2];
  }
  else
  {
    if (delta[1] > delta[2])
      max = delta[1];
    else
      max = delta[2];
  }


  /* apply scale */
  *xmin = center[0] - max * scale;
  *xmax = center[0] + max * scale;
  *ymin = center[1] - max * scale;
  *ymax = center[1] + max * scale;
  *zmin = center[2] - max * scale;
  *zmax = center[2] + max * scale;
}


/* -------------------------------------------------------------------
** Surf3DCheckValid - this routine makes geometrical checks to make sure
**                   that no edges of a candidate element cross the
**                   the current boundary.
*/
static int  Surf3DCheckValid
(
Surf3DBdryEdge     *edge,
int                node_indx,
int                found       /* the node was extract from the heap */
)
{
  int               i, j, type, edg_vts[3], *cur_elem;
  double            area, zero_area, len = edge->length*edge->length;
  double            tri_pts[3][3], n[3], vectors[2][3], base_matrix[3][3], *tmp_n;
  double            edge_nrm[3], dot_nrm, tol; /* markos - base edge and candidate node must point to the same side */
  double            xmax, xmin, ymax, ymin, zmax, zmin;
  double            xmx, xmn, ymx, ymn, zmx, zmn;
  Surf3DBdryEdge    *current = NULL;
  Surf3DAdjElem     *adj_elem = NULL; /* markos - new edges must not exist out of the front */
  Surf3DElem        *elem = NULL; /* markos - new edges must not exist out of the front */
  Surf3DAdjEdge     *adj_edge = NULL;
  int               edge_count[2] = {0,0}, edge_exists[2] = {0,0}; /* markos - if a new edge already exists */

  /*printf("markos - testing validity for edge %d %d, node %d\n", edge->verts[0]+1, edge->verts[1]+1, node_indx+1);*/

  /* first of all, check if the adjacent current front is a triangle */
  if (found)
  {
    if (Surf3DAdjEdgsFormTriangle (edge, node_indx))
      return 1;
  }

  /* edge vertices */
  edg_vts[0] = edge->verts[1];
  edg_vts[1] = node_indx;
  edg_vts[2] = edge->verts[0];

  #if PRINT_DEBUG
    fprintf (stdout, "    Funcao Surf3DCheckValid -> %d %d %d\n", edg_vts[0], edg_vts[1], edg_vts[2]);
  #endif

  /* check minimum area */
  zero_area = len * len / 500.0 / (edge->use+1);
  area = Surf3DAreaEdge_Node (edge, node_indx, node_list);
  if (area < zero_area && edge->use == 0)
  {
    /*printf("markos - debug almost zero area\n");*/
    #if PRINT_DEBUG
      fprintf (stdout, "    * Area quase zero  %f\n", area);
    #endif
     return 0;
  }

  /*if (area < 0.0) printf("markos - debug really zero area\n");*/

  /* really zero area */
  if (area < 0.0)
    return 0;

  /* markos - base edge and candidate node must have normals pointing to the same direction

     when a candidate node already existed, and the advancing front is on the geometric phase,
     the candidate node and the base edge must have normals pointing to the same direction,
     because a base edge cannot choose a node that's on the other side of a hill, for example.
    */
  if ((found > 0) && (edge->use == 0))
  {
      edge_nrm[0] = (node_list[edge->verts[0]].normal[0] + node_list[edge->verts[1]].normal[0]) * 0.5;
      edge_nrm[1] = (node_list[edge->verts[0]].normal[1] + node_list[edge->verts[1]].normal[1]) * 0.5;
      edge_nrm[2] = (node_list[edge->verts[0]].normal[2] + node_list[edge->verts[1]].normal[2]) * 0.5;

      dot_nrm = edge_nrm[0] * node_list[node_indx].normal[0] +
                edge_nrm[1] * node_list[node_indx].normal[1] +
                edge_nrm[2] * node_list[node_indx].normal[2];

      /*if (dot_nrm < 0.0) printf("markos - debug point on other side of the hill\n");*/

      if (dot_nrm < 0.0)
        return 0;
  }

  /* markos - new edges must not exist

     if one of the new edges already belongs to two elements, the new element
     cannot be added, because there would be an edge with 3 adjacent faces */
  for (adj_elem = node_list[node_indx].elems; adj_elem; adj_elem = adj_elem->next)
  {
    elem = &ptr_ElemList[adj_elem->id];

    for (i = 0; i < 3; i++)
    {
        if (elem->conn[i] == edge->verts[0]) edge_count[0]++;
        if (elem->conn[i] == edge->verts[1]) edge_count[1]++;
    }
  }

  if ((edge_count[0] >= 2) || (edge_count[1] >= 2))
  {
    /*printf("markos - debug edge %d %d, node %d\n", edge->verts[0]+1, edge->verts[1]+1, node_indx+1);*/
    return 0;
  }

  /* markos - checks if a new edge already exists; if it exists, it doesn't need to be checked
     for intersection with other existent edges */
  if (found > 0)
  {
    for (adj_edge = node_list[node_indx].edges; adj_edge; adj_edge = adj_edge->next)
    {
      if ((adj_edge->edge->verts[0] == edge->verts[0]) || (adj_edge->edge->verts[1] == edge->verts[0]))
        edge_exists[0] = 1; /* edge_exists[0] corresponds to edge->verts[0], which is also edg_vts[2] */
      if ((adj_edge->edge->verts[0] == edge->verts[1]) || (adj_edge->edge->verts[1] == edge->verts[1]))
        edge_exists[1] = 1; /* edge_exists[1] corresponds to edge->verts[1], which is also edg_vts[0] */

      if ((edge_exists[0]) && (edge_exists[1])) break;
    }
  }

  /* pts */
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
     tri_pts[j][i] = node_list[edg_vts[j]].coord[i];

    vectors[0][i] = tri_pts[1][i] - tri_pts[0][i];
    vectors[1][i] = tri_pts[2][i] - tri_pts[0][i];
  }

  /* normal */
  Surf3DCrossprodNorm (vectors[0], vectors[1], n);

  /* base matrix */
//  Surf3DGetMatrixTransfZ (n, base_matrix);
  base_matrix[2][0] = n[0];
  base_matrix[2][1] = n[1];
  base_matrix[2][2] = n[2];
  /* y axis */
  Surf3DCrossprodNorm (base_matrix[2], vectors[0], base_matrix[1]);
  /* x axis */
  Surf3DCrossprodNorm (base_matrix[1], base_matrix[2], base_matrix[0]);

  /* get boundbox and use a scale of amplification */
  Surf3DGetMaxMin (tri_pts, &xmax, &xmin, &ymax, &ymin, &zmax, &zmin, 2.0);

  /* find edges around element */
  RtreeInitSearchBox (Surf3DEdgeTree, xmin, xmax, ymin, ymax, zmin, zmax);

   /* loop through all the edges in the current boundbox and see
      if they intersect any of the edges of the proposed new element */
  while ((current = (Surf3DBdryEdge *)
          RtreeSearchBox(Surf3DEdgeTree, &xmn, &xmx, &ymn,
                                            &ymx, &zmn, &zmx)) != NULL)
  {
    int    cur_vts[2];
    double n_e[3], dot;
    /* markos - initializing type */
    type = 0;

    /* test if the edge exist */
    if (current == NULL)
      continue;

    /* test if the current edge. */
    if ( current == edge )
      continue;

    /* current vertices */
    cur_vts[0] = current->verts[0];
    cur_vts[1] = current->verts[1];

    /* internal edge */
    if (node_list[cur_vts[0]].active_flag == 0 &&
        node_list[cur_vts[1]].active_flag == 0 )
      continue;

    #if PRINT_DEBUG
      fprintf (stdout, "    * Testando edge %d %d\n", cur_vts[0], cur_vts[1]);
    #endif


    /* test if normal edges are compatible */
    for (i = 0; i < 3; i++)
      n_e[i] = (node_list[cur_vts[0]].normal[i] + node_list[cur_vts[1]].normal[i]) * 0.5;
    dot = n[0]*n_e[0] + n[1]*n_e[1] + n[2]*n_e[2];
    if (dot <= 0.0)
      continue;

    /* markos - if edge already existed, it doesn't need to be checked for intersection with other existent edges */
    if (!edge_exists[1])
    {
      /* test first edge */
      type = Surf3DCrossEdge2d (edg_vts[0], edg_vts[1], cur_vts[0], cur_vts[1], node_list, base_matrix);

      /*if (Surf3DCrossEdge2d (edg_vts[0], edg_vts[1], cur_vts[0], cur_vts[1], node_list, base_matrix) != Surf3DCrossEdge3d (edg_vts[0], edg_vts[1], cur_vts[0], cur_vts[1], node_list, base_matrix))
      {
        int cross = Surf3DCrossEdge2d (edg_vts[0], edg_vts[1], cur_vts[0], cur_vts[1], node_list, base_matrix);
        int cross3d = Surf3DCrossEdge3d (edg_vts[0], edg_vts[1], cur_vts[0], cur_vts[1], node_list, base_matrix);
        printf("markos - debug edge %d %d crossed edge %d %d A, cross %d, cross3d %d\n", edg_vts[0]+1, edg_vts[1]+1, cur_vts[0]+1, cur_vts[1]+1, cross, cross3d);
      }*/

      #if PRINT_DEBUG
        if (type == 1)
          fprintf (stdout, "    * Cruzou no  Surf3DCrossEdge 1\n");
      #endif
      /*if (type == 1) printf("markos - debug edge %d %d crossed edge %d %d A (2)\n", edg_vts[0]+1, edg_vts[1]+1, cur_vts[0]+1, cur_vts[1]+1);*/
      if (type == 1)
        return 0;
    }

    /* markos - if edge already existed, it doesn't need to be checked for intersection with other existent edges */
    if (!edge_exists[0])
    {
      /* test second edge */
      type = Surf3DCrossEdge2d (edg_vts[1], edg_vts[2], cur_vts[0], cur_vts[1], node_list, base_matrix);

      /*if (Surf3DCrossEdge2d (edg_vts[1], edg_vts[2], cur_vts[0], cur_vts[1], node_list, base_matrix) != Surf3DCrossEdge3d (edg_vts[1], edg_vts[2], cur_vts[0], cur_vts[1], node_list, base_matrix))
      {
        int cross = Surf3DCrossEdge2d (edg_vts[1], edg_vts[2], cur_vts[0], cur_vts[1], node_list, base_matrix);
        int cross3d = Surf3DCrossEdge3d (edg_vts[1], edg_vts[2], cur_vts[0], cur_vts[1], node_list, base_matrix);
        printf("markos - debug edge %d %d crossed edge %d %d B, cross %d, cross3d %d\n", edg_vts[1]+1, edg_vts[2]+1, cur_vts[0]+1, cur_vts[1]+1, cross, cross3d);
      }*/

      #if PRINT_DEBUG
        if (type == 1)
          fprintf (stdout, "    * Cruzou no  Surf3DCrossEdge 2\n");
      #endif
      /*if (type == 1) printf("markos - debug edge %d %d crossed edge %d %d B (2)\n", edg_vts[1]+1, edg_vts[2]+1, cur_vts[0]+1, cur_vts[1]+1);*/
      if (type == 1)
        return 0;
    }

    /* check if there is a point inside the new element */
    if (node_list[cur_vts[0]].active_flag)
    {
      tmp_n = node_list[cur_vts[0]].normal;
      dot = tmp_n[0]*n[0] + tmp_n[1]*n[1] + tmp_n[2]*n[2];
      if (dot > 0)
      {
        type = Surf3DInsideElem (tri_pts, node_list[cur_vts[0]].coord, base_matrix);

        /*if (Surf3DInsideElem (tri_pts, node_list[cur_vts[0]].coord, base_matrix) != Surf3DInsideElem3d (tri_pts, node_list[cur_vts[0]].coord, base_matrix))
        {
            int prox = Surf3DInsideElem (tri_pts, node_list[cur_vts[0]].coord, base_matrix);
            int prox3d = Surf3DInsideElem3d (tri_pts, node_list[cur_vts[0]].coord, base_matrix);
            printf("markos - debug triangle A %d %d %d, node %d, prox %d, prox3d %d\n", edge->verts[0]+1, edge->verts[1]+1, node_indx+1, cur_vts[0]+1, prox, prox3d);
        }*/

        /*if (type == 1) printf("markos - debug triangle A (2) %d %d %d, node %d\n", edge->verts[0]+1, edge->verts[1]+1, node_indx+1, cur_vts[0]+1);*/

        #if PRINT_DEBUG
          if( type == 1 )
            fprintf (stdout, "    * Point -%d- Inside of new Element\n", cur_vts[0]);
        #endif
        if (type == 1)
          return 0;
      }
    }


    /* check if there is a point inside the new element */
    if (node_list[cur_vts[1]].active_flag)
    {
      tmp_n = node_list[cur_vts[1]].normal;
      dot = tmp_n[0]*n[0] + tmp_n[1]*n[1] + tmp_n[2]*n[2];
      if (dot > 0)
      {
        type = Surf3DInsideElem (tri_pts, node_list[cur_vts[1]].coord, base_matrix);

        /*if (Surf3DInsideElem (tri_pts, node_list[cur_vts[1]].coord, base_matrix) != Surf3DInsideElem3d (tri_pts, node_list[cur_vts[1]].coord, base_matrix))
        {
            int prox = Surf3DInsideElem (tri_pts, node_list[cur_vts[1]].coord, base_matrix);
            int prox3d = Surf3DInsideElem3d (tri_pts, node_list[cur_vts[1]].coord, base_matrix);
            printf("markos - debug triangle B %d %d %d, node %d, prox %d, prox3d %d\n", edge->verts[0]+1, edge->verts[1]+1, node_indx+1, cur_vts[1]+1, prox, prox3d);
        }*/

        /*if (type == 1) printf("markos - debug triangle B (2) %d %d %d, node %d\n", edge->verts[0]+1, edge->verts[1]+1, node_indx+1, cur_vts[1]+1);*/

        #if PRINT_DEBUG
          if( type == 1 )
              fprintf (stdout, "    * Point -%d- Inside of new Element\n", cur_vts[1]);
        #endif
        if (type == 1)
          return 0;
      }
    }

    /* test to see if new point is near from current edge. If it's
       then this point should not be considered */
    if (found != 1)
    {
      type = Surf3DProxEdge (cur_vts[0], cur_vts[1], edg_vts[1], (edge->use*2 + 1)*1.0, base_matrix);

      /*if (Surf3DProxEdge (cur_vts[0], cur_vts[1], edg_vts[1], (edge->use*2 + 1)*1.0, base_matrix) != Surf3DProxEdge3d (cur_vts[0], cur_vts[1], edg_vts[1], (edge->use*2 + 1)*1.0, base_matrix))
      {
          int prox = Surf3DProxEdge (cur_vts[0], cur_vts[1], edg_vts[1], (edge->use*2 + 1)*1.0, base_matrix);
          int prox3d = Surf3DProxEdge3d (cur_vts[0], cur_vts[1], edg_vts[1], (edge->use*2 + 1)*1.0, base_matrix);
          printf("markos - debug edge A %d %d, node %d, prox %d, prox3d %d\n", cur_vts[0]+1, cur_vts[1]+1, edg_vts[1]+1, prox, prox3d);
      }*/

      /*if (type == 1) printf("markos - debug edge A (2) %d %d, node %d\n", cur_vts[0]+1, cur_vts[1]+1, edg_vts[1]+1);*/

      #if PRINT_DEBUG
        if (type == 1)
          fprintf (stdout, "    * Falhou no  Surf3DProxEdge\n");
      #endif

      if( type == 1 )
          return 0;
    }

    /* test to see if the new edges are near of the current edge */
    if (edge->use < 3) /* markos - new use for edge */
    {
      for (i = 0; i < 2; ++i)
      {
        /* markos - if edge already existed, it doesn't need to be checked for intersection with other existent edges */
        if (edge_exists[1-i]) continue;

        for (j = 0; j < 2; ++j)
        {
          type = Surf3DProxEdge (edg_vts[i+1], edg_vts[i], cur_vts[j], (edge->use*2 + 1)*1.0, base_matrix);

          /*if (Surf3DProxEdge (edg_vts[i+1], edg_vts[i], cur_vts[j], (edge->use*2 + 1)*1.0, base_matrix) != Surf3DProxEdge3d (edg_vts[i+1], edg_vts[i], cur_vts[j], (edge->use*2 + 1)*1.0, base_matrix))
          {
              int prox = Surf3DProxEdge (edg_vts[i+1], edg_vts[i], cur_vts[j], (edge->use*2 + 1)*1.0, base_matrix);
              int prox3d = Surf3DProxEdge3d (edg_vts[i+1], edg_vts[i], cur_vts[j], (edge->use*2 + 1)*1.0, base_matrix);
              printf("markos - debug edge B %d %d, node %d, prox %d, prox3d %d\n", edg_vts[i+1]+1, edg_vts[i]+1, cur_vts[j]+1, prox, prox3d);
          }*/

          /*if (type == 1) printf("markos - debug edge B (2) %d %d, node %d\n", edg_vts[i+1]+1, edg_vts[i]+1, cur_vts[j]+1);*/

          #if PRINT_DEBUG
            if (type == 1)
              fprintf (stdout, "    * Falhou no  Surf3DProxEdge\n");
          #endif

          if( type == 1 )
             return 0;
        }
      }
    }

  }

  /* find elements around the tryed element */
  RtreeInitSearchBox (Surf3DElemTree, xmin, xmax, ymin, ymax, zmin, zmax);
  while ((cur_elem = (int *)
          RtreeSearchBox(Surf3DElemTree, &xmn, &xmx, &ymn,
                                            &ymx, &zmn, &zmx)) != NULL)
  {
    double dot, pts[3][3], vec[3], normal[3], tri_base_matrix[3][3];
    int    idx[3];
    /* markos - initializing type */
    type = 0;

    /* markos - triangle x triangle intersection tests; they are performed
       for all existing triangles (in the search region), independently whether they
       have normals compatible with the new triangle */
    type = Surf3DCrossTri3d(edge->verts[0], edge->verts[1], node_indx,
        ptr_ElemList[*cur_elem].conn[0], ptr_ElemList[*cur_elem].conn[1], ptr_ElemList[*cur_elem].conn[2],
        node_list, n, ptr_ElemList[*cur_elem].normal, base_matrix);

    /*if (type == 1) printf("markos - debug triangle intersection %d %d %d\n", ptr_ElemList[*cur_elem].conn[0]+1, ptr_ElemList[*cur_elem].conn[1]+1, ptr_ElemList[*cur_elem].conn[2]+1);*/

    if (type == 1)
        return 0;

    /* markos - if node already existed, i.e., it is not a new node, it doesn't
       need to be checked for proximity of other triangles */
    if (found != 1)
    {
      dot = n[0] * ptr_ElemList[*cur_elem].normal[0] +
            n[1] * ptr_ElemList[*cur_elem].normal[1] +
            n[2] * ptr_ElemList[*cur_elem].normal[2];

      if (dot <= 0.0)
        continue;

      for (i = 0; i < 3; i++)
      {
        idx[i] = ptr_ElemList[*cur_elem].conn[i];
        for (j = 0; j < 3; j++)
          pts[i][j] = node_list[idx[i]].coord[j];
      }

      /* markos - using the base matrix of the existing triangle, not the new triangle (as was before) */
      for (i = 0; i < 3; i++)
      {
        vec[i] = pts[1][i] - pts[0][i];
        normal[i] = ptr_ElemList[*cur_elem].normal[i];
      }

      tri_base_matrix[2][0] = normal[0];
      tri_base_matrix[2][1] = normal[1];
      tri_base_matrix[2][2] = normal[2];
      Surf3DCrossprodNorm (tri_base_matrix[2], vec, tri_base_matrix[1]);
      Surf3DCrossprodNorm (tri_base_matrix[1], tri_base_matrix[2], tri_base_matrix[0]);

      type = Surf3DInsideElem (pts, node_list[node_indx].coord, tri_base_matrix);

      /*if (Surf3DInsideElem (pts, node_list[node_indx].coord, tri_base_matrix) != Surf3DInsideElem3d (pts, node_list[node_indx].coord, tri_base_matrix))
      {
        int prox = Surf3DInsideElem (pts, node_list[node_indx].coord, tri_base_matrix);
        int prox3d = Surf3DInsideElem3d (pts, node_list[node_indx].coord, tri_base_matrix);
        printf("markos - debug triangle C %d %d %d, node %d, prox %d, prox3d %d\n", idx[0]+1, idx[1]+1, idx[2]+1, node_indx+1, prox, prox3d);
      }*/

      /*if (type == 1) printf("markos - debug triangle C (2) %d %d %d, node %d\n", idx[0]+1, idx[1]+1, idx[2]+1, node_indx+1);*/

#if PRINT_DEBUG
    fprintf (stdout, "    * Creck Surf3DInsideElem %d-%d-%d\n", idx[0], idx[1], idx[2]);
    fprintf (stdout, "    * Creck Surf3DInsideElem => %d\n", type);
      if( type == 1 )
        fprintf (stdout, "    * New Point Inside of old Elements\n\n");
#endif
    if (type == 1)
      return 0;
    }
  }

  /* markos - collinearity test is being performed here now */
  if (found > 0)
  {
    type = colinear(node_list[edge->verts[0]].coord[0], node_list[edge->verts[0]].coord[1], node_list[edge->verts[0]].coord[2],
                    node_list[node_indx].coord[0],      node_list[node_indx].coord[1],      node_list[node_indx].coord[2],
                    node_list[edge->verts[1]].coord[0], node_list[edge->verts[1]].coord[1], node_list[edge->verts[1]].coord[2]);

    /* markos other type of collinearity test, disabled now */
    /*Surf3DCrossprod3d (vectors[0], vectors[1], n);

    len = fabs(n[0]) + fabs(n[1]) + fabs(n[2]);

    tol = MIN(Surf3D_toler, edge->length)*0.01;

    type = (len < tol);*/

    /*{
        Surf3DCrossprod3d (vectors[0], vectors[1], n);

        len = fabs(n[0]) + fabs(n[1]) + fabs(n[2]);

        tol = MIN(Surf3D_toler, edge->length)*0.01;

        if ((type == 1) && (len < tol))
        {
            printf("markos - debug collinear on both tests\n");
        }
        else if (type == 1)
        {
            printf("markos - debug collinear on colinear() test\n");
        }
        else if (len < tol)
        {
            printf("markos - debug collinear on MIN() test\n");
        }
    }*/

    if (type == 1)
      return 0;
  }

  /* if we get here, the new element is valid */
  return 1;
}

/* -------------------------------------------------------------------
** Surf3DAdjEdgsFormTriangle - check if the adjacent current front is a triangle.
*/
static int Surf3DAdjEdgsFormTriangle
(
Surf3DBdryEdge  *edge,
int             node_indx
)
{
  Surf3DAdjEdgeList *adj_edge;
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


/*************************** Surf3DCrossTri3d *********************************/
/* markos - triangle x triangle intersection tests (based on Msh3D) */
static int Surf3DCrossTri3d (int a, int b, int c, int u, int v, int w, Surf3DBdryNode *nodes, double *n1, double *n2, double base_matrix[3][3]) //markos
{
    /* test intersection against two given triangles and return if
       intersection that it's not coincident with any vertices is
       found */

    int i, j, type = 0;
    double box1[2], box2[2], len[2][3], d, dot, tol;

    int p[3] = {1, 0, 0};
    int q[3] = {0, 0, 0};

    int tri1[3] = {a, b, c};
    int tri2[3] = {u, v, w};

    for (i = 0; i < 3; i++)
    {
        box1[0] = box1[1] = nodes[tri1[0]].coord[i];
        box2[0] = box2[1] = nodes[tri2[0]].coord[i];

        for (j = 1; j < 3; j++)
        {
            box1[0] = MIN(box1[0], nodes[tri1[j]].coord[i]);
            box1[1] = MAX(box1[1], nodes[tri1[j]].coord[i]);

            box2[0] = MIN(box2[0], nodes[tri2[j]].coord[i]);
            box2[1] = MAX(box2[1], nodes[tri2[j]].coord[i]);
        }

        if ((box1[0] > box2[1]) || (box2[0] > box1[1]))
        {
            return 0;
        }
    }

    for (i = 0; i < 3; i++)
    {
        len[0][i] = 0.0;
        len[1][i] = 0.0;

        for (j = 0; j < 3; j++)
        {
            d = nodes[tri1[i]].coord[j] - nodes[tri1[(i+1)%3]].coord[j];
            len[0][i] += d*d;

            d = nodes[tri2[i]].coord[j] - nodes[tri2[(i+1)%3]].coord[j];
            len[1][i] += d*d;
        }

        len[0][i] = sqrt(len[0][i]);
        len[1][i] = sqrt(len[1][i]);
    }

    len[0][0] = MAX(len[0][0], MAX(len[0][1], len[0][2]));
    len[1][0] = MAX(len[1][0], MAX(len[1][1], len[1][2]));
    len[0][0] = MIN(len[0][0], len[1][0]);

    tol = 1.0e-08 * MIN(Surf3D_toler, len[0][0]/10.0);

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            if (tri1[i] == tri2[j])
            {
                p[i] = 1;
                p[(i+2)%3] = 1;

                q[j] = 1;
                q[(j+2)%3] = 1;
            }
        }
    }

    type = GeoTriIntersect (nodes[a].coord, nodes[b].coord, nodes[c].coord, p, nodes[u].coord, nodes[v].coord, nodes[w].coord, q, tol);

    if (type == 0) return 1;

    /* if we get here no intersection was found */

    return 0;
}


/*************************** Surf3DCrossEdge3d *********************************/
/* markos - same as Surf3DCrossEdge2d, but tests 3D bbox before it (currently unused) */
static int Surf3DCrossEdge3d
(
int                  I1,           /* Node pointer   */
int                  I2,           /* Node pointer   */
int                  J1,           /* Node pointer   */
int                  J2,           /* Node pointer   */
Surf3DBdryNode       *nodes,
double               BaseMatrix[3][3]
)
{
  /* locals variables */
  int    cross, i;
  double PI1[3], PI2[3], PJ1[3], PJ2[3];
  double orig[3] = {0.0, 0.0, 0.0};
  double v[3], d[2], tol;

  for (i = 0; i < 3; i++)
    orig[i] = nodes[I1].coord[i];

  Surf3DTransf_xyz (orig, BaseMatrix, nodes[I1].coord, PI1);
  Surf3DTransf_xyz (orig, BaseMatrix, nodes[I2].coord, PI2);
  Surf3DTransf_xyz (orig, BaseMatrix, nodes[J1].coord, PJ1);
  Surf3DTransf_xyz (orig, BaseMatrix, nodes[J2].coord, PJ2);


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

  //begin markos
  v[0] = PI2[0] - PI1[0];
  v[1] = PI2[1] - PI1[1];
  v[2] = PI2[2] - PI1[2];
  d[0] = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] = PJ2[0] - PJ1[0];
  v[1] = PJ2[1] - PJ1[1];
  v[2] = PJ2[2] - PJ1[2];
  d[1] = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  //d[0] = MAX(d[0], d[1]);
  d[0] = MIN(d[0], d[1]);
  //tol = MIN(d[0]/10.0, Surf3D_toler);
  tol = MIN(d[0]/10.0, 10.0*Surf3D_toler);
  //end markoss

  /*  First do some simple box checks to eliminate lines that are not
      close together */
  if (( (PI1[0] > PJ1[0] + tol) && (PI1[0] > PJ2[0] + tol) && (PI2[0] > PJ1[0] + tol) && (PI2[0] > PJ2[0] + tol) ) ||
      ( (PI1[0] < PJ1[0] - tol) && (PI1[0] < PJ2[0] - tol) && (PI2[0] < PJ1[0] - tol) && (PI2[0] < PJ2[0] - tol) ))
  {
    cross = FALSE;
    return(cross);
  }

  if (( (PI1[1] > PJ1[1] + tol) && (PI1[1] > PJ2[1] + tol) && (PI2[1] > PJ1[1] + tol) && (PI2[1] > PJ2[1] + tol) ) ||
      ( (PI1[1] < PJ1[1] - tol) && (PI1[1] < PJ2[1] - tol) && (PI2[1] < PJ1[1] - tol) && (PI2[1] < PJ2[1] - tol) ))
  {
    cross = FALSE;
    return(cross);
  }

  if (( (PI1[2] > PJ1[2] + tol) && (PI1[2] > PJ2[2] + tol) && (PI2[2] > PJ1[2] + tol) && (PI2[2] > PJ2[2] + tol) ) ||
      ( (PI1[2] < PJ1[2] - tol) && (PI1[2] < PJ2[2] - tol) && (PI2[2] < PJ1[2] - tol) && (PI2[2] < PJ2[2] - tol) ))
  {
    cross = FALSE;
    return(cross);
  }

  /*  Now cross line I with J1 and line I with J2, if have the same sign they
     cannot cross */
  if (Surf3DCrossprod2d (PI1, PI2, PJ1) * Surf3DCrossprod2d (PI1, PI2, PJ2) >= 0.0)
  {
    cross = FALSE;
    return(cross);
  }

  /*  Now cross line J with I1 and line J with I2, if have the same sign they
      cannot cross */
  if (Surf3DCrossprod2d (PJ1, PJ2, PI1) * Surf3DCrossprod2d (PJ1, PJ2, PI2) >= 0.0)
  {
    cross = FALSE;
    return(cross);
  }

  return(cross);
}


/*************************** Surf3DProxEdge *********************************/
static int Surf3DProxEdge
(
int    edge_id_i,
int    edge_id_j,
int    test_pts,
double factor,
double base[3][3]
)
{
  double   dot1, dot2,len, tol; /* markos - added tolerance */
  double   orig[3], pts_j[2], pts_test[2];
  int      i;

  /* to avoid test with the same nodes  */
  if (edge_id_i == test_pts || edge_id_j == test_pts)
    return 0;

  /* origin */
  for (i = 0; i < 3; i++)
    orig[i] = node_list[edge_id_i].coord[i];

  /* transform to 2D space */
  Surf3DTransf_xy (orig, base, node_list[edge_id_j].coord, pts_j);
  Surf3DTransf_xy (orig, base, node_list[test_pts].coord , pts_test);

  /* obs: because edge_id_i is origin, pts_j and pts_test are vector too */
  len = sqrt (pts_j[0]*pts_j[0] + pts_j[1]*pts_j[1]);

  if (len == 0.0)
    return 0;

  /* check if the point is in the domain of the line */
  dot1 = (pts_j[0]*pts_test[0] + pts_test[1]*pts_j[1]) / len;
  if (dot1 < -0.01*len || dot1 > 1.01*len)
    return 0;

  /* check the correct side */
  dot2 = (pts_j[0]*pts_test[1] - pts_test[0]*pts_j[1]) / len;
  if (dot2 < 0.0) /* opposite side */
    return 0;

  /* markos - added tolerance, which might depend on global tolerance */
  tol = MIN(len, 100.0*Surf3D_toler);

  if (dot2 < tol * BDY_FACTOR / factor) /* markos - changed len for tol */
    return 1;

 return 0;
}



/*************************** Surf3DProxEdge3d ******************************/
/* markos - same as Surf3DProxEdge, but tests 3D bbox before it (currently unused) */
static int Surf3DProxEdge3d
(
int    edge_id_i,
int    edge_id_j,
int    test_pts,
double factor,
double base[3][3]
)
{
  double   dot1, dot2,len, tol; /* markos - added tolerance */
  double   orig[3], pts_i[3], pts_j[3], pts_test[3];
  int      i;

  /* to avoid test with the same nodes  */
  if (edge_id_i == test_pts || edge_id_j == test_pts)
    return 0;

  /* origin */
  for (i = 0; i < 3; i++)
    orig[i] = node_list[edge_id_i].coord[i];

  /* transform to 2D space */
  Surf3DTransf_xyz (orig, base, node_list[edge_id_i].coord, pts_i); /* markos - changes the 3D coordinate base */
  Surf3DTransf_xyz (orig, base, node_list[edge_id_j].coord, pts_j); /* markos - changes the 3D coordinate base */
  Surf3DTransf_xyz (orig, base, node_list[test_pts].coord , pts_test); /* markos - changes the 3D coordinate base */

  /* obs: because edge_id_i is origin, pts_j and pts_test are vector too */
  len = sqrt (pts_j[0]*pts_j[0] + pts_j[1]*pts_j[1]);

  if (len == 0.0)
    return 0;

  /* markos - added tolerance, which might depend on global tolerance */
  tol = MIN(len, 100.0*Surf3D_toler);

  /* markos - tests if the tested point is outside the bounding box of the edge */
  for (i = 0; i < 3; i++)
  {
    if (((pts_test[i] > pts_i[i] + tol) && (pts_test[i] > pts_j[i] + tol)) ||
        ((pts_test[i] < pts_i[i] - tol) && (pts_test[i] < pts_j[i] - tol)))
    {
      return 0;
    }
  }

  /* check if the point is in the domain of the line */
  dot1 = (pts_j[0]*pts_test[0] + pts_test[1]*pts_j[1]) / len;
  if (dot1 < -0.01*len || dot1 > 1.01*len)
    return 0;

  /* check the correct side */
  dot2 = (pts_j[0]*pts_test[1] - pts_test[0]*pts_j[1]) / len;
  if (dot2 < 0.0) /* opposite side */
    return 0;

  if (dot2 < tol * BDY_FACTOR / factor) /* markos - changed len for tol */
    return 1;

 return 0;
}


/*************************** Surf3DInisideElem ******************************/
static int Surf3DInsideElem (double pts[3][3], double pts_idx[3], double base[3][3])
{
  int    i, count = 0;
  double orig[3], P[3][2], Pj[2];
  double *Pi, *Pk;
  double d1, d2, d3, d, tol; /* markos - added tolerance */

  /* origin */
  for (i = 0; i < 3; i++)
    orig[i] = pts[0][i];

  /* transform element vert to 2d */
  for (i = 0; i < 3; i++)
    Surf3DTransf_xy (orig, base, pts[i], P[i]);

  /* tranform node_indx to 2d */
  Surf3DTransf_xy (orig, base, pts_idx, Pj);

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

/*
  fprintf (stdout, "%lf %lf\n", P[0][0], P[0][1]);
  fprintf (stdout, "%lf %lf\n", P[1][0], P[1][1]);
  fprintf (stdout, "%lf %lf\n", P[2][0], P[2][1]);
  fprintf (stdout, "%lf %lf\n", Pj[0], Pj[1]);
*/

  /* markos - added tolerance, which depends on the triangle's edges' length */
  d1 = sqrt((P[0][0]-P[1][0])*(P[0][0]-P[1][0]) + (P[0][1]-P[1][1])*(P[0][1]-P[1][1]));
  d2 = sqrt((P[1][0]-P[2][0])*(P[1][0]-P[2][0]) + (P[1][1]-P[2][1])*(P[1][1]-P[2][1]));
  d3 = sqrt((P[2][0]-P[0][0])*(P[2][0]-P[0][0]) + (P[2][1]-P[0][1])*(P[2][1]-P[0][1]));
  d = MAX(d1, MAX(d2, d3));
  tol = MIN(d/10.0, Surf3D_toler);

  /* markos - changed Surf3D_toler for tol */
  if ((fabs (P[0][0] - Pj[0])) < tol &&
      (fabs (P[0][1] - Pj[1])) < tol )
    return 0;

  if ((fabs (P[1][0] - Pj[0])) < tol &&
      (fabs (P[1][1] - Pj[1])) < tol )
    return 0;

  if ((fabs (P[2][0] - Pj[0])) < tol &&
      (fabs (P[2][1] - Pj[1])) < tol )
    return 0;

  for (i = 0; i < 3; i++)
  {
    Pi = P[i];
    Pk = P[(i+1)%3];
    if (Surf3DCrossprod2d (Pi, Pj, Pk) >= 0.0)
      count++;
  }

  if (count == 3) return 1;

  return 0;
}



/*************************** Surf3DInisideElem3d ***************************/
/* markos - same as Surf3DInsideElem, but tests 3D bbox before it (currently unused) */
static int Surf3DInsideElem3d (double pts[3][3], double pts_idx[3], double base[3][3])
{
  int    i, count = 0;
  double orig[3], P[3][3], Pj[3];
  double *Pi, *Pk;
  double d1, d2, d3, d, tol; /* markos - added tolerance */

  /* origin */
  for (i = 0; i < 3; i++)
    orig[i] = pts[0][i];

  /* transform element vert to 2d */
  for (i = 0; i < 3; i++)
    Surf3DTransf_xyz (orig, base, pts[i], P[i]); /* markos - changes the 3D coordinate base */

  /* tranform node_indx to 2d */
  Surf3DTransf_xyz (orig, base, pts_idx, Pj); /* markos - changes the 3D coordinate base */

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

/*
  fprintf (stdout, "%lf %lf\n", P[0][0], P[0][1]);
  fprintf (stdout, "%lf %lf\n", P[1][0], P[1][1]);
  fprintf (stdout, "%lf %lf\n", P[2][0], P[2][1]);
  fprintf (stdout, "%lf %lf\n", Pj[0], Pj[1]);
*/

  /* markos - added tolerance, which depends on the triangle's edges' length */
  d1 = sqrt((P[0][0]-P[1][0])*(P[0][0]-P[1][0]) + (P[0][1]-P[1][1])*(P[0][1]-P[1][1]) + (P[0][2]-P[1][2])*(P[0][2]-P[1][2]));
  d2 = sqrt((P[1][0]-P[2][0])*(P[1][0]-P[2][0]) + (P[1][1]-P[2][1])*(P[1][1]-P[2][1]) + (P[1][2]-P[2][2])*(P[1][2]-P[2][2]));
  d3 = sqrt((P[2][0]-P[0][0])*(P[2][0]-P[0][0]) + (P[2][1]-P[0][1])*(P[2][1]-P[0][1]) + (P[2][2]-P[0][2])*(P[2][2]-P[0][2]));
  d = MAX(d1, MAX(d2, d3));
  tol = MIN(d/10.0, Surf3D_toler);

  /* markos - tests if the tested point is outside the bounding box of the triangle */
  for (i = 0; i < 3; i++)
  {
    if (((Pj[i] > P[0][i] + tol) && (Pj[i] > P[1][i] + tol) && (Pj[i] > P[2][i] + tol)) ||
        ((Pj[i] < P[0][i] - tol) && (Pj[i] < P[1][i] - tol) && (Pj[i] < P[2][i] - tol)))
    {
      return 0;
    }
  }

  /* markos - changed Surf3D_toler for tol */
  if ((fabs (P[0][0] - Pj[0])) < tol &&
      (fabs (P[0][1] - Pj[1])) < tol &&
      (fabs (P[0][2] - Pj[2])) < tol )
    return 0;

  if ((fabs (P[1][0] - Pj[0])) < tol &&
      (fabs (P[1][1] - Pj[1])) < tol &&
      (fabs (P[1][2] - Pj[2])) < tol )
    return 0;

  if ((fabs (P[2][0] - Pj[0])) < tol &&
      (fabs (P[2][1] - Pj[1])) < tol &&
      (fabs (P[2][2] - Pj[2])) < tol )
    return 0;

  for (i = 0; i < 3; i++)
  {
    Pi = P[i];
    Pk = P[(i+1)%3];
    if (Surf3DCrossprod2d (Pi, Pj, Pk) >= 0.0)
      count++;
  }

  if (count == 3) return 1;

  return 0;
}


/* -------------------------------------------------------------------
** Surf3DDetNode  - this routine determines the node that should be
**                 taken in case of two nodes choosen has the same
**                 coordinates.
*/

static int Surf3DDetNode
(
Surf3DBdryEdge      *edge,
int                node_indx,
int                next_indx
)
{
  int             num_adj;
  double          vector[3];
  Surf3DAdjIniEdge *cur ;
  double          med_nrm[3] = {0.0, 0.0, 0.0}, dot1, dot2;

  /* vector from center to node_indx */
  vector[0] = edge->center[0] - node_list[node_indx].coord[0];
  vector[1] = edge->center[1] - node_list[node_indx].coord[1];
  vector[2] = edge->center[2] - node_list[node_indx].coord[2];

  /* test the first node => node_indx */
  num_adj = 0;
  for( cur = node_list[node_indx].iedges ; cur ; cur = cur->next )
  {
    med_nrm[0] += cur->nrm[0];
    med_nrm[1] += cur->nrm[1];
    med_nrm[2] += cur->nrm[2];
    num_adj++;
  }
  med_nrm[0] /= num_adj;
  med_nrm[1] /= num_adj;
  med_nrm[2] /= num_adj;
  dot1 = vector[0]*med_nrm[0] + vector[1]*med_nrm[1] + vector[2]*med_nrm[2];


  /* test the second node => next_indx */
  num_adj = 0;
  med_nrm[0] = med_nrm[1] = med_nrm[2] = 0.0;
  for( cur = node_list[next_indx].iedges ; cur ; cur = cur->next )
  {
    med_nrm[0] += cur->nrm[0];
    med_nrm[1] += cur->nrm[1];
    med_nrm[2] += cur->nrm[2];
    num_adj++;
  }
  med_nrm[0] /= num_adj;
  med_nrm[1] /= num_adj;
  med_nrm[2] /= num_adj;
  dot2 = vector[0]*med_nrm[0] + vector[1]*med_nrm[1] + vector[2]*med_nrm[2];

  if (dot1 > dot2)
    return node_indx;
  else
    return next_indx;

  return -1 ;
}

/* -------------------------------------------------------------------
** Surf3DAddEdges - this routine updates the edge stack to include/
**                 delete the necessary boundary edges to include
**                 the new element.  It also updates the active
**                 flags in the node list.
*/

static void Surf3DAddEdges
(
Surf3DBdryEdge      *edge,
int                 node_indx
)
{
    int             v1, i, v[3];
    int             found ;
    Surf3DAdjEdge   *cur, **temp, *save ;
    Surf3DBdryEdge  *edge_to_delete, *_new ;
    int             layer[3]; /* markos - added parameter for new edges' layer */

    /* start with each of the base nodes.  See if they are adjacent
       to a edge that includes the cap node along with the previous
       vertex on the base.  If so, this edge is already part of the
       boundary, and it should be deleted. */

    v[0] = edge->verts[0];
    v[1] = node_indx;
    v[2] = edge->verts[1];

    /* markos - layer of the new edges depends only on the layer of the
       base edge (currently disabled) */
    /*layer[0] = layer[1] = layer[2] = edge->layer + 1;*/

    /*markos - layer of a new edge depends on the minimum layer of an edge
      adjacent to any of its nodes (default) */
    for (i = 0; i < 3; i++)
    {
        layer[i] = edge->layer;

        for (cur = node_list[v[i]].edges; cur; cur = cur->next)
        {
            if (cur->edge->layer < layer[i])
            {
                layer[i] = cur->edge->layer;
            }
        }

        ++layer[i];
    }

    layer[0] = MIN(layer[0], layer[1]);
    layer[1] = MIN(layer[1], layer[2]);

    /* markos - layer of a new edge depends on the minimum layer of the
       edges forming the new triangle */
    /*layer[0] = layer[1] = edge->layer;

    for (cur = node_list[v[1]].edges; cur; cur = cur->next)
    {
        if ((cur->edge->verts[0] == v[0]) || (cur->edge->verts[1] == v[0]))
        {
            layer[1] = MIN(layer[1], cur->edge->layer);
        }

        if ((cur->edge->verts[0] == v[2]) || (cur->edge->verts[1] == v[2]))
        {
            layer[0] = MIN(layer[0], cur->edge->layer);
        }
    }

    ++layer[0];
    ++layer[1];*/

    /* test if there is a existent edge */
    for (i = 0; i < 2 ; i++)
    {
      v1 = edge->verts[i];
      found = 0;
      for ( cur=node_list[v1].edges ; cur ; cur=cur->next )
      {
#if PRINT_DEBUG
  fprintf (stdout, "Test edges to delete => %d-%d e %d-%d\n", cur->edge->verts[1], cur->edge->verts[0], v[i], v[i+1]);
#endif
        if ( (cur->edge->verts[0] == v[i+1]) && (cur->edge->verts[1] == v[i]) )
        {
#if PRINT_DEBUG
  fprintf (stdout, "\t*****Delete Edge => %d-%d\n", v[i], v[i+1]);
#endif
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
        for ( temp = &(node_list[v1].edges) ; *temp ;
              temp = &((*temp)->next) )
        {
          if ( (*temp)->edge == edge_to_delete )
          {
            save = *temp ;
            *temp = (*temp)->next ;
            Surf3DAdjFree( save ) ;
            break ;
          }
        }

        for ( temp = &(node_list[node_indx].edges) ; *temp ;
              temp = &((*temp)->next) )
        {
          if ( (*temp)->edge == edge_to_delete )
          {
            save = *temp ;
            *temp = (*temp)->next ;
            Surf3DAdjFree( save ) ;
            break ;
          }
        }

        Surf3DBdryDelete( cur->edge ) ;
      }
      else
      {
         /* if we did not find the edge, we must create it and update
          all the adjacent lists */

        _new = Surf3DPushBdryEdge( node_list, v[i], v[i+1], 0, 1, layer[i]) ; /* markos - changed use for 0 (always); added parameter for layer */

        cur = Surf3DAdjEdgeAlloc() ;
        cur->edge = _new ;
        cur->next = node_list[v[i]].edges ;
        node_list[v[i]].edges = cur ;

        cur = Surf3DAdjEdgeAlloc() ;
        cur->edge = _new ;
        cur->next = node_list[v[i+1]].edges ;
        node_list[v[i+1]].edges = cur ;
      }
    }

    /* update all the adjacent edge lists so they no longer point
       to the base edge */

    for (i = 0; i < 2 ; i++)
    {
      for (temp = &(node_list[edge->verts[i]].edges) ; *temp; temp = &((*temp)->next))
      {
        if ( (*temp)->edge == edge )
        {
          save = *temp ;
          *temp = (*temp)->next ;
          Surf3DAdjFree( save ) ;
#if PRINT_DEBUG
    fprintf (stdout, "\t*****Release Adj Edge => %d-%d from node %d\n",
             edge->verts[0], edge->verts[1], edge->verts[i]);
#endif
          break ;
        }
      }
    }

    /* update the active flags.  Any node on the new element that
       no longer adjacent to at least on edge on the boundary
       becomes inactive */

    for ( i=0 ; i<2 ; i++ )
    {
      if ( node_list[edge->verts[i]].edges == NULL)
           node_list[edge->verts[i]].active_flag = 0 ;
    }
    if ( node_list[node_indx].edges == NULL)
     node_list[node_indx].active_flag = 0 ;
    else
     node_list[node_indx].active_flag = 1 ;

#if PRINT_DEBUG
    fprintf (stdout, "\t*****Delete Curr Edge => %d-%d\n", edge->verts[0], edge->verts[1]);
#endif

    /* markos - remove edge from rtree */
    RtreeDelete (Surf3DEdgeTree, edge, edge->min[0], edge->max[0],
                   edge->min[1], edge->max[1], edge->min[2], edge->max[2]);

    /* give up the memory associated with this edge */
    Surf3DBdryFree( edge ) ;
}

/* -------------------------------------------------------------------
** Surf3DAddElemRTree - add a new element to the RTree.
*/
static void Surf3DAddElemRTree (Surf3DElemList *elemlist, int id_elem)
{
  double xmin, xmax, ymin, ymax, zmin, zmax;
  int    i;
  int    *new_id = (int *) malloc (sizeof (int));


  xmin = xmax = node_list[elemlist->conn[0]].coord[0];
  ymin = ymax = node_list[elemlist->conn[0]].coord[1];
  zmin = zmax = node_list[elemlist->conn[0]].coord[2];
  for (i = 1; i < 3; i++)
  {
    if (node_list[elemlist->conn[i]].coord[0] < xmin)
      xmin = node_list[elemlist->conn[i]].coord[0];
    if (node_list[elemlist->conn[i]].coord[0] > xmax)
      xmax = node_list[elemlist->conn[i]].coord[0];

    if (node_list[elemlist->conn[i]].coord[1] < ymin)
      ymin = node_list[elemlist->conn[i]].coord[1];
    if (node_list[elemlist->conn[i]].coord[1] > ymax)
      ymax = node_list[elemlist->conn[i]].coord[1];

    if (node_list[elemlist->conn[i]].coord[2] < zmin)
      zmin = node_list[elemlist->conn[i]].coord[2];
    if (node_list[elemlist->conn[i]].coord[2] > zmax)
      zmax = node_list[elemlist->conn[i]].coord[2];

  }

  new_id[0] = id_elem;
  RtreeInsert (Surf3DElemTree, (void *)new_id, xmin, xmax, ymin, ymax, zmin, zmax);

}


/* -------------------------------------------------------------------
** Surf3DAddElem - add a new element to the element list.
*/

static void Surf3DAddElem (int *num_elems, int **elements, Surf3DElemList **elemlist,
                           Surf3DBdryEdge *edge, int node_indx)
{
  int             *slot, i;
  Surf3DElemList  *curr_elem;
  double          r[3], s[3], len;


  if ( *num_elems >= num_elem_alloced )
  {
    if ( !num_elem_alloced )
    {
      num_elem_alloced = Surf3D_ELEM_QUANTUM ;
      *elements = (int *)Surf3DMalloc (num_elem_alloced * 3 * sizeof(int)) ;
      *elemlist = (Surf3DElemList *) Surf3DMalloc (num_elem_alloced * sizeof(Surf3DElemList)) ;
    }
    else
    {
      num_elem_alloced += Surf3D_ELEM_QUANTUM ;
      *elements = (int *)Surf3DRealloc (*elements, num_elem_alloced * 3 * sizeof(int)) ;
      *elemlist = (Surf3DElemList *)Surf3DRealloc (*elemlist, num_elem_alloced * sizeof(Surf3DElemList)) ;
    }
  }

  /* store away the vertex numbers */
  slot = &(*elements)[3*(*num_elems)] ;
  slot[0] = edge->verts[0] ;
  slot[1] = edge->verts[1] ;
  slot[2] = node_indx ;

#if PRINT_DEBUG
      fprintf (stdout, "Create new triangle id %d = -%d-%d-%d- use %d\n",
               (*num_elems), edge->verts[0], node_indx, edge->verts[1], edge->use);
#endif

  /* store vertex numbers */
  curr_elem = &(*elemlist)[*num_elems];
  curr_elem->conn[0] = edge->verts[0] ;
  curr_elem->conn[1] = edge->verts[1] ;
  curr_elem->conn[2] = node_indx ;

  /* compute normal */
  for (i = 0; i < 3; i++)
  {
    r[i] = node_list[node_indx].coord[i]      - node_list[edge->verts[0]].coord[i];
    s[i] = node_list[edge->verts[1]].coord[i] - node_list[edge->verts[0]].coord[i];
  }

  curr_elem->normal[0] = s[1] * r[2] - r[1] * s[2] ;
  curr_elem->normal[1] = s[2] * r[0] - r[2] * s[0] ;
  curr_elem->normal[2] = s[0] * r[1] - r[0] * s[1] ;
  len = sqrt (curr_elem->normal[0] * curr_elem->normal[0] +
              curr_elem->normal[1] * curr_elem->normal[1] +
              curr_elem->normal[2] * curr_elem->normal[2] );

  //rgd
  if (len!=0)
  {
  for (i = 0; i < 3; i++)
    curr_elem->normal[i] /= len;
  }

  /* update adjacente elements to vertexs */
  Surf3DUpdateAdjNodeElem (edge->verts[0], *elemlist, *num_elems);
  Surf3DUpdateAdjNodeElem (edge->verts[1], *elemlist, *num_elems);
  Surf3DUpdateAdjNodeElem (node_indx,      *elemlist, *num_elems);

  Surf3DAddElemRTree (curr_elem, *num_elems);

  (*num_elems)++ ;

}

/* -------------------------------------------------------------------
** Surf3DUpdateAdjNodeElem - update adjacente element to node.
*/
void Surf3DUpdateAdjNodeElem (int node_indx, Surf3DElemList *elem_list, int new_elem)
{
  Surf3DAdjElem   *adj_elem;
  double          n[3], i;

  /* update adjacente elements to vertexs */
  adj_elem = (Surf3DAdjElem *) calloc (1, sizeof (Surf3DAdjElem));
  adj_elem->next = node_list[node_indx].elems;
  adj_elem->id   = new_elem;
  node_list[node_indx].elems = adj_elem;

  /* update normal vertexs */
  n[0] = n[1] = n[2] = 0.0;
  i = 0;
  while (adj_elem)
  {
    n[0] += elem_list[adj_elem->id].normal[0];
    n[1] += elem_list[adj_elem->id].normal[1];
    n[2] += elem_list[adj_elem->id].normal[2];
    i++;
    adj_elem = adj_elem->next;
  }
  node_list[node_indx].normal[0] = n[0]/i;
  node_list[node_indx].normal[1] = n[1]/i;
  node_list[node_indx].normal[2] = n[2]/i;
}

/* -------------------------------------------------------------------
** Surf3DCheckArea - Check the area of added element.
*/
static int Surf3DCheckArea (Surf3DBdryEdge *edge, int node_indx)
{
  double node[3], old[3], area = 0.0, min_area = 0.0 ;
  int    i;

  /* return if it's not the phase of ideal elements */
  if (edge->use != 0) return 1 ;

  /* compute the area of polygon */
  area = Surf3DAreaEdge_Node (edge, node_indx, node_list) ;

  /* compute the minimal node */
  Surf3DFindNode (edge, 10, node) ;

  /* save node_indx choosen */
  for (i = 0; i < 3; i++)
    old[i] = node_list[node_indx].coord[i] ;

  /* set new node_indx */
  for (i = 0; i < 3; i++)
    node_list[node_indx].coord[i] = node[i] ;

  /* compute the minimal area */
  min_area = Surf3DAreaEdge_Node (edge, node_indx, node_list) ;

  /* restore node_indx */
  for (i = 0; i < 3; i++)
    node_list[node_indx].coord[i] = old[i] ;

  /* Compare the two areas */
  if( fabs(area) < fabs(min_area) )
    return 0 ;
  else
    return 1 ;

}






/* -------------------------------------------------------------------*/
/* -------------------------------------------------------------------*/
/* -------------------------------------------------------------------*/
/* -------------------------------------------------------------------*/
/* -------------------------------------------------------------------*/
/* -------------------------------------------------------------------*/
/**   DATA ESTRUCTURES                                                */
/* -------------------------------------------------------------------*/
/* -------------------------------------------------------------------*/
/* -------------------------------------------------------------------*/
/* -------------------------------------------------------------------*/
/* -------------------------------------------------------------------*/
/* -------------------------------------------------------------------*/



/* -------------------------------------------------------------------
** Surf3DHeap - these routines manage a priorty queue using a heap
**             data structure.
*/

typedef struct _Surf3DHeapEntry
{
    double    value ;
    int       indx ;
} Surf3DHeapEntry, *Surf3DHeap ;

static Surf3DHeap the_heap = 0 ;
static int       alloced_size ;
static int       heap_last ;

/* -------------------------------------------------------------------*/
static void Surf3DHeapInit( int size )
{
    /* make sure we have enough room for the heap */

    if ( !the_heap )
    {
      alloced_size = size * 10;
      the_heap = (Surf3DHeap)malloc (alloced_size * sizeof(Surf3DHeapEntry) ) ;
      alloced_size = size ;
    }
    else if ( alloced_size*10 < size )
    {
      Surf3DFree( the_heap ) ;
      alloced_size = size * 10;
      the_heap = (Surf3DHeap)malloc (alloced_size * sizeof(Surf3DHeapEntry) ) ;
    }
    heap_last = -1 ;
}

/* -------------------------------------------------------------------*/
static void Surf3DHeapDelete(void)
{
    Surf3DFree( the_heap ) ;
    the_heap = 0 ;
}

/* -------------------------------------------------------------------*/
static void Surf3DHeapInsert
(
double    value,
int      indx
)
{
    int  cur, parent ;

    /* put the new values at the end of the heap */

    heap_last++ ;

    if (heap_last < 0)
      heap_last = 0;
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

/* -------------------------------------------------------------------*/
static int Surf3DHeapExtract( double  *value )
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

/* -------------------------------------------------------------------*/
static double Surf3DHeapAngle
(
int                vi,
int                vj,
int                indx
)
{
  double v1[3], v2[3];
  double angle;
  int    i;

  /* get the index vertex */
  for (i = 0; i < 3; i++)
  {
    v1[i] = node_list[vi].coord[i] - node_list[indx].coord[i];
    v2[i] = node_list[vj].coord[i] - node_list[indx].coord[i];
  }

  angle = Surf3DVectAngle (v1, v2);

  return(angle);
}




/* -------------------------------------------------------------------
** Surf3DEdgeStack - these routines manage a stack that is used to store
**                  the active boundary edge data.
*/

static Surf3DBdryEdge *Surf3DEdgeAlloc( )
{
    Surf3DBdryEdge  *new_block, *alloced ;
    int            i ;

    /* if the free pointer is null we need to allocate a new block
       of boundary nodes */

    if ( !bdry_free )
    {
     new_block = (Surf3DBdryEdge *)Surf3DMalloc(
      Surf3D_BDRY_FACE_QUANTUM * sizeof(Surf3DBdryEdgeRec) ) ;
     new_block[0].next = bdry_block_ptr ;
     bdry_block_ptr = new_block ;
     for ( i=1 ; i<(Surf3D_BDRY_FACE_QUANTUM-1) ; i++ )
     {
      new_block[i].next = &(new_block[i+1]) ;
     }
     new_block[Surf3D_BDRY_FACE_QUANTUM-1].next = 0 ;
     bdry_free = &(new_block[1]) ;
    }

    /* return the top thing on the free list */

    alloced = bdry_free ;
    bdry_free = bdry_free->next ;

    return( alloced ) ;
}

/* -------------------------------------------------------------------*/
static void Surf3DBdryFree( Surf3DBdryEdge  *edge )
{
    /* put this edge back on the free list */

    edge->next = bdry_free ;
    bdry_free = edge ;
}

/* -------------------------------------------------------------------*/
static void Surf3DEdgeFreeAll()
{
    Surf3DBdryEdge  *cur, *next ;

    /* free all blocks allocated to store edge information */

    if ( bdry_block_ptr ) cur = bdry_block_ptr ;
    else return ;

    while ( cur->next ) {
    next = cur->next ;
    Surf3DFree( cur ) ;
    cur = next ;
    }
    Surf3DFree( cur ) ;

    bdry_stack = 0 ;
    bdry_tail  = 0 ;
    bdry_free  = 0 ;
    bdry_try   = 0 ;
    bdry_cursor    = 0 ;
    bdry_block_ptr = 0 ;
}

/* -------------------------------------------------------------------*/
static void Surf3DBdryPush( Surf3DBdryEdge  *edge )
{
    /* push this edge on the end of the stack that is implemented
       as a doubly linked list */

    edge->prev = bdry_tail ;
    edge->next = 0 ;
    if ( bdry_tail ) bdry_tail->next = edge ;
    bdry_tail = edge ;
    if ( !bdry_stack ) bdry_stack = edge ;
}

/* -------------------------------------------------------------------*/
static void Surf3DBdryPushCorre( Surf3DBdryEdge  *edge )
{
  Surf3DBdryEdge  *bdry_cur = 0 ;

  /* find the bdry_try, that is, the first edge that went one time
      to the end of the stack */

  bdry_try = 0 ;
  bdry_cur = bdry_stack;
  do
  {
    // rgd
    if (bdry_cur==NULL)
     {bdry_cur = bdry_tail;continue;}
    if (bdry_cur->use > edge->use)
    {
      bdry_try = bdry_cur ;
      break ;
    }
    bdry_cur = bdry_cur->next ;
  } while (bdry_cur != bdry_tail && bdry_cur != NULL);

  /* modify the list to push the edge in its position considering the
      list only from bdry_stack to bdry_try, because edges from bdry_try
      until bdry_tail are edges already tested and pushed to the end of
      the list */

  if (bdry_try == 0)
  {
    edge->prev = bdry_tail ;
    edge->next = 0 ;
    if ( bdry_tail )
      bdry_tail->next = edge ;
    bdry_tail = edge ;
    if ( !bdry_stack )
      bdry_stack = edge ;
  }
  else
  {
    edge->prev = bdry_try->prev ;
    edge->next = bdry_try ;
    if ( bdry_try->prev )
      bdry_try->prev->next = edge ;
    bdry_try->prev = edge ;
    if ( bdry_try == bdry_stack )
      bdry_stack = edge ;
  }

}

/* -------------------------------------------------------------------*/
static void Surf3DBdryPushSmall( Surf3DBdryEdge  *edge )
{
  Surf3DBdryEdge  *bdry_cur = 0 ;
  Surf3DBdryEdge  *bdry_real_tail = 0 ;
  int            found_pos = 0 ;

  /* find the key value for the edge */
  edge->key = edge->length ;

  /* find the bdry_try, that is, the first edge that went one time
      to the end of the stack */

  bdry_try = 0 ;
  bdry_cur = bdry_stack;
  do
  {
    if( bdry_cur )
    {
      if (bdry_cur->use > 0)
      {
        bdry_try = bdry_cur ;
        break ;
      }
      bdry_cur = bdry_cur->next ;
      if( !bdry_cur )
        break ;
    }
  } while ( bdry_cur != bdry_tail ) ;

  /* modify the list to push the edge in its position considering the
      list only from bdry_stack to bdry_try, because edges from bdry_try
      until bdry_tail are edges already tested and pushed to the end of
      the list */

  if ( (bdry_try != 0) && (bdry_try != 0) )
  {
    if (bdry_real_tail == 0 && bdry_try->prev != 0)
    {
      bdry_real_tail = bdry_tail ;
      bdry_tail = bdry_try->prev ;
      bdry_tail->next = 0 ;
    }
  }

  /* push this edge in its right position in the list acoording
      wiht its length. Here is only considered the list from the
      bdry_stack until bdry_try */

  if ( !found_pos )
  {
    if ( (bdry_stack == 0) || (edge->key <= bdry_stack->key) )
    {
      edge->prev = 0 ;
      edge->next = bdry_stack ;
      if ( bdry_stack )
        bdry_stack->prev = edge ;
      bdry_stack = edge ;
      if ( !bdry_tail )
        bdry_tail = edge ;
      found_pos = 1 ;
    }
  }

  if ( !found_pos )
  {
    if ( (bdry_tail != 0) && (edge->key >= bdry_tail->key) )
    {
      edge->prev = bdry_tail ;
      edge->next = 0 ;
      if ( bdry_tail ) bdry_tail->next = edge ;
      bdry_tail = edge ;
      found_pos = 1 ;
    }
  }

  if ( !found_pos )
  {
    bdry_cur = bdry_stack;
    do
    {
      if ( edge->key == bdry_cur->key )
      {
        edge->prev = bdry_cur->prev ;
        edge->next = bdry_cur ;
        if ( bdry_cur->prev ) bdry_cur->prev->next = edge ;
        bdry_cur->prev = edge ;
        found_pos = 1 ;
        break ;
      }
      else if ( (edge->key > bdry_cur->key) && (edge->key < bdry_cur->next->key) )
      {
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

  if ( (bdry_try != 0) && (bdry_try->prev != 0) )
  {
    if (bdry_real_tail != 0)
    {
      bdry_tail->next = bdry_try ;
      bdry_try->prev = bdry_tail ;
      bdry_tail = bdry_real_tail ;
    }
  }
}

/* -------------------------------------------------------------------*/
/* markos - adding edge to the front based on use, layer and length, in that order */
static void Surf3DBdryPushLayer( Surf3DBdryEdge  *edge )
{
    Surf3DBdryEdge  *bdry_cur = NULL, *bdry_cur_back = NULL;
    int from_back = 0;

    /* find the key value for the edge */
    edge->key = edge->length;

    if (!bdry_stack)
    {
        edge->prev = NULL;
        edge->next = NULL;
        bdry_stack = edge;
        bdry_tail = edge;

        return;
    }

    for (bdry_cur = bdry_stack, bdry_cur_back = bdry_tail; bdry_cur; bdry_cur = bdry_cur->next, bdry_cur_back = bdry_cur_back->prev)
    {
        if ((bdry_cur_back->use < edge->use) ||
            ((bdry_cur_back->use == edge->use) && (bdry_cur_back->layer < edge->layer)) ||
            ((bdry_cur_back->use == edge->use) && (bdry_cur_back->layer == edge->layer) && (bdry_cur_back->key <= edge->key)))
        {
            from_back = 1;

            break;
        }

        if ((bdry_cur->use > edge->use) ||
            ((bdry_cur->use == edge->use) && (bdry_cur->layer > edge->layer)) ||
            ((bdry_cur->use == edge->use) && (bdry_cur->layer == edge->layer) && (bdry_cur->key > edge->key)))
        {
            break;
        }
    }

    if (from_back)
    {
        bdry_cur = (!bdry_cur_back) ? bdry_stack : bdry_cur_back->next;
    }

    if (!bdry_cur)
    {
        edge->prev = bdry_tail;
        edge->next = NULL;
        bdry_tail->next = edge;
        bdry_tail = edge;

        return;
    }

    edge->prev = bdry_cur->prev;
    edge->next = bdry_cur;
    if (edge->prev) edge->prev->next = edge;
    bdry_cur->prev = edge;
    if (bdry_cur == bdry_stack) bdry_stack = edge;
}

/* -------------------------------------------------------------------*/
static Surf3DBdryEdge *Surf3DBdryPop()
{
    Surf3DBdryEdge  *popped ;

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

/* -------------------------------------------------------------------*/
static void Surf3DBdryDelete( Surf3DBdryEdge  *edge )
{
    /* markos - remove edge from rtree */
    RtreeDelete (Surf3DEdgeTree, edge, edge->min[0], edge->max[0],
                   edge->min[1], edge->max[1], edge->min[2], edge->max[2]);
    /* delete this edge from the middle of the doubly linked list */
    if ( bdry_stack == edge ) bdry_stack = edge->next ;
    if ( bdry_tail  == edge ) bdry_tail  = edge->prev ;
    if ( edge->next ) edge->next->prev = edge->prev ;
    if ( edge->prev ) edge->prev->next = edge->next ;
    Surf3DBdryFree( edge ) ;
}

/* -------------------------------------------------------------------*/
static void Surf3DBdryReset(void)
{
    /* reset the cursor for scanning through the list */

    bdry_cursor = bdry_stack ;
}

/* -------------------------------------------------------------------*/
static Surf3DBdryEdge *Surf3DBdryNext(void)
{
    Surf3DBdryEdge *current ;

    /* return the edge the cursor is pointing to, and increment the
       cursor */

    current = bdry_cursor ;
    if ( bdry_cursor ) bdry_cursor = bdry_cursor->next ;
    return( current ) ;
}

#if 0
/* -------------------------------------------------------------------
** Surf3DAdjElemAlloc - these routines manage the allocation and freeing
**                     of adjacent elem list entries
*/

static Surf3DAdjElem *Surf3DAdjElemAlloc( )
{
    Surf3DAdjElem  *new_block, *alloced ;
    int           i ;

    /* if the free pointer is null we need to allocate a new block
       of elements */

    if ( !adj_elem_free )
    {
     new_block = (Surf3DAdjElem *)Surf3DMalloc(
                  Surf3D_ADJ_ELEM_QUANTUM * sizeof(Surf3DAdjElemRec) ) ;
     new_block[0].next = adj_elem_block_ptr ;
     adj_elem_block_ptr = new_block ;
     for ( i=1 ; i<(Surf3D_ADJ_ELEM_QUANTUM-1) ; i++ )
     {
       new_block[i].next = &(new_block[i+1]) ;
     }
      new_block[Surf3D_ADJ_ELEM_QUANTUM-1].next = 0 ;
      adj_elem_free = &(new_block[1]) ;
    }

    /* return the top thing on the free list */

    alloced = adj_elem_free ;
    adj_elem_free = adj_elem_free->next ;

    return( alloced ) ;
}
#endif


/* -------------------------------------------------------------------*/
static void Surf3DAdjElemFreeAll()
{
    Surf3DAdjElem  *cur, *next ;

    /* free all blocks allocated to store elem information */

    if ( adj_elem_block_ptr ) cur = adj_elem_block_ptr ;
    else return ;

    while ( cur->next )
    {
     next = cur->next ;
     Surf3DFree( cur ) ;
     cur = next ;
    }
    Surf3DFree( cur ) ;

    adj_elem_free = 0 ;
    adj_elem_block_ptr = 0 ;
}

/* -------------------------------------------------------------------
** Surf3DAdjNodeAlloc - these routines manage the allocation and freeing
**                     of node list entries
*/

static void Surf3DAdjNodeFreeAll()
{
    Surf3DBdryNodeList *cur ;

    /* free all blocks allocated to store node information */

    cur = node_list ;

    Surf3DFree( cur ) ;

    node_list = NULL ;
}

/* -------------------------------------------------------------------
** Surf3DAdjEdgeAlloc - these routines manage the allocation and freeing
**                     of adjacent edge list entries
*/

static Surf3DAdjEdge *Surf3DAdjEdgeAlloc( )
{
    Surf3DAdjEdge  *new_block, *alloced ;
    int           i ;

    /* if the free pointer is null we need to allocate a new block
       of boundary nodes */

    if ( adj_free == NULL)
    {
     new_block = (Surf3DAdjEdge *)calloc(Surf3D_ADJ_FACE_QUANTUM, sizeof(Surf3DAdjEdgeRec)) ;
     new_block[0].next = adj_block_ptr ;
     adj_block_ptr = new_block ;
     for ( i=1 ; i<(Surf3D_ADJ_FACE_QUANTUM-1) ; i++ )
     {
      new_block[i].next = &(new_block[i+1]) ;
     }
     new_block[Surf3D_ADJ_FACE_QUANTUM-1].next = 0 ;
     adj_free = &(new_block[1]) ;
    }

    /* return the top thing on the free list */

    alloced = adj_free ;
    adj_free = adj_free->next ;

    return( alloced ) ;
}

/* -------------------------------------------------------------------*/
static void Surf3DAdjFree( Surf3DAdjEdge  *edge )
{
    /* put this edge back on the free list */

    edge->next = adj_free ;
    adj_free = edge ;
}

/* -------------------------------------------------------------------*/
static void Surf3DAdjFreeAll()
{
    Surf3DAdjEdge  *cur, *next ;

    /* free all blocks allocated to store edge information */

    if ( adj_block_ptr ) cur = adj_block_ptr ;
    else return ;

    while ( cur->next )
    {
     next = cur->next ;
     Surf3DFree( cur ) ;
     cur = next ;
    }
    Surf3DFree( cur ) ;

    adj_free = 0 ;
    adj_block_ptr = 0 ;
}

/* -------------------------------------------------------------------
** Surf3DAdjIniEdgeAlloc - these routines manage the allocation and
**                        freeing of initial adjacent edge list entries
*/

static Surf3DAdjIniEdge *Surf3DAdjIniEdgeAlloc( )
{
    Surf3DAdjIniEdge  *new_block, *alloced ;
    int              i ;

    /* if the free pointer is null we need to allocate a new block
       of boundary nodes */

    if ( !adj_ini_free )
    {
     new_block = (Surf3DAdjIniEdge *)Surf3DMalloc(
       Surf3D_ADJ_INI_FACE_QUANTUM * sizeof(Surf3DAdjIniEdgeRec) ) ;
     new_block[0].next = adj_ini_block_ptr ;
     adj_ini_block_ptr = new_block ;
     for ( i=1 ; i<(Surf3D_ADJ_INI_FACE_QUANTUM-1) ; i++ )
     {
       new_block[i].next = &(new_block[i+1]) ;
     }
     new_block[Surf3D_ADJ_INI_FACE_QUANTUM-1].next = 0 ;
     adj_ini_free = &(new_block[1]) ;
    }

    /* return the top thing on the free list */

    alloced = adj_ini_free ;
    adj_ini_free = adj_ini_free->next ;

    return( alloced ) ;
}

/* -------------------------------------------------------------------*/
static void Surf3DAdjIniFreeAll()
{
    Surf3DAdjIniEdge  *cur, *next ;

    /* free all blocks allocated to store edge information */

    if ( adj_ini_block_ptr ) cur = adj_ini_block_ptr ;
    else return ;

    while ( cur->next )
    {
     next = cur->next ;
     Surf3DFree( cur ) ;
     cur = next ;
    }
    Surf3DFree( cur ) ;

    adj_ini_free = 0 ;
    adj_ini_block_ptr = 0 ;
}

/* -------------------------------------------------------------------
** Surf3DTestStack - these routines manage a stack that is used to store
**                  the boundary edge data already done.
*/

static void Surf3DTestFreeAll()
{
    Surf3DBdryEdge  *cur, *next ;

    /* free all blocks allocated to store edge information */

    if ( test_block_ptr ) cur = test_block_ptr ;
    else return ;

    while ( cur->next )
    {
     next = cur->next ;
     Surf3DFree( cur ) ;
     cur = next ;
    }
    Surf3DFree( cur ) ;

    test_tail  = 0 ;
    test_block_ptr = 0 ;
}


#if 0
/* -------------------------------------------------------------------
** Surf3DCrossProd - Compute the Cross Prod.
*/
static double Surf3DCrossProd
(
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

   kx = node_list[K].coord[0] - node_list[I].coord[0];
   ky = node_list[K].coord[1] - node_list[I].coord[1];
   jx = node_list[J].coord[0] - node_list[I].coord[0];
   jy = node_list[J].coord[1] - node_list[I].coord[1];

   return (kx*jy - ky*jx);
}


/* -------------------------------------------------------------------
** Surf3DCrossEdge
*/
static int Surf3DCrossEdge
(
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

  if ((node_list[I1].coord[0] > node_list[J1].coord[0]) &&
      (node_list[I1].coord[0] > node_list[J2].coord[0]) &&
      (node_list[I2].coord[0] > node_list[J1].coord[0]) &&
      (node_list[I2].coord[0] > node_list[J2].coord[0]))
  {
    cross = FALSE;
    return(cross);
  }

  if ((node_list[I1].coord[1] > node_list[J1].coord[1]) &&
      (node_list[I1].coord[1] > node_list[J2].coord[1]) &&
      (node_list[I2].coord[1] > node_list[J1].coord[1]) &&
      (node_list[I2].coord[1] > node_list[J2].coord[1]))
  {
    cross = FALSE;
    return(cross);
  }

/*  Now cross line I with J1 and line I with J2, if have the same sign they
   cannot cross */

  if (Surf3DCrossProd(I1,I2,J1) *
      Surf3DCrossProd(I1,I2,J2) >= 0.0)
  {
    cross = FALSE;
    return(cross);
  }

/*  Now cross line J with I1 and line J with I2, if have the same sign they
     cannot cross */

  if (Surf3DCrossProd(J1,J2,I1) *
      Surf3DCrossProd(J1,J2,I2) >= 0.0)
  {
    cross = FALSE;
    return(cross);
  }

  return(cross);
}
#endif


/* -------------------------------------------------------------------
** Surf3DDraw - Draw current entities from mesh.
*/
void Surf3DDraw (void *dummy)
{
#if SURF3DDRAW

  Surf3DDrawBegin ( );

  /* Curr edge and optimal node */
  if (draw_edge != NULL)
  {
/*    Surf3DDrawCurrEdge (node_list[draw_edge->verts[0]].coord,
                        node_list[draw_edge->verts[1]].coord, opt_node);
*/
    Surf3DDrawCurrEdge (node_list[draw_edge->verts[0]].coord,
                        node_list[draw_edge->verts[1]].coord, &try_pts[0]);
  }

  /* node */
  Surf3DDrawNodes (*ptr_nnodes, node_list);
  // Surf3DDrawNodeNumber (*ptr_nnodes, node_list);
  Surf3DDrawNodeNormal (*ptr_nnodes, node_list);

  /* Boundary  */
  Surf3DDrawBoundary (node_list, bdry_stack);
  Surf3DDrawBoundaryNormal (node_list, bdry_stack);

  /* elements */
  Surf3DDrawElements (node_list, *ptr_nelem, *ptr_elems);
  Surf3DDrawElementNormal (node_list, *ptr_nelem, *ptr_elems);

  Surf3DDrawFlush ( );

  Delay ( );
#endif
}

/* markos - functions that test the triangle x triangle intersection, based on Msh3D */

/*
** ----------------------------------------------------------------------
** Geo routines - These routines check if two faces, defined by their
**                nodes (a,b,c - first face and u,v,w - second face)
**                intersect each other. This package is formed by the
**                following routines:
**                      1) GeoTriIntersect ;
**                      2) GeoTriSegIntersect ;
**                      3) GeoTriSegIntersectCheck ;
**                      4) GeoSegSegIntersect ;
**                      5) GeoNorm ;
**                      6) GeoCrossProd ;
**                      7) GeoCross ;
**                      8) GeoCrossNorm ;
**                      9) GeoDot.
*/

static int  GeoTriIntersect(
                     double a[3], double b[3], double c[3], int p[3],
                     double u[3], double v[3], double w[3], int q[3],
                     double tol_inters )
{
    int                type ;

    /* test intersection against two given triangles and return if
       intersection that it's not coincident with any vertices is
       found */

    if (q[0]==0) {
     type = GeoTriSegIntersect(a,b,c,u,v,tol_inters) ;
     if (type == 0)  return type ;
    }
    if (q[1]==0) {
     type = GeoTriSegIntersect(a,b,c,v,w,tol_inters);
     if (type == 0)  return type ;
    }
    if (q[2]==0) {
     type = GeoTriSegIntersect(a,b,c,w,u,tol_inters);
     if (type == 0)  return type ;
    }
    if (p[0]==0) {
     type = GeoTriSegIntersect(u,v,w,a,b,tol_inters);
     if (type == 0)  return type ;
    }
    if (p[1]==0) {
     type = GeoTriSegIntersect(u,v,w,b,c,tol_inters);
     if (type == 0)  return type ;
    }
    if (p[2]==0) {
     type = GeoTriSegIntersect(u,v,w,c,a,tol_inters);
     if (type == 0)  return type ;
    }

   /* if we get here no intersection was found */

   return 1 ;
}

static int GeoTriSegIntersect(
    double             a[3],
    double             b[3],
    double             c[3],
    double             u[3],
    double             v[3],
    double             tol_inters )
{
    int    inter ;
    double d, nc;
    double v1, v2;
    double p[3], n[3];
    double nt[3];          /* auxiliar vector to store normal */

    /* compute triangle normal vector */

    GeoCross(a,b,c,n);
    GeoNorm(n);

    /* compute triangle plane independent term */

    d = -a[0]*n[0] - a[1]*n[1] - a[2]*n[2];

    /* compute distances of segment endpoints to the triangle plane */

    v1 = u[0]*n[0] + u[1]*n[1] + u[2]*n[2] + d;
    v2 = v[0]*n[0] + v[1]*n[1] + v[2]*n[2] + d;

    /* check if triangle and segment are coplanar */

    if (fabs(v1) < tol_inters && fabs(v2) < tol_inters)
    {
     /* check if endpoints are inside triangle */

     if (GeoDot(n,GeoCrossNorm(u,a,b,nt)) > 0.0 &&
         GeoDot(n,GeoCrossNorm(u,b,c,nt)) > 0.0 &&
         GeoDot(n,GeoCrossNorm(u,c,a,nt)) > 0.0
        ) {
      inter = 1 ;   /* intersection point is u */
      if (GeoTriSegIntersectCheck (u,a,tol_inters)) inter = 0 ;
      if (GeoTriSegIntersectCheck (u,b,tol_inters)) inter = 0 ;
      if (GeoTriSegIntersectCheck (u,c,tol_inters)) inter = 0 ;
      if (inter) return 0 ;
     }
     if (GeoDot(n,GeoCrossNorm(v,a,b,nt)) > 0.0 &&
         GeoDot(n,GeoCrossNorm(v,b,c,nt)) > 0.0 &&
         GeoDot(n,GeoCrossNorm(v,c,a,nt)) > 0.0
        ) {
      inter = 1 ;   /* intersection point is v */
      if (GeoTriSegIntersectCheck (v,a,tol_inters)) inter = 0 ;
      if (GeoTriSegIntersectCheck (v,b,tol_inters)) inter = 0 ;
      if (GeoTriSegIntersectCheck (v,c,tol_inters)) inter = 0 ;
      if (inter) return 0 ;
     }

     if (GeoSegSegIntersect(a,b,u,v,p,tol_inters))
     {
      inter = 1 ;   /* intersection point is p */
      if (GeoTriSegIntersectCheck (p,a,tol_inters)) inter = 0 ;
      if (GeoTriSegIntersectCheck (p,b,tol_inters)) inter = 0 ;
      if (GeoTriSegIntersectCheck (p,c,tol_inters)) inter = 0 ;
      if (inter) return 0 ;
     }
     if (GeoSegSegIntersect(b,c,u,v,p,tol_inters))
     {
      inter = 1 ;   /* intersection point is p */
      if (GeoTriSegIntersectCheck (p,a,tol_inters)) inter = 0 ;
      if (GeoTriSegIntersectCheck (p,b,tol_inters)) inter = 0 ;
      if (GeoTriSegIntersectCheck (p,c,tol_inters)) inter = 0 ;
      if (inter) return 0 ;
     }
     if (GeoSegSegIntersect(c,a,u,v,p,tol_inters))
     {
      inter = 1 ;   /* intersection point is p */
      if (GeoTriSegIntersectCheck (p,a,tol_inters)) inter = 0 ;
      if (GeoTriSegIntersectCheck (p,b,tol_inters)) inter = 0 ;
      if (GeoTriSegIntersectCheck (p,c,tol_inters)) inter = 0 ;
      if (inter) return 0 ;
     }
    }

    /* check if segment cross or touch the triangle plane */

    else if (v1*v2 < 0.0)
    {
     /* compute intersection between segment and triangle */

     double t = v1/(v1-v2);
     p[0] = u[0] + t*(v[0]-u[0]);
     p[1] = u[1] + t*(v[1]-u[1]);
     p[2] = u[2] + t*(v[2]-u[2]);

     /* check if intersection point is inside triangle */

     if (GeoDot(n,GeoCrossNorm(p,a,b,nt)) >= -tol_inters &&
         GeoDot(n,GeoCrossNorm(p,b,c,nt)) >= -tol_inters &&
         GeoDot(n,GeoCrossNorm(p,c,a,nt)) >= -tol_inters
        )
     {
      {
       inter = 1 ;   /* intersection point is p */
       if (GeoTriSegIntersectCheck (p,a,tol_inters)) inter = 0 ;
       if (GeoTriSegIntersectCheck (p,b,tol_inters)) inter = 0 ;
       if (GeoTriSegIntersectCheck (p,c,tol_inters)) inter = 0 ;
       if (inter) return 0 ;
      }
     }
    }

    /* if we get here no intersection was found */

    return 1 ;
}

static int GeoTriSegIntersectCheck( double p[3], double v[3], double tol_inters )
{
    /* verify if the intersection point and the vertex are the same */

    if( ABOUT_ZERO((p[0]-v[0]),tol_inters) &&
        ABOUT_ZERO((p[1]-v[1]),tol_inters) &&
        ABOUT_ZERO((p[2]-v[2]),tol_inters)  )
     return 1 ;
    else
     return 0 ;
}

static int GeoSegSegIntersect
(double p1[3], double q1[3], double p2[3], double q2[3], double p[3], double tol_inters)
{
 double det, t, l1, l2, c2, c[3];
 double v1[3], v2[3] ;

 v1[0] = q1[0]-p1[0];
 v1[1] = q1[1]-p1[1];
 v1[2] = q1[2]-p1[2];
 v2[0] = q2[0]-p2[0];
 v2[1] = q2[1]-p2[1];
 v2[2] = q2[2]-p2[2];
 l1 = GeoNorm(v1);
 l2 = GeoNorm(v2);
 c[0] = v1[1]*v2[2] - v1[2]*v2[1];
 c[1] = v1[2]*v2[0] - v1[0]*v2[2];
 c[2] = v1[0]*v2[1] - v1[1]*v2[0];
 c2 = (c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);

 /* test if segments are colinear */
 if (c2 < tol_inters)
  return 0;

 /* segments are not colinear */

 /* compute parameter related to second line */
 det = ((p2[0]-p1[0])*v1[1]*c[2] +
        (p2[1]-p1[1])*v1[2]*c[0] +
        (p2[2]-p1[2])*v1[0]*c[1]
       )
       -
       (
        (p2[2]-p1[2])*v1[1]*c[0] +
        (p2[0]-p1[0])*v1[2]*c[1] +
        (p2[1]-p1[1])*v1[0]*c[2]
       );
 t = det/c2;
 if (t < 0.0 || t > l2)
  return 0;

 /* compute parameter related to first line */
 det = ((p2[0]-p1[0])*v2[1]*c[2] +
        (p2[1]-p1[1])*v2[2]*c[0] +
        (p2[2]-p1[2])*v2[0]*c[1]
       )
       -
       (
        (p2[2]-p1[2])*v2[1]*c[0] +
        (p2[0]-p1[0])*v2[2]*c[1] +
        (p2[1]-p1[1])*v2[0]*c[2]
       );
 t = det/c2;
 if (t < 0.0 || t > l1)
  return 0;
 else
 {
  p[0] = p1[0] + t*v1[0];
  p[1] = p1[1] + t*v1[1];
  p[2] = p1[2] + t*v1[2];
  return 1;
 }
}

static double GeoNorm (double v[3])
{
 double lenght = (double)(sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
 /*if (jmesh_version == SMALL_VERSION) {*/
  if( !ABOUT_ZERO(lenght,1.0e-7) ) {
   v[0] /= lenght; v[1] /= lenght; v[2] /= lenght;
  }
 /*} else {
  if( !ABOUT_ZERO(lenght,1.0e-12) ) {
   v[0] /= lenght; v[1] /= lenght; v[2] /= lenght;
  }
 }*/
 return lenght;
}

static double *GeoCrossProd(double a[3], double b[3], double n[3])
{
 n[0] = (a[1]*b[2]) - (a[2]*b[1]) ;
 n[1] = (a[2]*b[0]) - (a[0]*b[2]) ;
 n[2] = (a[0]*b[1]) - (a[1]*b[0]) ;
 return n;
}

static double *GeoCross (double a[3], double b[3], double c[3], double n[3])
{
 n[0] = (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]);
 n[1] = (b[2]-a[2])*(c[0]-a[0]) - (b[0]-a[0])*(c[2]-a[2]);
 n[2] = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
 return n;
}

static double *GeoCrossNorm (double a[3], double b[3], double c[3], double n[3])
{
 n[0] = (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]);
 n[1] = (b[2]-a[2])*(c[0]-a[0]) - (b[0]-a[0])*(c[2]-a[2]);
 n[2] = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
 GeoNorm (n);
 return n;
}

static double GeoDot (double u[3], double v[3])
{
 return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

