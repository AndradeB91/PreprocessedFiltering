/*
** ---------------------------------------------------------------
**
** msh_quad2d.c:  Definitions and prototypes of package for
**                implementing usual 2D quadtree operations.
**
** ---------------------------------------------------------------
**
** Created:      xx-xxx-xx
** Modify:       01-Mar-98      Antonio C.O. Miranda (2D)
** ---------------------------------------------------------------
**
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <time.h>

#define  QTREE_INTP   1


#ifdef QUAD_DRAW
#ifdef USE_IUP
#include "iup.h"
#endif
#endif

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1.0e+06
#endif

#include "msh_quadsurf.h"
#include "msh_defsurf.h"

#define OLD_IMP 1

/* --------------------------------------------------------------
** Private definitons and data types:
*/

#define REF_FACTOR     0.9       /* changed by ref_factor */
#define NEW_NODES      10
#define INFPOS         1e+06
#ifndef PI
#define PI             acos(-1.0)
#endif
#define MAX(a,b)       (((a)>(b))?(a):(b))
#define MIN(a,b)       (((a)<(b))?(a):(b))
#define MSH2D_QUAD_TREE_QUANTUM 500

typedef struct _metric {
 double                     dev[6];       /* Derived */
 double                     normal[3];
 double                     center[2];
 int                        choose;
} Metric;

typedef struct _Msh2DQuadTreeRec {
  struct _Msh2DQuadTreeRec  *child[4] ;
  struct _Msh2DQuadTreeRec  *parent ;
  int                        flags ;
  int                        marks ;
  Metric                     *m[4];
  double                     size;
}
Msh2DQuadTreeRec, *Msh2DQuadTree ;

typedef struct _Box2d {
  double min[2], max[2], center[2] ;
}
Box2d ;


typedef struct ref_curvature
{
#if OLD_IMP
   double size;
   double mid[2];
#else
  struct _Msh2DQuadTreeRec  *tree;
  double                    newsize;
  double                    oldsize;
#endif
  struct ref_curvature      *next;
} Ref_curvature;
static Ref_curvature *init_curv = NULL;
static Ref_curvature *tail_curv = NULL;

/* --------------------------------------------------------------
** Global variables:
*/



/* --------------------------------------------------------------
** Private variables:
*/

static Msh2DQuadTree root = 0 ;
static Box2d  box ;
static double gmin[2], gmax[2] ;
static double bmin[2], bmax[2] ;
static double Distor_x, Distor_y;

static Msh2DQuadTree quad_free = 0 ;
static Msh2DQuadTree quad_block_ptr = 0 ;

double mshsurfCosAngle = (1.2457309396155173259666803366403/2.0);
double mshsurfAngle    = 3.1416/12.0;

/* --------------------------------------------------------------
** Private functions prototypes:
*/

static Msh2DQuadTree Msh2DBdryQuad( Msh2DQuadTree, Msh2DQuadTree, double,
                     double *, double *, double *, double, int, int *) ;
static void          Msh2DTransQuad( Msh2DQuadTree, int ) ;
static int           Msh2DQuadNeighbor( Msh2DQuadTree, int, int,
                     Msh2DQuadTree *, int *, int *) ;
static int           Msh2DQuadSize( Msh2DQuadTree, double *, double *,
                     double *, double *, int, double *, Msh2DQuadTree *, int *) ;
static Msh2DQuadTree Msh2DQuadAlloc() ;
static Msh2DQuadTree Msh2DRefQuad( Msh2DQuadTree, Msh2DQuadTree, int, int) ;

static double        SurfDist3D_00(double *,double *,
                     void (*f_surf)(double,double,double *));
static void          SurfQuadMetric( Msh2DQuadTree, double *, double *,
                     void (*f_surf)(double,double,double *));

static void          SurfQuadCurvature (Msh2DQuadTree, double *, double *, int);
#if OLD_IMP
static void          SurfQuadMetricCurvature (Msh2DQuadTree, int ,
                     Msh2DQuadTree, int, double);
#else


static void          SurfQuadMetricCurvature (Msh2DQuadTree, int , double,
                     Msh2DQuadTree, int, double);
#endif

static Msh2DQuadTree Msh2DRefSizeQuad ( Msh2DQuadTree, Msh2DQuadTree, int, int, double) ;

#if OLD_IMP

#else
static void          SurfQuadRefineCurvature (void);
#endif

#ifdef QUAD_DRAW
 static void         Msh2DQuadDraw( Msh2DQuadTree, double *, double *) ;
 static void         Msh2DQuadDrawM     ( Msh2DQuadTree, double *, double *, char *file) ;
#endif



/* -------------------------------------------------------------------
** Msh2Dk, m, InternalNodes - find the location of internal nodes.
*/

int MshSurfGenQuadTree
(
int     num_org_nodes,
int     num_org_edges,
double  original_nodes[][2],
int     original_edges[][2],
void (*f_surf)(double u,double v, double *f_dev)
)
{
    double  max[2], min[2], mid[2], span[2] ;
    double  coords[2][2];
    double  length, dist_par ;
    double  A[2], B[2];
    double  l_left, l_right, l_top, l_botton, tmp_lenght;
    double  p_left, p_right, p_top, p_botton;
    double  cx, cy, Rx, Ry;
    int     i, j, k, min_level;
    int     bound_min_level = (int) INFPOS ;
    double  ref_factor = REF_FACTOR;
    Msh2DQuadTree tree = 0;
    Ref_curvature *curr_curv;
    long int cpu_time;

    /* init angle */
   mshsurfAngle = MshSurfGetMaxAngle();
   mshsurfCosAngle = cos (mshsurfAngle);



	/* 1.1 Find the max and min x and y coordinates.  Determine the
           smallest square that will contain all the nodes, and center it
           on the domain */

    for ( i=0 ; i<2 ; i++ ) min[i] = max[i] = original_nodes[0][i] ;
    for ( j=0 ; j<num_org_nodes ; j++ )
    {
        for ( i=0 ; i<2 ; i++ )
        {
            max[i] = (max[i] > original_nodes[j][i]) ?
                     max[i] : original_nodes[j][i] ;
            min[i] = (min[i] < original_nodes[j][i]) ?
                     min[i] : original_nodes[j][i] ;
        }
    }
    for ( i=0 ; i<2 ; i++ )
    {
        bmin[i] = min[i] ;
        bmax[i] = max[i] ;
    }
    for ( i=0 ; i<2 ; i++ )
    {
        mid[i] = (max[i]+min[i]) / 2.0 ;
        span[i] = max[i] - min[i] ;
    }

    /* calcula uma QTree aproximada e as distorcoes em x e y */
    /* calulo o Box que engloba toda estrutura */
    l_left = l_right = l_top = l_botton = 0.0;
    p_left = p_right = p_top = p_botton = 0.0;
    for (j=0; j<num_org_nodes ; j++ )
    {
      /* distancia na horizontal */
      A[0] = mid[0];  A[1] = original_nodes[j][1];
      B[0] = original_nodes[j][0];  B[1] = original_nodes[j][1];
      tmp_lenght = SurfDist3D_00(A,B,f_surf);
      dist_par = sqrt((A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]));
      if ((tmp_lenght > l_right) && (B[0] > A[0]))  /* lado direito */
      {
        l_right = tmp_lenght;
        p_right = dist_par;
      }
      if ((tmp_lenght > l_left) && (B[0] < A[0]))  /* lado esquerdo */
      {
        l_left = tmp_lenght;
        p_left = dist_par;
      }

      /* distancia na vertical */
      A[0] = original_nodes[j][0];  A[1] = mid[1];
      B[0] = original_nodes[j][0];  B[1] = original_nodes[j][1];
      tmp_lenght = SurfDist3D_00(A,B,f_surf);
      dist_par = sqrt((A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]));
      if ((tmp_lenght > l_top) && (B[1] > A[1]))  /* lado de cima */
      {
        l_top = tmp_lenght;
        p_top = dist_par;
      }
      if ((tmp_lenght > l_botton) && (B[1] < A[1]))  /* lado de baixo */
      {
        l_botton = tmp_lenght;
        p_botton = dist_par;
      }
    }

    /* distorcoes aproximadas */
    span[0] = (l_right + l_left); /* span[0] = 2 * max(l_right,l_left), para garantir englobar tudo */
    Distor_x = Rx = span[0] / (p_right + p_left); /*  (l_right + l_left) / (p_right + p_left) */
    span[1] = (l_top + l_botton) ; /* span[1] = 2 * max(l_right,l_left), para garantir englobar tudo */
    Distor_y = Ry = span[1] / (p_botton + p_top); /*  (l_top + l_botton) / (p_botton + p_top) */

    if ( span[0] >= span[1] )
    {
        min[0] = mid[0] - (span[0]/2.0) ;
        max[0] = mid[0] + (span[0]/2.0) ;
        min[1] = mid[1] - (span[0]/2.0) ;
        max[1] = mid[1] + (span[0]/2.0) ;
    }
    else
    {
        min[0] = mid[0] - (span[1]/2.0) ;
        max[0] = mid[0] + (span[1]/2.0) ;
        min[1] = mid[1] - (span[1]/2.0) ;
        max[1] = mid[1] + (span[1]/2.0) ;
    }

    /* fill in box structure */
    for ( i=0 ; i<2 ; i++ )
    {
       box.min[i] = gmin[i] = min[i];
       box.max[i] = gmax[i] = max[i];
       box.center[i] = mid[i];
    }


    /* 1.2 For each edge in the boundary, find the mid-point of the edge
           and it's length.  Locate the mid-point in the quadtree and look
           at the cell size.  If the cell size is approximately equal
           to or smaller than the edge size, go on to the next edge.
           If not, subdivide the quadant. */

    for ( i=0 ; i<num_org_edges ; i++ )
    {
        for ( j=0 ; j<2 ; j++ )
        {
            for ( k=0 ; k<2 ; k++ )
            {
                coords[j][k] = original_nodes[original_edges[i][j]][k] ;
            }
        }

        /* compute the center */

        cx = (coords[0][0] + coords[1][0])/2.0;
        cy = (coords[0][1] + coords[1][1])/2.0;

        mid[0] = (cx - box.center[0])*Rx + box.center[0] ;
        mid[1] = (cy - box.center[1])*Ry + box.center[1] ;

        /* compute the length  */
        length = SurfDist3D_00(coords[0],coords[1],f_surf);

        min_level = -1 ;
        if (length > 0)
        {
          tree = Msh2DBdryQuad (tree, 0, length, mid, max, min, ref_factor,
                                1, &min_level ) ;

          if( min_level < bound_min_level )
            bound_min_level = min_level ;
        }
    }

    if ( !tree )
    {
        tree = Msh2DQuadAlloc() ;
        tree->parent = 0 ;
        tree->marks  = 0 ;
        tree->size   = 0.0;
    }

#ifdef QUAD_DRAW
    Msh2DQuadDrawM( tree, max, min, "quadtree_init.m" );
#ifdef USE_IUP
    IupMessage("Pause","QTreeINIT");
#endif
#endif


    /* 1.3 Find the largest cell size that contains a boundary edge.
           Refine the entire quadtree to at least this level. */
   if (MshSurfRefineByMaxEdgeSize_ask ())
   {
     tree = Msh2DRefQuad( tree, 0, 1, bound_min_level ) ;
   }

#ifdef QUAD_DRAW
   Msh2DQuadDrawM( tree, max, min, "quadtree_min_level_1.m" );
#ifdef USE_IUP
   IupMessage("Pause","QTreeINIT");
#endif
#endif

    /* 1.4 Refine the quadtree so that no adjacent cells differ by a
           refinement factor of more than one. */
    Msh2DTransQuad( tree, 0 ) ;

#ifdef QUAD_DRAW
    Msh2DQuadDrawM( tree, max, min, "quadtree_one_level_1.m" );
#ifdef USE_IUP
    IupMessage("Pause","QTreeINIT");
#endif
#endif

#if  QTREE_INTP
    /* 1.5 Compute the Riemanninan metric for each cells */
    cpu_time = clock( );
    SurfQuadMetric( tree, max, min, f_surf);
    cpu_time = (clock( ) - cpu_time)/CLOCKS_PER_SEC;
    printf("\n\t\tCPU time SurfQuadMetric........ %0.3f (s)\n", (double)cpu_time);
#endif


    /* Refine considering curvatures */
    if (MshSurfRefineCurvature_ask ( ))
    {
      root = tree;
      cpu_time = clock( );

      /* 1.6 Compute the curvature */
      init_curv = NULL;
      tail_curv = NULL;
      SurfQuadCurvature (tree, max, min, 0);

    /* and Refine */
#if OLD_IMP
      for (curr_curv = init_curv; curr_curv != NULL; curr_curv = curr_curv->next)
      {
         mid[0] = curr_curv->mid[0];
         mid[1] = curr_curv->mid[1];

         min_level = -1 ;
         if (curr_curv->size > 0)
         {
           tree = Msh2DBdryQuad (tree, 0, curr_curv->size/2.0, mid, max, min,
                                  ref_factor, 1, &min_level ) ;
           if( min_level < bound_min_level )
             bound_min_level = min_level ;
         }
      }
#else
      SurfQuadRefineCurvature ( );
#endif


      cpu_time = (clock( ) - cpu_time)/CLOCKS_PER_SEC;
      printf("\n\t\tCPU time SurfQuadCurvature........ %0.3f (s)", (double)cpu_time);

#ifdef QUAD_DRAW
      Msh2DQuadDrawM( tree, max, min, "quadtree_refine_curvature.m" );
#ifdef USE_IUP
      IupMessage("Pause","QTreeINIT");
#endif
#endif

      /* Release auxiliar memory */
      while (init_curv != NULL)
      {
        curr_curv = init_curv;
        init_curv = curr_curv->next;
        free (curr_curv);
      };

      /* 1.7 Refine the quadtree so that no adjacent cells differ by a
             refinement factor of more than one. */
      Msh2DTransQuad( tree, 0 ) ;
#ifdef QUAD_DRAW
      Msh2DQuadDrawM( tree, max, min, "quadtree_one_level_2.m" );
#ifdef USE_IUP
      IupMessage("Pause","QTreeINIT");
#endif
#endif

#if  QTREE_INTP
      /* 1.5 Compute the Riemanninan metric for each cells */
      cpu_time = clock( );
      SurfQuadMetric( tree, max, min, f_surf);
      cpu_time = (clock( ) - cpu_time)/CLOCKS_PER_SEC;
      printf("\n\t\tCPU time SurfQuadMetric........ %0.3f (s)", (double)cpu_time);
#endif

    }


    root = tree ;

    return(1) ;
}


/* -------------------------------------------------------------------
** MshSurfSetMaxQuadSize
*/
void MshSurfSetMaxQuadSize (double max_size)
{
  if( root == NULL )
  	return;

  /* Find the largest cell size that contains a boundary face.
     Refine the entire quadtree to at least this level. */
  if (max_size > 0)
  {
    int    level;
    double cur_size = box.max[0] - box.min[0];
    for (level = 0; cur_size > max_size; level++)
       cur_size *= 0.5;
    root = Msh2DRefSizeQuad (root, 0, 1, level, max_size);
  }

  /* Refine the quadtree so that no adjacent cells differ by a
     refinement factor of more than one. */
  Msh2DTransQuad (root, 0);
}


/* -------------------------------------------------------------------
** MshSurfQuadFreeAll - free quadtree
*/
void MshSurfQuadFreeAll()
{
    Msh2DQuadTree  cur, next ;

    /* free all blocks allocated to store edge information */

    if ( quad_block_ptr ) cur = quad_block_ptr ;
    else return ;

    while ( cur->child[0] ) {
        next = cur->child[0] ;
        free( cur ) ;
        cur = next ;
    }
    free( cur ) ;

    quad_free = 0 ;
    quad_block_ptr = 0 ;

	root = 0;
}


/* -------------------------------------------------------------------
** Msh2DOptimalNodes - find the size of the cell where the given edge
**                     center is to use for the location of optimal nodes.
*/

double MshSurfOptimalNodes
(
double  edge_center[2],
double  tree_center[2],
int   *tree_level
)
{
    int    level, index ;
    double  size ;
    double new_edge_center[2];
    Msh2DQuadTree  leaf;

    if( root == NULL ) return 0.0 ;

    new_edge_center[0] = (edge_center[0] - box.center[0])*Distor_x + box.center[0] ;
    new_edge_center[1] = (edge_center[1] - box.center[1])*Distor_y + box.center[1] ;

    /* 1.1 Find the size of the cell where the given point is */
    level = Msh2DQuadSize( root, gmax, gmin, new_edge_center,
                           tree_center, 1, &size, &leaf, &index ) ;
    *tree_level = level ;

    /* 1.2 Return the found size */

    return size ;
}

/*-------------------------------------------------------------------
** SurfQuadTreeDeriv - Find derivad of the cell where is the 2D point
*/

void MshSurfQuadTreeDeriv
(
double    u,
double    v,
double    *dev,
int       interpol
)
{
    double  new_coord[2];
    double  size;
    /*int     level;*/
    double  tree_center[2];
    Msh2DQuadTree  leaf;
    int     i, j, cont = 0, index;

    if( root == NULL ) return;

    new_coord[0] = (u - box.center[0])*Distor_x + box.center[0] ;
    new_coord[1] = (v - box.center[1])*Distor_y + box.center[1] ;

    /* 1.1 Find the cell where the given point is */
    /*level = */Msh2DQuadSize( root, gmax, gmin, new_coord,
                           tree_center, 1, &size, &leaf, &index);

    if (!interpol)
      for ( j = 0; j < 6; j++) dev[j] = leaf->m[index]->dev[j];
    else
    {
      for (i = 0; i < 4; i++)
      {
        if (leaf->m[i] != NULL)
        {
          for(j=0; j<6; j++) dev[j] = leaf->m[i]->dev[j];
          cont++;
        }
      }
      for ( i = 0; i < 6; i++) dev[i] = dev[i]/cont;
    }
}

/* -------------------------------------------------------------------
** Msh2DBdryQuad - subdivide the quad tree until the length of the sides
**                 of the quad containing the given point is about
**                 the same as the length associated with the point.
*/

static Msh2DQuadTree Msh2DBdryQuad
(
Msh2DQuadTree  tree,
Msh2DQuadTree  parent,
double         length,
double         coord[2],
double         max[2],
double         min[2],
double         factor,
int            this_level,
int           *min_level
)
{
    double         clength ;
    int            nindx, indx[2], i ;
    double         lmin[2][2], lmax[2][2], mid[2] ;

    /* find the size of the cell and the current length */

    clength = max[0] - min[0] ;

    if ( length > clength * factor )
    {
        if (( *min_level == -1 ) || ((this_level-1) < *min_level ))
            *min_level = this_level-1 ;

      if( parent != NULL )
        parent->size = length;

      return( tree ) ;
    }

    /* if we get here we descend the tree */

    if ( !tree )
    {
        tree = Msh2DQuadAlloc() ;
        tree->parent = parent ;
        tree->marks  = 0 ;
        tree->size   = 0.0;
    }

    for ( i=0 ; i<2 ; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

    if ( coord[0] == mid[0] || coord[1] == mid[1] )
      nindx = 2 ;
    else
      nindx = 1 ;
    indx[0] = indx[1] = 0 ;

    for ( i=0 ; i<2 ; i++ )
    {
        if ( coord[i] == mid[i] )
        {
            indx[0] |= (1<<i) ; lmin[0][i] = mid[i] ; lmax[0][i] = max[i] ;
            lmin[1][i] = min[i] ; lmax[1][i] = mid[i] ;
        }
        else if ( coord[i] > mid[i] )
        {
             indx[0] |= (1<<i) ;  lmin[0][i] = mid[i] ;  lmax[0][i] = max[i] ;
             indx[1] |= (1<<i) ;  lmin[1][i] = mid[i] ;  lmax[1][i] = max[i] ;
        }
        else
        {
            lmin[0][i] = min[i] ;  lmax[0][i] = mid[i] ;
            lmin[1][i] = min[i] ;  lmax[1][i] = mid[i] ;
        }

    }

    for (i = 0; i < nindx ; i++)
     tree->child[indx[i]] = Msh2DBdryQuad( tree->child[indx[i]], tree, length,
                             coord, lmax[i], lmin[i],
                             factor, this_level+1, min_level ) ;

    return( tree ) ;
}


/* -------------------------------------------------------------------
** Msh2DRefQuad - subdivide the quad tree to at least the specified level.
*/

static Msh2DQuadTree Msh2DRefQuad
(
Msh2DQuadTree  tree,
Msh2DQuadTree  parent,
int            this_level,
int            min_level
)
{
    int j ;

    if ( this_level > min_level ) return( tree ) ;
    if ( tree == NULL)
    {
        tree = Msh2DQuadAlloc() ;
        tree->parent = parent ;
        tree->size   = 0.0;
    }

    for ( j=0 ; j<4 ; j++ )
    {
        tree->child[j] = Msh2DRefQuad( tree->child[j], tree, this_level+1, min_level ) ;
    }
    return( tree ) ;
}




/* -------------------------------------------------------------------
** Msh2DTransQuad - subdivide the quad tree so that no two neighbor cells
**                  differ in level by more than one.
*/

static int direct_tbl[4][2] = {
    { 2, 3 },
    { 0, 3 },
    { 2, 1 },
    { 0, 1 } } ;

static void Msh2DTransQuad
(
Msh2DQuadTree  tree,
int            level
)
{
    int i, j,direction, clevel, nchild ;
    Msh2DQuadTree  neighbor ;

    for ( j=0 ; j<4 ; j++ )
    {
        if ( tree->child[j] )
        {
            Msh2DTransQuad( tree->child[j], level+1 ) ;
        }
        else
        {
            for ( i=0 ; i<2 ; i++ )
            {

                direction = direct_tbl[j][i] ;

                clevel = level - 1 ;
                Msh2DQuadNeighbor( tree, j, direction, &neighbor,
                                         &nchild, &clevel ) ;

                /* if the neighbors differ by more than one, subdivide */

                if ( neighbor && ((level-clevel) > 1) )
                {
                    neighbor->child[nchild] = Msh2DQuadAlloc() ;
                    neighbor->child[nchild]->size   = 0.0;
                    neighbor->child[nchild]->parent = neighbor ;

                    Msh2DTransQuad( neighbor->child[nchild], clevel+1 ) ;
                }
            }
        }
    }
}

/* -------------------------------------------------------------------
** Msh2DQuadNeighbor - given a node in the quadtree along with a search
**                     direction, this routine finds the neighbor cell
**                     in the quad tree.
*/

static int dir_tbl[4][4] = {
    { 0, 0, 1, 1 },
    { 1, 0, 0, 1 },
    { 0, 1, 1, 0 },
    { 1, 1, 0, 0 } } ;

static int mirror_tbl[4][4] = {
    { 1, 2, 1, 2 },
    { 0, 3, 0, 3 },
    { 3, 0, 3, 0 },
    { 2, 1, 2, 1 } } ;

static int Msh2DQuadNeighbor
(
Msh2DQuadTree  cell,
int            child,
int            direction,
Msh2DQuadTree  *cousin,
int            *cchild,
int            *level
)
{
    Msh2DQuadTree  parent ;
    int           indx, found = 0 ;

    /* first move up the tree until we find a common ancester
       we do this by finding the first time we find a decendent
       in the direction opposite the search direction */

    parent = cell->parent ;
    *level -= 1 ;
    if ( !parent )
    {
        *cousin = 0 ;
        return(1) ;
    }

    for ( indx=0 ; indx<4 ; indx++ )
    {
        if ( parent->child[indx] == cell ) break ;
    }

    if ( dir_tbl[indx][direction] )
    {
        found = Msh2DQuadNeighbor( parent, indx, direction, cousin, cchild, level ) ;
    }
    else
    {
        *cousin = parent ;
        *cchild = mirror_tbl[indx][direction] ;
        *level += 1 ;
    }

    /* decend the tree in the mirror direction */

    if ( !found )
    {
        child = mirror_tbl[indx][direction] ;
        if ( (*cousin)->child[child] )
        {
            *level += 1 ;
            *cchild = child ;
            *cousin = (*cousin)->child[child] ;
            return ( 0 ) ;
        }
        else
            *cchild = child ;
    }

    return ( 1 ) ;
}

/* -------------------------------------------------------------------
** Msh2DQuadSize - given the coordinates of a point, this routine finds
**                the terminal cell that contains the point, and returns
**                the size of this cell.
*/

static int Msh2DQuadSize
(
Msh2DQuadTree  tree,
double         max[2],
double         min[2],
double         edge_coord[2],
double         tree_coord[2],
int            tree_level,
double         *size,
Msh2DQuadTree  *leaf,
int            *index
)
{
    int            indx, i ;
    double         lmin[2], lmax[2], mid[2] ;

    /* find the mid point and the appropriate child cell */

    indx = 0 ;
    for ( i=0 ; i<2 ; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

    for ( i=0 ; i<2 ; i++ )
    {
        if ( edge_coord[i] > mid[i] )
        {
            indx |= (1<<i) ;  lmin[i] = mid[i] ;  lmax[i] = max[i] ;
        }
        else
        {
            lmin[i] = min[i] ;  lmax[i] = mid[i] ;
        }
    }

    /* check to see if the cell has children.  If so, descend the tree.
       If not, set it's flag */

    if ( tree->child[indx] )
    {
        Msh2DQuadSize( tree->child[indx], lmax, lmin, edge_coord,
                       tree_coord, tree_level+1, size, leaf, index ) ;
    }
    else
    {
        /* tree center and size */

        if ( !(tree->flags & (1<<indx)) )
        {
          if (tree->size != 0.0)
            *size = tree->size;
          else
         *size = lmax[0] - lmin[0] ;
         tree_coord[0] = (lmin[0]+lmax[0])/2.0 ;
         tree_coord[1] = (lmin[1]+lmax[1])/2.0 ;
         *leaf = tree;
         *index = indx;
        }
        else
        {
         *size = 0.0 ;
         tree_coord[0] = edge_coord[0] ;
         tree_coord[1] = edge_coord[1] ;
         *index = -1;
        }
        return tree_level ;
    }

    return 1;
}


static void get_normal (Msh2DQuadTree tree, int i)
{
  double *A = &(tree->m[i]->dev[0]);
  double *B = &(tree->m[i]->dev[3]);
  double C[3], lenC;


  C[0] = A[1]*B[2] - A[2]*B[1];
  C[1] = A[2]*B[0] - A[0]*B[2];
  C[2] = A[0]*B[1] - A[1]*B[0];

  lenC = sqrt(C[0]*C[0] + C[1]*C[1] + C[2]*C[2]);

  tree->m[i]->normal[0] = C[0] / lenC;
  tree->m[i]->normal[1] = C[1] / lenC;
  tree->m[i]->normal[2] = C[2] / lenC;

}




/* -------------------------------------------------------------------
** SurfQuadMetric - Evaluate the derived on center of quadtree leaf.
*/

static void SurfQuadMetric
(
Msh2DQuadTree tree,
double        max[2],
double        min[2],
void (*f_surf)(double u,double v,double *f_dev)
)
{
    int           indx, i;
    /* int           mark ; */
    double        lmin[2], lmax[2], mid[2];
    double        par_u, par_v;

    /* find the mid point of the cell */

    for ( i=0; i<2; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

    /* visit all the child cells */

    for ( indx=0; indx<4; indx++ )
    {
     tree->m[indx] = NULL;
     if ( tree->child[indx] )
     {
         for ( i=0 ; i<2 ; i++ )
         {
           if ( indx & (1<<i) )
           {
              lmin[i] = mid[i] ;  lmax[i] = max[i] ;
           }
           else
           {
              lmin[i] = min[i] ;  lmax[i] = mid[i] ;
           }
         }
         SurfQuadMetric( tree->child[indx], lmax, lmin, f_surf ) ;
     }
     else
     {
         /* initiate mark for this terminal cell */
        /* mark = 0 ; */

        /* find the limits of the terminal cell */
        for( i=0; i<2; i++ )
        {
           if ( indx & (1<<i) )
           {
              lmin[i] = mid[i] ; lmax[i] = max[i] ;
           }
           else
           {
              lmin[i] = min[i] ; lmax[i] = mid[i] ;
           }
        }

        if (tree->m[indx] == NULL)
        {
          tree->m[indx] = (Metric *) calloc (1, sizeof (Metric));

          tree->m[indx]->center[0] = (lmin[0]+lmax[0]) / 2.0;
          tree->m[indx]->center[1] = (lmin[1]+lmax[1]) / 2.0;
          par_u = (tree->m[indx]->center[0] - box.center[0])/Distor_x +
          box.center[0];
          par_v = (tree->m[indx]->center[1] - box.center[1])/Distor_y +
          box.center[1];
          (*f_surf) (par_u, par_v, tree->m[indx]->dev);
          get_normal (tree, indx);
          tree->m[indx]->choose = 0;
        }
    }
 }

}



#ifdef QUAD_DRAW

static double tt=0.0;
/* -------------------------------------------------------------------
** Msh2DQuadDraw -
*/

static void Msh2DQuadDraw
(
Msh2DQuadTree tree,
double        max[2],
double        min[2]
)
{
    int           indx, i;
    /* int           mark ; */
    double        lmin[2], lmax[2], mid[2];

    /* find the mid point of the cell */

    for ( i=0; i<2; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

    /* visit all the child cells */

    for ( indx=0; indx<4; indx++ )
    {
     if ( tree->child[indx] )
     {
         for ( i=0 ; i<2 ; i++ )
         {
           if ( indx & (1<<i) )
           {
              lmin[i] = mid[i] ;  lmax[i] = max[i] ;
           }
           else
           {
              lmin[i] = min[i] ;  lmax[i] = mid[i] ;
           }
         }
         Msh2DQuadDraw( tree->child[indx], lmax, lmin ) ;
     }
     else
     {
	    	 /* initiate mark for this terminal cell */
        /* mark = 0 ; */

	      /* find the limits of the terminal cell */
        for( i=0; i<2; i++ )
        {
           if ( indx & (1<<i) )
           {
              lmin[i] = mid[i] ; lmax[i] = max[i] ;
           }
           else
           {
              lmin[i] = min[i] ; lmax[i] = mid[i] ;
           }
        }

#if 0
        G3dForeground(G3D_YELLOW);
	G3dLine (lmin[0]+tt,lmin[1]+tt,0.0,lmin[0]+tt,lmax[1]-tt,0.0);
        G3dLine (lmin[0]+tt,lmax[1]-tt,0.0,lmax[0]-tt,lmax[1]-tt,0.0);
	G3dLine (lmax[0]-tt,lmax[1]-tt,0.0,lmax[0]-tt,lmin[1]+tt,0.0);
        G3dLine (lmax[0]-tt,lmin[1]+tt,0.0,lmin[0]+tt,lmin[1]+tt,0.0);
#endif
    }
 }
}


/* -------------------------------------------------------------------*/
static void Msh2DQuadDrawMCell
(
 Msh2DQuadTree tree,
 double        max[2],
 double        min[2],
 FILE          *ptr
 )
 {
   int           indx, i, mark ;
   double        lmin[2], lmax[2], mid[2];

   /* find the mid point of the cell */

   for ( i=0; i<2; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

   /* visit all the child cells */

   for ( indx=0; indx<4; indx++ )
   {
     if ( tree->child[indx] )
     {
       for ( i=0 ; i<2 ; i++ )
       {
         if ( indx & (1<<i) )
         {
           lmin[i] = mid[i] ;  lmax[i] = max[i] ;
         }
         else
         {
           lmin[i] = min[i] ;  lmax[i] = mid[i] ;
         }
       }
       Msh2DQuadDrawMCell( tree->child[indx], lmax, lmin, ptr) ;
     }
     else
     {
       /* initiate mark for this terminal cell */
       mark = 0 ;

       /* find the limits of the terminal cell */
       for( i=0; i<2; i++ )
       {
         if ( indx & (1<<i) )
         {
           lmin[i] = mid[i] ; lmax[i] = max[i] ;
         }
         else
         {
           lmin[i] = min[i] ; lmax[i] = mid[i] ;
         }
       }

       fprintf (ptr, "line ([%f, %f], [%f, %f])\n", lmin[0]+tt, lmin[0]+tt, lmin[1]+tt, lmax[1]-tt);
       fprintf (ptr, "line ([%f, %f], [%f, %f])\n", lmin[0]+tt, lmax[0]-tt, lmax[1]-tt, lmax[1]-tt);
       fprintf (ptr, "line ([%f, %f], [%f, %f])\n", lmax[0]-tt, lmax[0]-tt, lmax[1]-tt, lmin[1]+tt);
       fprintf (ptr, "line ([%f, %f], [%f, %f])\n", lmax[0]-tt, lmin[0]+tt, lmin[1]+tt, lmin[1]+tt);

    }
 }
}

/* -------------------------------------------------------------------*/
static void Msh2DQuadDrawM
(
 Msh2DQuadTree tree,
 double        max[2],
 double        min[2],
 char          *file
 )
 {
   FILE *ptr = fopen(file, "wt");
   Msh2DQuadDrawMCell (tree, max, min, ptr);
   fclose(ptr);
 }

#endif


/* -------------------------------------------------------------------
** Msh2DQuadAlloc - these routines manage the allocation and freeing
**                 quad tree entries.
*/


static Msh2DQuadTree Msh2DQuadAlloc( )
{
    Msh2DQuadTree  new_block, alloced ;
    int           i ;

    /* if the free pointer is null we need to allocate a new block
       of boundary nodes */

    if ( !quad_free )
    {
        new_block = (Msh2DQuadTree)malloc(
           MSH2D_QUAD_TREE_QUANTUM * sizeof(Msh2DQuadTreeRec) ) ;
/*        memset (new_block, 0, MSH2D_QUAD_TREE_QUANTUM * sizeof(Msh2DQuadTreeRec)); */

        new_block[0].child[0] = quad_block_ptr ;
        quad_block_ptr = new_block ;
        for ( i=1 ; i<(MSH2D_QUAD_TREE_QUANTUM-1) ; i++ )
        {
            new_block[i].child[0] = &(new_block[i+1]) ;
        }
        new_block[MSH2D_QUAD_TREE_QUANTUM-1].child[0] = 0 ;
        quad_free = &(new_block[1]) ;
    }

    /* return the top thing on the free list */

    alloced = quad_free ;
    quad_free = quad_free->child[0] ;

    for ( i=0 ; i<4 ; i++ ) alloced->child[i] = 0 ;
    alloced->parent = 0 ;
    alloced->flags = 0 ;
    alloced->m[0] = NULL;
    alloced->m[1] = NULL;
    alloced->m[2] = NULL;
    alloced->m[3] = NULL;
    return( alloced ) ;
}


static double SurfDist3D_00
(
double  *A,
double  *B,
void (*f_surf)(double u,double v,double *f_dev)
)
{
 double A_E, A_F, A_G;
 double B_E, B_F, B_G;
 double distAB, distBA;
 double dev[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
 double du = (A[0] - B[0]);
 double dv = (A[1] - B[1]);

 (*f_surf)(A[0],A[1],dev);
 A_E = dev[0]*dev[0] + dev[1]*dev[1] + dev[2]*dev[2];
 A_F = dev[0]*dev[3] + dev[1]*dev[4] + dev[2]*dev[5];
 A_G = dev[3]*dev[3] + dev[4]*dev[4] + dev[5]*dev[5];

 distAB = sqrt(A_E*du*du + 2*A_F*du*dv + A_G*dv*dv);

 (*f_surf)(B[0],B[1],dev);
 B_E = dev[0]*dev[0] + dev[1]*dev[1] + dev[2]*dev[2];
 B_F = dev[0]*dev[3] + dev[1]*dev[4] + dev[2]*dev[5];
 B_G = dev[3]*dev[3] + dev[4]*dev[4] + dev[5]*dev[5];

 distBA = sqrt(B_E*du*du + 2*B_F*du*dv + B_G*dv*dv);

 return((distAB+distBA)/2.0);
}





/**************************************************************************/
/**************************                    ****************************/
/**************************  C U R V A T U R A ****************************/
/**************************                    ****************************/
/**************************************************************************/


/* -------------------------------------------------------------------------
** SurfQuadNeighborCurvature -
*/

static void SurfQuadNeighborCurvature
(
Msh2DQuadTree tree,
int           indx,
double        size
)
{
    int            nchild ;
    Msh2DQuadTree  neighbor ;
    double         dir_x[2], dir_y[2], tree_center[2];
    /*int            level;*/

#if OLD_IMP
    dir_x[0] = tree->m[indx]->center[0] + size/2.0 + size/10.0;
    dir_x[1] = tree->m[indx]->center[1];
    dir_y[0] = tree->m[indx]->center[0];
    dir_y[1] = tree->m[indx]->center[1] + size/2.0 + size/10.0;

    /* 1 Find the neighbor cell in x diretion */
    /*level = */Msh2DQuadSize (root, gmax, gmin, dir_x, tree_center,
                           1, &size, &neighbor, &nchild);
    /* Check the curvature */
    SurfQuadMetricCurvature (tree, indx, neighbor, nchild, size);

    /* 2 Find the neighbor cell in y diretion */
    /*level = */Msh2DQuadSize (root, gmax, gmin, dir_y, tree_center,
                           1, &size, &neighbor, &nchild);
    /* Check the curvature */
    SurfQuadMetricCurvature (tree, indx, neighbor, nchild, size);

#else
    double         neigSize;

    dir_x[0] = tree->m[indx]->center[0] + size/2.0 + size/10.0;
    dir_x[1] = tree->m[indx]->center[1];
    dir_y[0] = tree->m[indx]->center[0];
    dir_y[1] = tree->m[indx]->center[1] + size/2.0 + size/10.0;

    /* 1 Find the neighbor cell in x direction */
    /*level = */Msh2DQuadSize (root, gmax, gmin, dir_x, tree_center,
                           1, &neigSize, &neighbor, &nchild);
    /* Check the curvature */
    SurfQuadMetricCurvature (tree, indx, size, neighbor, nchild, neigSize);

    /* 2 Find the neighbor cell in y direction */
    /*level = */Msh2DQuadSize (root, gmax, gmin, dir_y, tree_center,
                           1, &neigSize, &neighbor, &nchild);
    /* Check the curvature */
    SurfQuadMetricCurvature (tree, indx, size, neighbor, nchild, neigSize);
#endif
}



static void SurfQuadCurvature
(
Msh2DQuadTree  tree,
double        max[2],
double        min[2],
int           level
)
{
    int           i, j;
    double        lmin[2], lmax[2], mid[2];

    /* find the mid point of the cell */
    for ( i=0; i<2; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

    for ( j=0 ; j<4 ; j++ )
    {
       if ( tree->child[j] )
       {
          for ( i=0 ; i<2 ; i++ )
          {
             if ( j & (1<<i) )
             {
               lmin[i] = mid[i] ;  lmax[i] = max[i] ;
             }
             else
             {
               lmin[i] = min[i] ;  lmax[i] = mid[i] ;
             }
          }
          SurfQuadCurvature( tree->child[j], lmax, lmin, level+1) ;
       }
       else
       {
         /* find the limits of the terminal cell */
         for( i=0; i<2; i++ )
         {
           if ( j & (1<<i) )
           {
              lmin[i] = mid[i] ; lmax[i] = max[i] ;
           }
           else
           {
              lmin[i] = min[i] ; lmax[i] = mid[i] ;
           }
         }

#if 0
         G3dForeground(G3D_BLUE);
	 G3dMark ((lmin[0]+lmax[0])/2,(lmin[1]+lmax[1])/2,0.0);
         G3dForeground(G3D_WHITE);
	 G3dMark (tree->m[j]->center[0], tree->m[j]->center[1], 0.0);
#endif

         SurfQuadNeighborCurvature (tree, j, lmax[0]-lmin[0]);
       }
    }
}


/* -------------------------------------------------------------------------
** SurfQuadMetricCurvature - Evaluate the curvature cells of quadtree.
*/

#if OLD_IMP

static void SurfQuadMetricCurvature
(
Msh2DQuadTree tree1,
int           i,
Msh2DQuadTree tree2,
int           j,
double size
)
{
  double *A = tree1->m[i]->normal;
  double *B = tree2->m[j]->normal;
  double lenA = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
  double lenB = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
  double cos_AB;

  double A_E, A_F, A_G;
  double B_E, B_F, B_G;
  double distAB, distBA, dist_total;
  double *dev;
  double du, dv;
  double par_u1, par_u2, par_v1, par_v2;

/*  FILE *ptr_file = fopen ("SenosAB.dat", "a"); */

  cos_AB = (A[0]*B[0] + A[1]*B[1] + A[2]*B[2]) / (lenA*lenB);

  /* calculo da distancia */
  par_u1 = (tree1->m[i]->center[0] - box.center[0])/Distor_x + box.center[0];
  par_v1 = (tree1->m[i]->center[1] - box.center[1])/Distor_y + box.center[1];

  par_u2 = (tree2->m[j]->center[0] - box.center[0])/Distor_x + box.center[0];
  par_v2 = (tree2->m[j]->center[1] - box.center[1])/Distor_y + box.center[1];

  du = par_u1 - par_u2; dv = par_v1 - par_v2;

  dev = tree1->m[i]->dev;
  A_E = dev[0]*dev[0] + dev[1]*dev[1] + dev[2]*dev[2];
  A_F = dev[0]*dev[3] + dev[1]*dev[4] + dev[2]*dev[5];
  A_G = dev[3]*dev[3] + dev[4]*dev[4] + dev[5]*dev[5];

  dev = tree2->m[j]->dev;
  B_E = dev[0]*dev[0] + dev[1]*dev[1] + dev[2]*dev[2];
  B_F = dev[0]*dev[3] + dev[1]*dev[4] + dev[2]*dev[5];
  B_G = dev[3]*dev[3] + dev[4]*dev[4] + dev[5]*dev[5];

  distAB = sqrt(A_E*du*du + 2*A_F*du*dv + A_G*dv*dv);
  distBA = sqrt(B_E*du*du + 2*B_F*du*dv + B_G*dv*dv);

  dist_total = (distAB + distBA) / 2.0;


  /* calcula o tamanho ideal */
  if ( cos_AB < mshsurfCosAngle)
  {
    double r = dist_total / acos(cos_AB);
    if (tree1->m[i]->choose == 0)
    {
      Ref_curvature *new_curv = (Ref_curvature *) malloc (sizeof (Ref_curvature));

      new_curv->size = r * mshsurfAngle;
      new_curv->mid[0] = tree1->m[i]->center[0];
      new_curv->mid[1] = tree1->m[i]->center[1];

      if (tail_curv != NULL ) tail_curv->next = new_curv;
      tail_curv = new_curv;
      tail_curv->next = NULL;
      if( init_curv == NULL) init_curv = new_curv;
      tree1->m[i]->choose = 1;
    }
    if (tree2->m[j]->choose == 0)
    {
      Ref_curvature *new_curv = (Ref_curvature *) malloc (sizeof (Ref_curvature));

      new_curv->size = r * mshsurfAngle;
      new_curv->mid[0] = tree2->m[j]->center[0];
      new_curv->mid[1] = tree2->m[j]->center[1];

      if (tail_curv != NULL ) tail_curv->next = new_curv;
      tail_curv = new_curv;
      tail_curv->next = NULL;
      if( init_curv == NULL) init_curv = new_curv;
      tree2->m[j]->choose = 1;
    }
  }

}

#else

static void SurfQuadMetricCurvature
(
Msh2DQuadTree tree1,
int           i,
double        size1,
Msh2DQuadTree tree2,
int           j,
double        size2
)
{
  double *A = tree1->m[i]->normal;
  double *B = tree2->m[j]->normal;
  double lenA = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
  double lenB = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
  double C[3], lenC, sen_AB, cos_AB;

  double A_E, A_F, A_G;
  double B_E, B_F, B_G;
  double distAB, distBA, dist_total;
  double *dev;
  double du, dv;
  double par_u1, par_u2, par_v1, par_v2;

/*  FILE *ptr_file = fopen ("SenosAB.dat", "a"); */

  C[0] = A[1]*B[2] - A[2]*B[1];
  C[1] = A[2]*B[0] - A[0]*B[2];
  C[2] = A[0]*B[1] - A[1]*B[0];

  lenC = sqrt(C[0]*C[0] + C[1]*C[1] + C[2]*C[2]);

  cos_AB = (A[0]*B[0] + A[1]*B[1] + A[2]*B[2]) / (lenA*lenB);
  sen_AB = lenC / (lenA*lenB);


  /* calculo da distancia */
  par_u1 = (tree1->m[i]->center[0] - box.center[0])/Distor_x + box.center[0];
  par_v1 = (tree1->m[i]->center[1] - box.center[1])/Distor_y + box.center[1];

  par_u2 = (tree2->m[j]->center[0] - box.center[0])/Distor_x + box.center[0];
  par_v2 = (tree2->m[j]->center[1] - box.center[1])/Distor_y + box.center[1];

  du = par_u1 - par_u2; dv = par_v1 - par_v2;

  dev = tree1->m[i]->dev;
  A_E = dev[0]*dev[0] + dev[1]*dev[1] + dev[2]*dev[2];
  A_F = dev[0]*dev[3] + dev[1]*dev[4] + dev[2]*dev[5];
  A_G = dev[3]*dev[3] + dev[4]*dev[4] + dev[5]*dev[5];

  dev = tree2->m[j]->dev;
  B_E = dev[0]*dev[0] + dev[1]*dev[1] + dev[2]*dev[2];
  B_F = dev[0]*dev[3] + dev[1]*dev[4] + dev[2]*dev[5];
  B_G = dev[3]*dev[3] + dev[4]*dev[4] + dev[5]*dev[5];

  distAB = sqrt(A_E*du*du + 2*A_F*du*dv + A_G*dv*dv);
  distBA = sqrt(B_E*du*du + 2*B_F*du*dv + B_G*dv*dv);

  dist_total = (distAB + distBA) / 2.0;


  /* compute the ideal size */
  if ( cos_AB < mshsurfCosAngle)
  {
    double r = dist_total / acos(cos_AB);
    if (tree1->m[i]->choose == 0)
    {
      Ref_curvature *new_curv = (Ref_curvature *) malloc (sizeof (Ref_curvature));

      new_curv->newsize = r * mshsurfAngle;
      new_curv->tree = tree1;
      new_curv->oldsize = size1;

      if (tail_curv != NULL ) tail_curv->next = new_curv;
      tail_curv = new_curv;
      tail_curv->next = NULL;
      if( init_curv == NULL) init_curv = new_curv;
      tree1->m[i]->choose = 1;
    }
    if (tree2->m[j]->choose == 0)
    {
      Ref_curvature *new_curv = (Ref_curvature *) malloc (sizeof (Ref_curvature));

      new_curv->newsize = r * mshsurfAngle;
      new_curv->tree = tree2;
      new_curv->oldsize = size2;

      if (tail_curv != NULL ) tail_curv->next = new_curv;
      tail_curv = new_curv;
      tail_curv->next = NULL;
      if( init_curv == NULL) init_curv = new_curv;
      tree2->m[j]->choose = 1;
    }
  }

}

#endif

/* -------------------------------------------------------------------
** Msh2DRefSizeQuad - subdivide the oct tree to at least the specified level.
*/

static Msh2DQuadTree Msh2DRefSizeQuad (
Msh2DQuadTree  tree,
Msh2DQuadTree  parent,
int            this_level,
int            min_level,
double         size
)
{
  int j ;

  if (this_level > min_level)
  {
    if (parent != NULL)
    {
      if (parent->size == 0.0 || parent->size > size)
      {
        parent->size = size;
      }
    }
    return( tree ) ;
  }

  if ( !tree )
  {
    tree = Msh2DQuadAlloc() ;
    tree->parent = parent ;
    tree->size   = 0.0;
  }

  for (j = 0 ;j < 4 ;j++)
  {
    tree->child[j] = Msh2DRefSizeQuad (tree->child[j], tree, this_level+1, min_level, size);
  }
  return( tree ) ;
}



#if OLD_IMP

#else
/*
** SurfQuadRefineCurvature - Refine quadtree cell based on curvature
*/
static void SurfQuadRefineCurvature (void)
{
  Ref_curvature *curr_curv;
  int level;

  for (curr_curv = init_curv; curr_curv != NULL; curr_curv = curr_curv->next)
  {
    int min_level = -1 ;
    if (curr_curv->newsize > 0)
    {
      double cur_size = curr_curv->newsize;
      for (level = 0; cur_size > curr_curv->newsize; level++)
        cur_size *= 0.5;

      /* refine the current cell */
      root = Msh2DRefSizeQuad (curr_curv->tree, curr_curv->tree->parent,
                               1, level, curr_curv->newsize);
    }
  }

}

#endif
