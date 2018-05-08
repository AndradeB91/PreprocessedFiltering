/*
** -----------------------------------------------------------------
**
** surf3d_octree.c - Octree generation module.
**
** -----------------------------------------------------------------
**
** Description:
**  Package for implementing usual 3D octree operations based on
**  the boundary defined by triangular facets given as an input.
**
** Copyright:
**  (c) 1996 TECGRAF/PUC-Rio && CFG/Cornell University
**           All rights reserved
**
** History:
**  Created:        Oct-01-1996        Joaquim Bento Cavalcante Neto
**   Initial version. Created during my PhD Sandwich in Cornell.
**   It was based on Wash's implementation for octree operations.
**  Modified:        Oct-25-2001        Joaquim Bento Cavalcante Neto
**   Changed the use of variable fac to allow more control over the
**   degree of refinement of the mesh.
**  Modified:        Mar-26-2003        Joaquim Bento Cavalcante Neto
**   Created #defines starting with OCT_ for use of algorithms.
**   Also deleted #define BDR_OCTREE that was changed by OCT_FREE.
**  Modified:        Oct-23-2003        Joaquim Bento Cavalcante Neto
**   Modified variable bound_min_level to init as 1000000 instead
**   of INFPOS. Also changed function Msh3DOctSize to retun calling
**   (return Msh3DOctSize) instead of only (Msh3DOctSize) to avoid
**   compilers' warnings related to returns for all control paths.
**   Finally, modified to get rid of remaining compiler's warnings.
**  Modified:        28-Feb-2007        Antonio Carlos de O. Miranda
**   Modified API to use in Surf3dMesh.
** -----------------------------------------------------------------
**
*/

#include "surf3d_octree.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DRAW_OCTREE 0

#if DRAW_OCTREE
#include <GL/gl.h>
#include <GL/glu.h>
#endif
/* --------------------------------------------------------------
** Definitons for use of algorithms:
*/

#define OCT_REFI        1        /* Refine octree for largest cell   */
#define OCT_MARK        0        /* Mark cells near given boundary   */
#define OCT_CLAS        0        /* Classify all cells of the octree */
#define OCT_GENE        0        /* Generate nodes at unmarked cells */
#define OCT_DISP        0        /* Print octree cells for display   */
#define OCT_FREE        1        /* Release memory for octree cells  */

/* ----------------------------------------------------------------
** Private definitons and data types:
*/

#define REF_FACTOR     0.4       /* changed by ref_factor */
#define MULT_FACTOR    2.0       /* changed by ref_factor */
#define INT_OCTREE     1         /* changed by int_octree */
#define NEW_NODES      10
#define INFPOS         1e+06
#define INFNEG        -1e+06
#define PI             3.141592653589793324
#define MAX(a,b)       (((a)>(b))?(a):(b))
#define MIN(a,b)       (((a)<(b))?(a):(b))

typedef struct _Surf3DOctTreeRec {
  struct _Surf3DOctTreeRec  *child[8] ;
  struct _Surf3DOctTreeRec  *parent ;
  int                       flags ;
  int                       marks ;
  double                    size;
}

Surf3DOctTreeRec, *Surf3DOctTree ;

typedef struct _Box3d {
  double min[3], max[3] ;
}
Box3d ;

/* --------------------------------------------------------------
** Global variables:
*/

int    bound = 0;
double bound_min[3], bound_max[3] ;  /* boundary box */
double bmin[3], bmax[3] ;

/* --------------------------------------------------------------
** Private variables:
*/

static Surf3DOctTree root = NULL;
static Box3d  box ;
static int bound_min_level;

/* --------------------------------------------------------------
** Private functions prototypes:
*/

static Surf3DOctTree Msh3DBdryOct( Surf3DOctTree, Surf3DOctTree, double,
                                  double [3], double [3], double [3],
                                  double, int, int * ) ;
#if 0
static Surf3DOctTree Msh3DBdryOctFull( Surf3DOctTree, Surf3DOctTree, double,
                                     double [3], double [3], double [3],
                                      double, int, int *, int, int ) ;
#endif
static Surf3DOctTree Msh3DRefOct( Surf3DOctTree, Surf3DOctTree, int, int, double,
                                 double [3], double [3]) ;
static void         Msh3DTransOct( Surf3DOctTree, int ) ;
static int          Msh3DOctNeighbor( Surf3DOctTree, int, int, Surf3DOctTree *,
                                      int *, int * ) ;
#if 0
static void         Msh3DOctMark( Surf3DOctTree, double [3], double [3],
                                  double [3] ) ;
#endif
static int          Msh3DOctSize( Surf3DOctTree, double [3], double [3],
                                  double [3], double [3], int, double * ) ;
#if 0
static void         Msh3DOctClass( Surf3DOctTree, int, int, double **, int **,
                                   double [3], double [3] ) ;
static void         Msh3DOctClassFull( Surf3DOctTree, int, int, double **,
                                       int **, double [3], double [3] ) ;
static double       Msh3DOctClassSolidAngle( double **, int, int, int,
                                             double [3] ) ;
static double       Msh3DOctClassSolidEdge( double [3], int, int, double [3],
                                            double [3], double [3], double [3] ) ;
static void         Msh3DOctClassVolum( double [3][3], double [3], double [3],
                                        double *, double * ) ;
static double       Msh3DOctClassVolumFace( double [4][3], int, int [3] ) ;
static void         Msh3DOctGenNodes( Surf3DOctTree, double [3], double [3],
                                      int *, double ** );
static void         Msh3DOctGenNodesFull( Surf3DOctTree, double [3], double [3],
                                          int *, double ** );
#endif
static Surf3DOctTree Msh3DOctAlloc( void ) ;
static void         Msh3DOctFreeAll( void ) ;
#if 0
static void         Msh3DOctFree( Surf3DOctTree ) ;
static void         Msh3DOctDTree( Surf3DOctTree, double [3], double [3] ) ;
static void         Msh3DOctDBound( double [3], double [3] ) ;
static void         Msh3DOctDNodes( int num, double * ) ;
#endif
#if DRAW_OCTREE
static void         Msh3DOctreeDspBody (Surf3DOctTree tree, double min[3], double max[3]);
#endif
static void         Msh3DOctreeVisitLevels (Surf3DOctTree tree, Surf3DOctTree parent, double min[3],
                                            double max[3], void ((*func)(void *octree, void *parent,
                                            int leaf, double min[3], double max[3])));

static Surf3DOctTree Msh3DRefOctToLevel (Surf3DOctTree  *tree, Surf3DOctTree  parent,
                                                      int level);


/*
** ---------------------------------------------------------------
** Public functions:
*/


/* Surf3DOctreeInit          - init the octree with boundbox parameters
************************************************************************/
void* Surf3DOctreeInit (double xmax, double ymax, double zmax,
                        double xmin, double ymin, double zmin)
{
  double  mid[3],  span[3], octmax[3], octmin[3];
  int     i;
  Surf3DOctTree  tree;

  octmax[0] = xmax;
  octmax[1] = ymax;
  octmax[2] = zmax;
  octmin[0] = xmin;
  octmin[1] = ymin;
  octmin[2] = zmin;

  mid[0] = (xmax + xmin) / 2.0;
  mid[1] = (ymax + ymin) / 2.0;
  mid[2] = (zmax + zmin) / 2.0;

  span[0] = xmax - xmin;
  span[1] = ymax - ymin;
  span[2] = zmax - zmin;

  if (( span[0] >= span[1] ) && ( span[0] >= span[2] ))
  {
      octmin[1] = mid[1] - (span[0]/2.0) ;
      octmax[1] = mid[1] + (span[0]/2.0) ;
      octmin[2] = mid[2] - (span[0]/2.0) ;
      octmax[2] = mid[2] + (span[0]/2.0) ;
  }
  else if (( span[1] >= span[0] ) && ( span[1] >= span[2] ))
  {
      octmin[0] = mid[0] - (span[1]/2.0) ;
      octmax[0] = mid[0] + (span[1]/2.0) ;
      octmin[2] = mid[2] - (span[1]/2.0) ;
      octmax[2] = mid[2] + (span[1]/2.0) ;
  }
  else
  {
      octmin[0] = mid[0] - (span[2]/2.0) ;
      octmax[0] = mid[0] + (span[2]/2.0) ;
      octmin[1] = mid[1] - (span[2]/2.0) ;
      octmax[1] = mid[1] + (span[2]/2.0) ;
  }

  /* fill in box structure */

  for ( i=0 ; i<3 ; i++ )
  {
    box.min[i] = octmin[i];
    box.max[i] = octmax[i];
  }

  bound_min_level = 1000000;

  tree = Msh3DOctAlloc() ;
  tree->parent = NULL;
  tree->marks  = 0 ;
  tree->size   = 0.0;

  bound = 0;

  return (void *) tree;
}

/*
************************************************************************/
void Surf3DOctreeBound( double xmax, double ymax, double zmax,
                       double xmin, double ymin, double zmin )
{
  bound = 1;
  bound_max[0] = xmax;
  bound_max[1] = ymax;
  bound_max[2] = zmax;
  bound_min[0] = xmin;
  bound_min[1] = ymin;
  bound_min[2] = zmin;
}
/* Surf3DOctreeAddVertexSize - refine the octree using a size in a specific point
************************************************************************/
int Surf3DOctreeAddPointSize (void *octree, double x, double y, double z, double size)
{
  int min_level = -1;
  double mid[3] = {x, y, z};
  Surf3DOctTree  tree = (Surf3DOctTree) octree;

  if (size <= 0.0)
  {
    printf ("Error! Inserting size = 0.0  in octree!");
    return 0;
  }

  root = Msh3DBdryOct (tree, 0, size, mid, box.max, box.min, 1.0, 1, &min_level ) ;

  if (min_level < bound_min_level)
    bound_min_level = min_level ;

  // check bound box
  if (bound)
  {
    if (x > bound_max[0])
      bound_max[0] = x;
    if (y > bound_max[1])
      bound_max[1] = y;
    if (z > bound_max[2])
      bound_max[2] = z;

    if (x < bound_min[0])
      bound_min[0] = x;
    if (y < bound_min[1])
      bound_min[1] = y;
    if (z < bound_min[2])
      bound_min[2] = z;
  }

  if (root == NULL)
    return 0;

  return 1;
}

/* Surf3DOctreeEnd - refine the whole octree to avoid disparity of cell size
************************************************************************/
void Surf3DOctreeEnd (void *octree, double max_size, int adj_ref)
{
  Surf3DOctTree  tree = (Surf3DOctTree) octree;

#if OCT_REFI
  /* Find the largest cell size that contains a boundary face.
     Refine the entire octree to at least this level. */
  if (max_size > 0)
  {
    int    level = 0;
    double cur_size = box.max[0] - box.min[0];
   /* printf ("\n1 - %d %f %f - ", level, cur_size, max_size);*/
    for (level = 0; cur_size > max_size; level++)
       cur_size *= 0.5;
    /*printf ("%d %f %f - ", level, cur_size, max_size);*/
    root = Msh3DRefOct (tree, 0, 1, level, max_size, box.max, box.min);
    /*printf ("fim\n", level, cur_size, max_size);*/
  }
#endif

  /* Refine the octree so that no adjacent cells differ by a
     refinement factor of more than one. */
  Msh3DTransOct (tree, 0);

}

/* Surf3DOctreeLevelRefine
************************************************************************/
void Surf3DOctreeLevelRefine (void **tree, void *parent, int level)
{
  Surf3DOctTree* tree_c   = (Surf3DOctTree *) tree;
  Surf3DOctTree  parent_c = (Surf3DOctTree) parent;

  Msh3DRefOctToLevel (tree_c, parent_c, level);
}


/* Surf3DOctreeSize - obtain the size of octree in a point
************************************************************************/
double Surf3DOctreeSize (void *octree, double x, double y, double z)
{
  Surf3DOctTree  tree = (Surf3DOctTree) octree;
  double  face_center[3] = {x, y, z};
  double  tree_center[3], size;
  /*int     level;*/

  if (tree == NULL)
    return 0.0;

  /* Find the size of the cell where the given point is */
  /*level = */Msh3DOctSize (tree, box.max, box.min, face_center, tree_center, 1, &size);

  return (size);
}

/* Release memory of internal octree */
/************************************************************************/
void Surf3DOctreeRelease (void *octree)
{
  Msh3DOctFreeAll ( );
}


/* Surf3DOctreeDraw
 ************************************************************************/
void Surf3DOctreeDraw (void *octree)
{
//   Surf3DOctTree  tree = (Surf3DOctTree) octree;
#if DRAW_OCTREE
  Msh3DOctreeDspBody (tree, box.min, box.max);
#endif
}


/* Surf3DOctreeVisitLevels
 ************************************************************************/
void Surf3DOctreeVisitLevels (void *octree, void ((*func)(void *tree, void *parent,
                              int leaf, double min[3], double max[3])))
{
  Surf3DOctTree  tree = (Surf3DOctTree) octree;
  Msh3DOctreeVisitLevels (tree, NULL, box.min, box.max, func);
}

/*
** ---------------------------------------------------------------
** Private functions:
*/


/* -------------------------------------------------------------------
** Msh3DBdryOct - subdivide the oct tree until the area of the sides
**                of the octant containing the given point is about
**                the same as the area associated with the point.
*/

static Surf3DOctTree Msh3DBdryOct
(
Surf3DOctTree  tree,
Surf3DOctTree  parent,
double        size,
double        coord[3],
double        max[3],
double        min[3],
double        factor,
int           this_level,
int           *min_level
)
{
  int           nindx, indx[2], i ;
  double        len, lmin[2][3], lmax[2][3], mid[3] ;

  /* find the size of the cell and the current length */
  len = max[0] - min[0];

  if ( size > len * factor )
  {
    if (( *min_level == -1 ) || ((this_level-1) < *min_level ))
        *min_level = this_level-1 ;

    if( parent != NULL )
       parent->size = size;

    return( tree ) ;
  }

  /* if we get here we decend the tree */
  if ( !tree )
  {
    tree = Msh3DOctAlloc() ;
    tree->parent = parent ;
    tree->size   = 0.0;
  }

  for ( i=0 ; i<3 ; i++ )
    mid[i] = (max[i] + min[i]) / 2.0 ;

  if ( coord[0] == mid[0] || coord[1] == mid[1] || coord[2] == mid[2] )
    nindx = 2 ;
  else
    nindx = 1 ;
  indx[0] = indx[1] = 0 ;

  for ( i=0 ; i<3 ; i++ )
  {
    if ( coord[i] == mid[i] )
    {
      indx[0] |= (1<<i);
      lmin[0][i] = mid[i];
      lmax[0][i] = max[i];
      lmin[1][i] = min[i];
      lmax[1][i] = mid[i];
    }
    else
    {
      if ( coord[i] > mid[i] )
      {
        indx[0] |= (1<<i);
        lmin[0][i] = mid[i];
        lmax[0][i] = max[i];
        indx[1] |= (1<<i);
        lmin[1][i] = mid[i];
        lmax[1][i] = max[i];
      }
      else
      {
        lmin[0][i] = min[i];
        lmax[0][i] = mid[i];
        lmin[1][i] = min[i];
        lmax[1][i] = mid[i];
      }
    }
  }

  for (i = 0; i<nindx; i++)
  {
    tree->child[indx[i]] = Msh3DBdryOct (tree->child[indx[i]], tree, size, coord,
                                         lmax[i], lmin[i], factor, this_level+1, min_level) ;
  }

  return( tree ) ;
}

#if 0
static Surf3DOctTree Msh3DBdryOctFull(
    Surf3DOctTree  tree,
    Surf3DOctTree  parent,
    double        area,
    double        coord[3],
    double        max[3],
    double        min[3],
    double        factor,
    int           this_level,
    int           *min_level,
    int           this_child,
    int           indx_child )
{
    double        len, carea ;
    int           nindx, indx[2], i, this ;
    double        lmin[2][3], lmax[2][3], mid[3] ;

    /* alloc a cell of the tree if it's necessary */

    if ( !tree )
    {
        tree = Msh3DOctAlloc() ;
        tree->parent = parent ;
        tree->size   = 0.0;
    }

    /* return tree if this indx is not the child indx where coord is */

    if( this_child != indx_child )
     return( tree ) ;

    /* find the size of the cell and the current length */

    len = max[0] - min[0] ;  carea = len * len ;

    /* return tree if face's area is bigger than cell's area */

    if ( area > carea * factor ) {
        if (( *min_level == -1 ) || ((this_level-1) < *min_level ))
            *min_level = this_level-1 ;
        return( tree ) ;
    }

    /* if we get here we decend the tree */

    for( this = 0; this < 8; this++ )
    {

     for ( i=0 ; i<3 ; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

     if ( coord[0] == mid[0] || coord[1] == mid[1] || coord[2] == mid[2] )
      nindx = 2 ;
     else
      nindx = 1 ;
     indx[0] = indx[1] = 0 ;

     for ( i=0 ; i<3 ; i++ ) {
         if ( coord[i] == mid[i] ) {
             indx[0] |= (1<<i) ; lmin[0][i] = mid[i] ; lmax[0][i] = max[i] ;
             lmin[1][i] = min[i] ; lmax[1][i] = mid[i] ;
         }
         else {
          if ( coord[i] > mid[i] ) {
              indx[0] |= (1<<i) ;  lmin[0][i] = mid[i] ;  lmax[0][i] = max[i] ;
              indx[1] |= (1<<i) ;  lmin[1][i] = mid[i] ;  lmax[1][i] = max[i] ;
          }
          else {
              lmin[0][i] = min[i] ;  lmax[0][i] = mid[i] ;
              lmin[1][i] = min[i] ;  lmax[1][i] = mid[i] ;
          }
        }
     }

     for ( i=0; i<nindx ; i++ )
      tree->child[this] = Msh3DBdryOctFull( tree->child[this], tree, area,
                             coord, lmax[i], lmin[i],
                             factor, this_level+1, min_level, this, indx[i] ) ;
    }
    return( tree ) ;
}
#endif

/* -------------------------------------------------------------------
** Msh3DRefOct - subdivide the oct tree to at least the specified level.
*/

static Surf3DOctTree Msh3DRefOct(
Surf3DOctTree  tree,
Surf3DOctTree  parent,
int            this_level,
int            min_level,
double         size,
double         max[3],
double         min[3]
)
{
  int j, i;
  double  lmin[3], lmax[3], mid[3];


  /* check boundary */
  if (bound)
  {
    if( max[0] < bound_min[0] || min[0] > bound_max[0] ||
        max[1] < bound_min[1] || min[1] > bound_max[1] ||
        max[2] < bound_min[2] || min[2] > bound_max[2] )
      return( tree ) ;
  }


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
    tree = Msh3DOctAlloc() ;
    tree->parent = parent ;
    tree->size   = 0.0;
  }

  for ( i=0 ; i<3 ; i++ )
    mid[i] = (max[i] + min[i]) / 2.0 ;


  for (j = 0 ;j < 8 ;j++)
  {
    for (i=0 ; i < 3; i++)
    {
      if ( j & (1<<i) )
      {
        lmin[i] = mid[i] ;
        lmax[i] = max[i] ;
      }
      else
      {
        lmin[i] = min[i] ;
        lmax[i] = mid[i] ;
      }
    }

    tree->child[j] = Msh3DRefOct (tree->child[j], tree, this_level+1, min_level, size,
                                  lmax, lmin);
  }
  return( tree ) ;
}

/* -------------------------------------------------------------------
** Msh3DTransOct - subdivide the oct tree so that no two neighbor cells
**                 differ in level by more than one.
*/

static int direct_tbl[8][3] = {
    { 3, 4, 5 },
    { 0, 4, 5 },
    { 3, 1, 5 },
    { 0, 1, 5 },
    { 3, 4, 2 },
    { 0, 4, 2 },
    { 3, 1, 2 },
    { 0, 1, 2 },
    } ;

static void Msh3DTransOct( Surf3DOctTree tree, int level )
{
  int i, j, direction, clevel, nchild ;
  Surf3DOctTree  neighbor ;

  for ( j=0 ; j<8 ; j++ )
  {
    if ( tree->child[j] )
    {
        Msh3DTransOct( tree->child[j], level+1 ) ;
    }
    else
    {
      for ( i=0 ; i<3 ; i++ )
      {

        direction = direct_tbl[j][i] ;

        clevel = level - 1 ;
        Msh3DOctNeighbor (tree, j, direction, &neighbor, &nchild, &clevel ) ;

        /* if the neighbors differ by more than one, subdivide */
        if ( neighbor )
        {
          if ((level-clevel) > 1)
          {
            neighbor->child[nchild] = Msh3DOctAlloc() ;
            neighbor->child[nchild]->parent = neighbor ;
            neighbor->child[nchild]->size = tree->size * MULT_FACTOR;
            Msh3DTransOct( neighbor->child[nchild], clevel+1 ) ;
          }

          //if (neighbor->child[nchild]->size == 0.0)
          //  neighbor->child[nchild]->size = tree->size * MULT_FACTOR;

        }

      }
    }
  }
}

/* -------------------------------------------------------------------
** Msh3DOctNeighbor - given a node in the octree along with a search
**                    direction, this routine finds the neighbor cell
**                    in the oct tree.
**
** Hang on, this algorithm is not for the faint at heart.  It is an
** extenstion of the quadtree algorithm described by Samet in his
** comp surveys review.  Here, however, I play a bunch of bit tricks
** to make the code managable, and to eliminate a bunch of case
** statements.
**
** The basic idea is this, a cells index tells the location of a cell
** in its parent based on the bit pattern in the cell.  Each index
** is comprised of 3 bits, and each bit codes for a particular direction.
** The low order (right most) bit is for x, 2nd lowest for y, and
** third lowest for z.  If we assume a local coordinate system with
** its origin at the center of the parent cell, then a set bit implies
** that the coordinate values for the corresponding direction (x,y,z)
** have a positive value.  Likewise, a cleared bit indicates negative
** values.
**
** When searching for a neighbor, we need to indicate a search
** direction.  This is done with another set of bit patterns which
** allow us to do some simple bit operations.  The patterns are
** as follows:
**     +x -> 000 (1), -x -> 011 (6), +y -> 001 (2), -y -> 100 (5),
**     +z -> 010 (4), -z -> 101 (3)
*/

static int dir_tbl[8][6] = {
    { 0, 0, 0, 1, 1, 1 },
    { 1, 0, 0, 0, 1, 1 },
    { 0, 1, 0, 1, 0, 1 },
    { 1, 1, 0, 0, 0, 1 },
    { 0, 0, 1, 1, 1, 0 },
    { 1, 0, 1, 0, 1, 0 },
    { 0, 1, 1, 1, 0, 0 },
    { 1, 1, 1, 0, 0, 0 },
    } ;

static int mirror_tbl[8][6] = {
    { 1, 2, 4, 1, 2, 4 },
    { 0, 3, 5, 0, 3, 5 },
    { 3, 0, 6, 3, 0, 6 },
    { 2, 1, 7, 2, 1, 7 },
    { 5, 6, 0, 5, 6, 0 },
    { 4, 7, 1, 4, 7, 1 },
    { 7, 4, 2, 7, 4, 2 },
    { 6, 5, 3, 6, 5, 3 },
    } ;

static int Msh3DOctNeighbor(
    Surf3DOctTree  cell,
    int           child,
    int           direction,
    Surf3DOctTree  *cousin,
    int           *cchild,
    int           *level )
{
    Surf3DOctTree  parent ;
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
    for ( indx=0 ; indx<8 ; indx++ )
    {
        if ( parent->child[indx] == cell )
          break ;
    }
    if ( dir_tbl[indx][direction] )
    {
        found = Msh3DOctNeighbor( parent, indx, direction, cousin, cchild, level ) ;
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
** Msh3DOctMark - given the coordinates of a point, this routine finds
**                the terminal cell that contains the point, and marks
**                the appropriate bit in the flag word.
*/

#if 0
static void Msh3DOctMark( Surf3DOctTree tree, double max[3], double min[3],
                          double coord[3] )
{
    int           nindx, indx[2], i ;
    double        lmin[2][3], lmax[2][3], mid[3] ;

    /* find the mid point and the appropriate child cell */

    for ( i=0 ; i<3 ; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

    if ( coord[0] == mid[0] || coord[1] == mid[1] || coord[2] == mid[2] )
     nindx = 2 ;
    else
     nindx = 1 ;
    indx[0] = indx[1] = 0 ;

    for ( i=0 ; i<3 ; i++ )
    {
        if ( coord[i] == mid[i] )
        {
            indx[0] |= (1<<i) ; lmin[0][i] = mid[i] ; lmax[0][i] = max[i] ;
            lmin[1][i] = min[i] ; lmax[1][i] = mid[i] ;
        }
        else
        {
          if ( coord[i] > mid[i] )
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
    }

    /* check to see if the cell has children.  If so, decend the tree.
       If not, set it's flag */

    for ( i=0; i<nindx ; i++ )
    {
      if ( tree->child[indx[i]] )
      {
         Msh3DOctMark( tree->child[indx[i]], lmax[i], lmin[i], coord ) ;
      }
      else
      {
        tree->flags |= (1<<indx[i]) ;
      }
    }
}
#endif

#if 0
static void Msh3DOctMark( tree, max, min, coord )
    Surf3DOctTree  tree ;
    double        coord[3], max[3], min[3] ;
{
    int           indx, i ;
    double        lmin[3], lmax[3], mid[3] ;

    /* find the mid point and the appropriate child cell */

    indx = 0 ;
    for ( i=0 ; i<3 ; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

    for ( i=0 ; i<3 ; i++ ) {
        if ( coord[i] > mid[i] ) {
            indx |= (1<<i) ;  lmin[i] = mid[i] ;  lmax[i] = max[i] ;
        }
        else {
            lmin[i] = min[i] ;  lmax[i] = mid[i] ;
        }
    }

    /* check to see if the cell has children.  If so, decend the tree.
       If not, set it's flag */

    if ( tree->child[indx] ) {
        Msh3DOctMark( tree->child[indx], lmax, lmin, coord ) ;
    }
    else {
        tree->flags |= (1<<indx) ;
    }
}
#endif

/* -------------------------------------------------------------------
** Msh3DOctSize - given the coordinates of a point, this routine finds
**                the terminal cell that contains the point, and returns
**                the size of this cell.
*/

static int Msh3DOctSize(
    Surf3DOctTree  tree,
    double        max[3],
    double        min[3],
    double        face_coord[3],
    double        tree_coord[3],
    int           tree_level,
    double        *size )
{
    int           indx, i ;
    double        lmin[3], lmax[3], mid[3] ;

    /* find the mid point and the appropriate child cell */

    indx = 0 ;
    for ( i=0 ; i<3 ; i++ )
      mid[i] = (max[i] + min[i]) * 0.5 ;

    for ( i=0 ; i<3 ; i++ )
    {
        if ( face_coord[i] > mid[i] )
        {
            indx |= (1<<i) ;
            lmin[i] = mid[i] ;
            lmax[i] = max[i] ;
        }
        else
        {
            lmin[i] = min[i] ;
            lmax[i] = mid[i] ;
        }
    }

    /* check to see if the cell has children.  If so, decend the tree.
       If not, set it's flag */

    if ( tree->child[indx] )
    {
        return ( Msh3DOctSize( tree->child[indx], lmax, lmin, face_coord, tree_coord, tree_level+1, size ) );
    }
    else
    {
        /* tree center and size */
        if ( !(tree->flags & (1<<indx)) )
        {
          /* printf ("Sizes %lf - %lf\n", lmax[0] - lmin[0], tree->size); */
          if (tree->size != 0.0)
            *size = tree->size;
          else
            *size = lmax[0] - lmin[0] ;
          tree_coord[0] = (lmin[0]+lmax[0])/2.0 ;
          tree_coord[1] = (lmin[1]+lmax[1])/2.0 ;
          tree_coord[2] = (lmin[2]+lmax[2])/2.0 ;
        }
        else
        {
          *size = 0.0 ;
          tree_coord[0] = face_coord[0] ;
          tree_coord[1] = face_coord[1] ;
          tree_coord[2] = face_coord[2] ;
        }
        return tree_level ;
    }
}

#if 0
/* -------------------------------------------------------------------
** Msh3DOctClass - this routine classifies the tree. The cells are of
**                 the three types: cells that are in the boundary of
**                 the body, cells that are outside of the body   and
**                 cells inside of the body. This routine also  marks
**                 those ones that are in the boundary or outside  of
**                 the body.
*/

static void Msh3DOctClass(
    Surf3DOctTree  tree,
    int           num_org_nodes,
    int           num_org_faces,
    double        **original_nodes,
    int           **original_faces,
    double        max[3],
    double        min[3] )
{
    int           indx, i, j, k, mark ;
    double        lmin[3], lmax[3], mid[3], cmid[3] ;
    double        cen[3], vec[3], angle, dot ;
    double        coords[3][3], cross[3], pit[3], r[3], s[3] ;
    double        area, length, height, norm, d_pit, d_cen ;

    /* find the mid point of the cell */

    for ( i=0; i<3; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

    /* visit all the child cells */

    for ( indx=0; indx<8; indx++ )
    {
     if ( tree->child[indx] )
     {
         for ( i=0 ; i<3 ; i++ )
         {
             if ( indx & (1<<i) )
             {
                 lmin[i] = mid[i] ;  lmax[i] = max[i] ;
             }
             else {
                 lmin[i] = min[i] ;  lmax[i] = mid[i] ;
             }
        }
        Msh3DOctClass( tree->child[indx], num_org_nodes, num_org_faces,
              original_nodes, original_faces, lmax, lmin ) ;
      }
      else {

        /* initiate mark for this terminal cell */

        mark = 0 ;

        /* find the limits of the terminal cell */

        for( i=0; i<3; i++ ) {
           if ( indx & (1<<i) ) {
               lmin[i] = mid[i] ; lmax[i] = max[i] ;
           }
           else {
               lmin[i] = min[i] ; lmax[i] = mid[i] ;
           }
        }

        /* find the centroid of the terminal cell */

        cen[0] = (lmin[0]+lmax[0])/2.0 ;
        cen[1] = (lmin[1]+lmax[1])/2.0 ;
        cen[2] = (lmin[2]+lmax[2])/2.0 ;

        /* mark the appropiate bit in the flag word of the cell that is
           outside of the body */

        if( !mark ) {

         /* loop trough all boundary faces to increase the contribution to
             the solid angle to this centroid */

         angle = 0.0 ;
         for( i = 0; i < num_org_faces; i++ )
         {
          angle += Msh3DOctClassSolidAngle( original_nodes,original_faces[i][0],
           original_faces[i][1], original_faces[i][2], cen ) ;
         }

         /* mark the cell if the centroid is outside of the body */

         if( !((angle>2*PI)||(angle<-2*PI)) )
          mark = 1 ;

         /* if this cell was marked by the faces, mark the appropiate bit
            in the flag word */

         if( mark )
          tree->flags |= (1<<indx) ;
        }

        /* mark the appropiate bit in the flag word of the cell that is
           in the boundary of the body */

        if( !mark ) {

         /* loop trough all boundary faces to see if the volum of the
            tetrahedron formed by the face and the centroid of the cell is
            lower than the volum of the tetrahedron formed by the the face
            and its optimal point, what represent that the centroid will
            make a very poor tetrahedron if it's choosen by this face */

         for( i = 0; i < num_org_faces; i++ )
         {
          /* find the coords of the face */

          for( j = 0; j < 3; j++ )
           for( k = 0; k < 3; k++ )
            coords[j][k] = original_nodes[original_faces[i][j]][k] ;

          /* find the center, the normal and the area of the face */

          for( j = 0; j < 3; j++ ) {
           cmid[j] = 0.0 ;
           for ( k=0 ; k<3 ; k++ ) cmid[j] += coords[k][j] ;
           cmid[j] /= 3.0 ;
          }
          for ( j=0 ; j<3 ; j++ ) {
           r[j] = coords[1][j] - coords[0][j] ;
           s[j] = coords[2][j] - coords[0][j] ;
          }
          cross[0] = r[1] * s[2] - s[1] * r[2] ;
          cross[1] = r[2] * s[0] - s[2] * r[0] ;
          cross[2] = r[0] * s[1] - s[0] * r[1] ;
          area = sqrt( cross[0] * cross[0] +
                       cross[1] * cross[1] +
                       cross[2] * cross[2] ) / 2.0 ;

          /* find the optimal point (divided by a factor) of the face */

          length = sqrt( (4*area) / sqrt(3.0) ) ;
          height = (length * sqrt(13.0)) / 4.0 ;
          norm = sqrt( cross[0] * cross[0] +
                       cross[1] * cross[1] +
                       cross[2] * cross[2] ) ;
          cross[0] /= norm ;
          cross[1] /= norm ;
          cross[2] /= norm ;
          pit[0] = cmid[0] + (height/NEW_NODES) * cross[0] ;
          pit[1] = cmid[1] + (height/NEW_NODES) * cross[1] ;
          pit[2] = cmid[2] + (height/NEW_NODES) * cross[2] ;

          /* the test shouldn't be done is the cen isn't in the side of
             the face */

          vec[0] = cen[0] - cmid[0] ;
          vec[1] = cen[1] - cmid[1] ;
          vec[2] = cen[2] - cmid[2] ;
          dot = vec[0]*cross[0] + vec[1]*cross[1] + vec[2]*cross[2] ;
          if( dot <= 0.0 ) continue ;

          /* find the volums of the two tetrahedrons */

          Msh3DOctClassVolum( coords, pit, cen, &d_pit, &d_cen );

          /* mark the cell if the volum of the tetrahedron formed by the face
             and the centroid of the cell is less than the volum of the
             tetrahedron formed by the face and its optimal point */

          if( fabs(d_cen) < fabs(d_pit) )
           mark = 1 ;

          /* if this cell was marked by this face, mark the appropiate bit
             in the flag word and break the loop because it's not necessary
             test against other faces anymore */

          if( mark ) {
           tree->flags |= (1<<indx) ;
           break ;
          }
         }
        }
    }
 }
}

static void Msh3DOctClassFull(
    Surf3DOctTree  tree,
    int           num_org_nodes,
    int           num_org_faces,
    double        **original_nodes,
    int           **original_faces,
    double        max[3],
    double        min[3] )
{
    int           indx, i, j, k, found = 0 ;
    double        lmin[3], lmax[3], mid[3] ;
    double        fmin[3], fmax[3], cen[3], angle ;

    /* find the mid point of the cell */

    for ( i=0; i<3; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

    /* verify if the cell has child */

    for ( indx=0; indx<8; indx++ )
     if ( !tree->child[indx] )
      found++ ;

    /* visit all the child cells */

    if ( !found ) {
     for ( indx=0; indx<8; indx++ ) {
      if ( tree->child[indx] ) {
          for ( i=0 ; i<3 ; i++ ) {
              if ( indx & (1<<i) ) {
                      lmin[i] = mid[i] ;  lmax[i] = max[i] ;
              }
              else {
                  lmin[i] = min[i] ;  lmax[i] = mid[i] ;
              }
         }
         Msh3DOctClassFull( tree->child[indx], num_org_nodes, num_org_faces,
               original_nodes, original_faces, lmax, lmin ) ;
       }
      }
    }

    /* if this child is terminal, then classify the cell and mark those
       ones that are in the boundary or outside of the body */

    else {

     /* find the limits of the terminal cell */

     for( i=0; i<3; i++ ) {
       lmax[i] = max[i] ;
       lmin[i] = min[i] ;
     }

     /* find the centroid of the terminal cell */

     cen[0] = (lmin[0]+lmax[0])/2.0 ;
     cen[1] = (lmin[1]+lmax[1])/2.0 ;
     cen[2] = (lmin[2]+lmax[2])/2.0 ;

     /* mark the appropiate bit in the flag word of the cell that is
        in the boundary of the body */

     for( i = 0; i < num_org_faces; i++ )
     {
      /* find the cube related to the face */

      for( j=0; j<3; j++ ) {
       fmin[j] = INFPOS ;
       fmax[j] = INFNEG ;
      }
      for( j=0; j<3; j++ ) {
       for( k=0; k<3; k++ ) {
        fmax[k] = (fmax[k] > original_nodes[original_faces[i][j]][k]) ?
                   fmax[k] : original_nodes[original_faces[i][j]][k] ;
        fmin[k] = (fmin[k] < original_nodes[original_faces[i][j]][k]) ?
                   fmin[k] : original_nodes[original_faces[i][j]][k] ;
       }
      }

      /* mark the cell if it is in the boundary. This is done verifying
         if the cube related to the face intersects the cell */

      if( !( ( ((lmin[0]+(lmax[0]-lmin[0]))<=fmin[0]) ||
               ((lmin[1]+(lmax[1]-lmin[1]))<=fmin[1]) ||
               ((lmin[2]+(lmax[2]-lmin[2]))<=fmin[2])  ) ||
             ( ((fmin[0]+(fmax[0]-fmin[0]))<=lmin[0]) ||
               ((fmin[1]+(fmax[1]-fmin[1]))<=lmin[1]) ||
               ((fmin[2]+(fmax[2]-fmin[2]))<=lmin[2])  )  )  )
       tree->flags = 1 ;
     }

     /* mark the appropiate bit in the flag word of the cell that is
        outside of the body */

     angle = 0.0 ;
     for( i = 0; i < num_org_faces; i++ )
     {
      /* loop trough all boundary faces to increase the contribution to
         the solid angle to the centroid of the cell */

      angle += Msh3DOctClassSolidAngle( original_nodes,original_faces[i][0],
       original_faces[i][1], original_faces[i][2], cen ) ;
     }

     /* mark the cell if the centroid is outside of the body */

     if( !((angle>2*PI)||(angle<-2*PI)) )
      tree->flags = 1 ;
    }
}

static double Msh3DOctClassSolidAngle(
    double             **original_nodes,
    int                v0,
    int                       v1,
    int                v2,
    double             cen[3] )
{
    double              *c0, *c1, *c2 ;
    double              area = 0.0, mod, dplane, dot ;
    double              r[3], s[3], r_x_s[3], a[3], r1[3], b[3] ;

    c0 = original_nodes[v0] ;
    c1 = original_nodes[v1] ;
    c2 = original_nodes[v2] ;

    r[0] = c1[0]-c0[0] ; r[1] = c1[1]-c0[1] ; r[2] = c1[2]-c0[2] ;
    s[0] = c2[0]-c1[0] ; s[1] = c2[1]-c1[1] ; s[2] = c2[2]-c1[2] ;
    r_x_s[0] = r[1] * s[2] - s[1] * r[2] ;
    r_x_s[1] = r[2] * s[0] - s[2] * r[0] ;
    r_x_s[2] = r[0] * s[1] - s[0] * r[1] ;
    mod = sqrt( r_x_s[0]*r_x_s[0] + r_x_s[1]*r_x_s[1] + r_x_s[2]*r_x_s[2] ) ;
    r_x_s[0] /= mod ;
    r_x_s[1] /= mod ;
    r_x_s[2] /= mod ;
    dplane = -((r_x_s[0]*c0[0])+(r_x_s[1]*c0[1])+(r_x_s[2]*c0[2])) ;

    a[0] = c2[0]-c0[0] ; a[1] = c2[1]-c0[1] ; a[2] = c2[2]-c0[2] ;
    r1[0] = cen[0]-c0[0] ; r1[1] = cen[1]-c0[1] ; r1[2] = cen[2]-c0[2] ;
    b[0] = c1[0]-c0[0] ; b[1] = c1[1]-c0[1] ; b[2] = c1[2]-c0[2] ;
    area += Msh3DOctClassSolidEdge( cen, v0, v1, a, b, r1, r_x_s ) ;

    a[0] = -1.0*b[0] ; a[1] = -1.0*b[1] ; a[2] = -1.0*b[2] ;
    r1[0] = cen[0]-c1[0] ; r1[1] = cen[1]-c1[1] ; r1[2] = cen[2]-c1[2] ;
    b[0] = c2[0]-c1[0] ; b[1] = c2[1]-c1[1] ; b[2] = c2[2]-c1[2] ;
    area += Msh3DOctClassSolidEdge( cen, v1, v2, a, b, r1, r_x_s ) ;

    a[0] = -1.0*b[0] ; a[1] = -1.0*b[1] ; a[2] = -1.0*b[2] ;
    r1[0] = cen[0]-c2[0] ; r1[1] = cen[1]-c2[1] ; r1[2] = cen[2]-c2[2] ;
    b[0] = c0[0]-c2[0] ; b[1] = c0[1]-c2[1] ; b[2] = c0[2]-c2[2] ;
    area += Msh3DOctClassSolidEdge( cen, v2, v0, a, b, r1, r_x_s ) ;

    area -= PI ;
    dot = ( r_x_s[0] * r1[0] + r_x_s[1] * r1[1] + r_x_s[2] * r1[2] ) ;
    return (dot>0.0) ? (-area) : (area) ;
}

static double Msh3DOctClassSolidEdge(
    double             cen[3],
    int                vi,
    int                vj,
    double             a[3],
    double             b[3],
    double             r1[3],
    double             r_x_s[3] )
{
    double             n1[3], n2[3], l1, l2, s, ang, ba[3] ;

    n1[0] = a[1] * r1[2] - r1[1] * a[2] ;
    n1[1] = a[2] * r1[0] - r1[2] * a[0] ;
    n1[2] = a[0] * r1[1] - r1[0] * a[1] ;
    n2[0] = r1[1] * b[2] - b[1] * r1[2] ;
    n2[1] = r1[2] * b[0] - b[2] * r1[0] ;
    n2[2] = r1[0] * b[1] - b[0] * r1[1] ;
    l1 = sqrt( n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2] ) ;
    l2 = sqrt( n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2] ) ;
    s = (( n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2] ) / ( l1 * l2 )) ;
    ang = acos(MAX(-1.0,MIN(1.0,s))) ;
    ba[0] = b[1] * a[2] - a[1] * b[2] ;
    ba[1] = b[2] * a[0] - a[2] * b[0] ;
    ba[2] = b[0] * a[1] - a[0] * b[1] ;
    s = ( ba[0] * r_x_s[0] + ba[1] * r_x_s[1] + ba[2] * r_x_s[2] ) ;
    return (s>0.0) ? (PI-ang) : (PI+ang) ;
}

static double Msh3DOctClassVolumFace(
    double             coord_v[4][3],
    int                base_v,
    int                v[3] )
{
    int                i ;
    double             vol, area[3], cross[3] ;

    /* compute the area for the face */

    area[0] = area[1] = area[2] = 0.0 ;
    for( i = 0; i < 3; i++ ) {

     cross[0] = coord_v[v[i]][1] * coord_v[v[(i+1)%3]][2] -
                coord_v[v[i]][2] * coord_v[v[(i+1)%3]][1] ;
     cross[1] = coord_v[v[i]][2] * coord_v[v[(i+1)%3]][0] -
                coord_v[v[i]][0] * coord_v[v[(i+1)%3]][2] ;
     cross[2] = coord_v[v[i]][0] * coord_v[v[(i+1)%3]][1] -
                coord_v[v[i]][1] * coord_v[v[(i+1)%3]][0] ;

     area[0] += ((1/6.) * cross[0]) ;
     area[1] += ((1/6.) * cross[1]) ;
     area[2] += ((1/6.) * cross[2]) ;
    }

    /* compute the volume contribution for the face */

    vol = area[0] * coord_v[base_v][0] +
          area[1] * coord_v[base_v][1] +
          area[2] * coord_v[base_v][2] ;

    return vol ;
}
#endif

#if 0
static void Msh3DOctClassVolum(
    double             coords[3][3],
    double             pit[3],
    double             cen[3],
    double             *d_pit,
    double             *d_cen )
{
    int                i, j, verts[3] ;
    double             coord_verts[4][3] ;

    /* initialize volum for pit and cen */

    (*d_pit) = (*d_cen) = 0.0 ;

    /* fill vertexs */

    for( i = 0; i < 3; i++ )
     for( j = 0; j < 3; j++ )
      coord_verts[i][j] = coords[i][j] ;

    /* set pit as point out the plane face in vertexs */

    for( j = 0; j < 3; j++ )
     coord_verts[3][j] = pit[j] ;

    /* compute volum for face plus pit */

    verts[0] = 0 ;
    verts[1] = 1 ;
    verts[2] = 2 ;
    (*d_pit) += Msh3DOctClassVolumFace( coord_verts, verts[0], verts ) ;
    verts[0] = 1 ;
    verts[1] = 0 ;
    verts[2] = 3 ;
    (*d_pit) += Msh3DOctClassVolumFace( coord_verts, verts[0], verts ) ;
    verts[0] = 2 ;
    verts[1] = 1 ;
    verts[2] = 3 ;
    (*d_pit) += Msh3DOctClassVolumFace( coord_verts, verts[0], verts ) ;
    verts[0] = 0 ;
    verts[1] = 2 ;
    verts[2] = 3 ;
    (*d_pit) += Msh3DOctClassVolumFace( coord_verts, verts[0], verts ) ;

    /* set cen as point out the plane face in vertexs */

    for( j = 0; j < 3; j++ )
     coord_verts[3][j] = cen[j] ;

    /* compute volum for face plus pit */

    verts[0] = 0 ;
    verts[1] = 1 ;
    verts[2] = 2 ;
    (*d_cen) += Msh3DOctClassVolumFace( coord_verts, verts[0], verts ) ;
    verts[0] = 1 ;
    verts[1] = 0 ;
    verts[2] = 3 ;
    (*d_cen) += Msh3DOctClassVolumFace( coord_verts, verts[0], verts ) ;
    verts[0] = 2 ;
    verts[1] = 1 ;
    verts[2] = 3 ;
    (*d_cen) += Msh3DOctClassVolumFace( coord_verts, verts[0], verts ) ;
    verts[0] = 0 ;
    verts[1] = 2 ;
    verts[2] = 3 ;
    (*d_cen) += Msh3DOctClassVolumFace( coord_verts, verts[0], verts ) ;
}

/* -------------------------------------------------------------------
** Msh3DOctGenNodes - generate nodes at the center of all unmarked cells.
*/

#define MSH3D_NODE_QUANTUM 100
static int nodes_alloced = 0 ;

static void Msh3DOctGenNodes(
    Surf3DOctTree  tree,
    double        max[3],
    double        min[3],
    int           *num_nodes,
    double        **nodes )
{
    int           indx, i ;
    double        lmin[3], lmax[3], mid[3] ;

    /* find the mid point of the cell */

    for ( i=0 ; i<3 ; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

    /* visit all the child cells */

    for ( indx=0 ; indx<8 ; indx++ ) {
        if ( tree->child[indx] ) {
            for ( i=0 ; i<3 ; i++ ) {
                if ( indx & (1<<i) ) {
                    lmin[i] = mid[i] ;  lmax[i] = max[i] ;
                }
                else {
                    lmin[i] = min[i] ;  lmax[i] = mid[i] ;
                }
            }
            Msh3DOctGenNodes( tree->child[indx], lmax, lmin,
                              num_nodes, nodes ) ;
        }

        /* if this child is terminal, then check for a mark
           and generate a node */

        else {
            if ( !(tree->flags & (1<<indx)) ) {

                /* find the limits of the terminal cell */

                for ( i=0 ; i<3 ; i++ ) {
                    if ( indx & (1<<i) ) {
                        lmin[i] = mid[i] ;  lmax[i] = max[i] ;
                    }
                    else {
                        lmin[i] = min[i] ;  lmax[i] = mid[i] ;
                    }
                }

                /* create the node if the cell isn't out of the boundary */

                if ( !( (lmin[0] >= bmax[0]) ||
                        (lmin[1] >= bmax[1]) ||
                        (lmin[2] >= bmax[2]) ||
                        (lmax[0] <= bmin[0]) ||
                        (lmax[1] <= bmin[1]) ||
                        (lmax[2] <= bmin[2]) ) ) {

                /* allocate space if we need it */

                if ( *num_nodes >= nodes_alloced ) {
                    if ( !nodes_alloced ) {
                        nodes_alloced = MSH3D_NODE_QUANTUM ;
                        *nodes = (double *)malloc(
                            nodes_alloced * 3 * sizeof(double)) ;
                    }
                    else {
                        nodes_alloced += MSH3D_NODE_QUANTUM ;
                        *nodes = (double *)realloc( *nodes,
                            nodes_alloced * 3 * sizeof(double)) ;
                    }
                }

                /* create the node */


                for ( i=0 ; i<3 ; i++ )
                   (*nodes)[(*num_nodes)*3+i] = (lmin[i]+lmax[i])/2.0 ;
                *num_nodes += 1 ;
               }
            }
        }
    }
}

static void Msh3DOctGenNodesFull(
    Surf3DOctTree  tree,
    double        max[3],
    double        min[3],
    int           *num_nodes,
    double        **nodes )
{
    int           indx, i, found = 0 ;
    double        lmin[3], lmax[3], mid[3] ;

    /* find the mid point of the cell */

    for ( i=0 ; i<3 ; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

    /* verify if cell has child */

    for ( indx=0; indx<8; indx++ )
     if ( !tree->child[indx] )
      found++ ;

    /* visit all the child cells */

    if ( !found ) {
     for ( indx=0 ; indx<8 ; indx++ ) {
        if ( tree->child[indx] ) {
            for ( i=0 ; i<3 ; i++ ) {
                if ( indx & (1<<i) ) {
                    lmin[i] = mid[i] ;  lmax[i] = max[i] ;
                }
                else {
                    lmin[i] = min[i] ;  lmax[i] = mid[i] ;
                }
            }
            Msh3DOctGenNodesFull( tree->child[indx], lmax, lmin,
                              num_nodes, nodes ) ;
        }
      }
    }

    /* if this child is terminal, then check for a mark
       and generate a node */

    else {
      if ( !tree->marks /* !tree->flags */ ) {

      /* allocate space if we need it */

      if ( *num_nodes >= nodes_alloced ) {
        if ( !nodes_alloced ) {
          nodes_alloced = MSH3D_NODE_QUANTUM ;
          *nodes = (double *)malloc(
                 nodes_alloced * 3 * sizeof(double)) ;
         }
         else {
          nodes_alloced += MSH3D_NODE_QUANTUM ;
          *nodes = (double *)realloc( *nodes,
                 nodes_alloced * 3 * sizeof(double)) ;
         }
       }

       /* create the node */

       for ( i=0 ; i<3 ; i++ ) {
        lmax[i] = max[i] ;
        lmin[i] = min[i] ;

        (*nodes)[(*num_nodes)*3+i] = (lmin[i]+lmax[i])/2.0 ;
       }
       *num_nodes += 1 ;
      }
    }
}

/* -------------------------------------------------------------------
** Msh3DOctDBound - display the octreee boundary.
*/

static void Msh3DOctDBound( double min[3], double max[3] )
{
#if OCT_DISP
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], min[1], min[2], max[0], min[1], min[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], min[1], max[2], max[0], min[1], max[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], max[1], max[2], max[0], max[1], max[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], max[1], min[2], max[0], max[1], min[2] ) ;

    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], min[1], min[2], min[0], min[1], max[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            max[0], min[1], min[2], max[0], min[1], max[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            max[0], max[1], min[2], max[0], max[1], max[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], max[1], min[2], min[0], max[1], max[2] ) ;

    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], min[1], min[2], min[0], max[1], min[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            max[0], min[1], min[2], max[0], max[1], min[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            max[0], min[1], max[2], max[0], max[1], max[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], min[1], max[2], min[0], max[1], max[2] ) ;
#else
    return ;
#endif
}
#endif

#if 0
/* -------------------------------------------------------------------
** Msh3DOctDTree - display the octreee innards.
*/

static void Msh3DOctDTree( Surf3DOctTree tree, double min[3], double max[3] )
{
    int     i, j ;
    double  lmin[3], lmax[3], mid[3] ;

#if OCT_DISP
    if ( !tree ) return ;
    for ( i=0 ; i<3 ; i++ ) mid[i] = (max[i] + min[i]) / 2.0 ;

    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], mid[1], mid[2], max[0], mid[1], mid[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            mid[0], min[1], mid[2], mid[0], max[1], mid[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            mid[0], mid[1], min[2], mid[0], mid[1], max[2] ) ;

    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], mid[1], min[2], min[0], mid[1], max[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], min[1], mid[2], min[0], max[1], mid[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            max[0], mid[1], min[2], max[0], mid[1], max[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            max[0], min[1], mid[2], max[0], max[1], mid[2] ) ;

    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], min[1], mid[2], max[0], min[1], mid[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            mid[0], min[1], min[2], mid[0], min[1], max[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], max[1], mid[2], max[0], max[1], mid[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            mid[0], max[1], min[2], mid[0], max[1], max[2] ) ;

    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], mid[1], min[2], max[0], mid[1], min[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            mid[0], min[1], min[2], mid[0], max[1], min[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            min[0], mid[1], max[2], max[0], mid[1], max[2] ) ;
    printf( "define line = {{ %f %f %f } { %f %f %f }}\n",
            mid[0], min[1], max[2], mid[0], max[1], max[2] ) ;

    for ( j=0 ; j<8 ; j++ ) {
        if ( tree->child[j] ) {
            for ( i=0 ; i<3 ; i++ ) {
                if ( j & (1<<i) ) {
                    lmin[i] = mid[i] ;  lmax[i] = max[i] ;
                }
                else {
                    lmin[i] = min[i] ;  lmax[i] = mid[i] ;
                }
            }
            Msh3DOctDTree( tree->child[j], lmin, lmax ) ;
        }
    }
#else
    for (i = 0; i < 3; i++) {
     lmin[i] = 0 ;
     lmax[i] = 0 ;
    }
    for (j = 0; j < 3; j++) {
     mid[j] = 0 ;
    }
    return ;
#endif
}
#endif

/* -------------------------------------------------------------------
** Msh3DOctDNode - display the internal nodes.
*/

#if 0
static void Msh3DOctDNodes( int num, double *nodes )
{
    int     i ;

#if OCT_DISP
    for ( i=0 ; i<num ; i++ ) {
        printf( "define point = { %f %f %f }\n",
            nodes[i*3], nodes[i*3+1], nodes[i*3+2] ) ;
    }
#else
    i = 0 ;
    return ;
#endif
}
#endif

/* -------------------------------------------------------------------
** Msh3DOctAlloc - these routines manage the allocation and freeing
**                 oct tree entries.
*/

#define MSH3D_OCT_TREE_QUANTUM 1000
static Surf3DOctTree oct_free = 0 ;
static Surf3DOctTree oct_block_ptr = 0 ;

static Surf3DOctTree Msh3DOctAlloc(void)
{
    Surf3DOctTree  new_block, alloced ;
    int           i ;

    /* if the free pointer is null we need to allocate a new block
       of boundary nodes */

    if ( !oct_free )
    {
        new_block = (Surf3DOctTree) malloc(MSH3D_OCT_TREE_QUANTUM * sizeof(Surf3DOctTreeRec) ) ;
        new_block[0].child[0] = oct_block_ptr ;
        oct_block_ptr = new_block ;
        for ( i=1 ; i<(MSH3D_OCT_TREE_QUANTUM-1) ; i++ ) {
            new_block[i].child[0] = &(new_block[i+1]) ;
        }
        new_block[MSH3D_OCT_TREE_QUANTUM-1].child[0] = 0 ;
        oct_free = &(new_block[1]) ;
    }

    /* return the top thing on the free list */

    alloced = oct_free ;
    oct_free = oct_free->child[0] ;

    for ( i=0 ; i<8 ; i++ ) alloced->child[i] = 0 ;
    alloced->parent = 0 ;
    alloced->flags = 0 ;
    return( alloced ) ;
}

#if 0
static void Msh3DOctFree( Surf3DOctTree oct )
{
    /* Put this face back on the free list */

#if OCT_FREE
    if (oct) oct->child[0] = oct_free ;
    oct_free = oct ;
#else
    return ;
#endif
}
#endif

static void Msh3DOctFreeAll(void)
{
    Surf3DOctTree  cur, next ;

    /* free all blocks allocated to store face information */

#if OCT_FREE
    if ( oct_block_ptr ) cur = oct_block_ptr ;
    else return ;

    while ( cur->child[0] ) {
        next = cur->child[0] ;
        free( cur ) ;
        cur = next ;
    }
    free( cur ) ;

    oct_free = 0 ;
    oct_block_ptr = 0 ;
#else
    cur = next = 0 ;
    return ;
#endif
}

#if DRAW_OCTREE
static void Msh3DOctreeDspBody (Surf3DOctTree tree, double min[3], double max[3] )
{
#if DRAW_OCTREE
  int     i, j ;
  double   lmin[3], lmax[3], mid[3] ;

  for ( i=0 ; i<3 ; i++ )
    mid[i] = (max[i] + min[i]) / 2.0 ;

  for ( j=0 ; j<8 ; j++ )
  {
    if ( tree->child[j] )
    {
      for ( i=0 ; i<3 ; i++ )
      {
        if ( j & (1<<i) )
        {
          lmin[i] = mid[i] ;
          lmax[i] = max[i] ;
        }
        else
        {
          lmin[i] = min[i] ;
          lmax[i] = mid[i] ;
        }
      }

      Msh3DOctreeDspBody ( tree->child[j], lmin, lmax ) ;
    }
    else
    {
      for( i = 0; i < 3; i++ )
      {
        if( j & (1<<i) )
        {
	         lmin[i] = mid[i] ;
           lmax[i] = max[i] ;
        }
        else
        {
	         lmin[i] = min[i] ;
           lmax[i] = mid[i] ;
        }
      }

      /* First face of the tree node */
      glBegin( GL_LINES );
      glVertex3d( lmin[0], lmin[1], lmin[2] );
      glVertex3d( lmax[0], lmin[1], lmin[2] );
      glVertex3d( lmax[0], lmin[1], lmax[2] );
      glVertex3d( lmin[0], lmin[1], lmax[2] );
      glEnd( );

      /* Second face of the tree node */
      glBegin( GL_LINES );
      glVertex3d( lmin[0], lmax[1], lmax[2] );
      glVertex3d( lmax[0], lmax[1], lmax[2] );
      glVertex3d( lmax[0], lmax[1], lmin[2] );
      glVertex3d( lmin[0], lmax[1], lmin[2] );
      glEnd( );

      /* Third face of the tree node */
      glBegin( GL_LINES );
      glVertex3d( lmin[0], lmin[1], lmax[2] );
      glVertex3d( lmax[0], lmin[1], lmax[2] );
      glVertex3d( lmax[0], lmax[1], lmax[2] );
      glVertex3d( lmin[0], lmax[1], lmax[2] );
      glEnd( );

      /* Fourth face of the tree node */
      glBegin( GL_LINES );
      glVertex3d( lmin[0], lmax[1], lmin[2] );
      glVertex3d( lmax[0], lmax[1], lmin[2] );
      glVertex3d( lmax[0], lmin[1], lmin[2] );
      glVertex3d( lmin[0], lmin[1], lmin[2] );
      glEnd( );

      /* Fifth face of the tree node */
      glBegin( GL_LINES );
      glVertex3d( lmax[0], lmin[1], lmin[2] );
      glVertex3d( lmax[0], lmax[1], lmin[2] );
      glVertex3d( lmax[0], lmax[1], lmax[2] );
      glVertex3d( lmax[0], lmin[1], lmax[2] );
      glEnd( );

      /* Sixty face of the tree node */
      glBegin( GL_LINES );
      glVertex3d( lmin[0], lmax[1], lmax[2] );
      glVertex3d( lmin[0], lmax[1], lmin[2] );
      glVertex3d( lmin[0], lmin[1], lmin[2] );
      glVertex3d( lmin[0], lmin[1], lmax[2] );
      glEnd( );
    }
  }
#endif
}
#endif

static void Msh3DOctreeVisitLevels (Surf3DOctTree tree, Surf3DOctTree parent, double min[3],
                                    double max[3], void ((*func)(void *octree, void *parent,
                                    int leaf, double min[3], double max[3])))
{
  int     i, j ;
  double   lmin[3], lmax[3], mid[3] ;

  for (i = 0; i < 3 ; i++)
    mid[i] = (max[i] + min[i]) / 2.0;

  for (j = 0 ;j < 8 ;j++)
  {
    if (tree->child[j] != NULL)
    {
      for (i=0 ; i < 3; i++)
      {
        if ( j & (1<<i) )
        {
          lmin[i] = mid[i] ;
          lmax[i] = max[i] ;
        }
        else
        {
          lmin[i] = min[i] ;
          lmax[i] = mid[i] ;
        }
      }
      if (func != NULL)
        (func) (&(tree->child[j]), tree, 0, lmin, lmax);
      Msh3DOctreeVisitLevels (tree->child[j], tree, lmin, lmax, func);
    }
    else
    {
      for (i = 0; i < 3; i++)
      {
        if (j & (1<<i))
        {
	         lmin[i] = mid[i] ;
           lmax[i] = max[i] ;
        }
        else
        {
	         lmin[i] = min[i] ;
           lmax[i] = mid[i] ;
        }
      }
      if (func != NULL)
        (func) (&(tree->child[j]), tree, 1, lmin, lmax);
    }
  }
}


static Surf3DOctTree Msh3DRefOctToLevel (
    Surf3DOctTree  *tree,
    Surf3DOctTree  parent,
    int            level)
{
    int j ;

    if (level == 0)
      return (*tree);

    if (*tree == NULL)
    {
      *tree = Msh3DOctAlloc();
      (*tree)->parent = parent;
    }

    for (j = 0 ; j < 8; j++)
      (*tree)->child[j] = Msh3DRefOctToLevel (&((*tree)->child[j]), *tree, level-1) ;

    return (*tree) ;
}

