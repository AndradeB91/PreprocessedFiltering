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


#define  QUAD_DRAW  0 
#define  QUAD_REF   1  

#if QUAD_DRAW
#include "iup.h"
#include "cd.h"
#include "cdiup.h"
#include "wd.h"
#endif


#include "msh_def2d.h" 
#include "msh2d.h"


/* --------------------------------------------------------------
** Private definitons and data types:
*/

#define REF_FACTOR     1.0       /* changed by ref_factor */
#define NEW_NODES      10 
#define INFPOS         1e+06
#define PI             acos(-1.0)
#define MAX(a,b)       (((a)>(b))?(a):(b))
#define MIN(a,b)       (((a)<(b))?(a):(b))

typedef struct _Msh2DQuadTreeRec {
  struct _Msh2DQuadTreeRec  *child[4] ;
  struct _Msh2DQuadTreeRec  *parent ;
  int                        flags ;
  int                        marks ;
}
Msh2DQuadTreeRec, *Msh2DQuadTree ;

typedef struct _Box2d {
  double min[2], max[2] ;
}
Box2d ;

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


/* --------------------------------------------------------------
** Private functions prototypes: 
*/

static Msh2DQuadTree Msh2DBdryQuad( Msh2DQuadTree, Msh2DQuadTree, double,
                     double *, double *, double *, double, int, int *) ;
static void          Msh2DTransQuad( Msh2DQuadTree, int ) ;
static int           Msh2DQuadNeighbor( Msh2DQuadTree, int, int,
                     Msh2DQuadTree *, int *, int *) ;
static int           Msh2DQuadSize( Msh2DQuadTree, double *, double *,
                     double *, double *, int, double *) ;
static Msh2DQuadTree Msh2DQuadAlloc() ;

#if QUAD_DRAW
 static void         Msh2DQuadDraw( Msh2DQuadTree, int, int,
                     double original_nodes[][2], int original_edges[][2],
                     double *, double *) ;
#endif

#ifdef QUAD_REF
 static Msh2DQuadTree Msh2DRefQuad( Msh2DQuadTree, Msh2DQuadTree, int, int) ;
#endif


/* -------------------------------------------------------------------
** Msh2DInternalNodes - find the location of internal nodes.
*/

int Msh2DGenQuadTree
(
int     num_org_nodes, 
int     num_org_edges, 
double  original_nodes[][2],
int     original_edges[][2],
int     *num_gen_nodes,
double  **generated_nodes
)
{
    double  max[2], min[2], mid[2], span[2] ;
    double  coords[2][2];
    double  length ;
    int     i, j, k, min_level;
    int     bound_min_level = (int) INFPOS ;
    double  ref_factor; /* = REF_FACTOR;     */
    Msh2DQuadTree tree = 0;

    *num_gen_nodes =0;
    
    ref_factor = Msh2DGetRefFactor ();
    // printf ("Fator de refinamento = %f\n", ref_factor);
	
	
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

        for ( j=0 ; j<2 ; j++ )
        {
            mid[j] = 0.0 ;
            for ( k=0 ; k<2 ; k++ ) mid[j] += coords[k][j] ;
            mid[j] /= 2.0 ;
        }

		
        /* compute the length  */

        length = sqrt( 
          ((coords[1][0]-coords[0][0]) * (coords[1][0]-coords[0][0])) +
          ((coords[1][1]-coords[0][1]) * (coords[1][1]-coords[0][1])) ); 

        min_level = -1 ;

        if (length > 0.)
        {
          tree = Msh2DBdryQuad( tree, 0, length, mid, max, min, ref_factor,
            1, &min_level ) ;
          if( min_level < bound_min_level ) bound_min_level = min_level ;
        }
    }
    
    if ( !tree )
    {
        tree = Msh2DQuadAlloc() ;
        tree->parent = 0 ;
        tree->marks  = 0 ;
    }

#if QUAD_DRAW 
    Msh2DQuadDraw( tree, num_org_nodes, num_org_edges, original_nodes,
                   original_edges, max, min );
    IupMessage("Pause","QTree");
#endif


#ifdef QUAD_REF 
    /* 1.3 Find the largest cell size that contains a boundary edge.
           Refine the entire quadtree to at least this level. */
   tree = Msh2DRefQuad( tree, 0, 1, bound_min_level ) ;

#if  QUAD_DRAW 
    Msh2DQuadDraw( tree, num_org_nodes, num_org_edges, original_nodes,
                   original_edges, max, min );
    IupMessage("Pause","QTree");
#endif
  
#endif  

    /* 1.4 Refine the quadtree so that no adjacent cells differ by a
           refinement factor of more than one. */

    Msh2DTransQuad( tree, 0 ) ;


#if QUAD_DRAW
    Msh2DQuadDraw( tree, num_org_nodes, num_org_edges, original_nodes,
		   original_edges, max, min );
    IupMessage("Pause","QTree"); 
#endif


    root = tree ;

    return(1) ;
}

/* -------------------------------------------------------------------
** Msh2DOptimalNodes - find the size of the cell where the given edge 
**                     center is to use for the location of optimal nodes.
*/

double Msh2DOptimalNodes
(
double  edge_center[2], 
double  tree_center[2], 
int   *tree_level
) 
{
    int    level ;
    double  size ;

    if( root == NULL ) return 0.0 ;

    /* 1.1 Find the size of the cell where the given point is */
    
    level = Msh2DQuadSize( root, gmax, gmin, edge_center, tree_center, 1, &size ) ;
    *tree_level = level ;

    /* 1.2 Return the found size */

    return size ;
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
	
       return( tree ) ;
    }

    /* if we get here we decend the tree */

    if ( !tree )
    {
        tree = Msh2DQuadAlloc() ;
        tree->parent = parent ;
        tree->marks  = 0 ;
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

    for ( i=0; i<nindx ; i++ )
     tree->child[indx[i]] = Msh2DBdryQuad( tree->child[indx[i]], tree, length, 
                             coord, lmax[i], lmin[i],
                             factor, this_level+1, min_level ) ;

    return( tree ) ;
}


#ifdef QUAD_REF
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
    if ( !tree )
    {
        tree = Msh2DQuadAlloc() ;
        tree->parent = parent ;
    }

    for ( j=0 ; j<4 ; j++ )
    {
        tree->child[j] =
            Msh2DRefQuad( tree->child[j], tree, this_level+1, min_level ) ;
    }
    return( tree ) ;
}
#endif




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
double        *size
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

    /* check to see if the cell has children.  If so, decend the tree.
       If not, set it's flag */

    if ( tree->child[indx] ) 
    {
        Msh2DQuadSize( tree->child[indx], lmax, lmin, edge_coord, 
                       tree_coord, tree_level+1, size ) ;
    }
    else 
    {
        /* tree center and size */

        if ( !(tree->flags & (1<<indx)) ) 
        {
         *size = lmax[0] - lmin[0] ;
         tree_coord[0] = (lmin[0]+lmax[0])/2.0 ;
         tree_coord[1] = (lmin[1]+lmax[1])/2.0 ;
        } 
        else 
        {
         *size = 0.0 ;
         tree_coord[0] = edge_coord[0] ;
         tree_coord[1] = edge_coord[1] ;
        }
        return tree_level ;
    }

    return -1;
}


#if QUAD_DRAW

/* -------------------------------------------------------------------
** Msh2DQuadDraw - 
*/

static void Msh2DQuadDraw
(
Msh2DQuadTree tree,
int           num_org_nodes,
int           num_org_edges, 
double        original_nodes[][2],
int           original_edges[][2],
double        max[2],
double        min[2]
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
          Msh2DQuadDraw( tree->child[indx], num_org_nodes, num_org_edges, 
	                         original_nodes, original_edges, lmax, lmin ) ;
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
        
        cdForeground(CD_YELLOW);
	wdLine (lmin[0],lmin[1],lmin[0],lmax[1]);
        wdLine (lmin[0],lmax[1],lmax[0],lmax[1]);
	wdLine (lmax[0],lmax[1],lmax[0],lmin[1]);
        wdLine (lmax[0],lmin[1],lmin[0],lmin[1]); 

    }
 }
}

#endif


/* -------------------------------------------------------------------
** Msh2DQuadAlloc - these routines manage the allocation and freeing
**                 quad tree entries.
*/

#define MSH2D_QUAD_TREE_QUANTUM 100
static Msh2DQuadTree quad_free = 0 ;
static Msh2DQuadTree quad_block_ptr = 0 ;

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
    return( alloced ) ;
}

#if 0 /* NAO USADA */
static void Msh2DQuadFree( Msh2DQuadTree quad )
{
    /* Put this edge back on the free list */

    quad->child[0] = quad_free ;
    quad_free = quad ;
}
#endif

