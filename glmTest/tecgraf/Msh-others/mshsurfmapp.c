/*
**	mshmap2d - package for creating meshes based on transfinite mapping.
**		   It will be used the following meshing techniques:
**
**		 - MshSurfTransfin:
**		     Bilinear mapping for quadrilateral regions with sides 
**		     of any shape, i.e., the boundary may be formed by 
**		     lines, circles, beziers, and so on.  The number of 
**		     segments on opposite sides must be equal.
**		 - MshSurfLofting:
**		     Linear mapping between two opposite sides of a 
**		     quadrilateral region.  The other two sides must 
**		     be straight lines.  The number of segments on opposite
**		     sides must be equal.
**		 - MshSurfTrsfncoll
**		     Bilinear mapping for triangular regions which conceptually
**		     is formed by a quadrilateral region with one of the sides
**		     collapsed to a point.  The 3 sides can be of any shape and
**		     the two sides adjacent to the collapsed side must have
**		     the same number of segments.
**		 - MshSurfLoftcoll:
**		     Linear mapping between a point and a curve, which form
**		     a triangular region.  The two sides of the region which
**		     are adjacent to the point of lofting must be straight 
**		     lines and have the same number of segments.
**		 - MshSurfTrimap:
**		     Trilinear mapping for a triangular region with sides of
**		     any shape.  The number of segments on all sides must 
**		     be equal.
**
** ---------------------------------------------------------------------------
**
**	Remarks :
**	The connectivity structure is formed by a vector that follows the
**	rule:
**	- The first position is the number of element nodes,
**	- The other nth positions is the connectivity in clockwise order,
**	  where nth is the number of element nodes.
**	- Analogous to the other elements.
**	For quadratic elements, the midside nodes must be given in the 
**	boundary coordinate vector.
**	For Q8 element it's necessary to allocate a ficticious node in
**	the middle of the face, although this node will never be used
**	in the element connectivity.
**	It was necessary to include the number of element nodes in the 
**	structure, since the algorithm triangularizes some regions which
**	will be impossible to put a four-node or eigth-node element.
**
** ---------------------------------------------------------------------------
**
**	Entry points :
** ---------------------------------------------------------------------------
**
**	void MshSurfTransfin(
**
**	double *bdynodes,   Boundary coordinate vector data X,Y coordinates
**                          of points on boundary given in clockwise order  
**	int    nu,          Number of nodes in first direction              
**	int    nv,          Number of nodes in the other direction          
**	int    elemtype,    Element type T3 (3), T6 (6), Q4 (4), or Q8 (8)  
**	int    diagtype,    Triangle option (option for cell diagonal)
**                      = 1 --> diagonal oriented to right direction
**                      = 2 --> diagonal oriented to left direction
**                      = 3 --> union jack alternation
**                      = 4 --> optimum diagonal (smallest of two possible)   
**	double *gennodes,   Element coordinates	    (out)
**	int    *genelems    Element Connectivity    (out)
**
**	MshSurfTransfin ( bilinear mapping meshing with any boundary shape ).
**	Based on the given points coordinates on the boundary , the number of
**	nodes in the two directions, and the element type, the routine fills
**	an array of element nodes coordinates, and a structure composed by
**  the number of element nodes and the connectivity.
**	If a triangular element is given, the diagonal orientation is also
**	required (diagtype).
**
** ---------------------------------------------------------------------------
**
**	void MshSurfLofting(
**
**	double *bdynodes,   boundary coordinate vector data X,Y coordinates
**	                    of points on boundary given in clockwise order
**	int    nu,          Number of nodes in first direction
**	int    nv,          Number of nodes in the other direction
**	int    dir,         Direction of lofting :
**                      (0 => first direction, 1 => other)
**	double  weight,     Lofting weight
**	int     elemtype,   Element type  = T3 (3), T6 (6), Q4 (4), or Q8 (8)
**	int     diagtype,   Triangle option (option for cell diagonal)
**	                    = 1 --> diagonal oriented to right direction
**	                    = 2 --> diagonal oriented to left direction
**	                    = 3 --> union jack alternation
**	                    = 4 --> optimum diagonal (smallest of two possible)
**	double *gennodes,   Element coordinates	                        (out)	
**	int    *genelems    Element genelemsectivity                    (out)
**
**	MshSurfLofting ( bilinear mapping meshing with straigth boundary )
**	Based on the given points coordinates on the boundary , the number of
**	nodes in the two directions, the element type, the direction of
**	lofting and the lofting weigth the routine fills an array of element
**	nodes coordinates, and a structure composed by the number of element
**	nodes and the connectivity.
**	If a triangular element is given, the diagonal orientation is also
**	required (diagtype).
**
** ---------------------------------------------------------------------------
**
**	void MshSurfTrsfncoll(
**
**	double *bdynodes,  boundary coordinate vector data X,Y coordinates
**	                   of points on boundary given in clockwise order 
**	int    nu,         Number of nodes in the side opposite
**	                   to the collapsed one 
**	int    nv,         Number of nodes in the other direction 
**	int    elemtype,   Element type = T3 (3), T6 (6), Q4 (4), or Q8 (8)
**	int    diagtype,   Triangle option (option for cell diagonal)
**	                   = 1 --> diagonal oriented to right direction
**	                   = 2 --> diagonal oriented to left direction
**	                   = 3 --> union jack alternation
**	                   = 4 --> optimum diagonal (smallest of two possible)
**	double *gennodes,  Element coordinates	                       (out)
**	int    *genelems   Element Connectivity                        (out)
**
**	MshSurfTrsfncoll ( bilinear collapsed mapping meshing with any
**	boundary shape ).
**	Based on the given points coordinates on the boundary , the number of
**	nodes in the two directions, and the element type, the routine fills
**	an array of element nodes coordinates, and a structure composed by
**  the number of element nodes and the connectivity.
**	If a triangular element is given, the diagonal orientation is also
**	required (diagtype).
**
** ---------------------------------------------------------------------------
**
**	void MshSurfLoftcoll(
**	
**	double *bdynodes,  boundary coordinate vector data X,Y coordinates
**	                   of points on boundary given in clockwise order
**	int    nu,         Numb. of nodes in the side opposite to the
**	                   collapsed one
**	int    nv,         Number of nodes in the other direction 
**	double  weight,    Lofting weight 
**	int    elemtype,   Element type = T3 (3), T6 (6), Q4 (4), or Q8 (8) 
**	int    diagtype,   Triangle option (option for cell diagonal)
**	                   = 1 --> diagonal oriented to right direction
**	                   = 2 --> diagonal oriented to left direction
**	                   = 3 --> union jack alternation
**	                   = 4 --> optimum diagonal (smallest of two possible) 
**	double *gennodes,  Element coordinates                         (out)
**	int    *genelems   Element Connectivity                        (out)
**	
**	MshSurfLoftcoll ( bilinear collapsed mapping meshing with straigth
**	boundary shape ).
**	Based on the given points coordinates on the boundary , the number of
**	nodes in the two directions, and the element type, the routine fills
**	an array of element nodes coordinates, and a structure composed by
**  the number of element nodes and the connectivity.
**	If a triangular element is given, the diagonal orientation is also
**	required (diagtype).
**
** ---------------------------------------------------------------------------
**
**	void MshSurfTrimap(
**
**	double   *bdynodes,    boundary coordinate vector data X,Y coordinates
**	                       of points on boundary given in clockwise order
**	int      n,            Number of nodes on one boundary
**	int      elemtype,     Element type= T3 (3), T6 (6), Q4 (4), or Q8 (8) 
**	double   *gennodes,    Element coordinates          (out)
**	int  	 *genelems     Element Connectivity         (out)
**
**	MshSurfTrimap ( trilinear mapping meshing with any boundary shape ).
**	Based on the given points coordinates on the boundary , the number of
**	nodes in the three directions, and the element type, the routine fills
**	an array of element nodes coordinates, and a structure composed by
**  the number of element nodes and the connectivity.
**
** ---------------------------------------------------------------------------
*/
	
/* Global variables and symbols           */
#include 	<math.h>
#include <stdlib.h>
#include <string.h>
#include 	"mshsurfmapp.h"
#include "mshsurf_geo.h"

typedef struct _MshPnt3D3d
{
    double x, y, z;
} MshPnt3D ;



/* Local variables and symbols            */

#define  SQR(a) ((a)*(a))

#define id(i,j,m) (((i)*(m))+(j))


/* static functions
***************************************************************************/
static int    compInt ( const void *c1, const void *c2 );
static int    MshSurfGetMappDirection4_( int *bd, int nPtsLoop, int *m, int *n, int *pos);
static int    MshSurfGetMappDirection5_( int *bd, int nPtsLoop, int *m, int *n, int *pos);
static int    MshSurfGetMappDirection6_( int *bd, int nPtsLoop, int *m, int *n, int *pos);

static int    MshSurfGetMappTriDirection3_( int *bd, int nPtsLoop, int *n, int *pos);
static int    MshSurfGetMappTriDirection4_( int *bd, int nPtsLoop, int *n, int *pos);
static int    MshSurfGetMappTriDirection5_( int *bd, int nPtsLoop, int *n, int *pos);


/*=======================================================================*/
/*============================ MshSurfTransfin ============================*/
/*=======================================================================*/


void MshSurfTransfin(

double *bdynodes,   /* boundary coordinate vector data X,Y coordinates
                     of points on boundary given in clockwise order     */
int    nu,          /* Number of nodes in first direction               */
int    nv,          /* Number of nodes in the other direction           */
int    elemtype,    /* Element type T3 (3), T6 (6), Q4 (4), or Q8 (8)   */
int    diagtype,    /* Triangle option (option for cell diagonal)
                    = 1 --> diagonal oriented to right direction
                    = 2 --> diagonal oriented to left direction
                    = 3 --> union jack alternation
                    = 4 --> optimum diagonal (smallest of two possible) */
double   *gennodes, /* Element coordinates          (out)               */
int  	 *genelems  /* Element genelemsectivity     (out)               */
)

{

/*  Local variables: */

  MshPnt3D *def_nds;
  MshPnt3D *pt;

  int     in1,in2,in3,in4;
  int     node_inc;
  int     ip;
  int     ident = 0;
  double   dis1,dis2;
  int     i,j,k,l;

/* Processing begins here */

  def_nds = (MshPnt3D *)bdynodes;
  pt = (MshPnt3D *)gennodes;

/*
* Get node increment
*/
  if(elemtype == 6 || elemtype == 8)
  {
    node_inc = 2;
  }
  else
  {
    node_inc = 1;
  }
/*
* Make sure triangular option are set only for triangular elements
*/

  if(elemtype == 4  ||  elemtype == 8)
  {
    diagtype = 0;
  } 
/*
* Get the starting index (offset position) of each boundary on boundary array
*/
  in1 = 0;
  in2 = in1 + nu - 1 ;
  in3 = in2 + nv - 1;
  in4 = in3 + nu - 1 ;
/*
* Input boundary nodes
*/
  for (i = 0;i < (nu-1);i++)
  {
    ip = in1 + i;
    pt[id(0,i,nu)].x = def_nds[ip].x;
    pt[id(0,i,nu)].y = def_nds[ip].y;
    pt[id(0,i,nu)].z = def_nds[ip].z;
  }
  for (i = 0;i < (nv-1);i++)
  {
    ip = in2 + i;
    pt[id(i,nu-1,nu)].x = def_nds[ip].x;
    pt[id(i,nu-1,nu)].y = def_nds[ip].y;
    pt[id(i,nu-1,nu)].z = def_nds[ip].z;
  }
  for (i = 1;i < nu;i++)
  {
    ip = in3 + (nu-i-1);
    pt[id(nv-1,i,nu)].x = def_nds[ip].x;
    pt[id(nv-1,i,nu)].y = def_nds[ip].y;
    pt[id(nv-1,i,nu)].z = def_nds[ip].z;
  }
  for (i = 1;i < nv;i++)
  {
    ip = in4 + (nv-i-1);
    pt[id(i,0,nu)].x = def_nds[ip].x;
    pt[id(i,0,nu)].y = def_nds[ip].y;
    pt[id(i,0,nu)].z = def_nds[ip].z;
  }
/*
* Generate interior nodes
*/
  for (k = node_inc;k < nu-node_inc;k = k+node_inc)
  {
    for (j = node_inc;j < nv-node_inc;j = j + node_inc)
    {
      pt[id(j,k,nu)].x=
         pt[id(0,k,nu)].x*(nv-j-1)/(nv-1) +
         pt[id(nv-1,k,nu)].x*j/(nv-1) +
         pt[id(j,0,nu)].x*(nu-k-1)/(nu-1) +
         pt[id(j,nu-1,nu)].x*k/(nu-1) -
         pt[id(0,0,nu)].x*(nu-k-1)/(nu-1)*(nv-j-1)/(nv-1) -
         pt[id(nv-1,0,nu)].x*(nu-k-1)/(nu-1)*j/(nv-1) -
         pt[id(nv-1,nu-1,nu)].x*k/(nu-1)*j/(nv-1) -
         pt[id(0,nu-1,nu)].x*k/(nu-1)*(nv-j-1)/(nv-1);

      pt[id(j,k,nu)].y=
         pt[id(0,k,nu)].y*(nv-j-1)/(nv-1) +
         pt[id(nv-1,k,nu)].y*j/(nv-1) +
         pt[id(j,0,nu)].y*(nu-k-1)/(nu-1) +
         pt[id(j,nu-1,nu)].y*k/(nu-1) -
         pt[id(0,0,nu)].y*(nu-k-1)/(nu-1)*(nv-j-1)/(nv-1) -
         pt[id(nv-1,0,nu)].y*(nu-k-1)/(nu-1)*j/(nv-1) -
         pt[id(nv-1,nu-1,nu)].y*k/(nu-1)*j/(nv-1) -
         pt[id(0,nu-1,nu)].y*k/(nu-1)*(nv-j-1)/(nv-1);

      pt[id(j,k,nu)].z=
         pt[id(0,k,nu)].z*(nv-j-1)/(nv-1) +
         pt[id(nv-1,k,nu)].z*j/(nv-1) +
         pt[id(j,0,nu)].z*(nu-k-1)/(nu-1) +
         pt[id(j,nu-1,nu)].z*k/(nu-1) -
         pt[id(0,0,nu)].z*(nu-k-1)/(nu-1)*(nv-j-1)/(nv-1) -
         pt[id(nv-1,0,nu)].z*(nu-k-1)/(nu-1)*j/(nv-1) -
         pt[id(nv-1,nu-1,nu)].z*k/(nu-1)*j/(nv-1) -
         pt[id(0,nu-1,nu)].z*k/(nu-1)*(nv-j-1)/(nv-1);
    }
  }
/*
* Generate midside nodes
*/
  if (node_inc == 2)
  {
    for (k = 2;k < nu-2;k = k+2) 
    {
      for (j = 1;j < nv-1;j = j+2)
      {
        pt[id(j,k,nu)].x =
          (pt[id(j-1,k,nu)].x + pt[id(j+1,k,nu)].x) / 2.0;
        pt[id(j,k,nu)].y =
          (pt[id(j-1,k,nu)].y + pt[id(j+1,k,nu)].y) / 2.0;
        pt[id(j,k,nu)].z =
          (pt[id(j-1,k,nu)].z + pt[id(j+1,k,nu)].z) / 2.0;
      }
    }
    for (k = 1;k < nu-1;k = k+2)
    {
      for (j = 2;j < nv-2;j = j+2)
      {
        pt[id(j,k,nu)].x =
          (pt[id(j,k-1,nu)].x + pt[id(j,k+1,nu)].x) / 2.0;
        pt[id(j,k,nu)].y =
          (pt[id(j,k-1,nu)].y + pt[id(j,k+1,nu)].y) / 2.0;
        pt[id(j,k,nu)].z =
          (pt[id(j,k-1,nu)].z + pt[id(j,k+1,nu)].z) / 2.0;
      }
     }
  }
/*
* Get code for diagonal if we have simple cases of right or left
*/
  if(diagtype == 1)
  {
    ident = 1;
  }
  else
  {
    if(diagtype == 2)
    {
       ident = 2;
    }
  }
/*
* Start generating
*/

/*
* Initialize the genelemsectivity vector
*/
  l=0;

  for (k = 1;k < nu+1-node_inc;k = k+node_inc)
  {
    for (j = 1;j < nv+1-node_inc;j = j+node_inc)
    {
/*
*    Get code for diagonal if we have union jack option
*/
       if(diagtype == 3)
       {
         if( (((j+k)*(2/node_inc)) % 4) == 0)
         {
           ident = 1;
         }
         else
         {
           ident = 2;
         }
       }
/*
*    Get code for diagonal based on diagonal sizes if we have optimum option
*/
       else if(diagtype == 4)
       {
         if (node_inc == 2)
         {
           dis1 = sqrt(SQR((pt[id(j-1,k-1,nu)].x-pt[id(j+1,k+1,nu)].x)) +
                       SQR((pt[id(j-1,k-1,nu)].y-pt[id(j+1,k+1,nu)].y)) +
                       SQR((pt[id(j-1,k-1,nu)].z-pt[id(j+1,k+1,nu)].z)) );
           dis2 = sqrt(SQR((pt[id(j+1,k-1,nu)].x-pt[id(j-1,k+1,nu)].x)) +
                       SQR((pt[id(j+1,k-1,nu)].y-pt[id(j-1,k+1,nu)].y)) +
                       SQR((pt[id(j+1,k-1,nu)].z-pt[id(j-1,k+1,nu)].z)) );
         }
         else
         {
           dis1 = sqrt(SQR((pt[id(j-1,k-1,nu)].x-pt[id(j  ,k  ,nu)].x)) +
                       SQR((pt[id(j-1,k-1,nu)].y-pt[id(j  ,k  ,nu)].y)) + 
                       SQR((pt[id(j-1,k-1,nu)].z-pt[id(j  ,k  ,nu)].z)) );
           dis2 = sqrt(SQR((pt[id(j  ,k-1,nu)].x-pt[id(j-1,k  ,nu)].x)) +
                       SQR((pt[id(j  ,k-1,nu)].y-pt[id(j-1,k  ,nu)].y)) + 
                       SQR((pt[id(j  ,k-1,nu)].z-pt[id(j-1,k  ,nu)].z)) );
         }
			
         if(dis1<=dis2)
         {
           ident = 1;
         }
         else
         {
           ident = 2;
         }
       }

/*
*     Generate T3 triangular element
*/
      if(elemtype == 3)
      {
        if(ident == 1)
        {
          genelems[l++] = 3;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j  ,k-1,nu);
		
          genelems[l++] = 3;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j  ,k  ,nu);
        }
        else
        {
          genelems[l++] = 3;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j  ,k-1,nu);

          genelems[l++] = 3;
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j  ,k-1,nu);
        }
      }
/*
*     Generate Q4 quadrilateral element
*/
      else
      if(elemtype == 4)
      {
         genelems[l++] = 4;
         genelems[l++] = id(j-1,k-1,nu);
         genelems[l++] = id(j-1,k  ,nu);
         genelems[l++] = id(j  ,k  ,nu);
         genelems[l++] = id(j  ,k-1,nu);
      }

/*
*     Generate T6 triangular element
*/
      else
      if(elemtype == 6)
      {
        if(ident == 1)
        {
/*
*          Coordinates for node on inclined edge
*/
           pt[id(j,k,nu)].x =
            (pt[id(j-1,k-1,nu)].x + pt[id(j+1,k+1,nu)].x)/2.0;
           pt[id(j,k,nu)].y =
            (pt[id(j-1,k-1,nu)].y + pt[id(j+1,k+1,nu)].y)/2.0;
           pt[id(j,k,nu)].z =
            (pt[id(j-1,k-1,nu)].z + pt[id(j+1,k+1,nu)].z)/2.0;

           genelems[l++] = 6;
           genelems[l++] = id(j-1,k-1,nu);
           genelems[l++] = id(j  ,k  ,nu);
           genelems[l++] = id(j+1,k+1,nu);
           genelems[l++] = id(j+1,k  ,nu);
           genelems[l++] = id(j+1,k-1,nu);
           genelems[l++] = id(j  ,k-1,nu);

           genelems[l++] = 6;
           genelems[l++] = id(j-1,k-1,nu);
           genelems[l++] = id(j-1,k  ,nu);
           genelems[l++] = id(j-1,k+1,nu);
           genelems[l++] = id(j  ,k+1,nu);
           genelems[l++] = id(j+1,k+1,nu);
           genelems[l++] = id(j  ,k  ,nu);
        }
        else
        {
/*
*          Coordinates for node on inclined edge
*/
          pt[id(j,k,nu)].x =
            (pt[id(j+1,k-1,nu)].x + pt[id(j-1,k+1,nu)].x)/2.0;
          pt[id(j,k,nu)].y =
            (pt[id(j+1,k-1,nu)].y + pt[id(j-1,k+1,nu)].y)/2.0;
          pt[id(j,k,nu)].z =
            (pt[id(j+1,k-1,nu)].z + pt[id(j-1,k+1,nu)].z)/2.0;

          genelems[l++] = 6;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j-1,k+1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j+1,k-1,nu);

          genelems[l++] = id(j  ,k-1,nu);
          genelems[l++] = 6;
          genelems[l++] = id(j+1,k-1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j-1,k+1,nu);
          genelems[l++] = id(j  ,k+1,nu);
          genelems[l++] = id(j+1,k+1,nu);
          genelems[l++] = id(j+1,k  ,nu);
        }
      }
/*
*     Generate Q8 quadrilateral element
*/
      else
      if(elemtype == 8)
      {
         genelems[l++] = 8;
         genelems[l++] = id(j-1,k-1,nu);
         genelems[l++] = id(j-1,k  ,nu);
         genelems[l++] = id(j-1,k+1,nu);
         genelems[l++] = id(j  ,k+1,nu);
         genelems[l++] = id(j+1,k+1,nu);
         genelems[l++] = id(j+1,k  ,nu);
         genelems[l++] = id(j+1,k-1,nu);
         genelems[l++] = id(j,  k-1,nu);
      }

    }
  }
}



/*============================================================================*/
/*============================== MshSurfLofting ================================*/
/*============================================================================*/

void MshSurfLofting(

double * bdynodes,	/* boundary coordinate vector data X,Y coordinates
                      of points on boundary given in clockwise order          */
int    nu,          /* Number of nodes in first direction                     */
int    nv,          /* Number of nodes in the other direction                 */
int    dir,         /* Direction of lofting (0 => first direction, 1 => other)*/
double  weight,		/* Lofting weight                                         */
int    elemtype,	/* Element type  = T3 (3), T6 (6), Q4 (4), or Q8 (8)      */
int    diagtype,	/* Triangle option (option for cell diagonal)
                       = 1 --> diagonal oriented to right direction
                       = 2 --> diagonal oriented to left direction
                       = 3 --> union jack alternation
                       = 4 --> optimum diagonal (smallest of two possible)    */
double   *gennodes,	/* Element coordinates			(out)                     */
int  	 *genelems	/* Element genelemsectivity			(out)                     */
)

{
/*  Local variables: */
  MshPnt3D *def_nds;
  MshPnt3D *pt;

  int      loft_seg;
  int      in1,in2,in3,in4;
  int      node_inc;
  int      ip;
  int      ident = 0;
  double   dis1,dis2;
  double   aa,bb,cc,u,v;
  int      i,j,k,l;


/* Processing begins here */

  def_nds = (MshPnt3D *)bdynodes;
  pt = (MshPnt3D *)gennodes;


/*
* Get node increment
*/
  if(elemtype == 6  ||  elemtype == 8)
  {
    node_inc = 2;
  }
  else
  {
    node_inc = 1;
  }
/*
* Get number of segments for lofting based on the given direction
*/
  if (dir  ==  0)
  {
    loft_seg = (nu - 1) / node_inc;
  }
  else
  {
    loft_seg = (nv - 1) / node_inc;
  }

/*
* Make sure triangular option are set only for triangular elements
*/
  if(elemtype == 4  ||  elemtype == 8)
  {
    diagtype = 0;
  }
/*
* Get the starting index (offset position) of each boundary on boundary array
*/
  in1 = 0;
  in2 = in1 + (nu - 1);
  in3 = in2 + (nv - 1);
  in4 = in3 + (nu - 1);
/*
* Input boundary nodes
*/
  for (i = 0;i < nu-1; i++)
  {
    ip = in1 + i;
    pt[id(0,i,nu)].x = def_nds[ip].x;
    pt[id(0,i,nu)].y = def_nds[ip].y;
    pt[id(0,i,nu)].z = def_nds[ip].z;
  }
  for (i = 0;i < nv-1;i++)
  {
    ip = in2 + i;
    pt[id(i,nu-1,nu)].x = def_nds[ip].x;
    pt[id(i,nu-1,nu)].y = def_nds[ip].y;
    pt[id(i,nu-1,nu)].z = def_nds[ip].z;
  }
  for (i = 1;i < nu;i++) 
  {
    ip = in3 + (nu-i-1);
    pt[id(nv-1,i,nu)].x = def_nds[ip].x;
    pt[id(nv-1,i,nu)].y = def_nds[ip].y;
    pt[id(nv-1,i,nu)].z = def_nds[ip].z;
  }
  for (i= 1;i < nv;i++) 
  {
    ip = in4 + (nv-i-1);
    pt[id(i,0,nu)].x = def_nds[ip].x;
    pt[id(i,0,nu)].y = def_nds[ip].y;
    pt[id(i,0,nu)].z = def_nds[ip].z;
  }
/*
* Generate interior nodes
*
*   First compute auxiliar coefficients
*/
  if (loft_seg  >  1)
  { 
    aa = (2.0 * (weight)) / (((weight) + 1.0) * (double)(loft_seg));
    bb = (aa * (1.0 - (weight))) /
         (2.0 * (weight) * ((double)(loft_seg) - 1.0));

    if (dir  ==  0)
    {
      for (j = node_inc;j < nv-node_inc;j = j+node_inc) 
      {
        for (k = node_inc;k < nu-node_inc;k = k + node_inc)
        {

          cc = (double)(k) / (double)(node_inc);
          v = aa * cc + bb * cc * (cc - 1.0);
          u = 1.0 - v;

          pt[id(j,k,nu)].x = u * pt[id(j,0,nu)].x  +  v * pt[id(j,nu-1,nu)].x;
          pt[id(j,k,nu)].y = u * pt[id(j,0,nu)].y  +  v * pt[id(j,nu-1,nu)].y;
          pt[id(j,k,nu)].z = u * pt[id(j,0,nu)].z  +  v * pt[id(j,nu-1,nu)].z;

        }
      }
    }
    else
    {
      for (k = node_inc;k < nu-node_inc;k = k + node_inc)
      {
	    for (j = node_inc;j < nv-node_inc;j = j + node_inc)
        {

           cc = (double)(j) / (double)(node_inc);
           v = aa * cc + bb * cc * (cc - 1.0);
           u = 1.0 - v;

           pt[id(j,k,nu)].x = u * pt[id(0,k,nu)].x  +  v * pt[id(nv-1,k,nu)].x;
           pt[id(j,k,nu)].y = u * pt[id(0,k,nu)].y  +  v * pt[id(nv-1,k,nu)].y;
           pt[id(j,k,nu)].z = u * pt[id(0,k,nu)].z  +  v * pt[id(nv-1,k,nu)].z;

        }
      }

    }
  }
/*
* Generate midside nodes
*/
  if(node_inc == 2)
  {
    for (k = 2;k < nu-2;k = k+2)
    {
	  for (j = 1;j < nv-1;j = j+2) 
      {
         pt[id(j,k,nu)].x = (pt[id(j-1,k,nu)].x + pt[id(j+1,k,nu)].x) / 2.0;
         pt[id(j,k,nu)].y = (pt[id(j-1,k,nu)].y + pt[id(j+1,k,nu)].y) / 2.0;
         pt[id(j,k,nu)].z = (pt[id(j-1,k,nu)].z + pt[id(j+1,k,nu)].z) / 2.0;
      }
    }

    for (k = 1;k < nu-1;k = k+2)
    {
      for (j = 2;j < nv-2;j = j+2)
      {
         pt[id(j,k,nu)].x = (pt[id(j,k-1,nu)].x + pt[id(j,k+1,nu)].x) / 2.0;
         pt[id(j,k,nu)].y = (pt[id(j,k-1,nu)].y + pt[id(j,k+1,nu)].y) / 2.0;
         pt[id(j,k,nu)].z = (pt[id(j,k-1,nu)].z + pt[id(j,k+1,nu)].z) / 2.0;
      }
    }
  }
/*
* Get code for diagonal if we have simple cases of right or left
*/
  if(diagtype == 1)
  {
     ident = 1;
  }
  else if(diagtype == 2)
  {
     ident = 2;
  }
/*
* Start generating
*/

/* 
* Initialize the genelemsectivity vector
*/
  l=0;

  for (k = 1;k < nu+1-node_inc;k = k+node_inc)
  {
    for (j = 1;j < nv+1-node_inc;j = j+node_inc)
    {
/*
*    Get code for diagonal if we have union jack option
*/

      if (diagtype == 3)
      {
        if( (((j+k)*(2/node_inc)) % 4) == 0)
        {
           ident = 1;
        }
        else
        {
           ident = 2;
        }
      }
/*
*    Get code for diagonal based on diagonal sizes if we have optimum option
*/
      else if(diagtype == 4)
      {
        if (node_inc == 2)
        {
          dis1 =
            sqrt(SQR((pt[id(j-1,k-1,nu)].x-pt[id(j+1,k+1,nu)].x)) +
                 SQR((pt[id(j-1,k-1,nu)].y-pt[id(j+1,k+1,nu)].y)) +
                 SQR((pt[id(j-1,k-1,nu)].z-pt[id(j+1,k+1,nu)].z)) );
           dis2 = 
              sqrt(SQR((pt[id(j+1,k-1,nu)].x-pt[id(j-1,k+1,nu)].x)) +
                   SQR((pt[id(j+1,k-1,nu)].y-pt[id(j-1,k+1,nu)].y)) +
                   SQR((pt[id(j+1,k-1,nu)].z-pt[id(j-1,k+1,nu)].z)) );
        }
        else
        {
          dis1 = 
            sqrt(SQR((pt[id(j-1,k-1,nu)].x-pt[id(j  ,k  ,nu)].x)) +
                 SQR((pt[id(j-1,k-1,nu)].y-pt[id(j  ,k  ,nu)].y)) +
                 SQR((pt[id(j-1,k-1,nu)].z-pt[id(j  ,k  ,nu)].z)) );
          dis2 =
            sqrt(SQR((pt[id(j  ,k-1,nu)].x-pt[id(j-1,k  ,nu)].x)) +
                 SQR((pt[id(j  ,k-1,nu)].y-pt[id(j-1,k  ,nu)].y)) +
                 SQR((pt[id(j  ,k-1,nu)].z-pt[id(j-1,k  ,nu)].z)) );
        }

        if(dis1 <= dis2)
        {
           ident = 1;
        }
        else
        {
            ident = 2;
        }
      }

/*
*     Generate T3 triangular element
*/
      if(elemtype == 3)
      {
        if(ident == 1)
        {
          genelems[l++] = 3;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j  ,k-1,nu);

          genelems[l++] = 3;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j  ,k  ,nu);
        }
        else
        {
          genelems[l++] = 3;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j  ,k-1,nu);

          genelems[l++] = 3;
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j  ,k-1,nu);
        }
      }
/*
*     Generate Q4 quadrilateral element
*/
      else if(elemtype == 4)
      {
        genelems[l++] = 4;
        genelems[l++] = id(j-1,k-1,nu);
        genelems[l++] = id(j-1,k  ,nu);
        genelems[l++] = id(j  ,k  ,nu);
        genelems[l++] = id(j  ,k-1,nu);
      }
/*
*     Generate T6 triangular element
*/
      else if(elemtype == 6)
      {
        if(ident == 1)
        {
/*
*          Coordinates for node on inclined edge
*/
          pt[id(j,k,nu)].x = (pt[id(j-1,k-1,nu)].x + pt[id(j+1,k+1,nu)].x)/2.0;
          pt[id(j,k,nu)].y = (pt[id(j-1,k-1,nu)].y + pt[id(j+1,k+1,nu)].y)/2.0;
          pt[id(j,k,nu)].z = (pt[id(j-1,k-1,nu)].z + pt[id(j+1,k+1,nu)].z)/2.0;

          genelems[l++] = 6;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j+1,k+1,nu);
          genelems[l++] = id(j+1,k  ,nu);
          genelems[l++] = id(j+1,k-1,nu);
          genelems[l++] = id(j  ,k-1,nu);

          genelems[l++] = 6;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j-1,k+1,nu);
          genelems[l++] = id(j  ,k+1,nu);
          genelems[l++] = id(j+1,k+1,nu);
          genelems[l++] = id(j  ,k  ,nu);
        }
        else
        {
/*
*          Coordinates for node on inclined edge
*/
          pt[id(j,k,nu)].x = (pt[id(j+1,k-1,nu)].x + pt[id(j-1,k+1,nu)].x)/2.0;
          pt[id(j,k,nu)].y = (pt[id(j+1,k-1,nu)].y + pt[id(j-1,k+1,nu)].y)/2.0;
          pt[id(j,k,nu)].z = (pt[id(j+1,k-1,nu)].z + pt[id(j-1,k+1,nu)].z)/2.0;

          genelems[l++] = 6;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j-1,k+1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j+1,k-1,nu);
          genelems[l++] = id(j  ,k-1,nu);

          genelems[l++] = 6;
          genelems[l++] = id(j+1,k-1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j-1,k+1,nu);
          genelems[l++] = id(j  ,k+1,nu);
          genelems[l++] = id(j+1,k+1,nu);
          genelems[l++] = id(j+1,k  ,nu);
        }
      }
/*
*     Generate Q8 quadrilateral element
*/
      else if(elemtype == 8)
      {
        genelems[l++] = 8;
        genelems[l++] = id(j-1,k-1,nu);
        genelems[l++] = id(j-1,k  ,nu);
        genelems[l++] = id(j-1,k+1,nu);
        genelems[l++] = id(j  ,k+1,nu);
        genelems[l++] = id(j+1,k+1,nu);
        genelems[l++] = id(j+1,k  ,nu);
        genelems[l++] = id(j+1,k-1,nu);
        genelems[l++] = id(j,  k-1,nu);
      }
    }
  }
}


/*============================================================================*/
/* ============================ Mhs2DTrsfncoll ===============================*/
/*============================================================================*/

void MshSurfTrsfncoll(

 double * bdynodes, /* boundary coordinate vector data X,Y coordinates
                       of points on boundary given in clockwise order  */
 int    nu,         /* Number of nodes in the side opposite
                       to the collapsed one */
 int    nv,         /* Number of nodes in the other direction */
 int    elemtype,   /* Element type = T3 (3), T6 (6), Q4 (4), or Q8 (8) */
 int    diagtype,   /* Triangle option (option for cell diagonal)
                        = 1 --> diagonal oriented to right direction
                        = 2 --> diagonal oriented to left direction
                        = 3 --> union jack alternation
                        = 4 --> optimum diagonal (smallest of two possible) */
 double *gennodes,    /* Element coordinates(out)*/
 int    *genelems     /* Element Connectivity(out)*/
)

{
/*  Local variables:  */
  MshPnt3D *def_nds;
  MshPnt3D *pt;

  int     in1,in2,in3,in4;
  int     node_inc;
  int     ip;
  int     ident = 0;
  double  dis1,dis2;
  int     i,j,k,l;


/* Processing begins here */
  def_nds = (MshPnt3D *)bdynodes;
  pt = (MshPnt3D *)gennodes;

/*
* Get node increment
*/
  if((elemtype == 6) || (elemtype == 8))
  {
    node_inc = 2;
  }
  else
  {
    node_inc = 1;
  }
/*
* Make sure triangular option are set only for triangular elements
*/
  if((elemtype == 4) || (elemtype == 8))
  {
    diagtype = 0;
  }
/*
* Get the starting index (offset position) of each boundary on boundary array
*/
  in1 = 0;
  in2 = in1;
  in3 = in2 + (nv - 1);
  in4 = in3 + (nu - 1);
/*
* Input boundary nodes
*/
  for (i = 0;i < nu-1;i++)
  {
    ip = in1;
    pt[id(0,i,nu)].x = def_nds[ip].x;
    pt[id(0,i,nu)].y = def_nds[ip].y;
    pt[id(0,i,nu)].z = def_nds[ip].z;
  }
  for (i = 0;i < nv-1;i++)
  {
    ip = in2 + i;
    pt[id(i,nu-1,nu)].x = def_nds[ip].x;
    pt[id(i,nu-1,nu)].y = def_nds[ip].y;
    pt[id(i,nu-1,nu)].z = def_nds[ip].z;
  }
  for (i = 1;i < nu;i++)
  {
    ip = in3 + (nu-i-1);
    pt[id(nv-1,i,nu)].x = def_nds[ip].x;
    pt[id(nv-1,i,nu)].y = def_nds[ip].y;
    pt[id(nv-1,i,nu)].z = def_nds[ip].z;
  }
  for (i = 1;i < nv;i++)
  {
    ip = in4 + (nv-i-1);
    pt[id(i,0,nu)].x = def_nds[ip].x;
    pt[id(i,0,nu)].y = def_nds[ip].y;
    pt[id(i,0,nu)].z = def_nds[ip].z;
  }
/*
* Generate interior nodes
*/
  for (k = node_inc;k < nu-node_inc;k = k+node_inc)
  {
    for (j = node_inc;j < nv-node_inc;j = j+node_inc)
    {
      pt[id(j,k,nu)].x =
       pt[id(0,k,nu)].x*(nv-j-1)/(nv-1) +
       pt[id(nv-1,k,nu)].x*j/(nv-1) +
       pt[id(j,0,nu)].x*(nu-k-1)/(nu-1) +
       pt[id(j,nu-1,nu)].x*k/(nu-1) -
       pt[id(0,0,nu)].x*(nu-k-1)/(nu-1)*(nv-j-1)/(nv-1)  -
       pt[id(nv-1,0,nu)].x*(nu-k-1)/(nu-1)*j/(nv-1)-
       pt[id(nv-1,nu-1,nu)].x*k/(nu-1)*j/(nv-1)-
       pt[id(0,nu-1,nu)].x*k/(nu-1)*(nv-j-1)/(nv-1);

      pt[id(j,k,nu)].y =
       pt[id(0,k,nu)].y*(nv-j-1)/(nv-1) +
       pt[id(nv-1,k,nu)].y*j/(nv-1) +
       pt[id(j,0,nu)].y*(nu-k-1)/(nu-1) +
       pt[id(j,nu-1,nu)].y*k/(nu-1) -
       pt[id(0,0,nu)].y*(nu-k-1)/(nu-1)*(nv-j-1)/(nv-1)-
       pt[id(nv-1,0,nu)].y*(nu-k-1)/(nu-1)*j/(nv-1)-
       pt[id(nv-1,nu-1,nu)].y*k/(nu-1)*j/(nv-1)-
       pt[id(0,nu-1,nu)].y*k/(nu-1)*(nv-j-1)/(nv-1);

      pt[id(j,k,nu)].z =
       pt[id(0,k,nu)].z*(nv-j-1)/(nv-1) +
       pt[id(nv-1,k,nu)].z*j/(nv-1) +
       pt[id(j,0,nu)].z*(nu-k-1)/(nu-1) +
       pt[id(j,nu-1,nu)].z*k/(nu-1) -
       pt[id(0,0,nu)].z*(nu-k-1)/(nu-1)*(nv-j-1)/(nv-1)-
       pt[id(nv-1,0,nu)].z*(nu-k-1)/(nu-1)*j/(nv-1)-
       pt[id(nv-1,nu-1,nu)].z*k/(nu-1)*j/(nv-1)-
       pt[id(0,nu-1,nu)].z*k/(nu-1)*(nv-j-1)/(nv-1);
    }
  }

/*
* Generate midside nodes
*/
  if(node_inc == 2)
  {
    for (k = 2;k < nu-2;k = k+2)
    {
      for (j = 1;j < nv-1;j = j+2)
      {
        pt[id(j,k,nu)].x = (pt[id(j-1,k,nu)].x + pt[id(j+1,k,nu)].x) / 2.0;
        pt[id(j,k,nu)].y = (pt[id(j-1,k,nu)].y + pt[id(j+1,k,nu)].y) / 2.0;
        pt[id(j,k,nu)].z = (pt[id(j-1,k,nu)].z + pt[id(j+1,k,nu)].z) / 2.0;
      }
    }

    for (k = 1;k < nu-1;k = k+2)
    {
      for (j = 2;j < nv-2;j = j+2)
      {
        pt[id(j,k,nu)].x = (pt[id(j,k-1,nu)].x + pt[id(j,k+1,nu)].x) / 2.0;
        pt[id(j,k,nu)].y = (pt[id(j,k-1,nu)].y + pt[id(j,k+1,nu)].y) / 2.0;
        pt[id(j,k,nu)].z = (pt[id(j,k-1,nu)].z + pt[id(j,k+1,nu)].z) / 2.0;
      }
    }
  }
/*
* Get code for diagonal if we have simple cases of right or left
*/
  if(diagtype == 1)
  {
    ident = 2;
  }
  else if(diagtype == 2)
  {
    ident = 1;
  }
/*
* Start generating
*/

/* 
* Initialize the Connectivity vector
*/
  l=0;

  for (k = 1;k < nu+1-node_inc;k = k+node_inc)
  {
/*
*     The first element of this row is a triangular element close
*     to the collapsed side.
*/
    j = 1;
/*
*     Generate T3 triangular element or Q4 quadrilateral element
*/
    if((elemtype == 3) || (elemtype == 4))
    {
      genelems[l++] = 3;
      genelems[l++] = id(j-1,k-1,nu);
      genelems[l++] = id(j  ,k  ,nu);
      genelems[l++] = id(j  ,k-1,nu);
    }
/*
*     Generate T6 triangular element or Q8 quadrilateral element
*/
    else if((elemtype == 6) || (elemtype == 8))
    {
/*
*          Coordinates for node on inclined edge
*/
      if (k < nu-1)
      {
        pt[id(j,k,nu)].x = (pt[id(j-1,k-1,nu)].x + pt[id(j+1,k+1,nu)].x)/2.0;
        pt[id(j,k,nu)].y = (pt[id(j-1,k-1,nu)].y + pt[id(j+1,k+1,nu)].y)/2.0;
        pt[id(j,k,nu)].z = (pt[id(j-1,k-1,nu)].z + pt[id(j+1,k+1,nu)].z)/2.0;
      }
      else
      {
        pt[id(j,k,nu)].x = pt[id(2,nu,nu)].x;
        pt[id(j,k,nu)].y = pt[id(2,nu,nu)].y;
        pt[id(j,k,nu)].z = pt[id(2,nu,nu)].z;
      }

      genelems[l++] = 6;
      genelems[l++] = id(j-1,k-1,nu);
      genelems[l++] = id(j  ,k  ,nu);
      genelems[l++] = id(j+1,k+1,nu);
      genelems[l++] = id(j+1,k  ,nu);
      genelems[l++] = id(j+1,k-1,nu);
      genelems[l++] = id(j  ,k-1,nu);
    }
/*
*   Now the remaining elements of this row
*/
    for (j = 1+node_inc;j < nv+1-node_inc;j = j+node_inc)
    {
/*
*    Get code for diagonal if we have union jack option
*/
      if (diagtype == 3)
      {
        if( (((j+k)*(2/node_inc)) % 4) == 0)
        {
          ident = 1;
        }
        else
        {
          ident = 2;
        }
      }
/*
*    Get code for diagonal based on diagonal sizes if we have optimum option
*/
      else if(diagtype == 4)
      {
        if (node_inc == 2)
        {
          dis1 = sqrt(SQR((pt[id(j-1,k-1,nu)].x-pt[id(j+1,k+1,nu)].x)) +
                      SQR((pt[id(j-1,k-1,nu)].y-pt[id(j+1,k+1,nu)].y)) +
                      SQR((pt[id(j-1,k-1,nu)].z-pt[id(j+1,k+1,nu)].z)) );
          dis2 = sqrt(SQR((pt[id(j+1,k-1,nu)].x-pt[id(j-1,k+1,nu)].x)) +
                      SQR((pt[id(j+1,k-1,nu)].y-pt[id(j-1,k+1,nu)].y)) +
                      SQR((pt[id(j+1,k-1,nu)].z-pt[id(j-1,k+1,nu)].z)) );
        }
        else
        {
          dis1 = sqrt(SQR((pt[id(j-1,k-1,nu)].x-pt[id(j  ,k  ,nu)].x)) +
                      SQR((pt[id(j-1,k-1,nu)].y-pt[id(j  ,k  ,nu)].y)) + 
                      SQR((pt[id(j-1,k-1,nu)].z-pt[id(j  ,k  ,nu)].z)) );
          dis2 = sqrt(SQR((pt[id(j  ,k-1,nu)].x-pt[id(j-1,k  ,nu)].x)) +
                      SQR((pt[id(j  ,k-1,nu)].y-pt[id(j-1,k  ,nu)].y)) + 
                      SQR((pt[id(j  ,k-1,nu)].z-pt[id(j-1,k  ,nu)].z)) );
        }

        if(dis1 <= dis2)
        {
          ident = 1;
        }
        else
        {
          ident = 2;
        }
      }
/*
*     Generate T3 triangular element
*/
      if(elemtype == 3)
      {
        if(ident == 1)
        {
          genelems[l++] = 3;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j  ,k-1,nu);

          genelems[l++] = 3;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j  ,k  ,nu);
        }
        else
        {
          genelems[l++] = 3;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j  ,k-1,nu);

          genelems[l++] = 3;
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j  ,k-1,nu);
        }
      }
/*
*     Generate Q4 quadrilateral element
*/
      else if(elemtype == 4)
      {
        genelems[l++] = 4;
        genelems[l++] = id(j-1,k-1,nu);
        genelems[l++] = id(j-1,k  ,nu);
        genelems[l++] = id(j  ,k  ,nu);
        genelems[l++] = id(j  ,k-1,nu);
      }
/*
*     Generate T6 triangular element
*/
      else if(elemtype == 6)
      {
        if(ident == 1)
        {
/*
*          Coordinates for node on inclined edge
*/
          pt[id(j,k,nu)].x = (pt[id(j-1,k-1,nu)].x+pt[id(j+1,k+1,nu)].x)/2.0;
          pt[id(j,k,nu)].y = (pt[id(j-1,k-1,nu)].y+pt[id(j+1,k+1,nu)].y)/2.0;
          pt[id(j,k,nu)].z = (pt[id(j-1,k-1,nu)].z+pt[id(j+1,k+1,nu)].z)/2.0;

          genelems[l++] = 6;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j+1,k+1,nu);
          genelems[l++] = id(j+1,k  ,nu);
          genelems[l++] = id(j+1,k-1,nu);
          genelems[l++] = id(j  ,k-1,nu);

          genelems[l++] = 6;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j-1,k+1,nu);
          genelems[l++] = id(j  ,k+1,nu);
          genelems[l++] = id(j+1,k+1,nu);
          genelems[l++] = id(j  ,k  ,nu);
        }
        else
        {
/*
*          Coordinates for node on inclined edge
*/
          pt[id(j,k,nu)].x = (pt[id(j+1,k-1,nu)].x + pt[id(j-1,k+1,nu)].x)/2.0;
          pt[id(j,k,nu)].y = (pt[id(j+1,k-1,nu)].y + pt[id(j-1,k+1,nu)].y)/2.0;
          pt[id(j,k,nu)].z = (pt[id(j+1,k-1,nu)].z + pt[id(j-1,k+1,nu)].z)/2.0;

          genelems[l++] = 6;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j-1,k+1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j+1,k-1,nu);
          genelems[l++] = id(j  ,k-1,nu);

          genelems[l++] = 6;
          genelems[l++] = id(j+1,k-1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j-1,k+1,nu);
          genelems[l++] = id(j  ,k+1,nu);
          genelems[l++] = id(j+1,k+1,nu);
          genelems[l++] = id(j+1,k  ,nu);
        }
      }
/*
*     Generate Q8 quadrilateral element
*/
      else if(elemtype == 8)
      {
        genelems[l++] = 8;
        genelems[l++] = id(j-1,k-1,nu);
        genelems[l++] = id(j-1,k  ,nu);
        genelems[l++] = id(j-1,k+1,nu);
        genelems[l++] = id(j  ,k+1,nu);
        genelems[l++] = id(j+1,k+1,nu);
        genelems[l++] = id(j+1,k  ,nu);
        genelems[l++] = id(j+1,k-1,nu);
        genelems[l++] = id(j,  k-1,nu);
      }
    }
  }
}

/*============================================================================*/
/* =============================== MshSurfLoftcoll ============================ */
/*============================================================================*/

void MshSurfLoftcoll(

double *bdynodes,  /* boundary coordinate vector data X,Y coordinates
                      of points on boundary given in clockwise order */
int    nu,         /* Numb. of nodes in the side opposite to the collap. one */
int    nv,         /* Number of nodes in the other direction */
double weight,     /* Lofting weight */
int    elemtype,   /* Element type = T3 (3), T6 (6), Q4 (4), or Q8 (8) */
int    diagtype,   /* Triangle option (option for cell diagonal)
                     = 1 --> diagonal oriented to right direction
                     = 2 --> diagonal oriented to left direction
                     = 3 --> union jack alternation
                     = 4 --> optimum diagonal (smallest of two possible) */
double *gennodes, /* Element coordinates(out) */
int    *genelems  /* Element Connectivity(out) */

)

{
/*  Local variables: */
  MshPnt3D *def_nds;
  MshPnt3D *pt;

  int    loft_seg;
  int    in1,in2,in3,in4;
  int    node_inc;
  int    ip;
  int    ident = 0;
  double dis1,dis2;
  double aa,bb,cc,u,v;
  int    i,j,k,l;


/* Processing begins here */
  def_nds = (MshPnt3D *)bdynodes;
  pt = (MshPnt3D *)gennodes;

/*
* Get node increment
*/
  if(elemtype == 6 || elemtype == 8)
  {
    node_inc = 2;
  }
  else
  {
    node_inc = 1;
  }
/*
* Get number of segments for lofting
*/
  loft_seg = (nv - 1) / node_inc;
/*
* Make sure triangular option are set only for triangular elements
*/
  if(elemtype == 4 || elemtype == 8)
  {
    diagtype = 0;
  }
/*
* Get the starting index (offset position) of each boundary on boundary array
*/
  in1 = 0;
  in2 = in1;
  in3 = in2 + (nv - 1);
  in4 = in3 + (nu - 1);
/*
* Input boundary nodes
*/
  for (i = 0;i < nu-1;i++)
  {
    ip = in1;
    pt[id(0,i,nu)].x = def_nds[ip].x;
    pt[id(0,i,nu)].y = def_nds[ip].y;
    pt[id(0,i,nu)].z = def_nds[ip].z;
  }
  for (i = 0;i < nv-1;i++)
  {
    ip = in2 + i;
    pt[id(i,nu-1,nu)].x = def_nds[ip].x;
    pt[id(i,nu-1,nu)].y = def_nds[ip].y;
    pt[id(i,nu-1,nu)].z = def_nds[ip].z;
  }
  for (i = 1;i < nu;i++)
  {
    ip = in3 + (nu-i-1);
    pt[id(nv-1,i,nu)].x = def_nds[ip].x;
    pt[id(nv-1,i,nu)].y = def_nds[ip].y;
    pt[id(nv-1,i,nu)].z = def_nds[ip].z;
  }
  for (i = 1;i < nv;i++)
  {
    ip = in4 + (nv-i-1);
    pt[id(i,0,nu)].x = def_nds[ip].x;
    pt[id(i,0,nu)].y = def_nds[ip].y;
    pt[id(i,0,nu)].z = def_nds[ip].z;
  }
/*
* Generate interior nodes
*
*
*   First compute auxiliar coefficients
*/
  if (loft_seg > 1)
  {
    aa = (2.0 * (weight)) / (((weight) + 1.0) * (double)(loft_seg));
    bb = (aa * (1.0 - (weight))) /
         (2.0 * (weight) * ((double)(loft_seg) - 1.0));

    for (k = node_inc;k < nu-node_inc;k = k+node_inc)
    {
      for (j = node_inc;j <= nv-node_inc;j = j+node_inc)
      {
        cc = (double)(j) / (double)(node_inc);
        v = aa * cc + bb * cc * (cc - 1.0);
        u = 1.0 - v;

        pt[id(j,k,nu)].x = u * pt[id(0,k,nu)].x + v * pt[id(nv-1,k,nu)].x;
        pt[id(j,k,nu)].y = u * pt[id(0,k,nu)].y + v * pt[id(nv-1,k,nu)].y;
        pt[id(j,k,nu)].z = u * pt[id(0,k,nu)].z + v * pt[id(nv-1,k,nu)].z;
      }
    }
  }
/*
* Generate midside nodes
*/
  if(node_inc == 2)
  {
    for (k = 2;k < nu-2;k = k+2)
    {
      for (j = 1;j < nv-1;j = j+2)
      {
        pt[id(j,k,nu)].x = (pt[id(j-1,k,nu)].x + pt[id(j+1,k,nu)].x) / 2.0;
        pt[id(j,k,nu)].y = (pt[id(j-1,k,nu)].y + pt[id(j+1,k,nu)].y) / 2.0;
        pt[id(j,k,nu)].z = (pt[id(j-1,k,nu)].z + pt[id(j+1,k,nu)].z) / 2.0;
      }
    }

    for (k = 1;k < nu-1;k = k+2)
    {
      for (j = 2;j < nv-2;j = j+2)
      {
        pt[id(j,k,nu)].x = (pt[id(j,k-1,nu)].x + pt[id(j,k+1,nu)].x) / 2.0;
        pt[id(j,k,nu)].y = (pt[id(j,k-1,nu)].y + pt[id(j,k+1,nu)].y) / 2.0;
        pt[id(j,k,nu)].z = (pt[id(j,k-1,nu)].z + pt[id(j,k+1,nu)].z) / 2.0;
      }
    }
  }
/*
* Get code for diagonal if we have simple cases of right or left
*/
  if(diagtype == 1)
  {
    ident = 2;
  }
  else if(diagtype == 2)
  {
    ident = 1;
  }
/*
* Start generating
*/

/* 
* Initialize the Connectivity vector
*/
  l=0;

  for (k = 1;k < nu+1-node_inc;k = k+node_inc)
  {
/*
*     The first element of this row is a triangular element close
*     to the collapsed side.
*/
    j = 1;
/*
*     Generate T3 triangular element or Q4 quadrilateral element
*/
    if((elemtype == 3) || (elemtype == 4))
    {
      genelems[l++] = 3;
      genelems[l++] = id(j-1,k-1,nu);
      genelems[l++] = id(j  ,k  ,nu);
      genelems[l++] = id(j  ,k-1,nu);
    }
/*
*     Generate T6 triangular element or Q8 quadrilateral element
*/
    else if((elemtype == 6) || (elemtype == 8))
    {
/*
*          Coordinates for node on inclined edge
*/
      if (k < nu-1)
      {
        pt[id(j,k,nu)].x = (pt[id(j-1,k-1,nu)].x + pt[id(j+1,k+1,nu)].x)/2.0;
        pt[id(j,k,nu)].y = (pt[id(j-1,k-1,nu)].y + pt[id(j+1,k+1,nu)].y)/2.0;
        pt[id(j,k,nu)].z = (pt[id(j-1,k-1,nu)].z + pt[id(j+1,k+1,nu)].z)/2.0;
      }
      else
      {
        pt[id(j,k,nu)].x = pt[id(2,nu,nu)].x;
        pt[id(j,k,nu)].y = pt[id(2,nu,nu)].y;
        pt[id(j,k,nu)].z = pt[id(2,nu,nu)].z;
      }

      genelems[l++] = 6;
      genelems[l++] = id(j-1,k-1,nu);
      genelems[l++] = id(j  ,k  ,nu);
      genelems[l++] = id(j+1,k+1,nu);
      genelems[l++] = id(j+1,k  ,nu);
      genelems[l++] = id(j+1,k-1,nu);
      genelems[l++] = id(j  ,k-1,nu);
    }
/*
*   Now the remaining elements of this row
*/
    for (j = 1+node_inc;j < nv+1-node_inc;j = j+node_inc)
    {
/*
*    Get code for diagonal if we have union jack option
*/
      if(diagtype == 3)
      {
        if( (((j+k)*(2/node_inc)) % 4) == 0)
        {
          ident = 1;
        }
        else
        {
          ident = 2;
        }
      }
/*
*    Get code for diagonal based on diagonal sizes if we have optimum option
*/
      else if(diagtype == 4)
      {
        if (node_inc == 2)
        {
          dis1 = sqrt(SQR((pt[id(j-1,k-1,nu)].x-pt[id(j+1,k+1,nu)].x)) +
                      SQR((pt[id(j-1,k-1,nu)].y-pt[id(j+1,k+1,nu)].y)) +
                      SQR((pt[id(j-1,k-1,nu)].z-pt[id(j+1,k+1,nu)].z))
                      );
          dis2 = sqrt(SQR((pt[id(j+1,k-1,nu)].x-pt[id(j-1,k+1,nu)].x)) +
                      SQR((pt[id(j+1,k-1,nu)].y-pt[id(j-1,k+1,nu)].y)) +
                      SQR((pt[id(j+1,k-1,nu)].z-pt[id(j-1,k+1,nu)].z))
                      );
        }
        else
        {
          dis1 = sqrt(SQR((pt[id(j-1,k-1,nu)].x-pt[id(j  ,k  ,nu)].x)) +
                      SQR((pt[id(j-1,k-1,nu)].y-pt[id(j  ,k  ,nu)].y)) +
                      SQR((pt[id(j-1,k-1,nu)].z-pt[id(j  ,k  ,nu)].z))
                      );
          dis2 = sqrt(SQR((pt[id(j  ,k-1,nu)].x-pt[id(j-1,k  ,nu)].x)) +
                      SQR((pt[id(j  ,k-1,nu)].y-pt[id(j-1,k  ,nu)].y)) +
                      SQR((pt[id(j  ,k-1,nu)].z-pt[id(j-1,k  ,nu)].z))
                      );
        }

        if(dis1 < dis2)
        {
          ident = 1;
        }
        else
        {
          ident = 2;
        }
      }
/*
*     Generate T3 triangular element
*/
      if(elemtype == 3)
      {
        if(ident == 1)
        {
          genelems[l++] = 3;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j  ,k-1,nu);

          genelems[l++] = 3;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j  ,k  ,nu);
        }
        else
        {
          genelems[l++] = 3;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j  ,k-1,nu);

          genelems[l++] = 3;
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j  ,k-1,nu);
        }
      }
/*
*     Generate Q4 quadrilateral element
*/
      else if(elemtype == 4)
      {
        genelems[l++] = 4;
        genelems[l++] = id(j-1,k-1,nu);
        genelems[l++] = id(j-1,k  ,nu);
        genelems[l++] = id(j  ,k  ,nu);
        genelems[l++] = id(j  ,k-1,nu);
      }
/*
*     Generate T6 triangular element
*/
      else if(elemtype == 6)
      {
        if(ident == 1)
        {
/*
*          Coordinates for node on inclined edge
*/
          pt[id(j,k,nu)].x = (pt[id(j-1,k-1,nu)].x + pt[id(j+1,k+1,nu)].x)/2.0;
          pt[id(j,k,nu)].y = (pt[id(j-1,k-1,nu)].y + pt[id(j+1,k+1,nu)].y)/2.0;
          pt[id(j,k,nu)].z = (pt[id(j-1,k-1,nu)].z + pt[id(j+1,k+1,nu)].z)/2.0;

          genelems[l++] = 6;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j+1,k+1,nu);
          genelems[l++] = id(j+1,k  ,nu);
          genelems[l++] = id(j+1,k-1,nu);
          genelems[l++] = id(j  ,k-1,nu);

          genelems[l++] = 6;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j-1,k+1,nu);
          genelems[l++] = id(j  ,k+1,nu);
          genelems[l++] = id(j+1,k+1,nu);
          genelems[l++] = id(j  ,k  ,nu);
        }
        else
        {
/*
*          Coordinates for node on inclined edge
*/
          pt[id(j,k,nu)].x = (pt[id(j+1,k-1,nu)].x + pt[id(j-1,k+1,nu)].x)/2.0;
          pt[id(j,k,nu)].y = (pt[id(j+1,k-1,nu)].y + pt[id(j-1,k+1,nu)].y)/2.0;
          pt[id(j,k,nu)].z = (pt[id(j+1,k-1,nu)].z + pt[id(j-1,k+1,nu)].z)/2.0;

          genelems[l++] = 6;
          genelems[l++] = id(j-1,k-1,nu);
          genelems[l++] = id(j-1,k  ,nu);
          genelems[l++] = id(j-1,k+1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j+1,k-1,nu);
          genelems[l++] = id(j  ,k-1,nu);

          genelems[l++] = 6;
          genelems[l++] = id(j+1,k-1,nu);
          genelems[l++] = id(j  ,k  ,nu);
          genelems[l++] = id(j-1,k+1,nu);
          genelems[l++] = id(j  ,k+1,nu);
          genelems[l++] = id(j+1,k+1,nu);
          genelems[l++] = id(j+1,k  ,nu);
        }
      }
/*
*     Generate Q8 quadrilateral element
*/
      else if(elemtype == 8)
      {
        genelems[l++] = 8;
        genelems[l++] = id(j-1,k-1,nu);
        genelems[l++] = id(j-1,k  ,nu);
        genelems[l++] = id(j-1,k+1,nu);
        genelems[l++] = id(j  ,k+1,nu);
        genelems[l++] = id(j+1,k+1,nu);
        genelems[l++] = id(j+1,k  ,nu);
        genelems[l++] = id(j+1,k-1,nu);
        genelems[l++] = id(j,  k-1,nu);
      }
    }
  }
}

/*============================================================================*/
/* ================================ MshSurfTrimap ============================= */
/*============================================================================*/

void MshSurfTrimap(

 double   *bdynodes,  /* boundary coordinate vector data X,Y coordinates
                         of points on boundary given in clockwise order */
 int      n,          /* Number of nodes on one boundary     */
 int      elemtype,   /* Element type= T3 (3), T6 (6), Q4 (4), or Q8 (8)  */
 double   *gennodes,  /* Element coordinates (out) */
 int      *genelems   /* Element Connectivity (out) */

)

{

/*  Local variables: */
  MshPnt3D *def_nds;
  MshPnt3D *pt;

  int    in1,in2,in3;
  int    node_inc;
  int    ip;
  int    j,k,l;


/* Processing begins here */
  def_nds = (MshPnt3D *)bdynodes;
  pt = (MshPnt3D *)gennodes;

/*
* Get node increment
*/
  if ((elemtype==6) || (elemtype==8))
  {
    node_inc = 2;
  }
  else
  {
    node_inc = 1;
  }

/*
* Get the starting index (offset position) of each boundary on boundary array
*/
  in1 = 0;
  in2 = in1 + (n - 1);
  in3 = in2 + (n - 1);
/*
* Input coordinates of nodes on first (K) boundary
*/
  for (k = 0;k < n;k++)
  {
    ip = in1 + k;
    pt[id(0,k,n)].x = def_nds[ip].x;
    pt[id(0,k,n)].y = def_nds[ip].y;
    pt[id(0,k,n)].z = def_nds[ip].z;
  }
/*
* Input coordinates of nodes on second (J) boundary
*/
  for (j = 1;j < n;j++)
  {
    ip = in3 + (n-j-1);
    pt[id(j,0,n)].x = def_nds[ip].x;
    pt[id(j,0,n)].y = def_nds[ip].y;
    pt[id(j,0,n)].z = def_nds[ip].z;
  }

/*
* Input coordinates of nodes on third boundary
*/
  for (k = 1; k < n-1;k++)
  {
    ip = in2 + (n-k-1);
    pt[id(n-k-1,k,n)].x = def_nds[ip].x;
    pt[id(n-k-1,k,n)].y = def_nds[ip].y;
    pt[id(n-k-1,k,n)].z = def_nds[ip].z;
  }

/*
* Generate interior nodes
*/
  for (k = node_inc; k < n-node_inc-1;k = k+node_inc)
  {
    for (j = node_inc;j < n-k-node_inc+1;j = j+node_inc)
    {
      pt[id(j,k,n)].x =
        (pt[id(j,n-j-1,n)].x*k/(n-j-1)+
         pt[id(j,0,n)].x*(n-j-k-1)/(n-j-1)+
         pt[id(j+k,0,n)].x*j/(j+k)+
         pt[id(0,j+k,n)].x*k/(j+k)+
         pt[id(0,k,n)].x*(n-j-k-1)/(n-k-1)+
         pt[id(n-k-1,k,n)].x*j/(n-k-1)-
         pt[id(0,0,n)].x*(n-j-k-1)/(n-1)-
         pt[id(0,n-1,n)].x*k/(n-1)-
         pt[id(n-1,0,n)].x*j/(n-1))*0.5;

      pt[id(j,k,n)].y =
        (pt[id(j,n-j-1,n)].y*k/(n-j-1)+
         pt[id(j,0,n)].y*(n-j-k-1)/(n-j-1)+
         pt[id(j+k,0,n)].y*j/(j+k)+
         pt[id(0,j+k,n)].y*k/(j+k)+
         pt[id(0,k,n)].y*(n-j-k-1)/(n-k-1)+
         pt[id(n-k-1,k,n)].y*j/(n-k-1)-
         pt[id(0,0,n)].y*(n-j-k-1)/(n-1)-
         pt[id(0,n-1,n)].y*k/(n-1)-
         pt[id(n-1,0,n)].y*j/(n-1))*0.5;

      pt[id(j,k,n)].z =
        (pt[id(j,n-j-1,n)].z*k/(n-j-1)+
         pt[id(j,0,n)].z*(n-j-k-1)/(n-j-1)+
         pt[id(j+k,0,n)].z*j/(j+k)+
         pt[id(0,j+k,n)].z*k/(j+k)+
         pt[id(0,k,n)].z*(n-j-k-1)/(n-k-1)+
         pt[id(n-k-1,k,n)].z*j/(n-k-1)-
         pt[id(0,0,n)].z*(n-j-k-1)/(n-1)-
         pt[id(0,n-1,n)].z*k/(n-1)-
         pt[id(n-1,0,n)].z*j/(n-1))*0.5;
    }
  }

  if(node_inc==2)
  {
/*
* Generate midside nodes in J direction
*/
    for (k = 2;k < n-node_inc;k = k+node_inc)
    {
      for (j = node_inc-1;j < n-k;j = j+node_inc)
      {
        pt[id(j,k,n)].x = (pt[id(j-1,k,n)].x + pt[id(j+1,k,n)].x)/2.0;
        pt[id(j,k,n)].y = (pt[id(j-1,k,n)].y + pt[id(j+1,k,n)].y)/2.0;
        pt[id(j,k,n)].z = (pt[id(j-1,k,n)].z + pt[id(j+1,k,n)].z)/2.0;
      }
    }

/*
* Generate midside nodes in K direction
*/
    for (k = 1;k < n-1;k = k+node_inc)
    {
      for (j = 2;j < n-k;j = j+node_inc)
      {
        pt[id(j,k,n)].x = (pt[id(j,k-1,n)].x + pt[id(j,k+1,n)].x)/2.0;
        pt[id(j,k,n)].y = (pt[id(j,k-1,n)].y + pt[id(j,k+1,n)].y)/2.0;
        pt[id(j,k,n)].z = (pt[id(j,k-1,n)].z + pt[id(j,k+1,n)].z)/2.0;
      }
    }

/*
* Generate midside nodes in the third direction
*/
    for (k = 1;k < n-1;k = k+node_inc)
    {
      for (j = 1;j < n-k;j = j+node_inc)
      {
        pt[id(j,k,n)].x = (pt[id(j-1,k+1,n)].x + pt[id(j+1,k-1,n)].x)/2.0;
        pt[id(j,k,n)].y = (pt[id(j-1,k+1,n)].y + pt[id(j+1,k-1,n)].y)/2.0;
        pt[id(j,k,n)].z = (pt[id(j-1,k+1,n)].z + pt[id(j+1,k-1,n)].z)/2.0;
      }
    }
  }
/*
* Generate elements
*/

/* 
* Initialize the Connectivity vector
*/
  l=0;

  for (k = 0;k < n-node_inc;k = k+node_inc)
  {
    for (j = 0;j < n-k-node_inc-1;j = j+node_inc)
    {
/*
*     Generate T3 triangular element
*/
      if (elemtype == 3)
      {
        genelems[l++] = 3;
        genelems[l++] = id(j  ,k  ,n);
        genelems[l++] = id(j  ,k+1,n);
        genelems[l++] = id(j+1,k  ,n);

        genelems[l++] = 3;
        genelems[l++] = id(j+1,k  ,n);
        genelems[l++] = id(j  ,k+1,n);
        genelems[l++] = id(j+1,k+1,n);
      }
/*
*     Generate Q4 quadrilateral element
*/
      else if (elemtype == 4)
      {
        genelems[l++] = 4;
        genelems[l++] = id(j  ,k  ,n);
        genelems[l++] = id(j  ,k+1,n);
        genelems[l++] = id(j+1,k+1,n);
        genelems[l++] = id(j+1,k  ,n);
      }
/*
*     Generate T6 triangular element
*/
      else if (elemtype == 6)
      {
        genelems[l++] = 6;
        genelems[l++] = id(j  ,k  ,n);
        genelems[l++] = id(j  ,k+1,n);
        genelems[l++] = id(j  ,k+2,n);
        genelems[l++] = id(j+1,k+1,n);
        genelems[l++] = id(j+2,k  ,n);
        genelems[l++] = id(j+1,k  ,n);

        genelems[l++] = 6;
        genelems[l++] = id(j+2,k  ,n);
        genelems[l++] = id(j+1,k+1,n);
        genelems[l++] = id(j  ,k+2,n);
        genelems[l++] = id(j+1,k+2,n);
        genelems[l++] = id(j+2,k+2,n);
        genelems[l++] = id(j+2,k+1,n);
      }
/*
*     Generate Q8 quadrilateral element
*/
      else if (elemtype == 8)
      {
        genelems[l++] = 8;
        genelems[l++] = id(j  ,k  ,n);
        genelems[l++] = id(j  ,k+1,n);
        genelems[l++] = id(j  ,k+2,n);
        genelems[l++] = id(j+1,k+2,n);
        genelems[l++] = id(j+2,k+2,n);
        genelems[l++] = id(j+2,k+1,n);
        genelems[l++] = id(j+2,k  ,n);
        genelems[l++] = id(j+1,k  ,n);
      }
    }
/*
*     Generate triangular elements to complete the triangular patch
*/
    if(node_inc==1)
    {
      genelems[l++] = 3;
      genelems[l++] = id(n-k-2,k  ,n);
      genelems[l++] = id(n-k-2,k+1,n);
      genelems[l++] = id(n-k-1,k  ,n);
    }
    else
    {
      genelems[l++] = 6;
      genelems[l++] = id(n-k-3,k  ,n);
      genelems[l++] = id(n-k-3,k+1,n);
      genelems[l++] = id(n-k-3,k+2,n);
      genelems[l++] = id(n-k-2,k+1,n);
      genelems[l++] = id(n-k-1,k  ,n);
      genelems[l++] = id(n-k-2,k  ,n);
    }
  }
}





/*============================================================================*/
/* ===================== MshSurfAutomaticBilinearCornes ===================== */
/*============================================================================*/
int MshSurfAutomaticBilinearCornes (int np, double *pts, int *m, int *n)
{
  int nbordes, *borders, i, status, pos;
  double *tmp_pts;
	
  /* get border stones */
  MshSurfFindBorders(np, pts, &nbordes, &borders );

  if (nbordes > 6)
    nbordes = 6;

  /* order border stones id sequentially to loop */ 
  qsort (borders, nbordes, sizeof( int ), compInt );


  /* number of elements in each direction */
  switch(nbordes)
  {
  case 4:
    status = MshSurfGetMappDirection4_ (borders, np, m, n, &pos);
    break;

  case 5:
    status = MshSurfGetMappDirection5_ (borders, np, m, n, &pos);
    break;

  case 6:
    status = MshSurfGetMappDirection6_ (borders, np, m, n, &pos);
    break;


  default:
    return 0;
    break;
  }

  if (!status)
    return 0;


  /* align with first border stone */
  tmp_pts = (double *) calloc (np*3, sizeof (double));
  memcpy(tmp_pts, pts, np*3*sizeof (double));
  for (i = 0; i < np; ++i)
  {
    pts[i*3+0] = tmp_pts[((pos+i)%np)*3+0];
    pts[i*3+1] = tmp_pts[((pos+i)%np)*3+1];
    pts[i*3+2] = tmp_pts[((pos+i)%np)*3+2];
  }

  free (tmp_pts);
  free (borders);

  return 1;
}


/*============================================================================*/
/* ==================== MshSurfAutomaticTrilinearCornes ===================== */
/*============================================================================*/
int MshSurfAutomaticTrilinearCornes (int np, double *pts, int *n)
{
  int nbordes, *borders, i, status, pos;
  double *tmp_pts;

  /* get border stones */
  MshSurfFindBorders(np, pts, &nbordes, &borders );

  if (nbordes > 5)
    nbordes = 5;

  /* order border stones id sequentially to loop */
  qsort (borders, nbordes, sizeof( int ), compInt );


  /* number of elements in each direction */
  switch(nbordes)
  {
  case 3:
    status = MshSurfGetMappTriDirection3_(borders, np, n, &pos);
    break;

  case 4:
    status = MshSurfGetMappTriDirection4_ (borders, np, n, &pos);
    break;

  case 5:
    status = MshSurfGetMappTriDirection5_ (borders, np, n, &pos);
    break;

  default:
    return 0;
    break;
  }

  if (!status)
    return 0;


  /* align with first border stone */
  tmp_pts = (double *) calloc (np*3, sizeof (double));
  memcpy(tmp_pts, pts, np*3*sizeof (double));
  for (i = 0; i < np; ++i)
  {
    pts[i*3+0] = tmp_pts[((pos+i)%np)*3+0];
    pts[i*3+1] = tmp_pts[((pos+i)%np)*3+1];
    pts[i*3+2] = tmp_pts[((pos+i)%np)*3+2];
  }

  free (tmp_pts);
  free (borders);

  return 1;
}



/* --------------------------------------------------------------- */
static int compInt( const void *c1, const void *c2 )
{
  int *m1 = (int *) c1;
  int *m2 = (int *) c2;
  if      (*m1 < *m2) return -1;
  else if (*m1 > *m2) return +1;
  else                return  0;
}



/* --------------------------------------------------------------- */
void MshSurfAlignWithBS( int nNodes, int *idNodes, int nbordes, int *borders )
{
  int i, j;
  int *tmpnodes, reference;
  /* if nodes are aligned with the initial border stone */
  if (borders[0] != 0)
  {
    /* copy idNodes to tmpnodes */
    tmpnodes = (int *) calloc(nNodes, sizeof(int));
    memcpy (tmpnodes, idNodes, nNodes * sizeof(int));

    /* fill idNodes aligned with borders[0] */
    reference = borders[0];
    for (i = 0; i<nNodes; ++i)
    {
      j = (i+reference)%nNodes;
      idNodes[i] = tmpnodes[j];   
    }

    /* correct id border stones */
    for (i = 0; i<nbordes; ++i)
    {
      int new_id = borders[i] - reference;
      if (new_id < 0)
        borders[i] = new_id + nNodes;
      else
        borders[i] = new_id;
    }

    free(tmpnodes);
  }

  /* order border stones id sequentially to loop */
  qsort (borders, nbordes, sizeof( int ), compInt );

}


/* --------------------------------------------------------------- */
int MshSurfGetMappTriDirection3_( int *bd, int nPtsLoop, int *n, int *pos)
{
  int b[3] = {bd[1]-bd[0], bd[2]-bd[1], nPtsLoop-bd[2]+bd[0]};

  if (b[0] == b[1] && b[1] == b[2]) /* match */
  {
    *n = b[0] + 1;
    *pos = 0;
    return 1;
  }
  return 0;
}


/* --------------------------------------------------------------- */
int MshSurfGetMappTriDirection4_( int *bd, int nPtsLoop, int *n, int *pos)
{
  int b[4] = {bd[1]-bd[0], bd[2]-bd[1], bd[3]-bd[2], nPtsLoop-bd[3]+bd[0]};
  int i, div[3];

  for (i = 0; i < 4; ++i)
  {
    div[0] = b[i] + b[(i+1)%4];
    div[1] = b[(i+2)%4];
    div[2] = b[(i+3)%4];
    if (div[0] == div[1] && div[1] == div[2]) /* match */
    {
      *n   = div[0] + 1;
      *pos = i;
      return 1;
    }
  }

  return 0;
}


/* --------------------------------------------------------------- */
int MshSurfGetMappTriDirection5_( int *bd, int nPtsLoop, int *n, int *pos)
{
  int i, div[3];
  int b[5] = {bd[1]-bd[0], bd[2]-bd[1], bd[3]-bd[2], bd[4]-bd[3], nPtsLoop-bd[4]+bd[0]};

  /* first situation */
  for (i = 0; i < 5; ++i)
  {
    div[0] = b[i]       + b[(i+1)%5];
    div[1] = b[(i+2)%5] + b[(i+3)%5];
    div[2] = b[(i+4)%5];
    if (div[0] == div[1] && div[1] == div[2])  /* match */
    {
      *n = div[0] + 1;
      *pos = i;
      return 1;
    }
  }

  return 0;
}

/* --------------------------------------------------------------- */
int MshSurfGetMappDirection4_( int *bd, int nPtsLoop, int *m, int *n, int *pos)
{
  int b[4] = {bd[1]-bd[0], bd[2]-bd[1], bd[3]-bd[2], nPtsLoop-bd[3]+bd[0]};

  if (b[0] == b[2] && b[1] == b[3])  /* match */
  {
    *m = b[0] + 1;
    *n = b[1] + 1;
    *pos = 0;
    return 1;
  }
  return 0;
}

/* --------------------------------------------------------------- */
int MshSurfGetMappDirection5_( int *bd, int nPtsLoop, int *m, int *n, int *pos)
{
  int i, div[4];
  int b[5] = {bd[1]-bd[0], bd[2]-bd[1], bd[3]-bd[2], bd[4]-bd[3], nPtsLoop-bd[4]+bd[0]};

  for (i = 0; i < 5; ++i)
  {
    div[0] = b[i] + b[(i+1)%5];
    div[1] = b[(i+2)%5];
    div[2] = b[(i+3)%5];
    div[3] = b[(i+4)%5];
    if (div[0] == div[2] && div[1] == div[3])  /* match */
    {
      *m = div[i%2+0] + 1;
      *n = div[i%2+1] + 1;
      *pos = i;
      return 1;
    }
  }
  return 0;
}

/* --------------------------------------------------------------- */
int MshSurfGetMappDirection6_( int *bd, int nPtsLoop, int *m, int *n, int *pos)
{
  int i, div[4];
  int b[6] = {bd[1]-bd[0], bd[2]-bd[1], bd[3]-bd[2], bd[4]-bd[3], bd[5]-bd[4], nPtsLoop-bd[5]+bd[0]};
  int div1[5] = {0, 1, 0, -1, 1};

  /* first situation */
  for (i = 0; i < 5; ++i)
  {
    if (div1[i] == -1)
      continue;
    div[0] = b[i]       + b[(i+1)%6];
    div[1] = b[(i+2)%6] + b[(i+3)%6];
    div[2] = b[(i+4)%6];
    div[3] = b[(i+5)%6];
    if (div[0] == div[2] && div[1] == div[3])  /* match */
    {
      int mm = div1[i];
      int nn = !div1[i];
      *m = div[mm] + 1;
      *n = div[nn] + 1;
      *pos = i;
      return 1;
    }
  }

  /* second situation */
  for (i = 0; i < 6; ++i)
  {
    div[0] = b[i]       + b[(i+1)%6];
    div[1] = b[(i+2)%6];
    div[2] = b[(i+3)%6] + b[(i+4)%6];
    div[3] = b[(i+5)%6];
    if (div[0] == div[2] && div[1] == div[3])  /* match */
    {
      *m = div[i%2+0] + 1;
      *n = div[i%2+1] + 1;
      *pos = i;
      return 1;
    }
  }

  return 0;
}

