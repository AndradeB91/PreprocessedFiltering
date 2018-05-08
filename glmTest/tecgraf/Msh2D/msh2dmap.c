/*
**	mshmap2d - package for creating meshes based on transfinite mapping.
**		   It will be used the following meshing techniques:
**
**		 - Msh2DTransfin:
**		     Bilinear mapping for quadrilateral regions with sides 
**		     of any shape, i.e., the boundary may be formed by 
**		     lines, circles, beziers, and so on.  The number of 
**		     segments on opposite sides must be equal.
**		 - Msh2DLofting:
**		     Linear mapping between two opposite sides of a 
**		     quadrilateral region.  The other two sides must 
**		     be straight lines.  The number of segments on opposite
**		     sides must be equal.
**		 - Msh2DTrsfncoll
**		     Bilinear mapping for triangular regions which conceptually
**		     is formed by a quadrilateral region with one of the sides
**		     collapsed to a point.  The 3 sides can be of any shape and
**		     the two sides adjacent to the collapsed side must have
**		     the same number of segments.
**		 - Msh2DLoftcoll:
**		     Linear mapping between a point and a curve, which form
**		     a triangular region.  The two sides of the region which
**		     are adjacent to the point of lofting must be straight 
**		     lines and have the same number of segments.
**		 - Msh2DTrimap:
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
**	void Msh2DTransfin(
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
**	Msh2DTransfin ( bilinear mapping meshing with any boundary shape ).
**	Based on the given points coordinates on the boundary , the number of
**	nodes in the two directions, and the element type, the routine fills
**	an array of element nodes coordinates, and a structure composed by
**  the number of element nodes and the connectivity.
**	If a triangular element is given, the diagonal orientation is also
**	required (diagtype).
**
** ---------------------------------------------------------------------------
**
**	void Msh2DLofting(
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
**	Msh2DLofting ( bilinear mapping meshing with straigth boundary )
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
**	void Msh2DTrsfncoll(
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
**	Msh2DTrsfncoll ( bilinear collapsed mapping meshing with any
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
**	void Msh2DLoftcoll(
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
**	Msh2DLoftcoll ( bilinear collapsed mapping meshing with straigth
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
**	void Msh2DTrimap(
**
**	double   *bdynodes,    boundary coordinate vector data X,Y coordinates
**	                       of points on boundary given in clockwise order
**	int      n,            Number of nodes on one boundary
**	int      elemtype,     Element type= T3 (3), T6 (6), Q4 (4), or Q8 (8) 
**	double   *gennodes,    Element coordinates          (out)
**	int  	 *genelems     Element Connectivity         (out)
**
**	Msh2DTrimap ( trilinear mapping meshing with any boundary shape ).
**	Based on the given points coordinates on the boundary , the number of
**	nodes in the three directions, and the element type, the routine fills
**	an array of element nodes coordinates, and a structure composed by
**  the number of element nodes and the connectivity.
**
** ---------------------------------------------------------------------------
*/
	
/* Global variables and symbols           */
#include 	<math.h>
#include 	"msh2dmap.h" 

typedef struct _mshpnt
{
    double x, y;
} MshPnt ;

/* Local variables and symbols            */

#define  SQR(a) ((a)*(a))

#define id(i,j,m) (((i)*(m))+(j))


/*=======================================================================*/
/*============================ Msh2DTransfin ============================*/
/*=======================================================================*/


void Msh2DTransfin(

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

  MshPnt *def_nds;
  MshPnt *pt;

  int     in1,in2,in3,in4;
  int     node_inc;
  int     ip;
  int     ident = 0;
  double   dis1,dis2;
  int     i,j,k,l;

/* Processing begins here */

  def_nds = (MshPnt *)bdynodes;
  pt = (MshPnt *)gennodes;

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
  }
  for (i = 0;i < (nv-1);i++)
  {
    ip = in2 + i;
    pt[id(i,nu-1,nu)].x = def_nds[ip].x;
    pt[id(i,nu-1,nu)].y = def_nds[ip].y;
  }
  for (i = 1;i < nu;i++)
  {
    ip = in3 + (nu-i-1);
    pt[id(nv-1,i,nu)].x = def_nds[ip].x;
    pt[id(nv-1,i,nu)].y = def_nds[ip].y;
  }
  for (i = 1;i < nv;i++)
  {
    ip = in4 + (nv-i-1);
    pt[id(i,0,nu)].x = def_nds[ip].x;
    pt[id(i,0,nu)].y = def_nds[ip].y;
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
                  SQR((pt[id(j-1,k-1,nu)].y-pt[id(j+1,k+1,nu)].y)));
           dis2 = sqrt(SQR((pt[id(j+1,k-1,nu)].x-pt[id(j-1,k+1,nu)].x)) +
                  SQR((pt[id(j+1,k-1,nu)].y-pt[id(j-1,k+1,nu)].y)));
         }
         else
         {
           dis1 = sqrt(SQR((pt[id(j-1,k-1,nu)].x-pt[id(j  ,k  ,nu)].x)) +
                  SQR((pt[id(j-1,k-1,nu)].y-pt[id(j  ,k  ,nu)].y)));
           dis2 = sqrt(SQR((pt[id(j  ,k-1,nu)].x-pt[id(j-1,k  ,nu)].x)) +
                  SQR((pt[id(j  ,k-1,nu)].y-pt[id(j-1,k  ,nu)].y)));
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
/*============================== Msh2DLofting ================================*/
/*============================================================================*/

void Msh2DLofting(

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
  MshPnt *def_nds;
  MshPnt *pt;

  int      loft_seg;
  int      in1,in2,in3,in4;
  int      node_inc;
  int      ip;
  int      ident = 0;
  double   dis1,dis2;
  double   aa,bb,cc,u,v;
  int      i,j,k,l;


/* Processing begins here */

  def_nds = (MshPnt *)bdynodes;
  pt = (MshPnt *)gennodes;


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
  }
  for (i = 0;i < nv-1;i++)
  {
    ip = in2 + i;
    pt[id(i,nu-1,nu)].x = def_nds[ip].x;
    pt[id(i,nu-1,nu)].y = def_nds[ip].y;
  }
  for (i = 1;i < nu;i++) 
  {
    ip = in3 + (nu-i-1);
    pt[id(nv-1,i,nu)].x = def_nds[ip].x;
    pt[id(nv-1,i,nu)].y = def_nds[ip].y;
  }
  for (i= 1;i < nv;i++) 
  {
    ip = in4 + (nv-i-1);
    pt[id(i,0,nu)].x = def_nds[ip].x;
    pt[id(i,0,nu)].y = def_nds[ip].y;
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

          pt[id(j,k,nu)].x =
              u * pt[id(j,0,nu)].x  +  v * pt[id(j,nu-1,nu)].x;
          pt[id(j,k,nu)].y =
              u * pt[id(j,0,nu)].y  +  v * pt[id(j,nu-1,nu)].y;

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

           pt[id(j,k,nu)].x =
		     u * pt[id(0,k,nu)].x  +  v * pt[id(nv-1,k,nu)].x;
           pt[id(j,k,nu)].y =
		     u * pt[id(0,k,nu)].y  +  v * pt[id(nv-1,k,nu)].y;

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
         pt[id(j,k,nu)].x =
           (pt[id(j-1,k,nu)].x + pt[id(j+1,k,nu)].x) / 2.0;
         pt[id(j,k,nu)].y =
           (pt[id(j-1,k,nu)].y + pt[id(j+1,k,nu)].y) / 2.0;
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
                 SQR((pt[id(j-1,k-1,nu)].y-pt[id(j+1,k+1,nu)].y)));
           dis2 = 
              sqrt(SQR((pt[id(j+1,k-1,nu)].x-pt[id(j-1,k+1,nu)].x)) +
                   SQR((pt[id(j+1,k-1,nu)].y-pt[id(j-1,k+1,nu)].y)));
        }
        else
        {
          dis1 = 
            sqrt(SQR((pt[id(j-1,k-1,nu)].x-pt[id(j  ,k  ,nu)].x)) +
                 SQR((pt[id(j-1,k-1,nu)].y-pt[id(j  ,k  ,nu)].y)));
          dis2 =
            sqrt(SQR((pt[id(j  ,k-1,nu)].x-pt[id(j-1,k  ,nu)].x)) +
                 SQR((pt[id(j  ,k-1,nu)].y-pt[id(j-1,k  ,nu)].y)));
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
          pt[id(j,k,nu)].x =
            (pt[id(j-1,k-1,nu)].x + pt[id(j+1,k+1,nu)].x)/2.0;
          pt[id(j,k,nu)].y =
            (pt[id(j-1,k-1,nu)].y + pt[id(j+1,k+1,nu)].y)/2.0;

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

void Msh2DTrsfncoll(

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
  MshPnt *def_nds;
  MshPnt *pt;

  int     in1,in2,in3,in4;
  int     node_inc;
  int     ip;
  int     ident = 0;
  double  dis1,dis2;
  int     i,j,k,l;


/* Processing begins here */
  def_nds = (MshPnt *)bdynodes;
  pt = (MshPnt *)gennodes;

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
  }
  for (i = 0;i < nv-1;i++)
  {
    ip = in2 + i;
    pt[id(i,nu-1,nu)].x = def_nds[ip].x;
    pt[id(i,nu-1,nu)].y = def_nds[ip].y;
  }
  for (i = 1;i < nu;i++)
  {
    ip = in3 + (nu-i-1);
    pt[id(nv-1,i,nu)].x = def_nds[ip].x;
    pt[id(nv-1,i,nu)].y = def_nds[ip].y;
  }
  for (i = 1;i < nv;i++)
  {
    ip = in4 + (nv-i-1);
    pt[id(i,0,nu)].x = def_nds[ip].x;
    pt[id(i,0,nu)].y = def_nds[ip].y;
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
        pt[id(j,k,nu)].x =
          (pt[id(j-1,k,nu)].x + pt[id(j+1,k,nu)].x) / 2.0;
        pt[id(j,k,nu)].y =
          (pt[id(j-1,k,nu)].y + pt[id(j+1,k,nu)].y) / 2.0;
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
        pt[id(j,k,nu)].x =
          (pt[id(j-1,k-1,nu)].x + pt[id(j+1,k+1,nu)].x)/2.0;
        pt[id(j,k,nu)].y =
          (pt[id(j-1,k-1,nu)].y + pt[id(j+1,k+1,nu)].y)/2.0;
      }
      else
      {
        pt[id(j,k,nu)].x = pt[id(2,nu,nu)].x;
        pt[id(j,k,nu)].y = pt[id(2,nu,nu)].y;
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
                      SQR((pt[id(j-1,k-1,nu)].y-pt[id(j+1,k+1,nu)].y)));
          dis2 = sqrt(SQR((pt[id(j+1,k-1,nu)].x-pt[id(j-1,k+1,nu)].x)) +
                      SQR((pt[id(j+1,k-1,nu)].y-pt[id(j-1,k+1,nu)].y)));
        }
        else
        {
          dis1 = sqrt(SQR((pt[id(j-1,k-1,nu)].x-pt[id(j  ,k  ,nu)].x)) +
                      SQR((pt[id(j-1,k-1,nu)].y-pt[id(j  ,k  ,nu)].y)));
          dis2 = sqrt(SQR((pt[id(j  ,k-1,nu)].x-pt[id(j-1,k  ,nu)].x)) +
                      SQR((pt[id(j  ,k-1,nu)].y-pt[id(j-1,k  ,nu)].y)));
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
          pt[id(j,k,nu)].x =
            (pt[id(j-1,k-1,nu)].x+pt[id(j+1,k+1,nu)].x)/2.0;
          pt[id(j,k,nu)].y =
            (pt[id(j-1,k-1,nu)].y+pt[id(j+1,k+1,nu)].y)/2.0;

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
/* =============================== Msh2DLoftcoll ============================ */
/*============================================================================*/

void Msh2DLoftcoll(

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
  MshPnt *def_nds;
  MshPnt *pt;

  int    loft_seg;
  int    in1,in2,in3,in4;
  int    node_inc;
  int    ip;
  int    ident = 0;
  double dis1,dis2;
  double aa,bb,cc,u,v;
  int    i,j,k,l;


/* Processing begins here */
  def_nds = (MshPnt *)bdynodes;
  pt = (MshPnt *)gennodes;

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
  }
  for (i = 0;i < nv-1;i++)
  {
    ip = in2 + i;
    pt[id(i,nu-1,nu)].x = def_nds[ip].x;
    pt[id(i,nu-1,nu)].y = def_nds[ip].y;
  }
  for (i = 1;i < nu;i++)
  {
    ip = in3 + (nu-i-1);
    pt[id(nv-1,i,nu)].x = def_nds[ip].x;
    pt[id(nv-1,i,nu)].y = def_nds[ip].y;
  }
  for (i = 1;i < nv;i++)
  {
    ip = in4 + (nv-i-1);
    pt[id(i,0,nu)].x = def_nds[ip].x;
    pt[id(i,0,nu)].y = def_nds[ip].y;
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
        pt[id(j,k,nu)].x =
          (pt[id(j-1,k,nu)].x + pt[id(j+1,k,nu)].x) / 2.0;
        pt[id(j,k,nu)].y =
          (pt[id(j-1,k,nu)].y + pt[id(j+1,k,nu)].y) / 2.0;
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
        pt[id(j,k,nu)].x =
          (pt[id(j-1,k-1,nu)].x + pt[id(j+1,k+1,nu)].x)/2.0;
        pt[id(j,k,nu)].y =
          (pt[id(j-1,k-1,nu)].y + pt[id(j+1,k+1,nu)].y)/2.0;
      }
      else
      {
        pt[id(j,k,nu)].x = pt[id(2,nu,nu)].x;
        pt[id(j,k,nu)].y = pt[id(2,nu,nu)].y;
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
                      SQR((pt[id(j-1,k-1,nu)].y-pt[id(j+1,k+1,nu)].y)));
          dis2 = sqrt(SQR((pt[id(j+1,k-1,nu)].x-pt[id(j-1,k+1,nu)].x)) +
                      SQR((pt[id(j+1,k-1,nu)].y-pt[id(j-1,k+1,nu)].y)));
        }
        else
        {
          dis1 = sqrt(SQR((pt[id(j-1,k-1,nu)].x-pt[id(j  ,k  ,nu)].x)) +
                      SQR((pt[id(j-1,k-1,nu)].y-pt[id(j  ,k  ,nu)].y)));
          dis2 = sqrt(SQR((pt[id(j  ,k-1,nu)].x-pt[id(j-1,k  ,nu)].x)) +
                      SQR((pt[id(j  ,k-1,nu)].y-pt[id(j-1,k  ,nu)].y)));
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
          pt[id(j,k,nu)].x =
            (pt[id(j-1,k-1,nu)].x + pt[id(j+1,k+1,nu)].x)/2.0;
          pt[id(j,k,nu)].y =
            (pt[id(j-1,k-1,nu)].y + pt[id(j+1,k+1,nu)].y)/2.0;

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
/* ================================ Msh2DTrimap ============================= */
/*============================================================================*/

void Msh2DTrimap(

 double   *bdynodes,  /* boundary coordinate vector data X,Y coordinates
                         of points on boundary given in clockwise order */
 int      n,          /* Number of nodes on one boundary     */
 int      elemtype,   /* Element type= T3 (3), T6 (6), Q4 (4), or Q8 (8)  */
 double   *gennodes,  /* Element coordinates (out) */
 int      *genelems   /* Element Connectivity (out) */

)

{

/*  Local variables: */
  MshPnt *def_nds;
  MshPnt *pt;

  int    in1,in2,in3;
  int    node_inc;
  int    ip;
  int    j,k,l;


/* Processing begins here */
  def_nds = (MshPnt *)bdynodes;
  pt = (MshPnt *)gennodes;

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
  }
/*
* Input coordinates of nodes on second (J) boundary
*/
  for (j = 1;j < n;j++)
  {
    ip = in3 + (n-j-1);
    pt[id(j,0,n)].x = def_nds[ip].x;
    pt[id(j,0,n)].y = def_nds[ip].y;
  }

/*
* Input coordinates of nodes on third boundary
*/
  for (k = 1; k < n-1;k++)
  {
    ip = in2 + (n-k-1);
    pt[id(n-k-1,k,n)].x = def_nds[ip].x;
    pt[id(n-k-1,k,n)].y = def_nds[ip].y;
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
        pt[id(j,k,n)].x =
          (pt[id(j-1,k,n)].x + pt[id(j+1,k,n)].x)/2.0;
        pt[id(j,k,n)].y =
          (pt[id(j-1,k,n)].y + pt[id(j+1,k,n)].y)/2.0;
      }
    }

/*
* Generate midside nodes in K direction
*/
    for (k = 1;k < n-1;k = k+node_inc)
    {
      for (j = 2;j < n-k;j = j+node_inc)
      {
        pt[id(j,k,n)].x =
          (pt[id(j,k-1,n)].x + pt[id(j,k+1,n)].x)/2.0;
        pt[id(j,k,n)].y =
          (pt[id(j,k-1,n)].y + pt[id(j,k+1,n)].y)/2.0;
      }
    }

/*
* Generate midside nodes in the third direction
*/
    for (k = 1;k < n-1;k = k+node_inc)
    {
      for (j = 1;j < n-k;j = j+node_inc)
      {
        pt[id(j,k,n)].x =
          (pt[id(j-1,k+1,n)].x + pt[id(j+1,k-1,n)].x)/2.0;
        pt[id(j,k,n)].y =
          (pt[id(j-1,k+1,n)].y + pt[id(j+1,k-1,n)].y)/2.0;
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
