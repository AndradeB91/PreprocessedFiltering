/*
** ---------------------------------------------------------------
**
** msh_2dtosurf.c: Functions to create a 2d mesh a move to a 
**                 surface mesh using
**
** ---------------------------------------------------------------
**
** Created: July-2008      Antonio Miranda
**
** ---------------------------------------------------------------
**
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* #include <time.h> */
#include <string.h>

#include "mshsurf3d.h"
#include "topsurfbin.h"
#include "mshsurfmap.h"
#include "mshsurf_geo.h"
#include "msh_copysurf.h"

#include "msh2d.h"

/*   static functions
 *******************************************************************/
static void MshSurfGetMatrix     (int nNodes, int *idNodes, double mat[4][4]);
static void MshSurfGetMatrixPts  (int n, double *pts, double mat[4][4]);
static int  MshSurfGet2dBoundPts (int nNodes, int *idNodes, double mat[4][4], 
                                  double *pts2d);
static int  MshSurfGet2dPts      (int n, double *pts, double mat[4][4], 
                                  double *pts2d, double *h);
static int  MshSurfCreateMapping (int nPtsLoop, double *boundpts2d, int nb, int *bd, int invert, 
                                  int *n_node, double **pts2d, int *n_elem, int **Conn);

static void   MshSurf2dTo3d        (int n_node, double *p2d, int n_elem, int *conn, 
                                    int invert, double mat [4][4], double **p3d);
static void   MshSurfBordersStones (int nNodes, int *idNodes, int *nbordes, int **borders);
static void   MshSurfAlignWithBS   (int nNodes, int *idNodes, int nbordes, int *borders);

#if 0
static void   write_neutralfile    (char *fn, int n_node, double *Coords, int n_elem, int *conn);
#endif
static int    compmetr             ( const void *c1, const void *c2 );
static int    compInt              ( const void *c1, const void *c2 );

static int    MshSurfGetMappDirection4 (int *bd, int nPtsLoop, int *m, int *n);
static int    MshSurfGetMappDirection5 (int *bd, int nPtsLoop, int *m, int *n);
static int    MshSurfGetMappDirection6 (int *bd, int nPtsLoop, int *m, int *n);

static void   Msh2DInvertPoints        (int n_pts, double *boundpts2d);
static int    minDist2                 (double *min, double *p0, double *p1);

typedef struct
{
  int id;
  double value;
} metricBS;


/*******************************************************************/
/*******************************************************************/
/*******   2 D   M E S H    t o   a   S U R F A C E    *************/
/*******************************************************************/
/*******************************************************************/

/*   static variables
 *******************************************************************/
static void *support_mesh = NULL;

/* First Function - Set the support mesh */
/*******************************************************************/
void MshSurfSetSupport2d (
int     n_node,       /* # of pts in the mesh                    (in) */
double  *coords,      /* coordinate array of the mesh            (in) */
int     n_elem,       /* number of elements generated            (in) */
int     *Conn         /* elem.connectivity list of the mesh      (in) */
)
{
  support_mesh = SurfTopInsertMesh (n_node, coords, n_elem, Conn, NULL);
} 

/* Second Function - Create a new mesh using the current support mesh */
/*******************************************************************/
int MshSurfCreate2d (
int     alg_2d,        /* 2d algorithm                         (in)
                          0 - quadrilateral mapping
                          1 - boundary contraction
                          2 - quadrilateral with quadtree
                          3 - quadrilateral seam                    */
int     *n_node,       /* # of pts in the mesh                (out) */
double  **coords,      /* coordinate array of the mesh        (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn         /* elem.connectivity list of the mesh  (out) */
)
{
  int   nloops, nPtsLoop, *idPtsLoops;
  int   i, status, invert;
  double mat[4][4], *boundpts2d, *pts2d = NULL, base_pts[3], *normal;
  int *loop_segs, *edges, nbordes, *borders;
 
 if (support_mesh == NULL)
  return 0;

  if (alg_2d < 0 || alg_2d > 3)
    return 0;

  /* check number of loops */
  nloops = SurfTopGetNumLoops (support_mesh);
  if (nloops == 0)
    return 0;

  /* check for holes - the algorithm does not consider holes */
  if (nloops > 1) 
    return 0;

  /* get number of nodes on boundary */
  SurfTopGetLoop (support_mesh, 0, &nPtsLoop, &idPtsLoops);
  
  /* check for quadrilateral algorithms */
  if ((nPtsLoop%2) != 0  &&  /* even number boundary nodes */ 
       alg_2d != 1 /* boundary contraction*/ )
  {
    free(idPtsLoops);
    return 0;
  }

  /* get border stones */
  MshSurfBordersStones (nPtsLoop, idPtsLoops, &nbordes, &borders);

  /* Align idNodes with first border stone */
  MshSurfAlignWithBS (nPtsLoop, idPtsLoops, nbordes, borders);

  /* obtain matrix transformation to 2d */
  MshSurfGetMatrix   (nPtsLoop, idPtsLoops, mat);

  /* transform boundary pts to 2d                */
  /* obtain if the 2d boundary pts were inverted */ 
  boundpts2d = (double *) calloc(nPtsLoop*2, sizeof(double));
  invert = MshSurfGet2dBoundPts (nPtsLoop, idPtsLoops, mat, boundpts2d);

  if (invert == -1) /* zero area - impossible to create mesh */
  {
    free(boundpts2d);
    return 0;
  }

  /* create 2d mesh */
  status = 0;
  switch(alg_2d)
  {
  case 0:
    status = MshSurfCreateMapping (nPtsLoop, boundpts2d, nbordes, borders, invert, 
                                   n_node, &pts2d, n_elem, Conn);
   break;
  case 1:
    loop_segs = (int *) calloc(1, sizeof(int));
    loop_segs[0] = nPtsLoop;
    /* call mesh generation algorithm */
    status = Msh2DShape (1, loop_segs, boundpts2d, 3, 
                         n_node, &pts2d, n_elem, Conn);
    free(loop_segs);
   break;
  case 2:
    loop_segs = (int *) calloc(1, sizeof(int));
    loop_segs[0] = nPtsLoop;
    /* call mesh generation algorithm */
    status = Msh2DQuadBound (1, loop_segs, boundpts2d, 4, 1, 
                             n_node, &pts2d, n_elem, Conn);
    free(loop_segs);
      break;
  case 3:
    /* get edges */
    edges = (int *) calloc(nPtsLoop*2, sizeof(int));
    for (i = 0; i < nPtsLoop; ++i)
    {
      edges[2*i+0] = i;
      if (i == nPtsLoop-1)
        edges[2*i+1] = 0;
      else
       edges[2*i+1] = i+1;
    }
    /* call mesh generation algorithm */
    status = Msh2DQuadSeamEdge (nPtsLoop, boundpts2d, nPtsLoop, edges, 
                                n_node, &pts2d, n_elem, Conn);
    free(edges);
      break;
  default:
      break;
  }

#if 0
  if (status)
    write_neutralfile ("C:\\temp\\surf2d.nf", *n_node, pts2d, *n_elem, *Conn);
#endif

  /* check mesh generation */
  if (!status)
  {
    free(idPtsLoops);
    free(boundpts2d);
    free(borders);
    return 0;
  }
  
  /* transform to 3d */
  MshSurf2dTo3d (*n_node, pts2d, *n_elem, *Conn, invert, mat, coords);

  /* move mesh to surface */
  /* set target mesh */
  SurfTopGetCoordNode (support_mesh, idPtsLoops[0], base_pts);
  MshSurfSetTargetMesh (base_pts, support_mesh, 0, NULL, 0, NULL);
  /* set source mesh */
  MshSurfCoord_xyz(boundpts2d[0], boundpts2d[1], mat, base_pts);
  MshSurfSetSourceMesh (base_pts, NULL, *n_node, *coords, *n_elem, *Conn);

  /* Get the copy mesh */
  status = MshSurfGetCopyMesh (0, NULL, n_node, coords, &normal, n_elem, Conn);

  /* release mesh */
//   if (support_mesh != NULL)
//     SurfTopRelease (support_mesh);
  support_mesh = NULL;
  free (pts2d);
  free(idPtsLoops);
  free(boundpts2d);
  free(borders);

 return 1;
}

/* Second Function - Create a new mesh using the current support mesh */
/*                   This function use boundary from edges            */
int MshSurfCreate2dEdge (
int     alg_2d,        /* 2d algorithm                         (in)
                          0 - boundary contraction
                          1 - quadrilateral seam                    */
int     n_pts,         /* # of points                          (in) */
double  *bdry_pts,     /* coordinates of all points            (in) */
int     bound_edge,    /* # of boundary free edges             (in) */
int     *edges,        /* edge vector (i0,j0; i1, j1; ...)     (in) */
int     *n_node,       /* # of pts in the mesh                (out) */
double  **coords,      /* coordinate array of the mesh        (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn         /* elem.connectivity list of the mesh  (out) */
)
{
  int    status, i, id_i, id_j, invert = 0;
  double mat[4][4], *boundpts2d, *pts2d = NULL, *h, area = 0.0;
  double *reoder_pts;
  void   *mesh_2d;
  double _param = 1.0;

  if (support_mesh == NULL)
    return 0;

  if (alg_2d < 0 || alg_2d > 2)
    return 0;

  /* obtain matrix transformation to 2d */
  MshSurfLstSqrPlanMtx   (n_pts, bdry_pts, mat);

  /* reorder nodes to normalize */
  reoder_pts = (double *) calloc(n_pts*3, sizeof(double));
  for (i = 0; i < bound_edge; ++i)
  {
    reoder_pts[i*3+0] = bdry_pts[edges[i*2+0]*3+0];
    reoder_pts[i*3+1] = bdry_pts[edges[i*2+0]*3+1];
    reoder_pts[i*3+2] = bdry_pts[edges[i*2+0]*3+2];
    edges[i*2+0] = i;
    edges[i*2+1] = i+1;
  }
  edges[bound_edge*2-1] = edges[0];

  /* transform boundary pts to 2d                */
  /* obtain if the 2d boundary pts were inverted */ 
  boundpts2d = (double *) calloc(n_pts*2, sizeof(double));
  h          = (double *) calloc(n_pts, sizeof(double));
  MshSurfGet2dPts (n_pts, reoder_pts, mat, boundpts2d, h);

  // compute area
  area = MshSurfArea2DEdges (n_pts, boundpts2d, bound_edge, edges);

  // invert edge orientation if area > 0
  if (area > 0)
  {
    invert = 1;
    printf ("Inverte contorno\n");
    for (i = 0; i < bound_edge; ++i)
    {
      id_i  = edges[i*2+0];
      id_j  = edges[i*2+1];
      edges[i*2+0] = id_j;
      edges[i*2+1] = id_i;
    }
  }
  else if (area == 0)
    return 0;

  /* Msh2DSetRefFactor (1.0); */
  Msh2DEdgeParams (1, &_param); /* Generate intenal points */ 

  /* create 2d mesh */
  status = 0;
  switch(alg_2d)
  {
    case 0:
      /* call mesh generation algorithm */
      status = Msh2DEdge (n_pts, boundpts2d, bound_edge, 0, edges, 3, 
                          n_node, &pts2d, n_elem, Conn);
      break;
    case 1:


#if 0
      {
        FILE *ptr = fopen ("C:\\temp\\bound.m", "wt");
        fprintf (ptr, "hold on\n");
        for (i = 0; i < n_pts; ++i)
        {
          fprintf (ptr, "scatter(%lf,%lf)\n", boundpts2d[i*2+0], boundpts2d[i*2+1]);
          fprintf (ptr, "text(%lf,%lf,'%d')\n", boundpts2d[i*2+0], boundpts2d[i*2+1], i);
        }
        fclose (ptr);
      }
#endif
        /* call mesh generation algorithm */
        status = Msh2DQuadSeamEdge (n_pts, boundpts2d, bound_edge, edges, 
                                    n_node, &pts2d, n_elem, Conn);
      break;
    case 2:
      {
        // create boundary as loops
        //int n_loops = 1;
        int loops[1] = {n_pts};
        if (invert)
          Msh2DInvertPoints (n_pts, boundpts2d);
        
#if 0
        {
          FILE *ptr = fopen ("C:\\temp\\bound.dat", "wt");
          for (i = 0; i < n_pts; ++i)
            fprintf (ptr, "%lf %lf\n", boundpts2d[i*2+0], boundpts2d[i*2+1]);
          fclose (ptr);
        }
#endif

        status = Msh2DQuadBound (1, loops, boundpts2d, 4, 1, 
                                 n_node, &pts2d, n_elem, Conn);
      }
      break;
    default:
      break;
  }
  free (boundpts2d);

  if (!status || *n_node == 0 || *n_elem == 0)
  {
    return 0;
  }

#if 0
  if (status)
    write_neutralfile ("C:\\temp\\surf2d.nf", *n_node, pts2d, *n_elem, *Conn);
#endif


  /* transform to 3d */
  MshSurf2dTo3d (*n_node, pts2d, *n_elem, *Conn, invert, mat, coords);

  // release tmp points
  free (pts2d);

  // create a topology surface
  mesh_2d = SurfTopInsertMesh (*n_node, *coords, *n_elem, *Conn, NULL);

#if 0
  writeNF (mesh_2d, "C:\\temp\\surf3d.nf");
#endif

  if (mesh_2d == NULL)
  {
    *n_node = *n_elem = 0;
    return 0;
  }

  /* move boundary points to original boundary  */
  {
    int nloops = SurfTopGetNumLoops (mesh_2d);
    int npts_loop, *idpts_loop;
    int m_pos = 0;
    double min = -1.0;

    if (nloops != 1)
      return 0;
    
    SurfTopGetLoop (mesh_2d, 0, &npts_loop, &idpts_loop);
    if (npts_loop != n_pts)
    {
      free (reoder_pts);
      free (idpts_loop);
      return 0;
    }

    /* get first point position in loop */
    for (i = 0; i < npts_loop; ++i)
    {
      if (minDist2 (&min, &((*coords)[idpts_loop[i]*3]), reoder_pts))
        m_pos = i;
    }

    /* update boundary nodes */
    if (1) /* !invert) */
    {
      int k;
      for (i = 0; i < npts_loop; ++i)
      {
        k = (i+m_pos)%npts_loop;
        (*coords)[idpts_loop[k]*3+0] = reoder_pts[i*3+0];
        (*coords)[idpts_loop[k]*3+1] = reoder_pts[i*3+1];
        (*coords)[idpts_loop[k]*3+2] = reoder_pts[i*3+2];
      }
    }
    else
    {
      for (i = 0; i < npts_loop; ++i)
      {
        (*coords)[idpts_loop[i]*3+0] = reoder_pts[(npts_loop-i-1)*3+0];
        (*coords)[idpts_loop[i]*3+1] = reoder_pts[(npts_loop-i-1)*3+1];
        (*coords)[idpts_loop[i]*3+2] = reoder_pts[(npts_loop-i-1)*3+2];
      }
    }
    free (idpts_loop);
  }
  free (reoder_pts);

  // move domain points to target surface
  /* Copy the node position to target surface */
  MshSurfCpToTarget (mesh_2d, support_mesh, *coords, 1);

  /* release meshes */
  SurfTopRelease (mesh_2d);
  SurfTopRelease (support_mesh);

  return 1;
}





/*******************************************************************/
/*   static functions
 *******************************************************************/

/************************************************************************/
static void MshSurfGetMatrixPts  (int n, double *pts3d, double mat[4][4])
{
  int    i, j;
  double plan_eq[4], tmp_axis[3];

  /* get plane equation */
  MshSurfLstSqrPlanEqn (n, pts3d, plan_eq);
  free(pts3d);

  /* init transformation matrix */
  for (i = 0; i < 4; ++i)
  {
    for (j = 0; j < 4; ++j)
    {
      mat[i][j] = 0.0;
    }
  }
  mat[3][3] = 1;

  /* compute origin pts to translation operation */
  mat[3][0] = -1.0 * plan_eq[0] * plan_eq[3];
  mat[3][1] = -1.0 * plan_eq[1] * plan_eq[3];
  mat[3][2] = -1.0 * plan_eq[2] * plan_eq[3];

  /* try to set a local axis */
  for (i = 0; i < 3; i++)
  {
    double try_axis[3] = {0.0, 0.0, 0.0};
    try_axis[i] = 1.0;
    MshSurfCrossProd (plan_eq, try_axis, tmp_axis);
    if ( MshSurfSquareCompr (tmp_axis) > 0.0)
      break;
  }

  /* compute rotation matrix */
  /* Z */
  mat[2][0] = plan_eq[0];
  mat[2][1] = plan_eq[1];
  mat[2][2] = plan_eq[2];
  /* X, tmp_axis = y */
  MshSurfCrossProdNorm (tmp_axis, mat[2], mat[0]);
  /* y */
  MshSurfCrossProdNorm (mat[0], mat[2], mat[1]);
}

/************************************************************************/
void MshSurfGetMatrix( int nNodes, int *idNodes, double mat[4][4] )
{
  int    i;
  double *pts3d, coord[3];

  /* get pts */
  pts3d = (double *) calloc(nNodes*3, sizeof(double));
  for (i = 0; i < nNodes; ++i)
  {
    SurfTopGetCoordNode (support_mesh, idNodes[i], coord);
    pts3d[i*3+0] = coord[0];
    pts3d[i*3+1] = coord[1];
    pts3d[i*3+2] = coord[2];
  }

  MshSurfGetMatrixPts  (nNodes, pts3d, mat);
}

/************************************************************************/
static int MshSurfGet2dBoundPts( int nNodes, int *idNodes, double mat[4][4], double *pts2d )
{
  int i, j, invert = 0;
  double pts[3], area, h;

  for (i = 0; i < nNodes; ++i)
  {
    SurfTopGetCoordNode (support_mesh, idNodes[i], pts);
    MshSurfCoord_xy(pts[0], pts[1], pts[2], mat, &pts2d[i*2], &h);
  }

  /* compute area */
  area = MshSurfArea2D (nNodes, pts2d);
  if (area == 0.0)
    return -1;

  /* invert points orientation */
  if (area > 0.0)
  {
    double x, y;
    for (i = 1, j = nNodes-1; i < nNodes/2; ++i, --j)
    {
      x = pts2d[j*2+0];
      y = pts2d[j*2+1];
      pts2d[j*2+0] = pts2d[i*2+0];
      pts2d[j*2+1] = pts2d[i*2+1];
      pts2d[i*2+0] = x;
      pts2d[i*2+1] = y;
    }
    invert = 1;
  }

  return invert;
}

/************************************************************************/
static int MshSurfGet2dPts( int n, double *pts, double mat[4][4], 
                            double *pts2d, double *h )
{
  int i, invert = 0;

  for (i = 0; i < n; ++i)
    MshSurfCoord_xy(pts[i*3+0], pts[i*3+1], pts[i*3+2], mat, &pts2d[i*2], &h[i]);

#if 0
  double area = 0.0;

  /* compute area */
  area = MshSurfArea2D (n, pts2d);
  if (area == 0.0)
    return -1;

  /* invert points orientation */
  if (area > 0.0)
  {
    double x, y;
    for (i = 1, j = n-1; i < n/2; ++i, --j)
    {
      x = pts2d[j*2+0];
      y = pts2d[j*2+1];
      pts2d[j*2+0] = pts2d[i*2+0];
      pts2d[j*2+1] = pts2d[i*2+1];
      pts2d[i*2+0] = x;
      pts2d[i*2+1] = y;
    }
    invert = 1;
  }
#endif

  return invert;
}


/************************************************************************/
int MshSurfCreateMapping( int nPtsLoop, double *boundpts2d, int nb, int *bd, int invert,  
                         int *n_node, double **pts2d, int *n_elem, int **Conn )
{
  int status = 0, m, n;

  /* number of elements in each direction */
  switch(nb)
  {
  case 4:
    status = MshSurfGetMappDirection4 (bd, nPtsLoop, &m, &n);
   break;

  case 5:
    status = MshSurfGetMappDirection5 (bd, nPtsLoop, &m, &n);
   break;

  case 6:
    status = MshSurfGetMappDirection6 (bd, nPtsLoop, &m, &n);
      break;

  //case :
  //    break;

  default:
    return 0;
      break;
  }

  /* fail in finding number of elements in each direction */
  if (!status)
    return 0;
  
  /* invert orientation of boundary */
  if (invert)
  {
    int tmp = m;
    m = n;
    n = tmp;
  }

  /* mesh generation algorithm */
  status = Msh2DBilinear (boundpts2d, m, n, 4, 0, n_node, n_elem, pts2d, Conn);  

  return status;
}


#if 0
/* --------------------------------------------------------------- */
static void write_neutralfile( char *fn, int n_node, double *Coords, int n_elem, int *conn )
{

  int index =0, j, i;
  FILE *arquivo;
  int numq4 = 0, numt3 = 0, numq8 = 0, numt6 = 0;

  arquivo = fopen(fn,"w");

  fprintf(arquivo,"\n%%HEADER\nFile created by program 'Quebra2D' at - AMiranda\n");
  fprintf(arquivo,"\n%%HEADER.VERSION\n1.00\n");
  fprintf(arquivo,"\n%%HEADER.TITLE\n' untitled'\n");

  fprintf(arquivo,"%%NODE\n");
  fprintf(arquivo,"%d\n\n",n_node);
  fprintf(arquivo,"%%NODE.COORD\n");
  fprintf(arquivo,"%d\n\n",n_node);

  for(i=0; i<n_node; i++)
  {
    fprintf(arquivo,"%d    %f    %f    %f\n",i+1, Coords[i*2+0], Coords[i*2+1], 0.0);
  }

  fprintf(arquivo,"\n%%MATERIAL\n1\n");
  fprintf(arquivo,"\n%%MATERIAL.LABEL\n1\n1\t'mat1'\n");          
  fprintf(arquivo,"\n%%MATERIAL.ISOTROPIC\n1\n1\t10000\t0.2\n");  
  fprintf(arquivo,"\n%%THICKNESS\n1\n1\t1\n");                    
  fprintf(arquivo,"\n%%INTEGRATION.ORDER\n1\n1\t2\t2\t2\t2\t2\t2\n");


  index = 0;
  for(i=0; i<n_elem; i++)
  {
    if (conn[index] == 3)
      numt3++;
    if (conn[index] == 4)
      numq4++;
    if (conn[index] == 6)
      numt6++;
    if (conn[index] == 8)
      numq8++;

    index=index+conn[index]+1;
  }


  fprintf(arquivo,"%%ELEMENT\n");
  fprintf(arquivo,"%d\n\n",n_elem);

  index = 0;
  if (numt3 > 0)
  {
    fprintf(arquivo,"%%ELEMENT.T3\n");
    fprintf(arquivo,"%d\n", numt3);
    for(i=0; i<n_elem; i++)
    {
      if (conn[index] == 3)
      {
        fprintf(arquivo,"%d  1  1  1",i+1);
        for(j=0; j<conn[index]; j++)
          fprintf(arquivo,"  %d",conn[index+j+1]+1);
      }
      index=index+conn[index]+1;
      fprintf(arquivo,"\n");
    }
  }

  index = 0;
  if (numt6 > 0)
  {
    fprintf(arquivo,"%%ELEMENT.T6\n");
    fprintf(arquivo,"%d\n", numt6);
    for(i=0; i<n_elem; i++)
    {
      if (conn[index] == 6)
      {
        fprintf(arquivo,"%d  1  1  1",i+1);
        for(j=0; j<conn[index]; j++)
          fprintf(arquivo,"  %d",conn[index+j+1]+1);
      }
      index=index+conn[index]+1;
      fprintf(arquivo,"\n");
    }
  }

  index = 0;
  if (numq4 > 0)
  {
    fprintf(arquivo,"%%ELEMENT.Q4\n");
    fprintf(arquivo,"%d\n", numq4);
    for(i=0; i<n_elem; i++)
    {
      if (conn[index] == 4)
      {
        fprintf(arquivo,"%d  1  1  1",i+1);
        for(j=0; j<conn[index]; j++)
          fprintf(arquivo,"  %d",conn[index+j+1]+1);
      }
      index=index+conn[index]+1;
      fprintf(arquivo,"\n");
    }
  }

  index = 0;
  if (numq8 > 0)
  {
    fprintf(arquivo,"%%ELEMENT.Q8\n");
    fprintf(arquivo,"%d\n", numq8);
    for(i=0; i<n_elem; i++)
    {
      if (conn[index] == 8)
      {
        fprintf(arquivo,"%d  1  1  1",i+1);
        for(j=0; j<conn[index]; j++)
          fprintf(arquivo,"  %d",conn[index+j+1]+1);
      }
      index=index+conn[index]+1;
      fprintf(arquivo,"\n");
    }
  }

  fprintf(arquivo,"\n%%END\n");
  fclose(arquivo);
}
#endif


/* --------------------------------------------------------------- */
void MshSurf2dTo3d (int n_node, double *p2d, int n_elem, int *conn, 
                    int invert, double mat [4][4], double **p3d )
{
  int i, idtmp;
  double pt[3];

  /* move 2d mesh to 3d plane */
  *p3d = (double *) calloc(n_node*3, sizeof (double));
  for (i = 0; i < n_node; ++i)
  {
    MshSurfCoord_xyz(p2d[i*2+0], p2d[i*2+1], mat, pt);
    (*p3d)[i*3+0] = pt[0];
    (*p3d)[i*3+1] = pt[1];
    (*p3d)[i*3+2] = pt[2];
  }
  /* invert normal faces if it is necessary */
  if (invert)
  {
    int index = 0;
    for (i = 0; i < n_elem; ++i)
    {
      idtmp = conn[index+2];

      /* triangular element */
      if (conn[index] == 3)
      {
        conn[index+2] = conn[index+3];
        conn[index+3] = idtmp;
      }
      /* quadrilateral element */
      if (conn[index] == 4)
      {
        conn[index+2] = conn[index+4];
        conn[index+4] = idtmp;
      }

      index = index + conn[index] + 1;
    }
  }
}

/* --------------------------------------------------------------- */
static void MshSurfBordersStones( int nNodes, int *idNodes, 
                                  int *nbordes, int **borders )
{
  int i;
  double pprev[3], pcurr[3], pnext[3];

  /*fill metrics structure with current id node and angle*/
  metricBS *metrics = (metricBS *)calloc( nNodes, sizeof( metricBS ) );
  for( i = 0; i < nNodes; i++ )
  {
    int curr = idNodes[i];
    int next = idNodes[(i+1)%nNodes];
    int prev;

    if( i == 0 ) 
      prev = idNodes[nNodes-1];
    else         
      prev = idNodes[i-1];

    SurfTopGetCoordNode (support_mesh, prev, pprev);
    SurfTopGetCoordNode (support_mesh, curr, pcurr);
    SurfTopGetCoordNode (support_mesh, next, pnext);

    /* compute angle*/
    metrics[i].value = MshSurfAngleAboutPt( pcurr, pprev, pnext );
    metrics[i].id    = i;
    /*printf ("curr = %d , %f - %d-%d\n", curr, metrics[i].value, prev, next);*/
  }

  /* sort metrics*/
  qsort( metrics, nNodes, sizeof( metricBS ), compmetr );

  /* count number of border stones*/
  *nbordes = 0;
  for( i = 0; i < nNodes; ++i )
  {
    /* 2.356 = 3 / 4 * PI  = 135 degrees*/
    if( metrics[i].value < 2.356 )
      ++(*nbordes);
  }


  /* consider at least four border stones*/
  if (*nbordes < 4)
    *nbordes = 4;

  /* fill borders id vector*/
  *borders = (int *) calloc(*nbordes, sizeof(int));
  for( i = 0; i < *nbordes; ++i )
    (*borders)[i] = metrics[i].id;

  free( metrics );
}

/* --------------------------------------------------------------- */
static int compmetr( const void *c1, const void *c2 )
{
  metricBS *m1 = (metricBS *) c1;
  metricBS *m2 = (metricBS *) c2;
  if      (m1->value < m2->value) return -1;
  else if (m1->value > m2->value) return  1;
  else                            return  0;
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

    // correct id border stones
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

  // order border stones id sequentially to loop 
  qsort (borders, nbordes, sizeof( int ), compInt );

}

/* --------------------------------------------------------------- */
int MshSurfGetMappDirection5( int *bd, int nPtsLoop, int *m, int *n )
{
  int i, div[4];
  int b[5] = {bd[1]-bd[0], bd[2]-bd[1], bd[3]-bd[2], bd[4]-bd[3], nPtsLoop-bd[4]+bd[0]};

  for (i = 0; i < 4; ++i)
  {
    div[0] = b[i] + b[(i+1)%5];
    div[1] = b[(i+2)%5];
    div[2] = b[(i+3)%5];
    div[3] = b[(i+4)%5];
    if (div[0] == div[2] && div[1] == div[3]) // match
    {
      *m = div[i%2+0] + 1;
      *n = div[i%2+1] + 1;
      return 1;
    }
  }
  return 0;
}

/* --------------------------------------------------------------- */
int MshSurfGetMappDirection4( int *bd, int nPtsLoop, int *m, int *n )
{
  int b[4] = {bd[1]-bd[0], bd[2]-bd[1], bd[3]-bd[2], nPtsLoop-bd[3]+bd[0]};

  if (b[0] == b[2] && b[1] == b[3]) // match
  {
    *m = b[0] + 1;
    *n = b[1] + 1;
    return 1;
  }
  return 0;
}

/* --------------------------------------------------------------- */
int MshSurfGetMappDirection6( int *bd, int nPtsLoop, int *m, int *n )
{
  int i, div[4];
  int b[6] = {bd[1]-bd[0], bd[2]-bd[1], bd[3]-bd[2], bd[4]-bd[3], bd[5]-bd[4], nPtsLoop-bd[5]+bd[0]};

  /* first situation */
  for (i = 0; i < 4; ++i)
  {
    div[0] = b[i]       + b[(i+1)%6];
    div[1] = b[(i+2)%6] + b[(i+3)%6];
    div[2] = b[(i+4)%6];
    div[3] = b[(i+5)%6];
    if (div[0] == div[2] && div[1] == div[3]) // match
    {
      *m = div[i%2+0] + 1;
      *n = div[i%2+1] + 1;
      return 1;
    }
  }

  /* second situation */
  for (i = 0; i < 4; ++i)
  {
    div[0] = b[i]       + b[(i+1)%6];
    div[1] = b[(i+2)%6];
    div[2] = b[(i+3)%6] + b[(i+4)%6];
    div[3] = b[(i+5)%6];
    if (div[0] == div[2] && div[1] == div[3]) // match
    {
      *m = div[i%2+0] + 1;
      *n = div[i%2+1] + 1;
      return 1;
    }
  }

  return 0;
}

/*********************************************************************/
void Msh2DInvertPoints( int n, double *boundpts2d )
{
  int    i;
  double tmp;
  for (i = 0; i < n/2; ++i)
  {
    tmp = boundpts2d[i*2+0];
    boundpts2d[i*2+0] = boundpts2d[(n-1-i)*2+0];
    boundpts2d[(n-1-i)*2+0] = tmp;
    tmp = boundpts2d[i*2+1];
    boundpts2d[i*2+1] = boundpts2d[(n-1-i)*2+1];
    boundpts2d[(n-1-i)*2+1] = tmp;
  }
}

/*******************************************************************/
int minDist2( double *min, double *p0, double *p1 )
{
  double d[3] = {p0[0]-p1[0], p0[1]-p1[1], p0[2]-p1[2]};
  double dist2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
  if (*min == -1.0 || dist2 < *min)
  {
    *min = dist2;
    return 1;
  }
  return 0;
}
