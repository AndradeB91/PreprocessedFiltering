/*
** ---------------------------------------------------------------
**
** msh_copysurf.c: Functions to copy topological surface mesh using
** source to target aproach
**
** ---------------------------------------------------------------
**
** Created: May-2007      Antonio Miranda
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
#include "msh_smoothing.h"
#include "amr3bind.h"  /* Rtree */
#include "mshsurfmap.h"
#include "msh_copysurf.h"
#include "mshsurf_geo.h"

#define MSHSURF_LAPLACIAN_FACTOR 1.0
#define vector(a,b,c) \
	(a)[0] = (b)[0] - (c)[0];	\
	(a)[1] = (b)[1] - (c)[1];	\
	(a)[2] = (b)[2] - (c)[2];
#define innerProduct(v,q) \
		((v)[0] * (q)[0] + \
		(v)[1] * (q)[1] + \
		(v)[2] * (q)[2])	
#define crossProduct(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

/*   static functions
 *******************************************************************/
static void MshSurfCpGetBoundNodes (int nb, double *bound, int *n_node, double **coords);
static void LaplacianOperator      (double *curr_pos, int nadj, 
                                    double *adj_coords, double *v_lap);
static void GetLaplacianVectors    (void *surf, double *p_pos, double *v_lap);
#if 0
static void MshSurfCpSmooth        (void *surf, double *coords);
#endif
static void MshSurfBuildElemRtree      (void *s_mesh);
static void MshSurfInsertElemIntoRtree (void *surf, int *elem);
static int  rayIntersectsTriangle      (double *, double *, double *, double *, 
                                        double *, double *); 
static int  MshSurfSnapPoint  (void *curr_surf, double pts[3], double normal[3], 
                               double box, double snap_pts[3]);

#if 0
static int    GetNumAdjNodeNode (void *surf, int i);
static void   GetAdjNodeNode    (void *surf, int id, int *adj_nodes);
#endif

static int    minDist2          (double *min, double *p0, double *p1);
static void   MshSurfLSPLane    (int n, double *pts, double mt[4][4], int *i_, int *j, int *k);
static void   MshSurfGetXYScale (int n, double *p0, double mt0[4][4], 
								 double *p1, double mt1[4][4], double *sx, double *sy); 



/*******************************************************************/
/*******************************************************************/
/*******   C O P Y    M E S H    ***********************************/
/*******************************************************************/
/*******************************************************************/

/*   static variables
 *******************************************************************/
void *source_mesh = NULL;
void *target_mesh = NULL;
double base_target[3];
double base_source[3];
static int  source_base  = -1;
static int  target_base  = -1;
static void *elm3d_tree = NULL;


/*   external variables
 *******************************************************************/
extern int MshSurf_smoothOnProjection;

/* First Function - Set the source mesh */
/*******************************************************************/
void MshSurfSetSourceMesh (
double  base_pts[3],  /* point on boundary where inits surf copy (in) */
void    *surfmesh,    /* surface mesh using topology lib or the 
                      inserting a mesh as bellow              (in) */
int     n_node,       /* # of pts in the mesh                    (in) */
double  *coords,      /* coordinate array of the mesh            (in) */
int     n_elem,       /* number of elements generated            (in) */
int     *Conn         /* elem.connectivity list of the mesh      (in) */
)
{
  if (surfmesh != NULL)
  {
    source_mesh = surfmesh;
    source_base = SurfTopReturnNodeId (source_mesh, base_pts[0], base_pts[1], base_pts[2]);
  } 
  else
  {
    source_mesh = SurfTopInsertMesh (n_node, coords, n_elem, Conn, NULL);
    source_base = SurfTopReturnNodeId (source_mesh, base_pts[0], base_pts[1], base_pts[2]);
  }
  base_source[0] = base_pts[0];
  base_source[1] = base_pts[1];
  base_source[2] = base_pts[2];
} 

/* Second Function - Set the target mesh */
/*******************************************************************/
void MshSurfSetTargetMesh (
double  base_pts[3],  /* point on boundary where inits surf copy (in) */
void    *surfmesh,    /* surface mesh using topology lib or the 
                      inserting a mesh as bellow              (in) */
int     n_node,       /* # of pts in the mesh                    (in) */
double  *coords,      /* coordinate array of the mesh            (in) */
int     n_elem,       /* number of elements generated            (in) */
int     *Conn         /* elem.connectivity list of the mesh      (in) */
)
{
  if (surfmesh != NULL)
  {
    target_mesh = surfmesh;
    if (base_pts != NULL)
    {
      target_base = SurfTopReturnNodeId (target_mesh, base_pts[0], base_pts[1], base_pts[2]);
      base_target[0] = base_pts[0];
      base_target[1] = base_pts[1];
      base_target[2] = base_pts[2];
    }
  }
  else
  {
    target_mesh = SurfTopInsertMesh (n_node, coords, n_elem, Conn, NULL);
    if (base_pts != NULL)
    {
      target_base = SurfTopReturnNodeId (target_mesh, base_pts[0], base_pts[1], base_pts[2]);
      base_target[0] = base_pts[0];
      base_target[1] = base_pts[1];
      base_target[2] = base_pts[2];
    }
  }
} 

/* Third Function - Get the copy mesh */
/*******************************************************************/
int MshSurfGetCopyMesh( 
int     n_bound_nodes, /* # of boundary nodes (optional or 0)  (in) */
double  *bound_nodes,  /* boundary nodes (optional or NULL)    (in) */
int     *n_node,       /* # of pts in the mesh                (out) */
double  **coords,      /* coordinate array of the mesh        (out) */
double  **normal,      /* normal of the coordinate            (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn         /* elem.connectivity list of the mesh  (out) */
)
{
  int   n0, n1, l_n_boud;
  int   index, i, npts, *conn;
  void  *new_mesh;
  double *l_bound;
	
	if (source_mesh == NULL || target_mesh == NULL)
		return 0;

  if (n_bound_nodes == 0)
  {
    int *idsloop;
    /* get number of loops */
    n1 = SurfTopGetNumLoops (target_mesh);
    if (n1 > 1)
      return 0;

    SurfTopGetLoop (target_mesh, 0, &l_n_boud, &idsloop);
    l_bound = (double *) calloc (l_n_boud*3, sizeof(double));
    for (i = 0; i < l_n_boud; ++i)
      SurfTopGetCoordNode (target_mesh, idsloop[i], &(bound_nodes[i*3]));
  }
  else
  {
    l_n_boud = n_bound_nodes;
    l_bound  = bound_nodes;
  }

  /* compare number of boundary edges */
  n0 = SurfTopNumBoundEdge (source_mesh);
  if (n0 != l_n_boud)
    return 0;

  /* copy element connectivity */
  *n_elem = SurfTopNumElems (source_mesh);
  *Conn = (int *) calloc ((*n_elem) * 5, sizeof (int)); /* considering all quad elements */
  for (i = 0, index = 0; i < *n_elem; i++)
  {
    npts = SurfTopGetElemNNodes (source_mesh, i);
    conn = SurfTopGetElemConn (source_mesh, i);
    (*Conn)[index] = npts;
    memcpy (&((*Conn)[index+1]), conn, npts * sizeof (int));
    index = index + npts + 1;
  }

  /* obtain the position of boundary nodes */
  MshSurfCpGetBoundNodes (l_n_boud, l_bound, n_node, coords);

#if 0
  /* Smooth the node to a better position */
  for (i = 0; i < 100; i++)
    MshSurfCpSmooth (source_mesh, *coords);
#endif


  /* build topological structure to new mesh */
  new_mesh = SurfTopInsertMesh (*n_node, *coords, *n_elem, *Conn, NULL);


  /* Smooth the node to a better position */

   /* writeNF (new_mesh, "c:\\temp\\antes_smooth.nf"); */ 

  if (MshSurf_smoothOnProjection == 1)
  {
    for (i = 0; i < 1000; i++)
    SmoothLaplacian (new_mesh);
  }

  /* writeNF (new_mesh, "c:\\temp\\depois_smooth.nf"); */ 

  /* get current nodes */
  for (i = 0; i < *n_node; ++i)
    SurfTopGetCoordNode (new_mesh, i, &((*coords)[i*3]));
  

  /* Copy the node position to target surface */
  MshSurfCpToTarget (new_mesh, target_mesh, *coords, 1);

  /* writeNF (new_mesh, "c:\\temp\\superficie_smooth.nf"); */ 

  /* release meshes */
  SurfTopRelease (source_mesh);
  SurfTopRelease (target_mesh);
  SurfTopRelease (new_mesh);
  source_mesh = target_mesh = NULL;

  if (n_bound_nodes == 0)
    free (l_bound);


	return 1;
}



/* MshSurfCpToTarget - move coords to target surface 
 *******************************************************************/
void MshSurfCpToTarget (void *s_mesh, void *t_mesh, double *coords, int n_smooth)
{
  int    n_node, i, j;
  double try_coords[3], normal[3], snap_pts[3], *lap_vec, box;

  if (s_mesh == NULL)
    return;

  printf ("Smooth on projection = %d\n", MshSurf_smoothOnProjection);

  n_node = SurfTopNumNodes (s_mesh);
  lap_vec = (double *) calloc (n_node*3, sizeof (double));

  /* fill a range tree with bound box of target elements */
  MshSurfBuildElemRtree (t_mesh);

  /* smooth and snap to surface */
  for (j = 0; j < n_smooth; ++j)
  {
    /* Regular Laplacian Operator */
    GetLaplacianVectors (s_mesh, coords, lap_vec);

    /* update nodes coords */
    for (i = 0; i < n_node; i++)
    {
      if (!SurfTopIsBoundaryNode (s_mesh, i))
      {
        /* size box to search elements */
        box = MSHSURF_LAPLACIAN_FACTOR * fabs(lap_vec[i*3+0] + 
          lap_vec[i*3+1] + lap_vec[i*3+2]) / 3.0;
        if (box == 0.0)
          box = 0.001;

        /* try node coordinate */
#if 1
        if (MshSurf_smoothOnProjection == 0)
        {
          try_coords[0] = coords[i*3+0];
          try_coords[1] = coords[i*3+1];
          try_coords[2] = coords[i*3+2];
        }
        else
        {
          try_coords[0] = coords[i*3+0] + MSHSURF_LAPLACIAN_FACTOR * lap_vec[i*3+0];
          try_coords[1] = coords[i*3+1] + MSHSURF_LAPLACIAN_FACTOR * lap_vec[i*3+1];
          try_coords[2] = coords[i*3+2] + MSHSURF_LAPLACIAN_FACTOR * lap_vec[i*3+2];
        }
#else
        SurfTopGetCoordNode (s_mesh, i, try_coords);
#endif
        /* normal of current node */
        SurfTopGetNormalNode (s_mesh, i, normal);

        /* get to new node snapping to target surface */
        if (MshSurfSnapPoint (t_mesh, try_coords, normal, box, snap_pts))
        {
          /* update node coordinate */
          coords[i*3+0] = snap_pts[0];
          coords[i*3+1] = snap_pts[1];
          coords[i*3+2] = snap_pts[2];
          SurfTopSetCoordNode  (s_mesh, i, snap_pts);
        }

      }

    }

    SurfTopUdateNormalNodes (s_mesh);
  }


  free (lap_vec);
}


/************************************************************************/
int MshSurfCreate3dBil 
(
 int     np,            /* # of points                          (in) */
 double  *bry,          /* coordinates of all points            (in) */
 int     elem_type,     /* element type:  T3 (3), T6 (6), Q4 (4) ou  Q8 (8) (in) */
 int     diagtype,      /* Triangle option (option for cell diagonal) (in)
                        = 1 --> diagonal oriented to right direction
                        = 2 --> diagonal oriented to left direction
                        = 3 --> union jack alternation
                        = 4 --> optimum diagonal (smallest of two possible) */
int     *n_node,       /* # of pts in the mesh                (out) */
double  **coords,      /* coordinate array of the mesh        (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn         /* elem.connectivity list of the mesh  (out) */
)
{
  void *bil_mesh;

  if (target_mesh == NULL)
    return 0;

  // try to create a 3d bilinear
  if (!MshSurfTryBilinear (bry, np, elem_type, diagtype, n_node, n_elem, 
                           coords, Conn))
  {
    *n_node = *n_elem = 0;
    return 0;
  }

  // create a topology surface
  bil_mesh = SurfTopInsertMesh (*n_node, *coords, *n_elem, *Conn, NULL);
  if (bil_mesh == NULL)
  {
    *n_node = *n_elem = 0;
    return 0;
  }

  /* set all normal as plane */
  {
    int i;
    double plan_eqn[4], normal[3];
    MshSurfLstSqrPlanEqn (np, bry, plan_eqn);
    normal[0] = plan_eqn[0];
    normal[1] = plan_eqn[1];
    normal[2] = plan_eqn[2];
    for (i = 0; i < *n_node; ++i)
      SurfTopSetNormalNode (bil_mesh, i, normal);
  }


  /*if (MshSurf_smoothOnProjection == 1)
  {
    int i;
    for (i = 0; i < 1000; i++)
      SmoothLaplacian (bil_mesh);
  }

  writeNF (bil_mesh, "c:\\temp\\mapmesh.nf");*/


  /* move points to target surface */
  /* Copy the node position to target surface */
  MshSurfCpToTarget (bil_mesh, target_mesh, *coords, 1);

  /* release meshes */
  SurfTopRelease (bil_mesh);
  SurfTopRelease (target_mesh);
  target_mesh = NULL;
  bil_mesh    = NULL;

  return 1;
}

/*******************************************************************/
/*   static functions
 *******************************************************************/

/* MshSurfCpGetBoundNodes
 *******************************************************************/
static void MshSurfCpGetBoundNodes (int nb, double *bound, int *n_node, double **coords)
{
  int    n0, *pts0, i, axis_0, axis_j, axis_k;
  int    pos0 = 0, pos1 = 0;
  double pt0[3], /*delta[3],*/ min0 = -1.0, min1 = -1.0;
  double mt0[4][4], mt1[4][4];
  double *ord_bound0, *ord_bound1;
  double scale_x = 1.0, scale_y = 1.0;

  /* copy nodes */
  *n_node = SurfTopNumNodes (source_mesh);
  *coords = (double *) calloc ((*n_node) * 3, sizeof (double));
  ord_bound0 = (double *) calloc (nb * 3, sizeof (double));
  ord_bound1 = (double *) calloc (nb * 3, sizeof (double));

  /* get boundary loops */
  SurfTopGetLoop (source_mesh, 0, &n0, &pts0);

  /* get position of base points */
  for (i = 0; i < n0; i++)
  {
    SurfTopGetCoordNode (source_mesh, pts0[i], pt0);

    /* find base point in source mesh */
    if (minDist2 (&min0, pt0, base_source))
      pos0 = i;
    /* find base point in target mesh */
    if (minDist2 (&min1, &(bound[i*3]), base_target))
      pos1 = i;
  }

  /* initially, copy boundary points */
  for (i = 0; i < n0; i++)
  {
    int m = (pos0+i)%n0;  // source 
    int n = (pos1+i)%n0;  // target
    (*coords)[pts0[m]*3+0] = bound[n*3+0];
    (*coords)[pts0[m]*3+1] = bound[n*3+1];
    (*coords)[pts0[m]*3+2] = bound[n*3+2];

    // order boundaries
    SurfTopGetCoordNode (source_mesh, pts0[m], &(ord_bound0[i*3]));
    ord_bound1[i*3+0] = bound[n*3+0];
    ord_bound1[i*3+1] = bound[n*3+1];
    ord_bound1[i*3+2] = bound[n*3+2];
  }

  /* create matrix */
  axis_j = axis_k = -1;
  /* get local plane and the position of points to create local axis */
  MshSurfLSPLane (nb, ord_bound0, mt0, &axis_0, &axis_j, &axis_k);
  /* get local plane with point obtained before */
  MshSurfLSPLane (nb, ord_bound1, mt1, &axis_0, &axis_j, &axis_k);

  /* get local scale x and y */
  MshSurfGetXYScale (nb, ord_bound0, mt0, ord_bound1, mt1, &scale_x, &scale_y); 


  /* internal nodes */
  for (i = 0; i < (*n_node); i++)
  {
    double pt_tmp[3], pt1[3];

    if (SurfTopIsBoundaryNode (source_mesh, i))
      continue;

    SurfTopGetCoordNode (source_mesh, i, pt0);
    MshSurfTransf_xyx    (pt0, mt0, pt_tmp);
    
    /* apply scale */
    pt_tmp[0] *= scale_x; 
    pt_tmp[1] *= scale_y; 
    
    MshSurfInvTransf_xyz (pt_tmp, mt1, pt1);

    (*coords)[i*3+0] = pt1[0];
    (*coords)[i*3+1] = pt1[1];
    (*coords)[i*3+2] = pt1[2];
  }

}


/* LaplacianOperator
 *******************************************************************/
static void LaplacianOperator (double *curr_pos, int nadj, 
                               double *adj_coords, double *v_lap)
{
	int i;
  double dx, dy, dz;

  dx = dy = dz = 0.0;
  for (i = 0; i < nadj; i++)
  {
    dx += (adj_coords[i*3+0] - curr_pos[0]);
    dy += (adj_coords[i*3+1] - curr_pos[1]);
    dz += (adj_coords[i*3+2] - curr_pos[2]);
  }

  v_lap[0] = dx / nadj;
  v_lap[1] = dy / nadj;
  v_lap[2] = dz / nadj;
}

/* GetLaplacianVectors 
 *******************************************************************/
static void GetLaplacianVectors (void *surf, double *p_pos, double *v_lap)
{
  int    n_node, nadj_node, *adj_nodes, i, k;
  double coords[3], *adj_coords, old_coords[3];
  int    MAX_ADJ_EDGES = 50;
  
  if (surf == NULL)
    return;
  
  n_node = SurfTopNumNodes (surf);

  adj_coords = (double *) calloc (MAX_ADJ_EDGES*3, sizeof (double));
  adj_nodes  = (int *) calloc (MAX_ADJ_EDGES, sizeof (int));

  for (i = 0; i < n_node; i++)
  {

    if (SurfTopIsBoundaryNode (surf, i))
    {
      v_lap[i*3+0] = v_lap[i*3+1] = v_lap[i*3+2] = 0.0;
      continue;
    }

    /* number of adj nodes */
    nadj_node = SurfTopNumAdjNodeNode (surf, i);
    /* nadj_node = GetNumAdjNodeNode (surf, i); */

    /* old coordinate */
    if (p_pos == NULL)
    {
      SurfTopGetCoordNode   (surf, i, coords);
      memcpy (old_coords, coords, 3 * sizeof (double));
    }
    else
    {
      memcpy (old_coords, &(p_pos[i*3]), 3 * sizeof (double));
    }

    /* get adj ids of current node */
    SurfTopAdjNodeNode (surf, i, adj_nodes);
    /* GetAdjNodeNode (surf, i, adj_nodes); */

    /* get adj coords of current node */
    for (k = 0; k < nadj_node; k++)
    {
      if (p_pos == NULL)
      {
        SurfTopGetCoordNode   (surf, adj_nodes[k], coords);
        memcpy (&(adj_coords[k*3]), coords, 3 * sizeof (double));
      }
      else
      {
        memcpy (&(adj_coords[k*3]), &(p_pos[adj_nodes[k]*3]), 3 * sizeof (double));
      }
    }

    /* smooth */
    LaplacianOperator (old_coords, nadj_node, adj_coords, &(v_lap[i*3]));
  }

  free (adj_coords);
  free (adj_nodes);
}


#if 0
/* MshSurfCpSmooth - smooth the nodes
 *******************************************************************/
static void MshSurfCpSmooth (void *surf, double *coords)
{
  int    n_node, i, k;
  double old_coords[3], *lap_vec;
  
  if (surf == NULL)
    return;

  n_node = SurfTopNumNodes (surf);
  lap_vec = (double *) calloc (n_node*3, sizeof (double));

  /* Regular Laplacian Operator */
  GetLaplacianVectors (surf, coords, lap_vec);

  /* update nodes coords */
  for (i = 0; i < n_node; i++)
  {
    if (!SurfTopIsBoundaryNode (surf, i))
    {
      /* compute new coordinate */
      SurfTopGetCoordNode (surf, i, old_coords);
      for (k = 0; k < 3; k++)
        coords[i*3+k] = coords[i*3+k] + MSHSURF_LAPLACIAN_FACTOR * lap_vec[i*3+k];
    }
  }

  free (lap_vec);
}
#endif



static int *id_elem;
/* MshSurfBuildElemRtree
 *******************************************************************/
static void MshSurfBuildElemRtree (void *s_mesh)
{
  int i, ne;
  /* hold triangular element on rtree */
  elm3d_tree = RtreeCreate( );

  ne = SurfTopNumElems (s_mesh);
  id_elem = (int *) calloc (ne, sizeof (int));

  for (i = 0; i < ne; i++)
  {
    id_elem[i] = i;
    MshSurfInsertElemIntoRtree (s_mesh, &(id_elem[i]));
  }
}

/* MshSurfInsertElemIntoRtree
 *******************************************************************/
static void MshSurfInsertElemIntoRtree (void *surf, int *elem)
{
  int    *conn, nn, i;
  double pt[4], xmax, xmin, ymax, ymin, zmax, zmin;

  nn   = SurfTopGetElemNNodes (surf, *elem);
  conn = SurfTopGetElemConn  (surf, *elem);

  if (conn == NULL)
    return;

  if (nn >= 4)
    return;

  SurfTopGetCoordNode (surf, conn[0], pt);
  xmax = xmin = pt[0];
  ymax = ymin = pt[1];
  zmax = zmin = pt[2];

  for (i = 1; i < nn; i++)
  {
    SurfTopGetCoordNode (surf, conn[i], pt);

    if (pt[0] < xmin) xmin = pt[0];
    if (pt[1] < ymin) ymin = pt[1];
    if (pt[2] < zmin) zmin = pt[2];
    if (pt[0] > xmax) xmax = pt[0];
    if (pt[1] > ymax) ymax = pt[1];
    if (pt[2] > zmax) zmax = pt[2];

  }

  RtreeInsert (elm3d_tree, (void *)elem, xmin, xmax, ymin, ymax, zmin, zmax);
}



/* rayIntersectsTriangle
 *******************************************************************/
static int rayIntersectsTriangle 
(
  double *p, double *d, /* point and direction */
  double *v0, double *v1, double *v2, /* triangle */
  double *pts  /* return intersected point */
) 
{
	double e1[3], e2[3], h[3], s[3], q[3];
	double a, f, u, v, t;
	
	vector (e1, v1, v0);
	vector (e2, v2, v0);
	crossProduct (h, d, e2);
	a = innerProduct (e1, h);
	
	if (a > -0.00001 && a < 0.00001)
		return 0;
	
	f = 1.0 / a;
	vector (s,p,v0);
	u = f * (innerProduct (s, h));
	
	if (u < 0.0 || u > 1.0)
		return 0;
	
	crossProduct (q, s, e1);
	v = f * innerProduct (d, q);
	if (v < 0.0 || u + v > 1.0)
		return 0;
	// at this stage we can compute t to find out where 
	// the intersection point is on the line
	t = f * innerProduct (e2, q);
  pts[0] = p[0] + t * d[0];
  pts[1] = p[1] + t * d[1];
  pts[2] = p[2] + t * d[2];

  return 1;
}

/* MshSurfSnapPoint
 *******************************************************************/
static int MshSurfSnapPoint (void *curr_surf, double pts[3], double normal[3], 
                              double box, double snap_pts[3])
{
  int    k, i, *conn, n_nodes;
  double v1[3], v2[3], v3[3];

#if 0
  double box_fac = 1.0; 

  n = 0;
  do
  {
    /* find the nearest element */
    for (i = 0; i < 3; i++)
    {
      pmin[i] = pts[i] - box_fac * box * 0.5;
      pmax[i] = pts[i] + box_fac * box * 0.5;
    }

    RtreeInitSearchBox (elm3d_tree, pmin[0], pmax[0], pmin[1], pmax[1], pmin[2], pmax[2]);

    while ((cur_elem = (int *) RtreeSearchBox(elm3d_tree, &pmin[0], &pmax[0], &pmin[1], 
                                                &pmax[1], &pmin[2], &pmax[2])) != NULL)
    {
      n_nodes = SurfTopGetElemNNodes (curr_surf, *cur_elem);
      if (n_nodes < 3)
        continue;

      /* connectivity */
      conn = SurfTopGetElemConn (curr_surf, *cur_elem);
      SurfTopGetCoordNode (curr_surf, conn[0], v1);
      SurfTopGetCoordNode (curr_surf, conn[1], v2);
      SurfTopGetCoordNode (curr_surf, conn[2], v3);

      /* try to intersect with triangle */
      if (rayIntersectsTriangle (pts, normal, v1, v2, v3, snap_pts))
        return;

      /* in case of number of nodes in the element is four */
      if (n_nodes == 4)
        SurfTopGetCoordNode (curr_surf, conn[3], v1);
      if (rayIntersectsTriangle (pts, normal, v1, v2, v3, snap_pts))
        return;

    }

    box_fac *= 1.5; /* increase factor */

    ++n;
    if (n == 10)  /* try ten times */
      break;

  } while (1);

#else
  int ne = SurfTopNumElems (curr_surf);

  /* for all elements */
  for (i = 0; i < ne; ++i)
  {
    n_nodes = SurfTopGetElemNNodes (curr_surf, i);
    if (n_nodes < 3)
      continue;

    /* connectivity */

    conn = SurfTopGetElemConn (curr_surf, i);
    SurfTopGetCoordNode (curr_surf, conn[0], v1);
    SurfTopGetCoordNode (curr_surf, conn[1], v2);
    SurfTopGetCoordNode (curr_surf, conn[2], v3);

    /* try to intersect with triangle */
    if (rayIntersectsTriangle (pts, normal, v1, v2, v3, snap_pts))
      return 1;

    /* in case of number of nodes in the element is four */
    if (n_nodes == 4)
    {
      v2[0] = v3[0]; v2[1] = v3[1]; v2[2] = v3[2];
      SurfTopGetCoordNode (curr_surf, conn[3], v3);
    }
    if (rayIntersectsTriangle (pts, normal, v1, v2, v3, snap_pts))
      return 1;
  }


#endif

  return 0;

  /*  */
  for (k = 0; k < 3; k++)
     snap_pts[k] = pts[k];  
}

#if 0
/*******************************************************************/
static int GetNumAdjNodeNode( void *surf, int i )
{
#if 1
  return (SurfTopNumAdjNodeNode (surf, i));
#else
  int nnodes = 0;
  int j, ne, *elem;
  SurfTopAdjElemToNode (surf, i, &ne, &elem);

  for (j = 0; j < ne; ++j)
    nnodes += SurfTopGetElemNNodes (surf, elem[j]) - 2;

  free(elem);

  return nnodes*2;
#endif
}


/*******************************************************************/
static void GetAdjNodeNode( void *surf, int id, int *adj_nodes )
{
#if 1
  SurfTopAdjNodeNode (surf, id, adj_nodes);
#else
  int j, ne, *edges, ei, ej;

  SurfTopOppAdjEdgeToNode  (surf, id, &ne, &edges);
  for (j = 0; j < ne; ++j)
  {
    SurfTopAdjNodeToEdge (surf, edges[j], &ei, &ej);
    adj_nodes[j*2+0] = ei;
    adj_nodes[j*2+1] = ej;
  }
#endif
}
#endif


/*******************************************************************/
int minDist2( double *min, double *p0, double *p1 )
{
  double d[3] = {p0[0]-p1[0], p0[1]-p1[1], p0[2]-p1[2]};
  double dist2 = d[0]*d[0] + d[1]*d[1] + d[1]*d[1];
  if (*min == -1.0 || dist2 < *min)
  {
    *min = dist2;
    return 1;
  }
  return 0;
}

/*******************************************************************/
void MshSurfLSPLane( int n, double *pts, double mat[4][4], int *i_, int *j, int *k )
{
  int    i;
  double ltol = 0.0001, dist;
  double  seg1[3], seg2[3], cross[3];
  double  plan_eqn[4];

  /* set initial values to matrix */
  for(i = 0; i < 4; i++ )
    mat[i][0] = mat[i][1] = mat[i][2] = mat[i][3] = 0.0;
  mat[0][0] = mat[1][1] = mat[2][2] = mat[3][3] = 1.0;


  /* Find first two segments which are non-colinear. Use these to form the plane eqn. */
  if (*j == -1 && *k == -1)
  {
    /* get border stones */
    int nbordes, *borders;
    MshSurfFindBorders   (n, pts, &nbordes, &borders );

    /* look for better position to create a plane */
    for(i = 0; i < nbordes - 2; i++ )
    {
      MshSurfDiffVector( &pts[borders[i+1]*3], &pts[borders[0]*3], seg1 ) ;
      MshSurfDiffVector( &pts[borders[i+2]*3], &pts[borders[0]*3], seg2 ) ;
      MshSurfCrossProd( seg1, seg2, cross ) ;
      if (!(fabs(cross[0]) < ltol && fabs(cross[1]) < ltol && fabs(cross[2]) < ltol)) 
      {
        *i_ = borders[0];
        *j = borders[i+1];
        *k = borders[i+2];
        break;
      }
    }
  }

  /* compute z axis */
  MshSurfGenPlanEqn( &pts[(*i_)*3], &pts[(*j)*3], &pts[(*k)*3], plan_eqn);
  mat[2][0] = plan_eqn[0];
  mat[2][1] = plan_eqn[1];
  mat[2][2] = plan_eqn[2];

  /* set x axis  -- first segment */
  MshSurfDiffVector( &pts[(*j)*3], &pts[(*i_)*3], seg1 ) ;
  dist = MshSurfCompr (seg1);
  mat[0][0] = seg1[0] / dist;
  mat[0][1] = seg1[1] / dist;
  mat[0][2] = seg1[2] / dist;

  /* set y axis  -- first segment */
  MshSurfCrossProd (mat[0], mat[2], mat[1]);

  /* translate point */
  mat[3][0] = pts[(*i_)*3+0];
  mat[3][1] = pts[(*i_)*3+1];
  mat[3][2] = pts[(*i_)*3+2];
}

/*******************************************************************/
static void MshSurfGetXYScale (int n, double *p0, double mt0[4][4], 
                               double *p1, double mt1[4][4], double *sx, double *sy)
{
  double min0_x, min0_y, max0_x, max0_y;
  double min1_x, min1_y, max1_x, max1_y;
  double pt_tmp[3];
  int    i;

  /* init max and min values */
  MshSurfTransf_xyx (p0, mt0, pt_tmp);
  min0_x = max0_x = pt_tmp[0];
  min0_y = max0_y = pt_tmp[1];
  MshSurfTransf_xyx (p1, mt1, pt_tmp);
  min1_x = max1_x = pt_tmp[0];
  min1_y = max1_y = pt_tmp[1];

  /* get max and min */
  for (i = 1; i < n; ++i)
  {
    MshSurfTransf_xyx (&(p0[i*3]), mt0, pt_tmp);
    if (pt_tmp[0] < min0_x)
      min0_x = pt_tmp[0];
    if (pt_tmp[1] < min0_y)
      min0_y = pt_tmp[1];
    if (pt_tmp[0] > max0_x)
      max0_x = pt_tmp[0];
    if (pt_tmp[1] > max0_y)
      max0_y = pt_tmp[1];
    
    MshSurfTransf_xyx (&(p1[i*3]), mt1, pt_tmp);
    if (pt_tmp[0] < min1_x)
      min1_x = pt_tmp[0];
    if (pt_tmp[1] < min1_y)
      min1_y = pt_tmp[1];
    if (pt_tmp[0] > max1_x)
      max1_x = pt_tmp[0];
    if (pt_tmp[1] > max1_y)
      max1_y = pt_tmp[1];
  }

  *sx = (max1_x - min1_x) / (max0_x - min0_x);
  *sy = (max1_y - min1_y) / (max0_y - min0_y);

  printf ("Escalas = %lf, %lf\n", *sx, *sy);
}
