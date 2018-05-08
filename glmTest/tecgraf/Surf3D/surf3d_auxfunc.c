/*
**  surf3d_auxfunc.c - This file contains auxiliar functions (AuxFuncGetElemSize,
    AuxFuncIdealPoints, AuxFuncSnapPoint) used in advance front algorithm.
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

//#include <GL/gl.h>
//#include <GL/glu.h>
#include "surf3d_auxfunc.h"
#include "surf3d_octree.h"
#include "surf3d_geom.h"
#include "surf3d.h"
#include "topsurfbin.h"
#include "amr3bind.h"

#define  MSHSURF_CURV  1

#define  ANG_MAX       18   // degree
#define  MAX_ANG_TEST  80   // degree

/*
** ---------------------------------------------------------------
** Local definitions and variables:
*/

static void *CurrOctree = NULL;
static void *curr_surf  = NULL;
static void *elm3d_tree = NULL;
static int  *id_elem    = NULL;
static void *nodes_tree = NULL;
static int  *id_node    = NULL;

static double msh_max_size = 0;

/*
** ---------------------------------------------------------------
** Private functions:
*/

static void InsertElemIntoRtree  (void *surf, int *elem);
static void InsertNodeIntoRtree  (void *surf, int *node, double size_box);
static void GetMatrixTransf (double normal[3], double dir[3], double matrix[3][3]);
static void GetBoundBox     (double pts[3], double normal[3], double direction[3],
                             double length, double pmin[3], double pmax[3]);
static int  FilterElement   (double orig[3], double base[3][3], double cur_len,
                             int *conn, double idealpts[3]);
static int AuxFuncGetPointInsideTriang (double pts[3], double matrix[3][3],
                                        int cur_elem, double snap_pts[3]);

static void MathInver3x3 (double a[3][3], double b[3][3], double *det);
#if 0
static int  ComputeEigenvalues (int n, double *pts, double eigenval[3]);
static void AuxFuncRefCurvedSurfaces (void **tree, void *parent, int leaf,
                                      double min[3], double max[3]);
#endif
static void RefineHighCurvatureSurf ( );


/*
** ---------------------------------------------------------------
** Public functions:
*/


/* ====================== AuxFuncSetCurrentMesh ========================= */
int  AuxFuncSetCurrentMesh (void *surf, double max_size, int np, double *pts,
                            int ne, int *edges)
{
  double max[3], min[3], delta[3], mid[3], len, box_size;
  int    i, j, nn, id0, id1;
  double bound_max[3], bound_min[3];

  msh_max_size = max_size;

  /* checks */
  if (surf == NULL)
    return 0;

  curr_surf = surf;

  /* set boundbox of the surface */
  SurfTopBoundBox (surf, min, max);

  /* get boundary box */
  bound_max[0] = bound_min[0] = pts[0];
  bound_max[1] = bound_min[1] = pts[1];
  bound_max[2] = bound_min[2] = pts[2];
  for (i = 0; i < np; i++)
  {
    for (j = 0; j < 3; j++)
    {
      if (pts[i*3+j] > bound_max[j])
        bound_max[j] = pts[i*3+j];
      if (pts[i*3+j] < bound_min[j])
        bound_min[j] = pts[i*3+j];
    }
  }

  /* new octree */
  if (CurrOctree != NULL)
    Surf3DOctreeRelease (CurrOctree);
  CurrOctree = Surf3DOctreeInit (max[0], max[1], max[2], min[0], min[1], min[2]);
  if( CurrOctree == NULL )
    return 0;

  /* set boundary box */
  Surf3DOctreeBound (bound_max[0], bound_max[1], bound_max[2],
                     bound_min[0], bound_min[1], bound_min[2]);

  /* complete the octree based on edges sizes */
  for (i = 0; i < ne; i++)
  {
    id0 = edges[i*2+0];
    id1 = edges[i*2+1];

    for (j = 0; j < 3; j++)
    {
      delta[j] = pts[id0*3+j] - pts[id1*3+j];
      mid[j]   = (pts[id0*3+j] + pts[id1*3+j]) * 0.5;
    }
    len = sqrt (delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);

    /*printf ("len = %lf ", len);*/

    if (len > 0.0) /* get mid point and size from surface */
      Surf3DOctreeAddPointSize (CurrOctree, mid[0], mid[1], mid[2], len);
  }

  /* refine octree based on maximum size and one diference level */
  Surf3DOctreeEnd (CurrOctree, max_size, 1);

  /* hold nodes on rtree */
  nodes_tree = RtreeCreate( );
  nn = SurfTopNumNodes (surf);
  id_node = (int *) calloc (nn, sizeof (int));
  box_size = ((max[0]-min[0])+(max[1]-min[1])+(max[2]-min[2])) * 0.3 * 0.0001;
  for (i = 0; i < nn; i++)
  {
    id_node[i] = i;
    InsertNodeIntoRtree (surf, &(id_node[i]), box_size);
  }

  /* hold triangular element on rtree */
  elm3d_tree = RtreeCreate( );
  ne = SurfTopNumElems (surf);
  id_elem = (int *) calloc (ne, sizeof (int));
  for (i = 0; i < ne; i++)
  {
    id_elem[i] = i;
    InsertElemIntoRtree (surf, &(id_elem[i]));
  }

#if MSHSURF_CURV
  /* refine octree cell on high curvature */
  RefineHighCurvatureSurf ( );

  /* refine octree based on one difference level */
  Surf3DOctreeEnd (CurrOctree, 0.0, 1);

#endif

  return 1;
}

/* =========================== AuxFuncRelease =========================== */
int AuxFuncGetNormal (double *pts, double hint_size, double *normal)
{
  double pmin[3], pmax[3];
  int    i, *cur_elem, n;
  double triang_n[3];

	if (curr_surf == NULL)
		return 0;

  if (hint_size == 0)
  	hint_size = 0.000001;

  do
  {
    /* find the nearest element/face */
    for (i = 0; i < 3; i++)
    {
      pmin[i] = pts[i] - hint_size * 0.5;
      pmax[i] = pts[i] + hint_size * 0.5;
    }
    RtreeInitSearchBox (elm3d_tree, pmin[0], pmax[0], pmin[1], pmax[1], pmin[2], pmax[2]);

    n = 0;
    normal[0] = normal[1] = normal[2] = 0.0;
    while ((cur_elem = (int *) RtreeSearchBox(elm3d_tree, &pmin[0], &pmax[0], &pmin[1],
                                                &pmax[1], &pmin[2], &pmax[2])) != NULL)
    {
      SurfTopGetElemNorm  (curr_surf, *cur_elem, triang_n);
    	normal[0] += triang_n[0];
    	normal[1] += triang_n[1];
    	normal[2] += triang_n[2];
      n++;
    }

    if (n > 0)
    {
    	normal[0] /= n;
    	normal[1] /= n;
    	normal[2] /= n;
      return 1;
    }
    else
      hint_size *= 1.5; /* increase box */

  } while (1);

  return 0;
}


/* =========================== AuxFuncRelease =========================== */
void AuxFuncRelease (void)
{
  if (elm3d_tree != NULL)
    RtreeDestroy( elm3d_tree );
  elm3d_tree = NULL;

  if (nodes_tree != NULL)
    RtreeDestroy (nodes_tree);
  elm3d_tree = NULL;

  if (id_elem != NULL)
    free (id_elem);
  id_elem = NULL;

  /* octree */
  Surf3DOctreeRelease (CurrOctree);
  CurrOctree = NULL;

  if (id_node != NULL)
    free (id_node);
  id_node = NULL;
}


/* ======================= AuxFuncGetElemSize =========================== */
void AuxFuncGetElemSize (double x, double y, double z, double *size)
{
	*size = Surf3DOctreeSize (CurrOctree, x, y, z);
  if (*size > msh_max_size)
    *size = msh_max_size;
}


/* ======================= AuxFuncIdealPoints =========================== */
int  AuxFuncIdealPoints (double pts[3], double normal[3], double direction[3],
                         double length, double ipts[3])
{
  double pmin[3], pmax[3], base_matriz[3][3];
  double /*tol = 0.001,*/ dot, triang_n[3];
  int    cand, *cur_elem, *conn;

  /* obtem o sistema de coordenada local, sendo que x = direction
     z = normal do elemento */
  GetMatrixTransf (normal, direction, base_matriz);

  /* obtain a point to build a boundbox */
  GetBoundBox (pts, normal, direction, length, pmin, pmax);

  /* find the nearest element/face */
  RtreeInitSearchBox (elm3d_tree, pmin[0], pmax[0], pmin[1], pmax[1], pmin[2], pmax[2]);

  cand = 0;
  while ((cur_elem = (int *) RtreeSearchBox(elm3d_tree, &pmin[0], &pmax[0], &pmin[1],
                                               &pmax[1], &pmin[2], &pmax[2])) != NULL)
  {
    SurfTopGetElemNorm  (curr_surf, *cur_elem, triang_n);
    dot = triang_n[0] * normal[0] + triang_n[1] * normal[1] + triang_n[2] * normal[2];
    if (dot > 0.0)
    {
      conn = SurfTopGetElemConn  (curr_surf, *cur_elem);

      cand = FilterElement (pts, base_matriz, length, conn, ipts);
      if (cand)
         break;
    }
  }

  return cand;
}


/* ========================= AuxFuncSnapPoint =========================== */
void AuxFuncSnapPoint (double pts[3], double normal[3], double box,
                       double snap_pts[3])
{
  double pmin[3], pmax[3], box_fac = 1.0;
  int    k, i, *cur_elem, n, n_try = 0;
  double dot, triang_n[3];
  double matrix_base[3][3];

/*  fprintf (stdout, "Tentando ponto = %lf %lf %lf\n", pts[0], pts[1], pts[2]); */

  /* obtain transformation to 2D matrix */
  Surf3DGetMatrixTransfZ (normal, matrix_base);

  // try to find
  do
  {
    /* find the nearest element/face */
    for (i = 0; i < 3; i++)
    {
      pmin[i] = pts[i] - box_fac * box * 0.5;
      pmax[i] = pts[i] + box_fac * box * 0.5;
    }
    RtreeInitSearchBox (elm3d_tree, pmin[0], pmax[0], pmin[1], pmax[1], pmin[2], pmax[2]);

    n = 0;
    while ((cur_elem = (int *) RtreeSearchBox(elm3d_tree, &pmin[0], &pmax[0], &pmin[1],
                                                &pmax[1], &pmin[2], &pmax[2])) != NULL)
    {
      SurfTopGetElemNorm  (curr_surf, *cur_elem, triang_n);
      dot = triang_n[0] * normal[0] + triang_n[1] * normal[1] + triang_n[2] * normal[2];

      if (dot > 0.0)
      {
        if (AuxFuncGetPointInsideTriang (pts, matrix_base, *cur_elem, snap_pts))
        {
          /* fprintf (stdout, "Pegou ponto\n"); */
          return;
        }
        n++;
      }
    }

    box_fac *= 1.5; /* increase factor */
    n_try++;

  } while (n_try < 5);

/*  fprintf (stdout, "Nao pegou valor\n"); */
  for (k = 0; k < 3; k++)
     snap_pts[k] = pts[k];
}


/*===================  AuxFuncDrawOctree =======================*/
void AuxFuncDrawOctree (void)
{
  if (CurrOctree != NULL)
	  Surf3DOctreeDraw (CurrOctree);
}


/*
** ---------------------------------------------------------------
** Private functions:
*/


/*===================  InsertElemIntoRtree =======================*/
static void InsertElemIntoRtree (void *surf, int *elem)
{
  int    *conn, nn, i;
  double pt[3], xmax, xmin, ymax, ymin, zmax, zmin;

  nn   = SurfTopGetElemNNodes (surf, *elem);
  conn = SurfTopGetElemConn  (surf, *elem);

  if (conn == NULL)
    return;

  if (nn != 3)
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

/*=================== InsertNodeIntoRtree =======================*/
static void InsertNodeIntoRtree  (void *surf, int *node, double size_box)
{
  double pt[3], xmax, xmin, ymax, ymin, zmax, zmin;

  SurfTopGetCoordNode (surf, *node, pt);
  xmax = pt[0] + size_box * 0.5;
  xmin = pt[0] - size_box * 0.5;
  ymax = pt[1] + size_box * 0.5;
  ymin = pt[1] - size_box * 0.5;
  zmax = pt[2] + size_box * 0.5;
  zmin = pt[2] - size_box * 0.5;

  RtreeInsert (nodes_tree, (void *)node, xmin, xmax, ymin, ymax, zmin, zmax);
}


/*===================  GetMatrixTransf =======================*/
static void GetMatrixTransf (double normal[3], double dir[3], double matrix[3][3])
{
  double lendir;

  lendir = Surf3DVectLen (dir);
  matrix[0][0] = dir[0] / lendir;
  matrix[0][1] = dir[1] / lendir;
  matrix[0][2] = dir[2] / lendir;
  /* Z direction */
  Surf3DCrossprodNorm (matrix[0], normal, matrix[2]);
  /* Y direction */
  Surf3DCrossprodNorm (matrix[2], matrix[0], matrix[1]);
}

/*===================  GetBoundBox =======================*/
static void GetBoundBox (double pts[3], double normal[3], double direction[3],
                  double length, double pmin[3], double pmax[3])
{
  double corner[4][3];
  int    i, j;

  /* first */
  for (i = 0; i < 3; i++)
  {
    corner[0][i] = pts[i] + 0.5 * length * normal[i];
    corner[1][i] = pts[i] - 0.5 * length * normal[i];
    corner[2][i] = corner[0][i] + 1.5 * length * direction[i];
    corner[3][i] = corner[1][i] + 1.5 * length * direction[i];
    pmin[i] = pmax[i] = corner[0][i];
  }

  for (i = 1; i < 4; i++)
  {
    for (j = 0; j < 3; j++)
    {
      if (corner[i][j] < pmin[j])
        pmin[j] = corner[i][j];
      if (corner[i][j] > pmax[j])
        pmax[j] = corner[i][j];
    }
  }
}


/*===================  FilterElement  =======================*/
static int FilterElement (double orig[3], double base[3][3], double cur_len, int *conn,
                          double idealpts[3])
{
  double tmp_dist, dist[2], dx, dy,  dr2, D, delta, dmin, dmax, pts[3];
  double tri[3][3], p2d[3][2], p2dt[3][2], localp2d[2], tmp_p2d[2][2], dy_sign;
  int    i, n, flag, i_min, i_max;
  double v_unit[2], dt;

  /* pega as coordenadas dos pontos */
  for (i = 0; i < 3; i++)
  {
    SurfTopGetCoordNode (curr_surf, conn[i], pts);
    tri[i][0] = pts[0];
    tri[i][1] = pts[1];
    tri[i][2] = pts[2];
  }
  /* pega os dois pontos de intesecao, se existir */
  if (!Surf3DPlaneLineIntersec_xy (orig, base, tri, &n, p2dt))
    return 0;

  /* pega os pontos maximos e minimos para teste */
  flag  = 1;
  i_min = i_max = -1;
  for (i = 0; i < n; i++)
  {
    tmp_dist = p2dt[i][0];
    if (flag)
    {
      i_min = i_max = i;
      flag = 0;
    }
    if (tmp_dist < p2dt[i_min][0])
      i_min = i;

    if (tmp_dist > p2dt[i_max][0])
      i_max = i;

  }

  if (i_min == -1 || i_max == -1)
    return 0;

  /* atualiza pontos */
  p2d[0][0] = p2dt[i_min][0];
  p2d[0][1] = p2dt[i_min][1];
  p2d[1][0] = p2dt[i_max][0];
  p2d[1][1] = p2dt[i_max][1];

  /* unit vector from min to max */
  v_unit[0] = p2d[1][0] - p2d[0][0];
  v_unit[1] = p2d[1][1] - p2d[0][1];
  dt = sqrt (v_unit[0]*v_unit[0] + v_unit[1]*v_unit[1]);

  dist[0] = p2d[0][0]*p2d[0][0] + p2d[0][1]*p2d[0][1];
  dist[1] = p2d[1][0]*p2d[1][0] + p2d[1][1]*p2d[1][1];

  if (p2d[0][0] < 0.0 && p2d[1][0] < 0.0)
    return 0;


  /* situacao especial apenas um ponto toca */
  if (dt == 0.0)
  {
    if (sqrt(dist[0])*0.999 < cur_len && sqrt(dist[0])*1.001 > cur_len && p2d[0][0] > 0.0)
    {
      Surf3DTransf_XYZ (orig, base, p2d[0], idealpts);
      return 1;
    }
  }

  /* update unit vector */
  v_unit[0] /= dt;
  v_unit[1] /= dt;

  /* intersecao de um circulo com a reta */
  dx  = p2d[1][0] - p2d[0][0];
  dy  = p2d[1][1] - p2d[0][1];
  dr2 = dx*dx + dy*dy;
  D  = p2d[0][0]*p2d[1][1] - p2d[1][0]*p2d[0][1];
  delta = cur_len*cur_len * dr2 - D*D;

  if (p2d[0][0] < 0.0)
  {
    dmin = - sqrt (dist[0]);
    dmax =   sqrt (dist[1]);
  }
  else
  {
    dmin = sqrt (dist[0]);
    dmax = dmin + sqrt(dr2);
  }


  /* fora do intervalo */
  if (!(cur_len > dmin && cur_len < dmax*1.001) || delta < 0 || dr2 <= 0.0)
    return 0;


  //fprintf (stdout, "OK\n");

  if (dy < 0)
    dy_sign = -1.0;
  else
    dy_sign = +1.0;

  // compute both points of inteserction
  tmp_p2d[0][0] = ( 1.0 * (D*dy) + dy_sign*dx*sqrt(delta)) / dr2;
  tmp_p2d[0][1] = (-1.0 * (D*dx) + fabs(dy)*sqrt(delta)) / dr2;

  tmp_p2d[1][0] = ( 1.0 * (D*dy) - dy_sign*dx*sqrt(delta)) / dr2;
  tmp_p2d[1][1] = (-1.0 * (D*dx) - fabs(dy)*sqrt(delta)) / dr2;

  /* check the points on the line */
  for (i = 0; i < 2; i++)
  {
    double vec_r[2], dot;
    vec_r[0] = tmp_p2d[i][0] - p2d[0][0];
    vec_r[1] = tmp_p2d[i][1] - p2d[0][1];
    dot = v_unit[0]*vec_r[0] + v_unit[1]*vec_r[1];

    if (tmp_p2d[i][0] > 0.0 && (dot > -0.01*dt && dot < 1.01*dt))
  {
      localp2d[0] = tmp_p2d[i][0];
      localp2d[1] = tmp_p2d[i][1];
      break;
  }

  }

  Surf3DTransf_XYZ (orig, base, localp2d, idealpts);

  return 1;
}

#if 0
/*************************** Surf3DInisideElem ******************************/
static int Surf3DInsideElem (double P[3][2], double Pj[2])
{
  int    i, count = 0;
  double *Pi, *Pk;

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
#endif


/* ====================== AuxFuncGetPointInsideTriang ========================= */
static int AuxFuncGetPointInsideTriang (double pts[3], double matrix_base[3][3],
                                        int cur_elem, double snap_pts[3])
{
  int    *conn, k, j;
  double P[3][3], C[3][3], D[3][3], pts2d[2], pt[3], det, lambda[3];

  conn = SurfTopGetElemConn  (curr_surf, cur_elem);
  if (conn == NULL)
    return 0;

  /* get corner points of triangle */
  for (k = 0; k < 3; k++)
  {
    SurfTopGetCoordNode (curr_surf, conn[k], pt);
    for (j = 0; j < 3; j++)
      P[k][j] = pt[j];
/*    fprintf (stdout, " %lf %lf %lf\n", P[k][0], P[k][1], P[k][2]); */
  }

  if (!Surf3DPointInsideElem (P, pts, matrix_base))
  {
/*    fprintf (stdout, "Determinante igual a zero\n"); */
    return 0;
  }

  /* coordenadas baricentricas */
  /* primeiro ponto */
  C[0][0] = 0.0;
  C[1][0] = 0.0;
  C[2][0] = 1.0;
  Surf3DTransf_xy  (P[0], matrix_base, P[1], pts2d);
  C[0][1] = pts2d[0];
  C[1][1] = pts2d[1];
  C[2][1] = 1.0;
  Surf3DTransf_xy  (P[0], matrix_base, P[2], pts2d);
  C[0][2] = pts2d[0];
  C[1][2] = pts2d[1];
  C[2][2] = 1.0;

  Surf3DTransf_XYZ (P[0], matrix_base, pts2d, snap_pts);
  /* inverte a matrix para obter as coordenadas parametricas */
  MathInver3x3 (C, D, &det);
  if (det == 0.0)
  {
/*    fprintf (stdout, "Determinante igual a zero\n"); */
    return 0;
  }

  Surf3DTransf_xy  (P[0], matrix_base, pts, pts2d);

  lambda[0] = D[0][0]*pts2d[0] + D[0][1]*pts2d[1] + D[0][2];
  lambda[1] = D[1][0]*pts2d[0] + D[1][1]*pts2d[1] + D[1][2];
  lambda[2] = D[2][0]*pts2d[0] + D[2][1]*pts2d[1] + D[2][2];

  for (k = 0; k < 3; k++)
    snap_pts[k] = P[0][k]*lambda[0] + P[1][k]*lambda[1] + P[2][k]*lambda[2];

/*  fprintf (stdout, "Ponto gerado = %lf %lf %lf\n", snap_pts[0], snap_pts[1], snap_pts[2]); */
    return 1;

}


/* ============================= MathInver3x3 =========================== */
static void MathInver3x3( double a[3][3], double b[3][3], double *det )
{
  int i,j;

  /* Clear inverse matrix.
  */
   for (i = 0;i < 3; i++)
     for (j = 0;j < 3; j++) b[i][j] = 0.0;

  /* Compute matrix determinant.
  */
  (*det) = a[0][0] * ( a[1][1] * a[2][2] - a[1][2] * a[2][1] ) -
           a[0][1] * ( a[1][0] * a[2][2] - a[1][2] * a[2][0] ) +
           a[0][2] * ( a[1][0] * a[2][1] - a[1][1] * a[2][0] );

  if ( fabs(*det) != 0.0  )
  {
    b[0][0] =  ( a[1][1] * a[2][2] - a[1][2] * a[2][1] ) / (*det);
    b[0][1] = -( a[0][1] * a[2][2] - a[0][2] * a[2][1] ) / (*det);
    b[0][2] =  ( a[0][1] * a[1][2] - a[0][2] * a[1][1] ) / (*det);
    b[1][0] = -( a[1][0] * a[2][2] - a[1][2] * a[2][0] ) / (*det);
    b[1][1] =  ( a[0][0] * a[2][2] - a[0][2] * a[2][0] ) / (*det);
    b[1][2] = -( a[0][0] * a[1][2] - a[0][2] * a[1][0] ) / (*det);
    b[2][0] =  ( a[1][0] * a[2][1] - a[1][1] * a[2][0] ) / (*det);
    b[2][1] = -( a[0][0] * a[2][1] - a[0][1] * a[2][0] ) / (*det);
    b[2][2] =  ( a[0][0] * a[1][1] - a[0][1] * a[1][0] ) / (*det);
  }
  else
    return;

}  /* End of MathInver3x3  */


#if 0
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);
/* Computes all eigenvalues and eigenvectors of a real symmetric matrix
a[0..n-1][0..n-1]. On output, elements of a above the diagonal are destroyed.
d[0..n-1] returns the eigenvalues of a. v[0..n-1][0..n-1] is a matrix whose
columns contain, on output, the normalized eigenvectors of a. nrot returns the
number of Jacobi rotations that were required. */
static int Jacobi (double **a, int n, double d[], double **v, int *nrot)
{
	int    i, j, ip, iq;
	double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

	b = (double *) calloc (n, sizeof (double));
  z = (double *) calloc (n, sizeof (double));

	for (ip = 0; ip < n; ip++)
  {
		for (iq = 0; iq < n; iq++)
      v[ip][iq] = 0.0;
		v[ip][ip] = 1.0;
	}
	for (ip = 0; ip < n; ip++)
  {
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}
	nrot=0;
	for (i = 1; i <= 50; i++)
  {
		sm = 0.0;
		for (ip = 0; ip < n-1; ip++)
    {
			for (iq = ip+1; iq < n; iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0)
    {
      free (b);
      free (z);
			return 1;
    }
		if (i < 4)
			tresh = 0.2 * sm / (n*n);
		else
			tresh = 0.0;
		for (ip = 0; ip < n-1; ip++)
    {
			for (iq = ip+1; iq < n; iq++)
      {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip]) && (fabs(d[iq])+g) == fabs(d[iq]))
        {
					a[ip][iq]=0.0;
        }
				else if (fabs(a[ip][iq]) > tresh)
        {
					h = d[iq]-d[ip];
					if ((fabs(h)+g) == fabs(h))
          {
						t = (a[ip][iq]) / h;
          }
					else
          {
						theta=0.5*h/(a[ip][iq]);
						t = 1.0 / (fabs(theta) + sqrt(1.0+theta*theta));
						if (theta < 0.0)
              t = -t;
					}
					c = 1.0 / sqrt(1+t*t);
					s = t*c;
					tau = s/(1.0+c);
					h = t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq] = 0.0;
					for (j = 0; j < ip; j++)
          {
						ROTATE (a,j,ip,j,iq);
          }
					for (j = ip+1; j < iq; j++)
          {
						ROTATE (a,ip,j,j,iq);
          }
					for (j = iq+1; j < n; j++)
          {
						ROTATE (a,ip,j,iq,j);
          }
					for (j = 0; j < n; j++)
          {
						ROTATE (v,j,ip,j,iq);
          }
					++nrot;
				}
			}
		}
		for (ip = 0; ip < n; ip++)
    {
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}

	return 0;
}


/*===================  ComputeEigenvalues  =======================*/
static int ComputeEigenvalues (int n, double *pts, double eigenval[3])
{
  int    i, m, k, nrot;
  double pm[3] = {0.0, 0.0, 0.0};
  double **C = NULL, *eigval, **eigvec;

  if (C == NULL)
  {
    C = (double **) malloc (3 * sizeof (double *));
    eigvec = (double **) malloc (3 * sizeof (double *));
    eigval = (double *) malloc (3 * sizeof (double));

    for (i = 0; i < 3; i++)
    {
      C[i] = (double *) malloc (3 * sizeof (double));
      memset (C[i], 0, 3 * sizeof (double));
      eigvec[i] = (double *) malloc (3 * sizeof (double));
    }
  }

  /* average point */
  for (i = 0; i < n; i++)
  {
    for (m = 0; m < 3; m++)
      pm[m] += pts[i*3+m];

  }
  for (m = 0; m < 3; m++)
    pm[m] /= n;

  /* matrix C */
  for (k = 0; k < 3; k++)
  {
    for (m = k; m < 3; m++)
    {
      for (i = 0; i < n; i++)
        C[k][m] += ((pts[i*3+k] - pm[k]) * (pts[i*3+m] - pm[m]));
    }
  }
  /* simetrical part */
  C[1][0] = C[0][1];
  C[2][0] = C[0][1];
  C[2][1] = C[1][2];

  /* calcula os autovalores */
  if (!Jacobi (C, 3, eigval, eigvec, &nrot))
    return 0;

  /* classifica e preenche a variavel de retorno */
  for (k = 0; k < 3; k++)
  {
    for (m = k; m < 3; m++)
    {
      if (eigval[m] < eigval[k])
      {
        double swap = eigval[m];
        eigval[m]   = eigval[k];
        eigval[k]   = swap;
      }
    }
  }

  for (k = 0; k < 3; k++)
    eigenval[k] = eigval[k];

  return 1;
}

#define MAX_TEST_NODES  500

/* ==================== AuxFuncRefCurvedSurfaces ========================= */
static void AuxFuncRefCurvedSurfaces (void **tree, void *parent, int leaf,
                                      double min[3], double max[3])
{
  int    *curr_node, n, i;
  double pmin[3], pmax[3], newmin[3], newmax[3], p[3];
  double sigma, eigenval[3];
  static  double *test_nodes = NULL;

  if (parent == NULL)
    return;

  if (!leaf)
    return;

  if (test_nodes == NULL)
    test_nodes = (double *) calloc (MAX_TEST_NODES*3, sizeof (double));

  for (i = 0; i < 3; i++)
  {
    newmin[i] = (min[i] + max[i]) * 0.5 - (max[i] - min[i]) * 2.0;
    newmax[i] = (min[i] + max[i]) * 0.5 + (max[i] - min[i]) * 2.0;
  }

  RtreeInitSearchBox (nodes_tree, newmin[0], newmax[0], newmin[1],
                                     newmax[1], newmin[2], newmax[2]);
  n = 0;
  while ((curr_node = (int *) RtreeSearchBox(nodes_tree, &pmin[0], &pmax[0], &pmin[1],
                                                &pmax[1], &pmin[2], &pmax[2])) != NULL)
  {
    /* get coordinate node */
    SurfTopGetCoordNode (curr_surf, *curr_node, p);
    if (p == NULL)
      continue;

    test_nodes[n*3+0] = p[0];
    test_nodes[n*3+1] = p[1];
    test_nodes[n*3+2] = p[2];

    n++;
    if (n == MAX_TEST_NODES)
      break;
  }

  if (n < 5)
    return;

  if (!ComputeEigenvalues (n, test_nodes, eigenval))
    return;

  if ((eigenval[0] + eigenval[1] + eigenval[2]) == 0.0)
    return;

  /* curvature metric */
  sigma = eigenval[0] / (eigenval[0] + eigenval[1] + eigenval[2]);

  if (sigma > 0.05)
  {
/*    int level = (int) (sigma / 0.05) + 1;
    Surf3DOctreeLevelRefine (tree, parent, level);
*/
/*    printf ("Octree sigma (%d) = %lf\n", n, sigma); */
  }

  return;
}
#endif

/* ====================== RefineHighCurvatureSurf ========================= */
static void RefineHighCurvatureSurf ( )
{
  int    i, k, ne, *curr_elem;
  double curr_normal[3], normal[3], min_cos, size, new_size, dot;
  double pmin[3], pmax[3], newmin[3], newmax[3], curr_center[3], center[3];
  double ang_max = ANG_MAX * 3.1416 / 180.0;
//   double COS_MIN = cos (ang_max);
  double COS_MIN_TEST = cos (MAX_ANG_TEST);
  int flag = 0;

  ne = SurfTopNumElems (curr_surf);

  /* compute for all elements */
  for (i = 0; i < ne; i++)
  {
    flag = 0;
    min_cos = 1.0;

    /* obtain the center of element */
    if (!SurfTopGetElemCenter  (curr_surf, i, curr_center))
      continue;

    /* get the current size in octree using the center of element */
    size = Surf3DOctreeSize (CurrOctree, curr_center[0], curr_center[1], curr_center[2]);

    /* compute a bound box */
    for (k = 0; k < 3; k++)
    {
      newmin[k] = curr_center[k] - size * 0.5;
      newmax[k] = curr_center[k] + size * 0.5;
    }

    /* find element into bound box */
    RtreeInitSearchBox (elm3d_tree, newmin[0], newmax[0], newmin[1],
                                       newmax[1], newmin[2], newmax[2]);

    /* current normal of element */
    SurfTopGetElemNorm  (curr_surf, i, curr_normal);

    /*  */
    while ((curr_elem = (int *) RtreeSearchBox(elm3d_tree, &pmin[0], &pmax[0], &pmin[1],
                                                  &pmax[1], &pmin[2], &pmax[2])) != NULL)
    {
      /* normal of elements */
      SurfTopGetElemNorm  (curr_surf, *curr_elem, normal);
      dot = curr_normal[0]*normal[0] + curr_normal[1]*normal[1] + curr_normal[2]*normal[2];

      if (dot <= COS_MIN_TEST)
        continue;

      if ((dot < min_cos) && (SurfTopGetElemCenter (curr_surf, *curr_elem, center)))
      {
        min_cos = dot;
        flag = 1;
      }
    }

    /* compute new size if min_cos < COS_30 */
    if (flag) // min_cos < COS_MIN)
    {
      double dt[3] = {center[0]-curr_center[0], center[1]-curr_center[1], center[2]-curr_center[2]};
      double d     = sqrt (dt[0]*dt[0] + dt[1]*dt[1] + dt[2]*dt[2]);
      double r, sin_;

      if (min_cos > 1.0)
        min_cos = 1.0;

      sin_ = sin (acos(min_cos)/2.0);

      if (sin_ <= 0.0)
        continue;

      r = d * 0.5 / sin (acos(min_cos)/2.0);

      new_size = r * ang_max;

      if (new_size <= 0.0)
        continue;

      Surf3DOctreeAddPointSize (CurrOctree, curr_center[0], curr_center[1], curr_center[2], new_size);
      Surf3DOctreeAddPointSize (CurrOctree, center[0], center[1], center[2], new_size);
    }

  }

  Surf3DOctreeEnd (CurrOctree, 0, 1);
}
