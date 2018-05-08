/*
**  surf3d_auxfunc.c - This file contains auxiliar functions (AuxFuncGetElemSize,
    AuxFuncIdealPoints, AuxFuncSnapPoint) used in advance front algorithm.

	Author: Antonio Miranda - Nov/2007.

	Modified: Insert quadtree as background mesh.
	          By Antonio Miranda - March/2008.
*/


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

//#include <GL/gl.h>
//#include <GL/glu.h>
#include "surf3d_analt_auxfunc.h"
#include "msh_quadsurf.h"
/* #include "surf3d_octree.h" */
#include "surf3d.h"


#define  MSHSURF_CURV  0


/*
** ---------------------------------------------------------------
** Local definitions and variables:
*/

static double local_tol = 0.001;

/*
** ---------------------------------------------------------------
** Private functions:
*/

static void AuxGetGradVector (double u, double v, double *dev);


/*
** ---------------------------------------------------------------
** Public functions:
*/

static  int (*AuxFuncGetUVfromXYZ)   (double xyz[3], double tol, double uv[2]) = NULL;
static  int (*AuxFuncGetXYZfromUV)   (double uv[2], double xyz[3])  = NULL;
static  int (*AuxFuncGetGradVectors) (double uv[2], double Gu[3], double Gv[3]) = NULL;


/* AuxFuncSetCurrentSurface - define the current surface using feedback functions */
int  AuxFuncSetCurrentSurface
(
  int (*getUVfromXYZ)   (double xyz[3], double tol, double uv[2]),
  int (*getXYZfromUV)   (double uv[2], double xyz[3]),
  int (*getGradVectors) (double uv[2], double Gu[3], double Gv[3])
)
{
  if (getUVfromXYZ == NULL || getXYZfromUV == NULL || getGradVectors == NULL)
    return 0;

  AuxFuncGetUVfromXYZ   = getUVfromXYZ;
  AuxFuncGetXYZfromUV   = getXYZfromUV;
  AuxFuncGetGradVectors = getGradVectors;
  return 1;
}

/* ====================== AuxFuncSetCurrentMesh ========================= */
int  AuxFuncSetCurrentBoundary (double max_size, int np, double *pts,
                                int ne, int *edges)
{
  double delta[3], len, xyz[3], uv[2];
  int    i, j, id0, id1/*, status*/;

  double *par_pts; /* points in parametric space */

  par_pts = (double *) malloc (np * 2 * sizeof (double));

  /* find a snallest edge to compute a tolerance */
  for (i = 0; i < ne; i++)
  {
    id0 = edges[i*2+0];
    id1 = edges[i*2+1];
    for (j = 0; j < 3; j++)
    {
      delta[j] = pts[id0*3+j] - pts[id1*3+j];
    }
    len = sqrt (delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
    if ((local_tol == 0.0 || local_tol > len) && len > 0.0)
      local_tol = len;
  }
  local_tol /= 1000;

  /* get points in parametric space */
  for (i = 0; i < ne; i++)
  {
    xyz[0] = pts[i*3+0];
    xyz[1] = pts[i*3+1];
    xyz[2] = pts[i*3+2];
    AuxFuncGetUVfromXYZ (xyz, local_tol, uv);
    par_pts[i*2+0] = uv[0];
    par_pts[i*2+1] = uv[1];
  }

  /*status = */MshSurfGenQuadTree (np, ne, (double (*)[2]) par_pts, (int (*)[2]) edges, /* IN */
                               AuxGetGradVector); /* IN */

  // set maximum size of quadtree
  MshSurfSetMaxQuadSize (max_size);

  return 1;
}

/* =========================== AuxFuncRelease =========================== */
int AuxFuncGetNormal2 (double *pts, double hint_size, double *normal)
{
  int    n;
  double ptsv[3] = {pts[0], pts[1], pts[2]};
  double gu[3], gv[3];
  double uv[2], len;

  if (hint_size == 0)
  	hint_size = 0.000001;

  if (AuxFuncGetUVfromXYZ == NULL)
  	return 0;

  n = 0;
  do
  {
    /* find the uv pts */
    if (AuxFuncGetUVfromXYZ (ptsv, hint_size, uv))
    {
    	AuxFuncGetGradVectors (uv, gu, gv);
    	normal[0] = gv[1] * gu[2] - gu[1] * gv[2];
    	normal[1] = gv[2] * gu[0] - gu[2] * gv[0];
    	normal[2] = gv[0] * gu[1] - gu[0] * gv[1];
      len = sqrt (normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
      if (len > 0.0)
      {
        normal[0] /= len;
        normal[1] /= len;
        normal[2] /= len;
    	  return 1;
      }
    }
    else
      hint_size *= 2.0;

    ++n;

    if (n == 20)
    	return 0;

  } while (1);

  return 0;
}


/* =========================== AuxFuncRelease =========================== */
void AuxFuncRelease2 (void)
{
  MshSurfQuadFreeAll();
}


/* ======================= AuxFuncGetElemSize =========================== */
void AuxFuncGetElemSize2 (double x, double y, double z, double *size)
{
  double xyz[3] = {x, y, z};
  double uv[2], _new[2];
  int    level;

  AuxFuncGetUVfromXYZ (xyz, local_tol, uv);
	*size = MshSurfOptimalNodes (uv, _new, &level ) ;
}


/* ======================= AuxFuncIdealPoints =========================== */
int  AuxFuncIdealPoints2 (double pts[3], double normal[3], double direction[3],
                          double length, double ipts[3])
{
  double trypts[3], uv[2];
  int    n = 1;

  trypts[0] = pts[0] + length * direction[0];
  trypts[1] = pts[1] + length * direction[1];
  trypts[2] = pts[2] + length * direction[2];

  do
  {
    if (AuxFuncGetUVfromXYZ (trypts, length*0.1*n, uv))
      return AuxFuncGetXYZfromUV (uv, ipts);

    ++n;

    if (n == 10)
    	return 0;
  }
  while (1);

  return 0;
}


/* ========================= AuxFuncSnapPoint =========================== */
void AuxFuncSnapPoint2 (double pts[3], double normal[3], double box,
                        double snap_pts[3])
{
  double uv[2];
  int    n = 1;

  do
  {
    if (AuxFuncGetUVfromXYZ (pts, box*0.01*n, uv))
    {
      AuxFuncGetXYZfromUV (uv, snap_pts);
      return;
    }
    ++n;

    if (n == 10)
    	break;
  }
  while (1);

  snap_pts[0] = pts[0];
  snap_pts[1] = pts[1];
  snap_pts[2] = pts[2];
}



/*
** ---------------------------------------------------------------
** Private functions:
*/


/*****************************************************************/
static void AuxGetGradVector (double u, double v, double *dev)
    {
  double uv[2] = {u, v};
  double gu[3] = {0.0, 0.0, 0.0};
  double gv[3] = {0.0, 0.0, 0.0};
  AuxFuncGetGradVectors (uv, gu, gv);
  dev[0] = gu[0];
  dev[1] = gu[1];
  dev[2] = gu[2];
  dev[3] = gv[0];
  dev[4] = gv[1];
  dev[5] = gv[2];
}
