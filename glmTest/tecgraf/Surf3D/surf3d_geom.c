/*
** ----------------------------------------------------------------------
**
** surf3d_geom.c - Functions to geometric tests and transformations.  
**
** ----------------------------------------------------------------------
**
** Created:       21-Jan-2007    Antonio C.O. Miranda (surf3D)
**
** ----------------------------------------------------------------------
**
*/

#include <math.h>
#include "surf3d_geom.h"

#define TRUE                            1
#define FALSE                           0

/* --------------------------------------------------------------- */
double Surf3DVectLen2 (double *u)
{
  return (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
}

/* --------------------------------------------------------------- */
double Surf3DVectLen (double *u)
{
  return (sqrt (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]));
}

/* --------------------------------------------------------------- */
   /* Orientation:
 
   |v       
   |       / u
   |      /
   |     /
   |    /
   |   /
   |  /
   | /
   |/
   *
 */
void Surf3DCrossprod3d (double *u, double *v, double *w)
{
  w[0] = u[1]*v[2]-u[2]*v[1];
  w[1] = u[2]*v[0]-u[0]*v[2];
  w[2] = u[0]*v[1]-u[1]*v[0];
}

/* --------------------------------------------------------------- */
void Surf3DCrossprodNorm (double *u, double *v, double *w)
{
  double tmp[3], length;
  tmp[0] = u[1]*v[2]-u[2]*v[1];
  tmp[1] = u[2]*v[0]-u[0]*v[2];
  tmp[2] = u[0]*v[1]-u[1]*v[0];
  length = Surf3DVectLen (tmp);
  if (length > 0.0)
  {
    w[0] = tmp[0] / length;
    w[1] = tmp[1] / length;
    w[2] = tmp[2] / length;
  }
  else
  {
    w[0] = tmp[0]; w[1] = tmp[1]; w[2] = tmp[2];
  }
}


/* --------------------------------------------------------------- */
void Surf3DGetMatrixTransfX (double *xv, double matrix[3][3])
{
  double tmp_axis[3];
  int    i;

  /* try to set a local axis */
  for (i = 0; i < 3; i++)
  {
    double try_axis[3] = {0.0, 0.0, 0.0};
    try_axis[i] = 1.0;
    Surf3DCrossprod3d (xv, try_axis, tmp_axis);
    if (Surf3DVectLen2 (tmp_axis) > 0.0)
      break;
  }

  /* X axis */
  matrix[0][0] = xv[0];
  matrix[0][1] = xv[1];
  matrix[0][2] = xv[2];
  /* Y, tmp_axis = z */ 
  Surf3DCrossprodNorm (tmp_axis, matrix[0], matrix[1]);
  /* y */
  Surf3DCrossprodNorm (matrix[0], matrix[1], matrix[2]);
}

/* --------------------------------------------------------------- */
void Surf3DGetMatrixTransfZ (double *n, double matrix[3][3])
{
  double tmp_axis[3];
  int    i;

  /* try to set a local axis */
  for (i = 0; i < 3; i++)
  {
    double try_axis[3] = {0.0, 0.0, 0.0};
    try_axis[i] = 1.0;
    Surf3DCrossprod3d (n, try_axis, tmp_axis);
    if (Surf3DVectLen2 (tmp_axis) > 0.0)
      break;
  }

  /* Z axis */
  matrix[2][0] = n[0];
  matrix[2][1] = n[1];
  matrix[2][2] = n[2];
  /* X, tmp_axis = y */ 
  Surf3DCrossprodNorm (tmp_axis, matrix[2], matrix[0]);
  /* y */
  Surf3DCrossprodNorm (matrix[2], matrix[0], matrix[1]);
}

/* --------------------------------------------------------------- */
/* Surf3DTransf_xy - retorna os pontos x e y depois de uma transformacao
 * de translacao e rotacao 
 */
void Surf3DTransf_xy (double orig[3], double BaseMatrix[3][3], 
                      double p3d[3], double p2d[2])
{
  double dx = p3d[0] - orig[0];
  double dy = p3d[1] - orig[1];
  double dz = p3d[2] - orig[2];
  p2d[0] = BaseMatrix[0][0]*dx + BaseMatrix[0][1]*dy + BaseMatrix[0][2]*dz;
  p2d[1] = BaseMatrix[1][0]*dx + BaseMatrix[1][1]*dy + BaseMatrix[1][2]*dz;
} 

/*/ Surf3DTransf_XYZ
/////////////////////////////////////////////////////////////////////*/
void Surf3DTransf_XYZ (double orig[3], double BaseMatrix[3][3], 
                       double p2d[2],double p3d[3])
{
  double dx = BaseMatrix[0][0]*p2d[0] + BaseMatrix[1][0]*p2d[1];
  double dy = BaseMatrix[0][1]*p2d[0] + BaseMatrix[1][1]*p2d[1];
  double dz = BaseMatrix[0][2]*p2d[0] + BaseMatrix[1][2]*p2d[1];
  p3d[0] = orig[0] + dx;
  p3d[1] = orig[1] + dy;
  p3d[2] = orig[2] + dz;
}


/* --------------------------------------------------------------- */
/* Surf3DTransf_xy - retorna os pontos (x, y, z) locais depois de uma 
 * transformacao de translacao e rotacao 
 */
void Surf3DTransf_xyz (double orig[3], double BaseMatrix[3][3], 
                       double p3d[3], double p3d_l[3])
{
  double dx = p3d[0] - orig[0];
  double dy = p3d[1] - orig[1];
  double dz = p3d[2] - orig[2];
  p3d_l[0] = BaseMatrix[0][0]*dx + BaseMatrix[0][1]*dy + BaseMatrix[0][2]*dz;
  p3d_l[1] = BaseMatrix[1][0]*dx + BaseMatrix[1][1]*dy + BaseMatrix[1][2]*dz;
  p3d_l[2] = BaseMatrix[2][0]*dx + BaseMatrix[2][1]*dy + BaseMatrix[2][2]*dz;
}

/* --------------------------------------------------------------- */
/* Surf3DPlaneLineIntersec_xy - dado um plano BaseMatrix e orig, uma triangulo 
 * no espaco, tri[3][3], retorna os pontos (x, y) de intersecao locais. 
 */
int Surf3DPlaneLineIntersec_xy (double orig[3], double BaseMatrix[3][3],
                                double tri[3][3], int *n, double p2d[3][2])
{
  double p_xyz[3][3], z1_z0, t;
  int    i, j, k;

  /* obtem pontos no sistema local */
  for (i = 0; i < 3; i++)
    Surf3DTransf_xyz (orig, BaseMatrix, tri[i], p_xyz[i]);

  k = 0;
  for (i = 0; i < 3; i++)
  {
    j = (i+1)%3;
    /* avalia o valor da coodenada parametrica em z */
    z1_z0 = p_xyz[j][2] - p_xyz[i][2];

    /* z1_z0 == 0, significa linha toda no plano */
    if (z1_z0 == 0.0)
      continue;

    t = -1.0 * p_xyz[i][2] / z1_z0;

    /* intersecao fora da faixa da linha */
    if (t < -0.001 || t > 1.001)
      continue;

    /* pontos de intersecao */
    p2d[k][0] = p_xyz[i][0] + (p_xyz[j][0] - p_xyz[i][0]) * t;
    p2d[k][1] = p_xyz[i][1] + (p_xyz[j][1] - p_xyz[i][1]) * t;
    k++;

  }

  if (k < 2)
    return 0;

  *n = k;

  return 1;
}

/* --------------------------------------------------------------- */
double Surf3DVectArea (double *u, double *v)
{
  double w[3];
  Surf3DCrossprod3d (u, v, w);
  return ((w[0]*w[0] + w[1]*w[1] + w[2]*w[2]) * 0.5);
}

/* --------------------------------------------------------------- */
int Surf3DCrossEdge2d
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
  double PI1[2], PI2[2], PJ1[2], PJ2[2];
  double orig[3] = {0.0, 0.0, 0.0};

  for (i = 0; i < 3; i++)
    orig[i] = nodes[I1].coord[i];

  Surf3DTransf_xy (orig, BaseMatrix, nodes[I1].coord, PI1); 
  Surf3DTransf_xy (orig, BaseMatrix, nodes[I2].coord, PI2); 
  Surf3DTransf_xy (orig, BaseMatrix, nodes[J1].coord, PJ1); 
  Surf3DTransf_xy (orig, BaseMatrix, nodes[J2].coord, PJ2); 


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
  if ((PI1[0] > PJ1[0]) && (PI1[0] > PJ2[0]) && (PI2[0] > PJ1[0]) && (PI2[0] > PJ2[0]) )
  {
    cross = FALSE;
    return(cross);
  }

  if ((PI1[1] > PJ1[1]) && (PI1[1] > PJ2[1]) && (PI2[1] > PJ1[1]) && (PI2[1] > PJ2[1]) )
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


/* --------------------------------------------------------------- */
double Surf3DCrossprod2d (double I[2], double J[2], double K[2]) 
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

   kx = K[0] - I[0];
   ky = K[1] - I[1];
   jx = J[0] - I[0];
   jy = J[1] - I[1];
 
   return (kx*jy - ky*jx);
}


/* --------------------------------------------------------------- */
double Surf3DAreaEdge_Node (Surf3DBdryEdge *edge, int node_indx,
                                   Surf3DBdryNode *nodes)
{
  double u[3], v[3];
  int    i;

  for (i = 0; i < 3; i++)
  {
    u[i] = nodes[node_indx].coord[i]      - nodes[edge->verts[0]].coord[i];
    v[i] = nodes[edge->verts[1]].coord[i] - nodes[edge->verts[0]].coord[i];
  }

  return (Surf3DVectArea (u, v));
}

/* --------------------------------------------------------------- */
double Surf3DVectAngle (double *u, double *v)
{
  double angle = 0.0, lu, lv;

  lu = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
  lv = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  if (lu > 0.0 && lv > 0.0)
    angle = ((u[0]*v[0]) + (u[1]*v[1]) + (u[2]*v[2])) / (lu * lv);

  if (angle > 1.0) angle = 1.0;

  angle = acos (angle);

  return angle;
}

/* --------------------------------------------------------------- */
double Surf3DMinDist2Lines   (double L1i[3], double L1j[3], 
                              double L2i[3], double L2j[3])
{
  double n1[3], n2[3], n3[3], dist;
  int    i;

  for (i = 0; i < 3; i++)
  {
    n1[i] = L1j[i] - L1i[i];
    n2[i] = L2j[i] - L2i[i];
  }
  Surf3DCrossprodNorm (n1, n2, n3);
  dist = (L2i[0]-L1i[0])*n3[0] + (L2i[1]-L1i[1])*n3[1] + (L2i[2]-L1i[2])*n3[2];

  return dist;
}


#include <stdio.h>

/* --------------------------------------------------------------- */
int Surf3DPointInsideElem (double pts[3][3], double pts_idx[3], double base[3][3])
{
  int    i, count = 0;
  double orig[3], P[3][2], Pj[2];
  double *Pi, *Pk;
  
  /* origin */
  for (i = 0; i < 3; i++)
    orig[i] = pts[0][i];

  /* transform element vert to 2d */
  for (i = 0; i < 3; i++)
    Surf3DTransf_xy (orig, base, pts[i], P[i]);

  /* tranform node_indx to 2d */
  Surf3DTransf_xy (orig, base, pts_idx, Pj);
/*
  fprintf (stdout, "Surf3DPointInsideElem\n");
  fprintf (stdout, "%lf %lf\n", P[0][0], P[0][1]);
  fprintf (stdout, "%lf %lf\n", P[1][0], P[1][1]);
  fprintf (stdout, "%lf %lf\n", P[2][0], P[2][1]);
  fprintf (stdout, "%lf %lf\n", Pj[0], Pj[1]);
*/
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
    if (Surf3DCrossprod2d (Pi, Pj, Pk) > -0.0001)
      count++;
  }
  
  if (count == 3) return 1;

  return 0;
}
