/*
** ----------------------------------------------------------------------
**
** mshsurf_geo.c - File for auxiliar 3D meshing routine. 
**
** ----------------------------------------------------------------------
**
** Created:      April-2010      Antonio C.O. Miranda
**
** ----------------------------------------------------------------------
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "mshsurf_geo.h"


typedef struct
{
  int id;
  double value;
} metricBs3D;

static int compmetr( const void *c1, const void *c2 );


#define pi 3.141592654


/*
//////////////////////////////////////////////////////////////////////////
*/

/************************************************************************/
int  MshSurfLstSqrPlanEqn(
                          int    num_pts, /* number of points */
                          double *pts,     /* list of points */
                          double plan_eqn[4]
)
{
  double a11,a12,a13,a22,a23,a33 ; /* A matrix                */
  double b11,b12,b13,b22,b23,b33 ; /* A inverse matrix        */
  double sum_x,sum_y,sum_z ;      /* summed values           */
  double div,norm ;             /* auxiliary coefficient.  */
  int  i ;                      /* counter                */
  double ltol = 0.0001;



  /* loop through the points and calculate the sums and products  */

  a11 = 0.0 ;     a12 = 0.0 ;     a13 = 0.0 ;
  a22 = 0.0 ;     a23 = 0.0 ;     a33 = 0.0 ;

  sum_x = 0.0 ;   sum_y = 0.0 ;   sum_z = 0.0 ;

  for ( i = 0; i < num_pts; i++ ) {
    a11 = a11 + pts[i*3+0] * pts[i*3+0] ;
    a12 = a12 + pts[i*3+0] * pts[i*3+1] ;
    a13 = a13 + pts[i*3+0] * pts[i*3+2] ;
    a22 = a22 + pts[i*3+1] * pts[i*3+1] ;
    a23 = a23 + pts[i*3+1] * pts[i*3+2] ;
    a33 = a33 + pts[i*3+2] * pts[i*3+2] ;

    sum_x = sum_x + pts[i*3+0] ;
    sum_y = sum_y + pts[i*3+1] ;
    sum_z = sum_z + pts[i*3+2] ;
  }

  /* compute the inverse of A       */

  div = a11 * ( a22*a33 - a23*a23 ) + 
    a12 * ( a13*a23 - a12*a33 ) +
    a13 * ( a12*a23 - a13*a22 ) ;

  /* Check for possible trouble.  Either colinear points or
  all pts. lie near to a plane which passes through origin resulting
  in the A matrix to be singular, e.g., if only three pts. then
  third row is linear combo. of top two rows.
  */

  if ( fabs(div) < ltol )
  {

    /* Find first two segments which are non-colinear.
    Use these to form the plane eqn.
    */
    double  seg1[3], seg2[3], cross[3];

    for( i = 0; i < num_pts - 2; i++ )
    {
      MshSurfDiffVector( &pts[(i+1)*3], &pts[i*3], seg1 ) ;
      MshSurfDiffVector( &pts[(i+2)*3], &pts[(i+1)*3], seg2 ) ;
      MshSurfCrossProd( seg1, seg2, cross ) ;
      if (!(fabs(cross[0]) < ltol && fabs(cross[1]) < ltol && fabs(cross[2]) < ltol)) 
      {
        MshSurfGenPlanEqn( &pts[i*3], &pts[(i+1)*3], &pts[(i+2)*3], plan_eqn ) ;
        return 1;
      }
    }
    return 0;
  }
  else 
  {
    /* compute the inverse of A coefficients      */

    b11 = ( a22*a33 - a23*a23 ) / div ;
    b12 = ( a13*a23 - a12*a33 ) / div ;
    b13 = ( a12*a23 - a13*a22 ) / div ;
    b22 = ( a11*a33 - a13*a13 ) / div ;
    b23 = ( a12*a13 - a11*a23 ) / div ;
    b33 = ( a11*a22 - a12*a12 ) / div ;

    /* compute the equation coefficients      */

    plan_eqn[0] = b11*sum_x + b12*sum_y + b13*sum_z ;
    plan_eqn[1] = b12*sum_x + b22*sum_y + b23*sum_z ;
    plan_eqn[2] = b13*sum_x + b23*sum_y + b33*sum_z ;

    /* normalize the equation coefficients      */

    norm = sqrt( plan_eqn[0] * plan_eqn[0] + plan_eqn[1] * plan_eqn[1] + 
      plan_eqn[2] * plan_eqn[2] ) ;

    plan_eqn[0] /= norm ;
    plan_eqn[1] /= norm ;
    plan_eqn[2] /= norm ;
    plan_eqn[3] = - 1.0 / norm ;
  }

  return 1;
}

/* --------------------------------------------------------------- */
int MshSurfLstSqrPlanMtx( int num_pts, double *pts, double mat[4][4] )
{
  int i, j;
  double plan_eq[4], tmp_axis[3];

  if (num_pts < 3)
    return 0;

  // plane equation
  if (!MshSurfLstSqrPlanEqn (num_pts, pts, plan_eq))
    return 0;

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

  return 1;
}

/* --------------------------------------------------------------- */
void MshSurfGenPlanEqn (double pt0[3], double pt1[3], double pt2[3], double plan_eqn[4]) 
{
  double vec1[3], vec2[3], normal[3];
  double length ;

  vec1[0] =  pt1[0] - pt0[0] ;
  vec1[1] =  pt1[1] - pt0[1] ;
  vec1[2] =  pt1[2] - pt0[2] ;
  vec2[0] =  pt2[0] - pt0[0] ;
  vec2[1] =  pt2[1] - pt0[1] ;
  vec2[2] =  pt2[2] - pt0[2] ;

  normal[0] = (vec1[1] * vec2[2]) - (vec1[2] * vec2[1]) ; 
  normal[1] = (vec1[2] * vec2[0]) - (vec1[0] * vec2[2]) ; 
  normal[2] = (vec1[0] * vec2[1]) - (vec1[1] * vec2[0]) ; 

  length = sqrt (normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

  plan_eqn[0] = normal[0] / length ;
  plan_eqn[1] = normal[1] / length ;
  plan_eqn[2] = normal[2] / length ;
  plan_eqn[3] = (-plan_eqn[0]*pt0[0] - plan_eqn[1]*pt0[1] - plan_eqn[2]*pt0[2]) ;

  return;
}


/* --------------------------------------------------------------- */
void MshSurfCrossProd (double *u, double *v, double *w)
{
  w[0] = u[1]*v[2]-u[2]*v[1];
  w[1] = u[2]*v[0]-u[0]*v[2];
  w[2] = u[0]*v[1]-u[1]*v[0];
}

/* --------------------------------------------------------------- */
void MshSurfDiffVector (double *a, double *b, double *c)
{
  c[0] = a[0]-b[0];
  c[1] = a[1]-b[1];
  c[2] = a[2]-b[2];
}

/* --------------------------------------------------------------- */
void MshSurfCrossProdNorm (double *u, double *v, double *w)
{
  double tmp[3], length;
  tmp[0] = u[1]*v[2]-u[2]*v[1];
  tmp[1] = u[2]*v[0]-u[0]*v[2];
  tmp[2] = u[0]*v[1]-u[1]*v[0];
  length = MshSurfCompr (tmp);
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
double MshSurfCompr (double *u)
{
  return (sqrt (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]));
}

/* --------------------------------------------------------------- */
double MshSurfSquareCompr (double *u)
{
  return (sqrt (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]));
}

/* --------------------------------------------------------------- */
void MshSurfCoord_xy (double x, double y, double z,
                             double BaseMatrix[4][4], double *p, double *h)
{
  double dx = x - BaseMatrix[3][0];
  double dy = y - BaseMatrix[3][1];
  double dz = z - BaseMatrix[3][2];
  p[0] = BaseMatrix[0][0]*dx + BaseMatrix[0][1]*dy + BaseMatrix[0][2]*dz;
  p[1] = BaseMatrix[1][0]*dx + BaseMatrix[1][1]*dy + BaseMatrix[1][2]*dz;
  *h   = BaseMatrix[2][0]*dx + BaseMatrix[2][1]*dy + BaseMatrix[2][2]*dz;
}

/* --------------------------------------------------------------- */
void MshSurfCoord_xyz (double x, double y, double BaseMatrix[4][4],
                              double *p)
{
  double dx = BaseMatrix[0][0]*x + BaseMatrix[1][0]*y + BaseMatrix[2][0]*0.0;
  double dy = BaseMatrix[0][1]*x + BaseMatrix[1][1]*y + BaseMatrix[2][1]*0.0;
  double dz = BaseMatrix[0][2]*x + BaseMatrix[1][2]*y + BaseMatrix[2][2]*0.0;
  p[0] = BaseMatrix[3][0] + dx;
  p[1] = BaseMatrix[3][1] + dy;
  p[2] = BaseMatrix[3][2] + dz;
}


/* --------------------------------------------------------------- */
void MshSurfTransf_xyx( double *p_in, double BaseMatrix[4][4], double *p_out )
{
  double dx = p_in[0] - BaseMatrix[3][0];
  double dy = p_in[1] - BaseMatrix[3][1];
  double dz = p_in[2] - BaseMatrix[3][2];
  p_out[0] = BaseMatrix[0][0]*dx + BaseMatrix[0][1]*dy + BaseMatrix[0][2]*dz;
  p_out[1] = BaseMatrix[1][0]*dx + BaseMatrix[1][1]*dy + BaseMatrix[1][2]*dz;
  p_out[2] = BaseMatrix[2][0]*dx + BaseMatrix[2][1]*dy + BaseMatrix[2][2]*dz;
}

/* --------------------------------------------------------------- */
void MshSurfInvTransf_xyz( double *p_in, double BaseMatrix[4][4], double *p_out )
{
  double x = p_in[0];
  double y = p_in[1];
  double z = p_in[2];
  double dx = BaseMatrix[0][0]*x + BaseMatrix[1][0]*y + BaseMatrix[2][0]*z;
  double dy = BaseMatrix[0][1]*x + BaseMatrix[1][1]*y + BaseMatrix[2][1]*z;
  double dz = BaseMatrix[0][2]*x + BaseMatrix[1][2]*y + BaseMatrix[2][2]*z;
  p_out[0] = BaseMatrix[3][0] + dx;
  p_out[1] = BaseMatrix[3][1] + dy;
  p_out[2] = BaseMatrix[3][2] + dz;
}


/* --------------------------------------------------------------- */
double MshSurfArea2D( int nNodes, double *boundpts2d )
{
  int i, j;
  double area = 0.0;

  for (i = 0, j = 1; i < nNodes; ++i, ++j)
  {
    if (i == nNodes-1)
    {
      area = area + boundpts2d[i*2+0]*boundpts2d[1] - 
        boundpts2d[0]*boundpts2d[i*2+1];
    }
    else
    {
      area = area + boundpts2d[i*2+0]*boundpts2d[j*2+1] - 
        boundpts2d[j*2+0]*boundpts2d[i*2+1];
    }
  }

  return area;
}

/* --------------------------------------------------------------- */
double MshSurfArea2DEdges   (int nNodes, double *boundpts2d, 
                             int nedge, int *edges)
{
  int i, id_i, id_j;
  double area = 0.0;

  for (i = 0; i < nedge; ++i)
  {
    id_i  = edges[i*2+0];
    id_j  = edges[i*2+1];
    area += (boundpts2d[id_i*2+0]*boundpts2d[id_j*2+1] - 
             boundpts2d[id_j*2+0]*boundpts2d[id_i*2+1]);
  }

  return area;
}


/* --------------------------------------------------------------- */
double MshSurfAngleAboutPt( double a[3], double b[3], double c[3] )
{
  double v0[3], v1[3], l0, l1, cosang, angle;


  /* first the vectors from the pivot to the two other points.
  get the cross product of these vectors. */
  v0[0] = ( b[0] - a[0] );
  v0[1] = ( b[1] - a[1] );
  v0[2] = ( b[2] - a[2] );
  v1[0] = ( c[0] - a[0] );
  v1[1] = ( c[1] - a[1] );
  v1[2] = ( c[2] - a[2] );
  l0 =  MshSurfCompr (v0);
  l1 =  MshSurfCompr (v1);

  if (l0 == 0.0 || l1 == 0.0)
    return pi;

  cosang = (v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2]) / l0 / l1;

  if (cosang > 1.0)
    cosang = 1.0;
  if (cosang < -1.0)
    cosang = -1.0;


  angle = acos(cosang);

  if (angle < -1.0*pi || angle > pi)
	  return pi;

  return angle;
}

/* --------------------------------------------------------------- */
void MshSurfFindBorders(int nNodes, double *pts, int *nbordes, int **borders )
{
  int i, k;
  double pprev[3], pcurr[3], pnext[3];

  /*fill metrics structure with current id node and angle*/
  metricBs3D *metrics = (metricBs3D *)calloc( nNodes, sizeof( metricBs3D ) );
  for( i = 0; i < nNodes; i++ )
  {
    int curr = i;
    int next = (i+1)%nNodes;
    int prev;

    if( i == 0 ) 
      prev = nNodes-1;
    else         
      prev = i-1;

    for (k = 0; k < 3; ++k)
    {
      pprev[k] = pts[prev*3+k];
      pcurr[k] = pts[curr*3+k];
      pnext[k] = pts[next*3+k];
    }

    /* compute angle*/
    metrics[i].value = MshSurfAngleAboutPt( pcurr, pprev, pnext );
    metrics[i].id    = i;
    /*printf ("curr = %d , %f - %d-%d\n", curr, metrics[i].value, prev, next);*/
  }

  /* sort metrics*/
  qsort( metrics, nNodes, sizeof( metricBs3D ), compmetr );

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
static int 
compmetr( const void *c1, const void *c2 )
{
  metricBs3D *m1 = (metricBs3D *) c1;
  metricBs3D *m2 = (metricBs3D *) c2;
  if      (m1->value < m2->value) return -1;
  else if (m1->value > m2->value) return  1;
  else                            return  0;
}
