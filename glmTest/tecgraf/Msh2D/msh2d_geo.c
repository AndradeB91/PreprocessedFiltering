/*
** ----------------------------------------------------------------------
**
** msh2D_geo.c - File for auxiliar 2D meshing routine. 
**
** ----------------------------------------------------------------------
**
** Created:      June-2010      Antonio C.O. Miranda
**
** ----------------------------------------------------------------------
*/

#include <stdlib.h>
#include <string.h> // memcpy
#include <math.h>

#include "msh2d_geo.h"


typedef struct
{
  int id;
  double value;
} metricBs2D;

static int compmetr( const void *c1, const void *c2 );



/*
//////////////////////////////////////////////////////////////////////////
*/
static int    Msh2DGetMappBilDirection4_( int *bd, int nPtsLoop, int *m, int *n );
static int    Msh2DGetMappBilDirection5_( int *bd, int nPtsLoop, int *m, int *n );
static int    Msh2DGetMappBilDirection6_( int *bd, int nPtsLoop, int *m, int *n );
static int    Msh2DGetMappTriDirection3_( int *bd, int nPtsLoop, int *m);
static int    Msh2DGetMappTriDirection4_( int *bd, int nPtsLoop, int *m);
static int    Msh2DGetMappTriDirection5_( int *bd, int nPtsLoop, int *m);

static int    compInt ( const void *c1, const void *c2 );


/* --------------------------------------------------------------- */
double Msh2DCrossProd (double *u, double *v, double *w)
{
  return (u[0]*v[1]-u[1]*v[0]);
}

/* --------------------------------------------------------------- */
void Msh2DDiffVector (double *a, double *b, double *c)
{
  c[0] = a[0]-b[0];
  c[1] = a[1]-b[1];
}


/* --------------------------------------------------------------- */
double Msh2DCompr (double *u)
{
  return (sqrt (u[0]*u[0] + u[1]*u[1]));
}

/* --------------------------------------------------------------- */
double Msh2DSquareCompr (double *u)
{
  return (sqrt (u[0]*u[0] + u[1]*u[1]));
}


/* --------------------------------------------------------------- */
double Msh2DArea2D( int nNodes, double *boundpts2d )
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
double Msh2DAngleAboutPt( double a[2], double b[2], double c[2] )
{
  double v0[2], v1[2], l0, l1, cosang, angle;

  /* first the vectors from the pivot to the two other points.
  get the cross product of these vectors. */
  v0[0] = ( b[0] - a[0] );
  v0[1] = ( b[1] - a[1] );
  v1[0] = ( c[0] - a[0] );
  v1[1] = ( c[1] - a[1] );
  l0 =  Msh2DCompr (v0);
  l1 =  Msh2DCompr (v1);

  cosang = (v0[0]*v1[0] + v0[1]*v1[1]) / l0 / l1;

  if (cosang > 1.0)
    cosang = 1.0;
  if (cosang < -1.0)
    cosang = -1.0;


  angle = acos(cosang);
  /*printf("value = %f - %f\n", cosang, angle);*/

  return angle;
}

/* --------------------------------------------------------------- */
void Msh2DFindBorders(int nNodes, double *pts, int *nbordes, int **borders )
{
  int i, k;
  double pprev[2], pcurr[2], pnext[2];

  /*fill metrics structure with current id node and angle*/
  metricBs2D *metrics = (metricBs2D *)calloc( nNodes, sizeof( metricBs2D ) );
  for( i = 0; i < nNodes; i++ )
  {
    int curr = i;
    int next = (i+1)%nNodes;
    int prev;

    if( i == 0 ) 
      prev = nNodes-1;
    else         
      prev = i-1;

    for (k = 0; k < 2; ++k)
    {
      pprev[k] = pts[prev*2+k];
      pcurr[k] = pts[curr*2+k];
      pnext[k] = pts[next*2+k];
    }

    /* compute angle*/
    metrics[i].value = Msh2DAngleAboutPt( pcurr, pprev, pnext );
    metrics[i].id    = i;
    /*printf ("curr = %d , %f - %d-%d\n", curr, metrics[i].value, prev, next);*/
  }

  /* sort metrics*/
  qsort( metrics, nNodes, sizeof( metricBs2D ), compmetr );

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
  metricBs2D *m1 = (metricBs2D *) c1;
  metricBs2D *m2 = (metricBs2D *) c2;
  if      (m1->value < m2->value) return -1;
  else if (m1->value > m2->value) return  1;
  else                            return  0;
}


/* --------------------------------------------------------------- */
int Msh2DAutomaticBilinearCornes (int np, double *pts, int *m, int *n)
{
  int nbordes, *borders, i, status;
  double *tmp_pts;

  // get border stones
  Msh2DFindBorders(np, pts, &nbordes, &borders );

  if (nbordes > 6)
    nbordes = 6;

  // order border stones id sequentially to loop 
  qsort (borders, nbordes, sizeof( int ), compInt );


  /* number of elements in each direction */
  switch(nbordes)
  {
  case 4:
    status = Msh2DGetMappBilDirection4_ (borders, np, m, n);
    break;

  case 5:
    status = Msh2DGetMappBilDirection5_ (borders, np, m, n);
    break;

  case 6:
    status = Msh2DGetMappBilDirection6_ (borders, np, m, n);
    break;

    //case :
    //    break;

  default:
    return 0;
    break;
  }

  if (!status)
    return 0;


  // align with first border stone
  tmp_pts = (double *) calloc (np*2, sizeof (double));
  memcpy(tmp_pts, pts, np*2*sizeof (double));
  for (i = 0; i < np; ++i)
  {
    pts[i*2+0] = tmp_pts[((borders[0]+i)%np)*2+0];
    pts[i*2+1] = tmp_pts[((borders[0]+i)%np)*2+1];
  }

  free (tmp_pts);
  free (borders);

  return 1;
}


/* --------------------------------------------------------------- */
int Msh2DAutomaticTrilinearCornes (int np, double *pts, int *m)
{
  int nbordes, *borders, i, status;
  double *tmp_pts;

  // get border stones
  Msh2DFindBorders(np, pts, &nbordes, &borders );

  if (nbordes > 5)
    nbordes = 5;

  // order border stones id sequentially to loop 
  qsort (borders, nbordes, sizeof( int ), compInt );


  /* number of elements in each direction */
  switch(nbordes)
  {
  case 3:
    status = Msh2DGetMappTriDirection3_ (borders, np, m);
    break;

  case 4:
    status = Msh2DGetMappTriDirection4_ (borders, np, m);
    break;

  case 5:
    status = Msh2DGetMappTriDirection5_ (borders, np, m);
    break;

    //case :
    //    break;

  default:
    return 0;
    break;
  }

  if (!status)
    return 0;


  // align with first border stone
  tmp_pts = (double *) calloc (np*2, sizeof (double));
  memcpy(tmp_pts, pts, np*2*sizeof (double));
  for (i = 0; i < np; ++i)
  {
    pts[i*2+0] = tmp_pts[((borders[0]+i)%np)*2+0];
    pts[i*2+1] = tmp_pts[((borders[0]+i)%np)*2+1];
  }

  free (tmp_pts);
  free (borders);

  return 1;
}


/* --------------------------------------------------------------- */
int Msh2DGetMappBilDirection4_( int *bd, int nPtsLoop, int *m, int *n )
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
int Msh2DGetMappBilDirection5_( int *bd, int nPtsLoop, int *m, int *n )
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
int Msh2DGetMappBilDirection6_( int *bd, int nPtsLoop, int *m, int *n )
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


/* --------------------------------------------------------------- */
int Msh2DGetMappTriDirection3_( int *bd, int nPtsLoop, int *m)
{
  int b[3] = {bd[1]-bd[0], bd[2]-bd[1], nPtsLoop-bd[2]+bd[0]};

  if (b[0] == b[1] && b[0] == b[2]) // match
  {
    *m = b[0] + 1;
    return 1;
  }
  return 0;
}

/* --------------------------------------------------------------- */
int Msh2DGetMappTriDirection4_( int *bd, int nPtsLoop, int *m)
{
  int i, div[3];
  int b[4] = {bd[1]-bd[0], bd[2]-bd[1], bd[3]-bd[2], nPtsLoop-bd[3]+bd[0]};

  for (i = 0; i < 3; ++i)
  {
    div[0] = b[i] + b[(i+1)%4];
    div[1] = b[(i+2)%4];
    div[2] = b[(i+3)%4];
    if (div[0] == div[1] && div[0] == div[2]) // match
    {
      *m = div[0] + 1;
      return 1;
    }
  }
  return 0;
}

/* --------------------------------------------------------------- */
int Msh2DGetMappTriDirection5_( int *bd, int nPtsLoop, int *m)
{
  int i, div[4];
  int b[5] = {bd[1]-bd[0], bd[2]-bd[1], bd[3]-bd[2], bd[4]-bd[3], nPtsLoop-bd[4]+bd[0]};

  /* first situation */
  for (i = 0; i < 4; ++i)
  {
    div[0] = b[i]       + b[(i+1)%5];
    div[1] = b[(i+2)%5] + b[(i+3)%5];
    div[2] = b[(i+4)%5];
    if (div[0] == div[1] && div[0] == div[2]) // match
    {
      *m = div[0] + 1;
      return 1;
    }
  }

  /* second situation */
  for (i = 0; i < 4; ++i)
  {
    div[0] = b[i]       + b[(i+1)%5];
    div[1] = b[(i+2)%5];
    div[2] = b[(i+3)%5] + b[(i+4)%5];
    if (div[0] == div[1] && div[0] == div[2]) // match
    {
      *m = div[0] + 1;
      return 1;
    }
  }

  return 0;
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
