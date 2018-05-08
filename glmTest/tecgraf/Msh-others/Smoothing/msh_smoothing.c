/*
**  msh_smoothing.c
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "msh_smoothing.h"
#include "topsurfbin.h"
#include "mshsurf3d.h"


#define MAX_ADJ_EDGES  50
#define NUM_NEWTON_RAPHSON_INT  10

static double LAPLACIAN_FACTOR = 1.0;

/*
** ---------------------------------------------------------------
** Private functions:
*/
static void   SmoothLocalDenoisePoint (double *old_pos, double *normal, int nadj, 
                                       double *adj_coords, double *new_pos);
static void   LaplacianOperator       (double *old_pos, int nadj, 
                                       double *adj_coords, double *v_lap);
static void   LocalMeanCurvature      (double *old_pos, int nadj, 
                                       double *adj_coords, double *curvature);
static double StandardDeviation       (int n, double *x);
static double ComputeCotAngle         (double *p0, double *pi, double *pj);
#if 0
static double CoumputeArea            (double *p0, double *pi, double *pj);
#endif
static void   GetLaplacianVectors     (void *surf, double *p_pos, double *q_pos);
static void   GetCurvatureVectors     (void *surf, double *p_pos, double *v_curv);
static void   NormalizeVector         (double *vec, double *vnorm);
static double ComputeCosBetweenVector (double *vec1, double *vec2);
static void   getPositiveParall       (void *surf, int i, double coords[3]);
static void   getNegativeParall       (void *surf, int i, double coords[3]);

/* get gradient and Hessian of metric function */
static double getQuadMetric   (double p0[2], double p1[2], double p2[2], double p3[2]);
static double getTriangMetric (double p0[2], double p1[2], double p2[2]);


/*
** ---------------------------------------------------------------
** Public function:
*/

/* MshSurfSmooth
** --------------------------------------------------------------- */
int MshSurfSmooth (
int     n_node,       /* # of pts in the mesh                    (in) */
double  *coords,      /* coordinate array of the mesh        (in/out) */
int     n_elem,       /* number of elements generated            (in) */
int     *Conn,        /* elem.connectivity list of the mesh      (in) */
int     algorithm,    /* algorithm: 1-laplacian, 2 - mod. laplacian, 3 - taubin (in) */
int     num_steps,    /* number of steps                         (in) */
int     n_param,      /* number of parameters                    (in) */
double  *param        /* parameters                              (in) */
)
{
  int i;
  
  /* create internal mesh */
  void *surf = SurfTopInsertMesh (n_node, coords, n_elem, Conn, NULL);
  
  if (param == NULL)
    return 0;
  if (surf == NULL)
    return 0;
  if (n_param == 0)
    return 0;

  // 
  switch (algorithm)
  {
    case 1: /* Laplacian */
      SmoothSetLaplacianFactor (param[0]);
      for (i = 0; i < num_steps; i++)
        SmoothLaplacian (surf);
    break;

    case 2: /* Laplacian Using Origin Surface */
      SmoothSetLaplacianFactor (param[0]);
      /* for (i = 0; i < number_of_interactions; i++) */
        SmoothHC (surf, num_steps);
    break;

    case 3: /* Taubin */
      if (n_param != 2)
        return 0;
      
      for (i = 0; i < num_steps*2; ++i)
      {
        if (i % 2 == 0)
          SmoothSetLaplacianFactor (param[0]);
        else
          SmoothSetLaplacianFactor (-1.0 * param[1]);
        SmoothLaplacian (surf);
      }
    break;

#if 0
    case 4: /* Mean Curvature Laplacian */
      for (i = 0; i < number_of_interactions; i++) 
        SmoothMeanCurvFlow (surf);
    break;

    case 5: /* Bilateral Mesh Denoising */
      for (i = 0; i < number_of_interactions; i++) 
        SmoothDenoisePoint (surf);
    break;

    case 6: /* No-Interative Mesh Smoothing */
    break;

    case 7: /* Use local Bezier Surface */
      SmoothLocalBezierSurf (surf);
    break;

    case 8: /* Use local Spline Surface */
      SmoothLocalSplineSurf (surf);
    break;
#endif    
    
  }

  // update nodes coords.
  for (i = 0; i < n_node; i++)
  {
    double coord[3];
    SurfTopGetCoordNode (surf, i, coord);
    memcpy(&(coords[i*3]), coord, 3*sizeof(double));
  }

  return 1;
}



/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
void SmoothSetLaplacianFactor (double value)
{
  if (value <= 0.0)
    LAPLACIAN_FACTOR = 1.0;

  LAPLACIAN_FACTOR = value;
}


/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
void SmoothLaplacian (void *surf)
{
  int    n_node, i, k;
  double old_coords[3], new_coords[3], *lap_vec;
  
  if (surf == NULL)
    return;

  n_node = SurfTopNumNodes (surf);
  lap_vec = (double *) calloc (n_node*3, sizeof (double));

  /* Regular Laplacian Operator */
  GetLaplacianVectors (surf, NULL, lap_vec);

  /* update nodes coords */
  for (i = 0; i < n_node; i++)
  {
    if (!SurfTopIsBoundaryNode (surf, i))
    {
      /* compute new coodinate */
      SurfTopGetCoordNode (surf, i, old_coords);
      for (k = 0; k < 3; k++)
        new_coords[k] = old_coords[k] + LAPLACIAN_FACTOR * lap_vec[i*3+k];

      /* update in the topological surface */
      SurfTopSetCoordNode  (surf, i, new_coords);
    }
  }

  SurfTopUdateNormalNodes (surf);

  free (lap_vec);
}

/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
void SmoothBiLaplacian (void *surf)
{
  int    n_node, i, k;
  double old_coords[3], new_coords[3], *first_lap_vec, *second_lap_vec;
  
  if (surf == NULL)
    return;

  n_node = SurfTopNumNodes (surf);
  first_lap_vec  = (double *) calloc (n_node*3, sizeof (double));
  second_lap_vec = (double *) calloc (n_node*3, sizeof (double));

  /* First Regular Laplacian Operator */
  GetLaplacianVectors (surf, NULL, first_lap_vec);

  /* Second Regular Laplacian Operator */
  GetLaplacianVectors (surf, first_lap_vec, second_lap_vec);

  /* update nodes coords */
  for (i = 0; i < n_node; i++)
  {
    if (!SurfTopIsBoundaryNode (surf, i))
    {
      /* compute new coordinate */
      SurfTopGetCoordNode (surf, i, old_coords);
      for (k = 0; k < 3; k++)
        new_coords[k] = old_coords[k] + LAPLACIAN_FACTOR * second_lap_vec[i*3+k];

      /* update in the topological surface */
      SurfTopSetCoordNode  (surf, i, new_coords);
    }
  }

  SurfTopUdateNormalNodes (surf);

  free (first_lap_vec);
  free (second_lap_vec);
}


/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
void SmoothDenoisePoint (void *surf)
{
  int    n_node, nadj_node, *adj_nodes, i, k;
  double normal[3], coords[3], *adj_coords, *new_coords, old_coords[3];

  if (surf == NULL)
    return;

  n_node = SurfTopNumNodes (surf);

  adj_coords = (double *) calloc (MAX_ADJ_EDGES*3, sizeof (double));
  adj_nodes  = (int *) calloc (MAX_ADJ_EDGES, sizeof (int));
  new_coords = (double *) calloc (n_node*3, sizeof (double));

  for (i = 0; i < n_node; i++)
  {
    nadj_node = SurfTopNumAdjNodeNode (surf, i);
    SurfTopGetCoordNode   (surf, i, coords);
    memcpy (old_coords, coords, 3 * sizeof (double));
    SurfTopGetNormalNode  (surf, i, normal);
    /* get adj ids of current node */
    SurfTopAdjNodeNode (surf, i, adj_nodes);
    /* get adj coords of current node */
    for (k = 0; k < nadj_node; k++)
    {
      SurfTopGetCoordNode   (surf, adj_nodes[k], coords);
      memcpy (&(adj_coords[k*3]), coords, 3 * sizeof (double));
    }
    /* smooth */
    SmoothLocalDenoisePoint (old_coords, normal, nadj_node, adj_coords, &(new_coords[i*3+0]));
  }

  /* update nodes coords */
  for (i = 0; i < n_node; i++)
  {
    if (!SurfTopIsBoundaryNode (surf, i))
      SurfTopSetCoordNode  (surf, i, &(new_coords[i*3+0]));
  }

  SurfTopUdateNormalNodes (surf);
  
  free (adj_coords);
  free (adj_nodes);
  free (new_coords);
}


/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
void SmoothHC (void *surf, int num_interactions)
{
  int    n_node, i, k, m;
  double *bi_coords, *ori_coords, *lap_vec, *curr_coords;

  /* parameter used in this algorithm */
  double alpha = 0.2;
  double beta  = 0.6;

  if (surf == NULL)
    return;

  n_node = SurfTopNumNodes (surf);

  curr_coords = (double *) calloc (n_node*3, sizeof (double));
  lap_vec     = (double *) calloc (n_node*3, sizeof (double));
  bi_coords   = (double *) calloc (n_node*3, sizeof (double));
  ori_coords  = (double *) calloc (n_node*3, sizeof (double));

  /* get original nodes */
  for (i = 0; i < n_node; i++)
  {
    double coords[3];
    SurfTopGetCoordNode   (surf, i, coords);
    memcpy (&(ori_coords[i*3]), coords, 3 * sizeof (double));
    memcpy (&(curr_coords[i*3]), coords, 3 * sizeof (double));
  }

  /* until the number of interactions */
  for (m = 0; m < num_interactions; m++)
  {

    /* Regular Laplacian Operator */
    GetLaplacianVectors (surf, curr_coords, lap_vec);

    /* compute bi coordinates */
    for (i = 0; i < n_node; i++)
    {
      /* compute bi coordinates */
      for (k = 0; k < 3; k++)
      {
        bi_coords[i*3+k] = (curr_coords[i*3+k] + LAPLACIAN_FACTOR * lap_vec[i*3+k]) -
                           (alpha * ori_coords[i*3+k] + (1-alpha) * curr_coords[i*3+k]);
        curr_coords[i*3+k] = curr_coords[i*3+k] + LAPLACIAN_FACTOR * lap_vec[i*3+k];
      }
    }

    /* Regular Laplacian Operator */
    GetLaplacianVectors (surf, bi_coords, lap_vec);

    /* compute bi coordinates */
    for (i = 0; i < n_node; i++)
    {
      /* compute bi coordinates */
      for (k = 0; k < 3; k++)
        curr_coords[i*3+k] = curr_coords[i*3+k] - (beta * bi_coords[i*3+k] + 
                            (1-beta) * (bi_coords[i*3+k] + LAPLACIAN_FACTOR * lap_vec[i*3+k]));
    }

  }

  /* update coords  nodes */
  for (i = 0; i < n_node; i++)
  {
    if (!SurfTopIsBoundaryNode (surf, i))
      SurfTopSetCoordNode  (surf, i, &(curr_coords[i*3]));
  }

  /* update normal nodes */
  SurfTopUdateNormalNodes (surf);

  free (curr_coords);
  free (lap_vec);
  free (bi_coords);
  free (ori_coords);
}


/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
void SmoothMeanCurvFlow (void *surf)
{
  int    n_node, i, k;
  double new_coords[3], old_coords[3], res_vec[3];
  double *lap_vec, *cur_vec, m[3];
  double cos_a, eps = 0.1, H;
  
  if (surf == NULL)
    return;
  
  n_node = SurfTopNumNodes (surf);

  lap_vec = (double *) calloc (n_node*3, sizeof (double));
  cur_vec = (double *) calloc (n_node*3, sizeof (double));

  /* Regular Laplacian Operator */
  GetLaplacianVectors (surf, NULL, lap_vec);

  /* Curvature aproximation vector */
  GetCurvatureVectors (surf, NULL, cur_vec);


  /* update nodes coords */
  for (i = 0; i < n_node; i++)
  {
    if (!SurfTopIsBoundaryNode (surf, i))
    {
      /* compute the vector of moving */
      /* cos angle between lap_vec and curv_vec */
      NormalizeVector (&lap_vec[i*3], m);
      cos_a = ComputeCosBetweenVector (&lap_vec[i*3], &cur_vec[i*3]);
      H = sqrt (cur_vec[i*3+0]*cur_vec[i*3+0] + 
                cur_vec[i*3+1]*cur_vec[i*3+1] + 
                cur_vec[i*3+2]*cur_vec[i*3+2]);

      /* tests */
      if (cos_a > eps)
      {
        for (k = 0; k < 3; k++)
          res_vec[k] = H * m[k] / cos_a;
      }
      else if (cos_a < -1.0 * eps)
      {
        for (k = 0; k < 3; k++)
          res_vec[k] = 2.0*cur_vec[i*3+k] - H*m[k]/cos_a;
      }
      else 
      {
        res_vec[0] = res_vec[1] = res_vec[2] = 0.0;
      }


      /* compute new coodinate */
      SurfTopGetCoordNode (surf, i, old_coords);
      for (k = 0; k < 3; k++)
        new_coords[k] = old_coords[k] + LAPLACIAN_FACTOR * res_vec[k];

      /* update in the topological surface */
      SurfTopSetCoordNode  (surf, i, new_coords);
    }
  }

  SurfTopUdateNormalNodes (surf);

  free (lap_vec);
  free (cur_vec);
}


/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
void SmoothLocalBezierSurf (void *surf)
{
  int    n_node, i, size_net = 2, **net;
  double *lap_vec;
  
  if (surf == NULL)
    return;

  n_node = SurfTopNumNodes (surf);
  lap_vec = (double *) calloc (n_node*3, sizeof (double));

  net = (int **) calloc (size_net*2+1, sizeof (int *));
  for (i = 0; i < size_net*2+1; i++)
    net[i] = (int *) calloc (size_net*2+1, sizeof (int));

  /* update nodes coords */
  for (i = 0; i < n_node; i++)
  {

#if 0
    if (SurfTopGetNetNodes  (surf, i, size_net, net))
    {
      printf ("Net para id = %d\n", i);
      for (k = 0; k < size_net*2+1; k++)
      {
        for (m = 0; m < size_net*2+1; m++)
          printf ("%d\t", net[k][m]);
        printf ("\n");
      }

    }

    if (!SurfTopIsBoundaryNode (surf, i))
    {
      /* compute new coordinate */
      old_coords = SurfTopGetCoordNode (surf, i);
      for (k = 0; k < 3; k++)
        new_coords[k] = old_coords[k] + LAPLACIAN_FACTOR * lap_vec[i*3+k];

      /* update in the topological surface */
      SurfTopSetCoordNode  (surf, i, new_coords);
    }
#endif

}

  SurfTopUdateNormalNodes (surf);

  free (lap_vec);
}

/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
void SmoothLocalSplineSurf (void *surf)
{

}





/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
void SmoothParallelogram (void *surf, double f1, double f2)
{
  int    n_node, i;
  double new_coords[3], positive_coords[3], negative_coords[3], *smooth_coord;

  /* check if surface is valid */
  if (surf == NULL)
    return;

  /* check for quad elements using the first element */
  n_node = SurfTopNumNodes (surf);

  smooth_coord = (double *) calloc (n_node*3, sizeof (double));

  /* update nodes coords */
  for (i = 0; i < n_node; i++)
  {
    if (!SurfTopIsBoundaryNode (surf, i))
    {
      int nadj_node = SurfTopNumAdjNodeNode (surf, i);

      if ( nadj_node == 4)
      {
        /* get positive part of new coordinate */
        getPositiveParall (surf, i, positive_coords);
        /* get negative part of new coordinate */
        getNegativeParall (surf, i, negative_coords);

        /* compute new coordinate */
        smooth_coord[i*3+0] = f1 * positive_coords[0] - f2 * negative_coords[0];
        smooth_coord[i*3+1] = f1 * positive_coords[1] - f2 * negative_coords[1];
        smooth_coord[i*3+2] = f1 * positive_coords[2] - f2 * negative_coords[2];
      }
      else
      {
        int adj_nodes[5], k;
        double coords[3];

        new_coords[0] = new_coords[1] = new_coords[2] = 0.0;

        /* get adj ids of current node */

        SurfTopAdjNodeNode (surf, i, adj_nodes);

        /* get adj coords of current node */
        for (k = 0; k < nadj_node; k++)
        {
          SurfTopGetCoordNode   (surf, adj_nodes[k], coords);
          new_coords[0] += coords[0];
          new_coords[1] += coords[1];
          new_coords[2] += coords[2];
        }

        smooth_coord[i*3+0] = new_coords[0] / nadj_node;
        smooth_coord[i*3+1] = new_coords[1] / nadj_node;
        smooth_coord[i*3+2] = new_coords[2] / nadj_node;
      }
    }
  }

  /* update nodes coords */
  for (i = 0; i < n_node; i++)
  {
    if (!SurfTopIsBoundaryNode (surf, i))
    {
      /* update in the topological surface */
      SurfTopSetCoordNode  (surf, i, &(smooth_coord[i*3]));
    }
  }

  SurfTopUdateNormalNodes (surf);
  free (smooth_coord);
}


/*
** ---------------------------------------------------------------
** Private functions:
*/



/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
static void SmoothLocalDenoisePoint (double *old_pos, double *normal, int nadj, 
                                     double *adj_coords, double *new_pos)
{
  int i;
  double *t, *h, dx, dy, dz, sig_c = 0.0, sig_s;
  double  wc, ws, sum = 0.0, normalizer = 0.0;

  /* compute distance and offset */
  t = (double *) calloc (nadj, sizeof (double));
  h = (double *) calloc (nadj+1, sizeof (double));
  for (i = 0; i < nadj; i++)
  {
    dx = adj_coords[i*3+0] - old_pos[0];
    dy = adj_coords[i*3+1] - old_pos[1];
    dz = adj_coords[i*3+2] - old_pos[2];
    t[i] = sqrt (dx*dx + dy*dy + dz*dz);
    if (t[i] > sig_c)
      sig_c = t[i];
    h[i] = normal[0]*dx + normal[1]*dy + normal[2]*dz;
  }
  h[i] = 0.0; /* old node */
  sig_s = StandardDeviation (nadj+1, h);

  /* compute sum and normalizer */
  for (i = 0; i < nadj; i++)
  {
    wc = exp (-1.0*t[i]*t[i] / (2*sig_c*sig_c));
    if (h[i] == 0.0)
      ws = 1.0;
    else
      ws = exp (-1.0*h[i]*h[i] / (2*sig_s*sig_s));
    sum += (wc*ws*h[i]);
    normalizer += (wc*ws);
  }

  if (normalizer != 0.0)
  {
    for (i = 0; i < 3; i++)
      new_pos[i] = old_pos[i] + normal[i] * sum / normalizer;
  }
  else
  {
    for (i = 0; i < 3; i++)
      new_pos[i] = old_pos[i];
  }

}

/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
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


/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
static void LocalMeanCurvature (double *old_pos, int nadj, 
                                double *adj_coords, double *curvature)
{
  int    i, k, qi, qi_prev, qi_next;
  double area, cot_alpha, cot_beta;
  double Q_P[3] = {0.0, 0.0, 0.0};

  area = 0.0;
  for (i = 0; i < nadj; i++)
  {
    qi      = i;
    qi_next = (i+1)%nadj;
    qi_prev = (i+nadj-1)%nadj;
 
    cot_alpha = ComputeCotAngle (&adj_coords[qi_next*3], old_pos, &adj_coords[qi*3]);
    cot_beta  = ComputeCotAngle (&adj_coords[qi_prev*3], &adj_coords[qi*3], old_pos);
//    curr_area = CoumputeArea    (&adj_coords[qi_next*3], old_pos, &adj_coords[qi*3]);

    for (k = 0; k < 3; k++)
      Q_P[k] += ((cot_alpha + cot_beta) * (adj_coords[qi*3+k] - old_pos[k]));

//    area += curr_area; 
    area += (cot_alpha + cot_beta); 
  }

  /* divide by 4 times area */
  for (k = 0; k < 3; k++)
    curvature[k] = Q_P[k] / (1.0 * area);
//    Q_P[k] /= (4 * area);

}

/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
static double StandardDeviation (int n, double *x)
{
  int i;
  double mean, sum;

  if (n == 0)
    return 0;

  /* mean */
  mean = 0;
  for (i = 0; i < n; i++)
    mean += x[i];
  mean /= n;

  sum = 0;
  for (i = 0; i < n; i++)
    sum += ((x[i]-mean) * (x[i]-mean));

  return (sqrt(sum/(n-1)));
}


/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
static double ComputeCotAngle (double *p0, double *pi, double *pj)
{
  double vi[3], vj[3], di, dj, angle;
  int    k;

  for (k = 0; k < 3; k++)
  {
    vi[k] = pi[k] - p0[k];
    vj[k] = pj[k] - p0[k];
  }
  di = sqrt (vi[0]*vi[0] + vi[1]*vi[1] + vi[2]*vi[2]);
  dj = sqrt (vj[0]*vj[0] + vj[1]*vj[1] + vj[2]*vj[2]);

  if (di == 0.0 || dj == 0.0)
    return 0;

  angle = acos ((vi[0]*vj[0] + vi[1]*vj[1] + vi[2]*vj[2]) / (di * dj));

  if (angle == 0.0)
    return 0.0;

  if (angle == acos(-1.0) / 2.0)
    return 0.0;

  return (1.0 / tan (angle));
}

#if 0
/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
static double CoumputeArea (double *p0, double *pi, double *pj)
{
  double u[3], v[3], w[3], area;
  int    k;

  for (k = 0; k < 3; k++)
  {
    u[k] = pi[k] - p0[k];
    v[k] = pj[k] - p0[k];
  }
  w[0] = u[1]*v[2]-u[2]*v[1];
  w[1] = u[2]*v[0]-u[0]*v[2];
  w[2] = u[0]*v[1]-u[1]*v[0];
  area = (w[0]*w[0] + w[1]*w[1] + w[2]*w[2]) * 0.5;

  if (area < 0.0)
    area *= -1.0;

  return area;
}
#endif

/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
static void GetLaplacianVectors (void *surf, double *p_pos, double *v_lap)
{
  int    n_node, nadj_node, *adj_nodes, i, k;
  double coords[3], *adj_coords, old_coords[3];
  
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

    /* nunber of adj nodes */
    nadj_node = SurfTopNumAdjNodeNode (surf, i);

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


/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
static void GetCurvatureVectors (void *surf, double *p_pos, double *v_curv)
{
  int    n_node, nadj_node, *adj_nodes, i, k;
  double coords[3], *adj_coords, old_coords[3];
  
  if (surf == NULL)
    return;
  
  n_node = SurfTopNumNodes (surf);

  adj_coords = (double *) calloc (MAX_ADJ_EDGES*3, sizeof (double));
  adj_nodes  = (int *) calloc (MAX_ADJ_EDGES, sizeof (int));
  memset (v_curv, 0, n_node * 3 * sizeof (double));

  for (i = 0; i < n_node; i++)
  {

    if (SurfTopIsBoundaryNode (surf, i))
    {
      v_curv[i*3+0] = v_curv[i*3+1] = v_curv[i*3+2] = 0.0;
      continue;
    }

    /* nunber of adj nodes */
    nadj_node = SurfTopNumAdjNodeNode (surf, i);

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

#if 0
    /* get adj ids of current node */
    if (!SurfTopOrientAdjNodeNode (surf, i, adj_nodes))
      continue;
#endif

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
    LocalMeanCurvature (old_coords, nadj_node, adj_coords, &(v_curv[i*3]));
  }

  free (adj_coords);
  free (adj_nodes);
}

/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
static void NormalizeVector (double *vec, double *vnorm)
{
  double len;
  len = sqrt (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  if (len == 0.0)
  {
    vnorm [0] = vnorm [1] = vnorm [2] = 0.0;
    return;
  }
  vnorm [0] = vec[0] / len;
  vnorm [1] = vec[1] / len;
  vnorm [2] = vec[2] / len;
}

/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
static double ComputeCosBetweenVector (double *vec1, double *vec2)
{
  double v1n[3], v2n[3];
  NormalizeVector (vec1, v1n);
  NormalizeVector (vec2, v2n);
  return (v1n[0]*v2n[0] + v1n[1]*v2n[1] + v1n[2]*v2n[2]);
}

/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
void getNegativeParall( void *surf, int id, double coords[3] )
{
  int ne, *elem, *conn, i, j;
  int prev, curr, next;

  SurfTopAdjElemToNode (surf, id, &ne, &elem);

  /* checks */
  if (ne == 0)
    return;

  /* check for quad elements */
  for (i = 0; i < ne; i++)
  {
    if (SurfTopGetElemNNodes (surf, elem[i]) != 4)
    {
      free(elem);
      return;
    }
  }

  /* */
  coords[0] = coords[1] = coords[2] = 0.0;
  for (i = 0; i < ne; i++)
  {
    conn = SurfTopGetElemConn (surf, elem[i]);
    for (j = 0; j < 4; ++j)
    {
      prev = conn[(j+3)%4];
      curr = conn[j];
      next = conn[(j+1)%4];
      if (prev != id && curr != id && next != id)
      {
        double curr_coord[3];
        SurfTopGetCoordNode (surf, curr, curr_coord);
        coords[0] += curr_coord[0];
        coords[1] += curr_coord[1];
        coords[2] += curr_coord[2];
      }
    }
  }

  /* average */
  coords[0] /= ne;
  coords[1] /= ne;
  coords[2] /= ne;
 

  free (elem);

}

/* ------------------------------------------------------- */
/* ------------------------------------------------------- */
void getPositiveParall( void *surf, int id, double coords[3] )
{
  int nn, i, *adj_nodes;
  double div;

  nn = SurfTopNumAdjNodeNode (surf, id);
  /* checks */
  if (nn == 0)
    return;

  adj_nodes = (int *) calloc (nn, sizeof(int));
  SurfTopAdjNodeNode (surf, id, adj_nodes);


  coords[0] = coords[1] = coords[2] = 0.0;
  for (i = 0; i < nn; i++)
  {
    double curr_coord[3];
    SurfTopGetCoordNode (surf, adj_nodes[i], curr_coord);    
    coords[0] += curr_coord[0];
    coords[1] += curr_coord[1];
    coords[2] += curr_coord[2];
  }
  div = nn * 1.0;

  /* average */
  coords[0] /= div;
  coords[1] /= div;
  coords[2] /= div;
}

/***************************************************************************/
static int getG_H_Quad( double p[4][2], int pos, double G[2], double H[2][2] )
{
  double x1 = p[(pos+2)%4][0];
  double y1 = p[(pos+2)%4][1];
  double x2 = p[(pos+3)%4][0];
  double y2 = p[(pos+3)%4][1];
  double x3 = p[pos][0];
  double y3 = p[pos][1];
  double x4 = p[(pos+1)%4][0];
  double y4 = p[(pos+1)%4][1];
  double p_1 = (-x2+2*x3-x4);
  double p_2 = (-y2+y4);
  double p_3 = (-y2+2*y3-y4);
  double p_4 = (-x4+x2);
  double num, denom;

  num   = -x2*y1 - x1*y4 + y2*x1 + y1*x4 - x4*y3 - x3*y2 + y4*x3 + y3*x2;
  denom =  x2*x2 - x2*x1 + x1*x1 + y2*y2 - y2*y1 + y1*y1 + x4*x4 - x4*x1 + 
    y4*y4 - y4*y1 - x2*x3 + x3*x3 - y2*y3 + y3*y3 - x4*x3 - y4*y3;

  if (denom == 0.0)
  {
    G[0] = G[1] = 0;
    H[0][0] =  H[1][1] =  1.0;
    H[0][1] =  H[1][0] =  0.0;
    return 1;
  }

  G[0] += ((-y2+y4) / denom) - ((num * (-x2+2*x3-x4)) / (denom*denom));
  G[1] += ((-x4+x2) / denom) - ((num * (-y2+2*y3-y4)) / (denom*denom));


  H[0][0] += ((2 * num * p_1 * p_1) / (denom*denom*denom)) -
    ((2 * p_1 * p_2) / (denom*denom)) - ((2 * num) / (denom*denom));
  H[1][1] += ((2 * num * p_3 * p_3) / (denom*denom*denom)) -
    ((2 * p_3 * p_4) / (denom*denom)) - ((2 * num) / (denom*denom));
  H[0][1] += ((2 * num * p_1 * p_3) / (denom*denom*denom)) -
    ((p_2 * p_3) / (denom*denom)) - ((p_1 * p_4) / (denom*denom));
  H[1][0] = H[0][1];

  return 1;
}

/***************************************************************************/
static int getG_H_Triang( double p[4][2], int pos, double G[2], double H[2][2] )
{
  double x1 = p[(pos+1)%3][0];
  double y1 = p[(pos+1)%3][1];
  double x2 = p[(pos+2)%3][0];
  double y2 = p[(pos+2)%3][1];
  double x3 = p[pos][0];
  double y3 = p[pos][1];

  double A_ = y1 - y2;
  double B_ = x2 - x1;
  /*double C_ = y2*x1 - x2*y1;*/
  double E_ = -x1 - x2;
  double F_ = -y1 - y2;
  /*double G_ = x1*x1 + x2*x2 + y1*y1 + y2*y2 - x2*x1 - y2*y1;*/
  double p_x = 2*x3 + E_;
  double p_y = 2*y3 + F_;

  double num   = x2*y3-x2*y1-x1*y3-y2*x3+y2*x1+y1*x3;
  double denom = x2*x2-x2*x1+x1*x1+y2*y2-y2*y1+y1*y1+x3*x3-x3*x1+y3*y3-y3*y1-x3*x2-y3*y2;

  if (denom == 0.0)
  {
    G[0] = G[1] = 0;
    H[0][0] =  H[1][1] =  1.0;
    H[0][1] =  H[1][0] =  0.0;
    return 1;
  }

  G[0] += (A_ / denom) - ((num * p_x) / (denom*denom));
  G[1] += (B_ / denom) - ((num * p_y) / (denom*denom));


  H[0][0] += ((2 * num * p_x * p_x) / (denom*denom*denom)) -
    ((2 * A_ * p_x) / (denom*denom)) - ((2 * num) / (denom*denom));
  H[1][1] += ((2 * num * p_y * p_y) / (denom*denom*denom)) -
    ((2 * B_ * p_y) / (denom*denom)) - ((2 * num) / (denom*denom));
  H[0][1] += ((2 * num * p_x * p_y) / (denom*denom*denom)) -
    ((A_ * p_y) / (denom*denom)) - ((B_ * p_x) / (denom*denom));
  H[1][0] = H[0][1];

  return 1;
}

/***************************************************************************/
double getQuadMetric( double p0[2], double p1[2], double p2[2], double p3[2] )
{
  double x1 = p0[0];
  double y1 = p0[1];
  double x2 = p1[0];
  double y2 = p1[1];
  double x3 = p2[0];
  double y3 = p2[1];
  double x4 = p3[0];
  double y4 = p3[1];

  double num, denom;

  num   = -x2*y1 - x1*y4 + y2*x1 + y1*x4 - x4*y3 - x3*y2 + y4*x3 + y3*x2;
  denom =  x2*x2 - x2*x1 + x1*x1 + y2*y2 - y2*y1 + y1*y1 + x4*x4 - x4*x1 + 
           y4*y4 - y4*y1 - x2*x3 + x3*x3 - y2*y3 + y3*y3 - x4*x3 - y4*y3;

  if (denom == 0.0)
    return 0.0;

  return (num/denom);
}


/***************************************************************************/
double getTriangMetric( double p0[2], double p1[2], double p2[2] )
{
  double x1 = p0[0];
  double y1 = p0[1];
  double x2 = p1[0];
  double y2 = p1[1];
  double x3 = p2[0];
  double y3 = p2[1];

  double num, denom;

  num   = x2*y3-x2*y1-x1*y3-y2*x3+y2*x1+y1*x3;
  denom = x2*x2-x2*x1+x1*x1+y2*y2-y2*y1+y1*y1+x3*x3-x3*x1+y3*y3-y3*y1-x3*x2-y3*y2;

  if (denom == 0.0)
    return 0.0;

  return (num/denom);
}

//////////////////////////////////////////////////////////////////////////
static void invert2_2 (double M[2][2], double MI[2][2])
{
  double denom = M[0][0]*M[1][1]-M[1][0]*M[0][1];
  MI[0][0] =  M[1][1] / denom;
  MI[0][1] = -M[1][0] / denom;
  MI[1][0] = -M[0][1] / denom; 
  MI[1][1] =  M[0][0] / denom;
}

//////////////////////////////////////////////////////////////////////////
static void multiply_M_V (double M[2][2], double v[2], double R[2])
{
  R[0] = M[0][0] * v[0] + M[0][1] * v[1];
  R[1] = M[1][0] * v[0] + M[1][1] * v[1];
}

/***************************************************************************/
static double getAdjElemMetric (void *surf, int id, double initial_coords[3])
{
  int ne, *elem, *conn, i, j;
  int n_nodes;
  double metric, curr_coord[3], pts[4][2];

  SurfTopAdjElemToNode (surf, id, &ne, &elem);

  /* checks */
  if (ne == 0)
    return 0.0;

  /* get initial metric */
  metric = 0.0;
  for (i = 0; i < ne; i++)
  {
    conn = SurfTopGetElemConn (surf, elem[i]);
    n_nodes = SurfTopGetElemNNodes(surf, elem[i]);
    
    for (j = 0; j < n_nodes; ++j)
    {
      if (conn[j] == id) /* we have the coordinate, so it is not necessary to get from structure */
      {
        pts[j][0] = initial_coords[0];
        pts[j][1] = initial_coords[1];
      }
      else
      {
        SurfTopGetCoordNode (surf, conn[j], curr_coord);
        pts[j][0] = curr_coord[0];
        pts[j][1] = curr_coord[1];
      }
    }

    if (n_nodes == 4)
      metric += getQuadMetric(pts[0], pts[1], pts[2], pts[3]);
    else
      metric += getTriangMetric(pts[0], pts[1], pts[2]);
  }

  free (elem);

  return metric;
}

/***************************************************************************/
void newPtsMaxmetric (void *surf, int id, double initial_coords[3], double new_coords[3])
{
  int ne, *elem, *conn, i, j, k;
  int num_nodes;
  double initial_metric, final_metric, curr_coord[3], pts[4][2];
  double grad[2], H[2][2], inv_H[2][2], delta[2];


  SurfTopAdjElemToNode (surf, id, &ne, &elem);

  /* set new coordinate as initial coordinate */
  new_coords[0] = initial_coords[0];
  new_coords[1] = initial_coords[1];
  new_coords[2] = initial_coords[2];

  /* checks */
  if (ne == 0)
    return;

  /* get initial metric */
  initial_metric = getAdjElemMetric (surf, id, new_coords); 


  /* step of interaction on newton-raphson */
  for (k = 0; k < NUM_NEWTON_RAPHSON_INT; k++)
  {
    grad[0] = grad[1] = H [0][0] = H [1][1] = H [0][1] = H [1][0] = 0.0;

    for (i = 0; i < ne; i++)
    {
      int pos = 0;
      conn = SurfTopGetElemConn (surf, elem[i]);
      num_nodes = SurfTopGetElemNNodes (surf, elem[i]);

      /* get node position */
      for (j = 0; j < num_nodes; ++j)
      {
        if (conn[j] == id)
        {
          pos = j;
          pts[j][0] = new_coords[0];
          pts[j][1] = new_coords[1];
        }
        else
        {
          SurfTopGetCoordNode (surf, conn[j], curr_coord);
          pts[j][0] = curr_coord[0];
          pts[j][1] = curr_coord[1];
        }
      }
      
      /* get gradient and Hessian */
      if (num_nodes == 4)
        getG_H_Quad (pts, pos, grad, H);
      else
        getG_H_Triang(pts, pos, grad, H);
    }

    /* compute delta */
    invert2_2    (H, inv_H);
    multiply_M_V (inv_H, grad, delta);
#if 1
    new_coords[0] = new_coords[0] - delta[0];
    new_coords[1] = new_coords[1] - delta[1];
#else
    new_coords[0] = new_coords[0] + grad[0];
    new_coords[1] = new_coords[1] + grad[1];

    if (fabs(grad[0])<0.000001 && fabs(grad[1])<0.000001)
      break;
#endif
  }
  
#if 0
  printf ("Gradiente = %lf, %lf\n", grad[0], grad[1]);
#endif

  /* get final metric */
  final_metric = getAdjElemMetric (surf, id, new_coords); 

  if (final_metric < initial_metric)
  {
    printf ("Metrics = %lf, %lf\n", initial_metric, final_metric);
    new_coords[0] = initial_coords[0];
    new_coords[1] = initial_coords[1];
    new_coords[2] = initial_coords[2];
  }

  /* */

    free (elem);
}


/***************************************************************************/
void Smooth2D_Optim ( void *surf )
{
  int    n_node, i, k;
  double old_coords[3], initial_coords[3], *lap_vec;

  if (surf == NULL)
    return;

  n_node = SurfTopNumNodes (surf);
  lap_vec = (double *) calloc (n_node*3, sizeof (double));

  /* Regular Laplacian Operator */
  GetLaplacianVectors (surf, NULL, lap_vec);

  /* update nodes coords */
  for (i = 0; i < n_node; i++)
  {
    if (!SurfTopIsBoundaryNode (surf, i))
    {
      /* compute initial coordinate  */
      SurfTopGetCoordNode (surf, i, old_coords);
      for (k = 0; k < 3; k++)
        initial_coords[k] = old_coords[k]; /* + LAPLACIAN_FACTOR * lap_vec[i*3+k]; */

      newPtsMaxmetric (surf, i, initial_coords, &(lap_vec[i*3]));
    }
  }

  /* update nodes coords */
  for (i = 0; i < n_node; i++)
  {
    if (!SurfTopIsBoundaryNode (surf, i))
      SurfTopSetCoordNode  (surf, i, &(lap_vec[i*3]));
  }

  SurfTopUdateNormalNodes (surf);

  free (lap_vec);

}








