/*
** ----------------------------------------------------------------------
**
** mshsurf3dmain.c - Surface mesh generator directly in 3D.
**
** ----------------------------------------------------------------------
**
** Created:      18-Jun-2007      Antonio C.O. Miranda
** ----------------------------------------------------------------------
**
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* #include <time.h> */
#include <string.h>

#include "mshsurf3d.h"
#include "surf3d.h"
#include "surf3d_auxfunc.h"
#include "surf3d_analt_auxfunc.h"
#include "topsurfbin.h"
#include "mshsurf_geo.h"

static void *support_mesh = NULL;
static  int (*MshSurfGetUVfromXYZ)   (double xyz[3], double tol, double uv[2]) = NULL;
static  int (*MshSurfGetXYZfromUV)   (double uv[2], double xyz[3])  = NULL;
static  int (*MshSurfGetGradVectors) (double uv[2], double Gu[3], double Gv[3]) = NULL;


int MshSurf_smoothOnProjection = 1; /* true */

/* type of support surface
** 1 - a mesh
** 2 - an analytical surface
*/
int MshSurfType = -1;

#if 0
static void write_neutralfile_teste(
char *name_file_teste,
int n_elem,
int n_node,
double *Coords,
int    *conn
);
#endif

/*
** static functions - definition
** ----------------------------------------------------------------------
*/
int MshSurfGetNormal (int n, double *coords, double size, double *normal);



/*
** public functions
** ----------------------------------------------------------------------
*/

/* First Function - Set the support surface */

/* The support surface is a mesh */
void MshSurfSetSupportMesh (
void    *surfmesh,    /* surface mesh using topology lib or the
                         inserting a mesh as bellow              (in) */
int     n_node,       /* # of pts in the mesh                    (in) */
double  *coords,      /* coordinate array of the mesh            (in) */
int     n_elem,       /* number of elements generated            (in) */
int     *Conn         /* elem.connectivity list of the mesh      (in) */
)
{
  /* write_neutralfile_teste ("c:\\temp\\out.nf", n_elem, n_node, coords, Conn); */
	if (surfmesh != NULL)
	  support_mesh = surfmesh;
	else
	  support_mesh = SurfTopInsertMesh (n_node, coords, n_elem, Conn, NULL);
  MshSurfType = 1;
}



/* The support surface is an analytical surface */
void MshSurfSetSupportAnalytical (
  int (*getUVfromXYZ)   (double xyz[3], double tol, double uv[2]),
  int (*getXYZfromUV)   (double uv[2], double xyz[3]),
  int (*getGradVectors) (double uv[2], double Gu[3], double Gv[3])
)
{
  MshSurfGetUVfromXYZ   = getUVfromXYZ;
  MshSurfGetXYZfromUV   = getXYZfromUV;
  MshSurfGetGradVectors = getGradVectors;

  AuxFuncSetCurrentSurface (getUVfromXYZ, getXYZfromUV, getGradVectors);

  MshSurfType = 2;
}


/* Second Function - Create a new mesh using the currente support mesh */
int MshSurf3D (
int     n_pts,            /* # of points                          (in) */
double  *bdry_pts,        /* coordinates of all points            (in) */
int     bound_edge,       /* # of boundary edges                  (in) */
int     inter_edge,       /* # of internal free edges             (in) */
int     *edges,           /* vector with the edges                (in) */
double  max_elm_size,     /* maximum size of elements             (in) */
int     curvature,        /* compute curvature ?                  (in) */
Surf3DMessFunc *mes_func, /* message function                     (in) */
int     *n_node,          /* counts # of pts for meshing         (out) */
double  **coords,         /* coordinate array used for meshing   (out) */
int     *n_elem,          /* number of elements generated        (out) */
int     **Conn            /* elem.connectivity list from meshing (out) */
)
{
	int    i, j, status, num_edge, *new_edges;
	double *normal, smallest_edge = max_elm_size, dx, dy, dz, len;

 /*
  long double cpu_time;
*/

  if (MshSurfType == -1) /* the support surface is not defined */
    return 0;

  /* open memory */
  num_edge = ( bound_edge + (inter_edge * 2) );
  new_edges = (int *) calloc ((num_edge*2),sizeof(int));

#if 0
  for (i = 0; i < n_pts; ++i)
  {
   printf ("%f %f %f\n", bdry_pts[i*3+0], bdry_pts[i*3+1], bdry_pts[i*3+2]);
  }
  printf ("\n");
#endif

  /* boundary edges */
  for (i = 0; i < bound_edge; i++)
  {
    /* insert one side of edge */
    new_edges[i*2 + 0] = edges[i*2 + 0];
    new_edges[i*2 + 1] = edges[i*2 + 1];
    dx = bdry_pts[new_edges[i*2+0]*3+0] - bdry_pts[new_edges[i*2+1]*3+0];
    dy = bdry_pts[new_edges[i*2+0]*3+1] - bdry_pts[new_edges[i*2+1]*3+1];
    dz = bdry_pts[new_edges[i*2+0]*3+2] - bdry_pts[new_edges[i*2+1]*3+2];
    len = sqrt (dx*dx + dy*dy + dz*dz);
    if (len < smallest_edge)
      smallest_edge  = len;
  }
  /* get the edges on the domain */
  for (i = bound_edge; i < (inter_edge+bound_edge); i++)
  {
    j = inter_edge + i;

    /* insert one side of edge */
    new_edges[i*2 + 0] = edges[i*2 + 0];
    new_edges[i*2 + 1] = edges[i*2 + 1];
    dx = bdry_pts[new_edges[i*2+0]*3+0] - bdry_pts[new_edges[i*2+1]*3+0];
    dy = bdry_pts[new_edges[i*2+0]*3+1] - bdry_pts[new_edges[i*2+1]*3+1];
    dz = bdry_pts[new_edges[i*2+0]*3+2] - bdry_pts[new_edges[i*2+1]*3+2];
    len = sqrt (dx*dx + dy*dy + dz*dz);
    if (len < smallest_edge)
      smallest_edge  = len;

    /* insert another side of edge */
    new_edges[j*2 + 0] = edges[i*2 + 1];
    new_edges[j*2 + 1] = edges[i*2 + 0];
  }

	/* set auxiliar functions */
  if (mes_func != NULL)
    mes_func ("Creating background mesh...");

/*
  cpu_time = clock( );
*/
  if (MshSurfType == 1) /* support mesh */
  {
    if (!AuxFuncSetCurrentMesh (support_mesh, max_elm_size, n_pts,
                                bdry_pts, num_edge, new_edges))
      return 0;
  }
  else if (MshSurfType == 2) /* support analytical surface */
  {
    if (!AuxFuncSetCurrentBoundary (max_elm_size, n_pts, bdry_pts,
                                    num_edge, new_edges))
      return 0;
  }
  else
    return 0;

  /* printf("\tSurface Mesh Generator Direct in 3D Space, by Antonio Miranda\n\n"); */

/*
  cpu_time = (clock( ) - cpu_time)/CLOCKS_PER_SEC;
  printf("\n\n");
  printf("\t\tCPU time background mesh............. %0.3f (s)\n", (double)cpu_time);
*/


  /* get normal points */
  normal = (double *) calloc (n_pts * 3, sizeof (double));
  if (!MshSurfGetNormal (n_pts, bdry_pts, smallest_edge, normal))
  	return 0;


  /* set register functions and call surface mesh generator */
  if (MshSurfType == 1) /* support mesh */
  {
    Surf3DRegFunc  (AuxFuncGetElemSize, AuxFuncIdealPoints, AuxFuncSnapPoint);
  }
  else if (MshSurfType == 2) /* support analytical surface */
  {
    Surf3DRegFunc  (AuxFuncGetElemSize2, AuxFuncIdealPoints2, AuxFuncSnapPoint2);
  }


#if 0

  /* write boundary to matlab .m file */
  if (1)
  {
    FILE *ptr = fopen ("C:\\tecgraf\\prj\\mvgeo\\trunk\\mshsurf\\test\\mshsurf.m", "wt");
    if (ptr != NULL)
    {
      char text[10];
      int i;
      for(i = 0; i < num_edge; i++ )
      {
        int ei = new_edges[i*2+0];
        int ej = new_edges[i*2+1];
        sprintf (text, "%d", i);
        fprintf (ptr, "text(%f, %f, %f, '%s', 'Color', 'r')\n", (bdry_pts[3*ei+0] + bdry_pts[3*ej+0])* 0.5,
          (bdry_pts[3*ei+1] + bdry_pts[3*ej+1])* 0.5,
          (bdry_pts[3*ei+2] + bdry_pts[3*ej+2])* 0.5, text);

        fprintf (ptr, "line([%f, %f], [%f, %f],[%f, %f], 'Marker', 'o')\n",
          bdry_pts[3*ei+0], bdry_pts[3*ei+0] + (bdry_pts[3*ej+0] - bdry_pts[3*ei+0]) *0.7,
          bdry_pts[3*ei+1], bdry_pts[3*ei+1] + (bdry_pts[3*ej+1] - bdry_pts[3*ei+1]) *0.7,
          bdry_pts[3*ei+2], bdry_pts[3*ei+2] + (bdry_pts[3*ej+2] - bdry_pts[3*ei+2]) *0.7 );
      }

      // text
      for(i = 0; i < n_pts; i++ )
      {
        sprintf (text, "%d", i);
        fprintf (ptr, "text(%f, %f, %f, '%s')\n", bdry_pts[3*i+0], bdry_pts[3*i+1],
          bdry_pts[3*i+2], text);
      }

    }
    fclose(ptr);
  }

  // write boundary and support mesh
  if (0)
  {
    FILE *ptr = fopen ("C:\\tecgraf\\prj\\mvgeo\\trunk\\mshsurf\\test\\mshsurf_bdy.dat", "wt");

    if (ptr != NULL)
    {
      char text[10];
      int i;

      fprintf (ptr, "%d\n", n_pts);
      for(i = 0; i < n_pts; i++ )
      {
        fprintf (ptr, "%e, %e, %e\n", bdry_pts[3*i+0], bdry_pts[3*i+1],
          bdry_pts[3*i+2]);
      }

      for(i = 0; i < num_edge; i++ )
      {
        int ei = new_edges[i*2+0];
        int ej = new_edges[i*2+1];
        sprintf (text, "%d", i);
        fprintf (ptr, "text(%f, %f, %f, '%s', 'Color', 'r')\n", (bdry_pts[3*ei+0] + bdry_pts[3*ej+0])* 0.5,
          (bdry_pts[3*ei+1] + bdry_pts[3*ej+1])* 0.5,
          (bdry_pts[3*ei+2] + bdry_pts[3*ej+2])* 0.5, text);

        fprintf (ptr, "line([%f, %f], [%f, %f],[%f, %f], 'Marker', 'o')\n",
          bdry_pts[3*ei+0], bdry_pts[3*ei+0] + (bdry_pts[3*ej+0] - bdry_pts[3*ei+0]) *0.7,
          bdry_pts[3*ei+1], bdry_pts[3*ei+1] + (bdry_pts[3*ej+1] - bdry_pts[3*ei+1]) *0.7,
          bdry_pts[3*ei+2], bdry_pts[3*ei+2] + (bdry_pts[3*ej+2] - bdry_pts[3*ei+2]) *0.7 );
      }

      // text

    }
    fclose(ptr);
  }

#endif


  if (mes_func != NULL)
     mes_func ("Generating new triangles...");
  status = Surf3DAdvFront (n_pts, bdry_pts, normal, num_edge, new_edges, mes_func,
  	                       n_node, coords, n_elem, Conn);


#if 0
  write_neutralfile_teste ("c:\\temp\\out.nf", n_elem, n_node, coords, Conn);
#endif


  free (normal);
  free (new_edges);

  /* release memory */
  if (MshSurfType == 1) /* support mesh */
  {
    AuxFuncRelease ();
  }
  else if (MshSurfType == 2) /* support analytical surface */
  {
    AuxFuncRelease2 ();
  }

  MshSurfType = -1;

  return status;
}


/* Third Function - Release current suppor mesh */
void MshSurfReleaseSupportMesh (void)
{
  if (support_mesh != NULL)
	  SurfTopRelease (support_mesh);
	support_mesh = NULL;
}


/*
** static functions - implementation
** ----------------------------------------------------------------------
*/
int MshSurfGetNormal (int n, double *coords, double size, double *normal)
{
	int i;
	double hint_size;

  if (MshSurfType == -1) /* the support surface is not defined */
    return 0;

	hint_size = size * 0.01;

  if (MshSurfType == 1) /* support mesh */
  {
    for (i = 0; i < n; i++)
	    AuxFuncGetNormal (&(coords[i*3]), hint_size, &(normal[i*3]));
  }
  else if (MshSurfType == 2) /* support analytical surface */
  {
    for (i = 0; i < n; i++)
	    AuxFuncGetNormal2 (&(coords[i*3]), hint_size, &(normal[i*3]));
  }

	return 1;
}


#if 0
/**************************************************************/
/* write_neutralfile_teste
*/
static void write_neutralfile_teste(
char *name_file_teste,
int n_elem,
int n_node,
double *Coords,
int    *conn
)
{

 int index =0, j, i;
 FILE *arquivo;
 int numq4 = 0, numt3 = 0, numq8 = 0, numt6 = 0;

 arquivo = fopen(name_file_teste,"w");

 fprintf(arquivo,"\n%%HEADER\nFile created by program 'Quebra2D' at - AMiranda\n");
 fprintf(arquivo,"\n%%HEADER.VERSION\n1.00\n");
 fprintf(arquivo,"\n%%HEADER.TITLE\n' untitled'\n");

 fprintf(arquivo,"%%NODE\n");
 fprintf(arquivo,"%d\n\n",n_node);
 fprintf(arquivo,"%%NODE.COORD\n");
 fprintf(arquivo,"%d\n\n",n_node);

 for(i=0; i<n_node; i++)
 {
   fprintf(arquivo,"%d    %f    %f    %f\n",i+1, Coords[i*3+0], Coords[i*3+1], Coords[i*3+2]);
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


#include "msh2d.h"


/* Mesh Algorithm
   - Create surface mesh using Msh2DEdge, but the boundary points is used
     to create a plane. If the option to create internal points, these
     internal points is smoothed to 3D */
int MshSurfEdge2D (
int     n_pts,            /* # of points                          (in) */
double  *bdry_pts,        /* coordinates of all points            (in) */
int     bound_edge,       /* # of boundary edges                  (in) */
int     inter_edge,       /* # of internal free edges             (in) */
int     *edges,           /* vector with the edges                (in) */
int     internal_pts,     /* create internal points?               (in) */
int     *n_node,          /* counts # of pts for meshing         (out) */
double  **coords,         /* coordinate array used for meshing   (out) */
int     *n_elem,          /* number of elements generated        (out) */
int     **Conn            /* elem.connectivity list from meshing (out) */
)
{
  double BaseMatrix[4][4];
  double *pts2d, area, *coords2d, h;
  int    i, j, k;

  /* create base matrix that transform 3d points do 2d */
  MshSurfLstSqrPlanMtx (n_pts, bdry_pts, BaseMatrix);

  /* create 2d points */
  pts2d = (double *) calloc (n_pts*2, sizeof (double));
  for (i = 0; i < n_pts; ++i)
    MshSurfCoord_xy (bdry_pts[i*3+0], bdry_pts[i*3+1],bdry_pts[i*3+2],
                     BaseMatrix, &(pts2d[i*2]), &h);

  /* check area */
  /* compute area to check the orientation */
  /* area = x1*y2 - x2*y1  */
  area = 0.0;
  for (k = 0; k < bound_edge; k++)
  {
    i = edges[k*2 + 0];
    j = edges[k*2 + 1];
    area = area + pts2d[i*2+0]*pts2d[j*2+1] -
      pts2d[j*2+0]*pts2d[i*2+1];
  }

  /* invert edges if necessary */
  if (area > 0)
  {
    int tmp_edge;
    for (i = 0; i < bound_edge; i++)
    {
      tmp_edge = edges[i*2 + 0];
      edges[i*2 + 0] = edges[i*2 + 1];
      edges[i*2 + 1] = tmp_edge;
    }
  }

  /* set parameter */
  if (!internal_pts)
  {
    double gen_pts = 0.0;
    Msh2DEdgeParams (1, &gen_pts);
  }

  /* generate points */
  if (!Msh2DEdge (n_pts, pts2d, bound_edge, inter_edge, edges, 3,
                  n_node, &coords2d, n_elem, Conn))
    return 0;


  /* restore default value */
  if (!internal_pts)
  {
    double gen_pts = 1.0;
    Msh2DEdgeParams (1, &gen_pts);
  }


  /* move points to 3d */
  *coords = (double *) calloc ((*n_node)*3, sizeof(double));

  /* boundary points */
  for (i = 0; i < n_pts; ++i)
  {
    (*coords)[i*3+0] = bdry_pts[i*3+0];
    (*coords)[i*3+1] = bdry_pts[i*3+1];
    (*coords)[i*3+2] = bdry_pts[i*3+2];
  }


  for (i = n_pts; i < (*n_node); ++i)
    MshSurfCoord_xyz (coords2d[i*2+0], coords2d[i*2+1], BaseMatrix,
                      &((*coords)[i*3]));

  // smooth
  if (internal_pts)
    MshSurfSmooth (*n_node, *coords, *n_elem, *Conn, 1, 10, 0, NULL);

  free (pts2d);

  return 1;
}

/************************************************************************/
/* MshSurfSetParams                                                     */
/************************************************************************/
void MshSurfSetParams( int index, double value )
{
  switch (index)
  {
  case 0:
    MshSurf_smoothOnProjection = (int) value;
    break;
  }
}

/************************************************************************/
/* MshSurfSetParams                                                     */
/************************************************************************/
double MshSurfGetParams( int index )
{
  switch (index)
  {
  case 0:
    return (double) MshSurf_smoothOnProjection;
    break;
  }

  return -1.0;
}
