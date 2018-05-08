/*
** ---------------------------------------------------------------
**
** msh_gen2d.c: Main driver to generate 2D mesh.
**
** ---------------------------------------------------------------
**
** Created:      01-Oct-97      Joaquim B.C. Neto
**
** ---------------------------------------------------------------
**
*/

#define DISPLAY_CLOCK  1
#define SAVE_FILE      0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "msh_defsurf.h"
#include "msh_quadsurf.h"
#include "msh_bdrsurf.h"


#ifdef _UNIX_
#include <sys/time.h>
#else
#include <sys/timeb.h>
#endif

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1.0e+06
#endif

#define MSHSURF_MAX_PARAM   20

void MshSurfSetDefaultFlagsParams ( void );

/* --------------------------------------------------------------
** Static variables: 
*/
static int    MshSurfFlags  [MSHSURF_MAX_PARAM] = 
{1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
static double MshSurfParams [MSHSURF_MAX_PARAM] =
{-1.0, 3.1416/12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


/* --------------------------------------------------------------
** Static function: 
*/
#if SAVE_FILE
static void writeNF(char *fn, int n_elem, int n_node, double *Coords, int *conn);
#endif


/* --------------------------------------------------------------
** Public function: 
*/
int MshSurfEdge( 
int     n_pts,         /* # of points                          (in) */
double  *bdry_pts,     /* coordinates of all points            (in) */
int     bound_edge,    /* # of boundary free edges             (in) */
int     inter_edge,    /* # of internal free edges             (in) */
int     *edges,        /* edge vector (i0,j0; i1, j1; ...)     (in) */
int     type_mesh,     /* 3 -> T3;  6 -> T6                    (in) */ 
int     *n_node,       /* # of pts in the mesh                (out) */
double  **coords,      /* coordinate array of the mesh        (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn,        /* elem.connectivity list of the mesh  (out) */
void    (*f_surf)      (double, double, double *) 
                       /* function to compute the deviations   (in) */
)
{
 int status ;
 int i, j, k;

 int     num_node = 0, num_edge, *new_edges;
 double  area;

 int     num_gen_nodes ;
 double  *generated_nodes ;
 int     num_elems ;
 int     *elems ;
 double cpu_time;

#if SAVE_FILE

 FILE *ptr = fopen ("contorno.m", "wt");


 /* print points */
 if (n_pts > 0)
 {
   /* create matrix */
   fprintf (ptr,"pts_x = [ ...\n");
   for ( i = 0; i < n_pts; i++)
     fprintf (ptr, " %f ...\n", bdry_pts[i*2+0]);
   fprintf (ptr,"]\n");

   fprintf (ptr,"pts_y = [ ...\n");
   for ( i = 0; i < n_pts; i++)
     fprintf (ptr, " %f ...\n", bdry_pts[i*2+1]);
   fprintf (ptr,"]\n");

   fprintf (ptr,"plot (pts_x, pts_y, 'r+')\n");

   for ( i = 0; i < n_pts; i++)
     fprintf (ptr, "text (%f, %f, '%d')\n", bdry_pts[i*2+0], bdry_pts[i*2+1], i+1);
 }


 for ( i = 0; i < bound_edge; i++)
 { 
   int idi = edges[i*2+0];
   int idj = edges[i*2+1];
   fprintf (ptr, "line ([%f, %f], [%f, %f])\n", bdry_pts[idi*2+0], bdry_pts[idj*2+0],
            bdry_pts[idi*2+1], bdry_pts[idj*2+1]);
 }

 for ( i = bound_edge; i < bound_edge+inter_edge; i++)
 { 
   int idi = edges[i*2+0];
   int idj = edges[i*2+1];
   fprintf (ptr, "line ([%f, %f], [%f, %f])\n", bdry_pts[idi*2+0], bdry_pts[idj*2+0],
     bdry_pts[idi*2+1], bdry_pts[idj*2+1]);
 }

 fclose (ptr);

#endif


 /* compute area to check the orientation */
 /* area = x1*y2 - x2*y1  */
 area = 0.0;
 for (k = 0; k < bound_edge; k++)
 {
   i = edges[k*2 + 0];
   j = edges[k*2 + 1];
   area = area + bdry_pts[i*2+0]*bdry_pts[j*2+1] - 
                 bdry_pts[j*2+0]*bdry_pts[i*2+1];
 }

 /* open memory */
 num_node = n_pts;
 num_edge = ( bound_edge + (inter_edge * 2) );
 new_edges = (int *) calloc (num_edge*2, sizeof(int));

 /* get the edges from the boundary */
 if (area < 0)
 {
   for (i = 0; i < bound_edge; i++)
   {
     /* insert one side of edge */
     new_edges[i*2 + 0] = edges[i*2 + 0];
     new_edges[i*2 + 1] = edges[i*2 + 1];
   }
 }
 else
 {
   for (i = 0; i < bound_edge; i++)
   {
     /* insert one side of edge */
     new_edges[i*2 + 0] = edges[i*2 + 1];
     new_edges[i*2 + 1] = edges[i*2 + 0];
   }
 }
 

 /* get the edges on the domain */
 for (i = bound_edge; i < (inter_edge+bound_edge); i++)
 {
   j = inter_edge + i;

   /* insert one side of edge */
   new_edges[i*2 + 0] = edges[i*2 + 0];
   new_edges[i*2 + 1] = edges[i*2 + 1];

   /* insert another side of edge */
   new_edges[j*2 + 0] = edges[i*2 + 1];
   new_edges[j*2 + 1] = edges[i*2 + 0];
 }

 /* 0.2 Generate internal points if wanted by usual 2D quadtree
    operations */

 cpu_time = clock( );

#if DISPLAY_CLOCK
    printf("\n");
    printf("\t\tCPU Init......................");
#endif

 num_gen_nodes = 0;
 generated_nodes = NULL;

 status = MshSurfGenQuadTree (num_node, num_edge,  /* IN */
                           (double (*)[2]) bdry_pts, (int (*)[2]) new_edges, /* IN */
                           f_surf); /* IN */
 if (status == 0) return 0;

    cpu_time = (clock( ) - cpu_time)/(CLOCKS_PER_SEC*1.0);

#if DISPLAY_CLOCK
    printf("\n");
    printf("\t\tCPU time QuadTree............. %0.3f (s)\n", (double)cpu_time);
#endif

    cpu_time = clock( );

 /* 0.3 Generate Surf mesh by Boundary Contraction operations */
 
 status = SurfBdryContraction (num_node, num_edge,  /* IN */
                               (double (*)[2]) bdry_pts, (int (*)[2]) new_edges, /* IN */
                               &num_gen_nodes, &generated_nodes, &num_elems, &elems,   /* OUT */
                               f_surf ); /* IN */
 if (status == 0) return 0;

 cpu_time = (clock( ) - cpu_time)/(CLOCKS_PER_SEC*1.0);

#if DISPLAY_CLOCK
    printf("\t\tCPU time Contraction.......... %0.3f (s)\n", (double)cpu_time);
    printf("\t\tNumero de elementos..........  %d \n", num_elems);
    printf("\n");
#endif

 *n_node = num_gen_nodes;
 *coords = (double *)generated_nodes;

 *n_elem = num_elems;

 *Conn = (int *) calloc( 4*num_elems, sizeof(int) );
 for( i=0; i<num_elems; i++)
 {
   (*Conn)[i*4+0] = 3;
   (*Conn)[i*4+1] = elems[i*3+0];
   (*Conn)[i*4+2] = elems[i*3+1];
   (*Conn)[i*4+3] = elems[i*3+2];
 }

 /* if the orientation was changed */
 if (area > 0)
 {
   int  pts_swap;
   int  num_field = type_mesh + 1;
   for (i = 0; i < *n_elem; i++)
   {
     for (j = 0; j < type_mesh/2; j++)
     {
       pts_swap = (*Conn)[(i*num_field)+1+j];
       (*Conn)[(i*num_field)+1+j] = (*Conn)[(i*num_field)+type_mesh-j];
       (*Conn)[(i*num_field)+type_mesh-j] = pts_swap;
     }
   }
 }

 free(elems);
 free(new_edges);
 /* release quadtree created with previous function */
 MshSurfQuadFreeAll();

 /* set default parameters */
 MshSurfSetDefaultFlagsParams ( );

#if SAVE_FILE
 writeNF("out.nf", *n_elem,*n_node,*coords,*Conn);
#endif

 return status ;
}


/**************************************************************/
/* MshSurfRegFunc 
*/
void MshSurfRegFunc (MshSurfSizeElement *mshsurf_size)
{
  MshSurfBdryRegFunc (mshsurf_size);
}


/* This function set options to surface mesh generator
**
** nf               - number of flags 
** flags[0] - (0,1) - refine mesh considering the max size of edge.
**                    Default is 1.
** flags[1] - (0,1) - refine mesh considering the max size of element
**                    given in param variable. Default is 0.
** flags[2] - (0,1) - refine locations of high curvature given a max
**                    angle in param variable. Default is 1.
**
** np               - number of parameters
** param[0]         - max size of element. Default is the max size of edge.
** param[1]         - max angle between elements (radius). Default is PI/12.
**
*/
void MshSurfOptions (int nf, int *flags, int np, double *param)
{
  int i;

  /* set flags */
  if (nf > 0 && nf <= 3)
  {
    for (i = 0; i < nf; ++i)
      MshSurfFlags[i] = flags[i];
  }

  /* set parameters */
  if (np > 0 && np <= 2)
  {
    for (i = 0; i < np; ++i)
      MshSurfParams[i] = flags[i];
  }
}


/**************************************************************/
/* functions to get flags and parameters */
void MshSurfSetDefaultFlagsParams ( void )
{
  MshSurfFlags  [0] = 1;
  MshSurfFlags  [1] = 0;
  MshSurfFlags  [2] = 1;
  MshSurfFlags  [3] = 0;
  MshSurfParams [0] = -1.0;
  MshSurfParams [1] = 3.1416 / 12.0;
}


/*******************************************************************/
int MshSurfRefineByMaxEdgeSize_ask  ( void )
{
  return (MshSurfFlags [0]);
}

/*******************************************************************/
int MshSurfRefineByMaxSetSize_ask   ( void )
{
  return (MshSurfFlags [1]);
}

/*******************************************************************/
int MshSurfRefineCurvature_ask ( void )
{
  return (MshSurfFlags [2]);
}

/*******************************************************************/
int MshSurfRefineCurvature2_ask ( void )
{
  return (MshSurfFlags [3]);
}

/*******************************************************************/
double MshSurfGetMaxSizeElement ( void )
{
  return (MshSurfParams [0]);
}

/*******************************************************************/
double MshSurfGetMaxAngle ( void )
{
  return (MshSurfParams [1]);
}


#if SAVE_FILE
/**************************************************************/
/* writeNF
*/
static void writeNF(
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

