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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "msh_def2d.h"

#include "btree.h"   // mshaux library


/* --------------------------------------------------------------
** global variable:
*/
int Msh2DPointGeneration = 1;

/* --------------------------------------------------------------
** Local function:
*/
static void  QuadraticShape (Sh_pt **,int **,int *,int,int,int *,double *);
static void  QuadraticEdge  (Sh_pt **, int **, int *, int, int, double *, int, int *, int *);


/* --------------------------------------------------------------
** Public function:
*/
int Msh2DShape(
int     n_loops,       /* # of circuits (connect.bdry parts)   (in)*/
int     loop_segs[],   /* # of segs (elem.side) on circuits    (in)*/
double  *bdry_pts,     /* coordinates of points on boundary    (in)*/
int     type_mesh,     /********************************************/
int     *n_node,       /* counts # of pts for meshing         (out)*/
double  **coords,      /* coordinate array used for meshing   (out)*/
int     *n_elem,       /* number of elements generated        (out)*/
int     **Conn         /* elem.connectivity list from meshing (out)*/
)
{
    int status ;
    int i,j, first_node = 0;
    int     num_node = 0;
    int     num_edge ;
    double  *nodes ;
    int     *edges ;
    int     num_gen_nodes = 0 ;
    double  *generated_nodes ;
    int     num_elems ;
    int     *elems ;
    int     pass;

    if(type_mesh==3)
     pass=1;
    else if(type_mesh==6)
     pass=2;
    else return 0;

    for ( i=0; i<n_loops; i++)
     num_node = num_node + loop_segs[i];

    num_edge = num_node/pass;
    nodes = (double *) calloc( num_node*2/pass, sizeof(double) );
    edges = (int *) calloc( num_edge*2, sizeof(int) );

    /* get the coodinates of bonundary to nodes vector */
    num_node = 0;
    for( i=0; i<n_loops; i++)
    {
     for( j=0; j<loop_segs[i]/pass; j++ )
     {
      nodes[num_node*2] = bdry_pts[num_node*2*pass];
      nodes[(num_node*2)+1] = bdry_pts[num_node*2*pass+1];
      if( j!=0 )
      {
       edges[(num_node-1)*2+0] = num_node-1;
       edges[(num_node-1)*2+1] = num_node;
      }
      else
       first_node = num_node;
      num_node++;
     }
     edges[(num_node-1)*2+0] = num_node-1;
     edges[(num_node-1)*2+1] = first_node;
    }

    /* 0.1 Initialize global parameters to generate the mesh */


    /* 0.2 Generate 2D quadtree */
    if (Msh2DPointGeneration)
    {
      status = Msh2DGenQuadTree (num_node, num_edge,
        (double (*)[2]) nodes, (int (*)[2]) edges,
        &num_gen_nodes, &generated_nodes ) ;
    }



    /* 0.3 Generate 2D mesh by Boundary Contraction operations */

    status = Msh2DBdryContraction (num_node, num_edge,
             (double (*)[2]) nodes, (int (*)[2]) edges,
             &num_gen_nodes, &generated_nodes, &num_elems, &elems ) ;
    if( status == 0 )
    {
      *n_node = 0;
      *coords = NULL;
      *n_elem = 0;
      *Conn = NULL;
      return status;
    }

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

/*  Case type_mesh = 6 (T6) => generate quadratic elements    */
    if ( type_mesh == 6)
    {
       Sh_pt *Coords=(Sh_pt *)(*coords);
       QuadraticShape ( &Coords, Conn, n_node, *n_elem, n_loops, loop_segs, bdry_pts);
       *coords = (double *)Coords;
    }

    if (num_elems > 0) free(elems);
    free(edges);
    free(nodes);

    Msh2DPointGeneration = 1;

    return status ;
}



/* --------------------------------------------------------------
** Public function:
*/
int Msh2DEdge(
int     n_pts,         /* # of points                          (in) */
double  *bdry_pts,     /* coordinates of all points            (in) */
int     bound_edge,    /* # of boundary free edges             (in) */
int     inter_edge,    /* # of internal free edges             (in) */
int     *edges,        /* vector with thed edges               (in) */
int     type_mesh,     /*********************************************/
int     *n_node,       /* counts # of pts for meshing         (out) */
double  **coords,      /* coordinate array used for meshing   (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn         /* elem.connectivity list from meshing (out) */
)
{
 int status ;
 int i, j, id_i, id_j;

 int     num_node = 0, num_edge, *new_edges;
 double  *nodes ;

 int     num_gen_nodes ;
 double  *generated_nodes ;
 int     num_elems ;
 int     *elems ;
 int     pass, *id_nodes;

 if(type_mesh==3)
   pass=1;
 else if(type_mesh==6)
   pass=2;
 else return 0;

 num_node = n_pts; /* /pas; */
 num_edge = ( bound_edge + (inter_edge * 2) ) / pass;

 /* open memory */
 nodes = (double *) calloc (n_pts*2, sizeof(double));
 new_edges = (int *) calloc (num_edge*2, sizeof(int));

 id_nodes = (int *) calloc (n_pts, sizeof(int));
 memset (id_nodes, -1, sizeof (int) * n_pts);

 /* when element is linear, the procedure is simple. */
 if(type_mesh == 3)
 {
   memcpy (nodes, bdry_pts, n_pts*2*sizeof (double));
   memcpy (new_edges, edges, (bound_edge+inter_edge)*2*sizeof (int));

   /* duplicate edges */
   for (i = bound_edge; i < (bound_edge+inter_edge); i++)
   {
     new_edges[(i+inter_edge)*2 + 0] = edges[i*2+0];
     new_edges[(i+inter_edge)*2 + 1] = edges[i*2+1];
   }
 }

 /* get the edges from the boundary */
 if (type_mesh == 6)
 {
   num_node = 0;
   for (i = 0; i < bound_edge/pass; i++)
   {
     id_i = edges[i*2*pass + 0];
     id_j = edges[i*2*pass + (pass-1)*2 + 1];

     /* check if the node id_i already has included in nodal vector */
     if (id_i < id_j)
     {
       if (id_nodes[id_i] == -1)
       {
         nodes[num_node*2]     = bdry_pts[id_i*2];
         nodes[(num_node*2)+1] = bdry_pts[id_i*2+1];
         id_nodes[id_i] = num_node;
         id_i = num_node;
         num_node++;
       }
       else
         id_i = id_nodes[id_i];

       /* check if the node id_j already has included in nodal vector */
       if (id_nodes[id_j] == -1)
       {
         nodes[num_node*2]     = bdry_pts[id_j*2];
         nodes[(num_node*2)+1] = bdry_pts[id_j*2+1];
         id_nodes[id_j] = num_node;
         id_j = num_node;
         num_node++;
       }
       else
         id_j = id_nodes[id_j];
     }
     else
     {
       /* check if the node id_j already has included in nodal vector */
       if (id_nodes[id_j] == -1)
       {
         nodes[num_node*2]     = bdry_pts[id_j*2];
         nodes[(num_node*2)+1] = bdry_pts[id_j*2+1];
         id_nodes[id_j] = num_node;
         id_j = num_node;
         num_node++;
       }
       else
         id_j = id_nodes[id_j];

       if (id_nodes[id_i] == -1)
       {
         nodes[num_node*2]     = bdry_pts[id_i*2];
         nodes[(num_node*2)+1] = bdry_pts[id_i*2+1];
         id_nodes[id_i] = num_node;
         id_i = num_node;
         num_node++;
       }
       else
         id_i = id_nodes[id_i];
     }

     /* insert one side of edge */
     new_edges[i*2 + 0] = id_i;
     new_edges[i*2 + 1] = id_j;
   }

   /* get the edges on the domain */
   for (i = bound_edge/pass; i < (inter_edge+bound_edge)/pass; i++)
   {
     id_i = edges[i*2*pass + 0];
     id_j = edges[i*2*pass + (pass-1)*2 + 1];
     j = inter_edge/pass + i;

     /* check if the node id_i already has included in nodes vector */
     if (id_nodes[id_i] == -1)
     {
       nodes[num_node*2]     = bdry_pts[id_i*2];
       nodes[(num_node*2)+1] = bdry_pts[id_i*2+1];
       id_nodes[id_i] = num_node;
       id_i = num_node;
       num_node++;
     }
     else
       id_i = id_nodes[id_i];

     /* check if the node id_j already has included in nodes vector */
     if (id_nodes[id_j] == -1)
     {
       nodes[num_node*2]     = bdry_pts[id_j*2];
       nodes[(num_node*2)+1] = bdry_pts[id_j*2+1];
       id_nodes[id_j] = num_node;
       id_j = num_node;
       num_node++;
     }
     else
       id_j = id_nodes[id_j];

     /* insert one side of edge */
     new_edges[i*2 + 0] = id_i;
     new_edges[i*2 + 1] = id_j;

     /* insert another side of edge */
     new_edges[j*2 + 0] = id_j;
     new_edges[j*2 + 1] = id_i;
   }

   free (id_nodes);
 }


  /* 0.1 Initialize global parameters to generate the mesh */



 /* 0.2 Generate internal points if wanted by usual 2D quadtree
    operations */

 status = Msh2DGenQuadTree (num_node, num_edge,  /* IN */
                            (double (*)[2]) nodes, (int (*)[2]) new_edges, /* IN */
                            &num_gen_nodes, &generated_nodes ) ; /* OUT */

 /* 0.3 Generate 2D mesh by Boundary Contraction operations */

 status = Msh2DBdryContraction (num_node, num_edge,  /* IN */
                                (double (*)[2]) nodes, (int (*)[2]) new_edges, /* IN */
                                &num_gen_nodes, &generated_nodes, &num_elems, &elems ); /* OUT */

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

 /*  Case type_mesh = 6 (T6) => generate quadratic elements    */
 if ( type_mesh == 6)
 {
    Sh_pt *Coords=(Sh_pt *)(*coords);
    QuadraticEdge (&Coords, Conn, n_node, *n_elem,
                   n_pts, bdry_pts, bound_edge + inter_edge, edges, new_edges);
    *coords = (double *)Coords;
 }

 if (num_elems > 0)
   free(elems);
 free(new_edges);
 free(nodes);

 return status ;
}


/*  Msh2DRegFunc */
void Msh2DRegFunc (Msh2DSizeElement *msh2d_size)
{
  Msh2DSizeFunction = msh2d_size;
 }

/*
** Quadratic mesh generation.
**
**  This function updates the mesh structure from triangulation and quadtree
**  generation, to assure quadratic mesh. It is included quadratic mid nodes
**  in each element side, and in a binary tree to an eficient search.
**  it returns 1 if all quadratic generation it was done.
*/

static Tbtree *edgetree;

typedef struct
{
  long int i;
  long int j;

  long int mid;

  int elem;
  int ciclo;
  int flag_change;
} Tcvtedge;


static int fedgecomp(void *e1,void *e2)
{
 Tcvtedge *ef1=e1, *ef2=e2;
 if      (ef1->i < ef2->i)      return -1;
 else if (ef1->i > ef2->i)      return  1;
 else if (ef1->j < ef2->j)      return -1;
 else if (ef1->j > ef2->j)      return  1;
 else return  0;

}



/*=========================  QuadraticEdge ========================*/
static void QuadraticEdge (
Sh_pt **coords,     /* Coord. of the domain mesh        */
int **Conn,         /* Connectivity                     */
int *n_node,        /* number of nodes                  */
int n_elem,         /* number of elements               */

int n_pts,          /* nunber of points  */
double *Bdry_pts,    /* Coord               */
int n_edge,
int *edges,
int *new_edges
)
{
 int       i, j, index, id1, id2;
 int       indx_edge;
 int       *conn;
 Tcvtedge  *Edge, *TmpEdge, FindEdge;

 /****** init quadratic elements ********/
 conn = (int *) calloc(n_elem*7,sizeof(int));
 index=0;

 /****** realloc the middle coordenates *******/
 (*coords) = realloc ((*coords), (6*n_elem) * sizeof (Sh_pt));

 /* init Btree */
 edgetree = NULL;
 BtreeInit (fedgecomp);

 Edge = (Tcvtedge *) calloc (3*n_elem, sizeof(Tcvtedge));

 /* build edge list */
 index=0;
 indx_edge = 0;
 for (i = 0; i < n_elem; i++)  /* for all elements */
 {
   int type = (*Conn)[index];
   conn[i*7] = type * 2;
   for (j = 0; j < type; j++)
   {
      /* fill the estructure */
     id1 = (*Conn)[index+1+j];
     id2 = (*Conn)[index+1+((j+1)%type)];
     Edge[indx_edge].elem = i;
     Edge[indx_edge].ciclo = j;
     if (id1 < id2)
     {
       Edge[indx_edge].i = id1;  Edge[indx_edge].j = id2;
       Edge[indx_edge].flag_change = 0;
     }
     else
     {
       Edge[indx_edge].i = id2;  Edge[indx_edge].j = id1;
       Edge[indx_edge].flag_change = 1;
     }

     /* look for a edge */
     TmpEdge = NULL;
     TmpEdge = BtreeFind (edgetree, &Edge[indx_edge]);

     /* case the edge does't exist -> insert new node and update id */
     if (TmpEdge == NULL)
     {
       /* add one node */
       (*coords)[(*n_node) + indx_edge].x = ((*coords)[id1].x + (*coords)[id2].x)/2.0;
       (*coords)[(*n_node) + indx_edge].y = ((*coords)[id1].y + (*coords)[id2].y)/2.0;
       /* add Edge structure in Btree */
       Edge[indx_edge].mid = (*n_node) + indx_edge;
       Edge[indx_edge].elem = i;
       Edge[indx_edge].ciclo = j;
       edgetree = BtreeInsert (edgetree, &Edge[indx_edge]);

       TmpEdge = &Edge[indx_edge];
       indx_edge++;
     }

     /* update conn */
     conn[(i*7)+1+(j*2)]   = (*Conn)[index+1+j];
     conn[(i*7)+1+(j*2)+1] = TmpEdge->mid;
   }
   index=index+(*Conn)[index]+1;
 }

 free(*Conn);
 *Conn = conn;
 *n_node = *n_node + indx_edge;

#if 1

 /* update nodes with boundary and internal edges */
 for (i = 0; i < n_edge/2; i++)
 {
   int id_mid = edges[i*4 + 1];
   id1 = new_edges[i*2 + 0];
   id2 = new_edges[i*2 + 1];

   if (id1 < id2)
   {
     FindEdge.i = id1;  FindEdge.j = id2;
     FindEdge.flag_change = 0;
   }
   else
   {
     FindEdge.i = id2;  FindEdge.j = id1;
     FindEdge.flag_change = 1;
   }


   /* look for a edge */
   TmpEdge = NULL;
   TmpEdge = BtreeFind (edgetree, &FindEdge);
   if (TmpEdge != NULL)
   {
/*     printf ("mid = %d\n", TmpEdge->mid);
 */
#if 1

     (*coords)[TmpEdge->mid].x = Bdry_pts [id_mid*2 + 0];
     (*coords)[TmpEdge->mid].y = Bdry_pts [id_mid*2 + 1];

#else

{
 Tpoint pts_edge[2];
 Tpoint pts_1[2];

 pts_edge[0].x = (*coords)[id1].x;
 pts_edge[0].y = (*coords)[id1].y;

 pts_edge[1].x = (*coords)[id2].x;
 pts_edge[1].y = (*coords)[id2].y;

 pts_1[0].x = (*coords)[TmpEdge->mid].x;
 pts_1[0].y = (*coords)[TmpEdge->mid].y;


 GraSetLineColor (GRA_RED);

 GraPolyLine (2, pts_edge);

 GraSetMarkColor (GRA_BLACK);

 GraPolyMark (1, pts_1);



     (*coords)[TmpEdge->mid].x = Bdry_pts [id_mid*2 + 0];
     (*coords)[TmpEdge->mid].y = Bdry_pts [id_mid*2 + 1];


 pts_1[0].x = (*coords)[TmpEdge->mid].x;
 pts_1[0].y = (*coords)[TmpEdge->mid].y;

 GraSetMarkColor (GRA_BLUE);

 GraPolyMark (1, pts_1);

 IupMessage ("Points", "Contorno");

}

#endif


   }


 }

/* UpdateQuadraticEdges (*coords, Bdry_pts, n_edge, edges);
*/
#endif


/* BtreeRelease(edgetree); */
 free (Edge);
}

/*=========================  QuadraticShape ============================*/
static void QuadraticShape(
Sh_pt **coords,     /* Coord. of   domain mesh         */
int **Conn,         /* Connectivity                    */
int *n_node,        /* number of nodes                 */
int n_elem,         /* number of elements              */
int n_loops,        /* nunber of loops of boundary     */
int *loop_segs,     /* number of segments in each loop */
double *Bdry_pts    /* Coord. of boundary              */
)
{
 int       i, j, index, id1, id2;
 int       indx_edge, num_node;
 int       *conn;
 Tcvtedge  *Edge, *TmpEdge, FindEdge;

 /****** init quadratic elements ********/
 conn = (int *) calloc(n_elem*7,sizeof(int));
 index=0;

 /****** realloc the middle coordenates *******/
 (*coords) = realloc ((*coords), (6*n_elem) * sizeof (Sh_pt));

 /* init Btree */
 edgetree = NULL;
 BtreeInit (fedgecomp);

 Edge = (Tcvtedge *) calloc (3*n_elem, sizeof(Tcvtedge));

 /* build edge list */
 index=0;
 indx_edge = 0;
 for (i = 0; i < n_elem; i++)  /* for all elements */
 {
   int type = (*Conn)[index];
   conn[i*7] = type * 2;
   for (j = 0; j < type; j++)
   {
      /* fill the estructure */
     id1 = (*Conn)[index+1+j];
     id2 = (*Conn)[index+1+((j+1)%type)];
     Edge[indx_edge].elem = i;
     Edge[indx_edge].ciclo = j;
     if (id1 < id2)
     {
       Edge[indx_edge].i = id1;  Edge[indx_edge].j = id2;
       Edge[indx_edge].flag_change = 0;
     }
     else
     {
       Edge[indx_edge].i = id2;  Edge[indx_edge].j = id1;
       Edge[indx_edge].flag_change = 1;
     }

     /* look for a edge */
     TmpEdge = NULL;
     TmpEdge = BtreeFind (edgetree, &Edge[indx_edge]);

     /* case the edge does't exist -> insert new node and update id */
     if (TmpEdge == NULL)
     {
       /* add one node */
       (*coords)[(*n_node) + indx_edge].x = ((*coords)[id1].x + (*coords)[id2].x)/2.0;
       (*coords)[(*n_node) + indx_edge].y = ((*coords)[id1].y + (*coords)[id2].y)/2.0;
       /* add Edge structure in Btree */
       Edge[indx_edge].mid = (*n_node) + indx_edge;
       Edge[indx_edge].elem = i;
       Edge[indx_edge].ciclo = j;
       edgetree = BtreeInsert (edgetree, &Edge[indx_edge]);

       TmpEdge = &Edge[indx_edge];
       indx_edge++;
     }

     /* update conn */
     conn[(i*7)+1+(j*2)]   = (*Conn)[index+1+j];
     conn[(i*7)+1+(j*2)+1] = TmpEdge->mid;
   }
   index=index+(*Conn)[index]+1;
 }

 free(*Conn);
 *Conn = conn;
 *n_node = *n_node + indx_edge;

 /* update nodes with real boundary */
 num_node = 0;
 for (i = 0; i < n_loops; i++)
 {
   int pass = 2;
   int first_node = num_node;
   int id_mid, id1, id2;

   for (j = 0; j < loop_segs[i] / pass; j++ )
   {
     id_mid = num_node*pass+1;
     if (j == loop_segs[i] / pass - 1) /* last point */
     {
       id1 = num_node;
       id2 = first_node;
     }
     else
     {
       id1 = num_node;
       id2 = num_node+1;
     }
     num_node++;


     /* look for a edge */
     if (id1 < id2)
     {
       FindEdge.i = id1;  FindEdge.j = id2;
       FindEdge.flag_change = 0;
     }
     else
     {
       FindEdge.i = id2;  FindEdge.j = id1;
       FindEdge.flag_change = 1;
     }
     TmpEdge = NULL;
     TmpEdge = BtreeFind (edgetree, &FindEdge);

     /* update mid point */
     if (TmpEdge != NULL)
     {
       (*coords)[TmpEdge->mid].x = Bdry_pts[id_mid*2 + 0];
       (*coords)[TmpEdge->mid].y = Bdry_pts[id_mid*2 + 1];
     }
   }

 }


 /* BtreeRelease(edgetree); */
 free (Edge);

}

/*
************************************************************/
void Msh2DEdgeParams (int _n, double *_param)
{
  if (_n == 0 || _param == NULL)
    return;
  Msh2DPointGeneration = (int) _param[0];
}

