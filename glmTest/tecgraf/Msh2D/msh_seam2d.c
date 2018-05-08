/*
** ---------------------------------------------------------------------------
**
**  msh_seam2d.c - This file contains a bind to cpp for generate quadrilateral 
**                 elements from triangles (Seam technique).
** 
** ---------------------------------------------------------------------------
**
**  Public Mesh Generating Function :
**
** ---------------------------------------------------------------------------
**
** void Msh2DSeam ( n_loops, loop_segs[], bdry_pts, flag_mesh, n_node, 
**                  coords, n_elem, conn )
**
**   int     n_loops      - # of circuits (connect.bdry parts)         (in )
**   int     loop_segs[]  - # of segs (elem.side) on circuits          (in )
**   GeoPnt  *bdry_pts    - coordinates of points on boundary          (in )
**   int     flag_mesh    - flag for triangular or quadrilateral mesh  (in )
**   int     *n_node      - counts # of pts for meshing                (out)
**   GeoPnt  *coords      - coordinate array used for meshing          (out)
**   int     *n_elem      - number of elements generated               (out)
**   int     *conn        - elem.connectivity list from meshing        (out)
** 
** ---------------------------------------------------------------------------
**
** Private Tree Functions :
**
** ---------------------------------------------------------------------------
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mshq2d.h"
#include "msh2d.h"

#include "btree.h"			// mshaux library




/* --------------------------------------------------------------
** Local function:
*/
static void  QuadraticShape (double **,int **,int *,int,int,int *,double *);


/** 
    @param n_pts      (in) número de pontos de entrada 
    @param bdry_pts   (in) vetor de coordenadas
    @bound_edge       (in) número de aresta do contorno
    @inter_edge       (in) número de aresta internas (restricões)
    @edges            (in) arestas (id i, id j)   
    @param type_mesh  (in) Tipo de elemento: T3 (3) ou T6 (6)
    @param n_node     (out) número de nós gerados
    @param coords     (out) coordenadas dos nós da malha
    @param n_elem     (out) número de elementos gerados
    @param conn       (out) conectividade dos elementos
*/
int Msh2DQuadSeamEdge (int n_pts, double *bdry_pts, int bound_edge, int *edges, 
                       int *n_node, double  **coords, int *n_elem, int **Conn)
{
  int    ok, i, j, index_l, index_q;
  int    tmp_n_node, tmp_n_elem, *tmp_conn;
  double *tmp_coords, mid_coord[2];
  int     n_loops, *loop_segs; 
  
  if (bound_edge%2 != 0)
     return 0;
  
  n_loops = 1;
  loop_segs = calloc (1, sizeof (int));
  loop_segs[0] = bound_edge;

  if (edges[0] > edges[1])
  {
    double *tmp_bdry;
    tmp_bdry = (double *) calloc (n_pts * 2, sizeof (double));

    tmp_bdry[0] = bdry_pts[0];
    tmp_bdry[1] = bdry_pts[1];
    for (i = 1; i < n_pts; i++)
    {
      tmp_bdry[i*2+0] = bdry_pts[(n_pts-i)*2+0];
      tmp_bdry[i*2+1] = bdry_pts[(n_pts-i)*2+1];
    }

    /* calls mesh generation */
    ok = Msh2DQuadSeam (n_loops, loop_segs, tmp_bdry, 8,
                        &tmp_n_node, &tmp_coords, &tmp_n_elem, &tmp_conn);
    free (tmp_bdry);
  }
  else
  {
    /* calls mesh generation */
    ok = Msh2DQuadSeam (n_loops, loop_segs, bdry_pts, 8,
                        &tmp_n_node, &tmp_coords, &tmp_n_elem, &tmp_conn);
  }

  
  if (!ok)
    return 0;
    
  /* build all quads */
  *n_elem = index_q = 0;
  *n_node = tmp_n_node + tmp_n_elem;
  (*coords) = (double *) calloc ((tmp_n_node + tmp_n_elem)*2, sizeof (double));
  (*Conn)  = (int *) calloc (tmp_n_elem * 4 * 5, sizeof (int));

  /* copy node vector */
  for (i = 0; i < tmp_n_node*2; i++)
    (*coords)[i] = tmp_coords[i];
  
  for (i = 0; i < tmp_n_elem; i++)
  {
    index_l = i * 9 + 1;
    mid_coord[0] = 0.0;
    mid_coord[1] = 0.0;
    if (tmp_conn[index_l+0] == tmp_conn[index_l+7]) /* triangular elements */
    {
      for (j = 0; j < 6; j+=2)  /* obtain only corners */
      {
        mid_coord[0] += tmp_coords[(tmp_conn[index_l+j]*2)+0]; 
        mid_coord[1] += tmp_coords[(tmp_conn[index_l+j]*2)+1];
      }
      mid_coord[0] /= 3.0; 
      mid_coord[1] /= 3.0;
      (*coords)[(tmp_n_node+i)*2 + 0] = mid_coord[0];
      (*coords)[(tmp_n_node+i)*2 + 1] = mid_coord[1];
      
      (*Conn)[index_q+0] = 4;
      (*Conn)[index_q+1] = tmp_conn[index_l+0];
      (*Conn)[index_q+2] = tmp_conn[index_l+1];
      (*Conn)[index_q+3] = tmp_n_node+i;
      (*Conn)[index_q+4] = tmp_conn[index_l+5];
      index_q += 5;
      (*Conn)[index_q+0] = 4;
      (*Conn)[index_q+1] = tmp_conn[index_l+2];
      (*Conn)[index_q+2] = tmp_conn[index_l+3];
      (*Conn)[index_q+3] = tmp_n_node+i;
      (*Conn)[index_q+4] = tmp_conn[index_l+1];
      index_q += 5;
      (*Conn)[index_q+0] = 4;
      (*Conn)[index_q+1] = tmp_conn[index_l+4];
      (*Conn)[index_q+2] = tmp_conn[index_l+5];
      (*Conn)[index_q+3] = tmp_n_node+i;
      (*Conn)[index_q+4] = tmp_conn[index_l+3];
      index_q += 5;
      
      *n_elem += 3;
      
    }
    else  /* quadrilateral elements */
    {
      for (j = 0; j < 8; j+=2)  /* obtain only corners */
      {
        mid_coord[0] += tmp_coords[(tmp_conn[index_l+j]*2)+0]; 
        mid_coord[1] += tmp_coords[(tmp_conn[index_l+j]*2)+1];
      }
      mid_coord[0] /= 4.0; 
      mid_coord[1] /= 4.0;
      (*coords)[(tmp_n_node+i)*2 + 0] = mid_coord[0];
      (*coords)[(tmp_n_node+i)*2 + 1] = mid_coord[1];
      
      (*Conn)[index_q+0] = 4;
      (*Conn)[index_q+1] = tmp_conn[index_l+0];
      (*Conn)[index_q+2] = tmp_conn[index_l+1];
      (*Conn)[index_q+3] = tmp_n_node+i;
      (*Conn)[index_q+4] = tmp_conn[index_l+7];
      index_q += 5;
      (*Conn)[index_q+0] = 4;
      (*Conn)[index_q+1] = tmp_conn[index_l+2];
      (*Conn)[index_q+2] = tmp_conn[index_l+3];
      (*Conn)[index_q+3] = tmp_n_node+i;
      (*Conn)[index_q+4] = tmp_conn[index_l+1];
      index_q += 5;
      (*Conn)[index_q+0] = 4;
      (*Conn)[index_q+1] = tmp_conn[index_l+4];
      (*Conn)[index_q+2] = tmp_conn[index_l+5];
      (*Conn)[index_q+3] = tmp_n_node+i;
      (*Conn)[index_q+4] = tmp_conn[index_l+3];
      index_q += 5;
      (*Conn)[index_q+0] = 4;
      (*Conn)[index_q+1] = tmp_conn[index_l+6];
      (*Conn)[index_q+2] = tmp_conn[index_l+7];
      (*Conn)[index_q+3] = tmp_n_node+i;
      (*Conn)[index_q+4] = tmp_conn[index_l+5];
      index_q += 5;
      
      *n_elem += 4;

    }
     
  }
  
  free (tmp_conn);
  free (tmp_coords);
  

  return ok;
}
               


/** 
    @param n_loops    (in) número de circuitos
    @param loop_segs  (in) vetor com o número de segmentos em cada  circuito    
    @param bdry_pts   (in) vetor de coordenadas do contorno
    @param type_mesh  (in) Tipo de elemento: 4 (Q4) ou 8 (Q8)
    @param n_node     (out) número de nós gerados
    @param coords     (out) coordenadas dos nós da malha
    @param n_elem     (out) número de elementos gerados
    @param conn       (out) conectividade dos elementos
*/
int Msh2DQuadSeam (int n_loops, int *loop_segs, double *bdry_pts, int type_mesh,
                   int *n_node, double **coords, int *n_elem, int **conn)
{
  int    i, j, ok, pass, first_node = 0;
  int    *edges, inter_edge, num_edge, num_node;
  double *nodes;

  /* get the number of edges */
  if (type_mesh == 4)
    pass = 1;
  else if (type_mesh == 8)
    pass = 2;
  else return 0;

  /* get the number of nodes */
  num_node = 0;
  for (i = 0; i < n_loops; i++)
    num_node += loop_segs[i];


  num_edge   = num_node / pass;
  inter_edge = 0;
  nodes = (double *) calloc (num_node * 2 / pass, sizeof (double));
  edges = (int *)    calloc (num_edge * 2 , sizeof (int));

  /* get the coodinates of boundary to nodes vector */
  num_node = 0;
  for (i = 0; i < n_loops; i++)
  {
    for (j = 0; j < loop_segs[i] / pass; j++ )
    {
      nodes[num_node*2]     = bdry_pts[num_node*2*pass];
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

  /* Msh2DSetRefFactor (3.0); */

  /* calls mesh generation */
  ok = Msh2DSeamGeneration (num_node, nodes, num_edge, inter_edge, edges, 4,
                            n_node, coords, n_elem, conn);

  if (!ok)
    return 0;

  /*  Case type_mesh = 8 (Q8) => generate quadratic elements    */
  if ( type_mesh == 8)
     QuadraticShape (coords, conn, n_node, *n_elem, n_loops, loop_segs, bdry_pts);

  free (nodes);
  free (edges);

  return ok;
}



/*
** Quadratic mesh generation to Msh2dShape.
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


/*=========================  QuadraticShape ============================*/
static void QuadraticShape(
double **coords,    /* Coord. of   domain mesh         */
int **Conn,         /* Connectivity                    */
int *n_node,        /* number of nodes                 */
int n_elem,         /* number of elements              */
int n_loops,        /* nunber of loops of boundary     */
int *loop_segs,     /* number of segments in each loop */
double *Bdry_pts    /* Coord. of boundary              */
)
#if 1
{
 int       i, j, index, id1, id2;
 int       indx_edge, num_node;
 int       *conn;
 Tcvtedge  *Edge, *TmpEdge, FindEdge;

 /****** init quadratic elements ********/                    
 conn = (int *) calloc(n_elem*9,sizeof(int));                  
 index=0;

 /****** realloc the middle coordenates *******/               
 (*coords) = realloc ((*coords), (8*n_elem) * 2 * sizeof (double));

 /* init Btree */
 edgetree = NULL;
 BtreeInit (fedgecomp);

 Edge = (Tcvtedge *) calloc (4*n_elem, sizeof(Tcvtedge));

 /* build edge list */
 index=0;
 indx_edge = 0;
 for (i = 0; i < n_elem; i++)  /* for all elements */
 {
   int type = (*Conn)[index];
   conn[i*9] = type * 2;
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
     else if (id1 > id2)
     {
       Edge[indx_edge].i = id2;  Edge[indx_edge].j = id1;
       Edge[indx_edge].flag_change = 1;
     }
     else
     {
       conn[(i*9)+1+(j*2)]   = (*Conn)[index+1+j];
       conn[(i*9)+1+(j*2)+1] = (*Conn)[index+1+j];
       continue;
     }

     /* look for a edge */
     TmpEdge = NULL;
     TmpEdge = BtreeFind (edgetree, &Edge[indx_edge]);
  
     /* case the edge does't exist -> insert new node and update id */
     if (TmpEdge == NULL)
     { 
       /* add one node */
       (*coords)[((*n_node) + indx_edge) * 2 + 0] = ((*coords)[id1*2+0] + (*coords)[id2*2+0])/2.0;
       (*coords)[((*n_node) + indx_edge) * 2 + 1] = ((*coords)[id1*2+1] + (*coords)[id2*2+1])/2.0;
       /* add Edge structure in Btree */  
       Edge[indx_edge].mid = (*n_node) + indx_edge;
       Edge[indx_edge].elem = i;
       Edge[indx_edge].ciclo = j;
       edgetree = BtreeInsert (edgetree, &Edge[indx_edge]);

       TmpEdge = &Edge[indx_edge];
       indx_edge++;
     }

     /* update conn */
     conn[(i*9)+1+(j*2)]   = (*Conn)[index+1+j];
     conn[(i*9)+1+(j*2)+1] = TmpEdge->mid;
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
     else if (id1 > id2)
     {
       FindEdge.i = id2;  FindEdge.j = id1;
       FindEdge.flag_change = 1;
     }
     else
       continue;

     TmpEdge = NULL;
     TmpEdge = BtreeFind (edgetree, &FindEdge);

     /* update mid point */
     if (TmpEdge != NULL)
     {
       (*coords)[TmpEdge->mid * 2 + 0] = Bdry_pts[id_mid*2 + 0];
       (*coords)[TmpEdge->mid * 2 + 1] = Bdry_pts[id_mid*2 + 1];
     }
   }

 }


 /* BtreeRelease(edgetree); */
 free (Edge);

}


#else
{
  int        i, j, index;
  int        nedge;
  int        nn_mid, *conn;
  Seam_edge *Edge;
  int        num_node_bound = 0, mid_id;

  for ( i=0; i<n_loops; i++)
    num_node_bound = num_node_bound + loop_segs[i];

  nedge = 4*n_elem;
  Edge = (Seam_edge *) calloc (nedge, sizeof(Seam_edge));

  /* build edge list */
  index = 0;
  for (i = 0; i < n_elem; i++)
  {
    int type = (*Conn)[index];
    for ( j = 0; j < (*Conn)[index]; j++)
    {
      int id1 = (*Conn)[index+1+j];
      int id2 = (*Conn)[index+1+((j+1)%type)];
      Edge[i*4+j].elem  = i;
      Edge[i*4+j].ciclo = j;
      if (id1 < id2)
      {
        Edge[i*4+j].un.node.i = id1;
        Edge[i*4+j].un.node.j = id2;
        Edge[i*4+j].flag_change = 0;
      }
      else
      {
        Edge[i*4+j].un.node.i = id2;
        Edge[i*4+j].un.node.j = id1;
        Edge[i*4+j].flag_change = 1;
      }
    }
    index = index + (*Conn)[index] + 1;
  }

  /* sort edge list( using C quicksort routine ) */
  qsort (Edge, nedge, sizeof(Seam_edge), fcmp);

#if 0
  {
    FILE *fl = fopen ("d:\\temp\\qsort.dat", "wt");
    for (i = 0; i < nedge; i++)
      fprintf (fl, "%4d - %4d -> %d\n", Edge[i].un.node.i, Edge[i].un.node.j, Edge[i].un.code);
    fclose (fl);
  }
#endif

  /******  fill quadratic elements ********/
  conn = (int *) calloc(n_elem*9,sizeof(int));
  index=0;

  /* get the old Conn  */
  for (i = 0; i < n_elem; i++)
  {
    conn[i*9] = (*Conn)[index] * 2;
    for (j = 0; j < (*Conn)[index]; j++)
      conn[(i*9)+1+(j*2)] = (*Conn)[index+1+j];

    index=index+(*Conn)[index]+1;
  }
 
  nn_mid=0;

  /* count the number of midnodes and fill conn with new nodes */
  for (i = 0; i < nedge-1; i++)
  {
    if (Edge[i].un.code == Edge[i+1].un.code)
    {
      conn[(Edge[i+0].elem*9)+2+(Edge[i+0].ciclo*2)] = *n_node + nn_mid;  /* elem edge i+0 */
      conn[(Edge[i+1].elem*9)+2+(Edge[i+1].ciclo*2)] = *n_node + nn_mid;  /* elem edge i+1 */
      i++;
    }
    else if (Edge[i].un.node.i == Edge[i].un.node.j) /* colapsed quad. elements */
    {
      conn[(Edge[i+0].elem*9)+2+(Edge[i+0].ciclo*2)] = Edge[i].un.node.i;
      nn_mid--;
    }
    else
    {
      conn[(Edge[i+0].elem*9)+2+(Edge[i+0].ciclo*2)] = *n_node + nn_mid;
      if ((i+1) == (nedge-1))
      {
         conn[(Edge[i+1].elem*9)+2+(Edge[i+1].ciclo*2)] = *n_node + nn_mid;
         nn_mid++;
      }
    }
    nn_mid++;
  }
  /* special case */
  if (Edge[nedge-1].un.node.i == Edge[nedge-1].un.node.j) /* colapsed quad. elements */
  {
    conn[(Edge[nedge-1].elem*9)+2+(Edge[nedge-1].ciclo*2)] = Edge[nedge-1].un.node.i;
  }

  free(*Conn); *Conn = conn;

  /****** evaluate the midle coordenates *******/
  (*coords) = realloc ((*coords), ((*n_node)+nn_mid) * 2 * sizeof (double));
  nn_mid = 0;
  for(i = 0; i < nedge-1; i++)
  {
    if(Edge[i].un.code == Edge[i+1].un.code)
    {
      (*coords)[((*n_node)+nn_mid)*2+0] = ((*coords)[Edge[i].un.node.i*2+0] + 
                                           (*coords)[Edge[i].un.node.j*2+0]) / 2.0;
      (*coords)[((*n_node)+nn_mid)*2+1]  =((*coords)[Edge[i].un.node.i*2+1] + 
                                           (*coords)[Edge[i].un.node.j*2+1]) / 2.0;
    }
    else if (Edge[i].un.node.i != Edge[i].un.node.j)
    {
      if (!Edge[i].flag_change)
        mid_id = Edge[i].un.node.i*2 + 1;
      else
        mid_id = Edge[i].un.node.j*2 + 1;

      (*coords)[((*n_node)+nn_mid)*2+0] = Bdry_pts[mid_id*2];
      (*coords)[((*n_node)+nn_mid)*2+1] = Bdry_pts[mid_id*2+1];
    }

    if(Edge[i].un.code == Edge[i+1].un.code)
    {
      i++;
    }
    else if (Edge[i].un.node.i == Edge[i].un.node.j)
    {
      nn_mid--;
    }
    else if((i+1) == (nedge-1))
    {
      (*coords)[((*n_node)+nn_mid)*2+0] = ((*coords)[Edge[i+1].un.node.i*2+0] +
                                           (*coords)[Edge[i+1].un.node.j*2+0]) / 2.0;
      (*coords)[((*n_node)+nn_mid)*2+1] = ((*coords)[Edge[i+1].un.node.i*2+1] +
                                           (*coords)[Edge[i+1].un.node.j*2+1]) / 2.0;
      nn_mid++;
    }

    nn_mid++;
  }

  *n_node = *n_node + nn_mid;
 
  free (Edge);

}


#endif
               
