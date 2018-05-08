
//
// Test driver for region meshing.  This program reads a boundary
// description from stdin and generates a mesh
//
// the boundary description file has the following format:
//
// number_of_nodes
// node_id_0 node_x_coord_0 node_y_coord_0
// node_id_1 node_x_coord_1 node_y_coord_1
// ...
// number_of_edges
// edge_0_node_id_0 edge_0_node_id_1 [edge_0_node_id_2] <- only for quadratic
// edge_1_node_id_0 edge_1_node_id_1 [edge_1_node_id_2]
// ...

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


extern "C"
{
#include "mshq2d.h"
#include "msh2d.h"
}

#include "ArbMsh.hpp"
#include "ArbMshRegion2D.hpp"
#include "ArbHashTable.hpp"

#ifdef MEMDEBUG
#include "MemDbg.hpp"
#define new new(__FILE__,__LINE__)
#endif

extern int print_num ;
extern int quit_num ;
extern bool do_labels ;


void Msh2d_write_neutralfile_teste(
char *name_file_teste,
int n_elem,
int n_node,
double *Coords,
int    *conn
);


int Msh2DSeamGeneration (
int     n_pts,         /* # of points                          (in) */
double  *bdry_pts,     /* coordinates of all points            (in) */
int     bound_edge,    /* # of boundary edges                  (in) */
int     inter_edge,    /* # of internal free edges             (in) */
int     *edges_bound,  /* vector with the edges                (in) */
int     type_mesh,     /* Ps: only are linear elements (3, 4)
                          0..working     (in) */ 
int     *n_node,       /* counts # of pts for meshing         (out) */
double  **coords,      /* coordinate array used for meshing   (out) */
int     *n_elem,       /* number of elements generated        (out) */
int     **Conn         /* elem.connectivity list from meshing (out) */
)
{
    CArbMshRegion2D  *region;
    ArbMshNode       *nodes;
    ArbMshEdge       *edges;
    bool             quad_flg = false;
    bool             bf_flag = false;
    bool             sm_flag = true;
    bool             ql_flag = false;
    bool             ls_flag = false;
    bool             tc_flag = true;
    bool             qt_dbg_flag = false;
    double           bf ;
    int              i, status, num_edge;

    do_labels = false ;

    //  -ql     - do quadrilateral elements
    //  -quad   - quadratic order elements
    //  -bf val - set near boundary factor 
    //  -ns     - turns of all smoothing
    //  -pn     - start debugging print number
    //  -qn     - debugging quit number
    //  -labels - print debugging labels
    //  -ls     - force Laplace smoothing
    //  -ntc    - turn off topological cleanup
    //  -qtd    - turn on quad tree debugging
    //  -dm     - dumps the mesh to a .dmp file

    if (type_mesh == 6 || type_mesh == 8)
      quad_flg = true; // quadratic order elements
    else
      quad_flg = false;
    
    if (type_mesh == 4 || type_mesh == 8)
      ql_flag = true; // do quadrilateral elements
    else
      ql_flag = false;

    bf_flag = true; // set near boundary factor 
    //bf = 0.7; 
    bf = 0.3; 
 
    quit_num = false;

    do_labels = false; // print debugging labels

    ls_flag = true;  //  force Laplace smoothing

    tc_flag = true; // turn off topological cleanup

    qt_dbg_flag = false;  // turn on quad tree debugging

    ////////////////////////////////////////
    if (quad_flg)  /* quadratic elements */
    {
      if (ql_flag)  /* do quadrilateral elements */
          region = new CArbMshRegion2D(QUADRATIC,QUADRILATERAL) ;
      else         /* do triangular elements */
          region = new CArbMshRegion2D(QUADRATIC) ;
    } 
    else /* linear elements */
    {
      if (ql_flag) /* do quadrilateral elements */
          region = new CArbMshRegion2D(LINEAR,QUADRILATERAL) ;
      else         /* do triangular elements */
          region = new CArbMshRegion2D() ;
    }

    if (ls_flag) /* force Laplace smoothing */
      region->SetWinslowSmoothing(true) ;
    else
      region->SetWinslowSmoothing(false) ;
    
    if (tc_flag) /* turn off topological cleanup */
      region->SetTopoCleanup(true) ;
    else
      region->SetTopoCleanup(false) ;
    
    if (qt_dbg_flag) /* turn on quad tree debugging */
    {
      region->SetQuadTreeDebug() ;
    }

//    region->SetDebugDisplayFlags(CArbMshRegion2D::RmshAll) ;

#ifdef MEMDEBUG
    MemDbg *dbg_mem = GetMemDbg() ;
    dbg_mem->set_overwrite_check() ;
#endif

    // read the nodes
    nodes = new ArbMshNode[n_pts];
    for (i = 0; i < n_pts; ++i)
    {
        nodes[i].id = i;
        nodes[i].coord[0] = bdry_pts[i*2+0];
        nodes[i].coord[1] = bdry_pts[i*2+1];
        nodes[i].coord[2] = 0.0;
    }

    // read edges
    if (quad_flg == true)
    {
      num_edge = (bound_edge + inter_edge) / 2;
      edges = new ArbMshEdge[num_edge];
      for (i=0 ; i<num_edge*2 ; i+=2)
      {
          edges[i/2].node_id[0] = edges_bound[((i+1)*2)%(num_edge*4)+1];
          edges[i/2].node_id[1] = edges_bound[i*2+0];
          edges[i/2].node_id[2] = edges_bound[i*2+1];  // mid
      }
    }
    else
    {
      num_edge = bound_edge + inter_edge;
      edges = new ArbMshEdge[num_edge];
      for (i=0 ; i<num_edge ; ++i)
      {
          edges[i].node_id[0] = edges_bound[i*2+1];
          edges[i].node_id[1] = edges_bound[i*2+0];
          edges[i].node_id[2] = -1; // edges_bound[i*3+2];
      }
    }

    // Add node list to region
    // printf ("adicinou nos\n");
    status = region->AddNodeList (n_pts, nodes);
    // printf ("adicinou nos\n");
    switch (status) 
    {
        case ARB_NORMAL_STATUS:
            break ;

        case ARB_DUPLICATE_NODE_ID:
            fprintf (stderr, "Duplicate Node Exception\n");
            break ;

        case ARB_BAD_NODE_ID:
            fprintf (stderr, "Bad Node Exception\n");
            break ;
    }

    // Add edge list to region
    status = region->AddEdgeList (num_edge, edges);
    // printf ("adicinou arestas\n");
    switch (status)
    {
        case ARB_NORMAL_STATUS:
            break ;

        case ARB_DUPLICATE_EDGE:
            fprintf (stderr, "Duplicate Edge Exception\n");
            break ;
    }

    delete [] nodes;
    delete [] edges;

    if (bf_flag) region->SetNearBoundaryFactor (bf);
    if (!sm_flag) region->SetNodeSmoothing (false);
 

    // Generate Mesh
    ////////////////////////////////////////////////////
    if (quad_flg == false && type_mesh == 4)
    {
       int     tmp_n_node, tmp_n_elem, *tmp_conn;
       double  *tmp_coords;
       int     tmp_type_msh;
       
       tmp_type_msh = 3;
       
       Msh2DEdge (n_pts, bdry_pts, bound_edge, inter_edge, edges_bound, tmp_type_msh, 
                  &tmp_n_node, &tmp_coords,&tmp_n_elem, &tmp_conn);
       
       // printf ("Create Triangular Mesh.\n") ;
       // write_neutralfile_teste("d:\\temp\\triangulacao.dat", tmp_n_elem, tmp_n_node, tmp_coords, tmp_conn);
       
       region->SetTriangularMesh (tmp_n_node, tmp_coords, tmp_n_elem, tmp_conn);
       
       free (tmp_coords);
       free (tmp_conn);

    }

    status = region->GenerateMesh ( );
    ////////////////////////////////////////////////////
    printf ("Create Quadrilateral Mesh.\n") ;

    switch (status) 
    {
        case ARB_NORMAL_STATUS:
            break ;

        case ARB_ILLEGAL_BOUNDARY:
            fprintf(stderr,"Detected an illegal boundry.\n") ;
            fprintf(stderr,"Try again without checks\n") ;
            return 0;
            break ;
    }

    ArbMshElement2D *elems = region->GetElements() ;
    nodes = region->GetNodes() ;

    if (1)
    {
      int  nn;
      *n_node = region->NumNodes();
      *coords = (double *) calloc (region->NumNodes() * 2, sizeof (double));
      for (i=0 ; i<region->NumNodes() ; ++i)
      {
        (*coords)[nodes[i].id*2+0] = nodes[i].coord[0];
        (*coords)[nodes[i].id*2+1] = nodes[i].coord[1];
      }

      *n_elem = region->NumElements ();
       if (quad_flg == false)
         *Conn = (int *) calloc (region->NumElements() * 5, sizeof (int));
       else
         *Conn = (int *) calloc (region->NumElements() * 9, sizeof (int));

      for (i=0 ; i < region->NumElements(); i++)
      {
        nn = elems[i].num_nodes;
        if (quad_flg == false)
        {
          (*Conn)[i*5+0] = 4;
          if (nn == 4)
          {
            (*Conn)[i*5+1] = elems[i].nodes[3];
            (*Conn)[i*5+2] = elems[i].nodes[2];
            (*Conn)[i*5+3] = elems[i].nodes[1];
            (*Conn)[i*5+4] = elems[i].nodes[0];
          }
          else if (nn == 3)
          {
            (*Conn)[i*5+1] = elems[i].nodes[0];
            (*Conn)[i*5+2] = elems[i].nodes[2];
            (*Conn)[i*5+3] = elems[i].nodes[1];
            (*Conn)[i*5+4] = elems[i].nodes[0];
          }
        }
        else
        {
          (*Conn)[i*9+0] = 8;
          if (nn == 8)
          {
            (*Conn)[i*9+1] = elems[i].nodes[0];
            (*Conn)[i*9+2] = elems[i].nodes[7];
            (*Conn)[i*9+3] = elems[i].nodes[3];
            (*Conn)[i*9+4] = elems[i].nodes[6];
            (*Conn)[i*9+5] = elems[i].nodes[2];
            (*Conn)[i*9+6] = elems[i].nodes[5];
            (*Conn)[i*9+7] = elems[i].nodes[1];
            (*Conn)[i*9+8] = elems[i].nodes[4];
          }
          else if (nn == 6)
          {
            (*Conn)[i*9+1] = elems[i].nodes[0];
            (*Conn)[i*9+2] = elems[i].nodes[5];
            (*Conn)[i*9+3] = elems[i].nodes[2];
            (*Conn)[i*9+4] = elems[i].nodes[4];
            (*Conn)[i*9+5] = elems[i].nodes[1];
            (*Conn)[i*9+6] = elems[i].nodes[3];
            (*Conn)[i*9+7] = elems[i].nodes[0];
            (*Conn)[i*9+8] = elems[i].nodes[0];
          }
        }
      }

    }

#if 0
    fprintf(stderr,"Nodes : %d\n",region->NumNodes()) ;
    fprintf(stderr,"Elems : %d\n",region->NumElements()) ;

    ArbMshStats *stats = region->GetMeshStatistics() ;
    fprintf(stderr,"Mean Shape Measure    : %f\n",stats->mean_shape_measure) ;
    fprintf(stderr,"Minimum Shape Measure : %f\n",stats->minimum_shape_measure) ;
    
    for (i=1 ; i<11 ; ++i)
    {
        double x = double(i) / 10.0 ;
        fprintf(stderr,"num less than %g : %d\n", x, stats->num_less_than[i]) ;
    }
    delete [] stats ;

#endif


    delete [] nodes ;
    delete [] elems ;
    delete region ;

#ifdef MEMDEBUG
    int mem_cnt = dbg_mem->count_allocated() ;
    fprintf(stderr,"Allocated Memory: %d\n",mem_cnt) ;
    dbg_mem->print_allocated() ;
#endif

    // write_neutralfile_teste("D:\\temp\\quad.nf", *n_elem, *n_node, *coords, *Conn);

    return(1) ;
}




/* ------------------------------------------------------- */
void Msh2d_write_neutralfile_teste(
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
   fprintf(arquivo,"%d    %f    %f    0.0\n",i+1,Coords[i*2],Coords[i*2+1]);
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
       fprintf(arquivo,"\n");
     }
     index=index+conn[index]+1;
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
       fprintf(arquivo,"\n");
     }
     index=index+conn[index]+1;
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
       fprintf(arquivo,"\n");
     }
     index=index+conn[index]+1;
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
       fprintf(arquivo,"\n");
     }
     index=index+conn[index]+1;
   }
 }

 fprintf(arquivo,"\n%%END\n");
 fclose(arquivo);

}

