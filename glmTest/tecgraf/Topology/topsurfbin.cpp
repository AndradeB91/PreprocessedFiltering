/* topology.cpp: implementation of the topology class.
//
////////////////////////////////////////////////////////////////////*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

extern "C"
{
#include "topsurfbin.h"
}

#include "topology.hpp"


/* InsertMesh
////////////////////////////////////////////////////////////////////*/
void* SurfTopInsertMesh (int n_node, double *coords, int n_elem, int *conn,
                         void (*func) (const char *))
{
  cTopology* new_topo = new cTopology;

  new_topo->InsertMesh (n_node, n_elem, coords, conn, func);

  return (void *) new_topo;
}

/* SurfTopUdateNormalNodes 
////////////////////////////////////////////////////////////////////*/
void SurfTopUdateNormalNodes (void *surf)
{
  cTopology *curr_topol = (cTopology *) surf;
  curr_topol->UpdateNormals ( );
}

/* ReturnNodeId
////////////////////////////////////////////////////////////////////*/
int SurfTopReturnNodeId (void *surf, double x, double y, double z)
{
  cTopology *curr_topol = (cTopology *) surf;
  return curr_topol ->ReturnNodeId (x, y, z);
}


// SurfTopRelease 
//////////////////////////////////////////////////////////////////////
void SurfTopRelease (void *surf)
{
  cTopology *curr_topol = (cTopology *) surf;
  delete curr_topol;
}


// SurfTopBoundBox 
//////////////////////////////////////////////////////////////////////
void SurfTopBoundBox (void *surf, double min[3], double max[3])
{
  cTopology *topo = (cTopology *) surf;

  topo->BoundBox (min, max);
}

/* SurfTopGetPts
////////////////////////////////////////////////////////////////////*/
int SurfTopGetCoordNode (void *surf, int id, double coord[3])
{
  cTopology *topo = (cTopology *) surf;

  return (topo->GetCoordNode (id, coord));
}

/* SurfTopGetNormalNode
////////////////////////////////////////////////////////////////////*/
int SurfTopGetNormalNode (void *surf, int id, double normal[3])
{
  cTopology *topo = (cTopology *) surf;

  return (topo->GetNormalNode (id, normal));
}

/* SurfTopGetNormalNode
////////////////////////////////////////////////////////////////////*/
void SurfTopSetNormalNode (void *surf, int id, double normal[3])
{
  cTopology *topo = (cTopology *) surf;

  topo->SetNormalNode (id, normal);
}


/* SurfTopSetPtsInative 
////////////////////////////////////////////////////////////////////*/
void SurfTopSetPtsInative (void *surf, int id)
{
  cTopology *topo = (cTopology *) surf;

  topo->SetPtsInative (id);
}

/* SurfTopNumAdjNodeNode 
////////////////////////////////////////////////////////////////////*/
int SurfTopNumAdjNodeNode( void *surf, int id )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return 0;

  return topo->NumAdjNodeNode(id);
}

/* SurfTopAdjNodeNode
////////////////////////////////////////////////////////////////////*/
void SurfTopAdjNodeNode( void *surf, int id, int *adj_nodes )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return;
  topo->AdjNodeNode (id, adj_nodes);
}


/* AMTopOrientAdjNodeNode
////////////////////////////////////////////////////////////////////*/
int AMTopOrientAdjNodeNode (void *surf, int id, int *adj_nodes)
{
  cTopology *topo = (cTopology *) surf;

  return (topo->OrientAdjNodeNode (id, adj_nodes));
}


/* SurfTopIsBoundaryNode
////////////////////////////////////////////////////////////////////*/
int SurfTopIsBoundaryNode (void *surf, int id)
{
  cTopology *topo = (cTopology *) surf;

  return (topo->IsBoundaryNode (id));
}


/* SurfTopSetCoordNode
////////////////////////////////////////////////////////////////////*/
void SurfTopSetCoordNode  (void *surf, int id, double *new_coord)
{
  cTopology *topo = (cTopology *) surf;

  topo->SetCoordNode (id, new_coord);
}

/* SurfTopIsActNode
////////////////////////////////////////////////////////////////////*/
int SurfTopIsActNode (void *surf, int id)
{
  cTopology *topo = (cTopology *) surf;

  return (topo->IsActNode (id));
}





/* SurfTopGetEdge 
////////////////////////////////////////////////////////////////////*/
int SurfTopGetEdge (void *surf, int i, int *ei, int *ej)
{
  cTopology *topo = (cTopology *) surf;

  return topo->GetEdge (i, ei, ej);
}

/* SurfTopGetEdgeSize 
////////////////////////////////////////////////////////////////////*/
int SurfTopGetEdgeSize (void *surf, int i, double mid[3], double *size)
{
  cTopology *topo = (cTopology *) surf;

  return topo->GetEdgeSize (i, mid, size);
}

/* SurfTopGetBoundEdge 
////////////////////////////////////////////////////////////////////*/
int SurfTopGetBoundEdge (void *surf, int i, int *ei, int *ej)
{
  cTopology *topo = (cTopology *) surf;

  return topo->GetBoundEdge (i, ei, ej);
}


/* SurfTopGetBEdgeSize 
////////////////////////////////////////////////////////////////////*/
int SurfTopGetBEdgeSize (void *surf, int i, double mid[3], double *size)
{
  cTopology *topo = (cTopology *) surf;

  return topo->GetBEdgeSize (i, mid, size);
}

/* SurfTopGetIsEdgeValid 
////////////////////////////////////////////////////////////////////*/
int SurfTopGetIsEdgeValid (void *surf, int i)
{
  cTopology *topo = (cTopology *) surf;

  return topo->GetIsEdgeValid (i);
}


/* SurfTopGetSmalestEdge
////////////////////////////////////////////////////////////////////*/
double SurfTopGetSmalestEdge (void *surf)
{
  cTopology *topo = (cTopology *) surf;

  return topo->GetSmalestEdge ( );
}



/* SurfTopNumEdge 
////////////////////////////////////////////////////////////////////*/
int SurfTopNumEdge (void *surf)
{
  cTopology *topo = (cTopology *) surf;

  if (topo == NULL)
    return 0;
  return (topo->NumEdges ());
}

/* SurfTopNumBoundEdge
////////////////////////////////////////////////////////////////////*/
int SurfTopNumBoundEdge (void *surf)
{
  cTopology *topo = (cTopology *) surf;

  if (topo == NULL)
    return 0;
  return (topo->NumBoundEdges ());
}


// SurfTopNumNodes 
//////////////////////////////////////////////////////////////////////
int SurfTopNumNodes (void *surf)
{
  cTopology *topo = (cTopology *) surf;

  if (topo == NULL)
    return 0;
  return (topo->NumNodes ());
}

/* SurfTopNumElems
////////////////////////////////////////////////////////////////////*/
int SurfTopNumElems (void *surf)
{
  cTopology *topo = (cTopology *) surf;

  if (topo == NULL)
    return 0;
  return (topo->NumElems ());
}



/* SurfTopIsActElem
////////////////////////////////////////////////////////////////////*/
int SurfTopIsActElem (void *surf, int id)
{
  cTopology *topo = (cTopology *) surf;

  if (topo == NULL)
    return 0;

  return (topo->IsActElem (id));
}


// SurfTopGetElemNNodes 
//////////////////////////////////////////////////////////////////////
int SurfTopGetElemNNodes (void *surf, int id)
{
  cTopology *topo = (cTopology *) surf;

  if (topo == NULL)
    return 0;

  return (topo->GetElemNNodes (id));
}

/* SurfTopGetElemConn
/////////////////////////////////////////////////////////////////////*/
int* SurfTopGetElemConn  (void *surf, int id)
{
  cTopology *topo = (cTopology *) surf;

  if (topo == NULL)
    return NULL;

  return (topo->GetElemConn (id));
}

/* SurfTopGetElemNorm  
/////////////////////////////////////////////////////////////////////*/
int SurfTopGetElemNorm  (void *surf, int id, double nornal[3])
{
  cTopology *topo = (cTopology *) surf;

  if (topo == NULL)
    return 0;

  return (topo->GetElemNorm (id, nornal));
}

/* SurfTopGetElemIdEdge
/////////////////////////////////////////////////////////////////////*/
int SurfTopGetElemIdEdge (void *surf, int id, int pos_adj)
{
 cTopology *topo = (cTopology *) surf;

  if (topo == NULL)
    return -1;

  return (topo->GetElemIdEdge (id, pos_adj));
}

/* SurfTopGetIsElemValid
/////////////////////////////////////////////////////////////////////*/
int SurfTopGetIsElemValid (void *surf, int id)
{
 cTopology *topo = (cTopology *) surf;

  if (topo == NULL)
    return 0;

  return (topo->GetIsElemValid (id));
}

/* SurfTopGetElemCenter
/////////////////////////////////////////////////////////////////////*/
int SurfTopGetElemCenter (void *surf, int id, double center[3])
{
  cTopology *topo = (cTopology *) surf;
 
  if (topo == NULL)
    return 0;

  return topo->GetElemCenter (id, center);
}


/* SurfTopGetNumLoops 
/////////////////////////////////////////////////////////////////////*/
int SurfTopGetNumLoops (void *surf)
{
  cTopology *topo = (cTopology *) surf;

  if (topo == NULL)
    return 0;

  return (topo->GetNumLoops ());
}

/* SurfTopGetLoop
/////////////////////////////////////////////////////////////////////*/
int SurfTopGetLoop (void *surf, int id, int *npts, int **idpts)
{
  cTopology *topo = (cTopology *) surf;

  if (topo == NULL)
    return 0;

  return topo->GetLoop (id, npts, idpts);
}

#if 0

/* SurfTopFindEdgeNodeI
/////////////////////////////////////////////////////////////////////*/
/* pos = 0 => i, pos = 1 => j */ 
static EdgeBtree *SurfTopFindInElem_EdgeWithNodeId (Elm *curr_elem, int id, int pos)
{
  EdgeBtree  *curr_edge;
  int        k, n = curr_elem->type;

  /* find the oposite element */
  for (k = 0; k < n; k++)
  {
    curr_edge = curr_elem->edge[k];
    if (curr_edge->id[(curr_elem->ciclo[k]+pos)%2] == id)
      return curr_edge;
  }

  return NULL;
}

/* SurfTopGetAdjacentElemToEdge
/////////////////////////////////////////////////////////////////////*/
static Elm *SurfTopGetAdjacentElemToEdge (Elm *curr_elem, EdgeBtree *edge)
{
  if (edge->AdjElm[0] == curr_elem)
    return edge->AdjElm[1];
  else
    return edge->AdjElm[0];
}

/* SurfTopGeEdgePosInElem
/////////////////////////////////////////////////////////////////////*/
static int SurfTopGeEdgePosInElem (Elm *elem, EdgeBtree *edge)
{
  int i;
  for (i = 0; i < elem->type; i++)
  {
    if (elem->edge[i] == edge)
      return i;
  }

  return -1;
}

/* SurfTopGetNetNodes
/////////////////////////////////////////////////////////////////////*/
/*
   Net:

   *---------*---------*---------*---------*  \
   |         |         |         |         |  |
   |         |         |         |         |  |
   |         |         |         |         |  |
   |         |         |         |         |  |
   *---------*---------*---------*---------*   > size (=2)
   |         |         |         |         |  |
   |         |         |         |         |  |
   |         |         |         |         |  |
   |         |         |id       |         |  |
   *---------*---------O---------*---------*  /
   |         |         |         |         |
   |         |   /->   |         |         |
   |         |  |      |         |         |
   |         |   \_    |         |         |
   *---------*---------*---------*---------*
   |         |         |         |         |
   |         |         |         |         |
   |         |         |         |         |
   |         |         |         |         |
   o---------*---------*---------*---------*
   first node
*/

int SurfTopGetNetNodes (void *surf, int id, int size, int **net)
{
  cTopology *topo = (cTopology *) surf;
  Elm        *curr_elem, *base_elem;
  EdgeBtree  *curr_edge;
  int        i, k, curr_node, base_node, pos;


  if (topo == NULL)
    return 0;
  if (id > topo->NumNode-1)
    return 0;

  /* curr element */
  curr_elem = topo->NodeVector[id].AdjEdges[0]->AdjElm[0];
  if (!curr_elem->valid)
    return 0;

  if (topo->NodeVector[id].n_adj != 4)
    return 0;

  /* find the first node */
  curr_node = id;
  for (i = 0; i < size; i++)
  {
    /* find the opposite element */
    for (k = 0; k < 2; k++)
    {
      /* find a edge with node id with i (0) position */
      curr_edge = SurfTopFindInElem_EdgeWithNodeId (curr_elem, curr_node, 0);

      /* move to adjacent element */
      curr_elem = SurfTopGetAdjacentElemToEdge (curr_elem, curr_edge);

      if (!curr_elem->valid)
        return 0;
    }

    /* find a edge with node id with i (0) position */
    curr_edge = SurfTopFindInElem_EdgeWithNodeId (curr_elem, curr_node, 0);

    /* position of edge in curr_elem adj vector */
    pos = SurfTopGeEdgePosInElem (curr_elem, curr_edge);

    /* get the next edge */
    pos = (pos+1)%curr_elem->type;
    curr_edge = curr_elem->edge[pos];

    /* get the new corner node => j position*/
    curr_node = curr_edge->id[(curr_elem->ciclo[pos]+1)%2];
  }

  /* fill net ids */
  base_node = curr_node;
  base_elem = curr_elem;
  for (i = 0; i < 2*size; i++)
  {
    net[i][0] = curr_node = base_node;
    curr_elem = base_elem;

    for (k = 1; k <= 2*size; k++)
    {
      /* find a edge with node id with i (0) position */
      curr_edge = SurfTopFindInElem_EdgeWithNodeId (curr_elem, curr_node, 0);
      /* position of edge in curr_elem adj vector */
      pos = SurfTopGeEdgePosInElem (curr_elem, curr_edge);

      /* get the next edge */
      pos = (pos+1)%curr_elem->type;
      curr_edge = curr_elem->edge[pos];

      /* get the new corner node => i position*/
      net[i][k] = curr_node = curr_edge->id[curr_elem->ciclo[pos]];

      /* last column */
      if (i == 2*size-1)
        net[i+1][k] = curr_edge->id[(curr_elem->ciclo[pos]+1)%2];

      /* update current element => move to adjacent element */
      curr_elem = SurfTopGetAdjacentElemToEdge (curr_elem, curr_edge);
      if (!curr_elem->valid)
        return 0;
    }

    /* update base element and node */

    /* find a edge with node id with j (1) position */
    curr_edge = SurfTopFindInElem_EdgeWithNodeId (base_elem, base_node, 1);
    pos = SurfTopGeEdgePosInElem (base_elem, curr_edge);

    /* get i (first node) => base node */
    base_node = curr_edge->id[base_elem->ciclo[pos]];

    /* find a edge with node id with j (1) position */
    curr_edge = SurfTopFindInElem_EdgeWithNodeId (base_elem, base_node, 1);

    /* update current element => move to adjacent element */
    base_elem = SurfTopGetAdjacentElemToEdge (base_elem, curr_edge);

    /* last column */
    if (i == 2*size-1)
      net[i+1][0] = base_node;
  }

  /* check if internal nodes are valid */
  for (i = 1; i < 2*size; i++)
  {
    for (k = 1; k < 2*size; k++)
    {
      if (topo->NodeVector[net[i][k]].n_adj != 4)
        return 0;    
    }
  }

  return 1;
}

#endif

/* SurfTopOrientAdjNodeNode
/////////////////////////////////////////////////////////////////////*/
int SurfTopOrientAdjNodeNode( void *surf, int id, int *adj_nodes )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return 0;
  
  return topo->OrientAdjNodeNode (id, adj_nodes);
}

/* SurfTopNumAdjEdgeToNode
/////////////////////////////////////////////////////////////////////*/
int SurfTopNumAdjEdgeToNode( void *surf, int id )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return 0;
  return topo->NumAdjEdgeToNode (id);
}

/* SurfTopAdjEdgeToNode
/////////////////////////////////////////////////////////////////////*/
int SurfTopAdjEdgeToNode( void *surf, int id, int /*pos*/ )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return 0;
  return topo->NumAdjEdgeToNode (id);
}

/* SurfTopOppAdjEdgeToNode
/////////////////////////////////////////////////////////////////////*/
int SurfTopOppAdjEdgeToNode( void *surf, int id, int *ne, int **edges )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return 0;
  return topo->OppAdjEdgeToNode(id, ne, edges);
}

/* SurfTopAdjElemToNode
/////////////////////////////////////////////////////////////////////*/
int SurfTopAdjElemToNode( void *surf, int id, int *ne, int **elem )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return 0;
  return topo->AdjElemToNode (id, ne, elem);
}

/* SurfTopAdjNodeToEdge
/////////////////////////////////////////////////////////////////////*/
int SurfTopAdjNodeToEdge( void *surf, int id, int *ei, int *ej )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return 0;

  return topo->AdjNodeToEdge (id, ei, ej);
}

/* SurfTopNumAdjElemToEdge
/////////////////////////////////////////////////////////////////////*/
int SurfTopNumAdjElemToEdge( void *surf, int id )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return 0;

  return topo->NumAdjElemToEdge (id);
}

/* SurfTopAdjElemToEdge
/////////////////////////////////////////////////////////////////////*/
int SurfTopAdjElemToEdge( void *surf, int id, int pos )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return 0;

  return topo->AdjElemToEdge (id, pos);
}

/* SurfTopNumAdjEdgeToElem
/////////////////////////////////////////////////////////////////////*/
int SurfTopNumAdjEdgeToElem( void *surf, int id )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return 0;

  return topo->NumAdjEdgeToElem (id);
}

/* SurfTopAdjEdgeToElem
/////////////////////////////////////////////////////////////////////*/
int SurfTopAdjEdgeToElem( void *surf, int id, int pos )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return 0;

  return topo->AdjEdgeToElem (id, pos);
}

/* writeNF
/////////////////////////////////////////////////////////////////////*/
int writeNF( void *surf, char *filename )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return 0;

  return topo->writeNF(filename);
}

/////////////////////////////////////////////////////////////////////*/
void SurfTopSetElemInfo( void *surf, int id, int flag )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return;
  topo->SetElemInfo (id, flag);
}

int SurfTopGetElemInfo( void *surf, int id )
{
  cTopology *topo = (cTopology *) surf;
  if (topo == NULL)
    return -1;
  return topo->GetElemInfo (id);
}

