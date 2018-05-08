// topology.cpp: implementation of the topology class.
//
//////////////////////////////////////////////////////////////////////

/* #include <time.h>  */

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <map>
#include <cstring>

#ifdef _UNIX_
#include <sys/time.h>
#else
#include <sys/timeb.h>
#endif

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1.0e+06
#endif


#include "Vec3D.hpp"
#include "topology.hpp"
#include "mgIndexMap.hpp"

#include "topology_struc.hpp"

using namespace FTools;
using namespace std;


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

// cTopology
//////////////////////////////////////////////////////////////////////
cTopology::cTopology()
{
  NodeVector = NULL;
  NumNode = 0;

  Elem = NULL;
  NumElem = 0;

  EdgeVector = NULL;
  NumEdge = 0;
  NumSurfaces = 0;
  feedbackFunction = NULL;
  EdgeList = NULL;
}

// ~cTopology
//////////////////////////////////////////////////////////////////////
cTopology::~cTopology()
{
  Clear ( );
}

// Clean
//////////////////////////////////////////////////////////////////////
void cTopology::Clear ()
{
  if (NodeVector != NULL)
    delete []NodeVector;
  NodeVector = NULL;
  NumNode = 0;

  if (Elem != NULL)
    delete []Elem;
  Elem = NULL;
  NumElem = 0;

  int i;
  for (i = 0; i < NumEdge; ++i)
    delete EdgeVector[i];
  if (EdgeVector != 0)
    delete []EdgeVector;
  EdgeVector = NULL;
  NumEdge = 0;

  BoundList.clear ( );

  if (EdgeList != NULL)
    EdgeList->clear ( );

  NumSurfaces = 0;
  feedbackFunction = NULL;

  m_avgEdgeSize = 0.0;
  m_minEdgeSize = 0.0;
}


// InsertMesh
//////////////////////////////////////////////////////////////////////
int cTopology::InsertMesh (int n_node, int n_elem, double *coords, int *conn,
                           void (*func) (const char *))
{
  int        i, j, index;

  Clear ( );

  feedbackFunction = func;

  NumNode = n_node;
  NumElem = n_elem;

  if (feedbackFunction != NULL)
    feedbackFunction ("Filling node and element vectors...");

  if (n_node == 0 || n_elem == 0) return 0;

  // fill nodes
  NodeVector = new TopNode_[n_node];
  for (i = 0; i < n_node; ++i)
  {
    NodeVector[i].coords[0] = coords[i*3+0];
    NodeVector[i].coords[1] = coords[i*3+1];
    NodeVector[i].coords[2] = coords[i*3+2];
    NodeVector[i].active    = 1;
  }

  // boudbox
  for (j = 0; j < 3; j++)
    max[j] = min[j] = NodeVector[0].coords[j];
  for (i = 1; i < n_node; i++)
  {
    for (j = 0; j < 3; j++)
    {
      if (NodeVector[i].coords[j] < min[j])
        min[j] = NodeVector[i].coords[j];
      if (NodeVector[i].coords[j] > max[j])
        max[j] = NodeVector[i].coords[j];
    }
  }


  // fill elements
  Elem = new TopElem_[n_elem];
  index = 0;
  for (i = 0; i < n_elem; i++)
  {
    Elem[i].type = conn[index];
    // percorre os nos
    for (j = 0; j < conn[index]; j++)
      Elem[i].inc[j] = conn[index+1+j];

    Elem[i].active = 1;
    Elem[i].valid  = 1;
    Elem[i].id     = i;
    index = index + conn[index] + 1;
  }

  if (feedbackFunction != NULL)
    feedbackFunction ("Building and filling edge vector...");
  // Build the list of edges
  BuildEdgeList ( );

  if (feedbackFunction != NULL)
    feedbackFunction ("Building and filling adjacent edges to nodes...");
  // Adjacent Edges to Node
  BuildAdjEdgesToNode ( );

  if (feedbackFunction != NULL)
    feedbackFunction ("Identifing differente surfaces...");
  // Identify different surfaces
  //IdentifySurfaces ( );

  BuildOrientation ( );

  if (feedbackFunction != NULL)
    feedbackFunction ("Identifing boundary edges...");
  // Build the boundary edge list
  BuildEdgesBoundary ( );

  // SetValidFaces ( );

  if (feedbackFunction != NULL)
    feedbackFunction ("Computing normals...");
  // Build normals of elements and nodes
  UpdateNormals ( );


  // set local tolerance
  SetLocalToler ( );

  if (feedbackFunction != NULL)
    feedbackFunction (" ");

  return 1;
}


// BoundBox
//////////////////////////////////////////////////////////////////////
void cTopology::BoundBox (double bmin[3], double bmax[3])
{
  int i;
  for (i = 0; i < 3; i++)
  {
    bmin[i] = min[i];
    bmax[i] = max[i];
  }
}


// BuildNormals
//////////////////////////////////////////////////////////////////////
void cTopology::UpdateNormals (void)
{
  int    i, j, *n_adj, id_node, n;
  Vec3D  u, v, w, center;

  //normal of elements
  for (i = 0; i < NumElem; ++i)
  {
    if (Elem[i].type < 3)
    {
      continue;
    }
    else if (Elem[i].type == 3)
    {
    u = NodeVector[Elem[i].inc[1]].coords - NodeVector[Elem[i].inc[0]].coords;
    v = NodeVector[Elem[i].inc[2]].coords - NodeVector[Elem[i].inc[0]].coords;
    w = CrossProd (u, v);
    }
    else
    {
      // compute center
      center.x(0); center.y(0); center.z(0);
      n = Elem[i].type;
      for (j = 0; j < n; ++j)
      {
        id_node = Elem[i].inc[j];
        center += NodeVector[id_node].coords;
      }
      center /= n;

      // compute normal
      w.x(0.0); w.y(0.0); w.z(0.0);
      for (j = 0; j < n; ++j)
      {
        u = NodeVector[Elem[i].inc[j]].coords - center;
        v = NodeVector[Elem[i].inc[(j+1)%n]].coords - center;
        w += CrossProd (u, v);
      }
      w /= n;
    }

    Elem[i].normal = w.Normalize ( );
  }

  /* update normal to points */
  n_adj = new int[NumNode];
  memset (n_adj, 0, NumNode * sizeof (int));
  for (i = 0; i < NumNode; ++i)
  {
    n_adj[i] = 0;
    NodeVector[i].normal[0] = 0.0;
    NodeVector[i].normal[1] = 0.0;
    NodeVector[i].normal[2] = 0.0;
  }

  for (i = 0; i < NumElem; ++i)
  {
    if (Elem[i].type < 3)
      continue;

    for (j = 0; j < Elem[i].type; ++j)
    {
      id_node = Elem[i].inc[j];
      NodeVector[id_node].normal += Elem[i].normal;
      n_adj[id_node]++;
    }
  }

  /* average */
  for (i = 0; i < NumNode; i++)
    NodeVector[i].normal /= n_adj[i];

  delete []n_adj;

  // normal of nodes

}

//////////////////////////////////////////////////////////////////////
int cTopology::GetCoordNode (int id, double coord[3])
{
  if (id > NumNode-1)
    return 0;
  coord[0] = NodeVector[id].coords[0];
  coord[1] = NodeVector[id].coords[1];
  coord[2] = NodeVector[id].coords[2];
  return (1);
}

//////////////////////////////////////////////////////////////////////
void cTopology::SetCoordNode (int id, double *new_coord)
{
  if (id > NumNode-1)
    return;

  NodeVector[id].coords[0] = new_coord[0];
  NodeVector[id].coords[1] = new_coord[1];
  NodeVector[id].coords[2] = new_coord[2];
}

//////////////////////////////////////////////////////////////////////
int  cTopology::IsBoundaryNode (int id)
{
  if (id > NumNode-1)
    return 0;

  return NodeVector[id].bound;
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetNormalNode (int id, double normal[3])
{
  if (id > NumNode-1)
    return 0;

  normal[0] = NodeVector[id].normal[0];
  normal[1] = NodeVector[id].normal[1];
  normal[2] = NodeVector[id].normal[2];
  return 1;
}

//////////////////////////////////////////////////////////////////////
void cTopology::SetNormalNode (int id, double normal[3])
{
  if (id > NumNode-1)
    return;

  NodeVector[id].normal[0] = normal[0];
  NodeVector[id].normal[1] = normal[1];
  NodeVector[id].normal[2] = normal[2];
  return;
}


//////////////////////////////////////////////////////////////////////
void cTopology::SetPtsInative (int id)
{
  if (id > NumNode-1)
    return;

  NodeVector[id].active = 0;
}

//////////////////////////////////////////////////////////////////////
int cTopology::NumAdjNodeNode (int id)
{
  if (id > NumNode-1)
    return 0;
  return ((int) NodeVector[id].AdjEdges.size());
}

//////////////////////////////////////////////////////////////////////
void cTopology::AdjNodeNode (int id, int *adj_nodes)
{
  if (id > NumNode-1)
    return;

  int size = (int) NodeVector[id].AdjEdges.size();
  int i;
  for (i = 0; i < size; ++i)
  {
    if (NodeVector[id].AdjEdges[i]->id[0] != id)
      adj_nodes[i] = NodeVector[id].AdjEdges[i]->id[0];
    else
      adj_nodes[i] = NodeVector[id].AdjEdges[i]->id[1];
  }
}

//////////////////////////////////////////////////////////////////////
int cTopology::OrientAdjNodeNode (int id, int *adj_nodes)
{
  if (id > NumNode-1)
    return 0;

  int size = (int) NodeVector[id].AdjEdges.size();
  int i;

  // current edge
  int curr_edge = -1;
  if (NodeVector[id].bound)  // get las
  {
    for (int j = 0; j < (int) NodeVector[id].AdjEdges.size(); ++j)
    {
      if (NodeVector[id].AdjEdges[j]->id[0] == id)
      {
        curr_edge = NodeVector[id].AdjEdges[j]->i;
        break;
      }
    }
  }
  else  // any edge can be selected
  {
    curr_edge = NodeVector[id].AdjEdges[0]->i;
  }

  // find adjnodes
  for (i = 0; i < size; ++i)
  {
    if (EdgeVector[curr_edge]->id[0] != id)
      adj_nodes[i] = EdgeVector[curr_edge]->id[0];
    else
      adj_nodes[i] = EdgeVector[curr_edge]->id[1];

    if (i == size-1) // avoid last search
      break;

    curr_edge = getAdjEdgeNextNodeEdge (id, curr_edge);
  }


  return 1;
}


//////////////////////////////////////////////////////////////////////
int cTopology::OrientAdjEdgeToNode (int id, int *adj_edges)
{
  if (id > NumNode-1)
    return 0;

  int size = (int) NodeVector[id].AdjEdges.size();
  int i;

  // current edge
  int curr_edge = -1;
  if (NodeVector[id].bound)  // get las
  {
    TopEdge_ *ptr_edge;
    for (int j = 0; j < (int) NodeVector[id].AdjEdges.size(); ++j)
    {
      ptr_edge = NodeVector[id].AdjEdges[j];
      if (ptr_edge->id[0] == id && ptr_edge->type == 1)
      {
        curr_edge = ptr_edge->i;
        break;
      }
    }
  }
  else  // any edge can be selected
  {
    curr_edge = NodeVector[id].AdjEdges[0]->i;
  }

  // find adjnodes
  for (i = 0; i < size; ++i)
  {
    adj_edges[i] = curr_edge;

    if (i == size-1) // avoid last search
      break;

    curr_edge = getAdjEdgeNextNodeEdge (id, curr_edge);
  }


  return 1;
}

//////////////////////////////////////////////////////////////////////
int  cTopology::IsActNode  (int id)
{
  if (id > NumNode-1)
    return 0;

  return (NodeVector[id].active);
}


// ReturnNodeId
//////////////////////////////////////////////////////////////////////
int cTopology::ReturnNodeId (double x, double y, double z)
{
  Vec3D pts (x, y, z);

  int    i, idmin = -1;
  double dist, mindist;

  mindist = (NodeVector[0].coords - pts).Magnitude();
  idmin = 0;

  for (i = 1; i < NumNode; ++i)
  {
    dist = (NodeVector[i].coords - pts).Magnitude();
    if (dist < mindist)
    {
    	mindist = dist;
      idmin = i;
    }
  }

 return idmin;
}

//////////////////////////////////////////////////////////////////////
void cTopology::SetNodeInfo (int id, int flag)
{
  if (id > NumNode-1)
    return;
  NodeVector[id].extraInfo = flag;
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetNodeInfo (int id)
{
  if (id > NumNode-1)
    return 0;
  return (NodeVector[id].extraInfo);
}


//////////////////////////////////////////////////////////////////////
int cTopology::GetEdge (int i, int *ei, int *ej)
{
  if (i > NumEdge-1)
    return 0;

  *ei = EdgeVector[i]->id[0];
  *ej = EdgeVector[i]->id[1];
  return 1;
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetEdge( int ei, int ej )
{
  if (ei > NumNode-1)
    return -1;
  if (ej > NumNode-1)
    return -1;

  int n_adj_edge = NodeVector[ei].AdjEdges.size();
  for (int i = 0; i < n_adj_edge; ++i)
  {
    if (NodeVector[ei].AdjEdges[i]->id[0] == ej ||
        NodeVector[ei].AdjEdges[i]->id[1] == ej )
        return NodeVector[ei].AdjEdges[i]->i;
  }


  return -1;
}


//////////////////////////////////////////////////////////////////////
bool cTopology::IsBoundaryEdge( int id )
{
  if (id > NumEdge-1)
    return false;
  if (EdgeVector[id]->type == 1)
    return true;
  return false;
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetEdgeSize (int i, double mid[3], double *size)
{
  if (i > NumEdge-1)
    return 0;
  Vec3D  middle (0.0, 0.0, 0.0);

  middle = (NodeVector[EdgeVector[i]->id[0]].coords + NodeVector[EdgeVector[i]->id[1]].coords) * 0.5;
  mid[0] = middle[0];
  mid[1] = middle[1];
  mid[2] = middle[2];
  *size = EdgeVector[i]->size;
  return 1;
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetBoundEdge (int i, int *ei, int *ej)
{
  if (i > NumEdgeBound-1)
    return 0;

  *ei = BoundVector[i]->id[0];
  *ej = BoundVector[i]->id[1];
  return 1;
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetBEdgeSize (int i, double mid[3], double *size)
{
  if (i > NumEdgeBound-1)
    return 0;
  Vec3D  middle (0.0, 0.0, 0.0);

  middle = (NodeVector[BoundVector[i]->id[0]].coords +
            NodeVector[BoundVector[i]->id[1]].coords) * 0.5;
  mid[0] = middle[0];
  mid[1] = middle[1];
  mid[2] = middle[2];
  *size = BoundVector[i]->size;
  return 1;
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetIsEdgeValid (int i)
{
  if (i > NumEdge-1)
    return 0;
  return EdgeVector[i]->valid;
}

//////////////////////////////////////////////////////////////////////
double cTopology::GetSmalestEdge (void)
{
  return m_minEdgeSize;
}

//////////////////////////////////////////////////////////////////////
double cTopology::GetAverageEdge (void)
  {
  return m_avgEdgeSize;
  }

double cTopology::GetLongestEdge (void)
{
  return m_maxEdgeSize;
}



//////////////////////////////////////////////////////////////////////
void cTopology::SetEdgeInfo (int id, int flag)
{
  if (id > NumEdge-1)
    return;
  EdgeVector[id]->extraInfo = flag;
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetEdgeInfo (int id)
{
  if (id > NumEdge-1)
    return 0;
  return (EdgeVector[id]->extraInfo);
}


//////////////////////////////////////////////////////////////////////
int cTopology::GetElemNNodes (int id)
{
  if (id > NumElem-1)
    return 0;

  return Elem[id].type;
}

//////////////////////////////////////////////////////////////////////
int* cTopology::GetElemConn (int id)
{
  if (id > NumElem-1)
    return 0;

  return Elem[id].inc;
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetElemNorm (int id, double normal[3])
{
  if (id > NumElem-1)
    return 0;

  normal[0] = Elem[id].normal[0];
  normal[1] = Elem[id].normal[1];
  normal[2] = Elem[id].normal[2];

  return 1;
}

//
//////////////////////////////////////////////////////////////////////////
int cTopology::GetElemNormArea( int id, double normal[3] )
{
  int    j, n, id_node;

  if (id > NumElem-1)
    return 0;

  Vec3D  u, v, w, center;

  if (Elem[id].type == 3)
  {
  u = NodeVector[Elem[id].inc[1]].coords - NodeVector[Elem[id].inc[0]].coords;
  v = NodeVector[Elem[id].inc[2]].coords - NodeVector[Elem[id].inc[0]].coords;
  w = CrossProd (u, v);
  }
  else
  {
    // compute center
    center.x(0); center.y(0); center.z(0);
    n = Elem[id].type;
    for (j = 0; j < n; ++j)
    {
      id_node = Elem[id].inc[j];
      center += NodeVector[id_node].coords;
    }
    center /= n;

    // compute normal
    w.x(0.0); w.y(0.0); w.z(0.0);
    for (j = 0; j < n; ++j)
    {
      u = NodeVector[Elem[id].inc[j]].coords - center;
      v = NodeVector[Elem[id].inc[(j+1)%n]].coords - center;
      w += CrossProd (u, v);
    }
  }

  normal[0] = w[0];
  normal[1] = w[1];
  normal[2] = w[2];
  return 1;
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetElemIdEdge (int id, int pos_adj)
{
  if (id > NumElem-1)
    return 0;
  if (pos_adj >= Elem[id].type)
    return -1;

  return (Elem[id].AdjEdges[pos_adj]->i);
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetIsElemValid (int id)
{
  if (id > NumElem-1)
    return 0;

  return Elem[id].valid;
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetElemCenter (int id, double center[3])
{
  int i;
  if (id > NumElem-1)
    return 0;
  Vec3D  middle (0,0,0);

  for (i = 0; i < Elem[id].type; ++i)
  {
    middle += NodeVector[Elem[id].inc[i]].coords;
  }
  middle /= Elem[id].type;

  center[0] = middle[0];
  center[1] = middle[1];
  center[2] = middle[2];

  return 1;
}

//////////////////////////////////////////////////////////////////////
int cTopology::IsActElem (int id)
{
  if (id > NumElem-1)
    return 0;

  return Elem[id].active;
}


//////////////////////////////////////////////////////////////////////
void cTopology::SetElemInfo (int id, int flag)
{
  if (id > NumElem-1)
    return;
  Elem[id].extraInfo = flag;
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetElemInfo (int id)
{
  if (id > NumElem-1)
    return 0;
  return (Elem[id].extraInfo);
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetNumLoops (void)
{
  return (int) BoundList.size ();
}

//////////////////////////////////////////////////////////////////////
int cTopology::GetLoop (int id, int *npts, int **idpts)
{
  int i;
  int nloops = (int) BoundList.size ();

  if (id >= nloops)
    return 0;

  int  size = (int) BoundList[id]->Edges.size();

  *npts = size;
  *idpts = new int[size];

  for (i = 0; i < size; ++i)
    (*idpts)[i] = BoundList[id]->Edges[i]->id[0];

  return 1;
}


// NumAdjEdgeToNode
//////////////////////////////////////////////////////////////////////
int cTopology::NumAdjEdgeToNode (int id)
{
  if (id > NumNode-1)
    return 0;

  return (int) NodeVector[id].AdjEdges.size();
}

// AdjEdgeToNode
//////////////////////////////////////////////////////////////////////
int cTopology::AdjEdgeToNode (int id, int pos)
{
  if (id > NumNode-1)
    return -1;

  if ((int) NodeVector[id].AdjEdges.size() == 0 ||
      (int) NodeVector[id].AdjEdges.size() < pos+1 )
    return -1;

  return NodeVector[id].AdjEdges[pos]->i;
}

// OppAdjEdgeToNode
//////////////////////////////////////////////////////////////////////
int cTopology::OppAdjEdgeToNode (int id, int *ne, int **edges)
{
  int i, j, ei, ej;

  if (id > NumNode-1)
    return 0;

  if ((int) NodeVector[id].AdjEdges.size() == 0)
    return 0;

  int nelm, *elem;
  vector <int> tmpedges;
  AdjElemToNode (id, &nelm, &elem);
  for (i = 0; i < nelm; ++i)
  {
    int nadjedge = NumAdjEdgeToElem (elem[i]);
    for (j = 0; j < nadjedge; ++j)
    {
      int adjedge = AdjEdgeToElem (elem[i], j);
      AdjNodeToEdge (adjedge, &ei, &ej);
      if (!(ei == id || ej == id))
        tmpedges.push_back (adjedge);
    }
  }
  delete []elem;

  // get adjacent edges
  if (tmpedges.size () == 0)
  {
    *ne = 0;
    *edges = NULL;
    return 0;
  }

  *ne = (int) tmpedges.size ();
  *edges = new int[*ne];
  for (i = 0; i < *ne; ++i)
    (*edges)[i] = tmpedges[i];

  return 1;
}

// GetAdjNodeToEdge
//////////////////////////////////////////////////////////////////////
int cTopology::AdjNodeToEdge (int id, int *ei, int *ej)
{
  if (id > NumEdge-1)
    return 0;

  *ei = EdgeVector[id]->id[0];
  *ej = EdgeVector[id]->id[1];
  return 1;
}

// AdjElemToNode
//////////////////////////////////////////////////////////////////////
int cTopology::AdjElemToNode (int id, int *ne, int **elem)
{
  if (id > NumNode-1)
    return 0;

  map <int, int> adj_elem;

  // for all adjacent edges
  int i, j, nedges = NumAdjEdgeToNode (id);
  for (i = 0; i < nedges; ++i)
  {
    int adjedge  = AdjEdgeToNode (id, i);
    int nadjelem = NumAdjElemToEdge (adjedge);
    for (j = 0; j < nadjelem; ++j)
    {
      int adjelem = AdjElemToEdge (adjedge, j);
      if (adj_elem.find (adjelem) == adj_elem.end())
        adj_elem.insert (pair <int, int> (adjelem, adjelem));
    }
  }

  // any adj_elems
  if (adj_elem.size() == 0)
  {
    *ne = 0;
    *elem = NULL;
    return 0;
  }

  *ne = (int) adj_elem.size();
  *elem = new int[*ne];
  map <int,int>::iterator ii;
  for (ii = adj_elem.begin(), i = 0; ii != adj_elem.end(); ++i, ++ii)
  {
    (*elem)[i] = (*ii).second;
  }

  return 1;
}


// NumAdjElemToEdge
//////////////////////////////////////////////////////////////////////
int cTopology::NumAdjElemToEdge (int id)
{
  if (id > NumEdge-1)
    return 0;
  return (int) EdgeVector[id]->AdjElms.size();
}

// AdjElemToEdge
//////////////////////////////////////////////////////////////////////
int cTopology::AdjElemToEdge (int id, int pos)
{
  if (id > NumEdge-1)
    return -1;
  int size = (int) EdgeVector[id]->AdjElms.size();
  if (size == 0 || size < pos+1)
    return -1;

  return EdgeVector[id]->AdjElms[pos]->id;
}

// NumAdjEdgeToElem
//////////////////////////////////////////////////////////////////////
int cTopology::NumAdjEdgeToElem (int id)
{
  if (id > NumElem-1)
    return 0;
  return (int) Elem[id].AdjEdges.size ();
}

// AdjEdgeToElem
//////////////////////////////////////////////////////////////////////
int cTopology::AdjEdgeToElem (int id, int pos)
{
  if (id > NumElem-1)
    return -1;

  if ((int) Elem[id].AdjEdges.size () == 0 ||
      (int) Elem[id].AdjEdges.size () < pos+1 )
   return -1;

  return Elem[id].AdjEdges[pos]->i;
}

// GetVolume
//////////////////////////////////////////////////////////////////////
double cTopology::GetVolume ( )
{
  int i;
  Vec3D p0, p1, p2;
  double volume = 0.0;

  // loop for all edges in the elements
  for (i = 0; i < NumElem; i++)
  {
    p0 = NodeVector[Elem[i].inc[0]].coords;
    p1 = NodeVector[Elem[i].inc[1]].coords;
    p2 = NodeVector[Elem[i].inc[2]].coords;

    volume = volume + 1.0 / 6.0 * TripleProd (p0, p1, p2);

    if (Elem[i].type == 4)
    {
      p0 = NodeVector[Elem[i].inc[0]].coords;
      p1 = NodeVector[Elem[i].inc[2]].coords;
      p2 = NodeVector[Elem[i].inc[3]].coords;

      volume = volume + 1.0 / 6.0 * TripleProd (p0, p1, p2);
    }

  }

  return volume;
}



//////////////////////////////////////////////////////////////////////
// Private methods
//////////////////////////////////////////////////////////////////////

#include <set>

//////////////////////////////////////////////////////////////////////////
// local class to store adjacent nodes
class mshAdjNodes
{
public:
  mshAdjNodes(void) {};
  ~mshAdjNodes(void)
  {
    m_adjNodes.clear();
  }

  void insert (int node, int elem)
  {
    itr = m_adjNodes.find(node);
    if (itr == m_adjNodes.end())
    {
      vector <int> new_v;
      new_v.push_back (elem);
      m_adjNodes.insert(make_pair(node, new_v));
    }
    else
    {
      (*itr).second.push_back (elem);
    }
  }

  int size (void)
  {
    return (int) m_adjNodes.size();
  }

  map < int, vector <int> > m_adjNodes;

private:
  map < int, vector <int> >::iterator itr;
};


// BuildEdgeList
//////////////////////////////////////////////////////////////////////
int cTopology::BuildEdgeList (void)
{
  int      id1, id2, id3 = -1;
  int      i, j, curr_step = 0;
  TopEdge_  *currEdge;
  Vec3D    vEdge;
  double   steps = 0.05;
  char     message[126];

#if 0

  if (NumNode == 0)
    return 0;

  double cpu_time = clock( );

  //////////////////////////////////////////////////////////////////////////
  mshAdjNodes *adj_nodes = new mshAdjNodes[NumNode];

  // insert
  for (i = 0; i < NumElem; i++)
  {
    for (j = 0; j < Elem[i].type; j++)
    {
      id1 = Elem[i].inc[j];
      id2 = Elem[i].inc[(j+1)%Elem[i].type];
      if (id1 < id2)
        adj_nodes[id1].insert(id2);
      else
        adj_nodes[id2].insert(id1);
    }
  }

  // count the number of edges
  int n_edges = 0;
  for (i = 0; i < NumNode; i++)
    n_edges += adj_nodes[i].size();
  printf ("Numero de arestas = %d\n", n_edges);

  cpu_time = (clock( ) - cpu_time)/(CLOCKS_PER_SEC*1.0);
  printf("\t\tCPU time Edge 1 ............. %0.5f (s)\n", (double)cpu_time);




    cpu_time = clock( );


#endif

  // clear edge list
  EdgeList = new mshsurfIndexMap;
  EdgeList->clear ( );

  // loop for all edges in the elements
  for (i = 0; i < NumElem; i++)
  {
    // print to feedback function
    if (feedbackFunction != NULL)
    {
      if (i*1.0/NumElem > curr_step * steps)
      {
        sprintf (message, "Building and filling edge vector (%3.1f %%)...", curr_step * steps * 100.0);
        feedbackFunction (message);
        ++curr_step;
      }
    }

    // for all
    for (j = 0; j < Elem[i].type; j++)
    {
      id1 = Elem[i].inc[j];
      id2 = Elem[i].inc[(j+1)%Elem[i].type];

      currEdge = (TopEdge_ *) EdgeList->find (id1, id2, id3);

      // try to find a edge
      if (currEdge == NULL)  // new edge
      {
        TopEdge_  *newEdge = new TopEdge_;

        newEdge->id[0] = id1;
        newEdge->id[1] = id2;
        newEdge->valid = 1;
        newEdge->type  = 0;

        // size
        vEdge = NodeVector[id1].coords - NodeVector[id2].coords;
        newEdge->size = vEdge.Magnitude ( );

        EdgeList->insert (id1, id2, id3, (void *)newEdge);

        currEdge = newEdge;
      }
      // update type
      currEdge->type = (currEdge->type + 1);
//      cout << "Edge type = " << i << " " << currEdge->type << endl;


      // update adjacent elements
      currEdge->AdjElms.push_back (&(Elem[i]));

      // update adjacent edge
      Elem[i].AdjEdges.push_back (currEdge);
      if (id1 == currEdge->id[0])
        Elem[i].cycles.push_back (0);
      else
        Elem[i].cycles.push_back (1);
    }
  }

  //cpu_time = (clock( ) - cpu_time)/(CLOCKS_PER_SEC*1.0);
  //printf("\t\tCPU time Edge 2 ............. %0.5f (s)\n", (double)cpu_time);

  // fill edge vector
  NumEdge = EdgeList->size ();
  // printf ("Numero de arestas = %d\n", NumEdge);

  EdgeVector = new TopEdge_*[NumEdge];
  m_avgEdgeSize = m_minEdgeSize = m_maxEdgeSize = 0.0;
  for (currEdge = (TopEdge_ *) EdgeList->first (), i = 0; currEdge != NULL;
       currEdge = (TopEdge_ *) EdgeList->next (), ++i)
  {
    EdgeVector[i] = currEdge;
    currEdge->i = i;
    m_avgEdgeSize += currEdge->size;
    if ((m_minEdgeSize == 0.0 || currEdge->size < m_minEdgeSize) && currEdge->size > 0.0)
      m_minEdgeSize = currEdge->size;
    if (currEdge->size > m_maxEdgeSize)
      m_maxEdgeSize = currEdge->size;
  }

  // Average edge size
  if (NumEdge > 0)
    m_avgEdgeSize /= NumEdge;


  return 1;
}


// BuildAdjEdgesToNode
//////////////////////////////////////////////////////////////////////
int cTopology::BuildAdjEdgesToNode (void)
{
  int i, id0, id1;

  for (i = 0; i < NumEdge; ++i)
  {
    id0 = EdgeVector[i]->id[0];
    id1 = EdgeVector[i]->id[1];
    NodeVector[id0].AdjEdges.push_back (EdgeVector[i]);
    NodeVector[id1].AdjEdges.push_back (EdgeVector[i]);
    if (EdgeVector[i]->type == 1)
    {
      NodeVector[id0].bound = 1;
      NodeVector[id1].bound = 1;
    }
  }

  return 1;
}

// BuildAdjEdgesToNode
//////////////////////////////////////////////////////////////////////
int cTopology::BuildEdgesBoundary  (void)
{
  int i, j, k;

  // look for a boundary edge
  for (i = 0; i < NumEdge; i++)
  {
    if (EdgeVector[i]->flag == 0  &&  EdgeVector[i]->type == 1)
    {
      // cout << "Uma aresta de contorno" << endl;
      RunAlongEdges (EdgeVector[i]);
    }
  }

  // fill boundary edge vector
  int nloops = (int) BoundList.size ();
  NumEdgeBound = 0;
  for (i = 0; i < nloops; ++i)
    NumEdgeBound += (int) BoundList[i]->Edges.size();

  BoundVector = new TopEdge_*[NumEdgeBound];

   for (i = 0, k = 0; i < nloops; ++i)
   {
     for (j = 0; j < (int) BoundList[i]->Edges.size(); ++j, ++k)
        BoundVector[k] = BoundList[i]->Edges[j];
   }

  return 1;
}


// RunAlongEdges
//////////////////////////////////////////////////////////////////////
void cTopology::RunAlongEdges (TopEdge_ *initEdge)
{
  int k;

  TopEdge_ *currEdge, *adjEdge;
  TopBoundary *newBoundary = new TopBoundary;  // new bounday

  // Add initial edge
  newBoundary->Edges.push_back (initEdge);

  currEdge = initEdge;
  currEdge->flag = 1;

  // obtain the boundary from initial edge
  do
  {
    int id_j, n_edge_adj = 0;

    /* pega o no na posicao de j */
    id_j = currEdge->id[1];

    // count how many adjacent boundary edges to current point
    for (k = 0; k < (int) NodeVector[id_j].AdjEdges.size(); ++k)
    {
      adjEdge = NodeVector[id_j].AdjEdges[k];
      if (adjEdge->type == 1 && adjEdge != currEdge && adjEdge->flag == 0)
        ++n_edge_adj;
    }

    // it must find only one edge
    if (n_edge_adj != 1)
      break;

    // update the current edge
    for (k = 0; k < (int) NodeVector[id_j].AdjEdges.size(); ++k)
    {
      adjEdge = NodeVector[id_j].AdjEdges[k];
      if (adjEdge->type == 1 && adjEdge != currEdge  && adjEdge->flag == 0)
        currEdge = adjEdge;
    }

    currEdge->flag = 1;

    // add current edge to current boundary
    if (currEdge != initEdge)
      newBoundary->Edges.push_back (currEdge);

  } while (currEdge != initEdge);

  // hold current edge boundary
  BoundList.push_back (newBoundary);
}


// SetLocalToler
//////////////////////////////////////////////////////////////////////
void cTopology::SetLocalToler (void)
{
  int i;
  LocalTol = 0.0;

  for (i = 0; i < NumEdge; ++i)
  {
    if (EdgeVector[i]->size > 0.0 && (LocalTol == 0.0 || EdgeVector[i]->size < LocalTol))
      LocalTol = EdgeVector[i]->size;
  }

  LocalTol *= 0.01;
}


// IdentifySurfaces
//////////////////////////////////////////////////////////////////////
void cTopology::IdentifySurfaces ( )
{
  int i;

  NumSurfaces = 0;

  // reset all surfaces in elements
  for (i = 0; i < NumElem; ++i)
    Elem[i].surfId = -1;

  // reset all surfaces in elements
  for (i = 0; i < NumElem; ++i)
  {
    if (Elem[i].surfId == -1)
      SetIdSurface (NumSurfaces, &(Elem[i]));
    ++NumSurfaces;
  }
}

// SetIdSurface
//////////////////////////////////////////////////////////////////////
void cTopology::SetIdSurface (int idSurf, TopElem_ *initElem)
{
  TopEdge_ *currEdge;
  TopElem_ *currElem;
  deque<TopEdge_ *> dequeEdges;

  int i;

  // reset flags in the edges
  for (i = 0; i < (int) NumEdge; ++i)
    EdgeVector[i]->flag = 0;

  // add all edges in initial element
  for (i = 0; i < (int) initElem->AdjEdges.size(); ++i)
  {
    if (initElem->AdjEdges[i]->type == 2)
      dequeEdges.push_back (initElem->AdjEdges[i]);
  }

  currElem = initElem;
  currElem->surfId = idSurf;

  //
  while (!dequeEdges.empty ( ))
  {
    // get the first edge in deque
    currEdge = dequeEdges.front ( );

    // Set flag
    currEdge->flag = 1;

    // set adjacent elements to edge

    // add adjacent edges to deque
    currElem = currEdge->AdjElms[0];
    if (currElem->surfId == -1)
    {
      for (i = 0; i < (int) currElem->AdjEdges.size(); ++i)
      {
        if ((currElem->AdjEdges[i]->type == 2) &&
            (currElem->AdjEdges[i]->flag == 0) )
          dequeEdges.push_back (currElem->AdjEdges[i]);
      }
      currEdge->AdjElms[0]->surfId = idSurf;
    }

    currElem = currEdge->AdjElms[1];
    if (currElem->surfId == -1)
    {
      for (i = 0; i < (int) currElem->AdjEdges.size(); ++i)
      {
        if ((currElem->AdjEdges[i]->type == 2) &&
            (currElem->AdjEdges[i]->flag == 0) )
          dequeEdges.push_back (currElem->AdjEdges[i]);
      }
      currEdge->AdjElms[1]->surfId = idSurf;
    }

    // remove the first edge from deque
    dequeEdges.pop_front ( );
  }
}


// SetValidFaces
//////////////////////////////////////////////////////////////////////
void cTopology::SetValidFaces (void)
{
  int i, j;

  // reset all surfaces in elements
  for (i = 0; i < NumElem; ++i)
    Elem[i].valid = 1;

  // set all faces as valid
  for (i = 0; i < NumEdge; ++i)
  {
    if (EdgeVector[i]->AdjElms.size() > 2)
    {
      for (j = 0; j < (int) EdgeVector[i]->AdjElms.size(); ++j)
         EdgeVector[i]->AdjElms[j]->valid = 0;
    }
  }
}


// BuildOrientation
//////////////////////////////////////////////////////////////////////
void cTopology::BuildOrientation (void)
{
  int      i, k/*, ok = 1*/;
  TopElem_ *head_elem = NULL, *tail_elem;

  // set all element as valid
  for (i = 0; i < NumElem; ++i)
    Elem[i].valid = 1;

  // try to any inconsistence in the mesh
  for (i = 0; i < NumEdge; ++i)
  {
    if (EdgeVector[i]->AdjElms.size() > 2)
    {
      EdgeVector[i]->valid = 0;
      for (k = 0; k < (int) EdgeVector[i]->AdjElms.size(); ++k)
        EdgeVector[i]->AdjElms[k]->valid = 0;
    }
    else
      EdgeVector[i]->valid = 1;
  }

  //if (!ok)
  //return;


  /* set all element flags to zero */
  for (i = 0; i < NumElem; ++i)
    Elem[i].flag = 0;

  while (true)
  {
    // look for a element with flag = 0
    for (i = 0; i < NumElem; ++i)
    {
      if (Elem[i].flag == 0)
      {
        head_elem = tail_elem = &(Elem[i]);
        break;
      }
    }
    if (i == NumElem)
      break;

    /* from the head element, it orients the all elements */
    do
    {
      head_elem->flag = 1;

      for (i = 0; i < head_elem->type; i++)
        OrientAdjElement (head_elem->AdjEdges[i], head_elem->cycles[i], &tail_elem);

      head_elem = head_elem->next;

    } while (head_elem != NULL);

  }


}


/* OrientAdjElement
////////////////////////////////////////////////////////////////////*/
void cTopology::OrientAdjElement (TopEdge_ *edge, int ciclo, TopElem_ **tail_elem)
{
  TopElem_  *curr_elem;
  int       i, pos = -1, n;
  TopEdge_  *tmp_edge;

  if (edge->type == 1) /* boundary */
    return;

  /* both adjacent elements already were passed */
  if (edge->AdjElms[0]->flag && edge->AdjElms[1]->flag)
    return;

  /* get the element to be checked */
  if (edge->AdjElms[0]->flag == 0)
    curr_elem = edge->AdjElms[0];
  else
    curr_elem = edge->AdjElms[1];
  n = curr_elem->type;

  /* find edge position */
  for (i = 0; i < n; i++)
  {
    if (edge == curr_elem->AdjEdges[i])
    {
      pos = i;
      break;
    }
  }

  /* check ciclo */
  if (curr_elem->cycles[pos] == ciclo)
  {
    int tmp;
    /* invert orientation */
    for (i = 0; i < (n-1)/2; i++)
    {
      /* incidency */
      tmp                   = curr_elem->inc[i+1];
      curr_elem->inc[i+1]   = curr_elem->inc[n-i-1];
      curr_elem->inc[n-i-1] = tmp;
    }

     /* invert edge orientation */
    for (i = 0; i < n/2; i++)
    {
#if 1
      /* ciclo */
      tmp                      = curr_elem->cycles[i];
      curr_elem->cycles[i]     = curr_elem->cycles[n-i-1];
      curr_elem->cycles[n-i-1] = tmp;
      /* edge */
      tmp_edge                   = curr_elem->AdjEdges[i];
      curr_elem->AdjEdges[i]     = curr_elem->AdjEdges[n-i-1];
      curr_elem->AdjEdges[n-i-1] = tmp_edge;
#endif
    }

    for (i = 0; i < n; i++)
    {
      curr_elem->cycles[i] = !curr_elem->cycles[i]; /* invert ciclo */
      if (curr_elem->AdjEdges[i]->type == 1) /* invert id boundary */
      {
        tmp = curr_elem->AdjEdges[i]->id[0];
        curr_elem->AdjEdges[i]->id[0] = curr_elem->AdjEdges[i]->id[1];
        curr_elem->AdjEdges[i]->id[1] = tmp;
      }
    }
  }

  curr_elem->flag = 1;  /* set flag */
  /* orient adj elements */
  (*tail_elem)->next = curr_elem;
  (*tail_elem) = curr_elem;

}

//////////////////////////////////////////////////////////////////////////
int cTopology::getAdjEdgeNextNodeEdge (int node, int edge)
{
  int n_adj_elem = (int) EdgeVector[edge]->AdjElms.size();

  if (n_adj_elem > 2)
    return -1;

  // next node adjacent to node
  int next_node;
  if (EdgeVector[edge]->id[0] == node)
    next_node = EdgeVector[edge]->id[1];
  else
    next_node = EdgeVector[edge]->id[0];


  // find element

  for (int i = 0; i < n_adj_elem; ++i)
  {
    int adj_elem = EdgeVector[edge]->AdjElms[i]->id;
    // for all edges
    int n_edges = (int) Elem[adj_elem].AdjEdges.size();
    for (int j = 0; j < n_edges; ++j)
    {
      int prev = (j+n_edges-1)%n_edges;
      int next = (j+1)%n_edges;
      int curr = Elem[adj_elem].inc[j];
      int curr_next = Elem[adj_elem].inc[next];
      //int curr_prev = Elem[adj_elem].inc[prev];


      // find by the corner node and next node
      if (curr == node &&  curr_next == next_node)
        return (Elem[adj_elem].AdjEdges[prev]->i);

      //curr_prev = -1;
    }
  }
  return -1;
}


//////////////////////////////////////////////////////////////////////////
int cTopology::writeNF( char *filename )
{
  //int index =0, j, i;
  int j, i;
  FILE *arquivo;
  int numq4 = 0, numt3 = 0, numq8 = 0, numt6 = 0, num_gen = 0;

  arquivo = fopen(filename,"w");
  if (arquivo == NULL)
    return 0;

  fprintf(arquivo,"\n%%HEADER\nFile created by program 'topology' at - AMiranda\n");
  fprintf(arquivo,"\n%%HEADER.VERSION\n1.00\n");
  fprintf(arquivo,"\n%%HEADER.TITLE\n' untitled'\n");

  fprintf(arquivo,"%%NODE\n");
  fprintf(arquivo,"%d\n\n",NumNode);
  fprintf(arquivo,"%%NODE.COORD\n");
  fprintf(arquivo,"%d\n\n",NumNode);

  for(i=0; i<NumNode; i++)
  {
    fprintf(arquivo,"%d    %f    %f    %f\n",i+1, NodeVector[i].coords[0],
            NodeVector[i].coords[1], NodeVector[i].coords[2]);
  }

  fprintf(arquivo,"\n%%MATERIAL\n1\n");
  fprintf(arquivo,"\n%%MATERIAL.LABEL\n1\n1\t'mat1'\n");
  fprintf(arquivo,"\n%%MATERIAL.ISOTROPIC\n1\n1\t10000\t0.2\n");
  fprintf(arquivo,"\n%%THICKNESS\n1\n1\t1\n");
  fprintf(arquivo,"\n%%INTEGRATION.ORDER\n1\n1\t2\t2\t2\t2\t2\t2\n");


  //index = 0;
  for(i=0; i<NumElem; i++)
  {
    if (Elem[i].type == 3)
      numt3++;
    else if (Elem[i].type == 4)
      numq4++;
    else if (Elem[i].type == 6)
      numt6++;
    else if (Elem[i].type == 8)
      numq8++;
    else
      num_gen++;

  }


  fprintf(arquivo,"%%ELEMENT\n");
  fprintf(arquivo,"%d\n\n",NumElem);

  //index = 0;
  if (numt3 > 0)
  {
    fprintf(arquivo,"%%ELEMENT.T3\n");
    fprintf(arquivo,"%d\n", numt3);
    for(i=0; i<NumElem; i++)
    {
      if (Elem[i].type == 3)
      {
        fprintf(arquivo,"%d  1  1  1",i+1);
        for(j=0; j<Elem[i].type; j++)
          fprintf(arquivo,"  %d",Elem[i].inc[j] +1);
        fprintf(arquivo,"\n");
      }
    }
  }

  //index = 0;
  if (numt6 > 0)
  {
    fprintf(arquivo,"%%ELEMENT.T6\n");
    fprintf(arquivo,"%d\n", numt6);
    for(i=0; i<NumElem; i++)
    {
      if (Elem[i].type == 6)
      {
        fprintf(arquivo,"%d  1  1  1",i+1);
        for(j=0; j<Elem[i].type; j++)
          fprintf(arquivo,"  %d",Elem[i].inc[j]+1);
        fprintf(arquivo,"\n");
      }
    }
  }

  //index = 0;
  if (numq4 > 0)
  {
    fprintf(arquivo,"%%ELEMENT.Q4\n");
    fprintf(arquivo,"%d\n", numq4);
    for(i=0; i<NumElem; i++)
    {
      if (Elem[i].type == 4)
      {
        fprintf(arquivo,"%d  1  1  1",i+1);
        for(j=0; j<Elem[i].type; j++)
          fprintf(arquivo,"  %d",Elem[i].inc[j]+1);
        fprintf(arquivo,"\n");
      }
    }
  }

  //index = 0;
  if (numq8 > 0)
  {
    fprintf(arquivo,"%%ELEMENT.Q8\n");
    fprintf(arquivo,"%d\n", numq8);
    for(i=0; i<NumElem; i++)
    {
      if (Elem[i].type == 8)
      {
        fprintf(arquivo,"%d  1  1  1",i+1);
        for(j=0; j<Elem[i].type; j++)
          fprintf(arquivo,"  %d",Elem[i].inc[j]+1);
        fprintf(arquivo,"\n");
      }
      }
  }

  //index = 0;
  if (num_gen > 0)
  {
    fprintf(arquivo,"%%ELEMENT.GEN\n");
    fprintf(arquivo,"%d\n", num_gen);
    for(i=0; i<NumElem; i++)
    {
      if (Elem[i].type != 3 && Elem[i].type != 6 && Elem[i].type != 4 && Elem[i].type != 8)
      {
        fprintf(arquivo,"%d  1  1  1",i+1);
        for(j=0; j<Elem[i].type; j++)
          fprintf(arquivo,"  %d",Elem[i].inc[j]+1);
      fprintf(arquivo,"\n");
    }
  }
  }

  fprintf(arquivo,"\n%%END\n");
  fclose(arquivo);

  return 1;
}

//
//////////////////////////////////////////////////////////////////////////
int cTopology::writeMatlab( char *filename, char *open_type)
{
  //int index =0, i;
  int i;
  FILE *f;
  int num_gen = 0;

  f = fopen(filename, open_type);
  if (f == NULL)
    return 0;

  fprintf (f, "X = [ \n");
  for (i = 0; i < NumNode; ++i)
    fprintf (f, "     %10.5f\n", NodeVector[i].coords[0]);
  fprintf (f, "] \n\n");

  fprintf (f, "Y = [ \n");
  for (i = 0; i < NumNode; ++i)
    fprintf (f, "     %10.5f\n", NodeVector[i].coords[1]);
  fprintf (f, "] \n\n");

  fprintf (f, "Z = [ \n");
  for (i = 0; i < NumNode; ++i)
    fprintf (f, "     %10.5f\n", NodeVector[i].coords[2]);
  fprintf (f, "] \n\n");

  fprintf (f, "C = [ \n");


  //index = 0;
  for(i=0; i<NumElem; i++)
  {
    if (Elem[i].type == 3)
    {
      fprintf (f, "     %5d %5d %5d; \n", Elem[i].inc[0]+1, Elem[i].inc[1]+1, Elem[i].inc[2]+1);
    }
    else if (Elem[i].type == 4)
    {
      fprintf (f, "     %5d %5d %5d; \n", Elem[i].inc[0]+1, Elem[i].inc[1]+1, Elem[i].inc[2]+1);
      fprintf (f, "     %5d %5d %5d; \n", Elem[i].inc[0]+1, Elem[i].inc[2]+1, Elem[i].inc[3]+1);
    }
    else if (Elem[i].type == 6)
    {
      fprintf (f, "     %5d %5d %5d; \n", Elem[i].inc[0]+1, Elem[i].inc[1]+1, Elem[i].inc[5]+1);
      fprintf (f, "     %5d %5d %5d; \n", Elem[i].inc[1]+1, Elem[i].inc[2]+1, Elem[i].inc[3]+1);
      fprintf (f, "     %5d %5d %5d; \n", Elem[i].inc[3]+1, Elem[i].inc[4]+1, Elem[i].inc[5]+1);
      fprintf (f, "     %5d %5d %5d; \n", Elem[i].inc[1]+1, Elem[i].inc[3]+1, Elem[i].inc[5]+1);
    }
    else if (Elem[i].type == 8)
    {
      fprintf (f, "     %5d %5d %5d; \n", Elem[i].inc[0]+1, Elem[i].inc[1]+1, Elem[i].inc[7]+1);
      fprintf (f, "     %5d %5d %5d; \n", Elem[i].inc[1]+1, Elem[i].inc[2]+1, Elem[i].inc[3]+1);
      fprintf (f, "     %5d %5d %5d; \n", Elem[i].inc[3]+1, Elem[i].inc[4]+1, Elem[i].inc[5]+1);
      fprintf (f, "     %5d %5d %5d; \n", Elem[i].inc[5]+1, Elem[i].inc[6]+1, Elem[i].inc[7]+1);
      fprintf (f, "     %5d %5d %5d; \n", Elem[i].inc[1]+1, Elem[i].inc[3]+1, Elem[i].inc[7]+1);
      fprintf (f, "     %5d %5d %5d; \n", Elem[i].inc[3]+1, Elem[i].inc[5]+1, Elem[i].inc[7]+1);
    }
    else
      num_gen++;
  }

  fprintf (f, "] \n\n");

  fprintf (f, "trisurf (C, X, Y, Z)\n");

  fprintf (f, "hold on \n\n");

  fclose(f);

  return 1;
}

//
//////////////////////////////////////////////////////////////////////////
int cTopology::WriteVTK( char *filename )
{
  int j, i;
  FILE *arquivo = NULL;

  arquivo = fopen(filename,"w");
  if (arquivo == NULL)
    return 0;

  fprintf(arquivo,"# vtk DataFile Version 2.0\n");
  fprintf(arquivo,"Modelos com template 3D\n");
  fprintf(arquivo,"ASCII\n\n");

  fprintf(arquivo,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(arquivo,"POINTS %d float\n",NumNode);

  for(i=0; i<NumNode; i++)
  {
    fprintf(arquivo,"%f    %f    %f\n", NodeVector[i].coords[0],
      NodeVector[i].coords[1], NodeVector[i].coords[2]);
  }



  fprintf(arquivo,"CELLS %d %d\n", NumElem, NumElem*5);
  for(i=0; i<NumElem; i++)
  {
    fprintf(arquivo,"%d", Elem[i].type);
    for(j=0; j< Elem[i].type; j++)
      fprintf(arquivo,"  %d", Elem[i].inc[j]);
    fprintf(arquivo,"\n");
  }

  fprintf(arquivo,"CELL_TYPES %d\n", NumElem);
  for(i=0; i<NumElem; i++)
  {
    if (Elem[i].type == 3)
      fprintf(arquivo,"5\n");
    else if (Elem[i].type == 4)
      fprintf(arquivo,"9\n");
  }

  fclose(arquivo);

  return 1;
}


