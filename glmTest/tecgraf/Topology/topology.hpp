// topology.h: interHedFace for the topology class.
//
//////////////////////////////////////////////////////////////////////
#ifndef AM_SURFACE_TOPOLOGY
#define AM_SURFACE_TOPOLOGY

#include "deflib.h"

#include <vector>

using namespace std;

class TopNode_;
class TopEdge_;
class TopElem_;
class TopBoundary;
class mshsurfIndexMap;

#ifdef _WINDOWS
#pragma warning (disable : 4251 ) // retirado o aviso 4251 apenas para essa interface
// warning 4251:  class 'std::vector<_Ty>' needs to have dll-interface to be used by clients of class 'Spline'
// How to export an instantiation of a Standard Template Library (STL) class and a class that contains a data member that is an STL object
// http://support.microsoft.com/kb/q168958/
#endif


// cTopology
//////////////////////////////////////////////////////////////////////
class MSHSURF_API cTopology  
{
public:
           cTopology();
  virtual ~cTopology();
 
  void     Clear           (void);

  // main methods
  int      InsertMesh      (int n_node, int n_elem, double *coords, int *conn,
                            void (*func) (const char *));
  void     BoundBox        (double min[3], double max[3]);
  void     UpdateNormals   (void);

  int      writeNF         (char *filename);
  int      WriteVTK        (char *filename);
  int      writeMatlab     (char *filename, char *open_type);

  // number of nodes, edges, bound edges, and elements
  int  NumNodes          (void) { return NumNode; };
  int  NumEdges          (void) { return NumEdge; };
  int  NumBoundEdges     (void) { return NumEdgeBound; };
  int  NumElems          (void) { return NumElem; };


  // node information
  int     GetCoordNode      (int id, double coord[3]);
  void    SetCoordNode      (int id, double *new_coord);
  int     IsBoundaryNode    (int id);
  int     GetNormalNode     (int id, double normal[3]);
  void    SetNormalNode     (int id, double normal[3]);
  void    SetPtsInative     (int id);
  int     IsActNode         (int id);
  int     ReturnNodeId      (double x, double y, double z);
  void    SetNodeInfo       (int id, int flag);
  int     GetNodeInfo       (int id);

  // edge information
  int     GetEdge           (int id, int *ei, int *ej);
  int     GetEdge           (int ei, int ej);
  int     GetEdgeSize       (int i, double mid[3], double *size);
  int     GetBoundEdge      (int i, int *ei, int *ej);
  int     GetBEdgeSize      (int i, double mid[3], double *size);
  int     GetIsEdgeValid    (int i);
  double  GetSmalestEdge    (void);
  double  GetLongestEdge    (void);
  double  GetAverageEdge    (void);
  void    SetEdgeInfo       (int id, int flag);
  int     GetEdgeInfo       (int id);
  bool    IsBoundaryEdge    (int id);

  // element information
  int     GetElemNNodes     (int id); 
  int*    GetElemConn       (int id); 
  int     GetElemNorm       (int id, double normal[3]); 
  int     GetElemNormArea   (int id, double normal[3]); 
  int     GetElemIdEdge     (int id, int pos_adj); 
  int     GetIsElemValid    (int id);
  int     GetElemCenter     (int id, double center[3]);
  int     IsActElem         (int id);
  void    SetElemInfo       (int id, int flag);
  int     GetElemInfo       (int id);

  // loops
  int     GetNumLoops       (void);
  int     GetLoop           (int id, int *npts, int **idpts);

  // volume
  double  GetVolume ( );

  // adjacency information
  ///////////////////////////////////////////////////
  // node - nodes
  int     NumAdjNodeNode    (int id);
  void    AdjNodeNode       (int id, int *adj_nodes);
  int     OrientAdjNodeNode (int id, int *adj_nodes);
  // node - edges
  int     NumAdjEdgeToNode    (int id);
  int     AdjEdgeToNode       (int id, int pos);
  int     OppAdjEdgeToNode    (int id, int *ne, int **edges); // opposite adjacent edges to node
  int     OrientAdjEdgeToNode (int id, int *adj_edges);
  // node - elements
  int     AdjElemToNode     (int id, int *ne, int **elem);
  // edges - nodes
  int     AdjNodeToEdge     (int id, int *ei, int *ej);
  // edges - elements
  int     NumAdjElemToEdge  (int id);
  int     AdjElemToEdge     (int id, int pos);
  // elements - edges
  int     NumAdjEdgeToElem  (int id);
  int     AdjEdgeToElem     (int id, int pos);

private:

  // node information
  int         NumNode;
  TopNode_     *NodeVector;

  // element information
  int         NumElem;
  TopElem_     *Elem;

  // edge information
  mshsurfIndexMap      *EdgeList;
  vector<TopBoundary *> BoundList;
  TopEdge_              **EdgeVector;
  int                   NumEdge;
  TopEdge_              **BoundVector;
  int                   NumEdgeBound;

  // boundbox
  double min[3], max[3];

  // surfaces
  int     NumSurfaces;

  double     LocalTol;    // local tolerance 
  double     m_avgEdgeSize;
  double     m_minEdgeSize;
  double     m_maxEdgeSize;

  // Feedback function
  void (*feedbackFunction) (const char *);

  // methods
  int BuildEdgeList (void);

  int BuildAdjEdgesToNode (void);

  int BuildEdgesBoundary  (void);
  void RunAlongEdges (TopEdge_ *init_edge);

  void SetLocalToler (void);

  void IdentifySurfaces (void);
  void SetIdSurface     (int idSurf, TopElem_ *initElem);

  void SetValidFaces (void);

  void BuildOrientation (void);
  void OrientAdjElement (TopEdge_ *edge, int ciclo, TopElem_ **tail_elem);

  int  getAdjEdgeNextNodeEdge (int node, int edge);


};

#ifdef _WINDOWS
#pragma warning ( default : 4251 )  // voltando com os avisos 4251
// warning 4251:  class 'std::vector<_Ty>' needs to have dll-interface to be used by clients of class 'PolyLine'
// How to export an instantiation of a Standard Template Library (STL) class and a class that contains a data member that is an STL object
// http://support.microsoft.com/kb/q168958/
#endif


#endif
