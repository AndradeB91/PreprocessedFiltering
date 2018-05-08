//
// CArbMshRegion2D Class definition
//
// Description -
//   This class implements an object that will create a mesh
//   for an arbitrarily shaped region.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.39 $  $Date: 2002/09/06 14:57:57 $  $Author: wash $
//
// -------------------------------------

// #define DEBUG_DUMP

#include <stdio.h>
#include <math.h>
#include <assert.h>
/* #include <time.h> */

#include "ArbMshRegion2D.hpp"
#include "ArbQuadTree.hpp"
#include "ArbHeap.hpp"
#include "ArbMshTopo2D.hpp"
#include "ArbKDTree2D.hpp"
#include "ArbRectangle.hpp"
#include "ArbFeasableRegion.hpp"
#include "ArbMshSmooth2D.hpp"
#include "ArbSet.cpp"

#ifdef MEMDEBUG
#include "MemDbg.hpp"
#define new new(__FILE__,__LINE__)
#endif

#define DEFAULT_REFINE_BDRY_FACT 0.8  // fator para refinar a quadtree
//#define DEFAULT_NEAR_BDRY_FACT 0.6
#define DEFAULT_NEAR_BDRY_FACT 0.2
#define DEFAULT_MINIMUM_SHAPE_MEASURE 0.3

#define PI     3.141592654
#define TWO_PI 6.283185307

// ------------------------------------------------------------

static int CompareEdge(const ArbIntEdge &e1,const ArbIntEdge &e2)
{
    if (e1.node_id[0] > e2.node_id[0])
      return(1) ;
    else if (e1.node_id[0] < e2.node_id[0])
      return(-1) ;
    else
    {
      if (e1.node_id[1] > e2.node_id[1])
        return(1) ;
      else if (e1.node_id[1] < e2.node_id[1])
        return(-1) ;
    }
    return(0) ;
}


// %(CArbMshRegion2D::CArbMshRegion2D-constructor-|-ArbMshOrder-const|-ArbMshElemType-const|-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    CArbMshRegion2D - constructor method
**
**      CArbMshRegion2D(
**              const ArbMshOrder    iorder = LINEAR,
**              const ArbMshElemType itype = TRIANGLE,
**              const int            istart_id = 0,
**              const int            imat_id = 0)
**
**        iorder    - (in)  LINEAR or QUADRATIC
**        itype     - (in)  TRIANGLE or QUADRILATERAL
**        istart_id - (in)  starting generated node id
**        imat_id   - (in)  region's material id
**
**      Description: This is a constructor method for a CArbMshRegion2D
**          object.
**
**
** -- */

CArbMshRegion2D::CArbMshRegion2D (
const ArbMshOrder iorder,
const ArbMshElemType itype,
const int istart_id,
const int imat_id
) : Order(iorder),Type(itype),StartId(istart_id),MatId(imat_id)
{
    NodeGeneration       = true ;
    BoundaryChecks       = false ;
    DoSmoothNodes        = true ;
    DoCleanup            = true ;
    DoQuadAngleChecks    = false ;
    WinslowSmoothing     = true ;
    QuadTreeDebug        = false ;
    StartId              = 0 ;
    MaxId                = 0 ;
    StartElemId          = 0 ;
    RefineBoundaryFactor = DEFAULT_REFINE_BDRY_FACT ;
    NearBoundaryFactor   = DEFAULT_NEAR_BDRY_FACT ;
    MinimumShapeMeasure  = DEFAULT_MINIMUM_SHAPE_MEASURE ;
    pnode_table          = new CArbHashTable<int,ArbIntNode> ;
    pedge_table          = new CArbSet<ArbIntEdge>(CompareEdge) ;
    pelem_table          = new CArbHashTable<int,ArbMshElement2D> ;
    MateTable            = 0 ;
#ifdef DOING_QUADS
    CoordCache           = new CArbArray<CArbCoord2D>() ;
    BdryVtxCache         = new CArbArray<int>() ;
#endif
    DebugDisplayFlags = 0 ;

    nnode = nelem = 0;
    conn = NULL;
    coords = NULL;


}




// %(CArbMshRegion2D::CArbMshRegion2D-destructor-|~)
/* ++ ----------------------------------------------------------
**
**    CArbMshRegion2D - destructor
**
**      ~CArbMshRegion2D()
**
**      Description: This is a destructor for a CArbMshRegion2D object.
**
**
** -- */

CArbMshRegion2D::~CArbMshRegion2D()
{
    delete pnode_table ;
    delete pedge_table ;
    delete pelem_table ;
    if (MateTable != 0) delete MateTable ;
#ifdef DOING_QUADS
    delete CoordCache ;
    delete BdryVtxCache ;
#endif

    nnode = nelem = 0;
    if (conn != NULL)
      delete []conn;
    if (coords != NULL)
      delete []coords;
}




// %(CArbMshRegion2D::SetOrder-void-|-ArbMshOrder-const|)
/* ++ ----------------------------------------------------------
**
**    SetOrder - sets the polynomial order
**
**      void SetOrder(const ArbMshOrder iorder)
**
**        iorder - (in)  LINEAR or QUADRATIC
**
**      Description: This method sets the polynomial order for
**          generated elements.
**
**
** -- */

void CArbMshRegion2D::SetOrder(const ArbMshOrder iorder)
{
    Order = iorder ;
}




// %(CArbMshRegion2D::SetElemType-void-|-ArbMshElemType-const|)
/* ++ ----------------------------------------------------------
**
**    SetElemType - sets the element type
**
**      void SetElemType(const ArbMshElemType itype)
**
**        itype - (in)  TRIANGLE or QUADRILATERA
**
**      Description: This method sets the element type (shape) of the
**          generated elements.
**
**
** -- */

void CArbMshRegion2D::SetElemType(const ArbMshElemType itype)
{
    Type = itype ;
}




// %(CArbMshRegion2D::SetStartNodeID-void-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    SetStartNodeID - sets the first generated node id
**
**      void SetStartNodeID(const int is_id)
**
**        is_id - (in)  first node id
**
**      Description: This method sets the first id to be assigned to
**          generated nodes. Subsequent id's will increase from this.
**
**
** -- */

void CArbMshRegion2D::SetStartNodeID(const int is_id)
{
    StartId = is_id ;
}




// %(CArbMshRegion2D::SetStartElemID-void-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    SetStartElemID - sets the first generated element id
**
**      void SetStartElemID(const int is_id)
**
**        is_id - (in)  first element id
**
**      Description: This method sets the first id to be assigned to
**          generated elements. Subsequent id's will increase from
**          this.
**
**
** -- */

void CArbMshRegion2D::SetStartElemID(const int is_id)
{
    StartElemId = is_id ;
}




// %(CArbMshRegion2D::SetMaterialID-void-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    SetMaterialID - sets the assigned material id
**
**      void SetMaterialID(const int mat_id)
**
**        mat_id - (in)  material id
**
**      Description: This method sets the material id to be assigned to
**          generated elements.
**
**
** -- */

void CArbMshRegion2D::SetMaterialID(const int mat_id)
{
    MatId = mat_id ;
}




// %(CArbMshRegion2D::SetNodeGeneration-void-|-bool-|)
/* ++ ----------------------------------------------------------
**
**    SetNodeGeneration - set the node generation option
**
**      void SetNodeGeneration(bool flag = true)
**
**        flag - (in)  on or off flag
**
**      Description: This method turns on or off the generation of
**          internal nodes during meshing.
**
**
** -- */

void CArbMshRegion2D::SetNodeGeneration(bool flag)
{
    NodeGeneration = flag ;
}




// %(CArbMshRegion2D::SetNodeSmoothing-void-|-bool-|)
/* ++ ----------------------------------------------------------
**
**    SetNodeSmoothing - set the nodal smoothing option
**
**      void SetNodeSmoothing(bool flag = true)
**
**        flag - (in)  on or off flag
**
**      Description: This method turns on or off the option to smooth
**          internal node to improve element quality.
**
**
** -- */

void CArbMshRegion2D::SetNodeSmoothing(bool flag)
{
    DoSmoothNodes = flag ;
}




// %(CArbMshRegion2D::SetBoundaryValidityChecks-void-|-bool-|)
/* ++ ----------------------------------------------------------
**
**    SetBoundaryValidityChecks - set the boundary check option
**
**      void SetBoundaryValidityChecks(bool flag = true)
**
**        flag - (in)  on or off flag
**
**      Description: This method turns on or off the option to check
**          the validity of the input boundary before elements are
**          generated.
**
**
** -- */

void CArbMshRegion2D::SetBoundaryValidityChecks(bool flag)
{
    BoundaryChecks = flag ;
}




// %(CArbMshRegion2D::SetRefineBoundaryFactor-void-|-double-|)
/* ++ ----------------------------------------------------------
**
**    SetRefineBoundaryFactor - set the refine boundary factor
**
**      void SetRefineBoundaryFactor(double factor)
**
**        factor - (in)  boundary factor
**
**      Description: This method sets a factor used when generating
**          internal nodes. Higher values will generate denser meshes.
**
**
** -- */

void CArbMshRegion2D::SetRefineBoundaryFactor(double factor)
{
    RefineBoundaryFactor = factor ;
}




// %(CArbMshRegion2D::SetNearBoundaryFactor-void-|-double-|)
/* ++ ----------------------------------------------------------
**
**    SetNearBoundaryFactor - set the near boundary factor
**
**      void SetNearBoundaryFactor(double factor)
**
**        factor - (in)  boundary factor
**
**      Description: This method sets a factor used when generating
**          internal nodes. The factor is used to remove generated
**          nodes that lie too near a boundary.
**
**
** -- */

void CArbMshRegion2D::SetNearBoundaryFactor(double factor)
{
    NearBoundaryFactor = factor ;
}


/* ++ ----------------------------------------------------------
**
**    SetQuadAngleChecks - turns on checks on max and min quad angles
**
**      void SetQuadElemAngleChecks(double min_angle,
**                                  double max_angle,
**                                  AngleCheckLocation location)
**
**        min_angle - (in) minimum acceptable angle (radians)
**        max_angle - (in) maximum acceptable angle (radians)
**        location  - (in) specify if angles should be checked
**                         at the nodes or at the integration points.
**
**      Description: This turns on quadrilateral angle checking and
**          sets the check parameters.  If a quadrilateral element
**          is found to have angles that are out of bounds it is
**          replaced with two triangular elements.
**
** -- */

void CArbMshRegion2D::SetQuadElemAngleChecks(
                                double min_angle,
                                double max_angle,
                                AngleCheckLocation location)
{
    DoQuadAngleChecks = true ;
    MinQuadAngle = min_angle ;
    MaxQuadAngle = max_angle ;
    QuadAngleLocation = location ;
}


// %(CArbMshRegion2D::GetOrder-ArbMshOrder-|^const)
/* ++ ----------------------------------------------------------
**
**    GetOrder - get the polynomial order
**
**      ArbMshOrder GetOrder() const
**
**      Description: This method returns the polynomial order of
**          generated elements.
**
**      Return Value: LINEAR or QUADRATIC
**
**
** -- */

ArbMshOrder CArbMshRegion2D::GetOrder() const
{
    return(Order) ;
}




// %(CArbMshRegion2D::GetElemType-ArbMshElemType-|^const)
/* ++ ----------------------------------------------------------
**
**    GetElemType - get the element type
**
**      ArbMshElemType GetElemType() const
**
**      Description: This method returns the type (shape) of generated
**          elements
**
**      Return Value: TRIANGLE or QUADRILATERAL
**
**
** -- */

ArbMshElemType CArbMshRegion2D::GetElemType() const
{
    return(Type) ;
}




// %(CArbMshRegion2D::GetStartNodeID-int-|^const)
/* ++ ----------------------------------------------------------
**
**    GetStartNodeID - get the start node id
**
**      int GetStartNodeID() const
**
**      Description: This method returns the first id that will be
**          assigned to generated internal nodes.
**
**      Return Value: lowest id of generated nodes
**
**
** -- */

int CArbMshRegion2D::GetStartNodeID() const
{
    return(StartId) ;
}




// %(CArbMshRegion2D::SetQuadTreeDebug-void-|-bool-|)
/* ++ ----------------------------------------------------------
**
**    SetQuadTreeDebug - set a debug flag for the quadtree
**
**      void SetQuadTreeDebug(bool flag = true)
**
**        flag - (in)  on or off
**
**      Description: This method sets a quadtree debugging flag. If set
**          the debugging data printed.
**
**
** -- */

void CArbMshRegion2D::SetQuadTreeDebug(bool flag)
{
    QuadTreeDebug = flag ;
}




// %(CArbMshRegion2D::SetWinslowSmoothing-void-|-bool-|)
/* ++ ----------------------------------------------------------
**
**    SetWinslowSmoothing - sets the Winslow smoothing option
**
**      void SetWinslowSmoothing(bool flag = true)
**
**        flag - (in)  on or off
**
**      Description: This method turns on or off the option to us a
**          Winslow smoothing algorithm. If turned off a Laplace
**          algorithm is used.
**
**
** -- */

void CArbMshRegion2D::SetWinslowSmoothing(bool flag)
{
    WinslowSmoothing = flag ;
}




// %(CArbMshRegion2D::SetTopoCleanup-void-|-bool-|)
/* ++ ----------------------------------------------------------
**
**    SetTopoCleanup - set the topological cleanup option
**
**      void SetTopoCleanup(bool flag = true)
**
**        flag - (in)  on or off flag
**
**      Description: This method turns on or off the option to perform
**          "topological cleanup" to improve element quality.
**
**
** -- */

void CArbMshRegion2D::SetTopoCleanup(bool flag)
{
    DoCleanup = flag ;
}




// %(CArbMshRegion2D::AddNode-void-|-ArbMshNode-const|&-ArbMshNodeType-const|-ArbMshNodeMotion-const|)
/* ++ ----------------------------------------------------------
**
**    AddNode - add a node to the region
**
**      int AddNode(
**              const ArbMshNode       &inode,
**              const ArbMshNodeType   itype = BOUNDARY,
**              const ArbMshNodeMotion imotion = ARB_FLOATING)
**
**        inode   - (in)  node description
**        itype   - (in)  BOUNDARY or INTERIOR
**        imotion - (in)  ARB_FIXED or ARB_FLOATING
**
**      Description: This method adds a node to the boundary of of the
**          current region.
**
**      Returns:
**          ARB_DUPLICATE_NODE_ID - duplicate node id
**          ARB_NORMAL_STATUS     - normal return status
**
** -- */

int CArbMshRegion2D::AddNode(const ArbMshNode &inode,
                         const ArbMshNodeType itype,
                         const ArbMshNodeMotion imotion)
{
    ArbIntNode int_node ;
    int lnode = inode.id ;

    int_node.id = inode.id ;
    for (int i=0 ; i<2 ; ++i) int_node.coord[i] = inode.coord[i] ;
    int_node.type = itype ;
    int_node.motion = (itype == BOUNDARY) ? ARB_FIXED : imotion ;

    // If this is an interior node, set it as a corner node.
    // The boundary nodes are set to corner nodes in add egdge
    // if they are part of a boundary edge.

    int_node.corner = (itype == BOUNDARY) ? false : true ;

    if (pnode_table->Store(lnode,int_node) == false) {
        return (ARB_DUPLICATE_NODE_ID) ;
    }

    if (inode.id > MaxId) MaxId = inode.id ;

    return (ARB_NORMAL_STATUS) ;
}




// %(CArbMshRegion2D::AddNodeList-void-|-int-const|-ArbMshNode-const|*const-ArbMshNodeType-const|-ArbMshNodeMotion-const|)
/* ++ ----------------------------------------------------------
**
**    AddNodeList - add a list of nodes to the region.
**
**      int AddNodeList(
**              const int              num_nodes,
**              const ArbMshNode       *constnode_list,
**              const ArbMshNodeType   itype = BOUNDARY,
**              const ArbMshNodeMotion imotion = ARB_FLOATING)
**
**        num_nodes - (in)  number of nodes in list
**        node_list - (in)  list of nodes
**        itype     - (in)  BOUNDARY or INTERIOR
**        imotion   - (in)  ARB_FIXED or ARB_FLOATING
**
**      Description: This method adds a list of nodes to the boundary
**          of the current region.
**
**      Returns:
**          ARB_DUPLICATE_NODE_ID - duplicate node id
**          ARB_NORMAL_STATUS     - normal return status
**
**
** -- */

int CArbMshRegion2D::AddNodeList(const int num_nodes,
                         const ArbMshNode *const node_list,
                         const ArbMshNodeType itype,
                         const ArbMshNodeMotion imotion)
{
    for (int i=0 ; i<num_nodes ; ++i) {
        int status = AddNode(node_list[i],itype,imotion) ;
        if (status != ARB_NORMAL_STATUS) return(status) ;
    }
    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshRegion2D::AddEdge-void-|-ArbMshEdge-const|&)
/* ++ ----------------------------------------------------------
**
**    AddEdge - add an edge to the region
**
**      int AddEdge(const ArbMshEdge &edge)
**
**        edge - (in)  edge description
**
**      Description: This method adds an edge to the boundary of the
**          current region.
**
**      Returns:
**          ARB_BAD_NODE_ID    - invalid edge end node id
**          ARB_DUPLICATE_EDGE - duplicate edge
**          ARB_NORMAL_STATUS  - normal return status
**
**
** -- */

int CArbMshRegion2D::AddEdge(const ArbMshEdge &edge)
{
    ArbIntNode *pTmpNode ;
    int tmp_node ;

    // first check to see if the nodes are valid

    tmp_node = edge.node_id[0] ;
    if ((pTmpNode = pnode_table->Fetch(tmp_node)) == 0) {
        return(ARB_BAD_NODE_ID) ;
    } else {
        pTmpNode->corner = true ;
    }
    tmp_node = edge.node_id[1] ;
    if ((pTmpNode = pnode_table->Fetch(tmp_node)) == 0) {
        return(ARB_BAD_NODE_ID) ;
    } else {
        pTmpNode->corner = true ;
    }
    if (Order == QUADRATIC) {
        tmp_node = edge.node_id[2] ;
        if ((pTmpNode = pnode_table->Fetch(tmp_node)) == 0) {
            return(ARB_BAD_NODE_ID) ;
        }
    }

    // now build a internal edge structure, create a key
    // by adding the first two node numbers then try to
    // add them to the set of edges.

    ArbIntEdge tmp_edge ;
    for (int i=0 ; i<3 ; ++i) tmp_edge.node_id[i] = edge.node_id[i] ;
    tmp_edge.type = BOUNDARY ;

    if (pedge_table->Insert(tmp_edge) == false) {
        return(ARB_DUPLICATE_EDGE) ;
    }
    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshRegion2D::AddEdgeList-void-|-int-const|-ArbMshEdge-const|*const)
/* ++ ----------------------------------------------------------
**
**    AddEdgeList - add a list of edges to the region.
**
**      void AddEdgeList(
**              const int        num_edges,
**              const ArbMshEdge *constedge_list)
**
**        num_edges - (in)  number of edges in the list
**        edge_list - (in)  list of edge descriptions
**
**      Description: This method adds a list of edges to the boundary
**          of the current region.
**
**      Returns:
**          ARB_BAD_NODE_ID    - invalid edge end node id
**          ARB_DUPLICATE_EDGE - duplicate edge
**          ARB_NORMAL_STATUS  - normal return status
**
** -- */

int CArbMshRegion2D::AddEdgeList(const int num_edges,
                         const ArbMshEdge *const edge_list)
{
    for (int i=0 ; i<num_edges ; ++i) {
        int status = AddEdge(edge_list[i]) ;
        if (status != ARB_NORMAL_STATUS) return(status) ;
    }
    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshRegion2D::NumBoundaryNodes-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumBoundaryNodes - number of boundary nodes
**
**      int NumBoundaryNodes() const
**
**      Description: This method returns the number of nodes on the
**          regions's boundary.
**
**      Return Value: number of boundary nodes
**
**
** -- */

int CArbMshRegion2D::NumBoundaryNodes() const
{
    CArbHashTableIterator<int,ArbIntNode> iter(pnode_table) ;
    int num = 0 ;

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.Entry()->type == BOUNDARY) ++num ;
    }
    return(num) ;
}




// %(CArbMshRegion2D::NumInternalNodes-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumInternalNodes - number of internal nodes
**
**      int NumInternalNodes() const
**
**      Description: This method returns the number of nodes in the
**          regions's Interior.
**
**      Return Value: number of internal nodes
**
**
** -- */

int CArbMshRegion2D::NumInternalNodes() const
{
    CArbHashTableIterator<int,ArbIntNode> iter(pnode_table) ;
    int num = 0 ;

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.Entry()->type == INTERIOR) ++num ;
    }
    return(num) ;
}




// %(CArbMshRegion2D::NumNodes-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumNodes - number of nodes
**
**      int NumNodes() const
**
**      Description: This method returns the total number of nodes in
**          the region.
**
**      Return Value: total number of nodes
**
**
** -- */

int CArbMshRegion2D::NumNodes() const
{
    return(pnode_table->NumEntries()) ;
}




// %(CArbMshRegion2D::GetBoundaryNodes-ArbMshNode-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetBoundaryNodes - get the boundary nodes
**
**      ArbMshNode *GetBoundaryNodes() const
**
**      Description: This method returns a list of the nodes on the
**          boundary of the region. Ownership of this memory passes to
**          the client, which must eventually call delete [].
**
**      Return Value: a list of the nodes on the boundary
**
**
** -- */

ArbMshNode *CArbMshRegion2D::GetBoundaryNodes() const
{
    CArbHashTableIterator<int,ArbIntNode> iter(pnode_table) ;
    ArbMshNode *bound_entries ;
    int num_bound = NumBoundaryNodes() ;
    int cur = 0 ;

    bound_entries = new ArbMshNode[num_bound] ;

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.Entry()->type == BOUNDARY) {
            bound_entries[cur].id = iter.Entry()->id ;
            for (int j=0 ; j<2 ; ++j)
                bound_entries[cur].coord[j] = iter.Entry()->coord[j] ;
            bound_entries[cur].coord[2] = 0.0 ;
            ++cur ;
        }
    }
    return(bound_entries) ;
}




// %(CArbMshRegion2D::GetInternalNodes-ArbMshNode-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetInternalNodes - get the internal nodes
**
**      ArbMshNode *GetInternalNodes() const
**
**      Description: This method returns a list of the nodes in the
**          interior of the region. Ownership of this memory passes to
**          the client, which must eventually call delete [].
**
**      Return Value: a list of the nodes in the interior
**
**
** -- */

ArbMshNode *CArbMshRegion2D::GetInternalNodes() const
{
    CArbHashTableIterator<int,ArbIntNode> iter(pnode_table) ;
    ArbMshNode *intern_entries ;
    int num_intern = NumInternalNodes() ;
    int cur = 0 ;

    intern_entries = new ArbMshNode[num_intern] ;

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.Entry()->type == INTERIOR) {
            intern_entries[cur].id = iter.Entry()->id ;
            for (int j=0 ; j<2 ; ++j)
                intern_entries[cur].coord[j] = iter.Entry()->coord[j] ;
            intern_entries[cur].coord[2] = 0.0 ;
            ++cur ;
        }
    }
    return(intern_entries) ;
}




// %(CArbMshRegion2D::GetNodes-ArbMshNode-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetNodes - get all the nodes
**
**      ArbMshNode *GetNodes() const
**
**      Description: This method returns a list of all the nodes in the
**          region. Ownership of this memory passes to the client,
**          which must eventually call delete [].
**
**      Return Value: a list of all the nodes in the region
**
**
** -- */

ArbMshNode *CArbMshRegion2D::GetNodes() const
{
    CArbHashTableIterator<int,ArbIntNode> iter(pnode_table) ;
    ArbMshNode *local_entries ;
    int num_total = pnode_table->NumEntries() ;

    local_entries = new ArbMshNode[num_total] ;

    for (int i=0 ; i<num_total ; ++i,++iter) {
        local_entries[i].id = iter.Entry()->id ;
        local_entries[i].coord[0] = iter.Entry()->coord.x() ;
        local_entries[i].coord[1] = iter.Entry()->coord.y() ;
        local_entries[i].coord[2] = 0.0 ;
    }
    return(local_entries) ;
}




// %(CArbMshRegion2D::NumBoundaryEdges-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumBoundaryEdges - number of boundary edges
**
**      int NumBoundaryEdges() const
**
**      Description: This method retuns the number of boundary edges in
**          the region.
**
**      Return Value: number of boundary edges
**
**
** -- */

int CArbMshRegion2D::NumBoundaryEdges() const
{
    ArbIntEdge **entries ;
    int nume = pedge_table->NumEntries() ;
    int num = 0 ;

    entries = pedge_table->GetKeyList() ;
    for (int i=0 ; i<nume ; ++i) {
        if (entries[i]->type == BOUNDARY) ++num ;
    }
    delete [] entries ;
    return(num) ;
}




// %(CArbMshRegion2D::GetBoundaryEdges-ArbMshEdge-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetBoundaryEdges - get the boundary edges
**
**      ArbMshEdge *GetBoundaryEdges() const
**
**      Description: This method returns a list of all the edges the
**          make up the boundary of the region. Ownership of this
**          memory passes to the client, which must eventually call
**          delete [].
**
**      Return Value: a list of the boundary edges
**
**
** -- */

ArbMshEdge *CArbMshRegion2D::GetBoundaryEdges() const
{
    ArbIntEdge **full_entries ;
    ArbMshEdge *bound_entries ;
    int num_bound = NumBoundaryNodes() ;
    int num_total = pedge_table->NumEntries() ;
    int cur = 0 ;

    full_entries = pedge_table->GetKeyList() ;
    bound_entries = new ArbMshEdge[num_bound] ;

    for (int i=0 ; i<num_total ; ++i) {
        if (full_entries[i]->type == BOUNDARY) {
            for (int j=0 ; j<3 ; ++j)
                bound_entries[cur].node_id[j] =
                    full_entries[i]->node_id[j] ;
            ++cur ;
        }
    }
    delete [] full_entries ;
    return(bound_entries) ;
}




// %(CArbMshRegion2D::GetElements-ArbMshElement2D-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetElements - get the elements
**
**      ArbMshElement2D *GetElements() const
**
**      Description: This method returns a list of all the elements in
**          the region. Ownership of this memory passes to the client,
**          which must eventually call delete [].
**
**      Return Value: a list of the elements
**
**
** -- */

ArbMshElement2D *CArbMshRegion2D::GetElements() const
{
    ArbMshElement2D *local_entries ;
    int num_total = pelem_table->NumEntries() ;

    CArbHashTableIterator<int,ArbMshElement2D> eiter(pelem_table) ;
    local_entries = new ArbMshElement2D[num_total] ;

    int i=0 ;
    for (eiter.First() ; eiter.More() ; ++eiter) {
        local_entries[i] = *(eiter.Entry()) ;
        ++i ;
    }
    return(local_entries) ;
}




// %(CArbMshRegion2D::NumElements-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumElements - number of elements
**
**      int NumElements() const
**
**      Description: This method returns the number of elements in the
**          region.
**
**      Return Value: number of elements
**
**
** -- */

int CArbMshRegion2D::NumElements() const
{
    return(pelem_table->NumEntries()) ;
}

// -----------------------------------------------------------

/* ------------------------------------------------------------
    ArbGenerateIntNodes - callback function to generate nodes
                          at the center of each quadtree cell.

    Input:
        void *p_data - callback data which is realy a pointer
                       to the calling CArbMshRegion2D object.
        double origin_x - x coordinate of the center of this cell.
        double origin_y - y coordinate of the center of this cell.
        double half - half size of this cell.
        bool is_root - true if this is the root node.
        bool is_leaf - true if this is a leaf node.
*/

void ArbGenerateIntNodes(void *p_data,
                double origin_x,double origin_y,double half,
                bool is_root,bool is_leaf)
{
    CArbMshRegion2D *p_region = (CArbMshRegion2D *)p_data ;

    if (is_leaf) {

        // if this is a leaf node check to see if a node
        // generated at this point is inside the region
        // and that it is not too close to any boundary edge

        int num_edges = p_region->pedge_table->NumEntries() ;
        ArbIntEdge **edges = p_region->pedge_table->GetKeyList() ;
        bool insert = true ;
        int n_cross = 0 ;

        CArbCoord2D b(origin_x,origin_y) ;

        for (int i=0 ; i<num_edges ; ++i) {
            ArbIntNode *p_nd0, *p_nd1 ;
            double len, rr, dist_sqr ;

            p_nd0 = p_region->pnode_table->Fetch(edges[i]->node_id[0]) ;
            p_nd1 = p_region->pnode_table->Fetch(edges[i]->node_id[1]) ;

            // find the minimum distance from the point to the
            // current edge

            dist_sqr = p_region->DistSqr(b,p_nd0->coord,p_nd1->coord,&len) ;

            // find the radius from the edge within which nodes
            // will be deleted.  We square this so that we can
            // do comparions on the square of the length and not
            // have to take square roots.

            rr = len*p_region->NearBoundaryFactor ;
            rr = rr * rr ;
            if (dist_sqr < rr) {
                insert = false ;
                break ;
            }

            // check to see if a scan line from minus infinity
            // to this point crosses the edge

            if (p_region->ScanCross(b,p_nd0->coord,p_nd1->coord))
                ++n_cross ;
        }

        if (insert && ((n_cross % 2) == 1)) {
            (void)p_region->NewNode(origin_x,origin_y,INTERIOR,
                                    ARB_FLOATING,true) ;
        }
        delete [] edges ;
    }


    if (p_region->QuadTreeDebug) {
        double x0,y0,x1,y1 ;

        if (is_root) {
            x0 = origin_x - half ; y0 = origin_y - half ;
            x1 = origin_x + half ; y1 = origin_y - half ;
            printf("l %g %g %g %g\n",x0,y0,x1,y1) ;

            x0 = x1 ; y0 = y1 ;
            x1 = origin_x + half ; y1 = origin_y + half ;
            printf("l %g %g %g %g\n",x0,y0,x1,y1) ;

            x0 = x1 ; y0 = y1 ;
            x1 = origin_x - half ; y1 = origin_y + half ;
            printf("l %g %g %g %g\n",x0,y0,x1,y1) ;

            x0 = x1 ; y0 = y1 ;
            x1 = origin_x - half ; y1 = origin_y - half ;
            printf("l %g %g %g %g\n",x0,y0,x1,y1) ;
        }

        if (!is_leaf) {
            x0 = origin_x - half ; y0 = origin_y ;
            x1 = origin_x + half ; y1 = origin_y ;
            printf("l %g %g %g %g\n",x0,y0,x1,y1) ;

            x0 = origin_x ; y0 = origin_y + half ;
            x1 = origin_x ; y1 = origin_y - half ;
            printf("l %g %g %g %g\n",x0,y0,x1,y1) ;
        }
    }
}


// %(CArbMshRegion2D::DebugDumpRegion-void-|)
/* ++ ----------------------------------------------------------
**
**    DebugDumpRegion - debug routine to save region information
**
**      void DebugDumpRegion()
**
**      Description: This method provides debugging support. It stores
**          information about the region in the file "debug.rgn".
**
**
** -- */

void CArbMshRegion2D::DebugDumpRegion()
{
    FILE *fd = fopen("debug.rgn","w") ;
    int i ;
    if (fd != 0) {
        ArbIntNode **full_entries ;
        int num_total = pnode_table->NumEntries() ;
        fprintf(fd,"%d\n",num_total) ;
        full_entries = pnode_table->GetEntryList() ;
        for (i=0 ; i<num_total ; ++i) {
            fprintf(fd,"%d %18.10g %18.10g\n",full_entries[i]->id,
                                    full_entries[i]->coord[0],
                                    full_entries[i]->coord[1]) ;
        }
        delete full_entries ;

        ArbIntEdge **edge_entries ;
        int num_bound = NumBoundaryNodes() ;
        num_total = pedge_table->NumEntries() ;
        edge_entries = pedge_table->GetKeyList() ;
        fprintf(fd,"%d\n",num_bound) ;
        for (i=0 ; i<num_total ; ++i) {
            if (edge_entries[i]->type == BOUNDARY) {
                if (Order == QUADRATIC) {
                    fprintf(fd,"%d %d %d\n",
                        edge_entries[i]->node_id[0],
                        edge_entries[i]->node_id[1],
                        edge_entries[i]->node_id[2]) ;
                } else {
                    fprintf(fd,"%d %d\n",
                        edge_entries[i]->node_id[0],
                        edge_entries[i]->node_id[1]) ;
                }
            }
        }
        delete edge_entries ;
        fclose(fd) ;
    }
}


// %(CArbMshRegion2D::DebugDumpMesh-void-|)
/* ++ ----------------------------------------------------------
**
**    DebugDumpMesh - debug routine to save mesh information
**
**      void DebugDumpMesh()
**
**      Description: This method provides debugging support. It stores
**          information about the mesh in the files "debug.cds" and
**          "debug.msh".
**
**
** -- */

#ifdef DEBUG_DUMP

void CArbMshRegion2D::DebugDumpMesh()
{
  FILE *fd = fopen("d:\\temp\\debug.cds","w") ;
    if (fd != 0) {
        ArbIntNode **nodes ;
        int i, num_node = pnode_table->NumEntries() ;

        nodes = pnode_table->GetEntryList() ;
        fprintf(fd,"%d\n",num_node) ;
        for (i=0 ; i<num_node ; ++i) {
            fprintf(fd,"%d %g %g\n",nodes[i]->id,
                                    nodes[i]->coord[0],
                                    nodes[i]->coord[1]) ;
        }
        delete [] nodes ;
        fclose(fd) ;
    }

    fd = fopen("d:\\temp\\debug.msh","w") ;
    if (fd != 0) {
        ArbMshElement2D **elements ;
        int i, j, num_elem = pelem_table->NumEntries() ;

        elements = pelem_table->GetEntryList() ;
        fprintf(fd,"%d\n",num_elem) ;
        for (i=0 ; i<num_elem ; ++i) {
            fprintf(fd,"%d %d",elements[i]->elem_id,
                               elements[i]->num_nodes) ;
            for (j=0 ; j<elements[i]->num_nodes ; ++j)
                fprintf(fd," %d",elements[i]->nodes[j]) ;
            fprintf(fd,"\n") ;

        }
        delete [] elements ;
        fclose(fd) ;
    }
}

#endif





// %(CArbMshRegion2D::GenerateMesh-void-|)
/* ++ ----------------------------------------------------------
**
**    GenerateMesh - generate a mesh
**
**      void GenerateMesh()
**
**      Description: This method generates a mesh for the region.
**
**      Returns:
**          ARB_ILLEGAL_BDRY   - the boundary specified is illegal
**          ARB_NORMAL_STATUS  - normal return status
**
** -- */
int CArbMshRegion2D::GenerateMesh ( )
{
/*
    static long double cpu_time;
*/

#ifdef DEBUG_DUMP
    DebugDumpRegion ( ) ;
#endif

    if (DebugDisplayFlags & RmshBoundary)
        DisplayBoundary("Remesh Boundary") ;

/*
    cpu_time = clock( );
*/

    // if necessary, check the validity of the boundary
    ////////////////////////////////////////////////////
    if (BoundaryChecks)
    {
        int status = CheckBoundary() ;
        if (status != ARB_NORMAL_STATUS) return(status) ;
    }

/*
    cpu_time = (clock( ) - cpu_time)/CLOCKS_PER_SEC;
    printf("\n");
    printf("\t\tCPU MsqQuad - Boundary Checks............. %0.3f (s)\n", (double)cpu_time);
*/

    // generate internal nodes if required
    ////////////////////////////////////////////////////
    if (StartId == 0) StartId = MaxId+1 ;
    StartIdSave = StartId ;
    StartElemIdSave = StartElemId ;
    if (NodeGeneration) GenerateNodes ( ) ;


/*
    cpu_time = clock( );
*/

    // do the boundary contraction
    ////////////////////////////////////////////////////
    BoundaryContraction() ;

/*
    cpu_time = (clock( ) - cpu_time)/CLOCKS_PER_SEC;
    printf("\n");
    printf("\t\tCPU MsqQuad - Boundary Contraction............. %0.3f (s)\n", (double)cpu_time);
*/


#ifdef DEBUG_DUMP
    DebugDumpMesh() ;
#endif

    if (DebugDisplayFlags & RmshTriBeforeSmooth)
        DisplayMesh("Trimesh Before Smoothing") ;

/*
    cpu_time = clock( );
*/

    // smooth the internal nodes
    ////////////////////////////////////////////////////
    if (DoSmoothNodes)
    {
        SmoothNodes() ;

        // Sometimes it seems that winslow smoothing can do bad
        // things to triangular meshes with poor element size
        // distribution.  Therefore, we check for bad elements
        // and if found fall back to laplace smoothing
        if (WinslowSmoothing && !CheckValidTriangles())
        {
            CArbMshSmooth2D smooth(pnode_table,pelem_table) ;
            smooth.SmoothNodesLaplace() ;
        }
    }

/*
    cpu_time = (clock( ) - cpu_time)/CLOCKS_PER_SEC;
    printf("\n");
    printf("\t\tCPU MsqQuad - Smooth Nodes............. %0.3f (s)\n", (double)cpu_time);
*/

    if (DebugDisplayFlags & RmshTriAfterSmooth)
        DisplayMesh("Trimesh After Smoothing") ;

#ifdef DOING_QUADS
/*
    cpu_time = clock( );
*/

    // if quadrilateral elements are required, generate them
    ////////////////////////////////////////////////////
    if (Type == QUADRILATERAL)
    {
        if (!GenerateQuads())
         return (ARB_ILLEGAL_BOUNDARY);
        RenumberNodes() ;
    }

/*
    cpu_time = (clock( ) - cpu_time)/CLOCKS_PER_SEC;
    printf("\n");
    printf("\t\tCPU MsqQuad - Quadrilateral Generate............. %0.3f (s)\n", (double)cpu_time);
*/
#endif

    return(ARB_NORMAL_STATUS) ;
}

// CheckValidTriangles
//////////////////////////////////////////////////////////
bool CArbMshRegion2D::CheckValidTriangles ( )
{
    CArbHashTableIterator<int,ArbMshElement2D> iter(pelem_table) ;
    for (iter.First() ; iter.More() ; ++iter) {
        ArbIntNode *nd0, *nd1, *nd2 ;
        nd0 = pnode_table->Fetch(iter.Entry()->nodes[0]) ;
        nd1 = pnode_table->Fetch(iter.Entry()->nodes[1]) ;
        nd2 = pnode_table->Fetch(iter.Entry()->nodes[2]) ;
        double sm = TriShapeMeasure(nd0->coord,nd1->coord,nd2->coord) ;
        if (sm <= 0.0) return(false) ;
    }
    return(true) ;
}


// %(CArbMshRegion2D::CheckBoundary-void-|)
/* ++ ----------------------------------------------------------
**
**    CheckBoundary - perform a check of the boundary
**
**      void CheckBoundary()
**
**      Description: This method performs a check of the region's
**          boundary to see if it is valid.
**
**      Returns:
**          ARB_ILLEGAL_BDRY   - the boundary specified is illegal
**          ARB_NORMAL_STATUS  - normal return status
**
** -- */
int CArbMshRegion2D::CheckBoundary()
{
    int i ;
//    int j ;

    // get a local list of all the edges

    int num_edges = pedge_table->NumEntries() ;
    ArbIntEdge **edges = pedge_table->GetKeyList() ;

    // check to make sure that we have at least 3 edges and nodes

    if ((num_edges < 3) || (pnode_table->NumEntries() < 3))
    {
        printf ("Illegal Boundary 1\n");
        return(ARB_ILLEGAL_BOUNDARY) ;
    }

    // Check for "loop" edges (same starting and ending nodes)

    for (i=0 ; i<num_edges ; ++i) {
        if (edges[i]->node_id[0] == edges[i]->node_id[1])
         {
                printf ("Illegal Boundary 2\n");
                return(ARB_ILLEGAL_BOUNDARY) ;
        }
    }

    // Now check to make sure that the boundary edges form
    // properly oriented closed loops

    // The algorithm is as follows:
    // 1. Rehash all the edges so that the hash key the first
    //    node on an edge, and the hash valus is the second.
    // 2. Get one edge from the hash table.  Make the
    //    first node the "start" point and set the angle to zero.
    // 3. Use the second node as the key to find the next
    //    edge in chain.  Compute the included angle between the
    //    start point and edge.  Add this to the total angle.
    // 4. Repeat until we get back to the start node.  If
    //    the total angle is positive, this is an external
    //    boundary.  If negative, it is an internal boundary.
    // 5. Check to see if the hash table is empty.  If not,
    //    repeat steps 2-4.
    // 6. Check to make shure that we saw only one external
    //    boundary.

    CArbHashTable<int,int> *pedge ;
    pedge = new CArbHashTable<int,int>(true) ;
    bool have_outside = false ;

    for (i=0 ; i<num_edges ; ++i)
    {
        pedge->Store(edges[i]->node_id[0],edges[i]->node_id[1]) ;
    }

    while (pedge->NumEntries() > 0)
    {
        int *pKey, *pEntry ;
        int start, nd0, nd1 ;
        double angle = 0 ;

        // get the starting node, note that the angle between
        // the start node and the first and last edges is zero.

        pedge->FetchAny(&pKey,&pEntry) ;
        start = *pKey ;
        nd0 = *pEntry ;
        CArbCoord2D base = (pnode_table->Fetch(start))->coord ;
        pedge->RemoveEntry(start,nd0) ;

        // look through all nodes on the loop.
        if ((pEntry = pedge->Fetch(nd0)) == 0)
        {
             printf ("Illegal Boundary 3\n");
            return(ARB_ILLEGAL_BOUNDARY) ;
        }

        nd1 = *pEntry ;
        pedge->RemoveEntry(nd0,nd1) ;

        while (nd1 != start)
        {
            CArbCoord2D i = (pnode_table->Fetch(nd0))->coord ;
            CArbCoord2D j = (pnode_table->Fetch(nd1))->coord ;
            angle += Angle(base,i,j) ;
            nd0 = nd1 ;
            if ((pEntry = pedge->Fetch(nd0)) == 0)
            {
                  printf ("Illegal Boundary 4\n");
                  return(ARB_ILLEGAL_BOUNDARY) ;
            }
            nd1 = *pEntry ;
            pedge->RemoveEntry(nd0,nd1) ;
        }

        // check for outside loop

        if (angle > 0.0)
        {
            if (have_outside)
            {
                 printf ("Illegal Boundary 5\n");
                 return(ARB_ILLEGAL_BOUNDARY) ;
            }
            else
                have_outside = true ;
        }
        else if (angle == 0.0)
        {
             printf ("Illegal Boundary 6\n");
             return(ARB_ILLEGAL_BOUNDARY) ;
        }
    }

    if (!have_outside)
    {
         printf ("Illegal Boundary 7\n");
        return(ARB_ILLEGAL_BOUNDARY) ;
    }

    delete pedge ;
    delete [] edges ;
    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshRegion2D::DisplayBoundary-void-|)
/* ++ ----------------------------------------------------------
**
**    DisplayBoundary - debug routine to display a boundary
**
**      void DisplayBoundary()
**
**      Description: This method provides debugging support. It prints
**          the commands necessary to display the boundary of a region.
**
**
** -- */
void CArbMshRegion2D::DisplayBoundary(const char *label)
{
    int i ;

    int num_edges = pedge_table->NumEntries() ;
    ArbIntEdge **edges = pedge_table->GetKeyList() ;

    // check to make sure that none of the boundary
    // edges cross each other

    for (i=0 ; i<num_edges ; ++i) {

        CArbCoord2D i1, i2 ;

        i1 = (pnode_table->Fetch(edges[i]->node_id[0]))->coord ;
        i2 = (pnode_table->Fetch(edges[i]->node_id[1]))->coord ;

        printf("l %g %g %g %g\n",i1.x(),i1.y(),i2.x(),i2.y()) ;
        printf("t %g %g %d\n",i1.x(),i1.y(),edges[i]->node_id[0]) ;
        printf("# edge: %d %d\n",edges[i]->node_id[0],
                                edges[i]->node_id[1]) ;
    }

//    CArbHashTableIterator<int,ArbIntNode> iter(pnode_table) ;
//    for (iter.First() ; iter.More() ; ++iter) {
//        ArbIntNode *node = iter.Entry() ;
//        printf("t %g %g %d\n",node->coord[0],node->coord[1],node->id) ;
//    }

    printf("a %s\n",label) ;
    printf("n\n") ;
    fflush(stdout) ;
    delete [] edges ;
}




// %(CArbMshRegion2D::GenerateNodes-void-|)
/* ++ ----------------------------------------------------------
**
**    GenerateNodes - generate interior nodes
**
**      void GenerateNodes()
**
**      Description: This method uses a quadtree procedure to generate
**          interior nodes. This is done by:
**
**          1. Generate a quadtree based on the max and min extent of
**          the region. The tree is refined localy so that the region's
**          boundary edges fall in cells who's size about equal to the
**          boundary edge length (tuned by the refine boundary factor).
**
**          2. Refine the quadtree so that no cell is larger than the
**          largest cell on the boundary.
**
**          3. Refine the quadtree so that no two adjacent cells differ
**          in refinement level by more than one.
**
**          4. Attempt to generate a node at the center of each cell,
**          rejecting those that would fall too close to the boundary
**          (tuned by the near boundary factor).
**
**
** -- */
void CArbMshRegion2D::GenerateNodes ( )
{
    int i ;

  // it inserts nodes
  int num_bedges = pedge_table->NumEntries() ;
  if (nnode > num_bedges)
  {
    for (i = num_bedges; i < nnode; i++)
      NewNode (coords[i*2+0], coords[i*2+1], INTERIOR, ARB_FLOATING,true) ;
    return;
  }

    // loop through all the nodes and find the max and min
    // x and y values

    CArbHashTableIterator <int, ArbIntNode> iter(pnode_table);
    double minx,maxx,miny,maxy ;

    minx = maxx = iter.Entry()->coord[0];
    miny = maxy = iter.Entry()->coord[1];
    for (iter++ ; iter.More() ; ++iter)
    {
      if (iter.Entry()->coord[0] < minx) minx = iter.Entry()->coord[0] ;
      if (iter.Entry()->coord[0] > maxx) maxx = iter.Entry()->coord[0] ;
      if (iter.Entry()->coord[1] < miny) miny = iter.Entry()->coord[1] ;
      if (iter.Entry()->coord[1] > maxy) maxy = iter.Entry()->coord[1] ;
    }

    // Create a quadtree that is centered on the boundary and
    // large enough the contain all the nodes plus a little fudge
    // to avoid numerical tolerance troubles.

    double dx = maxx - minx ;
    double dy = maxy - miny ;
    double size = (dx > dy) ? dx : dy ;

    CArbQuadTree *pQTree = new CArbQuadTree((maxx+minx)/2, (maxy+miny)/2, size*1.01);

    // Locally refine the quad tree at the center of each
    // edge.  Refine to size proportional to the edge size.

    int num_edges = pedge_table->NumEntries ( );
    ArbIntEdge **edges = pedge_table->GetKeyList ( );
    double max_edge = 0.0 ;

    for (i=0 ; i<num_edges ; ++i)
    {
      ArbIntNode *p_nd0, *p_nd1 ;
      double dx, dy, cx, cy, len ;

      p_nd0 = pnode_table->Fetch (edges[i]->node_id[0]) ;
      p_nd1 = pnode_table->Fetch (edges[i]->node_id[1]) ;

      dx =  p_nd1->coord[0] - p_nd0->coord[0] ;
      dy =  p_nd1->coord[1] - p_nd0->coord[1] ;
      cx = (p_nd1->coord[0] + p_nd0->coord[0]) / 2 ;
      cy = (p_nd1->coord[1] + p_nd0->coord[1]) / 2 ;

      len = sqrt (dx*dx + dy*dy) ;
      if (len > max_edge) max_edge = len ;

      pQTree->RefineToSize (cx, cy, len*RefineBoundaryFactor);
    }
    delete [] edges ;

    // Now refine the quad tree so that no cell is larger than
    // the larges cell on the boundary.
    pQTree->UniformRefine (max_edge*RefineBoundaryFactor, minx, maxx, miny, maxy);

    // Now refine the quad tree so that no two adjacent cells
    // differ in level of refinement by more than one.
    pQTree->RefineOneLevelDiff ( );


    // Visit all the cells in the quad tree and generate nodes
    // in the center of all the leaf nodes.

    if (QuadTreeDebug)
      printf("# Start Quadtree\n") ;

    if (StartId == 0) StartId = MaxId + 1 ;
    pQTree->VisitLevels ((void *)this, ArbGenerateIntNodes);

    if (QuadTreeDebug)
      printf("# End Quadtree\n");

    delete pQTree ;
}



// %(CArbMshRegion2D::DisplayMesh-void-|)
/* ++ ----------------------------------------------------------
**
**    DisplayMesh - debug routine to display a mesh
**
**      void DisplayMesh()
**
**      Description: This method provides debugging support. It prints
**          the commands necessary to display the mesh.
**
**
** -- */

void CArbMshRegion2D::DisplayMesh(const char *label)
{
    CArbHashTableIterator<int,ArbMshElement2D> eiter(pelem_table) ;

#if 1
    bool do_labels = false ;
//    bool do_labels = true ;

    for (eiter.First() ; eiter.More() ; ++eiter) {
        ArbMshElement2D *elem = eiter.Entry() ;
        double x[8], y[8] ;
        int id[8] ;
        int nn ;
        if (elem->num_nodes == 8)
            nn = 4 ;
        else if (elem->num_nodes == 6)
            nn = 3 ;
        else
            nn = elem->num_nodes ;
        for (int j=0 ; j<elem->num_nodes ; ++j) {
            ArbIntNode *node = pnode_table->Fetch(elem->nodes[j]) ;
            x[j] = node->coord[0] ;
            y[j] = node->coord[1] ;
            id[j] = node->id ;
        }
        for (int jj=0 ; jj<nn ; ++jj) {
            int kk = (jj+1) % nn ;
            printf("l %g %g %g %g\n",x[jj],y[jj],x[kk],y[kk]) ;
        }
        if (do_labels) {
            for (int jj=0 ; jj<elem->num_nodes ; ++jj) {
                printf("t %g %g %d\n",x[jj],y[jj],id[jj]) ;
            }
        }
    }
    printf("a %s\n",label) ;
    printf("n\n") ;
    fflush(stdout) ;
#endif

#if 0
    for (eiter.First() ; eiter.More() ; ++eiter) {
        ArbMshElement2D *elem = eiter.Entry() ;
        double x[8], y[8] ;
        int id[8] ;
        int nn ;
        if (elem->num_nodes == 8)
            nn = 4 ;
        else if (elem->num_nodes == 6)
            nn = 3 ;
        else
            nn = elem->num_nodes ;

        printf("e %d\n",elem->num_nodes) ;

        for (int j=0 ; j<elem->num_nodes ; ++j) {
            ArbIntNode *node = pnode_table->Fetch(elem->nodes[j]) ;
            x[j] = node->coord[0] ;
            y[j] = node->coord[1] ;
            id[j] = node->id ;
            printf(" %d",elem->nodes[j]) ;
        }
        printf("\n") ;

        for (int jj=0 ; jj<nn ; ++jj) {
            int kk = (jj+1) % nn ;
            printf("%g %g\n",x[jj],y[jj]) ;
        }
    }
    printf("a %s\n",label) ;
    printf("n\n") ;
    fflush(stdout) ;
#endif
}


#define MATE_TOL 0.01

void CArbMshRegion2D::GetBdryInfo(
                   int num_nodes,ArbIntNode **nodes,
                   int num_edges,ArbIntEdge **edges,
                   int *node_map,
                   CArbHashTable<int,CArbCoord2D> &bdry_start,
                   CArbHashTable<int,CArbCoord2D> &bdry_stop,
                   CArbHashTable<int,CArbSmallSet<int,2> > &mate_table)
{
    int i,j ;
    double *edge_len = new double[num_nodes] ;
    int *edge_num = new int[num_nodes] ;

    for (i=0 ; i<num_nodes ; ++i) {
        edge_len[i] = 0.0 ;
        edge_num[i] = 0 ;
    }

    // for every node on the boundary we want to find the starting
    // and stopping vectors, between which valid elements can be
    // formed.  These tables are created to deal with the cases
    // where distinct portions of the boundary are coincident
    // e.g., crack faces.

    for (i=0 ; i<num_edges ; ++i) {
        int nd0 = node_map[edges[i]->node_id[0]] ;
        int nd1 = node_map[edges[i]->node_id[1]] ;

        CArbCoord2D v0 =
           (nodes[nd1]->coord-nodes[nd0]->coord).Normalize() ;

        bdry_start.Store(nodes[nd0]->id,v0) ;
        bdry_stop.Store(nodes[nd1]->id,-v0) ;

        double elen = (nodes[nd0]->coord - nodes[nd1]->coord).Magnitude() ;
        edge_len[nd0] += elen ;
        edge_len[nd1] += elen ;
        edge_num[nd0]++ ;
        edge_num[nd1]++ ;
    }

    // find the average edge length

    for (i=0 ; i<num_nodes ; ++i) {
        if (edge_num[i] != 0)
            edge_len[i] /= double(edge_num[i]) ;
    }

    // put all the boundary nodes into a kd tree

    int num_used = 0 ;
    ArbIntNode **used = new ArbIntNode*[num_nodes] ;
    for (i=0 ; i<num_nodes ; ++i) {
        if (nodes[i]->type == BOUNDARY) {
            used[num_used] = nodes[i] ;
            num_used++ ;
        }
    }

    CArbKDTree2D kd_tree(num_used,used) ;

    // now loop through the edges and look for mate nodes

    CArbArray<ArbIntNode*> *range_points = 0 ;
    for (i=0 ; i<num_nodes ; ++i) {

        // if this is a boundary node (average length set)
        // do a range search looking for the node

        if (edge_len[i] > 0.0) {
            ArbIntNode ll,ur ;
            double delt = edge_len[i] * 0.1 ;
            ll.coord[0] = nodes[i]->coord[0] - delt ;
            ll.coord[1] = nodes[i]->coord[1] - delt ;
            ur.coord[0] = nodes[i]->coord[0] + delt ;
            ur.coord[1] = nodes[i]->coord[1] + delt ;
            CArbRectangle rect(ll,ur,delt) ;
            range_points = kd_tree.RangeQuery(rect) ;

            // if one or more nodes were found in the search then
            // check to see if it is within the tolerance.  If so
            // then store the information in the mate table

            for (j=0 ; j<range_points->NumEntries() ; ++j) {
                ArbIntNode *p = (*range_points)[j] ;
                if (p->id == nodes[i]->id) continue ;
                double dist = (nodes[i]->coord - p->coord).Magnitude() ;
                if (dist < edge_len[i]*MATE_TOL) {
                    CArbSmallSet<int,2> *mate =
                            mate_table.Fetch(nodes[i]->id) ;
                    if (mate == 0) {
                        CArbSmallSet<int,2> ss ;
                        ss.Insert(p->id) ;
                        mate_table.Store(nodes[i]->id,ss) ;
                    } else {
                        mate->Insert(p->id) ;
                    }
                }
            }
            delete range_points ;
        }
    }

    delete [] edge_len ;
    delete [] edge_num ;
    delete [] used ;

//     CArbHashTableIterator<int,CArbSmallSet<int,2> > iter(&mate_table) ;
//     for (iter.First() ; iter.More() ; ++iter) {
//         int key = iter.Key() ;
//         CArbSmallSet<int,2> *ss = iter.Entry() ;
//         fprintf(stderr,"mate: %d %d %d\n",key,ss->NumElements(),
//                                           ss->Element(0)) ;
//     }
}


#define BETWEEN_TOL -0.000001

static bool Between(CArbCoord2D &b,CArbCoord2D &i,CArbCoord2D &j)
{
    if (CrossProduct(i,j) >= 0.0) {
        if ((CrossProduct(i,b) > BETWEEN_TOL) &&
            (CrossProduct(b,j) > BETWEEN_TOL)) return(true) ;
    } else {
        if (CrossProduct(i,b) >= 0) return(true) ;
        if (CrossProduct(b,j) > BETWEEN_TOL) return(true) ;
    }
    return(false) ;
}


#define EDGE_FACTOR 0.01

bool CArbMshRegion2D::CheckCross(
                   CArbHeap<ArbIntEdge> *pedge_heap,
                   ArbIntEdge **bdry,ArbIntEdge *entry,
                   ArbIntNode **nodes,int *node_map,
                   ArbIntNode *p,double blensqr,
                   ArbIntNode *p_nd0,ArbIntNode *p_nd1,
                   CArbCoord2D & /*enorm*/,CArbCoord2D & /*emid*/,
                   CArbHashTable<int,CArbCoord2D> &bdry_start,
                   CArbHashTable<int,CArbCoord2D> &bdry_stop,
                   int *lside_0,int *lside_1,
                   bool *nd_0_flg,bool *nd_1_flg,
                   CArbHashTable<int,CArbSmallSet<int,2> > &mate_table)
{
    // get the mates list for the edge

    CArbSmallSet<int,2> *mates0 = mate_table.Fetch(entry->node_id[0]) ;
    CArbSmallSet<int,2> *mates1 = mate_table.Fetch(entry->node_id[1]) ;

    for (int j=0 ; j<pedge_heap->NumEntries() ; ++j) {
        if ((bdry[j]->node_id[0] != entry->node_id[0]) ||
            (bdry[j]->node_id[1] != entry->node_id[1])) {

            // check for any mates

            if ((mates0 != 0) &&
                (mates0->HasElement(bdry[j]->node_id[0]) ||
                 mates0->HasElement(bdry[j]->node_id[1]))) continue ;

            if ((mates1 != 0) &&
                (mates1->HasElement(bdry[j]->node_id[0]) ||
                 mates1->HasElement(bdry[j]->node_id[1]))) continue ;

            ArbIntNode *p_nd2, *p_nd3 ;
            p_nd2 = nodes[node_map[bdry[j]->node_id[0]]] ;
            p_nd3 = nodes[node_map[bdry[j]->node_id[1]]] ;

            double dx, dy, len0sqr, len1sqr ;
            dx = p_nd2->coord[0] - p->coord[0] ;
            dy = p_nd2->coord[1] - p->coord[1] ;
            len0sqr = dx*dx + dy*dy ;
            dx = p_nd3->coord[0] - p->coord[0] ;
            dy = p_nd3->coord[1] - p->coord[1] ;
            len1sqr = dx*dx + dy*dy ;

            if ((len0sqr > blensqr) && (len1sqr > blensqr)) {
                dx = p_nd2->coord[0] - p_nd0->coord[0] ;
                dy = p_nd2->coord[1] - p_nd0->coord[1] ;
                len0sqr = dx*dx + dy*dy ;
                dx = p_nd3->coord[0] - p_nd0->coord[0] ;
                dy = p_nd3->coord[1] - p_nd0->coord[1] ;
                len1sqr = dx*dx + dy*dy ;
                if ((len0sqr > blensqr) && (len1sqr > blensqr) &&
                     Cross(p_nd0->coord,p->coord,
                           p_nd2->coord,p_nd3->coord)) {
                    double elen ;
                    double dist = DistSqr(p->coord,
                                          p_nd2->coord,
                                          p_nd3->coord,
                                          &elen) ;
                    if (sqrt(dist) < elen*EDGE_FACTOR) {
                        if (len0sqr < len1sqr) {
                            CArbCoord2D *start = bdry_start.Fetch(p_nd2->id) ;
                            if (start != 0) {
                                CArbCoord2D *stop = bdry_stop.Fetch(p_nd2->id) ;
                                CArbCoord2D v0 =
                                   (p_nd0->coord-p->coord).Normalize() ;
                                CArbCoord2D v1 =
                                   (p_nd1->coord-p->coord).Normalize() ;
                                if (Between(v0,*start,*stop) &&
                                    Between(v1,*start,*stop)) return(true) ;
                            } else {
                                return(true) ;
                            }
                        } else {
                            CArbCoord2D *start = bdry_start.Fetch(p_nd3->id) ;
                            if (start != 0) {
                                CArbCoord2D *stop = bdry_stop.Fetch(p_nd3->id) ;
                                CArbCoord2D v0 =
                                   (p_nd0->coord-p->coord).Normalize() ;
                                CArbCoord2D v1 =
                                   (p_nd1->coord-p->coord).Normalize() ;
                                if (Between(v0,*start,*stop) &&
                                    Between(v1,*start,*stop)) return(true) ;
                            } else {
                                return(true) ;
                            }
                        }
                    } else {
                        return(true) ;
                    }
                }

                dx = p_nd3->coord[0] - p_nd1->coord[0] ;
                dy = p_nd3->coord[1] - p_nd1->coord[1] ;
                len0sqr = dx*dx + dy*dy ;
                dx = p_nd2->coord[0] - p_nd1->coord[0] ;
                dy = p_nd2->coord[1] - p_nd1->coord[1] ;
                len1sqr = dx*dx + dy*dy ;
                if ((len0sqr > blensqr) && (len1sqr > blensqr) &&
                     Cross(p->coord,p_nd1->coord,
                          p_nd2->coord,p_nd3->coord)) {
                    double elen ;
                    double dist = DistSqr(p->coord,
                                          p_nd2->coord,
                                          p_nd3->coord,
                                          &elen) ;
                    if (sqrt(dist) < elen*EDGE_FACTOR) {
                        if (len0sqr < len1sqr) {
                            CArbCoord2D *start = bdry_start.Fetch(p_nd2->id) ;
                            if (start != 0) {
                                CArbCoord2D *stop = bdry_stop.Fetch(p_nd2->id) ;
                                CArbCoord2D v0 =
                                   (p_nd0->coord-p->coord).Normalize() ;
                                CArbCoord2D v1 =
                                   (p_nd1->coord-p->coord).Normalize() ;
                                if (Between(v0,*start,*stop) &&
                                    Between(v1,*start,*stop)) return(true) ;
                            } else {
                                return(true) ;
                            }
                        } else {
                            CArbCoord2D *start = bdry_start.Fetch(p_nd3->id) ;
                            if (start != 0) {
                                CArbCoord2D *stop = bdry_stop.Fetch(p_nd3->id) ;
                                CArbCoord2D v0 =
                                   (p_nd0->coord-p->coord).Normalize() ;
                                CArbCoord2D v1 =
                                   (p_nd1->coord-p->coord).Normalize() ;
                                if (Between(v0,*start,*stop) &&
                                    Between(v1,*start,*stop)) return(true) ;
                            } else {
                                return(true) ;
                            }
                        }
                    } else {
                        return(true) ;
                    }
                }
            }

            // check to see if this is a pre-existing
            // edge of the triangle.

            if (((p_nd0->id == p_nd2->id) &&
                 (p->id == p_nd3->id)) ||
                ((p->id == p_nd2->id) &&
                 (p_nd0->id == p_nd3->id))) {
                *nd_0_flg = true ;
                *lside_0 = bdry[j]->node_id[2] ;
            }

            if (((p_nd1->id == p_nd2->id) &&
                 (p->id == p_nd3->id)) ||
                ((p->id == p_nd2->id) &&
                 (p_nd1->id == p_nd3->id))) {
                *nd_1_flg = true ;
                *lside_1 = bdry[j]->node_id[2] ;
            }
        }
    }
    return(false) ;
}


// %(CArbMshRegion2D::BoundaryContraction-void-|)
/* ++ ----------------------------------------------------------
**
**    BoundaryContraction - do boundary contraction mesh generation
**
**      void BoundaryContraction()
**
**      Description: This method generates a triangular mesh for the
**          region by performing a boundary contraction (advancing
**          front) procedure. The procedure is to select an edge from
**          the current active boundary and find the node, which when
**          combined with this edge, forms the best triangle. The
**          boundary is updated and the procedure is repeated.
**
**
** -- */

static int ArbCmpEdge(const ArbIntEdge &ed1,const ArbIntEdge &ed2)
{
    if (ed1.length > ed2.length) return(1) ;
    if (ed1.length < ed2.length) return(-1) ;
    return(0) ;
}

#define RANGE_FACTOR_1 5.0
#define RANGE_FACTOR_2 10.0
#define RANGE_FACTOR_3 20.0

// this is the tolerance used to determine if two nodes are at
// the same location.

#define ANGLE_TOL 0.0017
//#define ANGLE_TOL 0.002
#define POS_TOL 0.005
//#define POS_TOL 0.0025
//#define POS_TOL 0.001
//#define POS_TOL 0.01
#define MIN_ANGLE_TOL 0.01

void CArbMshRegion2D::BoundaryContraction ( )
{
    ArbIntNode **nodes ;
    int        *node_map ;
    CArbArray<ArbIntNode*> *range_points = 0 ;
    int valid, retry=0 ;// flags for checking elements
    int i ;

    int tnum = 0 ;

    // get a list of all the nodes

    int num_nodes = pnode_table->NumEntries() ;
    nodes = pnode_table->GetEntryList() ;

    CArbHashTableIterator<int,ArbIntNode> iter(pnode_table) ;

    // Find the max node ID
    int max_id = 0 ;
    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.Entry()->id > max_id) max_id = iter.Entry()->id ;
    }

    node_map = new int[max_id+1];
    for (i=0,iter.First() ; iter.More() ; ++i,++iter) {
      node_map[iter.Entry()->id]=i;
    }

    // place all the boundary edges in a heap data
    // structure so that we can extract the smallest edge

    CArbHeap <ArbIntEdge> *pedge_heap ;
    pedge_heap = new CArbHeap <ArbIntEdge> (ArbCmpEdge);

    int num_edges = pedge_table->NumEntries ( );
    ArbIntEdge **edges = pedge_table->GetKeyList() ;

    for (i=0 ; i<num_edges ; ++i)
    {
        ArbIntNode *p_nd0, *p_nd1 ;
        double dx, dy ;

        p_nd0 = nodes[node_map[edges[i]->node_id[0]]] ;
        p_nd1 = nodes[node_map[edges[i]->node_id[1]]] ;

        dx = p_nd1->coord[0] - p_nd0->coord[0] ;
        dy = p_nd1->coord[1] - p_nd0->coord[1] ;
        edges[i]->length = sqrt(dx*dx + dy*dy) ;

        pedge_heap->Insert(*(edges[i])) ;
    }

    // here we want to preprocess things for the case where we
    // have two surfaces that are coincident but not connected.
    // This may happen in the case of crack faces or slide line
    // type contacts.  For each node that is adjacent to a boundary
    // we create an entry in a hash table indexed by the node's
    // id.  The value of the entry is the negative of the normal
    // to the adjacent edge.

    CArbHashTable <int, CArbCoord2D> bdry_start ;
    CArbHashTable <int, CArbCoord2D> bdry_stop ;
    MateTable = new CArbHashTable <int, CArbSmallSet <int, 2> > ;

    GetBdryInfo(num_nodes,nodes,num_edges,edges,
                node_map,bdry_start,bdry_stop,
                *MateTable) ;

    //
    // it inserts a external triangular mesh
    if (nelem > 0)
    {
       // adds elements
       for (i = 0; i < nelem; i++)
       {
         ArbMshElement2D elem;
         elem.elem_id = NewElemNum();
         elem.mat_id = MatId;
         elem.num_nodes = 3;
         elem.nodes[0] = conn[i*4+3];
         elem.nodes[1] = conn[i*4+2];
         elem.nodes[2] = conn[i*4+1];

         pelem_table->Store (elem.elem_id, elem);
       }

       delete [] nodes ;
       delete [] node_map ;
       delete [] edges ;
       delete pedge_heap ;
       return;
    }

    // get a list of all the nodes

    // build a TwoDTree

    CArbKDTree2D kd_tree(num_nodes,nodes) ;

    // loop until there are no more edges in the boundary

    while (pedge_heap->NumEntries() > 0)
    {
        double max_angle = 0.0 ;
        bool flg_0 = false ;
        bool flg_1 = false ;
        int tri_nodes[6] ;
        int side_0 = -1 ;
        int side_1 = -1 ;
        ArbIntEdge *entry, an_edge ;
        ArbIntNode *p_nd0, *p_nd1 ;
        double dx, dy, blen, blensqr ;

        // get the shortest edge in the boundary

        entry = pedge_heap->GetMin() ;
        p_nd0 = nodes[node_map[entry->node_id[0]]] ;
        p_nd1 = nodes[node_map[entry->node_id[1]]] ;
        dx = p_nd1->coord[0] - p_nd0->coord[0] ;
        dy = p_nd1->coord[1] - p_nd0->coord[1] ;
        blen = sqrt(dx*dx + dy*dy) * POS_TOL ;
        blensqr = blen * blen ;
        CArbCoord2D emid = CArbCoord2D (0.5*(p_nd1->coord[0]+p_nd0->coord[0]),
                                        0.5*(p_nd1->coord[1]+p_nd0->coord[1])) ;
        CArbCoord2D enorm = CArbCoord2D(-dy,dx) ;

        // define a range query rectangle around the edge
        // and get the points that fall in range
        if ( retry==0 )
        {
            CArbRectangle rect(*p_nd0, *p_nd1, RANGE_FACTOR_1*entry->length);
            range_points = kd_tree.RangeQuery (rect) ;

            // if there are no points in range other than the
            // edge end points, increase the range size
            if ( range_points->NumEntries() <= 3 )
            {
                delete range_points;
                CArbRectangle rect (*p_nd0, *p_nd1, RANGE_FACTOR_2*entry->length);
                range_points = kd_tree.RangeQuery(rect) ;
            }
        }
        else if ( retry==1 )
        {
            CArbRectangle rect (*p_nd0, *p_nd1, RANGE_FACTOR_3*entry->length);
            range_points = kd_tree.RangeQuery (rect) ;
        }
        else if ( retry==2 )
        {
            range_points = kd_tree.ReturnAll () ;
        }


        // now loop through all the nodes and find the node
        // that will make a triangle with the minimum included
        // angle

        valid = 0 ;
        int cur, num ;
        num = range_points->NumEntries() ;
        for (cur=0 ; cur < num ; ++cur)
        {
          ArbIntNode *p = (*range_points)[cur] ;

          if (p->corner)
          {

            // inserted to deal with nodes adjacent to an edge.

            CArbCoord2D *start = bdry_start.Fetch(p->id) ;
            if (start != 0)
            {
                CArbCoord2D *stop = bdry_stop.Fetch(p->id) ;
                CArbCoord2D v0 = (p_nd0->coord-p->coord).Normalize() ;
                if (!Between(v0,*start,*stop)) continue ;
                CArbCoord2D v1 = (p_nd1->coord-p->coord).Normalize() ;
                if (!Between(v1,*start,*stop)) continue ;
            }

            bool duplicate_node = false ;
            double angle = 0.0 ;
            double dx, dy, len0, len1 ;

            dx = p_nd0->coord[0] - p->coord[0] ;
            dy = p_nd0->coord[1] - p->coord[1] ;
            len0 = sqrt(dx*dx + dy*dy) ;

            dx = p_nd1->coord[0] - p->coord[0] ;
            dy = p_nd1->coord[1] - p->coord[1] ;
            len1 = sqrt(dx*dx + dy*dy) ;

            // reject the node if it one of the edge's end nodes
            // or if it is to the right of the edge

            if ((p->id != entry->node_id[0]) && (p->id != entry->node_id[1]) &&
                (len0 > blen) && (len1 > blen) &&
                (CrossProd(p->coord, p_nd0->coord,p_nd1->coord) > 0.0))
            {

                angle = Angle (p->coord, p_nd0->coord,p_nd1->coord) ;

                // check to see if this is the
                // biggest angle we've seen so far, check to
                // see if this is a valid triangle (i.e., it
                // does not cross the existing boundary

                if ((angle > max_angle) || (duplicate_node))
                {
                    ArbIntEdge **bdry = pedge_heap->GetEntryList() ;
                    bool cross = false ;
                    bool nd_0_flg = false ;
                    bool nd_1_flg = false ;
                    int lside_0 = -1 ;
                    int lside_1 = -1 ;

                    cross = CheckCross(pedge_heap,bdry,
                                entry,nodes,node_map,p,blensqr,
                                p_nd0,p_nd1,enorm,emid,bdry_start,
                                bdry_stop,&lside_0,&lside_1,
                                &nd_0_flg,&nd_1_flg,*MateTable) ;

                    // when we get here, either we have crossed an
                    // edge so we forget this node or we have run out
                    // edges

                    if ((!cross) && (angle > 1e-5))
                    {
                        max_angle = angle ;
                        tri_nodes[0] = p_nd0->id ;
                        tri_nodes[1] = p_nd1->id ;
                        tri_nodes[2] = p->id ;
                        tri_nodes[3] = entry->node_id[2] ;
                        flg_0 = nd_0_flg ;
                        flg_1 = nd_1_flg ;
                        side_0 = lside_0 ;
                        side_1 = lside_1 ;
                        valid = 1 ;
                    }
                    delete [] bdry ;
                }
            }
          }
        }

        // when we get here we have the best triangle we can form
        // so add it to the list and update the boundary
        if ( !valid )
        {
          if ( retry==2 ) return ;
          else
          {
            retry++;
            pedge_heap->Insert(*entry) ;
          }
        }
        else
        {
          retry = 0;
          ArbMshElement2D elem ;

          elem.elem_id = NewElemNum() ;
          elem.mat_id = MatId ;
          if (Order == LINEAR)
          {
            elem.num_nodes = 3 ;
            for (int i=0 ; i<3 ; ++i) elem.nodes[i] = tri_nodes[i] ;
          }
          else
          {
            if (side_0 == -1)
            {
                double x, y ;
                p_nd0 = nodes[node_map[tri_nodes[0]]] ;
                p_nd1 = nodes[node_map[tri_nodes[2]]] ;
                x = (p_nd1->coord[0] + p_nd0->coord[0]) / 2 ;
                y = (p_nd1->coord[1] + p_nd0->coord[1]) / 2 ;
                tri_nodes[5] = NewNode (x, y, INTERIOR, ARB_FLOATING, false) ;
            }
            else
            {
                tri_nodes[5] = (int)side_0 ;
            }

            if (side_1 == -1)
            {
                double x, y ;
                p_nd0 = nodes[node_map[tri_nodes[1]]] ;
                p_nd1 = nodes[node_map[tri_nodes[2]]] ;
                x = (p_nd1->coord[0] + p_nd0->coord[0]) / 2 ;
                y = (p_nd1->coord[1] + p_nd0->coord[1]) / 2 ;
                tri_nodes[4] = NewNode(x,y,INTERIOR,ARB_FLOATING,false) ;
            } else
            {
                tri_nodes[4] = (int)side_1 ;
            }

            elem.num_nodes = 6 ;
            for (int i=0 ; i<6 ; ++i) elem.nodes[i] = tri_nodes[i] ;
          }

          pelem_table->Store(elem.elem_id,elem) ;

          if (!flg_0)
          {
            double dx, dy ;
            an_edge.node_id[0] = tri_nodes[0] ;
            an_edge.node_id[1] = tri_nodes[2] ;
            an_edge.node_id[2] = tri_nodes[5] ;
            p_nd0 = nodes[node_map[an_edge.node_id[0]]] ;
            p_nd1 = nodes[node_map[an_edge.node_id[1]]] ;
            dx = p_nd1->coord[0] - p_nd0->coord[0] ;
            dy = p_nd1->coord[1] - p_nd0->coord[1] ;
            an_edge.length = sqrt(dx*dx + dy*dy) ;
            pedge_heap->Insert(an_edge) ;
          }
          else
          {
            an_edge.node_id[0] = tri_nodes[2] ;
            an_edge.node_id[1] = tri_nodes[0] ;
            pedge_heap->Remove(an_edge) ;
          }

          if (!flg_1)
          {
            double dx, dy ;
            an_edge.node_id[0] = tri_nodes[2] ;
            an_edge.node_id[1] = tri_nodes[1] ;
            an_edge.node_id[2] = tri_nodes[4] ;
            p_nd0 = nodes[node_map[an_edge.node_id[0]]] ;
            p_nd1 = nodes[node_map[an_edge.node_id[1]]] ;
            dx = p_nd1->coord[0] - p_nd0->coord[0] ;
            dy = p_nd1->coord[1] - p_nd0->coord[1] ;
            an_edge.length = sqrt(dx*dx + dy*dy) ;
            pedge_heap->Insert(an_edge) ;
          }
          else
          {
            an_edge.node_id[0] = tri_nodes[1] ;
            an_edge.node_id[1] = tri_nodes[2] ;
            pedge_heap->Remove(an_edge) ;
          }
        }

        delete range_points;
        ++tnum ;

    }

    delete [] nodes ;
    delete [] node_map ;
    delete [] edges ;
    delete pedge_heap ;
}


// %(CArbMshRegion2D::NewNode-int-|-double-|-double-|-ArbMshNodeType-|-ArbMshNodeMotion-|-bool-|)
/* ++ ----------------------------------------------------------
**
**    NewNode - generate a new node
**
**      int NewNode(
**              double           x,
**              double           y,
**              ArbMshNodeType   type,
**              ArbMshNodeMotion motion,
**              bool             corner)
**
**        x      - (in)  the node's x coordinate
**        y      - (in)  the node's y coordinate
**        type   - (in)  BOUNDARY or INTERIOR
**        motion - (in)  ARB_FIXED or ARB_FLOATING
**        corner - (in)  corner node flag
**
**      Description: This method generates a new node and assigns it a
**          unique node id.
**
**      Return Value: the new node id
**
**
** -- */

int CArbMshRegion2D::NewNode(double x,double y,
                                      ArbMshNodeType type,
                                      ArbMshNodeMotion motion,
                                      bool corner)
{
    int id ;
    ArbIntNode int_node ;

    id = int_node.id = StartId ;
    int_node.coord[0] = x ;
    int_node.coord[1] = y ;
    int_node.type = type ;
    int_node.motion = motion ;
    int_node.corner = corner ;

    while (pnode_table->Store(id,int_node) == false) {
        ++id ;
        ++int_node.id ;
    }
    StartId = id + 1 ;
    return(id) ;
}




// %(CArbMshRegion2D::DuplicateNode-int-|-int-|)
/* ++ ----------------------------------------------------------
**
**    DuplicateNode - generate a second node a coordinate location
**
**      int DuplicateNode(int id)
**
**        id - (in)  existing node id
**
**      Description: Given the id of an existing node, this method
**          generates a new node that will share the coordinates off
**          the existing node, but have a unique node id.
**
**      Return Value: the id of the generated node
**
**
** -- */

int CArbMshRegion2D::DuplicateNode(int id)
{
    ArbIntNode *node = pnode_table->Fetch(id) ;
    assert(node != 0) ;
    return(NewNode(node->coord[0],node->coord[1],node->type,
                   node->motion,node->corner)) ;
}




// %(CArbMshRegion2D::SmoothNodes-void-|)
/* ++ ----------------------------------------------------------
**
**    SmoothNodes - do nodal smoothing
**
**      void SmoothNodes()
**
**      Description: This method performs smoothing (repositioning) of
**          internal nodes to improve element quality. This can be done
**          using either a Laplace or a Winslow (default) algorithm.
**
**
** -- */

void CArbMshRegion2D::SmoothNodes()
{
    CArbMshSmooth2D smooth(pnode_table,pelem_table) ;
    if (WinslowSmoothing)
        smooth.SmoothNodesWinslow() ;
    else
        smooth.SmoothNodesLaplace() ;
}




// %(CArbMshRegion2D::GetMeshStatistics-ArbMshStats-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetMeshStatistics - return element shape statistics
**
**      ArbMshStats *GetMeshStatistics() const
**
**      Description: This method returns statistics regarding the
**          quality of the elements generated for the region. Ownership
**          of this structure passes to the client, which must
**          eventually call delete.
**
**      Return Value: statistics regarding element quality
**
**
** -- */

ArbMshStats *CArbMshRegion2D::GetMeshStatistics() const
{
    ArbMshStats *stats = new ArbMshStats ;

    stats->num_elements = pelem_table->NumEntries() ;
    stats->num_nodes = pnode_table->NumEntries() ;

    ArbMshElement2D **elems ;
    elems = pelem_table->GetEntryList() ;

    double sm, sm_sum = 0.0 ;
    stats->minimum_shape_measure = 100.0 ;

    int i ;
    for (i=0 ; i<11 ; ++i) stats->num_less_than[i] = 0 ;

    for (i=0 ; i<stats->num_elements ; ++i) {
        if ((elems[i]->num_nodes == 3) ||
            (elems[i]->num_nodes == 6)) {

            // compute the element shape measure.  We are
            // using the mean edge ratio, which is defined
            // as sm = (4*sqrt(3)*Area)/(l01^2+l12^2+l20^2)
            // where lij is the length of the edge from node
            // i to j.

            ArbIntNode *nd0, *nd1, *nd2 ;
            nd0 = pnode_table->Fetch(elems[i]->nodes[0]) ;
            nd1 = pnode_table->Fetch(elems[i]->nodes[1]) ;
            nd2 = pnode_table->Fetch(elems[i]->nodes[2]) ;
            sm = TriShapeMeasure(nd0->coord,nd1->coord,nd2->coord) ;
            if (sm <= 0.0) {
                fprintf(stderr,"Invalid Element: %d\n", elems[i]->elem_id) ;
            }
        } else {
            ArbIntNode *nd0, *nd1, *nd2, *nd3 ;
            nd0 = pnode_table->Fetch(elems[i]->nodes[0]) ;
            nd1 = pnode_table->Fetch(elems[i]->nodes[1]) ;
            nd2 = pnode_table->Fetch(elems[i]->nodes[2]) ;
            nd3 = pnode_table->Fetch(elems[i]->nodes[3]) ;
            sm = QuadMetric(nd0->coord,nd1->coord,
                            nd2->coord,nd3->coord) ;
            if (sm <= 0.0) {
                fprintf(stderr,"Invalid Element: %d\n",elems[i]->elem_id) ;
            }
        }

        sm_sum += sm ;
        if (sm < stats->minimum_shape_measure)
            stats->minimum_shape_measure = sm ;
        for (int i = 1 ; i < 11 ; ++i) {
            double x = double(i) / 10.0 ;
            if (sm < x) stats->num_less_than[i]++ ;
        }
    }

    stats->mean_shape_measure = sm_sum / stats->num_elements ;

    delete [] elems ;
    return(stats) ;
}




// %(CArbMshRegion2D::Cross-bool-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    Cross - check to see if lines cross
**
**      bool Cross(
**              CArbCoord2D i1,
**              CArbCoord2D i2,
**              CArbCoord2D j1,
**              CArbCoord2D j2) const
**
**        i1 - (in)  first node of first line
**        i2 - (in)  second node of first line
**        j1 - (in)  first node of second line
**        j2 - (in)  second node of second line
**
**      Description: This method check to see if two line segments
**          cross.
**
**      Return Value: true if they cross
**
**
** -- */

bool CArbMshRegion2D::Cross(const CArbCoord2D &i1,
                            const CArbCoord2D &i2,
                            const CArbCoord2D &j1,
                            const CArbCoord2D &j2) const
{
    /* This is the assumed node configuration:

           \I1   /J2
            \   /
             \ /
              \
             / \
            /   \
           /J1   \I2
    */

    // first do a simple min max test

    if ((i1[0] > j1[0]) && (i1[0] > j2[0]) &&
        (i2[0] > j1[0]) && (i2[0] > j2[0])) return false ;

    if ((i1[1] > j1[1]) && (i1[1] > j2[1]) &&
        (i2[1] > j1[1]) && (i2[1] > j2[1])) return false ;

    // now check to make sure that J1 and J2 are on opposite sides
    // of line I.

    CArbCoord2D delt = i2 - i1 ;
    double tol = -1e-10 * delt.Magnitude() ;

    int sj1 = 0 ;
    int sj2 = 0 ;

    double vj1 = CrossProd(i1,i2,j1) ;
    double vj2 = CrossProd(i1,i2,j2) ;

    if (vj1 >= tol) {
        sj1 = 1 ;
    } else if (vj1 <= -tol) {
        sj1 = -1 ;
    }

    if (vj2 >= tol) {
        sj2 = 1 ;
    } else if (vj2 <= -tol) {
        sj2 = -1 ;
    }

    if (sj1 * sj2 != -1) return false ;

    // now check to make sure that I1 and I2 are on opposite sides
    // of line J.

    int si1 = 0 ;
    int si2 = 0 ;

    double vi1 = CrossProd(j1,j2,i1) ;
    double vi2 = CrossProd(j1,j2,i2) ;

    if (vi1 >= tol) {
        si1 = 1 ;
    } else if (vi1 <= -tol) {
        si1 = -1 ;
    }

    if (vi2 >= tol) {
        si2 = 1 ;
    } else if (vi2 <= -tol) {
        si2 = -1 ;
    }

    return (si1 * si2 < 0) ;
}




// %(CArbMshRegion2D::CrossProd-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    CrossProd - compute a cross product
**
**      double CrossProd(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j) const
**
**        b - (in)  base vertex coordinate
**        i - (in)  end of first vector
**        j - (in)  end of second vector
**
**      Description: This method computes the cross product of two
**          vector that share a common starting vertex.
**
**      Return Value: the (2D) cross product value
**
**
** -- */

double CArbMshRegion2D::CrossProd(const CArbCoord2D &b,
                                  const CArbCoord2D &i,
                                  const CArbCoord2D &j) const
{
    double cross ;

    cross = ((i[0] - b[0]) * (j[1] - b[1])) -
            ((i[1] - b[1]) * (j[0] - b[0])) ;
    return cross ;
}




// %(CArbMshRegion2D::Angle-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    Angle - compute an angle
**
**      double Angle(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j) const
**
**        b - (in)  base vertex coordinate
**        i - (in)  end of first vector
**        j - (in)  end of second vector
**
**      Description: Compute the magnitude of the included angle
**          between two vectors that share a common starting vertex (in
**          the range -pi to pi).
**
**      Return Value: the magnitude of the angle (in radians)
**
**
** -- */

double CArbMshRegion2D::Angle(const CArbCoord2D &b,
                              const CArbCoord2D &i,
                              const CArbCoord2D &j) const
{
    CArbCoord2D bi = i - b ;
    CArbCoord2D bj = j - b ;

    double tmp = ((bi.x() * bj.x()) + (bi.y() * bj.y())) /
                 (sqrt(bi.x()*bi.x() + bi.y()*bi.y()) *
                  sqrt(bj.x()*bj.x() + bj.y()*bj.y())) ;
    if (tmp >  1.0) tmp =  1.0 ;
    if (tmp < -1.0) tmp = -1.0 ;

    double cross = CrossProd(b,i,j) ;

    if (cross == 0.0) {
        return((tmp > 0) ? 0.0 : acos(-1.0)) ;
    }
    return(acos(tmp) * (cross/fabs(cross))) ;
}




// %(CArbMshRegion2D::Angle2Pi-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    Angle2Pi - compute an angle
**
**      double Angle2Pi(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j) const
**
**        b - (in)  base vertex coordinate
**        i - (in)  end of first vector
**        j - (in)  end of second vector
**
**      Description: Compute the magnitude of the included angle
**          between two vectors that share a common starting vertex (in
**          the range 0 to 2*pi).
**
**      Return Value: the magnitude of the angle (in radians)
**
**
** -- */

double CArbMshRegion2D::Angle2Pi(const CArbCoord2D &b,
                                 const CArbCoord2D &i,
                                 const CArbCoord2D &j) const
{
    CArbCoord2D bi = i - b ;
    CArbCoord2D bj = j - b ;

    double tmp = ((bi.x() * bj.x()) + (bi.y() * bj.y())) /
                 (sqrt(bi.x()*bi.x() + bi.y()*bi.y()) *
                  sqrt(bj.x()*bj.x() + bj.y()*bj.y())) ;
    if (tmp >  1.0) tmp =  1.0 ;
    if (tmp < -1.0) tmp = -1.0 ;
    double cross = CrossProd(b,i,j) ;

    if (cross == 0.0) {
        return((tmp > 0) ? 0.0 : acos(-1.0)) ;
    }

    return((cross > 0.0) ? acos(tmp) : TWO_PI - acos(tmp)) ;
}




// %(CArbMshRegion2D::Area-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    Area - find the area of a triangle
**
**      double Area(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j) const
**
**        b - (in)  first vertex
**        i - (in)  second vertex
**        j - (in)  third vertex
**
**      Description: Compute the area of a triangle.
**
**      Return Value: triangle area
**
**
** -- */

double CArbMshRegion2D::Area(const CArbCoord2D &b,
                             const CArbCoord2D &i,
                             const CArbCoord2D &j) const
{
    double area ;

    area = (i[0]-b[0])*(j[1]-b[1]) - (j[0]-b[0])*(i[1]-b[1]) ;

    return 0.5*area ;
}




// %(CArbMshRegion2D::DistSqr-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|-double-|*)
/* ++ ----------------------------------------------------------
**
**    DistSqr - compute the distance from a point to a line segment
**
**      double DistSqr(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j,
**              double      *blen) const
**
**        b    - (in)  input vertex
**        i    - (in)  first line segment vertex
**        j    - (in)  second line segment vertex
**        blen - (out) distance from i to j
**
**      Description: This method computes the quare of the minimum
**          distance from a point to a line segment
**
**      Return Value: the square of the distance
**
**
** -- */

double CArbMshRegion2D::DistSqr(
        const CArbCoord2D &b,
        const CArbCoord2D &i,
        const CArbCoord2D &j,
        double *blen) const
{
  /* This routine finds the perpendicular minimum distance
     between a point (b) and a line segment (i,j).  It also
     returns the length of the triangle base (i,j) in *blen

     This is the picture:

                      b                       b
                    *                       *
        case 1     /|\           case 2    /
                  / | \   or              /
                 /  |  \                 /
                *---+---*       *-------*
               i    k    j     i        j,k

    the variables are defined as follows:

    blength = |ij|
    area    = area of triangle ijb
    pdist   = |bk|
    adist1  = |bi|
    adist2  = |bj|
  */

    double twice_area, pdist, adist1, adist2 ;

    CArbCoord2D base = j - i ;
    *blen = base.Magnitude() ;

    // check for case 1 or case 2

    if ((b-i)*base * (b-j)*base < 0.0) {

        twice_area = fabs((i[0]-b[0])*(j[1]-b[1]) -
                          (j[0]-b[0])*(i[1]-b[1])) ;
        pdist = twice_area / *blen ;
        pdist = pdist * pdist ;
        return(pdist) ;

    }

    // case 2, find distance to the ends of the interface

    adist1 = (b - i).Magnitude() ;
    adist2 = (b - j).Magnitude() ;

    return(adist1 < adist2 ? adist1 : adist2) ;
}




// %(CArbMshRegion2D::ScanCross-bool-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    ScanCross - check for a scan line crossing
**
**      bool ScanCross(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j) const
**
**        b - (in)  input point
**        i - (in)  first line segment vertex
**        j - (in)  second line segment vertex
**
**      Description: This method determines if a horizontal line drawn
**          from minus infinity to a given point crosses a given line
**
**      Return Value: true if crosses
**
**
** -- */

bool CArbMshRegion2D::ScanCross(const CArbCoord2D &b,
                                const CArbCoord2D &i,
                                const CArbCoord2D &j) const
{
  /*
      Given the coordinates of the vertices of a straight
      edge and a point location, this function returns a
      status to indicate if a scan line drawn from the
      point to minus infinity crosses the line segement.
  */

    double y_min,y_max,x_min ;

    // first reject if the line segment is to the right of
    // the point

    if ((i[0] > b[0]) && (j[0] > b[0])) return(false) ;

    // find the minimum and maximum values of the edge

    if (i[1] == j[1]) {                // horizontal line segment
        return(false) ;
    } else if (i[1] <= j[1]) {
        y_max = j[1] ;   y_min = i[1] ;   x_min = i[0] ;
    } else {
        y_max = i[1] ;   y_min = j[1] ;   x_min = j[0] ;
    }

    if (i[0] == j[0]) {                // vertical line segment
        if ((b[1] >= y_min) && (b[1] < y_max))
            return(true) ;
        else
            return(false) ;
    }

    if ((b[1] == y_min) && (b[0] > x_min)) return(true) ;

    // now reject if the line segment is above or below
    // the point

    if ((b[1] > y_max) || (b[1] < y_min)) return(false) ;

    // Find the x coordinates of the intersection and see
    // if this is to the left of the point

    double m = (i[1] - j[1]) / (i[0] - j[0]) ;
    double c = i[1] - m * i[0] ;
    double xint = (b[1] - c)/m ;

    if (xint < b[0]) return(true) ;

    return(false) ;
}




// %(CArbMshRegion2D::TriShapeMeasure-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    TriShapeMeasure - compute a shape measure
**
**      double TriShapeMeasure(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j) const
**
**        b - (in)  first vertex
**        i - (in)  second vertex
**        j - (in)  third vertex
**
**      Description: This method computes a shape measure for a
**          triangular element.
**
**      Return Value: triangle's shape measure
**
**
** -- */

double CArbMshRegion2D::TriShapeMeasure(const CArbCoord2D &i,
                                        const CArbCoord2D &j,
                                        const CArbCoord2D &k) const
{
    // compute the element shape measure.  We are using the mean edge ratio,
    // which is defined as sm = (4*sqrt(3)*Area)/(l01^2+l12^2+l20^2)
    // where lij is the length of the edge from node i to j.

    static const double factor = 6.92820323 ;
    double len_sqr_sum = 0 ;

    CArbCoord2D d = i - j ;
    len_sqr_sum += d.x()*d.x() + d.y()*d.y() ;
    d = j - k ;
    len_sqr_sum += d.x()*d.x() + d.y()*d.y() ;
    d = k - i ;
    len_sqr_sum += d.x()*d.x() + d.y()*d.y() ;

    return(factor * Area(i,j,k) / len_sqr_sum) ;
}




// %(CArbMshRegion2D::BisectNorm-CArbCoord2D-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    BisectNorm - find the bisector of an angle
**
**      CArbCoord2D BisectNorm(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j) const
**
**        b - (in)  base vertex coordinate
**        i - (in)  end of first vector
**        j - (in)  end of second vector
**
**      Description: This method finds the normalized direction of a
**          bisector of an angle formed by to vectors that share a
**          common starting vertex.
**
**      Return Value: a normalized vector in the bisector direction
**
**
** -- */

CArbCoord2D CArbMshRegion2D::BisectNorm(const CArbCoord2D &b,
                                        const CArbCoord2D &i,
                                        const CArbCoord2D &j) const
{
    /* This is the assumed node configuration:

            * i
             \
              \
               \
                \
                 *----------*
                b          j

       We want to find the normalized vector that is the bisector of
       angle ibj.  The procedure is as follows:

       1. normalize the vectors bi and bj.
       2. find the "right hand" normal to bj (that is, bj cross bj_normal
          is positive).
       3. find the "left hand" normal to bi
       4. Average the normals and renormalize.

       There is only one degenerate case, were the include angle is 0.
       for this case return the normal at 180.
    */

    CArbCoord2D n, d, bi, bj, nbi, nbj, nmean ;
    double len ;

    d = i - b ;
    len = d.Magnitude() ;
    bi = d / len ;
    nbi[0] = bi[1] ;  nbi[1] = -bi[0] ;

    d = j - b ;
    len = d.Magnitude() ;
    bj = d / len ;
    nbj[0] = -bj[1] ;  nbj[1] = bj[0] ;

    nmean = 0.5 * (nbi + nbj) ;
    len = nmean.x()*nmean.x() + nmean.y()*nmean.y() ;

    if (len < 0.0001) {   // degenerate case
        n = (bi + bj) / 2.0 ;
        n = n.Normalize() ;
        double angle = Angle(b,j,i) ;
        if (angle < 0.0) n = -n ;
    } else {
        len = sqrt(len) ;
        n = nmean / len ;
    }
    return(n) ;
}


// -----------------------------------------------------------

void CArbMshReg2DNodeIterator::First()
{
    iter.First() ;
    if (ntype != UNSPECIFIED) {
        while (iter.More() && (iter.Entry()->type != ntype))
            iter.Next() ;
    }
}

void CArbMshReg2DNodeIterator::Next()
{
    iter.Next() ;
    if (ntype != UNSPECIFIED) {
        while (iter.More() && (iter.Entry()->type != ntype))
            iter.Next() ;
    }
}


// -----------------------------------------------------------

static inline int id(int i,int j,int m)
{
    return ((i)*(m))+(j) ;
}

void CArbMshMapRegion2D::GenerateMesh(ArbMshNode bdry_nodes[],
                                      int nu,int nv)
{
    ArbIntNode **pts = new ArbIntNode*[nu*nv] ;
    int i,j,k ;

    // set the node step increment and diagonal type

    int node_inc = (Order == LINEAR) ? 1 : 2 ;

    // Get the starting index (offset position) of each boundary
    // in the boundary array

    int in1 = 0 ;
    int in2 = in1 + nu - 1 ;
    int in3 = in2 + nv - 1 ;
    int in4 = in3 + nu - 1 ;

    // move the boundary nodes into the output node array

    for (i=0 ; i<(nu-1) ; ++i) {
        int ip = in1 + i ;
        int nid = NewNode(bdry_nodes[ip].coord[0],
                         bdry_nodes[ip].coord[1],
                         BOUNDARY,ARB_FIXED,
                         (node_inc == 1) || ((i%2) == 0)) ;
        pts[id(0,i,nu)] = pnode_table->Fetch(nid) ;
    }
    for (i=0 ; i<(nv-1) ; ++i) {
        int ip = in2 + i ;
        int nid = NewNode(bdry_nodes[ip].coord[0],
                         bdry_nodes[ip].coord[1],
                         BOUNDARY,ARB_FIXED,
                         (node_inc == 1) || ((i%2) == 0)) ;
        pts[id(i,nu-1,nu)] = pnode_table->Fetch(nid) ;
    }
    for (i=1 ; i<nu ; ++i) {
        int ip = in3 + (nu-i-1) ;
        int nid = NewNode(bdry_nodes[ip].coord[0],
                         bdry_nodes[ip].coord[1],
                         BOUNDARY,ARB_FIXED,
                         (node_inc == 1) || ((i%2) == 0)) ;
        pts[id(nv-1,i,nu)] = pnode_table->Fetch(nid) ;
    }
    for (i=1 ; i<nv ; ++i) {
        int ip = in4 + (nv-i-1) ;
        int nid = NewNode(bdry_nodes[ip].coord[0],
                         bdry_nodes[ip].coord[1],
                         BOUNDARY,ARB_FIXED,
                         (node_inc == 1) || ((i%2) == 0)) ;
        pts[id(0,i,nu)] = pnode_table->Fetch(nid) ;
    }

    // Generate interior nodes

    for (k=node_inc ; k<(nu-node_inc) ; k+=node_inc) {
        for (j=node_inc ; j<(nv-node_inc) ; k+=node_inc) {
            double x, y ;
            x = (pts[id(0,k,nu)]->coord[0]*(nv-j-1)/(nv-1) +
                 pts[id(nv-1,k,nu)]->coord[0]*j/(nv-1) +
                 pts[id(j,0,nu)]->coord[0]*(nu-k-1)/(nu-1) +
                 pts[id(j,nu-1,nu)]->coord[0]*k/(nu-1) -
                 pts[id(0,0,nu)]->coord[0]*(nu-k-1)/(nu-1)*(nv-j-1)/(nv-1) -
                 pts[id(nv-1,0,nu)]->coord[0]*(nu-k-1)/(nu-1)*j/(nv-1) -
                 pts[id(nv-1,nu-1,nu)]->coord[0]*k/(nu-1)*j/(nv-1) -
                 pts[id(0,nu-1,nu)]->coord[0]*k/(nu-1)*(nv-j-1)/(nv-1)) ;

            y = (pts[id(0,k,nu)]->coord[1]*(nv-j-1)/(nv-1) +
                 pts[id(nv-1,k,nu)]->coord[1]*j/(nv-1) +
                 pts[id(j,0,nu)]->coord[1]*(nu-k-1)/(nu-1) +
                 pts[id(j,nu-1,nu)]->coord[1]*k/(nu-1) -
                 pts[id(0,0,nu)]->coord[1]*(nu-k-1)/(nu-1)*(nv-j-1)/(nv-1) -
                 pts[id(nv-1,0,nu)]->coord[1]*(nu-k-1)/(nu-1)*j/(nv-1) -
                 pts[id(nv-1,nu-1,nu)]->coord[1]*k/(nu-1)*j/(nv-1) -
                 pts[id(0,nu-1,nu)]->coord[1]*k/(nu-1)*(nv-j-1)/(nv-1)) ;

            int nid = NewNode(x,y,BOUNDARY,ARB_FIXED,true) ;
            pts[id(j,k,nu)] = pnode_table->Fetch(nid) ;
        }
    }

    // Generate midside nodes

    if (node_inc == 2) {
        for (k=2 ; k<(nu-2) ; k+=2) {
            for (j=1 ; j<(nv-1) ; j+=2) {
                CArbCoord2D coord = 0.5 *
                    (pts[id(j-1,k,nu)]->coord +
                     pts[id(j+1,k,nu)]->coord) ;
                int nid = NewNode(coord.x(),coord.y(),
                                 BOUNDARY,ARB_FIXED,false) ;
                pts[id(j,k,nu)] = pnode_table->Fetch(nid) ;
            }
        }
        for (k=1 ; k<(nu-1) ; k+=2) {
            for (j=2 ; j<(nv-2) ; j+=2) {
                CArbCoord2D coord = 0.5 *
                    (pts[id(j,k-1,nu)]->coord +
                     pts[id(j,k+1,nu)]->coord) ;
                int nid = NewNode(coord.x(),coord.y(),
                                 BOUNDARY,ARB_FIXED,false) ;
                pts[id(j,k,nu)] = pnode_table->Fetch(nid) ;
            }
        }
    }

    // Start generating

    ArbMshElement2D elem ;

    if ((Order == LINEAR) && (Type == TRIANGLE)) {
        for (k=1 ; k<(nu+1-node_inc) ; k+=node_inc) {
            for (j=1 ; j<(nv+1-node_inc) ; j+=node_inc) {
                elem.elem_id = NewElemNum() ;
                elem.mat_id = MatId ;
                elem.num_nodes = 3 ;
                elem.nodes[0] = pts[id(j-1,k-1,nu)]->id ;
                elem.nodes[1] = pts[id(j  ,k  ,nu)]->id ;
                elem.nodes[2] = pts[id(j  ,k-1,nu)]->id ;
                pelem_table->Store(elem.elem_id,elem) ;
                elem.elem_id = NewElemNum() ;
                elem.mat_id = MatId ;
                elem.num_nodes = 3 ;
                elem.nodes[0] = pts[id(j-1,k-1,nu)]->id ;
                elem.nodes[1] = pts[id(j-1,k  ,nu)]->id ;
                elem.nodes[2] = pts[id(j  ,k  ,nu)]->id ;
                pelem_table->Store(elem.elem_id,elem) ;
            }
        }
    } else if ((Order == LINEAR) && (Type == QUADRILATERAL)) {
        for (k=1 ; k<(nu+1-node_inc) ; k+=node_inc) {
            for (j=1 ; j<(nv+1-node_inc) ; j+=node_inc) {
                elem.elem_id = NewElemNum() ;
                elem.mat_id = MatId ;
                elem.num_nodes = 4 ;
                elem.nodes[0] = pts[id(j-1,k-1,nu)]->id ;
                elem.nodes[1] = pts[id(j-1,k  ,nu)]->id ;
                elem.nodes[2] = pts[id(j  ,k  ,nu)]->id ;
                elem.nodes[3] = pts[id(j  ,k-1,nu)]->id ;
                pelem_table->Store(elem.elem_id,elem) ;
            }
        }
    } else if ((Order == QUADRATIC) && (Type == TRIANGLE)) {
        for (k=1 ; k<(nu+1-node_inc) ; k+=node_inc) {
            for (j=1 ; j<(nv+1-node_inc) ; j+=node_inc) {
                CArbCoord2D coord = 0.5 *
                    (pts[id(j-1,k-1,nu)]->coord +
                     pts[id(j+1,k+1,nu)]->coord) ;
                int nid = NewNode(coord.x(),coord.y(),
                                  BOUNDARY,ARB_FIXED,false) ;
                pts[id(j,k,nu)] = pnode_table->Fetch(nid) ;

                elem.elem_id = NewElemNum() ;
                elem.mat_id = MatId ;
                elem.num_nodes = 6 ;
                elem.nodes[0] = pts[id(j-1,k-1,nu)]->id ;
                elem.nodes[1] = pts[id(j+1,k+1,nu)]->id ;
                elem.nodes[2] = pts[id(j+1,k-1,nu)]->id ;
                elem.nodes[3] = pts[id(j  ,k  ,nu)]->id ;
                elem.nodes[4] = pts[id(j+1,k  ,nu)]->id ;
                elem.nodes[5] = pts[id(j  ,k-1,nu)]->id ;
                pelem_table->Store(elem.elem_id,elem) ;

                elem.elem_id = NewElemNum() ;
                elem.mat_id = MatId ;
                elem.num_nodes = 6 ;
                elem.nodes[0] = pts[id(j-1,k-1,nu)]->id ;
                elem.nodes[1] = pts[id(j-1,k+1,nu)]->id ;
                elem.nodes[2] = pts[id(j+1,k+1,nu)]->id ;
                elem.nodes[3] = pts[id(j-1,k  ,nu)]->id ;
                elem.nodes[4] = pts[id(j  ,k+1,nu)]->id ;
                elem.nodes[5] = pts[id(j  ,k  ,nu)]->id ;
                pelem_table->Store(elem.elem_id,elem) ;
            }
        }
    } else if ((Order == QUADRATIC) && (Type == QUADRILATERAL)) {
        for (k=1 ; k<(nu+1-node_inc) ; k+=node_inc) {
            for (j=1 ; j<(nv+1-node_inc) ; j+=node_inc) {
                elem.elem_id = NewElemNum() ;
                elem.mat_id = MatId ;
                elem.num_nodes = 8 ;
                elem.nodes[0] = pts[id(j-1,k-1,nu)]->id ;
                elem.nodes[1] = pts[id(j-1,k+1,nu)]->id ;
                elem.nodes[2] = pts[id(j+1,k+1,nu)]->id ;
                elem.nodes[3] = pts[id(j+1,k-1,nu)]->id ;
                elem.nodes[4] = pts[id(j-1,k  ,nu)]->id ;
                elem.nodes[5] = pts[id(j  ,k+1,nu)]->id ;
                elem.nodes[6] = pts[id(j+1,k  ,nu)]->id ;
                elem.nodes[7] = pts[id(j,  k-1,nu)]->id ;
                pelem_table->Store(elem.elem_id,elem) ;
            }
        }
    }
}


// -----------------------------------------------------------
void CArbMshRegion2D:: SetTriangularMesh (int n_node, double *Coords,
                                          int n_elem, int *Conn)
{
  int i;

  // nodes
  nnode = n_node;
  coords = new double[n_node*2];
  for (i = 0; i < n_node*2; i++)
    coords[i] = Coords[i];
  // elements
  nelem = n_elem;
  conn = new int[n_elem*4];
  for (i = 0; i < n_elem*4; i++)
    conn[i] = Conn[i];


#if 0

  // it inserts nodes
  int i, index = 0;

  ArbMshNode node;
  int num_bedges = pedge_table->NumEntries() ;
  for (i = num_bedges; i < n_node; i++)
  {
    node.id = i;
    node.coord[0] = coords[i*2+0];
    node.coord[1] = coords[i*2+1];
    node.coord[2] = 0.0;
    NewNode (coords[i*2+0], coords[i*2+0], INTERIOR, ARB_FLOATING,true) ;
    //AddNode (node, INTERIOR, ARB_FLOATING);
  }


  // adds elements
  for (index = 0, i = 0; i < n_elem; i++)
  {
    ArbMshElement2D elem;
    elem.elem_id = NewElemNum();
    elem.mat_id = MatId;
    elem.num_nodes = 3;
    if (conn[index] == 3)
    {
      elem.nodes[0] = conn[index+3];
      elem.nodes[1] = conn[index+2];
      elem.nodes[2] = conn[index+1];
    }
    else if (conn[index] == 6)
    {
      elem.nodes[0] = conn[index+5];
      elem.nodes[1] = conn[index+3];
      elem.nodes[2] = conn[index+1];
    }

    pelem_table->Store (elem.elem_id, elem);

    index = index + conn[index]+1;
  }

  int num_nodes = pnode_table->NumEntries() ;
  ArbIntNode **nodes = pnode_table->GetEntryList() ;
  int num_edges = pedge_table->NumEntries ( );
  ArbIntEdge **edges = pedge_table->GetKeyList() ;


  int *node_map;
  CArbHashTableIterator<int,ArbIntNode> iter(pnode_table) ;

  // Find the max node ID
  int max_id = 0 ;
  for (iter.First() ; iter.More() ; ++iter)
  {
    if (iter.Entry()->id > max_id)
      max_id = iter.Entry()->id ;
  }

  node_map = new int[max_id+1];
  for (i=0,iter.First() ; iter.More() ; ++i,++iter)
    node_map[iter.Entry()->id]=i;


  CArbHashTable <int, CArbCoord2D> bdry_start ;
  CArbHashTable <int, CArbCoord2D> bdry_stop ;
  MateTable = new CArbHashTable <int, CArbSmallSet <int, 2> > ;

  GetBdryInfo(num_nodes,nodes,num_edges,edges,
              node_map,bdry_start,bdry_stop,
              *MateTable) ;
#endif

}






