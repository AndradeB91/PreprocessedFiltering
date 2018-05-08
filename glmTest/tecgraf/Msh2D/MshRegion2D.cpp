//
// MshRegion2D Class definition
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
//   $Revision: 1.58 $  $Date: 2004/11/19 17:57:01 $  $Author: wash $
//
// -------------------------------------

//#define DEBUG_DUMP

#include <cstdio>
#include <cmath>
#include <cassert>

#include "Heap.hpp"
#include "QuadTree.hpp"

#include "MshRegion2D.hpp"
#include "MshSmooth2D.hpp"
#include "KDTree2D.hpp"

//#include "ArbMshTopo2D.hpp"
//#include "ArbFeasableRegion.hpp"

using FTools::QuadTree ;

namespace Msh2D {

#ifdef MEMDEBUG
#include "MemDbg.hpp"
//#define new new(__FILE__,__LINE__)
#endif

#define DEFAULT_REFINE_BDRY_FACT 1.1
//#define DEFAULT_NEAR_BDRY_FACT 0.6
#define DEFAULT_NEAR_BDRY_FACT 0.2
#define DEFAULT_MINIMUM_SHAPE_MEASURE 0.3

#define PI     3.141592654 
#define TWO_PI 6.283185307

// ------------------------------------------------------------

int MshRegion2D::IntEdgePriority::Compare(
    const IntEdge &e1,
    const IntEdge &e2)
{
    if ((e1.length-e2.length) < e1.tol) return(0) ;
    if (e1.length > e2.length) return(1) ;
    if (e1.length < e2.length) return(-1) ;
    return(0) ;
}

int MshRegion2D::IntEdgeSet::Compare(
    const IntEdge &e1,
    const IntEdge &e2) const
{
    if (e1.node_id[0] > e2.node_id[0]) {
        return(1) ;
    } else if (e1.node_id[0] < e2.node_id[0]) {
        return(-1) ;
    } else {
        if (e1.node_id[1] > e2.node_id[1]) {
            return(1) ;
        } else if (e1.node_id[1] < e2.node_id[1]) {
            return(-1) ;
        }
    }
    return(0) ;
}

//static int DictHashIndex(int i) { return(static_cast<int>(i)) ; }

/* ++ ----------------------------------------------------------
**
**    MshRegion2D - constructor method 
**
**      MshRegion2D(
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
**      Description: This is a constructor method for a MshRegion2D 
**          object. 
**
**
** -- */

MshRegion2D::MshRegion2D(const MshOrder iorder,
                         const int istart_id,
                         const int imat_id) :
    Order(iorder),StartId(istart_id),MatId(imat_id)
{
    NodeGeneration = true ;
    BoundaryChecks = true ;
    DoSmoothNodes = true ;
    WinslowSmoothing = true ;
    QuadTreeDebug = false ;
//    QuadTreeDebug = true ;
    StartId = 0 ;
    MaxId = 0 ;
    StartElemId = 0 ;
    MaxNumElem = 0 ;
    RefineBoundaryFactor = DEFAULT_REFINE_BDRY_FACT ;
    NearBoundaryFactor = DEFAULT_NEAR_BDRY_FACT ;
    MinimumShapeMeasure = DEFAULT_MINIMUM_SHAPE_MEASURE ;
    pnode_table = new Dict<int,IntNode> ;
    pedge_table = new IntEdgeSet() ;
    pelem_table = new Dict<int,MshElement2D> ;
    MateTable = 0 ;
#ifdef DOING_QUADS
//    CoordCache = new List<Vec2D>() ;
//    CoordCacheFlags = new List<bool>() ;
//    BdryVtxCache = new List<int>() ;
#endif
    DebugDisplayFlags = 0 ;
    CheckAllPoints = false ;
    dbg_file = 0 ;
}


MshRegion2D::MshRegion2D(const MshRegion2D &other)
{
    NodeGeneration = other.NodeGeneration ;
    BoundaryChecks = other.BoundaryChecks ;
    DoSmoothNodes = other.DoSmoothNodes ;
    WinslowSmoothing = other.WinslowSmoothing ;
    QuadTreeDebug = other.QuadTreeDebug ;
    StartId = other.StartId ;
    MaxId = other.MaxId ;
    StartElemId = other.StartElemId ;
    MaxNumElem = other.MaxNumElem ;
    RefineBoundaryFactor = other.RefineBoundaryFactor ;
    NearBoundaryFactor = other.NearBoundaryFactor ;
    MinimumShapeMeasure = other.MinimumShapeMeasure ;

    pnode_table = new Dict<int,IntNode> ;
    *pnode_table = *other.pnode_table ;
    pedge_table = new IntEdgeSet() ;
    OrderedSet<IntEdge>::SetIterator iter = other.pedge_table->Iterator() ;
    for (iter.First() ; iter.More() ; ++iter) {
        pedge_table->Insert(*iter) ;
    }
    //*pedge_table = *other.pedge_table ;
    pelem_table = new Dict<int,MshElement2D> ;
    *pelem_table = *other.pelem_table ;
    MateTable = new Dict<int,SmallSet<int,2> > ;
    *MateTable = *other.MateTable ;

    DebugDisplayFlags = other.DebugDisplayFlags ;
    CheckAllPoints = other.CheckAllPoints ;
    dbg_file = other.dbg_file ;
}

MshRegion2D MshRegion2D::operator =
     (const MshRegion2D &other)
{
    NodeGeneration = other.NodeGeneration ;
    BoundaryChecks = other.BoundaryChecks ;
    DoSmoothNodes = other.DoSmoothNodes ;
    WinslowSmoothing = other.WinslowSmoothing ;
    QuadTreeDebug = other.QuadTreeDebug ;
    StartId = other.StartId ;
    MaxId = other.MaxId ;
    StartElemId = other.StartElemId ;
    MaxNumElem = other.MaxNumElem ;
    RefineBoundaryFactor = other.RefineBoundaryFactor ;
    NearBoundaryFactor = other.NearBoundaryFactor ;
    MinimumShapeMeasure = other.MinimumShapeMeasure ;

    pnode_table = new Dict<int,IntNode> ;
    *pnode_table = *other.pnode_table ;
    pedge_table = new IntEdgeSet() ;
    OrderedSet<IntEdge>::SetIterator iter = other.pedge_table->Iterator() ;
    for (iter.First() ; iter.More() ; ++iter) {
        pedge_table->Insert(*iter) ;
    }
    //*pedge_table = *other.pedge_table ;
    pelem_table = new Dict<int,MshElement2D> ;
    *pelem_table = *other.pelem_table ;
    MateTable = new Dict<int,SmallSet<int,2> > ;
    *MateTable = *other.MateTable ;

    DebugDisplayFlags = other.DebugDisplayFlags ;
    CheckAllPoints = other.CheckAllPoints ;
    dbg_file = other.dbg_file ;

    return(*this) ;
}


/* ++ ----------------------------------------------------------
**
**    MshRegion2D - destructor 
**
**      ~MshRegion2D()
**
**      Description: This is a destructor for a MshRegion2D object. 
**
**
** -- */

MshRegion2D::~MshRegion2D()
{
    delete pnode_table ;
    delete pedge_table ;
    delete pelem_table ;
    if (MateTable != 0) delete MateTable ;
#ifdef DOING_QUADS
//    delete CoordCache ;
//    delete CoordCacheFlags ;
//    delete BdryVtxCache ;
#endif
}



// %(MshRegion2D::AddNode-void-|-ArbMshNode-const|&-ArbMshNodeType-const|-ArbMshNodeMotion-const|)
/* ++ ----------------------------------------------------------
**
**    AddNode - add a node to the region 
**
**      int AddNode(
**              const ArbMshNode       &inode,
**              const ArbMshNodeType   itype = BOUNDARY,
**              const ArbMshNodeMotion imotion = MSH_FLOATING)
**
**        inode   - (in)  node description 
**        itype   - (in)  BOUNDARY or INTERIOR 
**        imotion - (in)  MSH_FIXED or MSH_FLOATING 
**
**      Description: This method adds a node to the boundary of of the 
**          current region. 
**
**      Returns:
**          MSH_DUPLICATE_NODE_ID - duplicate node id
**          MSH_NORMAL_STATUS     - normal return status
**
** -- */

void MshRegion2D::AddNode(
    int id,
    Vec2D coord,
    bool retained_flag,
    MshNodeType itype,
    MshNodeMotion imotion)
{
    if (pnode_table->HasKey(id)) throw DuplicateNodeError() ;

    IntNode int_node ;
    int_node.id = id ;
    int_node.coord = coord ;
    int_node.type = itype ;
    int_node.motion = (itype == BOUNDARY) ? MSH_FIXED : imotion ;

    // If this is an interior node, set it as a corner node.
    // The boundary nodes are set to corner nodes in add egdge
    // if they are part of a boundary edge.

    int_node.corner = (itype == BOUNDARY) ? false : true ;
    pnode_table->Store(id,int_node) ;

    if (id > MaxId) MaxId = id ;

    if (dbg_file) DebugAddNodeEcho(id,coord,retained_flag,itype,imotion) ;
}


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
**          MSH_BAD_NODE_ID    - invalid edge end node id
**          MSH_DUPLICATE_EDGE - duplicate edge
**          MSH_NORMAL_STATUS  - normal return status
**
**
** -- */

void MshRegion2D::AddEdge(int nd0,int nd1,int nd2)
{
    if (!pnode_table->HasKey(nd0)) throw NodeIdError() ;
    if (!pnode_table->HasKey(nd1)) throw NodeIdError() ;
    if ((nd2 >= 0) && !pnode_table->HasKey(nd2)) throw NodeIdError() ;

    IntNode *pTmpNode ;

    // first check to see if the nodes are valid

    pTmpNode = pnode_table->Get(nd0) ; 
    pTmpNode->corner = true ;

    pTmpNode = pnode_table->Get(nd1) ; 
    pTmpNode->corner = true ;

    // now build a internal edge structure, create a key
    // by adding the first two node numbers then try to
    // add them to the set of edges.

    IntEdge tmp_edge(nd0,nd1,nd2,BOUNDARY) ;
    if (pedge_table->Insert(tmp_edge) == false) throw DuplicateEdgeError() ;

    if (dbg_file) DebugAddEdgeEcho(nd0,nd1,nd2) ;
}




void MshRegion2D::AddRefinementPoint(const Vec2D& xy,double size)
{
    RefinePts.Append(xy) ;
    RefineSize.Append(size) ;
    if (dbg_file) DebugAddRefineEcho(xy,size) ;
}


// %(MshRegion2D::NumBoundaryNodes-int-|^const)
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

int MshRegion2D::NumBoundaryNodes() const
{
    Dict<int,IntNode>::DictIterator iter(pnode_table) ;
    int num = 0 ;

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.Entry().type == BOUNDARY) ++num ;
    }
    return(num) ;
}




// %(MshRegion2D::NumInternalNodes-int-|^const)
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

int MshRegion2D::NumInternalNodes() const
{
    Dict<int,IntNode>::DictIterator iter(pnode_table) ;
    int num = 0 ;

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.Entry().type == INTERIOR) ++num ;
    }
    return(num) ;
}


// %(MshRegion2D::NumBoundaryEdges-int-|^const)
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

int MshRegion2D::NumBoundaryEdges() const
{
    int num = 0 ;

    OrderedSet<IntEdge>::SetIterator iter = pedge_table->Iterator() ;
    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.Entry().type == BOUNDARY) ++num ;
    }
    return(num) ;
}


MshRegion2D::EdgeIterator MshRegion2D::GetBoundaryEdges() const
{
    List<UPair<int> > bound_entries ;
    OrderedSet<IntEdge>::SetIterator iter = pedge_table->Iterator() ;
    for (iter.First() ; iter.More() ; ++iter) {
        const IntEdge& edge = iter.Entry() ; 
        if (edge.type == BOUNDARY) {
            bound_entries.Append(UPair<int>(edge.node_id[0],
                                            edge.node_id[1])) ;
        }
    }
    return EdgeIterator(bound_entries) ;
}

// -----------------------------------------------------------

// This is a bad random number generator.  It is not suitable if
// you realy want a random number.  It is good enough for our
// puposes because we only want to use it add a little node to
// node coordinate in order to break symmetries and make the mesher
// behave deterministly on different systems

class BadRand {

    public:

        BadRand() {
            Seed = 287 ;
            M = 65536 ;
            A = 1664525 ;
            C = 1013904223 ;
        }

        ~BadRand() {}

        void Reset() {
            Seed = 287 ;
        }

        int RandInt() {
            int tmp = A * Seed + C  ;
            Seed = tmp % M ;
            return Seed ;
        }

        double Rand() {
            return double(RandInt()) / M ;
        }

    private:

        int M,A,C,Seed ;
} ;

static BadRand rand ;

/* ------------------------------------------------------------
    GenerateIntNodes - callback function to generate nodes
                          at the center of each quadtree cell.

    Input:
        void *p_data - callback data which is realy a pointer
                       to the calling MshRegion2D object.
        double origin_x - x coordinate of the center of this cell.
        double origin_y - y coordinate of the center of this cell.
        double half - half size of this cell.
        bool is_root - true if this is the root node.
        bool is_leaf - true if this is a leaf node.
*/

void GenerateIntNodes(void *p_data,
                double origin_x,double origin_y,double half,
                bool is_root,bool is_leaf)
{
    MshRegion2D *p_region = (MshRegion2D *)p_data ;

    if (is_leaf) {

        // if this is a leaf node check to see if a node
        // generated at this point is inside the region
        // and that it is not too close to any boundary edge

        OrderedSet<IntEdge>::SetIterator iter =
            p_region->pedge_table->Iterator() ;
        bool insert = true ;
        int n_cross = 0 ;

        Vec2D b(origin_x,origin_y) ;

        for (iter.First() ; iter.More() ; ++iter) {
            IntNode *p_nd0, *p_nd1 ;
            double len, rr, dist_sqr ;

            p_nd0 = p_region->pnode_table->Get(iter.Entry().node_id[0]) ;
            p_nd1 = p_region->pnode_table->Get(iter.Entry().node_id[1]) ;

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
//if (origin_y > 0.0 && origin_y < 10.0) {
//   printf("v %g %g 0 1 0 0\n",origin_x,origin_y) ;
//}
                break ;
            }

            // check to see if a scan line from minus infinity
            // to this point crosses the edge

            if (p_region->ScanCross(b,p_nd0->coord,p_nd1->coord))
                ++n_cross ;
        }

        if (p_region->QuadTreeDebug) {
            if (insert) printf("v %g %g 0\n",origin_x,origin_y) ;
        }

        if (insert && ((n_cross % 2) == 1)) {

            // add a little random "salt" to the node coordinate
            // to break any symmetries and make the mesher work
            // deterministicly on all systems

            double scale = 0.001 * half ;
            double x = origin_x + scale * (rand.Rand() - 0.5) ;
            double y = origin_y + scale * (rand.Rand() - 0.5) ;

            (void)p_region->NewNode(x,y,INTERIOR,MSH_FLOATING,true) ;
        }
    }


    if (p_region->QuadTreeDebug) {
        double x0,y0,x1,y1 ;

        if (is_root) {
            x0 = origin_x - half ; y0 = origin_y - half ;
            x1 = origin_x + half ; y1 = origin_y - half ;
            printf("e %g %g 0 %g %g 0\n",x0,y0,x1,y1) ;

            x0 = x1 ; y0 = y1 ;
            x1 = origin_x + half ; y1 = origin_y + half ;
            printf("e %g %g 0 %g %g 0\n",x0,y0,x1,y1) ;

            x0 = x1 ; y0 = y1 ;
            x1 = origin_x - half ; y1 = origin_y + half ;
            printf("e %g %g 0 %g %g 0\n",x0,y0,x1,y1) ;

            x0 = x1 ; y0 = y1 ;
            x1 = origin_x - half ; y1 = origin_y - half ;
            printf("e %g %g 0 %g %g 0\n",x0,y0,x1,y1) ;
        }

        if (!is_leaf) {
            x0 = origin_x - half ; y0 = origin_y ;
            x1 = origin_x + half ; y1 = origin_y ;
            printf("e %g %g 0 %g %g 0\n",x0,y0,x1,y1) ;

            x0 = origin_x ; y0 = origin_y + half ;
            x1 = origin_x ; y1 = origin_y - half ;
            printf("e %g %g 0 %g %g 0\n",x0,y0,x1,y1) ;
        }
    }
}

#if 0
static void DisplayQTree(void *p_data,
                double origin_x,double origin_y,double half,
                bool is_root,bool is_leaf)
{

        double x0,y0,x1,y1 ;

        if (is_root) {
            x0 = origin_x - half ; y0 = origin_y - half ;
            x1 = origin_x + half ; y1 = origin_y - half ;
            printf("e %g %g 0 %g %g 0\n",x0,y0,x1,y1) ;

            x0 = x1 ; y0 = y1 ;
            x1 = origin_x + half ; y1 = origin_y + half ;
            printf("e %g %g 0 %g %g 0\n",x0,y0,x1,y1) ;

            x0 = x1 ; y0 = y1 ;
            x1 = origin_x - half ; y1 = origin_y + half ;
            printf("e %g %g 0 %g %g 0\n",x0,y0,x1,y1) ;

            x0 = x1 ; y0 = y1 ;
            x1 = origin_x - half ; y1 = origin_y - half ;
            printf("e %g %g 0 %g %g 0\n",x0,y0,x1,y1) ;
        }

        if (!is_leaf) {
            x0 = origin_x - half ; y0 = origin_y ;
            x1 = origin_x + half ; y1 = origin_y ;
            printf("e %g %g 0 %g %g 0\n",x0,y0,x1,y1) ;

            x0 = origin_x ; y0 = origin_y + half ;
            x1 = origin_x ; y1 = origin_y - half ;
            printf("e %g %g 0 %g %g 0\n",x0,y0,x1,y1) ;
        }
}
#endif



// %(MshRegion2D::DebugDumpRegion-void-|)
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

void MshRegion2D::DebugDumpRegion()
{
    FILE *fd = fopen("debug.rgn","w") ;
    if (fd != 0) {
        int num_total = pnode_table->Len() ;
        fprintf(fd,"%d\n",num_total) ;
        Dict<int,IntNode>::DictIterator niter = pnode_table->Iterator() ;
        for (niter.First() ; niter.More() ; ++niter) {
            fprintf(fd,"%d %18.10g %18.10g\n",niter.Entry().id,
                    niter.Entry().coord[0],niter.Entry().coord[1]) ;
        }

        int num_bound = NumBoundaryNodes() ;

        OrderedSet<IntEdge>::SetIterator iter = pedge_table->Iterator() ;
        fprintf(fd,"%d\n",num_bound) ;
        for (iter.First() ; iter.More() ; ++iter) {
            if (iter.Entry().type == BOUNDARY) {
                if (Order == QUADRATIC) {
                    fprintf(fd,"%d %d %d\n",
                        iter.EntryPtr()->node_id[0],
                        iter.EntryPtr()->node_id[1],
                        iter.EntryPtr()->node_id[2]) ;
                } else {
                    fprintf(fd,"%d %d\n",
                        iter.EntryPtr()->node_id[0],
                        iter.EntryPtr()->node_id[1]) ;
                }
            }
        }
        fclose(fd) ;
    }
}


// %(MshRegion2D::DebugDumpMesh-void-|)
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

void MshRegion2D::DebugDumpMesh()
{
    FILE *fd = fopen("debug.cds","w") ;
    if (fd != 0) {
        int i, num_node = pnode_table->Len() ;
        fprintf(fd,"%d\n",num_node) ;

        Dict<int,IntNode>::DictIterator pn_iter = pnode_table->Iterator() ;
        for (pn_iter.First() ; pn_iter.More() ; ++pn_iter) {
            fprintf(fd,"%d %g %g\n",pn_iter.Entry().id,
                                    pn_iter.Entry().coord[0],
                                    pn_iter.Entry().coord[1]) ;
        }
        fclose(fd) ;
    }

    fd = fopen("debug.msh","w") ;
    if (fd != 0) {
        int i, j, num_elem = pelem_table->Len() ;
        fprintf(fd,"%d\n",num_elem) ;

        Dict<int,MshElement2D>::DictIterator pe_iter = pelem_table->Iterator() ;
        for (pe_iter.First() ; pe_iter.More() ; ++pe_iter) {
            fprintf(fd,"%d %d",pe_iter.Entry().elem_id,
                               pe_iter.Entry().num_nodes) ;
            for (j=0 ; j<pe_iter.Entry().num_nodes ; ++j)
                fprintf(fd," %d",pe_iter.Entry().nodes[j]) ;
            fprintf(fd,"\n") ;

        }
        fclose(fd) ;
    }
}

#endif




// %(MshRegion2D::GenerateMesh-void-|)
/* ++ ----------------------------------------------------------
**
**    GenerateMesh - generate a mesh 
**
**      void GenerateMesh()
**
**      Description: This method generates a mesh for the region. 
**
**      Returns:
**          MSH_ILLEGAL_BDRY   - the boundary specified is illegal
**          MSH_NORMAL_STATUS  - normal return status
**
** -- */

void MshRegion2D::GenerateMesh()
{
    if (dbg_file) {
        fprintf(dbg_file,"DOIT: \n") ;
        fflush(dbg_file) ;
    }


#ifdef DEBUG_DUMP
    DebugDumpRegion() ;
#endif
    if (DebugDisplayFlags & RmshBoundary)
        DisplayBoundary("Remesh Boundary") ;

    // if necessary, check the validity of the boundary

    if (BoundaryChecks) CheckBoundary() ;

    // generate internal nodes if required

    if (StartId == 0) StartId = MaxId+1 ;
    StartIdSave = StartId ; 
    StartElemIdSave = StartElemId ;
    rand.Reset() ;
    if (NodeGeneration && pedge_table->Len() > 3) GenerateNodes() ;

    // do the boundary contraction

    if (BoundaryContraction() != 0) throw CannotMeshError() ;

#ifdef DEBUG_DUMP
    DebugDumpMesh() ;
#endif

    if (DebugDisplayFlags & RmshTriBeforeSmooth)
        DisplayMesh("Trimesh Before Smoothing") ;

    // smooth the internal nodes

    if (DoSmoothNodes) {

        SmoothNodes() ;

        // Sometimes it seems that winslow smoothing can do bad
        // things to triangular meshes with poor element size
        // distribution.  Therefore, we check for bad elements
        // and if found fall back to laplace smoothing

        if (WinslowSmoothing && !CheckValidTriangles()) {
            MshSmooth2D smooth(pnode_table,pelem_table) ;
            smooth.SmoothNodesLaplace() ;
        }
    }

    if (DebugDisplayFlags & RmshTriAfterSmooth)
        DisplayMesh("Trimesh After Smoothing") ;

//     // if quadrilateral elements are required, generate them
// 
//     if (Type == QUADRILATERAL) {
//         Dict<int,MshElement2D> elem_save = *(pelem_table) ;
//         Dict<int,IntNode> node_save = *(pnode_table) ;
// 
//         if (GenerateQuads()) {
//             RenumberNodes() ;
//         } else {
//             // quad meshing failed so revert to
//             // the saved triangular mesh
//             *(pelem_table) = elem_save ;
//             *(pnode_table) = node_save ;
//         }
//     }
}


bool MshRegion2D::CheckValidTriangles()
{
    Dict<int,MshElement2D>::DictIterator iter(pelem_table) ;
    for (iter.First() ; iter.More() ; ++iter) {
        IntNode *nd0, *nd1, *nd2 ;
        nd0 = pnode_table->Get(iter.Entry().nodes[0]) ;
        nd1 = pnode_table->Get(iter.Entry().nodes[1]) ;
        nd2 = pnode_table->Get(iter.Entry().nodes[2]) ;
        double sm = TriShapeMeasure(nd0->coord,nd1->coord,nd2->coord) ;
        if (sm <= 0.0) return(false) ;
    }
    return(true) ;
}


// %(MshRegion2D::CheckBoundary-void-|)
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
**          MSH_ILLEGAL_BDRY   - the boundary specified is illegal
**          MSH_NORMAL_STATUS  - normal return status
**
** -- */

void MshRegion2D::CheckBoundary()
{
    double tol = 1e-12 ;
//    int j ;

    // get a local list of all the edges

    int num_edges = pedge_table->Len() ;
    OrderedSet<IntEdge>::SetIterator iter = pedge_table->Iterator() ;

    // check to make sure that we have at least 3 edges and nodes

    if ((num_edges < 3) || (pnode_table->Len() < 3)) {
        throw IllegalBoundaryError() ;
    }

    // Check for "loop" edges (same starting and ending nodes)

    for (iter.First() ; iter.More() ; ++iter) {
        if ((*iter).node_id[0] == (*iter).node_id[1]) {
            throw IllegalBoundaryError() ;
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

    Dict<int,int> *pedge ;
    pedge = new Dict<int,int>(true) ;
    bool have_outside = false ;

    for (iter.First() ; iter.More() ; ++iter) {
        pedge->Store(iter.Entry().node_id[0],iter.Entry().node_id[1]) ;
    }

    while (pedge->Len() > 0) {
        int *pEntry ;
        int start, nd0, nd1 ;
        double angle = 0 ;

        // get the starting node, note that the angle between
        // the start node and the first and last edges is zero.

        Dict<int,int>::DictIterator itmp = pedge->Iterator() ;
        itmp.First() ;
        start = itmp.Key() ;
        nd0 = itmp.Entry() ;
        Vec2D base = (pnode_table->Get(start))->coord ;
        pedge->DelEntry(start,nd0) ;

        // look through all nodes on the loop.

        if ((pEntry = pedge->Get(nd0)) == 0) {
            delete pedge ;
            throw IllegalBoundaryError() ;
        }
        nd1 = *pEntry ;
        pedge->DelEntry(nd0,nd1) ;
        while (nd1 != start) {
            Vec2D i = (pnode_table->Get(nd0))->coord ;
            Vec2D j = (pnode_table->Get(nd1))->coord ;
            angle += Angle(base,i,j) ;
            nd0 = nd1 ;
            if ((pEntry = pedge->Get(nd0)) == 0) {
                delete pedge ;
                throw IllegalBoundaryError() ;
            }
            nd1 = *pEntry ;
            pedge->DelEntry(nd0,nd1) ;
        }

        // check for outside loop

        if (angle > tol) {
            if (have_outside) {
                delete pedge ;
                throw IllegalBoundaryError() ;
            } else
                have_outside = true ;
        }

// This check is commented out because it will cause us to fail
// for internal cracks.
// else if (angle == 0.0) {
//            return(MSH_ILLEGAL_BOUNDARY) ;
//        }
    }

//std::cerr << "Past loops compute: " << have_outside << std::endl ;

    if (!have_outside) {
        delete pedge ;
        throw IllegalBoundaryError() ;
    }

    delete pedge ;
}




// %(MshRegion2D::DisplayBoundary-void-|)
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

void MshRegion2D::DisplayBoundary(const char *label)
{
    OrderedSet<IntEdge>::SetIterator iter = pedge_table->Iterator() ;

    // check to make sure that none of the boundary
    // edges cross each other

    for (iter.First() ; iter.More() ; ++iter) {

        Vec2D i1, i2 ;

        i1 = (pnode_table->Get(iter.Entry().node_id[0]))->coord ;
        i2 = (pnode_table->Get(iter.Entry().node_id[1]))->coord ;

        printf("e %g %g 0 %g %g 0\n",i1.x(),i1.y(),i2.x(),i2.y()) ;
        printf("t %g %g 0 %d\n",i1.x(),i1.y(),iter.Entry().node_id[0]) ;
        printf("# edge: %d %d\n",iter.Entry().node_id[0],
                                 iter.Entry().node_id[1]) ;
    }

//    Dict<int,IntNode>::DictIterator iter(pnode_table) ;
//    for (iter.First() ; iter.More() ; ++iter) {
//        IntNode *node = iter.Entry() ;
//        printf("t %g %g %d\n",node->coord[0],node->coord[1],node->id) ;
//    }

    printf("a %s\n",label) ;
    printf("f\n") ;
    fflush(stdout) ;
}




// %(MshRegion2D::GenerateNodes-void-|)
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

void MshRegion2D::GenerateNodes()
{
    int i ;

    // loop through all the nodes and find the max and min
    // x and y values

    Dict<int,IntNode>::DictIterator iter(pnode_table) ;
    double minx,maxx,miny,maxy ;

    minx = maxx = iter.Entry().coord[0] ;
    miny = maxy = iter.Entry().coord[1] ;
    for (iter++ ; iter.More() ; ++iter) {
        if (iter.Entry().coord[0] < minx) minx = iter.Entry().coord[0] ;
        if (iter.Entry().coord[0] > maxx) maxx = iter.Entry().coord[0] ;
        if (iter.Entry().coord[1] < miny) miny = iter.Entry().coord[1] ;
        if (iter.Entry().coord[1] > maxy) maxy = iter.Entry().coord[1] ;
    }

    // Create a quadtree that is centered on the boundary and
    // large enough the contain all the nodes plus a little fudge
    // to avoid numerical tolerance troubles.

    double dx = maxx - minx ;
    double dy = maxy - miny ;
    double size = (dx > dy) ? dx : dy ;

    QuadTree *pQTree = new QuadTree((maxx+minx)/2,
                                (maxy+miny)/2,size*1.01) ;

    // Locally refine the quad tree at the center of each
    // edge.  Refine to size proportional to the edge size.

    OrderedSet<IntEdge>::SetIterator siter = pedge_table->Iterator() ;
    double max_edge = 0.0 ;
    double max_cell_size = 0.0 ;

    for (siter.First() ; siter.More() ; ++siter) {
        IntNode *p_nd0, *p_nd1 ;
        double dx, dy, cx, cy, len ;

        p_nd0 = pnode_table->Get(siter.Entry().node_id[0]) ;
        p_nd1 = pnode_table->Get(siter.Entry().node_id[1]) ;

        dx = p_nd1->coord[0] - p_nd0->coord[0] ;
        dy = p_nd1->coord[1] - p_nd0->coord[1] ;
        cx = (p_nd1->coord[0] + p_nd0->coord[0]) / 2 ;
        cy = (p_nd1->coord[1] + p_nd0->coord[1]) / 2 ;

        len = sqrt(dx*dx + dy*dy) ;
        if (len > max_edge) max_edge = len ;

        pQTree->RefineToSize(cx,cy,len*RefineBoundaryFactor) ;
            // DisplayBoundary("Remesh Boundary") ;
            // if (QuadTreeDebug) printf("# Start Quadtree: %d\n",i) ;
            // pQTree->VisitLevels((void *)this,DisplayQTree) ;
            // if (QuadTreeDebug) printf("a after refine %d\n",i) ;
            // if (QuadTreeDebug) printf("# End Quadtree\nn\n") ;
    }

    // Now any additionally define refinement points

    for (i=0 ; i<RefinePts.Len() ; ++i) {
        pQTree->RefineToSize(RefinePts[i].x(),RefinePts[i].y(),
                             RefineSize[i]) ;
    }

    // look to find maximum boundary cell size

    for (siter.First() ; siter.More() ; ++siter) {
        IntNode *p_nd0, *p_nd1 ;
        double cx, cy, csize ;

        p_nd0 = pnode_table->Get((*siter).node_id[0]) ;
        p_nd1 = pnode_table->Get((*siter).node_id[1]) ;

        cx = (p_nd1->coord[0] + p_nd0->coord[0]) / 2 ;
        cy = (p_nd1->coord[1] + p_nd0->coord[1]) / 2 ;

        csize = pQTree->ContainingCellSize(cx,cy) ;
        if (csize > max_cell_size) max_cell_size = csize ;
    }

        // QuadTreeDebug = true ;
        // if (QuadTreeDebug) printf("# Start Quadtree: %d\n",i) ;
        // DisplayBoundary("Remesh Boundary") ;
        // pQTree->VisitLevels((void *)this,DisplayQTree) ;
        // if (QuadTreeDebug) printf("a after refine\n") ;
        // if (QuadTreeDebug) printf("# End Quadtree\nn\n") ;

    // Now refine the quad tree so that no cell is larger than
    // the larges cell on the boundary.

    pQTree->UniformRefine(1.99*max_cell_size,
                          minx,maxx,miny,maxy) ;
    // pQTree->UniformRefine(max_edge*RefineBoundaryFactor,
    //                       minx,maxx,miny,maxy) ;

        //   if (QuadTreeDebug) printf("# Start Quadtree\n") ;
        //   DisplayBoundary("Remesh Boundary") ;
        //   pQTree->VisitLevels((void *)this,DisplayQTree) ;
        //   if (QuadTreeDebug) printf("a after uniform\n") ;
        //   if (QuadTreeDebug) printf("# End Quadtree\nn\n") ;

    // Now refine the quad tree so that no two adjacent cells
    // differ in level of refinement by more than one.

    pQTree->RefineOneLevelDiff() ;

        //   if (QuadTreeDebug) printf("# Start Quadtree\n") ;
        //   pQTree->VisitLevels((void *)this,DisplayQTree) ;
        //   if (QuadTreeDebug) printf("a after one level\n") ;
        //   if (QuadTreeDebug) printf("# End Quadtree\nn\n") ;
        //   fflush(stdout) ;

    // Visit all the cells in the quad tree and generate nodes
    // in the center of all the leaf nodes.

    if (QuadTreeDebug) printf("# Start Quadtree\n") ;

    if (StartId == 0) StartId = MaxId + 1 ;
    pQTree->VisitLevels((void *)this,GenerateIntNodes) ;

    if (QuadTreeDebug) printf("# End Quadtree\n") ;

    delete pQTree ;
}




// %(MshRegion2D::DisplayMesh-void-|)
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

void MshRegion2D::DisplayMesh(const char *label,FILE *dfd)
{
    FILE *fd = (dfd == 0) ? stdout : dfd ;
    Dict<int,MshElement2D>::DictIterator eiter(pelem_table) ;

#if 1
//    bool do_labels = false ;
//    bool do_labels = true ;

    for (eiter.First() ; eiter.More() ; ++eiter) {
        MshElement2D& elem = eiter.Entry() ;
        double x[8], y[8] ;
        int id[8] ;
        int nn ;
        if (elem.num_nodes == 8)
            nn = 4 ;
        else if (elem.num_nodes == 6)
            nn = 3 ;
        else
            nn = elem.num_nodes ;
        for (int j=0 ; j<elem.num_nodes ; ++j) {
            IntNode *node = pnode_table->Get(elem.nodes[j]) ;
            x[j] = node->coord[0] ;
            y[j] = node->coord[1] ;
            id[j] = node->id ;
        }

/*
double xc = 0.0 ;
double yc = 0.0 ;
for (int jj=0 ; jj<nn ; ++jj) {
    xc += x[jj] ;
    yc += y[jj] ;
}
xc /= double(nn) ;
yc /= double(nn) ;

for (int jj=0 ; jj<nn ; ++jj) {
    x[jj] = xc + 0.75*(x[jj]-xc) ;
    y[jj] = yc + 0.75*(y[jj]-yc) ;
}
*/

        for (int jj=0 ; jj<nn ; ++jj) {
            int kk = (jj+1) % nn ;
            fprintf(fd,"e %g %g 0 %g %g 0\n",x[jj],y[jj],x[kk],y[kk]) ;
        }

        if (elem.elem_id > 0) {
            for (int jj=0 ; jj<elem.num_nodes ; ++jj) {
                fprintf(fd,"t %g %g 0 %d\n",x[jj],y[jj],id[jj]) ;
            }
        }

        // if (do_labels) {
        //     for (int jj=0 ; jj<elem->num_nodes ; ++jj) {
        //         fprintf(fd,"t %g %g 0 %d\n",x[jj],y[jj],id[jj]) ;
        //     }
        // }
    }
    fprintf(fd,"a %s\n",label) ;
    fprintf(fd,"f\n") ;
    fflush(fd) ;

#if 0
    for (eiter.First() ; eiter.More() ; ++eiter) {
        MshElement2D& elem = eiter.Entry() ;
        double x[8], y[8] ;
        int id[8] ;
        int nn ;
        if (elem.num_nodes == 8)
            nn = 4 ;
        else if (elem.num_nodes == 6)
            nn = 3 ;
        else
            nn = elem.num_nodes ;
        for (int j=0 ; j<elem.num_nodes ; ++j) {
            IntNode *node = pnode_table->Get(elem.nodes[j]) ;
            x[j] = node->coord[0] ;
            y[j] = node->coord[1] ;
            id[j] = node->id ;
        }

/*
double xc = 0.0 ;
double yc = 0.0 ;
for (int jj=0 ; jj<nn ; ++jj) {
    xc += x[jj] ;
    yc += y[jj] ;
}
xc /= double(nn) ;
yc /= double(nn) ;

for (int jj=0 ; jj<nn ; ++jj) {
    x[jj] = xc + 0.75*(x[jj]-xc) ;
    y[jj] = yc + 0.75*(y[jj]-yc) ;
}
*/

        for (int jj=0 ; jj<nn ; ++jj) {
            int kk = (jj+1) % nn ;
            fprintf(fd,"e %g %g 0 %g %g 0\n",x[jj],y[jj],x[kk],y[kk]) ;
        }

        if (elem->elem_id > 0) {
            for (int jj=0 ; jj<elem->num_nodes ; ++jj) {
                fprintf(fd,"t %g %g 0 %d\n",x[jj],y[jj],id[jj]) ;
            }
        }

        // if (do_labels) {
        //     for (int jj=0 ; jj<elem->num_nodes ; ++jj) {
        //         fprintf(fd,"t %g %g 0 %d\n",x[jj],y[jj],id[jj]) ;
        //     }
        // }
    }
    fprintf(fd,"a %s\n",label) ;
    fprintf(fd,"f\n") ;
    fflush(fd) ;
#endif

#endif

#if 0
    for (eiter.First() ; eiter.More() ; ++eiter) {
        MshElement2D& elem = eiter.Entry() ;
        double x[8], y[8] ;
        int id[8] ;
        int nn ;
        if (elem.num_nodes == 8)
            nn = 4 ;
        else if (elem.num_nodes == 6)
            nn = 3 ;
        else
            nn = elem.num_nodes ;

        printf("e %d\n",elem.num_nodes) ;

        for (int j=0 ; j<elem.num_nodes ; ++j) {
            IntNode *node = pnode_table->Get(elem.nodes[j]) ;
            x[j] = node->coord[0] ;
            y[j] = node->coord[1] ;
            id[j] = node->id ;
            printf(" %d",elem.nodes[j]) ;
        }
        printf("\n") ;

        for (int jj=0 ; jj<nn ; ++jj) {
            int kk = (jj+1) % nn ;
            printf("%g %g\n",x[jj],y[jj]) ;
        }
    }
    printf("a %s\n",label) ;
    printf("f\n") ;
    fflush(stdout) ;
#endif
}


#define MATE_TOL 0.005

void MshRegion2D::GetBdryInfo(
                   int num_nodes,
                   List<IntNode*>& nodes,
                   OrderedSet<IntEdge>::SetIterator& eiter,
                   int *node_map,
                   Dict<int,Vec2D> &bdry_start,
                   Dict<int,Vec2D> &bdry_stop,
                   Dict<int,SmallSet<int,2> > &mate_table)
{
    int i,j ;
    double *edge_len = new double[num_nodes] ;
    double *edge_min = new double[num_nodes] ;
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

    for (eiter.First() ; eiter.More() ; ++eiter) {
        int nd0 = node_map[eiter.Entry().node_id[0]] ;
        int nd1 = node_map[eiter.Entry().node_id[1]] ;

        Vec2D v0 =
           (nodes[nd1]->coord-nodes[nd0]->coord).Normalize() ;

        bdry_start.Store(nodes[nd0]->id,v0) ;
        bdry_stop.Store(nodes[nd1]->id,-v0) ;

        double elen = (nodes[nd0]->coord - nodes[nd1]->coord).Magnitude() ;
        edge_len[nd0] += elen ;
        edge_len[nd1] += elen ;
        if ((edge_num[nd0] == 0) || (elen < edge_min[nd0])) edge_min[nd0] = elen ;
        if ((edge_num[nd1] == 0) || (elen < edge_min[nd1])) edge_min[nd1] = elen ;
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
    IntNode **used = new IntNode*[num_nodes] ;
    for (i=0 ; i<num_nodes ; ++i) {
        if (nodes[i]->type == BOUNDARY) {
            used[num_used] = nodes[i] ;
            num_used++ ;
        }
    }

    KDTree2D<IntNode> kd_tree(num_used,used) ;

    // now loop through the edges and look for mate nodes

    List<IntNode*> *range_points = 0 ;
    for (i=0 ; i<num_nodes ; ++i) {

        // if this is a boundary node (average length set)
        // do a range search looking for the node

        if (edge_len[i] > 0.0) {
            IntNode ll,ur ;
            double delt = edge_len[i] * 0.1 ;
            ll.coord[0] = nodes[i]->coord[0] - delt ;
            ll.coord[1] = nodes[i]->coord[1] - delt ;
            ur.coord[0] = nodes[i]->coord[0] + delt ;
            ur.coord[1] = nodes[i]->coord[1] + delt ;
            KDTree2D<IntNode>::Rectangle rect(ll,ur,delt) ;
            range_points = kd_tree.RangeQuery(rect) ;

            // if one or more nodes were found in the search then
            // check to see if it is within the tolerance.  If so
            // then store the information in the mate table

            for (j=0 ; j<range_points->Len() ; ++j) {
                IntNode *p = (*range_points)[j] ;
                if (p->id == nodes[i]->id) continue ;
                double dist = (nodes[i]->coord - p->coord).Magnitude() ;
//                if ((dist < edge_len[i]*MATE_TOL) &&
//                    (dist < edge_min[i]*0.99)) {         // NEW
                if (dist < edge_min[i]*MATE_TOL) {         // NEW
                    SmallSet<int,2> *mate =
                            mate_table.Get(nodes[i]->id) ;
                    if (mate == 0) {
                        SmallSet<int,2> ss ;
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
    delete [] edge_min ;
    delete [] edge_num ;
    delete [] used ;

//     Dict<int,SmallSet<int,2>::DictIterator > iter(&mate_table) ;
//     for (iter.First() ; iter.More() ; ++iter) {
//         int key = iter.Key() ;
//         SmallSet<int,2> *ss = iter.Entry() ;
//         fprintf(stderr,"mate: %d %d %d\n",key,ss->NumElements(),
//                                           ss->Element(0)) ;
//     }
}


#define BETWEEN_TOL -0.000001

static bool Between(Vec2D &b,Vec2D &i,Vec2D &j)
{
    double cross = CrossProd(i,j) ;
    if (fabs(cross) < 1e-14 && (i*j > 0.0)) {
        return(true) ;
    } else if (cross > 0.0) {
        if ((CrossProd(i,b) > BETWEEN_TOL) &&
            (CrossProd(b,j) > BETWEEN_TOL)) return(true) ;
    } else {
        //if (CrossProduct(i,b) >= 0) return(true) ;
        if (CrossProd(i,b) > BETWEEN_TOL) return(true) ;
        if (CrossProd(b,j) > BETWEEN_TOL) return(true) ;
    }
    return(false) ;
}


//#define EDGE_FACTOR 0.01
#define EDGE_FACTOR 0.001

bool MshRegion2D::CheckCross(
                   IntEdgePriority *pedge_heap,
                   IntEdge **bdry,IntEdge *entry,
                   List<IntNode*>& nodes,
                   int *node_map,
                   IntNode *p,double blensqr,
                   IntNode *p_nd0,IntNode *p_nd1,
                   Vec2D &enorm,Vec2D &emid,
                   Dict<int,Vec2D> &bdry_start,
                   Dict<int,Vec2D> &bdry_stop,
                   int *lside_0,int *lside_1,
                   bool *nd_0_flg,bool *nd_1_flg,
                   Dict<int,SmallSet<int,2> > &mate_table)
{
    // get the mates list for the edge

    SmallSet<int,2> *mates0 = mate_table.Get(entry->node_id[0]) ;
    SmallSet<int,2> *mates1 = mate_table.Get(entry->node_id[1]) ;

    // first check for a mate

    if (((mates0 != 0) && mates0->HasElement(p->id)) ||
        ((mates1 != 0) && mates1->HasElement(p->id))) return(true) ; 

    // now check for a cross

    for (int j=0 ; j<pedge_heap->Len() ; ++j) {
        if ((bdry[j]->node_id[0] != entry->node_id[0]) ||
            (bdry[j]->node_id[1] != entry->node_id[1])) {

            // check for any mates

            if ((mates0 != 0) && 
                (mates0->HasElement(bdry[j]->node_id[0]) ||
                 mates0->HasElement(bdry[j]->node_id[1]))) continue ;          

            if ((mates1 != 0) && 
                (mates1->HasElement(bdry[j]->node_id[0]) ||
                 mates1->HasElement(bdry[j]->node_id[1]))) continue ;          

            IntNode *p_nd2, *p_nd3 ;
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
                            Vec2D *start = bdry_start.Get(p_nd2->id) ;
                            if (start != 0) {
                                Vec2D *stop = bdry_stop.Get(p_nd2->id) ;
                                Vec2D v0 = 
                                   (p_nd0->coord-p->coord).Normalize() ;
                                Vec2D v1 = 
                                   (p_nd1->coord-p->coord).Normalize() ;
                                if (Between(v0,*start,*stop) &&
                                    Between(v1,*start,*stop)) return(true) ;
                            } else {
                                return(true) ;
                            }
                        } else {
                            Vec2D *start = bdry_start.Get(p_nd3->id) ;
                            if (start != 0) {
                                Vec2D *stop = bdry_stop.Get(p_nd3->id) ;
                                Vec2D v0 = 
                                   (p_nd0->coord-p->coord).Normalize() ;
                                Vec2D v1 = 
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
                            Vec2D *start = bdry_start.Get(p_nd2->id) ;
                            if (start != 0) {
                                Vec2D *stop = bdry_stop.Get(p_nd2->id) ;
                                Vec2D v0 = 
                                   (p_nd0->coord-p->coord).Normalize() ;
                                Vec2D v1 = 
                                   (p_nd1->coord-p->coord).Normalize() ;
                                if (Between(v0,*start,*stop) &&
                                    Between(v1,*start,*stop)) return(true) ;
                            } else {
                                return(true) ;
                            }
                        } else {
                            Vec2D *start = bdry_start.Get(p_nd3->id) ;
                            if (start != 0) {
                                Vec2D *stop = bdry_stop.Get(p_nd3->id) ;
                                Vec2D v0 = 
                                   (p_nd0->coord-p->coord).Normalize() ;
                                Vec2D v1 = 
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


// %(MshRegion2D::BoundaryContraction-void-|)
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

// static int ArbCmpEdge(const MshRegion2D::IntEdge &ed1,const MshRegion2D::IntEdge &ed2)
// {
//     if (ed1.length > ed2.length) return(1) ;
//     if (ed1.length < ed2.length) return(-1) ;
//     return(0) ;
// }

#define RANGE_FACTOR_1 5.0
#define RANGE_FACTOR_2 10.0
#define RANGE_FACTOR_3 20.0

// this is the tolerance used to determine if two nodes are at
// the same location.

#define ANGLE_TOL 0.0017
//#define ANGLE_TOL 0.002
#define POS_TOL 0.0075    // set for debugFlaw81.dmp
//#define POS_TOL 0.005
//#define POS_TOL 0.0025
//#define POS_TOL 0.001
//#define POS_TOL 0.01
#define MIN_ANGLE_TOL 0.01

#define TRIAL_FACTOR 10

int MshRegion2D::BoundaryContraction()
{
    //IntNode **nodes ;
    int        *node_map ;
    List<IntNode*> *range_points = 0 ;
    int valid, retry=0 ;// flags for checking elements
    int i ;

    int tnum = 0 ;

    // get a list of all the nodes

    int num_nodes = pnode_table->Len() ;
    //nodes = pnode_table->GetEntryList() ;
    List<IntNode*> nodes ;
        // for (i=0 ; i<num_nodes ; ++i) {
        //     fprintf(stderr,"%d %g %g\n",nodes[i]->id,
        //             nodes[i]->coord[0],nodes[i]->coord[1]) ;
        // }

    Dict<int,IntNode>::DictIterator iter = pnode_table->Iterator() ;

    // Find the max node ID

    int max_id = 0 ;
    for (iter.First() ; iter.More() ; ++iter) {
        nodes.Append(pnode_table->Get(iter.Key())) ;
        if (iter.Entry().id > max_id) max_id = iter.Entry().id ;
    }

    node_map = new int[max_id+1];
    for (i=0,iter.First() ; iter.More() ; ++i,++iter) {
      node_map[iter.Entry().id]=i;
    }

    // place all the boundary edges in a heap data
    // structure so that we can extract the smallest edge

    IntEdgePriority *pedge_heap = new IntEdgePriority ;

    OrderedSet<IntEdge>::SetIterator siter = pedge_table->Iterator() ;

    for (siter.First() ; siter.More() ; ++siter) {
        IntEdge *edge = siter.EntryPtr() ;
        IntNode *p_nd0, *p_nd1 ;
        double dx, dy ;

        p_nd0 = nodes[node_map[edge->node_id[0]]] ;
        p_nd1 = nodes[node_map[edge->node_id[1]]] ;

        // fprintf(stderr,"edge: %d %d\n",edges[i]->node_id[0],
        //                                edges[i]->node_id[1]) ;

        dx = p_nd1->coord[0] - p_nd0->coord[0] ;
        dy = p_nd1->coord[1] - p_nd0->coord[1] ;
        edge->length = sqrt(dx*dx + dy*dy) ;
        edge->tol = 1e-12 * edge->length ;

        
        pedge_heap->Insert(siter.Entry()) ;
    }

    // here we want to preprocess things for the case where we
    // have two surfaces that are coincident but not connected.
    // This may happen in the case of crack faces or slide line
    // type contacts.  For each node that is adjacent to a boundary
    // we create an entry in a hash table indexed by the node's
    // id.  The value of the entry is the negative of the normal
    // to the adjacent edge.

    Dict<int,Vec2D> bdry_start ;
    Dict<int,Vec2D> bdry_stop ;
    MateTable = new Dict<int,SmallSet<int,2> > ;

    GetBdryInfo(num_nodes,nodes,siter,
                node_map,bdry_start,bdry_stop,
                *MateTable) ;

    // get a list of all the nodes

    // build a TwoDTree 

    KDTree2D<IntNode> kd_tree(nodes) ;

    // loop until there are no more edges in the boundary

    while (pedge_heap->Len() > 0) {
        double max_angle = 0.0 ;
        bool flg_0 = false ;
        bool flg_1 = false ;
        int tri_nodes[6] ;
        int side_0 = -1 ;
        int side_1 = -1 ;
        IntEdge *entry, an_edge ;
        IntNode *p_nd0, *p_nd1 ;
        double dx, dy, blen, blensqr ;

        // get the shortest edge in the boundary

        entry = pedge_heap->GetMin() ;
        p_nd0 = nodes[node_map[entry->node_id[0]]] ;
        p_nd1 = nodes[node_map[entry->node_id[1]]] ;
        dx = p_nd1->coord[0] - p_nd0->coord[0] ;
        dy = p_nd1->coord[1] - p_nd0->coord[1] ;
        blen = sqrt(dx*dx + dy*dy) * POS_TOL ;
        blensqr = blen * blen ;
        Vec2D emid =
            Vec2D(0.5*(p_nd1->coord[0]+p_nd0->coord[0]),
                        0.5*(p_nd1->coord[1]+p_nd0->coord[1])) ;
        Vec2D enorm = Vec2D(-dy,dx) ;

// if ((p_nd0->id == 13506) && (p_nd1->id == 13504)) {
//     fprintf(stderr,"BoundaryContract: bug\n") ;
// }
    
        // define a range query rectangle around the edge
        // and get the points that fall in range
        
        if ((num_nodes < 10) || (CheckAllPoints)) {
            range_points = kd_tree.ReturnAll() ;
        } else if ( retry==0 ) {
            KDTree2D<IntNode>::Rectangle rect(*p_nd0, *p_nd1,
                           RANGE_FACTOR_1*entry->length);
            range_points = kd_tree.RangeQuery(rect) ;

            //fprintf(stderr,"Range points (1): %d\n",range_points->length());

            // if there are no points in range other than the
            // edge end points, increase the range size
            if ( range_points->Len() <= 3 ) {
                delete range_points;
                KDTree2D<IntNode>::Rectangle rect(*p_nd0, *p_nd1,
                               RANGE_FACTOR_2*entry->length);
                range_points = kd_tree.RangeQuery(rect) ;
                //fprintf(stderr,"Range points (1b): %d\n",range_points->length());
            }
        } else if ( retry==1 ) {
            KDTree2D<IntNode>::Rectangle rect(*p_nd0, *p_nd1,
                           RANGE_FACTOR_3*entry->length);
            range_points = kd_tree.RangeQuery(rect) ;
            //fprintf(stderr,"Range points (2): %d\n",range_points->length());
        } else if ( retry==2 ) {
            range_points = kd_tree.ReturnAll() ;
            //fprintf(stderr,"Range points (3): %d\n",range_points->length());
        }

        // now loop through all the nodes and find the node
        // that will make a triangle with the minimum included
        // angle

        valid = 0 ;
        int cur, num ;
        num = range_points->Len() ;
        for (cur=0 ; cur < num ; ++cur) {
          IntNode *p = (*range_points)[cur] ;

          if (p->corner) {

            // inserted to deal with nodes adjacent to an edge.

            Vec2D *start = bdry_start.Get(p->id) ;
            if (start != 0) {
                Vec2D *stop = bdry_stop.Get(p->id) ;
                Vec2D v0 = (p_nd0->coord-p->coord).Normalize() ;
                if (!Between(v0,*start,*stop)) continue ;
                Vec2D v1 = (p_nd1->coord-p->coord).Normalize() ;
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

            if ((p->id != entry->node_id[0]) &&
                (p->id != entry->node_id[1]) &&
                (len0 > blen) && (len1 > blen) &&
                (CrossProd(p->coord,
                           p_nd0->coord,p_nd1->coord) > 0.0)) {

                angle = Angle(p->coord,
                              p_nd0->coord,p_nd1->coord) ;

                if ((retry < 2) && (angle < 1.0e-6)) continue ;
//                if ((retry < 2) && (angle < 5.0e-5)) continue ;

                // check to see if this is the
                // biggest angle we've seen so far, check to
                // see if this is a valid triangle (i.e., it
                // does not cross the existing boundary

                if ((angle > max_angle) || (duplicate_node)) {
                    IntEdge **bdry = pedge_heap->GetEntryList() ;
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

//                    if ((!cross) && (angle > 1e-5)) {
//                    if ((!cross) && (angle > 1e-6)) {
                    if (!cross) {
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

        if ( !valid ) {
          //fprintf(stderr,"p_nd0: %d: %f, %f\n", 
          //        p_nd0->id, p_nd0->coord[0],p_nd0->coord[1]);
          //fprintf(stderr,"p_nd1: %d: %f, %f\n", 
          //        p_nd1->id, p_nd1->coord[0],p_nd1->coord[1]);
          if ( retry==2 ) {
              delete range_points;
              delete [] node_map ;
              delete pedge_heap ;

//fprintf(stderr,"Meshing failure\n") ;
//assert(0) ;
// 
// fprintf(stdout,"f\n") ;
// 
// Dict<int,MshElement2D>::DictIterator eiter(pelem_table) ;
// for (eiter.First() ; eiter.More() ; ++eiter) {
//     MshElement2D *elem = eiter.EntryPtr() ;
//     p_nd0 = nodes[node_map[elem->nodes[0]]] ;
//     p_nd1 = nodes[node_map[elem->nodes[1]]] ;
//     fprintf(stdout,"e %g %g 0 %g %g 0\n",
//             p_nd0->coord[0],p_nd0->coord[1],
//             p_nd1->coord[0],p_nd1->coord[1]) ;
//     p_nd0 = nodes[node_map[elem->nodes[1]]] ;
//     p_nd1 = nodes[node_map[elem->nodes[2]]] ;
//     fprintf(stdout,"e %g %g 0 %g %g 0\n",
//             p_nd0->coord[0],p_nd0->coord[1],
//             p_nd1->coord[0],p_nd1->coord[1]) ;
//     p_nd0 = nodes[node_map[elem->nodes[2]]] ;
//     p_nd1 = nodes[node_map[elem->nodes[0]]] ;
//     fprintf(stdout,"e %g %g 0 %g %g 0\n",
//             p_nd0->coord[0],p_nd0->coord[1],
//             p_nd1->coord[0],p_nd1->coord[1]) ;
// }
// exit(0) ;

return(1) ;

        } else {
            retry++;
            pedge_heap->Insert(*entry) ;
          }
        }
        
        else {
          retry = 0;
          MshElement2D elem ;

// if ((tri_nodes[0] == 181535) ||
//     (tri_nodes[1] == 181535) || 
//     (tri_nodes[2] == 181535)) {
//     std::cerr << tnum << ' ' << tri_nodes[0] << ' ' << tri_nodes[1] << ' ' << tri_nodes[2] << ' ' << std::endl ;
// }

// if ((tri_nodes[0] == 181529) ||
//     (tri_nodes[1] == 181529) || 
//     (tri_nodes[2] == 181529)) {
//     std::cerr << tnum << ' '  << tri_nodes[0] << ' ' << tri_nodes[1] << ' ' << tri_nodes[2] << ' ' << std::endl ;
// }

          elem.elem_id = NewElemNum() ;
          elem.mat_id = MatId ;
          if (Order == LINEAR) {
            elem.num_nodes = 3 ;
            for (int i=0 ; i<3 ; ++i) elem.nodes[i] = tri_nodes[i] ;
          } else {
            if (side_0 == -1) {
                double x, y ;
                p_nd0 = nodes[node_map[tri_nodes[0]]] ;
                p_nd1 = nodes[node_map[tri_nodes[2]]] ;
                x = (p_nd1->coord[0] + p_nd0->coord[0]) / 2 ;
                y = (p_nd1->coord[1] + p_nd0->coord[1]) / 2 ;
                tri_nodes[5] = NewNode(x,y,INTERIOR,MSH_FLOATING,false) ;
            } else {
                tri_nodes[5] = (int)side_0 ;
            }

            if (side_1 == -1) {
                double x, y ;
                p_nd0 = nodes[node_map[tri_nodes[1]]] ;
                p_nd1 = nodes[node_map[tri_nodes[2]]] ;
                x = (p_nd1->coord[0] + p_nd0->coord[0]) / 2 ;
                y = (p_nd1->coord[1] + p_nd0->coord[1]) / 2 ;
                tri_nodes[4] = NewNode(x,y,INTERIOR,MSH_FLOATING,false) ;
            } else {
                tri_nodes[4] = (int)side_1 ;
            }

            elem.num_nodes = 6 ;
            for (int i=0 ; i<6 ; ++i) elem.nodes[i] = tri_nodes[i] ;
          }

          pelem_table->Store(elem.elem_id,elem) ;
//          fprintf(stderr,"Elem : %d %d %d : %d\n",tri_nodes[0],
//                tri_nodes[1],tri_nodes[2],tnum) ;
          if (MaxNumElem && (pelem_table->Len() > MaxNumElem)) {
              delete range_points;
              delete [] node_map ;
              delete pedge_heap ;
              throw MaxElements() ;
          }

          if (!flg_0) {
            double dx, dy ;
            an_edge.node_id[0] = tri_nodes[0] ;
            an_edge.node_id[1] = tri_nodes[2] ;
            an_edge.node_id[2] = tri_nodes[5] ;
            p_nd0 = nodes[node_map[an_edge.node_id[0]]] ;
            p_nd1 = nodes[node_map[an_edge.node_id[1]]] ;
            dx = p_nd1->coord[0] - p_nd0->coord[0] ;
            dy = p_nd1->coord[1] - p_nd0->coord[1] ;
            an_edge.length = sqrt(dx*dx + dy*dy) ;
            an_edge.tol = 1e-12 * an_edge.length ;
// if (an_edge.node_id[0] == 181529 && an_edge.node_id[0] == 181086) {
//    std::cerr << "Found:" << std::endl ;
// }
            pedge_heap->Insert(an_edge) ;
          } else {
            an_edge.node_id[0] = tri_nodes[2] ;
            an_edge.node_id[1] = tri_nodes[0] ;
            pedge_heap->Remove(an_edge) ;
          }

          if (!flg_1) {
            double dx, dy ;
            an_edge.node_id[0] = tri_nodes[2] ;
            an_edge.node_id[1] = tri_nodes[1] ;
            an_edge.node_id[2] = tri_nodes[4] ;
            p_nd0 = nodes[node_map[an_edge.node_id[0]]] ;
            p_nd1 = nodes[node_map[an_edge.node_id[1]]] ;
            dx = p_nd1->coord[0] - p_nd0->coord[0] ;
            dy = p_nd1->coord[1] - p_nd0->coord[1] ;
            an_edge.length = sqrt(dx*dx + dy*dy) ;
            an_edge.tol = 1e-12 * an_edge.length ;
// if (an_edge.node_id[0] == 181529 && an_edge.node_id[0] == 181086) {
//    std::cerr << "Found:" << std::endl ;
// }
            pedge_heap->Insert(an_edge) ;
          } else {
            an_edge.node_id[0] = tri_nodes[1] ;
            an_edge.node_id[1] = tri_nodes[2] ;
            pedge_heap->Remove(an_edge) ;
          }
        }

        delete range_points;
// #ifdef DEBUG
//        DisplayMesh("elem") ;
//        fprintf(stderr,"Elem : %d %d %d : %d\n",tri_nodes[0],
//                tri_nodes[1],tri_nodes[2],tnum) ;
// #endif
        ++tnum ;
        if (MaxNumElem && (tnum*TRIAL_FACTOR > MaxNumElem)) {
              delete [] node_map ;
              delete pedge_heap ;
              throw MaxElements() ;
        }

//         if (tnum >= 150) {
//             DisplayMesh("elem") ;
//            // exit(1) ;
//         }

//         if (tnum > 10) {
//             char buff[40] ;
//             sprintf(buff,"tnum_%d",tnum) ;
//             DisplayMesh(buff) ;
//             fprintf(stderr,"tnum Display %d\n", tnum) ;
// //            exit(1) ;
//         }
    }

    delete [] node_map ;
    delete pedge_heap ;
    return(0) ;
} 


// %(MshRegion2D::NewNode-int-|-double-|-double-|-ArbMshNodeType-|-ArbMshNodeMotion-|-bool-|)
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
**        motion - (in)  MSH_FIXED or MSH_FLOATING 
**        corner - (in)  corner node flag 
**
**      Description: This method generates a new node and assigns it a 
**          unique node id. 
**
**      Return Value: the new node id 
**
**
** -- */

int MshRegion2D::NewNode(double x,double y,
                         MshNodeType type,
                         MshNodeMotion motion,
                         bool corner)
{
    int id ;
    IntNode int_node ;

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




// %(MshRegion2D::DuplicateNode-int-|-int-|)
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

int MshRegion2D::DuplicateNode(int id)
{
    IntNode *node = pnode_table->Get(id) ;

    if (node == 0) throw NodeIdError() ;

    return(NewNode(node->coord[0],node->coord[1],node->type,
                   node->motion,node->corner)) ;
}




// %(MshRegion2D::SmoothNodes-void-|)
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

void MshRegion2D::SmoothNodes()
{
   MshSmooth2D smooth(pnode_table,pelem_table) ;
   if (WinslowSmoothing)
       smooth.SmoothNodesWinslow() ;
   else
       smooth.SmoothNodesLaplace() ;
}




// %(MshRegion2D::GetMeshStatistics-ArbMshStats-|*^const)
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

MshRegion2D::MshStats MshRegion2D::GetMeshStatistics() const
{
    MshStats stats ; 

    stats.num_elements = pelem_table->Len() ;
    stats.num_nodes = pnode_table->Len() ;

    Dict<int,MshElement2D>::DictIterator eiter = pelem_table->Iterator() ;

    double sm, sm_sum = 0.0 ;
    stats.minimum_shape_measure = 100.0 ;

    int i ;
    for (i=0 ; i<11 ; ++i) stats.num_less_than[i] = 0 ;

    for (eiter.First() ; eiter.More() ; ++eiter) {
        if ((eiter.Entry().num_nodes == 3) ||
            (eiter.Entry().num_nodes == 6)) {

            // compute the element shape measure.  We are
            // using the mean edge ratio, which is defined
            // as sm = (4*sqrt(3)*Area)/(l01^2+l12^2+l20^2)
            // where lij is the length of the edge from node
            // i to j.

            IntNode *nd0, *nd1, *nd2 ;
            nd0 = pnode_table->Get(eiter.Entry().nodes[0]) ;
            nd1 = pnode_table->Get(eiter.Entry().nodes[1]) ;
            nd2 = pnode_table->Get(eiter.Entry().nodes[2]) ;
            sm = TriShapeMeasure(nd0->coord,nd1->coord,nd2->coord) ;
            if (sm <= 0.0) {
                fprintf(stderr,"Invalid Element: %d\n", eiter.Entry().elem_id) ;
            }
        } else {
            IntNode *nd0, *nd1, *nd2, *nd3 ;
            nd0 = pnode_table->Get(eiter.Entry().nodes[0]) ;
            nd1 = pnode_table->Get(eiter.Entry().nodes[1]) ;
            nd2 = pnode_table->Get(eiter.Entry().nodes[2]) ;
            nd3 = pnode_table->Get(eiter.Entry().nodes[3]) ;
            sm = QuadMetric(nd0->coord,nd1->coord,
                            nd2->coord,nd3->coord) ;
            if (sm <= 0.0) {
                fprintf(stderr,"Invalid Element: %d\n",eiter.Entry().elem_id) ;
            }
        }

        sm_sum += sm ;
        if (sm < stats.minimum_shape_measure)
            stats.minimum_shape_measure = sm ;
        for (int i = 1 ; i < 11 ; ++i) {
            double x = double(i) / 10.0 ;
            if (sm < x) stats.num_less_than[i]++ ;
        }
    }

    stats.mean_shape_measure = sm_sum / stats.num_elements ;

    return(stats) ;    
}




// %(MshRegion2D::Cross-bool-|^const-Vec2D-|-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    Cross - check to see if lines cross 
**
**      bool Cross(
**              Vec2D i1,
**              Vec2D i2,
**              Vec2D j1,
**              Vec2D j2) const
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

bool MshRegion2D::Cross(Vec2D &i1,
                        Vec2D &i2,
                        Vec2D &j1,
                        Vec2D &j2) const
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

    if ((i1[0] < j1[0]) && (i1[0] < j2[0]) &&
        (i2[0] < j1[0]) && (i2[0] < j2[0])) return false ;

    if ((i1[1] > j1[1]) && (i1[1] > j2[1]) &&
        (i2[1] > j1[1]) && (i2[1] > j2[1])) return false ;

    if ((i1[1] < j1[1]) && (i1[1] < j2[1]) &&
        (i2[1] < j1[1]) && (i2[1] < j2[1])) return false ;

    // now check to make sure that J1 and J2 are on opposite sides
    // of line I.

    Vec2D delt = i2 - i1 ;
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




// %(MshRegion2D::CrossProd-double-|^const-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    CrossProd - compute a cross product 
**
**      double CrossProd(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j) const
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

double MshRegion2D::CrossProd(Vec2D &b,
                              Vec2D &i,
                              Vec2D &j) const
{
    double cross ;

    cross = ((i[0] - b[0]) * (j[1] - b[1])) -
            ((i[1] - b[1]) * (j[0] - b[0])) ;
    return cross ;
}




// %(MshRegion2D::Angle-double-|^const-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    Angle - compute an angle 
**
**      double Angle(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j) const
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

double MshRegion2D::Angle(Vec2D &b,
                          Vec2D &i,
                          Vec2D &j) const
{
    Vec2D bi = i - b ;
    Vec2D bj = j - b ;

    double tmp0 = bi.x()*bi.x() + bi.y()*bi.y() ;
    if (tmp0 <= 0) return(0.0) ;

    double tmp1 = bj.x()*bj.x() + bj.y()*bj.y() ;
    if (tmp1 <= 0) return(0.0) ;

    double tmp = ((bi.x() * bj.x()) + (bi.y() * bj.y())) /
                 (sqrt(tmp0) * sqrt(tmp1)) ;

    if (tmp >  1.0) tmp =  1.0 ;
    if (tmp < -1.0) tmp = -1.0 ;

    double cross = CrossProd(b,i,j) ;

//std::cerr << "ang: " << tmp << ' ' << cross << std::endl ;

    if (cross == 0.0) {
        return((tmp > 0) ? 0.0 : acos(-1.0)) ;
    }
    return(acos(tmp) * (cross/fabs(cross))) ;
}




// %(MshRegion2D::Angle2Pi-double-|^const-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    Angle2Pi - compute an angle 
**
**      double Angle2Pi(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j) const
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

double MshRegion2D::Angle2Pi(Vec2D &b,
                             Vec2D &i,
                             Vec2D &j) const
{
    Vec2D bi = i - b ;
    Vec2D bj = j - b ;

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




// %(MshRegion2D::Area-double-|^const-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    Area - find the area of a triangle 
**
**      double Area(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j) const
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

double MshRegion2D::Area(Vec2D &b,
                         Vec2D &i,
                         Vec2D &j) const
{
    double area ;

    area = (i[0]-b[0])*(j[1]-b[1]) - (j[0]-b[0])*(i[1]-b[1]) ;

    return 0.5*area ;
}




// %(MshRegion2D::DistSqr-double-|^const-Vec2D-|-Vec2D-|-Vec2D-|-double-|*)
/* ++ ----------------------------------------------------------
**
**    DistSqr - compute the distance from a point to a line segment 
**
**      double DistSqr(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j,
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

double MshRegion2D::DistSqr(
        Vec2D &b,
        Vec2D &i,
        Vec2D &j,
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

    Vec2D base = j - i ;
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

    return(adist1 < adist2 ? adist1*adist1 : adist2*adist2) ;
}




// %(MshRegion2D::ScanCross-bool-|^const-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    ScanCross - check for a scan line crossing 
**
**      bool ScanCross(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j) const
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

bool MshRegion2D::ScanCross(Vec2D &b,
                            Vec2D &i,
                            Vec2D &j) const
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




// %(MshRegion2D::TriShapeMeasure-double-|^const-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    TriShapeMeasure - compute a shape measure 
**
**      double TriShapeMeasure(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j) const
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

double MshRegion2D::TriShapeMeasure(Vec2D &i,
                                    Vec2D &j,
                                    Vec2D &k) const
{
    // compute the element shape measure.  We are using the mean edge ratio,
    // which is defined as sm = (4*sqrt(3)*Area)/(l01^2+l12^2+l20^2)
    // where lij is the length of the edge from node i to j.

    static const double factor = 6.92820323 ;
    double len_sqr_sum = 0 ;

    Vec2D d = i - j ;
    len_sqr_sum += d.x()*d.x() + d.y()*d.y() ;
    d = j - k ;
    len_sqr_sum += d.x()*d.x() + d.y()*d.y() ;
    d = k - i ;
    len_sqr_sum += d.x()*d.x() + d.y()*d.y() ;

    return(factor * Area(i,j,k) / len_sqr_sum) ;
}




// %(MshRegion2D::BisectNorm-Vec2D-|^const-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    BisectNorm - find the bisector of an angle 
**
**      Vec2D BisectNorm(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j) const
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

Vec2D MshRegion2D::BisectNorm(Vec2D &b,
                              Vec2D &i,
                              Vec2D &j) const
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

    Vec2D n, d, bi, bj, nbi, nbj, nmean ;
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


/* ------------------------------------------------------------
    LenSqr - computes the square of the distance between to points
*/

inline double LenSqr(Vec2D i,Vec2D j)
{
    Vec2D delta = i - j ;
    return(delta.x()*delta.x() + delta.y()*delta.y()) ;
}



// %(MshQuad2D::TriMetric-double-|^const-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    TriMetric - compute a shape metric 
**
**      double TriMetric(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j) const
**
**        b - (in)  first vertex 
**        i - (in)  second vertex 
**        j - (in)  third vertex 
**
**      Description: This method computes a shape metric for a 
**          triangular element. A zero or negative metric indicates an 
**          invalid element. The maximum possible metric value is 1.0 
**          for an equilateral triangle. 
**
**      Return Value: the shape metric 
**
**
** -- */

double MshRegion2D::TriMetric(
        Vec2D i,
        Vec2D j,
        Vec2D k) const
{
    /* This is the assumed node configuration:
    
              * k
             / \
            /   \
         i *-----* j
    */

    double cross = CrossProd(k,i,j) ;
    if (fabs(cross) < 1.0e-12) cross = 0 ;
    double sum = LenSqr(i,j) + LenSqr(j,k) + LenSqr(k,i) ;
    static double factor = 3.464101615 ;  // 2 * sqrt(3)
    return(factor * cross / sum) ;
}




// %(MshQuad2D::QuadMetric-double-|^const-Vec2D-|-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    QuadMetric - compute a shape metric 
**
**      double QuadMetric(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j,
**              Vec2D l) const
**
**        b - (in)  first vertex 
**        i - (in)  second vertex 
**        j - (in)  third vertex 
**        l - (in)  forth vertex 
**
**      Description: This method computes a shape metric for a 
**          quadrilateral element. A zero or negative metric indicates 
**          an invalid element. The maximum possible metric value is 
**          1.0 for a square. 
**
**      Return Value: the shape metric 
**
**
** -- */

double MshRegion2D::QuadMetric(
        Vec2D i,
        Vec2D j,
        Vec2D k,
        Vec2D l) const
{
    /* This is the assumed node configuration:
    
         l *-----* k
           |     |
           |     |
         i *-----* j
    */

    // find the triangle metrics

    double a[4] ;
    a[0] = TriMetric(i,j,l) ;
    a[1] = TriMetric(j,k,i) ;
    a[2] = TriMetric(k,l,j) ;
    a[3] = TriMetric(l,i,k) ;

    // find the minimum

    int cur, min = 0 ;
    for (cur=1 ; cur<4 ; ++cur) {
        if (a[cur] < a[min]) min = cur ;
    }

    // figure out the negative value

    double negval = 0 ;
    for (cur=0 ; cur<4 ; ++cur) {
        if (a[cur] < 0) ++negval ;
    }

    if (negval >= 2) --negval ;

    // check for minimum angles

    static double min_angle = 0.10471975 ;  // 6 degrees
    if (negval == 0) {
       if (Angle(i,j,l) < min_angle) negval = 1 ;
       if (Angle(j,k,i) < min_angle) negval = 1 ;
       if (Angle(k,l,j) < min_angle) negval = 1 ;
       if (Angle(l,i,k) < min_angle) negval = 1 ;
    }

    a[min] /= 0.866025 ;
    return(a[min] - double(negval)) ;
}

double MshRegion2D::QuadMetricAreaOnly(
        Vec2D i,
        Vec2D j,
        Vec2D k,
        Vec2D l) const
{
    /* This is the assumed node configuration:
    
         l *-----* k
           |     |
           |     |
         i *-----* j
    */

    // find the triangle metrics

    double a[4] ;
    a[0] = TriMetric(i,j,l) ;
    a[1] = TriMetric(j,k,i) ;
    a[2] = TriMetric(k,l,j) ;
    a[3] = TriMetric(l,i,k) ;

    // find the minimum

    int cur, min = 0 ;
    for (cur=1 ; cur<4 ; ++cur) {
        if (a[cur] < a[min]) min = cur ;
    }

    // figure out the negative value

    double negval = 0 ;
    for (cur=0 ; cur<4 ; ++cur) {
        if (a[cur] < 0) ++negval ;
    }

    if ((negval == 1) && (a[min] > -0.05)) {
        negval = 0 ;
        a[min] = fabs(a[min]) ;
    }

    a[min] /= 0.866025 ;
    return(a[min] - double(negval)) ;
}

/* ++ ----------------------------------------------------------
**
**    FlatnessMetric - Compute a measure of the "flatness" of a quad
**
**      double QuadMetric(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j,
**              Vec2D l) const
**
**        b - (in)  first vertex 
**        i - (in)  second vertex 
**        j - (in)  third vertex 
**        l - (in)  forth vertex 
**
**      Description: This method computes a measure of the "flatness"
**          of a quad that will vary from 1.0 for a square to 0.0 for
**          flat quad with no area.  
**      
**      Return Value: the flatness metric 
**
**
** -- */

double MshRegion2D::FlatnessMetric(
    Vec2D &i,
    Vec2D &j,
    Vec2D &k,
    Vec2D &l) const
{
    /* This is the assumed node configuration:
    
         l *-----* k
           |     |
           |     |
         i *-----* j
    */

    // first consider the triangles ijl and klj, with
    // ij and lk being bases and compute the heights
    // and the ratios between hight and the base

    double r0 = 0.0 ;
    double b0 = (j-i).Magnitude() ;
    double b1 = (l-k).Magnitude() ;
    if (b0 != 0.0) {
        double h = 2.0*Area(i,j,l)/b0 ;
        r0 = h / b0 ;
    }
    if (b1 != 0.0) {
        double h = 2.0*Area(k,l,j)/b1 ;
        double r = h / b1 ;
        if (r > r0) r0 = r ;
    }

    // now consider the triangles lik and jki, with
    // li and jk being bases and compute the heights
    // and the ratios between hight and the base

    double r1 = 0.0 ;
    b0 = (j-i).Magnitude() ;
    b1 = (l-k).Magnitude() ;
    if (b0 != 0.0) {
        double h = 2.0*Area(l,i,k)/b0 ;
        r1 = h / b0 ;
    }
    if (b1 != 0.0) {
        double h = 2.0*Area(j,k,i)/b1 ;
        double r = h / b1 ;
        if (r > r1) r1 = r ;
    }

    return((r0<r1) ? r0 : r1) ;
}

// -----------------------------------------------------------

void MshRegion2D::NodeIterator::First()
{
    iter.First() ;
    if (ntype != UNSPECIFIED) {
        while (iter.More() && (iter.Entry().type != ntype))
            iter.Next() ;
    }
}

void MshRegion2D::NodeIterator::Next()
{
    iter.Next() ;
    if (ntype != UNSPECIFIED) {
        while (iter.More() && (iter.Entry().type != ntype))
            iter.Next() ;
    }
}

// -------------------------------------------------------------


void MshRegion2D::DebugAddNodeEcho(
    int id,
    Vec2D coord,
    bool retained_flag,
    MshNodeType itype,
    MshNodeMotion imotion)
{
    fprintf(dbg_file,"NODE: %d %23.16f %23.16f %d %d %d\n",
            id,coord[0],coord[1],(int)retained_flag,(int)itype,(int)imotion) ;
}

void MshRegion2D::DebugAddEdgeEcho(int nd0,int nd1,int nd2)
{
    fprintf(dbg_file,"EDGE: %d %d %d\n",nd0,nd1,nd2) ;
}

void MshRegion2D::DebugAddRefineEcho(const Vec2D& xy,double size)
{
    fprintf(dbg_file,"REFINE: %23.16f %23.16f %23.16f\n",
            xy[0],xy[1],size) ;
}


} // namespace
