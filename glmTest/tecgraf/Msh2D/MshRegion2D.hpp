//
// MshRegion2D Class header file
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
//   $Revision: 1.40 $  $Date: 2004/06/23 14:49:03 $  $Author: wash $
//

#ifndef MshRegion2D_h
#define MshRegion2D_h

#include <cstdlib>
#include <cstdio>
#include <stdexcept>

#include "Dict.hpp"
#include "Geom2DMixIn.hpp"
#include "Heap.hpp"
#include "List.hpp"
#include "Pair.hpp"
#include "Queue.hpp"
#include "OrderedSet.hpp"
#include "SmallSet.hpp"
#include "Vec2D.hpp"

#include "Msh2DTypes.hpp"
#include "MshTopo2D.hpp"
#include "MshEdgeList.hpp"

using FTools::Dict ;
using FTools::List ;
using FTools::Queue ;
using FTools::OrderedSet ;
using FTools::SmallSet ;
using FTools::UPair ;
using FTools::Vec2D ;


namespace Msh2D {

#define DOING_QUADS 0

//
// Object descritpion
//

class MshRegion2D : public FTools::Geom2DMixIn {

    public:

        class NodeIterator ;
        class ElemIterator ;
        class EdgeIterator ;

        class DuplicateNodeError : public std::runtime_error {
            public:
                DuplicateNodeError() : runtime_error("duplicate node id specified") {}
        } ;

        class DuplicateEdgeError : public std::runtime_error {
            public:
                DuplicateEdgeError() : runtime_error("duplicate edge specified") {}
        } ;

        class CannotMeshError : public std::runtime_error {
            public:
                CannotMeshError() : runtime_error("cannot mesh the region") {}
        } ;

        class CleanupMeshError : public std::runtime_error {
            public:
                CleanupMeshError() : runtime_error("bad mesh in the cleanup phase") {}
        } ;

        class IllegalBoundaryError : public std::runtime_error {
            public:
                IllegalBoundaryError() : runtime_error("illegal region boundary") {}
        } ;

        class MaxElements : public std::runtime_error {
            public:
                MaxElements() : runtime_error("max allowed elements exceeded") {}
        } ;

        class QuadConversionError : public std::runtime_error {
            public:
                QuadConversionError() : runtime_error("error convertng to quads") {}
        } ;

        class NodeIdError : public std::runtime_error {
            public:
                NodeIdError() : runtime_error("illegal node id") {}
        } ;

//         struct MshNode {
//             int id ;               // unique node id
//             Vec2D coord ;          // nodal coordinates
//             bool retained_flag ;   // retain this node if on a boundary
//         } ;
// 
//         struct MshEdge {
//             int node_id[3] ;       // id's of nodes defining edge
//         } ;

        struct MshStats {
            int num_elements ;
            int num_nodes ;
            double mean_shape_measure ;
            double minimum_shape_measure ;
            int num_less_than[11] ;
        } ;

        struct EdgeKey {
            int id_0 ;
            int id_1 ;

            EdgeKey() : id_0(0), id_1(0) {} ;
            EdgeKey(const int id0,const int id1) :
                id_0(id0), id_1(id1) {} ;

            int operator == (const EdgeKey &op) {
                return((id_0 == op.id_0) && (id_1 == op.id_1)) ;
            }
        } ;

        void DebugDumpRegion() ;

    // constructors and destructors

        MshRegion2D(const MshOrder iorder = LINEAR,
                    const int istart_id = 0,
                    const int imat_id = 0) ;

        MshRegion2D(const MshRegion2D &other) ;
        MshRegion2D operator = (const MshRegion2D &other) ;

        virtual ~MshRegion2D() ;

    // set and query configuration

        void SetOrder(const MshOrder iorder) { Order = iorder ; }
//        void SetElemType(const MshElemType itype) { Type = itype ; }
        void SetStartNodeID(const int is_id) { StartId = is_id ; }
        void SetStartElemID(const int is_id) { StartElemId = is_id ; }
        void SetMaterialID(const int mat_id) { MatId = mat_id ; }
        void SetNodeGeneration(bool flag = true) { NodeGeneration = flag ; }
        void SetBoundaryValidityChecks(bool flag = true) { BoundaryChecks = flag ; }
        void SetNodeSmoothing(bool flag = true) { DoSmoothNodes = flag ; }
        void SetRefineBoundaryFactor(double factor) { RefineBoundaryFactor = factor ; }
        void SetNearBoundaryFactor(double factor) { NearBoundaryFactor = factor ; }
        void SetQuadTreeDebug(bool flag = true) { QuadTreeDebug = flag ; }
        void SetWinslowSmoothing(bool flag = true) { WinslowSmoothing = flag ; }
        void SetMaxNumElements(int max_elem) { MaxNumElem = max_elem ; } ;
        void SetCheckAllPoints() { CheckAllPoints = true ; } ;

        MshOrder GetOrder() const { return(Order) ; }
//        MshElemType GetElemType() const { return(Type) ; }
        int GetStartNodeID() const { return(StartId) ; }

    // routines to add nodes and edges that define the boundary of
    // of the region

        void AddNode(
            int id,
            Vec2D coord,
            bool retained_flag,
            MshNodeType itype = BOUNDARY,
            MshNodeMotion imotion = MSH_FLOATING) ;

        void AddEdge(int nd0,int nd1,int nd2 = -1) ;

        void AddRefinementPoint(const Vec2D& xy,double size) ;

    // routines to check the boundary and generate a mesh

        virtual void GenerateMesh() ;
        virtual void CheckBoundary() ;

    // routines to query info about the generated mesh

        int NumBoundaryNodes() const ;
        int NumInternalNodes() const ;
        int NumNodes() const { return(pnode_table->Len()) ; }
        int NumBoundaryEdges() const ;
        int NumElements() const { return(pelem_table->Len()) ; }
 
        NodeIterator GetBoundaryNodes() const {
            return NodeIterator(*this,BOUNDARY) ;
        }
//         ArbMshNode *GetInternalNodes() const ;
//         ArbMshNode *GetNodes() const ;
        EdgeIterator GetBoundaryEdges() const ;
//         ArbMshElement2D *GetElements() const ;
// 
        NodeIterator GetNodeIterator() const { return NodeIterator(*this,UNSPECIFIED) ; }
        ElemIterator GetElemIterator() const { return ElemIterator(*this) ; }

        MshStats GetMeshStatistics() const ;

#ifdef DEBUG_DUMP
        void DebugDumpMesh() ;
#endif
        void DisplayMesh(const char *label,FILE *dfd=0) ;
        void DisplayBoundary(const char *label) ;
        void DisplayTopo(const char *label,MshTopo2D *topo) ;

#ifdef DOING_QUADS

        int NewElemNum() { return(StartElemId++) ; } ;
        int DuplicateNode(int id) ;

        struct SeamElemCache {
            int elem ;
            int num_nodes ;
            int nodes[4] ;
            SeamElemCache() : elem(0), num_nodes(0)
                { nodes[0]=nodes[1]=nodes[2]=nodes[3] = 0 ; } ;
            SeamElemCache(const int ielem,
                          const int inum_nodes,
                          const int *inodes) :
                elem(ielem),num_nodes(inum_nodes)
                { for (int i=0 ; i<inum_nodes ; ++i)
                      nodes[i] = inodes[i] ; } ;
        } ;

        struct SeamCoordCache {
            int num ;
            Vec2D coord[4] ;
            bool boundary ;
            bool ignore ;
        } ;

#endif

        enum DDFlags { RmshBoundary=1,
                       RmshTriBeforeSmooth=2,
                       RmshTriAfterSmooth=4,
                       RmshQuadBeforeCleanup=8,
                       RmshQuadBeforeSmooth=16,
                       RmshQuadAfterSmooth=32,
                       RmshQuadElemNumbers=64,
                       RmshAll=127 } ;

        int DebugDisplayFlags ;

        void SetDebugDisplayFlags(DDFlags flags) {
            DebugDisplayFlags |= flags ; } ;

        void SetDebugFile(FILE* fd) { dbg_file = fd ; }

    protected:

        class IntEdgeSet : public FTools::OrderedSet<IntEdge> {
            public:
                IntEdgeSet() : FTools::OrderedSet<IntEdge>() {}
                int Compare(const IntEdge& op0,
                            const IntEdge& op1) const ;
        } ;


        MshOrder        Order ;
//        MshElemType     Type ;
        int       StartId, StartIdSave ;
        int       StartElemId, StartElemIdSave ;
        int       MaxId ;
        int       MatId ;
        int                MaxNumElem ;
        bool               DoSmoothNodes ;
        bool               WinslowSmoothing ;

        Dict<int,IntNode> *pnode_table ;
        IntEdgeSet *pedge_table ;
        //Dict<EdgeKey,IntEdge> *pedge_table ;
        Dict<int,MshElement2D> *pelem_table ;
        Dict<int,SmallSet<int,2> > *MateTable ;
        List<Vec2D> RefinePts ;
        List<double> RefineSize ;

        int NewNode(double x,double y,MshNodeType type,
                    MshNodeMotion motion,bool corner) ;

        double Area(Vec2D &b,
                    Vec2D &i,
                    Vec2D &j) const ;

        double Angle(Vec2D &b,
                     Vec2D &i,
                     Vec2D &j) const ;

        double Angle2Pi(Vec2D &b,
                        Vec2D &i,
                        Vec2D &j) const ;

        Vec2D BisectNorm(Vec2D &b,
                         Vec2D &i,
                         Vec2D &j) const ;

        bool   Cross(Vec2D &i1,Vec2D &i2,
                     Vec2D &j1,Vec2D &j2) const ;

        double CrossProd(Vec2D &b,
                         Vec2D &i,
                         Vec2D &j) const ;

        double DistSqr(Vec2D &b,
                       Vec2D &i,
                       Vec2D &j,double *blen) const ;

        double TriMetric(Vec2D i,Vec2D j,Vec2D k) const ;

        double QuadMetric(Vec2D i,Vec2D j,Vec2D k,Vec2D l) const ;

        double QuadMetricAreaOnly(Vec2D i,Vec2D j,
                                  Vec2D k,Vec2D l) const ;

        double FlatnessMetric(Vec2D &i,Vec2D &j,
                              Vec2D &k,Vec2D &l) const ;

    private:

        class IntEdgePriority : public FTools::Heap<IntEdge> {
            public:
                IntEdgePriority() : FTools::Heap<IntEdge>() {}
                int Compare(const IntEdge& op0,
                            const IntEdge& op1) ;
        } ;

        bool               NodeGeneration ;
        bool               BoundaryChecks ;
        double             RefineBoundaryFactor ;
        double             NearBoundaryFactor ;
        double             MinimumShapeMeasure ;

        bool               QuadTreeDebug ;

        bool               CheckAllPoints ;

        FILE*  dbg_file ;

    // private member functions

        void   GenerateNodes() ;
        int    BoundaryContraction() ;
        void   SmoothNodes() ;

        void GetBdryInfo(int num_nodes,
                         List<IntNode*>& nodes,
                         OrderedSet<IntEdge>::SetIterator& eiter,
                         int *node_map,
                         Dict<int,Vec2D> &bdry_start,
                         Dict<int,Vec2D> &bdry_stop,
                         Dict<int,SmallSet<int,2> > &mate_tbl) ;

        bool CheckCross(IntEdgePriority *pedge_heap,
                        IntEdge **bdry,IntEdge *entry,
                        List<IntNode*> &nodes,
                        int *node_map,
                        IntNode *p,double blensqr,
                        IntNode *p_nd0,IntNode *p_nd1,
                        Vec2D &enorm,Vec2D &emid,
                        Dict<int,Vec2D> &bdry_start,
                        Dict<int,Vec2D> &bdry_stop,
                        int *lside_0,int *lside_1,
                        bool *nd_0_flg,bool *nd_1_flg,
                        Dict<int,SmallSet<int,2> > &mate_table) ;

#ifdef DOING_QUADS
        double AdjAspect(MshEdgeList::NewQuadData *qdata,
                         MshTopo2D *msh_topo) const ;

        double AdjNearBoundary(MshEdgeList::NewQuadData *qdata,
                               MshTopo2D *msh_topo,
                               MshEdgeList *edge_list) const ;

        bool   CloseSimplePoly(
                     int *current_num,
                     MshEdgeList::QdEdge *qd_edge,
                     MshTopo2D *msh_topo,
                     MshTopo2D *quad_topo,
                     MshEdgeList *edge_list) ;

        bool CloseSeam(
                     bool left,
                     int *current_num,
                     MshEdgeList::QdEdge *qd_edge,
                     MshTopo2D *msh_topo,
                     MshTopo2D *quad_topo,
                     MshEdgeList *edge_list,
                     bool *transition_flag,
                     bool *modified_flg) ;

        double CornerAspect(MshEdgeList::NewQuadData *qdata,
                            MshTopo2D *msh_topo) const ;

        bool ExtractOneTriangle(
            int *current_num,
            MshEdgeList::QdEdge *qd_edge,
            MshTopo2D *msh_topo,
            MshTopo2D *quad_topo,
            MshEdgeList *edge_list,
            int right_node,
            int left_node,
            int next_right,
            int next_left) ;

        bool CheckSurroundQuad(
            int *current_num,
            int *poly_nodes,
            MshTopo2D *msh_topo,
            MshTopo2D *quad_topo,
            MshEdgeList *edge_list) ;

        bool TransitionSplit(
                     int *current_num,
                     MshEdgeList::NewQuadData *qdata,
                     MshEdgeList::QdEdge *qd_edge,
                     MshTopo2D *msh_topo,
                     MshTopo2D *quad_topo,
                     MshEdgeList *edge_list,
                     bool *transition_flg) ;

        void UpdateForQuadMesh(MshTopo2D *quad_topo) ;

        bool   GenerateQuads() ;

        void   RenumberNodes() ;

        void InitializeQuadGen(
                MshTopo2D *msh_topo,
                MshEdgeList *edge_list) ;

        int FindElementSide(MshTopo2D *msh_topo,
                            IntNode *this_node,
                            IntNode *prev_node,
                            IntNode *next_node,
                            double base_length,
                            bool prev_edge) ;

        bool CheckAllValid(
            double threshold) const ;

        int CheckForEdge(MshTopo2D *msh_topo,
                         IntNode *this_node,
                         IntNode *prev_node,
                         IntNode *next_node,
                         bool prev_edge,
                         Vec2D normal,
                         int *vtx_before,
                         int *vtx_after) const ;

        bool CheckValidTriangles() ;

        Queue<MshTopo2D::Edge> *FindCrossedEdges(
                              int start_id,
                              int stop_id,
                              MshTopo2D *msh_topo,
                              MshEdgeList *edge_list,
                              List<int> &moved,
                              List<Vec2D> &orig) ;

        bool RecoverTopEdge(const int start_id,
                            const int stop_id,
                            MshTopo2D *msh_topo,
                            MshEdgeList *edge_list) ;

        void ClearRegion(const int num_corners,
                         const int *boundary,
                         MshTopo2D *msh_topo,
                         const MshTopo2D::Edge *edge) ;

        void SmoothFront(const int qcase,
                         const int num_frt_nodes,
                         const int *frt_nodes,
                         MshTopo2D *msh_topo,
                         MshTopo2D *quad_topo,
                         MshEdgeList *edge_list) ;

        void SmoothSeamFront(
                         MshEdgeList::NewSeamData *sdata,
                         MshTopo2D *msh_topo,
                         MshTopo2D *quad_topo,
                         MshEdgeList *edge_list) ;

        void SmoothOneFront(const int node_id,
                            MshTopo2D *msh_topo,
                            MshTopo2D *quad_topo,
                            MshEdgeList *edge_list) ;

        Vec2D LaplaceSmoothFrnt(const int vtx,
                               MshTopo2D *msh_topo,
                               MshTopo2D *quad_topo) ;

        Vec2D IsoParametricSmooth(IntNode *node,
                                 MshTopo2D *quad_topo) ;

        Vec2D LaplaceSmooth(IntNode *node,
                                  MshTopo2D *topo) ;

        void CheckSmoothCoords(IntNode *node,
                               Vec2D delta,
                               MshTopo2D *msh_topo,
                               MshTopo2D *quad_topo) ;

        void SmoothAdjacent(const int num_frt_nodes,
                            const int *frt_nodes,
                            MshTopo2D *topo,
                            const bool look_for_bdry,
                            const bool triangle,
                            MshEdgeList *edge_list) ;

        void SmoothOneAdj(const int vtx,
                          const int prev,
                          const int next,
                          MshTopo2D *topo,
                          const bool look_for_bdry,
                          const bool triangle,
                          MshEdgeList *edge_list) ;

        void FindQuadData(const MshEdgeList::QdEdge *edge,
                          MshEdgeList::NewQuadData *qdata,
                          MshEdgeList *edge_list,
                          MshTopo2D *msh_topo) ;

        bool FindSeamData(const MshEdgeList::QdEdge *edge,
                          MshEdgeList::NewSeamData *sdata,
                          MshEdgeList *edge_list,
                          MshTopo2D *msh_topo,
                          MshTopo2D *quad_topo) ;

        bool DoSeam(MshEdgeList::NewSeamData *sdata,
                    MshTopo2D *msh_topo,
                    MshTopo2D *quad_topo,
                    MshEdgeList *edge_list) ;

        bool DoTransitionSeam(MshEdgeList::QdEdge *edge,
                              MshEdgeList::NewSeamData *sdata,
                              MshTopo2D *msh_topo,
                              MshTopo2D *quad_topo,
                              MshEdgeList *edge_list) ;

        void DoTransitionSplit(MshEdgeList::QdEdge *edge,
                               MshEdgeList::NewQuadData *qdata,
                               MshTopo2D *msh_topo,
                               MshTopo2D *quad_topo,
                               MshEdgeList *edge_list) ;

        void DoTemplateSplit(MshEdgeList::NewQuadData *qdata,
                             MshTopo2D *msh_topo,
                             MshTopo2D *quad_topo,
                             MshEdgeList *edge_list,
                             const bool base) ;

        bool TweakNode(IntNode *nd,MshTopo2D *msh_topo,
                       List<int> &moved,
                       List<Vec2D> &orig) ;

        bool CheckValidCoord(IntNode *node,
                             Vec2D delta,
                             MshTopo2D *topo) const ;

        bool CheckValidTriangleList(const Vec2D coord,
                   const Vec2D delta,
                   const List<SeamCoordCache> *triangles) const ;

        bool CheckValidQuadList(const Vec2D coord,
                   const Vec2D delta,
                   const List<SeamCoordCache> *quads) const ;

        void DoPrint(
            MshTopo2D *msh_topo,
            FILE *fd = 0,
            bool label = 0) ;

        void QuadAngleChecks() ;

        bool QuadsInPolygon(int num,
                            const int *ids,
                            const Vec2D *vts,
                            const MshTopo2D *quad_topo,
                            const MshEdgeList *edge_list) const ;
#endif

        bool   ScanCross(Vec2D &b,
                         Vec2D &i,
                         Vec2D &j) const ;
        double TriShapeMeasure(Vec2D &i,
                               Vec2D &j,
                               Vec2D &k) const ;

#ifdef DOING_QUADS
        Vec2D IntersectLines(Vec2D i1,
                                   Vec2D i2,
                                   Vec2D j1,
                                   Vec2D j2) const ;
        List<Vec2D> *CoordCache ;
        //List<bool> *CoordCacheFlags ;
        //List<int> *BdryVtxCache ;

#endif

        void DebugAddNodeEcho(int id,
                              Vec2D coord,
                              bool retained_flag,
                              MshNodeType itype,
                              MshNodeMotion imotion) ;

        void DebugAddEdgeEcho(int nd0,int nd1,int nd2) ;

        void DebugAddRefineEcho(const Vec2D& xy,double size) ;


    public:

        class NodeIterator {
            public:
                NodeIterator(const MshRegion2D &reg,
                             const MshNodeType type=UNSPECIFIED) :
                    iter(reg.pnode_table),ntype(type) {} ;
                NodeIterator() {}
                NodeIterator(const NodeIterator& other) { Copy(other) ; }
                NodeIterator operator = (const NodeIterator& other) {
                    Copy(other) ;
                    return *this ;
                }
                void First() ;
                void Next() ;
                bool More()  { return(iter.More()) ; } ;
                int  Id()    { return(iter.Key()) ; } ;
                double NCoordI(int i) { return(iter.Entry().coord[i]) ; } ;
                Vec2D Coord2D() { return(iter.Entry().coord) ; } ;

                void operator ++ ()      { Next() ; } ;
                void operator ++ (int i) { Next() ; } ;

            private:
                Dict<int,IntNode>::DictIterator iter ;
                MshNodeType ntype ;

                void Copy(const NodeIterator& other) {
                    iter = other.iter ;
                    ntype = other.ntype ;
                }
        } ;

        class ElemIterator {
            public:
                ElemIterator(const MshRegion2D &reg) : iter(reg.pelem_table) {} ;
                ElemIterator() {} ;
                ElemIterator(const ElemIterator& other) { Copy(other) ; }
                ElemIterator operator = (const ElemIterator& other) {
                    Copy(other) ;
                    return *this ;
                }
                void First() { iter.First() ; };
                void Next()  { iter.Next() ; };
                bool More()  { return(iter.More()) ; } ;
                int  Id()    { return(iter.Entry().elem_id) ; } ;
                int  Mat()   { return(iter.Entry().mat_id) ; } ;
                int  NumNodes() { return(iter.Entry().num_nodes) ; } ;
                int  Node(int i) { return(iter.Entry().nodes[i]) ; } ;
                int* AllNodes()  { return(iter.Entry().nodes) ; } ;

                void operator ++ ()      { Next() ; } ;
                void operator ++ (int i) { Next() ; } ;

            private:
                Dict<int,MshElement2D>::DictIterator iter ;
                MshElement2D loc_elem ;

                void Copy(const ElemIterator& other) {
                    iter = other.iter ;
                    loc_elem = other.loc_elem ;
                }
        } ;

        class EdgeIterator {
            public:
                EdgeIterator(const List<UPair<int> >& iedges) :
                    edges(iedges),cur(0) {}
                EdgeIterator() {}
                EdgeIterator(const EdgeIterator& other) { Copy(other) ; }
                EdgeIterator operator = (const EdgeIterator& other) {
                    Copy(other) ;
                    return *this ;
                }
                void First()             { cur = 0 ; }
                void Next()              { ++cur ; }
                bool More()              { return cur < edges.Len() ; }
                int  Id0()               { return edges[cur][0] ; }
                int  Id1()               { return edges[cur][1] ; }
                void operator ++ ()      { Next() ; }
                void operator ++ (int i) { Next() ; }

            private:
                List<UPair<int> > edges ;
                int cur ;

                void Copy(const EdgeIterator& other) {
                    edges = other.edges ;
                    cur = other.cur ;
                }
        } ;

    // friends

    friend int DictHashIndex(const EdgeKey& key) ;

    friend void GenerateIntNodes(void *p_data,
                double origin_x,double origin_y,double half,
                bool is_root,bool is_leaf) ;

} ;

inline int DictHashIndex(const MshRegion2D::EdgeKey& key) {
    return key.id_0 ;
}


/*
CLASS MshRegion2D

  This object generates a 2D mesh of triangles or quadrilaterals for an 
  arbitrarily shaped region. 

  The normal procedure for using this object is as follows: 

  1. Create an instance of the object. 

  2. Add nodes to the region with AddNode and/or AddNodeList. 

  3. Add edges to the boundary of the region with AddEdge and/or 
  AddEdgeList. The edges must make up one or more loops that define 
  portions of the boundary of the region. There must be at least one 
  outside loop specified with a counter clockwise orientation. In 
  addition, there may be any number of internal loops specified with a 
  clockwise orientation. None of the edges defining the loops may 
  cross. 

  4. GenerateMesh is called the generate a mesh for the region. 

  5. GetNodes and GetElements are called to retrieve information about 
  the generated mesh. 


PUBLIC INTERFACE

  Public Data Structures:

    struct SeamElemCache

      This data structure is used to cache information about the 
      elements near a seam while the local topology is being updated. 

      Member Functions:

        SeamElemCache - no argument constructor 

          SeamElemCache()

          Description: This is a no argument constructor for a 
              SeamElemCache structure. 


        SeamElemCache - constructor 

          SeamElemCache(
                  const int ielem,
                  const int inum_nodes,
                  const int *inodes)

            ielem      - (in)  element id 
            inum_nodes - (in)  number of nodes 
            inodes     - (in)  node list 

          Description: This is a constructor for a SeamElemCache 
              structure. 


      Member Variables:

        int elem - element id 

        int num_nodes - number of nodes 

        int nodes[4] - node id's 


    struct SeamCoordCache

      This data structure is used to cache information about nodes near 
      a seam while the local topology is being updated. 

      Member Variables:

        Vec2D coord[4] - element coordinates 

        bool boundary - boundary flag 


  Public Member Functions:

    MshRegion2D - constructor method 

      MshRegion2D(
              const ArbMshOrder    iorder = LINEAR,
              const ArbMshElemType itype = TRIANGLE,
              const int            istart_id = 0,
              const int            imat_id = 0)

        iorder    - (in)  LINEAR or QUADRATIC 
        itype     - (in)  TRIANGLE or QUADRILATERAL 
        istart_id - (in)  starting generated node id 
        imat_id   - (in)  region's material id 

      Description: This is a constructor method for a MshRegion2D 
          object. 


    MshRegion2D - destructor 

      ~MshRegion2D()

      Description: This is a destructor for a MshRegion2D object. 


    SetOrder - sets the polynomial order 

      void SetOrder(const ArbMshOrder iorder)

        iorder - (in)  LINEAR or QUADRATIC 

      Description: This method sets the polynomial order for 
          generated elements. 


    SetElemType - sets the element type 

      void SetElemType(const ArbMshElemType itype)

        itype - (in)  TRIANGLE or QUADRILATERA 

      Description: This method sets the element type (shape) of the 
          generated elements. 


    SetStartNodeID - sets the first generated node id 

      void SetStartNodeID(const int is_id)

        is_id - (in)  first node id 

      Description: This method sets the first id to be assigned to 
          generated nodes. Subsequent id's will increase from this. 


    SetStartElemID - sets the first generated element id 

      void SetStartElemID(const int is_id)

        is_id - (in)  first element id 

      Description: This method sets the first id to be assigned to 
          generated elements. Subsequent id's will increase from 
          this. 


    SetMaterialID - sets the assigned material id 

      void SetMaterialID(const int mat_id)

        mat_id - (in)  material id 

      Description: This method sets the material id to be assigned to 
          generated elements. 


    SetNodeGeneration - set the node generation option 

      void SetNodeGeneration(bool flag = true)

        flag - (in)  on or off flag 

      Description: This method turns on or off the generation of 
          internal nodes during meshing. 


    SetBoundaryValidityChecks - set the boundary check option 

      void SetBoundaryValidityChecks(bool flag = true)

        flag - (in)  on or off flag 

      Description: This method turns on or off the option to check 
          the validity of the input boundary before elements are 
          generated. 


    SetNodeSmoothing - set the nodal smoothing option 

      void SetNodeSmoothing(bool flag = true)

        flag - (in)  on or off flag 

      Description: This method turns on or off the option to smooth 
          internal node to improve element quality. 


    SetTopoCleanup - set the topological cleanup option 

      void SetTopoCleanup(bool flag = true)

        flag - (in)  on or off flag 

      Description: This method turns on or off the option to perform 
          "topological cleanup" to improve element quality. 


    SetRefineBoundaryFactor - set the refine boundary factor 

      void SetRefineBoundaryFactor(double factor)

        factor - (in)  boundary factor 

      Description: This method sets a factor used when generating 
          internal nodes. Higher values will generate denser meshes. 


    SetNearBoundaryFactor - set the near boundary factor 

      void SetNearBoundaryFactor(double factor)

        factor - (in)  boundary factor 

      Description: This method sets a factor used when generating 
          internal nodes. The factor is used to remove generated 
          nodes that lie too near a boundary. 


    SetQuadTreeDebug - set a debug flag for the quadtree 

      void SetQuadTreeDebug(bool flag = true)

        flag - (in)  on or off 

      Description: This method sets a quadtree debugging flag. If set 
          the debugging data printed. 


    SetWinslowSmoothing - sets the Winslow smoothing option 

      void SetWinslowSmoothing(bool flag = true)

        flag - (in)  on or off 

      Description: This method turns on or off the option to us a 
          Winslow smoothing algorithm. If turned off a Laplace 
          algorithm is used. 


    GetOrder - get the polynomial order 

      ArbMshOrder GetOrder() const

      Description: This method returns the polynomial order of 
          generated elements. 

      Return Value: LINEAR or QUADRATIC 


    GetElemType - get the element type 

      ArbMshElemType GetElemType() const

      Description: This method returns the type (shape) of generated 
          elements 

      Return Value: TRIANGLE or QUADRILATERAL 


    GetStartNodeID - get the start node id 

      int GetStartNodeID() const

      Description: This method returns the first id that will be 
          assigned to generated internal nodes. 

      Return Value: lowest id of generated nodes 


    AddNode - add a node to the region 

      void AddNode(
              const ArbMshNode       &inode,
              const ArbMshNodeType   itype = BOUNDARY,
              const ArbMshNodeMotion imotion = ARB_FLOATING)

        inode   - (in)  node description 
        itype   - (in)  BOUNDARY or INTERIOR 
        imotion - (in)  ARB_FIXED or ARB_FLOATING 

      Description: This method adds a node to the boundary of of the 
          current region. 

      Exceptions:
          CArbMshDuplNode - duplicate node id


    AddNodeList - add a list of nodes to the region. 

      void AddNodeList(
              const int              num_nodes,
              const ArbMshNode       *constnode_list,
              const ArbMshNodeType   itype = BOUNDARY,
              const ArbMshNodeMotion imotion = ARB_FLOATING)

        num_nodes - (in)  number of nodes in list 
        node_list - (in)  list of nodes 
        itype     - (in)  BOUNDARY or INTERIOR 
        imotion   - (in)  ARB_FIXED or ARB_FLOATING 

      Description: This method adds a list of nodes to the boundary 
          of the current region. 

      Exceptions:
          CArbMshDuplNode - duplicate node id


    AddEdge - add an edge to the region 

      void AddEdge(const ArbMshEdge &edge)

        edge - (in)  edge description 

      Description: This method adds an edge to the boundary of the 
          current region. 

      Exceptions:
          CArbMshBadNode - bad node id
          CArbMshDuplEdge - duplicate edge description


    AddEdgeList - add a list of edges to the region. 

      void AddEdgeList(
              const int        num_edges,
              const ArbMshEdge *constedge_list)

        num_edges - (in)  number of edges in the list 
        edge_list - (in)  list of edge descriptions 

      Description: This method adds a list of edges to the boundary 
          of the current region. 

      Exceptions:
          CArbMshBadNode - bad node id
          CArbMshDuplEdge - duplicate edge description


    GenerateMesh - generate a mesh 

      void GenerateMesh()

      Description: This method generates a mesh for the region. 

      Exceptions:
          CArbMshCrossingBdry - the boundary crosses itself
          CArbMshIllegalBdry - the boundary is invalid


    CheckBoundary - perform a check of the boundary 

      void CheckBoundary()

      Description: This method performs a check of the region's 
          boundary to see if it is valid. 

      Exceptions:
          CArbMshCrossingBdry - the boundary crosses itself
          CArbMshIllegalBdry - the boundary is invalid


    NumBoundaryNodes - number of boundary nodes 

      int NumBoundaryNodes() const

      Description: This method returns the number of nodes on the 
          regions's boundary. 

      Return Value: number of boundary nodes 


    NumInternalNodes - number of internal nodes 

      int NumInternalNodes() const

      Description: This method returns the number of nodes in the 
          regions's Interior. 

      Return Value: number of internal nodes 


    NumNodes - number of nodes 

      int NumNodes() const

      Description: This method returns the total number of nodes in 
          the region. 

      Return Value: total number of nodes 


    NumBoundaryEdges - number of boundary edges 

      int NumBoundaryEdges() const

      Description: This method retuns the number of boundary edges in 
          the region. 

      Return Value: number of boundary edges 


    NumElements - number of elements 

      int NumElements() const

      Description: This method returns the number of elements in the 
          region. 

      Return Value: number of elements 


    GetBoundaryNodes - get the boundary nodes 

      ArbMshNode *GetBoundaryNodes() const

      Description: This method returns a list of the nodes on the 
          boundary of the region. Ownership of this memory passes to 
          the client, which must eventually call delete []. 

      Return Value: a list of the nodes on the boundary 


    GetInternalNodes - get the internal nodes 

      ArbMshNode *GetInternalNodes() const

      Description: This method returns a list of the nodes in the 
          interior of the region. Ownership of this memory passes to 
          the client, which must eventually call delete []. 

      Return Value: a list of the nodes in the interior 


    GetNodes - get all the nodes 

      ArbMshNode *GetNodes() const

      Description: This method returns a list of all the nodes in the 
          region. Ownership of this memory passes to the client, 
          which must eventually call delete []. 

      Return Value: a list of all the nodes in the region 


    GetBoundaryEdges - get the boundary edges 

      ArbMshEdge *GetBoundaryEdges() const

      Description: This method returns a list of all the edges the 
          make up the boundary of the region. Ownership of this 
          memory passes to the client, which must eventually call 
          delete []. 

      Return Value: a list of the boundary edges 


    GetElements - get the elements 

      ArbMshElement2D *GetElements() const

      Description: This method returns a list of all the elements in 
          the region. Ownership of this memory passes to the client, 
          which must eventually call delete []. 

      Return Value: a list of the elements 


    GetMeshStatistics - return element shape statistics 

      ArbMshStats *GetMeshStatistics() const

      Description: This method returns statistics regarding the 
          quality of the elements generated for the region. Ownership 
          of this structure passes to the client, which must 
          eventually call delete. 

      Return Value: statistics regarding element quality 


    DebugDumpRegion - debug routine to save region information 

      void DebugDumpRegion()

      Description: This method provides debugging support. It stores 
          information about the region in the file "debug.rgn". 


    DebugDumpMesh - debug routine to save mesh information 

      void DebugDumpMesh()

      Description: This method provides debugging support. It stores 
          information about the mesh in the files "debug.cds" and 
          "debug.msh". 


    DisplayMesh - debug routine to display a mesh 

      void DisplayMesh()

      Description: This method provides debugging support. It prints 
          the commands necessary to display the mesh. 


    DisplayBoundary - debug routine to display a boundary 

      void DisplayBoundary()

      Description: This method provides debugging support. It prints 
          the commands necessary to display the boundary of a region. 


    NewElemNum - generate a new element number 

      int NewElemNum()

      Description: This method generates a new and unique element id. 

      Return Value: a unique element id 


    DuplicateNode - generate a second node a coordinate location 

      int DuplicateNode(int id)

        id - (in)  existing node id 

      Description: Given the id of an existing node, this method 
          generates a new node that will share the coordinates off 
          the existing node, but have a unique node id. 

      Return Value: the id of the generated node 


PRIVATE INTERFACE

  Private Member Functions:

    GenerateNodes - generate interior nodes 

      void GenerateNodes()

      Description: This method uses a quadtree procedure to generate 
          interior nodes. This is done by: 

          1. Generate a quadtree based on the max and min extent of 
          the region. The tree is refined localy so that the region's 
          boundary edges fall in cells who's size about equal to the 
          boundary edge length (tuned by the refine boundary factor). 

          2. Refine the quadtree so that no cell is larger than the 
          largest cell on the boundary. 

          3. Refine the quadtree so that no two adjacent cells differ 
          in refinement level by more than one. 

          4. Attempt to generate a node at the center of each cell, 
          rejecting those that would fall too close to the boundary 
          (tuned by the near boundary factor). 


    BoundaryContraction - do boundary contraction mesh generation 

      void BoundaryContraction()

      Description: This method generates a triangular mesh for the 
          region by performing a boundary contraction (advancing 
          front) procedure. The procedure is to select an edge from 
          the current active boundary and find the node, which when 
          combined with this edge, forms the best triangle. The 
          boundary is updated and the procedure is repeated. 


    SmoothNodes - do nodal smoothing 

      void SmoothNodes()

      Description: This method performs smoothing (repositioning) of 
          internal nodes to improve element quality. This can be done 
          using either a Laplace or a Winslow (default) algorithm. 


    CloseSimplePoly - turn a simple polygonal region into an element 

      bool CloseSimplePoly(
              int                        *current_num,
              MshEdgeList::QdEdge *qd_edge,
              MshTopo2D              *msh_topo,
              MshTopo2D              *quad_topo,
              MshEdgeList            *edge_list)

        current_num - (i/o) current element number 
        qd_edge     - (in)  edge description 
        msh_topo    - (i/o) current triangular mesh description 
        quad_topo   - (i/o) current quadrilateral mesh description 
        edge_list   - (i/o) current active boundary 

      Description: This method is used during quadrilateral mesh 
          generation to turn a simple 3 or 4 sided region into an 
          element. 

      Return Value: true if an element was created 


    CloseSeam - performs a seam closing 

      bool CloseSeam(
              bool                       left,
              int                        *current_num,
              MshEdgeList::QdEdge *qd_edge,
              MshTopo2D              *msh_topo,
              MshTopo2D              *quad_topo,
              MshEdgeList            *edge_list,
              bool                       *transition_flag)

        left            - (in)  true if seam end is on the left 
        current_num     - (i/o) current element number 
        qd_edge         - (in)  edge description 
        msh_topo        - (i/o) current triangular mesh description 
        quad_topo       - (i/o) current quadrilateral mesh 
                                description 
        edge_list       - (i/o) current active boundary 
        transition_flag - (out) true if a transition seam was 
                                performed 

      Description: This method checks to see if a seaming operation 
          can be performed and if so does it. 

      Return Value: true if a seaming option was performed. 


    TransitionSplit - performs a transition split 

      bool TransitionSplit(
              int                             *current_num,
              MshEdgeList::NewQuadData *qdata,
              MshEdgeList::QdEdge      *qd_edge,
              MshTopo2D                   *msh_topo,
              MshTopo2D                   *quad_topo,
              MshEdgeList                 *edge_list,
              bool                            *transition_flg)

        current_num    - (in)  current element number 
        qdata          - (i/o) proposed quad elem data 
        qd_edge        - (i/o) edge description 
        msh_topo       - (i/o) current triangular mesh description 
        quad_topo      - (i/o) current quadrilateral mesh description 
        edge_list      - (i/o) current active boundary 
        transition_flg - (out) true if a split was performed 
                               (redundant with function value) 

      Description: This method checks to see if a transion split 
          operation can be performed and if so does it. 

      Return Value: true if a transition split was performed 


    UpdateForQuadMesh - update the region with quad mesh information 

      void UpdateForQuadMesh(MshTopo2D *quad_topo)

        quad_topo - (in)  quadrilateral mesh data 

      Description: This method updates the instance to replace the 
          triangular mesh information with quad mesh information. 


    GenerateQuads - driver to generate quad elements 

      void GenerateQuads()

      Description: This method is a driver for the Q-mesh 
          quadrilateral mesh generation algorithm. 


    RenumberNodes - renumber nodes so that unused numbers are removed 

      void RenumberNodes()

      Description: This method renumbers the generated nodes so that 
          nodes generated during triangulation but not used when 
          generating quads are eliminated. 


    InitializeQuadGen - initialization for quad generation 

      void InitializeQuadGen(
              MshTopo2D   *msh_topo,
              MshEdgeList *edge_list)

        msh_topo  - (in)  triangular mesh 
        edge_list - (i/o) active boundary list 

      Description: This method initializes the adjacent vertex and 
          edge list structures used for quadrilateral element 
          generation 


    FindElementSide - find a candidate element side 

      int FindElementSide(
              MshTopo2D *msh_topo,
              IntNode    *this_node,
              IntNode    *prev_node,
              IntNode    *next_node,
              bool          prev_edge)

        msh_topo  - (in)  triangular mesh 
        this_node - (in)  near node on the side 
        prev_node - (in)  adjacent node on boundary 
        next_node - (in)  adjacent node on boundary 
        prev_edge - (in)  flag that tells if we are looking for a 
                          right or left side 

      Description: Given nodes on the active boundary this method 
          finds a candidate for a side for a new element. This can be 
          either an existing edge, an edge generated by "swapping", 
          or a new edge generated by splitting one of the existing 
          edges. 

      Return Value: id of the far node on the side 


    CheckForEdge - find an edge in the search direction 

      int CheckForEdge(
              MshTopo2D *msh_topo,
              IntNode    *this_node,
              IntNode    *prev_node,
              IntNode    *next_node,
              bool          prev_edge,
              Vec2D   normal,
              int           *vtx_before,
              int           *vtx_after) const

        msh_topo   - (in)  triangular mesh 
        this_node  - (in)  near node on this side 
        prev_node  - (in)  adjacent node on boundary 
        next_node  - (in)  adjacent node on boundary 
        prev_edge  - (in)  flag that tells if we are looking for a 
                           right or left side 
        normal     - (in)  search normal direction 
        vtx_before - (out) bracketing vertex 
        vtx_after  - (out) bracketing vertex 

      Description: This method looks at the nodes adjacent to the 
          input node to see if there is one close to the normal 
          direction. If this is the case, the node is returned. If 
          not, the routine returns the nodes that bracket the search 
          direction. 

      Return Value: non-zero is the far node on a found edge. Zero 
          indicates that the bracketing edges are being returned. 


    FindCrossedEdges - find crossed edges 

      CArbQueue <MshTopo2D::Edge>*FindCrossedEdges(
              int             start_id,
              int             stop_id,
              MshTopo2D   *msh_topo,
              MshEdgeList *edge_list) const

        start_id  - (in)  start node id 
        stop_id   - (in)  end node id 
        msh_topo  - (in)  triangular mesh 
        edge_list - (in)  active boundary 

      Description: This method returns a list of triangle edges 
          crossed by a vector drawn between two nodes 

      Return Value: a ArbQueue containing a list of the edges crossed 


    RecoverTopEdge - do edge swapping 

      bool RecoverTopEdge(
              const int       start_id,
              const int       stop_id,
              MshTopo2D   *msh_topo,
              MshEdgeList *edge_list) const

        start_id  - (in)  start node id 
        stop_id   - (in)  end node id 
        msh_topo  - (i/o) triangular mesh 
        edge_list - (in)  active boundary 

      Description: This method performs "edge swapping" to recover 
          the top edge of a new element. 

      Return Value: true if succesful 


    ClearRegion - clear a quadrilateral region 

      void ClearRegion(
              const int                    num_corners,
              const int                    *boundary,
              MshTopo2D                *msh_topo,
              const MshTopo2D::Edge *edge)

        num_corners - (in)  3 or 4 
        boundary    - (in)  corner node id's 
        msh_topo    - (i/o) triangle mesh 
        edge        - (in)  one boundary edge 

      Description: Once the four corners of a quadrilateral or a 
          triangular region have been identified, this routine 
          deletes all the triangles that currently cover the region. 


    SmoothFront - smooth advancing front nodes 

      void SmoothFront(
              const int       qcase,
              const int       num_frt_nodes,
              const int       *frt_nodes,
              MshTopo2D   *msh_topo,
              MshTopo2D   *quad_topo,
              MshEdgeList *edge_list)

        qcase         - (in)  new quad case 
        num_frt_nodes - (in)  number of front nodes 
        frt_nodes     - (in)  list of front nodes 
        msh_topo      - (in)  triangle mesh 
        quad_topo     - (in)  quadrilateral mesh 
        edge_list     - (in)  active boundary 

      Description: Perform smoothing of advancing front nodes 


    SmoothSeamFront - smooth advancing front nodes 

      void SmoothSeamFront(
              MshEdgeList::NewSeamData *sdata,
              MshTopo2D                   *msh_topo,
              MshTopo2D                   *quad_topo,
              MshEdgeList                 *edge_list)

        sdata     - (in)  seaming data 
        msh_topo  - (in)  triangle mesh 
        quad_topo - (in)  quadrilateral mesh 
        edge_list - (in)  active boundary 

      Description: Perform smoothing of advancing front nodes after a 
          seaming operation. 


    SmoothOneFront - smooth one advancing front node 

      void SmoothOneFront(
              const int       node_id,
              MshTopo2D   *msh_topo,
              MshTopo2D   *quad_topo,
              MshEdgeList *edge_list)

        node_id   - (in)  id of node to smooth 
        msh_topo  - (in)  triangle mesh 
        quad_topo - (in)  quadrilateral mesh 
        edge_list - (in)  active boundary 

      Description: Perform smoothing for one advancing front node. 


    LaplaceSmoothFrnt - smooth one advancing front node 

      Vec2D LaplaceSmoothFrnt(
              const int     vtx,
              MshTopo2D *msh_topo,
              MshTopo2D *quad_topo)

        vtx       - (in)  id of node to smooth 
        msh_topo  - (in)  triangle mesh 
        quad_topo - (in)  quadrilateral mesh 

      Description: Performs smoothing of one advancing front node 
          using a Laplace algorithm. 

      Return Value: new nodal coordinates 


    IsoParametricSmooth - smooth one advancing front node 

      Vec2D IsoParametricSmooth(
              IntNode    *node,
              MshTopo2D *quad_topo)

        node      - (in)  id of node to smooth 
        quad_topo - (in)  quadrilateral mesh 

      Description: Performs smoothing of one advancing front node 
          using an isoparametric algorithm. 

      Return Value: new nodal coordinates 


    LaplaceSmooth - smooth one internal node 

      Vec2D LaplaceSmooth(
              IntNode    *node,
              MshTopo2D *topo)

        node - (in)  id of node to smooth 
        topo - (in)  mesh topology 

      Description: Performs smoothing of one internal node. 

      Return Value: new nodal coordinates 


    CheckSmoothCoords - check a proposed smoothed node location 

      void CheckSmoothCoords(
              IntNode    *node,
              Vec2D   delta,
              MshTopo2D *msh_topo,
              MshTopo2D *quad_topo)

        node      - (in)  node id 
        delta     - (in)  proposed coordinate update 
        msh_topo  - (in)  triangle mesh 
        quad_topo - (in)  quadrilateral mesh 

      Description: Check a proposed new node location to make sure 
          that all of the adjacent elements will have a valid shape 
          (not inverted). If the elements are valid the nodal 
          coordinates are updated. 


    SmoothAdjacent - smooth adjacent nodes 

      void SmoothAdjacent(
              const int       num_frt_nodes,
              const int       *frt_nodes,
              MshTopo2D   *topo,
              const bool      look_for_bdry,
              const bool      triangle,
              MshEdgeList *edge_list)

        num_frt_nodes - (in)  number of front nodes 
        frt_nodes     - (in)  list of front nodes 
        topo          - (in)  mesh topology 
        look_for_bdry - (in)  ignored 
        triangle      - (in)  true means a triangle mesh 
        edge_list     - (in)  active boundary 

      Description: Smooth all internal nodes adjacent to active front 
          nodes. 


    SmoothOneAdj - smooth one adjacent node 

      void SmoothOneAdj(
              const int       vtx,
              const int       prev,
              const int       next,
              MshTopo2D   *topo,
              const bool      look_for_bdry,
              const bool      triangle,
              MshEdgeList *edge_list)

        vtx           - (in)  input vertex id 
        prev          - (in)  adjacent vertex on boundary 
        next          - (in)  adjacent vertex on boundary 
        topo          - (in)  mesh topology 
        look_for_bdry - (in)  ignored 
        triangle      - (in)  true means triangle mesh 
        edge_list     - (in)  active boundary 

      Description: Smooth all internal nodes adjacent to one active 
          front node. 

      Return Value: ? 


    FindQuadData - find a candidate quadrilateral 

      void FindQuadData(
              const MshEdgeList::QdEdge *edge,
              MshEdgeList::NewQuadData  *qdata,
              MshEdgeList                  *edge_list,
              MshTopo2D                    *msh_topo)

        edge      - (in)  base edge for the quad 
        qdata     - (out) quadrilateral data 
        edge_list - (in)  active boundary 
        msh_topo  - (in)  triangle mesh 

      Description: This method identifies and generates information 
          about a candidate new quadrilateral element. 


    FindSeamData - find data for a seam operation 

      bool FindSeamData(
              const MshEdgeList::QdEdge *edge,
              MshEdgeList::NewSeamData  *sdata,
              MshEdgeList                  *edge_list,
              MshTopo2D                    *msh_topo,
              MshTopo2D                    *quad_topo)

        edge      - (in)  base edge 
        sdata     - (out) seam data 
        edge_list - (in)  active boundary 
        msh_topo  - (in)  triangle mesh 
        quad_topo - (in)  quadrilateral mesh 

      Description: This method generates information needed for a 
          seam operation. 

      Return Value: true if valid seam data is returned 


    DoSeam - do a seam update 

      bool DoSeam(
              MshEdgeList::NewSeamData *sdata,
              MshTopo2D                   *msh_topo,
              MshTopo2D                   *quad_topo,
              MshEdgeList                 *edge_list)

        sdata     - (in)  seam data 
        msh_topo  - (i/o) triangle mesh 
        quad_topo - (i/o) quadrilateral mesh 
        edge_list - (i/o) active boundary 

      Description: This method performs a seam update 

      Return Value: true if the update was successful 


    DoTransitionSeam - do a transition seam 

      bool DoTransitionSeam(
              MshEdgeList::QdEdge      *edge,
              MshEdgeList::NewSeamData *sdata,
              MshTopo2D                   *msh_topo,
              MshTopo2D                   *quad_topo,
              MshEdgeList                 *edge_list)

        edge      - (in)  base edge 
        sdata     - (in)  seam data 
        msh_topo  - (i/o) triangle mesh 
        quad_topo - (i/o) quadrilateral mesh 
        edge_list - (i/o) active boundary 

      Description: This method performs a transition seam update 

      Return Value: true if the update was successful 


    DoTransitionSplit - do a transition split 

      void DoTransitionSplit(
              MshEdgeList::QdEdge      *edge,
              MshEdgeList::NewQuadData *qdata,
              MshTopo2D                   *msh_topo,
              MshTopo2D                   *quad_topo,
              MshEdgeList                 *edge_list)

        edge      - (in)  base edge 
        qdata     - (in)  quad data 
        msh_topo  - (i/o) triangle mesh 
        quad_topo - (i/o) quadrilateral mesh 
        edge_list - (i/o) active boundary 

      Description: This method performs a transition split update 


    DoTemplateSplit - do a template split 

      void DoTemplateSplit(
              MshEdgeList::NewQuadData *qdata,
              MshTopo2D                   *msh_topo,
              MshTopo2D                   *quad_topo,
              MshEdgeList                 *edge_list,
              const bool                      base)

        qdata     - (in)  quad data 
        msh_topo  - (i/o) triangle mesh 
        quad_topo - (i/o) quadrilateral mesh 
        edge_list - (i/o) active boundary 
        base      - (in)  true means split the base 

      Description: This method performs a template split update 


    CheckValidCoord - check for valid triangles 

      bool CheckValidCoord(
              IntNode    *node,
              Vec2D   delta,
              MshTopo2D *topo) const

        node  - (in)  node to move 
        delta - (in)  proposed update 
        topo  - (in)  mesh topology 

      Description: Given a proposed nodal coordinate update this 
          method checks to see if all the adjacent triangular 
          elements will still be valid. 

      Return Value: true if the update is valid 


    CheckValidTriangleList - check a list of triangles 

      bool CheckValidTriangleList(
              const Vec2D                coord,
              const Vec2D                delta,
              const List<SeamCoordCache>* triangles) const

        coord     - (in)  nodal coord 
        delta     - (in)  proposed update 
        triangles - (in)  list of triangles 

      Description: This method check a list of triangular element to 
          see if they are all valid. 

      Return Value: true if all are valid 


    CheckValidQuadList - check a list of quadrilaterals 

      bool CheckValidQuadList(
              const Vec2D                coord,
              const Vec2D                delta,
              const List<SeamCoordCache>* triangles) const

        coord     - (in)  nodal coord 
        delta     - (in)  proposed update 
        triangles - (in)  list of quads 

      Description: This method check a list of quadrilateral element 
          to see if they are all valid. 

      Return Value: true if all are valid 


    DoPrint - debug routine to display a mesh 

      void DoPrint(MshTopo2D *msh_topo)

        msh_topo - (in)  mesh topology 

      Description: This method provides debugging support. It prints 
          information that can be used to display a mesh. 


    QuadCleanup - do topological cleanups 

      void QuadCleanup(MshTopo2D *quad_topo)

        quad_topo - (in)  quad mesh topology 

      Description: This method performs topological cleanups on a 
          quadrilateral mesh. 


    NewNode - generate a new node 

      int NewNode(
              double           x,
              double           y,
              ArbMshNodeType   type,
              ArbMshNodeMotion motion,
              bool             corner)

        x      - (in)  the node's x coordinate 
        y      - (in)  the node's y coordinate 
        type   - (in)  BOUNDARY or INTERIOR 
        motion - (in)  ARB_FIXED or ARB_FLOATING 
        corner - (in)  corner node flag 

      Description: This method generates a new node and assigns it a 
          unique node id. 

      Return Value: the new node id 


    Cross - check to see if lines cross 

      bool Cross(
              Vec2D i1,
              Vec2D i2,
              Vec2D j1,
              Vec2D j2) const

        i1 - (in)  first node of first line 
        i2 - (in)  second node of first line 
        j1 - (in)  first node of second line 
        j2 - (in)  second node of second line 

      Description: This method check to see if two line segments 
          cross. 

      Return Value: true if they cross 


    CrossProd - compute a cross product 

      double CrossProd(
              Vec2D b,
              Vec2D i,
              Vec2D j) const

        b - (in)  base vertex coordinate 
        i - (in)  end of first vector 
        j - (in)  end of second vector 

      Description: This method computes the cross product of two 
          vector that share a common starting vertex. 

      Return Value: the (2D) cross product value 


    Angle - compute an angle 

      double Angle(
              Vec2D b,
              Vec2D i,
              Vec2D j) const

        b - (in)  base vertex coordinate 
        i - (in)  end of first vector 
        j - (in)  end of second vector 

      Description: Compute the magnitude of the included angle 
          between two vectors that share a common starting vertex (in 
          the range -pi to pi). 

      Return Value: the magnitude of the angle (in radians) 


    Angle2Pi - compute an angle 

      double Angle2Pi(
              Vec2D b,
              Vec2D i,
              Vec2D j) const

        b - (in)  base vertex coordinate 
        i - (in)  end of first vector 
        j - (in)  end of second vector 

      Description: Compute the magnitude of the included angle 
          between two vectors that share a common starting vertex (in 
          the range 0 to 2*pi). 

      Return Value: the magnitude of the angle (in radians) 


    Area - find the area of a triangle 

      double Area(
              Vec2D b,
              Vec2D i,
              Vec2D j) const

        b - (in)  first vertex 
        i - (in)  second vertex 
        j - (in)  third vertex 

      Description: Compute the area of a triangle. 

      Return Value: triangle area 


    DistSqr - compute the distance from a point to a line segment 

      double DistSqr(
              Vec2D b,
              Vec2D i,
              Vec2D j,
              double      *blen) const

        b    - (in)  input vertex 
        i    - (in)  first line segment vertex 
        j    - (in)  second line segment vertex 
        blen - (out) distance from i to j 

      Description: This method computes the quare of the minimum 
          distance from a point to a line segment 

      Return Value: the square of the distance 


    ScanCross - check for a scan line crossing 

      bool ScanCross(
              Vec2D b,
              Vec2D i,
              Vec2D j) const

        b - (in)  input point 
        i - (in)  first line segment vertex 
        j - (in)  second line segment vertex 

      Description: This method determines if a horizontal line drawn 
          from minus infinity to a given point crosses a given line 

      Return Value: true if crosses 


    TriShapeMeasure - compute a shape measure 

      double TriShapeMeasure(
              Vec2D b,
              Vec2D i,
              Vec2D j) const

        b - (in)  first vertex 
        i - (in)  second vertex 
        j - (in)  third vertex 

      Description: This method computes a shape measure for a 
          triangular element. 

      Return Value: triangle's shape measure 


    IntersectLines - find the intersection point of two lines 

      Vec2D IntersectLines(
              Vec2D i1,
              Vec2D i2,
              Vec2D j1,
              Vec2D j2) const

        i1 - (in)  first node of first line 
        i2 - (in)  second node of first line 
        j1 - (in)  first node of second line 
        j2 - (in)  second node of second line 

      Description: This method finds the intersection point of two 
          line segments. 

      Return Value: the intersection coordinates 


    TriMetric - compute a shape metric 

      double TriMetric(
              Vec2D b,
              Vec2D i,
              Vec2D j) const

        b - (in)  first vertex 
        i - (in)  second vertex 
        j - (in)  third vertex 

      Description: This method computes a shape metric for a 
          triangular element. A zero or negative metric indicates an 
          invalid element. The maximum possible metric value is 1.0 
          for an equilateral triangle. 

      Return Value: the shape metric 


    QuadMetric - compute a shape metric 

      double QuadMetric(
              Vec2D b,
              Vec2D i,
              Vec2D j,
              Vec2D l) const

        b - (in)  first vertex 
        i - (in)  second vertex 
        j - (in)  third vertex 
        l - (in)  forth vertex 

      Description: This method computes a shape metric for a 
          quadrilateral element. A zero or negative metric indicates 
          an invalid element. The maximum possible metric value is 
          1.0 for a square. 

      Return Value: the shape metric 


  Private Member Variables:

    ArbMshOrder Order - LINEAR or QUADRATIC 

    ArbMshElemType Type - TRIANGLE or QUADRILATERAL 

    int StartId - next node id to assign 

    int StartIdSave - temp starting node id 

    int StartElemId - next elem id to assign 

    int StartElemIdSave - temp starting elem id 

    int MaxId - maximum assigned node id 

    int MatId - material id 

    bool NodeGeneration - flag for internal node generation 

    bool BoundaryChecks - flag to perform boundary checs 

    bool DoSmoothNodes - flag for node smoothing 

    bool DoCleanup - flag for topological cleanups 

    double RefineBoundaryFactor - refine boundary factor 

    double NearBoundaryFactor - near boundary factor 

    double MinimumShapeMeasure - minimum acceptable shape measure 

    bool WinslowSmoothing - flag for Winslow smoothing 

    bool QuadTreeDebug - flag for quad tree debugging 

    Dict<int,IntNode>* pnode_table - node table 

    CArbSet<ArbIntEdge>* pedge_table - edge table 

    CArbSet<ArbMshElement2D>* pelem_table - elem_table 

    List<Vec2D>* CoordCache - temp coordinate storage 

    List<int>* BdryVtxCache - temp vertex storage 

*/

} // namespace

#endif
