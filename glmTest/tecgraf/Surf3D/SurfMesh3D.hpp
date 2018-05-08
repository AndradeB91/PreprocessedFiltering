
#ifndef SurfMesh_hpp
#define SurfMesh_hpp

#include "CubicBezSurf.hpp"
#include "Dict.hpp"
#include "RangeTree.hpp"
#include "SmallSet.hpp"
#include "Vec3D.hpp"
#include "Vec2D.hpp"

#include "MshTopo2D.hpp"

using FTools::CubicBezSurf ;
using FTools::Dict ;
using FTools::RangeTree ;
using FTools::SmallSet ;
using FTools::Vec2D ;
using FTools::Vec3D ;

using Msh2D::MshTopo2D ;

namespace Msh3D {

class SurfMesh3D {

    public:

        enum {NORMAL_STATUS, DUPLICATE_NODE_ID, BAD_NODE_ID, CANNOT_MESH} ;

        enum NodeType {BOUNDARY, INTERIOR, UNSPECIFIED} ;

        enum NodeMotion {M_FIXED, M_FLOATING} ;

        struct MshNode {
            int id ;
            Vec3D   coord ;
            Vec2D uv_coord ;
            Vec3D   normal ;
            NodeType    type ; 
            NodeMotion  motion ;
        } ;

        struct MshEdge {
            int node_id[2] ;
            double length ;
            Vec3D  normal ;
            bool bdry ;
        } ;

        struct MshElement {
            int elem_id ;
            int mat_id ;
            int num_nodes ;
            int nodes[4] ;
        } ;

        struct EdgeKey {
            int id_0 ;
            int id_1 ;
            EdgeKey() : id_0(0), id_1(0) {} ;
            EdgeKey(int id0,int id1) {
                if (id0 < id1) {
                   id_0 = id0 ;  id_1 = id1 ;
                } else {  
                   id_0 = id1 ;  id_1 = id0 ;
                }
            }
        } ;

        class NodeIterator ;
        class ElemIterator ;

        SurfMesh3D() ;

        ~SurfMesh3D() {
            if (PatchRangeTree != 0) delete PatchRangeTree ; }

        int AddEdge(
            int node0,
            int node1) ;

        int AddNode(
            int id,
            const Vec3D& coord,
            const Vec2D& uv_coord,
            const NodeType type = BOUNDARY,
            const NodeMotion motion = M_FLOATING) ;

        int AddPatch(
            int id,
            const Vec3D* cntrl_points) ;

        int AddSurface(
            int patch_id,
            const Vec2D uv[3],
            const Vec2D rs[3]) ;

        int GenerateMesh() ;

        void SetCheckAllPoints() { CheckAllPoints = true ; }
        void SetMaxNumElements(int max_elem) { MaxNumElem = max_elem ; }
        void SetStartNodeID(const int is_id) { StartNodeId = is_id ; }
        void SetStartElemID(const int is_id) { StartElemId = is_id ; }
        void SetMaterialID(const int mat_id) { MatId = mat_id ; }
        void SetNodeGeneration(bool flag = true) { NodeGeneration = flag ; }
        void SetNodeSmoothing(bool flag = true) { DoSmoothNodes = flag ; }
        void SetRefineBoundaryFactor(double factor) { RefineBoundaryFactor = factor ; }
        void SetNearBoundaryFactor(double factor) { NearBoundaryFactor = factor ; }
        void SetMaxSmoothIterations(int maxits) { Max_it = maxits ; }
        void SetStopTolerance(double tol) { It_tol = tol ; }
        void SetDebugPlot() { DebugPlot = true ; }

        NodeIterator GetNodeIterator() { return(NodeIterator(&NodeTable)) ; }

        ElemIterator GetElemIterator() { return(ElemIterator(&ElemTable)) ; }

        class MaxElements {} ;

#ifdef DEBUG_FAC
        void DebugDump() ;
#endif

    private:

        class FaceSurfaceRef {
            public:
                FaceSurfaceRef(int patch_id,
                               const Vec2D uv[3],
                               const Vec2D rs[3]) ;
                FaceSurfaceRef() {} ;

                bool Contains(double u,double v) const ;
                Vec2D GetNatural(double u,double v) const ;

                int PatchId ;
                Vec2D CornerUV[3] ;
                Vec2D CornerRS[3] ;
                Vec2D MinUV,MaxUV ;
                double Tolerance ;
        } ;

        Dict<int,MshNode> NodeTable ;
        Dict<int,MshEdge> EdgeTable ;
        Dict<int,MshElement> ElemTable ;

        bool   CheckAllPoints ;
        int    MaxNumElem ;
        int    StartNodeId ;
        int    StartElemId ;
        int    MatId ;
        bool   NodeGeneration ;
        bool   DoSmoothNodes ;
        double RefineBoundaryFactor ;
        double NearBoundaryFactor ;
        int    MaxId ;
        int    EdgeId ;
        int    SurfId ;
        int    Num_it ;
        int    Max_it ;
        double It_tol ;
        bool   DebugPlot ;

        Dict<int,CubicBezSurf> Patches ;
        Dict<int,FaceSurfaceRef> Surfaces ;
        RangeTree* PatchRangeTree ;

        int BoundaryContraction() ;

        double Angle(
            const Vec3D &b,
            const Vec3D &i,
            const Vec3D &j) const ;

        bool Between(
            const Vec3D& b,
            const Vec3D& n,
            const Vec3D& p,
            const Vec3D& i,
            const Vec3D& j) const ;

        MshTopo2D* BuildMeshTopo() ;

        bool CheckCross(
            const MshEdge* entry,
            const MshNode* p,
            const MshNode* p_nd0,
            const MshNode* p_nd1,
            double blen,
            int num_edges,
            MshEdge** bdry,
            Dict<int,Vec3D> &bdry_start,
            Dict<int,Vec3D> &bdry_stop,
            Dict<int,SmallSet<int,2> > &mate_table,
            bool& nd_0_flg,
            bool& nd_1_flg) ;

        bool Cross(
            const Vec3D& i1,
            const Vec3D& i2,
            const Vec3D& j1,
            const Vec3D& j2) const ;

        bool CrossInPlane(
            const Vec3D& i1,
            const Vec3D& i2,
            const Vec3D& j1,
            const Vec3D& j2,
            const Vec3D& norm) const ;

        double CrossProdMag(
            const Vec3D& b,
            const Vec3D& n,
            const Vec3D& i,
            const Vec3D& j) const ;

        void ComputeNodeNormals() ;

        double Dist(
            const Vec3D& b,
            const Vec3D& i,
            const Vec3D& j,
            double& blen) const ;

        double Dist2D(
            const Vec2D& b,
            const Vec2D& i,
            const Vec2D& j) const ;

        void DoSmoothNodesLaplace(
            MshTopo2D *msh_topo,
            int max_it,
            double adapt_tol,
            bool skip_checks=false) ;

        void FindBoundaryNodes() ;

        double FindMinCharEdgeLength(
            MshTopo2D *msh_topo) ;

        bool FindClosestSurfacePoint(
            const Vec3D& trial,
            Vec3D& surf_pt,
            int& surf_id,
            int surf_hint = -1) ;

        void GenerateNodes() ;
        void GenerateNodesFromPatches() ;

        void GetBdryInfo(
            Dict<int,Vec3D>& bdry_start,
            Dict<int,Vec3D>& bdry_stop,
            Dict<int,SmallSet<int,2> >& mate_table) const ;

        Vec3D LaplaceUpdate(
            int num_elem,
            List<int> *num_adj_nodes,
            List<Vec3D> *adj_coords,
            bool skip_checks) ;

        int NewElemNum() { return(StartElemId++) ; }

        int NewNode(
            const Vec3D& coord,
            const Vec2D& uv_coord,
            const Vec3D& norm,
            NodeType type,
            NodeMotion motion) ;

        double QuadMetric(
            Vec3D i,
            Vec3D j,
            Vec3D k,
            Vec3D l) const ;

       bool ScanCross(
            const Vec2D& b,
            const Vec2D& i,
            const Vec2D& j) const ;

        void SmoothNodesLaplace() ;

        void UVtoXYZ(
            const Vec2D& uv,
            Vec3D& xyz,
            Vec3D& norm) const ;

        double TriMetric(
            Vec3D i,
            Vec3D j,
            Vec3D k) const ;

    public:

        class NodeIterator {
            public:
                NodeIterator(
                   Dict<int,MshNode>* tbl) : iter(tbl) {}
                void First()      { iter.First() ; }
                void Next()       { iter.Next() ; }
                bool More()       { return(iter.More()) ; }
                int  Id()         { return(iter.Key()) ; }
                Vec3D Coord() { return(iter.Entry().coord) ; }
 
                void operator ++ ()      { Next() ; }
                void operator ++ (int i) { Next() ; }

            private:
                Dict<int,MshNode>::DictIterator iter ;
        } ;

        class ElemIterator {
            public:
                ElemIterator(
                   Dict<int,MshElement>* tbl) : iter(tbl) {}
                void First()       { iter.First() ; }
                void Next()        { iter.Next() ; }
                bool More()        { return(iter.More()) ; }
                int  Id()          { return(iter.Key()) ; }
                int  MatId()       { return(iter.Entry().mat_id) ; }
                int  NumNodes()    { return(iter.Entry().num_nodes) ; }
                int  NodeId(int i) { return(iter.Entry().nodes[i]) ; }
 
                void operator ++ ()      { Next() ; }
                void operator ++ (int i) { Next() ; }

            private:
                Dict<int,MshElement>::DictIterator iter ;
        } ;


    friend void SurfGenIntNodes(void *p_data,
                    double origin_x,double origin_y,double half,
                    bool is_root,bool is_leaf) ;
} ;

inline bool operator == (
    const SurfMesh3D::MshEdge& op1,
    const SurfMesh3D::MshEdge& op2)
{
    return((op1.node_id[0] == op2.node_id[0]) &&
           (op1.node_id[1] == op2.node_id[1])) ;
}

inline bool operator == (
    const SurfMesh3D::EdgeKey &op1,
    const SurfMesh3D::EdgeKey &op2)
{
    return((op1.id_0 == op2.id_0) && (op1.id_1 == op2.id_1)) ;
}

inline int DictHashIndex(
    const SurfMesh3D::EdgeKey &edge)
{
    return(edge.id_0 + edge.id_1) ;
}

} // namespace

#endif
