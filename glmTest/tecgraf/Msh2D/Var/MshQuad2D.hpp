//
// MshQuad2D Class header file
//
// Description -
//   This class implements an object that will create a quadrialteral mesh
//   for an arbitrarily shaped region.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2007
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//

#ifndef MshQuad2D_h
#define MshQuad2D_h

#include "MshRegion2D.hpp"

namespace Msh2D {

#define DOING_QUADS 0

//
// Object descritpion
//

class MshQuad2D : public MshRegion2D {

    public:

        enum AngleCheckLocation { AT_NODES, AT_INT_PTS } ;

        MshQuad2D(const MshOrder iorder = LINEAR,
                  const MshElemType itype = TRIANGLE,
                  const int istart_id = 0,
                  const int imat_id = 0) :
            MshRegion2D(iorder,istart_id,imat_id),Type(itype),
            DoCleanup(true),DoQuadAngleChecks(false),
            CoordCache(new List<Vec2D>),
            CoordCacheFlags(new List<bool>),
            BdryVtxCache(new List<int>) {}

        virtual ~MshQuad2D() { delete CoordCache ;
                               delete CoordCacheFlags ;
                               delete BdryVtxCache ; }

        void SetElemType(const MshElemType itype) { Type = itype ; }
        void SetTopoCleanup(bool flag = true) { DoCleanup = flag ; }
        void SetQuadElemAngleChecks(double min_angle,
                                    double max_angle,
                                    AngleCheckLocation location) ;

        MshElemType GetElemType() const { return(Type) ; }

        int NewElemNum() { return(StartElemId++) ; } ;
        int DuplicateNode(int id) ;

        void GenerateMesh() ;

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


    private:

        MshElemType        Type ;
        bool               DoCleanup ;
        bool               DoQuadAngleChecks ;
        double             MinQuadAngle ;
        double             MaxQuadAngle ;
        AngleCheckLocation QuadAngleLocation ;

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

        void QuadCleanup(MshTopo2D *quad_topo) ;

        void QuadAngleChecks() ;

        bool QuadsInPolygon(int num,
                            const int *ids,
                            const Vec2D *vts,
                            const MshTopo2D *quad_topo,
                            const MshEdgeList *edge_list) const ;

        Vec2D IntersectLines(Vec2D i1,
                                   Vec2D i2,
                                   Vec2D j1,
                                   Vec2D j2) const ;

        List<Vec2D> *CoordCache ;
        List<bool> *CoordCacheFlags ;
        List<int> *BdryVtxCache ;
} ;

} // namespace

#endif
