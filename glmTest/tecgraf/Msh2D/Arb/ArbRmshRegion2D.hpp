//
// ArbRmshRegion2D header file
//
// Description -
//   This is the header file for the ArbRmshRegion2D objects.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2000,2001,2002
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.16 $  $Date: 2002/09/04 20:51:17 $  $Author: wash $
//
 
#ifndef ArbRmshRegion2D_h
#define ArbRmshRegion2D_h

#define DEFAULT_ID 100000000

#include "ArbGeom2DMixIn.hpp"
#include "ArbMshCrackRegion2D.hpp"
//#include "ArbSet.hpp"
#include "ArbBdrySeg.hpp"
#include "ArbSmallSet.hpp"

typedef CArbHashTable<int,int> CornerNodesSet ;

class CArbRmshRegion2D : public CArbGeom2DMixIn {

    public:

        enum AngleCheckLocation { AT_NODES, AT_INT_PTS } ;

        // constructors and destructor

        CArbRmshRegion2D(CArbMshCrackRegion2D *owner) ;
        ~CArbRmshRegion2D() ;

        // routines that support crack insertion/growth

        void AddFlaw(const int num_tips,
                     const int *tip_ids,
                     const int *tip_indices,
                     const int num_points,
                     const CArbMshCrackRegion2D::CrackPt *const pts) ;

        void SetTriDeleteFactor(double f)  { TriDeleteFactor = f ; } ;
        void SetQuadDeleteFactor(double f) { QuadDeleteFactor = f ; } ;

        void SetQuadElemAngleChecks(double min_angle,
                                    double max_angle,
                                    AngleCheckLocation location)
            { DoQuadAngleChecks = true ;
              MinQuadAngle = min_angle ;
              MaxQuadAngle = max_angle ;
              QuadAngleLocation = location ; } ;

        double GetTriDeleteFactor()  { return(TriDeleteFactor) ; } ;
        double GetQuadDeleteFactor() { return(QuadDeleteFactor) ; } ;

        CArbArray<int> *GetMouthNodes() { return(MouthNodes) ; } ;

        void SetAnisotropicRatio(double val)
            { AnisotropicRatio = val ; } ;

        void SetOldCrackTip(CArbCoord2D &tip_coord) ;

        // debug display options

        enum DDFlags { RmshAll=63 } ;

        int DebugDisplayFlags ;

        void SetDebugDisplayFlags(DDFlags flags) {
            DebugDisplayFlags |= flags ; } ;

        // these data structures should be private, but VC++ does
        // not deal with friend access to private structures

        class EdgeKey {
            public:
                EdgeKey() : id0(0),id1(0) {} ;
                EdgeKey(int iid0,int iid1) {
                    if (iid0 < iid1) {
                        id0 = iid0 ;  id1 = iid1 ;
                    } else {
                        id0 = iid1 ;  id1 = iid0 ;
                    }
                }
                int id0, id1 ;
        } ;

        class EdgeFaceKey {
            public:
                EdgeFaceKey() : vtx(0),face(0) {} ;
                EdgeFaceKey(int ivtx,int iface) :
                    vtx(ivtx),face(iface) {} ;
                int vtx, face ;
        } ;

    private:

        // private data structures and types

        enum PointType {INSIDE, OUTSIDE, BOUNDARY, INTERFACE} ;

        struct LocCrackPt {
            CArbCoord2D coord ;
            double char_elem_size ;
            double boundary_tol ;
            bool respect_char_size ;
            int node_id ;
            int elem ;
            PointType type ;
            bool crack_tip ;
            int tip_id ;
            int tip_indx ;
            CArbSmallSet<int,4> parent_splines ;
        } ;

        struct IntscPoint {
            CArbBdrySeg *spline ;
            double spline_par, line_par ;
            int pt_indx ;
        } ;

        struct DeleteSegData {
            CArbCoord2D end_coords[2] ;
            double end_radii[2] ;
            double end_radii_sqr[2] ;
            double line_coefs[4][3] ;
            bool degenerate_case ;
        } ;

        // private member variables

        CArbMshCrackRegion2D *Owner ;

        CArbArray<int> *CreatedElems ;
        // CArbHashTable<ArbEdgeKey,CArbMshCrackRegion2D::BoundaryData> *OldBdryTable ;
        // CArbHashTable<ArbEdgeKey,CArbMshCrackRegion2D::BoundaryData> *NewBdryTable ;

        CornerNodesSet *CornerNodes ;
        CArbMshTopo2D  *Boundary ; 
        MaterialHash   *MatTable ;
        int NumPoints ;
        int OldCrackTipNode ;
        int NumOldRings ;

        CArbArray<CArbMshCrackRegion2D::CrackTipData> *TipData ;

        CArbArray<int> *MouthNodes ;

        double *InsituSize ;

        bool               DoQuadAngleChecks ;
        double             MinQuadAngle ;
        double             MaxQuadAngle ;
        AngleCheckLocation QuadAngleLocation ;
        double             AnisotropicRatio ;

        void Initialize(const int num_tips,
                        const int *tip_ids,
                        const int *tip_indices,
                        const int num_points,
                        const CArbMshCrackRegion2D::CrackPt *const pts) ;
        void DeleteElements() ;
        void SetupRegions() ;
        void RemeshRegions() ;
        void Finalize() ;

        double CharElemSize(
                 const int elem,
                 CArbMshTopo2D *topo) const ;

        DeleteSegData *BuildDeleteData(const double *insitu_size) const ;

        CArbMshTopo2D *DoDeleteElements(int num_segs,
                                        DeleteSegData *del_segs,
                                        MaterialHash *mat_table,
                                        CornerNodesSet *corner_nodes,
                                        bool closed) ;

        void GenerateCrackSegment(
                 CArbMshCrackRegion2D::CrackPt *start_pt,
                 CArbMshCrackRegion2D::CrackPt *end_pt,
                 bool start_crack_flag,
                 double template_length_s,
                 double template_projection_s,
                 bool end_crack_flag,
                 double template_length_e,
                 double template_projection_e,
                 CArbArray<ArbMshNode> *nodes) ;

        void FillCharElemSize(
                 int num_pts,
                 CArbMshCrackRegion2D::CrackPt *pts,
                 double *insitu_size) const ;

//         void FillCrackMouthElemSize(
//                  int num_pts,
//                  CArbMshCrackRegion2D::CrackPt *pts,
//                  CArbMshCrackRegion2D::CrackTipData &tip_data) const ;

        int FindContainingElem(
                 const CArbCoord2D &point,
                 CArbMshTopo2D *topo) const ;

//         int FindMouthNode(const CArbCoord2D mouth,
//                           const int elem,
//                           CArbMshTopo2D *topo) const ;
// 
        int FindNearestNode(CArbCoord2D point) const ;

        int FindNumRings(int tip_id,
                         double *size) const ;

        void FillTipData(
                 int crack_tip_id,
                 CArbMshCrackRegion2D::CrackTipData &tip_data) const ;

        void GenerateCrackTipBdry(int num_tip_elems,
                    CArbMshCrackRegion2D::CrackTipData &tip_data,
                    CArbArray<int> &bdry) ;

        void GenerateCrackTipData(CArbCoord2D tip,
                                  CArbCoord2D norm,
                                  double tip_angle,
                                  CArbMshCrackRegion2D::CrackTipData &tip_data,
                                  int mat_id,
                                  CArbMshTopo2D *crack_bdry,
                                  CArbMshRegion2D *region,
                                  CArbArray<ArbMshNode> &temp_nodes,
                                  CArbArray<int> &temp_bdry,
                                  int *tip_node_id) ;

        void GenerateCrackTipElems(int num_tip_elems,
                       CArbMshCrackRegion2D::CrackTipData &tip_data,
                       CArbArray<ArbMshNode> &nodes,
                       int mat_id,
                       CArbArray<ArbMshElement2D> &elems) ;

//         ArbMshNode *GenerateCrackSide(
//                  bool internal_crack,
//                  int num_points,
//                  CArbMshCrackRegion2D::CrackPt *local_pts,
//                  CArbMshCrackRegion2D::CrackTipData &start_tip_data,
//                  CArbMshCrackRegion2D::CrackTipData &end_tip_data,
//                  int *num_nodes) ;

        void GenerateCrackTipNodes(
                         CArbCoord2D tip,
                         CArbCoord2D norm,
                         int num_tip_elems,
                         CArbMshCrackRegion2D::CrackTipData &tip_data,
                         double tip_angle,
                         CArbArray<ArbMshNode> &nodes) ;

        double GetAdjCharSize(int node) const ;

//         void RefineBoundarySeg(CArbRmshLoopIterator rmsh_iter,
//                                        CArbMshTopo2D *boundary,
//                                        int elem, 
//                                        bool two_points, 
//                                        CArbCoord2D close_pnt0,
//                                        CArbCoord2D close_pnt1,
//                                        double char_size0,
//                                        double char_size1,
//                                        double len_first,
//                                        double len_last,
//                                        int prev,
//                                        int next,
//                                        bool two_sides,
//                                        int oelem,
//                                        int oprev,
//                                        int onext,
//                                        int *node0,
//                                        int *node1) ;

        void RemoveTemplateNodes(int num_rings,
                                 int tip_id,
                                 CArbMshTopo2D *boundary) ;

        void RemoveTipNodes(int tip_id,
                            CArbMshTopo2D *boundary) ;

        bool SetTipData(
                 const int tip_id,
                 CArbMshCrackRegion2D::CrackPt *crack_pnt,
                 CArbMshCrackRegion2D::CrackTipData &tip_data) const ;

        bool ScanCross(const CArbCoord2D &b,const CArbCoord2D &i,
                       const CArbCoord2D &j) const ;


        // data structures

        // delete element factor

        double TriDeleteFactor ;
        double QuadDeleteFactor ;


        void AddCrackTemplates() ;
        void AddCrossingPoints() ;
        void AddFlawToBoundary() ;
        void AddNodesToBoundary() ;
        void ClassifyPoints() ;
        void FindAllIntersections() ;
        void FindCoincidentSplines() ;
        void FindParentSplines() ;
        void MeshRegions() ;
        void MovePointsToBoundary() ;
        void DuplicateBoundaryNodes() ;
        void HandleInterfaces() ;
        void RebuildRemeshTopology() ;
        void ResizeCrackTemplates() ;
        void SetFlawCharSizes() ;
        void PlotUpdate() ;
        bool CheckForIntersection(CArbBdrySeg *spline,
                  CArbCoord2D &pt0,CArbCoord2D &pt1) const ;
        void UpdateBdrySegCharLengths(bool tip_only) ;
        void UpdateCrackSegCharLengths() ;

        void PrintBoundary() ;


//        PointType *CrkPointTypes ;

        CArbHashTable<EdgeKey,CArbBdrySeg*> *BoundaryReps ;

        CArbArray<LocCrackPt> *LocalCrackPoints ;

        CArbHashTable<int,int> *IgnoreElems ;      

        double TemplateScaleUpFactor ;
        double TemplateScaleDownFactor ;
        double TemplateCharSizeFactor ;

        CArbCoord2D AdjacentTangent(int vtx,int adj) const ;

        bool CrossingInterface(const LocCrackPt &this_pt,
                               const LocCrackPt &prev_pt,
                               const LocCrackPt &next_pt) const ;

        int FindAdjElem(int node_id,
                        const CArbCoord2D &pt0,
                        const CArbCoord2D &pt1) const ;

        int FindCWvtx(int node_id,
                      const CArbCoord2D &pt0,
                      const CArbCoord2D &pt1) ;

        int FindDplElem(int node_id,
                        const CArbCoord2D &cur,
                        const CArbCoord2D &prev,
                        const CArbCoord2D &next,
                        CArbCoord2D *perp) const ;

        int AddIntersectionPoint(
                          CArbBdrySeg *spline,
                          double spline_par,double line_par,
                          const LocCrackPt &prev,
                          const LocCrackPt &next,
                          CArbArray<LocCrackPt> *NewPoints) ;

        PointType ClassifyEdge(int nd0,int nd1,int *other_elem) const ;

        double CharLength(int node,int other) const ;

        void DivideSplineEdge(int start,int stop,
                              CArbBdrySeg *spline) ;

        int SplitFlawNode(const LocCrackPt &this_pt,
                          int start_elem,int stop_elem,
                          const CArbCoord2D &tang) ;

        void SplitSpline(CArbBdrySeg *spline,
                         int new_node,double split_par,
                         double new_char_len,
                         CArbMshTopo2D  *Boundary) ;

        void UpdateBRepsForTemplates(
                 int tip_id,int rmt_id,int new_id,
                 double radius,double new_char_len) ;

        void FindClosestSplineData(
                 CArbCoord2D &lcoord,double *min_dist,
                 double *min_u,CArbCoord2D *min_coord,
                 CArbBdrySeg **min_spline,int *min_end0,
                 int *min_end1,int *min_pspl) const ;

        int DuplicateInterfaceNode(int existing_idx) ;
        void DuplicateSplineEndNode(LocCrackPt &existing,int far_node,
                            int new_node) ;


//         CArbRmshLoop *FindBoundarySegments(
//                              const int elem,
//                              const CArbMshTopo2D *rmsh,
//                              const CArbMshTopo2D *orig,
//                              const CArbElemSet *elems) const ;

        int *FindTipNodeList(
                       const int tip_node_id,
                       int *num_nodes,
                       CArbMshTopo2D *topo) const ;

        void GenerateNewBdryData(
                 CArbMshTopo2D *boundary,
                 CornerNodesSet *corner_node) ;

        void GenTipNodesHelp(int ring,int num,
                             int a0,int a1,int a2,
                             int b0,int b1,int b2, 
                             bool qtrpt_flag,
                             CArbArray<ArbMshNode> &nodes) ;

        CArbCoord2D PtUpdate(CArbCoord2D pt,
                             CArbCoord2D dir,
                             double mag,
                             CArbArray<ArbMshNode> &nodes) ;

        void RebuildNodeTable() ;

        void UpdateBoundary(CArbMshTopo2D *boundary,
                            const int loop_id,
                            const int nnode,
                            const int *nodes,
                            CArbHashTable<int,int> *mat_table) ;


        void DebugPrintBoundaryInfo() ;

} ;

inline int operator == (const CArbRmshRegion2D::EdgeKey &edge0,
                        const CArbRmshRegion2D::EdgeKey &edge1)
{
    return((edge0.id0 == edge1.id0) && (edge0.id1 == edge1.id1)) ;
}

inline int ArbHashIndex(const CArbRmshRegion2D::EdgeKey &edge)
{
    return(edge.id0 + edge.id1) ;
} ;
 
inline int ArbHashIndex(const CArbRmshRegion2D::EdgeFaceKey &key)
{
    return(key.vtx + key.face) ;
} ;
 
inline int operator == (const CArbRmshRegion2D::EdgeFaceKey &key0,
                        const CArbRmshRegion2D::EdgeFaceKey &key1)
{
    return((key0.vtx == key1.vtx) && (key0.face == key1.face)) ;
}

#endif
