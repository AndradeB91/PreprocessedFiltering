//
// ArbRmshRegion2D implementation file
//
// Description -
//   This is the implementation file for the
//   ArbRmshRegion2D objects.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2000,2001,2002
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.34 $  $Date: 2002/09/04 20:51:17 $  $Author: wash $
//

// TipZoneList - gives the zone id for each tip
// TipIndicies - gives the index in the crack points for each tip
// LocalTips - local array of crack tip coordinates
// TipData - array of crack tip data


#include "ArbRmshRegion2D.hpp"
//#include "ArbMshCrackRegion2D.hpp"
#include "ArbMshRegion2D.hpp"
//#include "ArbArray.hpp"
//#include "ArbBSpline2D.hpp"
#include "ArbSet.cpp"
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifdef MEMDEBUG
#include "MemDbg.hpp"
#define new new(__FILE__,__LINE__)
#endif

#define PI      3.14159265359
#define HALF_PI 1.57079632680
#define TWO_PI  6.28318530718

//#define NO_NODE 100000000


double LinearSpacing(const double segment,
                     const double length,
                     const int num,
                     const double ratio) ;

double LinearPosition(const double segment,
                      const double length,
                      const int num,
                      const double ratio) ;

int LinearNumber(const double length,
                 const double length_first,
                 const double length_last) ;

//#define TRI_DELETE_FACTOR  2.25
//#define QUAD_DELETE_FACTOR 3.0
#define DEFAULT_TRI_DELETE_FACTOR  2.25
//#define QUAD_DELETE_FACTOR 3.0
#define DEFAULT_QUAD_DELETE_FACTOR 3.25


// This factor is used when adding crack templates.  If the
// template radius is greater than this percentage of the
// segment length the template size will be scaled up to
// be the full segment length.

#define DEFAULT_TEMPLATE_SCALE_UP_FACTOR 0.9

// This factor is used when adding crack templates.  If the
// template radius is greater than this percentage of the
// segment length the temple size will be scaled down for a
// more even distribution of elements in the segment.

#define DEFAULT_TEMPLATE_SCALE_DOWN_FACTOR 0.7

// This factor is used when adding crack templates.  When a
// template is added a characteristic element length must be
// set for the points on the crack face at the template radius.
// The characteristic length of this element is set to be some
// value between the template radius and the length of the element
// in the outer most ring of the template.  A value of 0.0 sets
// this to be the size of the template radius.  A value of 1.0
// sets it to be the outermost element size.

#define DEFAULT_TEMPLATE_CHAR_SIZE_FACTOR 0.5

/* ------------------------------------------------------------

    Static helper functions

*/

#define LIN_TRI_COL_TRI  0
#define LIN_TRI_COL_QUAD 1
#define LIN_TRI_UNCOL    2
#define LIN_QUAD         3
#define QDC_TRI_COL_TRI  4
#define QDC_TRI_COL_QUAD 5
#define QDC_TRI_UNCOL    6
#define QDC_QUAD         7

static int ClassifyTipType(ArbMshOrder order,
                           CArbMshCrackRegion2D::ElemShape shape,
                           CArbMshCrackRegion2D::ElemType type,
                           bool constrained)
{
    if (order == LINEAR) {
        if (shape == CArbMshCrackRegion2D::S_TRIANGLE)
            if (constrained)
                if (type == CArbMshCrackRegion2D::T_TRIANGLE)
                    return(LIN_TRI_COL_TRI) ;
                else
                    return(LIN_TRI_COL_QUAD) ;
            else
                return(LIN_TRI_UNCOL) ;
        else
            return(LIN_QUAD) ;
    } else {
        if (shape == CArbMshCrackRegion2D::S_TRIANGLE)
            if (constrained)
                if (type == CArbMshCrackRegion2D::T_TRIANGLE)
                    return(QDC_TRI_COL_TRI) ;
                else
                    return(QDC_TRI_COL_QUAD) ;
            else
                return(QDC_TRI_UNCOL) ;
        else
            return(QDC_QUAD) ;
    }
}

static void FillElem(int id,int mat_id,
                     int num, int nd0,
                     int nd1, int nd2,
                     int nd3, ArbMshElement2D *elem)
{
    elem->elem_id = id ;
    elem->mat_id = mat_id ;
    elem->num_nodes = num ;
    elem->nodes[0] = nd0 ;
    elem->nodes[1] = nd1 ;
    elem->nodes[2] = nd2 ;
    elem->nodes[3] = nd3 ;
}

static void FillQuadElem(int id,int mat_id,
                         int num, int nd0,
                         int nd1, int nd2,
                         int nd3, int nd4,
                         int nd5, int nd6,
                         int nd7, ArbMshElement2D *elem)
{
    elem->elem_id = id ;
    elem->mat_id = mat_id ;
    elem->num_nodes = num ;
    elem->nodes[0] = nd0 ;
    elem->nodes[1] = nd1 ;
    elem->nodes[2] = nd2 ;
    elem->nodes[3] = nd3 ;
    elem->nodes[4] = nd4 ;
    elem->nodes[5] = nd5 ;
    elem->nodes[6] = nd6 ;
    elem->nodes[7] = nd7 ;
}

static void DoLinQuad(bool side,int id,int mat,
                      CArbArray<ArbMshNode> &nodes,
                      CArbArray<ArbMshElement2D> &elem,
                      int *ocur,int *icur)
{
    ArbMshElement2D lelem ;

    if (side) {                                 // side element
        FillElem(id,mat,4,nodes[*icur].id,
                 nodes[*ocur].id,nodes[*ocur+1].id,
                 nodes[*icur+1].id,&lelem) ;
        elem.InsertAtEnd(lelem) ;
        ++(*ocur) ;  ++(*icur) ;
    } else {                                    // corner element
        FillElem(id,mat,4,nodes[*icur].id,
                 nodes[*ocur].id,nodes[*ocur+1].id,
                 nodes[*ocur+2].id,&lelem) ;
        elem.InsertAtEnd(lelem) ;
        *ocur += 2 ;
    }
}

static void DoQdcQuad(bool side,int id,int mat,
                      CArbArray<ArbMshNode> &nodes,
                      CArbArray<ArbMshElement2D> &elem,
                      int *ocur,
                      int *icur,int *s1cur,
                      int *s2cur,int *s3cur)
{
    ArbMshElement2D lelem ;

    if (side) {                                 // side element
        FillQuadElem(id,mat,8,
                     nodes[*icur].id,nodes[*ocur].id,
                     nodes[*ocur+1].id,nodes[*icur+1].id,
                     nodes[*s1cur].id,nodes[*s2cur].id,
                     nodes[*s1cur+1].id,nodes[*s3cur].id,
                     &lelem) ;
        elem.InsertAtEnd(lelem) ;
        ++(*ocur) ;  ++(*icur) ;
        ++(*s1cur) ;  ++(*s2cur) ;  ++(*s3cur) ;
    } else {                                    // corner element
        FillQuadElem(id,mat,8,
                     nodes[*icur].id,nodes[*ocur].id,
                     nodes[*ocur+1].id,nodes[*ocur+2].id,
                     nodes[*s1cur].id,nodes[*s2cur].id,
                     nodes[*s2cur+1].id,nodes[*s1cur+1].id,
                     &lelem) ;
        elem.InsertAtEnd(lelem) ;
        *ocur += 2 ;
        ++(*s1cur) ;  *s2cur += 2 ;  // ++(*s3cur) ;
    }
}

static bool HasCommonParent(CArbSmallSet<int,4> &s0,
                            CArbSmallSet<int,4> &s1)
{
    for (int i=0 ; i<s0.NumElements() ; ++i) {
        if (s1.HasElement(s0.Element(i))) return(true) ;
    }
    return(false) ;
}

extern void display_topo(CArbMshTopo2D *topo,
             CArbHashTable<int,CArbCoord2D> *NodeTable) ;


void display_topo_edge(CArbMshTopo2D *topo,
             CArbHashTable<int,CArbCoord2D> *NodeTable)
{
    int num_nodes = topo->NumNodes() ;
    int *nodes = topo->GetNodeList() ;

    for (int i=0 ; i<num_nodes ; ++i) {
        CArbTopoAdjVtxIterator iter(topo,nodes[i]) ;
        CArbCoord2D *nd0 = NodeTable->Fetch(nodes[i]) ;
//        printf("t %g %g %d\n",(*nd0)[0],(*nd0)[1],nodes[i]) ;

        while (iter.More()) {
            CArbCoord2D *nd1 = NodeTable->Fetch(iter.AdjVtx()) ;
//   if (iter.CcwElem() == 128) {
            printf("l %g %g %g %g\n",
                   (*nd0)[0],(*nd0)[1],(*nd1)[0],(*nd1)[1]) ;

            printf("t %g %g %d\n",
                   (*nd0)[0],(*nd0)[1],nodes[i]) ;
            printf("t %g %g %d\n",
                   (*nd1)[0],(*nd1)[1],iter.AdjVtx()) ;


//             double x = 0.5 * ((*nd0)[0] + (*nd1)[0]) ;
//             double y = 0.5 * ((*nd0)[1] + (*nd1)[1]) ;
//             if (iter.CcwElem() != NO_ELEM) {
//                 printf("t %g %g %d\n",x,y,iter.CcwElem()) ;
//             } else {
//                 printf("t %g %g     =\n",x,y) ;
//             }

//            printf("t %g %g %d\n",(*nd1)[0],(*nd1)[1],iter.AdjVtx()) ;
//            if (iter.CcwElem() != NO_ELEM)
//                fprintf(stderr,"Elem: %d\n",iter.CcwElem()) ;
//    }
            ++iter ;
        }
    }
    printf("n\n") ;
    fflush(stdout) ;

    delete [] nodes ;
}

#if 0
static void DoPrint(CArbMshTopo2D *msh_topo,bool /*do_labels*/,
             CArbHashTable<int,CArbCoord2D> *NodeTable)
{
    int num_node = msh_topo->NumNodes() ;
    int *nodes = msh_topo->GetNodeList() ;

    for (int i=0 ; i<num_node ; ++i) {

        CArbTopoAdjVtxIterator iter(msh_topo,nodes[i]) ;

        for (iter.First() ; iter.More() ; ++iter) {
            if (iter.AdjVtx() > nodes[i]) {
                CArbCoord2D *p_nd0 = NodeTable->Fetch(iter.AdjVtx()) ;
                CArbCoord2D *p_nd1 = NodeTable->Fetch(nodes[i]) ;
                printf("l %g %g %g %g\n",
                       (*p_nd0)[0],(*p_nd0)[1],
                       (*p_nd1)[0],(*p_nd1)[1]) ;
                printf("t %g %g %d\n",(*p_nd0)[0],(*p_nd0)[1],iter.AdjVtx()) ;
                printf("t %g %g %d\n",(*p_nd1)[0],(*p_nd1)[1],nodes[i]) ;
            }
        }
        // if (do_labels) {
        //     CArbCoord2D *p_nd0 = NodeTable->Fetch(nodes[i]) ;
        //     printf("t %g %g %d\n",(*p_nd0)[0],(*p_nd0)[1],nodes[i]) ;
        // }
    }
    printf("n\n") ;
    fflush(stdout) ;
    delete [] nodes ;
}
#endif

// ----------------------------------------------------------------
// ----------------------------------------------------------------
// ----------------------------------------------------------------
//  These routines are used to compute a distribution of nodes
//  along an edge so that there is a linear progression in the
//  size of the element edges.
//
//  Given an edge of total length "length", and number of segments
//  "num", and the ratio of the last segment length divided by the
//  the first segment length "ratio", then the length of the ith
//  segment (i=1..n) is l_seg = F*(i*(ratio-1)/(num-1) + 1.0).
//  The length from the beginning to the end of the jth segment
//  (j=0..n-1) is Len = F*((j+1) + j*(j+1)(ratio-1)/2(num-1)).
//  in both cases F = length/(num * (1.0 + (ratio-1)/2)).

// LinearSpacing gives the length of the ith segment, where i = 1..n

double LinearSpacing(const double segment,
                            const double length,
                            const int num,
                            const double ratio)
{
    if (num == 1) return(length) ;
    double r = 1.0 / ratio ;
    double fact = length / (num + 0.5*num*(r-1.0)) ;
    return (fact * (1.0 + (segment-1)*(r-1.0)/(num-1.0))) ;
}

// LinearPosition gives the total length from the beginning of
// the edge to the end of segment i, where i = 0..n-1.

double LinearPosition(const double segment,
                             const double length,
                             const int num,
                             const double ratio)
{
    if (num == 1) return(length) ;
    double r = 1.0 / ratio ;
    double fact = length / (num + 0.5*num*(r-1.0)) ;
    return (fact * ((segment+1) + 0.5*segment*(segment+1)*
                    (r-1.0)/(num-1.0))) ;
}

// Given the total length of and an edge and the lengths of the
// first and last segments, this routine computes the approximate
// number of segments along the edge.

#define ROUND_FACTOR 0.1
int LinearNumber(const double length,
                                 const double length_first,
                                 const double length_last)
{
    double num = 2.0*length/(length_first+length_last) ;
    double intpart ;
    double fracpart = modf(num,&intpart) ;
    return(fracpart < ROUND_FACTOR ? int(intpart) : int(intpart)+1) ;
}

// If we know the
// total length of an edge, the length of the first element, an
// the ratio between the size of the first element and the second,
// then we can compute (due to some algebraic manipulations in
// Maple) the approximate number of element edge segments that will
// be placed along that edge.  Furthermore, once we know the the
// number of segments (appropriately rounded) we can compute the
// ratio between the first and last segment.

static double ApproxNumSegs(double first_len,double first_ratio,
                            double total_len)
{
    double t1 = first_len * (3.0*first_ratio - 1.0) ;
    double t2 = first_len*first_len * (9.0*first_ratio*first_ratio -
                6.0*first_ratio + 1.0) +
                8.0*total_len*first_ratio*first_len*(1.0-first_ratio) ;
    double t3 = 2.0*first_len*(first_ratio-1.0) ;
    return((t1-sqrt(t2))/t3) ;
}


static double ApproxRatio(int num,double first_ratio)
{
    return(-(-2.0*first_ratio + num*first_ratio - num + 1.0)/
           first_ratio) ;
}



// ----------------------------------------------------------------
// ----------------------------------------------------------------

// %(CArbRmshRegion2D::CArbRmshRegion2D-constructor-|-CArbMshCrackRegion2D-|*)
/* ++ ----------------------------------------------------------
**
**    CArbRmshRegion2D - constructor
**
**      CArbRmshRegion2D(CArbMshCrackRegion2D *owner)
**
**        owner - (in)  associated ArbMshCrackRegion2D object
**
**      Description: This is a constructor for a RmshRegion2D
**
**
** -- */

CArbRmshRegion2D::CArbRmshRegion2D(CArbMshCrackRegion2D *owner)
{
    Owner = owner ;

    // create some data structures used during the remesh

    CreatedElems  = new CArbArray<int> ;
    TipData       = new CArbArray<CArbMshCrackRegion2D::CrackTipData> ;
    MatTable      = new MaterialHash ;

    TriDeleteFactor = DEFAULT_TRI_DELETE_FACTOR ;
    QuadDeleteFactor = DEFAULT_QUAD_DELETE_FACTOR ;

    InsituSize = 0 ;
    CornerNodes = 0 ;
    Boundary = 0 ;

    MouthNodes = 0 ;
    BoundaryReps = 0 ;
    LocalCrackPoints = 0 ;
    IgnoreElems = 0 ;
    OldCrackTipNode = -1 ;

    TemplateScaleUpFactor = DEFAULT_TEMPLATE_SCALE_UP_FACTOR ;
    TemplateScaleDownFactor = DEFAULT_TEMPLATE_SCALE_DOWN_FACTOR ;
    TemplateCharSizeFactor = DEFAULT_TEMPLATE_CHAR_SIZE_FACTOR ;

    DoQuadAngleChecks = false ;

    DebugDisplayFlags = 0 ;
}




// %(CArbRmshRegion2D::CArbRmshRegion2D-destructor-virtual|~)
/* ++ ----------------------------------------------------------
**
**    CArbRmshRegion2D - destructor
**
**      ~CArbRmshRegion2D()
**
**      Description: This is a destructor for a RmshRegion2D.
**
**
** -- */

CArbRmshRegion2D::~CArbRmshRegion2D()
{
    // delete some data structures

    delete CreatedElems ;
    delete TipData ;
    delete MatTable ;

    if (InsituSize != 0)  delete [] InsituSize ;
    if (CornerNodes != 0) delete CornerNodes ;
    if (Boundary != 0)    delete Boundary ;
    if (MouthNodes != 0)  delete MouthNodes ;
    if (LocalCrackPoints != 0) delete LocalCrackPoints ;
    if (IgnoreElems != 0) delete IgnoreElems ;
    if (BoundaryReps != 0) {
        CArbHashTableIterator<EdgeKey,CArbBdrySeg*> iter(BoundaryReps) ;
        for (iter.First() ; iter.More() ; ++iter) {
            CArbBdrySeg *tmp = *(iter.Entry()) ;
            delete tmp ;
        }
        delete BoundaryReps ;
    }
}

// -------------------------------------------------------------

void CArbRmshRegion2D::AddFlaw(
                    const int num_tips,
                    const int *tip_ids,
                    const int *tip_indices,
                    const int num_points,
                    const CArbMshCrackRegion2D::CrackPt *const pts)
{
    Initialize(num_tips,tip_ids,tip_indices,num_points,pts) ;
    DeleteElements() ;
        // PlotUpdate() ;
        // DoPrint(Boundary,true,Owner->NodeTable) ;
        // fflush(stdout) ;
    ClassifyPoints() ;
    RebuildRemeshTopology() ;
        // PlotUpdate() ;
        // printf("n\n") ;
        // fflush(stdout) ;
    FindParentSplines() ;
        // PlotUpdate() ;
        // printf("n\n") ;
        // fflush(stdout) ;
    MovePointsToBoundary() ;
        // PlotUpdate() ;
        // printf("n\n") ;
        // fflush(stdout) ;
    DuplicateBoundaryNodes() ;
        // PlotUpdate() ;
        // printf("n\n") ;
        // fflush(stdout) ;
    HandleInterfaces() ;
        // PlotUpdate() ;
        // printf("n\n") ;
        // fflush(stdout) ;
    FindAllIntersections() ;
    AddFlawToBoundary() ;
        // PlotUpdate() ;
        // printf("n\n") ;
        // fflush(stdout) ;
    ResizeCrackTemplates() ;
    UpdateCrackSegCharLengths() ;
        // PlotUpdate() ;
        // printf("n\n") ;
        // fflush(stdout) ;
    UpdateBdrySegCharLengths(true) ;
    UpdateBdrySegCharLengths(false) ;
        // PlotUpdate() ;
        // printf("n\n") ;
        // fflush(stdout) ;
    AddCrackTemplates() ;
    AddNodesToBoundary() ;
    MeshRegions() ;
    Finalize() ;
}


void CArbRmshRegion2D::Initialize(
                    const int num_tips,
                    const int *tip_ids,
                    const int *tip_indices,
                    const int num_points,
                    const CArbMshCrackRegion2D::CrackPt *const pts)
{
    NumPoints = num_points ;
    CArbMshCrackRegion2D::CrackPt *local_pts =
        new CArbMshCrackRegion2D::CrackPt[num_points] ;
    for (int j=0 ; j<num_points ; ++j) local_pts[j] = pts[j] ;
    InsituSize = new double[num_points] ;

    // rebuild the node table so that we do not consider any
    // notes not attached to an element

    RebuildNodeTable() ;

    // update stuff for the crack tips

    for (int i=0 ; i<num_tips ; ++i) {

        CArbMshCrackRegion2D::CrackTipData tip_data ;
        SetTipData(tip_ids[i],&local_pts[tip_indices[i]],
                   tip_data) ;
        TipData->InsertAtEnd(tip_data) ;

        // Go through the owner's TemplateElement array and
        // remove elements associated with this tip

        CArbHashTableIterator<int,int> iter(Owner->TemplateElems) ;
        while (iter.More()) {
            if (*(iter.Entry()) == tip_ids[i]) {
                Owner->TemplateElems->Remove(iter.Key()) ;
                iter.First() ;
            } else {
                ++iter ;
            }
        }
    }

    // if we have an existing crack tip then determine the
    // number of rings and set the characteristic size

    if (OldCrackTipNode != -1) {
        double csize ;
        NumOldRings = FindNumRings(OldCrackTipNode,&csize) ;
        if (!local_pts[0].has_char_size) {
            local_pts[0].char_elem_size = csize ;
            local_pts[0].has_char_size = true ;
        }
    }

    FillCharElemSize(num_points,local_pts,InsituSize) ;

    // get a local copy of the crack points

    LocalCrackPoints = new CArbArray<LocCrackPt> ;
    for (int l=0 ; l<num_points ; ++l) {
        LocCrackPt loc_pt ;
        loc_pt.coord = CArbCoord2D(local_pts[l].coord[0],
                                   local_pts[l].coord[1]) ;
        loc_pt.char_elem_size = local_pts[l].char_elem_size ;
        loc_pt.boundary_tol = -1.0 ;
        loc_pt.respect_char_size = local_pts[l].has_char_size ;
        loc_pt.node_id = -1 ;
        loc_pt.elem = -1 ;
        loc_pt.type = INSIDE ;
        loc_pt.crack_tip = false ;
        loc_pt.tip_id = -1 ;
        loc_pt.tip_indx = -1 ;
        LocalCrackPoints->InsertAtEnd(loc_pt) ;
    }

    // store the crack tip data for the local points

    for (int k=0 ; k<num_tips ; ++k) {
        LocalCrackPoints->At(tip_indices[k]).crack_tip = true ;
        LocalCrackPoints->At(tip_indices[k]).tip_id = tip_ids[k] ;
        LocalCrackPoints->At(tip_indices[k]).tip_indx = k ;
    }

    delete [] local_pts ;
}


// %(CArbRmshRegion2D::Finalize-void-virtual|)
/* ++ ----------------------------------------------------------
**
**    Finalize - finalize the remeshing
**
**      virtual void Finalize()
**
**      Description: This function rebuilds the table of used nodes and
**          generates information about the new element id's adjacent
**          to remeshed boundies.
**
**
** -- */

void CArbRmshRegion2D::Finalize()
{
    int i,j ;

    // First check to see if any of the crack tips are not
    // constrained.  If this is the case then we need to
    // delete all the crack tip nodes and add quad nodes
    // at the tip

    for (i=0 ; i<NumPoints ; ++i) {
        if (!LocalCrackPoints->At(i).crack_tip) continue ;

        int tip_id = LocalCrackPoints->At(i).node_id ;
        CArbMshCrackRegion2D::CrackTipData &tip_data =
                TipData->At(LocalCrackPoints->At(i).tip_indx) ;
        if (tip_data.const_tip_flag) continue ;

        CArbCoord2D *tip_coord = Owner->NodeTable->Fetch(tip_id) ;

        // here we have an unconstrained crack tip
        // loop through the edges adjacent to the crack
        // tip and determine which nodes are on elements

        // get memory to store the new element connectivities

        int num_adj = Owner->MshTopo->NumAdjElems(tip_id) ;
        int *tip_conn = new int[num_adj*10] ;
        int cur = 0 ;

        int first_node_id = tip_id ;

        CArbTopoAdjVtxIterator iter(Owner->MshTopo,tip_id) ;
        for (iter.First() ; iter.More() ; ++iter) {
            if (iter.CcwElem() == NO_ELEM) break ;

            int num,nodes[8] ;
            Owner->MshTopo->GetElemNodes(iter.CcwElem(),
                                         tip_id,&num,nodes) ;

            // if this is not a triangular element the bail

            if ((num != 3) && (num != 6)) continue ;

            // create the new crack tip node (and side node if
            // necessary and update the element list

            ArbMshNode new_node = Owner->NewNode(tip_coord->x(),
                                                 tip_coord->y()) ;
            Owner->NodeTable->Store(new_node.id,new_node.coord) ;
            CornerNodes->Store(new_node.id,1) ;

            if (num == 3) {
                nodes[0] = first_node_id ;
                nodes[3] = new_node.id ;
                num = 4 ;
            } else {
                ArbMshNode new_snode = Owner->NewNode(tip_coord->x(),
                                                      tip_coord->y()) ;
                Owner->NodeTable->Store(new_snode.id,new_snode.coord) ;
                nodes[0] = first_node_id ;
                nodes[6] = new_node.id ;
                nodes[7] = new_snode.id ;
                num = 8 ;
            }

            first_node_id = new_node.id ;

            tip_conn[cur*10] = iter.CcwElem() ;
            tip_conn[cur*10+1] = num ;
            for (int j=0 ; j<num ; ++j) {
                tip_conn[cur*10+j+2] = nodes[j] ;
            }

            // get a pointer to this element in the element table
            // and update the element for the new nodes

            ArbMshElement2D *elem =
                Owner->ElemTable->Fetch(iter.CcwElem()) ;
            elem->num_nodes = num ;
            if (num == 4) {
                for (int j=0 ; j<num ; ++j) {
                    elem->nodes[j] = nodes[j] ;
                }
            } else {
                elem->nodes[0] = nodes[0] ;
                elem->nodes[1] = nodes[2] ;
                elem->nodes[2] = nodes[4] ;
                elem->nodes[3] = nodes[6] ;
                elem->nodes[4] = nodes[1] ;
                elem->nodes[5] = nodes[3] ;
                elem->nodes[6] = nodes[5] ;
                elem->nodes[7] = nodes[7] ;
            }
            ++cur ;
        }

        // delete the old elements from the topology

        for (j=0 ; j<cur ; ++j) {
            Owner->MshTopo->DeleteElement(tip_conn[j*10]) ;
        }

        // add the updated to the topology

        for (j=0 ; j<cur ; ++j) {
            Owner->MshTopo->InsertCollapsedElement(
                               tip_conn[j*10],tip_conn[j*10+1],
                               &tip_conn[j*10+2]) ;
        }

        delete tip_conn ;
    }

    // redo the node table boundary

    RebuildNodeTable() ;
    GenerateNewBdryData(Boundary,CornerNodes) ;
}


// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------
//  ArbRmshGeneralFlaw2D
// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------


void CArbRmshRegion2D::DeleteElements()
{
    // build the array of data associated with the
    // delete regions

    DeleteSegData *delete_data = BuildDeleteData(InsituSize) ;

    // delete appropriate nodes and elements

    CornerNodes = new CornerNodesSet ;
    Boundary = DoDeleteElements(NumPoints,delete_data,
                                MatTable,CornerNodes,true) ;
    delete [] delete_data ;

    // if this is crack growth we need to deal with the nodes
    // associated with the old tip

    if (OldCrackTipNode != -1) {

        // look at the old crack tip and delete extra nodes that
        // may have been inserted if this was a collapsed but
        // not constrained crack tip

        RemoveTipNodes(OldCrackTipNode,Boundary) ;

        // look at the crack tip and delete extra nodes that
        // were inserted as part of the crack tip template

        RemoveTemplateNodes(NumOldRings,OldCrackTipNode,Boundary) ;
    }
}


void CArbRmshRegion2D::ClassifyPoints()
{
    // look through all the input crack points and classify
    // each as being Inside the delete region, Outside the
    // delete region, on the external boundary of the region
    // or on a bi-material interface in the boundary

    for (int i=0 ; i<NumPoints ; ++i) {

        CArbCoord2D cpoint = LocalCrackPoints->At(i).coord ;
        bool classified = false ;

        // find the distance from this point to the previous
        // and next crack points.  We store the smaller of these
        // which we will used as check to see if we should snap
        // to a boundary.

        CArbCoord2D next = LocalCrackPoints->At((i+1)%NumPoints).coord ;
        double next_dist = (cpoint-next).Magnitude() ;

        int j = i==0 ? NumPoints-1 : i-1 ;
        CArbCoord2D prev = LocalCrackPoints->At(j).coord ;
        double prev_dist = (cpoint-prev).Magnitude() ;

        double cdist = prev_dist < next_dist ? prev_dist : next_dist ;
        double min_b_tol = 0.5 * cdist ;

        CArbTopoElemIterator elem_iter(Boundary) ;
        for (elem_iter.First() ; elem_iter.More() ; ++elem_iter) {

            int num_cross = 0 ;
            CArbTopoEdgeOnElemIterator edge_iter(Boundary,*elem_iter) ;
            for (edge_iter.First() ; edge_iter.More() ; ++edge_iter) {

                CArbCoord2D *node0 = Owner->NodeTable->Fetch(edge_iter[0]) ;
                CArbCoord2D *node1 = Owner->NodeTable->Fetch(edge_iter[1]) ;

                // compute the perpendicular distance to this edge
                // to see if we are within the tolerance

                double elen ;
                double dist = PerpendicularDist(cpoint,*node0,
                                                *node1,&elen) ;

                double check = cdist < elen ? cdist : elen ;

                if ((dist >= 0) &&
                    (dist < check*Owner->BoundaryTolerance)) {

                    // we are in a boundary so determine if this is
                    // an external boundary or a bi-material interface

                    int elem_0, elem_1 ;
                    Boundary->GetElemsAboutEdge(edge_iter[0],edge_iter[1],
                                                &elem_0,&elem_1) ;

                    if ((elem_0 != NO_ELEM) && (elem_1 != NO_ELEM)) {

                         // if we are at the point were an interface
                         // meets a boundary, then make sure we
                         // classify the point as a boundary

                         if (!classified)
                             LocalCrackPoints->At(i).type = INTERFACE ;
                    } else {
                         LocalCrackPoints->At(i).type = BOUNDARY ;
                    }
                    LocalCrackPoints->At(i).boundary_tol =
                        (min_b_tol < elen*Owner->BoundaryTolerance) ?
                        min_b_tol : elen*Owner->BoundaryTolerance ;
                    classified = true ;
                } else {
                    if (ScanCross(cpoint,*node0,*node1)) {
                        ++num_cross ;
                    }
                }
            }

            if ((classified) &&
                ((LocalCrackPoints->At(i).type == INTERFACE) ||
                 (LocalCrackPoints->At(i).type == BOUNDARY))) break ;

            // check to see if there are an odd number of crossing
            // for this polygon.  If so we are inside

            if ((num_cross%2) == 1) {
                LocalCrackPoints->At(i).type = INSIDE ;
                LocalCrackPoints->At(i).elem = *elem_iter ;
                classified = true ;
            }
        }

        // if not classified as something else we must be outside

        if (!classified) LocalCrackPoints->At(i).type = OUTSIDE ;
    }
}

void CArbRmshRegion2D::FindClosestSplineData(
                 CArbCoord2D &lcoord,double *min_dist,
                 double *min_u,CArbCoord2D *min_coord,
                 CArbBdrySeg **min_spline,int *min_end0,
                 int *min_end1,int *min_pspl) const
{
    bool first = true ;

    CArbTopoElemIterator elem_iter(Boundary) ;
    for (elem_iter.First() ; elem_iter.More() ; ++elem_iter) {

        CArbTopoEdgeOnElemIterator edge_iter(Boundary,*elem_iter) ;
        for (edge_iter.First() ; edge_iter.More() ; ++edge_iter) {

            // check to see if we have a spline rep or not

            CArbBdrySeg **ptr = BoundaryReps->Fetch(
                EdgeKey(edge_iter[0],edge_iter[1])) ;

            if (ptr != 0) {
                double u ;
                CArbBdrySeg *spline = *ptr ;
                if (spline->ClosestPoint(lcoord,&u)) {
                    CArbCoord2D cpoint = spline->Evaluate(u) ;
                    double dist = (lcoord-cpoint).Magnitude() ;
                    if (first || (dist < *min_dist)) {
                        first = false ;
                        *min_dist = dist ;
                        *min_coord = cpoint ;
                        *min_spline = spline ;
                        *min_u = u ;
                        *min_end0 = edge_iter[0] ;
                        *min_end1 = edge_iter[1] ;
                        *min_pspl = spline->GetParentId() ;
                    }
                }
            }
        }
    }
}



void CArbRmshRegion2D::FindParentSplines()
{
    CArbArray<int> p_ids ;

    // get a tolerance based on the distance between the crack pts

    double tol = (LocalCrackPoints->At(1).coord -
                  LocalCrackPoints->At(0).coord).Magnitude() * 0.000001 ;

    // loop through all the boundary and interface points and
    // find the id of the spline on which they lie

    for (int cur=0 ; cur<NumPoints ; ++cur) {
        if ((LocalCrackPoints->At(cur).type == BOUNDARY) ||
            (LocalCrackPoints->At(cur).type == INTERFACE)) {

            // search through the boundary splines and find the
            // intersection points

            bool first = true ;
            double min_dist = .0f;
            CArbCoord2D &lcoord = LocalCrackPoints->At(cur).coord ;

            CArbTopoElemIterator elem_iter(Boundary) ;
            for (elem_iter.First() ; elem_iter.More() ; ++elem_iter) {

                CArbTopoEdgeOnElemIterator edge_iter(Boundary,*elem_iter) ;
                for (edge_iter.First() ; edge_iter.More() ; ++edge_iter) {

                // check to see if we have a spline rep or not

                    CArbBdrySeg **ptr = BoundaryReps->Fetch(
                        EdgeKey(edge_iter[0],edge_iter[1])) ;

                    if (ptr != 0) {
                        double u ;
                        CArbBdrySeg *spline = *ptr ;
                        if (spline->ClosestPoint(lcoord,&u)) {
                            CArbCoord2D cpoint = spline->Evaluate(u) ;
                            double dist = (lcoord-cpoint).Magnitude() ;

                            if (first || (dist <= min_dist-tol)) {
                                first = false ;
                                p_ids.Clear() ;
                                p_ids.InsertAtEnd(spline->GetParentId()) ;
                                min_dist = dist ;
                            } else if (dist <= min_dist+tol) {
                                p_ids.InsertAtEnd(spline->GetParentId()) ;
                            }
                        }
                    }
                }
            }

            for (int i=0 ; i<p_ids.NumEntries() ; ++i)
                LocalCrackPoints->At(cur).parent_splines.Insert(p_ids[i]);
        }
    }
}



#define CLOSE_FACTOR 0.001

static int Coincident(CArbCoord2D *pt0_0,CArbCoord2D *pt0_1,
                       double char0_0,double char0_1,
                       CArbCoord2D *pt1_0,CArbCoord2D *pt1_1,
                       double char1_0,double char1_1)
{
    double tol = char0_0 < char1_0 ? char0_0 : char1_0 ;
    if (((*pt0_0)-(*pt1_0)).Magnitude() < tol*CLOSE_FACTOR) {
        double tol = char0_1 < char1_1 ? char0_1 : char1_1 ;
        if (((*pt0_1)-(*pt1_1)).Magnitude() < tol*CLOSE_FACTOR) {
            return(1) ;
        } else {
            return(0) ;
        }
    }

    tol = char0_0 < char1_1 ? char0_0 : char1_1 ;
    if (((*pt0_0)-(*pt1_1)).Magnitude() < tol*CLOSE_FACTOR) {
        double tol = char0_1 < char1_0 ? char0_1 : char1_0 ;
        if (((*pt0_1)-(*pt1_0)).Magnitude() < tol*CLOSE_FACTOR) {
            return(2) ;
        } else {
            return(0) ;
        }
    }

    return(0) ;
}


void CArbRmshRegion2D::FindCoincidentSplines()
{

    // it may be that some of the splines making up the boundary
    // are coincident.  Probably part of another crack.  We need
    // to identify these pares and make sure that they are
    // oriented in the same direction geometrically so that later
    // when we add nodes they will be in identical locations.

    // Also, we keep track of the mate splines so that if one
    // side of the spline wants to be refined the other will
    // be refined also.

    CArbArray<CArbBdrySeg *> coinc ;
    CArbArray<int> orient ;

    // this is an N**2 comparison.  Let's hope that N is small

    CArbTopoElemIterator elem_iter0(Boundary) ;
    for (elem_iter0.First() ; elem_iter0.More() ; ++elem_iter0) {

        CArbTopoEdgeOnElemIterator edge_iter0(Boundary,*elem_iter0) ;
        for (edge_iter0.First() ; edge_iter0.More() ; ++edge_iter0) {

            // check to see if we have a spline rep or not

            CArbBdrySeg **ptr0 = BoundaryReps->Fetch(
                        EdgeKey(edge_iter0[0],edge_iter0[1])) ;

            if (ptr0 != 0) {

                CArbCoord2D *pt0_0 = Owner->NodeTable->Fetch(
                                           (*ptr0)->GetStartNode()) ;
                CArbCoord2D *pt0_1 = Owner->NodeTable->Fetch(
                                           (*ptr0)->GetStopNode()) ;
                double char0_0 = (*ptr0)->GetStartCharLen() ;
                double char0_1 = (*ptr0)->GetStopCharLen() ;

                // now loop through again looking for coincident eges

                CArbTopoElemIterator elem_iter1(Boundary) ;
                for (elem_iter1.First() ; elem_iter1.More() ; ++elem_iter1) {

                    CArbTopoEdgeOnElemIterator edge_iter1(Boundary,*elem_iter1) ;
                    for (edge_iter1.First() ; edge_iter1.More() ;
                         ++edge_iter1) {

                        // check to see if we have a spline rep or not

                        CArbBdrySeg **ptr1 = BoundaryReps->Fetch(
                             EdgeKey(edge_iter1[0],edge_iter1[1])) ;

                        if ((ptr1 != 0) && (ptr1 != ptr0)) {

                            CArbCoord2D *pt1_0 = Owner->NodeTable->Fetch(
                                                 (*ptr1)->GetStartNode()) ;
                            CArbCoord2D *pt1_1 = Owner->NodeTable->Fetch(
                                                 (*ptr1)->GetStopNode()) ;
                            double char1_0 = (*ptr1)->GetStartCharLen() ;
                            double char1_1 = (*ptr1)->GetStopCharLen() ;

                            int val = Coincident(
                                          pt0_0,pt0_1,char0_0,char0_1,
                                          pt1_0,pt1_1,char1_0,char1_1) ;
                            if (val > 0) {
                                coinc.InsertAtEnd(*ptr0) ;
                                coinc.InsertAtEnd(*ptr1) ;
                                orient.InsertAtEnd(val) ;
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i=0 ; i<orient.NumEntries() ; ++i) {
        CArbBdrySeg *ptr0 = coinc[i*2] ;
        CArbBdrySeg *ptr1 = coinc[i*2+1] ;
        if (ptr0->GetMate() == 0) {
            ptr0->SetMate(ptr1) ;
            ptr1->SetMate(ptr0) ;
            if (orient[i] == 2) ptr1->Reverse() ;
        }
    }
}



#define CORNER_TOL 0.01

void CArbRmshRegion2D::MovePointsToBoundary()
{
    // At this point we have classified all the crack points
    // and have rebuilt the remesh region boundary with a
    // combination of splines and straight line segments.

    // first loop through the crack points and snap them to
    // the boundary segments, splitting if necessary

    int cur ;

    for (cur=0 ; cur<NumPoints ; ++cur) {
        if ((LocalCrackPoints->At(cur).type == BOUNDARY) ||
            (LocalCrackPoints->At(cur).type == INTERFACE)) {

            // search through the boundary splines and find the
            // intersection point

            double min_dist, min_u ;
            CArbCoord2D min_coord ;
            CArbBdrySeg *min_spline ;
            int min_end0,min_end1,min_pspl ;

            FindClosestSplineData(
                LocalCrackPoints->At(cur).coord,
                &min_dist,&min_u,&min_coord,&min_spline,
                &min_end0,&min_end1,&min_pspl) ;

            // now add this node to the boundary

            int new_node = -1 ;
            bool do_split = true ;

            // if we are at an end of the spline, snap to the end
            // as long as the end point is not an adjacent crack point

            double dsdu = min_spline->TangentMag(min_u) ;
            double corner_tol = dsdu != 0.0 ?
                LocalCrackPoints->At(cur).boundary_tol / dsdu :
                CORNER_TOL ;

            if (min_u < corner_tol) {
                int end_node = min_spline->GetStartNode() ;
                if (cur == NumPoints-1) {
                    do_split =
                      (LocalCrackPoints->At(cur-1).node_id == end_node) ||
                      (LocalCrackPoints->At(0).node_id == end_node) ;
                } else {
                    do_split =
                      LocalCrackPoints->At(cur-1).node_id == end_node ;
                }
                if (!do_split) {
                    min_u = 0.0 ;
                    min_coord = min_spline->Evaluate(0.0) ;
                    new_node = min_spline->GetStartNode() ;
                }
            }

            if (min_u > (1.0-corner_tol)) {
                int end_node = min_spline->GetStopNode() ;
                if (cur == NumPoints-1) {
                    do_split =
                      (LocalCrackPoints->At(cur-1).node_id == end_node) ||
                      (LocalCrackPoints->At(0).node_id == end_node) ;
                } else {
                    do_split =
                      LocalCrackPoints->At(cur-1).node_id == end_node ;
                }
                if (!do_split) {
                    do_split = false ;
                    min_u = 1.0 ;
                    min_coord = min_spline->Evaluate(1.0) ;
                    new_node = min_spline->GetStopNode() ;
                }
            }

            if (new_node < 0) new_node = Owner->NewNodeId() ;

            // update the point information

            LocalCrackPoints->At(cur).coord = min_coord ;
            LocalCrackPoints->At(cur).node_id = new_node ;
            Owner->NodeTable->Store(new_node,min_coord) ;
            CornerNodes->Store(new_node,1) ;

            // if necessary split the spline

            if (do_split) {
                SplitSpline(min_spline,new_node,min_u,
                       LocalCrackPoints->At(cur).char_elem_size,
                       Boundary) ;
            }
        }
    }
}


void CArbRmshRegion2D::DuplicateBoundaryNodes()
{
    int cur, prev, next ;

    // go through the points a second time and duplicate
    // nodes at the boundary if necessary

    CArbArray<LocCrackPt> *NewPoints = new CArbArray<LocCrackPt> ;

    for (cur=0 ; cur<NumPoints ; ++cur) {

        LocCrackPt &this_pt = LocalCrackPoints->At(cur) ;

        if ((this_pt.type == BOUNDARY) || (this_pt.type == INTERFACE)) {

            bool duplicate = false ;

            // if this is a boundary point check to see if we
            // need to duplicate it

            if (this_pt.type == BOUNDARY) {

                prev = (cur == 0) ? LocalCrackPoints->NumEntries()-1 : cur-1 ;
                next = (cur+1) % LocalCrackPoints->NumEntries() ;

                LocCrackPt &prev_pt = LocalCrackPoints->At(prev) ;
                LocCrackPt &next_pt = LocalCrackPoints->At(next) ;

                // ignore the case where we move from outside to
                // the boundary and back out again

                if ((prev_pt.type == OUTSIDE) &&
                    (next_pt.type == OUTSIDE)) {
                    if (prev != next)
                        continue ;
                    else
                        duplicate = true ;

                // if we move from inside to the boundary and back
                // inside we duplicate the node

                } else if ((prev_pt.type == INSIDE) &&
                           (next_pt.type == INSIDE)) {
                    duplicate = true ;

                // if we move from an interface to the boundary and back
                // to the interface we duplicate

                } else if ((prev_pt.type == INTERFACE) &&
                           (next_pt.type == INTERFACE)) {
                    duplicate = true ;

                // if we move from a boundary to the boundary and back
                // to a boundary we then we duplicate (note if there is
                // only one boundary involved, the duplicated node will
                // be removed in the next phase

                } else if ((prev_pt.type == BOUNDARY) &&
                           (next_pt.type == BOUNDARY)) {
                    if ((!HasCommonParent(prev_pt.parent_splines,
                                          this_pt.parent_splines) ||
                         !HasCommonParent(this_pt.parent_splines,
                                          next_pt.parent_splines)) &&
                         (BoundaryReps->Fetch(
                           EdgeKey(prev_pt.node_id,this_pt.node_id)) == 0)) {
                        duplicate = true ;
                    }
                }
            }

            NewPoints->InsertAtEnd(this_pt) ;

            // if this point is outside, make sure it is not a crack tip

            if (this_pt.type == OUTSIDE)
                NewPoints->At(-1).crack_tip = false ;


            if (duplicate) {

                // if we are duplicating a crack  tip, make it
                // not a crack tip

                NewPoints->At(-1).crack_tip = false ;
                this_pt.crack_tip = false ;

                // determine if the in comming and out going flaw
                // segments are in bimaterial interfaces

                prev = (cur == 0) ? LocalCrackPoints->NumEntries()-1 : cur-1 ;
                LocCrackPt &prev_pt = LocalCrackPoints->At(prev) ;

                bool p_intfc = false ;
                CArbBdrySeg *p_seg ;
                CArbBdrySeg **ptr = BoundaryReps->Fetch(
                      EdgeKey(this_pt.node_id,prev_pt.node_id)) ;

                if (ptr != 0) {
                    p_seg = *ptr ;
                    if (((p_seg->GetStartNode() == this_pt.node_id) &&
                         (p_seg->GetStopNode() == prev_pt.node_id)) ||
                        ((p_seg->GetStopNode() == this_pt.node_id) &&
                         (p_seg->GetStartNode() == prev_pt.node_id)))
                        p_intfc = true ;
                }

                next = (cur+1) % LocalCrackPoints->NumEntries() ;
                LocCrackPt &next_pt = LocalCrackPoints->At(next) ;

                bool n_intfc = false ;
                CArbBdrySeg *n_seg ;
                ptr = BoundaryReps->Fetch(
                          EdgeKey(this_pt.node_id,next_pt.node_id)) ;

                if (ptr != 0) {
                    n_seg = *ptr ;
                    if (((n_seg->GetStartNode() == this_pt.node_id) &&
                         (n_seg->GetStopNode() == next_pt.node_id)) ||
                        ((n_seg->GetStopNode() == this_pt.node_id) &&
                         (n_seg->GetStartNode() == next_pt.node_id)))
                        n_intfc = true ;
                }

                // now determine where to do the split

                CArbCoord2D tang ;
                bool after = true ;
                int start_elem, stop_elem ;

                if (p_intfc && n_intfc &&
                    (prev_pt.node_id == next_pt.node_id)) {
                    start_elem = Boundary->GetCCWElem(this_pt.node_id,
                                                      next_pt.node_id) ;
                    stop_elem = NO_ELEM ;
                    after = false ;

                    CArbCoord2D delt = prev_pt.coord - this_pt.coord ;
                    tang = CArbCoord2D(-delt.y(),delt.x()) ;

                } else {
                    start_elem = NO_ELEM ;
                    stop_elem = FindDplElem(this_pt.node_id,
                                             this_pt.coord,
                                             prev_pt.coord,
                                             next_pt.coord,
                                             &tang) ;
                }

                // now duplicate the node

                int new_node = SplitFlawNode(this_pt,start_elem,
                                             stop_elem,tang) ;

                if (MouthNodes == 0) MouthNodes = new CArbArray<int> ;
                MouthNodes->InsertAtEnd(this_pt.node_id) ;
                MouthNodes->InsertAtEnd(new_node) ;

                NewPoints->InsertAtEnd(this_pt) ;

                if (after) {
                    NewPoints->At(-1).node_id = new_node ;
                } else {
                    NewPoints->At(-2).node_id = new_node ;
                }
            }
        } else {
            NewPoints->InsertAtEnd(this_pt) ;
        }
    }

    delete LocalCrackPoints ;
    LocalCrackPoints = NewPoints ;
    NumPoints = LocalCrackPoints->NumEntries() ;
}


int CArbRmshRegion2D::SplitFlawNode(const LocCrackPt &this_pt,
                                    int start_elem,int stop_elem,
                                    const CArbCoord2D &tang)
{
    // This routine splits a node that makes up part of a flaw
    // and adds a new edge between the new and old nodes

    int new_node = Owner->NewNodeId() ;
    Owner->NodeTable->Store(new_node,this_pt.coord) ;
    CornerNodes->Store(new_node,1) ;

    // update boundary segments

    CArbTopoAdjVtxCyclicIterator iter(Boundary,this_pt.node_id) ;
    while (iter.CcwElem() != start_elem) ++iter ;
    ++iter ;

    while (1) {
        CArbBdrySeg *seg = BoundaryReps->FetchValue(
                      EdgeKey(this_pt.node_id,iter.AdjVtx())) ;
        if (seg->GetStartNode() == this_pt.node_id) {
            seg->SetStartNode(new_node) ;
        } else {
            seg->SetStopNode(new_node) ;
        }
        BoundaryReps->Remove( EdgeKey(this_pt.node_id,iter.AdjVtx())) ;
        BoundaryReps->Store(EdgeKey(new_node,iter.AdjVtx()),seg) ;
        if (iter.CcwElem() == stop_elem) {
            break ;
        }
        ++iter ;
    }

    Boundary->SplitVertexMakeEdge(this_pt.node_id,new_node,
                                  start_elem,stop_elem) ;

    CArbBdrySeg *nseg = new CArbBdrySeg(
                        this_pt.coord,tang,
                        this_pt.node_id,this_pt.char_elem_size,
                        new_node,this_pt.char_elem_size,
                        this_pt.parent_splines.Element(0),
                        (this_pt.type == INTERFACE)) ;

    BoundaryReps->Store(EdgeKey(this_pt.node_id,new_node),nseg) ;
    return(new_node) ;
}


void CArbRmshRegion2D::HandleInterfaces()
{
    int i, cur, prev, next ;

    // go through the points and look for bi-material
    // interface crossings or cracks moving along
    // bi-material interfaces.

    CArbArray<LocCrackPt> *NewPoints = new CArbArray<LocCrackPt> ;
    CArbArray<int> interface_seg ;

    for (cur=0 ; cur<NumPoints ; ++cur) {

        LocCrackPt &this_pt = LocalCrackPoints->At(cur) ;
        NewPoints->InsertAtEnd(this_pt) ;

        if (this_pt.type == INTERFACE) {

            // check to see if we have seen this node already

            bool found = false ;
            for (i=0 ; i<cur ; ++i) {
                if (this_pt.node_id ==
                    LocalCrackPoints->At(i).node_id) {
                    found = true ;
                    break ;
                }
            }

            // if not found go on to the next point

            if (!found) continue ;

            // determine if we need to split this node

            prev = (cur == 0) ? LocalCrackPoints->NumEntries()-1 : cur-1 ;
            next = (cur+1) % LocalCrackPoints->NumEntries() ;

            LocCrackPt &prev_pt = LocalCrackPoints->At(prev) ;
            LocCrackPt &next_pt = LocalCrackPoints->At(next) ;

            bool duplicate = true ;
            if ((prev_pt.type == INTERFACE) &&
                HasCommonParent(prev_pt.parent_splines,
                                this_pt.parent_splines)) {
                duplicate = false ;
                interface_seg.InsertAtEnd(
                      NewPoints->NumEntries()-1) ;
                interface_seg.InsertAtEnd(prev_pt.node_id) ;
                interface_seg.InsertAtEnd(this_pt.node_id) ;
            } else {
                duplicate = CrossingInterface(this_pt,prev_pt,next_pt) ;
            }

            if (duplicate) {
                int new_node = DuplicateInterfaceNode(i) ;
                NewPoints->At(-1).node_id = new_node ;
                Owner->NodeTable->Store(new_node,
                           NewPoints->At(-1).coord) ;
                CornerNodes->Store(new_node,1) ;
            }
        }
    }

    // if necessary add segments for cracks that are in bimaterial
    // interfaces.

    if (interface_seg.NumEntries() > 0) {

        // add the first segment

        int indx = interface_seg[0] ;
        prev = interface_seg[1] ;
        cur  = interface_seg[2] ;
        int new_node = Owner->NewNodeId() ;
        NewPoints->At(indx).node_id = new_node ;
        Owner->NodeTable->Store(new_node,NewPoints->At(indx).coord) ;
        CornerNodes->Store(new_node,1) ;

        int CW_vtx0 = Boundary->GetCWNode(prev,cur) ;
        Boundary->AddEdgeVertex(prev,CW_vtx0,new_node) ;

        CArbBdrySeg *old_bdry =
            BoundaryReps->FetchValue(EdgeKey(prev,cur)) ;
        CArbBdrySeg *new_bdry = new CArbBdrySeg(*old_bdry) ;
        if (old_bdry->GetStartNode() == prev) {
            new_bdry->SetStartNode(prev) ;
            new_bdry->SetStopNode(new_node) ;
        } else {
            new_bdry->SetStartNode(new_node) ;
            new_bdry->SetStopNode(prev) ;
        }

        BoundaryReps->Store(EdgeKey(prev,new_node),new_bdry) ;
        int last_cur = new_node ;
        int last_prev = prev ;

        // now loop through and add the middle segments

        for (i=3 ; i<interface_seg.NumEntries() ; i+=3) {
            indx = interface_seg[i] ;
            prev = interface_seg[i+1] ;
            cur  = interface_seg[i+2] ;

            int new_node = Owner->NewNodeId() ;
            NewPoints->At(indx).node_id = new_node ;
            Owner->NodeTable->Store(new_node,NewPoints->At(indx).coord) ;
            CornerNodes->Store(new_node,1) ;

            Boundary->AddEdgeVertex(last_cur,last_prev,new_node) ;

            CArbBdrySeg *old_bdry =
                BoundaryReps->FetchValue(EdgeKey(prev,cur)) ;
            CArbBdrySeg *new_bdry = new CArbBdrySeg(*old_bdry) ;
            if (old_bdry->GetStartNode() == prev) {
                new_bdry->SetStartNode(last_cur) ;
                new_bdry->SetStopNode(new_node) ;
            } else {
                new_bdry->SetStartNode(new_node) ;
                new_bdry->SetStopNode(last_cur) ;
            }

            BoundaryReps->Store(EdgeKey(last_cur,new_node),new_bdry) ;
            last_prev = last_cur ;
            last_cur = new_node ;
        }

        // here we need to add the last segment to close off
        // the face.  First find the point that is on the
        // other side of the flaw from the current point

        for (i=0 ; i<NewPoints->NumEntries() ; ++i) {
            if (NewPoints->At(i).node_id == cur) break ;
        }
        int CW_vtx = NewPoints->At(i).node_id ;

        // move clockwise from this point

        int next_0 = (i==0) ? NewPoints->NumEntries()-1 : i-1 ;
        int next_bdry = NewPoints->At(next_0).node_id ;

        // we have two cases, the flaw may start at a boundary,
        // in which case the start node will be duplicated.
        // the second case is the flaw starting at a crack tip.
        //
        // Move clockwise one more point.  If the two clockwise
        // points have the same coordinates then we are at a
        // duplicated boundary and we want the second point.

        int next_1 = (next_0==0) ? NewPoints->NumEntries()-1 : next_0-1 ;

        int next_indx ;
        if (NewPoints->At(next_0).coord == NewPoints->At(next_1).coord) {
            next_indx = next_1 ;
            CW_vtx = NewPoints->At(next_0).node_id ;
        } else {
            next_indx = next_0 ;
        }

        int next = NewPoints->At(next_indx).node_id ;

        old_bdry = BoundaryReps->FetchValue(EdgeKey(cur,next_bdry)) ;
        new_bdry = new CArbBdrySeg(*old_bdry) ;
        if (old_bdry->GetStartNode() == cur) {
            new_bdry->SetStartNode(last_cur) ;
            new_bdry->SetStopNode(next) ;
        } else {
            new_bdry->SetStartNode(next) ;
            new_bdry->SetStopNode(last_cur) ;
        }

        int elem_id = Owner->NewElemId() ;
        int existing = Boundary->AddFaceEdge(last_cur,next,
                                  last_prev,CW_vtx,elem_id) ;
        BoundaryReps->Store(EdgeKey(last_cur,next),new_bdry) ;
        MatTable->Store(elem_id,MatTable->FetchValue(existing)) ;
    }

    delete LocalCrackPoints ;
    LocalCrackPoints = NewPoints ;
    NumPoints = LocalCrackPoints->NumEntries() ;
}


// this routine checks to see if a the crack is crossing a
// bi-material interface

bool CArbRmshRegion2D::CrossingInterface(
                           const LocCrackPt &this_pt,
                           const LocCrackPt &prev_pt,
                           const LocCrackPt &next_pt) const
{
    // This function does the geometry check necessary to
    // determine if the crack path is crossing a bi-material
    // interface.

    // The idea is to rank all the edges coming from this
    // point with respect to the incoming flaw edge vector.
    // the edges with the maximum and minimum rank bracket
    // this edge.  We then check to see if the outgoing
    // flaw edge is in the same sector.  If not we are crossing

    // first get a normalized version of the new edge

    CArbCoord2D norm = (prev_pt.coord-this_pt.coord).Normalize() ;

    // loop through adjacent edges in the boundary

    double min_score = .0f, max_score = .0f;
    bool first = true ;

    CArbTopoAdjVtxIterator iter(Boundary,this_pt.node_id) ;
    for (iter.First() ; iter.More() ; ++iter) {

        // get the tangent at this point

        CArbCoord2D adj =
             AdjacentTangent(this_pt.node_id,iter.AdjVtx()) ;

        double score = ScoreRelativeAngle(norm,adj) ;
        if (first || (score < min_score)) {
            min_score = score ;
        }
        if (first || (score > max_score)) {
            first = false ;
            max_score = score ;
        }
    }

    // check the outgoing

    CArbCoord2D adj = (next_pt.coord-this_pt.coord).Normalize() ;
    double score = ScoreRelativeAngle(norm,adj) ;

    return((score < max_score) && (score > min_score)) ;
}


void CArbRmshRegion2D::DuplicateSplineEndNode(
            LocCrackPt &existing,int far_node,int new_node)
{
    CArbBdrySeg *spline = BoundaryReps->FetchValue(
                      EdgeKey(existing.node_id,far_node)) ;
    BoundaryReps->Remove(EdgeKey(existing.node_id,far_node)) ;

    if (spline->GetStartNode() == existing.node_id) {
        spline->SetStartNode(new_node) ;
    } else {
        spline->SetStopNode(new_node) ;
    }

    BoundaryReps->Store(EdgeKey(new_node,far_node),spline) ;
    Boundary->SplitEdge(existing.node_id,far_node,new_node) ;

    // create a boundary segment for the duplicated node

    CArbCoord2D tang = spline->GetStartNode() == new_node ?
                 spline->Tangent(0.0) : spline->Tangent(1.0) ;

    CArbBdrySeg *seg = new CArbBdrySeg(existing.coord,tang,
                    existing.node_id,existing.char_elem_size,
                    new_node,existing.char_elem_size,
                    spline->GetParentId(),spline->IsBiMatInterface()) ;

    BoundaryReps->Store(EdgeKey(existing.node_id,new_node),seg) ;
}


int CArbRmshRegion2D::DuplicateInterfaceNode(int existing_indx)
{
    LocCrackPt &existing = LocalCrackPoints->At(existing_indx) ;
    LocCrackPt &next = LocalCrackPoints->At(existing_indx+1) ;

    // First we need to do geometry checks to determine on
    // which "side" of the existing node we need to do the
    // split.  Here we rely on the fact that all the flaw
    // segments are straight lines.

    int cw_vtx = FindCWvtx(existing.node_id,existing.coord,next.coord) ;

    // now loop ccw around this node to get the far vertex of the
    // edge we want to split

    int far_n = Boundary->GetCWNode(existing.node_id,cw_vtx) ;

    int new_node = Owner->NewNodeId() ;
    DuplicateSplineEndNode(existing,far_n,new_node) ;
    return(new_node) ;
}



void CArbRmshRegion2D::AddFlawToBoundary()
{
    int i,j,k ;
    CArbArray<LocCrackPt> *NewPoints = new CArbArray<LocCrackPt> ;

    // at this point we have found all intersections that the
    // the flaw makes with external boundaries and bimaterial
    // interfaces and have updated the topology of the remesh
    // boundary to include these points.  The list of points
    // defining the flaw has also been augmented to include
    // these points.

    // This routine steps through the flaw points and does
    // two things.  First, it updates the remesh boundary
    // topology to include portions of the flaw boundary
    // that are not already part of the boundary.  Second,
    // it create new spline representations for these segments.

    // for the logic below to work properly we want to start
    // with a point that is a boundary or bi-material interface
    // point.  If all points are inside it does not matter
    // where we start

    int start = 0 ;
    for (i=0 ; i<LocalCrackPoints->NumEntries() ; ++i) {
        if ((LocalCrackPoints->At(i).type == BOUNDARY) ||
            (LocalCrackPoints->At(i).type == INTERFACE)) {
            start = i ;
            break ;
        }
    }

    // If the first node is on the inside then we need to
    // create an id for it

    if (LocalCrackPoints->At(start).type == INSIDE) {
        int new_node = Owner->NewNodeId() ;
        LocalCrackPoints->At(start).node_id = new_node ;
        Owner->NodeTable->Store(new_node,
                                LocalCrackPoints->At(start).coord) ;
        CornerNodes->Store(new_node,1) ;
    }

    // now do it

    for (int cur=0 ; cur<LocalCrackPoints->NumEntries() ; ++cur) {

        i = (cur+start) % LocalCrackPoints->NumEntries() ;
        j = (i+1) % LocalCrackPoints->NumEntries() ;
        k = (i==0) ? LocalCrackPoints->NumEntries()-1 : i-1 ;

        LocCrackPt &this_pt = LocalCrackPoints->At(i) ;
        LocCrackPt &next_pt = LocalCrackPoints->At(j) ;
        LocCrackPt &prev_pt = LocalCrackPoints->At(k) ;

        if (this_pt.type != OUTSIDE) NewPoints->InsertAtEnd(this_pt) ;

        // if this point is outside or the next point is
        // outside then we ignore it

        if ((this_pt.type == OUTSIDE) || (next_pt.type == OUTSIDE)) continue ;

        // if the next point is inside the region then we
        // add a new edge and vertex unless this is the last point,
        // in that case add a face and edge

        if (next_pt.type == INSIDE) {

            if ((start == 0) && (j == 0)) {
                int CW_vtx0 = FindCWvtx(this_pt.node_id,
                                        this_pt.coord,
                                        next_pt.coord) ;

                int CW_vtx1 = FindCWvtx(next_pt.node_id,
                                        next_pt.coord,
                                        this_pt.coord) ;

                int elem_id = Owner->NewElemId() ;
                int existing = Boundary->AddFaceEdge(
                                      this_pt.node_id,
                                      next_pt.node_id,
                                      CW_vtx0,CW_vtx1,elem_id) ;
                MatTable->Store(elem_id,MatTable->FetchValue(existing)) ;
            } else {
                int new_node = Owner->NewNodeId() ;
                Owner->NodeTable->Store(new_node,next_pt.coord) ;
                CornerNodes->Store(new_node,1) ;
                LocalCrackPoints->At(j).node_id = new_node ;
                int CW_vtx = FindCWvtx(this_pt.node_id,
                                       this_pt.coord,
                                       next_pt.coord) ;

                if (CW_vtx != -1) {
                    Boundary->AddEdgeVertex(
                                    this_pt.node_id,
                                    CW_vtx,new_node) ;
                } else {
                    Boundary->AddEdgeTwoVertex(
                                    this_pt.node_id,
                                    next_pt.node_id,
                                    this_pt.elem) ;
                }
            }

        // else if this point and the next point are part of
        // the boundary or a bi-material interface we need
        // to determine if there is an edge connecting them.
        // if so we ignore them

        } else {
            if (((this_pt.type == BOUNDARY) ||
                 (this_pt.type == INTERFACE)) &&
                ((next_pt.type == BOUNDARY) ||
                 (next_pt.type == INTERFACE))) {

                // check for the special case of only two crack
                // points, both crack tips, and both in the
                // same bi-material interface.  For this we
                // need to add a new edge and face for the crack.

                if ((LocalCrackPoints->NumEntries() == 2) &&
                    HasCommonParent(this_pt.parent_splines,
                                    next_pt.parent_splines) &&
                    (this_pt.type == INTERFACE) &&
                    (next_pt.type == INTERFACE)) {

                    CArbBdrySeg *seg = BoundaryReps->FetchValue(
                         EdgeKey(this_pt.node_id,next_pt.node_id)) ;

                    // only add one side, the other is there already

                    if (seg->GetStartNode() == this_pt.node_id) {
                        int elem_id = Owner->NewElemId() ;
                        Boundary->TearEdgeAddFace(
                                          this_pt.node_id,
                                          next_pt.node_id,
                                          elem_id) ;
                    } else {
                        continue ;
                    }
                }

                CArbBdrySeg **ptr = BoundaryReps->Fetch(
                        EdgeKey(this_pt.node_id,next_pt.node_id)) ;

                if (ptr != 0) continue ;
            }

            // else we either moving from the inside to a boundary
            // or from one unconnected boundary to the other.  In
            // either case we add a new edge and element to the topology

            if (LocalCrackPoints->NumEntries() == 2) {

                int elem_id = Owner->NewElemId() ;
                Boundary->TearEdgeAddFace(this_pt.node_id,
                                          next_pt.node_id,
                                          elem_id) ;

            } else {

                int CW_vtx0 = FindCWvtx(this_pt.node_id,
                                        this_pt.coord,
                                        next_pt.coord) ;

                if (CW_vtx0 == -2) {
                    CW_vtx0 = Boundary->GetCWNode(this_pt.node_id,
                                                  prev_pt.node_id) ;
                }

                int CW_vtx1 = FindCWvtx(next_pt.node_id,
                                        next_pt.coord,
                                        this_pt.coord) ;

                if (CW_vtx1 == -2) {
                    CW_vtx1 = Boundary->GetCWNode(next_pt.node_id,
                                                  this_pt.node_id) ;
                }

                if (CW_vtx0 != -1) {
                    int celem = Boundary->FindCrossedFace(
                                          this_pt.node_id,
                                          next_pt.node_id,
                                          CW_vtx0,CW_vtx1) ;
                    if (celem == NO_ELEM) continue ;

                    int elem_id = Owner->NewElemId() ;
                    int existing = Boundary->AddFaceEdge(
                                          this_pt.node_id,
                                          next_pt.node_id,
                                          CW_vtx0,CW_vtx1,elem_id) ;
                    MatTable->Store(elem_id,MatTable->FetchValue(existing)) ;
                } else {
                    Boundary->AddEdgeVertex(next_pt.node_id,CW_vtx1,
                                            this_pt.node_id) ;
                }
            }
        }

        // Check for the special case of a surface crack along
        // bi-material interface.  For this we case we need to
        // duplicate the spline representation for the edge.
        // Otherwise create a (linear) spline representation for
        // this edge and add it to the edge rep store.

        if ((LocalCrackPoints->NumEntries() == 3) &&
            (next_pt.coord == prev_pt.coord)) {
            CArbBdrySeg **b_spline = BoundaryReps->Fetch(
                 EdgeKey(prev_pt.node_id,this_pt.node_id)) ;
            if (b_spline != 0) {
                CArbBdrySeg *n_spline = new CArbBdrySeg(**b_spline) ;
                if (this_pt.node_id == n_spline->GetStartNode()) {
                    n_spline->SetStopNode(next_pt.node_id) ;
                } else {
                    n_spline->SetStartNode(next_pt.node_id) ;
                }
                CArbBdrySeg **tmp = BoundaryReps->Fetch(
                      EdgeKey(next_pt.node_id,this_pt.node_id)) ;
                if (tmp != 0) delete *tmp ;
                BoundaryReps->Store(EdgeKey(next_pt.node_id,
                                    this_pt.node_id),
                                    n_spline) ;
                continue ;
            }
        }

        CArbCoord2D pts[2] ;
        pts[0] = next_pt.coord ;
        pts[1] = this_pt.coord ;
        CArbBdrySeg *b_spline =
                new CArbBdrySeg(2,pts,
                          next_pt.node_id,
                          next_pt.char_elem_size,
                          this_pt.node_id,
                          this_pt.char_elem_size,
                          0,false) ;
        CArbBdrySeg **tmp = BoundaryReps->Fetch(
             EdgeKey(next_pt.node_id,this_pt.node_id)) ;
        if (tmp != 0) delete *tmp ;
        BoundaryReps->Store(EdgeKey(next_pt.node_id,
                                    this_pt.node_id),
                                    b_spline) ;
    }

    // To determine which faces we can ignore we walk around
    // the flaw boundary and keep a list of the elements that
    // are to the "left" of these edges.

    IgnoreElems = new CArbHashTable<int,int> ;

    IgnoreElems->Store(NO_ELEM,1) ;
    for (i=0 ; i<LocalCrackPoints->NumEntries() ; ++i) {
        j = (i+1) % LocalCrackPoints->NumEntries() ;

        if ((LocalCrackPoints->At(i).type != OUTSIDE) &&
            (LocalCrackPoints->At(j).type != OUTSIDE)) {

            if ((LocalCrackPoints->At(i).type == BOUNDARY) &&
                (LocalCrackPoints->At(j).type == BOUNDARY) &&
                HasCommonParent(
                    LocalCrackPoints->At(i).parent_splines,
                    LocalCrackPoints->At(j).parent_splines)) continue ;

            IgnoreElems->Store(
                Boundary->GetCCWElem(
                    LocalCrackPoints->At(i).node_id,
                    LocalCrackPoints->At(j).node_id),1) ;
        }
    }

    delete LocalCrackPoints ;
    LocalCrackPoints = NewPoints ;
    NumPoints = LocalCrackPoints->NumEntries() ;
}


void CArbRmshRegion2D::ResizeCrackTemplates()
{
    int i, j, k, l, m ;

    // This routine resizes the crack-tip templates so that they
    // do not intersect each other or any of the boundaries

    // Here we check for three resize cases:
    //
    // 1. If the template radius is greater than the length
    //    of one of the adjacent segments then we scale the
    //    template size so that it is the size of the adjacent
    //    segment.
    //
    // 2. If by adding the template to an ajacent crack segment
    //    the remaining element edge is small relative to the
    //    template size we either scale up the template size to
    //    be the size of the segment or we scale down template
    //    size so that we have two more reasonably sized segments
    //
    // 3. If the previous or next points are crack-tips then
    //    make sure that the templates do not overlap

    for (i=0 ; i<NumPoints ; ++i) {

        if (LocalCrackPoints->At(i).crack_tip) {

            j = (i+1) % LocalCrackPoints->NumEntries() ;

            LocCrackPt &this_pt = LocalCrackPoints->At(i) ;
            LocCrackPt &next_pt = LocalCrackPoints->At(j) ;

            int tip_id = this_pt.node_id ;

            // get the crack tip data

            CArbMshCrackRegion2D::CrackTipData &tip_data =
                TipData->At(this_pt.tip_indx) ;

            // loop through all the segments adjacent to the crack

            CArbTopoAdjVtxIterator iter(Boundary,tip_id) ;
            for (iter.First() ; iter.More() ; ++iter) {

                // get the length of the segment

                CArbBdrySeg *spline =
                    BoundaryReps->FetchValue(EdgeKey(tip_id,iter.AdjVtx())) ;

                double len = spline->ApproxLength() ;

                // check for the length being too long

                if (tip_data.template_radius >= len) {
                    tip_data.template_radius = len ;

                // check to see if we need to expand the template

                } else if (tip_data.template_radius/len >=
                           TemplateScaleUpFactor) {
                    tip_data.template_radius = len ;

                // check to see if we need to scale down the template

                } else if (tip_data.template_radius/len >
                           TemplateScaleDownFactor) {
                    tip_data.template_radius =
                        TemplateScaleDownFactor*len ;
                }
            }

            // check for next crack tip and overlap

            if (next_pt.crack_tip) {

                CArbMshCrackRegion2D::CrackTipData &next_data =
                    TipData->At(next_pt.tip_indx) ;

                double len = (next_pt.coord - this_pt.coord).Magnitude() ;
                double left = len - tip_data.template_radius -
                              next_data.template_radius ;

                // if the two templates overlap or if the remaining
                // element would be small relative to the segment
                // length then we scale the templates so that the
                // ratio between them stays the same but they fill
                // the full segment.

                if (left < (1.0-TemplateScaleUpFactor)*len) {
                    double ratio = tip_data.template_radius /
                                   next_data.template_radius ;

                    tip_data.template_radius =
                        len * (1.0 - 1.0/(ratio + 1.0)) ;

                    next_data.template_radius =
                        len / (ratio + 1.0) ;
                }
            }
        }
    }

    // Now we loop through a second time and check to see if
    // the template will overlap or be too close to an existing
    // boundary edge.

    for (i=0 ; i<NumPoints ; ++i) {

        if (LocalCrackPoints->At(i).crack_tip) {

            j = (i+1) % LocalCrackPoints->NumEntries() ;
            k = (i==0) ? LocalCrackPoints->NumEntries()-1 : i-1 ;
            l = (j+1) % LocalCrackPoints->NumEntries() ;
            m = (k==0) ? LocalCrackPoints->NumEntries()-1 : k-1 ;

            LocCrackPt &this_pt = LocalCrackPoints->At(i) ;
            LocCrackPt &next_pt = LocalCrackPoints->At(j) ;
            LocCrackPt &prev_pt = LocalCrackPoints->At(k) ;
            LocCrackPt &next_next_pt = LocalCrackPoints->At(l) ;
            LocCrackPt &prev_prev_pt = LocalCrackPoints->At(m) ;

            CArbMshCrackRegion2D::CrackTipData &tip_data =
                TipData->At(this_pt.tip_indx) ;

            // loop through all the edges

            double min_dist = .0f;
            bool first = true ;

            CArbTopoEdgeIterator edge_iter(*Boundary) ;
            for (edge_iter.First() ; edge_iter.More() ; ++edge_iter) {

                if ((this_pt.node_id == edge_iter[0]) ||
                    (this_pt.node_id == edge_iter[1])) continue ;

                CArbBdrySeg **ptr = BoundaryReps->Fetch(
                          EdgeKey(edge_iter[0],edge_iter[1])) ;

                double dist,bedge ;

                // two cases, first for a spline edge, second
                // for a straight line case

                if ((ptr != 0) && ((*ptr)->GetSpline() != 0)) {
                    double u ;
                    CArbBdrySeg *p_seg = *ptr ;
                    p_seg->GetSpline()->ClosestPoint(this_pt.coord,&u) ;
                    dist = (p_seg->GetSpline()->Evaluate(u) -
                            this_pt.coord).Magnitude() ;
                    bedge = 0.5 ;
                } else {
                    dist = PerpendicularDist(this_pt.coord,
                            Owner->NodeTable->FetchValue(edge_iter[0]),
                            Owner->NodeTable->FetchValue(edge_iter[1]),
                            &bedge) ;
                }

                if (((edge_iter[0] == prev_pt.node_id) &&
                     (edge_iter[1] == prev_prev_pt.node_id)) ||
                    ((edge_iter[0] == next_pt.node_id) &&
                     (edge_iter[1] == next_next_pt.node_id)) ||
                    ((edge_iter[1] == prev_pt.node_id) &&
                     (edge_iter[0] == prev_prev_pt.node_id)) ||
                    ((edge_iter[1] == next_pt.node_id) &&
                     (edge_iter[0] == next_next_pt.node_id))) continue ;

                if (bedge != 0.0) {
                    if (first || (dist < min_dist)) {
                        first = false ;
                        min_dist = dist ;
                    }
                }
            }

            // if necessary we resize the template

            if (!first && (tip_data.template_radius/min_dist >
                          Owner->TemplateBoundaryFactor)) {
                    tip_data.template_radius =
                        Owner->TemplateBoundaryFactor*min_dist ;
            }
        }
    }
}


#define TEMP_FACTOR .75
#define SMALL_FACTOR 0.0001

void CArbRmshRegion2D::UpdateBdrySegCharLengths(bool tip_only)
{
    int i ;

    // This routines looks through all the segments on the
    // boundary and determines if the characteristic element
    // size should be reduced because there is a flaw point
    // close by.

    CArbTopoEdgeIterator edge_iter(*Boundary) ;
    for (edge_iter.First() ; edge_iter.More() ; ++edge_iter) {

        // if the segment does not have a BdrySeg
        // that means that it is an existing element
        // edge and we cannot resize it, so we ignore it.

        CArbBdrySeg **ptr = BoundaryReps->Fetch(
                  EdgeKey(edge_iter[0],edge_iter[1])) ;

        if ((ptr == 0) || ((*ptr)->GetSpline() == 0)) continue ;

        // find the flaw point closest to this edge

        // min_dist and min_u are arrays of size 3.  This is because
        // if we find a segment that needs to be refined we want to
        // check for refinement on it's start and end nodes in addition
        // to just the mid point.  The 0 index data is the mid
        // information, 1 index data for u = 0 and 2 index data
        // for u = 1

        bool first = true ;
        double min_dist[3] ;
        double min_u[3] ;
        int min_pt = 0;
        CArbBdrySeg *min_seg[3] ;
        for (i=0 ; i<NumPoints ; ++i) {

            LocCrackPt &tmp_pt = LocalCrackPoints->At(i) ;
            if (tip_only && !tmp_pt.crack_tip) continue ;
            double char_size ;
            if (tmp_pt.crack_tip) {
                CArbMshCrackRegion2D::CrackTipData &tip_data =
                    TipData->At(tmp_pt.tip_indx) ;
                char_size = tip_data.template_radius ;
            } else {
                char_size = tmp_pt.char_elem_size ;
            }

            // determine the shortest distance from the
            // crack tip to the segment

            double dist, u ;
            CArbBdrySeg *p_seg = *ptr ;
            p_seg->GetSpline()->ClosestPoint(tmp_pt.coord,&u) ;
            dist = (p_seg->GetSpline()->Evaluate(u) -
                        tmp_pt.coord).Magnitude() ;

            if (first ||
                ((dist < min_dist[0]) && (dist > 0.01*char_size))) {
                first = false ;
                min_pt = i ;
                min_u[0] = u ;
                min_dist[0] = dist ;
                min_seg[0] = p_seg ;
            }
        }
        if (first) continue ;

        // now check to see if we should refine the edge

        LocCrackPt &this_pt = LocalCrackPoints->At(min_pt) ;

        // if this point is one of the boundary edge points
        // then we go on to the next one

        if ((this_pt.node_id == edge_iter[0]) ||
            (this_pt.node_id == edge_iter[1])) continue ;

        // Also, if the distance is very small we continue
        // because the refinement will be done elsewere

        if (min_dist[0] < SMALL_FACTOR*this_pt.char_elem_size) continue ;

        // fill in the data for the segment end points

        min_u[1] = 0.0 ;
        min_u[2] = 1.0 ;
        CArbBdrySeg *p_seg = *ptr ;
        min_dist[1] = (p_seg->GetSpline()->Evaluate(0.0) -
                       this_pt.coord).Magnitude() ;
        min_dist[2] = (p_seg->GetSpline()->Evaluate(1.0) -
                       this_pt.coord).Magnitude() ;
        min_seg[1] = min_seg[0] ;
        min_seg[2] = min_seg[0] ;

        // get the characteristic size for this point

        double char_size ;
        if (this_pt.crack_tip) {
            CArbMshCrackRegion2D::CrackTipData &tip_data =
                TipData->At(this_pt.tip_indx) ;
            char_size = tip_data.template_radius ;
        } else {
            char_size = this_pt.char_elem_size ;
        }

        // If necessary split the spline by adding a point
        // with a smaller characteristic element size

        for (i=0 ; i<3 ; ++i) {
            if ((this_pt.crack_tip && (min_dist[i]/char_size > 1.001)) ||
                (min_dist[i] < char_size)) {

                double delta = this_pt.crack_tip ?
                               min_dist[i] - char_size : min_dist[i] ;
                double new_clen = .0f;
                bool update = false ;
                double csize = min_seg[i]->GetStartCharLen() +
                                min_u[i]*(min_seg[i]->GetStopCharLen() -
                                min_seg[i]->GetStartCharLen()) ;

                if (delta < char_size) {
                    if (csize/delta > Owner->BoundarySpacingFactor) {
                        update = true ;
                        new_clen = delta * Owner->BoundarySpacingFactor ;
                    }
                } else {

                    // this check is based on the assumption that
                    // as we move away from the template the
                    // characteristic element size should grow at
                    // a ratio no greater than 1.5.  That is,
                    // l_i+1 / l_i <= 1.5 and the total distance
                    // to element i is:
                    //    L_i = l_0*(1+1.5+1.5^2+...+1.5^i)
                    //
                    // We don't want to sum this series for each
                    // check so I have done an eyeball fit of this
                    // equation.  The approximate function is
                    //    L_i = l_0*(1+(0.6*i)^2.85)
                    // or
                    //    i = 1.667*((L_i-l_0)/l_0)^(1/2.85)
                    //
                    // with the corresponding maximum local size:
                    //    l_i = l_0 * 1.5^i
                    //
                    // l_0 is the effective template radius we
                    // set this to be approximately the size
                    // of one of the element edges on the outside
                    // of the template.  For 8 elements around the
                    // the tip this is approximately 0.75 times
                    // the radius.

                    double eff_radius = TEMP_FACTOR * char_size ;

                    double edist = 1.667 * pow((min_dist[i] -
                        eff_radius) / eff_radius,0.3509) ;
                    edist = floor(edist) ;
                    double loc_clen = pow(1.5,edist) *
                                     eff_radius ;
                    if (csize > loc_clen) {
                        update = true ;
                        new_clen = loc_clen ;
                    }
                }

                if (update) {
                    CArbBdrySeg *mate_seg = min_seg[i]->GetMate() ;
                    CArbCoord2D new_crd = min_seg[i]->Evaluate(min_u[i]) ;
                    if ((min_u[i] < 0.05) ||
                       ((min_seg[i]->Evaluate(0.0)-new_crd).Magnitude() < new_clen)) {
                        min_seg[i]->SetStartCharLen(new_clen) ;
                        if (mate_seg != 0)
                            mate_seg->SetStartCharLen(new_clen) ;
                    } else if ((min_u[i] > 0.95) ||
                       ((min_seg[i]->Evaluate(1.0)-new_crd).Magnitude() < new_clen)) {
                        min_seg[i]->SetStopCharLen(new_clen) ;
                        if (mate_seg != 0)
                            mate_seg->SetStopCharLen(new_clen) ;
                    } else {
                        int new_node = Owner->NewNodeId() ;
                        Owner->NodeTable->Store(new_node,
                                    min_seg[i]->Evaluate(min_u[i])) ;
                        CornerNodes->Store(new_node,1) ;
                        int start = min_seg[i]->GetStartNode() ;
                        int stop  = min_seg[i]->GetStopNode() ;
                        SplitSpline(min_seg[i],new_node,min_u[i],
                                    new_clen,Boundary) ;
                        CArbBdrySeg **ptr = BoundaryReps->Fetch(
                                     EdgeKey(start,new_node)) ;
                        min_seg[1] = *ptr ;
                        ptr = BoundaryReps->Fetch(
                                     EdgeKey(new_node,stop)) ;
                        min_seg[2] = *ptr ;
                        if (mate_seg != 0) {
                            int new_node = Owner->NewNodeId() ;
                            Owner->NodeTable->Store(new_node,
                                    mate_seg->Evaluate(min_u[i])) ;
                            CornerNodes->Store(new_node,1) ;
                            SplitSpline(mate_seg,new_node,min_u[i],
                                        new_clen,Boundary) ;
                        }
                    }
                }
            }
        }
    }
}


#define PROGRESSION_RATIO 0.6667

void CArbRmshRegion2D::UpdateCrackSegCharLengths()
{
    int i ;

    // This routines looks through all the segments defined for
    // the crack and sets the characteristic element sizes for
    // any points that are not crack tips or the size was not
    // specified on input.

    // The basic approach is that we start moving around the crack
    // points until we find a point that is either a crack tip or
    // has the char size specified.  We store this point and keep
    // moving and storing until we find a second point with a know
    // characteristic size.  Then, if there are any points between
    // these we start moving from the two end points towards the
    // the middle assigning characteristic element sizes along
    // the way.  If there are an odd number of points we average
    // in the middle.

    // find a starting point

    bool found = false ;
    int first_indx = 0, cur_indx = 0;
    for (i=0 ; i<NumPoints ; ++i) {
        LocCrackPt &this_pt = LocalCrackPoints->At(i) ;
        if (this_pt.crack_tip || this_pt.respect_char_size) {
            found = true ;
            first_indx = i ;
            cur_indx = i ;
            break ;
        }
    }

    // If we've gone all the way around and not found a point with
    // a specified characteristic length, we stick with the default
    // lengths that are a function of the insitu element sizes.

    if (!found) return ;

    // loop through all the points, and for any crack tips set
    // the characteristic element length field.

    for (i=0 ; i<NumPoints ; ++i) {
        LocCrackPt &this_pt = LocalCrackPoints->At(i) ;
        if (this_pt.crack_tip) {
            CArbMshCrackRegion2D::CrackTipData &tip_data =
                TipData->At(this_pt.tip_indx) ;
            this_pt.char_elem_size = tip_data.template_radius ;
            this_pt.respect_char_size = true ;
        }
    }

    // now go around the crack points again

    CArbArray<int> points ;
    bool first = true ;

    while (first || (first_indx != cur_indx)) {
        first = false ;
        points.Clear() ;

        do {
            points.InsertAtEnd(cur_indx) ;
            cur_indx = (cur_indx+1) % LocalCrackPoints->NumEntries() ;
        } while (!LocalCrackPoints->At(cur_indx).crack_tip &&
                 !LocalCrackPoints->At(cur_indx).respect_char_size) ;
        points.InsertAtEnd(cur_indx) ;

        // if there are more than two point in the array then
        // we need to update the middle points

        if (points.NumEntries() > 2) {
            for (i=1 ; i<(points.NumEntries()/2)-1 ; ++i) {
                LocCrackPt &first_pt = LocalCrackPoints->At(points[i-1]) ;
                LocCrackPt &last_pt = LocalCrackPoints->At(points[-i]) ;
                LocCrackPt &pt0 = LocalCrackPoints->At(points[i]) ;
                LocCrackPt &pt1 = LocalCrackPoints->At(points[-(i+1)]) ;

                double pt0_len, pt1_len ;
                double fnum,int_part,frac_part,aratio ;
                int num ;

                // determine the characteristic sizes at the first point

                CArbBdrySeg *b_seg = BoundaryReps->FetchValue(
                               EdgeKey(first_pt.node_id,pt0.node_id)) ;
                double len = b_seg->ApproxLength() ;

                if (len == 0.0) {
                    pt0_len = first_pt.char_elem_size ;
                } else {
                    fnum = ApproxNumSegs(first_pt.char_elem_size,
                                         PROGRESSION_RATIO,len) ;
                    frac_part = modf(fnum,&int_part) ;
                    num = (frac_part > 0.1) ? int(fnum+1.0) : int(fnum) ;
                    aratio = ApproxRatio(num,PROGRESSION_RATIO) ;
                    pt0_len = aratio * pt0.char_elem_size ;
                }

                // and the far point

                b_seg = BoundaryReps->FetchValue(
                           EdgeKey(pt1.node_id,last_pt.node_id)) ;
                len = b_seg->ApproxLength() ;

                if (len == 0.0) {
                    pt1_len = last_pt.char_elem_size ;
                } else {
                    fnum = ApproxNumSegs(last_pt.char_elem_size,
                                         PROGRESSION_RATIO,len) ;
                    frac_part = modf(fnum,&int_part) ;
                    num = (frac_part > 0.1) ? int(fnum+1.0) : int(fnum) ;
                    aratio = ApproxRatio(num,PROGRESSION_RATIO) ;
                    pt1_len = aratio * pt1.char_elem_size ;
                }

                // assign these values

                if (pt0.node_id == pt1.node_id) {
                    double avg = 0.5*(pt0_len+pt1_len) ;
                    if (avg <  pt0.char_elem_size)
                        pt0.char_elem_size = avg ;
                } else {
                    if (pt0_len < pt0.char_elem_size)
                        pt0.char_elem_size = pt0_len ;
                    if (pt1_len < pt1.char_elem_size)
                        pt1.char_elem_size = pt1_len ;
                }
            }
        }
    }
}


void CArbRmshRegion2D::AddCrackTemplates()
{
    int i, j, k ;

    // This routine adds the crack-tip element templates to the
    // remesh boundary

    for (i=0 ; i<NumPoints ; ++i) {

        if (LocalCrackPoints->At(i).crack_tip) {

            int tip_id = LocalCrackPoints->At(i).node_id ;
            Owner->CrackTable->Store(LocalCrackPoints->At(i).tip_id,
                                     tip_id) ;

            // get the crack tip data

            CArbMshCrackRegion2D::CrackTipData &tip_data =
                TipData->At(LocalCrackPoints->At(i).tip_indx) ;
            CArbCoord2D &tip = LocalCrackPoints->At(i).coord ;

            // the sector data array stores information for each
            // sector that needs to be meshed.  There are 3
            // integers for each sector, the prev remote id,
            // the next remote id, and the material id.

            CArbArray<int> sdata ;

            // loop through the edges adjacent to the crack
            // tip and determine which sectors need to be meshed

            CArbTopoAdjVtxCyclicIterator iter(Boundary,tip_id) ;

            int prev_id = iter.AdjVtx() ;
            int prev_elem = iter.CcwElem() ;

            ++iter ;
            int stop_vtx = iter.AdjVtx() ;
            int stop_elem = iter.CcwElem() ;
            bool first = true ;

            while (first || (iter.AdjVtx() != stop_vtx ||
                             iter.CcwElem() != stop_elem)) {

                first = false ;
                int next_id = iter.AdjVtx() ;
                int next_elem = iter.CcwElem() ;

                // check to see if this sector should be meshed

                if (IgnoreElems->Fetch(prev_elem) == 0) {
                    sdata.InsertAtEnd(prev_id) ;
                    sdata.InsertAtEnd(next_id) ;
                    sdata.InsertAtEnd(MatTable->FetchValue(prev_elem)) ;
                    sdata.InsertAtEnd(prev_elem) ;
                    sdata.InsertAtEnd(next_elem) ;
                }

                prev_elem = next_elem ;
                prev_id = next_id ;
                ++iter ;
            }

            // this is a hash table that keeps track of the
            // template nodes on the boundary of the segment
            // so that the nodes on the edges of the segments
            // do not get duplicated.  The hash key is the id
            // of the remote node.  The value is an array of
            // the nodes starting from the edge of the template
            // and going towards (but excluding) the template

            CArbHashTable<EdgeFaceKey,CArbArray<int> > bdry_nodes ;

            // now loop through the segments that need a mesh
            // and add it

            for (j=0 ; j<sdata.NumEntries() ; j+=5) {

                // find the crack normal, the sector angle
                // and the number of crack tip elements

                CArbCoord2D prev_tan =
                    AdjacentTangent(tip_id,sdata[j]) ;
                CArbCoord2D next_tan =
                    AdjacentTangent(tip_id,sdata[j+1]) ;

                CArbCoord2D norm = BisectNormVect(
                                     prev_tan,next_tan) ;

                double tip_angle = Angle2PiVect(prev_tan,next_tan) ;
                if (tip_angle < 0.001) tip_angle = TWO_PI ;

                double int_part,frac_part ;
                frac_part = modf(tip_angle * tip_data.num_tip_elems /
                                 TWO_PI,&int_part) ;
                int num_tip_elems = (frac_part > 0.2) ?
                                    int(int_part+1.0) : int(int_part) ;

                // generate the crack template nodes and boundary

                CArbArray<ArbMshNode> nodes ;
                CArbArray<int> bdry ;

                GenerateCrackTipNodes(tip,norm,num_tip_elems,
                                      tip_data,TWO_PI-tip_angle,nodes) ;
                GenerateCrackTipBdry(num_tip_elems,tip_data,bdry) ;

                // overwrite the crack tip node id

                nodes[0].id = tip_id ;

                // find the length of the last element in the template
                // and compute the characteristic element size for the
                // the template mouth node

                CArbCoord2D delt(nodes[bdry[0]].coord[0] -
                                 nodes[bdry[1]].coord[0],
                                 nodes[bdry[0]].coord[1] -
                                 nodes[bdry[1]].coord[1]) ;

                double new_char_len = tip_data.template_radius +
                    TemplateCharSizeFactor * (delt.Magnitude() -
                        tip_data.template_radius) ;

                // check to see if we have processed these sector
                // edges yet, first the previous edge

                CArbArray<int> *bnodes = bdry_nodes.Fetch(
                         EdgeFaceKey(sdata[j],sdata[j+3])) ;
                if (bnodes == 0) {

                    // check to see if the template takes up
                    // the full segment

                    double len = (Owner->NodeTable->FetchValue(sdata[j]) -
                                  tip).Magnitude() ;
                    double ratio = tip_data.template_radius / len ;

                    CArbArray<int> tmp ;
                    if (fabs(1.0-ratio) < 0.01) {
                        nodes[bdry[0]].id = sdata[j] ;
                        int last = sdata[j] ;
                        tmp.InsertAtEnd(last) ;
                        for (k=1 ; k<bdry.NumEntries()/2 ; ++k) {
                            Boundary->SplitEdgeElem(tip_id,last,
                                       nodes[bdry[k]].id,
                                       sdata[j+3],true) ;
                            last = nodes[bdry[k]].id ;
                            tmp.InsertAtEnd(last) ;
                        }
                    } else {
                        UpdateBRepsForTemplates(tip_id,sdata[j],
                                            nodes[bdry[0]].id,
                                            tip_data.template_radius,
                                            new_char_len) ;
                        int last = sdata[j] ;
                        for (k=0 ; k<bdry.NumEntries()/2 ; ++k) {
                            Boundary->SplitEdgeElem(tip_id,last,
                                       nodes[bdry[k]].id,
                                       sdata[j+3],true) ;
                            last = nodes[bdry[k]].id ;
                            tmp.InsertAtEnd(last) ;
                        }
                    }
                    bdry_nodes.Store(
                        EdgeFaceKey(sdata[j],sdata[j+3]),tmp) ;
                } else {
                    for (k=0 ; k<bdry.NumEntries()/2 ; ++k) {
                        nodes[bdry[k]].id = bnodes->At(k) ;
                    }
                }

                // then the next edge

                bnodes = bdry_nodes.Fetch(
                         EdgeFaceKey(sdata[j+1],sdata[j+4])) ;
                if (bnodes == 0) {

                    // check to see if the template takes up
                    // the full segment

                    double len = (Owner->NodeTable->FetchValue(sdata[j+1]) -
                                  tip).Magnitude() ;
                    double ratio = tip_data.template_radius / len ;
                    CArbArray<int> tmp ;

                    if (fabs(1.0-ratio) < 0.01) {
                        nodes[bdry[-1]].id = sdata[j+1] ;
                        int last = sdata[j+1] ;
                        tmp.InsertAtEnd(last) ;
                        for (k=1 ; k<bdry.NumEntries()/2 ; ++k) {
                            Boundary->SplitEdgeElem(tip_id,last,
                                       nodes[bdry[-(k+1)]].id,
                                       sdata[j+3],false) ;
                            last = nodes[bdry[-(k+1)]].id ;
                            tmp.InsertAtEnd(last) ;
                        }
                    } else {
                        UpdateBRepsForTemplates(tip_id,sdata[j+1],
                                            nodes[bdry[-1]].id,
                                            tip_data.template_radius,
                                            new_char_len) ;
                        int last = sdata[j+1] ;
                        for (k=0 ; k<bdry.NumEntries()/2 ; ++k) {
                            Boundary->SplitEdgeElem(tip_id,last,
                                       nodes[bdry[-(k+1)]].id,
                                       sdata[j+3],false) ;
                            last = nodes[bdry[-(k+1)]].id ;
                            tmp.InsertAtEnd(last) ;
                        }
                    }
                    bdry_nodes.Store(
                        EdgeFaceKey(sdata[j+1],sdata[j+4]),tmp) ;
                } else {
                    for (k=0 ; k<bdry.NumEntries()/2 ; ++k) {
                        nodes[bdry[-(k+1)]].id = bnodes->At(k) ;
                    }
                }

                // generate the elements

                CArbArray<ArbMshElement2D> elems ;
                GenerateCrackTipElems(num_tip_elems,
                                      tip_data,nodes,
                                      sdata[j+2],elems) ;

                // insert the new nodes into the Node table

                for (k=0 ; k<nodes.NumEntries() ; ++k) {
                    Owner->NodeTable->Store(nodes[k].id,
                       CArbCoord2D(nodes[k].coord[0],nodes[k].coord[1])) ;
                }

                // insert the new elements into the element table
                // and into the remesh boundary

                for (k=0 ; k<elems.NumEntries() ; ++k) {
                    Owner->AddElem(elems[k]) ;
                    IgnoreElems->Store(elems[k].elem_id,1) ;
                    int lnum = (elems[k].num_nodes == 3) ||
                               (elems[k].num_nodes == 6) ? 3 : 4 ;
                    for (int i=0 ; i<lnum ; ++i)
                        CornerNodes->Store(elems[k].nodes[i],1) ;
                }
                if (Owner->Order == LINEAR) {
                    for (int k=0 ; k<elems.NumEntries() ; ++k) {
                        Boundary->InsertCollapsedElement(
                                       elems[k].elem_id,
                                       elems[k].num_nodes,
                                       elems[k].nodes,sdata[j+3]) ;
                    }
                } else {  // if qaudratic reorder the nodes
                    int nodes[8], jj ;
                    for (int k=0 ; k<elems.NumEntries() ; ++k) {
                        int cur = 0 ;
                        for (jj=0 ; jj<elems[k].num_nodes/2 ; ++jj) {
                            nodes[jj*2] = elems[k].nodes[jj] ;
                            ++cur ;
                        }
                        for (jj=0 ; jj<elems[k].num_nodes/2 ; ++jj) {
                            nodes[jj*2+1] = elems[k].nodes[cur] ;
                            ++cur ;
                        }
                        Boundary->InsertCollapsedElement(
                                       elems[k].elem_id,
                                       elems[k].num_nodes,
                                       nodes,sdata[j+3]) ;
                    }
                }
            }
        }
    }

    // adding the template elements may have buggered up the
    // the element list in the topo object so we rebuild it

    Boundary->RebuildElemList() ;
}

void CArbRmshRegion2D::UpdateBRepsForTemplates(
                 int tip_id,int rmt_id,int new_id,
                 double radius,double new_char_len)
{

    // we need to update the boundary representations
    // of the crack faces.  Split the existing reps
    // so that they go to the edge of the template and
    // not all the way to the tip.  Also, remove and
    // delete the old reps

    EdgeKey old_key = EdgeKey(tip_id,rmt_id) ;

    CArbBdrySeg *spline = BoundaryReps->FetchValue(old_key) ;

    if (!Boundary->HasDblEdge(tip_id,rmt_id) &&
        !Boundary->HasDblEdge(rmt_id,tip_id))
        BoundaryReps->Remove(old_key) ;

    CArbBdrySeg *p_sub ;

    double split_par = radius / spline->ApproxLength() ;

    if (tip_id == spline->GetStartNode()) {

        p_sub = spline->SubSegment(
                                split_par,1.0,
                                new_id,new_char_len,
                                spline->GetStopNode(),
                                spline->GetStopCharLen(),
                                spline->IsBiMatInterface()) ;

    } else {

        p_sub = spline->SubSegment(
                                0.0,1.0-split_par,
                                spline->GetStartNode(),
                                spline->GetStartCharLen(),
                                new_id,new_char_len,
                                spline->IsBiMatInterface()) ;

    }

    BoundaryReps->Store(EdgeKey(new_id,rmt_id),p_sub) ;
//    delete spline ;
}

void CArbRmshRegion2D::AddNodesToBoundary()
{
    // This routine adds nodes to all edges in the remesh
    // boundary excluding those edges that fall inside the
    // the flaw.

    // build a list of all the edges that we may need to
    // subdivide.

    CArbHashTable<EdgeKey,CArbBdrySeg *> sub_div_edges ;

    // loop through all the nodes

    CArbTopoVtxIterator viter(*Boundary) ;
    for (viter.First() ; viter.More() ; ++viter) {

        // loop through all edges adjacent to this node

        CArbTopoAdjVtxIterator aiter(Boundary,*viter) ;
        for (aiter.First() ; aiter.More() ; ++aiter) {

            // check for no face or a noremesh face

            if ((aiter.CcwElem() == NO_ELEM) ||
                (IgnoreElems->Fetch(aiter.CcwElem()) != 0)) continue ;

            // check to see if there is a spline rep for this edge

            CArbBdrySeg **ptr = BoundaryReps->Fetch(
                    EdgeKey(*viter,aiter.AdjVtx())) ;

            if (ptr != 0) {
                sub_div_edges.Store(
                     EdgeKey(*viter,aiter.AdjVtx()),*ptr) ;
            }
        }
    }

    // For crack face segments we want to insure that the
    // characteristic lengths are the same on both sides of
    // the crack so that the node spacing will be the same.
    // To do this we search through the subdivide edge list
    // and do comparisons on the nodes.  The comparison seach
    // is dump (n^2), but there should not be many entries
    // in the hash table.

    CArbHashTableIterator<EdgeKey,CArbBdrySeg *>
            iteri(&sub_div_edges) ;
    for (iteri.First() ; iteri.More() ; ++iteri) {

        CArbCoord2D *cdi0 = Owner->NodeTable->Fetch(iteri.Key().id0) ;
        CArbCoord2D *cdi1 = Owner->NodeTable->Fetch(iteri.Key().id1) ;
        double tol = 0.001 * (*cdi0 - *cdi1).Magnitude() ;
        CArbBdrySeg *segi = *(iteri.Entry()) ;
        bool iorient = iteri.Key().id0 == segi->GetStartNode() ;

        CArbHashTableIterator<EdgeKey,CArbBdrySeg *> iterj(iteri) ;
        ++iterj ;
        while (iterj.More()) {

            CArbCoord2D *cdj0 = Owner->NodeTable->Fetch(iterj.Key().id0) ;
            CArbCoord2D *cdj1 = Owner->NodeTable->Fetch(iterj.Key().id1) ;

            if (((*cdi0 - *cdj0).Magnitude() < tol) &&
                ((*cdi1 - *cdj1).Magnitude() < tol)) {

                CArbBdrySeg *segj = *(iterj.Entry()) ;
                bool jorient = iterj.Key().id0 == segj->GetStartNode() ;
                double clen0i, clen1i, clen0j, clen1j ;

                if (iorient) {
                    clen0i = segi->GetStartCharLen() ;
                    clen1i = segi->GetStopCharLen() ;
                } else {
                    clen0i = segi->GetStopCharLen() ;
                    clen1i = segi->GetStartCharLen() ;
                }
                if (jorient) {
                    clen0j = segj->GetStartCharLen() ;
                    clen1j = segj->GetStopCharLen() ;
                } else {
                    clen0j = segj->GetStopCharLen() ;
                    clen1j = segj->GetStartCharLen() ;
                }

                double avg0 = 0.5 * (clen0i + clen0j) ;
                double avg1 = 0.5 * (clen1i + clen1j) ;

                if (iorient) {
                    segi->SetStartCharLen(avg0) ;
                    segi->SetStopCharLen(avg1) ;
                } else {
                    segi->SetStopCharLen(avg0) ;
                    segi->SetStartCharLen(avg1) ;
                }
                if (jorient) {
                    segj->SetStartCharLen(avg0) ;
                    segj->SetStopCharLen(avg1) ;
                } else {
                    segj->SetStopCharLen(avg0) ;
                    segj->SetStartCharLen(avg1) ;
                }

            } else if (((*cdi0 - *cdj1).Magnitude() < tol) &&
                       ((*cdi1 - *cdj0).Magnitude() < tol)) {

                CArbBdrySeg *segj = *(iterj.Entry()) ;
                bool jorient = iterj.Key().id0 == segj->GetStartNode() ;
                double clen0i, clen1i, clen0j, clen1j ;

                if (iorient) {
                    clen0i = segi->GetStartCharLen() ;
                    clen1i = segi->GetStopCharLen() ;
                } else {
                    clen0i = segi->GetStopCharLen() ;
                    clen1i = segi->GetStartCharLen() ;
                }
                if (jorient) {
                    clen0j = segj->GetStartCharLen() ;
                    clen1j = segj->GetStopCharLen() ;
                } else {
                    clen0j = segj->GetStopCharLen() ;
                    clen1j = segj->GetStartCharLen() ;
                }

                double avg0 = 0.5 * (clen0i + clen1j) ;
                double avg1 = 0.5 * (clen1i + clen0j) ;

                if (iorient) {
                    segi->SetStartCharLen(avg0) ;
                    segi->SetStopCharLen(avg1) ;
                } else {
                    segi->SetStopCharLen(avg0) ;
                    segi->SetStartCharLen(avg1) ;
                }
                if (jorient) {
                    segj->SetStartCharLen(avg1) ;
                    segj->SetStopCharLen(avg0) ;
                } else {
                    segj->SetStopCharLen(avg1) ;
                    segj->SetStartCharLen(avg0) ;
                }
            }

            ++iterj ;
        }
    }

    // now divide the edges

    CArbHashTableIterator<EdgeKey,CArbBdrySeg *>
            hiter(&sub_div_edges) ;
    for (hiter.First() ; hiter.More() ; ++hiter) {
        DivideSplineEdge(hiter.Key().id0,
                         hiter.Key().id1,
                         *(hiter.Entry())) ;
    }
}


void CArbRmshRegion2D::DivideSplineEdge(
                int start,int /*stop*/,
                CArbBdrySeg *spline)
{
    int i ;

    // add the nodes after checking
    // for the edge orientation

    int last,final ;
    double clen0,clen1,ratio ;

    if (start == spline->GetStartNode()) {
        clen0 = spline->GetStartCharLen() ;
        clen1 = spline->GetStopCharLen() ;
    } else {
        clen0 = spline->GetStopCharLen() ;
        clen1 = spline->GetStartCharLen() ;
    }

    ratio = spline->GetStartCharLen() /
            spline->GetStopCharLen() ;
    last = spline->GetStartNode() ;
    final = spline->GetStopNode() ;

    double sp_len = spline->ApproxLength() ;
    int num = LinearNumber(sp_len,clen0,clen1) ;

    // generate the nodes.  If this is a quadratic order mesh
    // add mid-side nodes at the average coordinate of the
    // adjacent nodes

    CArbCoord2D cd0, cd1 ;
    cd0 = spline->Evaluate(0.0) ;

    for (i=0 ; i<num-1 ; ++i) {
        double u = LinearPosition(i,1.0,num,ratio) ;
        cd1 = spline->Evaluate(u) ;
        int new_node = Owner->NewNodeId() ;
        Owner->NodeTable->Store(new_node,cd1) ;
        CornerNodes->Store(new_node,1) ;
        Boundary->SplitEdge(last,final,new_node) ;

        if (Owner->Order == QUADRATIC) {
            CArbCoord2D mid = 0.5 * (cd0 + cd1) ;
            int mid_node = Owner->NewNodeId() ;
            Owner->NodeTable->Store(mid_node,mid) ;
            Boundary->SplitEdge(new_node,last,mid_node) ;
        }

        last = new_node ;
        cd0 = cd1 ;
    }

    if (Owner->Order == QUADRATIC) {
        cd1 = spline->Evaluate(1.0) ;
        CArbCoord2D mid = 0.5 * (cd0 + cd1) ;
        int mid_node = Owner->NewNodeId() ;
        Owner->NodeTable->Store(mid_node,mid) ;
        Boundary->SplitEdge(last,final,mid_node) ;
    }
}


void CArbRmshRegion2D::DebugPrintBoundaryInfo()
{
    PlotUpdate() ;
    fflush(stdout) ;

    // print a list of the elements in the boundary

    CArbTopoElemIterator elem_iter(Boundary) ;
    int num = 0 ;
    for (elem_iter.First() ; elem_iter.More() ; ++elem_iter) {
        num++ ;
        bool ignore = false ;
        if (IgnoreElems->Fetch(*elem_iter) != 0) ignore = true ;
        fprintf(stderr,"%d %s\n",*elem_iter,(ignore?"I":" ")) ;
    }
    fprintf(stderr,"Num: %d\n\n",num) ;

    // for each element list the nodes

    for (elem_iter.First() ; elem_iter.More() ; ++elem_iter) {
        bool ignore = false ;
        if (IgnoreElems->Fetch(*elem_iter) != 0) ignore = true ;
        fprintf(stderr,"%d %s\n",*elem_iter,(ignore?"I":" ")) ;
        if (IgnoreElems->Fetch(*elem_iter) != 0) continue ;

        CArbHashTable<int,int> processed_nodes ;

        // loop through all the vertices in the boundary

        CArbTopoVtxIterator viter(*Boundary) ;
        for (viter.First() ; viter.More() ; ++viter) {

            // if we've seen this node already move on

            if (processed_nodes.Fetch(*viter) != 0) continue ;

            // if this is a quadratic order model and this
            // is not a corner node, then move on

            if ((Owner->Order == QUADRATIC) &&
                (CornerNodes->Fetch(*viter) == 0)) continue ;

            // look to see if this node is adjacent
            // to region we are currently meshing

            CArbTopoAdjVtxIterator aiter(Boundary,*viter) ;
            for (aiter.First() ; aiter.More() ; ++aiter) {

                if (aiter.CcwElem() == *elem_iter) {

                    // we are adjacent to the region so add
                    // this loop to the remesh zone. First nodes

                    int start = *viter ;
                    int cur = start ;
                    int adj = aiter.AdjVtx() ;
                    double first = true ;

                     while (first || (cur != start)) {
                        first = false ;
                        processed_nodes.Store(adj,1) ;
                        CArbCoord2D *ndc =
                            Owner->NodeTable->Fetch(adj) ;
                        fprintf(stderr,"    %d %g %g\n",
                                adj,ndc->x(),ndc->y()) ;
                        int nxt = Boundary->GetCCWNode(adj,cur) ;
                        adj = cur ;
                        cur = nxt ;
                    }

                }
            }
            processed_nodes.Store(*viter,1) ;
        }
    }
}


extern int print_num ;
static int curr = -1 ;

void CArbRmshRegion2D::MeshRegions()
{
    // DebugPrintBoundaryInfo() ;

    // mesh the remesh regions

    CArbTopoElemIterator elem_iter(Boundary) ;
    for (elem_iter.First() ; elem_iter.More() ; ++elem_iter) {

        ++curr ;
        //print_num = (curr == 26) ? 128 : 10000000 ;
// if (curr == 38) print_num = 0 ;
// if (curr == 53) print_num = 10 ;

        if (IgnoreElems->Fetch(*elem_iter) != 0) continue ;

        ArbMshElemType type =
            (Owner->ElemShapeType == CArbMshCrackRegion2D::S_TRIANGLE) ?
                              TRIANGLE : QUADRILATERAL ;
        CArbMshRegion2D *region =
            new CArbMshRegion2D(Owner->Order,type) ;

// if (curr == 3) region->SetElemType(TRIANGLE) ;
// if (curr == 4) region->SetElemType(TRIANGLE) ;
// if (curr == 12) region->SetElemType(TRIANGLE) ;
// if (curr == 13) region->SetElemType(TRIANGLE) ;

        region->SetDebugDisplayFlags(
            CArbMshRegion2D::DDFlags(DebugDisplayFlags)) ;

        region->SetBoundaryValidityChecks(false) ;
        region->SetMaterialID(MatTable->FetchValue(*elem_iter)) ;
        if (DoQuadAngleChecks) {
            CArbMshRegion2D::AngleCheckLocation loc =
                QuadAngleLocation == AT_NODES ?
                CArbMshRegion2D::AT_NODES :
                CArbMshRegion2D::AT_INT_PTS ;
            region->SetQuadElemAngleChecks(MinQuadAngle,MaxQuadAngle,loc) ;
        }

        // add the boundary information

        // this remesh region may have more than one loop.
        // Unfortunately, the current version of the topology
        // manager does not deal with multiply connected regions.
        // This means that we need to go hunt for the loops.

        CArbHashTable<int,int> processed_nodes ;

        // loop through all the vertices in the boundary

        CArbTopoVtxIterator viter(*Boundary) ;
        for (viter.First() ; viter.More() ; ++viter) {

            // if we've seen this node already move on

            if (processed_nodes.Fetch(*viter) != 0) continue ;

            // if this is a quadratic order model and this
            // is not a corner node, then move on

            if ((Owner->Order == QUADRATIC) &&
                (CornerNodes->Fetch(*viter) == 0)) continue ;

            // look to see if this node is adjacent
            // to region we are currently meshing

            CArbTopoAdjVtxIterator aiter(Boundary,*viter) ;
            for (aiter.First() ; aiter.More() ; ++aiter) {

                if (aiter.CcwElem() == *elem_iter) {

                    // we are adjacent to the region so add
                    // this loop to the remesh zone. First nodes

                    int start = *viter ;
                    int cur = start ;
                    int adj = aiter.AdjVtx() ;
                    double first = true ;

                     while (first || (cur != start)) {
                        first = false ;
                        processed_nodes.Store(adj,1) ;
                        CArbCoord2D *ndc =
                            Owner->NodeTable->Fetch(adj) ;
                        ArbMshNode node =
                            {adj,{ndc->x(),ndc->y(),0.0},false} ;
                        region->AddNode(node) ;
                        int nxt = Boundary->GetCCWNode(adj,cur) ;
                        adj = cur ;
                        cur = nxt ;
                    }

                    // now add the edges (move around the loop
                    // in clockwise order)

                    if (Owner->Order == LINEAR) {

                        start = *viter ;
                        cur = start ;
                        adj = aiter.AdjVtx() ;
                        first = true ;

                        while (first || (cur != start)) {
                            first = false ;
                            ArbMshEdge edge = {{cur,adj,-1}} ;
                            region->AddEdge(edge) ;
                            int nxt = Boundary->GetCCWNode(adj,cur) ;
                            adj = cur ;
                            cur = nxt ;
                        }

                    } else {

                        start = *viter ;
                        cur = start ;
                        int mid = Boundary->GetCCWNode(
                                     aiter.AdjVtx(),start) ;
                        adj = Boundary->GetCCWNode(cur,mid) ;
                        first = true ;

                        while (first || (cur != start)) {
                            first = false ;
                            ArbMshEdge edge = {{adj,cur,mid}} ;
                            region->AddEdge(edge) ;
                            mid = Boundary->GetCCWNode(mid,adj) ;
                            cur = adj ;
                            adj = Boundary->GetCCWNode(cur,mid) ;
                        }

                    }
                }
            }
            processed_nodes.Store(*viter,1) ;
        }

        // generate the mesh

        region->SetStartElemID(Owner->MaxElemId+1) ;
        region->SetStartNodeID(Owner->MaxNodeId+1) ;
        region->GenerateMesh() ;

        // extract the new mesh information and insert it into
        // the existing mesh

        CArbMshReg2DNodeIterator niter(*region,INTERIOR) ;
        for (niter.First() ; niter.More() ; ++niter) {
            ArbMshNode tmp ;
            tmp.id = niter.Id() ;
            tmp.coord[0] = niter.NCoordI(0) ;
            tmp.coord[1] = niter.NCoordI(1) ;
            tmp.coord[2] = 0.0 ;
            Owner->AddNode(tmp) ;
        }

        CArbMshReg2DElemIterator eiter(*region) ;
        for (eiter.First() ; eiter.More() ; ++eiter) {
            Owner->AddElem(*eiter.MshElement()) ;
            CreatedElems->InsertAtEnd(eiter.Id()) ;
            int lnum = (eiter.NumNodes() == 3) ||
                       (eiter.NumNodes() == 6) ? 3 : 4 ;
            for (int i=0 ; i<lnum ; ++i)
                CornerNodes->Store(eiter.Node(i),1) ;
        }

        delete region ;
    }
}


int CArbRmshRegion2D::FindCWvtx(int node_id,
                                     const CArbCoord2D &pt0,
                                     const CArbCoord2D &pt1)
{
    // This function does the geometry check necessary to
    // determine which boundary edge will be clockwise to
    // an new straight edge to be added to the boundary.  The
    // new straight edge starts at pt0 and ends at pt1.

    // first get a normalized version of the new edge

    CArbCoord2D norm = (pt1-pt0).Normalize() ;

    // loop through adjacent edges in the boundary

    int min_node = -1 ;
    double min_score = .0f;
    bool first = true ;
    bool parallel = false ;
    int num = 0 ;

    CArbTopoAdjVtxIterator iter(Boundary,node_id) ;
    for (iter.First() ; iter.More() ; ++iter) {

        // get the tangent at this point

        CArbCoord2D adj = AdjacentTangent(node_id,iter.AdjVtx()) ;

        // rank this vector and check to see if this is the smallest
        // value we have seen so far

        double score = ScoreRelativeAngle(norm,adj) ;

        // check for the special case of an edge along the check
        // direction

        if ((score < -.9999999) || (score > 0.9999999)) {
            parallel = true ;
        }

        if (first || (score < min_score)) {
            first = false ;
            min_node = iter.AdjVtx() ;
            min_score = score ;
        }
        ++num ;
    }

    if (parallel) return(num == 1 ? min_node : -2) ;
    return(min_node) ;
}


int CArbRmshRegion2D::FindDplElem(int node_id,
                                  const CArbCoord2D &cur,
                                  const CArbCoord2D &prev,
                                  const CArbCoord2D &next,
                                  CArbCoord2D *perp) const
{
    // This function does the geometry check necessary when
    // duplicating a node on the boundary because a flaw touches
    // the boundary at one point.  It finds the boundary edges
    // before and after the flaw vectors.  It returns the element
    // that is ccw to the before vector.  This element id is
    // needed when the vertex is split.  It also determines the
    // bisector to the angle between the before and after vectors.
    // The direction perpendicular to this is returned as the
    // the direction of the edge inserted in the split vertex

    // first get a normalized version of the incoming and outgoing
    // vectors

    CArbCoord2D incoming = (prev-cur).Normalize() ;
    CArbCoord2D outgoing = (next-cur).Normalize() ;

    // loop through adjacent edges in the boundary

    int before_elem = 0;
    double min_score = .0f, max_score = .0f;
    CArbCoord2D before_dir,after_dir ;
    bool first = true ;
    int num = 0 ;

    CArbTopoAdjVtxIterator iter(Boundary,node_id) ;
    for (iter.First() ; iter.More() ; ++iter) {

        // get the tangent at this point

        CArbCoord2D adj = AdjacentTangent(node_id,iter.AdjVtx()) ;

        // rank this vector and check to see if this is the
        // smallest value we have seen so far

        double before_score = ScoreRelativeAngle(outgoing,adj) ;
        double after_score = ScoreRelativeAngle(incoming,adj) ;

        if (first) {

            first = false ;
            if (iter.CcwElem() != NO_ELEM)
                before_elem = iter.CcwElem() ;
            min_score = before_score ;
            max_score = after_score ;
            before_dir = adj ;
            after_dir = adj ;

        } else {

            if (before_score < min_score) {
                if (iter.CcwElem() != NO_ELEM)
                    before_elem = iter.CcwElem() ;
                min_score = before_score ;
                before_dir = adj ;
            }

            if (after_score > max_score) {
                max_score = after_score ;
                after_dir = adj ;
            }
        }
        ++num ;
    }

    // now find the tangent vector for the split vertex edge.  First
    // the degenerate case of two coincident edges

    if (fabs(max_score-min_score) < 0.00001) {
        (*perp)[0] = -before_dir[1] ;
        (*perp)[1] = before_dir[0] ;
    } else {
        CArbCoord2D tmp ;
        tmp = BisectNormVect(before_dir,after_dir) ;
        (*perp)[0] = tmp[1] ;
        (*perp)[1] = -tmp[0] ;
    }

    return(before_elem) ;
}


int CArbRmshRegion2D::FindAdjElem(int node_id,
                                  const CArbCoord2D &pt0,
                                  const CArbCoord2D &pt1) const
{
    // Given a node id and a ray eminating from the point,
    // this function does the geometry check necessary to
    // determine which element (face) the ray crosses

    CArbCoord2D norm = (pt1-pt0).Normalize() ;

    // loop through adjacent edges in the boundary

    int min_node = -1 ;
    double min_score = .0f;
    int min_elem = 0;
    bool first = true ;
    bool parallel = false ;
    int num = 0 ;

    CArbTopoAdjVtxIterator iter(Boundary,node_id) ;
    for (iter.First() ; iter.More() ; ++iter) {

        // get the tangent at this point

        CArbCoord2D adj = AdjacentTangent(node_id,iter.AdjVtx()) ;

        // rank this vector and check to see if this is the smallest
        // value we have seen so far

        double score = ScoreRelativeAngle(norm,adj) ;

        // check for the special case of an edge along the check
        // direction

        if ((score < -.9999999) || (score > 0.9999999)) {
            parallel = true ;
        }

        if (first || (score < min_score)) {
            first = false ;
            min_node = iter.AdjVtx() ;
            min_elem = iter.CcwElem() ;
            min_score = score ;
        }
        ++num ;
    }

    if (parallel) return(num == 1 ? min_node : -2) ;
    return(min_elem) ;
}


CArbCoord2D CArbRmshRegion2D::AdjacentTangent(int vtx,int adj) const
{
    // this routine returns the tangent of the edge joining
    // the vertices vtx and adj evaluated at vtx

    CArbBdrySeg **ptr = BoundaryReps->Fetch(EdgeKey(vtx,adj)) ;

    if (ptr != 0) {
        CArbBdrySeg *spline = *ptr ;
        if (vtx == spline->GetStartNode()) {
            return(spline->Tangent(0.0)) ;
        } else {
            return(spline->Tangent(1.0)) ;
        }
    } else {
        CArbCoord2D pt1 = Owner->NodeTable->FetchValue(vtx) ;
        CArbCoord2D pt2 = Owner->NodeTable->FetchValue(adj) ;
        return((pt2-pt1).Normalize()) ;
    }
}



int CArbRmshRegion2D::AddIntersectionPoint(
                          CArbBdrySeg *spline,
                          double spline_par,double line_par,
                          const LocCrackPt &prev,const LocCrackPt &next,
                          CArbArray<LocCrackPt> *NewPoints)
{
    // this function adds a new point to the remesh boundary at
    // the intersection of a crack edge and an external boundary
    // or bi-material interface

    // get a new node ID

    int new_node = Owner->NewNodeId() ;

    // create a new crack point and add it to the list

    LocCrackPt loc_pt ;
    loc_pt.coord = spline->Evaluate(spline_par) ;
    loc_pt.char_elem_size = prev.char_elem_size + line_par *
                    (next.char_elem_size - prev.char_elem_size) ;
    loc_pt.node_id = new_node ;
    loc_pt.type = spline->IsBiMatInterface() ? INTERFACE : BOUNDARY ;
    loc_pt.crack_tip = false ;
    loc_pt.parent_splines.Insert(spline->GetParentId()) ;

    NewPoints->InsertAtEnd(loc_pt) ;
    Owner->NodeTable->Store(new_node,loc_pt.coord) ;
    CornerNodes->Store(new_node,1) ;
    SplitSpline(spline,new_node,spline_par,
                loc_pt.char_elem_size,Boundary) ;
    return(new_node) ;
}

void CArbRmshRegion2D::SplitSpline(
                          CArbBdrySeg *spline,
                          int new_node,double split_par,
                          double new_char_len,
                          CArbMshTopo2D *boundary)
{
    if (boundary != 0) {
        boundary->SplitEdge(spline->GetStartNode(),
                            spline->GetStopNode(),new_node) ;
    }

    CArbBdrySeg *first_part = spline->SubSegment(
                                        0.0,split_par,
                                        spline->GetStartNode(),
                                        spline->GetStartCharLen(),
                                        new_node,new_char_len,
                                        spline->IsBiMatInterface()) ;

    CArbBdrySeg *scnd_part = spline->SubSegment(
                                        split_par,1.0,
                                        new_node,new_char_len,
                                        spline->GetStopNode(),
                                        spline->GetStopCharLen(),
                                        spline->IsBiMatInterface()) ;

    BoundaryReps->Remove(EdgeKey(spline->GetStartNode(),
                                 spline->GetStopNode())) ;

    BoundaryReps->Store(EdgeKey(first_part->GetStartNode(),
                                first_part->GetStopNode()),first_part) ;

    BoundaryReps->Store(EdgeKey(scnd_part->GetStartNode(),
                                scnd_part->GetStopNode()),scnd_part) ;

    delete spline ;
}


void CArbRmshRegion2D::FindAllIntersections()
{
    int cur, next, k, l ;

    // What we want to do here is to search for all possible
    // intersections of the (straight) flaw edges and the
    // remesh boundary.  If intersections are found then
    // new points are added to the crack point list

    CArbArray<LocCrackPt> *NewPoints = new CArbArray<LocCrackPt> ;
    CArbArray<IntscPoint> intsc_list ;

    for (cur=0 ; cur<LocalCrackPoints->NumEntries() ; ++cur) {
        next = (cur+1) % LocalCrackPoints->NumEntries() ;
        NewPoints->InsertAtEnd(LocalCrackPoints->At(cur)) ;

        // loop through all the boundary edges looking for an
        // intersection.  If we find any add them to a list
        // that we can then sort so that the new points are
        // added in the proper order

        intsc_list.Clear() ;

        CArbTopoElemIterator elem_iter(Boundary) ;
        for (elem_iter.First() ; elem_iter.More() ; ++elem_iter) {

            CArbTopoEdgeOnElemIterator edge_iter(Boundary,*elem_iter) ;
            for (edge_iter.First() ; edge_iter.More() ; ++edge_iter) {

                // check to see if we have a spline rep or not

                CArbBdrySeg **ptr = BoundaryReps->Fetch(
                        EdgeKey(edge_iter[0],edge_iter[1])) ;

                if (ptr != 0) {
                    CArbBdrySeg *spline = *ptr ;

                    // if this is a bi-material interface we will
                    // see the edge twice, so ignore it one of the
                    // times

                    if (spline->IsBiMatInterface()) {
                        if (edge_iter[0] > edge_iter[1]) continue ;
                    }

                    // otherwise look for intersections

                    CArbArray<double> *intsct =
                       spline->IntersectLine(
                                 LocalCrackPoints->At(cur).coord,
                                 LocalCrackPoints->At(next).coord) ;

                    // if this flaw edge intersects this boundary
                    // edge (at least once) store this information
                    // so that we can process the information later

                    if (intsct != 0) {
                        for (int k=0 ; k<intsct->NumEntries() ; k+=2) {
                            if ((intsct->At(k) > 0.01) &&
                                (intsct->At(k) < 0.99) &&
                                (intsct->At(k+1) > 0.01) &&
                                (intsct->At(k+1) < 0.99)) {

                                IntscPoint tmp ;
                                tmp.spline = spline ;
                                tmp.spline_par = intsct->At(k) ;
                                tmp.line_par = intsct->At(k+1) ;
                                tmp.pt_indx = -1 ;
                                intsc_list.InsertAtEnd(tmp) ;
                            }
                        }
                        delete intsct ;
                    }
                }
            }
        }

        // Loop through the flaw points and see if we pass through
        // any of those

        for (k=0 ; k<NewPoints->NumEntries() ; ++k) {
            double param ;
            if (IsPointOnLine(NewPoints->At(k).coord,
                              LocalCrackPoints->At(cur).coord,
                              LocalCrackPoints->At(next).coord,
                              &param)) {
                if ((param > 0.01) && (param < 0.99)) {
                    IntscPoint tmp ;
                    tmp.spline = 0 ;
                    tmp.spline_par = 0 ;
                    tmp.line_par = param ;
                    tmp.pt_indx = k ;
                    intsc_list.InsertAtEnd(tmp) ;
                }
            }
        }

        // Sort the intersections along the flaw edge vector
        // so that they can be added in the proper order.  We use
        // a bubble sort here because we expect this list to be
        // very short

        for (k=0 ; k<intsc_list.NumEntries()-1 ; ++k) {
            for (l=k ; l<intsc_list.NumEntries() ; ++l) {
                if (intsc_list[k].line_par > intsc_list[l].line_par) {
                    IntscPoint tmp = intsc_list[k] ;
                    intsc_list[k] = intsc_list[l] ;
                    intsc_list[l] = tmp ;
                }
            }
        }

        // now do the intersection

        int new_node = 0, first = 0, second = 0 ;
        for (l=0 ; l<intsc_list.NumEntries() ; ++l) {

            if (intsc_list[l].spline != 0) {

                // if this is not the first time through then this may
                // be the second intersection of this spline.  In that
                // case we need to determine which segment of the
                // spline to split

                CArbBdrySeg *spline ;
                CArbArray<double> *intsct ;

                if ((l > 0) &&
                    (intsc_list[l].spline == intsc_list[l-1].spline)) {

                    CArbBdrySeg **ptr ;

                    if (intsc_list[l].spline_par <
                        intsc_list[l-1].spline_par) {
                        ptr = BoundaryReps->Fetch(EdgeKey(first,new_node)) ;
                    } else {
                        ptr = BoundaryReps->Fetch(EdgeKey(second,new_node)) ;
                    }

                    spline = *ptr ;
                    intsct = spline->IntersectLine(
                                 LocalCrackPoints->At(cur).coord,
                                 LocalCrackPoints->At(next).coord) ;
                } else {
                    spline = intsc_list[l].spline ;
                    intsct = spline->IntersectLine(
                                 LocalCrackPoints->At(cur).coord,
                                 LocalCrackPoints->At(next).coord) ;
                }

                first  = intsc_list[l].spline->GetStartNode() ;
                second = intsc_list[l].spline->GetStopNode() ;

                for (int k=0 ; k<intsct->NumEntries() ; k+=2) {
                    if ((intsct->At(k) > 0.01) &&
                        (intsct->At(k) < 0.99) &&
                        (intsct->At(k+1) > 0.01) &&
                        (intsct->At(k+1) < 0.99)) {

                        new_node = AddIntersectionPoint(
                                spline,intsct->At(k),intsct->At(k+1),
                                LocalCrackPoints->At(cur),
                                LocalCrackPoints->At(next),
                                NewPoints) ;
                        break ;
                    }
                }

                if (intsct != 0) delete intsct ;

            } else {

                // Find the face that contains the previous part
                // of the flaw leaving the point

                int indx = intsc_list[l].pt_indx ;
                LocCrackPt &this_pt = NewPoints->At(indx) ;

                int start_elem = FindAdjElem(this_pt.node_id,
                                       NewPoints->At(indx).coord,
                                       NewPoints->At(indx+1).coord) ;
                int stop_elem = FindAdjElem(this_pt.node_id,
                                       NewPoints->At(indx).coord,
                                       NewPoints->At(indx-1).coord) ;

                CArbCoord2D delt = (NewPoints->At(indx+1).coord -
                                    NewPoints->At(indx).coord).Normalize() ;
                CArbCoord2D tang(-delt.y(),delt.x()) ;

                // duplicate the node

                int new_node = SplitFlawNode(this_pt,start_elem,
                                             stop_elem,tang) ;

                // add this point

                LocCrackPt tmp ;
                tmp = this_pt ;
                tmp.node_id = new_node ;
                NewPoints->InsertAtEnd(tmp) ;

            }
        }
    }
    delete LocalCrackPoints ;
    LocalCrackPoints = NewPoints ;
    NumPoints = LocalCrackPoints->NumEntries() ;
}


#define KINK_TOL 0.349

void CArbRmshRegion2D::RebuildRemeshTopology()
{
    int prev_id, cur_id, next_id, nchk_id ;
    CArbCoord2D prev, cur, next ;
    double min_angle = PI - KINK_TOL ;
    double max_angle = PI + KINK_TOL ;
    int spline_id = 1 ;

// display_topo_edge(Boundary,Owner->NodeTable) ;

    // The current remesh boundary consists of current and
    // former element edges that are adjacent to the remesh
    // region.  These edges fall into one of three categories,
    // those still attached to an element, those that are part
    // of the external boundary, and those that are part of
    // a bi-material interfaoe.

    // This routine transforms the current remesh boundary to
    // one where all the boundary and interfaces segments are
    // collapsed into logical segments that span the regions
    // between "kink" angles.  These segments will be represented
    // with splines.

    CArbMshTopo2D *NewBoundary = new CArbMshTopo2D ;
    BoundaryReps = new CArbHashTable<EdgeKey,CArbBdrySeg*> ;

    CArbTopoElemIterator elem_iter(Boundary) ;
    for (elem_iter.First() ; elem_iter.More() ; ++elem_iter) {

        CArbTopoLoopOnElemIterator loop_iter(Boundary,*elem_iter) ;
        for (loop_iter.First() ; loop_iter.More() ; ++loop_iter) {

            int last[2] = {-1,-1} ;

            CArbArray<int> corners ;
            CArbTopoEdgeOnElemCyclicIterator
                edge_iter(Boundary,*elem_iter,*loop_iter) ;

            // get the first three nodes on the boundary

            if ((Owner->Order == QUADRATIC) &&
                (CornerNodes->Fetch(edge_iter[0]) == 0)) ++edge_iter ;

            prev_id = edge_iter[0] ;
            if (Owner->Order == QUADRATIC) {
                ++edge_iter ;
                cur_id = edge_iter[1] ;
            } else {
                cur_id = edge_iter[1] ;
            }
            prev = *(Owner->NodeTable->Fetch(prev_id)) ;
            cur  = *(Owner->NodeTable->Fetch(cur_id)) ;
            ++edge_iter ;
            if (Owner->Order == QUADRATIC) {
                nchk_id = edge_iter[1] ;
                ++edge_iter ;
                next_id  = edge_iter[1] ;
            } else {
                nchk_id = next_id  = edge_iter[1] ;
            }
            next = *(Owner->NodeTable->Fetch(next_id)) ;

            int other_elem ;
            bool done = false ;
            bool have_first_corner = false ;
            int first_node = cur_id ;

            // loop through the current boundary segments building
            // up an array of the corner node id's

            while (!done) {

                if (have_first_corner) corners.InsertAtEnd(cur_id) ;
                PointType edge_type = ClassifyEdge(cur_id,nchk_id,
                                               &other_elem) ;

                if (edge_type == INSIDE) {

                    // if this segment is the edge of an existing
                    // element then add it to the new boundary

                    prev_id = cur_id ;  cur_id  = next_id ;
                    prev = cur ;  cur = next ;
                    ++edge_iter ;
                    if (Owner->Order == QUADRATIC) {
                        if (have_first_corner)
                            corners.InsertAtEnd(nchk_id) ;
                        nchk_id = edge_iter[1] ;
                        ++edge_iter ;
                        next_id  = edge_iter[1] ;
                    } else {
                        nchk_id = next_id = edge_iter[1] ;
                    }
                    next = *(Owner->NodeTable->Fetch(next_id)) ;

                } else if (edge_type == BOUNDARY) {

                    // if this is part of an external boundary then
                    // move until we get to a kink or a change in type,
                    // keeping track of the node points along the way

                    int start_node,stop_node ;
                    CArbArray<CArbCoord2D> coords ;
                    CArbArray<int> nodeid ;
                    coords.InsertAtEnd(cur) ;
                    nodeid.InsertAtEnd(cur_id) ;
                    coords.InsertAtEnd(next) ;
                    nodeid.InsertAtEnd(next_id) ;

                    start_node = cur_id ;
                    stop_node = next_id ;
                    prev_id = cur_id ;  cur_id  = next_id ;
                    prev = cur ;  cur = next ;
                    ++edge_iter ;
                    if (Owner->Order == QUADRATIC) {
                        nchk_id = edge_iter[1] ;
                        ++edge_iter ;
                    } else {
                        nchk_id = edge_iter[1] ;
                    }

                    while (1) {
                        next_id = edge_iter[1] ;
                        PointType edge_type = ClassifyEdge(cur_id,nchk_id,
                                                   &other_elem) ;
                        next = *(Owner->NodeTable->Fetch(next_id)) ;

                        if ((edge_type == INSIDE) ||
                            (edge_type == INTERFACE)) break ;

                        double angle = Angle2Pi(cur,next,prev) ;
                        if ((angle < min_angle) || (angle > max_angle)) break ;

                        // deal with the special case of a closed
                        // loop boundary with no kinks

                        if ((!have_first_corner) &&
                            (next_id == first_node)) {
                            corners.InsertAtEnd(next_id) ;
                            coords.InsertAtEnd(next) ;
                            nodeid.InsertAtEnd(next_id) ;
                            start_node = next_id ;
                            stop_node = next_id ;
                            have_first_corner = true ;
                            done = true ;
                            break ;
                        }

                        coords.InsertAtEnd(next) ;
                        nodeid.InsertAtEnd(next_id) ;

                        stop_node = next_id ;
                        prev_id = cur_id ;  cur_id  = next_id ;
                        prev = cur ;  cur = next ;
                        ++edge_iter ;
                        if (Owner->Order == QUADRATIC) {
                            nchk_id = edge_iter[1] ;
                            ++edge_iter ;
                        } else {
                            nchk_id = edge_iter[1] ;
                        }
                    }

                    // Build the spline representation for this edge
                    // and store it away

                    if (have_first_corner) {
                        CArbCoord2D *pts = coords.AsVector() ;
                        CArbBdrySeg *b_spline = new CArbBdrySeg(
                                      coords.NumEntries(),pts,
                                      start_node,
                                      CharLength(start_node,nodeid[1]),
                                      stop_node,
                                      CharLength(stop_node,nodeid[-2]),
                                      spline_id,false) ;
                        CArbBdrySeg **tmp = BoundaryReps->Fetch(
                               EdgeKey(start_node,stop_node)) ;
                        if (tmp != 0) delete *tmp ;
                        BoundaryReps->Store(
                            EdgeKey(start_node,stop_node),b_spline) ;
                        ++spline_id ;
                        delete [] pts ;
                    }

                } else if (edge_type == INTERFACE) {

                    // if this is part of a bi-material interface then
                    // move until we get to a kink, change in type, or
                    // change of element on the other side of the interface

                    int start_node,stop_node ;
                    CArbArray<CArbCoord2D> coords ;
                    CArbArray<int> nodeid ;
                    coords.InsertAtEnd(cur) ;
                    nodeid.InsertAtEnd(cur_id) ;
                    coords.InsertAtEnd(next) ;
                    nodeid.InsertAtEnd(next_id) ;

                    start_node = cur_id ;
                    stop_node = next_id ;
                    prev_id = cur_id ;  cur_id  = next_id ;
                    prev = cur ;  cur = next ;
                    ++edge_iter ;
                    if (Owner->Order == QUADRATIC) {
                        nchk_id = edge_iter[1] ;
                        ++edge_iter ;
                    } else {
                        nchk_id = edge_iter[1] ;
                    }

                    while (1) {
                        next_id = edge_iter[1] ;
                        int new_other_elem ;
                        PointType edge_type = ClassifyEdge(cur_id,nchk_id,
                                                       &new_other_elem) ;
                        next = *(Owner->NodeTable->Fetch(next_id)) ;

                        if ((edge_type == INSIDE) ||
                            (edge_type == BOUNDARY)) break ;

                        double angle = Angle2Pi(cur,next,prev) ;
                        if ((angle < min_angle) || (angle > max_angle)) break ;

                        if (other_elem != new_other_elem) break ;

                        coords.InsertAtEnd(next) ;
                        nodeid.InsertAtEnd(next_id) ;

                        stop_node = next_id ;
                        prev_id = cur_id ;  cur_id  = next_id ;
                        prev = cur ;  cur = next ;
                        ++edge_iter ;
                        if (Owner->Order == QUADRATIC) {
                            nchk_id = edge_iter[1] ;
                            ++edge_iter ;
                        } else {
                            nchk_id = edge_iter[1] ;
                        }
                    }

                    // Build the spline representation for this edge
                    // and store it away

                    if (have_first_corner) {
                        CArbCoord2D *pts = coords.AsVector() ;
                        CArbBdrySeg *b_spline =
                            new CArbBdrySeg(coords.NumEntries(),pts,
                                  start_node,
                                  CharLength(start_node,nodeid[1]),
                                  stop_node,
                                  CharLength(stop_node,nodeid[-2]),
                                  spline_id,true) ;

                        // deal with the special case of an internal
                        // void with only two nodes (and two edges).
                        // If this is the case then we split both
                        // splines

                        CArbBdrySeg **tmp = BoundaryReps->Fetch(
                               EdgeKey(start_node,stop_node)) ;
                        if ((start_node == last[1]) &&
                            (stop_node == last[0])) {

                            CArbBdrySeg *old_spline = *tmp ;

                            int new_node = Owner->NewNodeId() ;
                            Owner->NodeTable->Store(new_node,
                                old_spline->Evaluate(0.5)) ;
                            CornerNodes->Store(new_node,1) ;

                            SplitSpline(old_spline,new_node,0.5,
                                0.5*(old_spline->GetStartCharLen() +
                                     old_spline->GetStopCharLen()),0) ;

                            corners.InsertAtEnd(corners[1]) ;
                            corners[1] = new_node ;

                            BoundaryReps->Store(
                                EdgeKey(start_node,stop_node),
                                b_spline) ;

                            new_node = Owner->NewNodeId() ;
                            Owner->NodeTable->Store(new_node,
                                b_spline->Evaluate(0.5)) ;
                            CornerNodes->Store(new_node,1) ;

                            SplitSpline(b_spline,new_node,0.5,
                                0.5*(b_spline->GetStartCharLen() +
                                     b_spline->GetStopCharLen()),0) ;

                            corners.InsertAtEnd(new_node) ;
                        } else {
                            if (tmp != 0) delete *tmp ;
                            BoundaryReps->Store(
                                EdgeKey(start_node,stop_node),
                                b_spline) ;
                        }
                        ++spline_id ;
                        delete [] pts ;
                        last[0] = start_node ;
                        last[1] = stop_node ;
                    }
                }

                if (!have_first_corner) {
                    first_node = cur_id ;
                    have_first_corner = true ;
                } else {
                    if (cur_id == first_node) done = true ;
                }
            }

            // insert segments into the new boundary

            int *corner_array = corners.AsVector() ;
            NewBoundary->InsertElement(*elem_iter,
                    corners.NumEntries(),corner_array) ;
            delete [] corner_array ;
        }
    }

    delete Boundary ;
    Boundary = NewBoundary ;

// #ifdef DEBUG
//     display_topo_edge(Boundary,Owner->NodeTable) ;
// #endif
}

void CArbRmshRegion2D::SetFlawCharSizes() {}


CArbRmshRegion2D::PointType
CArbRmshRegion2D::ClassifyEdge(
                   int nd0,int nd1,int *other_elem) const
{
    // first check to see if there is an element on the other
    // side of this edge in the remesh topology.  If so, this
    // is an inteface edge.

    *other_elem = Boundary->GetCCWElem(nd1,nd0) ;
    if (*other_elem != NO_ELEM) return(INTERFACE) ;

    // now check to see if there is an element on the other
    // side of this edge in the mesh.  If so, this is an
    // "INSIDE" edge.  Otherwise it is a boundary edge

    *other_elem = Owner->MshTopo->GetCCWElem(nd1,nd0) ;
    if (*other_elem != NO_ELEM) return(INSIDE) ;

    // check to see if both of these nodes are retained nodes

    if ((Owner->RetainedNodes->Fetch(nd0) != 0) &&
        (Owner->RetainedNodes->Fetch(nd1) != 0)) return(INSIDE) ;

    // else we are on the boundary

    return(BOUNDARY) ;
}


// The MIN_ANISOTROPIC_RATIO parameter is used when doing variable
// characteristic length.  If the length of the current old element
// edge length divided by the mimium old element length is greater
// than this value, then anisotropic characteristic lengths are
// considered.

#define MIN_ANISOTROPIC_RATIO 4.0

// The MAX_ANISOTROPIC_ANGLE parameter is used when doing variable
// characteristic length.  An anisotropic characteristic value
// is considered if the angle between the current segment and the
// segment with the minium value is less than this angle.

#define MAX_ANISOTROPIC_ANGLE 2.3562       // 3*PI/4


double CArbRmshRegion2D::CharLength(int node,int other) const
{
    // find the characteristic length for the element
    // edges at a node

    double clen = 0.0 ;

    if (Owner->SegEndCharSizeType == CArbMshCrackRegion2D::S_MINIMUM) {

        double dist,min = 0.0 ;
        CArbCoord2D *coord = Owner->NodeTable->Fetch(node) ;
        CArbTopoAdjVtxIterator iter(Boundary,node) ;

        for (iter.First() ; iter.More() ; ++iter) {
            CArbCoord2D *acoord =
                Owner->NodeTable->Fetch(iter.AdjVtx()) ;
            dist = (*acoord - *coord).Magnitude() ;
            if ((min == 0.0) || ((dist < min) && (dist > 0.0))) {
                min = dist ;
            }
        }

        clen = min ;

    } else if (Owner->SegEndCharSizeType == CArbMshCrackRegion2D::S_AVERAGE) {

        double sum = 0.0 ;
        int num = 0 ;
        CArbCoord2D *coord = Owner->NodeTable->Fetch(node) ;
        CArbTopoAdjVtxIterator iter(Boundary,node) ;

        for (iter.First() ; iter.More() ; ++iter) {
            CArbCoord2D *acoord =
                Owner->NodeTable->Fetch(iter.AdjVtx()) ;
            sum += (*acoord - *coord).Magnitude() ;
            ++num ;
        }

        clen = num > 0 ? sum/num : 0.0 ;

    } else if (Owner->SegEndCharSizeType == CArbMshCrackRegion2D::S_ANISOTROPIC) {

        CArbCoord2D *coord = Owner->NodeTable->Fetch(node) ;
        CArbTopoAdjVtxIterator iter(Boundary,node) ;
        bool found = false ;

        for (iter.First() ; iter.More() ; ++iter) {
            if (Owner->Order == LINEAR) {
                if (iter.AdjVtx() == other) {
                    CArbCoord2D *acoord =
                        Owner->NodeTable->Fetch(iter.AdjVtx()) ;
                    clen = (*acoord - *coord).Magnitude() ;
                    break ;
                }
            } else {
                CArbTopoAdjVtxIterator niter(Boundary,iter.AdjVtx()) ;
                for (niter.First() ; niter.More() ; ++niter) {
                    if (niter.AdjVtx() == other) {
                        CArbCoord2D *acoord =
                            Owner->NodeTable->Fetch(niter.AdjVtx()) ;
                        clen = 0.5 * (*acoord - *coord).Magnitude() ;
                        found = true ;
                        break ;
                    }
                }
                if (found) break ;
            }
        }

    } else if (Owner->SegEndCharSizeType ==
               CArbMshCrackRegion2D::S_USER_ANISOTROPIC) {

        double min = 0.0, cur = 0.0 ;
        CArbCoord2D *coord = Owner->NodeTable->Fetch(node) ;
        CArbTopoAdjVtxIterator iter(Boundary,node) ;
        bool found = false ;
        int num = 0 ;

        for (iter.First() ; iter.More() ; ++iter) {
            ++num ;
            CArbCoord2D *acoord =
                Owner->NodeTable->Fetch(iter.AdjVtx()) ;

            // check for the minimum adjacent edge length

            CArbCoord2D dir = *acoord - *coord ;
            double dist = dir.Magnitude() ;
            if ((min == 0.0) || ((dist < min) && (dist > 0.0))) {
                min = dist ;
            }

            // if we have not yet found "this" edge check to see
            // if we have it

            if (found) continue ;

            if (iter.AdjVtx() == other) {
                cur = dist ;
                found = true ;
            }
        }

        // determine a characteristic length

        if (cur == min) {
            clen = min ;
        } else {
            clen = AnisotropicRatio * min ;
        }

    } else if (Owner->SegEndCharSizeType == CArbMshCrackRegion2D::S_VARIABLE) {

        // The basic approach with variable is that in general
        // we wast the minimum element size except for the case
        // where the angle between this edge and the minimum
        // edge is large and the ratio between the segments is
        // "large"

        double min = 0.0, cur = 0.0 ;
        CArbCoord2D *coord = Owner->NodeTable->Fetch(node) ;
        CArbCoord2D min_dir,cur_dir ;
        CArbTopoAdjVtxIterator iter(Boundary,node) ;
        bool found = false ;

        for (iter.First() ; iter.More() ; ++iter) {
            CArbCoord2D *acoord =
                Owner->NodeTable->Fetch(iter.AdjVtx()) ;
            CArbCoord2D dir = *acoord - *coord ;
            double dist = dir.Magnitude() ;
            if ((min == 0.0) || ((dist < min) && (dist > 0.0))) {
                min = dist ;
                min_dir = dir ;
            }

            // if we have not yet found "this" edge check to see
            // if we have it

            if (found) continue ;

            if (Owner->Order == LINEAR) {
                if (iter.AdjVtx() == other) {
                    cur = dist ;
                    cur_dir = dir ;
                    found = true ;
                }
            } else {   // quadratic, move to next node
                CArbTopoAdjVtxIterator niter(Boundary,iter.AdjVtx()) ;
                for (niter.First() ; niter.More() ; ++niter) {
                    if (niter.AdjVtx() == other) {
                        cur = 0.5 * dist ;
                        cur_dir = dir ;
                        found = true ;
                        break ;
                    }
                }
            }
        }

        // determine a characteristic length

        if (cur/min >= MIN_ANISOTROPIC_RATIO) {
            if (fabs(AngleVect(cur_dir,min_dir)) <
                      MAX_ANISOTROPIC_ANGLE) {
                clen = cur ;
            } else {
                clen = min ;
            }
        } else {
            clen = min ;
        }
    }

    clen *= Owner->SegEndCharSizeFactor ;
    return(Owner->Order == QUADRATIC ? 2.0*clen : clen) ;
}



void CArbRmshRegion2D::PrintBoundary()
{
    CArbTopoElemIterator elem_iter(Boundary) ;
    for (elem_iter.First() ; elem_iter.More() ; ++elem_iter) {

        fprintf(stderr,"\nElement: %d\n",*elem_iter) ;

        CArbTopoEdgeOnElemIterator edge_iter(Boundary,*elem_iter) ;
        for (edge_iter.First() ; edge_iter.More() ; ++edge_iter) {

            fprintf(stderr,"    %d %d\n",edge_iter[0],edge_iter[1]) ;

        }
    }
}


void CArbRmshRegion2D::PlotUpdate()
{
    for (int i=0 ; i<NumPoints ; ++i) {
        int j = (i+1) % NumPoints ;
        printf("t %g %g %d\n",LocalCrackPoints->At(i).coord.x(),
                              LocalCrackPoints->At(i).coord.y(),i) ;
        printf("l %g %g %g %g\n",LocalCrackPoints->At(i).coord.x(),
                                 LocalCrackPoints->At(i).coord.y(),
                                 LocalCrackPoints->At(j).coord.x(),
                                 LocalCrackPoints->At(j).coord.y()) ;
    }


    CArbTopoElemIterator elem_iter(Boundary) ;
    for (elem_iter.First() ; elem_iter.More() ; ++elem_iter) {

        CArbTopoLoopOnElemIterator loop_iter(Boundary,*elem_iter) ;
        for (loop_iter.First() ; loop_iter.More() ; ++loop_iter) {

            CArbTopoEdgeOnElemIterator
            edge_iter(Boundary,*elem_iter,*loop_iter) ;
            for (edge_iter.First() ; edge_iter.More() ; ++edge_iter) {

                // check to see if we have a spline rep or not

                CArbBdrySeg **ptr = 0 ;

                if (BoundaryReps != 0) {
                    ptr = BoundaryReps->Fetch(
                          EdgeKey(edge_iter[0],edge_iter[1])) ;
                }

                if (ptr == 0) {
                    CArbCoord2D *node0 =
                         Owner->NodeTable->Fetch(edge_iter[0]) ;
                    CArbCoord2D *node1 =
                         Owner->NodeTable->Fetch(edge_iter[1]) ;
                    printf("t %g %g %d\n",node0->x(),node0->y(),edge_iter[0]) ;
                    printf("l %g %g %g %g\n",node0->x(),node0->y(),
                                             node1->x(),node1->y()) ;
                } else {
                    CArbBdrySeg *spline = *ptr ;
                    double du = 0.05 ;
                    CArbCoord2D pt0 = spline->Evaluate(0.0) ;

                    printf("t %g %g %d\n",pt0.x(),pt0.y(),
                           spline->GetStartNode()) ;

                    for (double u = du ; u < 1.0 ; u += du) {
                        CArbCoord2D pt1 = spline->Evaluate(u) ;
                        printf("l %g %g %g %g\n",pt0.x(),pt0.y(),
                                             pt1.x(),pt1.y()) ;
                        pt0 = pt1 ;
                    }
                    CArbCoord2D pt1 = spline->Evaluate(1.0) ;
                    printf("l %g %g %g %g\n",pt0.x(),pt0.y(),
                                             pt1.x(),pt1.y()) ;
                }
            }
        }
    }
}


// --------------------------------------------------------------
// --------------------------------------------------------------
// --------------------------------------------------------------


// --------------------------------------------------------------
//  Support routines:
//

static void FindLineCoefs(const CArbCoord2D &A,
                          const CArbCoord2D &B,
                          double *coefs)
{
    CArbCoord2D delta = B - A ;
    coefs[0] = -delta[1] ;
    coefs[1] =  delta[0] ;
    coefs[2] = -(coefs[0]*A[0] + coefs[1]*A[1]) ;
}

inline double EvlLineEqn(const CArbCoord2D &pt,double *coefs)
{
    return(pt[0]*coefs[0] + pt[1]*coefs[1] + coefs[2]) ;
}




// %(CArbRmshRegion2D::BuildDeleteData-CArbMshCrackRegion2D::DeleteSegData-|*^const-int-const|-CArbMshCrackRegion2D::CrackPt-const|*-double-const|*)
/* ++ ----------------------------------------------------------
**
**    BuildDeleteData - compute element delete segments
**
**      CArbMshCrackRegion2D::DeleteSegData *BuildDeleteData(
**              const int                           num_pts,
**              const CArbMshCrackRegion2D::CrackPt *pts,
**              const double                        *insitu_size) const
**
**        num_pts     - (in)  number of crack points
**        pts         - (in)  crack-point list
**        insitu_size - (in)  array of characteristic sizes of existing
**                            elements
**
**      Description: This function compiles an array of descriptions
**          for delete segments for remeshing. This array will have
**          size num_points-1. This contains information needed to
**          delete elements near the crack. Ownership of this array
**          passes to the client, which must eventually call delete []
**
**      Return Value: an array of delete element segments
**
**
** -- */

CArbRmshRegion2D::DeleteSegData *CArbRmshRegion2D::BuildDeleteData(
        const double *insitu_size) const
{
    double factor = (Owner->ElemShapeType ==
                     CArbMshCrackRegion2D::S_TRIANGLE) ?
                    TriDeleteFactor : QuadDeleteFactor ;

    DeleteSegData *data = new DeleteSegData[NumPoints] ;

    for (int i=0 ; i<NumPoints ; ++i) {
        int j = (i+1) % NumPoints ;

        LocCrackPt &this_pt = LocalCrackPoints->At(i) ;
        LocCrackPt &next_pt = LocalCrackPoints->At(j) ;

        // first the end point and radius data

        double char_size0 = (this_pt.char_elem_size > insitu_size[i]) ?
                             this_pt.char_elem_size : insitu_size[i] ;
        double char_size1 = (next_pt.char_elem_size > insitu_size[j]) ?
                             next_pt.char_elem_size : insitu_size[j] ;
        data[i].end_coords[0][0] = this_pt.coord[0] ;
        data[i].end_coords[0][1] = this_pt.coord[1] ;
        data[i].end_coords[1][0] = next_pt.coord[0] ;
        data[i].end_coords[1][1] = next_pt.coord[1] ;
        data[i].end_radii[0] = char_size0 * factor ;
        data[i].end_radii[1] = char_size1 * factor ;
        data[i].end_radii_sqr[0] = data[i].end_radii[0] *
                                   data[i].end_radii[0] ;
        data[i].end_radii_sqr[1] = data[i].end_radii[1] *
                                   data[i].end_radii[1] ;

        // now if we think of  circles of radius "end_radii"
        // around each of the end points, we want to find
        // the two lines that are mutually tangent to these
        // two circles.

        double dx = data[i].end_coords[1][0]-data[i].end_coords[0][0] ;
        double dy = data[i].end_coords[1][1]-data[i].end_coords[0][1] ;
        double dr = data[i].end_radii[0] - data[i].end_radii[1] ;

        double alpha = atan2(dy,dx) ;
        double len = sqrt(dx*dx + dy*dy) ;
        double ratio = fabs(dr) / len ;
        double theta ;

        if (ratio > 1.0) {
            data[i].degenerate_case = true ;
            if (dr < 0) {
                data[i].end_coords[0][0] = next_pt.coord[0] ;
                data[i].end_coords[0][1] = next_pt.coord[1] ;
            }
        } else {
            data[i].degenerate_case = false ;
            if (dr >= 0) {
                theta = acos(ratio) ;
            } else {
                theta = PI - acos(ratio) ;
            }

            CArbCoord2D G, Gp, H, Hp ;

            G[0] = data[i].end_coords[0][0] +
                   data[i].end_radii[0] * cos(alpha+theta) ;
            G[1] = data[i].end_coords[0][1] +
                   data[i].end_radii[0] * sin(alpha+theta) ;
            Gp[0] = data[i].end_coords[0][0] +
                    data[i].end_radii[0] * cos(alpha-theta) ;
            Gp[1] = data[i].end_coords[0][1] +
                    data[i].end_radii[0] * sin(alpha-theta) ;
            H[0] = data[i].end_coords[1][0] +
                   data[i].end_radii[1] * cos(alpha+theta) ;
            H[1] = data[i].end_coords[1][1] +
                   data[i].end_radii[1] * sin(alpha+theta) ;
            Hp[0] = data[i].end_coords[1][0] +
                    data[i].end_radii[1] * cos(alpha-theta) ;
            Hp[1] = data[i].end_coords[1][1] +
                    data[i].end_radii[1] * sin(alpha-theta) ;

            // now determine the equations for the
            // lines G - Gp, Gp - Hp, Hp - H, H - G

            FindLineCoefs(G,Gp,data[i].line_coefs[0]) ;
            FindLineCoefs(Gp,Hp,data[i].line_coefs[1]) ;
            FindLineCoefs(Hp,H,data[i].line_coefs[2]) ;
            FindLineCoefs(H,G,data[i].line_coefs[3]) ;
        }
    }
    return(data) ;
}




// %(CArbRmshRegion2D::CharElemSize-double-|^const-int-const|-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    CharElemSize - compute a characteristic element size
**
**      double CharElemSize(
**              const int     elem,
**              CArbMshTopo2D *topo) const
**
**        elem - (in)  element id
**        topo - (in)  current mesh topology
**
**      Description: This function computes a characteristic size for
**          an element. Currently this is the average of the lengths of
**          the element's sides.
**
**      Return Value: a characteristic element size
**
**
** -- */

double CArbRmshRegion2D::CharElemSize(
           const int elem,
           CArbMshTopo2D *topo) const
{
    CArbTopoEdgeOnElemIterator iter(topo,elem) ;
    CArbCoord2D *cd0 = Owner->NodeTable->Fetch(iter[0]) ;
    CArbCoord2D *cd1 = Owner->NodeTable->Fetch(iter[1]) ;
    CArbCoord2D delta = *cd1 - *cd0 ;
    double sum = delta.Magnitude() ;
    int num = 1 ;

    ++iter ;
    while (iter.More()) {
        cd0 = cd1 ;
        cd1 = Owner->NodeTable->Fetch(iter[1]) ;
        delta = *cd1 - *cd0 ;
        sum += delta.Magnitude() ;
        ++num ;
        ++iter ;
    }
    return((Owner->Order == LINEAR) ? sum/num : 2.0*sum/num) ;
}




// %(CArbRmshRegion2D::DoDeleteElements-CArbMshTopo2D-|*-int-|-CArbMshCrackRegion2D::DeleteSegData-|*-MaterialHash-|*-CornerNodesSet-|*)
/* ++ ----------------------------------------------------------
**
**    DoDeleteElements - delete elements for remeshing
**
**      CArbMshTopo2D *DoDeleteElements(
**              int                                 num_segs,
**              CArbMshCrackRegion2D::DeleteSegData *del_segs,
**              MaterialHash                        *mat_table,
**              CornerNodesSet                      *corner_nodes)
**
**        num_segs     - (in)  number of delete segments
**        del_segs     - (in)  array of delete segments
**        mat_table    - (out) hash table containing the material id's
**                             for the remesh regions
**        corner_nodes - (out) set of the id's of the corner nodes
**
**      Description: This function deletes all elements that are
**          adjacent to any node that falls in one of the delete
**          segments.
**
**      Return Value: a MshTopo2D object containing the boundary(s) of
**          the deleted element region(s)
**
**
** -- */

CArbMshTopo2D *CArbRmshRegion2D::DoDeleteElements(
                            int num_segs,
                            DeleteSegData *del_segs,
                            MaterialHash *mat_table,
                            CornerNodesSet *corner_nodes,
                            bool /*closed*/)
{
    CArbHashTableIterator<int,CArbCoord2D> iter(Owner->NodeTable) ;
    CArbMshTopo2D *boundary = new CArbMshTopo2D() ;
    int ii ;

    if (Owner->DeletedElems != 0) delete Owner->DeletedElems ;
    Owner->DeletedElems = new CArbArray<int> ;

    // loop through all the nodes and see if the fall in
    // any of the delete segments.  Also if this is a closed
    // boundary count the number of boundary crossings so that
    // we can delete elements that are completely contained
    // inside the boundary

    bool DeleteAll = true ;

    for (iter.First() ; iter.More() ; ++iter) {

        bool do_delete = false ;
        int num_cross = 0 ;

        if (DeleteAll) {
            do_delete = true ;
        } else {

            CArbCoord2D ncoord = *(iter.Entry()) ;

            for (int j=0 ; j<num_segs ; j++) {

                // first check to see if we are close to
                // either of the end points.  Note that we
                // only check the first point for the first
                // segment because after that we will have
                // already check for the previous segment.

                if (j == 0) {
                    CArbCoord2D da = ncoord - del_segs[j].end_coords[0] ;
                    double la_sqr = da[0]*da[0] + da[1]*da[1] ;

                    if (la_sqr <= del_segs[j].end_radii_sqr[0]) {
                        do_delete = true ;
                        break ;
                    }
                }

                CArbCoord2D db = ncoord - del_segs[j].end_coords[1] ;
                double lb_sqr = db[0]*db[0] + db[1]*db[1] ;

                if (lb_sqr <= del_segs[j].end_radii_sqr[1]) {
                    do_delete = true ;
                    break ;
                }

                if (!del_segs[j].degenerate_case) {

                    // now check to see if the node is inside the
                    // region bounded by the four defining lines

                    if ((EvlLineEqn(ncoord,del_segs[j].line_coefs[0]) >= 0.0) &&
                        (EvlLineEqn(ncoord,del_segs[j].line_coefs[1]) >= 0.0) &&
                        (EvlLineEqn(ncoord,del_segs[j].line_coefs[2]) >= 0.0) &&
                        (EvlLineEqn(ncoord,del_segs[j].line_coefs[3]) >= 0.0)) {
                        do_delete = true ;
                        break ;
                    }
                }

                if (ScanCross(ncoord,del_segs[j].end_coords[0],
                              del_segs[j].end_coords[1])) ++num_cross ;
            }
        }

        if (do_delete || ((num_cross%2) == 1)) {

            // build a list of the elements we want to delete

            CArbArray<int> elems ;
            CArbTopoAdjVtxIterator aiter(Owner->MshTopo,iter.Key()) ;

            for (aiter.First() ; aiter.More() ; ++aiter) {
                if (aiter.CcwElem() != NO_ELEM) {
                    ArbMshElement2D *eptr =
                        Owner->ElemTable->Fetch(aiter.CcwElem()) ;
                    if ((Owner->TemplateElems->Fetch(eptr->elem_id) == 0) &&
                        !Owner->NoRmshMats->HasKey(eptr->mat_id))
                        elems.InsertAtEnd(aiter.CcwElem()) ;
                }
            }

            for (int k=0 ; k<elems.NumEntries() ; ++k) {
                ArbMshElement2D *eptr = Owner->ElemTable->Fetch(elems[k]) ;
                mat_table->Store(eptr->elem_id,eptr->mat_id) ;
                Owner->ElemTable->Remove(elems[k]) ;
                int nnode, nodes[8] ;
                Owner->MshTopo->GetElemNodes(elems[k],&nnode,nodes) ;
                if (Owner->Order == LINEAR)
                    for (ii=0 ; ii<nnode ; ++ii)
                        corner_nodes->Store(nodes[ii],1) ;
                else
                    for (ii=0 ; ii<nnode ; ii+=2)
                        corner_nodes->Store(nodes[ii],1) ;

                UpdateBoundary(boundary,elems[k],nnode,nodes,mat_table) ;

                Owner->DeletedElems->InsertAtEnd(elems[k]) ;
                Owner->MshTopo->DeleteElement(nnode,nodes) ;
            }
        }
    }

    boundary->RebuildElemList() ;

    return(boundary) ;
}




// %(CArbRmshRegion2D::FillCharElemSize-void-|^const-int-|-CArbMshCrackRegion2D::CrackPt-|*-double-|*)
/* ++ ----------------------------------------------------------
**
**    FillCharElemSize - determines crack-tip element sizes
**
**      void FillCharElemSize(
**              int                           num_pts,
**              CArbMshCrackRegion2D::CrackPt *pts,
**              double                        *insitu_size) const
**
**        num_pts     - (in)  number of points in the list
**        pts         - (i/o) crack point description list
**        insitu_size - (in)  list of existing local element sizes
**
**      Description: This function computes and fills in any
**          characteristic element sizes not provided in the crack
**          point list
**
**
** -- */

void CArbRmshRegion2D::FillCharElemSize(
         int num_pts,
         CArbMshCrackRegion2D::CrackPt *pts,
         double *insitu_size) const
{
    // find the characteristic size of all the elements in
    // which the crack points fall.  The first and last points
    // may be on or near the surface of the object.  However,
    // we cannot guarantee that a search will find the element
    // that contains them. Therefore we fist find the nearest
    // node, and then find the characteristic size of the elements
    // adjacent to that node.

    for (int j=0 ; j < num_pts ; ++j) {
        int near_node = FindNearestNode(
                CArbCoord2D(pts[j].coord[0],pts[j].coord[1])) ;
        insitu_size[j] = GetAdjCharSize(near_node) ;
    }

    // if we have no better information, set the
    // element sizes based on the insitu element sizes

    for (int k=0 ; k < num_pts ; ++k) {
        if (!pts[k].has_char_size) {
            pts[k].char_elem_size = insitu_size[k] ;
        }
    }
}


// %(CArbRmshRegion2D::FindNumRings-int-|^const-int-|-double-|*)
/* ++ ----------------------------------------------------------
**
**    FindNumRings - find the number of rings in a crack-tip template
**
**      int FindNumRings(
**              int    tip_id,
**              double *size) const
**
**        tip_id - (in)  crack-tip id
**        size   - (out) template size
**
**      Description: This function attempts to determine the number of
**          template rings that were inserted at this crack tip during
**          the last remeshing.
**
**      Return Value: number of rings
**
**
** -- */

#define DIST_TOL 1e-10

int CArbRmshRegion2D::FindNumRings(int tip_id,double *size) const
{
 /* What we want to do here is to try to determine the number
    of template rings that were inserted at this crack tip
    during the last remeshing.

    The procedure is as follows:

    1. For one of the crack-tip elements adjacent to the crack,
       move down the side not on the crack, and compute the
       length of the side.

    2. Now go to the element that would be the next ring and
       move along it's side.  First check to see if the far
       node is on a straight line from crack tip through the
       near node.  If so, find the length of the element side.

    3. Move to the next element and check to see if the far node
       is on a straight line and also check to see if we are
       getting the expected progression the node ring radius.
       The equation for this is r_m = (m/(m-2))*(r_n-1 - r_1).

    4. Continue this until the conditions are not satisfied. */

    int prev, num_rings = 1 ;

    /* we start at the given crack tip.  Since this may be
       an unconstrained tip and we may have not started at
       at one end crack-tip chain of nodes, move CCW until
       we find an edge with nonzero length */

    CArbTopoAdjVtxCyclicIterator iter(Owner->MshTopo,tip_id) ;

    while (iter.CcwElem() != NO_ELEM) --iter ;
    ++iter ;

    int tipv = tip_id ;
    int nearv = iter.AdjVtx() ;
    CArbCoord2D tipc = *(Owner->NodeTable->Fetch(tipv)) ;
    CArbCoord2D nearc = *(Owner->NodeTable->Fetch(nearv)) ;
    CArbCoord2D diff0 = nearc - tipc ;
    double r_tip = diff0.Magnitude() ;

    while (r_tip <= DIST_TOL*(fabs(tipc[0])+fabs(tipc[1]))) {
        tipv = nearv ;
        iter.NewVtx(tipv) ;
        while (iter.CcwElem() != NO_ELEM) --iter ;
        ++iter ;
        nearv = iter.AdjVtx() ;
        tipc = *(Owner->NodeTable->Fetch(tipv)) ;
        nearc = *(Owner->NodeTable->Fetch(nearv)) ;
        diff0 = nearc - tipc ;
        r_tip = diff0.Magnitude() ;
    }

    /* now if this an unconstrained tip we need to move CW
       one or two nodes to find the element edge we want
       to move along.  We check to see if this is unconstrained
       by looking to see if the distance to the CW node is zero */

    iter.NewVtx(tipv) ;
    while (iter.CcwElem() != NO_ELEM) ++iter ;
    nearv = iter.AdjVtx() ;
    nearc = *(Owner->NodeTable->Fetch(nearv)) ;
    diff0 = nearc - tipc ;
    r_tip = diff0.Magnitude() ;

    if (r_tip <= DIST_TOL*(fabs(tipc[0])+fabs(tipc[1]))) {
        tipv = nearv ;
        if (Owner->Order == QUADRATIC) {
            iter.NewVtx(tipv) ;
            while (iter.CcwElem() != NO_ELEM) ++iter ;
            tipv = iter.AdjVtx() ;
        }
    }

    /* now we have a tip node.  Find the vertex on the crack
       flank, and then go one more vertex to find an edge
       that we can move along to count the rings */

    iter.NewVtx(tipv) ;
    while (iter.CcwElem() == NO_ELEM) ++iter ;
    ++iter ;
    if (Owner->Order == QUADRATIC) {
        prev = iter.AdjVtx() ;
        iter.NewVtx(iter.AdjVtx()) ;
        while (iter.AdjVtx() != tipv) ++iter ;
        ++iter ;
    } else {
        prev = tipv ;
    }
    nearv = iter.AdjVtx() ;
    nearc = *(Owner->NodeTable->Fetch(nearv)) ;
    diff0 = nearc - tipc ;
    double r_0 = diff0.Magnitude() ;
    *size = r_0 ;

    /* now loop around the "next" vertex until we find
       the previous vertex.  Once we find this, we
       go two more adjacent vertices to find the far
       node on the next ring of elements */

    iter.NewVtx(nearv) ;
    while (iter.AdjVtx() != prev) ++iter ;
    ++iter ;  ++iter ;
    if (Owner->Order == QUADRATIC) {
        prev = iter.AdjVtx() ;
        iter.NewVtx(iter.AdjVtx()) ;
        while (iter.AdjVtx() != nearv) ++iter ;
        ++iter ;
    } else {
        prev = nearv ;
    }

    /* now check to see that the tip, near, and far nodes
       all fall on a straight line */

    int farv = iter.AdjVtx() ;
    CArbCoord2D farc = *(Owner->NodeTable->Fetch(farv)) ;
    CArbCoord2D diff1 = farc - nearc ;
    CArbCoord2D diff2 = farc - tipc ;
    double r_1 = diff1.Magnitude() ;
    double r_2 = diff2.Magnitude() ;

    double dot = diff0 * diff1 ;
    double cosang = dot / (r_0 * r_1) ;

    if ((cosang < 0.9999) || (cosang > 1.0001)) return(num_rings) ;

    ++num_rings ;
    *size = r_2 ;

    /* now loop for the rest of the rings */

    double r_expected ;
    double r_tol = 0.001*r_tip ;

    while (1) {
        iter.NewVtx(farv) ;
        while (iter.AdjVtx() != prev) ++iter ;
        ++iter ;  ++iter ;
        if (Owner->Order == QUADRATIC) {
            prev = iter.AdjVtx() ;
            iter.NewVtx(iter.AdjVtx()) ;
            while (iter.AdjVtx() != farv) ++iter ;
            ++iter ;
        } else {
            prev = nearv ;
        }
        nearv = farv ;
        nearc = farc ;

        farv = iter.AdjVtx() ;
        farc = *(Owner->NodeTable->Fetch(farv)) ;
        diff1 = nearc - tipc ;
        diff2 = farc - tipc ;
        double r_i = diff2.Magnitude() ;

        dot = diff1 * diff2 ;
        cosang = dot / (diff1.Magnitude() * diff2.Magnitude()) ;

        if ((cosang < 0.995) || (cosang > 1.005)) return(num_rings) ;

        r_expected = 0.5 * (num_rings+1) *
            (4*r_0 + (num_rings+1)*r_2 - 2*(num_rings+1)*r_0 - r_2) ;

        if (((r_i+r_tol) < r_expected) || ((r_i-r_tol) > r_expected))
            return(num_rings) ;

        ++num_rings ;
        *size = r_i ;
    }
}




// %(CArbRmshRegion2D::FillTipData-void-|^const-int-|-CArbMshCrackRegion2D::CrackTipData-|*)
/* ++ ----------------------------------------------------------
**
**    FillTipData - get a local copy of crack-tip data
**
**      void FillTipData(
**              int                                crack_tip_id,
**              CArbMshCrackRegion2D::CrackTipData *tip_data) const
**
**        crack_tip_id - (in)  crack-tip id
**        tip_data     - (out) local copy of crack-tip parameters
**
**      Description: This function fills a local copy of a structure
**          containing parameters for a crack tip.
**
**
** -- */

void CArbRmshRegion2D::FillTipData(
              int crack_tip_id,
              CArbMshCrackRegion2D::CrackTipData &tip_data) const
{
    tip_data = Owner->DefaultTipData ;
    tip_data.tip_id = crack_tip_id ;

    if ((crack_tip_id != DEFAULT_ID) && (Owner->TipData != 0)) {
        CArbMshCrackRegion2D::CrackTipData *data =
           Owner->TipData->Fetch(crack_tip_id) ;
        if (data == 0) return ;
        if (data->elem_shape_set)
            tip_data.tip_elem_shape = data->tip_elem_shape ;
        if (data->qrtr_pt_set)
            tip_data.qrtr_pt_flag = data->qrtr_pt_flag ;
        if (data->const_tip_set)
            tip_data.const_tip_flag = data->const_tip_flag ;
        if (data->num_tip_set)
            tip_data.num_tip_elems = data->num_tip_elems ;
        if (data->temp_radius_set)
            tip_data.template_radius = data->template_radius ;
        if (data->prog_ratio_set)
            tip_data.progression_ratio = data->progression_ratio ;
        if (data->num_rings_set)
            tip_data.number_of_rings = data->number_of_rings ;
    }
}




// %(CArbRmshRegion2D::FindContainingElem-int-|^const-CArbCoord2D-const|&-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    FindContainingElem - find the containg element
**
**      int FindContainingElem(
**              const CArbCoord2D &point,
**              CArbMshTopo2D     *topo) const
**
**        point - (in)  point to check
**        topo  - (in)  mesh topology to check
**
**      Description: Given an input point, this function determines
**          which element contains the point.
**
**      Return Value: the containing element id
**
**
** -- */

int CArbRmshRegion2D::FindContainingElem(
                                      const CArbCoord2D &point,
                                      CArbMshTopo2D *topo) const
{
    int num_elems = topo->NumElements() ;
    int *elems = topo->GetElemList() ;

    // for each element, loop through the boundary edges
    // and see if it is crossed by a scan line from minus
    // infinity to the point.  When all done count the number
    // of crossings.  If odd, this is the element so return it.

    for (int i=0 ; i<num_elems ; ++i) {
        CArbTopoEdgeOnElemIterator iter(topo,elems[i]) ;
        CArbCoord2D *cd0 = Owner->NodeTable->Fetch(iter[0]) ;
        CArbCoord2D *cd1 = Owner->NodeTable->Fetch(iter[1]) ;
        int num_cross = ScanCross(point,*cd0,*cd1) ? 1 : 0 ;
        ++iter ;

        while (iter.More()) {
            cd0 = cd1 ;
            cd1 = Owner->NodeTable->Fetch(iter[1]) ;
            if (ScanCross(point,*cd0,*cd1)) ++num_cross ;
            ++iter ;
        }

        if ((num_cross % 2) == 1) {
            int found = elems[i] ;
            delete [] elems ;
            return(found) ;
        }
    }
    delete [] elems ;
    return(NO_ELEM) ;
}




// %(CArbRmshRegion2D::FindNearestNode-int-|^const-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    FindNearestNode - find nearest node
**
**      int FindNearestNode(CArbCoord2D point) const
**
**        point - (in)  search point
**
**      Description: This function determines the node closest to a
**          given point.
**
**      Return Value: nearest node id
**
**
** -- */

int CArbRmshRegion2D::FindNearestNode(CArbCoord2D point) const
{
    CArbHashTableIterator<int,CArbCoord2D> iter(Owner->NodeTable) ;

    CArbCoord2D delta = point - *(iter.Entry()) ;
    double min_dist = delta.Magnitude() ;
    int closest = iter.Key() ;

    for (++iter ; iter.More() ; ++iter) {
        delta = point - *(iter.Entry()) ;
        double dist = delta.Magnitude() ;
        if (dist < min_dist) {
            min_dist = dist ;
            closest = iter.Key() ;
        }
    }
    return(closest) ;
}




// %(CArbRmshRegion2D::FindTipNodeList-int-|*^const-int-const|-int-|*-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    FindTipNodeList - find the id's of nodes at the crack-tip
**
**      int *FindTipNodeList(
**              const int     tip_node_id,
**              int           *num_nodes,
**              CArbMshTopo2D *topo) const
**
**        tip_node_id - (in)  crack-tip node id.
**        num_nodes   - (out) number of crack-tip nodes
**        topo        - (out) mesh topology
**
**      Description: This function returns an array of the crack-tip
**          node id's for the given crack tip. If this is a constrained
**          crack-tip node there will only be one element in the list.
**          There will be more than one node for unconstrained
**          crack-tip nodes.
**
**      Return Value: List of crack-tip nodes. Ownership of this memory
**          passes to the client, which should eventually call delete
**          [].
**
**
** -- */

int *CArbRmshRegion2D::FindTipNodeList(
                       const int tip_node_id,
                       int *num_nodes,
                       CArbMshTopo2D *topo) const
{
    int loc_tip = tip_node_id ;

    // get the crack tip coordinates

    CArbCoord2D *tipc = Owner->NodeTable->Fetch(loc_tip) ;
    int num = 1 ;

    // first loop CW from this node finding all nodes with
    // the same coordinates

    CArbTopoAdjVtxCyclicIterator citer(topo,loc_tip) ;
    while(1) {
        while (citer.CcwElem() != NO_ELEM) ++citer ;
        ++citer ;
        CArbCoord2D *next = Owner->NodeTable->Fetch(citer.AdjVtx()) ;
        if (*next == *tipc)
            loc_tip = citer.AdjVtx() ;
        else
            break ;
        citer.NewVtx(citer.AdjVtx()) ;
    }


    // Now loop CCW from this node finding all nodes with
    // the same coordinates

    CArbTopoAdjVtxIterator iter(topo,loc_tip) ;

    while(1) {
        while (iter.CcwElem() != NO_ELEM) ++iter ;
        CArbCoord2D *next = Owner->NodeTable->Fetch(iter.AdjVtx()) ;
        if (*next == *tipc)
            ++num ;
        else
            break ;
        iter.NewVtx(iter.AdjVtx()) ;
    }

    // allocate memory

    int *nodes = new int[num] ;
    nodes[0] = loc_tip ;
    int cur = 1 ;

    // Now loop CCW from this node finding all nodes with
    // the same coordinates

    iter.NewVtx(loc_tip) ;
    while(cur < num) {
        while (iter.CcwElem() != NO_ELEM) ++iter ;
        nodes[cur] = iter.AdjVtx() ;
        ++cur ;
        iter.NewVtx(iter.AdjVtx()) ;
    }
    *num_nodes = num ;
    return(nodes) ;
}




// %(CArbRmshRegion2D::GetAdjCharSize-double-|^const-int-|)
/* ++ ----------------------------------------------------------
**
**    GetAdjCharSize - find the mean adjacent element size
**
**      double GetAdjCharSize(int node) const
**
**        node - (in)  node id
**
**      Description: This function finds the mean characteristic size
**          of all the elements adjacent to a node.
**
**      Return Value: mean adjacent element size
**
**
** -- */

double CArbRmshRegion2D::GetAdjCharSize(int node) const
{
    double sum = 0.0 ;
    int num = 0 ;

    // loop for all the elements adjacent to a node

    CArbTopoAdjVtxIterator iter(Owner->MshTopo,node) ;

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.CcwElem() != NO_ELEM) {
            sum += CharElemSize(iter.CcwElem(),Owner->MshTopo) ;
            ++num ;
        }
    }
    return((num == 0) ? 0 : sum/num) ;
}




// %(CArbRmshRegion2D::GenerateCrackSegment-void-|-CArbMshCrackRegion2D::CrackPt-|*-CArbMshCrackRegion2D::CrackPt-|*-bool-|-double-|-double-|-bool-|-double-|-double-|-CArbArray-|<ArbMshNode>*)
/* ++ ----------------------------------------------------------
**
**    GenerateCrackSegment - determine nodal spacing on a crack segment
**
**      void GenerateCrackSegment(
**              CArbMshCrackRegion2D::CrackPt *start_pt,
**              CArbMshCrackRegion2D::CrackPt *end_pt,
**              bool                          start_crack_flag,
**              double                        template_length_s,
**              double                        template_projection_s,
**              bool                          end_crack_flag,
**              double                        template_length_e,
**              double                        template_projection_e,
**              CArbArray<ArbMshNode>*        nodes)
**
**        start_pt              - (in)  starting point
**        end_pt                - (in)  ending point
**        start_crack_flag      - (in)  true if the start point is a
**                                      crack tip
**        template_length_s     - (in)  start point crack-template size
**        template_projection_s - (in)  size of the first element
**                                      outside of the template based
**                                      on the template grading
**                                      projection.
**        end_crack_flag        - (in)  true if the end point is a
**                                      crack tip
**        template_length_e     - (in)  end point crack-template size
**        template_projection_e - (in)  size of the first element
**                                      outside of the template based
**                                      on the template grading
**                                      projection.
**        nodes                 - (out) generated nodes
**
**      Description: This function determines node spacings and
**          locations along a crack segment. The nodes are positiond so
**          that they are graded in size along the segment.
**
**
** -- */

void CArbRmshRegion2D::GenerateCrackSegment(
            CArbMshCrackRegion2D::CrackPt *start_pt,
            CArbMshCrackRegion2D::CrackPt *end_pt,
            bool start_crack_flag,
            double template_length_s,
            double template_projection_s,
            bool end_crack_flag,
            double template_length_e,
            double template_projection_e,
            CArbArray<ArbMshNode> *nodes)
{
    /* Here's a map of the crack segment cases

         start                                        end

    no tip o-------+---------------------------+-------o no tip
              ^                                     ^
              |                                     |
    char_len -+                                     +- char_len

                  \
    tip    o-------*---+-----------------------+-------o no tip
              ^   /  ^                              ^
              |      |                              |
              |      +- template_projection_s       +- char_len
              +- template_length_s

                                                /
    no tip o-------+-----------------------+---*-------o tip
              ^                              ^  \   ^
              |                              |      |
    char_len -+       template_projection_e -+      |
                                 template_length_e -+

                  \                             /
    tip    o-------*---+-------------------+---*-------o tip
              ^   / ^                        ^  \   ^
              |     |                        |      |
              |     +- template_projection_s |      |
              +- template_length_s           |      |
                      template_projection_e -+      |
                                 template_length_e -+
    */

    int k ;

    // first make points for the start and end of the segment
    // and find the normailized vector pointing frome the start
    // to the end.

    CArbCoord2D pt0(start_pt->coord[0],start_pt->coord[1]) ;
    CArbCoord2D pt1(end_pt->coord[0],end_pt->coord[1]) ;
    CArbCoord2D dir = pt1 - pt0 ;
    CArbCoord2D norm = dir.Normalize() ;

    // now adjust the positions of the end if there are crack tips

    double start_len, end_len ;

    if (start_crack_flag) {
        pt0 += template_length_s * norm ;
        start_len = template_projection_s ;
    } else {
        start_len = start_pt->char_elem_size ;
    }

    if (end_crack_flag) {
        pt1 -= template_length_e * norm ;
        end_len = template_projection_e ;
    } else {
        end_len = end_pt->char_elem_size ;
    }

    CArbCoord2D delta = pt1 - pt0 ;
    double len = delta.Magnitude() ;
    double rat = end_len / start_len ;
    int num = unsigned((len * rat) /
           (end_len*(1.0+0.5*(rat-1.0)))) ;

    if (num < 1) num = 1 ;

    double fact = len / (num + 0.5*num*(rat-1.0)) ;
    CArbCoord2D pt_a, pt_b ;
    double dist ;

    pt_a = pt0 ;
    pt_b = (num == 1) ? pt1 : pt0 + norm * fact ;
    if (start_crack_flag) {
    } else {
        nodes->InsertAtEnd(Owner->NewNode(pt_a.x(),pt_a.y())) ;
    }

    for (k=0 ; k<num-1 ; ++k) {
        if (Owner->Order == QUADRATIC) {
            CArbCoord2D pt_m = 0.5 * (pt_a + pt_b) ;
            nodes->InsertAtEnd(Owner->NewNode(pt_m.x(),pt_m.y())) ;
        }
        nodes->InsertAtEnd(Owner->NewNode(pt_b.x(),pt_b.y())) ;
        pt_a = pt_b ;
        dist = fact*((k+2) + 0.5*(k+1)*(k+2)*(rat-1.0)/(num-1)) ;
        pt_b = pt0 + norm * dist ;
    }

    if (Owner->Order == QUADRATIC) {
        CArbCoord2D pt_m = 0.5 * (pt_a + pt_b) ;
        nodes->InsertAtEnd(Owner->NewNode(pt_m.x(),pt_m.y())) ;
    }
}




// %(CArbRmshRegion2D::GenerateCrackTipBdry-int-|*-CArbMshCrackRegion2D::CrackTipData-|*-int-|*)
/* ++ ----------------------------------------------------------
**
**    GenerateCrackTipBdry - find the crack-tip boundary
**
**      int *GenerateCrackTipBdry(
**              CArbMshCrackRegion2D::CrackTipData *tip_data,
**              int                                *num_bdry)
**
**        tip_data - (in)  crack-tip parameters
**        num_bdry - (out) number of boundary nodes
**
**      Description: This function finds a portion of the boundary of
**          the inside of a crack in the crack-tip template.
**
**      Return Value: A list of the crack-tip boundary node id's.
**          Ownership of this memory passes to the client, which should
**          eventually call delete [].
**
**
** -- */

void CArbRmshRegion2D::GenerateCrackTipBdry(int num_tip_elems,
              CArbMshCrackRegion2D::CrackTipData &tip_data,
              CArbArray<int> &bdry)
{
    int tip_type = ClassifyTipType(Owner->Order,
                                   tip_data.tip_elem_shape,
                                   tip_data.tip_elem_type,
                                   tip_data.const_tip_flag) ;
    int i ;

    switch (tip_type) {
        case LIN_TRI_COL_TRI:
        case LIN_TRI_COL_QUAD:
        case LIN_TRI_UNCOL:         // do the multiple tip nodes later
            for (i=tip_data.number_of_rings-1 ; i>=0 ; --i)
                bdry.InsertAtEnd(i * (num_tip_elems+1) + 1) ;
            bdry.InsertAtEnd(0) ;
            for (i=1 ; i <= int(tip_data.number_of_rings) ; ++i)
                bdry.InsertAtEnd(i * (num_tip_elems+1)) ;
            break ;

        case LIN_QUAD:
            for (i=tip_data.number_of_rings ; i>0 ; --i)
                bdry.InsertAtEnd(4 * i * (i-1) + i) ;
            bdry.InsertAtEnd(0) ;
            for (i=1 ; i <= int(tip_data.number_of_rings) ; ++i)
                bdry.InsertAtEnd(4 * i * (i+1) + i) ;
            break ;

        case QDC_TRI_COL_TRI:
        case QDC_TRI_COL_QUAD:
        case QDC_TRI_UNCOL:     // do the multiple crack tips later
            for (i=tip_data.number_of_rings-1 ; i>=0 ; --i) {
                bdry.InsertAtEnd(i * (3*num_tip_elems+2) + 1) ;
                bdry.InsertAtEnd(i * (3*num_tip_elems+2) +
                            (num_tip_elems+2)) ;
            }
            bdry.InsertAtEnd(0) ;
            for (i=0 ; i < int(tip_data.number_of_rings) ; ++i) {
                bdry.InsertAtEnd(i * (3*num_tip_elems+2) +
                                 2*(num_tip_elems+1)) ;
                bdry.InsertAtEnd(i * (3*num_tip_elems+2) +
                                 (num_tip_elems+1)) ;
            }
            break ;

        case QDC_QUAD:
            for (i=tip_data.number_of_rings-1 ; i>=0 ; --i) {
                bdry.InsertAtEnd(12*i*(i+1) - 2*i + 1) ;
                bdry.InsertAtEnd(12*i*(i+1) + 6*i + 10) ;
            }
            bdry.InsertAtEnd(0) ;
            for (i=0 ; i < int(tip_data.number_of_rings) ; ++i) {
                bdry.InsertAtEnd(12*i*(i+1) + 14*i + 14) ;
                bdry.InsertAtEnd(12*i*(i+1) + 6*i + 9) ;
            }
            break ;
    }
}




// %(CArbRmshRegion2D::GenerateCrackTipData-void-|-CArbCoord2D-|-CArbCoord2D-|-double-|-CArbMshCrackRegion2D::CrackTipData-|*-int-|-CArbMshTopo2D-|*-CArbMshRegion2D-|*-int-|*-ArbMshNode-|**-int-|*-int-|**-int-|*)
/* ++ ----------------------------------------------------------
**
**    GenerateCrackTipData - generate crack-tip template nodes and
**                           elements
**
**      void GenerateCrackTipData(
**              CArbCoord2D                        tip,
**              CArbCoord2D                        norm,
**              double                             tip_angle,
**              CArbMshCrackRegion2D::CrackTipData *tip_data,
**              int                                mat_id,
**              CArbMshTopo2D                      *crack_bdry,
**              CArbMshRegion2D                    *region,
**              int                                *num_temp_nodes,
**              ArbMshNode                         **temp_nodes,
**              int                                *num_temp_bdry,
**              int                                **temp_bdry,
**              int                                *tip_node_id)
**
**        tip            - (in)  crack-tip coordinates
**        norm           - (in)  normalized crack direction
**        tip_angle      - (in)  included angle of the crack-tip
**        tip_data       - (in)  crack-tip parameters for this tip
**        mat_id         - (in)  material to assign to new elements
**        crack_bdry     - (i/o) description of the crack boundary
**                               topology
**        region         - (out) region to which nodes and elements are
**                               added
**        num_temp_nodes - (out) number of template nodes
**        temp_nodes     - (out) generated nodes
**        num_temp_bdry  - (out) number of template boundary nodes
**        temp_bdry      - (out) template boundary nodes
**        tip_node_id    - (out) generated crack tip node id
**
**      Description: This function generates the nodes and elements for
**          a crack-tip template.
**
**
** -- */

void CArbRmshRegion2D::GenerateCrackTipData(
                             CArbCoord2D tip, CArbCoord2D norm,
                             double tip_angle,
                             CArbMshCrackRegion2D::CrackTipData &tip_data,
                             int mat_id,
                             CArbMshTopo2D *crack_bdry,
                             CArbMshRegion2D *region,
                             CArbArray<ArbMshNode> &temp_nodes,
                             CArbArray<int> &temp_bdry,
                             int *tip_id)
{
    int i ;
    GenerateCrackTipNodes(tip,norm,tip_data.num_tip_elems,
                          tip_data,tip_angle,temp_nodes) ;

    CArbArray<ArbMshElement2D> temp_elems ;
    GenerateCrackTipElems(tip_data.num_tip_elems,
                          tip_data,temp_nodes,mat_id,temp_elems) ;

    GenerateCrackTipBdry(tip_data.num_tip_elems,tip_data,temp_bdry) ;

    // insert the tip elements into the mesh, and into
    // the current boundary

    for (i=0 ; i<temp_nodes.NumEntries() ; ++i) {
        Owner->AddNode(temp_nodes[i]) ;
        region->AddNode(temp_nodes[i]) ;
    }
    for (i=0 ; i<temp_elems.NumEntries() ; ++i) {
        Owner->AddElem(temp_elems[i]) ;
        int lnum = (temp_elems[i].num_nodes == 3) ||
                   (temp_elems[i].num_nodes == 6) ? 3 : 4 ;
        for (int i=0 ; i<lnum ; ++i)
            CornerNodes->Store(temp_elems[i].nodes[i],1) ;
    }
    if (Owner->Order == LINEAR) {
        for (int j=0 ; j<temp_elems.NumEntries() ; ++j) {
            crack_bdry->InsertCollapsedElement(
                                       temp_elems[j].elem_id,
                                       temp_elems[j].num_nodes,
                                       temp_elems[j].nodes) ;
        }
    } else {
        int nodes[8], jj ;
        for (int j=0 ; j<temp_elems.NumEntries() ; ++j) {
            int cur = 0 ;
            for (jj=0 ; jj<temp_elems[j].num_nodes/2 ; ++jj) {
                nodes[jj*2] = temp_elems[j].nodes[jj] ;
                ++cur ;
            }
            for (jj=0 ; jj<temp_elems[j].num_nodes/2 ; ++jj) {
                nodes[jj*2+1] = temp_elems[j].nodes[cur] ;
                ++cur ;
            }
            crack_bdry->InsertCollapsedElement(
                                       temp_elems[j].elem_id,
                                       temp_elems[j].num_nodes,
                                       nodes) ;
        }
    }
    *tip_id = temp_nodes[0].id ;
}




// %(CArbRmshRegion2D::GenerateCrackTipElems-ArbMshElement2D-|*-CArbMshCrackRegion2D::CrackTipData-|*-ArbMshNode-|*-int-|*-int-|)
/* ++ ----------------------------------------------------------
**
**    GenerateCrackTipElems - generate crack-template elements
**
**      ArbMshElement2D *GenerateCrackTipElems(
**              CArbMshCrackRegion2D::CrackTipData *tip_data,
**              ArbMshNode                         *nodes,
**              int                                *num_elems,
**              int                                mat_id)
**
**        tip_data  - (in)  crack-tip parameters
**        nodes     - (in)  crack-tip template nodes
**        num_elems - (out) number of elements generated
**        mat_id    - (in)  material id assigned to new elements
**
**      Description: This function generates all the elements for a
**          crack-tip template.
**
**      Return Value: A list of all descriptions for the generated
**          elements. Ownership of this list passes to the client,
**          which should eventually call delete [].
**
**
** -- */

void CArbRmshRegion2D::GenerateCrackTipElems(int num_tip_elems,
                       CArbMshCrackRegion2D::CrackTipData &tip_data,
                       CArbArray<ArbMshNode> &nodes,
                       int mat_id,
                       CArbArray<ArbMshElement2D> &elem)
{
    ArbMshElement2D lelem ;
    int tip_type = ClassifyTipType(Owner->Order,
                                   tip_data.tip_elem_shape,
                                   tip_data.tip_elem_type,
                                   tip_data.const_tip_flag) ;

    if (tip_data.number_of_rings == 1) tip_data.progression_ratio = 1.0 ;

    if (tip_data.tip_elem_shape == CArbMshCrackRegion2D::S_TRIANGLE) {

        // create the crack tip elements

        int outer_start = 0 ;
        int i, j, k, l, m, n, o ;
        switch (tip_type) {
            case LIN_TRI_COL_TRI:
            case LIN_TRI_UNCOL:      // multipl tip nodes later
                for (j=0 ; j<num_tip_elems ; ++j) {
                    FillElem(Owner->NewElemId(),mat_id,3,nodes[0].id,
                             nodes[j+1].id,nodes[j+2].id,0,&lelem) ;
                    elem.InsertAtEnd(lelem) ;
                    Owner->TemplateElems->Store(lelem.elem_id,
                                                tip_data.tip_id) ;
                }
                outer_start = 1 ;
                break ;

            case LIN_TRI_COL_QUAD:
                for (j=0 ; j<num_tip_elems ; ++j) {
                    FillElem(Owner->NewElemId(),mat_id,4,nodes[0].id,
                             nodes[j+1].id,nodes[j+2].id,
                             nodes[0].id,&lelem) ;
                    elem.InsertAtEnd(lelem) ;
                    Owner->TemplateElems->Store(lelem.elem_id,
                                                tip_data.tip_id) ;
                }
                outer_start = 1 ;
                break ;

            case QDC_TRI_COL_TRI:
            case QDC_TRI_UNCOL:    // multiple tip-nodes later
                k = num_tip_elems + 2 ;
                l = k + num_tip_elems + 1 ;
                for (j=0 ; j<num_tip_elems ; ++j) {
                    FillQuadElem(Owner->NewElemId(),mat_id,6,nodes[0].id,
                                 nodes[j+1].id,nodes[j+2].id,
                                 nodes[j+k].id,nodes[j+l].id,
                                 nodes[j+k+1].id,0,0,&lelem) ;
                    elem.InsertAtEnd(lelem) ;
                    Owner->TemplateElems->Store(lelem.elem_id,
                                                tip_data.tip_id) ;
                }
                outer_start = 1 ;
                break ;

            case QDC_TRI_COL_QUAD:
                k = num_tip_elems + 2 ;
                l = k + num_tip_elems + 1 ;
                for (j=0 ; j<num_tip_elems ; ++j) {
                    FillQuadElem(Owner->NewElemId(),mat_id,8,nodes[0].id,
                                 nodes[j+1].id,nodes[j+2].id,
                                 nodes[0].id,nodes[j+k].id,
                                 nodes[j+l].id,nodes[j+k+1].id,
                                 nodes[0].id,&lelem) ;
                    elem.InsertAtEnd(lelem) ;
                    Owner->TemplateElems->Store(lelem.elem_id,
                                                tip_data.tip_id) ;
                }
                outer_start = 1 ;
                break ;
        }

        // now the outer rings of elements

        if (Owner->Order == LINEAR) {
            for (i=1 ; i<tip_data.number_of_rings ; ++i) {
                k = outer_start ;
                l = outer_start + num_tip_elems + 1 ;
                for (j=0 ; j<num_tip_elems ; ++j) {
                    FillElem(Owner->NewElemId(),mat_id,4,
                             nodes[j+k].id,nodes[j+l].id,
                             nodes[j+l+1].id,nodes[j+k+1].id,
                             &lelem) ;
                    elem.InsertAtEnd(lelem) ;
                    Owner->TemplateElems->Store(lelem.elem_id,
                                                tip_data.tip_id) ;
                }
                outer_start += num_tip_elems+1 ;
            }
        } else {
            for (i=1 ; i<tip_data.number_of_rings ; ++i) {
                k = outer_start ;
                if (i == 1)
                    if (outer_start == 1)
                        m = k + 2*(num_tip_elems+1) ;
                    else
                        m = outer_start + 2 +
                            3*num_tip_elems ;
                else
                    m = k + 2*(num_tip_elems + 1) ;
                l = m + num_tip_elems ;
                n = m + 2*num_tip_elems + 1 ;
                o = n + num_tip_elems + 1 ;

                for (j=0 ; j<num_tip_elems ; ++j) {
                    FillQuadElem(Owner->NewElemId(),mat_id,8,
                                 nodes[j+k].id,nodes[j+l].id,
                                 nodes[j+l+1].id,nodes[j+k+1].id,
                                 nodes[j+n].id,nodes[j+o].id,
                                 nodes[j+n+1].id,nodes[j+m].id,
                                 &lelem) ;
                    elem.InsertAtEnd(lelem) ;
                    Owner->TemplateElems->Store(lelem.elem_id,
                                                tip_data.tip_id) ;
                }
                outer_start = l ;
            }
        }
    } else {

        if (Owner->Order == LINEAR) {
            int r, i ;
            for (r=0 ; r<tip_data.number_of_rings ; ++r) {
                int icur = 4*r*r - 3*r ;
                int ocur = 4*r*r + 5*r + 1 ;

                for (i=0 ; i<r ; ++i) {
                    DoLinQuad(true,Owner->NewElemId(),
                              mat_id,nodes,elem,&ocur,&icur) ;
                    Owner->TemplateElems->Store(elem[-1].elem_id,
                                                tip_data.tip_id) ;
                }
                DoLinQuad(false,Owner->NewElemId(),mat_id,nodes,elem,
                          &ocur,&icur) ;
                Owner->TemplateElems->Store(elem[-1].elem_id,
                                            tip_data.tip_id) ;
                for (i=0 ; i<2*r ; ++i) {
                    DoLinQuad(true,Owner->NewElemId(),
                              mat_id,nodes,elem,&ocur,&icur) ;
                    Owner->TemplateElems->Store(elem[-1].elem_id,
                                                tip_data.tip_id) ;
                }
                DoLinQuad(false,Owner->NewElemId(),mat_id,nodes,elem,
                          &ocur,&icur) ;
                Owner->TemplateElems->Store(elem[-1].elem_id,
                                            tip_data.tip_id) ;
                for (i=0 ; i<2*r ; ++i) {
                    DoLinQuad(true,Owner->NewElemId(),
                              mat_id,nodes,elem,&ocur,&icur) ;
                    Owner->TemplateElems->Store(elem[-1].elem_id,
                                                tip_data.tip_id) ;
                }
                DoLinQuad(false,Owner->NewElemId(),mat_id,nodes,elem,
                          &ocur,&icur) ;
                Owner->TemplateElems->Store(elem[-1].elem_id,
                                            tip_data.tip_id) ;
                for (i=0 ; i<2*r ; ++i) {
                    DoLinQuad(true,Owner->NewElemId(),
                              mat_id,nodes,elem,&ocur,&icur) ;
                    Owner->TemplateElems->Store(elem[-1].elem_id,
                                                tip_data.tip_id) ;
                }
                DoLinQuad(false,Owner->NewElemId(),mat_id,nodes,elem,
                          &ocur,&icur) ;
                Owner->TemplateElems->Store(elem[-1].elem_id,
                                            tip_data.tip_id) ;
                for (i=0 ; i<r ; ++i) {
                    DoLinQuad(true,Owner->NewElemId(),
                              mat_id,nodes,elem,&ocur,&icur) ;
                    Owner->TemplateElems->Store(elem[-1].elem_id,
                                                tip_data.tip_id) ;
                }
            }
        } else {
            int r, i ;
            int ocur, icur, s1cur, s2cur, s3cur ;
            for (r=0 ; r<tip_data.number_of_rings ; ++r) {
                icur = (r == 0) ? 0 : 12*r*r - 14*r + 3 ;
                ocur = 12*r*r + 10*r + 1 ;
                s1cur = 12*r*r + 18*r + 10 ;
                s2cur = 12*r*r + 26*r + 15 ;
                s3cur = (r == 0) ? 11 : 12*r*r + 2*r + 1 ;

                for (i=0 ; i<r ; ++i) {
                    DoQdcQuad(true,Owner->NewElemId(),
                              mat_id,nodes,elem,&ocur,&icur,
                              &s1cur,&s2cur,&s3cur) ;
                    Owner->TemplateElems->Store(elem[-1].elem_id,
                                                tip_data.tip_id) ;
                }
                DoQdcQuad(false,Owner->NewElemId(),mat_id,nodes,elem,
                          &ocur,&icur,&s1cur,&s2cur,&s3cur) ;
                Owner->TemplateElems->Store(elem[-1].elem_id,
                                            tip_data.tip_id) ;
                for (i=0 ; i<2*r ; ++i) {
                    DoQdcQuad(true,Owner->NewElemId(),
                              mat_id,nodes,elem,&ocur,&icur,
                              &s1cur,&s2cur,&s3cur) ;
                    Owner->TemplateElems->Store(elem[-1].elem_id,
                                                tip_data.tip_id) ;
                }
                DoQdcQuad(false,Owner->NewElemId(),mat_id,nodes,elem,
                          &ocur,&icur,&s1cur,&s2cur,&s3cur) ;
                Owner->TemplateElems->Store(elem[-1].elem_id,
                                            tip_data.tip_id) ;
                for (i=0 ; i<2*r ; ++i) {
                    DoQdcQuad(true,Owner->NewElemId(),
                              mat_id,nodes,elem,&ocur,&icur,
                              &s1cur,&s2cur,&s3cur) ;
                    Owner->TemplateElems->Store(elem[-1].elem_id,
                                                tip_data.tip_id) ;
                }
                DoQdcQuad(false,Owner->NewElemId(),mat_id,nodes,elem,
                          &ocur,&icur,&s1cur,&s2cur,&s3cur) ;
                Owner->TemplateElems->Store(elem[-1].elem_id,
                                            tip_data.tip_id) ;
                for (i=0 ; i<2*r ; ++i) {
                    DoQdcQuad(true,Owner->NewElemId(),
                              mat_id,nodes,elem,&ocur,&icur,
                              &s1cur,&s2cur,&s3cur) ;
                    Owner->TemplateElems->Store(elem[-1].elem_id,
                                                tip_data.tip_id) ;
                }
                DoQdcQuad(false,Owner->NewElemId(),mat_id,nodes,elem,
                          &ocur,&icur,&s1cur,&s2cur,&s3cur) ;
                Owner->TemplateElems->Store(elem[-1].elem_id,
                                            tip_data.tip_id) ;
                for (i=0 ; i<r ; ++i) {
                    DoQdcQuad(true,Owner->NewElemId(),
                              mat_id,nodes,elem,&ocur,&icur,
                              &s1cur,&s2cur,&s3cur) ;
                    Owner->TemplateElems->Store(elem[-1].elem_id,
                                                tip_data.tip_id) ;
                }
            }
        }
    }
}




// %(CArbRmshRegion2D::GenerateCrackTipNodes-ArbMshNode-|*-CArbCoord2D-|-CArbCoord2D-|-CArbMshCrackRegion2D::CrackTipData-|*-int-|*-double-|)
/* ++ ----------------------------------------------------------
**
**    GenerateCrackTipNodes - generate crack-template nodes
**
**      ArbMshNode *GenerateCrackTipNodes(
**              CArbCoord2D                        tip,
**              CArbCoord2D                        norm,
**              CArbMshCrackRegion2D::CrackTipData *tip_data,
**              int                                *num_nodes,
**              double                             tip_angle)
**
**        tip       - (in)  crack-tip coordinates
**        norm      - (in)  normalized crack direction
**        tip_data  - (in)  crack-tip parameters
**        num_nodes - (out) number of generated nodes
**        tip_angle - (in)  included angle at the crack tip
**
**      Description: This function generates all the nodes for a
**          crack-tip template.
**
**      Return Value: A list of the generated template nodes. Ownership
**          of this list passes to the client, which should eventually
**          call delete [].
**
**
** -- */

void CArbRmshRegion2D::GenerateCrackTipNodes(
                  CArbCoord2D tip,
                  CArbCoord2D norm,
                  int num_tip_elems,
                  CArbMshCrackRegion2D::CrackTipData &tip_data,
                  double tip_angle,
                  CArbArray<ArbMshNode> &nodes)
{
    int i, j ;
    double tip_size = tip_data.template_radius ;

    CArbCoord2D rnorm = -norm ;
    double theta = atan2(rnorm[1],rnorm[0]) ;

    if (tip_data.tip_elem_shape == CArbMshCrackRegion2D::S_TRIANGLE) {

        // we now replicate the crack tip nodes elsewhere

        nodes.InsertAtEnd(Owner->NewNode(tip[0],tip[1])) ;

        // create the rest of the nodes for the crack-tip element

        double b = 1.0 / tip_data.progression_ratio ;
        double fact = tip_size / (tip_data.number_of_rings +
                      0.5*tip_data.number_of_rings*(b-1.0)) ;
        double radius = (tip_data.number_of_rings != 1) ? fact : tip_size ;
        double elem_angle = (2*PI - tip_angle) / num_tip_elems ;

        for (j=0 ; j<=num_tip_elems ; ++j) {
            double angle = j*elem_angle + tip_angle/2.0 ;
            CArbCoord2D cnorm(cos(theta+angle),sin(theta+angle)) ;
            CArbCoord2D pt = tip + cnorm * radius ;
            nodes.InsertAtEnd(Owner->NewNode(pt[0],pt[1])) ;
        }

        if (Owner->Order == QUADRATIC) {

            // side nodes, possibly quarter-points

            for (j=0 ; j<=num_tip_elems ; ++j) {
                double angle = j*elem_angle + tip_angle/2.0 ;
                CArbCoord2D cnorm(cos(theta+angle),sin(theta+angle)) ;
                CArbCoord2D pt ;
                if (tip_data.qrtr_pt_flag)
                    pt = tip + cnorm * radius / 4.0 ;
                else
                    pt = tip + cnorm * radius / 2.0 ;
                nodes.InsertAtEnd(Owner->NewNode(pt[0],pt[1])) ;
            }

            // now side nodes opposite the crack tip

            for (j=0 ; j<num_tip_elems ; ++j) {
                double angle0 = j*elem_angle + tip_angle/2.0 ;
                double angle1 = (j+1)*elem_angle + tip_angle/2.0 ;
                CArbCoord2D cnorm0(cos(theta+angle0),sin(theta+angle0)) ;
                CArbCoord2D cnorm1(cos(theta+angle1),sin(theta+angle1)) ;
                CArbCoord2D pt0 = tip + cnorm0 * radius ;
                CArbCoord2D pt1 = tip + cnorm1 * radius ;
                nodes.InsertAtEnd(Owner->NewNode(0.5*(pt0[0]+pt1[0]),
                                            0.5*(pt0[1]+pt1[1]))) ;
            }
        }

        // now create the nodes for the rest of the crack-tip region

        for (i=2 ; i<=tip_data.number_of_rings ; ++i) {
            double radius = (tip_data.number_of_rings != 1) ?
                            (fact*(i + 0.5*i*(i-1.0)*(b-1.0) /
                             (tip_data.number_of_rings-1.0))) :
                            tip_size ;

            double elem_angle = (2*PI - tip_angle) / num_tip_elems ;
            for (j=0 ; j<=num_tip_elems ; ++j) {
                double angle = j*elem_angle + tip_angle/2.0 ;
                CArbCoord2D cnorm(cos(theta+angle),sin(theta+angle)) ;
                CArbCoord2D pt = tip + cnorm * radius ;
                nodes.InsertAtEnd(Owner->NewNode(pt[0],pt[1])) ;
            }

            if (Owner->Order == QUADRATIC) {
                double prev_radius = (tip_data.number_of_rings != 1) ?
                            (fact*((i-1) + 0.5*(i-1)*(i-2.0)*(b-1.0) /
                             (tip_data.number_of_rings-1.0))) : 0 ;
                double side_rad = 0.5*(radius + prev_radius) ;
                for (j=0 ; j<=num_tip_elems ; ++j) {
                    double angle = j*elem_angle + tip_angle/2.0 ;
                    CArbCoord2D cnorm(cos(theta+angle),sin(theta+angle)) ;
                    CArbCoord2D pt = tip + cnorm * side_rad ;
                    nodes.InsertAtEnd(Owner->NewNode(pt[0],pt[1])) ;
                }
                for (j=0 ; j<num_tip_elems ; ++j) {
                    double angle0 = j*elem_angle + tip_angle/2.0 ;
                    double angle1 = (j+1)*elem_angle + tip_angle/2.0 ;
                    CArbCoord2D cnorm0(cos(theta+angle0),sin(theta+angle0)) ;
                    CArbCoord2D pt0 = tip + cnorm0 * radius ;
                    CArbCoord2D cnorm1(cos(theta+angle1),sin(theta+angle1)) ;
                    CArbCoord2D pt1 = tip + cnorm1 * radius ;
                    nodes.InsertAtEnd(Owner->NewNode(0.5*(pt0[0]+pt1[0]),
                                                0.5*(pt0[1]+pt1[1]))) ;
                }
            }
        }
    } else {

        nodes.InsertAtEnd(Owner->NewNode(tip[0],tip[1])) ;

        double b = 1.0 / tip_data.progression_ratio ;
        double fact = tip_size / (tip_data.number_of_rings +
                      0.5*tip_data.number_of_rings*(b-1.0)) ;
        int r, i ;
        CArbCoord2D perp(-norm[1],norm[0]) ;
        for (r=1 ; r<=tip_data.number_of_rings ; ++r) {
            double radius = (tip_data.number_of_rings != 1) ?
                            (fact*(r + 0.5*r*(r-1.0)*(b-1.0) /
                             (tip_data.number_of_rings-1.0))) :
                            tip_size ;
            double leng = radius / r ;
            double fleng = (radius+radius*sin(tip_angle/2.0)) / r ;

            CArbCoord2D flank(-norm[0]*cos(tip_angle/2.0) +
                               norm[1]*sin(tip_angle/2.0),
                              -norm[0]*sin(tip_angle/2.0) +
                              -norm[1]*cos(tip_angle/2.0)) ;
            CArbCoord2D pt = tip + flank * radius ;
            nodes.InsertAtEnd(Owner->NewNode(pt[0],pt[1])) ;

            for (i=0 ; i<r ; ++i)   pt = PtUpdate(pt,-perp,leng,nodes) ;
            for (i=0 ; i<r*2 ; ++i) pt = PtUpdate(pt,norm,leng,nodes) ;
            for (i=0 ; i<r*2 ; ++i) pt = PtUpdate(pt,perp,fleng,nodes) ;
            for (i=0 ; i<r*2 ; ++i) pt = PtUpdate(pt,-norm,leng,nodes) ;
            for (i=0 ; i<r ; ++i)   pt = PtUpdate(pt,-perp,leng,nodes) ;

            if (Owner->Order == QUADRATIC) {
                GenTipNodesHelp(r,r,29,-38,12,3,-14,12,
                                tip_data.qrtr_pt_flag,nodes) ;
                GenTipNodesHelp(r,(r*2)-1,28,-37,12,4,-13,12,
                                tip_data.qrtr_pt_flag,nodes) ;
                GenTipNodesHelp(r,(r*2)-1,26,-35,12,4,-11,12,
                                tip_data.qrtr_pt_flag,nodes) ;
                GenTipNodesHelp(r,(r*2)-1,24,-33,12,4,-9,12,
                                tip_data.qrtr_pt_flag,nodes) ;
                GenTipNodesHelp(r,r,22,-31,12,4,-7,12,
                                tip_data.qrtr_pt_flag,nodes) ;
                j = 12*r*r -14*r + 3 ;
                for (i=j ; i<(j + r*8) ; ++i) {
                    pt[0] = 0.5*(nodes[i].coord[0]+nodes[i+1].coord[0]) ;
                    pt[1] = 0.5*(nodes[i].coord[1]+nodes[i+1].coord[1]) ;
                    nodes.InsertAtEnd(Owner->NewNode(pt[0],pt[1])) ;
                }
            }
        }
    }
}




// %(CArbRmshRegion2D::GenTipNodesHelp-void-|-int-|-int-|-int-|-int-|-int-|-int-|-int-|-int-|-bool-|-ArbMshNode-|*-int-|*)
/* ++ ----------------------------------------------------------
**
**    GenTipNodesHelp - support routine for tip-node generation
**
**      void GenTipNodesHelp(
**              int        ring,
**              int        num,
**              int        a0,
**              int        a1,
**              int        a2,
**              int        b0,
**              int        b1,
**              int        b2,
**              bool       qtrpt_flag,
**              ArbMshNode *nodes,
**              int        *cur)
**
**        ring       - (in)  template ring number
**        num        - (in)  number of nodes to generate
**        a0         - (in)  parameter for node index equation
**        a1         - (in)  parameter for node index equation
**        a2         - (in)  parameter for node index equation
**        b0         - (in)  parameter for node index equation
**        b1         - (in)  parameter for node index equation
**        b2         - (in)  parameter for node index equation
**        qtrpt_flag - (in)  true if this is a quarter-point element
**        nodes      - (i/o) crack template nodes
**        cur        - (i/o) current node number
**
**      Description: This function is a support routine for crack-tip
**          node generation. It computes side-node coordinates.
**
**
** -- */

void CArbRmshRegion2D::GenTipNodesHelp(
                 int ring,int num,
                 int a0,int a1,int a2,
                 int b0,int b1,int b2,
                 bool qtrpt_flag,
                 CArbArray<ArbMshNode> &nodes)
{
    int i, j, k ;
    double x, y ;

    j = (ring == 1) ? 0 : a2*ring*ring + a1*ring + a0 ;
    k = b2*ring*ring + b1*ring + b0 ;

    for (i=0 ; i<num ; ++i) {
        if ((qtrpt_flag) && (ring==1) && (i==0)) {
            x = 0.25*(3.0*nodes[0].coord[0]+nodes[k].coord[0]) ;
            y = 0.25*(3.0*nodes[0].coord[1]+nodes[k].coord[1]) ;
            nodes.InsertAtEnd(Owner->NewNode(x,y)) ;
        } else {
            x = 0.5*(nodes[j+i].coord[0]+nodes[k+i].coord[0]) ;
            y = 0.5*(nodes[j+i].coord[1]+nodes[k+i].coord[1]) ;
            nodes.InsertAtEnd(Owner->NewNode(x,y)) ;
        }
    }
}






void CArbRmshRegion2D::GenerateNewBdryData(CArbMshTopo2D * /*boundary*/,
                                           CornerNodesSet *corner_node)
{
    // loop through all the edges on all the newly generated
    // elements and check to see if they form part of the external
    // boundary or part of a bi-material interface.  If so, then
    // add the edge to the new boundary list

    CArbHashTableIterator<int,ArbMshElement2D> iter(Owner->ElemTable) ;

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.Entry()->elem_id < Owner->FirstNewElem) continue ;

        // loop through the element's edges

        CArbTopoEdgeOnElemCyclicIterator e_iter(Owner->MshTopo,
                                                iter.Entry()->elem_id) ;

        while (corner_node->Fetch(e_iter[0]) == 0) ++e_iter ;

        int first = e_iter[0] ;
        do {
            int nd0, nd1, nd2, nxt ;
            nd0 = e_iter[0] ;
            if (Owner->Order == LINEAR) {
                nd1 = e_iter[1] ;
                nd2 = 0 ;
                nxt = nd1 ;
            } else {
                nd2 = e_iter[1] ;
                ++e_iter ;
                nd1 = e_iter[1] ;
                nxt = nd2 ;
            }

            // find the edge adjacent to the first vertex pointing
            // to the second, and get the ccw element's material type

            CArbTopoAdjVtxIterator v_iter0(Owner->MshTopo,nd0) ;
            while (v_iter0.AdjVtx() != nxt) ++v_iter0 ;

            ArbMshElement2D *elem =
                Owner->ElemTable->Fetch(v_iter0.CcwElem()) ;
            int m_id0 = elem->mat_id ;

            // find the edge adjacent to the second vertex pointing
            // to the first, and see if there is an element, and if
            // so get its material type

            CArbTopoAdjVtxIterator v_iter1(Owner->MshTopo,nxt) ;
            while (v_iter1.AdjVtx() != nd0) ++v_iter1 ;

            if (v_iter1.CcwElem() == NO_ELEM) {
                CArbMshCrackRegion2D::BoundaryData
                             b_data = {v_iter0.CcwElem(),
                                      {nd0,nd1,nd2},false} ;
                Owner->NewBdryTable->Store(ArbEdgeKey(nd0,nd1),b_data) ;
            } else {
                ArbMshElement2D *elem =
                    Owner->ElemTable->Fetch(v_iter1.CcwElem()) ;
                int m_id1 = elem->mat_id ;
                if (m_id0 != m_id1) {
                    CArbMshCrackRegion2D::BoundaryData
                                 b_data = {v_iter0.CcwElem(),
                                          {nd0,nd1,nd2},true} ;
                    Owner->NewBdryTable->Store(ArbEdgeKey(nd0,nd1),
                                               b_data) ;
                }
            }
            ++e_iter ;
        } while (e_iter[0] != first) ;
    }
}




// %(CArbRmshRegion2D::PtUpdate-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|-double-|-ArbMshNode-|*-int-|*)
/* ++ ----------------------------------------------------------
**
**    PtUpdate - support routine for tip-node generation
**
**      CArbCoord2D PtUpdate(
**              CArbCoord2D pt,
**              CArbCoord2D dir,
**              double      mag,
**              ArbMshNode  *nodes,
**              int         *cur)
**
**        pt    - (in)  existing node coordinate
**        dir   - (in)  direction of the offset
**        mag   - (in)  magnitude of the offset
**        nodes - (i/o) crack template nodes
**        cur   - (i/o) current node number
**
**      Description: This function is a support routine for crack tip
**          node generation. It takes a node coordinate and direction
**          and generates a new node that is offset in that direction.
**
**      Return Value: coordinates of the new node
**
**
** -- */

CArbCoord2D CArbRmshRegion2D::PtUpdate(CArbCoord2D pt,
                                       CArbCoord2D dir,
                                       double mag,
                                       CArbArray<ArbMshNode> &nodes)
{
    CArbCoord2D new_pt = pt + dir * mag ;
    nodes.InsertAtEnd(Owner->NewNode(new_pt[0],new_pt[1])) ;
    return(new_pt) ;
}




// %(CArbRmshRegion2D::RebuildNodeTable-void-|)
/* ++ ----------------------------------------------------------
**
**    RebuildNodeTable - rebuild the node table
**
**      void RebuildNodeTable()
**
**      Description: This function rebuilds the node table so that only
**          nodes attacehed to an element are retained.
**
**
** -- */

void CArbRmshRegion2D::RebuildNodeTable()
{
    CArbHashTable<int,CArbCoord2D> *tmp =
        new CArbHashTable<int,CArbCoord2D>() ;

    CArbHashTableIterator<int,ArbMshElement2D> iter(Owner->ElemTable) ;
    int  j ;

    for (iter.First() ; iter.More() ; ++iter) {
        for (j=0 ; j<iter.Entry()->num_nodes ; ++j) {
            CArbCoord2D *cd =
                Owner->NodeTable->Fetch(iter.Entry()->nodes[j]) ;
            tmp->Store(iter.Entry()->nodes[j],*cd) ;
        }
    }
    delete Owner->NodeTable ;
    Owner->NodeTable = tmp ;
}




// %(CArbRmshRegion2D::RemoveTemplateNodes-void-|-int-|-int-|-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    RemoveTemplateNodes - delete template nodes
**
**      void RemoveTemplateNodes(
**              int           num_rings,
**              int           tip_id,
**              CArbMshTopo2D *boundary)
**
**        num_rings - (in)  number of template rings
**        tip_id    - (in)  crack-tip id
**        boundary  - (i/o) mesh topology
**
**      Description: This function deletes the nodes on the crack flank
**          due to template refinement for the previous crack step.
**
**
** -- */

void CArbRmshRegion2D::RemoveTemplateNodes(int num_rings,
                                           int tip_id,
                                           CArbMshTopo2D *boundary)
{
    int i ;

    // start at the old crack tip node.  The first node
    // retrieved by the iterator should be on the crack face

    CArbTopoAdjVtxCyclicIterator iter(boundary,tip_id) ;

    int rmv, next, next_next, num_rmv ;
    int elem = iter.CcwElem() ;

    if (Owner->Order == LINEAR) {
        num_rmv = num_rings-1 ;
        rmv = iter.AdjVtx() ;
        next = tip_id ;
        ++iter ;
        next_next = iter.AdjVtx() ;
    } else {
        num_rmv = 2 * (num_rings-1) ;
        next = iter.AdjVtx() ;
        iter.NewVtx(next) ;
        while (iter.AdjVtx() != tip_id) ++iter ;
        ++iter ;
        rmv = iter.AdjVtx() ;
        ++iter ;
        next_next = iter.AdjVtx() ;
    }

    iter.NewVtx(rmv) ;
    while (iter.AdjVtx() != next) ++iter ;
    ++iter ;
    int prev = iter.AdjVtx() ;

    for (i=0 ; i<num_rmv ; ++i) {
        iter.NewVtx(prev) ;
        while (iter.AdjVtx() != rmv) ++iter ;
        ++iter ;
        int prev_prev = iter.AdjVtx() ;
        boundary->DeleteAngle(next,rmv,next_next) ;
        boundary->DeleteAngle(rmv,prev,next) ;
        boundary->DeleteAngle(prev,prev_prev,rmv) ;
        boundary->InsertAngle(elem,next,prev,next_next) ;
        boundary->InsertAngle(elem,prev,prev_prev,next) ;
        rmv = prev ;
        prev = prev_prev ;
    }

    // If this is a quadratic mesh, then we need to reposition
    // the side node

    if (Owner->Order == QUADRATIC) {
        CArbCoord2D tipc = *(Owner->NodeTable->Fetch(tip_id)) ;
        CArbCoord2D tipr = *(Owner->NodeTable->Fetch(rmv)) ;
        CArbCoord2D tipn = 0.5 * (tipc + tipr) ;
        CArbCoord2D *mid = Owner->NodeTable->Fetch(next) ;
        *mid = tipn ;
    }

    // now we need to go around to the other side of the
    // crack face.  If we have a collapsed crack tip then
    // we loop until we find a node with a different coordinate.

    CArbCoord2D tipc = *(Owner->NodeTable->Fetch(tip_id)) ;
    CArbCoord2D tipr = *(Owner->NodeTable->Fetch(rmv)) ;
    CArbCoord2D diff = tipc - tipr ;
    double tip_tol = 0.001 * diff.Magnitude() ;

    iter.NewVtx(tip_id) ;
    --iter ;
    rmv = iter.AdjVtx() ;
    int tip_n = tip_id ;

    while (1) {
        CArbCoord2D tipp = *(Owner->NodeTable->Fetch(rmv)) ;
        diff = tipc - tipp ;
        double dist = diff.Magnitude() ;

        if (dist > tip_tol) break ;

        tip_n = prev ;
        iter.NewVtx(tip_n) ;
        --iter ;
        rmv = iter.AdjVtx() ;
    }

    // now delete the nodes on this side of the crack

    if (Owner->Order == LINEAR) {
        next = tip_n ;
        iter.NewVtx(next) ;
        while (iter.AdjVtx() != rmv) ++iter ;
        --iter ;
        next_next = iter.AdjVtx() ;
        elem = iter.CcwElem() ;
    } else {
        next = rmv ;
        iter.NewVtx(next) ;
        while (iter.AdjVtx() != tip_n) ++iter ;
        --iter ;
        rmv = iter.AdjVtx() ;
        --iter ;
        next_next = iter.AdjVtx()  ;
        elem = iter.CcwElem() ;
    }
    iter.NewVtx(rmv) ;
    while (iter.AdjVtx() != next) ++iter ;
    --iter ;
    prev = iter.AdjVtx() ;

    for (i=0 ; i<num_rmv ; ++i) {
        iter.NewVtx(prev) ;
        while (iter.AdjVtx() != rmv) ++iter ;
        --iter ;
        int prev_prev = iter.AdjVtx() ;
        boundary->DeleteAngle(next,next_next,rmv) ;
        boundary->DeleteAngle(rmv,next,prev) ;
        boundary->DeleteAngle(prev,rmv,prev_prev) ;
        boundary->InsertAngle(elem,next,next_next,prev) ;
        boundary->InsertAngle(elem,prev,next,prev_prev) ;
        rmv = prev ;
        prev = prev_prev ;
    }

    // If this is a quadratic mesh, then we need to reposition
    // the side node

    if (Owner->Order == QUADRATIC) {
        CArbCoord2D tipc = *(Owner->NodeTable->Fetch(tip_n)) ;
        CArbCoord2D tipr = *(Owner->NodeTable->Fetch(rmv)) ;
        CArbCoord2D tipn = 0.5 * (tipc + tipr) ;
        CArbCoord2D *mid =Owner-> NodeTable->Fetch(next) ;
        *mid = tipn ;
    }
}




// %(CArbRmshRegion2D::RemoveTipNodes-void-|-int-|-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    RemoveTipNodes - delete crack-tip nodes
**
**      void RemoveTipNodes(
**              int           tip_id,
**              CArbMshTopo2D *boundary)
**
**        tip_id   - (in)  crack-tip id
**        boundary - (i/o) mesh topology
**
**      Description: This function deletes multiple nodes at a crack
**          tip due to collapsed but not constrained elements.
**
**
** -- */

void CArbRmshRegion2D::RemoveTipNodes(int tip_id,
                                      CArbMshTopo2D *boundary)
{
    // Note, the current implementation assumes that the
    // crack tip is not at a bi-material interface.

    int i,*tip_nodes,num_tip_nodes ;

    // get a list of the nodes at the old crack tip.  If
    // there is only one node in the list then just return

    tip_nodes = FindTipNodeList(tip_id,&num_tip_nodes,boundary) ;

    if (num_tip_nodes < 2) {
        delete [] tip_nodes ;
        return ;
    }

    int prev = tip_nodes[0] ;
    CArbTopoAdjVtxCyclicIterator iter(boundary,prev) ;
    int prev_prev = iter.AdjVtx() ;

    for (i=1 ; i<num_tip_nodes ; ++i) {
        int rmv = tip_nodes[i] ;
        iter.NewVtx(rmv) ;
        while (iter.AdjVtx() != prev) ++iter ;
        int elem = iter.CcwElem() ;
        --iter ;
        int next = iter.AdjVtx() ;
        iter.NewVtx(next) ;
        while (iter.AdjVtx() != rmv) ++iter ;
        --iter ;
        int next_next = iter.AdjVtx() ;

        boundary->DeleteAngle(prev,prev_prev,rmv) ;
        boundary->DeleteAngle(rmv,prev,next) ;
        boundary->DeleteAngle(next,rmv,next_next) ;
        boundary->InsertAngle(elem,prev,prev_prev,next) ;
        boundary->InsertAngle(elem,next,prev,next_next) ;
    }

    delete [] tip_nodes ;
    return ;
}




// %(CArbRmshRegion2D::ScanLineCross-int-|^const-CArbCoord2D-const|*-CArbCoord2D-const|*-CArbCoord2D-const|*)
/* ++ ----------------------------------------------------------
**
**    ScanLineCross - check for scan line crossing
**
**      int ScanLineCross(
**              const CArbCoord2D *b,
**              const CArbCoord2D *i,
**              const CArbCoord2D *j) const
**
**        b - (in)  point to check
**        i - (in)  first point on line
**        j - (in)  second point on line
**
**      Description: This function determines if a horizontal line
**          drawn from minus infinity to a given point crosses a given
**          line
**
**      Return Value: one if it crosses, zero otherwise
**
**
** -- */

bool CArbRmshRegion2D::ScanCross(const CArbCoord2D &b,
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

    if ((b[1] >= y_max) || (b[1] < y_min)) return(false) ;

    // Find the x coordinates of the intersection and see
    // if this is to the left of the point

    double m = (i[1] - j[1]) / (i[0] - j[0]) ;
    double c = i[1] - m * i[0] ;
    double xint = (b[1] - c)/m ;

    if (xint < b[0]) return(true) ;

    return(false) ;
}


void CArbRmshRegion2D::SetOldCrackTip(CArbCoord2D &tip_coord)
{
    OldCrackTipNode = FindNearestNode(tip_coord) ;
}


// %(CArbRmshRegion2D::SetTipData-bool-|^const-int-const|-CArbMshCrackRegion2D::CrackPt-|*-CArbMshCrackRegion2D::CrackTipData-|*)
/* ++ ----------------------------------------------------------
**
**    SetTipData - set crack-tip parameters
**
**      bool SetTipData(
**              const int                          tip_id,
**              CArbMshCrackRegion2D::CrackPt      *crack_pnt,
**              CArbMshCrackRegion2D::CrackTipData *tip_data) const
**
**        tip_id    - (in)  crack-tip id
**        crack_pnt - (i/o) crack-tip point
**        tip_data  - (i/o) crack-tip parameters
**
**      Description: This function updates missing crack-tip parameters
**          and crack-point data
**
**      Return Value: false if no element contains this crack tip
**
**
** -- */

bool CArbRmshRegion2D::SetTipData(
         const int tip_id,
         CArbMshCrackRegion2D::CrackPt *crack_pnt,
         CArbMshCrackRegion2D::CrackTipData &tip_data) const
{
    FillTipData(tip_id,tip_data) ;
    if (tip_data.number_of_rings == 1) tip_data.progression_ratio = 1.0 ;
    if (!tip_data.temp_radius_set || tip_data.template_radius == 0.0) {
        if (crack_pnt->has_char_size) {
            tip_data.template_radius = crack_pnt->char_elem_size ;
        } else {
            CArbCoord2D tip(crack_pnt->coord[0],
                            crack_pnt->coord[1]) ;
            int tip_elem = FindContainingElem(tip,Owner->MshTopo) ;

            // if we did not find an element check to see if
            // we are close to a boundary
            if (tip_elem == NO_ELEM) {
                int tip_node = FindNearestNode(tip) ;
                tip_data.template_radius =
                    GetAdjCharSize(tip_node) ;
            } else {
                tip_data.template_radius =
                    CharElemSize(tip_elem,Owner->MshTopo) ;
            }
            crack_pnt->has_char_size = true ;
            crack_pnt->char_elem_size = tip_data.template_radius ;
        }
    } else {
        crack_pnt->has_char_size = true ;
        crack_pnt->char_elem_size = tip_data.template_radius ;
    }
    return(true) ;
}




// %(CArbRmshRegion2D::UpdateBoundaries-void-virtual|-CArbCoord2D-|-CArbMshCrackRegion2D::CrackTipData-|&-int-|-CArbElemSet-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateBoundaries - update boundaries near crack-tips
**
**      virtual void UpdateBoundaries(
**              CArbCoord2D                        tip,
**              CArbMshCrackRegion2D::CrackTipData &tip_data,
**              int                                cur_zone,
**              CArbElemSet                        *zones)
**
**        tip      - (in)  crack-tip coordinates
**        tip_data - (in)  crack-tip parameters
**        cur_zone - (in)  current remesh zone number
**        zones    - (in)  set of all remesh zones
**
**      Description: This function refines a boundary if it is close to
**          a crack tip.
**
**
** -- */

void CArbRmshRegion2D::UpdateBoundary(CArbMshTopo2D *boundary,
                const int loop_id,
                const int nnode,
                const int *nodes,
                CArbHashTable<int,int> *mat_table)
{
    // first add this element to the topology

    boundary->InsertElement(loop_id,nnode,nodes) ;

    // now check to see if we need to remove any
    // of the element edges

    for (int i=0 ; i<nnode ; ++i) {
        int j = (i+1) % nnode ;

        // look at the nodes adjacent to the "next"
        // vertex about the element, and find the
        // record where it points back to this vertex.
        // Once we find this, check to see if the
        // the next counter-clockwise element is
        // a real element.  If so we need to
        // update the ccw element id's for all edges
        // adjacent to the new element, and delete the
        // edge.

        CArbTopoAdjVtxIterator iter0(boundary,nodes[j]) ;

        for (iter0.First() ; iter0.More() ; ++iter0 ) {
            if (iter0.AdjVtx() == nodes[i]) {
                if (iter0.CcwElem() != NO_ELEM) {

                    // find edge from the "this" node to the next

                    CArbTopoAdjVtxIterator iter1(boundary,nodes[i]) ;

                    for (iter1.First() ; iter1.More() ; ++iter1) {
                        if (iter1.AdjVtx() == nodes[j]) break ;
                    }

                    int eid = iter0.CcwElem() ;
                    int emid = *(mat_table->Fetch(eid)) ;
                    int curid = iter1.CcwElem() ;
                    int curmid = *(mat_table->Fetch(curid)) ;

                    if (emid == curmid) {

                        // update all the element pointers

                        if (iter1.CcwElem() != iter0.CcwElem()) {

                            // look through all edges and delete
                            // any refereces to this element

                            CArbTopoVtxIterator vtx_iter(*boundary) ;
                            for (vtx_iter.First() ; vtx_iter.More() ; ++vtx_iter) {
                                CArbTopoAdjVtxIterator
                                    iter2(boundary,*vtx_iter) ;
                                for (iter2.First() ; iter2.More() ; ++iter2) {
                                    if (iter2.CcwElem() == curid)
                                        iter2.CcwElemRef() = eid ;
                                }
                            }

                            boundary->DeleteElementRef(curid) ;
                        }

                        // determine the local topology around the
                        // edge that we want to delete.  The local
                        // nodes are labeled as follows:
                        //
                        //      4      1      5
                        //      o------o------o
                        //             |
                        //             |
                        //             |
                        //      o------o------o
                        //      2      0      3

                        // fill the data structure that describes
                        // this edges

                        int local[6] ;
                        local[0] = nodes[i] ;
                        local[1] = nodes[j] ;
                        local[2] = boundary->GetCCWNode(nodes[j],nodes[i]) ;
                        local[3] = boundary->GetCWNode(nodes[i],nodes[j]) ;
                        local[4] = boundary->GetCWNode(nodes[j],nodes[i]) ;
                        local[5] = boundary->GetCCWNode(nodes[i],nodes[j]) ;

                        // classify the ends and update the angles

                        if (local[2] != local[3]) {
                            boundary->DeleteAngle(local[0],local[1],local[2]) ;
                            boundary->DeleteAngle(local[0],local[3],local[1]) ;
                            boundary->InsertAngle(eid,local[0],
                                                  local[3],local[2]) ;
                        } else if (local[2] != local[1]) {
                            boundary->DeleteAngle(local[0],local[1],local[2]) ;
                            boundary->DeleteAngle(local[0],local[3],local[1]) ;
                            boundary->InsertAngle(eid,local[0],
                                                  local[3],local[2]) ;

                        } else {
                            boundary->DeleteAngle(local[0],local[1],local[2]) ;
                        }

                        if (local[4] != local[5]) {
                            boundary->DeleteAngle(local[1],local[4],local[0]) ;
                            boundary->DeleteAngle(local[1],local[0],local[5]) ;
                            boundary->InsertAngle(eid,local[1],
                                                  local[4],local[5]) ;
                        } else if (local[4] != local[0]) {
                            boundary->DeleteAngle(local[1],local[4],local[0]) ;
                            boundary->DeleteAngle(local[1],local[0],local[5]) ;
                            boundary->InsertAngle(eid,local[1],
                                                  local[4],local[5]) ;

                        } else {
                            boundary->DeleteAngle(local[1],local[4],local[0]) ;
                        }
                    } else {

                        // If we get here we have a bi-material interface.
                        // If this edge starts at a corner then save
                        // info about this edge

                        if (Owner->Order == LINEAR) {
                            CArbMshCrackRegion2D::BoundaryData b_data = {
                                loop_id,{nodes[i],nodes[j],0},true} ;
                            Owner->OldBdryTable->Store(
                                ArbEdgeKey(nodes[i],nodes[j]),b_data) ;
                        } else if ((i%2) == 0) {
                            int k = (i+2) % nnode ;
                            CArbMshCrackRegion2D::BoundaryData b_data = {
                                loop_id,{nodes[i],nodes[k],nodes[j]},true} ;
                            Owner->OldBdryTable->Store(
                               ArbEdgeKey(nodes[i],nodes[k]),b_data) ;
                        }
                    }

                } else {

                    // Here we need to check to see if this edge
                    // is part of the external boundary

                    if ((Owner->Order == LINEAR) || ((i%2) == 0)) {

                        CArbTopoAdjVtxIterator biter(Owner->MshTopo,nodes[j]) ;

                        for (biter.First() ; biter.More() ; ++biter) {
                            if (biter.AdjVtx() == nodes[i]) {
                                if (biter.CcwElem() == NO_ELEM) {
                                    if (Owner->Order == LINEAR) {
                                       CArbMshCrackRegion2D::BoundaryData
                                          b_data = {loop_id,
                                          {nodes[i],nodes[j],0},false} ;
                                       Owner->OldBdryTable->Store(
                                          ArbEdgeKey(nodes[i],nodes[j]),
                                          b_data) ;
                                    } else {
                                        int k = (i+2) % nnode ;
                                       CArbMshCrackRegion2D::BoundaryData
                                           b_data = {loop_id,
                                           {nodes[i],nodes[k],nodes[j]},false} ;
                                       Owner->OldBdryTable->Store(
                                           ArbEdgeKey(nodes[i],nodes[k]),
                                           b_data) ;
                                    }
                                }
                            }
                        }
                    }
                }
                break ;
            }
        }
    }
}

