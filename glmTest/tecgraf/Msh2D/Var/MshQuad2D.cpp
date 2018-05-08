//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000
//   All rights reserved
//
// Revision -
//   $Revision: 1.69 $  $Date: 2004/08/10 13:37:00 $  $Author: wash $
//

#include <cstdio>
#include <cmath>
#include <cassert>

#include "List.hpp"
#include "Vec2D.hpp"

#include "MshQuad2D.hpp"
#include "MshSmooth2D.hpp"
#include "FeasableRegion.hpp"

using FTools::List ;
using FTools::Vec2D ;

namespace Msh2D {

#ifdef MEMDEBUG
#include "MemDbg.hpp"
//#define new new(__FILE__,__LINE__)
#endif

//#define DEBUG_LOG 1
#ifdef DEBUG_LOG
FILE *lfd ;
#endif


/* ++ ----------------------------------------------------------
**
**    GenerateQuads - driver to generate quad elements 
**
**      void GenerateQuads()
**
**      Description: This method is a driver for the Q-mesh 
**          quadrilateral mesh generation algorithm. 
**
**
** -- */

#define QUAD_EDGE_ANGLE_TOLERANCE  2.356194490    // 3*pi/4
#define QUAD_NEAR_ANGLE_TOLERANCE  0.523598776    // pi/6
#define SEAM_ANGLE_BIG_TOLERANCE   0.785398163    // pi/4
#define SEAM_ANGLE_SMALL_TOLERANCE 0.785398163    // pi/4
//#define SEAM_ANGLE_SMALL_TOLERANCE 0.523598776    // pi/6
#define QUAD_EDGE_NORMAL_TOLERANCE 3.926990817    // 5*pi/4
//#define QUAD_EDGE_NORMAL_TOLERANCE 3.5342917    // 9*pi/8
//#define QUAD_EDGE_NORMAL_TOLERANCE 3.3379422    // 17*pi/16

#define MAX_ALLOWABLE_ANGLE 2.356194490


#define SMOOTH_BLACKER_TOLERANCE 2.5
//#define SMOOTH_OWEN_TOLERANCE 20.0
#define SMOOTH_OWEN_TOLERANCE 50.0

#define TRANSITION_SEAM_TOLERANCE 2.5
//#define TRANSITION_SEAM_TOLERANCE 2
#define TRANSITION_MIN_FM_TOLERANCE 0.1
#define TRANSITION_SPLIT_TOLERANCE 2.5
#define TRANSITION_TEMPLATE_TOLERANCE 3.0

#define TRIANGLE_SPLIT_TOLERANCE 3.0

#define TRANSITION_ADJACENT_TOLERANCE 20.0
#define TRANSITION_NEAR_BOUNDARY_TOLERANCE 2.0

#define FINAL_ACCEPTABLE_QUALITY_METRIC 0.05

#define SEAM_MASK 12
#define NO_SEAM_MASK 3

#define TWO_PI  6.283185307
#define HALF_PI 1.570796327
#define TRIAL_FACTOR 10

int print_num = 1000000 ;
//int print_num = 0 ;
int quit_num  = 1000000 ;
//bool do_labels = false ;
bool do_labels = true ;

static int debug_num = 0 ;

//int debug_visit = 0 ;
//FILE *dbfd = 0 ;

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

void MshQuad2D::GenerateMesh()
{
    MshRegion2D::GenerateMesh() ;

    // generate quads

    Dict<int,MshElement2D> elem_save = *(pelem_table) ;
    Dict<int,IntNode> node_save = *(pnode_table) ;

    if (GenerateQuads()) {
        RenumberNodes() ;
    } else {
        // quad meshing failed so revert to
        // the saved triangular mesh
        *(pelem_table) = elem_save ;
        *(pnode_table) = node_save ;
    }
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

void MshQuad2D::SetQuadElemAngleChecks(
    double min_angle,
    double max_angle,
    AngleCheckLocation location)
{
    DoQuadAngleChecks = true ;
    MinQuadAngle = min_angle ;
    MaxQuadAngle = max_angle ;
    QuadAngleLocation = location ;
}




bool MshQuad2D::GenerateQuads()
{
    MshTopo2D msh_topo ;
    MshTopo2D quad_topo ;
    MshEdgeList edge_list(pnode_table) ;

//    debug_visit++ ;
//    printf("debug_visit: %d\n",debug_visit) ;

    // Initialize the data structures

    InitializeQuadGen(&msh_topo,&edge_list) ;

    // process the edges

    MshEdgeList::QdEdge *qd_edge ;
    bool saved_flag = false ;
    int saved_0 = 0, saved_1 = 0 ;
    bool modified = true ;
    int stagnate_count = 0 ;

    int num = 0 ;
    int trial = -1 ;
//int chk = 100000 ;
//FILE *fd = 0 ;

    while ((qd_edge = edge_list.GetNextEdge()) != 0) {

        if (MaxNumElem && (num > MaxNumElem)) throw MaxElements() ;

        ++trial ;
        if (MaxNumElem && (trial*TRIAL_FACTOR > MaxNumElem)) throw MaxElements() ;

        debug_num = num ;

        bool transition_flg = false ;

if (num >= print_num) {
    DoPrint(&msh_topo,stdout,true) ;
    DoPrint(&quad_topo,stdout,true) ;
    fprintf(stdout,"a Start_Step_%d\n",num) ;
    fprintf(stdout,"f\n") ;
}

// if (num >= 571) {
//     if (fd == 0) fd = fopen("debug.msh","w") ;
//     DoPrint(&msh_topo,fd) ;
//     DoPrint(&quad_topo,fd) ;
//     fprintf(fd,"a Start Step %d\n",num) ;
//     fprintf(fd,"n\n") ;
//     fflush(fd) ;
//  
//   if (num == 582) {
//     fclose(fd) ;
//     exit(1) ;
//   }
// }


// if ((debug_visit == 53) && (num >= 62)) {
//     if (dbfd == 0) dbfd = fopen("debug.txt","w") ;
//     DoPrint(&msh_topo,dbfd) ;
//     DoPrint(&quad_topo,dbfd) ;
//     fprintf(dbfd,"a Start Step %d\n",num) ;
//     fprintf(dbfd,"n\n") ;
//     fflush(dbfd) ;
//     DoPrint(&msh_topo,dbfd,true) ;
//     DoPrint(&quad_topo,dbfd,true) ;
//     fprintf(dbfd,"n\n") ;
//     fflush(dbfd) ;
// }


        // here is were we try to prevent stagnation.  If
        // the mesh was modified on the last step we set
        // a counter equal to the number of edges in the edge
        // list plus a little fudge.  Each time through,
        // if there was no modifications on the last attempt
        // we decrement the counter.  When we get to zero we
        // have looked at all the edges with no updates so
        // we accept a triangle

        if (modified) {
            modified = false ;
            stagnate_count = int(edge_list.ListLength() * 1.5) ;
        } else {
            stagnate_count-- ;
        }

        // check to see if we are stagnating.  If so, accept
        // one triangle.

        if (stagnate_count <= 0) {
            if (stagnate_count > -(edge_list.ListLength()+1)) {
                int right_node, left_node, right_next, left_next ;
                right_node = edge_list.GetCCWNode(
                             qd_edge->id_0,
                             qd_edge->id_1) ;
                left_node = edge_list.GetCWNode(
                             qd_edge->id_0,
                             qd_edge->id_1) ;
                right_next = edge_list.GetCCWNode(
                             qd_edge->id_1,
                             right_node) ;
                left_next = edge_list.GetCWNode(
                             left_node,
                             qd_edge->id_0) ;

                if (ExtractOneTriangle(&num,qd_edge,&msh_topo,
                                       &quad_topo,&edge_list,
                                       right_node,left_node,
                                       right_next,left_next)) {
                    saved_flag = false ;
                    modified = true ;
                    continue ;
                }
            } else {
                // we have completely stalled out so we just accept
                // all remaining triangles
                break ;
            }
        }

        // check for the special case of closing a triangle
        // or a quad or a pentagon


        if ((qd_edge->end_code & NO_SEAM_MASK) == 3) {
            if (CloseSimplePoly(&num,qd_edge,&msh_topo,
                                &quad_topo,&edge_list)) {
                saved_flag = false ;
                modified = true ;
                continue ;
            }
        }

        // check to see if we need to do any "seaming"

        if (qd_edge->angle_0 < SEAM_ANGLE_BIG_TOLERANCE) {
            if (CloseSeam(true,&num,qd_edge,&msh_topo,&quad_topo,
                          &edge_list,&transition_flg,&modified)) {
                saved_flag = false ;
                modified = true ;
                continue ;
            } else {

                // see if we should form a triangle or other
                // simple polygon close a simple polygon

                if (qd_edge->level == 0) {
                    int top = edge_list.GetCWNode(
                               qd_edge->id_0,qd_edge->id_1) ;
                    if (edge_list.GetEdgeLevel(top,qd_edge->id_0) == 0) {
                        MshEdgeList::NewQuadData qdata ;
                        int top = edge_list.GetCWNode(
                                   qd_edge->id_0,qd_edge->id_1) ;
                        qdata.qcase = 6 ;
                        qdata.tl = qd_edge->id_1 ;
                        qdata.bl = top ;
                        qdata.br = qd_edge->id_0 ;
                        qdata.tr = qd_edge->id_1 ;
                        qdata.front_nodes[0] =
                            edge_list.GetCWNode(top,qd_edge->id_0) ;
                        qdata.front_nodes[1] = qdata.bl ;
                        qdata.front_nodes[2] = qdata.tr ;
                        qdata.front_nodes[3] = edge_list.GetCCWNode(
                            qd_edge->id_0,qd_edge->id_1) ;
                        qdata.num_front_nodes = 4 ;
                        if (!RecoverTopEdge(qdata.tr,qdata.bl,
                                            &msh_topo,&edge_list)) {
                            edge_list.PushToBack(qd_edge) ;
                            continue ;
                        }
                        edge_list.UpdateQuad(&qdata) ;
                        edge_list.UpdateGeom(&qdata) ;
                        int boundary[3] ;
                        boundary[0] = qd_edge->id_0 ;
                        boundary[1] = qd_edge->id_1 ;
                        boundary[2] = top ;
                        MshTopo2D::Edge edge ;
                        edge.nd0 = qd_edge->id_0 ;
                        edge.nd1 = qd_edge->id_1 ;
                        edge.elem0 = msh_topo.GetCCWElem(edge.nd0,edge.nd1) ;
                        edge.elem1 = NO_ELEM ;
                        if (edge.elem0 != NO_ELEM)
                            ClearRegion(3,boundary,&msh_topo,&edge) ;
                        quad_topo.InsertElement(num,3,boundary) ;
                        saved_flag = false ;
                        modified = true ;
                        num++ ;
                        continue ;
                    }
                }
                if (CloseSimplePoly(&num,qd_edge,&msh_topo,
                                    &quad_topo,&edge_list)) {
                    saved_flag = false ;
                    modified = true ;
                    continue ;
                }
            }
        }

        if (qd_edge->angle_1 < SEAM_ANGLE_BIG_TOLERANCE) {
            if (CloseSeam(true,&num,qd_edge,&msh_topo,&quad_topo,
                          &edge_list,&transition_flg,&modified)) {
                saved_flag = false ;
                modified = true ;
                continue ;
            } else {

                // see if we should form a triangle or other
                // simple polygon close a simple polygon

                if (qd_edge->level == 0) {
                    int top = edge_list.GetCCWNode(
                               qd_edge->id_0,qd_edge->id_1) ;
                    if (edge_list.GetEdgeLevel(qd_edge->id_1,top) == 0) {

                        MshEdgeList::NewQuadData qdata ;
                        int top = edge_list.GetCCWNode(
                                   qd_edge->id_0,qd_edge->id_1) ;
                        qdata.qcase = 6 ;
                        qdata.tl = top ;
                        qdata.bl = qd_edge->id_0 ;
                        qdata.br = qd_edge->id_1 ;
                        qdata.tr = top ;
                        qdata.front_nodes[0] = edge_list.GetCWNode(
                            qd_edge->id_0,qd_edge->id_1) ;
                        qdata.front_nodes[1] = qdata.bl ;
                        qdata.front_nodes[2] = qdata.tr ;
                        qdata.front_nodes[3] = edge_list.GetCCWNode(
                            qd_edge->id_1,top) ;
                        qdata.num_front_nodes = 4 ;
                        if (!RecoverTopEdge(qdata.tr,qdata.bl,
                                            &msh_topo,&edge_list)) {
                            edge_list.PushToBack(qd_edge) ;
                            continue ;
                        }
                        edge_list.UpdateQuad(&qdata) ;
                        edge_list.UpdateGeom(&qdata) ;
                        int boundary[3] ;
                        boundary[0] = qd_edge->id_0 ;
                        boundary[1] = qd_edge->id_1 ;
                        boundary[2] = top ;
                        MshTopo2D::Edge edge ;
                        edge.nd0 = qd_edge->id_0 ;
                        edge.nd1 = qd_edge->id_1 ;
                        edge.elem0 = msh_topo.GetCCWElem(edge.nd0,edge.nd1) ;
                        edge.elem1 = NO_ELEM ;
                        if (edge.elem0 != NO_ELEM)
                            ClearRegion(3,boundary,&msh_topo,&edge) ;
                        quad_topo.InsertElement(num,3,boundary) ;
                        saved_flag = false ;
                        modified = true ;
                        num++ ;
                        continue ;
                    }
                }
                if (CloseSimplePoly(&num,qd_edge,&msh_topo,
                                    &quad_topo,&edge_list)) {
                    saved_flag = false ;
                    modified = true ;
                    continue ;
                }
            }
        }

        // if we get here we are doint a quad, so get the
        // appropriate information

        MshEdgeList::NewQuadData qdata ;
        FindQuadData(qd_edge,&qdata,&edge_list,&msh_topo) ;

        // check for the case the right and left nodes are
        // the same.  If this happens, we push the edge back
        // on the edge list

        if ((qdata.qcase == 6) || (qdata.qcase == -1)) {
            edge_list.PushToBack(qd_edge) ;
            ++num ;
            continue ;
        }

        // check to see is we want to do a transition split

        if (TransitionSplit(&num,&qdata,qd_edge,&msh_topo,&quad_topo,
                            &edge_list,&transition_flg)) {
            saved_flag = false ;
            modified = true ;
            continue ;
        }
        if (transition_flg) modified = true ;

        // if we are not currently doing a transion, compute
        // a metric for the proposed element.  If the metric is very
        // negative (invalid), and we have not seen this edge
        // before, then go on to a new edge (saving information
        // about this edge).  If we have seen this edge before
        // then go ahead and use it, and hope that the smoothing
        // process forms a valid quad.

        IntNode *bl = pnode_table->Get(qdata.bl) ;
        IntNode *br = pnode_table->Get(qdata.br) ;
        IntNode *tr = pnode_table->Get(qdata.tr) ;
        IntNode *tl = pnode_table->Get(qdata.tl) ;
        bool constrained = ((bl->motion == MSH_FIXED) &&
                            (br->motion == MSH_FIXED) &&
                            (tr->motion == MSH_FIXED) &&
                            (tl->motion == MSH_FIXED)) ;

        if (!transition_flg) {
//            IntNode *bl = pnode_table->Get(qdata.bl) ;
//            IntNode *br = pnode_table->Get(qdata.br) ;
//            IntNode *tr = pnode_table->Get(qdata.tr) ;
//            IntNode *tl = pnode_table->Get(qdata.tl) ;
            double metric = QuadMetricAreaOnly(bl->coord,br->coord,
                                   tr->coord,tl->coord) ;

            // check to see if any of the existing quads
            // fall inside the proposed element

            int ids[4] = {qdata.bl,qdata.br,qdata.tr,qdata.tl} ;
            Vec2D vts[4] ;
            vts[0] = bl->coord ;
            vts[1] = br->coord ;
            vts[2] = tr->coord ;
            vts[3] = tl->coord ;

//            if ((metric <= 0.0001) || (fmetric <= -1.0) ||

            if ((metric <= 0.0001) ||
                QuadsInPolygon(4,ids,vts,&quad_topo,&edge_list)) {
                edge_list.PushToBack(qd_edge) ;
                if ((saved_0 != qd_edge->id_0) ||
                    (saved_1 != qd_edge->id_1)) {
                    if (!saved_flag) {
                        saved_0 = qd_edge->id_0 ;
                        saved_1 = qd_edge->id_1 ;
                        saved_flag = true ;
                    }
                }
                continue ;
            }
        }

        // now do a series of swaps to recover the top edge
        // if possible.

        if (qdata.qcase == 3) {
            if (!RecoverTopEdge(qdata.tr,qdata.tl,
                                &msh_topo,&edge_list)) {
                edge_list.PushToBack(qd_edge) ;
                continue ;
            }
        } else {
            if (!RecoverTopEdge(qdata.otr,qdata.otl,
                                &msh_topo,&edge_list)) {
                edge_list.PushToBack(qd_edge) ;
                continue ;
            }
        }

        // clear out the triangles that currently cover the
        // new quad

        int boundary[4] ;
        boundary[0] = qdata.bl ;
        boundary[1] = qdata.br ;
        boundary[2] = qdata.tr ;
        boundary[3] = qdata.tl ;

        MshTopo2D::Edge edge ;
        edge.nd0 = qdata.bl ;
        edge.nd1 = qdata.br ;
        edge.elem0 = msh_topo.GetCCWElem(edge.nd0,edge.nd1) ;
        edge.elem1 = NO_ELEM ;
        ClearRegion(4,boundary,&msh_topo,&edge) ;

        // update the edge list to add the quad

        edge_list.UpdateQuad(&qdata) ;
        saved_flag = false ;
        modified = true ;

        // add this quad to the topology

        quad_topo.InsertElement(num,4,boundary) ;
        if (qdata.qcase == 4) {
            num++ ;
            continue ;
        }

        // check for template transition if this new quad is
        // too distorted or the aspect ratio of the adjacent
        // triangle has a very poor aspect ratio

        bool template_tr = false ;
        if ((qdata.qcase == 3) && (qdata.edge_level > 1) &&
             !transition_flg &&
            (qdata.top_template_ratio > TRANSITION_TEMPLATE_TOLERANCE)) {
            DoTemplateSplit(&qdata,&msh_topo,&quad_topo,&edge_list,false) ;
            template_tr = true ;
        }

        // do smoothing of front nodes (three times through is
        // abitrary)

        for (int j=0 ; j<3 ; ++j) {
            SmoothFront(qdata.qcase,
                        qdata.num_front_nodes,
                        qdata.front_nodes,
                        &msh_topo,&quad_topo,&edge_list) ;
            SmoothAdjacent(qdata.num_front_nodes,
                           qdata.front_nodes,&msh_topo,
                           true,true,&edge_list) ;
            SmoothAdjacent(qdata.num_front_nodes,
                           qdata.front_nodes,&quad_topo,
                           false,false,&edge_list) ;

        }
        SmoothFront(qdata.qcase,
                    qdata.num_front_nodes,
                    qdata.front_nodes,
                    &msh_topo,&quad_topo,&edge_list) ;


        if ((qdata.qcase == 3) && !transition_flg && 
            !template_tr && !constrained) {
            if ((AdjAspect(&qdata,&msh_topo) > TRANSITION_ADJACENT_TOLERANCE) ||
                (AdjNearBoundary(&qdata,&msh_topo,&edge_list) > 
                                 TRANSITION_NEAR_BOUNDARY_TOLERANCE) ||
                CornerAspect(&qdata,&msh_topo)) {
                DoTemplateSplit(&qdata,&msh_topo,&quad_topo,
                                &edge_list,false) ;
                for (int j=0 ; j<3 ; ++j) {
                    SmoothFront(qdata.qcase,
                        qdata.num_front_nodes,
                        qdata.front_nodes,
                        &msh_topo,&quad_topo,&edge_list) ;
                    SmoothAdjacent(qdata.num_front_nodes,
                           qdata.front_nodes,&msh_topo,
                           true,true,&edge_list) ;
                    SmoothAdjacent(qdata.num_front_nodes,
                           qdata.front_nodes,&quad_topo,
                           false,false,&edge_list) ;
                }
                SmoothFront(qdata.qcase,
                    qdata.num_front_nodes,
                    qdata.front_nodes,
                    &msh_topo,&quad_topo,&edge_list) ;
            }
        }

        // update the classification of edges on the boundary

        edge_list.UpdateGeom(&qdata) ;
        for (int ii=0 ; ii<BdryVtxCache->Len() ; ++ii) {
            edge_list.UpdateBdryGeom((*BdryVtxCache)[ii]) ;
        }
        BdryVtxCache->Clear() ;

        // get a list of the edges and display

      ++num ;
    }

    // When we get here we are all done adding quads. 
    // if there are any triangles left in the  triangle mesh
    // add these to the quad mesh

    TopoElemIterator titer(&msh_topo) ;
    for (titer.First() ; titer.More() ; ++titer) {
        int nnode,nodes[3] ;
        msh_topo.GetElemNodes(*titer,&nnode,nodes) ;

        if (nnode > 3) throw QuadConversionError() ;

        quad_topo.InsertTriangle(num,nodes[0],nodes[1],nodes[2]) ;
        ++num ;
    }

    // Now we do topological improvements.

//    if (DebugDisplayFlags & RmshQuadBeforeCleanup)
//        DisplayTopo("Quad Mesh Before Cleanup",&quad_topo) ;
    if (DoCleanup) QuadCleanup(&quad_topo) ;

    // now we replace the triangular elements with the quad
    // elements and do global smoothing

    UpdateForQuadMesh(&quad_topo) ;

    // smooth the mesh

    if (DebugDisplayFlags & RmshQuadBeforeSmooth)
        DisplayMesh("Quad Mesh Before Smooth") ;
    
// if (debug_visit == 53) {
//     if (dbfd == 0) dbfd = fopen("debug.txt","w") ;
//     DisplayMesh("Quad Mesh Before Smooth",dbfd) ;
// }
	
	if (DoSmoothNodes) {
        MshSmooth2D smooth(pnode_table,pelem_table) ;
        if (WinslowSmoothing)
            smooth.SmoothNodesWinslow() ;
        else
            smooth.SmoothNodesLaplace() ;
        if (DebugDisplayFlags & RmshQuadAfterSmooth)
            DisplayMesh("Quad Mesh After Wins") ;

        smooth.SmoothNodesConsLaplace() ;
    }
	
    if (DebugDisplayFlags & RmshQuadAfterSmooth)
        DisplayMesh("Quad Mesh After Smooth") ;

    // do angle checks if required
    if (DoQuadAngleChecks) QuadAngleChecks() ;

    return(CheckAllValid(FINAL_ACCEPTABLE_QUALITY_METRIC)) ;
}



/* ++ ----------------------------------------------------------
**
**    CheckAllValid - checks to see if all elements are valid 
**
**      bool *CheckAllValid(double threshold) const
**
**      Description: This check the shape metrics of all the
**          elements. 
**
**      Return Value: true if all shape metrics are >= threshold,
**          false otherwise.
**
** -- */

bool MshQuad2D::CheckAllValid(double threshold) const
{
    // loop through all elements and find the nodes

    //bool valid = true ;
    Dict<int,MshElement2D>::DictIterator elems(pelem_table) ;
    for (elems.First() ; elems.More() ; ++elems) {

        int i ;
        MshElement2D& elem = elems.Entry() ;
        int num_node = ((elem.num_nodes == 3) ||
                        (elem.num_nodes == 6)) ? 3 : 4 ;

        // get the corner node coordinates

        IntNode *nd[4] ;
        for (i=0 ; i<num_node ; ++i) {
            nd[i] = pnode_table->Get(elem.nodes[i]) ;
        } 

        double metric ;
        if (num_node == 3) {
            metric = TriMetric(nd[0]->coord,
                               nd[1]->coord,
                               nd[2]->coord) ;
        } else {
            metric = QuadMetric(nd[0]->coord,nd[1]->coord,
                                nd[2]->coord,nd[3]->coord) ;
        }

        if (metric < threshold) return(false) ;
    }
    return(true) ;
}


double MshQuad2D::AdjAspect(
                             MshEdgeList::NewQuadData *qdata,
                             MshTopo2D *msh_topo) const
{
    // find the coords of the top of the quadrilateral

    IntNode *tr = pnode_table->Get(qdata->tr) ;
    IntNode *tl = pnode_table->Get(qdata->tl) ;

    // get the coordinates of the opposite node on the adjacent
    // triangle element

    int num_nodes,tri_nodes[3] ;
    int tri_elem =
            msh_topo->GetCCWElem(qdata->tl,qdata->tr) ;
    msh_topo->GetElemNodes(tri_elem,qdata->tr,
                           &num_nodes,tri_nodes) ;

    if (num_nodes > 3) throw QuadConversionError() ;

    IntNode *opp = pnode_table->Get(tri_nodes[1]) ;

    // compute the height of the triangle and return the ratio
    // of the base to the height

    double base_len ;
    double height = DistSqr(opp->coord,tl->coord,
                            tr->coord,&base_len) ;

    return(base_len/sqrt(height)) ;
}


double MshQuad2D::AdjNearBoundary(
                             MshEdgeList::NewQuadData *qdata,
                             MshTopo2D *msh_topo,
                             MshEdgeList *edge_list) const
{
    // find the coords of the top of the quadrilateral

    IntNode *tr = pnode_table->Get(qdata->tr) ;
    IntNode *tl = pnode_table->Get(qdata->tl) ;

    // get the coordinates of the opposite node on the adjacent
    // triangle element

    int num_nodes,tri_nodes[3] ;
    int tri_elem =
            msh_topo->GetCCWElem(qdata->tl,qdata->tr) ;
    msh_topo->GetElemNodes(tri_elem,qdata->tr,
                           &num_nodes,tri_nodes) ;

    if (num_nodes > 3) throw QuadConversionError() ;

    // check to see if this opposite node is part of a boundary

    if (!edge_list->ContainsNode(tri_nodes[1])) return(0.0) ;

    IntNode *opp = pnode_table->Get(tri_nodes[1]) ;

    // compute the height of the triangle

    double base_len ;
    double height = DistSqr(opp->coord,tl->coord,
                            tr->coord,&base_len) ;

    if (sqrt(height) > base_len) return(0.0) ;

    // find the average edge length at this node

    double char_len = edge_list->GetCharNodeLength(tri_nodes[1]) ;

    return(base_len/char_len) ;
}


double MshQuad2D::CornerAspect(
                             MshEdgeList::NewQuadData *qdata,
                             MshTopo2D *msh_topo) const
{
    IntNode *br = pnode_table->Get(qdata->br) ;

    if (br->motion != MSH_FIXED) return(false) ;

    IntNode *bl = pnode_table->Get(qdata->bl) ;

    if (bl->motion != MSH_FIXED) return(false) ;

    IntNode *tr = pnode_table->Get(qdata->tr) ;
    IntNode *tl = pnode_table->Get(qdata->tl) ;

    double rlen = (tr->coord - br->coord).Magnitude() ;
    double llen = (tl->coord - bl->coord).Magnitude() ;

    // check to see if we need to do a transition

    double ratio = rlen/llen ;
    if ((ratio > TRANSITION_NEAR_BOUNDARY_TOLERANCE) ||
        (1.0/ratio > TRANSITION_NEAR_BOUNDARY_TOLERANCE)) return(true) ;
    return(false) ;
}



// %(MshQuad2D::CloseSimplePoly-bool-|-int-|*-MshEdgeList::QdEdge-|*-MshTopo2D-|*-MshTopo2D-|*-MshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    CloseSimplePoly - turn a simple polygonal region into an element 
**
**      bool CloseSimplePoly(
**              int                        *current_num,
**              MshEdgeList::QdEdge *qd_edge,
**              MshTopo2D              *msh_topo,
**              MshTopo2D              *quad_topo,
**              MshEdgeList            *edge_list)
**
**        current_num - (i/o) current element number 
**        qd_edge     - (in)  edge description 
**        msh_topo    - (i/o) current triangular mesh description 
**        quad_topo   - (i/o) current quadrilateral mesh description 
**        edge_list   - (i/o) current active boundary 
**
**      Description: This method is used during quadrilateral mesh 
**          generation to turn a simple 3 or 4 sided region into an 
**          element. 
**
**      Return Value: true if an element was created 
**
**
** -- */

bool MshQuad2D::CloseSimplePoly(
                     int *current_num,
                     MshEdgeList::QdEdge *qd_edge,
                     MshTopo2D *msh_topo,
                     MshTopo2D *quad_topo,
                     MshEdgeList *edge_list)
{
    int right_node, left_node ;
    right_node = edge_list->GetCCWNode(qd_edge->id_0,
                                      qd_edge->id_1) ;
    left_node = edge_list->GetCWNode(qd_edge->id_0,
                                     qd_edge->id_1) ;

    // if the left and right nodes are the same, then the
    // remaining void forms a triangle, so make it an element

    if (right_node == left_node) {

        // build a list of the vertex coordinates so that we
        // can see if any quads are inside

        Vec2D vts[3] ;
        vts[0] = (pnode_table->Get(qd_edge->id_0))->coord ;
        vts[1] = (pnode_table->Get(qd_edge->id_1))->coord ;
        vts[2] = (pnode_table->Get(right_node))->coord ;
        int ids[3] ;
        ids[0] = qd_edge->id_0 ;
        ids[1] = qd_edge->id_1 ;
        ids[2] = right_node ;

        // check for anything inside the polygon

        if (QuadsInPolygon(3,ids,vts,quad_topo,edge_list)) return(false) ;

        edge_list->UpdateTriangle(qd_edge->id_0,
                                 qd_edge->id_1,
                                 right_node) ;
        int boundary[3] ;
        boundary[0] = qd_edge->id_0 ;
        boundary[1] = qd_edge->id_1 ;
        boundary[2] = right_node ;

        MshTopo2D::Edge edge ;
        edge.nd0 = qd_edge->id_0 ;
        edge.nd1 = qd_edge->id_1 ;
        edge.elem0 = msh_topo->GetCCWElem(edge.nd0,edge.nd1) ;
        edge.elem1 = NO_ELEM ;
        if (edge.elem0 != NO_ELEM)
            ClearRegion(3,boundary,msh_topo,&edge) ;
        quad_topo->InsertTriangle(*current_num,qd_edge->id_0,
                                 qd_edge->id_1,right_node) ;
        ++(*current_num) ;
        return(true) ;

    // otherwise, we look to see if the nodes adjacent to the
    // left and right are the same.  If this is the case, the
    // remaining void forms a quadrilateral, so make it one
    // or two elements.

    } else {
        int next_node ;
        next_node = edge_list->GetCCWNode(qd_edge->id_1,
                                         right_node) ;
        if (next_node == left_node) {

            MshEdgeList::NewQuadData qdata ;
            qdata.qcase = 4 ;
            qdata.bl = qd_edge->id_0 ;
            qdata.br = qd_edge->id_1 ;
            qdata.tr = right_node ;
            qdata.tl = left_node ;

            int boundary[4] ;
            boundary[0] = qdata.bl ;
            boundary[1] = qdata.br ;
            boundary[2] = qdata.tr ;
            boundary[3] = qdata.tl ;

            IntNode *nodes[4] ;
            nodes[0] = pnode_table->Get(qdata.bl) ;
            nodes[1] = pnode_table->Get(qdata.br) ;
            nodes[2] = pnode_table->Get(qdata.tr) ;
            nodes[3] = pnode_table->Get(qdata.tl) ;

            // check to make sure that there is nothing inside
            // the quad

            Vec2D vts[4] ;
            vts[0] = nodes[0]->coord ;
            vts[1] = nodes[1]->coord ;
            vts[2] = nodes[2]->coord ;
            vts[3] = nodes[3]->coord ;
            int ids[4] ;
            ids[0] = nodes[0]->id ;
            ids[1] = nodes[1]->id ;
            ids[2] = nodes[2]->id ;
            ids[3] = nodes[3]->id ;
            if (QuadsInPolygon(4,ids,vts,quad_topo,edge_list)) return(false) ;

            MshTopo2D::Edge edge ;
            edge.nd0 = qdata.bl ;
            edge.nd1 = qdata.br ;
            edge.elem0 = msh_topo->GetCCWElem(edge.nd0,edge.nd1) ;
            edge.elem1 = NO_ELEM ;
            ClearRegion(4,boundary,msh_topo,&edge) ;

            edge_list->UpdateQuad(&qdata) ;

            // here we need to check for the case where 3
            // of the 4 nodes are fixed, and just adding
            // this as a quad will form a bad element.
            // In this case add two triangles instead.

            int i = 0, j = 1, k = 2, l ;
            int num_fixed = 0 ;

            for (i=0 ; i<4 ; ++i)
                if (nodes[i]->motion == MSH_FIXED) ++num_fixed ;

            if (num_fixed == 3) {
                for (i=0 ; i<4 ; ++i) {
                    j = (i+1) % 4 ;
                    k = (i+2) % 4 ;
                    if ((nodes[i]->motion == MSH_FIXED) &&
                        (nodes[j]->motion == MSH_FIXED) &&
                        (nodes[k]->motion == MSH_FIXED)) break ;
                }
                double angle = Angle(nodes[j]->coord,
                                     nodes[k]->coord,
                                     nodes[i]->coord) ;
                if (angle > MAX_ALLOWABLE_ANGLE) {
                    l = (i+3) % 4 ;
                    quad_topo->InsertTriangle(*current_num,
                                             nodes[i]->id,
                                             nodes[j]->id,
                                             nodes[l]->id) ;
                    ++(*current_num) ;
                    quad_topo->InsertTriangle(*current_num,
                                             nodes[j]->id,
                                             nodes[k]->id,
                                             nodes[l]->id) ;
                    ++(*current_num) ;
                } else {
                    quad_topo->InsertElement(*current_num,4,boundary) ;
                    ++(*current_num) ;
                }
            } else {
                quad_topo->InsertElement(*current_num,4,boundary) ;
                ++(*current_num) ;
            }

            return(true) ;
        }

        // now check for the case of only 5 or 6 edges left on this
        // portion of the boundary.  Adding one edge makes two
        // quads, or a triangle and a quad so we add the edge that
        // makes the quads with the best aspect ratio.

        int poly_nodes[6] ;
        poly_nodes[3] = next_node ;
        poly_nodes[4] = edge_list->GetCCWNode(right_node,
                                             poly_nodes[3]) ;
        poly_nodes[5] = edge_list->GetCCWNode(poly_nodes[3],
                                             poly_nodes[4]) ;

        // first 5 nodes

        if (poly_nodes[4] == left_node) {
            poly_nodes[0] = qd_edge->id_0 ;
            poly_nodes[1] = qd_edge->id_1 ;
            poly_nodes[2] = right_node ;

            Vec2D vts[5] ;
            int ids[5] ;
            IntNode *nodes[5] ;
            for (int i=0 ; i<5 ; ++i) {
                nodes[i] = pnode_table->Get(poly_nodes[i]) ;
                vts[i] = nodes[i]->coord ;
                ids[i] = nodes[i]->id ;
            }

            // check to make sure that there is nothing inside
            // the quad

            if (QuadsInPolygon(5,ids,vts,quad_topo,edge_list)) return(false) ;

            static int cases[5][2][4] = {{{0,1,2,-1},{0,2,3,4}},
                                         {{1,2,3,-1},{1,3,4,0}},
                                         {{2,3,4,-1},{2,4,0,1}},
                                         {{3,4,0,-1},{3,0,1,2}},
                                         {{4,0,1,-1},{4,1,2,3}}} ;
                                         
            // check to see which configuration gives the
            // more favorable element shapes

            double min_metric = -1000.0 ;
            int j,do_case = 0 ;

            for (j=0 ; j<5 ; ++j) {
                double metric0,metric1,tmp ;

                metric0 = TriMetric(nodes[cases[j][0][0]]->coord,
                                    nodes[cases[j][0][1]]->coord,
                                    nodes[cases[j][0][2]]->coord) ;
                if (metric0 < 0) metric0 -= 2.0 ;
                metric1 = QuadMetric(nodes[cases[j][1][0]]->coord,
                                     nodes[cases[j][1][1]]->coord,
                                     nodes[cases[j][1][2]]->coord,
                                     nodes[cases[j][1][3]]->coord) ;
                if (metric1 < 0) metric1 += 1.0 ;
                tmp = metric0 < metric1 ? metric0 : metric1 ;

                if (tmp > min_metric) {
                    min_metric = tmp ;
                    do_case = j ;
                }
            }


            if (min_metric < -0.2) {

                // check for the special
                // case where 3 of the five edges surround one
                // quad.  In this case we delete the quad
                // and form one triange.

                if (CheckSurroundQuad(current_num,poly_nodes,
                        msh_topo,quad_topo,edge_list)) return(true) ;

                // if we cannot find a configuration with valid quads
                // then we try to extract a triangle, otherwise we may
                // get into a loop.
              
                return(ExtractOneTriangle(current_num,qd_edge,msh_topo,
                                   quad_topo,edge_list,
                                   right_node,left_node,
                                   poly_nodes[3],poly_nodes[3])) ;
            }

            if (min_metric < 0.0) {
                for (j=0 ; j<5 ; ++j) {
                     if (nodes[j]->motion == MSH_FIXED) return(
                         ExtractOneTriangle(current_num,qd_edge,msh_topo,
                                   quad_topo,edge_list,
                                   right_node,left_node,
                                   poly_nodes[3],poly_nodes[3])) ;

                }
            }

            // clear the region

            MshTopo2D::Edge edge ;
            edge.nd0 = poly_nodes[0] ;
            edge.nd1 = poly_nodes[1] ;
            edge.elem0 = msh_topo->GetCCWElem(edge.nd0,edge.nd1) ;
            edge.elem1 = NO_ELEM ;
            ClearRegion(5,poly_nodes,msh_topo,&edge) ;

            // add the new elements

            MshEdgeList::NewQuadData qdata ;
            qdata.qcase = 6 ;
            qdata.tl = nodes[cases[do_case][0][2]]->id ;
            qdata.bl = nodes[cases[do_case][0][0]]->id ;
            qdata.br = nodes[cases[do_case][0][1]]->id ;
            qdata.tr = nodes[cases[do_case][0][2]]->id ;
            qdata.front_nodes[0] = nodes[cases[do_case][1][3]]->id ;
            qdata.front_nodes[1] = qdata.bl ;
            qdata.front_nodes[2] = qdata.tr ;
            qdata.front_nodes[3] = nodes[cases[do_case][1][2]]->id ;
            qdata.num_front_nodes = 4 ;
            edge_list->UpdateQuad(&qdata) ;
            edge_list->UpdateGeom(&qdata) ;

            int boundary[4] ;
            for (j=0 ; j<3 ; ++j)
                boundary[j] = nodes[cases[do_case][0][j]]->id ;
            quad_topo->InsertElement(*current_num,3,boundary) ;
            ++(*current_num) ;

            qdata.qcase = 4 ;
            qdata.tl = nodes[cases[do_case][1][0]]->id ;
            qdata.bl = nodes[cases[do_case][1][1]]->id ;
            qdata.br = nodes[cases[do_case][1][2]]->id ;
            qdata.tr = nodes[cases[do_case][1][3]]->id ;
            edge_list->UpdateQuad(&qdata) ;

            for (j=0 ; j<4 ; ++j)
                boundary[j] = nodes[cases[do_case][1][j]]->id ;
            quad_topo->InsertElement(*current_num,4,boundary) ;
            ++(*current_num) ;

            return(true) ;
        }

        // now the 6 sided polygon

        if (poly_nodes[5] == left_node) {
            poly_nodes[0] = qd_edge->id_0 ;
            poly_nodes[1] = qd_edge->id_1 ;
            poly_nodes[2] = right_node ;

            IntNode *nodes[6] ;
            int ids[6] ;
            Vec2D vts[6] ;
            for (int i=0 ; i<6 ; ++i) {
                nodes[i] = pnode_table->Get(poly_nodes[i]) ;
                vts[i] = nodes[i]->coord ;
                ids[i] = nodes[i]->id ;
            }

            // check to make sure that there is nothing inside
            // the quad

            if (QuadsInPolygon(6,ids,vts,quad_topo,edge_list)) return(false) ;

            static int cases[3][2][4] = {{{0,1,2,3},{3,4,5,0}},
                                         {{1,2,3,4},{4,5,0,1}},
                                         {{2,3,4,5},{5,0,1,2}}} ;
                                         
            // check to see which configuration gives the
            // more favorable element shapes

            double min_metric = -1000.0 ;
            int j,do_case = 0 ;

            for (j=0 ; j<3 ; ++j) {
                double metric0,metric1,tmp ;

                metric0 = QuadMetric(nodes[cases[j][0][0]]->coord,
                                     nodes[cases[j][0][1]]->coord,
                                     nodes[cases[j][0][2]]->coord,
                                     nodes[cases[j][0][3]]->coord) ;
                metric1 = QuadMetric(nodes[cases[j][1][0]]->coord,
                                     nodes[cases[j][1][1]]->coord,
                                     nodes[cases[j][1][2]]->coord,
                                     nodes[cases[j][1][3]]->coord) ;
                tmp = metric0 < metric1 ? metric0 : metric1 ;

                if (tmp > min_metric) {
                    min_metric = tmp ;
                    do_case = j ;
                }
            }

            // if we cannot find a configuration with valid quads
            // then we try to extract a triangle, otherwise we may
            // get into a loop.

            if (min_metric < 0) return(
                ExtractOneTriangle(current_num,qd_edge,msh_topo,
                                   quad_topo,edge_list,
                                   right_node,left_node,
                                   poly_nodes[3],poly_nodes[4])) ;

            // clear the region

            MshTopo2D::Edge edge ;
            edge.nd0 = poly_nodes[0] ;
            edge.nd1 = poly_nodes[1] ;
            edge.elem0 = msh_topo->GetCCWElem(edge.nd0,edge.nd1) ;
            edge.elem1 = NO_ELEM ;
            ClearRegion(6,poly_nodes,msh_topo,&edge) ;

            // add the new elements

            MshEdgeList::NewQuadData qdata ;
            qdata.qcase = 3 ;
            qdata.tl = nodes[cases[do_case][0][0]]->id ;
            qdata.bl = nodes[cases[do_case][0][1]]->id ;
            qdata.br = nodes[cases[do_case][0][2]]->id ;
            qdata.tr = nodes[cases[do_case][0][3]]->id ;
            qdata.front_nodes[0] = nodes[cases[do_case][1][2]]->id ;
            qdata.front_nodes[1] = qdata.tl ;
            qdata.front_nodes[2] = qdata.tr ;
            qdata.front_nodes[3] = nodes[cases[do_case][1][1]]->id ;
            qdata.num_front_nodes = 4 ;
            edge_list->UpdateQuad(&qdata) ;
            edge_list->UpdateGeom(&qdata) ;

            int boundary[4] ;
            for (j=0 ; j<4 ; ++j)
                boundary[j] = nodes[cases[do_case][0][j]]->id ;
            quad_topo->InsertElement(*current_num,4,boundary) ;
            ++(*current_num) ;

            qdata.qcase = 4 ;
            qdata.tl = nodes[cases[do_case][1][0]]->id ;
            qdata.bl = nodes[cases[do_case][1][1]]->id ;
            qdata.br = nodes[cases[do_case][1][2]]->id ;
            qdata.tr = nodes[cases[do_case][1][3]]->id ;
            edge_list->UpdateQuad(&qdata) ;

            for (j=0 ; j<4 ; ++j)
                boundary[j] = nodes[cases[do_case][1][j]]->id ;
            quad_topo->InsertElement(*current_num,4,boundary) ;
            ++(*current_num) ;

            return(true) ;
        }
    }
    return(false) ;
}


bool MshQuad2D::CheckSurroundQuad(
    int *current_num,
    int *poly_nodes,
    MshTopo2D *msh_topo,
    MshTopo2D *quad_topo,
    MshEdgeList *edge_list) {

    // check to see if 3 of the edges are adjacent
    // to the same element

    int i,j,adj_elems[5] ;
    for (i=0 ; i<5 ; ++i) {
        j = (i+1) % 5 ;
        adj_elems[i] = quad_topo->GetCCWElem(
            poly_nodes[j],poly_nodes[i]) ;
    }

    int adj_elem = NO_ELEM ;
    for (i=0 ; i<3 ; ++i) {
        int num = 1 ;
        for (j=i+1 ; j<5 ; ++j) {
            if (adj_elems[i] == adj_elems[j]) ++num ;
        }
        if (num == 3) {
            adj_elem = adj_elems[i] ;
            break ;
        }
    }

    if (adj_elem == NO_ELEM) return(false) ;

    // otherwise we delete the quad element

    int num_nodes ;
    int nodes[4] ;
    quad_topo->GetElemNodes(adj_elem,
                            &num_nodes,nodes) ;
    quad_topo->DeleteElement(num_nodes,nodes) ;

    // now delete the triangles and fix up the
    // edge list

    for (i=0 ; i<5 ; ++i) {
        if (adj_elems[i] == adj_elem) {
            j = (i+1) % 5 ;
            int elem = msh_topo->GetCCWElem(
                poly_nodes[i],poly_nodes[j]) ;
            msh_topo->GetElemNodes(elem,
                            &num_nodes,nodes) ;
            msh_topo->DeleteElement(num_nodes,nodes) ;
        }
    }

    // remove the quad

    int bl_indx = 0 ;
    for (i=0 ; i<5 ; ++i) {
        j = (i==0) ? 4 : i-1 ;
        if ((adj_elems[j] != adj_elem) &&
            (adj_elems[i] == adj_elem)) {
            bl_indx = i ;
            break ;
        }
    } 

    MshEdgeList::NewQuadData qdata ;
    qdata.qcase = 1 ;
    qdata.bl = poly_nodes[bl_indx] ;
    qdata.tl = poly_nodes[(bl_indx+1)%5] ;
    qdata.tr = poly_nodes[(bl_indx+2)%5] ;
    qdata.br = poly_nodes[(bl_indx+3)%5] ;
    qdata.num_front_nodes = 6 ;
    qdata.front_nodes[0] = poly_nodes[(bl_indx+4)%5] ;
    qdata.front_nodes[1] = poly_nodes[bl_indx] ;
    qdata.front_nodes[2] = poly_nodes[(bl_indx+1)%5] ;
    qdata.front_nodes[3] = poly_nodes[(bl_indx+2)%5] ;
    qdata.front_nodes[4] = poly_nodes[(bl_indx+3)%5] ;
    qdata.front_nodes[5] = poly_nodes[(bl_indx+4)%5] ;

    edge_list->RemoveQuad(&qdata) ;

    // add the triangle

    int tnodes[3] ;
    tnodes[0] = poly_nodes[bl_indx] ;
    tnodes[1] = poly_nodes[(bl_indx+3)%5] ;
    tnodes[2] = poly_nodes[(bl_indx+4)%5] ;

    edge_list->UpdateTriangle(tnodes[0],
                   tnodes[1],tnodes[2]) ;
    quad_topo->InsertTriangle(*current_num,
                   tnodes[0],tnodes[1],tnodes[2]) ;
    ++(*current_num) ;
    return(true) ;
}


bool MshQuad2D::ExtractOneTriangle(
    int *current_num,
    MshEdgeList::QdEdge *qd_edge,
    MshTopo2D *msh_topo,
    MshTopo2D *quad_topo,
    MshEdgeList *edge_list,
    int right_node,
    int left_node,
    int next_right,
    int next_left)
{
    // find the triangular element adjacent to the front edge

    int num_nodes, ids[4] ;
    int tri_elem = msh_topo->GetCCWElem(qd_edge->id_0,qd_edge->id_1) ;
    msh_topo->GetElemNodes(tri_elem,&num_nodes,ids) ;

    // build a list of the vertex coordinates so that we
    // can see if any quads are inside

    Vec2D vts[3] ;
    vts[0] = (pnode_table->Get(ids[0]))->coord ;
    vts[1] = (pnode_table->Get(ids[1]))->coord ;
    vts[2] = (pnode_table->Get(ids[2]))->coord ;

    // check for anything inside the polygon

    if (QuadsInPolygon(3,ids,vts,quad_topo,edge_list)) return(false) ;

    // update the front, first find the vertex not
    // on the edge.

    int off_vt = 0 ;
    for (int i=0 ; i<3 ; ++i) {
        if ((ids[i] != qd_edge->id_0) && (ids[i] != qd_edge->id_1)) {
            off_vt = ids[i] ;
            break ;
        }
    }

    if (off_vt == right_node) {

        MshEdgeList::NewQuadData qdata ;
        qdata.qcase = 6 ;
        qdata.bl = qd_edge->id_0 ;
        qdata.br = qd_edge->id_1 ;
        qdata.tl = off_vt ;
        qdata.tr = off_vt ;
        qdata.num_front_nodes = 4 ;
        qdata.front_nodes[0] = left_node ;
        qdata.front_nodes[1] = qdata.bl ;
        qdata.front_nodes[2] = qdata.tr ;
        qdata.front_nodes[3] = next_right ;
        edge_list->UpdateQuad(&qdata) ;
        edge_list->UpdateGeom(&qdata) ;

    } else if (off_vt == left_node) {

        MshEdgeList::NewQuadData qdata ;
        qdata.qcase = 6 ;
        qdata.bl = off_vt ;
        qdata.br = qd_edge->id_0 ;
        qdata.tl = qd_edge->id_1 ;
        qdata.tr = qd_edge->id_1 ;
        qdata.num_front_nodes = 4 ;
        qdata.front_nodes[0] = next_left ;
        qdata.front_nodes[1] = qdata.bl ;
        qdata.front_nodes[2] = qdata.tr ;
        qdata.front_nodes[3] = right_node ;
        edge_list->UpdateQuad(&qdata) ;
        edge_list->UpdateGeom(&qdata) ;

    } else {

        return(false) ;

    }

    // update the triangular mesh

    MshTopo2D::Edge edge ;
    edge.nd0 = qd_edge->id_0 ;
    edge.nd1 = qd_edge->id_1 ;
    edge.elem0 = tri_elem ;
    edge.elem1 = NO_ELEM ;
    ClearRegion(3,ids,msh_topo,&edge) ;
    quad_topo->InsertTriangle(*current_num,qd_edge->id_0,
                              qd_edge->id_1,off_vt) ;
    ++(*current_num) ;
    return(true) ;
}

// %(MshQuad2D::CloseSeam-bool-|-bool-|-int-|*-MshEdgeList::QdEdge-|*-MshTopo2D-|*-MshTopo2D-|*-MshEdgeList-|*-bool-|*)
/* ++ ----------------------------------------------------------
**
**    CloseSeam - performs a seam closing 
**
**      bool CloseSeam(
**              bool                       left,
**              int                        *current_num,
**              MshEdgeList::QdEdge *qd_edge,
**              MshTopo2D              *msh_topo,
**              MshTopo2D              *quad_topo,
**              MshEdgeList            *edge_list,
**              bool                       *transition_flag)
**
**        left            - (in)  true if seam end is on the left 
**        current_num     - (i/o) current element number 
**        qd_edge         - (in)  edge description 
**        msh_topo        - (i/o) current triangular mesh description 
**        quad_topo       - (i/o) current quadrilateral mesh 
**                                description 
**        edge_list       - (i/o) current active boundary 
**        transition_flag - (out) true if a transition seam was 
**                                performed 
**
**      Description: This method checks to see if a seaming operation 
**          can be performed and if so does it. 
**
**      Return Value: true if a seaming option was performed. 
**
**
** -- */

bool MshQuad2D::CloseSeam(bool left,
                     int *current_num,
                     MshEdgeList::QdEdge *qd_edge,
                     MshTopo2D *msh_topo,
                     MshTopo2D *quad_topo,
                     MshEdgeList *edge_list,
                     bool *transition_flg,
                     bool *modified_flg)
{
    // Find the number of elements adjacent to the seam point.
    // We only close the seam if there is more than one and less
    // than five adjacent elements

    int num_adj ;
    *modified_flg = true ;

    if (left) {
        num_adj = quad_topo->HasVtx(qd_edge->id_0) ?
                    quad_topo->NumAdjElems(qd_edge->id_0) : 2 ;
    } else {
        num_adj = quad_topo->HasVtx(qd_edge->id_1) ?
                    quad_topo->NumAdjElems(qd_edge->id_1) : 2 ;
    }
    double angle = left ? qd_edge->angle_0 : qd_edge->angle_1 ;

    if (num_adj > 1) {
        if ((num_adj < 5) ||
            (angle < SEAM_ANGLE_SMALL_TOLERANCE)) {

            MshEdgeList::NewSeamData sdata ;
            bool status = FindSeamData(qd_edge,&sdata,edge_list,
                                       msh_topo,quad_topo) ;
            *modified_flg = false ;
            if (!status && (sdata.cw != sdata.ccw)) return(false) ;

            // if either of the legs of the seam are not movable, enforce
            // a smaller seam angle tolerance

            IntNode *tid0 = pnode_table->Get(sdata.id0) ;
            IntNode *tid2 = pnode_table->Get(sdata.id2) ;
            if ((tid0->motion == MSH_FIXED) ||
                (tid2->motion == MSH_FIXED)) {
                if (angle > QUAD_NEAR_ANGLE_TOLERANCE) {
                    *modified_flg = false ;
                    return(false) ;
                }
            }

            // check to see if the region inside the seam is
            // only a triangle.  If this is the case, then
            // all we can do is make a triangular element

            if ((sdata.id0 == sdata.ccw) &&
                (sdata.id2 == sdata.cw)) {

                edge_list->UpdateTriangle(sdata.id0,
                                         sdata.id1,
                                         sdata.id2) ;
                int boundary[3] ;
                boundary[0] = sdata.id0 ;
                boundary[1] = sdata.id1 ;
                boundary[2] = sdata.id2 ;

                MshTopo2D::Edge edge ;
                edge.nd0 = sdata.id0 ;
                edge.nd1 = sdata.id1 ;
                edge.elem0 = msh_topo->GetCCWElem(edge.nd0,edge.nd1) ;
                edge.elem1 = NO_ELEM ;
                if (edge.elem0 != NO_ELEM)
                    ClearRegion(3,boundary,msh_topo,&edge) ;
                quad_topo->InsertTriangle(*current_num,sdata.id0,
                                          sdata.id1,sdata.id2) ;
                ++(*current_num) ;
                return(true) ;

            // check for the case where we want to insert a quad

            } else if (!status && (sdata.cw == sdata.ccw)) {
                MshEdgeList::NewQuadData qdata ;
                qdata.obl = qdata.bl = sdata.id1 ;
                qdata.obr = qdata.br = sdata.id2 ;
                qdata.otr = qdata.tr = sdata.cw ;
                qdata.otl = qdata.tl = sdata.id0 ;
                qdata.qcase = 4 ;

                int boundary[4] ;
                boundary[0] = qdata.obl ;
                boundary[1] = qdata.obr ;
                boundary[2] = qdata.otr ;
                boundary[3] = qdata.otl ;

                MshTopo2D::Edge edge ;
                edge.nd0 = qdata.obl ;
                edge.nd1 = qdata.obr ;
                edge.elem0 = msh_topo->GetCCWElem(edge.nd0,edge.nd1) ;
                edge.elem1 = NO_ELEM ;
                ClearRegion(4,boundary,msh_topo,&edge) ;

                edge_list->UpdateQuad(&qdata) ;
                quad_topo->InsertElement(*current_num,4,boundary) ;
                ++(*current_num) ;
                return(true) ;

            // check the for the case where we want to close 
            // the seam and merge the two adjacent quads  The
            // actually merging will get done during the
            // topological cleanup phase.

            } else if (sdata.merge_quads) {

                if (!DoSeam(&sdata,msh_topo, quad_topo,edge_list)) {
                    edge_list->PushToBack(qd_edge) ;
                    *modified_flg = false ;
                } else {
                    if (sdata.cw != sdata.ccw)
                        edge_list->UpdateSeamGeom(&sdata) ;
                }
                ++(*current_num) ;
                return(true) ;

            // else check to see if we want to apply a
            // seam closing template

            } else if (sdata.do_template) {
                MshEdgeList::NewQuadData tqdata ;
                if (sdata.template_cw) {
                    int elem = !quad_topo->HasVtx(sdata.id1) ? NO_ELEM : 
                               quad_topo->GetCCWElem(sdata.id1,sdata.id0) ;
                    if (elem != NO_ELEM) {
                        int nnum, nodes[4] ;
                        quad_topo->GetElemNodes(elem,sdata.id1,&nnum,nodes) ;

                        if (nnum > 4) throw QuadConversionError() ;
                        tqdata.bl = nodes[2] ;
                        tqdata.br = nodes[3] ;
                        tqdata.tl = nodes[1] ;
                        tqdata.tr = nodes[0] ;
                        tqdata.num_front_nodes = 4 ;
                        tqdata.front_nodes[0] =
                            edge_list->GetCWNode(tqdata.tl,tqdata.tr) ;
                        tqdata.front_nodes[1] = tqdata.tl ;
                        tqdata.front_nodes[2] = tqdata.tr ;
                        tqdata.front_nodes[3] =
                            edge_list->GetCCWNode(tqdata.tl,tqdata.tr) ;

                        DoTemplateSplit(&tqdata,msh_topo,quad_topo,
                                        edge_list,true) ;
                        ++(*current_num) ;
                        return(true) ;
                    }
                } else {
                    int elem = !quad_topo->HasVtx(sdata.id2) ? NO_ELEM :
                               quad_topo->GetCCWElem(sdata.id2,sdata.id1) ;
                    if (elem != NO_ELEM) {
                        int nnum, nodes[4] ;
                        quad_topo->GetElemNodes(elem,sdata.id1,&nnum,nodes) ;

                        if (nnum > 4) throw QuadConversionError() ;

                        tqdata.bl = nodes[1] ;
                        tqdata.br = nodes[2] ;
                        tqdata.tl = nodes[0] ;
                        tqdata.tr = nodes[3] ;
                        tqdata.num_front_nodes = 4 ;
                        tqdata.front_nodes[0] =
                            edge_list->GetCWNode(tqdata.tl,tqdata.tr) ;
                        tqdata.front_nodes[1] = tqdata.tl ;
                        tqdata.front_nodes[2] = tqdata.tr ;
                        tqdata.front_nodes[3] =
                            edge_list->GetCCWNode(tqdata.tl,tqdata.tr) ;

                        DoTemplateSplit(&tqdata,msh_topo,quad_topo,
                                        edge_list,true) ;
                        ++(*current_num) ;
                        return(true) ;
                    }
                }

            // else check to see if we should do a transition seam

            } else if (((sdata.ratio > TRANSITION_SEAM_TOLERANCE) ||  
                        (sdata.ratio < 1.0/TRANSITION_SEAM_TOLERANCE)) &&
                        (sdata.min_fm > TRANSITION_MIN_FM_TOLERANCE)) { 

//            } else if ((sdata.ratio > TRANSITION_SEAM_TOLERANCE) ||  
//                       (sdata.ratio < 1.0/TRANSITION_SEAM_TOLERANCE)) { 

if (sdata.min_fm < TRANSITION_MIN_FM_TOLERANCE) {
    fprintf(stderr,"Rejecting trans\n") ;
}

                if (DoTransitionSeam(qd_edge,&sdata,msh_topo,
                                     quad_topo,edge_list)) {
                    *transition_flg = true ;
                } else {
                    edge_list->PushToBack(qd_edge) ;
                    *modified_flg = false ;
                    return(false) ;
                }

            // else a "regular" seam

            } else if (status) {

                if (!DoSeam(&sdata,msh_topo, quad_topo,edge_list)) {
                    edge_list->PushToBack(qd_edge) ;
                    *modified_flg = false ;
                } else {
                    SmoothSeamFront(&sdata,msh_topo,
                                    quad_topo,edge_list) ;
                    if (sdata.cw != sdata.ccw)
                        edge_list->UpdateSeamGeom(&sdata) ;
                }
                ++(*current_num) ;
                return(true) ;

            } else {
                edge_list->PushToBack(qd_edge) ;
                *modified_flg = false ;
                return(false) ;
            }
        }
    }
    *modified_flg = false ;
    return(false) ;
}




// %(MshQuad2D::TransitionSplit-bool-|-int-|*-MshEdgeList::NewQuadData-|*-MshEdgeList::QdEdge-|*-MshTopo2D-|*-MshTopo2D-|*-MshEdgeList-|*-bool-|*)
/* ++ ----------------------------------------------------------
**
**    TransitionSplit - performs a transition split 
**
**      bool TransitionSplit(
**              int                             *current_num,
**              MshEdgeList::NewQuadData *qdata,
**              MshEdgeList::QdEdge      *qd_edge,
**              MshTopo2D                   *msh_topo,
**              MshTopo2D                   *quad_topo,
**              MshEdgeList                 *edge_list,
**              bool                            *transition_flg)
**
**        current_num    - (in)  current element number 
**        qdata          - (i/o) proposed quad elem data 
**        qd_edge        - (i/o) edge description 
**        msh_topo       - (i/o) current triangular mesh description 
**        quad_topo      - (i/o) current quadrilateral mesh description 
**        edge_list      - (i/o) current active boundary 
**        transition_flg - (out) true if a split was performed 
**                               (redundant with function value) 
**
**      Description: This method checks to see if a transion split 
**          operation can be performed and if so does it. 
**
**      Return Value: true if a transition split was performed 
**
**
** -- */

bool MshQuad2D::TransitionSplit(
                     int *current_num,
                     MshEdgeList::NewQuadData *qdata,
                     MshEdgeList::QdEdge *qd_edge,
                     MshTopo2D *msh_topo,
                     MshTopo2D *quad_topo,
                     MshEdgeList *edge_list,
                     bool *transition_flg)
{
    // In this routine we check to see if we can do a
    // a transition split, and if so do it.  first check for
    // quad case 2 (adding two new edges to form the quad) where
    // we are within the split tolerances and this is a valid
    // split

    if ((qdata->qcase == 2) && (qdata->edge_level > 1) &&
        ((qdata->ratio > TRANSITION_SPLIT_TOLERANCE) ||  
         (qdata->ratio < 1.0/TRANSITION_SPLIT_TOLERANCE)) &&
          qdata->valid_transition_split) {

        DoTransitionSplit(qd_edge,qdata,msh_topo,
                          quad_topo,edge_list) ;
        SmoothFront(qdata->qcase,
                    3,&qdata->front_nodes[1],
                    msh_topo,quad_topo,edge_list) ;
        edge_list->UpdateGeomNodes(3,&qdata->front_nodes[1],
                                   qdata->edge_level) ;
        FindQuadData(qd_edge,qdata,edge_list,msh_topo) ;
        *transition_flg = true ;

    // If we have quad case 3 (adding one new edge to form the
    // the quad).  If there is a big disparity between the
    // lengths of the proposed elements sides, or top and bottom
    // we perform a template split.  This takes the element
    // adjacent to the side or bottom and splits it into 4
    // elements, adding 2 new nodes along an edge.

    // for the cases where we split a side, we fix up the
    // the qdata structure and "fall" through to add the
    // quad. 

    // for the case where we split the bottom, we return
    // a true flag which means that no element will be formed
    // but we will go back to the top of the quad building loop
    // and get a new edge from the active boundary.

    } else if (!(*transition_flg) && (qdata->qcase == 3) &&
                (qdata->edge_level > 1)) {

        // check for spliting the left side

        if (qdata->ratio > TRANSITION_TEMPLATE_TOLERANCE) {
            MshEdgeList::NewQuadData tqdata ;
            int elem = !quad_topo->HasVtx(qdata->bl) ? NO_ELEM :
                       quad_topo->GetCCWElem(qdata->bl,qdata->tl) ;
            if (elem != NO_ELEM) {
                int num, nodes[4] ;
                quad_topo->GetElemNodes(elem,qdata->bl,&num,nodes) ;
                if (num != 4) return(false) ;
                tqdata.bl = nodes[2] ;
                tqdata.br = nodes[3] ;
                tqdata.tl = nodes[1] ;
                tqdata.tr = nodes[0] ;
                tqdata.num_front_nodes = 4 ;
                tqdata.front_nodes[0] =
                    edge_list->GetCWNode(tqdata.tl,tqdata.tr) ;
                tqdata.front_nodes[1] = tqdata.tl ;
                tqdata.front_nodes[2] = tqdata.tr ;
                tqdata.front_nodes[3] =
                    edge_list->GetCCWNode(tqdata.tl,tqdata.tr) ;
                DoTemplateSplit(&tqdata,msh_topo,quad_topo,
                                edge_list,true) ;
                int tls = qdata->tl ;
                qdata->tl = edge_list->GetCWNode(qdata->bl,qdata->br) ;
                if (qdata->otl == tls) qdata->otl = qdata->tl ;
                if (qdata->otr == tls) qdata->otr = qdata->tl ;
                qdata->front_nodes[1] = qdata->tl ;
                qdata->front_nodes[0] =
                    edge_list->GetCWNode(qdata->tl,qdata->bl) ;
                *transition_flg = true ;
            }

        // check for spliting the right side

        } else if (qdata->ratio < 1.0/TRANSITION_TEMPLATE_TOLERANCE) {
            MshEdgeList::NewQuadData tqdata ;
            int elem = !quad_topo->HasVtx(qdata->tr) ? NO_ELEM :
                       quad_topo->GetCCWElem(qdata->tr,qdata->br) ;
            if (elem != NO_ELEM) {
                int num, nodes[4] ;
                quad_topo->GetElemNodes(elem,qdata->br,&num,nodes) ;
                if (num != 4) return(false) ;
                tqdata.bl = nodes[1] ;
                tqdata.br = nodes[2] ;
                tqdata.tl = nodes[0] ;
                tqdata.tr = nodes[3] ;
                tqdata.num_front_nodes = 4 ;
                tqdata.front_nodes[0] =
                    edge_list->GetCWNode(tqdata.tl,tqdata.tr) ;
                tqdata.front_nodes[1] = tqdata.tl ;
                tqdata.front_nodes[2] = tqdata.tr ;
                tqdata.front_nodes[3] =
                    edge_list->GetCCWNode(tqdata.tl,tqdata.tr) ;
                DoTemplateSplit(&tqdata,msh_topo,
                                quad_topo,edge_list,true) ;
                int trs = qdata->tr ;
                qdata->tr = edge_list->GetCCWNode(qdata->bl,qdata->br) ;
                if (qdata->otl == trs) qdata->otl = qdata->tr ;
                if (qdata->otr == trs) qdata->otr = qdata->tr ;
                qdata->front_nodes[2] = qdata->tr ;
                qdata->front_nodes[3] =
                    edge_list->GetCCWNode(qdata->br,qdata->tr) ;
                *transition_flg = true ;
            }

        // check for spliting the bottom

        } else if (qdata->base_template_ratio > TRANSITION_TEMPLATE_TOLERANCE) {
            MshEdgeList::NewQuadData tqdata ;
            int elem = !quad_topo->HasVtx(qdata->br) ? NO_ELEM :
                       quad_topo->GetCCWElem(qdata->br,qdata->bl) ;
            if (elem == NO_ELEM) return(false) ;
            int num_node, nodes[4] ;
            quad_topo->GetElemNodes(elem,qdata->br,&num_node,nodes) ;
            if (num_node != 4) return(false) ;
            tqdata.bl = nodes[2] ;
            tqdata.br = nodes[3] ;
            tqdata.tl = nodes[1] ;
            tqdata.tr = nodes[0] ;
            tqdata.num_front_nodes = 4 ;
            tqdata.front_nodes[0] = qdata->tl ;
            tqdata.front_nodes[1] = tqdata.tl ;
            tqdata.front_nodes[2] = tqdata.tr ;
            tqdata.front_nodes[3] = qdata->tr ;
            DoTemplateSplit(&tqdata,msh_topo,quad_topo,edge_list,true) ;
            *transition_flg = true ;
            ++(*current_num) ;
            return(true) ;
        }
    }

    // Here we check for a special case for doing transitions
    // starting from a boundary where we use a different tolerance

    if (!(*transition_flg) && (qdata->qcase == 3) &&
                (qdata->edge_level == 0)) {
        IntNode *tl = pnode_table->Get(qdata->tl) ;
        IntNode *bl = pnode_table->Get(qdata->bl) ;
        IntNode *br = pnode_table->Get(qdata->br) ;
        double blen = (bl->coord - br->coord).Magnitude() ;
        double llen = (bl->coord - tl->coord).Magnitude() ;
        if (tl->motion == MSH_FIXED) {
            double ratio = llen / blen ;
            if (ratio > TRANSITION_NEAR_BOUNDARY_TOLERANCE) {
                MshEdgeList::NewQuadData tqdata ;
                int elem = !quad_topo->HasVtx(qdata->bl) ? NO_ELEM :
                           quad_topo->GetCCWElem(qdata->bl,qdata->tl) ;
                if (elem != NO_ELEM) {
                    int num, nodes[4] ;
                    quad_topo->GetElemNodes(elem,qdata->bl,&num,nodes) ;
                    if (num != 4) return(false) ;
                    tqdata.bl = nodes[2] ;
                    tqdata.br = nodes[3] ;
                    tqdata.tl = nodes[1] ;
                    tqdata.tr = nodes[0] ;
                    tqdata.num_front_nodes = 4 ;
                    tqdata.front_nodes[0] =
                        edge_list->GetCWNode(tqdata.tl,tqdata.tr) ;
                    tqdata.front_nodes[1] = tqdata.tl ;
                    tqdata.front_nodes[2] = tqdata.tr ;
                    tqdata.front_nodes[3] =
                        edge_list->GetCCWNode(tqdata.tl,tqdata.tr) ;
                    DoTemplateSplit(&tqdata,msh_topo,quad_topo,
                                    edge_list,true) ;
                    int tls = qdata->tl ;
                    qdata->tl = edge_list->GetCWNode(qdata->bl,qdata->br) ;
                    if (qdata->otl == tls) qdata->otl = qdata->tl ;
                    if (qdata->otr == tls) qdata->otr = qdata->tl ;
                    qdata->front_nodes[1] = qdata->tl ;
                    qdata->front_nodes[0] =
                        edge_list->GetCWNode(qdata->tl,qdata->bl) ;
                    *transition_flg = true ;
                }
            }
        }
        IntNode *tr = pnode_table->Get(qdata->tr) ;
        double rlen = (bl->coord - tr->coord).Magnitude() ;
        if (tr->motion == MSH_FIXED) {
            double ratio = rlen / blen ;
            if (ratio > TRANSITION_NEAR_BOUNDARY_TOLERANCE) {
                MshEdgeList::NewQuadData tqdata ;
                int elem = !quad_topo->HasVtx(qdata->tr) ? NO_ELEM :
                       quad_topo->GetCCWElem(qdata->tr,qdata->br) ;
                if (elem != NO_ELEM) {
                    int num, nodes[4] ;
                    quad_topo->GetElemNodes(elem,qdata->br,&num,nodes) ;
                    if (num != 4) return(false) ;
                    tqdata.bl = nodes[1] ;
                    tqdata.br = nodes[2] ;
                    tqdata.tl = nodes[0] ;
                    tqdata.tr = nodes[3] ;
                    tqdata.num_front_nodes = 4 ;
                    tqdata.front_nodes[0] =
                        edge_list->GetCWNode(tqdata.tl,tqdata.tr) ;
                    tqdata.front_nodes[1] = tqdata.tl ;
                    tqdata.front_nodes[2] = tqdata.tr ;
                    tqdata.front_nodes[3] =
                        edge_list->GetCCWNode(tqdata.tl,tqdata.tr) ;
                    DoTemplateSplit(&tqdata,msh_topo,
                                    quad_topo,edge_list,true) ;
                    int trs = qdata->tr ;
                    qdata->tr = edge_list->GetCCWNode(qdata->bl,qdata->br) ;
                    if (qdata->otl == trs) qdata->otl = qdata->tr ;
                    if (qdata->otr == trs) qdata->otr = qdata->tr ;
                    qdata->front_nodes[2] = qdata->tr ;
                    qdata->front_nodes[3] =
                        edge_list->GetCCWNode(qdata->br,qdata->tr) ;
                    *transition_flg = true ;
                }
            }
        }
    }

    return(false) ;
}




// %(MshQuad2D::UpdateForQuadMesh-void-|-MshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateForQuadMesh - update the region with quad mesh information 
**
**      void UpdateForQuadMesh(MshTopo2D *quad_topo)
**
**        quad_topo - (in)  quadrilateral mesh data 
**
**      Description: This method updates the instance to replace the 
**          triangular mesh information with quad mesh information. 
**
**
** -- */

void MshQuad2D::UpdateForQuadMesh(MshTopo2D *quad_topo)
{
    // If these are quadratic order elements then we need to
    // build a data structure to keep track of the side nodes.

    Dict<EdgeKey,int> side_nodes ;
    if (Order == QUADRATIC) {
        Dict<int,MshElement2D>::DictIterator eiter(pelem_table) ;
        int nn, j, nd0, nd1, nd2 ;
        for (eiter.First() ; eiter.More() ; ++eiter) {
            MshElement2D& elem = eiter.Entry() ;
            nn = elem.num_nodes/2 ;
            for (j=0 ; j<nn ; ++j) {
                nd0 = elem.nodes[j] ;
                nd1 = elem.nodes[(j+1) % nn] ;
                nd2 = elem.nodes[j + nn] ;
                if (nd0 > nd1)
                    side_nodes.Store(EdgeKey(nd0,nd1),nd2) ;
                else
                    side_nodes.Store(EdgeKey(nd1,nd0),nd2) ;
            }
        }
    }

    // now delete the triangular mesh and build a quad mesh

    delete pelem_table ;
    pelem_table = new Dict<int,MshElement2D> ;

    int *elems = quad_topo->GetElemList() ;

    for (int jj=0 ; jj<quad_topo->NumElements() ; ++jj) {
        MshElement2D elem ;
        int num_nodes, nodes[4] ;

        quad_topo->GetElemNodes(elems[jj],&num_nodes,nodes) ;

        if (num_nodes > 4) throw QuadConversionError() ;

        elem.elem_id = StartElemIdSave + jj ;
        elem.mat_id = MatId ;

        if (Order == QUADRATIC) {
            elem.num_nodes = num_nodes * 2 ;
            for (int ii=0 ; ii<num_nodes ; ++ii) {
                int kk = (ii+1) % num_nodes ;
                elem.nodes[ii] = nodes[ii] ;

                // check to see this edge and thus the side
                // node already exist

                int *side ;
                if (nodes[ii] > nodes[kk])
                    side = side_nodes.Get(EdgeKey(nodes[ii],nodes[kk])) ;
                else
                    side = side_nodes.Get(EdgeKey(nodes[kk],nodes[ii])) ;

                if (side != 0) {
                    elem.nodes[ii+num_nodes] = *side ;
                } else {
                   
                    // new node, don't worry about the coords
                    // because they will be udated during smoothing

                    elem.nodes[ii+num_nodes] =
                        NewNode(0.0,0.0,INTERIOR,MSH_FLOATING,false) ; 
                    if (nodes[ii] > nodes[kk])
                        side_nodes.Store(
                                   EdgeKey(nodes[ii],nodes[kk]),
                                   elem.nodes[ii+num_nodes]) ;
                    else
                        side_nodes.Store(
                                   EdgeKey(nodes[kk],nodes[ii]),
                                   elem.nodes[ii+num_nodes]) ;
                }
            }
        } else {
            elem.num_nodes = num_nodes ;
            for (int ii=0 ; ii<num_nodes ; ++ii)
                elem.nodes[ii] = nodes[ii] ;
        }

        pelem_table->Store(elem.elem_id,elem) ;
    }

    delete [] elems ;
}




// %(MshQuad2D::DoPrint-void-|-MshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    DoPrint - debug routine to display a mesh 
**
**      void DoPrint(MshTopo2D *msh_topo)
**
**        msh_topo - (in)  mesh topology 
**
**      Description: This method provides debugging support. It prints 
**          information that can be used to display a mesh. 
**
**
** -- */

void MshQuad2D::DoPrint(MshTopo2D *msh_topo,FILE *fd,bool labels)
{
    FILE *lfd = (fd == 0) ? stdout : fd ;

    Queue<MshTopo2D::Edge> *elist ;
    elist = msh_topo->EdgeList() ;
    MshTopo2D::Edge dedge ;
    while (elist->RemoveFromFront(&dedge)) {
        IntNode *p_nd0, *p_nd1 ;
        p_nd0 = pnode_table->Get(dedge.nd0) ;
        p_nd1 = pnode_table->Get(dedge.nd1) ;
//        fprintf(lfd,"# %d %d\n",dedge.nd0,dedge.nd1) ;
//        fprintf(lfd,"l %g %g %g %g\n",
        fprintf(lfd,"e %20.13g %20.13g 0 %20.13g %20.13g 0\n",
               p_nd0->coord[0],p_nd0->coord[1],
               p_nd1->coord[0],p_nd1->coord[1]) ;
        if (labels) {
//            fprintf(lfd,"t %g %g %d\n",
            fprintf(lfd,"t %20.13g %20.13g 0 %d\n",
                   p_nd1->coord[0],p_nd1->coord[1],p_nd1->id) ;
//            fprintf(lfd,"t %g %g %d\n",
            fprintf(lfd,"t %20.13g %20.13g 0 %d\n",
                   p_nd0->coord[0],p_nd0->coord[1],p_nd0->id) ;
        }
    }
    fflush(lfd) ;
    delete elist ;
}

// void MshQuad2D::DisplayTopo(char *label,MshTopo2D *msh_topo)
// {
//     Queue<MshTopo2D::Edge> *elist ;
//     elist = msh_topo->EdgeList() ;
//     MshTopo2D::Edge dedge ;
//     while (elist->RemoveFromFront(&dedge)) {
//         IntNode *p_nd0, *p_nd1 ;
//         p_nd0 = pnode_table->Get(dedge.nd0) ;
//         p_nd1 = pnode_table->Get(dedge.nd1) ;
// //        printf("# %d %d\n",dedge.nd0,dedge.nd1) ;
//         printf("l %g %g %g %g\n",
//                p_nd0->coord[0],p_nd0->coord[1],
//                p_nd1->coord[0],p_nd1->coord[1]) ;
//         if (do_labels) {
//             printf("t %g %g %d\n",
//                    p_nd1->coord[0],p_nd1->coord[1],p_nd1->id) ;
//             printf("t %g %g %d\n",
//                    p_nd0->coord[0],p_nd0->coord[1],p_nd0->id) ;
//         }
//     }
//     printf("a %s\nn\n",label) ;
//     fflush(stdout) ;
//     delete elist ;
// }



// %(MshQuad2D::InitializeQuadGen-void-|-MshTopo2D-|*-MshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    InitializeQuadGen - initialization for quad generation 
**
**      void InitializeQuadGen(
**              MshTopo2D   *msh_topo,
**              MshEdgeList *edge_list)
**
**        msh_topo  - (in)  triangular mesh 
**        edge_list - (i/o) active boundary list 
**
**      Description: This method initializes the adjacent vertex and 
**          edge list structures used for quadrilateral element 
**          generation 
**
**
** -- */

void MshQuad2D::InitializeQuadGen(MshTopo2D *msh_topo,
                                        MshEdgeList *edge_list)
{
    // go through all the elements and build the adjacent
    // vertex table

    Dict<int,MshElement2D>::DictIterator eiter(pelem_table) ;
    for (eiter.First() ; eiter.More() ; ++eiter) {
        MshElement2D& elem = eiter.Entry() ;
        int nn ;
        if (elem.num_nodes == 8)
            nn = 4 ;
        else if (elem.num_nodes == 6)
            nn = 3 ;
        else
            nn = elem.num_nodes ;
        msh_topo->InsertElement(elem.elem_id,nn,elem.nodes) ;
    }

    // now loop through the boundary edges and
    // add them to the edge heap.

    OrderedSet<IntEdge>::SetIterator iter = pedge_table->Iterator() ;
    for (iter.First() ; iter.More() ; ++iter) {
        int cw_id = msh_topo->GetCWBdryNode(
                                 iter.Entry().node_id[0],
                                 iter.Entry().node_id[1]) ;
        int ccw_id = msh_topo->GetCCWBdryNode(
                                 iter.Entry().node_id[0],
                                 iter.Entry().node_id[1]) ;
        edge_list->InsertEdge(iter.Entry().node_id[0],
                              iter.Entry().node_id[1],
                              cw_id,ccw_id) ;
    }
    edge_list->InitializeCodes() ;
}




// %(MshQuad2D::FindElementSide-int-|-MshTopo2D-|*-IntNode-|*-IntNode-|*-IntNode-|*-bool-|)
/* ++ ----------------------------------------------------------
**
**    FindElementSide - find a candidate element side 
**
**      int FindElementSide(
**              MshTopo2D *msh_topo,
**              IntNode    *this_node,
**              IntNode    *prev_node,
**              IntNode    *next_node,
**              bool          prev_edge)
**
**        msh_topo  - (in)  triangular mesh 
**        this_node - (in)  near node on the side 
**        prev_node - (in)  adjacent node on boundary 
**        next_node - (in)  adjacent node on boundary 
**        prev_edge - (in)  flag that tells if we are looking for a 
**                          right or left side 
**
**      Description: Given nodes on the active boundary this method 
**          finds a candidate for a side for a new element. This can be 
**          either an existing edge, an edge generated by "swapping", 
**          or a new edge generated by splitting one of the existing 
**          edges. 
**
**      Return Value: id of the far node on the side 
**
**
** -- */

int MshQuad2D::FindElementSide(MshTopo2D *msh_topo,
                                              IntNode *this_node,
                                              IntNode *prev_node,
                                              IntNode *next_node,
                                              double base_length,
                                              bool prev_edge)
{
    Vec2D normal ;
    double tol = 0.01 ;

    // Find the "ideal" normal, which is the direction of the
    // bisector of the included angle, then turn this into
    // a coordinate point relative to the base node.

    double angle = Angle2Pi(this_node->coord,next_node->coord,
                            prev_node->coord) ;
    if (angle < QUAD_EDGE_NORMAL_TOLERANCE) {
        normal = BisectNorm(this_node->coord,prev_node->coord,
                            next_node->coord) ;
    } else {
        if (prev_edge) {
            normal[0] = prev_node->coord[1] - this_node->coord[1] ;
            normal[1] = this_node->coord[0] - prev_node->coord[0] ;
        } else {
            normal[0] = this_node->coord[1] - next_node->coord[1] ;
            normal[1] = next_node->coord[0] - this_node->coord[0] ;
        }
    }

    normal += this_node->coord ;

    // check to see if there is an existing edge that is
    // close to the normal direction.

    int use_vtx = 0, vtx_before = 0, vtx_after = 0 ;
    IntNode *swap_trial ;

    use_vtx = CheckForEdge(msh_topo,this_node,prev_node,
                           next_node,prev_edge,
                           normal,&vtx_before,&vtx_after) ;

    if (use_vtx != 0) return(use_vtx) ;

    // if we did not find a suitable edge then we look to
    // see if we can do edge swapping.  To do this we
    // start with the "after" vertex and look through it's
    // adjacent verticies for the "before" vertex.  Once
    // we find this, the next adjacent vertex is the 
    // candidate.

    TopoAdjVtxIterator iter0(msh_topo,this_node->id) ; 
    int elem_1, elem_2 ;

    while (iter0.AdjVtx() != vtx_before) ++iter0 ;
    elem_1 = iter0.CcwElem() ;

    TopoAdjVtxCyclicIterator iter1(msh_topo,vtx_after) ; 
    while (1) {
        if (iter1.AdjVtx() == vtx_before) break ;
        ++iter1 ;
    }

    elem_2 = iter1.CcwElem() ;
    ++iter1 ;

    // if elem_2 is NO_ELEM, then we have the case that the
    // edge joining the before and after verticies is part
    // of the current boundary.  The best we can do is to
    // return either the before or after vertex depending
    // on which makes a smaller angle with the normal

    if (elem_2 == NO_ELEM) {
        IntNode *before = pnode_table->Get(vtx_before) ;
        IntNode *after = pnode_table->Get(vtx_after) ;
        if (Angle(this_node->coord,before->coord,normal) <
            Angle(this_node->coord,normal,after->coord))
            return(vtx_before) ;
        else
            return(vtx_after) ;
    }

    // now we have a candidate, check the angle

    swap_trial = pnode_table->Get(iter1.AdjVtx()) ;
    angle = Angle(this_node->coord,swap_trial->coord,normal) ;

    if (fabs(angle) < QUAD_NEAR_ANGLE_TOLERANCE) {

        // check to make sure that the resulting edge
        // is not too long

        double tlen = (this_node->coord-swap_trial->coord).Magnitude() ;
        double e1len = (this_node->coord-prev_node->coord).Magnitude() ;
        double e2len = (this_node->coord-next_node->coord).Magnitude() ;

        if (tlen < 0.866025*(e1len + e2len)) {
            int thisid = this_node->id ;
            int before = vtx_before ;
            int after = vtx_after ;
            int swap = swap_trial->id ;
            int rtn_vtx = iter1.AdjVtx() ;

            // fix up the adjacent vertex list to reflect
            // the swap

            msh_topo->DeleteTriangle(thisid,before,after) ;
            msh_topo->DeleteTriangle(swap,after,before) ;

            msh_topo->InsertTriangle(elem_1,thisid,before,swap) ;
            msh_topo->InsertTriangle(elem_2,thisid,swap,after) ;

            return(rtn_vtx) ;
        }
    }
    
    // if we could not swap an edge, then we need to split an edge.
    // find the intersection of the "ideal" vector and the
    // line between the before and after verticies

    IntNode *before = pnode_table->Get(vtx_before) ;
    IntNode *after = pnode_table->Get(vtx_after) ;
    Vec2D intsc = IntersectLines(this_node->coord,normal,
                                       before->coord,after->coord) ;

    // reject if the resulting edge will be too small

    double dist = (intsc-this_node->coord).Magnitude() ;
    if (dist < tol*base_length) return(-1) ;

    // define the new node

    int new_id = NewNode(intsc.x(),intsc.y(),
                         INTERIOR,MSH_FLOATING,true) ;

    // update the vertex adjacency tables

    int thisid = this_node->id ;
    int beforeid = vtx_before ;
    int afterid = vtx_after ;
    int swapid = swap_trial->id ;

    int elem_3 = NewElemNum() ;
    int elem_4 = NewElemNum() ;

    msh_topo->DeleteTriangle(thisid,beforeid,afterid) ;
    msh_topo->DeleteTriangle(swapid,afterid,beforeid) ;

    msh_topo->InsertTriangle(elem_1,new_id,thisid,beforeid) ;
    msh_topo->InsertTriangle(elem_2,new_id,beforeid,swapid) ;
    msh_topo->InsertTriangle(elem_3,new_id,swapid,afterid) ;
    msh_topo->InsertTriangle(elem_4,new_id,afterid,thisid) ;
    return(new_id) ;
}    




// %(MshQuad2D::CheckForEdge-int-|^const-MshTopo2D-|*-IntNode-|*-IntNode-|*-IntNode-|*-bool-|-Vec2D-|-int-|*-int-|*)
/* ++ ----------------------------------------------------------
**
**    CheckForEdge - find an edge in the search direction 
**
**      int CheckForEdge(
**              MshTopo2D *msh_topo,
**              IntNode    *this_node,
**              IntNode    *prev_node,
**              IntNode    *next_node,
**              bool          prev_edge,
**              Vec2D   normal,
**              int           *vtx_before,
**              int           *vtx_after) const
**
**        msh_topo   - (in)  triangular mesh 
**        this_node  - (in)  near node on this side 
**        prev_node  - (in)  adjacent node on boundary 
**        next_node  - (in)  adjacent node on boundary 
**        prev_edge  - (in)  flag that tells if we are looking for a 
**                           right or left side 
**        normal     - (in)  search normal direction 
**        vtx_before - (out) bracketing vertex 
**        vtx_after  - (out) bracketing vertex 
**
**      Description: This method looks at the nodes adjacent to the 
**          input node to see if there is one close to the normal 
**          direction. If this is the case, the node is returned. If 
**          not, the routine returns the nodes that bracket the search 
**          direction. 
**
**      Return Value: non-zero is the far node on a found edge. Zero 
**          indicates that the bracketing edges are being returned. 
**
**
** -- */

int MshQuad2D::CheckForEdge(
                        MshTopo2D *msh_topo,
                        IntNode *this_node,
                        IntNode *prev_node,
                        IntNode *next_node,
                        bool prev_edge,
                        Vec2D normal,
                        int *vtx_before,
                        int *vtx_after) const
{
    // this routine looks at the nodes adjacent to the input
    // node to see if there is one close to the normal direction.
    // If this is the case, the node is returned.  If not, the
    // routine returns the nodes that bracket the search direction.

    // first, however, look around the node and determine if
    // it is adjacent to only one NO_ELEMENT region (the normal
    // case), or if it is adjacent to two or more.

    int use_vtx = 0 ;
    int vtx_before_next = -1 ;
    int num_free = 0 ;

    *vtx_before = *vtx_after = 0 ;

    TopoAdjVtxIterator iter(msh_topo,this_node->id) ;
    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.CcwElem() == NO_ELEM) ++num_free ;
    }

    if (num_free == 1) {
        double angle_min = 6.2831 ;
        double angle_before = 6.2831 ;
        double angle_after = 6.2831 ;
        *vtx_after = 0 ;
        for (iter.First() ; iter.More() ; ++iter) {
            IntNode *trial = pnode_table->Get(iter.AdjVtx()) ;
            double angle = Angle(this_node->coord,trial->coord,normal) ;

            if ((fabs(angle) < QUAD_NEAR_ANGLE_TOLERANCE) &&
                (fabs(angle) < angle_min)) {
                angle_min = fabs(angle) ;
                use_vtx = iter.AdjVtx() ;
            }

            if (angle > 0.0) {
                if ((*vtx_before != 0) &&
                    (vtx_before_next == iter.AdjVtx())) {
                    *vtx_before = iter.AdjVtx() ;
                    angle_before = fabs(angle) ;
                } else {
                    if (fabs(angle) <= angle_before) {
                        *vtx_before = iter.AdjVtx() ;
                        angle_before = fabs(angle) ;
                        TopoAdjVtxIterator biter = iter ;
                        ++biter ;
                        vtx_before_next = biter.AdjVtx() ;
                    }
                }
            } else {
                if (fabs(angle) < angle_after) {
                    *vtx_after = iter.AdjVtx() ;
                    angle_after = fabs(angle) ;
                    if ((*vtx_before != 0) &&
                        (vtx_before_next == *vtx_after)) break ;
                }
            }
        }
        return(use_vtx) ;

        // here deal with the nasty special case of more than
        // one NO_ELEMENT region about a vertex

    } else {

        TopoAdjVtxCyclicIterator citer(msh_topo,this_node->id) ;
        int start_vtx = 0, stop_vtx = 0 ;

        if (prev_edge) {

            // for this case, the base edge is from the "prev_node"
            // to "this_node".  Move ccw from the next node
            // to find the "start_vtx".  This is the first vertex
            // we would find adjacent to a NO_ELEM region if we
            // moved clockwise from the previous node

            // search until we find the "prev_node" about the vertex

            while (citer.AdjVtx() != prev_node->id) ++citer ;

            // search from this until we get back to the prev
            // vertex

            ++citer ;
            bool in_void = false ;

            while (citer.AdjVtx() != prev_node->id) {
                if (citer.CcwElem() == NO_ELEM) {
                    in_void = true ;
                } else if (in_void) {
                    start_vtx = citer.AdjVtx() ;
                    in_void = false ;
                }

                ++citer ;
            }

            // now search through the region and see if
            // if the normal direction is inside

            while (citer.AdjVtx() != start_vtx) ++citer ;
            double angle_min = 6.2831 ;
            double angle_before = 6.2831 ;
            double angle_after = 6.2831 ;
            while (citer.AdjVtx() != stop_vtx) {
                IntNode *trial = pnode_table->Get(citer.AdjVtx()) ;
                double angle = Angle(this_node->coord,trial->coord,normal) ;

                if ((fabs(angle) < QUAD_NEAR_ANGLE_TOLERANCE) &&
                    (fabs(angle) < angle_min)) {
                    angle_min = fabs(angle) ;
                    use_vtx = citer.AdjVtx() ;
                }

                if (angle > 0.0) {
                    if (fabs(angle) <= angle_before) {
                        *vtx_before = citer.AdjVtx() ;
                        angle_before = fabs(angle) ;
                        TopoAdjVtxCyclicIterator biter = citer ;
                        ++biter ;
                        vtx_before_next = biter.AdjVtx() ;
                    }
                } else {
                    if (fabs(angle) < angle_after) {
                        *vtx_after = citer.AdjVtx() ;
                        angle_after = fabs(angle) ;
                        if ((*vtx_before != 0) &&
                            (vtx_before_next == *vtx_after)) break ;
                    }
                }
                ++citer ;
            }
            if (use_vtx != 0) {
                return(use_vtx) ;
            } else if (*vtx_before != 0) {
                if (*vtx_after == 0) *vtx_after = stop_vtx ;
                return(0) ;
            } else {
                return(start_vtx) ;
            }

        } else {

            // for this case, the base edge is from "this_node"
            // to the "next_node".  Move ccw from the next node
            // to find the first node with a ccw_element pointer
            // to NO_ELEM.  This brackests the region in which
            // to search for valid nodes

            // search until we find the "next_node" about the vertex

            citer.First() ;
            while (citer.AdjVtx() != next_node->id) ++citer ;
            start_vtx = citer.AdjVtx() ;
   
            // continue around the vertex

            while (citer.CcwElem() != NO_ELEM) ++citer ;
            stop_vtx = citer.AdjVtx() ;

            // now search through the region and see if
            // if the normal direction is inside

            while (citer.AdjVtx() != start_vtx) ++citer ;
            double angle_min = 6.2831 ;
            double angle_before = 6.2831 ;
            double angle_after = 6.2831 ;
            while (citer.AdjVtx() != stop_vtx) {
                IntNode *trial = pnode_table->Get(citer.AdjVtx()) ;
                double angle = Angle(this_node->coord,trial->coord,normal) ;

                if ((fabs(angle) < QUAD_NEAR_ANGLE_TOLERANCE) &&
                    (fabs(angle) < angle_min)) {
                    angle_min = fabs(angle) ;
                    use_vtx = citer.AdjVtx() ;
                }

                if (angle > 0.0) {
                    if (fabs(angle) <= angle_before) {
                        *vtx_before = citer.AdjVtx() ;
                        angle_before = fabs(angle) ;
                        TopoAdjVtxCyclicIterator biter = citer ;
                        ++biter ;
                        vtx_before_next = biter.AdjVtx() ;
                    }
                } else {
                    if (fabs(angle) < angle_after) {
                        *vtx_after = citer.AdjVtx() ;
                        angle_after = fabs(angle) ;
                        if ((*vtx_before != 0) &&
                            (vtx_before_next == *vtx_after)) break ;
                    }
                }
                ++citer ;
            }
            if (use_vtx != 0) {
                return(use_vtx) ;
            } else if (*vtx_after != 0) {
                return(0) ;
            } else {
                return(stop_vtx) ;
            }
        }
    }
}


bool MshQuad2D::TweakNode(IntNode *nd,MshTopo2D *msh_topo,
                      List<int> &moved,
                      List<Vec2D> &orig)
{
    // if this node is not movable we cannot do anything

    if (nd->motion == MSH_FIXED) return(false) ;

    // loop through the adjacent nodes to compute a tolerance

    int num = 0 ;
    double tol = 0.0 ;
    TopoAdjVtxIterator iter(msh_topo,nd->id) ;
    for (iter.First() ; iter.More() ; ++iter) {
        int adj_id = iter.AdjVtx() ;
        IntNode *adj = pnode_table->Get(adj_id) ;
        tol += (nd->coord - adj->coord).Magnitude() ;
        ++num ;
    }

    tol = 0.000001 * tol/num ;

    // find a feasable region where we can place the point

    FeasableRegion feasable(tol) ;

    for (iter.First() ; iter.More() ; ++iter) {

        int elem_id = iter.CcwElem() ;
        MshTopo2D::Edge edge ;

        edge = msh_topo->OppositeEdge(elem_id,nd->id) ;
        IntNode *nd0 = pnode_table->Get(edge.nd0) ;
        IntNode *nd1 = pnode_table->Get(edge.nd1) ;

        feasable.AddConstraint(nd0->coord,
                    nd1->coord - nd0->coord) ;
    }

    // get the corners of the feasable region and find
    // the centroid

    List<Vec2D> *poly = feasable.GetVertices() ;

    // find the average of these points

    if (poly->Len() != 0) {
        int j, k ;
        Vec2D center ;
        for (j=0 ; j<poly->Len() ; ++j) {
            center += poly->At(j) ;
        }
        center /= poly->Len() ;

        // find the centroid

        double sum_area = 0.0 ;
        double sum_prod_x = 0.0 ;
        double sum_prod_y = 0.0 ;

        for (j=0 ; j<poly->Len() ; ++j) {
            k = (j+1) % poly->Len() ;
            Vec2D tcen = (poly->At(j)+poly->At(k)+center)/3 ;
            double area = 0.5 * fabs(
                poly->At(j).x()*poly->At(k).y() - 
                poly->At(j).y()*poly->At(k).x() +
                poly->At(k).x()*center.y() -
                poly->At(k).y()*center.x() +
                center.x()*poly->At(j).y() -
                center.y()*poly->At(j).x()) ;

            sum_area += area ;
            sum_prod_x += area * (tcen.x()-center.x()) ;
            sum_prod_y += area * (tcen.y()-center.y()) ;
        } 
        delete poly ;
        
        moved.Append(nd->id) ;
        orig.Append(nd->coord) ;
        nd->coord = Vec2D(center.x() + sum_prod_x/sum_area,
                                center.y() + sum_prod_y/sum_area) ;
        return (true) ;
    }
    return (false) ;
}


#if 0
bool MshQuad2D::TweakNode(IntNode *nd,MshTopo2D *msh_topo,
                      List<int> &moved,
                      List<Vec2D> &orig)
{
    // if this node is not movable we cannot do anything

    if (nd->motion == MSH_FIXED) return(false) ;

    // find the opposite edges on adjacent triangles and find
    // the closest

    bool first = true ;
    double min_dist = 0.0 ;
    Vec2D move_dir ;

    TopoAdjVtxIterator iter(msh_topo,nd->id) ;
    for (iter.First() ; iter.More() ; ++iter) {

        int elem_id = iter.CcwElem() ;
        MshTopo2D::Edge edge ;

        edge = msh_topo->OppositeEdge(elem_id,nd->id) ;
        IntNode *nd0 = pnode_table->Get(edge.nd0) ;
        IntNode *nd1 = pnode_table->Get(edge.nd1) ;

        double blen ;
        double dist = DistSqr(nd->coord,
                              nd0->coord,nd1->coord,&blen) ;

        if (first || (dist<min_dist)) {
            first = false ;
            min_dist = dist ;
            Vec2D delta = nd1->coord - nd0->coord ;
            move_dir = Vec2D(-delta.y(),delta.x()) ;
        }
    }

    // move the node one half the distance to the closest
    // edge in the opposite direction

    moved.Append(nd->id) ;
    orig.Append(nd->coord) ;

    nd->coord += (0.5*sqrt(min_dist)) * move_dir.Normalize() ;
    return(true) ;
}
#endif


// %(MshQuad2D::FindCrossedEdges-CArbQueue-|<MshTopo2D::Edge>*^const-int-|-int-|-MshTopo2D-|*-MshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    FindCrossedEdges - find crossed edges 
**
**      CArbQueue <MshTopo2D::Edge>*FindCrossedEdges(
**              int             start_id,
**              int             stop_id,
**              MshTopo2D   *msh_topo,
**              MshEdgeList *edge_list) const
**
**        start_id  - (in)  start node id 
**        stop_id   - (in)  end node id 
**        msh_topo  - (in)  triangular mesh 
**        edge_list - (in)  active boundary 
**
**      Description: This method returns a list of triangle edges 
**          crossed by a vector drawn between two nodes 
**
**      Return Value: a ArbQueue containing a list of the edges crossed 
**
**
** -- */

Queue<MshTopo2D::Edge>
    *MshQuad2D::FindCrossedEdges(int start_id,
                      int stop_id,
                      MshTopo2D *msh_topo,
                      MshEdgeList *bdry_list,
                      List<int> &moved,
                      List<Vec2D> &orig)
{
    Queue<MshTopo2D::Edge> *edge_list = 
        new Queue<MshTopo2D::Edge>() ;

    moved.Clear() ;
    orig.Clear() ;
   
    IntNode *nd_start = pnode_table->Get(start_id) ;
    IntNode *nd_stop  = pnode_table->Get(stop_id) ;

    // look around the start node and find the two edges that
    // are before and after the search vector (saved in adj_0
    // and adj_1).

    IntNode *adj_1 ;
    double cross_0 = 0, cross_1 = 0 ;
    int elem_id = 0, elem_id_1 = 0 ;
    bool found = false ;
    bool first_trip = true ;
    double tol = 1.0e-4 * (nd_start->coord-nd_stop->coord).Magnitude() ; 

    TopoAdjVtxIterator iter(msh_topo,start_id) ;
    if (!iter.More()) return(edge_list) ;

    // loop through and check to see if we are adjacent to the
    // stop id

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.AdjVtx() == stop_id) return(edge_list) ;
    }

    // otherwise do a more detailed search

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.AdjVtx() == stop_id) return(edge_list) ;
        adj_1 = pnode_table->Get(iter.AdjVtx()) ;
        elem_id_1 = iter.CcwElem() ;
        cross_1 = CrossProd(nd_start->coord,adj_1->coord,nd_stop->coord) ;

        if ((fabs(cross_1) < tol) && 
            (fabs(Angle(nd_start->coord,adj_1->coord,nd_stop->coord)) < 1.57) &&
            (adj_1->motion != MSH_FIXED)) {

            // if we get here then there is an edge that is along
            // the search direction.  The strategy is to move this
            // node off the search line.

            if (!TweakNode(adj_1,msh_topo,moved,orig)) {
                delete edge_list ;
                return(0) ;
            }

            // first loop around the node and find the length of
            // the shortest adjacent edge

            // double min_length = -1.0 ;
            // TopoAdjVtxIterator miter(msh_topo,adj_1->id) ;
            // for (miter.First() ; miter.More() ; ++miter) {
            //     IntNode *adj_t = pnode_table->Get(miter.AdjVtx()) ;
            //     double length = (adj_1->coord - adj_t->coord).Magnitude() ;
            //     if ((min_length < 0.0) ||
            //         (length < min_length)) min_length = length ;
            // }

            // // first loop around the edge and find edge adjacent
            // // to the start node

            // TopoAdjVtxCyclicIterator citer(msh_topo,adj_1->id) ;
            // while (1) {
            //     if (citer.AdjVtx() == start_id) break ;
            //     ++citer ;
            // }
            // ++citer ;
            // IntNode *adj_t = pnode_table->Get(citer.AdjVtx()) ;
            // Vec2D dir = (adj_t->coord - adj_1->coord).Normalize() ;
            // Vec2D mid = adj_1->coord + 0.5*min_length*dir ;

            // check to see if the new node position will cross
            // the current advancing front

            // if (bdry_list->DoesThisCrossBoundary(nd_start->coord,mid)) {
            //     delete edge_list ;
            //     return(0) ;
            // }

            // moved.Append(iter.AdjVtx()) ;
            // orig.Append(adj_1->coord) ;

            // adj_1->coord = mid ;
            cross_1 = CrossProd(nd_start->coord,adj_1->coord,
                                nd_stop->coord) ;
        }

        if (first_trip || (cross_1 > tol)) {
            first_trip = false ;
            cross_0 = cross_1 ;
            elem_id = elem_id_1 ;
        } else {
            if ((cross_0 > tol) && (cross_1 < -tol)) {
                found = true ;
                break ;
            }
        }
    }
    if (!found && elem_id_1 == NO_ELEM) {
        delete edge_list ;
        return(0) ;
    }

    // find the edge on this element opposite the start node

    MshTopo2D::Edge edge, next_edge ;
    edge = msh_topo->OppositeEdge(elem_id,start_id) ;

    // if the edge is not on the front, then add it to the
    // list of crossed edges

    if (bdry_list->ContainsEdge(edge.nd1,edge.nd0)) {
        delete edge_list ;
        return(0) ;
    }
    edge_list->AddToBack(edge) ;

    // Now that we have the first edge, we continue to move
    // towards the stop vertex

    while (1) {

        // find the element on the opposite side of the edge

        if (elem_id == edge.elem0)
            elem_id = edge.elem1 ;
        else
            elem_id = edge.elem0 ;

        // check for a boundary

        if (elem_id == NO_ELEM) {
            delete edge_list ;
            return(0) ;
        }

        // check to see if the stop node in adjacent to this
        // element.  If so we are done.

        if (msh_topo->ElemHasNode(elem_id,stop_id)) break ;

        // get the node opposite the current edge on the
        // element and find the vector from start node
        // to this node

        MshTopo2D::Edge rev_edge ;
        rev_edge.nd0 = edge.nd1 ;
        rev_edge.nd1 = edge.nd0 ;
        rev_edge.elem0 = edge.elem1 ;
        rev_edge.elem1 = edge.elem0 ;
        int opp_id = msh_topo->OppositeNode(elem_id,&rev_edge) ;
        IntNode *nd_opp = pnode_table->Get(opp_id) ;

        // check to see if the opposite node falls on the search
        // path.  If so move it off

        double blen ;
        double dist = sqrt(DistSqr(nd_opp->coord,
                        nd_start->coord,nd_stop->coord,
                        &blen)) ;

        if (dist < tol*100.0) {
            if (!TweakNode(nd_opp,msh_topo,moved,orig)) {
                delete edge_list ;
                return(0) ;
            }
        }

        // look at the cross product of the current vector
        // and the search vector to determine which of the
        // triangle's edges should be added to the list 

        if (CrossProd(nd_start->coord,nd_stop->coord,nd_opp->coord) > 0.0)
            next_edge = msh_topo->NextCCWEdge(elem_id,&rev_edge) ;
        else
            next_edge = msh_topo->NextCWEdge(elem_id,&rev_edge) ;

        edge = next_edge ;

        // if the edge is not on the front, then add it to the
        // list of crossed edges

        if (bdry_list->ContainsEdge(edge.nd1,edge.nd0)) {
            delete edge_list ;
            return(0) ;
        }
        edge_list->AddToBack(edge) ;
    }

    return(edge_list) ;
}




// %(MshQuad2D::RecoverTopEdge-bool-|^const-int-const|-int-const|-MshTopo2D-|*-MshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    RecoverTopEdge - do edge swapping 
**
**      bool RecoverTopEdge(
**              const int       start_id,
**              const int       stop_id,
**              MshTopo2D   *msh_topo,
**              MshEdgeList *edge_list) const
**
**        start_id  - (in)  start node id 
**        stop_id   - (in)  end node id 
**        msh_topo  - (i/o) triangular mesh 
**        edge_list - (in)  active boundary 
**
**      Description: This method performs "edge swapping" to recover 
**          the top edge of a new element. 
**
**      Return Value: true if succesful 
**
**
** -- */

// set a limit on the number of edges that can be crossed when
// trying to recover a top edge.  The limit is because if
// too many edges are crossed then we are transitioning from
// big elements to very small elements and it's likely to cause
// problems if this element is accepted.  Also, with lots of
// relatively small edges we may run into tolerancing problems.

#define MAX_CROSSED_EDGES 20

bool MshQuad2D::RecoverTopEdge(
             int start_id,
             int stop_id,
             MshTopo2D *msh_topo,
             MshEdgeList *edge_list)
{
    // first find the coordinates of the start and stop nodes

    IntNode *nd_start = pnode_table->Get(start_id) ;
    IntNode *nd_stop  = pnode_table->Get(stop_id) ;

    // now find get a list of edges crossing the vector between
    // the start and stop nodes

    Queue<MshTopo2D::Edge> *crossed ;
    List<int> moved ;
    List<Vec2D> orig ;

    crossed = FindCrossedEdges(start_id,stop_id,
                               msh_topo,edge_list,moved,orig) ;
    if (crossed == 0) return(false) ;
    if (crossed->Len() > MAX_CROSSED_EDGES) {
       delete crossed ;
       return(false) ;
    }

    int max_iters = crossed->Len() <= 10 ? 100 :
                    crossed->Len() * crossed->Len() ;

    MshTopo2D::Edge edge ;
    int saved_nd[2] ;
    bool saved = false ;
    int iter = 0 ;

//     if ((int)debug_num >= print_num) {
//         DoPrint(msh_topo) ;
//         //DoPrint(&quad_topo) ;
//         printf("n\n") ;
//     }

    // this is a list used to store the id's of the swapped
    // triangle nodes.  If we discover that the edge cannot
    // be recovered we run backwards through the list to
    // restore the original topology.

    List<int> swap_ids ;

    // loop while there are edges

    while (crossed->RemoveFromFront(&edge)) {

        iter++ ;
        if ((saved && (edge.nd0 == saved_nd[0]) &&
                      (edge.nd1 == saved_nd[1])) ||
            (iter > max_iters)) {

            // unwind the changes that have been made to
            // the mesh
            int cur = swap_ids.Len() - 6 ;
            while (cur >= 0) {
                msh_topo->DeleteTriangle(swap_ids[cur],
                                         swap_ids[cur+3],
                                         swap_ids[cur+2]) ;
                msh_topo->DeleteTriangle(swap_ids[cur+3],
                                         swap_ids[cur+1],
                                         swap_ids[cur+2]) ;
                msh_topo->InsertTriangle(swap_ids[cur+4],
                                         swap_ids[cur],
                                         swap_ids[cur+1],
                                         swap_ids[cur+2]) ;
                msh_topo->InsertTriangle(swap_ids[cur+5],
                                         swap_ids[cur],
                                         swap_ids[cur+3],
                                         swap_ids[cur+1]) ;
                cur -= 6 ;
            }

            for (int i=moved.Len()-1 ; i >= 0 ; --i) {
                IntNode *nd = pnode_table->Get(moved[i]) ;
                nd->coord = orig[i] ;
            }

            delete crossed ;
            return(false) ;
        }

        /*
           With reference to the figures below, the crossed
           edge is nd_0 to nd_1.  Check to see if the triangles
           in the swapped configuration have positive area.
           If so, do the swap.  If the edge nd_2 to nd_3 in
           the swapped configuration intersects the the line
           segment from the start to the end node, then add
           it to the crossed edge list.


               initial config     swapped config
    
                     * nd_2         nd_2 *
                    / \                 /|\
                   /   \               / | \
             nd_0 *-----* nd_1   nd_0 *  |  * nd_1
                   \   /               \ | /
                    \ /                 \|/
                     * nd_3         nd_3 *

        */

        int id_0 = edge.nd0 ;
        int id_1 = edge.nd1 ;

        IntNode *nd_0 = pnode_table->Get(id_0) ;
        IntNode *nd_1 = pnode_table->Get(id_1) ;

        edge.elem0 = msh_topo->GetCCWElem(id_0,id_1) ;
        edge.elem1 = msh_topo->GetCCWElem(id_1,id_0) ;

        int id_2 = msh_topo->OppositeNode(edge.elem0,&edge) ;
        IntNode *nd_2 = pnode_table->Get(id_2) ;

        MshTopo2D::Edge rev_edge ;
        rev_edge.nd0 = edge.nd1 ;
        rev_edge.nd1 = edge.nd0 ;
        rev_edge.elem0 = edge.elem1 ;
        rev_edge.elem1 = edge.elem0 ;
        int id_3 = msh_topo->OppositeNode(edge.elem1,&rev_edge) ;
        IntNode *nd_3 = pnode_table->Get(id_3) ;

        double base = (nd_3->coord - nd_2->coord).Magnitude() ;
        const double tol = 0.005 ;
        double check = tol * base * base ;

        // check the areas

        if ((Area(nd_0->coord,nd_3->coord,nd_2->coord) > check) &&
            (Area(nd_1->coord,nd_2->coord,nd_3->coord) > check)) {

            int eid_0 = edge.elem0 ;
            int eid_1 = edge.elem1 ;

            msh_topo->DeleteTriangle(id_0,id_1,id_2) ;
            msh_topo->DeleteTriangle(id_0,id_3,id_1) ;

            msh_topo->InsertTriangle(eid_0,id_0,id_3,id_2) ;
            msh_topo->InsertTriangle(eid_1,id_3,id_1,id_2) ;

//        if ((int)debug_num >= print_num) {
//            DoPrint(msh_topo) ;
//            //DoPrint(&quad_topo) ;
//            printf("n\n") ;
//        }

            swap_ids.Append(id_0) ;
            swap_ids.Append(id_1) ;
            swap_ids.Append(id_2) ;
            swap_ids.Append(id_3) ;
            swap_ids.Append(eid_0) ;
            swap_ids.Append(eid_1) ;

            // check to see if the new edge from nd_2 to nd_3
            // crosses the vector from the start to end node.
            // if so, add it to the edge list

            if (Cross(nd_start->coord,nd_stop->coord,
                      nd_2->coord,nd_3->coord)) {
                edge.nd0 = id_2 ;
                edge.nd1 = id_3 ;
                edge.elem0 = eid_0 ;
                edge.elem1 = eid_1 ;
                crossed->AddToBack(edge) ;
            }
            saved = false ;
        } else {

            // else push this edge on the back of the list

            crossed->AddToBack(edge) ;
            if (!saved) {
                saved = true ;
                saved_nd[0] = edge.nd0 ;
                saved_nd[1] = edge.nd1 ;
            }
        }
    }

    // when we get here do a sanity check to make sure that
    // things worked properly.  I have seen one case where
    // tolerancing problems allowed the code to get here
    // with no edge between the start and end verts.

    bool found = false ;
    TopoAdjVtxIterator viter(msh_topo,start_id) ;
    for (viter.First() ; viter.More() ; ++viter) {
        if (viter.AdjVtx() == stop_id) {
            found = true ;
            break ;
        }
    }

    delete crossed ;
    return(found) ;
}




// %(MshQuad2D::ClearRegion-void-|-int-const|-int-const|*-MshTopo2D-|*-MshTopo2D::Edge-const|*)
/* ++ ----------------------------------------------------------
**
**    ClearRegion - clear a quadrilateral region 
**
**      void ClearRegion(
**              const int                    num_corners,
**              const int                    *boundary,
**              MshTopo2D                *msh_topo,
**              const MshTopo2D::Edge *edge)
**
**        num_corners - (in)  3 or 4 
**        boundary    - (in)  corner node id's 
**        msh_topo    - (i/o) triangle mesh 
**        edge        - (in)  one boundary edge 
**
**      Description: Once the four corners of a quadrilateral or a 
**          triangular region have been identified, this routine 
**          deletes all the triangles that currently cover the region. 
**
**
** -- */

void MshQuad2D::ClearRegion(const int num_corners,
                                  const int *boundary,
                                  MshTopo2D *msh_topo,
                                  const MshTopo2D::Edge *edge)
{
    int i ;
    bool bound ;

    // store the next cw and ccw edges

    MshTopo2D::Edge cw_edge = 
         msh_topo->NextCWEdge(edge->elem0,edge) ;
    MshTopo2D::Edge ccw_edge = 
         msh_topo->NextCCWEdge(edge->elem0,edge) ;

    // delete this element from the topology

    msh_topo->DeleteTriangle(edge->nd0,edge->nd1,ccw_edge.nd1) ;

    // check to make sure that the ccw edge is not part of the
    // the boundary

    bound = false ;
    for (i=0 ; i<num_corners ; ++i) {
        if ((ccw_edge.nd0 == boundary[i]) &&
            (ccw_edge.nd1 == boundary[(i+1)%num_corners])) {
            bound = true ;
            break ;
        }
    }

    // check to make sure that there is an element on the
    // other side of the edge, and if so build the reverse
    // of the edge and call ourself recursively

    if (!bound) {
        if (msh_topo->GetCWElem(ccw_edge.nd0,ccw_edge.nd1) != NO_ELEM) {
            MshTopo2D::Edge rev_edge ;
            rev_edge.nd0 = ccw_edge.nd1 ; 
            rev_edge.nd1 = ccw_edge.nd0 ; 
            rev_edge.elem0 = ccw_edge.elem1 ; 
            rev_edge.elem1 = ccw_edge.elem0 ; 
            ClearRegion(num_corners,boundary,msh_topo,&rev_edge) ;
        }
    }

    // check to make sure that the cw edge is not part of the
    // the boundary

    bound = false ;
    for (i=0 ; i<num_corners ; ++i) {
        if ((cw_edge.nd0 == boundary[i]) &&
            (cw_edge.nd1 == boundary[(i+1)%num_corners])) {
            bound = true ;
            break ;
        }
    }

    // check to make sure that there is an element on the
    // other side of the edge, and if so build the reverse
    // of the edge and call ourself recursively

    if (!bound) {
        if (msh_topo->GetCWElem(cw_edge.nd0,cw_edge.nd1) != NO_ELEM) {
            MshTopo2D::Edge rev_edge ;
            rev_edge.nd0 = cw_edge.nd1 ; 
            rev_edge.nd1 = cw_edge.nd0 ; 
            rev_edge.elem0 = cw_edge.elem1 ; 
            rev_edge.elem1 = cw_edge.elem0 ; 
            ClearRegion(num_corners,boundary,msh_topo,&rev_edge) ;
        }
    }
}




// %(MshQuad2D::SmoothFront-void-|-int-const|-int-const|-int-const|*-MshTopo2D-|*-MshTopo2D-|*-MshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    SmoothFront - smooth advancing front nodes 
**
**      void SmoothFront(
**              const int       qcase,
**              const int       num_frt_nodes,
**              const int       *frt_nodes,
**              MshTopo2D   *msh_topo,
**              MshTopo2D   *quad_topo,
**              MshEdgeList *edge_list)
**
**        qcase         - (in)  new quad case 
**        num_frt_nodes - (in)  number of front nodes 
**        frt_nodes     - (in)  list of front nodes 
**        msh_topo      - (in)  triangle mesh 
**        quad_topo     - (in)  quadrilateral mesh 
**        edge_list     - (in)  active boundary 
**
**      Description: Perform smoothing of advancing front nodes 
**
**
** -- */

Vec2D MshQuad2D::LaplaceSmoothFrnt(const int vtx,
                                        MshTopo2D *msh_topo,
                                        MshTopo2D *quad_topo)
{
    int num_adj = quad_topo->NumAdjElems(vtx) ;
    IntNode *node = pnode_table->Get(vtx) ;
    IntNode *adj ;
    Vec2D sum ;
    int num = 0 ;

    // from the mesh topology get the coordinates
    // of the adjacent nodes

    TopoAdjVtxIterator iter0(msh_topo,vtx) ;
    int last = 0 ;

    if (num_adj < 3) {
        for (iter0.First() ; iter0.More() ; ++iter0) {
             adj = pnode_table->Get(iter0.AdjVtx()) ;
             sum += adj->coord - node->coord ;
             last = iter0.AdjVtx() ;
             ++num ;
        }
    }

    // now loop through the quad topology, ignoring
    // the first and last nodes, which should have
    // been dealt with in the triangle mesh

    TopoAdjVtxIterator iter1(quad_topo,vtx) ;

    for (iter1.First() ; iter1.More() ; ++iter1) {
        if (iter1.AdjVtx() == last) break ;
        adj = pnode_table->Get(iter1.AdjVtx()) ;
        sum += adj->coord - node->coord ;
        ++num ;
    }

    if (num == 0) num = 1 ;
    return(sum/num) ;
}




// %(MshQuad2D::IsoParametricSmooth-Vec2D-|-IntNode-|*-MshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    IsoParametricSmooth - smooth one advancing front node 
**
**      Vec2D IsoParametricSmooth(
**              IntNode    *node,
**              MshTopo2D *quad_topo)
**
**        node      - (in)  id of node to smooth 
**        quad_topo - (in)  quadrilateral mesh 
**
**      Description: Performs smoothing of one advancing front node 
**          using an isoparametric algorithm. 
**
**      Return Value: new nodal coordinates 
**
**
** -- */

Vec2D MshQuad2D::IsoParametricSmooth(IntNode *node,
                                          MshTopo2D *quad_topo)
{
    TopoAdjVtxIterator iter(quad_topo,node->id) ;
    int num_nodes, nodes[4] ;
    Vec2D sum(0.0,0.0) ;
    int num_adj = 0 ;
    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.CcwElem() != NO_ELEM) {
            quad_topo->GetElemNodes(iter.CcwElem(),node->id,
                                    &num_nodes,nodes) ;

            if (num_nodes > 4) throw QuadConversionError() ;

            sum += pnode_table->Get(nodes[1])->coord ;
            sum -= pnode_table->Get(nodes[2])->coord ;
            sum += pnode_table->Get(nodes[3])->coord ;
            ++num_adj ;
        }
    }
    Vec2D new_pos = sum / num_adj ;
    return(new_pos-node->coord) ;
}




// %(MshQuad2D::LaplaceSmooth-Vec2D-|-IntNode-|*-MshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    LaplaceSmooth - smooth one internal node 
**
**      Vec2D LaplaceSmooth(
**              IntNode    *node,
**              MshTopo2D *topo)
**
**        node - (in)  id of node to smooth 
**        topo - (in)  mesh topology 
**
**      Description: Performs smoothing of one internal node. 
**
**      Return Value: new nodal coordinates 
**
**
** -- */

Vec2D MshQuad2D::LaplaceSmooth(IntNode *node,
                                           MshTopo2D *topo)
{
    Vec2D delt_sum(0.0,0.0) ;
    IntNode *adj ;
    TopoAdjVtxIterator iter(topo,node->id) ;
    int num_adj = 0 ;
    for (iter.First() ; iter.More() ; ++iter) {
        adj = pnode_table->Get(iter.AdjVtx()) ;
        delt_sum += adj->coord - node->coord ;
        ++num_adj ;
    }
    if (num_adj == 0) num_adj = 1 ;
    return(delt_sum / num_adj) ;
}




// %(MshQuad2D::CheckValidCoord-bool-|^const-IntNode-|*-Vec2D-|-MshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    CheckValidCoord - check for valid triangles 
**
**      bool CheckValidCoord(
**              IntNode    *node,
**              Vec2D   delta,
**              MshTopo2D *topo) const
**
**        node  - (in)  node to move 
**        delta - (in)  proposed update 
**        topo  - (in)  mesh topology 
**
**      Description: Given a proposed nodal coordinate update this 
**          method checks to see if all the adjacent triangular 
**          elements will still be valid. 
**
**      Return Value: true if the update is valid 
**
**
** -- */

bool MshQuad2D::CheckValidCoord(IntNode *node,
                                      Vec2D delta,
                                      MshTopo2D *topo) const
{
    // here we chech to make sure that all new triangles 
    // have valid metrics.

    Vec2D ncoord, prev ;
    ncoord = node->coord + delta ;
    TopoAdjVtxIterator iter(topo,node->id) ;
    IntNode *adj ;

    if (!iter.More()) return(true) ;
    //++iter ;
    adj = pnode_table->Get(iter.AdjVtx()) ;
    prev = adj->coord ;

    for (++iter ; iter.More() ; ++iter) {
        adj = pnode_table->Get(iter.AdjVtx()) ;
        if (Area(ncoord,prev,adj->coord) <= 0.0) {
            return(false) ;
        }
        prev = adj->coord ;
    }
    return(true) ;
}




// %(MshQuad2D::CheckValidTriangleList-bool-|^const-Vec2D-const|-Vec2D-const|-List-const|<SeamCoordCache>*)
/* ++ ----------------------------------------------------------
**
**    CheckValidTriangleList - check a list of triangles 
**
**      bool CheckValidTriangleList(
**              const Vec2D                coord,
**              const Vec2D                delta,
**              const List<SeamCoordCache>* triangles) const
**
**        coord     - (in)  nodal coord 
**        delta     - (in)  proposed update 
**        triangles - (in)  list of triangles 
**
**      Description: This method check a list of triangular element to 
**          see if they are all valid. 
**
**      Return Value: true if all are valid 
**
**
** -- */

#define MINIMUM_BOUNDARY_METRIC 0.10

bool MshQuad2D::CheckValidTriangleList(const Vec2D coord,
                  const Vec2D delta,
                  const List<SeamCoordCache> *triangles) const
{
    Vec2D ncoord ;
    ncoord = coord + delta ;
    for (int i=0 ; i<triangles->Len() ; ++i) {
        if (triangles->At(i).ignore) continue ;
        double metric = TriMetric(ncoord,triangles->At(i).coord[1],
                                  triangles->At(i).coord[2]) ;
        if (triangles->At(i).boundary) {
            if (metric <= MINIMUM_BOUNDARY_METRIC) return(false) ;
        } else {
            if (metric <= 0.0) return(false) ;
        }
    }
    return(true) ;
}




// %(MshQuad2D::CheckValidQuadList-bool-|^const-Vec2D-const|-Vec2D-const|-List-const|<SeamCoordCache>*)
/* ++ ----------------------------------------------------------
**
**    CheckValidQuadList - check a list of quadrilaterals 
**
**      bool CheckValidQuadList(
**              const Vec2D                coord,
**              const Vec2D                delta,
**              const List<SeamCoordCache>* triangles) const
**
**        coord     - (in)  nodal coord 
**        delta     - (in)  proposed update 
**        triangles - (in)  list of quads 
**
**      Description: This method check a list of quadrilateral element 
**          to see if they are all valid. 
**
**      Return Value: true if all are valid 
**
**
** -- */

bool MshQuad2D::CheckValidQuadList(const Vec2D coord,
                  const Vec2D delta,
                  const List<SeamCoordCache> *quads) const
{
    Vec2D ncoord ;
    ncoord = coord + delta ;
    for (int j=0 ; j<quads->Len() ; ++j) {
        if (quads->At(j).ignore) continue ;
        double metric ;
        if (quads->At(j).num == 4) {
            metric = QuadMetric(ncoord,
                                quads->At(j).coord[1],
                                quads->At(j).coord[2],
                                quads->At(j).coord[3]) ;
        } else {
            metric = TriMetric(ncoord,
                               quads->At(j).coord[1],
                               quads->At(j).coord[2]) ;
        }
        if (quads->At(j).boundary) {
            if (metric <= MINIMUM_BOUNDARY_METRIC) return(false) ;
        } else {
            if (metric <= 0.0) return(false) ;
        }
    }
    return(true) ;
}




// %(MshQuad2D::CheckSmoothCoords-void-|-IntNode-|*-Vec2D-|-MshTopo2D-|*-MshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    CheckSmoothCoords - check a proposed smoothed node location 
**
**      void CheckSmoothCoords(
**              IntNode    *node,
**              Vec2D   delta,
**              MshTopo2D *msh_topo,
**              MshTopo2D *quad_topo)
**
**        node      - (in)  node id 
**        delta     - (in)  proposed coordinate update 
**        msh_topo  - (in)  triangle mesh 
**        quad_topo - (in)  quadrilateral mesh 
**
**      Description: Check a proposed new node location to make sure 
**          that all of the adjacent elements will have a valid shape 
**          (not inverted). If the elements are valid the nodal 
**          coordinates are updated. 
**
**
** -- */

#define VALID_METRIC_TOL 0.05

void MshQuad2D::CheckSmoothCoords(IntNode *node,
                                        Vec2D delta,
                                        MshTopo2D *msh_topo,
                                        MshTopo2D *quad_topo)
{
    int i,j ;

    // here we chech to make sure that all new triangles and
    // quads have valid metrics.  To do this
    // we start with the full delta.  If we detect an
    // invalid element we then try 0.75, 0.5, and 0.25
    // times the delta.

    // first save the coordinates of the triangle

    CoordCache->Clear() ;
    CoordCacheFlags->Clear() ;
    TopoAdjVtxIterator iter(msh_topo,node->id) ;
    for (iter.First() ; iter.More() ; ++iter) {
        IntNode *adj = pnode_table->Get(iter.AdjVtx()) ;
        CoordCache->Append(adj->coord) ;
        CoordCacheFlags->Append(
             iter.CcwElem()==NO_ELEM ? false : true) ; 
    }

    // compute the current minimum shape metric and the check
    // tolerance.  We allow existing bad elements to get up to
    // 10% worse.

    double min_metric = 1000.0 ;
    for (i=0 ; i<CoordCache->Len()-1 ; ++i) {
        if (!CoordCacheFlags->At(i)) continue ;
        double existing = TriMetric(node->coord,CoordCache->At(i),
                                     CoordCache->At(i+1)) ;
        if (existing < min_metric) {
            min_metric = existing ;
        }
    }
    double check_tol = min_metric < VALID_METRIC_TOL ?
                       0.9*min_metric : VALID_METRIC_TOL ;

    // now check the areas

// 8-jun-02 switch check from a positive area to a tolerance
// on the shape metric

    bool valid = true ;
    double tri_factor ;
    for (tri_factor = 1.0 ; tri_factor > 0.0 ; tri_factor -= 0.25) {
        Vec2D ncoord ;
        valid = true ;
        ncoord = node->coord + tri_factor * delta ;
        for (i=0 ; i<CoordCache->Len()-1 ; ++i) {
            if (!CoordCacheFlags->At(i)) continue ;
//            double value = TriMetric(ncoord,CoordCache->At(i),
//                          CoordCache->At(i+1)) ;
            if (TriMetric(ncoord,CoordCache->At(i),
                          CoordCache->At(i+1)) <= check_tol) {
                valid = false ;
                break ;
            }
        }
        if (valid) break ;
    }
    if (!valid) return ;

    // now look at the quads

    CoordCache->Clear() ;
    TopoAdjVtxIterator iter1(quad_topo,node->id) ;
    for (iter1.First() ; iter1.More() ; ++iter1) {
        if (iter1.CcwElem() != NO_ELEM) {
            int num_nodes, nodes[4] ;
            quad_topo->GetElemNodes(iter1.CcwElem(),node->id,
                                    &num_nodes,nodes) ;

            if (num_nodes > 4) throw QuadConversionError() ;
            if (num_nodes < 4) continue ;

            for (int ii=1 ; ii<num_nodes ; ++ii) {
                IntNode *adj = pnode_table->Get(nodes[ii]) ;
                CoordCache->Append(adj->coord) ;
            }
        }
    }

    // look for the minimum existing metric

    min_metric = 1000.0 ;
    for (j=0 ; j<CoordCache->Len() ; j += 3) {
        double existing = QuadMetric(node->coord,
                                   CoordCache->At(j),
                                   CoordCache->At(j+1),
                                   CoordCache->At(j+2)) ;
        if (existing < min_metric) {
            min_metric = existing ;
        }
    }
    if (min_metric < 0.0) {
        check_tol = min_metric ;
    } else {
        check_tol = min_metric < VALID_METRIC_TOL ?
                    0.9*min_metric : VALID_METRIC_TOL ;
    }

    // now check the metrics

    double quad_factor = 1.0 ;
    valid = true ;
    for (quad_factor = 1.0 ; quad_factor > 0.0 ; quad_factor -= 0.25) {
        valid = true ;
        Vec2D ncoord = node->coord + quad_factor * delta ;
        for (j=0 ; j<CoordCache->Len() ; j += 3) {
            double metric = QuadMetric(ncoord,
                                       CoordCache->At(j),
                                       CoordCache->At(j+1),
                                       CoordCache->At(j+2)) ;
            if (metric <= check_tol) {
                valid = false ;
                break ;
            }
        }
        if (valid) break ;
    }

    if (!valid) return ;

    // update the coordinates

    if (tri_factor < quad_factor) {
        node->coord = node->coord + tri_factor * delta ;
    } else {
        node->coord = node->coord + quad_factor * delta ;
    }
}




// %(MshQuad2D::SmoothFront-void-|-int-const|-int-const|-int-const|*-MshTopo2D-|*-MshTopo2D-|*-MshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    SmoothFront - smooth advancing front nodes 
**
**      void SmoothFront(
**              const int       qcase,
**              const int       num_frt_nodes,
**              const int       *frt_nodes,
**              MshTopo2D   *msh_topo,
**              MshTopo2D   *quad_topo,
**              MshEdgeList *edge_list)
**
**        qcase         - (in)  new quad case 
**        num_frt_nodes - (in)  number of front nodes 
**        frt_nodes     - (in)  list of front nodes 
**        msh_topo      - (in)  triangle mesh 
**        quad_topo     - (in)  quadrilateral mesh 
**        edge_list     - (in)  active boundary 
**
**      Description: Perform smoothing of advancing front nodes 
**
**
** -- */

void MshQuad2D::SmoothFront(const int qcase,
                                  const int num_frt_nodes,
                                  const int *frt_nodes,
                                  MshTopo2D *msh_topo,
                                  MshTopo2D *quad_topo,
                                  MshEdgeList *edge_list)
{
    int i ;

    // we smooth the new nodes first, then the adjacent nodes

    switch (qcase) {
        case 1:
            SmoothOneFront(frt_nodes[2],msh_topo,quad_topo,edge_list) ;
            SmoothOneFront(frt_nodes[3],msh_topo,quad_topo,edge_list) ;
            for (i=0 ; i<num_frt_nodes ; ++i) {
                if ((i < 2) || (i > 3))
                    SmoothOneFront(frt_nodes[i],msh_topo,
                                   quad_topo,edge_list) ;
            }
            break ;

        case 2:
            SmoothOneFront(frt_nodes[2],msh_topo,quad_topo,edge_list) ;
            for (i=0 ; i<num_frt_nodes ; ++i) {
                if (i != 2) SmoothOneFront(frt_nodes[i],msh_topo,
                                           quad_topo,edge_list) ;
            }
            break ;

        case 3:
        case 5:
            for (i=0 ; i<num_frt_nodes ; ++i) {
                SmoothOneFront(frt_nodes[i],msh_topo,
                               quad_topo,edge_list) ;
            }
            break ;
    }
}




// %(MshQuad2D::SmoothSeamFront-void-|-MshEdgeList::NewSeamData-|*-MshTopo2D-|*-MshTopo2D-|*-MshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    SmoothSeamFront - smooth advancing front nodes 
**
**      void SmoothSeamFront(
**              MshEdgeList::NewSeamData *sdata,
**              MshTopo2D                   *msh_topo,
**              MshTopo2D                   *quad_topo,
**              MshEdgeList                 *edge_list)
**
**        sdata     - (in)  seaming data 
**        msh_topo  - (in)  triangle mesh 
**        quad_topo - (in)  quadrilateral mesh 
**        edge_list - (in)  active boundary 
**
**      Description: Perform smoothing of advancing front nodes after a 
**          seaming operation. 
**
**
** -- */

void MshQuad2D::SmoothSeamFront(
                                  MshEdgeList::NewSeamData *sdata,
                                  MshTopo2D *msh_topo,
                                  MshTopo2D *quad_topo,
                                  MshEdgeList *edge_list)
{
    SmoothOneFront(sdata->id0,msh_topo,quad_topo,edge_list) ;
    SmoothOneFront(sdata->cw,msh_topo,quad_topo,edge_list) ;
    SmoothOneFront(sdata->ccw,msh_topo,quad_topo,edge_list) ;
}




// %(MshQuad2D::SmoothOneFront-void-|-int-const|-MshTopo2D-|*-MshTopo2D-|*-MshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    SmoothOneFront - smooth one advancing front node 
**
**      void SmoothOneFront(
**              const int       node_id,
**              MshTopo2D   *msh_topo,
**              MshTopo2D   *quad_topo,
**              MshEdgeList *edge_list)
**
**        node_id   - (in)  id of node to smooth 
**        msh_topo  - (in)  triangle mesh 
**        quad_topo - (in)  quadrilateral mesh 
**        edge_list - (in)  active boundary 
**
**      Description: Perform smoothing for one advancing front node. 
**
**
** -- */

void MshQuad2D::SmoothOneFront(const int node_id,
                                     MshTopo2D *msh_topo,
                                     MshTopo2D *quad_topo,
                                     MshEdgeList *edge_list)
{
    IntNode *node = pnode_table->Get(node_id) ;
    if (node->motion == MSH_FLOATING) {

        int num_total = quad_topo->NumAdjElems(node_id) ;
        int num_adj = 
            quad_topo->NumConsecutiveAdjElems(node_id) ;

        Vec2D delta_a, delta_b, delta_c ;

        if ((num_total == 2) && (num_adj == 2)) {

            // cache information about the two elements.  Node
            // coordinates will be stored in the cache in
            // the following
            // order.                 node
            //                 1 *------*------* 4
            //                   |      |      |
            //                   |      |      |
            //                 2 *------*------* 3
            //                          0

            CoordCache->Clear() ;
            IntNode *adj ;
            int num_nodes, nodes[4] ;
            int prev, next ;

            TopoAdjVtxIterator iter(quad_topo,node->id) ;
            quad_topo->GetElemNodes(iter.CcwElem(),node->id,
                                    &num_nodes,nodes) ;

            if (num_nodes < 4) return ;

            adj = pnode_table->Get(nodes[3]) ;   
            CoordCache->Append(adj->coord) ;
            adj = pnode_table->Get(nodes[1]) ;   
            CoordCache->Append(adj->coord) ;
            adj = pnode_table->Get(nodes[2]) ;   
            CoordCache->Append(adj->coord) ;
            prev = nodes[1] ;

            ++iter ;
            quad_topo->GetElemNodes(iter.CcwElem(),node->id,
                                    &num_nodes,nodes) ;

            if (num_nodes < 4) return ;

            adj = pnode_table->Get(nodes[2]) ;   
            CoordCache->Append(adj->coord) ;
            adj = pnode_table->Get(nodes[3]) ;   
            CoordCache->Append(adj->coord) ;
            next = nodes[3] ;

            // now determine new desired edge length

            double ld, la, lq ;
            double ratio = edge_list->GetEdgeLengthRatio() ;
            if (ratio < SMOOTH_BLACKER_TOLERANCE) {
                double d1 = (CoordCache->At(0) -
                             CoordCache->At(2)).Magnitude() ;
                double d2 = (CoordCache->At(0) -
                             CoordCache->At(3)).Magnitude() ;
                double angle = Angle(CoordCache->At(0),
                                     CoordCache->At(3),
                                     CoordCache->At(2)) ;
                if (angle < 0.0) angle += TWO_PI ;
                ld = 0.5 * (d1 + d2) / sin(0.5*angle) ;           
            } else if (ratio < SMOOTH_OWEN_TOLERANCE) {
                double sum = 0 ;
                sum += (CoordCache->At(1) -
                        CoordCache->At(2)).Magnitude() ;
                sum += (CoordCache->At(2) -
                        CoordCache->At(0)).Magnitude() ;
                sum += (CoordCache->At(4) -
                        CoordCache->At(3)).Magnitude() ;
                sum += (CoordCache->At(3) -
                        CoordCache->At(0)).Magnitude() ;

                TopoAdjVtxIterator iter(msh_topo,node_id) ;
                int num = 0 ;

                for (iter.First() ; iter.More() ; ++iter) {
                    if ((iter.AdjVtx() != prev) &&
                        (iter.AdjVtx() != next)) {
                        adj = pnode_table->Get(iter.AdjVtx()) ;
                        sum += (adj->coord - node->coord).Magnitude() ;
                        ++num ;
                    }
                }
                ld = sum / (4 + num) ;
            } else {
                // Laplace Smooth

                double sum = 0 ;
                int num = 0 ;
                TopoAdjVtxIterator iter1(quad_topo,node_id) ;

                for (iter1.First() ; iter1.More() ; ++iter1) {
                    adj = pnode_table->Get(iter1.AdjVtx()) ;
                    sum += (adj->coord - node->coord).Magnitude() ;
                    ++num ;
                }

                TopoAdjVtxIterator iter(msh_topo,node_id) ;
                for (iter.First() ; iter.More() ; ++iter) {
                    if ((iter.AdjVtx() != prev) &&
                        (iter.AdjVtx() != next)) {
                        adj = pnode_table->Get(iter.AdjVtx()) ;
                        sum += (adj->coord - node->coord).Magnitude() ;
                        ++num ;
                    }
                }

                ld = sum / num ;
            }

            // find the length of the vector to the
            // isoparametric point

            delta_a = IsoParametricSmooth(node,quad_topo) ;
            Vec2D v = node->coord - CoordCache->At(0) ;

            la = (v + delta_a).Magnitude() ;

            // find the proposed update

            double lratio = ld / la ;
            delta_b = -v + (delta_a + v) * lratio ;

            // now for the angles, first find bisector
            // of the angle between nodes 3 - 0 - 2
            // (with reference to the figure above)

            Vec2D pb1, pb2, q, tmp ;
            pb1 = BisectNorm(CoordCache->At(0),
                             CoordCache->At(2),
                             CoordCache->At(3)) ;
            pb1 += CoordCache->At(0) ;

            // find the bisector between this and the
            // edge 0 - node and then the point that
            // is the intersection of the new vector
            // and the vector from points 1 - 4

            if (pb1 == node->coord) {
                pb2 = pb1 - CoordCache->At(0) ;
            } else {
                pb2 = BisectNorm(CoordCache->At(0),pb1,node->coord) ;
                if (Angle(CoordCache->At(0),node->coord,pb1) < 0.0) {
                    pb2 = -pb2 ;
                }
            }
            tmp = pb2 + CoordCache->At(0) ;
            q = IntersectLines(node->coord,tmp,CoordCache->At(4),
                               CoordCache->At(1)) ;

            // now the distance to the intersection point

            lq = (q - CoordCache->At(0)).Magnitude() ;

            // finaly find the angle update

            if (ld > lq) ld = (ld + lq) / 2.0 ;

            pb2 = CoordCache->At(0) + pb2*ld ;
            delta_c = pb2 - node->coord ;
            delta_a = (delta_b + delta_c) / 2.0 ;
        } else {
            delta_a = LaplaceSmoothFrnt(node_id,msh_topo,quad_topo) ;
        }
        CheckSmoothCoords(node,delta_a,msh_topo,quad_topo) ;
    }
}




// %(MshQuad2D::SmoothAdjacent-void-|-int-const|-int-const|*-MshTopo2D-|*-bool-const|-bool-const|-MshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    SmoothAdjacent - smooth adjacent nodes 
**
**      void SmoothAdjacent(
**              const int       num_frt_nodes,
**              const int       *frt_nodes,
**              MshTopo2D   *topo,
**              const bool      look_for_bdry,
**              const bool      triangle,
**              MshEdgeList *edge_list)
**
**        num_frt_nodes - (in)  number of front nodes 
**        frt_nodes     - (in)  list of front nodes 
**        topo          - (in)  mesh topology 
**        look_for_bdry - (in)  ignored 
**        triangle      - (in)  true means a triangle mesh 
**        edge_list     - (in)  active boundary 
**
**      Description: Smooth all internal nodes adjacent to active front 
**          nodes. 
**
**
** -- */

void MshQuad2D::SmoothAdjacent(
                                  const int num_frt_nodes,
                                  const int *frt_nodes,
                                  MshTopo2D *topo,
                                  const bool look_for_bdry,
                                  const bool triangle,
                                  MshEdgeList *edge_list)
{
    for (int i=0 ; i<num_frt_nodes-2 ; ++i) {
        SmoothOneAdj(frt_nodes[i+1],frt_nodes[i],
                     frt_nodes[i+2],topo,look_for_bdry,
                     triangle,edge_list) ;
    }
}




// %(MshQuad2D::SmoothOneAdj-void-|-int-const|-int-const|-int-const|-MshTopo2D-|*-bool-const|-bool-const|-MshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    SmoothOneAdj - smooth one adjacent node 
**
**      void SmoothOneAdj(
**              const int       vtx,
**              const int       prev,
**              const int       next,
**              MshTopo2D   *topo,
**              const bool      look_for_bdry,
**              const bool      triangle,
**              MshEdgeList *edge_list)
**
**        vtx           - (in)  input vertex id 
**        prev          - (in)  adjacent vertex on boundary 
**        next          - (in)  adjacent vertex on boundary 
**        topo          - (in)  mesh topology 
**        look_for_bdry - (in)  ignored 
**        triangle      - (in)  true means triangle mesh 
**        edge_list     - (in)  active boundary 
**
**      Description: Smooth all internal nodes adjacent to one active 
**          front node. 
**
**      Return Value: ? 
**
**
** -- */

void MshQuad2D::SmoothOneAdj(const int vtx,
                                   const int prev,
                                   const int next,
                                   MshTopo2D *topo,
                                   const bool look_for_bdry,
                                   const bool triangles,
                                   MshEdgeList *edge_list)
{
    TopoAdjVtxIterator iter(topo,vtx) ;

    for (iter.First() ; iter.More() ; ++iter) {
        if ((iter.AdjVtx() != prev) &&
            (iter.AdjVtx() != next) &&
            !edge_list->ContainsNode(iter.AdjVtx())) {
            IntNode *node = pnode_table->Get(iter.AdjVtx()) ;
            if (node->motion == MSH_FLOATING) {

                Vec2D sum(0.0,0.0) ;
                int num = 0 ;
                TopoAdjVtxIterator aiter(topo,iter.AdjVtx()) ;

                for (aiter.First() ; aiter.More() ; ++aiter) {
                    IntNode *adj =
                        pnode_table->Get(aiter.AdjVtx()) ;
                    sum += adj->coord - node->coord ;
                    ++num ;
                }

                if (num == 0) num = 1 ;
                sum /= num ;

                if (triangles) {
                    CoordCache->Clear() ;
                    TopoAdjVtxIterator aiter(topo,iter.AdjVtx()) ;
                    for (aiter.First() ; aiter.More() ; ++aiter) {
                        IntNode *adj =
                            pnode_table->Get(aiter.AdjVtx()) ;
                        CoordCache->Append(adj->coord) ;
                    }
                    aiter.First() ;
                    IntNode *adj =
                        pnode_table->Get(aiter.AdjVtx()) ;
                    CoordCache->Append(adj->coord) ;

                    // now check the areas

                    bool valid = true ;
                    double tri_factor ;
                    for (tri_factor = 1.0 ; tri_factor > 0.0 ; tri_factor -= 0.25) {
                        valid = true ;
                        Vec2D ncoord = node->coord + tri_factor * sum ;
                        for (int i=0 ; i<CoordCache->Len()-1 ; ++i) {
                            if (Area(ncoord,CoordCache->At(i),
                                     CoordCache->At(i+1)) <= 0.0) {
                                valid = false ;
                                break ;
                            }
                        }
                        if (valid) break ;
                    }
                    if (!valid) return ;
                    node->coord = node->coord + tri_factor * sum ;
                } else {
                    CoordCache->Clear() ;
                    TopoAdjVtxIterator aiter(topo,iter.AdjVtx()) ;
                    for (aiter.First() ; aiter.More() ; ++aiter) {
                        if (aiter.CcwElem() != NO_ELEM) {
                            int num_nodes, nodes[4] ;
                            topo->GetElemNodes(aiter.CcwElem(),
                                    iter.AdjVtx(),
                                    &num_nodes,nodes) ;

                            if (num_nodes > 4) throw QuadConversionError() ;
                            if (num_nodes < 4) continue ;

                            for (int ii=1 ; ii<num_nodes ; ++ii) {
                                IntNode *adj = pnode_table->Get(nodes[ii]) ;
                                CoordCache->Append(adj->coord) ;
                            }
                        }
                    }

                    // now check the metrics

                    double quad_factor = 1.0 ;
                    bool valid = true ;
                    for (quad_factor = 1.0 ; quad_factor > 0.0 ; quad_factor -= 0.25) {
                        valid = true ;
                        Vec2D ncoord = node->coord + quad_factor * sum ;
                        for (int j=0 ; j<CoordCache->Len() ; j += 3) {
                            double metric = QuadMetric(ncoord,
                                       CoordCache->At(j),
                                       CoordCache->At(j+1),
                                       CoordCache->At(j+2)) ;
                            if (metric <= 0.0) {
                                valid = false ;
                                break ;
                            }
                        }
                        if (valid) break ;
                     }

                    if (!valid) return ;
                    node->coord = node->coord + quad_factor * sum ;
                }
            }
        }
    }
}




// %(MshQuad2D::FindQuadData-void-|-MshEdgeList::QdEdge-const|*-MshEdgeList::NewQuadData-|*-MshEdgeList-|*-MshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    FindQuadData - find a candidate quadrilateral 
**
**      void FindQuadData(
**              const MshEdgeList::QdEdge *edge,
**              MshEdgeList::NewQuadData  *qdata,
**              MshEdgeList                  *edge_list,
**              MshTopo2D                    *msh_topo)
**
**        edge      - (in)  base edge for the quad 
**        qdata     - (out) quadrilateral data 
**        edge_list - (in)  active boundary 
**        msh_topo  - (in)  triangle mesh 
**
**      Description: This method identifies and generates information 
**          about a candidate new quadrilateral element. 
**
**
** -- */

void MshQuad2D::FindQuadData(const MshEdgeList::QdEdge *edge,
                                   MshEdgeList::NewQuadData *qdata,
                                   MshEdgeList *edge_list,
                                   MshTopo2D *msh_topo)
{
    qdata->obl = edge->id_0 ;
    qdata->obr = edge->id_1 ;
    qdata->edge_level = edge->level ;

    // determine the nodes that will be at the ends of the
    // side edges.

    int end_code = edge->end_code & NO_SEAM_MASK ;

    if ((end_code == 0) || (end_code == 2)) {

        //  We have either code 0 or code 2, so define the
        //  edge at the right node.
        //
        //             code 0            *    code 2
        //                               |
        //      *-----*=====*-----*      *=====*-----*

        int next_nd_id ;
        IntNode *prev_node, *this_node, *next_node ;

        prev_node = pnode_table->Get(edge->id_0) ;
        this_node = pnode_table->Get(edge->id_1) ;
        next_nd_id = edge_list->GetCCWNode(edge->id_0,edge->id_1) ;
        next_node = pnode_table->Get(next_nd_id) ;
        double base_len = (prev_node->coord-this_node->coord).Magnitude() ;
        qdata->otr = FindElementSide(msh_topo,
                                    this_node,prev_node,
                                    next_node,base_len,true) ;
        if (qdata->otr < 0) {
            qdata->qcase = -1 ;
            return ;
        }

        // check to see if we want to split the triangle
        // edge

        if (!edge_list->ContainsEdge(qdata->obr,qdata->otr)) {
            double base_length = 
                   edge_list->GetEdgeLength(edge->id_0,edge->id_1) ;
            IntNode *found = pnode_table->Get(qdata->otr) ;
            double side_length =
                (this_node->coord - found->coord).Magnitude() ;

            bool split = false ;
            if (side_length/base_length > TRIANGLE_SPLIT_TOLERANCE) {
                split = true ;

                // check for a special case where the otl and otr are
                // the same and the triangle has very little area.  In
                // this case we do nothing

                if (qdata->otl == qdata->otr) {
                    double area = Area(this_node->coord,
                                       next_node->coord,
                                       found->coord) ;
                    double max = side_length > base_length ?
                                 side_length : base_length ;
                    double eps = 2.0*area/max ;
                    if (eps < 0.05*max) split = false ;
                }

                // if the found vertex is on the boundary then
                // don't split

                if (edge_list->ContainsNode(found->id)) split = false ;   
            }

            if (split) {
                double nx = 0.5 * (this_node->coord[0] + found->coord[0]) ;
                double ny = 0.5 * (this_node->coord[1] + found->coord[1]) ;
                int new_id = NewNode(nx,ny,INTERIOR,
                                              MSH_FLOATING,true) ;
                int tri_elem, num_nodes, tri_nodes[3] ;

                tri_elem = msh_topo->GetCCWElem(qdata->obr,qdata->otr) ;
                msh_topo->GetElemNodes(tri_elem,qdata->obr,
                                       &num_nodes,tri_nodes) ;

                if (num_nodes > 3) throw QuadConversionError() ;

                msh_topo->DeleteElement(num_nodes,tri_nodes) ;

                int ccw_id = tri_nodes[2] ;
                int elem_0 = tri_elem ;
                int elem_1 = NewElemNum() ;

                tri_elem = msh_topo->GetCWElem(qdata->obr,qdata->otr) ;
                msh_topo->GetElemNodes(tri_elem,qdata->otr,
                                       &num_nodes,tri_nodes) ;

                if (num_nodes > 3) throw QuadConversionError() ;

                msh_topo->DeleteElement(num_nodes,tri_nodes) ;

                int cw_id = tri_nodes[2] ;
                int elem_2 = tri_elem ;
                int elem_3 = NewElemNum() ;

                // add the new triangles

                int new_tri0[3] = {qdata->obr,new_id,ccw_id} ;
                int new_tri1[3] = {new_id,qdata->otr,ccw_id} ;
    
                msh_topo->InsertElement(elem_0,3,new_tri0) ;
                msh_topo->InsertElement(elem_1,3,new_tri1) ;

                // add the new triangles

                int new_tri2[3] = {qdata->otr,new_id,cw_id} ;
                int new_tri3[3] = {new_id,qdata->obr,cw_id} ;
 
                msh_topo->InsertElement(elem_2,3,new_tri2) ;
                msh_topo->InsertElement(elem_3,3,new_tri3) ;

                qdata->otr = new_id ;
            }
        }

        // now we want to check to see if by adding this
        // edge we close off a region with an odd number
        // edges.  If so, split the edge and return the
        // new vertex.

        // first look to see if any of the edges about the
        // node in the mesh topology is also part of the
        // boundary topology

        if (edge_list->ContainsNode(qdata->otr)) {

            // now count the number of edges until we get back
            // to the start location

            int num = 0 ;
            int vtx0 = qdata->obr ;
            int vtx1 =
                edge_list->GetCCWNode(qdata->obl,qdata->obr) ;

            while (vtx1 != qdata->otr) {
                int tmp = vtx1 ;
                vtx1 = edge_list->GetCCWNode(vtx0,vtx1) ;
                if (vtx1 == qdata->obr) {
                    num = 2 ;
                    break ;
                }
                vtx0 = tmp ;
                num++ ;
            }

            // check for odd number and split the edge.

            if ((num > 5) && ((num % 2) == 1)) {
                IntNode *node_0 = pnode_table->Get(qdata->otr) ;
                IntNode *node_1 = pnode_table->Get(qdata->obr) ;
                double nx = 0.5 * (node_0->coord[0] + node_1->coord[0]) ;
                double ny = 0.5 * (node_0->coord[1] + node_1->coord[1]) ;
                int new_id = NewNode(nx,ny,INTERIOR,
                                              MSH_FLOATING,true) ;

                int tri_elem, num_nodes, tri_nodes[3] ;

                tri_elem = msh_topo->GetCCWElem(qdata->obr,qdata->otr) ;
                msh_topo->GetElemNodes(tri_elem,qdata->obr,
                                       &num_nodes,tri_nodes) ;
                if (num_nodes > 3) throw QuadConversionError() ;

                msh_topo->DeleteElement(num_nodes,tri_nodes) ;

                int ccw_id = tri_nodes[2] ;
                int elem_0 = tri_elem ;
                int elem_1 = NewElemNum() ;

                tri_elem = msh_topo->GetCWElem(qdata->obr,qdata->otr) ;
                msh_topo->GetElemNodes(tri_elem,qdata->otr,
                                       &num_nodes,tri_nodes) ;

                if (num_nodes > 3) throw QuadConversionError() ;

                msh_topo->DeleteElement(num_nodes,tri_nodes) ;

                int cw_id = tri_nodes[2] ;
                int elem_2 = tri_elem ;
                int elem_3 = NewElemNum() ;

                // add the new triangles

                int new_tri0[3] = {qdata->obr,new_id,ccw_id} ;
                int new_tri1[3] = {new_id,qdata->otr,ccw_id} ;

                msh_topo->InsertElement(elem_0,3,new_tri0) ;
                msh_topo->InsertElement(elem_1,3,new_tri1) ;

                // add the new triangles

                int new_tri2[3] = {qdata->otr,new_id,cw_id} ;
                int new_tri3[3] = {new_id,qdata->obr,cw_id} ;
 
                msh_topo->InsertElement(elem_2,3,new_tri2) ;
                msh_topo->InsertElement(elem_3,3,new_tri3) ;

                qdata->otr = new_id ;
            }
        }
    } else {
        qdata->otr = edge_list->GetCCWNode(edge->id_0,edge->id_1) ;
    }

    if ((end_code == 0) || (end_code == 1)) {

        //  We have either code 0 or code 2, so define the
        //  edge at the right node.
        //
        //             code 0              code 1   *
        //                                          |
        //      *-----*=====*-----*     *-----*=====*

        int prev_nd_id ;
        IntNode *prev_node, *this_node, *next_node ;

        this_node = pnode_table->Get(edge->id_0) ;
        next_node = pnode_table->Get(edge->id_1) ;
        prev_nd_id = edge_list->GetCWNode(edge->id_0,edge->id_1) ;
        prev_node = pnode_table->Get(prev_nd_id) ;
        double base_len = (next_node->coord-this_node->coord).Magnitude() ;
        qdata->otl = FindElementSide(msh_topo,
                                    this_node,prev_node,
                                    next_node,base_len,false) ;
        if (qdata->otl < 0) {
            qdata->qcase = -1 ;
            return ;
        }

        // check to see if we want to split the triangle
        // edge

        if (!edge_list->ContainsEdge(qdata->otl,qdata->obl)) {
            double base_length = 
                   edge_list->GetEdgeLength(edge->id_0,edge->id_1) ;
            IntNode *found = pnode_table->Get(qdata->otl) ;
            double side_length =
                (this_node->coord - found->coord).Magnitude() ;

            bool split = false ;
            if (side_length/base_length > TRIANGLE_SPLIT_TOLERANCE) {
                split = true ;

                // check for a special case where the otl and otr are
                // the same and the triangle has very little area.  In
                // this case we do nothing

                if (qdata->otl == qdata->otr) {
                    double area = Area(this_node->coord,
                                       next_node->coord,
                                       found->coord) ;
                    double max = side_length > base_length ?
                                 side_length : base_length ;
                    double eps = 2.0*area/max ;
                    if (eps < 0.05*max) split = false ;
                }

                // if the found vertex is on the boundary then
                // don't split

                if (edge_list->ContainsNode(found->id)) split = false ;   
            }

            if (split) {
                double nx = 0.5 * (this_node->coord[0] + found->coord[0]) ;
                double ny = 0.5 * (this_node->coord[1] + found->coord[1]) ;
                int new_id = NewNode(nx,ny,INTERIOR,
                                              MSH_FLOATING,true) ;

                int tri_elem, num_nodes, tri_nodes[3] ;

                tri_elem = msh_topo->GetCCWElem(qdata->obl,qdata->otl) ;
                msh_topo->GetElemNodes(tri_elem,qdata->obl,
                                       &num_nodes,tri_nodes) ;

                if (num_nodes > 3) throw QuadConversionError() ;

                msh_topo->DeleteElement(num_nodes,tri_nodes) ;

                int ccw_id = tri_nodes[2] ;
                int elem_0 = tri_elem ;
                int elem_1 = NewElemNum() ;

                tri_elem = msh_topo->GetCWElem(qdata->obl,qdata->otl) ;
                msh_topo->GetElemNodes(tri_elem,qdata->otl,
                                       &num_nodes,tri_nodes) ;

                if (num_nodes > 3) throw QuadConversionError() ;

                msh_topo->DeleteElement(num_nodes,tri_nodes) ;

                int cw_id = tri_nodes[2] ;
                int elem_2 = tri_elem ;
                int elem_3 = NewElemNum() ;

                // add the new triangles

                int new_tri0[3] = {qdata->obl,new_id,ccw_id} ;
                int new_tri1[3] = {new_id,qdata->otl,ccw_id} ;

                msh_topo->InsertElement(elem_0,3,new_tri0) ;
                msh_topo->InsertElement(elem_1,3,new_tri1) ;

                // add the new triangles

                int new_tri2[3] = {qdata->otl,new_id,cw_id} ;
                int new_tri3[3] = {new_id,qdata->obl,cw_id} ;

                msh_topo->InsertElement(elem_2,3,new_tri2) ;
                msh_topo->InsertElement(elem_3,3,new_tri3) ;

                qdata->otl = new_id ;
            }
        }


        // now we want to check to see if by adding this
        // edge we close off a region with an odd number
        // edges.  If so, split the edge and return the
        // new vertex.

        // first look to see if any of the edges about the
        // node in the mesh topology is also part of the
        // boundary topology

        if (edge_list->ContainsNode(qdata->otl)) {

            TopoAdjVtxIterator iter(msh_topo,qdata->otl) ;
            for (iter.First() ; iter.More() ; ++iter) {
                if (edge_list->ContainsEdge(qdata->otl,iter.AdjVtx())) {
                    break ;
                }
            }

            // now count the number of edges until we get back
            // to the start location

            int num = 0 ;
            int vtx0 = qdata->otl ;
            int vtx1 = iter.AdjVtx() ;

            while (vtx1 != qdata->obl) {
                int tmp = vtx1 ;
                vtx1 = edge_list->GetCCWNode(vtx0,vtx1) ;
                if (vtx1 == qdata->otl) {
                    num = 2 ;
                    break ;
                }
                vtx0 = tmp ;
                num++ ;
            }

            // check for odd number and split the edge.

            if ((num > 5) && ((num % 2) == 1)) {
                IntNode *node_0 = pnode_table->Get(qdata->otl) ;
                IntNode *node_1 = pnode_table->Get(qdata->obl) ;
                double nx = 0.5 * (node_0->coord[0] + node_1->coord[0]) ;
                double ny = 0.5 * (node_0->coord[1] + node_1->coord[1]) ;
                int new_id = NewNode(nx,ny,INTERIOR,
                                              MSH_FLOATING,true) ;

                int tri_elem, num_nodes, tri_nodes[3] ;

                tri_elem = msh_topo->GetCCWElem(qdata->obl,qdata->otl) ;
                msh_topo->GetElemNodes(tri_elem,qdata->obl,
                                       &num_nodes,tri_nodes) ;
                if (num_nodes > 3) throw QuadConversionError() ;

                msh_topo->DeleteElement(num_nodes,tri_nodes) ;

                int ccw_id = tri_nodes[2] ;
                int elem_0 = tri_elem ;
                int elem_1 = NewElemNum() ;

                tri_elem = msh_topo->GetCWElem(qdata->obl,qdata->otl) ;
                msh_topo->GetElemNodes(tri_elem,qdata->otl,
                                       &num_nodes,tri_nodes) ;

                if (num_nodes > 3) throw QuadConversionError() ;

                msh_topo->DeleteElement(num_nodes,tri_nodes) ;

                int cw_id = tri_nodes[2] ;
                int elem_2 = tri_elem ;
                int elem_3 = NewElemNum() ;

                // add the new triangles

                int new_tri0[3] = {qdata->obl,new_id,ccw_id} ;
                int new_tri1[3] = {new_id,qdata->otl,ccw_id} ;

                msh_topo->InsertElement(elem_0,3,new_tri0) ;
                msh_topo->InsertElement(elem_1,3,new_tri1) ;

                // add the new triangles

                int new_tri2[3] = {qdata->otl,new_id,cw_id} ;
                int new_tri3[3] = {new_id,qdata->obl,cw_id} ;

                msh_topo->InsertElement(elem_2,3,new_tri2) ;
                msh_topo->InsertElement(elem_3,3,new_tri3) ;

                qdata->otl = new_id ;
            }
        }
    } else {
        qdata->otl = edge_list->GetCWNode(edge->id_0,edge->id_1) ;
    }

    edge_list->ClassifyQuad(qdata->obl,qdata->obr,qdata->otr,qdata->otl,
                            qdata) ;

    if (qdata->qcase == 2) {
        qdata->ratio = edge_list->GetEdgeLength(qdata->bl,qdata->br) /
                       edge_list->GetEdgeLength(qdata->br,qdata->tr) ;
        IntNode *node_0 ;
        IntNode *node_1 ;
        if (qdata->ratio > 1.0) {
            node_0 = pnode_table->Get(qdata->bl) ;
            node_1 = pnode_table->Get(qdata->br) ;
        } else {
            node_0 = pnode_table->Get(qdata->br) ;
            node_1 = pnode_table->Get(qdata->tr) ;
        }
        bool bdry_edge = (node_0->type == BOUNDARY) &&
                         (node_1->type == BOUNDARY) ;
        if (bdry_edge) {
            qdata->valid_transition_split = false ;
        } else {
            qdata->valid_transition_split = true ;
        }
    } else if (qdata->qcase == 3) {
        qdata->ratio = edge_list->GetEdgeLength(qdata->tl,qdata->bl) /
                       edge_list->GetEdgeLength(qdata->br,qdata->tr) ;
        IntNode *node_0 ;
        IntNode *node_1 ;
        if (qdata->ratio > 1.0) {
            node_0 = pnode_table->Get(qdata->bl) ;
            node_1 = pnode_table->Get(qdata->br) ;
        } else {
            node_0 = pnode_table->Get(qdata->br) ;
            node_1 = pnode_table->Get(qdata->tr) ;
        }
        bool bdry_edge = (node_0->type == BOUNDARY) &&
                         (node_1->type == BOUNDARY) ;
        if (bdry_edge) {
            qdata->valid_transition_split = false ;
        } else {
            qdata->valid_transition_split = true ;
        }

        node_0 = pnode_table->Get(qdata->tl) ;
        node_1 = pnode_table->Get(qdata->tr) ;
        double length = (node_1->coord - node_0->coord).Magnitude() ;
        double ratio_cw = length /
                   edge_list->GetEdgeLength(qdata->front_nodes[0],
                                            qdata->front_nodes[1]) ;
        double ratio_ccw = length /
                   edge_list->GetEdgeLength(qdata->front_nodes[2],
                                            qdata->front_nodes[3]) ;
        qdata->top_template_ratio = ratio_cw > ratio_ccw ?
                                    ratio_cw : ratio_ccw ;

        double ratio_tb = length /
                   edge_list->GetEdgeLength(qdata->bl,qdata->br) ;
        if (ratio_tb > qdata->top_template_ratio)
            qdata->top_template_ratio = ratio_tb ;

        length = edge_list->GetEdgeLength(qdata->bl,qdata->br) ;
        ratio_cw = length /
                   edge_list->GetEdgeLength(qdata->br,qdata->tr) ;
        ratio_ccw = length /
                    edge_list->GetEdgeLength(qdata->tl,qdata->bl) ;
        qdata->base_template_ratio = ratio_cw > ratio_ccw ?
                                     ratio_cw : ratio_ccw ;
    }
}




// %(MshQuad2D::FindSeamData-bool-|-MshEdgeList::QdEdge-const|*-MshEdgeList::NewSeamData-|*-MshEdgeList-|*-MshTopo2D-|*-MshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    FindSeamData - find data for a seam operation 
**
**      bool FindSeamData(
**              const MshEdgeList::QdEdge *edge,
**              MshEdgeList::NewSeamData  *sdata,
**              MshEdgeList                  *edge_list,
**              MshTopo2D                    *msh_topo,
**              MshTopo2D                    *quad_topo)
**
**        edge      - (in)  base edge 
**        sdata     - (out) seam data 
**        edge_list - (in)  active boundary 
**        msh_topo  - (in)  triangle mesh 
**        quad_topo - (in)  quadrilateral mesh 
**
**      Description: This method generates information needed for a 
**          seam operation. 
**
**      Return Value: true if valid seam data is returned 
**
**
** -- */

#define PLACE_TOL 0.005

bool MshQuad2D::FindSeamData(const MshEdgeList::QdEdge *edge,
                                   MshEdgeList::NewSeamData *sdata,
                                   MshEdgeList *edge_list,
                                   MshTopo2D *msh_topo,
                                   MshTopo2D *quad_topo)
{
    double small_ang ;

    if (edge->angle_0 < edge->angle_1) {
        sdata->id0 = edge_list->GetCWNode(edge->id_0,edge->id_1) ;
        sdata->id1 = edge->id_0 ;
        sdata->id2 = edge->id_1 ;
        sdata->cw_side = true ;
        small_ang = edge->angle_0 ;
    } else {
        sdata->id0 = edge->id_0 ;
        sdata->id1 = edge->id_1 ;
        sdata->id2 = edge_list->GetCCWNode(edge->id_0,edge->id_1) ;
        sdata->cw_side = false ;
        small_ang = edge->angle_1 ;
    }
    sdata->cw = edge_list->GetCWNode(sdata->id0,sdata->id1) ;
    sdata->ccw = edge_list->GetCCWNode(sdata->id1,sdata->id2) ;
    sdata->ratio = edge_list->GetEdgeLength(sdata->id0,sdata->id1) /
                   edge_list->GetEdgeLength(sdata->id1,sdata->id2) ;
    sdata->do_template = false ;
    int level_0 = edge_list->GetEdgeLevel(sdata->id0,sdata->id1) ;
    int level_1 = edge_list->GetEdgeLevel(sdata->id1,sdata->id2) ;
    sdata->min_level = (level_0 < level_1) ? level_0 : level_1 ;

    // check to make sure that the end nodes are movable

    IntNode *tid0 = pnode_table->Get(sdata->id0) ;
    if (tid0->motion == MSH_FIXED) {
        IntNode *tid2 = pnode_table->Get(sdata->id2) ;
        if (tid2->motion == MSH_FIXED) {
            sdata->cw = 0 ;
            sdata->ccw = 1 ;
            return(false) ;
        }
    }

    // check for the case where we want to do a template
    // transition before we do the seam.

    double len_cw = edge_list->GetEdgeLength(sdata->id0,sdata->id1) ;
    double len_ccw = edge_list->GetEdgeLength(sdata->id1,sdata->id2) ;
    IntNode *id1 = pnode_table->Get(sdata->id1) ;

    if (len_cw > len_ccw) {
        IntNode *ccw = pnode_table->Get(sdata->ccw) ;
        IntNode *id0 = pnode_table->Get(sdata->id0) ;
        double angle = Angle(id1->coord,ccw->coord,id0->coord) ;
        if (angle < small_ang) {
            double len = (ccw->coord - id1->coord).Magnitude() ;
            sdata->do_template = (len < len_cw) ? true : false ;
            sdata->template_cw = true ;
        }
        int elem = !quad_topo->HasVtx(sdata->id1) ? NO_ELEM : 
            quad_topo->GetCCWElem(sdata->id1,sdata->id0) ;
        if (elem != NO_ELEM) {
            int nnum, nodes[4] ;
            quad_topo->GetElemNodes(elem,sdata->id1,&nnum,nodes) ;
            if (nnum < 4) sdata->do_template = false ;
        }
    } else {
        IntNode *cw = pnode_table->Get(sdata->cw) ;
        IntNode *id2 = pnode_table->Get(sdata->id2) ;
        double angle = Angle(id1->coord,id2->coord,cw->coord) ;
        if (angle < small_ang) {
            double len = (cw->coord - id1->coord).Magnitude() ;
            sdata->do_template = (len < len_ccw) ? true : false ;
            sdata->template_cw = false ;
        }
//        int elem = !quad_topo->HasVtx(sdata->id1) ? NO_ELEM : 
//            quad_topo->GetCCWElem(sdata->id2,sdata->id1) ;
        int elem = !quad_topo->HasVtx(sdata->id2) ? NO_ELEM : 
            quad_topo->GetCCWElem(sdata->id2,sdata->id1) ;
        if (elem != NO_ELEM) {
            int nnum, nodes[4] ;
            quad_topo->GetElemNodes(elem,sdata->id1,&nnum,nodes) ;
            if (nnum < 4) sdata->do_template = false ;
        }
    }

    // try to recover an edge between id0 and id2 and clear the triangle

    if (!RecoverTopEdge(sdata->id2,sdata->id0,msh_topo,edge_list))
        return(false) ;

    int boundary[3] ;
    boundary[0] = sdata->id1 ;
    boundary[1] = sdata->id2 ;
    boundary[2] = sdata->id0 ;

    MshTopo2D::Edge bedge ;
    bedge.nd0 = sdata->id1 ;
    bedge.nd1 = sdata->id2 ;
    bedge.elem0 = msh_topo->GetCCWElem(bedge.nd0,bedge.nd1) ;
    bedge.elem1 = NO_ELEM ;
    ClearRegion(3,boundary,msh_topo,&bedge) ;
    msh_topo->InsertTriangle(bedge.elem0,sdata->id1,sdata->id2,
                             sdata->id0) ;

    // check for the special case where there are only two quad
    // elements about the mid point of the seam.  In this case
    // we can close the seam and turn the two quads into one

    sdata->merge_quads = false ;
    if (msh_topo->NumAdjElems(sdata->id1) == 1) {
        int num, triangle[3] ;
        TopoAdjVtxIterator titerq(msh_topo,sdata->id1) ;
        msh_topo->GetElemNodes(titerq.CcwElem(),sdata->id1,
                               &num,triangle) ;

        if (num > 3) throw QuadConversionError() ;

        if (quad_topo->NumAdjElems(sdata->id1) == 2) {
            TopoAdjVtxIterator iterq(quad_topo,sdata->id1) ;
            int num, nodes[4] ;
            num = 0 ;
            for (iterq.First() ; iterq.More() ; ++iterq) ++num ;
            if (num == 3) {
                iterq.First() ;
                quad_topo->GetElemNodes(iterq.CcwElem(),
                                        sdata->id1,&num,nodes) ;
                if (num == 4) {
                    int save0 = nodes[1] ;
                    int save1 = nodes[3] ;
                    ++iterq ;
                    quad_topo->GetElemNodes(iterq.CcwElem(),
                                            sdata->id1,&num,nodes) ;
                    if (num == 4) {
                        if ((triangle[2] == save0) &&
                            (save1 == nodes[1]) &&
                            (nodes[3] == triangle[1]))
                            sdata->merge_quads = true ;
                    }
                }
            }   
        }
    }

    // check to see if we can place the collapsed node if not, then
    // this is not a valid operation.

    // cache the coordinates of all the nodes for the triangular
    // elements about the first node

    List<SeamCoordCache> nd0_triangles ;
    List<SeamCoordCache> nd0_quads ;
    TopoAdjVtxIterator iterm(msh_topo,sdata->id0) ;
    int i ;
    
    for (iterm.First() ; iterm.More() ; ++iterm) {
        if (iterm.CcwElem() != NO_ELEM) {
            IntNode *adj ;
            SeamCoordCache cache ;
            int num, nodes[3], num_fix = 0 ;
            msh_topo->GetElemNodes(iterm.CcwElem(),
                                   sdata->id0,&num,nodes) ;
            cache.num = 3 ;
            cache.boundary = false ;
            for (i=0 ; i<num ; ++i) {
                adj = pnode_table->Get(nodes[i]) ;
                if (adj->motion == MSH_FIXED) ++num_fix ;
                cache.coord[i] = adj->coord ;
            }
            if (num_fix >= 2) cache.boundary = true ;
            cache.ignore = ((nodes[1] == sdata->id2) ||
                            (nodes[num-1] == sdata->id2)) ;
            nd0_triangles.Append(cache) ;
        }
    }

    // do the same for the quad elements

    sdata->min_fm = 2.0 ;
    TopoAdjVtxIterator iterq(quad_topo,sdata->id0) ;
    for (iterq.First() ; iterq.More() ; ++iterq) {
        if (iterq.CcwElem() != NO_ELEM) {
            IntNode *adj ;
            SeamCoordCache cache ;
            int num, nodes[4], num_fix = 0 ;
            quad_topo->GetElemNodes(iterq.CcwElem(),
                                    sdata->id0,
                                    &num,nodes) ;     
            cache.num = num ;
            cache.boundary = false ;
            cache.ignore = false ;
            for (i=0 ; i<num ; ++i) {
                adj = pnode_table->Get(nodes[i]) ;
                if (adj->motion == MSH_FIXED) ++num_fix ;
                cache.coord[i] = adj->coord ;
            }
            if (num_fix >= 2) cache.boundary = true ;
            if (sdata->merge_quads) {
                cache.ignore = ((nodes[1] == sdata->id1) ||
                                (nodes[num-1] == sdata->id1)) ;
            }
            nd0_quads.Append(cache) ;
        }
        SeamCoordCache cache = nd0_quads[-1] ;
        if (cache.num == 4) {
            double fm = FlatnessMetric(cache.coord[0],
                                       cache.coord[1],
                                       cache.coord[2],
                                       cache.coord[3]) ;
            if (fm < sdata->min_fm) sdata->min_fm = fm ;
        }
    }

    // look for the range of positions were we can place this
    // node on the vector between this in the second node (normalized
    // between 0 and 1) and still have all the adjacent elements
    // have valid shapes.

    double range_0, range_1 ;
    Vec2D vect ;
    IntNode *id0 = pnode_table->Get(sdata->id0) ;
    IntNode *id2 = pnode_table->Get(sdata->id2) ;

    vect = id2->coord - id0->coord ;

    if (id0->motion == MSH_FIXED) {
        Vec2D trial(0.0,0.0) ;
        if (!CheckValidTriangleList(id0->coord,trial,&nd0_triangles) ||
            !CheckValidQuadList(id0->coord,trial,&nd0_quads)) return(false) ;
        range_0 = 0.0 ;
    } else {
        Vec2D trial = vect ;
        if (CheckValidTriangleList(id0->coord,trial,&nd0_triangles) &&
            CheckValidQuadList(id0->coord,trial,&nd0_quads)) {
            range_0 = 1.0 ;
        } else {
            double low = 0.0 ;
            double high = 1.0 ;
            double mid = 0.5 ;
            while ((high-low) > PLACE_TOL) {
                trial = mid * vect ;
                if (CheckValidTriangleList(id0->coord,trial,&nd0_triangles) &&
                    CheckValidQuadList(id0->coord,trial,&nd0_quads))
                    low = mid ;
                else
                    high = mid ;
                mid = 0.5 * (low + high) ;
            }
            range_0 = mid ;
        }
    }


    // cache the coordinates of all the nodes for the triangular
    // elements about the second

    List<SeamCoordCache> nd2_quads ;
    List<SeamCoordCache> nd2_triangles ;

    for (iterm.NewVtx(sdata->id2) ; iterm.More() ; ++iterm) {
        if (iterm.CcwElem() != NO_ELEM) {
            IntNode *adj ;
            int num, nodes[3], num_fix = 0 ;
            SeamCoordCache cache ;
            msh_topo->GetElemNodes(iterm.CcwElem(),
                                   sdata->id2,&num,nodes) ;
            cache.num = 3 ;
            cache.boundary = false ;
            for (i=0 ; i<num ; ++i) {
                adj = pnode_table->Get(nodes[i]) ;
                if (adj->motion == MSH_FIXED) ++num_fix ;
                cache.coord[i] = adj->coord ;
            }
            if (num_fix >= 2) cache.boundary = true ;
            cache.ignore = ((nodes[1] == sdata->id0) ||
                            (nodes[num-1] == sdata->id0)) ;
            nd2_triangles.Append(cache) ;
        }
    }

    // do the same for the quad elements

    for (iterq.NewVtx(sdata->id2) ; iterq.More() ; ++iterq) {
        if (iterq.CcwElem() != NO_ELEM) {
            IntNode *adj ;
            int num, nodes[4], num_fix = 0 ;
            SeamCoordCache cache ;
            quad_topo->GetElemNodes(iterq.CcwElem(),
                                    sdata->id2,
                                    &num,nodes) ;     
            cache.num = num ;
            cache.boundary = false ;
            cache.ignore = false ;
            for (i=0 ; i<num ; ++i) {
                adj = pnode_table->Get(nodes[i]) ;
                if (adj->motion == MSH_FIXED) ++num_fix ;
                cache.coord[i] = adj->coord ;
            }
            if (num_fix >= 2) cache.boundary = true ;
            if (sdata->merge_quads) {
                cache.ignore = ((nodes[1] == sdata->id1) ||
                                (nodes[num-1] == sdata->id1)) ;
            }
            nd2_quads.Append(cache) ;
        }
        SeamCoordCache cache = nd2_quads[0] ;
        if (cache.num == 4) {
            double fm = FlatnessMetric(cache.coord[0],
                                       cache.coord[1],
                                       cache.coord[2],
                                       cache.coord[3]) ;
            if (fm < sdata->min_fm) sdata->min_fm = fm ;
        }
    }


    // look for the range of positions were we can place this
    // node on the vector between this in the first node (normalized
    // between 0 and 1) and still have all the adjacent elements
    // have valid shapes.

    if (id2->motion == MSH_FIXED) {
        Vec2D trial(0.0,0.0) ;
        if (!CheckValidTriangleList(id2->coord,trial,&nd2_triangles) ||
            !CheckValidQuadList(id2->coord,trial,&nd2_quads)) return(false) ;
        range_1 = 1.0 ;
    } else {
        Vec2D trial = -vect ;
        if (CheckValidTriangleList(id2->coord,trial,&nd2_triangles) &&
            CheckValidQuadList(id2->coord,trial,&nd2_quads)) {
            range_1 = 0.0 ;
        } else {
            double low = 0.0 ;
            double high = 1.0 ;
            double mid = 0.5 ;
            while ((high-low) > PLACE_TOL) {
                trial = mid * -vect ;
                if (CheckValidTriangleList(id2->coord,trial,&nd2_triangles) &&
                    CheckValidQuadList(id2->coord,trial,&nd2_quads))
                    low = mid ;
                else
                    high = mid ;
                mid = 0.5 * (low + high) ;
            }
            range_1 = 1.0 - mid ;
        }
    }

    // find the location to place the collapsed node and still
    // have valid elements

    if (id0->motion == MSH_FIXED) {
        if (range_1 == 0.0) {
            sdata->coord = id0->coord ;
        } else {
            return(false) ;
        }
    } else if (id2->motion == MSH_FIXED) {
        if (range_0 == 1.0) {
            sdata->coord = id2->coord ;
        } else {
            return(false) ;
        }
    } else if ((range_1 - range_0) >= 0.0){
        return(false) ;
    } else {
        double fact = 0.5 * (range_0 + range_1) ;
        sdata->coord = id0->coord + fact*vect ;
    }
    return(true) ;
}




// %(MshQuad2D::DoSeam-bool-|-MshEdgeList::NewSeamData-|*-MshTopo2D-|*-MshTopo2D-|*-MshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    DoSeam - do a seam update 
**
**      bool DoSeam(
**              MshEdgeList::NewSeamData *sdata,
**              MshTopo2D                   *msh_topo,
**              MshTopo2D                   *quad_topo,
**              MshEdgeList                 *edge_list)
**
**        sdata     - (in)  seam data 
**        msh_topo  - (i/o) triangle mesh 
**        quad_topo - (i/o) quadrilateral mesh 
**        edge_list - (i/o) active boundary 
**
**      Description: This method performs a seam update 
**
**      Return Value: true if the update was successful 
**
**
** -- */

bool MshQuad2D::DoSeam(MshEdgeList::NewSeamData *sdata,
                             MshTopo2D *msh_topo,
                             MshTopo2D *quad_topo,
                             MshEdgeList *edge_list)
{
    if (!RecoverTopEdge(sdata->id2,sdata->id0,msh_topo,edge_list))
        return(false) ;

    // now we want to find node that forms a triangle with
    // node 1 and the cw node.  To do this we loop around
    // the cw node in the ccw direction.  Once we find
    // node 1, the node we want is the next one.

    TopoAdjVtxIterator iter0(msh_topo,sdata->id0) ;

    while (iter0.AdjVtx() != sdata->id2) ++iter0 ;
    ++iter0 ;
    if (!iter0.More()) iter0.First() ;

    int opp = iter0.AdjVtx() ;

    // check for special case of a triangle

    if (sdata->id1 == opp) {
        edge_list->UpdateTriangle(sdata->id0,sdata->id1,sdata->id2) ;
        return(true) ;
    }

    // now clear the quad

    int boundary[4] ;
    MshTopo2D::Edge bedge ;
    if (sdata->cw_side) {
        boundary[0] = sdata->id1 ;
        boundary[1] = sdata->id2 ;
        boundary[2] = opp ;
        boundary[3] = sdata->id0 ;

        bedge.nd0 = sdata->id1 ;
        bedge.nd1 = sdata->id2 ;
        bedge.elem0 = msh_topo->GetCCWElem(bedge.nd0,bedge.nd1) ;
        bedge.elem1 = NO_ELEM ;
    } else {
        boundary[0] = sdata->id0 ;
        boundary[1] = sdata->id1 ;
        boundary[2] = sdata->id2 ;
        boundary[3] = opp ;

        bedge.nd0 = sdata->id0 ;
        bedge.nd1 = sdata->id1 ;
        bedge.elem0 = msh_topo->GetCCWElem(bedge.nd0,bedge.nd1) ;
        bedge.elem1 = NO_ELEM ;
    }

    ClearRegion(4,boundary,msh_topo,&bedge) ;

    // update the edge list

    int save_id0 = sdata->id0 ;
    int save_id2 = sdata->id2 ;
    int keep = edge_list->UpdateSeam(sdata) ;
    int remove = keep == save_id0 ? save_id2 : save_id0 ;

    // check for the special case where doing the seam either
    // collapsed or inverted two adjacent triangles

    int elem0 = msh_topo->GetCCWElem(opp,save_id2) ;
    int elem1 = msh_topo->GetCCWElem(save_id0,opp) ;

    if ((elem0 != NO_ELEM) && (elem1 != NO_ELEM)) {
        MshTopo2D::Edge tmp_edge0 = { opp, save_id2, 0, 0 } ;
        int vtx0 = msh_topo->OppositeNode(elem0,&tmp_edge0) ;

        MshTopo2D::Edge tmp_edge1 = { save_id0, opp, 0, 0 } ;
        int vtx1 = msh_topo->OppositeNode(elem1,&tmp_edge1) ;

        if (vtx0 == vtx1) {
            msh_topo->DeleteTriangle(opp,save_id2,vtx0) ;
            msh_topo->DeleteTriangle(save_id0,opp,vtx0) ;
        }
    }

    // update the triangular mesh topology

    int i ;
    List<SeamElemCache> elem_cache ;

    for (iter0.NewVtx(remove) ; iter0.More() ; ++iter0) {
        if (iter0.CcwElem() != NO_ELEM) {
            SeamElemCache edata ;
            edata.elem = iter0.CcwElem() ;
            msh_topo->GetElemNodes(iter0.CcwElem(),
                    remove,&edata.num_nodes,edata.nodes) ;
            elem_cache.Append(edata) ;
        }
    }
    for (i=0 ; i<elem_cache.Len() ; ++i) {
        msh_topo->DeleteElement(elem_cache[i].num_nodes,
                                elem_cache[i].nodes) ;
    }
    for (i=0 ; i<elem_cache.Len() ; ++i) {
        if ((elem_cache[i].nodes[1] != keep) &&
            (elem_cache[i].nodes[2] != keep)) {
            elem_cache[i].nodes[0] = keep ;
            msh_topo->InsertElement(elem_cache[i].elem,
                                    elem_cache[i].num_nodes,
                                    elem_cache[i].nodes) ;
        }
    }

    // update the quad mesh topology

    elem_cache.Clear() ;
    TopoAdjVtxIterator iter1(quad_topo,remove) ;
    for (iter1.First() ; iter1.More() ; ++iter1) {
        if (iter1.CcwElem() != NO_ELEM) {
            SeamElemCache edata ;
            edata.elem = iter1.CcwElem() ;
            quad_topo->GetElemNodes(iter1.CcwElem(),
                remove,&edata.num_nodes,edata.nodes) ;
            elem_cache.Append(edata) ;
        }
    }
    for (i=0 ; i<elem_cache.Len() ; ++i) {
        quad_topo->DeleteElement(elem_cache[i].num_nodes,
                                 elem_cache[i].nodes) ;
    }
    for (i=0 ; i<elem_cache.Len() ; ++i) {
        elem_cache[i].nodes[0] = keep ;
        quad_topo->InsertElement(elem_cache[i].elem,
                                 elem_cache[i].num_nodes,
                                 elem_cache[i].nodes) ;
    }

    // place the node

    IntNode *node_0 = pnode_table->Get(keep) ;
    node_0->coord = sdata->coord ;

    // update the quad and triangle topology

    SmoothOneAdj(sdata->id2,sdata->id1,sdata->ccw,quad_topo,
                 false,false,edge_list) ;
    SmoothOneAdj(sdata->id0,sdata->id1,sdata->cw,quad_topo,
                 false,false,edge_list) ;

    return(true) ;
}



// %(MshQuad2D::DoTransitionSeam-bool-|-MshEdgeList::QdEdge-|*-MshEdgeList::NewSeamData-|*-MshTopo2D-|*-MshTopo2D-|*-MshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    DoTransitionSeam - do a transition seam 
**
**      bool DoTransitionSeam(
**              MshEdgeList::QdEdge      *edge,
**              MshEdgeList::NewSeamData *sdata,
**              MshTopo2D                   *msh_topo,
**              MshTopo2D                   *quad_topo,
**              MshEdgeList                 *edge_list)
**
**        edge      - (in)  base edge 
**        sdata     - (in)  seam data 
**        msh_topo  - (i/o) triangle mesh 
**        quad_topo - (i/o) quadrilateral mesh 
**        edge_list - (i/o) active boundary 
**
**      Description: This method performs a transition seam update 
**
**      Return Value: true if the update was successful 
**
**
** -- */

bool MshQuad2D::DoTransitionSeam(
                             MshEdgeList::QdEdge *edge,
                             MshEdgeList::NewSeamData *sdata,
                             MshTopo2D *msh_topo,
                             MshTopo2D *quad_topo,
                             MshEdgeList *edge_list)
{
    if (sdata->ratio > 1.0) {

        IntNode *node_0 = pnode_table->Get(sdata->id0) ;
        IntNode *node_1 = pnode_table->Get(sdata->id1) ;
        if ((node_0->type == BOUNDARY) &&
            (node_1->type == BOUNDARY)) return(false) ;

        // delete the quad adjacent to the split edge

        int num_nodes ;
        int quad_nodes[4] ;

        int quad_elem =
            quad_topo->GetCCWElem(sdata->id1,sdata->id0) ;
        quad_topo->GetElemNodes(quad_elem,sdata->id0,
                                &num_nodes,quad_nodes) ;
        if (num_nodes != 4) return(false) ;

        if (edge_list->ContainsEdge(sdata->id1,quad_nodes[2]))
            return(false) ;

        // check to make shure that this element has a valid shape

        IntNode *nd0 = pnode_table->Get(quad_nodes[0]) ;
        IntNode *nd1 = pnode_table->Get(quad_nodes[1]) ;
        IntNode *nd2 = pnode_table->Get(quad_nodes[2]) ;
        IntNode *nd3 = pnode_table->Get(quad_nodes[3]) ;
        double metric = QuadMetric(nd0->coord,nd1->coord,
                                   nd2->coord,nd3->coord) ;
        if (metric <= 0.0) return(false) ;
        quad_topo->DeleteElement(num_nodes,quad_nodes) ;

        // create the new node

        double nx = 0.5 * (node_0->coord[0] + node_1->coord[0]) ;
        double ny = 0.5 * (node_0->coord[1] + node_1->coord[1]) ;
        int new_id = NewNode(nx,ny,INTERIOR,MSH_FLOATING,true) ;

        // add the new quad

        quad_nodes[3] = new_id ;
        quad_topo->InsertElement(quad_elem,4,quad_nodes) ;

        // delete the triangle adjacent to the split edge

        int tri_nodes[3] ;

        int tri_elem =
            msh_topo->GetCCWElem(sdata->id0,sdata->id1) ;
        msh_topo->GetElemNodes(tri_elem,sdata->id1,
                               &num_nodes,tri_nodes) ;

        if (num_nodes > 3) throw QuadConversionError() ;

        msh_topo->DeleteElement(num_nodes,tri_nodes) ;

        // add the new triangles

        int new_tri0[3] = {new_id,quad_nodes[2],sdata->id1} ;
        int new_tri1[3] = {new_id,sdata->id1,tri_nodes[1]} ;
        int new_tri2[3] = {new_id,tri_nodes[1],sdata->id0} ;
        int elem_1 = NewElemNum() ;
        int elem_2 = NewElemNum() ;

        msh_topo->InsertElement(tri_elem,3,new_tri0) ;
        msh_topo->InsertElement(elem_1,3,new_tri1) ;
        msh_topo->InsertElement(elem_2,3,new_tri2) ;

        // update the edge list

        edge_list->UpdateTransSeam(sdata,new_id,quad_nodes[2]) ;

        // now fix up the edge data structure so that we
        // add the correct new quad

        edge->id_0 = quad_nodes[2] ;
        edge->id_1 = sdata->id1 ;
        edge->end_code = 3 ;
        edge->angle_0 = HALF_PI ;
        edge->angle_1 = HALF_PI ;

    } else {

        IntNode *node_1 = pnode_table->Get(sdata->id1) ;
        IntNode *node_2 = pnode_table->Get(sdata->id2) ;
        if ((node_1->type == BOUNDARY) &&
            (node_2->type == BOUNDARY)) return(false) ;

        // find the nodes on the quad that we are going to delete

        int num_nodes ;
        int quad_nodes[4] ;

        int quad_elem =
            quad_topo->GetCCWElem(sdata->id2,sdata->id1) ;
        quad_topo->GetElemNodes(quad_elem,sdata->id1,
                                &num_nodes,quad_nodes) ;
        if (num_nodes != 4) return(false) ;

        if (edge_list->ContainsEdge(quad_nodes[1],sdata->id1))
            return(false) ;

        // check to make shure that this element has a valid shape

        IntNode *nd0 = pnode_table->Get(quad_nodes[0]) ;
        IntNode *nd1 = pnode_table->Get(quad_nodes[1]) ;
        IntNode *nd2 = pnode_table->Get(quad_nodes[2]) ;
        IntNode *nd3 = pnode_table->Get(quad_nodes[3]) ;
        double metric = QuadMetric(nd0->coord,nd1->coord,
                                   nd2->coord,nd3->coord) ;
        if (metric <= 0.0) return(false) ;
        quad_topo->DeleteElement(num_nodes,quad_nodes) ;

        // create the new node

        double nx = 0.5 * (node_1->coord[0] + node_2->coord[0]) ;
        double ny = 0.5 * (node_1->coord[1] + node_2->coord[1]) ;
        int new_id = NewNode(nx,ny,INTERIOR,MSH_FLOATING,true) ;

        // add the new quad

        quad_nodes[0] = new_id ;
        quad_topo->InsertElement(quad_elem,num_nodes,quad_nodes) ;

        // delete the triangle adjacent to the split edge

        int tri_nodes[3] ;

        int tri_elem =
            msh_topo->GetCCWElem(sdata->id1,sdata->id2) ;
        msh_topo->GetElemNodes(tri_elem,sdata->id2,
                               &num_nodes,tri_nodes) ;

        if (num_nodes > 3) throw QuadConversionError() ;

        msh_topo->DeleteElement(num_nodes,tri_nodes) ;

        // add the new triangles

        int new_tri0[3] = {new_id,sdata->id1,quad_nodes[1]} ;
        int new_tri1[3] = {new_id,tri_nodes[1],sdata->id1} ;
        int new_tri2[3] = {new_id,sdata->id2,tri_nodes[1]} ;
        int elem_1 = NewElemNum() ;
        int elem_2 = NewElemNum() ;

        msh_topo->InsertElement(tri_elem,3,new_tri0) ;
        msh_topo->InsertElement(elem_1,3,new_tri1) ;
        msh_topo->InsertElement(elem_2,3,new_tri2) ;

        edge_list->UpdateTransSeam(sdata,new_id,quad_nodes[1]) ;

        // now fix up the edge data structure so that we
        // add the correct new quad

        edge->id_0 = sdata->id1 ;
        edge->id_1 = quad_nodes[1] ;
        edge->end_code = 3 ;
        edge->angle_0 = HALF_PI ;
        edge->angle_1 = HALF_PI ;
    }
    return(true) ;
}




// %(MshQuad2D::DoTransitionSplit-void-|-MshEdgeList::QdEdge-|*-MshEdgeList::NewQuadData-|*-MshTopo2D-|*-MshTopo2D-|*-MshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    DoTransitionSplit - do a transition split 
**
**      void DoTransitionSplit(
**              MshEdgeList::QdEdge      *edge,
**              MshEdgeList::NewQuadData *qdata,
**              MshTopo2D                   *msh_topo,
**              MshTopo2D                   *quad_topo,
**              MshEdgeList                 *edge_list)
**
**        edge      - (in)  base edge 
**        qdata     - (in)  quad data 
**        msh_topo  - (i/o) triangle mesh 
**        quad_topo - (i/o) quadrilateral mesh 
**        edge_list - (i/o) active boundary 
**
**      Description: This method performs a transition split update 
**
**
** -- */

void MshQuad2D::DoTransitionSplit(
                             MshEdgeList::QdEdge *edge,
                             MshEdgeList::NewQuadData *qdata,
                             MshTopo2D *msh_topo,
                             MshTopo2D *quad_topo,
                             MshEdgeList *edge_list)
{
    if (qdata->ratio > 1.0) {

        // delete the quad adjacent to the split edge

        int num_nodes ;
        int quad_nodes[4] ;

        int quad_elem =
            quad_topo->GetCCWElem(qdata->br,qdata->bl) ;
        quad_topo->GetElemNodes(quad_elem,qdata->bl,
                                &num_nodes,quad_nodes) ;

        if (num_nodes > 4) throw QuadConversionError() ;

        quad_topo->DeleteElement(num_nodes,quad_nodes) ;

        // create the two new nodes

        IntNode *node_0 = pnode_table->Get(qdata->br) ;
        IntNode *node_1 = pnode_table->Get(qdata->bl) ;
        double nx = 0.5 * (node_0->coord[0] + node_1->coord[0]) ;
        double ny = 0.5 * (node_0->coord[1] + node_1->coord[1]) ;
        int new_id0 = NewNode(nx,ny,INTERIOR,MSH_FLOATING,true) ;

        node_1 = pnode_table->Get(quad_nodes[1]) ;
        nx = 0.5 * (node_0->coord[0] + node_1->coord[0]) ;
        ny = 0.5 * (node_0->coord[1] + node_1->coord[1]) ;
        int new_id1 = NewNode(nx,ny,INTERIOR,MSH_FLOATING,true) ;

        // add the new quads

        int nqn0[4] = { quad_nodes[0], quad_nodes[1],
                                 new_id1, new_id0 } ;
        quad_topo->InsertElement(quad_elem,4,nqn0) ;

        int nqn1[4] = { quad_nodes[1], quad_nodes[2],
                                 quad_nodes[3], new_id1 } ;
        int qelem1 = NewElemNum() ;
        quad_topo->InsertElement(qelem1,4,nqn1) ;

        // delete the triangle adjacent to the split edge

        int tri_nodes[3] ;

        int tri_elem =
            msh_topo->GetCCWElem(qdata->bl,qdata->br) ;
        msh_topo->GetElemNodes(tri_elem,qdata->br,
                               &num_nodes,tri_nodes) ;

        if (num_nodes > 3) throw QuadConversionError() ;

        msh_topo->DeleteElement(num_nodes,tri_nodes) ;

        // add the new triangles

        int new_tri0[3] = {new_id0,tri_nodes[1],qdata->bl} ;
        int new_tri1[3] = {new_id0,qdata->br,tri_nodes[1]} ;
        int new_tri2[3] = {new_id0,new_id1,qdata->br} ;
        int elem_1 = NewElemNum() ;
        int elem_2 = NewElemNum() ;

        msh_topo->InsertElement(tri_elem,3,new_tri0) ;
        msh_topo->InsertElement(elem_1,3,new_tri1) ;
        msh_topo->InsertElement(elem_2,3,new_tri2) ;

        // update the edge list

        edge_list->UpdateTransSplit(qdata,new_id0,new_id1) ;

        // now fix up the edge data structure so that we
        // add the correct new quad

        edge->id_0 = new_id1 ;
        edge->id_1 = qdata->br ;
        edge->end_code = 3 ;
        edge->angle_0 = HALF_PI ;
        edge->angle_1 = HALF_PI ;

        // update the front list

        qdata->num_front_nodes = 5 ;
        qdata->front_nodes[0] = qdata->bl ;
        qdata->front_nodes[1] = new_id0 ;
        qdata->front_nodes[2] = new_id1 ;
        qdata->front_nodes[3] = qdata->br ;
        qdata->front_nodes[4] = qdata->tr ;

    } else {

        // find the nodes on the quad that we are going to delete

        int num_nodes ;
        int quad_nodes[4] ;

        int quad_elem =
            quad_topo->GetCCWElem(qdata->tr,qdata->br) ;
        quad_topo->GetElemNodes(quad_elem,qdata->br,
                                &num_nodes,quad_nodes) ;

        if (num_nodes > 4) throw QuadConversionError() ;

        quad_topo->DeleteElement(num_nodes,quad_nodes) ;

        // create the new node

        IntNode *node_1 = pnode_table->Get(qdata->br) ;
        IntNode *node_0 = pnode_table->Get(qdata->tr) ;
        double nx = 0.5 * (node_1->coord[0] + node_0->coord[0]) ;
        double ny = 0.5 * (node_1->coord[1] + node_0->coord[1]) ;
        int new_id0 = NewNode(nx,ny,INTERIOR,MSH_FLOATING,true) ;

        node_0 = pnode_table->Get(quad_nodes[2]) ;
        nx = 0.5 * (node_1->coord[0] + node_0->coord[0]) ;
        ny = 0.5 * (node_1->coord[1] + node_0->coord[1]) ;
        int new_id1 = NewNode(nx,ny,INTERIOR,MSH_FLOATING,true) ;

        // add the new quads

        int nqn0[4] = { quad_nodes[0], quad_nodes[1],
                                 quad_nodes[2], new_id1 } ;
        quad_topo->InsertElement(quad_elem,4,nqn0) ;

        int nqn1[4] = { quad_nodes[2], quad_nodes[3],
                                 new_id0, new_id1 } ;
        int qelem1 = NewElemNum() ;
        quad_topo->InsertElement(qelem1,4,nqn1) ;

        // delete the triangle adjacent to the split edge

        int tri_nodes[3] ;

        int tri_elem =
            msh_topo->GetCCWElem(qdata->br,qdata->tr) ;
        msh_topo->GetElemNodes(tri_elem,qdata->tr,
                               &num_nodes,tri_nodes) ;

        if (num_nodes > 3) throw QuadConversionError() ;

        msh_topo->DeleteElement(num_nodes,tri_nodes) ;

        // add the new triangles

        int new_tri0[3] = {new_id0,qdata->tr,tri_nodes[1]} ;
        int new_tri1[3] = {new_id0,tri_nodes[1],qdata->br} ;
        int new_tri2[3] = {new_id0,qdata->br,new_id1} ;
        int elem_1 = NewElemNum() ;
        int elem_2 = NewElemNum() ;

        msh_topo->InsertElement(tri_elem,3,new_tri0) ;
        msh_topo->InsertElement(elem_1,3,new_tri1) ;
        msh_topo->InsertElement(elem_2,3,new_tri2) ;

        // update the edge list

        edge_list->UpdateTransSplit(qdata,new_id0,new_id1) ;

        // now fix up the edge data structure so that we
        // add the correct new quad

        edge->id_0 = qdata->br ;
        edge->id_1 = new_id1 ;
        edge->end_code = 3 ;
        edge->angle_0 = HALF_PI ;
        edge->angle_1 = HALF_PI ;

        // update the front list

        qdata->num_front_nodes = 5 ;
        qdata->front_nodes[0] = qdata->bl ;
        qdata->front_nodes[1] = qdata->br ;
        qdata->front_nodes[2] = new_id1 ;
        qdata->front_nodes[3] = new_id0 ;
        qdata->front_nodes[4] = qdata->tr ;
    }
}




// %(MshQuad2D::DoTemplateSplit-void-|-MshEdgeList::NewQuadData-|*-MshTopo2D-|*-MshTopo2D-|*-MshEdgeList-|*-bool-const|)
/* ++ ----------------------------------------------------------
**
**    DoTemplateSplit - do a template split 
**
**      void DoTemplateSplit(
**              MshEdgeList::NewQuadData *qdata,
**              MshTopo2D                   *msh_topo,
**              MshTopo2D                   *quad_topo,
**              MshEdgeList                 *edge_list,
**              const bool                      base)
**
**        qdata     - (in)  quad data 
**        msh_topo  - (i/o) triangle mesh 
**        quad_topo - (i/o) quadrilateral mesh 
**        edge_list - (i/o) active boundary 
**        base      - (in)  true means split the base 
**
**      Description: This method performs a template split update 
**
**
** -- */

void MshQuad2D::DoTemplateSplit(
                             MshEdgeList::NewQuadData *qdata,
                             MshTopo2D *msh_topo,
                             MshTopo2D *quad_topo,
                             MshEdgeList *edge_list,
                             bool base)
{
    //  This function performs a template transition, which
    //  takes a quad and breaks it into 4 smaller quads, adding
    //  two nodes along the top edge, as shown below:
    //
    //      *--------*      *--*--*--*
    //      |        |      |  |  |  |
    //      |        |      |  |  |  |
    //      |        |      |  *--*  |
    //      |        |      | /    \ |
    //      |        |      |/      \|
    //      *--------*      *--------*
    //
    //  It is assumed that the top edge is attached to the
    //  current boundary.  The triangle that was attached
    //  to the top boundary is split into three in the
    //  obvious way.

    int num_nodes ;
    int quad_nodes[4] ;

    int quad_elem = quad_topo->GetCCWElem(qdata->tr,qdata->tl) ;
    quad_topo->GetElemNodes(quad_elem,qdata->bl,&num_nodes,quad_nodes) ;

    if (num_nodes > 4) throw QuadConversionError() ;

    quad_topo->DeleteElement(num_nodes,quad_nodes) ;

    // create the four new nodes using the bilinear shape
    // functions to determine the new coordinates

    IntNode *br = pnode_table->Get(qdata->br) ;
    IntNode *bl = pnode_table->Get(qdata->bl) ;
    IntNode *tr = pnode_table->Get(qdata->tr) ;
    IntNode *tl = pnode_table->Get(qdata->tl) ;

    double nx, ny ;

    nx = (2.0*tr->coord[0] + tl->coord[0]) / 3.0 ; 
    ny = (2.0*tr->coord[1] + tl->coord[1]) / 3.0 ; 
    int new_id0 = NewNode(nx,ny,INTERIOR,MSH_FLOATING,true) ;

    nx = (tr->coord[0] + 2.0*tl->coord[0]) / 3.0 ; 
    ny = (tr->coord[1] + 2.0*tl->coord[1]) / 3.0 ; 
    int new_id1 = NewNode(nx,ny,INTERIOR,MSH_FLOATING,true) ;

    nx = (2.0*bl->coord[0] + br->coord[0] +
          tr->coord[0] + 2.0*tl->coord[0]) / 6.0 ; 
    ny = (2.0*bl->coord[1] + br->coord[1] +
          tr->coord[1] + 2.0*tl->coord[1]) / 6.0 ; 
    int new_id2 = NewNode(nx,ny,INTERIOR,MSH_FLOATING,true) ;

    nx = (bl->coord[0] + 2.0*br->coord[0] +
          2.0*tr->coord[0] + tl->coord[0]) / 6.0 ; 
    ny = (bl->coord[1] + 2.0*br->coord[1] +
          2.0*tr->coord[1] + tl->coord[1]) / 6.0 ; 
    int new_id3 = NewNode(nx,ny,INTERIOR,MSH_FLOATING,true) ;

    // add the new quads

    int nqn0[4] = { qdata->tl, qdata->bl, new_id2, new_id1 } ;
    int nqn1[4] = { qdata->bl, qdata->br, new_id3, new_id2 } ;
    int nqn2[4] = { qdata->br, qdata->tr, new_id0, new_id3 } ;
    int nqn3[4] = { new_id0, new_id1, new_id2, new_id3 } ;
    int qelem_0 = NewElemNum() ;
    int qelem_1 = NewElemNum() ;
    int qelem_2 = NewElemNum() ;

    quad_topo->InsertElement(quad_elem,4,nqn0) ;
    quad_topo->InsertElement(qelem_0,4,nqn1) ;
    quad_topo->InsertElement(qelem_1,4,nqn2) ;
    quad_topo->InsertElement(qelem_2,4,nqn3) ;

    // delete the triangle adjacent to the split edge

    int tri_nodes[3] ;
    int tri_elem =
            msh_topo->GetCCWElem(qdata->tl,qdata->tr) ;
    msh_topo->GetElemNodes(tri_elem,qdata->tr,
                           &num_nodes,tri_nodes) ;

    if (num_nodes > 3) throw QuadConversionError() ;

    msh_topo->DeleteElement(num_nodes,tri_nodes) ;

    // add the new triangles

    int new_tri0[3] = {tri_nodes[1],new_id0,qdata->tr} ;
    int new_tri1[3] = {tri_nodes[1],new_id1,new_id0} ;
    int new_tri2[3] = {tri_nodes[1],qdata->tl,new_id1} ;
    int elem_1 = NewElemNum() ;
    int elem_2 = NewElemNum() ;

    msh_topo->InsertElement(tri_elem,3,new_tri0) ;
    msh_topo->InsertElement(elem_1,3,new_tri1) ;
    msh_topo->InsertElement(elem_2,3,new_tri2) ;

    // update the edge list

    edge_list->UpdateTempSplit(qdata,new_id0,new_id1,base) ;

    // if we are splitting the top of the new element, update
    // the front list so that things get smoothed properly

    if (!base) {
        qdata->num_front_nodes = 6 ;
        qdata->front_nodes[5] = qdata->front_nodes[3] ;
        qdata->front_nodes[4] = qdata->front_nodes[2] ;
        qdata->front_nodes[3] = new_id0 ;
        qdata->front_nodes[2] = new_id1 ;
    }
}


bool MshQuad2D::QuadsInPolygon(int num,
                        const int *ids,
                        const Vec2D *vts,
                        const MshTopo2D *quad_topo,
                        const MshEdgeList *edge_list) const
{
    int i ;

    // here we set up data structures to look for mate nodes.  This
    // is an after the fact fix so the logic is a bit strange.  We
    // assume that the number of verts is <= 6.  If this is not
    // the case then we create more space

    SmallSet<int,2> *sets[6] = {0,0,0,0,0,0} ;
    SmallSet<int,2> **setptr = &sets[0] ;
    bool have_sets = false ;

    if (num > 6) {
        setptr = new SmallSet<int,2> *[num] ;
    }

    // determine a tolerance and store the max and min
    // extents

    Vec2D pmin = vts[0] ;
    Vec2D pmax = vts[0] ;

    double tol = 0 ;
    for (i=0 ; i<num ; ++i) {
        int j = (i+1) % num ;
        tol += (vts[j]-vts[i]).Magnitude() ;
        pmin = Minv(vts[i],pmin) ;
        pmax = Maxv(vts[i],pmax) ;
        setptr[i] = MateTable->Get(ids[i]) ;
        if (setptr[i] != 0) {
            have_sets = true ;
        }
    }
    tol = 0.001 * tol/double(num) ;

    // loop through all verts in the mesh.  If they are not
    // part of the polygon check to see if they are inside.

//    TopoVtxIterator iter(*quad_topo) ;

    EdgeListNodeIterator iter(edge_list) ;
    for (iter.First() ; iter.More() ; ++iter) {

//        if (quad_topo->NumAdjElems(*iter) == 0) continue ;
        if (!quad_topo->HasVtx(*iter) ||
            quad_topo->NumAdjElems(*iter) == 0) continue ;

        bool check_this = true ;
        for (i=0 ; i<num ; ++i) {
            if (*iter == ids[i]) {
                check_this = false ;
                break ;
            }
            if (have_sets && (setptr[i] != 0)) {
                if (setptr[i]->HasElement(*iter)) {
                    check_this = false ;
                    break ;
                }
            }
        }

        if (!check_this) continue ;

        Vec2D pt = (pnode_table->Get(*iter))->coord ;
        if ((pt.x() < pmin.x()) || (pt.x() > pmax.x()) ||
            (pt.y() < pmin.y()) || (pt.y() > pmax.y())) continue ;

        if (IsPointInPoly (pt,num,vts,tol)) {
            if (setptr != &sets[0]) delete [] setptr ;
            return(true) ;
        }
    }

    if (setptr != &sets[0]) delete [] setptr ;
    return(false) ;
}


// %(MshQuad2D::RenumberNodes-void-|)
/* ++ ----------------------------------------------------------
**
**    RenumberNodes - renumber nodes so that unused numbers are removed 
**
**      void RenumberNodes()
**
**      Description: This method renumbers the generated nodes so that 
**          nodes generated during triangulation but not used when 
**          generating quads are eliminated. 
**
**
** -- */

void MshQuad2D::RenumberNodes()
{
    Dict<int,int> node_map ; 

    Dict<int,MshElement2D>::DictIterator eiter(pelem_table) ;

    int j ;

    // first build a hash table that maps from the old number
    // to the new number for the nodes that are actually used.

    for (eiter.First() ; eiter.More() ; ++eiter) {
        MshElement2D& elem = eiter.Entry() ;
        for (j=0 ; j<elem.num_nodes ; ++j) {
            if (elem.nodes[j] > MaxId) {
                if (node_map.Get(elem.nodes[j]) == 0) {
                    node_map.Store(elem.nodes[j],StartIdSave) ;
                    StartIdSave++ ;
                }
            }
        }
    }

    // now loop again and update the node numbers in all
    // elements

    for (eiter.First() ; eiter.More() ; ++eiter) {
        MshElement2D& elem = eiter.Entry() ;
        for (j=0 ; j<elem.num_nodes ; ++j) {
            if (elem.nodes[j] > MaxId) {
                elem.nodes[j] = 
                    *(node_map.Get(elem.nodes[j])) ;
            }
        }
    }

    // build a new node list

    Dict<int,IntNode> *tmp_node = 
        new Dict<int,IntNode>() ;

    Dict<int,IntNode>::DictIterator iter(pnode_table) ;

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.Entry().id <= MaxId) {
            tmp_node->Store(iter.Entry().id,iter.Entry()) ;
        } else {
            int *nid = node_map.Get(iter.Entry().id) ;
            if (nid != 0) { 
                IntNode& tmp = iter.Entry() ;
                tmp.id = *nid ;
                tmp_node->Store(tmp.id,tmp) ;
            }
        }
    }
    delete pnode_table ;
    pnode_table = tmp_node ;
}




// %(MshQuad2D::IntersectLines-Vec2D-|^const-Vec2D-|-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    IntersectLines - find the intersection point of two lines 
**
**      Vec2D IntersectLines(
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
**      Description: This method finds the intersection point of two 
**          line segments. 
**
**      Return Value: the intersection coordinates 
**
**
** -- */

Vec2D MshQuad2D::IntersectLines(
         Vec2D i1,
         Vec2D i2,
         Vec2D j1,
         Vec2D j2) const
{
    Vec2D intsc ;

    /* This is the assumed node configuration:
    
           \I1   /J2
            \   /
             \ /
              \
             / \
            /   \
           /J1   \I2
    */

    Vec2D bi = i2 - i1 ;
    Vec2D bj = j2 - j1 ;

    // check for the degenerate case of parallel lines

    double denom = bi[0]*bj[1] - bi[1]*bj[0] ;
    double leni = bi.Magnitude() ;
    double lenj = bj.Magnitude() ;
    double len = leni < lenj ? leni : lenj ;

    if (fabs(denom) < (0.000001*len)) {

        // if the lines are parallel or coincident, just average all the
        // nodes.  This is foolish for some applications, but makes sense
        // for others.

        intsc = 0.25 * (i1 + i2 + j1 + j2) ;
    } else {
        double par = -(bj[0]*(j1[1]-i1[1]) + bj[1]*(i1[0]-j1[0])) / denom ;
        intsc = i1 + par * (i2 - i1) ;
    }
    return(intsc) ;
}

} // namespace
