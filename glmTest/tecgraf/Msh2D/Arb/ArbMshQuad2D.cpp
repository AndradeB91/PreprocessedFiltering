//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000
//   All rights reserved
//
// Revision -
//   $Revision: 1.59 $  $Date: 2002/09/12 18:32:05 $  $Author: wash $
//

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "ArbMshRegion2D.hpp"
#include "ArbMshEdgeList.hpp"
#include "ArbMshTopo2D.hpp"
#include "ArbArray.hpp"
#include "ArbCoord2D.hpp"
#include "ArbMshSmooth2D.hpp"
#include "ArbSet.cpp"

#ifdef MEMDEBUG
#include "MemDbg.hpp"
#define new new(__FILE__,__LINE__)
#endif

// #define DEBUG_LOG
#ifdef DEBUG_LOG
FILE *lfd ;
#endif


// static int CompareElement(const ArbMshElement2D &e1,
//                           const ArbMshElement2D &e2)
// {
//     if (e1.elem_id > e2.elem_id) {
//         return(1) ;
//     } else if (e1.elem_id < e2.elem_id) {
//         return(-1) ;
//     }
//     return(0) ;
// }


// %(CArbMshRegion2D::GenerateQuads-void-|)
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
#define QUAD_EDGE_NORMAL_TOLERANCE 3.926990817    // 5*pi/4

#define MAX_ALLOWABLE_ANGLE 2.356194490


#define SMOOTH_BLACKER_TOLERANCE 2.5
#define SMOOTH_OWEN_TOLERANCE 50.0

#define TRANSITION_SEAM_TOLERANCE 2.5
#define TRANSITION_SPLIT_TOLERANCE 2.5
#define TRANSITION_TEMPLATE_TOLERANCE 3.0

#define TRIANGLE_SPLIT_TOLERANCE 3.0

#define TRANSITION_ADJACENT_TOLERANCE 20.0
#define TRANSITION_NEAR_BOUNDARY_TOLERANCE 2.0

#define SEAM_MASK 12
#define NO_SEAM_MASK 3

#define TWO_PI  6.283185307
#define HALF_PI 1.570796327

int print_num = 1000000 ;
int quit_num  = 1000000 ;
bool do_labels = true ;

static int debug_num = 0 ;


// GenQuadFromMesh
//////////////////////////////////////////////////////////////
void CArbMshRegion2D::GenQuadFromMesh ()
{
  GenerateQuads() ;
  RenumberNodes() ;
}

// GenerateQuads
//////////////////////////////////////////////////////////////
int CArbMshRegion2D::GenerateQuads ()
{
  CArbMshTopo2D msh_topo ;
  CArbMshTopo2D quad_topo ;
  CArbMshEdgeList edge_list(pnode_table) ;

#ifdef DEBUG_LOG
  lfd = fopen("d:\\temp\\debugQuad.log","w") ;
#endif

  // Initialize the data structures
  //////////////////////////////////////////////////////////////
  if (!InitializeQuadGen (&msh_topo,&edge_list))
   return 0;

  // process the edges
  //////////////////////////////////////////////////////////////
  CArbMshEdgeList::ArbQdEdge *qd_edge ;
  bool saved_flag = false ;
  int saved_0 = 0, saved_1 = 0 ;
  bool modified = true ;
  int stagnate_check[2] ;
  int stagnate_count = 0 ;

  int num = 0 ;

  while ((qd_edge = edge_list.GetNextEdge()) != 0)
  {
    debug_num = num ;
    //printf("num %i\n",num);
    if (num == nelem*3)
     return 0;
    bool transition_flg = false ;

    if ((int)num >= print_num)
    {
        DoPrint(&msh_topo) ;
        DoPrint(&quad_topo) ;
        printf("a Start Step %d\n",num) ;
        printf("n\n") ;
        fflush(stdout) ;
    }

#ifdef DEBUG_LOG
fprintf(lfd,"\nNum : Num = %d, Edge = %5d %5d, NumElem = %3d, %1d %1d\n",num,
      qd_edge->id_0,qd_edge->id_1,quad_topo.NumElements(),
      (modified?1:0),stagnate_count) ;
fflush(lfd) ;
#endif

    if (modified)
    {
      stagnate_check[0] = qd_edge->id_0 ;
      stagnate_check[1] = qd_edge->id_1 ;
      modified = false ;
      stagnate_count = 0 ;
    }
    else
    {
      if ((stagnate_check[0] == qd_edge->id_0) && (stagnate_check[1] == qd_edge->id_1))
      {
        // it looks like we may be stagnating.  Look through
        // the edge list and see if any edges have not been
        // touched since the last sweep

        bool status = edge_list.StagnationSweep() ;
        if (status)
        {
          if (stagnate_count > 1) break ;  // stagnate
          ++stagnate_count ;
        }
      }
    }

    // check for the special case of closing a triangle
    // or a quad or a pentagon
    if ((qd_edge->end_code & NO_SEAM_MASK) == 3)
    {
      if (CloseSimplePoly(&num,qd_edge,&msh_topo, &quad_topo,&edge_list))
      {
        saved_flag = false ;
        modified = true ;
        continue ;
      }
    }

    //printf("1\n");
    // check to see if we need to do any "seaming"
    if (qd_edge->angle_0 < SEAM_ANGLE_BIG_TOLERANCE)
    {
      if (CloseSeam(true,&num,qd_edge,&msh_topo,&quad_topo, &edge_list,&transition_flg,&modified))
      {
        saved_flag = false ;
        continue ;
      }
      else
      {
        // see if we should form a triangle or other
        // simple polygon close a simple polygon
        if (qd_edge->level == 0)
        {
          int top = edge_list.GetCWNode (qd_edge->id_0,qd_edge->id_1);
          if (edge_list.GetEdgeLevel (top,qd_edge->id_0) == 0)
          {
            CArbMshEdgeList::ArbNewQuadData qdata ;
            int top = edge_list.GetCWNode (qd_edge->id_0,qd_edge->id_1);
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
            if (!RecoverTopEdge(qdata.tr,qdata.bl, &msh_topo,&edge_list))
            {
              edge_list.PushToBack(qd_edge) ;
              continue ;
            }
            edge_list.UpdateQuad(&qdata) ;
            if (!edge_list.UpdateGeom(&qdata))
             return 0;
            int boundary[3] ;
            boundary[0] = qd_edge->id_0 ;
            boundary[1] = qd_edge->id_1 ;
            boundary[2] = top ;
            CArbMshTopo2D::ArbEdge edge ;
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
        if (CloseSimplePoly(&num,qd_edge,&msh_topo, &quad_topo,&edge_list))
        {
          saved_flag = false ;
          modified = true ;
          continue ;
        }
      }
    }

    //printf("2\n");
    if (qd_edge->angle_1 < SEAM_ANGLE_BIG_TOLERANCE)
    {
      if (CloseSeam(true,&num,qd_edge,&msh_topo,&quad_topo,&edge_list,&transition_flg,&modified))
      {
          saved_flag = false ;
          continue ;
      }
      else
      {
        // see if we should form a triangle or other
        // simple polygon close a simple polygon

        if (qd_edge->level == 0)
        {
          int top = edge_list.GetCCWNode (qd_edge->id_0,qd_edge->id_1) ;
          if (edge_list.GetEdgeLevel(qd_edge->id_1,top) == 0)
          {
            CArbMshEdgeList::ArbNewQuadData qdata ;
            int top = edge_list.GetCCWNode (qd_edge->id_0,qd_edge->id_1) ;
            qdata.qcase = 6 ;
            qdata.tl = top ;
            qdata.bl = qd_edge->id_0 ;
            qdata.br = qd_edge->id_1 ;
            qdata.tr = top ;
            qdata.front_nodes[0] = edge_list.GetCWNode(
                qd_edge->id_0,qd_edge->id_1) ;
            if (qdata.front_nodes[0] == -1)
             return 0;
            qdata.front_nodes[1] = qdata.bl ;
            qdata.front_nodes[2] = qdata.tr ;
            qdata.front_nodes[3] = edge_list.GetCCWNode(
                qd_edge->id_1,top) ;
            if (qdata.front_nodes[3] == -1)
             return 0;
            qdata.num_front_nodes = 4 ;
            if (!RecoverTopEdge (qdata.tr,qdata.bl, &msh_topo,&edge_list))
            {
              edge_list.PushToBack(qd_edge) ;
              continue ;
            }
            edge_list.UpdateQuad(&qdata) ;
            if (!edge_list.UpdateGeom(&qdata))
             return 0;
            int boundary[3] ;
            boundary[0] = qd_edge->id_0 ;
            boundary[1] = qd_edge->id_1 ;
            boundary[2] = top ;
            CArbMshTopo2D::ArbEdge edge ;
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
        if (CloseSimplePoly (&num,qd_edge,&msh_topo, &quad_topo,&edge_list))
        {
          saved_flag = false ;
          modified = true ;
          continue ;
        }
      }
    }

    // if we get here we are doing a quad, so get the
    // appropriate information

    CArbMshEdgeList::ArbNewQuadData qdata ;

    //printf("3\n");

    if (!FindQuadData(qd_edge,&qdata,&edge_list,&msh_topo))
     return 0;

    //printf("4\n");

#ifdef DEBUG_LOG
if (qdata.qcase != 6) {
  fprintf(lfd,"QuadData : %1d %5d %5d %5d %5d\n",
      qdata.qcase,qdata.bl,qdata.br,qdata.tr,qdata.tl) ;
  fprintf(lfd,"             %5d %5d %5d %5d\n",
      qdata.obl,qdata.obr,qdata.otr,qdata.otl) ;
  fprintf(lfd,"           %1d : ",qdata.num_front_nodes) ;
  for (int di=0 ; di<qdata.num_front_nodes ; ++di) {
      fprintf(lfd," %5d",qdata.front_nodes[di]) ;
  }
  fprintf(lfd,"\n") ;
  fprintf(lfd,"           %5d %g\n",qdata.edge_level,qdata.ratio) ;
  fflush(lfd) ;
}
#endif

    // check for the case the right and left nodes are
    // the same.  If this happens, we push the edge back
    // on the edge list

    if ((qdata.qcase == 6) || (qdata.qcase == -1))
    {
      edge_list.PushToBack(qd_edge) ;
      ++num ;
      continue ;
    }

    // check to see is we want to do a transition split
    int op = TransitionSplit(&num,&qdata,qd_edge,&msh_topo,&quad_topo, &edge_list,&transition_flg);
    if (op==1)
    {
      saved_flag = false ;
      modified = true ;
      continue ;
    }
    if (op == -1)
     return 0;

    // if we are not currently doing a transion, compute
    // a metric for the proposed element.  If the metric is very
    // negative (invalid), and we have not seen this edge
    // before, then go on to a new edge (saving information
    // about this edge).  If we have seen this edge before
    // then go ahead and use it, and hope that the smoothing
    // process forms a valid quad.

    ArbIntNode *bl = pnode_table->Fetch(qdata.bl);
    ArbIntNode *br = pnode_table->Fetch(qdata.br);
    ArbIntNode *tr = pnode_table->Fetch(qdata.tr);
    ArbIntNode *tl = pnode_table->Fetch(qdata.tl);

    if (bl == 0)
     return 0;

    if (br == 0)
     return 0;

    if (tr == 0)
     return 0;

    if (tl == 0)
     return 0;

    bool constrained = ((bl->motion == ARB_FIXED) && (br->motion == ARB_FIXED) &&
                        (tr->motion == ARB_FIXED) && (tl->motion == ARB_FIXED) );

    if (!transition_flg)
    {
      double metric = QuadMetricAreaOnly(bl->coord,br->coord, tr->coord,tl->coord) ;

      // check to see if any of the existing quads
      // fall inside the proposed element
      int ids[4] = {qdata.bl,qdata.br,qdata.tr,qdata.tl} ;
      CArbCoord2D vts[4] ;
      vts[0] = bl->coord ;
      vts[1] = br->coord ;
      vts[2] = tr->coord ;
      vts[3] = tl->coord ;

      if ((metric <= 0.0001) || QuadsInPolygon(4,ids,vts,&quad_topo,&edge_list))
      {
        edge_list.PushToBack(qd_edge) ;
        if ((saved_0 != qd_edge->id_0) || (saved_1 != qd_edge->id_1))
        {
          if (!saved_flag)
          {
            saved_0 = qd_edge->id_0 ;
            saved_1 = qd_edge->id_1 ;
            saved_flag = true ;
          }
        }
        continue;
      }
    }


    // now do a series of swaps to recover the top edge
    // if possible.

    if (qdata.qcase == 3)
    {
      if (!RecoverTopEdge(qdata.tr,qdata.tl, &msh_topo,&edge_list))
      {
        edge_list.PushToBack(qd_edge) ;
        continue ;
      }
    }
    else
    {
      if (!RecoverTopEdge(qdata.otr,qdata.otl, &msh_topo,&edge_list))
      {
        edge_list.PushToBack(qd_edge) ;
        continue ;
      }
    }

    if ((int)debug_num >= print_num)
    {
      DoPrint(&msh_topo) ;
      DoPrint(&quad_topo) ;
      printf("a After recover top %d\n",num) ;
      printf("n\n") ;
      fflush(stdout) ;
    }

    // clear out the triangles that currently cover the
    // new quad

    int boundary[4] ;
    boundary[0] = qdata.bl ;
    boundary[1] = qdata.br ;
    boundary[2] = qdata.tr ;
    boundary[3] = qdata.tl ;

    CArbMshTopo2D::ArbEdge edge ;
    edge.nd0 = qdata.bl ;
    edge.nd1 = qdata.br ;
    edge.elem0 = msh_topo.GetCCWElem(edge.nd0,edge.nd1) ;
    edge.elem1 = NO_ELEM ;
    ClearRegion(4,boundary,&msh_topo,&edge) ;

    // update the edge list to add the quad

#ifdef DEBUG_LOG
fprintf(lfd,"Adding Quad : %5d %5d %5d %5d\n",
      boundary[0],boundary[1],boundary[2],boundary[3]) ;
fflush(lfd) ;
#endif

    edge_list.UpdateQuad(&qdata);
    if (qdata.edge_level == -1)
     return 0;
    saved_flag = false ;
    modified = true ;

    // add this quad to the topology

    quad_topo.InsertElement(num,4,boundary);
    if (qdata.qcase == 4)
    {
      num++ ;
      continue ;
    }

    // check for template transition if this new quad is
    // too distorted or the aspect ratio of the adjacent
    // triangle has a very poor aspect ratio

    bool template_tr = false ;
    if ((qdata.qcase == 3) && (qdata.edge_level > 1) && !transition_flg &&
        (qdata.top_template_ratio > TRANSITION_TEMPLATE_TOLERANCE))
    {
        if (!DoTemplateSplit(&qdata,&msh_topo,&quad_topo,&edge_list,false))
         return 0;
        template_tr = true ;
    }

    // do smoothing of front nodes (three times through is abitrary)

    for (int j=0 ; j<3 ; ++j)
    {
      SmoothFront (qdata.qcase, qdata.num_front_nodes, qdata.front_nodes,
                    &msh_topo,&quad_topo,&edge_list);
      SmoothAdjacent (qdata.num_front_nodes, qdata.front_nodes,&msh_topo,
                      true,true,&edge_list) ;
      SmoothAdjacent (qdata.num_front_nodes, qdata.front_nodes,&quad_topo,
                      false,false,&edge_list) ;
    }

    SmoothFront(qdata.qcase, qdata.num_front_nodes, qdata.front_nodes,
                &msh_topo,&quad_topo,&edge_list) ;

    if ((qdata.qcase == 3) && !transition_flg && !template_tr && !constrained)
    {
      if ((AdjAspect(&qdata,&msh_topo) > TRANSITION_ADJACENT_TOLERANCE) ||
          (AdjNearBoundary(&qdata,&msh_topo,&edge_list) > TRANSITION_NEAR_BOUNDARY_TOLERANCE) ||
          CornerAspect(&qdata,&msh_topo) )
      {
        if (!DoTemplateSplit(&qdata,&msh_topo,&quad_topo, &edge_list,false))
         return 0;

        for (int j=0 ; j<3 ; ++j)
        {
          SmoothFront (qdata.qcase, qdata.num_front_nodes, qdata.front_nodes,
                        &msh_topo,&quad_topo,&edge_list) ;
          SmoothAdjacent (qdata.num_front_nodes, qdata.front_nodes,&msh_topo,
                          true,true,&edge_list) ;
          SmoothAdjacent (qdata.num_front_nodes, qdata.front_nodes,&quad_topo,
                          false,false,&edge_list) ;
        }

        SmoothFront (qdata.qcase, qdata.num_front_nodes,qdata.front_nodes,
                      &msh_topo,&quad_topo,&edge_list) ;
      }
    }

    // update the classification of edges on the boundary

    if (!edge_list.UpdateGeom(&qdata))
     return 0;
    for (int ii=0 ; ii<BdryVtxCache->NumEntries() ; ++ii)
      edge_list.UpdateBdryGeom((*BdryVtxCache)[ii]) ;

    BdryVtxCache->Clear() ;

    // get a list of the edges and display

    ++num ;
  }

  // When we get here we are all done adding quads.
  // if there are any triangles left in the  triangle mesh
  // add these to the quad mesh
  /////////////////////////////////////////////////////////

  CArbTopoElemIterator titer(&msh_topo) ;
  for (titer.First() ; titer.More() ; ++titer)
  {
    int nnode,nodes[3] ;
    msh_topo.GetElemNodes(*titer,&nnode,nodes) ;
    quad_topo.InsertTriangle(num,nodes[0],nodes[1],nodes[2]) ;
    ++num ;
    if (num > nelem*3)
     return 0;
  }

  // Now we do topological improvements.
  /////////////////////////////////////////////////////////
  if (DebugDisplayFlags & RmshQuadBeforeCleanup)
      DisplayTopo("Quad Mesh Before Cleanup",&quad_topo) ;

  if (DoCleanup) QuadCleanup(&quad_topo) ;

  // now we replace the triangular elements with the quad
  // elements and do global smoothing
  /////////////////////////////////////////////////////////
  UpdateForQuadMesh(&quad_topo) ;

  // smooth the mesh
  /////////////////////////////////////////////////////////
  if (DebugDisplayFlags & RmshQuadBeforeSmooth)
      DisplayMesh("Quad Mesh Before Smooth") ;

  if (DoSmoothNodes)
  {
    CArbMshSmooth2D smooth(pnode_table,pelem_table) ;
    if (WinslowSmoothing)
    {
     if (!smooth.SmoothNodesWinslow())
      return 0;
    }
    else
    {
     smooth.SmoothNodesLaplace();
    }
    if (DebugDisplayFlags & RmshQuadAfterSmooth)
        DisplayMesh("Quad Mesh After Wins") ;
    smooth.SmoothNodesConsLaplace() ;
  }

  if (DebugDisplayFlags & RmshQuadAfterSmooth)
      DisplayMesh("Quad Mesh After Smooth") ;

  // do angle checks if required
  /////////////////////////////////////////////////////////
  if (DoQuadAngleChecks) QuadAngleChecks() ;

  return 1;
}


double CArbMshRegion2D::AdjAspect(
                             CArbMshEdgeList::ArbNewQuadData *qdata,
                             CArbMshTopo2D *msh_topo) const
{
    // find the coords of the top of the quadrilateral

    ArbIntNode *tr = pnode_table->Fetch(qdata->tr) ;
    ArbIntNode *tl = pnode_table->Fetch(qdata->tl) ;

    // get the coordinates of the opposite node on the adjacent
    // triangle element

    int num_nodes,tri_nodes[3] ;
    int tri_elem =
            msh_topo->GetCCWElem(qdata->tl,qdata->tr) ;
    msh_topo->GetElemNodes(tri_elem,qdata->tr,
                           &num_nodes,tri_nodes) ;

    ArbIntNode *opp = pnode_table->Fetch(tri_nodes[1]) ;

    // compute the height of the triangle and return the ratio
    // of the base to the height

    double base_len ;
    double height = DistSqr(opp->coord,tl->coord,
                            tr->coord,&base_len) ;

    return(base_len/sqrt(height)) ;
}


double CArbMshRegion2D::AdjNearBoundary(
                             CArbMshEdgeList::ArbNewQuadData *qdata,
                             CArbMshTopo2D *msh_topo,
                             CArbMshEdgeList *edge_list) const
{
    // find the coords of the top of the quadrilateral

    ArbIntNode *tr = pnode_table->Fetch(qdata->tr) ;
    ArbIntNode *tl = pnode_table->Fetch(qdata->tl) ;

    // get the coordinates of the opposite node on the adjacent
    // triangle element

    int num_nodes,tri_nodes[3] ;
    int tri_elem =
            msh_topo->GetCCWElem(qdata->tl,qdata->tr) ;
    msh_topo->GetElemNodes(tri_elem,qdata->tr,
                           &num_nodes,tri_nodes) ;

    // check to see if this opposite node is part of a boundary

    if (!edge_list->ContainsNode(tri_nodes[1])) return(0.0) ;

    ArbIntNode *opp = pnode_table->Fetch(tri_nodes[1]) ;

    // compute the height of the triangle

    double base_len ;
    double height = DistSqr(opp->coord,tl->coord,
                            tr->coord,&base_len) ;

    if (sqrt(height) > base_len) return(0.0) ;

    // find the average edge length at this node

    double char_len = edge_list->GetCharNodeLength(tri_nodes[1]) ;

    return(base_len/char_len) ;
}


double CArbMshRegion2D::CornerAspect(
                             CArbMshEdgeList::ArbNewQuadData *qdata,
                             CArbMshTopo2D * /*msh_topo*/) const
{
    ArbIntNode *br = pnode_table->Fetch(qdata->br) ;

    if (br->motion != ARB_FIXED) return(false) ;

    ArbIntNode *bl = pnode_table->Fetch(qdata->bl) ;

    if (bl->motion != ARB_FIXED) return(false) ;

    ArbIntNode *tr = pnode_table->Fetch(qdata->tr) ;
    ArbIntNode *tl = pnode_table->Fetch(qdata->tl) ;

    double rlen = (tr->coord - br->coord).Magnitude() ;
    double llen = (tl->coord - bl->coord).Magnitude() ;

    // check to see if we need to do a transition

    double ratio = rlen/llen ;
    if ((ratio > TRANSITION_NEAR_BOUNDARY_TOLERANCE) ||
        (1.0/ratio > TRANSITION_NEAR_BOUNDARY_TOLERANCE)) return(true) ;
    return(false) ;
}



// %(CArbMshRegion2D::CloseSimplePoly-bool-|-int-|*-CArbMshEdgeList::ArbQdEdge-|*-CArbMshTopo2D-|*-CArbMshTopo2D-|*-CArbMshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    CloseSimplePoly - turn a simple polygonal region into an element
**
**      bool CloseSimplePoly(
**              int                        *current_num,
**              CArbMshEdgeList::ArbQdEdge *qd_edge,
**              CArbMshTopo2D              *msh_topo,
**              CArbMshTopo2D              *quad_topo,
**              CArbMshEdgeList            *edge_list)
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

bool CArbMshRegion2D::CloseSimplePoly(
                     int *current_num,
                     CArbMshEdgeList::ArbQdEdge *qd_edge,
                     CArbMshTopo2D *msh_topo,
                     CArbMshTopo2D *quad_topo,
                     CArbMshEdgeList *edge_list)
{
    int right_node, left_node ;
    right_node = edge_list->GetCCWNode(qd_edge->id_0,
                                      qd_edge->id_1) ;
    left_node = edge_list->GetCWNode(qd_edge->id_0,
                                     qd_edge->id_1) ;
    if (right_node == -1)
     return false;
    if (left_node == -1)
     return false;

    // if the left and right nodes are the same, then the
    // remaining void forms a triangle, so make it an element

    if (right_node == left_node) {

        // build a list of the vertex coordinates so that we
        // can see if any quads are inside

        CArbCoord2D vts[3] ;
        vts[0] = (pnode_table->Fetch(qd_edge->id_0))->coord ;
        vts[1] = (pnode_table->Fetch(qd_edge->id_1))->coord ;
        vts[2] = (pnode_table->Fetch(right_node))->coord ;
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

        CArbMshTopo2D::ArbEdge edge ;
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

            CArbMshEdgeList::ArbNewQuadData qdata ;
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

            ArbIntNode *nodes[4] ;
            nodes[0] = pnode_table->Fetch(qdata.bl) ;
            nodes[1] = pnode_table->Fetch(qdata.br) ;
            nodes[2] = pnode_table->Fetch(qdata.tr) ;
            nodes[3] = pnode_table->Fetch(qdata.tl) ;

            // check to make sure that there is nothing inside
            // the quad

            CArbCoord2D vts[4] ;
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

            CArbMshTopo2D::ArbEdge edge ;
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
                if (nodes[i]->motion == ARB_FIXED) ++num_fixed ;

            if (num_fixed == 3) {
                for (i=0 ; i<4 ; ++i) {
                    j = (i+1) % 4 ;
                    k = (i+2) % 4 ;
                    if ((nodes[i]->motion == ARB_FIXED) &&
                        (nodes[j]->motion == ARB_FIXED) &&
                        (nodes[k]->motion == ARB_FIXED)) break ;
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

            CArbCoord2D vts[5] ;
            int ids[5] ;
            ArbIntNode *nodes[5] ;
            for (int i=0 ; i<5 ; ++i) {
                nodes[i] = pnode_table->Fetch(poly_nodes[i]) ;
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
            if (min_metric < -0.2) return(false) ;
            if (min_metric < 0.0) {
                for (j=0 ; j<5 ; ++j) {
                     if (nodes[j]->motion == ARB_FIXED) return(false) ;
                }
            }

            // clear the region

            CArbMshTopo2D::ArbEdge edge ;
            edge.nd0 = poly_nodes[0] ;
            edge.nd1 = poly_nodes[1] ;
            edge.elem0 = msh_topo->GetCCWElem(edge.nd0,edge.nd1) ;
            edge.elem1 = NO_ELEM ;
            ClearRegion(5,poly_nodes,msh_topo,&edge) ;

            // add the new elements

            CArbMshEdgeList::ArbNewQuadData qdata ;
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

            ArbIntNode *nodes[6] ;
            int ids[6] ;
            CArbCoord2D vts[6] ;
            for (int i=0 ; i<6 ; ++i) {
                nodes[i] = pnode_table->Fetch(poly_nodes[i]) ;
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
            if (min_metric < 0) return(false) ;

            // clear the region

            CArbMshTopo2D::ArbEdge edge ;
            edge.nd0 = poly_nodes[0] ;
            edge.nd1 = poly_nodes[1] ;
            edge.elem0 = msh_topo->GetCCWElem(edge.nd0,edge.nd1) ;
            edge.elem1 = NO_ELEM ;
            ClearRegion(6,poly_nodes,msh_topo,&edge) ;

            // add the new elements

            CArbMshEdgeList::ArbNewQuadData qdata ;
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




// %(CArbMshRegion2D::CloseSeam-bool-|-bool-|-int-|*-CArbMshEdgeList::ArbQdEdge-|*-CArbMshTopo2D-|*-CArbMshTopo2D-|*-CArbMshEdgeList-|*-bool-|*)
/* ++ ----------------------------------------------------------
**
**    CloseSeam - performs a seam closing
**
**      bool CloseSeam(
**              bool                       left,
**              int                        *current_num,
**              CArbMshEdgeList::ArbQdEdge *qd_edge,
**              CArbMshTopo2D              *msh_topo,
**              CArbMshTopo2D              *quad_topo,
**              CArbMshEdgeList            *edge_list,
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

bool CArbMshRegion2D::CloseSeam(bool left,
                     int *current_num,
                     CArbMshEdgeList::ArbQdEdge *qd_edge,
                     CArbMshTopo2D *msh_topo,
                     CArbMshTopo2D *quad_topo,
                     CArbMshEdgeList *edge_list,
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

            CArbMshEdgeList::ArbNewSeamData sdata ;
            bool status = FindSeamData(qd_edge,&sdata,edge_list,
                                       msh_topo,quad_topo) ;
            *modified_flg = false ;
            if (!status && (sdata.cw != sdata.ccw)) return(false) ;

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

                CArbMshTopo2D::ArbEdge edge ;
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
                CArbMshEdgeList::ArbNewQuadData qdata ;
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

                CArbMshTopo2D::ArbEdge edge ;
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
                CArbMshEdgeList::ArbNewQuadData tqdata ;
                if (sdata.template_cw) {
                    int elem = !quad_topo->HasVtx(sdata.id1) ? NO_ELEM :
                               quad_topo->GetCCWElem(sdata.id1,sdata.id0) ;
                    if (elem != NO_ELEM) {
                        int nnum, nodes[4] ;
                        quad_topo->GetElemNodes(elem,sdata.id1,&nnum,nodes) ;

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

            } else if ((sdata.ratio > TRANSITION_SEAM_TOLERANCE) ||
                       (sdata.ratio < 1.0/TRANSITION_SEAM_TOLERANCE)) {
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




// %(CArbMshRegion2D::TransitionSplit-bool-|-int-|*-CArbMshEdgeList::ArbNewQuadData-|*-CArbMshEdgeList::ArbQdEdge-|*-CArbMshTopo2D-|*-CArbMshTopo2D-|*-CArbMshEdgeList-|*-bool-|*)
/* ++ ----------------------------------------------------------
**
**    TransitionSplit - performs a transition split
**
**      bool TransitionSplit(
**              int                             *current_num,
**              CArbMshEdgeList::ArbNewQuadData *qdata,
**              CArbMshEdgeList::ArbQdEdge      *qd_edge,
**              CArbMshTopo2D                   *msh_topo,
**              CArbMshTopo2D                   *quad_topo,
**              CArbMshEdgeList                 *edge_list,
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

int CArbMshRegion2D::TransitionSplit(
                     int *current_num,
                     CArbMshEdgeList::ArbNewQuadData *qdata,
                     CArbMshEdgeList::ArbQdEdge *qd_edge,
                     CArbMshTopo2D *msh_topo,
                     CArbMshTopo2D *quad_topo,
                     CArbMshEdgeList *edge_list,
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

        if (!DoTransitionSplit(qd_edge,qdata,msh_topo,
                          quad_topo,edge_list))
          return 0;
        SmoothFront(qdata->qcase,
                    3,&qdata->front_nodes[1],
                    msh_topo,quad_topo,edge_list) ;
        if (!FindQuadData(qd_edge,qdata,edge_list,msh_topo))
         return -1;
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
            CArbMshEdgeList::ArbNewQuadData tqdata ;
            int elem = !quad_topo->HasVtx(qdata->bl) ? NO_ELEM :
                       quad_topo->GetCCWElem(qdata->bl,qdata->tl) ;
            if (elem != NO_ELEM) {
                int num, nodes[4] ;
                quad_topo->GetElemNodes(elem,qdata->bl,&num,nodes) ;
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
                if (!DoTemplateSplit(&tqdata,msh_topo,quad_topo,
                                edge_list,true))
                   return -1;
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
            CArbMshEdgeList::ArbNewQuadData tqdata ;
            int elem = !quad_topo->HasVtx(qdata->tr) ? NO_ELEM :
                       quad_topo->GetCCWElem(qdata->tr,qdata->br) ;
            if (elem != NO_ELEM) {
                int num, nodes[4] ;
                quad_topo->GetElemNodes(elem,qdata->br,&num,nodes) ;
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
                if (!DoTemplateSplit(&tqdata,msh_topo,
                                quad_topo,edge_list,true))
                  return -1;
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
            CArbMshEdgeList::ArbNewQuadData tqdata ;
            int elem = !quad_topo->HasVtx(qdata->br) ? NO_ELEM :
                       quad_topo->GetCCWElem(qdata->br,qdata->bl) ;
            if (elem == NO_ELEM) return(0) ;
            int num_node, nodes[4] ;
            quad_topo->GetElemNodes(elem,qdata->br,&num_node,nodes) ;
            tqdata.bl = nodes[2] ;
            tqdata.br = nodes[3] ;
            tqdata.tl = nodes[1] ;
            tqdata.tr = nodes[0] ;
            tqdata.num_front_nodes = 4 ;
            tqdata.front_nodes[0] = qdata->tl ;
            tqdata.front_nodes[1] = tqdata.tl ;
            tqdata.front_nodes[2] = tqdata.tr ;
            tqdata.front_nodes[3] = qdata->tr ;
            if (!DoTemplateSplit(&tqdata,msh_topo,quad_topo,edge_list,true))
             return -1;
            *transition_flg = true ;
            ++(*current_num) ;
            return(1) ;
        }
    }

    // Here we check for a special case for doing transitions
    // starting from a boundary where we use a different tolerance

    if (!(*transition_flg) && (qdata->qcase == 3) &&
                (qdata->edge_level == 0)) {
        ArbIntNode *tl = pnode_table->Fetch(qdata->tl) ;
        ArbIntNode *bl = pnode_table->Fetch(qdata->bl) ;
        ArbIntNode *br = pnode_table->Fetch(qdata->br) ;
        double blen = (bl->coord - br->coord).Magnitude() ;
        double llen = (bl->coord - tl->coord).Magnitude() ;
        if (tl->motion == ARB_FIXED) {
            double ratio = llen / blen ;
            if (ratio > TRANSITION_NEAR_BOUNDARY_TOLERANCE) {
                CArbMshEdgeList::ArbNewQuadData tqdata ;
                int elem = !quad_topo->HasVtx(qdata->bl) ? NO_ELEM :
                           quad_topo->GetCCWElem(qdata->bl,qdata->tl) ;
                if (elem != NO_ELEM) {
                    int num, nodes[4] ;
                    quad_topo->GetElemNodes(elem,qdata->bl,&num,nodes) ;
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
                    if (!DoTemplateSplit(&tqdata,msh_topo,quad_topo,
                                    edge_list,true))
                       return -1;
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
        ArbIntNode *tr = pnode_table->Fetch(qdata->tr) ;
        double rlen = (bl->coord - tr->coord).Magnitude() ;
        if (tr->motion == ARB_FIXED) {
            double ratio = rlen / blen ;
            if (ratio > TRANSITION_NEAR_BOUNDARY_TOLERANCE) {
                CArbMshEdgeList::ArbNewQuadData tqdata ;
                int elem = !quad_topo->HasVtx(qdata->tr) ? NO_ELEM :
                       quad_topo->GetCCWElem(qdata->tr,qdata->br) ;
                if (elem != NO_ELEM) {
                    int num, nodes[4] ;
                    quad_topo->GetElemNodes(elem,qdata->br,&num,nodes) ;
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
                    if (!DoTemplateSplit(&tqdata,msh_topo,
                                    quad_topo,edge_list,true))
                       return -1;
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

    return(0) ;
}




// %(CArbMshRegion2D::UpdateForQuadMesh-void-|-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateForQuadMesh - update the region with quad mesh information
**
**      void UpdateForQuadMesh(CArbMshTopo2D *quad_topo)
**
**        quad_topo - (in)  quadrilateral mesh data
**
**      Description: This method updates the instance to replace the
**          triangular mesh information with quad mesh information.
**
**
** -- */

void CArbMshRegion2D::UpdateForQuadMesh(CArbMshTopo2D *quad_topo)
{
    // If these are quadratic order elements then we need to
    // build a data structure to keep track of the side nodes.

    CArbHashTable<ArbEdgeKey,int> side_nodes ;
    if (Order == QUADRATIC) {
        CArbHashTableIterator<int,ArbMshElement2D> eiter(pelem_table) ;
        int nn, j, nd0, nd1, nd2 ;
        for (eiter.First() ; eiter.More() ; ++eiter) {
            ArbMshElement2D *elem = eiter.Entry() ;
            nn = elem->num_nodes/2 ;
            for (j=0 ; j<nn ; ++j) {
                nd0 = elem->nodes[j] ;
                nd1 = elem->nodes[(j+1) % nn] ;
                nd2 = elem->nodes[j + nn] ;
                if (nd0 > nd1)
                    side_nodes.Store(ArbEdgeKey(nd0,nd1),nd2) ;
                else
                    side_nodes.Store(ArbEdgeKey(nd1,nd0),nd2) ;
            }
        }
    }

    // now delete the triangular mesh and build a quad mesh

    delete pelem_table ;
    pelem_table = new CArbHashTable<int,ArbMshElement2D> ;

    int *elems = quad_topo->GetElemList() ;

    for (int jj=0 ; jj<quad_topo->NumElements() ; ++jj) {
        ArbMshElement2D elem ;
        int num_nodes, nodes[4] ;

        quad_topo->GetElemNodes(elems[jj],&num_nodes,nodes) ;
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
                    side = side_nodes.Fetch(ArbEdgeKey(nodes[ii],nodes[kk])) ;
                else
                    side = side_nodes.Fetch(ArbEdgeKey(nodes[kk],nodes[ii])) ;

                if (side != 0) {
                    elem.nodes[ii+num_nodes] = *side ;
                } else {

                    // new node, don't worry about the coords
                    // because they will be udated during smoothing

                    elem.nodes[ii+num_nodes] =
                        NewNode(0.0,0.0,INTERIOR,ARB_FLOATING,false) ;
                    if (nodes[ii] > nodes[kk])
                        side_nodes.Store(
                                   ArbEdgeKey(nodes[ii],nodes[kk]),
                                   elem.nodes[ii+num_nodes]) ;
                    else
                        side_nodes.Store(
                                   ArbEdgeKey(nodes[kk],nodes[ii]),
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




// %(CArbMshRegion2D::DoPrint-void-|-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    DoPrint - debug routine to display a mesh
**
**      void DoPrint(CArbMshTopo2D *msh_topo)
**
**        msh_topo - (in)  mesh topology
**
**      Description: This method provides debugging support. It prints
**          information that can be used to display a mesh.
**
**
** -- */

void CArbMshRegion2D::DoPrint(CArbMshTopo2D *msh_topo)
{
    CArbQueue<CArbMshTopo2D::ArbEdge> *elist ;
    elist = msh_topo->EdgeList() ;
    CArbMshTopo2D::ArbEdge dedge ;
    while (elist->RemoveFromFront(&dedge)) {
        ArbIntNode *p_nd0, *p_nd1 ;
        p_nd0 = pnode_table->Fetch(dedge.nd0) ;
        p_nd1 = pnode_table->Fetch(dedge.nd1) ;
//        printf("# %d %d\n",dedge.nd0,dedge.nd1) ;
//        printf("l %g %g %g %g\n",
        printf("l %20.13g %20.13g %20.13g %20.13g\n",
               p_nd0->coord[0],p_nd0->coord[1],
               p_nd1->coord[0],p_nd1->coord[1]) ;
        if (do_labels) {
//            printf("t %g %g %d\n",
            printf("t %20.13g %20.13g %d\n",
                   p_nd1->coord[0],p_nd1->coord[1],p_nd1->id) ;
//            printf("t %g %g %d\n",
            printf("t %20.13g %20.13g %d\n",
                   p_nd0->coord[0],p_nd0->coord[1],p_nd0->id) ;
        }
    }
    fflush(stdout) ;
    delete elist ;
}

void CArbMshRegion2D::DisplayTopo(const char *label,CArbMshTopo2D *msh_topo)
{
    CArbQueue<CArbMshTopo2D::ArbEdge> *elist ;
    elist = msh_topo->EdgeList() ;
    CArbMshTopo2D::ArbEdge dedge ;
    while (elist->RemoveFromFront(&dedge)) {
        ArbIntNode *p_nd0, *p_nd1 ;
        p_nd0 = pnode_table->Fetch(dedge.nd0) ;
        p_nd1 = pnode_table->Fetch(dedge.nd1) ;
//        printf("# %d %d\n",dedge.nd0,dedge.nd1) ;
        printf("l %g %g %g %g\n",
               p_nd0->coord[0],p_nd0->coord[1],
               p_nd1->coord[0],p_nd1->coord[1]) ;
        if (do_labels) {
            printf("t %g %g %d\n",
                   p_nd1->coord[0],p_nd1->coord[1],p_nd1->id) ;
            printf("t %g %g %d\n",
                   p_nd0->coord[0],p_nd0->coord[1],p_nd0->id) ;
        }
    }
    printf("a %s\nn\n",label) ;
    fflush(stdout) ;
    delete elist ;
}



// %(CArbMshRegion2D::InitializeQuadGen-void-|-CArbMshTopo2D-|*-CArbMshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    InitializeQuadGen - initialization for quad generation
**
**      void InitializeQuadGen(
**              CArbMshTopo2D   *msh_topo,
**              CArbMshEdgeList *edge_list)
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

int CArbMshRegion2D::InitializeQuadGen (CArbMshTopo2D *msh_topo, CArbMshEdgeList *edge_list)
{
  // go through all the elements and build the adjacent
  // vertex table
  int i;

  CArbHashTableIterator<int,ArbMshElement2D> eiter(pelem_table) ;
  for (eiter.First() ; eiter.More() ; ++eiter)
  {
    ArbMshElement2D *elem = eiter.Entry() ;
    int nn ;
    if (elem->num_nodes == 8)
      nn = 4 ;
    else if (elem->num_nodes == 6)
      nn = 3 ;
    else
      nn = elem->num_nodes ;
    msh_topo->InsertElement (elem->elem_id, nn, elem->nodes) ;
  }

  // now loop through the boundary edges and
  // add them to the edge heap.

  ArbIntEdge **b_edges = pedge_table->GetKeyList() ;
  int num_bedges = pedge_table->NumEntries() ;

  for (i=0 ; i<num_bedges ; ++i)
  {
    int cw_id = msh_topo->GetCWBdryNode (b_edges[i]->node_id[0], b_edges[i]->node_id[1]);
    int ccw_id = msh_topo->GetCCWBdryNode (b_edges[i]->node_id[0], b_edges[i]->node_id[1]);

    //assert((cw_id != 0) || (ccw_id != 0));
    //if ((cw_id == 0) || (ccw_id == 0))
    //return 0;

    edge_list->InsertEdge (b_edges[i]->node_id[0], b_edges[i]->node_id[1], cw_id,ccw_id) ;
  }

  if (!edge_list->InitializeCodes())
   return 0;
  delete [] b_edges ;
  return 1;
}




// %(CArbMshRegion2D::FindElementSide-int-|-CArbMshTopo2D-|*-ArbIntNode-|*-ArbIntNode-|*-ArbIntNode-|*-bool-|)
/* ++ ----------------------------------------------------------
**
**    FindElementSide - find a candidate element side
**
**      int FindElementSide(
**              CArbMshTopo2D *msh_topo,
**              ArbIntNode    *this_node,
**              ArbIntNode    *prev_node,
**              ArbIntNode    *next_node,
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

int CArbMshRegion2D::FindElementSide(CArbMshTopo2D *msh_topo,
                                              ArbIntNode *this_node,
                                              ArbIntNode *prev_node,
                                              ArbIntNode *next_node,
                                              double base_length,
                                              bool prev_edge)
{
    CArbCoord2D normal ;
    double tol = 0.01 ;

    // Find the "ideal" normal, which is the direction of the
    // bisector of the included angle, then turn this into
    // a coordinate point relative to the base node.
    //printf("ES0\n");
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
    ArbIntNode *swap_trial ;

    //printf("ES01\n");
    use_vtx = CheckForEdge(msh_topo,this_node,prev_node,
                           next_node,prev_edge,
                           normal,&vtx_before,&vtx_after) ;

    //printf("ES02\n");

    if (use_vtx != 0) return(use_vtx) ;

    // if we did not find a suitable edge then we look to
    // see if we can do edge swapping.  To do this we
    // start with the "after" vertex and look through it's
    // adjacent verticies for the "before" vertex.  Once
    // we find this, the next adjacent vertex is the
    // candidate.

    CArbTopoAdjVtxIterator iter0(msh_topo,this_node->id) ;
    int elem_1, elem_2 ;

    //printf("ES1\n");
    while (iter0.AdjVtx() != vtx_before)
    {
     if (iter0.AdjVtx() == -1)
      return NO_NODE;
     ++iter0 ;
    }
    elem_1 = iter0.CcwElem() ;

    //printf("ES2\n");

    CArbTopoAdjVtxCyclicIterator iter1(msh_topo,vtx_after) ;

    // rgd
    if ((vtx_before == 0) && (vtx_after == 0))
     return NO_NODE;

    //rgd08_1
    int cont = 0;
    //printf("ES3\n");
    while (1)
    {
        if (cont == nelem*3)
         return NO_NODE;
        if (iter1.AdjVtx() == -1)
         return NO_NODE;
        if (iter1.AdjVtx() == vtx_before)
         break ;
        ++iter1 ;
        ++cont;
    }
    //printf("ES4\n");
    elem_2 = iter1.CcwElem() ;
    ++iter1 ;

    // if elem_2 is NO_ELEM, then we have the case that the
    // edge joining the before and after verticies is part
    // of the current boundary.  The best we can do is to
    // return either the before or after vertex depending
    // on which makes a smaller angle with the normal

    if (elem_2 == NO_ELEM) {
        ArbIntNode *before = pnode_table->Fetch(vtx_before) ;
        ArbIntNode *after = pnode_table->Fetch(vtx_after) ;
        if (Angle(this_node->coord,before->coord,normal) <
            Angle(this_node->coord,normal,after->coord))
            return(vtx_before) ;
        else
            return(vtx_after) ;
    }

    // now we have a candidate, check the angle

    swap_trial = pnode_table->Fetch(iter1.AdjVtx()) ;
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

    ArbIntNode *before = pnode_table->Fetch(vtx_before) ;
    ArbIntNode *after = pnode_table->Fetch(vtx_after) ;
    CArbCoord2D intsc = IntersectLines(this_node->coord,normal,
                                       before->coord,after->coord) ;

    // reject if the resulting edge will be too small

    double dist = (intsc-this_node->coord).Magnitude() ;
    if (dist < tol*base_length) return(-1) ;

    // define the new node

    int new_id = NewNode(intsc.x(),intsc.y(),
                         INTERIOR,ARB_FLOATING,true) ;

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




// %(CArbMshRegion2D::CheckForEdge-int-|^const-CArbMshTopo2D-|*-ArbIntNode-|*-ArbIntNode-|*-ArbIntNode-|*-bool-|-CArbCoord2D-|-int-|*-int-|*)
/* ++ ----------------------------------------------------------
**
**    CheckForEdge - find an edge in the search direction
**
**      int CheckForEdge(
**              CArbMshTopo2D *msh_topo,
**              ArbIntNode    *this_node,
**              ArbIntNode    *prev_node,
**              ArbIntNode    *next_node,
**              bool          prev_edge,
**              CArbCoord2D   normal,
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

int CArbMshRegion2D::CheckForEdge(
                        CArbMshTopo2D *msh_topo,
                        ArbIntNode *this_node,
                        ArbIntNode *prev_node,
                        ArbIntNode *next_node,
                        bool prev_edge,
                        CArbCoord2D normal,
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

    //printf("g0\n");

    int use_vtx = 0 ;
    int vtx_before_next = -1 ;
    int num_free = 0 ;

    *vtx_before = *vtx_after = 0 ;

    CArbTopoAdjVtxIterator iter(msh_topo,this_node->id) ;
    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.CcwElem() == NO_ELEM) ++num_free ;
    }

    //printf("g1\n");
    //assert(num_free > 0) ;
    if (num_free == 0)
     return -1;

    if (num_free == 1) {
        double angle_min = 6.2831 ;
        double angle_before = 6.2831 ;
        double angle_after = 6.2831 ;
        *vtx_after = 0 ;
        for (iter.First() ; iter.More() ; ++iter)
        {
            ArbIntNode *trial = pnode_table->Fetch(iter.AdjVtx()) ;
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
                        CArbTopoAdjVtxIterator biter = iter ;
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
       //printf("g1\n");
    } else {
       //printf("g2\n");
        CArbTopoAdjVtxCyclicIterator citer(msh_topo,this_node->id) ;
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

            //printf("g3\n");
            while (citer.AdjVtx() != prev_node->id)
            {
                if (citer.CcwElem() == NO_ELEM) {
                    in_void = true ;
                } else if (in_void) {
                    start_vtx = citer.AdjVtx() ;
                    in_void = false ;
                }

                ++citer ;
            }
            //printf("g4\n");
            // now search through the region and see if
            // if the normal direction is inside

            //rgd08_1
            int cont = 0;
            while (citer.AdjVtx() != start_vtx)
            {
             cont++;
             if (cont == nelem*3)
              return -1;
             ++citer ;
            }

            double angle_min = 6.2831 ;
            double angle_before = 6.2831 ;
            double angle_after = 6.2831 ;
            cont=0;

            //printf("g5\n");

            while (citer.AdjVtx() != stop_vtx)
            {
                if (citer.AdjVtx() == -1)
                 return -1;
                if (cont == nelem*3)
                 return -1;
                ArbIntNode *trial = pnode_table->Fetch(citer.AdjVtx()) ;
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
                        CArbTopoAdjVtxCyclicIterator biter = citer ;
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
                ++cont;
            }
            //printf("g6\n");
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

            //rgd08_1
            int cont=0;
            while (citer.CcwElem() != NO_ELEM)
            {
             cont++;
             if (cont == nelem*3)
              return -1;
             ++citer ;
            }
            stop_vtx = citer.AdjVtx() ;

            // now search through the region and see if
            // if the normal direction is inside

            // rgd08_1
            cont = 0;
            while (citer.AdjVtx() != start_vtx)
            {
             cont++;
             if (cont == nelem*3)
              return -1;
             ++citer ;
            }
            double angle_min = 6.2831 ;
            double angle_before = 6.2831 ;
            double angle_after = 6.2831 ;
            //rgd08_1
            cont=0;
            while (citer.AdjVtx() != stop_vtx) {
                if (cont == nelem*3)
                 return -1;
                ArbIntNode *trial = pnode_table->Fetch(citer.AdjVtx()) ;
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
                        CArbTopoAdjVtxCyclicIterator biter = citer ;
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
                ++cont;
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



bool CArbMshRegion2D::TweakNode(ArbIntNode *nd,CArbMshTopo2D *msh_topo,
                      CArbArray<int> &moved,
                      CArbArray<CArbCoord2D> &orig)
{
    // if this node is not movable we cannot do anything

    if (nd->motion == ARB_FIXED) return(false) ;

    // find the opposite edges on adjacent triangles and find
    // the closest

    bool first = true ;
    double min_dist = 0.0 ;
    CArbCoord2D move_dir ;

    CArbTopoAdjVtxIterator iter(msh_topo,nd->id) ;
    for (iter.First() ; iter.More() ; ++iter) {

        int elem_id = iter.CcwElem() ;
        CArbMshTopo2D::ArbEdge edge ;

        edge = msh_topo->OppositeEdge(elem_id,nd->id) ;
        ArbIntNode *nd0 = pnode_table->Fetch(edge.nd0) ;
        ArbIntNode *nd1 = pnode_table->Fetch(edge.nd1) ;

        if (nd0 == 0 || nd1 == 0)
         return false;
        double blen ;
        double dist = DistSqr(nd->coord,
                              nd0->coord,nd1->coord,&blen) ;

        if (first || (dist<min_dist)) {
            first = false ;
            min_dist = dist ;
            CArbCoord2D delta = nd1->coord - nd0->coord ;
            move_dir = CArbCoord2D(-delta.y(),delta.x()) ;
        }
    }

    // move the node one half the distance to the closest
    // edge in the opposite direction

    moved.InsertAtEnd(nd->id) ;
    orig.InsertAtEnd(nd->coord) ;

    nd->coord += (0.5*sqrt(min_dist)) * move_dir.Normalize() ;
    return(true) ;
}


// %(CArbMshRegion2D::FindCrossedEdges-CArbQueue-|<CArbMshTopo2D::ArbEdge>*^const-int-|-int-|-CArbMshTopo2D-|*-CArbMshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    FindCrossedEdges - find crossed edges
**
**      CArbQueue <CArbMshTopo2D::ArbEdge>*FindCrossedEdges(
**              int             start_id,
**              int             stop_id,
**              CArbMshTopo2D   *msh_topo,
**              CArbMshEdgeList *edge_list) const
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

CArbQueue<CArbMshTopo2D::ArbEdge>
*CArbMshRegion2D::FindCrossedEdges(int start_id,
                      int stop_id,
                      CArbMshTopo2D *msh_topo,
                      CArbMshEdgeList *bdry_list,
                      CArbArray<int> &moved,
                      CArbArray<CArbCoord2D> &orig)
{
    CArbQueue<CArbMshTopo2D::ArbEdge> *edge_list =
        new CArbQueue<CArbMshTopo2D::ArbEdge>() ;

    moved.Clear() ;
    orig.Clear() ;

    ArbIntNode *nd_start = pnode_table->Fetch(start_id) ;
    ArbIntNode *nd_stop  = pnode_table->Fetch(stop_id) ;

    // look around the start node and find the two edges that
    // are before and after the search vector (saved in adj_0
    // and adj_1).

    ArbIntNode *adj_1 ;
    double cross_0 = .0f, cross_1 ;
    int elem_id = 0, elem_id_1 = 0;
    bool found = false ;
    bool first_trip = true ;
    double tol = 1.0e-4 * (nd_start->coord-nd_stop->coord).Magnitude() ;

    CArbTopoAdjVtxIterator iter(msh_topo,start_id) ;
    if (!iter.More() || (iter.AdjVtx() == stop_id)) return(edge_list) ;

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.AdjVtx() == stop_id) return(edge_list) ;
        adj_1 = pnode_table->Fetch(iter.AdjVtx()) ;
        if (adj_1 == 0)
         return 0;
        elem_id_1 = iter.CcwElem() ;
        cross_1 = CrossProd(nd_start->coord,adj_1->coord,nd_stop->coord) ;

        if ((fabs(cross_1) < tol) &&
            (fabs(Angle(nd_start->coord,adj_1->coord,nd_stop->coord)) < 1.57) &&
            (adj_1->motion != ARB_FIXED)) {

            // if we get here then there is an edge that is along
            // the search direction.  The strategy is to move this
            // node off the search line.

            // first loop around the edge and find edge adjacent
            // to the start node

            CArbTopoAdjVtxCyclicIterator citer(msh_topo,adj_1->id) ;
            while (1) {
                if (citer.AdjVtx() == start_id) break ;
                ++citer ;
            }
            ++citer ;
            ArbIntNode *adj_t = pnode_table->Fetch(citer.AdjVtx()) ;
            if (adj_t == 0)
             return 0;
            CArbCoord2D mid = 0.5 * (adj_1->coord + adj_t->coord) ;

            // check to see if the new node position will cross
            // the current advancing front

            if (bdry_list->DoesThisCrossBoundary(nd_start->coord,mid)) {
                delete edge_list ;
                return(0) ;
            }

            moved.InsertAtEnd(iter.AdjVtx()) ;
            orig.InsertAtEnd(adj_1->coord) ;

            adj_1->coord = mid ;
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

    CArbMshTopo2D::ArbEdge edge, next_edge ;
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

    int cont=0;
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

        CArbMshTopo2D::ArbEdge rev_edge ;
        rev_edge.nd0 = edge.nd1 ;
        rev_edge.nd1 = edge.nd0 ;
        rev_edge.elem0 = edge.elem1 ;
        rev_edge.elem1 = edge.elem0 ;
        int opp_id = msh_topo->OppositeNode(elem_id,&rev_edge) ;
        ArbIntNode *nd_opp = pnode_table->Fetch(opp_id) ;
        if (nd_opp == 0)
         return 0;

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
        cont++;
        //printf("c%i\n",cont);
        // rgd
        if (cont == nelem*3)
         return(0);
    }

    return(edge_list) ;
}




// %(CArbMshRegion2D::RecoverTopEdge-bool-|^const-int-const|-int-const|-CArbMshTopo2D-|*-CArbMshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    RecoverTopEdge - do edge swapping
**
**      bool RecoverTopEdge(
**              const int       start_id,
**              const int       stop_id,
**              CArbMshTopo2D   *msh_topo,
**              CArbMshEdgeList *edge_list) const
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

bool CArbMshRegion2D::RecoverTopEdge(
             int start_id,
             int stop_id,
             CArbMshTopo2D *msh_topo,
             CArbMshEdgeList *edge_list)
{
    // first find the coordinates of the start and stop nodes

    ArbIntNode *nd_start = pnode_table->Fetch(start_id) ;
    ArbIntNode *nd_stop  = pnode_table->Fetch(stop_id) ;

    if (nd_start == 0)
     return false;
    if (nd_stop == 0)
     return false;
    // now find get a list of edges crossing the vector between
    // the start and stop nodes

    CArbQueue<CArbMshTopo2D::ArbEdge> *crossed ;
    CArbArray<int> moved ;
    CArbArray<CArbCoord2D> orig ;

    crossed = FindCrossedEdges(start_id,stop_id,
                               msh_topo,edge_list,moved,orig) ;
    if (crossed == 0) return(false) ;
    if (crossed->NumEntries() > MAX_CROSSED_EDGES) {
       delete crossed ;
       return(false) ;
    }

    int max_iters = crossed->NumEntries() <= 10 ? 100 :
                    crossed->NumEntries() * crossed->NumEntries() ;

    CArbMshTopo2D::ArbEdge edge ;
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

    CArbArray<int> swap_ids ;

    // loop while there are edges

    while (crossed->RemoveFromFront(&edge)) {

        iter++ ;
        if ((saved && (edge.nd0 == saved_nd[0]) &&
                      (edge.nd1 == saved_nd[1])) ||
            (iter > max_iters)) {

            // unwind the changes that have been made to
            // the mesh
            int cur = swap_ids.NumEntries() - 6 ;
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

            for (int i=moved.NumEntries()-1 ; i >= 0 ; --i) {
                ArbIntNode *nd = pnode_table->Fetch(moved[i]) ;
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

        ArbIntNode *nd_0 = pnode_table->Fetch(id_0) ;
        ArbIntNode *nd_1 = pnode_table->Fetch(id_1) ;

        edge.elem0 = msh_topo->GetCCWElem(id_0,id_1) ;
        edge.elem1 = msh_topo->GetCCWElem(id_1,id_0) ;

        int id_2 = msh_topo->OppositeNode(edge.elem0,&edge) ;
        if (id_2 == -1)
         return false;
        ArbIntNode *nd_2 = pnode_table->Fetch(id_2) ;

        if (nd_2 == 0)
         return false;

        CArbMshTopo2D::ArbEdge rev_edge ;
        rev_edge.nd0 = edge.nd1 ;
        rev_edge.nd1 = edge.nd0 ;
        rev_edge.elem0 = edge.elem1 ;
        rev_edge.elem1 = edge.elem0 ;
        int id_3 = msh_topo->OppositeNode(edge.elem1,&rev_edge) ;
        ArbIntNode *nd_3 = pnode_table->Fetch(id_3) ;

        if (nd_3 == 0)
         return false;

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

            swap_ids.InsertAtEnd(id_0) ;
            swap_ids.InsertAtEnd(id_1) ;
            swap_ids.InsertAtEnd(id_2) ;
            swap_ids.InsertAtEnd(id_3) ;
            swap_ids.InsertAtEnd(eid_0) ;
            swap_ids.InsertAtEnd(eid_1) ;

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
    CArbTopoAdjVtxIterator viter(msh_topo,start_id) ;
    for (viter.First() ; viter.More() ; ++viter) {
        if (viter.AdjVtx() == stop_id) {
            found = true ;
            break ;
        }
    }

    delete crossed ;
    return(found) ;
}




// %(CArbMshRegion2D::ClearRegion-void-|-int-const|-int-const|*-CArbMshTopo2D-|*-CArbMshTopo2D::ArbEdge-const|*)
/* ++ ----------------------------------------------------------
**
**    ClearRegion - clear a quadrilateral region
**
**      void ClearRegion(
**              const int                    num_corners,
**              const int                    *boundary,
**              CArbMshTopo2D                *msh_topo,
**              const CArbMshTopo2D::ArbEdge *edge)
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

int CArbMshRegion2D::ClearRegion(const int num_corners,
                                  const int *boundary,
                                  CArbMshTopo2D *msh_topo,
                                  const CArbMshTopo2D::ArbEdge *edge)
{
    int i ;
    bool bound ;

    // store the next cw and ccw edges

    CArbMshTopo2D::ArbEdge cw_edge =
         msh_topo->NextCWEdge(edge->elem0,edge) ;
    if (cw_edge.nd0 == -1 && cw_edge.nd1 == -1 && cw_edge.elem0 == -1 && cw_edge.elem1 == -1)
     return 0;
    CArbMshTopo2D::ArbEdge ccw_edge =
         msh_topo->NextCCWEdge(edge->elem0,edge) ;
    if (ccw_edge.nd0 == -1 && ccw_edge.nd1 == -1 && ccw_edge.elem0 == -1 && ccw_edge.elem1 == -1)
     return 0;

    // delete this element from the topology

    if (!msh_topo->DeleteTriangle(edge->nd0,edge->nd1,ccw_edge.nd1))
     return 0;

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
            CArbMshTopo2D::ArbEdge rev_edge ;
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
            CArbMshTopo2D::ArbEdge rev_edge ;
            rev_edge.nd0 = cw_edge.nd1 ;
            rev_edge.nd1 = cw_edge.nd0 ;
            rev_edge.elem0 = cw_edge.elem1 ;
            rev_edge.elem1 = cw_edge.elem0 ;
            ClearRegion(num_corners,boundary,msh_topo,&rev_edge) ;
        }
    }
    return 1;
}




// %(CArbMshRegion2D::SmoothFront-void-|-int-const|-int-const|-int-const|*-CArbMshTopo2D-|*-CArbMshTopo2D-|*-CArbMshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    SmoothFront - smooth advancing front nodes
**
**      void SmoothFront(
**              const int       qcase,
**              const int       num_frt_nodes,
**              const int       *frt_nodes,
**              CArbMshTopo2D   *msh_topo,
**              CArbMshTopo2D   *quad_topo,
**              CArbMshEdgeList *edge_list)
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

CArbCoord2D CArbMshRegion2D::LaplaceSmoothFrnt(const int vtx,
                                        CArbMshTopo2D *msh_topo,
                                        CArbMshTopo2D *quad_topo)
{
    int num_adj = quad_topo->NumAdjElems(vtx) ;
    ArbIntNode *node = pnode_table->Fetch(vtx) ;
    ArbIntNode *adj ;
    CArbCoord2D sum ;
    int num = 0 ;

    // from the mesh topology get the coordinates
    // of the adjacent nodes

    CArbTopoAdjVtxIterator iter0(msh_topo,vtx) ;
    int last = 0 ;

    if (num_adj < 3) {
        for (iter0.First() ; iter0.More() ; ++iter0)
        {
             adj = pnode_table->Fetch(iter0.AdjVtx()) ;
             if (adj == 0)
             {
              CArbCoord2D sum1(-1,-1);
              return sum1;
             }
             sum += adj->coord - node->coord ;
             last = iter0.AdjVtx() ;
             ++num ;
        }
    }

    // now loop through the quad topology, ignoring
    // the first and last nodes, which should have
    // been dealt with in the triangle mesh

    CArbTopoAdjVtxIterator iter1(quad_topo,vtx) ;

    for (iter1.First() ; iter1.More() ; ++iter1)
    {
        if (iter1.AdjVtx() == last) break ;
        adj = pnode_table->Fetch(iter1.AdjVtx()) ;
        if (adj == 0)
        {
         CArbCoord2D sum1(-1,-1);
         return sum1;
        }
        sum += adj->coord - node->coord ;
        ++num ;
    }

    if (num == 0) num = 1 ;
    return(sum/num) ;
}




// %(CArbMshRegion2D::IsoParametricSmooth-CArbCoord2D-|-ArbIntNode-|*-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    IsoParametricSmooth - smooth one advancing front node
**
**      CArbCoord2D IsoParametricSmooth(
**              ArbIntNode    *node,
**              CArbMshTopo2D *quad_topo)
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

CArbCoord2D CArbMshRegion2D::IsoParametricSmooth(ArbIntNode *node,
                                          CArbMshTopo2D *quad_topo)
{
    CArbTopoAdjVtxIterator iter(quad_topo,node->id) ;
    int num_nodes, nodes[4] ;
    CArbCoord2D sum(0.0,0.0) ;
    int num_adj = 0 ;
    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.CcwElem() != NO_ELEM) {
            if (!quad_topo->GetElemNodes(iter.CcwElem(),node->id,
                                    &num_nodes,nodes))
              return 0;
            sum += pnode_table->Fetch(nodes[1])->coord ;
            sum -= pnode_table->Fetch(nodes[2])->coord ;
            sum += pnode_table->Fetch(nodes[3])->coord ;
            ++num_adj ;
        }
    }
    CArbCoord2D new_pos = sum / num_adj ;
    return(new_pos-node->coord) ;
}




// %(CArbMshRegion2D::LaplaceSmooth-CArbCoord2D-|-ArbIntNode-|*-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    LaplaceSmooth - smooth one internal node
**
**      CArbCoord2D LaplaceSmooth(
**              ArbIntNode    *node,
**              CArbMshTopo2D *topo)
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

CArbCoord2D CArbMshRegion2D::LaplaceSmooth(ArbIntNode *node,
                                           CArbMshTopo2D *topo)
{
    CArbCoord2D delt_sum(0.0,0.0) ;
    ArbIntNode *adj ;
    CArbTopoAdjVtxIterator iter(topo,node->id) ;
    int num_adj = 0 ;
    for (iter.First() ; iter.More() ; ++iter) {
        adj = pnode_table->Fetch(iter.AdjVtx()) ;
        delt_sum += adj->coord - node->coord ;
        ++num_adj ;
    }
    if (num_adj == 0) num_adj = 1 ;
    return(delt_sum / num_adj) ;
}




// %(CArbMshRegion2D::CheckValidCoord-bool-|^const-ArbIntNode-|*-CArbCoord2D-|-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    CheckValidCoord - check for valid triangles
**
**      bool CheckValidCoord(
**              ArbIntNode    *node,
**              CArbCoord2D   delta,
**              CArbMshTopo2D *topo) const
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

bool CArbMshRegion2D::CheckValidCoord(ArbIntNode *node,
                                      CArbCoord2D delta,
                                      CArbMshTopo2D *topo) const
{
    // here we chech to make sure that all new triangles
    // have valid metrics.

    CArbCoord2D ncoord, prev ;
    ncoord = node->coord + delta ;
    CArbTopoAdjVtxIterator iter(topo,node->id) ;
    ArbIntNode *adj ;

    if (!iter.More()) return(true) ;
    //++iter ;
    adj = pnode_table->Fetch(iter.AdjVtx()) ;
    prev = adj->coord ;

    for (++iter ; iter.More() ; ++iter) {
        adj = pnode_table->Fetch(iter.AdjVtx()) ;
        if (Area(ncoord,prev,adj->coord) <= 0.0) {
            return(false) ;
        }
        prev = adj->coord ;
    }
    return(true) ;
}




// %(CArbMshRegion2D::CheckValidTriangleList-bool-|^const-CArbCoord2D-const|-CArbCoord2D-const|-CArbArray-const|<SeamCoordCache>*)
/* ++ ----------------------------------------------------------
**
**    CheckValidTriangleList - check a list of triangles
**
**      bool CheckValidTriangleList(
**              const CArbCoord2D                coord,
**              const CArbCoord2D                delta,
**              const CArbArray<SeamCoordCache>* triangles) const
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

bool CArbMshRegion2D::CheckValidTriangleList(const CArbCoord2D coord,
                  const CArbCoord2D delta,
                  const CArbArray<SeamCoordCache> *triangles) const
{
    CArbCoord2D ncoord ;
    ncoord = coord + delta ;
    for (int i=0 ; i<triangles->NumEntries() ; ++i) {
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




// %(CArbMshRegion2D::CheckValidQuadList-bool-|^const-CArbCoord2D-const|-CArbCoord2D-const|-CArbArray-const|<SeamCoordCache>*)
/* ++ ----------------------------------------------------------
**
**    CheckValidQuadList - check a list of quadrilaterals
**
**      bool CheckValidQuadList(
**              const CArbCoord2D                coord,
**              const CArbCoord2D                delta,
**              const CArbArray<SeamCoordCache>* triangles) const
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

bool CArbMshRegion2D::CheckValidQuadList(const CArbCoord2D coord,
                  const CArbCoord2D delta,
                  const CArbArray<SeamCoordCache> *quads) const
{
    CArbCoord2D ncoord ;
    ncoord = coord + delta ;
    for (int j=0 ; j<quads->NumEntries() ; ++j) {
        if (quads->At(j).ignore) continue ;
        double metric = QuadMetric(ncoord,quads->At(j).coord[1],
                                          quads->At(j).coord[2],
                                          quads->At(j).coord[3]) ;
        if (quads->At(j).boundary) {
            if (metric <= MINIMUM_BOUNDARY_METRIC) return(false) ;
        } else {
            if (metric <= 0.0) return(false) ;
        }
    }
    return(true) ;
}




// %(CArbMshRegion2D::CheckSmoothCoords-void-|-ArbIntNode-|*-CArbCoord2D-|-CArbMshTopo2D-|*-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    CheckSmoothCoords - check a proposed smoothed node location
**
**      void CheckSmoothCoords(
**              ArbIntNode    *node,
**              CArbCoord2D   delta,
**              CArbMshTopo2D *msh_topo,
**              CArbMshTopo2D *quad_topo)
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

int CArbMshRegion2D::CheckSmoothCoords(ArbIntNode *node,
                                        CArbCoord2D delta,
                                        CArbMshTopo2D *msh_topo,
                                        CArbMshTopo2D *quad_topo)
{
    int i,j ;

    // here we chech to make sure that all new triangles and
    // quads have valid metrics.  To do this
    // we start with the full delta.  If we detect an
    // invalid element we then try 0.75, 0.5, and 0.25
    // times the delta.

    // first save the coordinates of the triangle

    CoordCache->Clear() ;
    CArbTopoAdjVtxIterator iter(msh_topo,node->id) ;
    for (iter.First() ; iter.More() ; ++iter)
    {
        if ((iter.AdjVtx() == -1) || (iter.AdjVtx() < 0))
         return 0;
        ArbIntNode *adj = pnode_table->Fetch(iter.AdjVtx()) ;
        CoordCache->InsertAtEnd(adj->coord) ;
    }

    // compute the current minimum shape metric and the check
    // tolerance.  We allow existing bad elements to get up to
    // 10% worse.

    double min_metric = 1000.0 ;
    for (i=0 ; i<CoordCache->NumEntries()-1 ; ++i) {
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
        CArbCoord2D ncoord ;
        valid = true ;
        ncoord = node->coord + tri_factor * delta ;
        for (i=0 ; i<CoordCache->NumEntries()-1 ; ++i) {
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
    if (!valid) return 1;

    // now look at the quads

    CoordCache->Clear() ;
    CArbTopoAdjVtxIterator iter1(quad_topo,node->id) ;
    for (iter1.First() ; iter1.More() ; ++iter1) {
        if (iter1.CcwElem() != NO_ELEM) {
            int num_nodes, nodes[4] ;
            quad_topo->GetElemNodes(iter1.CcwElem(),node->id,
                                    &num_nodes,nodes) ;
            for (int ii=1 ; ii<num_nodes ; ++ii) {
                ArbIntNode *adj = pnode_table->Fetch(nodes[ii]) ;
                CoordCache->InsertAtEnd(adj->coord) ;
            }
        }
    }

    // look for the minimum existing metric

    min_metric = 1000.0 ;
    for (j=0 ; j<CoordCache->NumEntries() ; j += 3) {
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
        CArbCoord2D ncoord = node->coord + quad_factor * delta ;
        for (j=0 ; j<CoordCache->NumEntries() ; j += 3) {
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

    if (!valid) return 1;

    // update the coordinates

    if (tri_factor < quad_factor) {
        node->coord = node->coord + tri_factor * delta ;
    } else {
        node->coord = node->coord + quad_factor * delta ;
    }
    return 1;
}




// %(CArbMshRegion2D::SmoothFront-void-|-int-const|-int-const|-int-const|*-CArbMshTopo2D-|*-CArbMshTopo2D-|*-CArbMshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    SmoothFront - smooth advancing front nodes
**
**      void SmoothFront(
**              const int       qcase,
**              const int       num_frt_nodes,
**              const int       *frt_nodes,
**              CArbMshTopo2D   *msh_topo,
**              CArbMshTopo2D   *quad_topo,
**              CArbMshEdgeList *edge_list)
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

int CArbMshRegion2D::SmoothFront(const int qcase,
                                  const int num_frt_nodes,
                                  const int *frt_nodes,
                                  CArbMshTopo2D *msh_topo,
                                  CArbMshTopo2D *quad_topo,
                                  CArbMshEdgeList *edge_list)
{
    int i ;

    // we smooth the new nodes first, then the adjacent nodes

    switch (qcase) {
        case 1:
            if (!SmoothOneFront(frt_nodes[2],msh_topo,quad_topo,edge_list))
             return 0;
            if (!SmoothOneFront(frt_nodes[3],msh_topo,quad_topo,edge_list))
             return 0;
            for (i=0 ; i<num_frt_nodes ; ++i)
            {
                if ((i < 2) || (i > 3))
                    if (!SmoothOneFront(frt_nodes[i],msh_topo,
                                   quad_topo,edge_list))
                      return 0;
            }
            break ;

        case 2:
            if (!SmoothOneFront(frt_nodes[2],msh_topo,quad_topo,edge_list))
             return 0;
            for (i=0 ; i<num_frt_nodes ; ++i)
            {
                if (i != 2)
                 if (!SmoothOneFront(frt_nodes[i],msh_topo,
                                           quad_topo,edge_list))
                   return 0;
            }
            break ;

        case 3:
        case 5:
            for (i=0 ; i<num_frt_nodes ; ++i)
            {
                if (!SmoothOneFront(frt_nodes[i],msh_topo,
                               quad_topo,edge_list))
                  return 0;
            }
            break ;
    }
    return 1;
}




// %(CArbMshRegion2D::SmoothSeamFront-void-|-CArbMshEdgeList::ArbNewSeamData-|*-CArbMshTopo2D-|*-CArbMshTopo2D-|*-CArbMshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    SmoothSeamFront - smooth advancing front nodes
**
**      void SmoothSeamFront(
**              CArbMshEdgeList::ArbNewSeamData *sdata,
**              CArbMshTopo2D                   *msh_topo,
**              CArbMshTopo2D                   *quad_topo,
**              CArbMshEdgeList                 *edge_list)
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

void CArbMshRegion2D::SmoothSeamFront(
                                  CArbMshEdgeList::ArbNewSeamData *sdata,
                                  CArbMshTopo2D *msh_topo,
                                  CArbMshTopo2D *quad_topo,
                                  CArbMshEdgeList *edge_list)
{
    SmoothOneFront(sdata->id0,msh_topo,quad_topo,edge_list) ;
    SmoothOneFront(sdata->cw,msh_topo,quad_topo,edge_list) ;
    SmoothOneFront(sdata->ccw,msh_topo,quad_topo,edge_list) ;
}




// %(CArbMshRegion2D::SmoothOneFront-void-|-int-const|-CArbMshTopo2D-|*-CArbMshTopo2D-|*-CArbMshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    SmoothOneFront - smooth one advancing front node
**
**      void SmoothOneFront(
**              const int       node_id,
**              CArbMshTopo2D   *msh_topo,
**              CArbMshTopo2D   *quad_topo,
**              CArbMshEdgeList *edge_list)
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

int CArbMshRegion2D::SmoothOneFront(const int node_id,
                                     CArbMshTopo2D *msh_topo,
                                     CArbMshTopo2D *quad_topo,
                                     CArbMshEdgeList *edge_list)
{
    ArbIntNode *node = pnode_table->Fetch(node_id) ;
    if (node->motion == ARB_FLOATING) {

        int num_total = quad_topo->NumAdjElems(node_id) ;
        int num_adj =
            quad_topo->NumConsecutiveAdjElems(node_id) ;

        CArbCoord2D delta_a, delta_b, delta_c ;

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
            ArbIntNode *adj ;
            int num_nodes, nodes[4] ;
            int prev, next ;

            CArbTopoAdjVtxIterator iter(quad_topo,node->id) ;
            if (!quad_topo->GetElemNodes(iter.CcwElem(),node->id,
                                    &num_nodes,nodes))
              return 0;

            adj = pnode_table->Fetch(nodes[3]) ;
            if (adj == 0)
             return 0;
            CoordCache->InsertAtEnd(adj->coord) ;
            adj = pnode_table->Fetch(nodes[1]) ;
            if (adj == 0)
             return 0;
            CoordCache->InsertAtEnd(adj->coord) ;
            adj = pnode_table->Fetch(nodes[2]) ;
            if (adj == 0)
             return 0;
            CoordCache->InsertAtEnd(adj->coord) ;
            prev = nodes[1] ;

            ++iter ;
            if (!quad_topo->GetElemNodes(iter.CcwElem(),node->id,
                                    &num_nodes,nodes))
             return 0;

            adj = pnode_table->Fetch(nodes[2]) ;
            if (adj == 0)
             return 0;
            CoordCache->InsertAtEnd(adj->coord) ;
            adj = pnode_table->Fetch(nodes[3]) ;
            if (adj == 0)
             return 0;
            CoordCache->InsertAtEnd(adj->coord) ;
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

                CArbTopoAdjVtxIterator iter(msh_topo,node_id) ;
                int num = 0 ;

                for (iter.First() ; iter.More() ; ++iter) {
                    if ((iter.AdjVtx() != prev) &&
                        (iter.AdjVtx() != next)) {
                        adj = pnode_table->Fetch(iter.AdjVtx()) ;
                        sum += (adj->coord - node->coord).Magnitude() ;
                        ++num ;
                    }
                }
                ld = sum / (4 + num) ;
            } else {
                // Laplace Smooth

                double sum = 0 ;
                int num = 0 ;
                CArbTopoAdjVtxIterator iter1(quad_topo,node_id) ;

                for (iter1.First() ; iter1.More() ; ++iter1) {
                    adj = pnode_table->Fetch(iter1.AdjVtx()) ;
                    sum += (adj->coord - node->coord).Magnitude() ;
                    ++num ;
                }

                CArbTopoAdjVtxIterator iter(msh_topo,node_id) ;
                for (iter.First() ; iter.More() ; ++iter) {
                    if ((iter.AdjVtx() != prev) &&
                        (iter.AdjVtx() != next)) {
                        adj = pnode_table->Fetch(iter.AdjVtx()) ;
                        sum += (adj->coord - node->coord).Magnitude() ;
                        ++num ;
                    }
                }

                ld = sum / num ;
            }

            // find the length of the vector to the
            // isoparametric point

            delta_a = IsoParametricSmooth(node,quad_topo) ;
            CArbCoord2D v = node->coord - CoordCache->At(0) ;

            la = (v + delta_a).Magnitude() ;

            // find the proposed update

            double lratio = ld / la ;
            delta_b = -v + (delta_a + v) * lratio ;

            // now for the angles, first find bisector
            // of the angle between nodes 3 - 0 - 2
            // (with reference to the figure above)

            CArbCoord2D pb1, pb2, q, tmp ;
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
        } else
        {
            delta_a = LaplaceSmoothFrnt(node_id,msh_topo,quad_topo) ;
            if (delta_a.x() == -1 && delta_a.y() == -1)
             return 0;
        }
        if (!CheckSmoothCoords(node,delta_a,msh_topo,quad_topo))
         return 0;
    }
    return 1;
}




// %(CArbMshRegion2D::SmoothAdjacent-void-|-int-const|-int-const|*-CArbMshTopo2D-|*-bool-const|-bool-const|-CArbMshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    SmoothAdjacent - smooth adjacent nodes
**
**      void SmoothAdjacent(
**              const int       num_frt_nodes,
**              const int       *frt_nodes,
**              CArbMshTopo2D   *topo,
**              const bool      look_for_bdry,
**              const bool      triangle,
**              CArbMshEdgeList *edge_list)
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

void CArbMshRegion2D::SmoothAdjacent(
                                  const int num_frt_nodes,
                                  const int *frt_nodes,
                                  CArbMshTopo2D *topo,
                                  const bool look_for_bdry,
                                  const bool triangle,
                                  CArbMshEdgeList *edge_list)
{
    for (int i=0 ; i<num_frt_nodes-2 ; ++i) {
        SmoothOneAdj(frt_nodes[i+1],frt_nodes[i],
                     frt_nodes[i+2],topo,look_for_bdry,
                     triangle,edge_list) ;
    }
}




// %(CArbMshRegion2D::SmoothOneAdj-void-|-int-const|-int-const|-int-const|-CArbMshTopo2D-|*-bool-const|-bool-const|-CArbMshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    SmoothOneAdj - smooth one adjacent node
**
**      void SmoothOneAdj(
**              const int       vtx,
**              const int       prev,
**              const int       next,
**              CArbMshTopo2D   *topo,
**              const bool      look_for_bdry,
**              const bool      triangle,
**              CArbMshEdgeList *edge_list)
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

void CArbMshRegion2D::SmoothOneAdj(const int vtx,
                                   const int prev,
                                   const int next,
                                   CArbMshTopo2D *topo,
                                   const bool /*look_for_bdry*/,
                                   const bool triangles,
                                   CArbMshEdgeList *edge_list)
{
    CArbTopoAdjVtxIterator iter(topo,vtx) ;

    for (iter.First() ; iter.More() ; ++iter) {
        if ((iter.AdjVtx() != prev) &&
            (iter.AdjVtx() != next) &&
            !edge_list->ContainsNode(iter.AdjVtx())) {
            ArbIntNode *node = pnode_table->Fetch(iter.AdjVtx()) ;
            if (node->motion == ARB_FLOATING) {

                CArbCoord2D sum(0.0,0.0) ;
                int num = 0 ;
                CArbTopoAdjVtxIterator aiter(topo,iter.AdjVtx()) ;

                for (aiter.First() ; aiter.More() ; ++aiter) {
                    ArbIntNode *adj =
                        pnode_table->Fetch(aiter.AdjVtx()) ;
                    sum += adj->coord - node->coord ;
                    ++num ;
                }

                if (num == 0) num = 1 ;
                sum /= num ;

                if (triangles) {
                    CoordCache->Clear() ;
                    CArbTopoAdjVtxIterator aiter(topo,iter.AdjVtx()) ;
                    for (aiter.First() ; aiter.More() ; ++aiter) {
                        ArbIntNode *adj =
                            pnode_table->Fetch(aiter.AdjVtx()) ;
                        CoordCache->InsertAtEnd(adj->coord) ;
                    }
                    aiter.First() ;
                    ArbIntNode *adj =
                        pnode_table->Fetch(aiter.AdjVtx()) ;
                    CoordCache->InsertAtEnd(adj->coord) ;

                    // now check the areas

                    bool valid = true ;
                    double tri_factor ;
                    for (tri_factor = 1.0 ; tri_factor > 0.0 ; tri_factor -= 0.25) {
                        valid = true ;
                        CArbCoord2D ncoord = node->coord + tri_factor * sum ;
                        for (int i=0 ; i<CoordCache->NumEntries()-1 ; ++i) {
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
                    CArbTopoAdjVtxIterator aiter(topo,iter.AdjVtx()) ;
                    for (aiter.First() ; aiter.More() ; ++aiter) {
                        if (aiter.CcwElem() != NO_ELEM) {
                            int num_nodes, nodes[4] ;
                            topo->GetElemNodes(aiter.CcwElem(),
                                    iter.AdjVtx(),
                                    &num_nodes,nodes) ;
                            for (int ii=1 ; ii<num_nodes ; ++ii) {
                                ArbIntNode *adj = pnode_table->Fetch(nodes[ii]) ;
                                CoordCache->InsertAtEnd(adj->coord) ;
                            }
                        }
                    }

                    // now check the metrics

                    double quad_factor = 1.0 ;
                    bool valid = true ;
                    for (quad_factor = 1.0 ; quad_factor > 0.0 ; quad_factor -= 0.25) {
                        valid = true ;
                        CArbCoord2D ncoord = node->coord + quad_factor * sum ;
                        for (int j=0 ; j<CoordCache->NumEntries() ; j += 3) {
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




// %(CArbMshRegion2D::FindQuadData-void-|-CArbMshEdgeList::ArbQdEdge-const|*-CArbMshEdgeList::ArbNewQuadData-|*-CArbMshEdgeList-|*-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    FindQuadData - find a candidate quadrilateral
**
**      void FindQuadData(
**              const CArbMshEdgeList::ArbQdEdge *edge,
**              CArbMshEdgeList::ArbNewQuadData  *qdata,
**              CArbMshEdgeList                  *edge_list,
**              CArbMshTopo2D                    *msh_topo)
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

int CArbMshRegion2D::FindQuadData(const CArbMshEdgeList::ArbQdEdge *edge,
                                   CArbMshEdgeList::ArbNewQuadData *qdata,
                                   CArbMshEdgeList *edge_list,
                                   CArbMshTopo2D *msh_topo)
{
    //printf("F0\n");
    qdata->obl = edge->id_0 ;
    qdata->obr = edge->id_1 ;
    qdata->edge_level = edge->level ;

    // determine the nodes that will be at the ends of the
    // side edges.

    int end_code = edge->end_code & NO_SEAM_MASK ;

    //printf("F1\n");
    if ((end_code == 0) || (end_code == 2)) {

        //  We have either code 0 or code 2, so define the
        //  edge at the right node.
        //
        //             code 0            *    code 2
        //                               |
        //      *-----*=====*-----*      *=====*-----*

        int next_nd_id ;
        ArbIntNode *prev_node, *this_node, *next_node ;

        prev_node = pnode_table->Fetch(edge->id_0) ;
        if (prev_node == 0)
         return 0;
        this_node = pnode_table->Fetch(edge->id_1) ;
        if (this_node == 0)
         return 0;
        next_nd_id = edge_list->GetCCWNode(edge->id_0,edge->id_1) ;
        next_node = pnode_table->Fetch(next_nd_id) ;
        if (next_node == 0)
         return 0;
        double base_len = (prev_node->coord-this_node->coord).Magnitude() ;
        //printf("FE1\n");
        qdata->otr = FindElementSide(msh_topo,
                                    this_node,prev_node,
                                    next_node,base_len,true) ;
        //printf("FE2\n");
        if (qdata->otr == NO_NODE)
         return 0;

        if (qdata->otr < 0) {
            qdata->qcase = -1 ;
            return 1;
        }

        // check to see if we want to split the triangle
        // edge

        if (!edge_list->ContainsEdge(qdata->obr,qdata->otr)) {
            double base_length =
                   edge_list->GetEdgeLength(edge->id_0,edge->id_1) ;
            ArbIntNode *found = pnode_table->Fetch(qdata->otr) ;
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
                                              ARB_FLOATING,true) ;
                int tri_elem, num_nodes, tri_nodes[3] ;

                tri_elem = msh_topo->GetCCWElem(qdata->obr,qdata->otr) ;
                if (tri_elem == NO_NODE)
                 return 0;
                msh_topo->GetElemNodes(tri_elem,qdata->obr,
                                       &num_nodes,tri_nodes) ;
                if (!msh_topo->DeleteElement(num_nodes,tri_nodes))
                 return 0;

                int ccw_id = tri_nodes[2] ;
                int elem_0 = tri_elem ;
                int elem_1 = NewElemNum() ;

                tri_elem = msh_topo->GetCWElem(qdata->obr,qdata->otr) ;
                if (tri_elem == NO_NODE)
                 return 0;
                msh_topo->GetElemNodes(tri_elem,qdata->otr,
                                       &num_nodes,tri_nodes) ;
                if (!msh_topo->DeleteElement(num_nodes,tri_nodes))
                 return 0;

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

            int cont = 0;
            while (vtx1 != qdata->otr)
            {
                int tmp = vtx1 ;
                if (cont == nelem*3)
                 return 0;
                vtx1 = edge_list->GetCCWNode(vtx0,vtx1) ;
                if (vtx1 == qdata->obr)
                {
                    num = 2 ;
                    break ;
                }
                vtx0 = tmp ;
                cont++;
                num++ ;
            }

            // check for odd number and split the edge.

            if ((num > 5) && ((num % 2) == 1)) {
                ArbIntNode *node_0 = pnode_table->Fetch(qdata->otr) ;
                ArbIntNode *node_1 = pnode_table->Fetch(qdata->obr) ;
                double nx = 0.5 * (node_0->coord[0] + node_1->coord[0]) ;
                double ny = 0.5 * (node_0->coord[1] + node_1->coord[1]) ;
                int new_id = NewNode(nx,ny,INTERIOR,
                                              ARB_FLOATING,true) ;

                int tri_elem, num_nodes, tri_nodes[3] ;

                tri_elem = msh_topo->GetCCWElem(qdata->obr,qdata->otr) ;
                if (tri_elem == NO_NODE)
                 return 0;
                msh_topo->GetElemNodes(tri_elem,qdata->obr,
                                       &num_nodes,tri_nodes) ;
                if (!msh_topo->DeleteElement(num_nodes,tri_nodes))
                 return 0;

                int ccw_id = tri_nodes[2] ;
                int elem_0 = tri_elem ;
                int elem_1 = NewElemNum() ;

                tri_elem = msh_topo->GetCWElem(qdata->obr,qdata->otr) ;
                if (tri_elem == NO_NODE)
                 return 0;
                msh_topo->GetElemNodes(tri_elem,qdata->otr,
                                       &num_nodes,tri_nodes) ;
                if (!msh_topo->DeleteElement(num_nodes,tri_nodes))
                 return 0;

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

    //printf("F2\n");
    if ((end_code == 0) || (end_code == 1)) {

        //  We have either code 0 or code 2, so define the
        //  edge at the right node.
        //
        //             code 0              code 1   *
        //                                          |
        //      *-----*=====*-----*     *-----*=====*

        int prev_nd_id ;
        ArbIntNode *prev_node, *this_node, *next_node ;

        this_node = pnode_table->Fetch(edge->id_0) ;
        if (this_node == 0)
         return 0;
        next_node = pnode_table->Fetch(edge->id_1) ;
        if (next_node == 0)
         return 0;
        prev_nd_id = edge_list->GetCWNode(edge->id_0,edge->id_1) ;
        prev_node = pnode_table->Fetch(prev_nd_id) ;
        if (prev_node == 0)
         return 0;
        double base_len = (next_node->coord-this_node->coord).Magnitude() ;
        qdata->otl = FindElementSide(msh_topo,
                                    this_node,prev_node,
                                    next_node,base_len,false) ;
        if (qdata->otl == NO_NODE)
         return 0;

        if (qdata->otl < 0) {
            qdata->qcase = -1 ;
            return 1;
        }

        // check to see if we want to split the triangle
        // edge

        if (!edge_list->ContainsEdge(qdata->otl,qdata->obl)) {
            double base_length =
                   edge_list->GetEdgeLength(edge->id_0,edge->id_1) ;
            ArbIntNode *found = pnode_table->Fetch(qdata->otl) ;
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
                                              ARB_FLOATING,true) ;

                int tri_elem, num_nodes, tri_nodes[3] ;

                tri_elem = msh_topo->GetCCWElem(qdata->obl,qdata->otl) ;
                if (tri_elem == NO_NODE)
                 return 0;
                msh_topo->GetElemNodes(tri_elem,qdata->obl,
                                       &num_nodes,tri_nodes) ;
                if (!msh_topo->DeleteElement(num_nodes,tri_nodes))
                 return 0;

                int ccw_id = tri_nodes[2] ;
                int elem_0 = tri_elem ;
                int elem_1 = NewElemNum() ;

                tri_elem = msh_topo->GetCWElem(qdata->obl,qdata->otl) ;
                if (tri_elem == NO_NODE)
                 return 0;
                msh_topo->GetElemNodes(tri_elem,qdata->otl,
                                       &num_nodes,tri_nodes) ;
                if (!msh_topo->DeleteElement(num_nodes,tri_nodes))
                 return 0;

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

            CArbTopoAdjVtxIterator iter(msh_topo,qdata->otl) ;
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

            int cont = 0;
            while (vtx1 != qdata->obl)
            {
                int tmp = vtx1 ;
                if (cont == nelem*3)
                 return 0;
                vtx1 = edge_list->GetCCWNode(vtx0,vtx1) ;
                if (vtx1 == qdata->otl)
                {
                    num = 2 ;
                    break ;
                }
                vtx0 = tmp ;
                num++ ;
                cont++;
            }

            // check for odd number and split the edge.

            if ((num > 5) && ((num % 2) == 1)) {
                ArbIntNode *node_0 = pnode_table->Fetch(qdata->otl) ;
                ArbIntNode *node_1 = pnode_table->Fetch(qdata->obl) ;
                double nx = 0.5 * (node_0->coord[0] + node_1->coord[0]) ;
                double ny = 0.5 * (node_0->coord[1] + node_1->coord[1]) ;
                int new_id = NewNode(nx,ny,INTERIOR,
                                              ARB_FLOATING,true) ;

                int tri_elem, num_nodes, tri_nodes[3] ;

                tri_elem = msh_topo->GetCCWElem(qdata->obl,qdata->otl) ;
                if (tri_elem == NO_NODE)
                 return 0;
                msh_topo->GetElemNodes(tri_elem,qdata->obl,
                                       &num_nodes,tri_nodes) ;
                if (!msh_topo->DeleteElement(num_nodes,tri_nodes))
                 return 0;

                int ccw_id = tri_nodes[2] ;
                int elem_0 = tri_elem ;
                int elem_1 = NewElemNum() ;

                tri_elem = msh_topo->GetCWElem(qdata->obl,qdata->otl) ;
                if (tri_elem == NO_NODE)
                 return 0;
                if (!msh_topo->GetElemNodes(tri_elem,qdata->otl,
                                       &num_nodes,tri_nodes))
                  return 0;
                if (!msh_topo->DeleteElement(num_nodes,tri_nodes))
                 return 0;

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

    //printf("F3\n");
    if (!edge_list->ClassifyQuad(qdata->obl,qdata->obr,qdata->otr,qdata->otl,
                            qdata))
     return 0;

    //printf("F4\n");

    if (qdata->qcase == 2) {
        qdata->ratio = edge_list->GetEdgeLength(qdata->bl,qdata->br) /
                       edge_list->GetEdgeLength(qdata->br,qdata->tr) ;
        ArbIntNode *node_0 ;
        ArbIntNode *node_1 ;
        if (qdata->ratio > 1.0) {
            node_0 = pnode_table->Fetch(qdata->bl) ;
            node_1 = pnode_table->Fetch(qdata->br) ;
        } else {
            node_0 = pnode_table->Fetch(qdata->br) ;
            node_1 = pnode_table->Fetch(qdata->tr) ;
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
        ArbIntNode *node_0 ;
        ArbIntNode *node_1 ;
        if (qdata->ratio > 1.0) {
            node_0 = pnode_table->Fetch(qdata->bl) ;
            node_1 = pnode_table->Fetch(qdata->br) ;
        } else {
            node_0 = pnode_table->Fetch(qdata->br) ;
            node_1 = pnode_table->Fetch(qdata->tr) ;
        }
        bool bdry_edge = (node_0->type == BOUNDARY) &&
                         (node_1->type == BOUNDARY) ;
        if (bdry_edge) {
            qdata->valid_transition_split = false ;
        } else {
            qdata->valid_transition_split = true ;
        }

        node_0 = pnode_table->Fetch(qdata->tl) ;
        node_1 = pnode_table->Fetch(qdata->tr) ;
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
    return 1;
}




// %(CArbMshRegion2D::FindSeamData-bool-|-CArbMshEdgeList::ArbQdEdge-const|*-CArbMshEdgeList::ArbNewSeamData-|*-CArbMshEdgeList-|*-CArbMshTopo2D-|*-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    FindSeamData - find data for a seam operation
**
**      bool FindSeamData(
**              const CArbMshEdgeList::ArbQdEdge *edge,
**              CArbMshEdgeList::ArbNewSeamData  *sdata,
**              CArbMshEdgeList                  *edge_list,
**              CArbMshTopo2D                    *msh_topo,
**              CArbMshTopo2D                    *quad_topo)
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

bool CArbMshRegion2D::FindSeamData(const CArbMshEdgeList::ArbQdEdge *edge,
                                   CArbMshEdgeList::ArbNewSeamData *sdata,
                                   CArbMshEdgeList *edge_list,
                                   CArbMshTopo2D *msh_topo,
                                   CArbMshTopo2D *quad_topo)
{
    double small_ang ;

    if (edge->angle_0 < edge->angle_1) {
        sdata->id0 = edge_list->GetCWNode(edge->id_0,edge->id_1) ;
        if (sdata->id0 == -1)
         return false;
        sdata->id1 = edge->id_0 ;
        sdata->id2 = edge->id_1 ;
        sdata->cw_side = true ;
        small_ang = edge->angle_0 ;
    } else {
        sdata->id0 = edge->id_0 ;
        sdata->id1 = edge->id_1 ;
        sdata->id2 = edge_list->GetCCWNode(edge->id_0,edge->id_1) ;
        if (sdata->id2 == -1)
         return false;
        sdata->cw_side = false ;
        small_ang = edge->angle_1 ;
    }
    sdata->cw = edge_list->GetCWNode(sdata->id0,sdata->id1) ;
    if (sdata->cw == -1)
     return false;
    sdata->ccw = edge_list->GetCCWNode(sdata->id1,sdata->id2) ;
    if (sdata->ccw == -1)
     return false;
    sdata->ratio = edge_list->GetEdgeLength(sdata->id0,sdata->id1) /
                   edge_list->GetEdgeLength(sdata->id1,sdata->id2) ;
    sdata->do_template = false ;
    int level_0 = edge_list->GetEdgeLevel(sdata->id0,sdata->id1) ;
    int level_1 = edge_list->GetEdgeLevel(sdata->id1,sdata->id2) ;
    sdata->min_level = (level_0 < level_1) ? level_0 : level_1 ;

    // check to make sure that the end nodes are movable

    ArbIntNode *tid0 = pnode_table->Fetch(sdata->id0) ;
    if (tid0->motion == ARB_FIXED) {
        ArbIntNode *tid2 = pnode_table->Fetch(sdata->id2) ;
        if (tid2->motion == ARB_FIXED) {
            sdata->cw = 0 ;
            sdata->ccw = 1 ;
            return(false) ;
        }
    }

    // check for the case where we want to do a template
    // transition before we do the seam.

    double len_cw = edge_list->GetEdgeLength(sdata->id0,sdata->id1) ;
    double len_ccw = edge_list->GetEdgeLength(sdata->id1,sdata->id2) ;
    ArbIntNode *id1 = pnode_table->Fetch(sdata->id1) ;

    if (len_cw > len_ccw) {
        ArbIntNode *ccw = pnode_table->Fetch(sdata->ccw) ;
        ArbIntNode *id0 = pnode_table->Fetch(sdata->id0) ;
        double angle = Angle(id1->coord,ccw->coord,id0->coord) ;
        if (angle < small_ang) {
            double len = (ccw->coord - id1->coord).Magnitude() ;
            sdata->do_template = (len < len_cw) ? true : false ;
            sdata->template_cw = true ;
        }
    } else {
        ArbIntNode *cw = pnode_table->Fetch(sdata->cw) ;
        ArbIntNode *id2 = pnode_table->Fetch(sdata->id2) ;
        double angle = Angle(id1->coord,id2->coord,cw->coord) ;
        if (angle < small_ang) {
            double len = (cw->coord - id1->coord).Magnitude() ;
            sdata->do_template = (len < len_ccw) ? true : false ;
            sdata->template_cw = false ;
        }
    }

    // try to recover an edge between id0 and id2 and clear the triangle

    if (!RecoverTopEdge(sdata->id2,sdata->id0,msh_topo,edge_list))
        return(false) ;

    int boundary[3] ;
    boundary[0] = sdata->id1 ;
    boundary[1] = sdata->id2 ;
    boundary[2] = sdata->id0 ;

    CArbMshTopo2D::ArbEdge bedge ;
    bedge.nd0 = sdata->id1 ;
    bedge.nd1 = sdata->id2 ;
    bedge.elem0 = msh_topo->GetCCWElem(bedge.nd0,bedge.nd1) ;
    if (bedge.elem0 == NO_NODE)
     return false;
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
        CArbTopoAdjVtxIterator titerq(msh_topo,sdata->id1) ;
        msh_topo->GetElemNodes(titerq.CcwElem(),sdata->id1,
                               &num,triangle) ;
        if (quad_topo->NumAdjElems(sdata->id1) == 2) {
            CArbTopoAdjVtxIterator iterq(quad_topo,sdata->id1) ;
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

    CArbArray<SeamCoordCache> nd0_triangles ;
    CArbArray<SeamCoordCache> nd0_quads ;
    CArbTopoAdjVtxIterator iterm(msh_topo,sdata->id0) ;
    int i ;

    for (iterm.First() ; iterm.More() ; ++iterm) {
        if (iterm.CcwElem() != NO_ELEM) {
            ArbIntNode *adj ;
            SeamCoordCache cache ;
            int num, nodes[3], num_fix = 0 ;
            msh_topo->GetElemNodes(iterm.CcwElem(),
                                   sdata->id0,&num,nodes) ;
            cache.boundary = false ;
            for (i=0 ; i<num ; ++i) {
                adj = pnode_table->Fetch(nodes[i]) ;
                if (adj->motion == ARB_FIXED) ++num_fix ;
                cache.coord[i] = adj->coord ;
            }
            if (num_fix >= 2) cache.boundary = true ;
            cache.ignore = ((nodes[1] == sdata->id2) ||
                            (nodes[num-1] == sdata->id2)) ;
            nd0_triangles.InsertAtEnd(cache) ;
        }
    }

    // do the same for the quad elements

    CArbTopoAdjVtxIterator iterq(quad_topo,sdata->id0) ;
    for (iterq.First() ; iterq.More() ; ++iterq) {
        if (iterq.CcwElem() != NO_ELEM) {
            ArbIntNode *adj ;
            SeamCoordCache cache ;
            int num, nodes[4], num_fix = 0 ;
            quad_topo->GetElemNodes(iterq.CcwElem(),
                                    sdata->id0,
                                    &num,nodes) ;
            cache.boundary = false ;
            cache.ignore = false ;
            for (i=0 ; i<num ; ++i) {
                adj = pnode_table->Fetch(nodes[i]) ;
                if (adj->motion == ARB_FIXED) ++num_fix ;
                cache.coord[i] = adj->coord ;
            }
            if (num_fix >= 2) cache.boundary = true ;
            if (sdata->merge_quads) {
                cache.ignore = ((nodes[1] == sdata->id1) ||
                                (nodes[num-1] == sdata->id1)) ;
            }
            nd0_quads.InsertAtEnd(cache) ;
        }
    }

    // look for the range of positions were we can place this
    // node on the vector between this in the second node (normalized
    // between 0 and 1) and still have all the adjacent elements
    // have valid shapes.

    double range_0, range_1 ;
    CArbCoord2D vect ;
    ArbIntNode *id0 = pnode_table->Fetch(sdata->id0) ;
    ArbIntNode *id2 = pnode_table->Fetch(sdata->id2) ;

    vect = id2->coord - id0->coord ;

    if (id0->motion == ARB_FIXED) {
        CArbCoord2D trial(0.0,0.0) ;
        if (!CheckValidTriangleList(id0->coord,trial,&nd0_triangles) ||
            !CheckValidQuadList(id0->coord,trial,&nd0_quads)) return(false) ;
        range_0 = 0.0 ;
    } else {
        CArbCoord2D trial = vect ;
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

    CArbArray<SeamCoordCache> nd2_quads ;
    CArbArray<SeamCoordCache> nd2_triangles ;

    for (iterm.NewVtx(sdata->id2) ; iterm.More() ; ++iterm) {
        if (iterm.CcwElem() != NO_ELEM) {
            ArbIntNode *adj ;
            int num, nodes[3], num_fix = 0 ;
            SeamCoordCache cache ;
            msh_topo->GetElemNodes(iterm.CcwElem(),
                                   sdata->id2,&num,nodes) ;
            cache.boundary = false ;
            for (i=0 ; i<num ; ++i) {
                adj = pnode_table->Fetch(nodes[i]) ;
                if (adj->motion == ARB_FIXED) ++num_fix ;
                cache.coord[i] = adj->coord ;
            }
            if (num_fix >= 2) cache.boundary = true ;
            cache.ignore = ((nodes[1] == sdata->id0) ||
                            (nodes[num-1] == sdata->id0)) ;
            nd2_triangles.InsertAtEnd(cache) ;
        }
    }

    // do the same for the quad elements

    for (iterq.NewVtx(sdata->id2) ; iterq.More() ; ++iterq) {
        if (iterq.CcwElem() != NO_ELEM) {
            ArbIntNode *adj ;
            int num, nodes[4], num_fix = 0 ;
            SeamCoordCache cache ;
            quad_topo->GetElemNodes(iterq.CcwElem(),
                                    sdata->id2,
                                    &num,nodes) ;
            cache.boundary = false ;
            cache.ignore = false ;
            for (i=0 ; i<num ; ++i) {
                adj = pnode_table->Fetch(nodes[i]) ;
                if (adj->motion == ARB_FIXED) ++num_fix ;
                cache.coord[i] = adj->coord ;
            }
            if (num_fix >= 2) cache.boundary = true ;
            if (sdata->merge_quads) {
                cache.ignore = ((nodes[1] == sdata->id1) ||
                                (nodes[num-1] == sdata->id1)) ;
            }
            nd2_quads.InsertAtEnd(cache) ;
        }
    }


    // look for the range of positions were we can place this
    // node on the vector between this in the first node (normalized
    // between 0 and 1) and still have all the adjacent elements
    // have valid shapes.

    if (id2->motion == ARB_FIXED) {
        CArbCoord2D trial(0.0,0.0) ;
        if (!CheckValidTriangleList(id2->coord,trial,&nd2_triangles) ||
            !CheckValidQuadList(id2->coord,trial,&nd2_quads)) return(false) ;
        range_1 = 1.0 ;
    } else {
        CArbCoord2D trial = -vect ;
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

    if (id0->motion == ARB_FIXED) {
        if (range_1 == 0.0) {
            sdata->coord = id0->coord ;
        } else {
            return(false) ;
        }
    } else if (id2->motion == ARB_FIXED) {
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




// %(CArbMshRegion2D::DoSeam-bool-|-CArbMshEdgeList::ArbNewSeamData-|*-CArbMshTopo2D-|*-CArbMshTopo2D-|*-CArbMshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    DoSeam - do a seam update
**
**      bool DoSeam(
**              CArbMshEdgeList::ArbNewSeamData *sdata,
**              CArbMshTopo2D                   *msh_topo,
**              CArbMshTopo2D                   *quad_topo,
**              CArbMshEdgeList                 *edge_list)
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

bool CArbMshRegion2D::DoSeam(CArbMshEdgeList::ArbNewSeamData *sdata,
                             CArbMshTopo2D *msh_topo,
                             CArbMshTopo2D *quad_topo,
                             CArbMshEdgeList *edge_list)
{
    if (!RecoverTopEdge(sdata->id2,sdata->id0,msh_topo,edge_list))
        return(false) ;

    // now we want to find node that forms a triangle with
    // node 1 and the cw node.  To do this we loop around
    // the cw node in the ccw direction.  Once we find
    // node 1, the node we want is the next one.

    CArbTopoAdjVtxIterator iter0(msh_topo,sdata->id0) ;

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
    CArbMshTopo2D::ArbEdge bedge ;
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
        CArbMshTopo2D::ArbEdge tmp_edge0 = { opp, save_id2, 0, 0 } ;
        int vtx0 = msh_topo->OppositeNode(elem0,&tmp_edge0) ;

        CArbMshTopo2D::ArbEdge tmp_edge1 = { save_id0, opp, 0, 0 } ;
        int vtx1 = msh_topo->OppositeNode(elem1,&tmp_edge1) ;

        if (vtx0 == vtx1) {
            msh_topo->DeleteTriangle(opp,save_id2,vtx0) ;
            msh_topo->DeleteTriangle(save_id0,opp,vtx0) ;
        }
    }

    // update the triangular mesh topology

    int i ;
    CArbArray<SeamElemCache> elem_cache ;

    for (iter0.NewVtx(remove) ; iter0.More() ; ++iter0) {
        if (iter0.CcwElem() != NO_ELEM) {
            SeamElemCache edata ;
            edata.elem = iter0.CcwElem() ;
            msh_topo->GetElemNodes(iter0.CcwElem(),
                    remove,&edata.num_nodes,edata.nodes) ;
            elem_cache.InsertAtEnd(edata) ;
        }
    }
    for (i=0 ; i<elem_cache.NumEntries() ; ++i) {
        msh_topo->DeleteElement(elem_cache[i].num_nodes,
                                elem_cache[i].nodes) ;
    }
    for (i=0 ; i<elem_cache.NumEntries() ; ++i) {
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
    CArbTopoAdjVtxIterator iter1(quad_topo,remove) ;
    for (iter1.First() ; iter1.More() ; ++iter1) {
        if (iter1.CcwElem() != NO_ELEM) {
            SeamElemCache edata ;
            edata.elem = iter1.CcwElem() ;
            quad_topo->GetElemNodes(iter1.CcwElem(),
                remove,&edata.num_nodes,edata.nodes) ;
            elem_cache.InsertAtEnd(edata) ;
        }
    }
    for (i=0 ; i<elem_cache.NumEntries() ; ++i) {
        quad_topo->DeleteElement(elem_cache[i].num_nodes,
                                 elem_cache[i].nodes) ;
    }
    for (i=0 ; i<elem_cache.NumEntries() ; ++i) {
        elem_cache[i].nodes[0] = keep ;
        quad_topo->InsertElement(elem_cache[i].elem,
                                 elem_cache[i].num_nodes,
                                 elem_cache[i].nodes) ;
    }

    // place the node

    ArbIntNode *node_0 = pnode_table->Fetch(keep) ;
    node_0->coord = sdata->coord ;

    // update the quad and triangle topology

    SmoothOneAdj(sdata->id2,sdata->id1,sdata->ccw,quad_topo,
                 false,false,edge_list) ;
    SmoothOneAdj(sdata->id0,sdata->id1,sdata->cw,quad_topo,
                 false,false,edge_list) ;

    return(true) ;
}



// %(CArbMshRegion2D::DoTransitionSeam-bool-|-CArbMshEdgeList::ArbQdEdge-|*-CArbMshEdgeList::ArbNewSeamData-|*-CArbMshTopo2D-|*-CArbMshTopo2D-|*-CArbMshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    DoTransitionSeam - do a transition seam
**
**      bool DoTransitionSeam(
**              CArbMshEdgeList::ArbQdEdge      *edge,
**              CArbMshEdgeList::ArbNewSeamData *sdata,
**              CArbMshTopo2D                   *msh_topo,
**              CArbMshTopo2D                   *quad_topo,
**              CArbMshEdgeList                 *edge_list)
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

bool CArbMshRegion2D::DoTransitionSeam(
                             CArbMshEdgeList::ArbQdEdge *edge,
                             CArbMshEdgeList::ArbNewSeamData *sdata,
                             CArbMshTopo2D *msh_topo,
                             CArbMshTopo2D *quad_topo,
                             CArbMshEdgeList *edge_list)
{
    if (sdata->ratio > 1.0) {

        ArbIntNode *node_0 = pnode_table->Fetch(sdata->id0) ;
        ArbIntNode *node_1 = pnode_table->Fetch(sdata->id1) ;
        if ((node_0->type == BOUNDARY) &&
            (node_1->type == BOUNDARY)) return(false) ;

        // delete the quad adjacent to the split edge

        int num_nodes ;
        int quad_nodes[4] ;

        int quad_elem =
            quad_topo->GetCCWElem(sdata->id1,sdata->id0) ;
        quad_topo->GetElemNodes(quad_elem,sdata->id0,
                                &num_nodes,quad_nodes) ;

        if (edge_list->ContainsEdge(sdata->id1,quad_nodes[2]))
            return(false) ;

        // check to make shure that this element has a valid shape

        ArbIntNode *nd0 = pnode_table->Fetch(quad_nodes[0]) ;
        ArbIntNode *nd1 = pnode_table->Fetch(quad_nodes[1]) ;
        ArbIntNode *nd2 = pnode_table->Fetch(quad_nodes[2]) ;
        ArbIntNode *nd3 = pnode_table->Fetch(quad_nodes[3]) ;
        double metric = QuadMetric(nd0->coord,nd1->coord,
                                   nd2->coord,nd3->coord) ;
        if (metric <= 0.0) return(false) ;
        quad_topo->DeleteElement(num_nodes,quad_nodes) ;

        // create the new node

        double nx = 0.5 * (node_0->coord[0] + node_1->coord[0]) ;
        double ny = 0.5 * (node_0->coord[1] + node_1->coord[1]) ;
        int new_id = NewNode(nx,ny,INTERIOR,ARB_FLOATING,true) ;

        // add the new quad

        quad_nodes[3] = new_id ;
        quad_topo->InsertElement(quad_elem,4,quad_nodes) ;

        // delete the triangle adjacent to the split edge

        int tri_nodes[3] ;

        int tri_elem =
            msh_topo->GetCCWElem(sdata->id0,sdata->id1) ;
        msh_topo->GetElemNodes(tri_elem,sdata->id1,
                               &num_nodes,tri_nodes) ;
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

        ArbIntNode *node_1 = pnode_table->Fetch(sdata->id1) ;
        ArbIntNode *node_2 = pnode_table->Fetch(sdata->id2) ;
        if ((node_1->type == BOUNDARY) &&
            (node_2->type == BOUNDARY)) return(false) ;

        // find the nodes on the quad that we are going to delete

        int num_nodes ;
        int quad_nodes[4] ;

        int quad_elem =
            quad_topo->GetCCWElem(sdata->id2,sdata->id1) ;
        quad_topo->GetElemNodes(quad_elem,sdata->id1,
                                &num_nodes,quad_nodes) ;

        if (edge_list->ContainsEdge(quad_nodes[1],sdata->id1))
            return(false) ;

        // check to make shure that this element has a valid shape

        ArbIntNode *nd0 = pnode_table->Fetch(quad_nodes[0]) ;
        ArbIntNode *nd1 = pnode_table->Fetch(quad_nodes[1]) ;
        ArbIntNode *nd2 = pnode_table->Fetch(quad_nodes[2]) ;
        ArbIntNode *nd3 = pnode_table->Fetch(quad_nodes[3]) ;
        double metric = QuadMetric(nd0->coord,nd1->coord,
                                   nd2->coord,nd3->coord) ;
        if (metric <= 0.0) return(false) ;
        quad_topo->DeleteElement(num_nodes,quad_nodes) ;

        // create the new node

        double nx = 0.5 * (node_1->coord[0] + node_2->coord[0]) ;
        double ny = 0.5 * (node_1->coord[1] + node_2->coord[1]) ;
        int new_id = NewNode(nx,ny,INTERIOR,ARB_FLOATING,true) ;

        // add the new quad

        quad_nodes[0] = new_id ;
        quad_topo->InsertElement(quad_elem,num_nodes,quad_nodes) ;

        // delete the triangle adjacent to the split edge

        int tri_nodes[3] ;

        int tri_elem =
            msh_topo->GetCCWElem(sdata->id1,sdata->id2) ;
        msh_topo->GetElemNodes(tri_elem,sdata->id2,
                               &num_nodes,tri_nodes) ;
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




// %(CArbMshRegion2D::DoTransitionSplit-void-|-CArbMshEdgeList::ArbQdEdge-|*-CArbMshEdgeList::ArbNewQuadData-|*-CArbMshTopo2D-|*-CArbMshTopo2D-|*-CArbMshEdgeList-|*)
/* ++ ----------------------------------------------------------
**
**    DoTransitionSplit - do a transition split
**
**      void DoTransitionSplit(
**              CArbMshEdgeList::ArbQdEdge      *edge,
**              CArbMshEdgeList::ArbNewQuadData *qdata,
**              CArbMshTopo2D                   *msh_topo,
**              CArbMshTopo2D                   *quad_topo,
**              CArbMshEdgeList                 *edge_list)
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

int CArbMshRegion2D::DoTransitionSplit(
                             CArbMshEdgeList::ArbQdEdge *edge,
                             CArbMshEdgeList::ArbNewQuadData *qdata,
                             CArbMshTopo2D *msh_topo,
                             CArbMshTopo2D *quad_topo,
                             CArbMshEdgeList *edge_list)
{
    if (qdata->ratio > 1.0) {

        // delete the quad adjacent to the split edge

        int num_nodes ;
        int quad_nodes[4] ;

        int quad_elem =
            quad_topo->GetCCWElem(qdata->br,qdata->bl) ;
        if (quad_elem == NO_NODE)
         return 0;
        quad_topo->GetElemNodes(quad_elem,qdata->bl,
                                &num_nodes,quad_nodes) ;
        quad_topo->DeleteElement(num_nodes,quad_nodes) ;

        // create the two new nodes

        ArbIntNode *node_0 = pnode_table->Fetch(qdata->br) ;
        ArbIntNode *node_1 = pnode_table->Fetch(qdata->bl) ;
        double nx = 0.5 * (node_0->coord[0] + node_1->coord[0]) ;
        double ny = 0.5 * (node_0->coord[1] + node_1->coord[1]) ;
        int new_id0 = NewNode(nx,ny,INTERIOR,ARB_FLOATING,true) ;

        node_1 = pnode_table->Fetch(quad_nodes[1]) ;
        nx = 0.5 * (node_0->coord[0] + node_1->coord[0]) ;
        ny = 0.5 * (node_0->coord[1] + node_1->coord[1]) ;
        int new_id1 = NewNode(nx,ny,INTERIOR,ARB_FLOATING,true) ;

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
        if (tri_elem == NO_NODE)
         return 0;
        msh_topo->GetElemNodes(tri_elem,qdata->br,
                               &num_nodes,tri_nodes) ;
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
        if (quad_elem == NO_NODE)
         return 0;
        quad_topo->GetElemNodes(quad_elem,qdata->br,
                                &num_nodes,quad_nodes) ;
        quad_topo->DeleteElement(num_nodes,quad_nodes) ;

        // create the new node

        ArbIntNode *node_1 = pnode_table->Fetch(qdata->br) ;
        ArbIntNode *node_0 = pnode_table->Fetch(qdata->tr) ;
        double nx = 0.5 * (node_1->coord[0] + node_0->coord[0]) ;
        double ny = 0.5 * (node_1->coord[1] + node_0->coord[1]) ;
        int new_id0 = NewNode(nx,ny,INTERIOR,ARB_FLOATING,true) ;

        node_0 = pnode_table->Fetch(quad_nodes[2]) ;
        nx = 0.5 * (node_1->coord[0] + node_0->coord[0]) ;
        ny = 0.5 * (node_1->coord[1] + node_0->coord[1]) ;
        int new_id1 = NewNode(nx,ny,INTERIOR,ARB_FLOATING,true) ;

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
        if (tri_elem == NO_NODE)
         return 0;
        msh_topo->GetElemNodes(tri_elem,qdata->tr,
                               &num_nodes,tri_nodes) ;
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
    return 1;
}




// %(CArbMshRegion2D::DoTemplateSplit-void-|-CArbMshEdgeList::ArbNewQuadData-|*-CArbMshTopo2D-|*-CArbMshTopo2D-|*-CArbMshEdgeList-|*-bool-const|)
/* ++ ----------------------------------------------------------
**
**    DoTemplateSplit - do a template split
**
**      void DoTemplateSplit(
**              CArbMshEdgeList::ArbNewQuadData *qdata,
**              CArbMshTopo2D                   *msh_topo,
**              CArbMshTopo2D                   *quad_topo,
**              CArbMshEdgeList                 *edge_list,
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

int CArbMshRegion2D::DoTemplateSplit(
                             CArbMshEdgeList::ArbNewQuadData *qdata,
                             CArbMshTopo2D *msh_topo,
                             CArbMshTopo2D *quad_topo,
                             CArbMshEdgeList *edge_list,
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
    if (quad_elem == -1)
     return 0;
    if (!quad_topo->GetElemNodes(quad_elem,qdata->bl,&num_nodes,quad_nodes))
     return 0;
    quad_topo->DeleteElement(num_nodes,quad_nodes) ;

    // create the four new nodes using the bilinear shape
    // functions to determine the new coordinates

    ArbIntNode *br = pnode_table->Fetch(qdata->br) ;
    ArbIntNode *bl = pnode_table->Fetch(qdata->bl) ;
    ArbIntNode *tr = pnode_table->Fetch(qdata->tr) ;
    ArbIntNode *tl = pnode_table->Fetch(qdata->tl) ;

    if (br == 0 || bl == 0 || tr == 0 || tl == 0)
     return 0;

    double nx, ny ;

    nx = (2.0*tr->coord[0] + tl->coord[0]) / 3.0 ;
    ny = (2.0*tr->coord[1] + tl->coord[1]) / 3.0 ;
    int new_id0 = NewNode(nx,ny,INTERIOR,ARB_FLOATING,true) ;

    nx = (tr->coord[0] + 2.0*tl->coord[0]) / 3.0 ;
    ny = (tr->coord[1] + 2.0*tl->coord[1]) / 3.0 ;
    int new_id1 = NewNode(nx,ny,INTERIOR,ARB_FLOATING,true) ;

    nx = (2.0*bl->coord[0] + br->coord[0] +
          tr->coord[0] + 2.0*tl->coord[0]) / 6.0 ;
    ny = (2.0*bl->coord[1] + br->coord[1] +
          tr->coord[1] + 2.0*tl->coord[1]) / 6.0 ;
    int new_id2 = NewNode(nx,ny,INTERIOR,ARB_FLOATING,true) ;

    nx = (bl->coord[0] + 2.0*br->coord[0] +
          2.0*tr->coord[0] + tl->coord[0]) / 6.0 ;
    ny = (bl->coord[1] + 2.0*br->coord[1] +
          2.0*tr->coord[1] + tl->coord[1]) / 6.0 ;
    int new_id3 = NewNode(nx,ny,INTERIOR,ARB_FLOATING,true) ;

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
    int tri_elem = msh_topo->GetCCWElem(qdata->tl,qdata->tr) ;
    if (tri_elem == -1)
     return 0;
    if (!msh_topo->GetElemNodes(tri_elem,qdata->tr,&num_nodes,tri_nodes))
     return 0;
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
    return 1;
}


bool CArbMshRegion2D::QuadsInPolygon(int num,
                        const int *ids,
                        const CArbCoord2D *vts,
                        const CArbMshTopo2D *quad_topo,
                        const CArbMshEdgeList *edge_list) const
{
    int i ;

    // here we set up data structures to look for mate nodes.  This
    // is an after the fact fix so the logic is a bit strange.  We
    // assume that the number of verts is <= 6.  If this is not
    // the case then we create more space

    CArbSmallSet<int,2> *sets[6] = {0,0,0,0,0,0} ;
    CArbSmallSet<int,2> **setptr = &sets[0] ;
    bool have_sets = false ;

    if (num > 6) {
        setptr = new CArbSmallSet<int,2> *[num] ;
    }

    // determine a tolerance and store the max and min
    // extents

    CArbCoord2D pmin = vts[0] ;
    CArbCoord2D pmax = vts[0] ;

    double tol = 0 ;
    for (i=0 ; i<num ; ++i) {
        int j = (i+1) % num ;
        tol += (vts[j]-vts[i]).Magnitude() ;
        pmin = min(vts[i],pmin) ;
        pmax = max(vts[i],pmax) ;
        setptr[i] = MateTable->Fetch(ids[i]) ;
        if (setptr[i] != 0) {
            have_sets = true ;
        }
    }
    tol = 0.001 * tol/double(num) ;

    // loop through all verts in the mesh.  If they are not
    // part of the polygon check to see if they are inside.

//    CArbTopoVtxIterator iter(*quad_topo) ;

    CArbEdgeListNodeIterator iter(edge_list) ;
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

        CArbCoord2D pt = (pnode_table->Fetch(*iter))->coord ;
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


// %(CArbMshRegion2D::RenumberNodes-void-|)
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

void CArbMshRegion2D::RenumberNodes()
{
    CArbHashTable<int,int> node_map ;

    CArbHashTableIterator<int,ArbMshElement2D> eiter(pelem_table) ;

    int j ;

    // first build a hash table that maps from the old number
    // to the new number for the nodes that are actually used.

    for (eiter.First() ; eiter.More() ; ++eiter) {
        ArbMshElement2D *elem = eiter.Entry() ;
        for (j=0 ; j<elem->num_nodes ; ++j) {
            if (elem->nodes[j] > MaxId) {
                if (node_map.Fetch(elem->nodes[j]) == 0) {
                    node_map.Store(elem->nodes[j],StartIdSave) ;
                    StartIdSave++ ;
                }
            }
        }
    }

    // now loop again and update the node numbers in all
    // elements

    for (eiter.First() ; eiter.More() ; ++eiter) {
        ArbMshElement2D *elem = eiter.Entry() ;
        for (j=0 ; j<elem->num_nodes ; ++j) {
            if (elem->nodes[j] > MaxId) {
                elem->nodes[j] =
                    *(node_map.Fetch(elem->nodes[j])) ;
            }
        }
    }

    // build a new node list

    CArbHashTable<int,ArbIntNode> *tmp_node =
        new CArbHashTable<int,ArbIntNode>() ;

    CArbHashTableIterator<int,ArbIntNode> iter(pnode_table) ;

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.Entry()->id <= MaxId) {
            tmp_node->Store(iter.Entry()->id,*(iter.Entry())) ;
        } else {
            int *nid = node_map.Fetch(iter.Entry()->id) ;
            if (nid != 0) {
                ArbIntNode tmp = *(iter.Entry()) ;
                tmp.id = *nid ;
                tmp_node->Store(tmp.id,tmp) ;
            }
        }
    }
    delete pnode_table ;
    pnode_table = tmp_node ;
}




// %(CArbMshRegion2D::IntersectLines-CArbCoord2D-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    IntersectLines - find the intersection point of two lines
**
**      CArbCoord2D IntersectLines(
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
**      Description: This method finds the intersection point of two
**          line segments.
**
**      Return Value: the intersection coordinates
**
**
** -- */

CArbCoord2D CArbMshRegion2D::IntersectLines(
         CArbCoord2D i1,
         CArbCoord2D i2,
         CArbCoord2D j1,
         CArbCoord2D j2) const
{
    CArbCoord2D intsc ;

    /* This is the assumed node configuration:

           \I1   /J2
            \   /
             \ /
              \
             / \
            /   \
           /J1   \I2
    */

    CArbCoord2D bi = i2 - i1 ;
    CArbCoord2D bj = j2 - j1 ;

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

/* ------------------------------------------------------------
    LenSqr - computes the square of the distance between to points
*/

inline double LenSqr(CArbCoord2D i,CArbCoord2D j)
{
    CArbCoord2D delta = i - j ;
    return(delta.x()*delta.x() + delta.y()*delta.y()) ;
}




// %(CArbMshRegion2D::TriMetric-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    TriMetric - compute a shape metric
**
**      double TriMetric(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j) const
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

double CArbMshRegion2D::TriMetric(
        CArbCoord2D i,
        CArbCoord2D j,
        CArbCoord2D k) const
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




// %(CArbMshRegion2D::QuadMetric-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    QuadMetric - compute a shape metric
**
**      double QuadMetric(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j,
**              CArbCoord2D l) const
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

double CArbMshRegion2D::QuadMetric(
        CArbCoord2D i,
        CArbCoord2D j,
        CArbCoord2D k,
        CArbCoord2D l) const
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

double CArbMshRegion2D::QuadMetricAreaOnly(
        CArbCoord2D i,
        CArbCoord2D j,
        CArbCoord2D k,
        CArbCoord2D l) const
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
