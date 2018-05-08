//
// CArbEdgeList Class definition
//
// Description -
//   This class maintains the edge list used in the Q-morph
//   algorithm
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.20 $  $Date: 2002/08/22 19:30:23 $  $Author: wash $
//

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "ArbMshEdgeList.hpp"
#include "ArbHashTable.hpp"
#include "ArbSet.cpp"

#ifdef MEMDEBUG
#include "MemDbg.hpp"
#define new new(__FILE__,__LINE__)
#endif

#define QUAD_EDGE_ANGLE_TOLERANCE 2.356194490             // 3*pi/4
#define QUAD_BOUNDARY_EDGE_ANGLE_TOLERANCE 1.963495408    // 5*pi/8
#define QUAD_EDGE_SEAM_TOLERANCE  0.785398163             // pi/4
//#define ASSUME_OVERLAP_TOLERANCE  6.108652382           // 350 degrees
//#define OVERLAP_ANGLE_VALUE       0.174532925           // -10 degrees
#define ASSUME_OVERLAP_TOLERANCE  6.195918845             // 355 degrees
#define OVERLAP_ANGLE_VALUE       0.087266463             // -5 degrees
#define QUAD_NEAR_ANGLE_TOLERANCE 0.523598776             // pi/6

#define SEAM_ANGLE_OVERRIDE_BDRY_TOL 0.035  // 2 degrees

#define LEFT_EDGE_MASK  10
#define RIGHT_EDGE_MASK  5
#define ANGLE_CODE 1
#define SEAM_CODE  5

#define PI      3.141592654
#define TWO_PI  6.283185307
#define HALF_PI 1.570796327

int operator == (const ArbEdgeKey &op1,const ArbEdgeKey &op2)
{
    return((op1.id_0 == op2.id_0) && (op1.id_1 == op2.id_1)) ;
}


/* ------------------------------------------------------------
    CArbCmpQdEdge - comparison function used to order edges
                    in the list so that the best candidate for
                    creating a new quad can be selected
*/

static int ArbCmpQdEdge(const CArbMshEdgeList::ArbIntQdEdgeDesc &edg1,
                        const CArbMshEdgeList::ArbIntQdEdgeDesc &edg2)
{
    // if we have a realy small seam then do this first

    if (((edg1.angle_0 < SEAM_ANGLE_OVERRIDE_BDRY_TOL) ||
         (edg1.angle_1 < SEAM_ANGLE_OVERRIDE_BDRY_TOL)) &&
         (edg1.level < edg2.level+3)) return(-1) ;
    if (((edg2.angle_0 < SEAM_ANGLE_OVERRIDE_BDRY_TOL) ||
         (edg2.angle_1 < SEAM_ANGLE_OVERRIDE_BDRY_TOL)) &&
         (edg2.level < edg1.level+3)) return(1) ;

    // otherwise do the original boundary first, because these
    // nodes don't float, we want to get away from them as
    // soon as possible

    if ((edg1.level > 0) && (edg2.level == 0)) return(1) ;
    if ((edg1.level == 0) && (edg2.level > 0)) return(-1) ;

//    if (edg1.level+1 > edg2.level) return(1) ;
//    if (edg1.level+1 < edg2.level) return(-1) ;
//    if (edg1.level+2 > edg2.level) return(1) ;
//    if (edg1.level+2 < edg2.level) return(-1) ;

//    int ldiff = edg1.level-edg2.level ;   // GNU compiler bug
//    if (ldiff < 0) ldiff = -ldiff ;            //
//    if (ldiff >= 2) {
    if (abs(edg1.level-edg2.level) >= 2) {
        if (edg1.level > edg2.level) return(1) ;
        if (edg1.level < edg2.level) return(-1) ;
    }

    // seams first

    if ((edg1.end_code & 0xc) < (edg2.end_code & 0xc)) return(1) ;
    if ((edg1.end_code & 0xc) > (edg2.end_code & 0xc)) return(-1) ;

//    if (edg1.level+1 > edg2.level) return(1) ;
//    if (edg1.level+1 < edg2.level) return(-1) ;
    if (edg1.level > edg2.level) return(1) ;
    if (edg1.level < edg2.level) return(-1) ;

    // the codes, higher code means lower value

    if (edg1.end_code < edg2.end_code) return(1) ;
    if (edg1.end_code > edg2.end_code) return(-1) ;

    // look at the level in the mesh

    if (edg1.level > edg2.level) return(1) ;
    if (edg1.level < edg2.level) return(-1) ;

    // the levels of adjacent elements

    if (edg1.adj_level_0 > edg2.adj_level_0) return(1) ;
    if (edg1.adj_level_0 < edg2.adj_level_0) return(-1) ;

    if (edg1.adj_level_1 > edg2.adj_level_1) return(1) ;
    if (edg1.adj_level_1 < edg2.adj_level_1) return(-1) ;

    // then the the edge lengths ;

    if (edg1.length > edg2.length) return(1) ;
    if (edg1.length < edg2.length) return(-1) ;
    return(0) ;
}




// %(CArbMshEdgeList::CArbMshEdgeList-constructor-|-CArbHashTable-|<int,ArbIntNode>*)
/* ++ ----------------------------------------------------------
**
**    CArbMshEdgeList - edge list constructor
**
**      CArbMshEdgeList(CArbHashTable<int,ArbIntNode>* inode_table)
**
**        inode_table - (in)  node hash table
**
**      Description: This is a constructor for an edge list. As an
**          argument it takes a pointer to a hash table that maps node
**          id's to ArbIntNodes.
**
**
** -- */

CArbMshEdgeList::CArbMshEdgeList(
    CArbHashTable<int,ArbIntNode> *inode_table) :
    node_table(inode_table)
{
    edge_table = new CArbHashTable<ArbEdgeKey,
              CArbHeap<ArbIntQdEdgeDesc>::EntryHandle>(false) ;
    order_heap = new CArbHeap<ArbIntQdEdgeDesc>(ArbCmpQdEdge) ;
    edge_lengths = new CArbSet<double>(ArbCmpDouble,true) ;
    front_nodes = new CArbHashTable<int,int>(true) ;
    bdry_topo = new CArbMshTopo2D() ;
}




// %(CArbMshEdgeList::CArbMshEdgeList-destructor-|~)
/* ++ ----------------------------------------------------------
**
**    CArbMshEdgeList - edge list destructor
**
**      ~CArbMshEdgeList()
**
**      Description: This is a desctructor for an edge list.
**
**
** -- */

CArbMshEdgeList::~CArbMshEdgeList()
{
    delete edge_table ;
    delete order_heap ;
    delete edge_lengths ;
    delete front_nodes ;
    delete bdry_topo ;
}




// %(CArbMshEdgeList::InsertEdge-void-|-int-const|-int-const|-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    InsertEdge - insert an edge into the list
**
**      void InsertEdge(
**              const int id_0,
**              const int id_1,
**              const int cw_id,
**              const int ccw_id)
**
**        id_0   - (in)  first node id
**        id_1   - (in)  second node id
**        cw_id  - (in)  id of node cw from id_0
**        ccw_id - (in)  if of node ccw from id_1
**
**      Description: This function inserts a new edge into the edge
**          list.
**
**
** -- */

void CArbMshEdgeList::InsertEdge(const int id_0,
                                 const int id_1,
                                 const int /*cw_id*/,
                                 const int ccw_id)
{
    // insert this information into the edge and back
    // tables.

    edge_table->Store(ArbEdgeKey(id_0,id_1),0) ;
    bdry_topo->InsertAngle(1,id_1,ccw_id,id_0) ;
    front_nodes->Store(id_0,id_1) ;
}




// %(CArbMshEdgeList::InitializeCodes-void-|)
/* ++ ----------------------------------------------------------
**
**    InitializeCodes - initialize the end codes for all edges
**
**      void InitializeCodes()
**
**      Description: This function goes through all edges currently in
**          the list and computes the angles they make with neighbors
**          and store this information in the heap that will be used to
**          select an edge from which an element will be constructed.
**
**
** -- */

int CArbMshEdgeList::InitializeCodes()
{
    CArbHashTableIterator<ArbEdgeKey,EdgeHandle> iter(edge_table) ;

    // go through all the edges.  Compute the angle codes
    // and length and insert this information in the selection
    // heap

    for (iter.First() ; iter.More() ; ++iter) {
        ArbIntQdEdgeDesc edge_desc ;
        ArbIntNode *p_nd0, *p_nd1, *adj ;
        double dx, dy ;

        // get pointers to the node info and initialize the
        // description

        p_nd0 = node_table->Fetch(iter.Key().id_0) ;
        p_nd1 = node_table->Fetch(iter.Key().id_1) ;

        if (p_nd0 == 0 || p_nd1 == 0)
         return 0;

        edge_desc.id_0 = iter.Key().id_0 ;
        edge_desc.id_1 = iter.Key().id_1 ;
        edge_desc.level = 0 ;
        edge_desc.end_code = 0 ;
        edge_desc.adj_level_0 = 0 ;
        edge_desc.adj_level_1 = 0 ;
        edge_desc.visited = false ;

        dx = p_nd1->coord[0] - p_nd0->coord[0] ;
        dy = p_nd1->coord[1] - p_nd0->coord[1] ;
        edge_desc.length = sqrt(dx*dx + dy*dy) ;

        // find the included angle at node 0 and classify the
        // the vertex

        adj = node_table->Fetch(bdry_topo->GetCWBdryNode(iter.Key().id_0,
                                                         iter.Key().id_1)) ;
        if (adj == 0)
         return 0;

        edge_desc.angle_0 = Angle2Pi(p_nd0->coord,p_nd1->coord,adj->coord) ;
        if (edge_desc.angle_0 < QUAD_EDGE_ANGLE_TOLERANCE)
            edge_desc.end_code = 2 ;

        // find the include angle at node 1 and classify the
        // the vertex

        adj = node_table->Fetch(bdry_topo->GetCCWBdryNode(iter.Key().id_0,
                                                          iter.Key().id_1)) ;
        if (adj == 0)
         return 0;
        edge_desc.angle_1 = Angle2Pi(p_nd1->coord,adj->coord,p_nd0->coord) ;
        if (edge_desc.angle_1 < QUAD_EDGE_ANGLE_TOLERANCE)
            edge_desc.end_code += 1 ;

        *(iter.Entry()) = order_heap->InsertWithHandle(edge_desc) ;
        edge_lengths->Insert(edge_desc.length) ;
    }
    return 1;
}

/* ------------------------------------------------------------
    Rotate - rotate nodes of a new element into the cannonical
             orientation
*/

static void Rotate(const int rotation,
                   const int in0,
                   const int in1,
                   const int in2,
                   const int in3,
                   int *out0,
                   int *out1,
                   int *out2,
                   int *out3)
{
    switch (rotation) {
        case 0:
            *out0 = in0 ;  *out1 = in1 ;
            *out2 = in2 ;  *out3 = in3 ;
            break ;
        case 1:
            *out0 = in1 ;  *out1 = in2 ;
            *out2 = in3 ;  *out3 = in0 ;
            break ;
         case 2:
            *out0 = in2 ;  *out1 = in3 ;
            *out2 = in0 ;  *out3 = in1 ;
            break ;
         case 3:
            *out0 = in3 ;  *out1 = in0 ;
            *out2 = in1 ;  *out3 = in2 ;
            break ;
    }
}




// %(CArbMshEdgeList::ClassifyQuad-void-|-int-const|-int-const|-int-const|-int-const|-ArbNewQuadData-|*)
/* ++ ----------------------------------------------------------
**
**    ClassifyQuad - find information for a new quadrilateral
**
**      void ClassifyQuad(
**              const int      id_0,
**              const int      id_1,
**              const int      id_2,
**              const int      id_3,
**              ArbNewQuadData *qdata)
**
**        id_0  - (out) first node
**        id_1  - (out) second node
**        id_2  - (out) third node
**        id_3  - (out) forth node
**        qdata - (out) quadrilateral information
**
**      Description: Given the id's for the four corner nodes for a new
**          quadrilatral, this function fills the qdata argument with
**          information about the new quadrilateral.
**
**
** -- */

int CArbMshEdgeList::ClassifyQuad(const int id0,
                                   const int id1,
                                   const int id2,
                                   const int id3,
                                   ArbNewQuadData *qdata)
{
    /*
        We have one of the following 3 cases, where "+"
        indicates new edges:

                                     |
                                     |
        *++++++*              *++++++*     ----*++++++*----
        +      +     	      +      |	       |      |
        +      +              +      |	       |      |
        +      +      	      +      |	       |      |
    ----*------*----  	  ----*------*         *------*

         case 1 	      case 2	          case 3

        The cononical node       top_left *----* top_right
        ordering is shown           (tl)  |    |  (tr)
        at the right                      |    |
                                base_left *----* base_right
                                    (bl)           (br)


         or case five:

    ----*------*----
        +      +
        +      +
        +      +
    ----*------*----


    */

    // do the special case of a triangle

    if ((id0 == id1) || (id1 == id2) || (id2 == id3) || (id3 == id0)) {
        qdata->qcase = 6 ;
        qdata->bl = id0 ;
        qdata->br = id1 ;
        qdata->tl = id2 ;
        qdata->tr = id3 ;
        return 1;
    }

    // first determine which of the edges exist so that
    // we know which case to deal with

    bool exists[4] ;
    int ii ;

    for (ii=0 ; ii<4 ; ++ii) exists[ii] = false ;

    if (edge_table->Fetch(ArbEdgeKey(id0,id1)) != 0) exists[0] = true ;
    if (edge_table->Fetch(ArbEdgeKey(id1,id2)) != 0) exists[1] = true ;
    if (edge_table->Fetch(ArbEdgeKey(id2,id3)) != 0) exists[2] = true ;
    if (edge_table->Fetch(ArbEdgeKey(id3,id0)) != 0) exists[3] = true ;

    qdata->qcase = 0 ;
    for (ii=0 ; ii<4 ; ++ii) {
        if (exists[ii] == true) qdata->qcase++ ;
    }

    // now that we know the case, swap the vertices to get
    // the elements into a cannonical orientation as shown
    // in the picture above

    switch (qdata->qcase) {
        case 1:
            for (ii=0 ; ii<4 ; ++ii) {
                if (exists[ii]) {
                    Rotate(ii,id0,id1,id2,id3,&qdata->bl,
                           &qdata->br,&qdata->tr,&qdata->tl) ;
                    break ;
                }
            }
            qdata->num_front_nodes = 6 ;
            qdata->front_nodes[0] =
                bdry_topo->GetCWBdryNode(qdata->bl,qdata->br) ;
            if (qdata->front_nodes[0] == -1)
             return 0;
            qdata->front_nodes[1] = qdata->bl ;
            qdata->front_nodes[2] = qdata->tl ;
            qdata->front_nodes[3] = qdata->tr ;
            qdata->front_nodes[4] = qdata->br ;
            qdata->front_nodes[5] =
                bdry_topo->GetCCWBdryNode(qdata->bl,qdata->br) ;
            if (qdata->front_nodes[5] == -1)
             return 0;
            break ;

        case 2:
            if (exists[0] != exists[2]) {
                for (ii=0 ; ii<4 ; ++ii) {
                    if (exists[ii] && exists[(ii+1)%4]) {
                        Rotate(ii,id0,id1,id2,id3,&qdata->bl,
                               &qdata->br,&qdata->tr,&qdata->tl) ;
                        break ;
                    }
                }
                qdata->num_front_nodes = 5 ;
                qdata->front_nodes[0] =
                    bdry_topo->GetCWBdryNode(qdata->bl,qdata->br) ;
                if (qdata->front_nodes[0] == -1)
                 return 0;
                qdata->front_nodes[1] = qdata->bl ;
                qdata->front_nodes[2] = qdata->tl ;
                qdata->front_nodes[3] = qdata->tr ;
                qdata->front_nodes[4] =
                    bdry_topo->GetCCWBdryNode(qdata->br,qdata->tr) ;
                if (qdata->front_nodes[4] == -1)
                 return 0;
            } else {
                qdata->qcase = 5 ;
                for (ii=0 ; ii<4 ; ++ii) {
                    if (exists[ii] && exists[(ii+2)%4]) {
                        Rotate(ii,id0,id1,id2,id3,&qdata->bl,
                               &qdata->br,&qdata->tr,&qdata->tl) ;
                        break ;
                    }
                }
                qdata->num_front_nodes = 8 ;
                qdata->front_nodes[0] =
                    bdry_topo->GetCWBdryNode(qdata->bl,qdata->br) ;
                if (qdata->front_nodes[0] == -1)
                 return 0;
                qdata->front_nodes[1] = qdata->bl ;
                qdata->front_nodes[2] = qdata->tl ;
                qdata->front_nodes[3] =
                    bdry_topo->GetCCWBdryNode(qdata->tr,qdata->tl) ;
                if (qdata->front_nodes[3] == -1)
                 return 0;

                qdata->front_nodes[4] =
                    bdry_topo->GetCWBdryNode(qdata->tr,qdata->tl) ;
                if (qdata->front_nodes[4] == -1)
                 return 0;
                qdata->front_nodes[5] = qdata->tr ;
                qdata->front_nodes[6] = qdata->br ;
                qdata->front_nodes[7] =
                    bdry_topo->GetCCWBdryNode(qdata->bl,qdata->br) ;
                if (qdata->front_nodes[7] == -1)
                 return 0;

            }
            break ;

        case 3:
            for (ii=0 ; ii<4 ; ++ii) {
                if (exists[ii] && exists[(ii+1)%4] && exists[(ii+2)%4]) {
                    Rotate((ii+1)%4,id0,id1,id2,id3,&qdata->bl,
                           &qdata->br,&qdata->tr,&qdata->tl) ;
                    break ;
                }
            }
            qdata->num_front_nodes = 4 ;
            qdata->front_nodes[0] =
                bdry_topo->GetCWBdryNode(qdata->tl,qdata->bl) ;
            if (qdata->front_nodes[0] == -1)
             return 0;

            qdata->front_nodes[1] = qdata->tl ;
            qdata->front_nodes[2] = qdata->tr ;
            qdata->front_nodes[3] =
                bdry_topo->GetCCWBdryNode(qdata->br,qdata->tr) ;
            if (qdata->front_nodes[3] == -1)
             return 0;

            break ;

        case 4:
            Rotate(0,id0,id1,id2,id3,&qdata->bl,
                   &qdata->br,&qdata->tr,&qdata->tl) ;
            qdata->num_front_nodes = 0 ;
            break ;
    }
    return 1;
}

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

static double CrossProd(const CArbCoord2D &b,const CArbCoord2D &i,
                        const CArbCoord2D &j)
{
    double cross ;

    cross = ((i[0] - b[0]) * (j[1] - b[1])) -
            ((i[1] - b[1]) * (j[0] - b[0])) ;
    return cross ;
}

static bool Cross(const CArbCoord2D &i1,const CArbCoord2D &i2,
                  const CArbCoord2D &j1,const CArbCoord2D &j2)
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

    // Now compute the cross product of line I with J1 and line I
    // with J2.  If have the same sign the lines cannot cross

    CArbCoord2D delt = i2 - i1 ;
    double tol = -1e-10 * delt.Magnitude() ;

    if ((CrossProd(i1,i2,j1) * CrossProd(i1,i2,j2)) >= tol)
        return false ;

    // find the cross product of line J with I1 and line J with I2.
    // If they have the same sign, the lines cannot cross

    if ((CrossProd(j1,j2,i1) * CrossProd(j1,j2,i2)) >= tol)
        return false ;

    return true ;
}

bool CArbMshEdgeList::DoesThisCrossBoundary(
                           const CArbCoord2D &pt0,
                           const CArbCoord2D &pt1) const
{
    // loop through all the edges in the front

    bool crossed = false ;

    CArbConstHashTableIterator<ArbEdgeKey,EdgeHandle> iter(edge_table) ;
    for (iter.First() ; iter.More() ; ++iter) {
        ArbEdgeKey edge = iter.Key() ;
        ArbIntNode *nd0 = node_table->Fetch(edge.id_0) ;
        ArbIntNode *nd1 = node_table->Fetch(edge.id_1) ;
        if ((nd0->coord != pt0) && (nd1->coord != pt1)) {
            if (Cross(nd0->coord,nd1->coord,pt0,pt1)) {
                crossed = true ;
                break ;
            }
        }
    }
    return(crossed) ;
}


// %(CArbMshEdgeList::UpdateQuad-void-|-ArbNewQuadData-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateQuad - update the list to add a quad element
**
**      void UpdateQuad(ArbNewQuadData *qdata)
**
**        qdata - (in)  quad description
**
**      Description: This function is called to update the edge list to
**          add a new quadrilateral to the mesh.
**
**
** -- */

void CArbMshEdgeList::UpdateQuad(ArbNewQuadData *qdata)
{

    switch (qdata->qcase) {
        case 1: qdata->edge_level = UpdateCase1(qdata) ; break ;
        case 2: qdata->edge_level = UpdateCase2(qdata) ; break ;
        case 3: qdata->edge_level = UpdateCase3(qdata) ; break ;
        case 4: qdata->edge_level = UpdateCase4(qdata) ; break ;
        case 5: qdata->edge_level = UpdateCase5(qdata) ; break ;
        case 6: qdata->edge_level = UpdateCase6(qdata) ; break ;
    }
}




// %(CArbMshEdgeList::FindAdjacent-void-|-int-const|-int-const|-int-|*-int-|*)
/* ++ ----------------------------------------------------------
**
**    FindAdjacent - find previous and next adjacent verticies
**
**      void FindAdjacent(
**              const int vtx,
**              const int new_vtx,
**              int       *prev_vtx,
**              int       *next_vtx)
**
**        vtx      - (in)  existing vertex
**        new_vtx  - (in)  new vertex
**        prev_vtx - (out) previous existing vertex
**        next_vtx - (out) next existing vertex
**
**      Description: This function looks around a vertex and determines
**          which currently adjacent verticies would come before and
**          after a new vertex.
**
**
** -- */

int CArbMshEdgeList::FindAdjacent(const int vtx,
                                   const int new_vtx,
                                   int *prev_vtx,
                                   int *next_vtx)
{
    ArbIntNode *node  = node_table->Fetch(vtx) ;
    ArbIntNode *nnode = node_table->Fetch(new_vtx) ;

    CArbTopoAdjVtxIterator iter(bdry_topo,vtx) ;
    int tprev_vtx ;

    ArbIntNode *prev = node_table->Fetch(iter.AdjVtx()) ;
    if (prev == 0)
     return 0;

    tprev_vtx = iter.AdjVtx() ;

    while (iter.More()) {
        ++iter ;
        ArbIntNode *next = node_table->Fetch(iter.AdjVtx()) ;

        if (Angle2Pi(node->coord,prev->coord,nnode->coord) <
            Angle2Pi(node->coord,prev->coord,next->coord)) {
            *prev_vtx = tprev_vtx ;
            *next_vtx = iter.AdjVtx() ;
            return 1;
        }

        ++iter ;
        prev = node_table->Fetch(iter.AdjVtx()) ;
        if (prev == 0)
         return 0;
        tprev_vtx = iter.AdjVtx() ;
    }
    return 1;
}




// %(CArbMshEdgeList::UpdateCase1-int-|-ArbNewQuadData-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateCase1 - update the edge list for a case 1 quad
**
**      int UpdateCase1(ArbNewQuadData *qdata)
**
**        qdata - (in)  quadrilateral data
**
**      Description: This function updates the edge list to add a case
**          1 quadrilateral. A case 1 quad has one edge on the front
**          and adds three new edges.
**
**      Return Value: The edge level of the new edges.
**
**
** -- */

int CArbMshEdgeList::UpdateCase1(ArbNewQuadData *qdata)
{
    int level ;
    double length ;

    // determine the level to assign to the new edges

    EdgeHandle *base =
        edge_table->Fetch(ArbEdgeKey(qdata->bl,qdata->br)) ;

    if (base == 0)
     return -1;

    ArbIntQdEdgeDesc *desc_ptr = order_heap->ViewWithHandle(*base) ;

    level = desc_ptr->level + 1 ;
    length = desc_ptr->length ;

    // delete the base information from the tables

    order_heap->RemoveWithHandle(*base) ;
    edge_lengths->Remove(length) ;
    edge_table->Remove(ArbEdgeKey(qdata->bl,qdata->br)) ;

    // create the new edges

    edge_table->Store(ArbEdgeKey(qdata->bl,qdata->tl),0) ;
    edge_table->Store(ArbEdgeKey(qdata->tl,qdata->tr),0) ;
    edge_table->Store(ArbEdgeKey(qdata->tr,qdata->br),0) ;

    front_nodes->Store(qdata->bl,qdata->tl) ;
    front_nodes->Store(qdata->tl,qdata->tr) ;
    front_nodes->Store(qdata->tr,qdata->br) ;

    // update the edges before and after the edge

    bdry_topo->DeleteAngle(qdata->bl,qdata->br,qdata->front_nodes[0]) ;
    bdry_topo->InsertAngle(1,qdata->bl,qdata->tl,qdata->front_nodes[0]) ;

    if (bdry_topo->HasVtx(qdata->tl)) {
        int prev_vtx, next_vtx ;
        if (!FindAdjacent(qdata->tl,qdata->bl,&prev_vtx,&next_vtx))
         return -1;
        bdry_topo->DeleteAngle(qdata->tl,prev_vtx,next_vtx) ;
        bdry_topo->InsertAngle(1,qdata->tl,prev_vtx,qdata->bl) ;
        bdry_topo->InsertAngle(1,qdata->tl,qdata->tr,next_vtx) ;
    } else {
        bdry_topo->InsertAngle(1,qdata->tl,qdata->tr,qdata->bl) ;
    }

    if (bdry_topo->HasVtx(qdata->tr)) {
        int prev_vtx, next_vtx ;
        if (!FindAdjacent(qdata->tr,qdata->tl,&prev_vtx,&next_vtx))
         return -1;
        bdry_topo->DeleteAngle(qdata->tr,prev_vtx,next_vtx) ;
        bdry_topo->InsertAngle(1,qdata->tr,prev_vtx,qdata->tl) ;
        bdry_topo->InsertAngle(1,qdata->tr,qdata->br,next_vtx) ;
    } else {
        bdry_topo->InsertAngle(1,qdata->tr,qdata->br,qdata->tl) ;
    }

    bdry_topo->DeleteAngle(qdata->br,qdata->front_nodes[5],qdata->bl) ;
    bdry_topo->InsertAngle(1,qdata->br,qdata->front_nodes[5],qdata->tr) ;

    return(level) ;
}




// %(CArbMshEdgeList::UpdateCase2-int-|-ArbNewQuadData-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateCase2 - update the edge list for a case 2 quad
**
**      int UpdateCase2(ArbNewQuadData *qdata)
**
**        qdata - (in)  quadrilateral data
**
**      Description: This function updates the edge list to add a case
**          2 quadrilateral. A case 2 quad has two edge on the front
**          and adds two new edges. The two edges on the front are
**          adjacent.
**
**      Return Value: The edge level of the new edges.
**
**
** -- */

int CArbMshEdgeList::UpdateCase2(ArbNewQuadData *qdata)
{
    int level ;

    // determine the level to assign to the new edges

    EdgeHandle *base0 =
        edge_table->Fetch(ArbEdgeKey(qdata->bl,qdata->br)) ;
    EdgeHandle *base1 =
        edge_table->Fetch(ArbEdgeKey(qdata->br,qdata->tr)) ;

    if (base0 == 0)
     return -1;
    if (base1 == 0)
     return -1;

    ArbIntQdEdgeDesc *desc_ptr0 = order_heap->ViewWithHandle(*base0) ;
    ArbIntQdEdgeDesc *desc_ptr1 = order_heap->ViewWithHandle(*base1) ;

    level = ((desc_ptr0->level < desc_ptr1->level) ?
              desc_ptr0->level : desc_ptr1->level) + 1 ;

    double length0 = desc_ptr0->length ;
    double length1 = desc_ptr1->length ;

    // delete the base information from the tables

    order_heap->RemoveWithHandle(*base0) ;
    order_heap->RemoveWithHandle(*base1) ;
    edge_lengths->Remove(length0) ;
    edge_lengths->Remove(length1) ;
    edge_table->Remove(ArbEdgeKey(qdata->bl,qdata->br)) ;
    edge_table->Remove(ArbEdgeKey(qdata->br,qdata->tr)) ;
    front_nodes->Remove(qdata->br) ;

    bdry_topo->DeleteAngle(qdata->br,qdata->tr,qdata->bl) ;

    edge_table->Store(ArbEdgeKey(qdata->bl,qdata->tl),0) ;
    edge_table->Store(ArbEdgeKey(qdata->tl,qdata->tr),0) ;

    front_nodes->Store(qdata->bl,qdata->tl) ;
    front_nodes->Store(qdata->tl,qdata->tr) ;

    // update the edges before and after the edge

    bdry_topo->DeleteAngle(qdata->bl,qdata->br,qdata->front_nodes[0]) ;
    bdry_topo->InsertAngle(1,qdata->bl,qdata->tl,qdata->front_nodes[0]) ;

    if (bdry_topo->HasVtx(qdata->tl)) {
        int prev_vtx, next_vtx ;
        if (!FindAdjacent(qdata->tl,qdata->bl,&prev_vtx,&next_vtx))
         return -1;
        bdry_topo->DeleteAngle(qdata->tl,prev_vtx,next_vtx) ;
        bdry_topo->InsertAngle(1,qdata->tl,prev_vtx,qdata->bl) ;
        bdry_topo->InsertAngle(1,qdata->tl,qdata->tr,next_vtx) ;
    } else {
        bdry_topo->InsertAngle(1,qdata->tl,qdata->tr,qdata->bl) ;
    }

    bdry_topo->DeleteAngle(qdata->tr,qdata->front_nodes[4],qdata->br) ;
    bdry_topo->InsertAngle(1,qdata->tr,qdata->front_nodes[4],qdata->tl) ;

    return(level) ;
}




// %(CArbMshEdgeList::UpdateCase3-int-|-ArbNewQuadData-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateCase3 - update the edge list for a case 3 quad
**
**      int UpdateCase3(ArbNewQuadData *qdata)
**
**        qdata - (in)  quadrilateral data
**
**      Description: This function updates the edge list to add a case
**          3 quadrilateral. A case 3 quad has three edge on the front
**          and adds one new edges.
**
**      Return Value: The edge level of the new edges.
**
**
** -- */

int CArbMshEdgeList::UpdateCase3(ArbNewQuadData *qdata)
{
    int level ;

    // determine the level to assign to the new edges

    EdgeHandle *base0 =
        edge_table->Fetch(ArbEdgeKey(qdata->tl,qdata->bl)) ;
    EdgeHandle *base1 =
        edge_table->Fetch(ArbEdgeKey(qdata->bl,qdata->br)) ;
    EdgeHandle *base2 =
        edge_table->Fetch(ArbEdgeKey(qdata->br,qdata->tr)) ;

    ArbIntQdEdgeDesc *desc_ptr0 = order_heap->ViewWithHandle(*base0) ;
    ArbIntQdEdgeDesc *desc_ptr1 = order_heap->ViewWithHandle(*base1) ;
    ArbIntQdEdgeDesc *desc_ptr2 = order_heap->ViewWithHandle(*base2) ;

    double length0 = desc_ptr0->length ;
    double length1 = desc_ptr1->length ;
    double length2 = desc_ptr2->length ;

    if (desc_ptr0->level == desc_ptr2->level)
        level = desc_ptr0->level ;
    else
        level = (desc_ptr0->level > desc_ptr2->level) ?
                 desc_ptr0->level : desc_ptr2->level ;

    // delete the base information from the tables

    order_heap->RemoveWithHandle(*base0) ;
    order_heap->RemoveWithHandle(*base1) ;
    order_heap->RemoveWithHandle(*base2) ;
    edge_lengths->Remove(length0) ;
    edge_lengths->Remove(length1) ;
    edge_lengths->Remove(length2) ;
    edge_table->Remove(ArbEdgeKey(qdata->tl,qdata->bl)) ;
    edge_table->Remove(ArbEdgeKey(qdata->bl,qdata->br)) ;
    edge_table->Remove(ArbEdgeKey(qdata->br,qdata->tr)) ;
    front_nodes->Remove(qdata->bl) ;
    front_nodes->Remove(qdata->br) ;

    bdry_topo->DeleteAngle(qdata->bl,qdata->br,qdata->tl) ;
    bdry_topo->DeleteAngle(qdata->br,qdata->tr,qdata->bl) ;

    // create the new top edge

    edge_table->Store(ArbEdgeKey(qdata->tl,qdata->tr),0) ;

    front_nodes->Store(qdata->tl,qdata->tr) ;

    // update the edges before and after the edge

    bdry_topo->DeleteAngle(qdata->tl,qdata->bl,qdata->front_nodes[0]) ;
    bdry_topo->InsertAngle(1,qdata->tl,qdata->tr,qdata->front_nodes[0]) ;

    bdry_topo->DeleteAngle(qdata->tr,qdata->front_nodes[3],qdata->br) ;
    bdry_topo->InsertAngle(1,qdata->tr,qdata->front_nodes[3],qdata->tl) ;

    return(level) ;
}




// %(CArbMshEdgeList::UpdateCase4-int-|-ArbNewQuadData-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateCase4 - update the edge list for a case 4 quad
**
**      int UpdateCase4(ArbNewQuadData *qdata)
**
**        qdata - (in)  quadrilateral data
**
**      Description: This function updates the edge list to add a case
**          4 quadrilateral. A case 4 quad has all for edge on the
**          front and adds no new edges.
**
**      Return Value: The edge level of the new edges.
**
**
** -- */

int CArbMshEdgeList::UpdateCase4(ArbNewQuadData *qdata)
{
    // determine the level to assign to the new edges

    EdgeHandle *base0 =
        edge_table->Fetch(ArbEdgeKey(qdata->tl,qdata->bl)) ;
    if (base0 == 0)
     return -1;
    EdgeHandle *base1 =
        edge_table->Fetch(ArbEdgeKey(qdata->bl,qdata->br)) ;
    if (base1 == 0)
     return -1;
    EdgeHandle *base2 =
        edge_table->Fetch(ArbEdgeKey(qdata->br,qdata->tr)) ;
    if (base2 == 0)
     return -1;
    EdgeHandle *base3 =
        edge_table->Fetch(ArbEdgeKey(qdata->tr,qdata->tl)) ;
    if (base3 == 0)
     return -1;

    ArbIntQdEdgeDesc *desc_ptr0 = order_heap->ViewWithHandle(*base0) ;
    ArbIntQdEdgeDesc *desc_ptr1 = order_heap->ViewWithHandle(*base1) ;
    ArbIntQdEdgeDesc *desc_ptr2 = order_heap->ViewWithHandle(*base2) ;
    ArbIntQdEdgeDesc *desc_ptr3 = order_heap->ViewWithHandle(*base3) ;

    double length0 = desc_ptr0->length ;
    double length1 = desc_ptr1->length ;
    double length2 = desc_ptr2->length ;
    double length3 = desc_ptr3->length ;

    // delete the base information from the tables

    order_heap->RemoveWithHandle(*base0) ;
    order_heap->RemoveWithHandle(*base1) ;
    order_heap->RemoveWithHandle(*base2) ;
    order_heap->RemoveWithHandle(*base3) ;
    edge_lengths->Remove(length0) ;
    edge_lengths->Remove(length1) ;
    edge_lengths->Remove(length2) ;
    edge_lengths->Remove(length3) ;
    edge_table->Remove(ArbEdgeKey(qdata->tl,qdata->bl)) ;
    edge_table->Remove(ArbEdgeKey(qdata->bl,qdata->br)) ;
    edge_table->Remove(ArbEdgeKey(qdata->br,qdata->tr)) ;
    edge_table->Remove(ArbEdgeKey(qdata->tr,qdata->tl)) ;
    front_nodes->Remove(qdata->bl) ;
    front_nodes->Remove(qdata->br) ;
    front_nodes->Remove(qdata->tl) ;
    front_nodes->Remove(qdata->tr) ;

    bdry_topo->DeleteAngle(qdata->bl,qdata->br,qdata->tl) ;
    bdry_topo->DeleteAngle(qdata->br,qdata->tr,qdata->bl) ;
    bdry_topo->DeleteAngle(qdata->tr,qdata->tl,qdata->br) ;
    bdry_topo->DeleteAngle(qdata->tl,qdata->bl,qdata->tr) ;

    return(0) ;
}




// %(CArbMshEdgeList::UpdateCase5-int-|-ArbNewQuadData-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateCase5 - update the edge list for a case 4 quad
**
**      int UpdateCase5(ArbNewQuadData *qdata)
**
**        qdata - (in)  quadrilateral data
**
**      Description: This function updates the edge list to add a case
**          2 quadrilateral. A case 2 quad has two edge on the front
**          and adds two new edges. The two edges on the front are not
**          adjacent.
**
**      Return Value: The edge level of the new edges.
**
**
** -- */

int CArbMshEdgeList::UpdateCase5(ArbNewQuadData *qdata)
{
    int level ;

    // determine the level to assign to the new edges

    EdgeHandle *base0 =
        edge_table->Fetch(ArbEdgeKey(qdata->bl,qdata->br)) ;
    EdgeHandle *base1 =
        edge_table->Fetch(ArbEdgeKey(qdata->tr,qdata->tl)) ;

    ArbIntQdEdgeDesc *desc_ptr0 = order_heap->ViewWithHandle(*base0) ;
    ArbIntQdEdgeDesc *desc_ptr1 = order_heap->ViewWithHandle(*base1) ;

    level = ((desc_ptr0->level < desc_ptr1->level) ?
              desc_ptr0->level : desc_ptr1->level) + 1 ;

    double length0 = desc_ptr0->length ;
    double length1 = desc_ptr1->length ;

    // delete the base information from the tables

    order_heap->RemoveWithHandle(*base0) ;
    order_heap->RemoveWithHandle(*base1) ;
    edge_lengths->Remove(length0) ;
    edge_lengths->Remove(length1) ;
    edge_table->Remove(ArbEdgeKey(qdata->bl,qdata->br)) ;
    edge_table->Remove(ArbEdgeKey(qdata->tr,qdata->tl)) ;

    edge_table->Store(ArbEdgeKey(qdata->bl,qdata->tl),0) ;
    edge_table->Store(ArbEdgeKey(qdata->tr,qdata->br),0) ;

    front_nodes->Store(qdata->bl,qdata->tl) ;
    front_nodes->Store(qdata->tr,qdata->br) ;

    // update the edges before and after the edge

    bdry_topo->DeleteAngle(qdata->bl,qdata->br,qdata->front_nodes[0]) ;
    bdry_topo->InsertAngle(1,qdata->bl,qdata->tl,qdata->front_nodes[0]) ;

    bdry_topo->DeleteAngle(qdata->tl,qdata->front_nodes[3],qdata->tr) ;
    bdry_topo->InsertAngle(1,qdata->tl,qdata->front_nodes[3],qdata->bl) ;

    bdry_topo->DeleteAngle(qdata->tr,qdata->tl,qdata->front_nodes[4]) ;
    bdry_topo->InsertAngle(1,qdata->tr,qdata->br,qdata->front_nodes[4]) ;

    bdry_topo->DeleteAngle(qdata->br,qdata->front_nodes[7],qdata->bl) ;
    bdry_topo->InsertAngle(1,qdata->br,qdata->front_nodes[7],qdata->tr) ;

    return(level) ;
}


// %(CArbMshEdgeList::UpdateCase6-int-|-ArbNewQuadData-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateCase6 - update the edge list for a case 6 quad
**                   (triangle)
**
**      int UpdateCase6(ArbNewQuadData *qdata)
**
**        qdata - (in)  quadrilateral data
**
**      Description: This function updates the edge list to add a case
**          6 quadrilateral. A case 6 quad is realy an triangular
**          update.  Unlike the UpdateTriangle method, this method
**          updates the edge list properly
**
**      Return Value: The edge level of the new edges.
**
**
** -- */

int CArbMshEdgeList::UpdateCase6(ArbNewQuadData *qdata)
{
    int level ;

    /*  This routine assumes the following configuration


      font[0]   bl     tl = tr   front[3]
        o--------*       *---------o
                  \     /
                   \   /
                    \ /
                     o br

    */

    // determine the level to assign to the new edges

    EdgeHandle *base0 =
        edge_table->Fetch(ArbEdgeKey(qdata->bl,qdata->br)) ;
    EdgeHandle *base1 =
        edge_table->Fetch(ArbEdgeKey(qdata->br,qdata->tr)) ;

    ArbIntQdEdgeDesc *desc_ptr0 = order_heap->ViewWithHandle(*base0) ;
    ArbIntQdEdgeDesc *desc_ptr1 = order_heap->ViewWithHandle(*base1) ;

    double length0 = desc_ptr0->length ;
    double length1 = desc_ptr1->length ;

    level = (desc_ptr0->level > desc_ptr1->level) ?
             desc_ptr0->level : desc_ptr1->level ;

    // delete the base information from the tables

    order_heap->RemoveWithHandle(*base0) ;
    order_heap->RemoveWithHandle(*base1) ;
    edge_lengths->Remove(length0) ;
    edge_lengths->Remove(length1) ;
    edge_table->Remove(ArbEdgeKey(qdata->bl,qdata->br)) ;
    edge_table->Remove(ArbEdgeKey(qdata->br,qdata->tr)) ;
    front_nodes->Remove(qdata->br) ;

    bdry_topo->DeleteAngle(qdata->br,qdata->tr,qdata->bl) ;

    // create the new top edge

    edge_table->Store(ArbEdgeKey(qdata->bl,qdata->tr),0) ;

    front_nodes->Store(qdata->bl,qdata->tr) ;

    // update the edges before and after the edge

    bdry_topo->DeleteAngle(qdata->bl,qdata->br,qdata->front_nodes[0]) ;
    bdry_topo->InsertAngle(1,qdata->bl,qdata->tr,qdata->front_nodes[0]) ;

    bdry_topo->DeleteAngle(qdata->tr,qdata->front_nodes[3],qdata->br) ;
    bdry_topo->InsertAngle(1,qdata->tr,qdata->front_nodes[3],qdata->bl) ;

    return(level) ;
}



// %(CArbMshEdgeList::UpdateTriangle-void-|-int-const|-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    UpdateTriangle - update the list to add a quad element
**
**      void UpdateTriangle(
**              const int id_0,
**              const int id_1,
**              const int id_2)
**
**        id_0 - (in)  first id
**        id_1 - (in)  second id
**        id_2 - (in)  third id
**
**      Description: This function is called to update the edge list to
**          add a new triangular element to the mesh. This is done only
**          to close off a triangular void.
**
**
** -- */

void CArbMshEdgeList::UpdateTriangle(const int id0,
                                     const int id1,
                                     const int id2)
{
    // delete the base information from the tables

    EdgeHandle *base0 = edge_table->Fetch(ArbEdgeKey(id0,id1)) ;
    EdgeHandle *base1 = edge_table->Fetch(ArbEdgeKey(id1,id2)) ;
    EdgeHandle *base2 = edge_table->Fetch(ArbEdgeKey(id2,id0)) ;

    ArbIntQdEdgeDesc *desc_ptr0 = order_heap->ViewWithHandle(*base0) ;
    ArbIntQdEdgeDesc *desc_ptr1 = order_heap->ViewWithHandle(*base1) ;
    ArbIntQdEdgeDesc *desc_ptr2 = order_heap->ViewWithHandle(*base2) ;

    double length0 = desc_ptr0->length ;
    double length1 = desc_ptr1->length ;
    double length2 = desc_ptr2->length ;

    order_heap->RemoveWithHandle(*base0) ;
    order_heap->RemoveWithHandle(*base1) ;
    order_heap->RemoveWithHandle(*base2) ;

    edge_lengths->Remove(length0) ;
    edge_lengths->Remove(length1) ;
    edge_lengths->Remove(length2) ;

    edge_table->Remove(ArbEdgeKey(id0,id1)) ;
    edge_table->Remove(ArbEdgeKey(id1,id2)) ;
    edge_table->Remove(ArbEdgeKey(id2,id0)) ;

    front_nodes->Remove(id0) ;
    front_nodes->Remove(id1) ;
    front_nodes->Remove(id2) ;

    bdry_topo->DeleteTriangle(id0,id1,id2) ;
}




// %(CArbMshEdgeList::UpdateSeam-int-|-ArbNewSeamData-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateSeam - update the front for a seam operation
**
**      int UpdateSeam(ArbNewSeamData *sdata)
**
**        sdata - (in)  seam description data
**
**      Description: This function updates the edge list for a seam
**          operation.
**
**      Return Value: The function returns the id of the node that will
**          be kept after the seam operation.
**
**
** -- */

int CArbMshEdgeList::UpdateSeam(ArbNewSeamData *sdata)
{
    /*
        Topologically the seaming operation looks like
        the following:

           cw *
              |                    cw *
              |                       |
          id0 *               =>      |
              |                   id0 *-----*
              |                            ccw
          id1 *-----*-----*
                   id2   ccw
    */

    // first determine the coordinates of the joined node
    // location.  This will be the average of the location
    // of nodes 0 and 2 unless one of the nodes is fixed


    ArbIntNode *node_0 = node_table->Fetch(sdata->id0) ;
    ArbIntNode *node_2 = node_table->Fetch(sdata->id2) ;
    int keep ;
    if (node_0->motion == ARB_FIXED) {
        keep = sdata->id0 ;
    } else if (node_2->motion == ARB_FIXED) {
        keep = sdata->id2 ;
    } else {
        keep = sdata->id0 ;
    }

    // now delete edges 0-1 and 1-2 from the edge list

    EdgeHandle *handle ;
    ArbIntQdEdgeDesc *desc_ptr ;

    handle = edge_table->Fetch(ArbEdgeKey(sdata->id0,sdata->id1)) ;
    desc_ptr = order_heap->ViewWithHandle(*handle) ;
    edge_lengths->Remove(desc_ptr->length) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_table->Remove(ArbEdgeKey(sdata->id0,sdata->id1)) ;

    handle = edge_table->Fetch(ArbEdgeKey(sdata->id1,sdata->id2)) ;
    desc_ptr = order_heap->ViewWithHandle(*handle) ;
    edge_lengths->Remove(desc_ptr->length) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_table->Remove(ArbEdgeKey(sdata->id1,sdata->id2)) ;

    bdry_topo->DeleteAngle(sdata->id0,sdata->id1,sdata->cw) ;
    bdry_topo->DeleteAngle(sdata->id1,sdata->id2,sdata->id0) ;
    bdry_topo->DeleteAngle(sdata->id2,sdata->ccw,sdata->id1) ;
    front_nodes->Remove(sdata->id1) ;

    if (sdata->cw == sdata->ccw) {

        // deal with the special case where cw == cww.  If this
        // is the case we have a "closing" seam, so we delete
        // all four edges.

        handle = edge_table->Fetch(ArbEdgeKey(sdata->cw,sdata->id0)) ;
        desc_ptr = order_heap->ViewWithHandle(*handle) ;
        edge_lengths->Remove(desc_ptr->length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Remove(ArbEdgeKey(sdata->cw,sdata->id0)) ;

        handle = edge_table->Fetch(ArbEdgeKey(sdata->id2,sdata->ccw)) ;
        desc_ptr = order_heap->ViewWithHandle(*handle) ;
        edge_lengths->Remove(desc_ptr->length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Remove(ArbEdgeKey(sdata->id2,sdata->ccw)) ;

        bdry_topo->DeleteAngle(sdata->cw,sdata->id0,sdata->id2) ;

        front_nodes->Remove(sdata->id0) ;
        front_nodes->Remove(sdata->id2) ;
        front_nodes->Remove(sdata->cw) ;

        // check to see if there are any edges attached to
        // id2 that need to be updated

        CArbTopoAdjVtxIterator iter(bdry_topo,sdata->id2) ;

        while(iter.More()) {
            int prev = iter.AdjVtx() ;
            int prev_elem = iter.CcwElem() ;
            ++iter ;
            if ((prev_elem == 1) && (prev != sdata->ccw)) {
                int next = iter.AdjVtx() ;

                EdgeHandle *handle, value ;
                handle = edge_table->Fetch(ArbEdgeKey(sdata->id2,prev)) ;
                value = *handle ;
                desc_ptr = order_heap->ViewWithHandle(value) ;
                desc_ptr->id_0 = sdata->id0 ;

                edge_table->Remove(ArbEdgeKey(sdata->id2,prev)) ;
                edge_table->Store(ArbEdgeKey(sdata->id0,prev),value) ;

                handle = edge_table->Fetch(ArbEdgeKey(next,sdata->id2)) ;
                value = *handle ;
                desc_ptr = order_heap->ViewWithHandle(value) ;
                desc_ptr->id_1 = sdata->id0 ;

                edge_table->Remove(ArbEdgeKey(next,sdata->id2)) ;
                edge_table->Store(ArbEdgeKey(next,sdata->id0),value) ;

                int ccw =
                    bdry_topo->GetCCWBdryNode(sdata->id2,prev) ;
                int cw =
                    bdry_topo->GetCWBdryNode(next,sdata->id2) ;

                bdry_topo->DeleteAngle(sdata->id2,prev,next) ;
                bdry_topo->InsertAngle(1,sdata->id0,prev,next) ;

                bdry_topo->DeleteAngle(prev,ccw,sdata->id2) ;
                bdry_topo->InsertAngle(1,prev,ccw,sdata->id0) ;

                bdry_topo->DeleteAngle(next,sdata->id2,cw) ;
                bdry_topo->InsertAngle(1,next,sdata->id0,cw) ;

                front_nodes->Remove(sdata->id2) ;
                front_nodes->Store(sdata->id0,prev) ;
                front_nodes->Store(next,sdata->id0) ;

                iter.NewVtx(sdata->id2) ;
             }
        }
    } else {

        // update the edges

        if (keep == sdata->id0) {
            bdry_topo->InsertAngle(1,sdata->id0,sdata->ccw,sdata->cw) ;

            EdgeHandle *handle, value ;
            handle = edge_table->Fetch(ArbEdgeKey(sdata->id2,sdata->ccw)) ;
            value = *handle ;
            desc_ptr = order_heap->ViewWithHandle(value) ;
            desc_ptr->id_0 = sdata->id0 ;

            edge_table->Remove(ArbEdgeKey(sdata->id2,sdata->ccw)) ;
            edge_table->Store(ArbEdgeKey(sdata->id0,sdata->ccw),value) ;

            front_nodes->Store(sdata->id0,sdata->ccw) ;

            int nccw = bdry_topo->GetCCWBdryNode(sdata->id2,sdata->ccw) ;
            bdry_topo->DeleteAngle(sdata->ccw,nccw,sdata->id2) ;
            bdry_topo->InsertAngle(1,sdata->ccw,nccw,sdata->id0) ;

            // check to see if there are any edges attached to
            // id2 that need to be updated

            CArbTopoAdjVtxIterator iter(bdry_topo,sdata->id2) ;

            while(iter.More()) {
                int prev = iter.AdjVtx() ;
                int prev_elem = iter.CcwElem() ;
                ++iter ;
                if ((prev_elem == 1) && (prev != sdata->ccw)) {
                    int next = iter.AdjVtx() ;

                    handle = edge_table->Fetch(ArbEdgeKey(sdata->id2,prev)) ;
                    value = *handle ;
                    desc_ptr = order_heap->ViewWithHandle(value) ;
                    desc_ptr->id_0 = sdata->id0 ;

                    edge_table->Remove(ArbEdgeKey(sdata->id2,prev)) ;
                    edge_table->Store(ArbEdgeKey(sdata->id0,prev),value) ;

                    handle = edge_table->Fetch(ArbEdgeKey(next,sdata->id2)) ;
                    value = *handle ;
                    desc_ptr = order_heap->ViewWithHandle(value) ;
                    desc_ptr->id_1 = sdata->id0 ;

                    edge_table->Remove(ArbEdgeKey(next,sdata->id2)) ;
                    edge_table->Store(ArbEdgeKey(next,sdata->id0),value) ;

                    int ccw =
                        bdry_topo->GetCCWBdryNode(sdata->id2,prev) ;
                    int cw =
                        bdry_topo->GetCWBdryNode(next,sdata->id2) ;

                    bdry_topo->DeleteAngle(sdata->id2,prev,next) ;
                    bdry_topo->InsertAngle(1,sdata->id0,prev,next) ;

                    bdry_topo->DeleteAngle(prev,ccw,sdata->id2) ;
                    bdry_topo->InsertAngle(1,prev,ccw,sdata->id0) ;

                    bdry_topo->DeleteAngle(next,sdata->id2,cw) ;
                    bdry_topo->InsertAngle(1,next,sdata->id0,cw) ;

                    front_nodes->Remove(sdata->id2) ;
                    front_nodes->Store(sdata->id0,prev) ;

                    iter.NewVtx(sdata->id2) ;
                }
            }
        } else {
            bdry_topo->InsertAngle(1,sdata->id2,sdata->ccw,sdata->cw) ;

            EdgeHandle *handle, value ;
            handle = edge_table->Fetch(ArbEdgeKey(sdata->cw,sdata->id0)) ;
            value = *handle ;
            desc_ptr = order_heap->ViewWithHandle(value) ;
            desc_ptr->id_1 = sdata->id2 ;

            edge_table->Remove(ArbEdgeKey(sdata->cw,sdata->id0)) ;
            edge_table->Store(ArbEdgeKey(sdata->cw,sdata->id2),value) ;

            front_nodes->Store(sdata->cw,sdata->id2) ;

            int ncw = bdry_topo->GetCWBdryNode(sdata->cw,sdata->id0) ;
            bdry_topo->DeleteAngle(sdata->cw,sdata->id0,ncw) ;
            bdry_topo->InsertAngle(1,sdata->cw,sdata->id2,ncw) ;

            sdata->id2 = sdata->id0 ;
            sdata->id0 = keep ;

            // check to see if there are any edges attached to
            // id0 that need to be updated

            CArbTopoAdjVtxIterator iter(bdry_topo,sdata->id2) ;

            while(iter.More()) {
                int prev = iter.AdjVtx() ;
                int prev_elem = iter.CcwElem() ;
                ++iter ;
                if ((prev_elem == 1) && (prev != sdata->ccw)) {
                    int next = iter.AdjVtx() ;

                    handle = edge_table->Fetch(ArbEdgeKey(sdata->id0,prev)) ;
                    value = *handle ;
                    desc_ptr = order_heap->ViewWithHandle(value) ;
                    desc_ptr->id_0 = sdata->id2 ;

                    edge_table->Remove(ArbEdgeKey(sdata->id0,prev)) ;
                    edge_table->Store(ArbEdgeKey(sdata->id2,prev),value) ;

                    handle = edge_table->Fetch(ArbEdgeKey(next,sdata->id0)) ;
                    value = *handle ;
                    desc_ptr = order_heap->ViewWithHandle(value) ;
                    desc_ptr->id_1 = sdata->id2 ;

                    edge_table->Remove(ArbEdgeKey(next,sdata->id0)) ;
                    edge_table->Store(ArbEdgeKey(next,sdata->id2),value) ;

                    int ccw =
                        bdry_topo->GetCCWBdryNode(sdata->id0,prev) ;
                    int cw =
                        bdry_topo->GetCWBdryNode(next,sdata->id0) ;

                    bdry_topo->DeleteAngle(sdata->id0,prev,next) ;
                    bdry_topo->InsertAngle(1,sdata->id2,prev,next) ;

                    bdry_topo->DeleteAngle(prev,ccw,sdata->id0) ;
                    bdry_topo->InsertAngle(1,prev,ccw,sdata->id2) ;

                    bdry_topo->DeleteAngle(next,sdata->id0,cw) ;
                    bdry_topo->InsertAngle(1,next,sdata->id2,cw) ;

                    front_nodes->Remove(sdata->id0) ;
                    front_nodes->Store(sdata->id2,prev) ;

                    iter.NewVtx(sdata->id2) ;
                }
            }
        }
    }
    return(keep) ;
}




// %(CArbMshEdgeList::UpdateTransSeam-void-|-ArbNewSeamData-const|*-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    UpdateTransSeam - update the front for a transition seam
**
**      void UpdateTransSeam(
**              const ArbNewSeamData *sdata,
**              const int            mid_id,
**              const int            far_id)
**
**        sdata  - (in)  seam description data
**        mid_id - (in)  id of the mid point for the transition
**        far_id - (in)  id of the far point for the transition
**
**      Description: This function updates the edge list for a
**          transition seam operation.
**
**
** -- */

void CArbMshEdgeList::UpdateTransSeam(
                 const ArbNewSeamData *sdata,
                 const int mid_id,
                 const int far_id)
{
    if (sdata->ratio > 1.0) {

        // delete the existing edge

        EdgeHandle *handle, new_handle ;
        ArbIntQdEdgeDesc desc ;

        handle = edge_table->Fetch(ArbEdgeKey(sdata->id0,sdata->id1)) ;
        desc = *(order_heap->ViewWithHandle(*handle)) ;
        double edge_length = desc.length ;
        edge_lengths->Remove(edge_length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Remove(ArbEdgeKey(sdata->id0,sdata->id1)) ;

        bdry_topo->DeleteAngle(sdata->id0,sdata->id1,sdata->cw) ;
        bdry_topo->DeleteAngle(sdata->id1,sdata->id2,sdata->id0) ;

        // now add the information for the new edges,
        // we make the descriptor information only partially
        // correct because it will be updated later

        desc.length = 0.5 * edge_length ;
        desc.id_1 = mid_id ;
        desc.angle_1 = HALF_PI ;
        desc.end_code &= LEFT_EDGE_MASK ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(sdata->id0,mid_id),new_handle) ;

        desc.length = edge_length ;
        desc.id_0 = mid_id ;
        desc.id_1 = far_id ;
        desc.angle_0 = desc.angle_1 = HALF_PI ;
        desc.end_code = 0 ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(mid_id,far_id),new_handle) ;

        desc.id_0 = far_id ;
        desc.id_1 = sdata->id1 ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(far_id,sdata->id1),new_handle) ;

        front_nodes->Store(sdata->id0,mid_id) ;
        front_nodes->Store(mid_id,far_id) ;
        front_nodes->Store(far_id,sdata->id1) ;

        // update the angles

        bdry_topo->InsertAngle(1,sdata->id0,mid_id,sdata->cw) ;
        bdry_topo->InsertAngle(1,mid_id,far_id,sdata->id0) ;
        bdry_topo->InsertAngle(1,far_id,sdata->id1,mid_id) ;
        bdry_topo->InsertAngle(1,sdata->id1,sdata->id2,far_id) ;

    } else {

        // delete the existing edge

        EdgeHandle *handle, new_handle ;
        ArbIntQdEdgeDesc desc ;

        handle = edge_table->Fetch(ArbEdgeKey(sdata->id1,sdata->id2)) ;
        desc = *(order_heap->ViewWithHandle(*handle)) ;
        double edge_length = desc.length ;
        edge_lengths->Remove(edge_length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Remove(ArbEdgeKey(sdata->id1,sdata->id2)) ;

        bdry_topo->DeleteAngle(sdata->id1,sdata->id2,sdata->id0) ;
        bdry_topo->DeleteAngle(sdata->id2,sdata->ccw,sdata->id1) ;

        // now add the information for the new edges,
        // we make the descriptor information only partially
        // correct because it will be updated later

        desc.length = 0.5 * edge_length ;
        desc.id_0 = mid_id ;
        desc.angle_0 = HALF_PI ;
        desc.end_code &= RIGHT_EDGE_MASK ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(mid_id,sdata->id2),new_handle) ;

        desc.length = edge_length ;
        desc.id_0 = far_id ;
        desc.id_1 = mid_id ;
        desc.angle_0 = desc.angle_1 = HALF_PI ;
        desc.end_code = 0 ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(far_id,mid_id),new_handle) ;

        desc.id_0 = sdata->id1 ;
        desc.id_1 = far_id ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(sdata->id1,far_id),new_handle) ;

        front_nodes->Store(mid_id,sdata->id2) ;
        front_nodes->Store(far_id,mid_id) ;
        front_nodes->Store(sdata->id1,far_id) ;

        // update the angles

        bdry_topo->InsertAngle(1,sdata->id1,far_id,sdata->id0) ;
        bdry_topo->InsertAngle(1,far_id,mid_id,sdata->id1) ;
        bdry_topo->InsertAngle(1,mid_id,sdata->id2,far_id) ;
        bdry_topo->InsertAngle(1,sdata->id2,sdata->ccw,mid_id) ;
    }
}




// %(CArbMshEdgeList::UpdateTransSplit-void-|-ArbNewQuadData-const|*-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    UpdateTransSplit - update the front for a transition split
**
**      void UpdateTransSplit(
**              const ArbNewQuadData *qdata,
**              const int            mid_id0,
**              const int            mid_id1)
**
**        qdata   - (in)  quadrilateral description
**        mid_id0 - (in)  id of the first mid point vertex
**        mid_id1 - (in)  id of the second mid point vertex
**
**      Description: This function updates the edge list for a
**          transition split operation.
**
**
** -- */

void CArbMshEdgeList::UpdateTransSplit(
                 const ArbNewQuadData *qdata,
                 const int mid_id0,
                 const int mid_id1)
{
    if (qdata->ratio > 1.0) {

        // delete the existing edge

        EdgeHandle *handle, new_handle ;
        ArbIntQdEdgeDesc desc ;

        handle = edge_table->Fetch(ArbEdgeKey(qdata->bl,qdata->br)) ;
        desc = *(order_heap->ViewWithHandle(*handle)) ;
        double edge_length = desc.length ;
        edge_lengths->Remove(edge_length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Remove(ArbEdgeKey(qdata->bl,qdata->br)) ;

        bdry_topo->DeleteAngle(qdata->bl,qdata->br,qdata->front_nodes[0]) ;
        bdry_topo->DeleteAngle(qdata->br,qdata->tr,qdata->bl) ;

        // now add the information for the new edges,
        // we make the descriptor information only partially
        // correct because it will be updated later

        desc.length = 0.5 * edge_length ;
        desc.id_1 = mid_id0 ;
        desc.angle_1 = HALF_PI ;
        desc.end_code &= LEFT_EDGE_MASK ;
        desc.visited = true ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(qdata->bl,mid_id0),new_handle) ;

        desc.id_0 = mid_id0 ;
        desc.id_1 = mid_id1 ;
        desc.angle_0 = desc.angle_1 = HALF_PI ;
        desc.end_code = 0 ;
        desc.visited = true ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(mid_id0,mid_id1),new_handle) ;

        desc.id_0 = mid_id1 ;
        desc.id_1 = qdata->br ;
        desc.visited = true ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(mid_id1,qdata->br),new_handle) ;

        front_nodes->Store(qdata->bl,mid_id0) ;
        front_nodes->Store(mid_id0,mid_id1) ;
        front_nodes->Store(mid_id1,qdata->br) ;

        // update the angles

        bdry_topo->InsertAngle(1,qdata->bl,mid_id0,qdata->front_nodes[0]) ;
        bdry_topo->InsertAngle(1,mid_id0,mid_id1,qdata->bl) ;
        bdry_topo->InsertAngle(1,mid_id1,qdata->br,mid_id0) ;
        bdry_topo->InsertAngle(1,qdata->br,qdata->tr,mid_id1) ;

    } else {

        // delete the existing edge

        EdgeHandle *handle, new_handle ;
        ArbIntQdEdgeDesc desc ;

        handle = edge_table->Fetch(ArbEdgeKey(qdata->br,qdata->tr)) ;
        desc = *(order_heap->ViewWithHandle(*handle)) ;
        double edge_length = desc.length ;
        edge_lengths->Remove(edge_length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Remove(ArbEdgeKey(qdata->tr,qdata->tr)) ;

        bdry_topo->DeleteAngle(qdata->br,qdata->tr,qdata->bl) ;
        bdry_topo->DeleteAngle(qdata->tr,qdata->front_nodes[4],qdata->br) ;

        // now add the information for the new edges,
        // we make the descriptor information only partially
        // correct because it will be updated later

        desc.length = 0.5 * edge_length ;
        desc.id_0 = mid_id0 ;
        desc.angle_0 = HALF_PI ;
        desc.end_code &= RIGHT_EDGE_MASK ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(mid_id0,qdata->tr),new_handle) ;

        desc.id_0 = mid_id1 ;
        desc.id_1 = mid_id0 ;
        desc.angle_0 = desc.angle_1 = HALF_PI ;
        desc.end_code = 0 ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(mid_id1,mid_id0),new_handle) ;

        desc.id_0 = qdata->br ;
        desc.id_1 = mid_id1 ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(qdata->br,mid_id1),new_handle) ;

        front_nodes->Store(mid_id0,qdata->tr) ;
        front_nodes->Store(mid_id1,mid_id0) ;
        front_nodes->Store(qdata->br,mid_id1) ;

        // update the angles

        bdry_topo->InsertAngle(1,qdata->br,mid_id1,qdata->bl) ;
        bdry_topo->InsertAngle(1,mid_id1,mid_id0,qdata->br) ;
        bdry_topo->InsertAngle(1,mid_id0,qdata->tr,mid_id1) ;
        bdry_topo->InsertAngle(1,qdata->tr,qdata->front_nodes[4],mid_id0) ;
    }
}




// %(CArbMshEdgeList::UpdateTempSplit-void-|-ArbNewQuadData-const|*-int-const|-int-const|-bool-const|)
/* ++ ----------------------------------------------------------
**
**    UpdateTempSplit - update the front for a template split
**
**      void UpdateTempSplit(
**              const ArbNewQuadData *qdata,
**              const int            mid_id0,
**              const int            mid_id1,
**              const bool           base)
**
**        qdata   - (in)  quadrilateral data
**        mid_id0 - (in)  first mid point id
**        mid_id1 - (in)  second mid point id
**        base    - (in)  true means split the base edge
**
**      Description: This function updates the edge list for a template
**          split operation.
**
**
** -- */

void CArbMshEdgeList::UpdateTempSplit(
                 const ArbNewQuadData *qdata,
                 const int mid_id0,
                 const int mid_id1,
                 bool base)
{
    // update the angles

    bdry_topo->DeleteAngle(qdata->tl,qdata->tr,qdata->front_nodes[0]) ;
    bdry_topo->DeleteAngle(qdata->tr,qdata->front_nodes[3],qdata->tl) ;

    front_nodes->Store(mid_id0,qdata->tr) ;
    front_nodes->Store(mid_id1,mid_id0) ;

    // update the angles

    bdry_topo->InsertAngle(1,qdata->tl,mid_id1,qdata->front_nodes[0]) ;
    bdry_topo->InsertAngle(1,mid_id1,mid_id0,qdata->tl) ;
    bdry_topo->InsertAngle(1,mid_id0,qdata->tr,mid_id1) ;
    bdry_topo->InsertAngle(1,qdata->tr,qdata->front_nodes[3],mid_id0) ;

    if (base) {
        EdgeHandle *handle, new_handle ;
        ArbIntQdEdgeDesc desc, orig_desc ;

        handle = edge_table->Fetch(ArbEdgeKey(qdata->tl,qdata->tr)) ;
        orig_desc = *(order_heap->ViewWithHandle(*handle)) ;
        double edge_length = orig_desc.length ;
        edge_lengths->Remove(edge_length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Remove(ArbEdgeKey(qdata->tl,qdata->tr)) ;

        // now add the information for the new edges.

        desc = orig_desc ;
        desc.length = edge_length / 3.0 ;
        desc.angle_1 = PI ;
        desc.id_1 = mid_id1 ;
        desc.end_code &= LEFT_EDGE_MASK ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(qdata->tl,mid_id1),new_handle) ;

        desc.angle_0 = PI ;
        desc.id_0 = mid_id1 ;
        desc.id_1 = mid_id0 ;
        desc.end_code = 0 ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(mid_id1,mid_id0),new_handle) ;

        desc.angle_1 = orig_desc.angle_1 ;
        desc.id_0 = mid_id0 ;
        desc.id_1 = qdata->tr ;
        desc.end_code = orig_desc.end_code & RIGHT_EDGE_MASK ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(ArbEdgeKey(mid_id0,qdata->tr),new_handle) ;
    } else {
        edge_table->Remove(ArbEdgeKey(qdata->tl,qdata->tr)) ;

        edge_table->Store(ArbEdgeKey(qdata->tl,mid_id1),0) ;
        edge_table->Store(ArbEdgeKey(mid_id1,mid_id0),0) ;
        edge_table->Store(ArbEdgeKey(mid_id0,qdata->tr),0) ;
    }

    front_nodes->Store(qdata->tl,mid_id1) ;
    front_nodes->Store(mid_id1,mid_id0) ;
    front_nodes->Store(mid_id0,qdata->tr) ;
}




// %(CArbMshEdgeList::UpdateGeom-void-|-ArbNewQuadData-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateGeom - update stored geometry parameters
**
**      void UpdateGeom(ArbNewQuadData *qdata)
**
**        qdata - (in)  quadrilateral description
**
**      Description: This function updates the stored mesh front edge
**          lengths and the angles associated with a new quad. This
**          function is called after nodal smoothing.
**
**
** -- */

int CArbMshEdgeList::UpdateGeom(ArbNewQuadData *qdata)
{
    if ((qdata->qcase < 5) || (qdata->qcase == 6)) {
        if (!UpdateGeomHelp(qdata->num_front_nodes,qdata->front_nodes,
                       qdata->edge_level))
                       return 0;
    } else {
        if (!UpdateGeomHelp(4,&qdata->front_nodes[0],qdata->edge_level))
         return 0;
        if (!UpdateGeomHelp(4,&qdata->front_nodes[4],qdata->edge_level))
         return 0;
    }

    if ((qdata->qcase == 3) || (qdata->qcase == 4))  {
        if (qdata->bl != qdata->obl) {
            int front_nodes[3] ;
            if (UpdateGeomOne(qdata->bl,front_nodes))
                if (!UpdateGeomHelp(3,front_nodes,qdata->edge_level))
                 return 0;
        }
    }
    if ((qdata->qcase >= 2) && (qdata->qcase <= 4)) {
        if (qdata->br != qdata->obr) {
            int front_nodes[3] ;
            if (UpdateGeomOne(qdata->br,front_nodes))
                if (!UpdateGeomHelp(3,front_nodes,qdata->edge_level))
                 return 0;
        }
    }
    return 1;
}




// %(CArbMshEdgeList::UpdateGeomHelp-void-|-int-const|-int-const|*-int-const|)
/* ++ ----------------------------------------------------------
**
**    UpdateGeomHelp - support for UpdateGeom
**
**      void UpdateGeomHelp(
**              const int num_front_nodes,
**              const int *front_nodes,
**              const int edge_level)
**
**        num_front_nodes - (in)  number of front nodes associated with
**                                this quad
**        front_nodes     - (in)  list of front nodes
**        edge_level      - (in)  edge level to assign to updated edges
**
**      Description: This function is a support routine for the
**          UpdateGeom routine. It updates an edges end codes and edge
**          lengths.
**
**
** -- */

static int EndCode(double *angle, int /*level*/,
                   int /*adj_level*/, bool fixed_flag=false)
{
    if (*angle < QUAD_EDGE_SEAM_TOLERANCE){
        return(SEAM_CODE) ;
    } else if (*angle > ASSUME_OVERLAP_TOLERANCE) {
        *angle = OVERLAP_ANGLE_VALUE ;
        return(SEAM_CODE) ;
    } else {
        if (fixed_flag) {
            if (*angle < QUAD_BOUNDARY_EDGE_ANGLE_TOLERANCE) {
                return(ANGLE_CODE) ;
            }
        } else {
            if (*angle < QUAD_EDGE_ANGLE_TOLERANCE) {
                return(ANGLE_CODE) ;
            }
        }
    }
    return(0) ;
}

int CArbMshEdgeList::UpdateGeomHelp(const int num_front_nodes,
                                     const int *front_nodes,
                                     const int edge_level)
{
   ArbIntNode *nd_cw, *nd_ccw, *nd_0, *nd_1 ;
    int vtx_cw, vtx_ccw, vtx_0, vtx_1 ;
    ArbIntNode *adj_nd_0, *adj_nd_1 ;
    int adj_vtx_0, adj_vtx_1 ;
    double dx, dy, angle ;
    EdgeHandle *handle ;
    short adj_level, end_code ;
    ArbIntQdEdgeDesc desc, adj_desc/*, *desc_ptr*/ ;

    // save edge level of the ccw edge

    handle = edge_table->Fetch(ArbEdgeKey(
               front_nodes[num_front_nodes-2],
               front_nodes[num_front_nodes-1])) ;
    if (handle == 0)
     return 0;
    if (*handle == 0)
     return 0;
    /*desc_ptr = */order_heap->ViewWithHandle(*handle) ;

    // first look at the edge that is twice clockwise from
    // where we are adding the new element.  For this edge
    // we update the angle on the right end.

    vtx_0 = front_nodes[0] ;
    vtx_1 = front_nodes[1] ;
    vtx_cw = bdry_topo->GetCWBdryNode(vtx_0,vtx_1) ;

    if (vtx_cw == -1)
     return 0;

    nd_0 = node_table->Fetch(vtx_0) ;
    nd_1 = node_table->Fetch(vtx_1) ;
    nd_cw = node_table->Fetch(vtx_cw) ;

    if (nd_0 == 0 || nd_1 == 0 || nd_cw == 0)
     return 0;

    handle = edge_table->Fetch(ArbEdgeKey(vtx_cw,vtx_0)) ;

    if (handle == 0)
     return 0;
    if (*handle == 0)
     return 0;
    desc = *(order_heap->ViewWithHandle(*handle)) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_lengths->Remove(desc.length) ;

    dx = nd_cw->coord[0] - nd_0->coord[0] ;
    dy = nd_cw->coord[1] - nd_0->coord[1] ;
    desc.length = sqrt(dx*dx + dy*dy) ;

    angle = Angle2Pi(nd_0->coord,nd_1->coord,nd_cw->coord) ;
    end_code = EndCode(&angle,desc.level,desc.adj_level_1) ;

    desc.angle_1 = angle ;
    desc.end_code &= LEFT_EDGE_MASK ;
    desc.end_code |= end_code ;
    desc.visited = false ;

    *handle = order_heap->InsertWithHandle(desc) ;
    edge_lengths->Insert(desc.length) ;

    // now look at the edge that is clockwise from
    // where we are adding the new element.  For this
    // edge we need to update it's length, and reclassify
    // the angle on it's right

    vtx_ccw = front_nodes[2] ;

    handle = edge_table->Fetch(ArbEdgeKey(vtx_0,vtx_1)) ;

    if (handle == 0)
     return 0;
    if (*handle == 0)
     return 0;

    desc = *(order_heap->ViewWithHandle(*handle)) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_lengths->Remove(desc.length) ;

    desc.angle_0 = angle ;
    desc.end_code &= RIGHT_EDGE_MASK ;
    desc.end_code |= end_code << 1 ;

    nd_ccw = node_table->Fetch(vtx_ccw) ;

    if (nd_ccw == 0)
     return 0;

    dx = nd_0->coord[0] - nd_1->coord[0] ;
    dy = nd_0->coord[1] - nd_1->coord[1] ;
    desc.length = sqrt(dx*dx + dy*dy) ;
    angle = Angle2Pi(nd_1->coord,nd_ccw->coord,nd_0->coord) ;

    desc.end_code &= LEFT_EDGE_MASK ;
    bool fix_flg = (nd_1->motion == ARB_FIXED) &&
                   (nd_0->motion == ARB_FIXED) &&
                   (nd_ccw->motion == ARB_FIXED) ;
    end_code = EndCode(&angle,desc.level,desc.adj_level_1,fix_flg) ;
    desc.end_code |= end_code ;
    desc.angle_1 = angle ;
    desc.adj_level_1 = edge_level ;
    adj_level = desc.level ;
    desc.visited = false ;

    *handle = order_heap->InsertWithHandle(desc) ;
    edge_lengths->Insert(desc.length) ;

    // Now update the data for the new edges that
    // were added with this element

    for (int i=0 ; i<num_front_nodes-3 ; ++i) {
        vtx_0 = vtx_1 ;
        vtx_1 = front_nodes[i+2] ;
        nd_0 = nd_1 ;
        nd_1 = node_table->Fetch(vtx_1) ;

        if (nd_1 == 0)
         return 0;

        // update the edge stuff

        desc.id_0 = vtx_0 ;
        desc.id_1 = vtx_1 ;
        dx = nd_0->coord[0] - nd_1->coord[0] ;
        dy = nd_0->coord[1] - nd_1->coord[1] ;
        desc.length = sqrt(dx*dx + dy*dy) ;
        desc.end_code = 0 ;
        desc.level = edge_level ;
        desc.adj_level_0 = adj_level ;
        desc.adj_level_1 = edge_level ;

        // now determine if we can use the previous angle
        // information

        if (vtx_ccw != vtx_1) {

            vtx_cw = bdry_topo->GetCWBdryNode(vtx_0,vtx_1) ;
            if (vtx_cw == -1)
             return 0;
            nd_cw = node_table->Fetch(vtx_cw) ;
            if (nd_cw == 0)
             return 0;
            angle = Angle2Pi(nd_0->coord,nd_1->coord,nd_cw->coord) ;

            // update the adjacent edge

            adj_vtx_0 = vtx_cw ;
            adj_vtx_1 = vtx_0 ;

            adj_nd_0 = node_table->Fetch(adj_vtx_0) ;
            if (adj_nd_0 == 0)
             return 0;
            adj_nd_1 = nd_0 ;

            handle = edge_table->Fetch(ArbEdgeKey(adj_vtx_0,adj_vtx_1)) ;
            if (handle == 0)
             return 0;
            if (*handle == 0)
             return 0;
            adj_desc = *(order_heap->ViewWithHandle(*handle)) ;
            order_heap->RemoveWithHandle(*handle) ;
            edge_lengths->Remove(adj_desc.length) ;

            dx = adj_nd_0->coord[0] - adj_nd_1->coord[0] ;
            dy = adj_nd_0->coord[1] - adj_nd_1->coord[1] ;
            adj_desc.length = sqrt(dx*dx + dy*dy) ;

            end_code = EndCode(&angle,adj_desc.level,adj_desc.adj_level_1) ;

            adj_desc.angle_1 = angle ;
            adj_desc.adj_level_1 = edge_level ;
            adj_desc.end_code &= LEFT_EDGE_MASK ;
            adj_desc.end_code |= end_code ;
            adj_desc.visited = false ;

            *handle = order_heap->InsertWithHandle(adj_desc) ;
            edge_lengths->Insert(adj_desc.length) ;
        }

        desc.angle_0 = angle ;
        desc.end_code |= end_code << 1 ;

        // compute the angle on the right side of the edge

        vtx_ccw = bdry_topo->GetCCWBdryNode(vtx_0,vtx_1) ;
        if (vtx_ccw == -1)
         return 0;
        nd_ccw = node_table->Fetch(vtx_ccw) ;
        if (nd_ccw == 0)
         return 0;

        angle = Angle2Pi(nd_1->coord,nd_ccw->coord,nd_0->coord) ;
        bool fix_flg = (nd_1->motion == ARB_FIXED) &&
                       (nd_0->motion == ARB_FIXED) &&
                       (nd_ccw->motion == ARB_FIXED) ;
        end_code = EndCode(&angle,desc.level,desc.adj_level_1,fix_flg) ;
        desc.angle_1 = angle ;
        desc.end_code |= end_code ;
        desc.visited = false ;

        handle = edge_table->Fetch(ArbEdgeKey(vtx_0,vtx_1)) ;
        if (handle == 0)
         return 0;
        *handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;

        // if necessary, update the adjacent edge on the right

        if (vtx_ccw != front_nodes[i+3]) {
            adj_vtx_0 = vtx_1 ;
            adj_vtx_1 = vtx_ccw ;

            adj_nd_0 = nd_1 ;
            adj_nd_1 = node_table->Fetch(adj_vtx_1) ;
            if (adj_nd_1 == 0)
             return 0;

            handle = edge_table->Fetch(ArbEdgeKey(adj_vtx_0,adj_vtx_1)) ;
            if (handle == 0)
             return 0;
            if (*handle == 0)
             return 0;
            adj_desc = *(order_heap->ViewWithHandle(*handle)) ;
            order_heap->RemoveWithHandle(*handle) ;
            edge_lengths->Remove(adj_desc.length) ;

            dx = adj_nd_0->coord[0] - adj_nd_1->coord[0] ;
            dy = adj_nd_0->coord[1] - adj_nd_1->coord[1] ;
            adj_desc.length = sqrt(dx*dx + dy*dy) ;

            adj_desc.angle_0 = angle ;
            adj_desc.adj_level_0 = edge_level ;
            adj_desc.end_code &= RIGHT_EDGE_MASK ;
            adj_desc.end_code |= end_code << 1 ;
            adj_desc.visited = false ;

            *handle = order_heap->InsertWithHandle(adj_desc) ;
            edge_lengths->Insert(adj_desc.length) ;
        }
    }

    // look at the edge that is counter clockwise from
    // where we are adding the new element.  For this
    // edge we need to update it's length, and reclassify
    // the angle on it's left end

    vtx_0 = vtx_1 ;
    vtx_1 = front_nodes[num_front_nodes-1] ;
    nd_0 = nd_1 ;
    nd_1 = node_table->Fetch(vtx_1) ;
    if (nd_1 == 0)
     return 0;

    handle = edge_table->Fetch(ArbEdgeKey(vtx_0,vtx_1)) ;
    if (handle == 0)
     return 0;
    if (*handle == 0)
     return 0;
    desc = *(order_heap->ViewWithHandle(*handle)) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_lengths->Remove(desc.length) ;

    dx = nd_0->coord[0] - nd_1->coord[0] ;
    dy = nd_0->coord[1] - nd_1->coord[1] ;
    desc.length = sqrt(dx*dx + dy*dy) ;

    if (vtx_ccw != vtx_1) {
        vtx_cw = bdry_topo->GetCWBdryNode(vtx_0,vtx_1) ;
        if (vtx_cw == -1)
         return 0;
        nd_cw = node_table->Fetch(vtx_cw) ;
        if (nd_cw == 0)
         return 0;
        angle = Angle2Pi(nd_0->coord,nd_1->coord,nd_cw->coord) ;
        end_code = EndCode(&angle,desc.level,desc.adj_level_0) ;
    }

    desc.angle_0 = angle ;
    desc.adj_level_0 = edge_level ;
    desc.end_code &= RIGHT_EDGE_MASK ;
    desc.end_code |= end_code << 1 ;

    // now update the angle on the right end

    vtx_ccw = bdry_topo->GetCCWBdryNode(vtx_0,vtx_1) ;
    if (vtx_ccw == -1)
     return 0;
    nd_ccw = node_table->Fetch(vtx_ccw) ;
    if (nd_ccw == 0)
     return 0;

    angle = Angle2Pi(nd_1->coord,nd_ccw->coord,nd_0->coord) ;
    end_code = EndCode(&angle,desc.level,desc.adj_level_1) ;
    desc.angle_1 = angle ;
    desc.end_code &= LEFT_EDGE_MASK ;
    desc.end_code |= end_code ;
    desc.visited = false ;

    *handle = order_heap->InsertWithHandle(desc) ;
    edge_lengths->Insert(desc.length) ;

    // now we update the angle for the next ccw edge

    handle = edge_table->Fetch(ArbEdgeKey(vtx_1,vtx_ccw)) ;
    if (handle == 0)
     return 0;
    if (*handle == 0)
     return 0;
    desc = *(order_heap->ViewWithHandle(*handle)) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_lengths->Remove(desc.length) ;

    dx = nd_1->coord[0] - nd_ccw->coord[0] ;
    dy = nd_1->coord[1] - nd_ccw->coord[1] ;
    desc.length = sqrt(dx*dx + dy*dy) ;

    desc.angle_0 = angle ;
    desc.end_code &= RIGHT_EDGE_MASK ;
    desc.end_code |= end_code << 1 ;
    desc.visited = false ;

    *handle = order_heap->InsertWithHandle(desc) ;
    edge_lengths->Insert(desc.length) ;
    return 1;
}


bool CArbMshEdgeList::UpdateGeomOne(const int node_id,
                                    int front_nodes[])
{
    // if this node is not on the boundary then return

    ArbIntNode *node = node_table->Fetch(node_id) ;
    if (node == 0) return(false) ;

    CArbTopoAdjVtxIterator iter(bdry_topo,node_id) ;
    if (iter.More()) {
        front_nodes[1] = node_id ;
        front_nodes[2] = iter.AdjVtx() ;
        ++iter ;
        front_nodes[0] = iter.AdjVtx() ;
        return(true) ;
    }
    return(false) ;
}


// %(CArbMshEdgeList::UpdateBdryGeom-void-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    UpdateBdryGeom - update stored geometry parameters
**
**      void UpdateBdryGeom(const int vtx)
**
**        vtx - (in)  vertex id
**
**      Description: This function updates the stored mesh front
**          lengths and angle associated with one vertex. This function
**          is called after nodal smoothing.
**
**
** -- */

void CArbMshEdgeList::UpdateBdryGeom(const int vtx)
{
    ArbIntNode *nd_cw, *nd_ccw, *nd_0, *nd_1 ;
    int vtx_cw, vtx_ccw, vtx_0, vtx_1 ;
    double dx, dy, angle ;
    EdgeHandle *handle ;
    short end_code ;
    ArbIntQdEdgeDesc desc ;

    // first look at the edge that is counter clockwise from
    // the vertex.

    vtx_0 = vtx ;

    CArbTopoAdjVtxIterator iter(bdry_topo,vtx) ;

    if (!iter.More()) return ;
    while (iter.More() && (iter.CcwElem() != 1)) ++iter ;

    vtx_1 = iter.AdjVtx() ;
    vtx_ccw = bdry_topo->GetCCWBdryNode(vtx_0,vtx_1) ;

    handle = edge_table->Fetch(ArbEdgeKey(vtx_0,vtx_1)) ;
    desc = *(order_heap->ViewWithHandle(*handle)) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_lengths->Remove(desc.length) ;

    nd_0 = node_table->Fetch(vtx_0) ;
    nd_1 = node_table->Fetch(vtx_1) ;
    nd_ccw = node_table->Fetch(vtx_ccw) ;

    dx = nd_0->coord[0] - nd_1->coord[0] ;
    dy = nd_0->coord[1] - nd_1->coord[1] ;
    desc.length = sqrt(dx*dx + dy*dy) ;
    angle = Angle2Pi(nd_1->coord,nd_ccw->coord,nd_0->coord) ;

    desc.end_code = 0 ;
    end_code = EndCode(&angle,desc.level,desc.adj_level_1) ;
    desc.end_code |= end_code ;
    desc.angle_1 = angle ;

    // now update the angle right at vertex

    vtx_cw = bdry_topo->GetCWBdryNode(vtx_0,vtx_1) ;
    nd_cw = node_table->Fetch(vtx_cw) ;
    angle = Angle2Pi(nd_0->coord,nd_1->coord,nd_cw->coord) ;
    end_code = EndCode(&angle,desc.level,desc.adj_level_0) ;
    desc.angle_0 = angle ;
    desc.end_code |= end_code << 1 ;
    desc.visited = false ;

    *handle = order_heap->InsertWithHandle(desc) ;
    edge_lengths->Insert(desc.length) ;

    // now move to the cw edge

    vtx_1 = vtx_0 ;
    vtx_0 = vtx_cw ;
    vtx_cw = bdry_topo->GetCWBdryNode(vtx_0,vtx_1) ;

    handle = edge_table->Fetch(ArbEdgeKey(vtx_0,vtx_1)) ;
    desc = *(order_heap->ViewWithHandle(*handle)) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_lengths->Remove(desc.length) ;

    nd_1 = nd_0 ;
    nd_0 = nd_cw ;
    nd_cw = node_table->Fetch(vtx_cw) ;

    dx = nd_0->coord[0] - nd_1->coord[0] ;
    dy = nd_0->coord[1] - nd_1->coord[1] ;
    desc.length = sqrt(dx*dx + dy*dy) ;
    desc.end_code |= end_code ;
    desc.angle_1 = angle ;

    angle = Angle2Pi(nd_0->coord,nd_1->coord,nd_cw->coord) ;
    end_code = EndCode(&angle,desc.level,desc.adj_level_0) ;
    desc.angle_0 = angle ;
    desc.end_code |= end_code << 1 ;
    desc.visited = false ;

    *handle = order_heap->InsertWithHandle(desc) ;
    edge_lengths->Insert(desc.length) ;
}




// %(CArbMshEdgeList::UpdateSeamGeom-void-|-ArbNewSeamData-const|*)
/* ++ ----------------------------------------------------------
**
**    UpdateSeamGeom - update stored geometry parameters
**
**      void UpdateSeamGeom(const ArbNewSeamData *sdata)
**
**        sdata - (in)  seam data
**
**      Description: This function updates the stored mesh front
**          lengths and angle. This function is called after nodal
**          smoothing.
**
**
** -- */

void CArbMshEdgeList::UpdateSeamGeom(const ArbNewSeamData *sdata)
{
    int front_nodes[3] ;

    front_nodes[0] = sdata->cw ;
    front_nodes[1] = sdata->id0 ;
    front_nodes[2] = sdata->ccw ;
    UpdateGeomHelp(3,front_nodes,1) ;
}




// %(CArbMshEdgeList::GetCWNode-int-|^const-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetCWNode - find the clockwise node
**
**      int GetCWNode(
**              const int id_0,
**              const int id_1) const
**
**        id_0 - (in)  first edge node
**        id_1 - (in)  second edge node
**
**      Description: This function returns the id of the next node on
**          the crack front that is clockwise from the edge defined by
**          the arguments.
**
**      Return Value: The id of the node clockwise from edge.
**
**
** -- */

int CArbMshEdgeList::GetCWNode(const int id_0,
                               const int id_1) const
{
    return(bdry_topo->GetCWBdryNode(id_0,id_1)) ;
}




// %(CArbMshEdgeList::GetCCWNode-int-|^const-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetCCWNode - find counter-clockwise node
**
**      int GetCCWNode(
**              const int id_0,
**              const int id_1) const
**
**        id_0 - (in)  first edge node
**        id_1 - (in)  second edge node
**
**      Description: This function returns the id of the next node on
**          the crack front that is counter- clockwise from the edge
**          defined by the arguments.
**
**      Return Value: The id of the node counter-clockwise from edge.
**
**
** -- */

int CArbMshEdgeList::GetCCWNode(const int id_0,
                                const int id_1) const
{
    return(bdry_topo->GetCCWBdryNode(id_0,id_1)) ;
}




// %(CArbMshEdgeList::GetEdgeLength-double-|^const-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetEdgeLength - return the length of an edge
**
**      double GetEdgeLength(
**              const int id_0,
**              const int id_1) const
**
**        id_0 - (in)  first edge node
**        id_1 - (in)  second edge node
**
**      Description: This function returns the length of the edge
**          defined by the arguments.
**
**      Return Value: The length of the edge.
**
**
** -- */

double CArbMshEdgeList::GetEdgeLength(const int id_0,
                                      const int id_1) const
{
    EdgeHandle *handle ;
    ArbIntQdEdgeDesc *desc_ptr ;
    handle = edge_table->Fetch(ArbEdgeKey(id_0,id_1)) ;
    desc_ptr = order_heap->ViewWithHandle(*handle) ;
    return(desc_ptr->length) ;
}




// %(CArbMshEdgeList::GetEdgeLevel-int-|^const-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetEdgeLevel - return the level of an edge
**
**      int GetEdgeLevel(
**              const int id_0,
**              const int id_1) const
**
**        id_0 - (in)  first edge node
**        id_1 - (in)  second edge node
**
**      Description: This function returns the level (approximate
**          number of steps required to get to the original boundary)
**          of the edge defined by the arguments.
**
**      Return Value: The level of the edge.
**
**
** -- */

int CArbMshEdgeList::GetEdgeLevel(const int id_0,
                                           const int id_1) const
{
    EdgeHandle *handle ;
    ArbIntQdEdgeDesc *desc_ptr ;
    handle = edge_table->Fetch(ArbEdgeKey(id_0,id_1)) ;
    desc_ptr = order_heap->ViewWithHandle(*handle) ;
    return(desc_ptr->level) ;
}




// %(CArbMshEdgeList::ContainsEdge-bool-|^const-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    ContainsEdge - check to see if an edge is in the list
**
**      bool ContainsEdge(
**              const int id_0,
**              const int id_1) const
**
**        id_0 - (in)  first edge node
**        id_1 - (in)  second edge node
**
**      Description: This function checks to see if the edge define by
**          the arguments is in the list.
**
**      Return Value: True if the edge is in the list, false otherwise.
**
**
** -- */

bool CArbMshEdgeList::ContainsEdge(const int id_0,
                                   const int id_1) const
{
    EdgeHandle *tmp = edge_table->Fetch(ArbEdgeKey(id_0,id_1)) ;
    return((tmp != 0) ? true : false) ;
}




// %(CArbMshEdgeList::ContainsNode-bool-|^const-int-const|)
/* ++ ----------------------------------------------------------
**
**    ContainsNode - checks to see if a node is in the list
**
**      bool ContainsNode(const int id) const
**
**        id - (in)  node id
**
**      Description: This function checks to see if the node specified
**          by the argument is contained in in the edge list.
**
**      Return Value: True if the node is found, false, otherwise
**
**
** -- */

bool CArbMshEdgeList::ContainsNode(const int id) const
{
    return((front_nodes->Fetch(id) != 0) ? true : false) ;
}




// %(CArbMshEdgeList::GetNextEdge-ArbQdEdge-|*)
/* ++ ----------------------------------------------------------
**
**    GetNextEdge - get next edge to use to advance the front
**
**      ArbQdEdge *GetNextEdge()
**
**      Description: This function returns information about the edge
**          that is ranked highest in the edge list's priority queue
**          for use in advancing the mesh front.
**
**      Return Value: An ArbQdEdge structure containing information the
**          highest ranked edge.
**
**
** -- */

CArbMshEdgeList::ArbQdEdge *CArbMshEdgeList::GetNextEdge()
{
    static ArbQdEdge edge ;

    ArbIntQdEdgeDesc *next = order_heap->ViewMin() ;
    if (next == 0) return(0) ;
    next->visited = true ;
    edge.id_0 = next->id_0 ;
    edge.id_1 = next->id_1 ;
    edge.end_code = next->end_code ;
    edge.angle_0 = next->angle_0 ;
    edge.angle_1 = next->angle_1 ;
    edge.level = next->level ;

// For Debug

//     CArbHeap<ArbIntQdEdgeDesc> tmp_heap = *order_heap ;
//     ArbIntQdEdgeDesc *tmp = tmp_heap.GetMin() ;
//     int i = 0 ;
//     while (tmp != 0) {
//         fprintf(stderr,"%d %d %d %d %g %g %d %g\n",i,
//                 tmp->id_0,tmp->id_1,tmp->end_code,
//                 tmp->angle_0,tmp->angle_1,tmp->level,
//                 tmp->length) ;
//         tmp = tmp_heap.GetMin() ;
//         i++ ;
//     }

    return(&edge) ;
}




// %(CArbMshEdgeList::GetEdgeLengthRatio-double-|^const)
/* ++ ----------------------------------------------------------
**
**    GetEdgeLengthRatio - get the edge length ratio
**
**      double GetEdgeLengthRatio() const
**
**      Description: This function returns the ratio of the longest to
**          shortest edges currently in the edge list.
**
**      Return Value: The ratio of the longest to shortest edge.
**
**
** -- */

double CArbMshEdgeList::GetEdgeLengthRatio() const
{
    const double *small, *big ;
    small = edge_lengths->GetSmallest() ;
    //assert(small > 0.0) ;
    big = edge_lengths->GetLargest() ;
    return((*big)/(*small)) ;
}




// %(CArbMshEdgeList::PushToBack-void-|-ArbQdEdge-const|*)
/* ++ ----------------------------------------------------------
**
**    PushToBack - put an edge back on the edge list
**
**      void PushToBack(const ArbQdEdge *edge)
**
**        edge - (in)  edge to place back on the list
**
**      Description: This function places an edge back on the edge list
**          priority queue, and flags it so that it will be ranked
**          lowest in the queue.
**
**
** -- */

void CArbMshEdgeList::PushToBack(const ArbQdEdge *edge)
{
    EdgeHandle *tedge =
        edge_table->Fetch(ArbEdgeKey(edge->id_0,edge->id_1)) ;

    if (tedge == 0) return ;

    ArbIntQdEdgeDesc desc = *(order_heap->ViewWithHandle(*tedge)) ;

    order_heap->RemoveWithHandle(*tedge) ;
    desc.level++ ;
    *tedge = order_heap->InsertWithHandle(desc) ;
}


bool CArbMshEdgeList::StagnationSweep()
{
    // get a list of the entries in the heap

    ArbIntQdEdgeDesc **list = order_heap->GetEntryList() ;

    // look for an unvisited edge

    bool stagnant = true ;
    for (int i=0 ; i<order_heap->NumEntries() ; ++i) {
         if (!list[i]->visited) stagnant = false ;
         list[i]->visited = false ;
    }
    delete [] list ;
    return(stagnant) ;
}

double CArbMshEdgeList::GetCharNodeLength(const int id_0) const
{
    int *ptr = front_nodes->Fetch(id_0) ;
    if (ptr == 0) return(0.0) ;
    int id_1 = *ptr ;

    if (!ContainsEdge(id_0,id_1)) return(1.0) ;

    ArbIntNode *node_0 = node_table->Fetch(id_0) ;
    ArbIntNode *node_1 = node_table->Fetch(id_1) ;

    double dist_0 = (node_0->coord - node_1->coord).Magnitude() ;

    int id_2 = GetCWNode(id_0,id_1) ;
    ArbIntNode *node_2 = node_table->Fetch(id_2) ;

    double dist_1 = (node_0->coord - node_2->coord).Magnitude() ;

    return(0.5*(dist_0+dist_1)) ;
}


// %(CArbMshEdgeList::Angle-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    Angle - compute the included angle
**
**      double Angle(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j) const
**
**        b - (in)  center point
**        i - (in)  first point
**        j - (in)  second point
**
**      Description: This function computes the included angle between
**          three points.
**
**      Return Value: Included angle, -PI < angle < PI
**
**
** -- */

double CArbMshEdgeList::Angle(CArbCoord2D b,CArbCoord2D i,
                              CArbCoord2D j) const
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




// %(CArbMshEdgeList::Angle2Pi-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    Angle2Pi - compute the included angle
**
**      double Angle2Pi(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j) const
**
**        b - (in)  center point
**        i - (in)  first point
**        j - (in)  second point
**
**      Description: This function computes the included angle between
**          three points.
**
**      Return Value: Included angle, 0 < angle < 2*PI
**
**
** -- */

double CArbMshEdgeList::Angle2Pi(CArbCoord2D b,CArbCoord2D i,
                                 CArbCoord2D j) const
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




// %(CArbMshEdgeList::CrossProd-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    CrossProd - cross product of two vectors
**
**      double CrossProd(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j) const
**
**        b - (in)  common start point
**        i - (in)  first end point
**        j - (in)  second end point
**
**      Description: This function computes the cross product of two
**          vectors. The vectors are defined by one common start point
**          and two separate end points.
**
**      Return Value: the cross product
**
**
** -- */

double CArbMshEdgeList::CrossProd(CArbCoord2D b,CArbCoord2D i,
                                  CArbCoord2D j) const
{
    double cross ;

    cross = ((i[0] - b[0]) * (j[1] - b[1])) -
            ((i[1] - b[1]) * (j[0] - b[0])) ;
    return cross ;
}
