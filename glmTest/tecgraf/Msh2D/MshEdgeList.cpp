//
// EdgeList Class definition
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
//   $Revision: 1.23 $  $Date: 2004/06/01 15:46:45 $  $Author: wash $
//

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include "Dict.hpp"

#include "MshEdgeList.hpp"

using FTools::Dict ;

namespace Msh2D {

#ifdef MEMDEBUG
#include "MemDbg.hpp"
//#define new new(__FILE__,__LINE__)
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


/* ------------------------------------------------------------
    CArbCmpQdEdge - comparison function used to order edges
                    in the list so that the best candidate for
                    creating a new quad can be selected
*/

int MshEdgeList::EdgePriority::Compare(
    const MshEdgeList::IntQdEdgeDesc &edg1,
    const MshEdgeList::IntQdEdgeDesc &edg2)
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




// %(MshEdgeList::MshEdgeList-constructor-|-Dict-|<int,IntNode>*) 
/* ++ ----------------------------------------------------------
**
**    MshEdgeList - edge list constructor 
**
**      MshEdgeList(Dict<int,IntNode>* inode_table)
**
**        inode_table - (in)  node hash table 
**
**      Description: This is a constructor for an edge list. As an 
**          argument it takes a pointer to a hash table that maps node 
**          id's to IntNodes. 
**
**
** -- */

MshEdgeList::MshEdgeList(
    Dict<int,IntNode> *inode_table) :
    node_table(inode_table)
{
    edge_table = new Dict<EdgeKey,
                          EdgePriority::EntryHandle>(false) ;
    order_heap = new EdgePriority() ;
    edge_lengths = new EdgeLengthSet() ;
    front_nodes = new Dict<int,int>(true) ;
    bdry_topo = new MshTopo2D() ;
}




// %(MshEdgeList::MshEdgeList-destructor-|~) 
/* ++ ----------------------------------------------------------
**
**    MshEdgeList - edge list destructor 
**
**      ~MshEdgeList()
**
**      Description: This is a desctructor for an edge list. 
**
**
** -- */

MshEdgeList::~MshEdgeList()
{
    delete edge_table ;
    delete order_heap ;
    delete edge_lengths ;
    delete front_nodes ;
    delete bdry_topo ;
}




// %(MshEdgeList::InsertEdge-void-|-int-const|-int-const|-int-const|-int-const|) 
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

void MshEdgeList::InsertEdge(const int id_0,
                                 const int id_1,
                                 const int cw_id,
                                 const int ccw_id)
{
    // insert this information into the edge and back
    // tables.

    edge_table->Store(EdgeKey(id_0,id_1),0) ;
    bdry_topo->InsertAngle(1,id_1,ccw_id,id_0) ;
    front_nodes->Store(id_0,id_1) ;
}




// %(MshEdgeList::InitializeCodes-void-|)
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

void MshEdgeList::InitializeCodes()
{
    Dict<EdgeKey,EdgePriority::EntryHandle>::DictIterator iter(edge_table) ;

    // go through all the edges.  Compute the angle codes
    // and length and insert this information in the selection
    // heap 

    for (iter.First() ; iter.More() ; ++iter) {
        IntQdEdgeDesc edge_desc ;
        IntNode *p_nd0, *p_nd1, *adj ;
        double dx, dy ;

        // get pointers to the node info and initialize the
        // description

        p_nd0 = node_table->Get(iter.Key().id_0) ;
        p_nd1 = node_table->Get(iter.Key().id_1) ;

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

        adj = node_table->Get(bdry_topo->GetCWBdryNode(iter.Key().id_0,
                                                         iter.Key().id_1)) ;
        edge_desc.angle_0 = Angle2Pi(p_nd0->coord,p_nd1->coord,adj->coord) ;
        if (edge_desc.angle_0 < QUAD_EDGE_ANGLE_TOLERANCE)
            edge_desc.end_code = 2 ;

        // find the include angle at node 1 and classify the
        // the vertex

        adj = node_table->Get(bdry_topo->GetCCWBdryNode(iter.Key().id_0,
                                                          iter.Key().id_1)) ;
        edge_desc.angle_1 = Angle2Pi(p_nd1->coord,adj->coord,p_nd0->coord) ;
        if (edge_desc.angle_1 < QUAD_EDGE_ANGLE_TOLERANCE)
            edge_desc.end_code += 1 ;

        iter.Entry() = order_heap->InsertWithHandle(edge_desc) ;
        edge_lengths->Insert(edge_desc.length) ;
    }
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




// %(MshEdgeList::ClassifyQuad-void-|-int-const|-int-const|-int-const|-int-const|-NewQuadData-|*)
/* ++ ----------------------------------------------------------
**
**    ClassifyQuad - find information for a new quadrilateral 
**
**      void ClassifyQuad(
**              const int      id_0,
**              const int      id_1,
**              const int      id_2,
**              const int      id_3,
**              NewQuadData *qdata)
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

void MshEdgeList::ClassifyQuad(const int id0,
                                   const int id1,
                                   const int id2,
                                   const int id3,
                                   NewQuadData *qdata)
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
        return ;
    }

    // first determine which of the edges exist so that
    // we know which case to deal with

    bool exists[4] ;
    int ii ;

    for (ii=0 ; ii<4 ; ++ii) exists[ii] = false ;

    if (edge_table->Get(EdgeKey(id0,id1)) != 0) exists[0] = true ;
    if (edge_table->Get(EdgeKey(id1,id2)) != 0) exists[1] = true ;
    if (edge_table->Get(EdgeKey(id2,id3)) != 0) exists[2] = true ;
    if (edge_table->Get(EdgeKey(id3,id0)) != 0) exists[3] = true ;

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
            qdata->front_nodes[1] = qdata->bl ;
            qdata->front_nodes[2] = qdata->tl ;
            qdata->front_nodes[3] = qdata->tr ;
            qdata->front_nodes[4] = qdata->br ;
            qdata->front_nodes[5] = 
                bdry_topo->GetCCWBdryNode(qdata->bl,qdata->br) ;
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
                qdata->front_nodes[1] = qdata->bl ;
                qdata->front_nodes[2] = qdata->tl ;
                qdata->front_nodes[3] = qdata->tr ;
                qdata->front_nodes[4] = 
                    bdry_topo->GetCCWBdryNode(qdata->br,qdata->tr) ;
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
                qdata->front_nodes[1] = qdata->bl ;
                qdata->front_nodes[2] = qdata->tl ;
                qdata->front_nodes[3] = 
                    bdry_topo->GetCCWBdryNode(qdata->tr,qdata->tl) ;

                qdata->front_nodes[4] = 
                    bdry_topo->GetCWBdryNode(qdata->tr,qdata->tl) ;
                qdata->front_nodes[5] = qdata->tr ;
                qdata->front_nodes[6] = qdata->br ;
                qdata->front_nodes[7] = 
                    bdry_topo->GetCCWBdryNode(qdata->bl,qdata->br) ;

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
            qdata->front_nodes[1] = qdata->tl ;
            qdata->front_nodes[2] = qdata->tr ;
            qdata->front_nodes[3] = 
                bdry_topo->GetCCWBdryNode(qdata->br,qdata->tr) ;
            break ;

        case 4:
            Rotate(0,id0,id1,id2,id3,&qdata->bl,
                   &qdata->br,&qdata->tr,&qdata->tl) ;
            qdata->num_front_nodes = 0 ;
            break ;
    }
}

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

static double CrossProd(const Vec2D &b,const Vec2D &i,
                        const Vec2D &j)
{
    double cross ;

    cross = ((i[0] - b[0]) * (j[1] - b[1])) -
            ((i[1] - b[1]) * (j[0] - b[0])) ;
    return cross ;
}

static bool Cross(const Vec2D &i1,const Vec2D &i2,
                  const Vec2D &j1,const Vec2D &j2)
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

    // Now compute the cross product of line I with J1 and line I
    // with J2.  If have the same sign the lines cannot cross

    Vec2D delt = i2 - i1 ;
    double tol = -1e-10 * delt.Magnitude() ;

    if ((CrossProd(i1,i2,j1) * CrossProd(i1,i2,j2)) >= tol)
        return false ;

    // find the cross product of line J with I1 and line J with I2.
    // If they have the same sign, the lines cannot cross

    if ((CrossProd(j1,j2,i1) * CrossProd(j1,j2,i2)) >= tol)
        return false ;

    return true ;
}

bool MshEdgeList::DoesThisCrossBoundary(
                           const Vec2D &pt0,
                           const Vec2D &pt1) const
{
    // loop through all the edges in the front

    bool crossed = false ;

    Dict<EdgeKey,EdgePriority::EntryHandle>::ConstDictIterator iter(edge_table) ;
    for (iter.First() ; iter.More() ; ++iter) {
        EdgeKey edge = iter.Key() ;
        IntNode *nd0 = node_table->Get(edge.id_0) ;
        IntNode *nd1 = node_table->Get(edge.id_1) ;
        if ((nd0->coord != pt0) && (nd1->coord != pt1)) {
            if (Cross(nd0->coord,nd1->coord,pt0,pt1)) {
                crossed = true ;
                break ;
            }
        }
    }
    return(crossed) ;
}


// %(MshEdgeList::UpdateQuad-void-|-NewQuadData-|*) 
/* ++ ----------------------------------------------------------
**
**    UpdateQuad - update the list to add a quad element 
**
**      void UpdateQuad(NewQuadData *qdata)
**
**        qdata - (in)  quad description 
**
**      Description: This function is called to update the edge list to 
**          add a new quadrilateral to the mesh. 
**
**
** -- */

void MshEdgeList::UpdateQuad(NewQuadData *qdata)
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


/* ++ ----------------------------------------------------------
**
**    RemoveQuad - update the list to remove a quad element 
**
**      void RemoveQuad(NewQuadData *qdata)
**
**        qdata - (in)  quad description 
**
**      Description: This function is called to update the edge list to 
**          remove quadrilateral to the mesh. 
**
**
** -- */

void MshEdgeList::RemoveQuad(NewQuadData *qdata)
{

    // currently only case 1 is implemented

    RemoveCase1(qdata) ;
}


// %(MshEdgeList::FindAdjacent-void-|-int-const|-int-const|-int-|*-int-|*) 
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

void MshEdgeList::FindAdjacent(const int vtx,
                                   const int new_vtx,
                                   int *prev_vtx, 
                                   int *next_vtx)
{
    IntNode *node  = node_table->Get(vtx) ;
    IntNode *nnode = node_table->Get(new_vtx) ;

    TopoAdjVtxIterator iter(bdry_topo,vtx) ;
    int tprev_vtx ;

    IntNode *prev = node_table->Get(iter.AdjVtx()) ;
    tprev_vtx = iter.AdjVtx() ;

    while (iter.More()) {
        ++iter ;
        IntNode *next = node_table->Get(iter.AdjVtx()) ;

        if (Angle2Pi(node->coord,prev->coord,nnode->coord) <
            Angle2Pi(node->coord,prev->coord,next->coord)) {
            *prev_vtx = tprev_vtx ;
            *next_vtx = iter.AdjVtx() ;
            return ;
        }

        ++iter ;
        prev = node_table->Get(iter.AdjVtx()) ;
        tprev_vtx = iter.AdjVtx() ;
    }
}




// %(MshEdgeList::UpdateCase1-int-|-NewQuadData-|*) 
/* ++ ----------------------------------------------------------
**
**    UpdateCase1 - update the edge list for a case 1 quad 
**
**      int UpdateCase1(NewQuadData *qdata)
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

int MshEdgeList::UpdateCase1(NewQuadData *qdata)
{
    int level ;
    double length ;

    // determine the level to assign to the new edges

    EdgePriority::EntryHandle *base = 
        edge_table->Get(EdgeKey(qdata->bl,qdata->br)) ;

    IntQdEdgeDesc *desc_ptr = order_heap->ViewWithHandle(*base) ;

    level = desc_ptr->level + 1 ;
    length = desc_ptr->length ;

    // delete the base information from the tables

    order_heap->RemoveWithHandle(*base) ;
    edge_lengths->Remove(length) ;
    edge_table->Del(EdgeKey(qdata->bl,qdata->br)) ;
    
    // create the new edges

    edge_table->Store(EdgeKey(qdata->bl,qdata->tl),0) ;
    edge_table->Store(EdgeKey(qdata->tl,qdata->tr),0) ;
    edge_table->Store(EdgeKey(qdata->tr,qdata->br),0) ;

    front_nodes->Store(qdata->bl,qdata->tl) ;
    front_nodes->Store(qdata->tl,qdata->tr) ;
    front_nodes->Store(qdata->tr,qdata->br) ;

    // update the edges before and after the edge

    bdry_topo->DeleteAngle(qdata->bl,qdata->br,qdata->front_nodes[0]) ;
    bdry_topo->InsertAngle(1,qdata->bl,qdata->tl,qdata->front_nodes[0]) ;

    if (bdry_topo->HasVtx(qdata->tl)) {
        int prev_vtx, next_vtx ;
        FindAdjacent(qdata->tl,qdata->bl,&prev_vtx,&next_vtx) ;
        bdry_topo->DeleteAngle(qdata->tl,prev_vtx,next_vtx) ;
        bdry_topo->InsertAngle(1,qdata->tl,prev_vtx,qdata->bl) ;
        bdry_topo->InsertAngle(1,qdata->tl,qdata->tr,next_vtx) ;
    } else {
        bdry_topo->InsertAngle(1,qdata->tl,qdata->tr,qdata->bl) ;
    }

    if (bdry_topo->HasVtx(qdata->tr)) {
        int prev_vtx, next_vtx ;
        FindAdjacent(qdata->tr,qdata->tl,&prev_vtx,&next_vtx) ;
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

/* ++ ----------------------------------------------------------
**
**    RemoveCase1 - update the edge list for a case 1 quad 
**
**      int RemobeCase1(NewQuadData *qdata)
**
**        qdata - (in)  quadrilateral data 
**
**      Description: This function updates the edge list to remove a case 
**          1 quadrilateral. A case 1 quad has one edge on the front. 
**
**
** -- */

void MshEdgeList::RemoveCase1(NewQuadData *qdata)
{
    int level ;
    double length[3] ;

    // determine the level to assign to the new edges

    EdgePriority::EntryHandle *base0 = 
        edge_table->Get(EdgeKey(qdata->bl,qdata->tl)) ;
    EdgePriority::EntryHandle *base1 = 
        edge_table->Get(EdgeKey(qdata->tl,qdata->tr)) ;
    EdgePriority::EntryHandle *base2 = 
        edge_table->Get(EdgeKey(qdata->tr,qdata->br)) ;

    IntQdEdgeDesc *desc_ptr0 = order_heap->ViewWithHandle(*base0) ;
    IntQdEdgeDesc *desc_ptr1 = order_heap->ViewWithHandle(*base1) ;
    IntQdEdgeDesc *desc_ptr2 = order_heap->ViewWithHandle(*base2) ;

    level = desc_ptr0->level - 1 ;
    length[0] = desc_ptr0->length ;
    length[1] = desc_ptr1->length ;
    length[2] = desc_ptr2->length ;

    // delete the base information from the tables

    order_heap->RemoveWithHandle(*base0) ;
    order_heap->RemoveWithHandle(*base1) ;
    order_heap->RemoveWithHandle(*base2) ;
    edge_lengths->Remove(length[0]) ;
    edge_lengths->Remove(length[1]) ;
    edge_lengths->Remove(length[2]) ;
    edge_table->Del(EdgeKey(qdata->bl,qdata->tl)) ;
    edge_table->Del(EdgeKey(qdata->tl,qdata->tr)) ;
    edge_table->Del(EdgeKey(qdata->tr,qdata->br)) ;
    
    // create the new edge

    IntQdEdgeDesc desc ;

    IntNode *nd_0 = node_table->Get(qdata->bl) ;
    IntNode *nd_1 = node_table->Get(qdata->br) ;

    double dx = nd_1->coord[0] - nd_0->coord[0] ;
    double dy = nd_1->coord[1] - nd_0->coord[1] ;
    desc.length = sqrt(dx*dx + dy*dy) ;
    desc.angle_0 = PI ;
    desc.angle_1 = PI ;
    desc.id_0 = qdata->bl ;
    desc.id_1 = qdata->br ;
    desc.end_code &= LEFT_EDGE_MASK ;
    desc.visited = false ;
    EdgePriority::EntryHandle new_handle = order_heap->InsertWithHandle(desc) ;
    edge_lengths->Insert(desc.length) ;

    edge_table->Store(EdgeKey(qdata->bl,qdata->br),new_handle) ;

    front_nodes->Store(qdata->bl,qdata->br) ;

    // update the edges before and after the edge

    bdry_topo->DeleteAngle(qdata->bl,qdata->tl,qdata->front_nodes[0]) ;
    bdry_topo->InsertAngle(1,qdata->bl,qdata->br,qdata->front_nodes[0]) ;

    bdry_topo->DeleteAngle(qdata->tl,qdata->tr,qdata->bl) ;
    bdry_topo->DeleteAngle(qdata->tr,qdata->br,qdata->tl) ;

    bdry_topo->DeleteAngle(qdata->br,qdata->front_nodes[5],qdata->tr) ;
    bdry_topo->InsertAngle(1,qdata->br,qdata->front_nodes[5],qdata->bl) ;

    // update the geometry

    // int fnodes[4] ;
    // fnodes[0] = qdata->front_nodes[0] ;
    // fnodes[1] = qdata->front_nodes[1] ;
    // fnodes[2] = qdata->front_nodes[4] ;
    // fnodes[2] = qdata->front_nodes[5] ;
    // UpdateGeomHelp(2,fnodes,level) ;
}



// %(MshEdgeList::UpdateCase2-int-|-NewQuadData-|*) 
/* ++ ----------------------------------------------------------
**
**    UpdateCase2 - update the edge list for a case 2 quad 
**
**      int UpdateCase2(NewQuadData *qdata)
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

int MshEdgeList::UpdateCase2(NewQuadData *qdata)
{
    int level ;

    // determine the level to assign to the new edges

    EdgePriority::EntryHandle *base0 = 
        edge_table->Get(EdgeKey(qdata->bl,qdata->br)) ;
    EdgePriority::EntryHandle *base1 =
        edge_table->Get(EdgeKey(qdata->br,qdata->tr)) ;

    IntQdEdgeDesc *desc_ptr0 = order_heap->ViewWithHandle(*base0) ;
    IntQdEdgeDesc *desc_ptr1 = order_heap->ViewWithHandle(*base1) ;

    level = ((desc_ptr0->level < desc_ptr1->level) ?
              desc_ptr0->level : desc_ptr1->level) + 1 ;

    double length0 = desc_ptr0->length ;
    double length1 = desc_ptr1->length ;

    // delete the base information from the tables

    order_heap->RemoveWithHandle(*base0) ;
    order_heap->RemoveWithHandle(*base1) ;
    edge_lengths->Remove(length0) ;
    edge_lengths->Remove(length1) ;
    edge_table->Del(EdgeKey(qdata->bl,qdata->br)) ;
    edge_table->Del(EdgeKey(qdata->br,qdata->tr)) ;
    front_nodes->Del(qdata->br) ;

    bdry_topo->DeleteAngle(qdata->br,qdata->tr,qdata->bl) ;
    
    edge_table->Store(EdgeKey(qdata->bl,qdata->tl),0) ;
    edge_table->Store(EdgeKey(qdata->tl,qdata->tr),0) ;

    front_nodes->Store(qdata->bl,qdata->tl) ;
    front_nodes->Store(qdata->tl,qdata->tr) ;

    // update the edges before and after the edge

    bdry_topo->DeleteAngle(qdata->bl,qdata->br,qdata->front_nodes[0]) ;
    bdry_topo->InsertAngle(1,qdata->bl,qdata->tl,qdata->front_nodes[0]) ;

    if (bdry_topo->HasVtx(qdata->tl)) {
        int prev_vtx, next_vtx ;
        FindAdjacent(qdata->tl,qdata->bl,&prev_vtx,&next_vtx) ;
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




// %(MshEdgeList::UpdateCase3-int-|-NewQuadData-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateCase3 - update the edge list for a case 3 quad 
**
**      int UpdateCase3(NewQuadData *qdata)
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

int MshEdgeList::UpdateCase3(NewQuadData *qdata)
{
    int level ;

    // determine the level to assign to the new edges

    EdgePriority::EntryHandle *base0 =
        edge_table->Get(EdgeKey(qdata->tl,qdata->bl)) ;
    EdgePriority::EntryHandle *base1 =
        edge_table->Get(EdgeKey(qdata->bl,qdata->br)) ;
    EdgePriority::EntryHandle *base2 =
        edge_table->Get(EdgeKey(qdata->br,qdata->tr)) ;

    IntQdEdgeDesc *desc_ptr0 = order_heap->ViewWithHandle(*base0) ;
    IntQdEdgeDesc *desc_ptr1 = order_heap->ViewWithHandle(*base1) ;
    IntQdEdgeDesc *desc_ptr2 = order_heap->ViewWithHandle(*base2) ;

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
    edge_table->Del(EdgeKey(qdata->tl,qdata->bl)) ;
    edge_table->Del(EdgeKey(qdata->bl,qdata->br)) ;
    edge_table->Del(EdgeKey(qdata->br,qdata->tr)) ;
    front_nodes->Del(qdata->bl) ;
    front_nodes->Del(qdata->br) ;

    bdry_topo->DeleteAngle(qdata->bl,qdata->br,qdata->tl) ;
    bdry_topo->DeleteAngle(qdata->br,qdata->tr,qdata->bl) ;
    
    // create the new top edge

    edge_table->Store(EdgeKey(qdata->tl,qdata->tr),0) ;

    front_nodes->Store(qdata->tl,qdata->tr) ;

    // update the edges before and after the edge

    bdry_topo->DeleteAngle(qdata->tl,qdata->bl,qdata->front_nodes[0]) ;
    bdry_topo->InsertAngle(1,qdata->tl,qdata->tr,qdata->front_nodes[0]) ;

    bdry_topo->DeleteAngle(qdata->tr,qdata->front_nodes[3],qdata->br) ;
    bdry_topo->InsertAngle(1,qdata->tr,qdata->front_nodes[3],qdata->tl) ;

    return(level) ;
}




// %(MshEdgeList::UpdateCase4-int-|-NewQuadData-|*) 
/* ++ ----------------------------------------------------------
**
**    UpdateCase4 - update the edge list for a case 4 quad 
**
**      int UpdateCase4(NewQuadData *qdata)
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

int MshEdgeList::UpdateCase4(NewQuadData *qdata)
{
    // determine the level to assign to the new edges

    EdgePriority::EntryHandle *base0 =
        edge_table->Get(EdgeKey(qdata->tl,qdata->bl)) ;
    EdgePriority::EntryHandle *base1 =
        edge_table->Get(EdgeKey(qdata->bl,qdata->br)) ;
    EdgePriority::EntryHandle *base2 =
        edge_table->Get(EdgeKey(qdata->br,qdata->tr)) ;
    EdgePriority::EntryHandle *base3 =
        edge_table->Get(EdgeKey(qdata->tr,qdata->tl)) ;

    IntQdEdgeDesc *desc_ptr0 = order_heap->ViewWithHandle(*base0) ;
    IntQdEdgeDesc *desc_ptr1 = order_heap->ViewWithHandle(*base1) ;
    IntQdEdgeDesc *desc_ptr2 = order_heap->ViewWithHandle(*base2) ;
    IntQdEdgeDesc *desc_ptr3 = order_heap->ViewWithHandle(*base3) ;

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
    edge_table->Del(EdgeKey(qdata->tl,qdata->bl)) ;
    edge_table->Del(EdgeKey(qdata->bl,qdata->br)) ;
    edge_table->Del(EdgeKey(qdata->br,qdata->tr)) ;
    edge_table->Del(EdgeKey(qdata->tr,qdata->tl)) ;
    front_nodes->Del(qdata->bl) ;
    front_nodes->Del(qdata->br) ;
    front_nodes->Del(qdata->tl) ;
    front_nodes->Del(qdata->tr) ;

    bdry_topo->DeleteAngle(qdata->bl,qdata->br,qdata->tl) ;
    bdry_topo->DeleteAngle(qdata->br,qdata->tr,qdata->bl) ;
    bdry_topo->DeleteAngle(qdata->tr,qdata->tl,qdata->br) ;
    bdry_topo->DeleteAngle(qdata->tl,qdata->bl,qdata->tr) ;

    return(0) ;
}




// %(MshEdgeList::UpdateCase5-int-|-NewQuadData-|*) 
/* ++ ----------------------------------------------------------
**
**    UpdateCase5 - update the edge list for a case 4 quad 
**
**      int UpdateCase5(NewQuadData *qdata)
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

int MshEdgeList::UpdateCase5(NewQuadData *qdata)
{
    int level ;

    // determine the level to assign to the new edges

    EdgePriority::EntryHandle *base0 = 
        edge_table->Get(EdgeKey(qdata->bl,qdata->br)) ;
    EdgePriority::EntryHandle *base1 =
        edge_table->Get(EdgeKey(qdata->tr,qdata->tl)) ;

    IntQdEdgeDesc *desc_ptr0 = order_heap->ViewWithHandle(*base0) ;
    IntQdEdgeDesc *desc_ptr1 = order_heap->ViewWithHandle(*base1) ;

    level = ((desc_ptr0->level < desc_ptr1->level) ?
              desc_ptr0->level : desc_ptr1->level) + 1 ;

    double length0 = desc_ptr0->length ;
    double length1 = desc_ptr1->length ;

    // delete the base information from the tables

    order_heap->RemoveWithHandle(*base0) ;
    order_heap->RemoveWithHandle(*base1) ;
    edge_lengths->Remove(length0) ;
    edge_lengths->Remove(length1) ;
    edge_table->Del(EdgeKey(qdata->bl,qdata->br)) ;
    edge_table->Del(EdgeKey(qdata->tr,qdata->tl)) ;

    edge_table->Store(EdgeKey(qdata->bl,qdata->tl),0) ;
    edge_table->Store(EdgeKey(qdata->tr,qdata->br),0) ;

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


// %(MshEdgeList::UpdateCase6-int-|-NewQuadData-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateCase6 - update the edge list for a case 6 quad
**                   (triangle)
**
**      int UpdateCase6(NewQuadData *qdata)
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

int MshEdgeList::UpdateCase6(NewQuadData *qdata)
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

    EdgePriority::EntryHandle *base0 =
        edge_table->Get(EdgeKey(qdata->bl,qdata->br)) ;
    EdgePriority::EntryHandle *base1 =
        edge_table->Get(EdgeKey(qdata->br,qdata->tr)) ;

    IntQdEdgeDesc *desc_ptr0 = order_heap->ViewWithHandle(*base0) ;
    IntQdEdgeDesc *desc_ptr1 = order_heap->ViewWithHandle(*base1) ;

    double length0 = desc_ptr0->length ;
    double length1 = desc_ptr1->length ;

    level = (desc_ptr0->level > desc_ptr1->level) ?
             desc_ptr0->level : desc_ptr1->level ;

    // delete the base information from the tables

    order_heap->RemoveWithHandle(*base0) ;
    order_heap->RemoveWithHandle(*base1) ;
    edge_lengths->Remove(length0) ;
    edge_lengths->Remove(length1) ;
    edge_table->Del(EdgeKey(qdata->bl,qdata->br)) ;
    edge_table->Del(EdgeKey(qdata->br,qdata->tr)) ;
    front_nodes->Del(qdata->br) ;

    bdry_topo->DeleteAngle(qdata->br,qdata->tr,qdata->bl) ;
    
    // create the new top edge

    edge_table->Store(EdgeKey(qdata->bl,qdata->tr),0) ;

    front_nodes->Store(qdata->bl,qdata->tr) ;

    // update the edges before and after the edge

    bdry_topo->DeleteAngle(qdata->bl,qdata->br,qdata->front_nodes[0]) ;
    bdry_topo->InsertAngle(1,qdata->bl,qdata->tr,qdata->front_nodes[0]) ;

    bdry_topo->DeleteAngle(qdata->tr,qdata->front_nodes[3],qdata->br) ;
    bdry_topo->InsertAngle(1,qdata->tr,qdata->front_nodes[3],qdata->bl) ;

    return(level) ;
}



// %(MshEdgeList::UpdateTriangle-void-|-int-const|-int-const|-int-const|)
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

void MshEdgeList::UpdateTriangle(const int id0,
                                     const int id1,
                                     const int id2)
{
    // delete the base information from the tables

    EdgePriority::EntryHandle *base0 = edge_table->Get(EdgeKey(id0,id1)) ;
    EdgePriority::EntryHandle *base1 = edge_table->Get(EdgeKey(id1,id2)) ;
    EdgePriority::EntryHandle *base2 = edge_table->Get(EdgeKey(id2,id0)) ;

    IntQdEdgeDesc *desc_ptr0 = order_heap->ViewWithHandle(*base0) ;
    IntQdEdgeDesc *desc_ptr1 = order_heap->ViewWithHandle(*base1) ;
    IntQdEdgeDesc *desc_ptr2 = order_heap->ViewWithHandle(*base2) ;

    double length0 = desc_ptr0->length ;
    double length1 = desc_ptr1->length ;
    double length2 = desc_ptr2->length ;

    order_heap->RemoveWithHandle(*base0) ;
    order_heap->RemoveWithHandle(*base1) ;
    order_heap->RemoveWithHandle(*base2) ;

    edge_lengths->Remove(length0) ;
    edge_lengths->Remove(length1) ;
    edge_lengths->Remove(length2) ;

    edge_table->Del(EdgeKey(id0,id1)) ;
    edge_table->Del(EdgeKey(id1,id2)) ;
    edge_table->Del(EdgeKey(id2,id0)) ;

    front_nodes->Del(id0) ;
    front_nodes->Del(id1) ;
    front_nodes->Del(id2) ;

    bdry_topo->DeleteTriangle(id0,id1,id2) ;
}




// %(MshEdgeList::UpdateSeam-int-|-NewSeamData-|*)
/* ++ ----------------------------------------------------------
**
**    UpdateSeam - update the front for a seam operation 
**
**      int UpdateSeam(NewSeamData *sdata)
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

int MshEdgeList::UpdateSeam(NewSeamData *sdata)
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


    IntNode *node_0 = node_table->Get(sdata->id0) ;
    IntNode *node_2 = node_table->Get(sdata->id2) ;
    int keep ;
    if (node_0->motion == MSH_FIXED) {
        keep = sdata->id0 ;
    } else if (node_2->motion == MSH_FIXED) {
        keep = sdata->id2 ;
    } else {
        keep = sdata->id0 ;
    }

    // now delete edges 0-1 and 1-2 from the edge list

    EdgePriority::EntryHandle *handle ;
    IntQdEdgeDesc *desc_ptr ;

    handle = edge_table->Get(EdgeKey(sdata->id0,sdata->id1)) ;
    desc_ptr = order_heap->ViewWithHandle(*handle) ;
    edge_lengths->Remove(desc_ptr->length) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_table->Del(EdgeKey(sdata->id0,sdata->id1)) ;

    handle = edge_table->Get(EdgeKey(sdata->id1,sdata->id2)) ;
    desc_ptr = order_heap->ViewWithHandle(*handle) ;
    edge_lengths->Remove(desc_ptr->length) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_table->Del(EdgeKey(sdata->id1,sdata->id2)) ;

    bdry_topo->DeleteAngle(sdata->id0,sdata->id1,sdata->cw) ;
    bdry_topo->DeleteAngle(sdata->id1,sdata->id2,sdata->id0) ;
    bdry_topo->DeleteAngle(sdata->id2,sdata->ccw,sdata->id1) ;
    front_nodes->Del(sdata->id1) ;

    if (sdata->cw == sdata->ccw) {

        // deal with the special case where cw == cww.  If this
        // is the case we have a "closing" seam, so we delete
        // all four edges.

        handle = edge_table->Get(EdgeKey(sdata->cw,sdata->id0)) ;
        desc_ptr = order_heap->ViewWithHandle(*handle) ;
        edge_lengths->Remove(desc_ptr->length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Del(EdgeKey(sdata->cw,sdata->id0)) ;

        handle = edge_table->Get(EdgeKey(sdata->id2,sdata->ccw)) ;
        desc_ptr = order_heap->ViewWithHandle(*handle) ;
        edge_lengths->Remove(desc_ptr->length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Del(EdgeKey(sdata->id2,sdata->ccw)) ;

        bdry_topo->DeleteAngle(sdata->cw,sdata->id0,sdata->id2) ;

        front_nodes->Del(sdata->id0) ;
        front_nodes->Del(sdata->id2) ;
        front_nodes->Del(sdata->cw) ;

        if (keep == sdata->id0) {

            // check to see if there are any edges attached to
            // id2 that need to be updated

            TopoAdjVtxIterator iter(bdry_topo,sdata->id2) ;

            while(iter.More()) {
                int prev = iter.AdjVtx() ;
                int prev_elem = iter.CcwElem() ;
                ++iter ;
                if ((prev_elem == 1) && (prev != sdata->ccw)) {
                    int next = iter.AdjVtx() ;

                    EdgePriority::EntryHandle *handle, value ;
                    handle = edge_table->Get(EdgeKey(sdata->id2,prev)) ;
                    value = *handle ;
                    desc_ptr = order_heap->ViewWithHandle(value) ;
                    desc_ptr->id_0 = sdata->id0 ;

                    edge_table->Del(EdgeKey(sdata->id2,prev)) ;
                    edge_table->Store(EdgeKey(sdata->id0,prev),value) ;

                    handle = edge_table->Get(EdgeKey(next,sdata->id2)) ;
                    value = *handle ;
                    desc_ptr = order_heap->ViewWithHandle(value) ;
                    desc_ptr->id_1 = sdata->id0 ;

                    edge_table->Del(EdgeKey(next,sdata->id2)) ;
                    edge_table->Store(EdgeKey(next,sdata->id0),value) ;

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

                    front_nodes->Del(sdata->id2) ;
                    front_nodes->Store(sdata->id0,prev) ;
                    front_nodes->Store(next,sdata->id0) ;

                    iter.NewVtx(sdata->id2) ;
                }
            }
        } else {

            // check to see if there are any edges attached to
            // id0 that need to be updated

            TopoAdjVtxIterator iter(bdry_topo,sdata->id0) ;

            while(iter.More()) {
                int prev = iter.AdjVtx() ;
                int prev_elem = iter.CcwElem() ;
                ++iter ;
                if ((prev_elem == 1) && (prev != sdata->cw)) {
                    int next = iter.AdjVtx() ;

                    EdgePriority::EntryHandle *handle, value ;
                    handle = edge_table->Get(EdgeKey(sdata->id0,prev)) ;
                    value = *handle ;
                    desc_ptr = order_heap->ViewWithHandle(value) ;
                    desc_ptr->id_0 = sdata->id2 ;

                    edge_table->Del(EdgeKey(sdata->id0,prev)) ;
                    edge_table->Store(EdgeKey(sdata->id2,prev),value) ;

                    handle = edge_table->Get(EdgeKey(next,sdata->id0)) ;
                    value = *handle ;
                    desc_ptr = order_heap->ViewWithHandle(value) ;
                    desc_ptr->id_1 = sdata->id2 ;

                    edge_table->Del(EdgeKey(next,sdata->id0)) ;
                    edge_table->Store(EdgeKey(next,sdata->id2),value) ;

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

                    front_nodes->Del(sdata->id0) ;
                    front_nodes->Store(sdata->id2,prev) ;
                    front_nodes->Store(next,sdata->id2) ;

                    iter.NewVtx(sdata->id0) ;
                }
            }
        }

    } else {

        // update the edges

        if (keep == sdata->id0) {
            bdry_topo->InsertAngle(1,sdata->id0,sdata->ccw,sdata->cw) ;

            EdgePriority::EntryHandle *handle, value ;
            handle = edge_table->Get(EdgeKey(sdata->id2,sdata->ccw)) ;
            value = *handle ;
            desc_ptr = order_heap->ViewWithHandle(value) ;
            desc_ptr->id_0 = sdata->id0 ;

            edge_table->Del(EdgeKey(sdata->id2,sdata->ccw)) ;
            edge_table->Store(EdgeKey(sdata->id0,sdata->ccw),value) ;

            front_nodes->Store(sdata->id0,sdata->ccw) ;

            int nccw = bdry_topo->GetCCWBdryNode(sdata->id2,sdata->ccw) ;
            bdry_topo->DeleteAngle(sdata->ccw,nccw,sdata->id2) ;
            bdry_topo->InsertAngle(1,sdata->ccw,nccw,sdata->id0) ;

            // check to see if there are any edges attached to
            // id2 that need to be updated

            TopoAdjVtxIterator iter(bdry_topo,sdata->id2) ;

            while(iter.More()) {
                int prev = iter.AdjVtx() ;
                int prev_elem = iter.CcwElem() ;
                ++iter ;
                if ((prev_elem == 1) && (prev != sdata->ccw)) {
                    int next = iter.AdjVtx() ;

                    handle = edge_table->Get(EdgeKey(sdata->id2,prev)) ;
                    value = *handle ;
                    desc_ptr = order_heap->ViewWithHandle(value) ;
                    desc_ptr->id_0 = sdata->id0 ;

                    edge_table->Del(EdgeKey(sdata->id2,prev)) ;
                    edge_table->Store(EdgeKey(sdata->id0,prev),value) ;

                    handle = edge_table->Get(EdgeKey(next,sdata->id2)) ;
                    value = *handle ;
                    desc_ptr = order_heap->ViewWithHandle(value) ;
                    desc_ptr->id_1 = sdata->id0 ;

                    edge_table->Del(EdgeKey(next,sdata->id2)) ;
                    edge_table->Store(EdgeKey(next,sdata->id0),value) ;

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

                    front_nodes->Del(sdata->id2) ;
                    front_nodes->Store(sdata->id0,prev) ;

                    iter.NewVtx(sdata->id2) ;
                }
            }
        } else {
            bdry_topo->InsertAngle(1,sdata->id2,sdata->ccw,sdata->cw) ;

            EdgePriority::EntryHandle *handle, value ;
            handle = edge_table->Get(EdgeKey(sdata->cw,sdata->id0)) ;
            value = *handle ;
            desc_ptr = order_heap->ViewWithHandle(value) ;
            desc_ptr->id_1 = sdata->id2 ;

            edge_table->Del(EdgeKey(sdata->cw,sdata->id0)) ;
            edge_table->Store(EdgeKey(sdata->cw,sdata->id2),value) ;

            front_nodes->Store(sdata->cw,sdata->id2) ;

            int ncw = bdry_topo->GetCWBdryNode(sdata->cw,sdata->id0) ;
            bdry_topo->DeleteAngle(sdata->cw,sdata->id0,ncw) ;
            bdry_topo->InsertAngle(1,sdata->cw,sdata->id2,ncw) ;

            sdata->id2 = sdata->id0 ;
            sdata->id0 = keep ;

            // check to see if there are any edges attached to
            // id0 that need to be updated

            TopoAdjVtxIterator iter(bdry_topo,sdata->id2) ;

            while(iter.More()) {
                int prev = iter.AdjVtx() ;
                int prev_elem = iter.CcwElem() ;
                ++iter ;
                if ((prev_elem == 1) && (prev != sdata->ccw)) {
                    int next = iter.AdjVtx() ;

                    handle = edge_table->Get(EdgeKey(sdata->id0,prev)) ;
                    value = *handle ;
                    desc_ptr = order_heap->ViewWithHandle(value) ;
                    desc_ptr->id_0 = sdata->id2 ;

                    edge_table->Del(EdgeKey(sdata->id0,prev)) ;
                    edge_table->Store(EdgeKey(sdata->id2,prev),value) ;

                    handle = edge_table->Get(EdgeKey(next,sdata->id0)) ;
                    value = *handle ;
                    desc_ptr = order_heap->ViewWithHandle(value) ;
                    desc_ptr->id_1 = sdata->id2 ;

                    edge_table->Del(EdgeKey(next,sdata->id0)) ;
                    edge_table->Store(EdgeKey(next,sdata->id2),value) ;

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

                    front_nodes->Del(sdata->id0) ;
                    front_nodes->Store(sdata->id2,prev) ;

                    iter.NewVtx(sdata->id2) ;
                }
            }
        }
    }
    return(keep) ;
}




// %(MshEdgeList::UpdateTransSeam-void-|-NewSeamData-const|*-int-const|-int-const|) 
/* ++ ----------------------------------------------------------
**
**    UpdateTransSeam - update the front for a transition seam 
**
**      void UpdateTransSeam(
**              const NewSeamData *sdata,
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

void MshEdgeList::UpdateTransSeam(
                 const NewSeamData *sdata,
                 const int mid_id,
                 const int far_id)
{
    if (sdata->ratio > 1.0) {

        // delete the existing edge

        EdgePriority::EntryHandle *handle, new_handle ;
        IntQdEdgeDesc desc ;

        handle = edge_table->Get(EdgeKey(sdata->id0,sdata->id1)) ;
        desc = *(order_heap->ViewWithHandle(*handle)) ;
        double edge_length = desc.length ;
        edge_lengths->Remove(edge_length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Del(EdgeKey(sdata->id0,sdata->id1)) ;

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
        edge_table->Store(EdgeKey(sdata->id0,mid_id),new_handle) ;

        desc.length = edge_length ;
        desc.id_0 = mid_id ;
        desc.id_1 = far_id ;
        desc.angle_0 = desc.angle_1 = HALF_PI ;
        desc.end_code = 0 ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(EdgeKey(mid_id,far_id),new_handle) ;

        desc.id_0 = far_id ;
        desc.id_1 = sdata->id1 ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(EdgeKey(far_id,sdata->id1),new_handle) ;

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

        EdgePriority::EntryHandle *handle, new_handle ;
        IntQdEdgeDesc desc ;

        handle = edge_table->Get(EdgeKey(sdata->id1,sdata->id2)) ;
        desc = *(order_heap->ViewWithHandle(*handle)) ;
        double edge_length = desc.length ;
        edge_lengths->Remove(edge_length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Del(EdgeKey(sdata->id1,sdata->id2)) ;

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
        edge_table->Store(EdgeKey(mid_id,sdata->id2),new_handle) ;

        desc.length = edge_length ;
        desc.id_0 = far_id ;
        desc.id_1 = mid_id ;
        desc.angle_0 = desc.angle_1 = HALF_PI ;
        desc.end_code = 0 ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(EdgeKey(far_id,mid_id),new_handle) ;

        desc.id_0 = sdata->id1 ;
        desc.id_1 = far_id ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(EdgeKey(sdata->id1,far_id),new_handle) ;

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




// %(MshEdgeList::UpdateTransSplit-void-|-NewQuadData-const|*-int-const|-int-const|) 
/* ++ ----------------------------------------------------------
**
**    UpdateTransSplit - update the front for a transition split 
**
**      void UpdateTransSplit(
**              const NewQuadData *qdata,
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

void MshEdgeList::UpdateTransSplit(
                 const NewQuadData *qdata,
                 const int mid_id0,
                 const int mid_id1)
{
    if (qdata->ratio > 1.0) {

        // delete the existing edge

        EdgePriority::EntryHandle *handle, new_handle ;
        IntQdEdgeDesc desc ;

        handle = edge_table->Get(EdgeKey(qdata->bl,qdata->br)) ;
        desc = *(order_heap->ViewWithHandle(*handle)) ;
        double edge_length = desc.length ;
        edge_lengths->Remove(edge_length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Del(EdgeKey(qdata->bl,qdata->br)) ;

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
        edge_table->Store(EdgeKey(qdata->bl,mid_id0),new_handle) ;

        desc.id_0 = mid_id0 ;
        desc.id_1 = mid_id1 ;
        desc.angle_0 = desc.angle_1 = HALF_PI ;
        desc.end_code = 0 ;
        desc.visited = true ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(EdgeKey(mid_id0,mid_id1),new_handle) ;

        desc.id_0 = mid_id1 ;
        desc.id_1 = qdata->br ;
        desc.visited = true ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(EdgeKey(mid_id1,qdata->br),new_handle) ;

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

        EdgePriority::EntryHandle *handle, new_handle ;
        IntQdEdgeDesc desc ;

        handle = edge_table->Get(EdgeKey(qdata->br,qdata->tr)) ;
        desc = *(order_heap->ViewWithHandle(*handle)) ;
        double edge_length = desc.length ;
        edge_lengths->Remove(edge_length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Del(EdgeKey(qdata->br,qdata->tr)) ;

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
        edge_table->Store(EdgeKey(mid_id0,qdata->tr),new_handle) ;

        desc.id_0 = mid_id1 ;
        desc.id_1 = mid_id0 ;
        desc.angle_0 = desc.angle_1 = HALF_PI ;
        desc.end_code = 0 ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(EdgeKey(mid_id1,mid_id0),new_handle) ;

        desc.id_0 = qdata->br ;
        desc.id_1 = mid_id1 ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(EdgeKey(qdata->br,mid_id1),new_handle) ;

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




// %(MshEdgeList::UpdateTempSplit-void-|-NewQuadData-const|*-int-const|-int-const|-bool-const|) 
/* ++ ----------------------------------------------------------
**
**    UpdateTempSplit - update the front for a template split 
**
**      void UpdateTempSplit(
**              const NewQuadData *qdata,
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

void MshEdgeList::UpdateTempSplit(
                 const NewQuadData *qdata,
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
        EdgePriority::EntryHandle *handle, new_handle ;
        IntQdEdgeDesc desc, orig_desc ;

        handle = edge_table->Get(EdgeKey(qdata->tl,qdata->tr)) ;
        orig_desc = *(order_heap->ViewWithHandle(*handle)) ;
        double edge_length = orig_desc.length ;
        edge_lengths->Remove(edge_length) ;
        order_heap->RemoveWithHandle(*handle) ;
        edge_table->Del(EdgeKey(qdata->tl,qdata->tr)) ;

        // now add the information for the new edges.

        desc = orig_desc ;
        desc.length = edge_length / 3.0 ;
        desc.angle_1 = PI ;
        desc.id_1 = mid_id1 ;
        desc.end_code &= LEFT_EDGE_MASK ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(EdgeKey(qdata->tl,mid_id1),new_handle) ;

        desc.angle_0 = PI ;
        desc.id_0 = mid_id1 ;
        desc.id_1 = mid_id0 ;
        desc.end_code = 0 ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(EdgeKey(mid_id1,mid_id0),new_handle) ;

        desc.angle_1 = orig_desc.angle_1 ;
        desc.id_0 = mid_id0 ;
        desc.id_1 = qdata->tr ;
        desc.end_code = orig_desc.end_code & RIGHT_EDGE_MASK ;
        desc.visited = false ;
        new_handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;
        edge_table->Store(EdgeKey(mid_id0,qdata->tr),new_handle) ;
    } else {
        edge_table->Del(EdgeKey(qdata->tl,qdata->tr)) ;

        edge_table->Store(EdgeKey(qdata->tl,mid_id1),0) ;
        edge_table->Store(EdgeKey(mid_id1,mid_id0),0) ;
        edge_table->Store(EdgeKey(mid_id0,qdata->tr),0) ;
    }

    front_nodes->Store(qdata->tl,mid_id1) ;
    front_nodes->Store(mid_id1,mid_id0) ;
    front_nodes->Store(mid_id0,qdata->tr) ;
}




// %(MshEdgeList::UpdateGeom-void-|-NewQuadData-|*) 
/* ++ ----------------------------------------------------------
**
**    UpdateGeom - update stored geometry parameters 
**
**      void UpdateGeom(NewQuadData *qdata)
**
**        qdata - (in)  quadrilateral description 
**
**      Description: This function updates the stored mesh front edge 
**          lengths and the angles associated with a new quad. This 
**          function is called after nodal smoothing. 
**
**
** -- */

void MshEdgeList::UpdateGeom(NewQuadData *qdata)
{
    if ((qdata->qcase < 5) || (qdata->qcase == 6)) {
        UpdateGeomHelp(qdata->num_front_nodes,qdata->front_nodes,
                       qdata->edge_level) ;
    } else {
        UpdateGeomHelp(4,&qdata->front_nodes[0],qdata->edge_level) ;
        UpdateGeomHelp(4,&qdata->front_nodes[4],qdata->edge_level) ;
    }

    if ((qdata->qcase == 3) || (qdata->qcase == 4))  {
        if (qdata->bl != qdata->obl) {
            int front_nodes[3] ;
            if (UpdateGeomOne(qdata->bl,front_nodes))
                UpdateGeomHelp(3,front_nodes,qdata->edge_level) ;
        }
    }
    if ((qdata->qcase >= 2) || (qdata->qcase <= 4)) {
        if (qdata->br != qdata->obr) {
            int front_nodes[3] ;
            if (UpdateGeomOne(qdata->br,front_nodes))
                UpdateGeomHelp(3,front_nodes,qdata->edge_level) ;
        }
    }
}




// %(MshEdgeList::UpdateGeomHelp-void-|-int-const|-int-const|*-int-const|)
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

static int EndCode(double *angle,int level,
                   int adj_level,bool fixed_flag=false)
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

void MshEdgeList::UpdateGeomHelp(const int num_front_nodes,
                                     const int *front_nodes,
                                     const int edge_level)
{
   IntNode *nd_cw, *nd_ccw, *nd_0, *nd_1 ;
    int vtx_cw, vtx_ccw, vtx_0, vtx_1 ;
    IntNode *adj_nd_0, *adj_nd_1 ;
    int adj_vtx_0, adj_vtx_1 ;
    double dx, dy, angle ;
    EdgePriority::EntryHandle *handle ;
    short adj_level, end_code ;
    IntQdEdgeDesc desc, adj_desc, *desc_ptr ;

    // save edge level of the ccw edge

    handle = edge_table->Get(EdgeKey(
               front_nodes[num_front_nodes-2],
               front_nodes[num_front_nodes-1])) ;
    desc_ptr = order_heap->ViewWithHandle(*handle) ;

    // first look at the edge that is twice clockwise from
    // where we are adding the new element.  For this edge
    // we update the angle on the right end.

    vtx_0 = front_nodes[0] ;
    vtx_1 = front_nodes[1] ;
    vtx_cw = bdry_topo->GetCWBdryNode(vtx_0,vtx_1) ;

    nd_0 = node_table->Get(vtx_0) ;
    nd_1 = node_table->Get(vtx_1) ;
    nd_cw = node_table->Get(vtx_cw) ;

    handle = edge_table->Get(EdgeKey(vtx_cw,vtx_0)) ;
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

    handle = edge_table->Get(EdgeKey(vtx_0,vtx_1)) ;
    desc = *(order_heap->ViewWithHandle(*handle)) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_lengths->Remove(desc.length) ;

    desc.angle_0 = angle ;
    desc.end_code &= RIGHT_EDGE_MASK ;
    desc.end_code |= end_code << 1 ;

    nd_ccw = node_table->Get(vtx_ccw) ;

    dx = nd_0->coord[0] - nd_1->coord[0] ;
    dy = nd_0->coord[1] - nd_1->coord[1] ;
    desc.length = sqrt(dx*dx + dy*dy) ;
    angle = Angle2Pi(nd_1->coord,nd_ccw->coord,nd_0->coord) ;

    desc.end_code &= LEFT_EDGE_MASK ;
    bool fix_flg = (nd_1->motion == MSH_FIXED) &&
                   (nd_0->motion == MSH_FIXED) &&
                   (nd_ccw->motion == MSH_FIXED) ;
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
        nd_1 = node_table->Get(vtx_1) ;

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
            nd_cw = node_table->Get(vtx_cw) ;
            angle = Angle2Pi(nd_0->coord,nd_1->coord,nd_cw->coord) ;

            // update the adjacent edge

            adj_vtx_0 = vtx_cw ;
            adj_vtx_1 = vtx_0 ;

            adj_nd_0 = node_table->Get(adj_vtx_0) ;
            adj_nd_1 = nd_0 ;

            handle = edge_table->Get(EdgeKey(adj_vtx_0,adj_vtx_1)) ;
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
        nd_ccw = node_table->Get(vtx_ccw) ;

        angle = Angle2Pi(nd_1->coord,nd_ccw->coord,nd_0->coord) ;
        bool fix_flg = (nd_1->motion == MSH_FIXED) &&
                       (nd_0->motion == MSH_FIXED) &&
                       (nd_ccw->motion == MSH_FIXED) ;
        end_code = EndCode(&angle,desc.level,desc.adj_level_1,fix_flg) ;
        desc.angle_1 = angle ;
        desc.end_code |= end_code ;
        desc.visited = false ;

        handle = edge_table->Get(EdgeKey(vtx_0,vtx_1)) ;
        *handle = order_heap->InsertWithHandle(desc) ;
        edge_lengths->Insert(desc.length) ;

        // if necessary, update the adjacent edge on the right

        if (vtx_ccw != front_nodes[i+3]) {
            adj_vtx_0 = vtx_1 ;
            adj_vtx_1 = vtx_ccw ;

            adj_nd_0 = nd_1 ;
            adj_nd_1 = node_table->Get(adj_vtx_1) ;

            handle = edge_table->Get(EdgeKey(adj_vtx_0,adj_vtx_1)) ;
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
    nd_1 = node_table->Get(vtx_1) ;

    handle = edge_table->Get(EdgeKey(vtx_0,vtx_1)) ;
    desc = *(order_heap->ViewWithHandle(*handle)) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_lengths->Remove(desc.length) ;

    dx = nd_0->coord[0] - nd_1->coord[0] ;
    dy = nd_0->coord[1] - nd_1->coord[1] ;
    desc.length = sqrt(dx*dx + dy*dy) ;

    if (vtx_ccw != vtx_1) {
        vtx_cw = bdry_topo->GetCWBdryNode(vtx_0,vtx_1) ;
        nd_cw = node_table->Get(vtx_cw) ;
        angle = Angle2Pi(nd_0->coord,nd_1->coord,nd_cw->coord) ;
        end_code = EndCode(&angle,desc.level,desc.adj_level_0) ;
    }

    desc.angle_0 = angle ;
    desc.adj_level_0 = edge_level ;
    desc.end_code &= RIGHT_EDGE_MASK ;
    desc.end_code |= end_code << 1 ;

    // now update the angle on the right end

    vtx_ccw = bdry_topo->GetCCWBdryNode(vtx_0,vtx_1) ;
    nd_ccw = node_table->Get(vtx_ccw) ;

    angle = Angle2Pi(nd_1->coord,nd_ccw->coord,nd_0->coord) ;
    end_code = EndCode(&angle,desc.level,desc.adj_level_1) ;
    desc.angle_1 = angle ;
    desc.end_code &= LEFT_EDGE_MASK ;
    desc.end_code |= end_code ;
    desc.visited = false ;

    *handle = order_heap->InsertWithHandle(desc) ;
    edge_lengths->Insert(desc.length) ;

    // now we update the angle for the next ccw edge

    handle = edge_table->Get(EdgeKey(vtx_1,vtx_ccw)) ;
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
}


bool MshEdgeList::UpdateGeomOne(const int node_id,
                                    int front_nodes[])
{
    // if this node is not on the boundary then return

    IntNode *node = node_table->Get(node_id) ;
    if (node == 0) return(false) ;

    TopoAdjVtxIterator iter(bdry_topo,node_id) ;
    if (iter.More()) {
        front_nodes[1] = node_id ;
        front_nodes[2] = iter.AdjVtx() ;
        ++iter ;
        front_nodes[0] = iter.AdjVtx() ;
        return(true) ;
    }
    return(false) ;
}


// %(MshEdgeList::UpdateBdryGeom-void-|-int-const|) 
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

void MshEdgeList::UpdateBdryGeom(const int vtx)
{
    IntNode *nd_cw, *nd_ccw, *nd_0, *nd_1 ;
    int vtx_cw, vtx_ccw, vtx_0, vtx_1 ;
    double dx, dy, angle ;
    EdgePriority::EntryHandle *handle ;
    short end_code ;
    IntQdEdgeDesc desc ;

    // first look at the edge that is counter clockwise from
    // the vertex.

    vtx_0 = vtx ;

    TopoAdjVtxIterator iter(bdry_topo,vtx) ;

    if (!iter.More()) return ;
    while (iter.More() && (iter.CcwElem() != 1)) ++iter ;

    vtx_1 = iter.AdjVtx() ;
    vtx_ccw = bdry_topo->GetCCWBdryNode(vtx_0,vtx_1) ;

    handle = edge_table->Get(EdgeKey(vtx_0,vtx_1)) ;
    desc = *(order_heap->ViewWithHandle(*handle)) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_lengths->Remove(desc.length) ;

    nd_0 = node_table->Get(vtx_0) ;
    nd_1 = node_table->Get(vtx_1) ;
    nd_ccw = node_table->Get(vtx_ccw) ;

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
    nd_cw = node_table->Get(vtx_cw) ;
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

    handle = edge_table->Get(EdgeKey(vtx_0,vtx_1)) ;
    desc = *(order_heap->ViewWithHandle(*handle)) ;
    order_heap->RemoveWithHandle(*handle) ;
    edge_lengths->Remove(desc.length) ;

    nd_1 = nd_0 ;
    nd_0 = nd_cw ;
    nd_cw = node_table->Get(vtx_cw) ;

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




// %(MshEdgeList::UpdateSeamGeom-void-|-NewSeamData-const|*)
/* ++ ----------------------------------------------------------
**
**    UpdateSeamGeom - update stored geometry parameters 
**
**      void UpdateSeamGeom(const NewSeamData *sdata)
**
**        sdata - (in)  seam data 
**
**      Description: This function updates the stored mesh front 
**          lengths and angle. This function is called after nodal 
**          smoothing. 
**
**
** -- */

void MshEdgeList::UpdateSeamGeom(const NewSeamData *sdata)
{
    int front_nodes[3] ;

    front_nodes[0] = sdata->cw ;
    front_nodes[1] = sdata->id0 ;
    front_nodes[2] = sdata->ccw ;
    UpdateGeomHelp(3,front_nodes,1) ;
}




// %(MshEdgeList::GetCWNode-int-|^const-int-const|-int-const|)
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

int MshEdgeList::GetCWNode(const int id_0,
                               const int id_1) const
{
    return(bdry_topo->GetCWBdryNode(id_0,id_1)) ;
}




// %(MshEdgeList::GetCCWNode-int-|^const-int-const|-int-const|)
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

int MshEdgeList::GetCCWNode(const int id_0,
                                const int id_1) const
{
    return(bdry_topo->GetCCWBdryNode(id_0,id_1)) ;
}




// %(MshEdgeList::GetEdgeLength-double-|^const-int-const|-int-const|)
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

double MshEdgeList::GetEdgeLength(const int id_0,
                                      const int id_1) const
{
    EdgePriority::EntryHandle *handle ;
    IntQdEdgeDesc *desc_ptr ;
    handle = edge_table->Get(EdgeKey(id_0,id_1)) ;
    desc_ptr = order_heap->ViewWithHandle(*handle) ;
    return(desc_ptr->length) ;
}




// %(MshEdgeList::GetEdgeLevel-int-|^const-int-const|-int-const|) 
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

int MshEdgeList::GetEdgeLevel(const int id_0,
                                           const int id_1) const
{
    EdgePriority::EntryHandle *handle ;
    IntQdEdgeDesc *desc_ptr ;
    handle = edge_table->Get(EdgeKey(id_0,id_1)) ;
    desc_ptr = order_heap->ViewWithHandle(*handle) ;
    return(desc_ptr->level) ;
}




// %(MshEdgeList::ContainsEdge-bool-|^const-int-const|-int-const|)
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

bool MshEdgeList::ContainsEdge(const int id_0,
                                   const int id_1) const
{
    EdgePriority::EntryHandle *tmp = edge_table->Get(EdgeKey(id_0,id_1)) ;
    return((tmp != 0) ? true : false) ;
}




// %(MshEdgeList::ContainsNode-bool-|^const-int-const|)
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

bool MshEdgeList::ContainsNode(const int id) const
{
    return((front_nodes->Get(id) != 0) ? true : false) ;
}




// %(MshEdgeList::GetNextEdge-QdEdge-|*)
/* ++ ----------------------------------------------------------
**
**    GetNextEdge - get next edge to use to advance the front 
**
**      QdEdge *GetNextEdge()
**
**      Description: This function returns information about the edge 
**          that is ranked highest in the edge list's priority queue 
**          for use in advancing the mesh front. 
**
**      Return Value: An QdEdge structure containing information the 
**          highest ranked edge. 
**
**
** -- */

MshEdgeList::QdEdge *MshEdgeList::GetNextEdge()
{
    static QdEdge edge ;

    IntQdEdgeDesc *next = order_heap->ViewMin() ;
    if (next == 0) return(0) ;
    next->visited = true ;
    edge.id_0 = next->id_0 ;
    edge.id_1 = next->id_1 ;
    edge.end_code = next->end_code ;
    edge.angle_0 = next->angle_0 ;
    edge.angle_1 = next->angle_1 ;
    edge.level = next->level ;

// For Debug

//     CArbHeap<IntQdEdgeDesc> tmp_heap = *order_heap ;
//     IntQdEdgeDesc *tmp = tmp_heap.GetMin() ;
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




// %(MshEdgeList::GetEdgeLengthRatio-double-|^const)
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

double MshEdgeList::GetEdgeLengthRatio() const
{
    double small, big ;
    small = *(edge_lengths->GetSmallest()) ;
    big = *(edge_lengths->GetLargest()) ;
    return(big/small) ;
}




// %(MshEdgeList::PushToBack-void-|-QdEdge-const|*) 
/* ++ ----------------------------------------------------------
**
**    PushToBack - put an edge back on the edge list 
**
**      void PushToBack(const QdEdge *edge)
**
**        edge - (in)  edge to place back on the list 
**
**      Description: This function places an edge back on the edge list 
**          priority queue, and flags it so that it will be ranked 
**          lowest in the queue. 
**
**
** -- */

void MshEdgeList::PushToBack(const QdEdge *edge)
{
    EdgePriority::EntryHandle *tedge =
        edge_table->Get(EdgeKey(edge->id_0,edge->id_1)) ;

    if (tedge == 0) return ;

    IntQdEdgeDesc desc = *(order_heap->ViewWithHandle(*tedge)) ;

    order_heap->RemoveWithHandle(*tedge) ;
    desc.level++ ;
    *tedge = order_heap->InsertWithHandle(desc) ;
}


bool MshEdgeList::StagnationSweep()
{
    // get a list of the entries in the heap

    IntQdEdgeDesc **list = order_heap->GetEntryList() ;

    // look for an unvisited edge

    bool stagnant = true ;
    for (int i=0 ; i<order_heap->Len() ; ++i) {
         if (!list[i]->visited) stagnant = false ;
         list[i]->visited = false ;
    }
    delete [] list ;
    return(stagnant) ;
}

double MshEdgeList::GetCharNodeLength(const int id_0) const
{
    int *ptr = front_nodes->Get(id_0) ;
    if (ptr == 0) return(0.0) ;
    int id_1 = *ptr ;

    if (!ContainsEdge(id_0,id_1)) return(1.0) ;

    IntNode *node_0 = node_table->Get(id_0) ;
    IntNode *node_1 = node_table->Get(id_1) ;

    double dist_0 = (node_0->coord - node_1->coord).Magnitude() ;

    int id_2 = GetCWNode(id_0,id_1) ;
    IntNode *node_2 = node_table->Get(id_2) ;

    double dist_1 = (node_0->coord - node_2->coord).Magnitude() ;

    return(0.5*(dist_0+dist_1)) ;
}


// %(MshEdgeList::Angle-double-|^const-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    Angle - compute the included angle 
**
**      double Angle(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j) const
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

double MshEdgeList::Angle(Vec2D b,Vec2D i,
                              Vec2D j) const
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
    return(acos(tmp) * (cross/fabs(cross))) ;
}




// %(MshEdgeList::Angle2Pi-double-|^const-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    Angle2Pi - compute the included angle 
**
**      double Angle2Pi(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j) const
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

double MshEdgeList::Angle2Pi(Vec2D b,Vec2D i,
                                 Vec2D j) const
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




// %(MshEdgeList::CrossProd-double-|^const-Vec2D-|-Vec2D-|-Vec2D-|)
/* ++ ----------------------------------------------------------
**
**    CrossProd - cross product of two vectors 
**
**      double CrossProd(
**              Vec2D b,
**              Vec2D i,
**              Vec2D j) const
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

double MshEdgeList::CrossProd(Vec2D b,Vec2D i,
                                  Vec2D j) const
{
    double cross ;

    cross = ((i[0] - b[0]) * (j[1] - b[1])) -
            ((i[1] - b[1]) * (j[0] - b[0])) ;
    return cross ;
}

} // namespace
