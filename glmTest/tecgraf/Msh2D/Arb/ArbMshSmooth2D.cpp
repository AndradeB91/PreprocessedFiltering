//
// ArbMshSmooth2D implementation file
//
// Description -
//   This is the implementation file for the
//   ArbMshSmooth2D objects.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2001
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.19 $  $Date: 2002/07/16 14:18:39 $  $Author: wash $
//

#include "ArbMshSmooth2D.hpp"
#include "ArbFeasableRegion.hpp"
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifdef MEMDEBUG
#include "MemDbg.hpp"
#define new new(__FILE__,__LINE__)
#endif

#define PI      3.14159265359
/*
#define HALF_PI 1.57079632680
#define TWO_PI  6.28318530718

#define NO_NODE 100000000
*/

#define DEFAULT_IT_TOL 0.01
#define DEFAULT_MAX_IT 25



// %(CArbMshSmooth2D::CArbMshSmooth2D-constructor-|)
/* ++ ----------------------------------------------------------
**
**    CArbMshSmooth2D - no argument constructor 
**
**      CArbMshSmooth2D()
**
**      Description: This is a no argument constructor for a 
**          CArbMshSmooth2D object. 
**
**
** -- */

CArbMshSmooth2D::CArbMshSmooth2D()
{
    NodeTable = new NodeHash ;
    ElemTable = new ElemHash ;
    ZCoords   = new ZCoordHash ;

    FirstElem    = true ;
    DeleteTables = true ;
    FindBound    = true ;
    Order        = LINEAR ;
    It_tol       = DEFAULT_IT_TOL ;
    Max_it       = DEFAULT_MAX_IT ;
}




// %(CArbMshSmooth2D::CArbMshSmooth2D-constructor-|-NodeHash-|*-ElemHash-|*)
/* ++ ----------------------------------------------------------
**
**    CArbMshSmooth2D - constructor 
**
**      CArbMshSmooth2D(
**              NodeHash *node_tbl,
**              ElemHash *elem_tbl)
**
**        node_tbl - (in)  node table 
**        elem_tbl - (in)  element table 
**
**      Description: This is a constructor for a CArbMshSmooth2D object 
**          with the node and element tables passed as arguments. 
**
**
** -- */
 
CArbMshSmooth2D::CArbMshSmooth2D(NodeHash *node_tbl,ElemHash *elem_tbl)
{
    NodeTable    = node_tbl ;
    ElemTable    = elem_tbl ;
    ZCoords      = 0 ;
    FirstElem    = false ;
    DeleteTables = false ;
    FindBound    = false ;
    It_tol       = DEFAULT_IT_TOL ;
    Max_it       = DEFAULT_MAX_IT ;

    ElemIter iter(ElemTable) ;
    ArbMshElement2D *elem = iter.Entry() ;
    if ((elem->num_nodes == 6) || (elem->num_nodes == 8))
        Order = QUADRATIC ;
    else
        Order = LINEAR ;
}




// %(CArbMshSmooth2D::CArbMshSmooth2D-destructor-|~)
/* ++ ----------------------------------------------------------
**
**    CArbMshSmooth2D - destructor 
**
**      ~CArbMshSmooth2D()
**
**      Description: This is a no destructor for a CArbMshSmooth2D 
**          object. 
**
**
** -- */

CArbMshSmooth2D::~CArbMshSmooth2D()
{
    if (DeleteTables) {
        delete NodeTable ;
        delete ElemTable ;
        delete ZCoords ;
    }
}




// %(CArbMshSmooth2D::AddNode-void-|-ArbMshNode-const|&)
/* ++ ----------------------------------------------------------
**
**    AddNode - add a node to the mesh 
**
**      void AddNode(const ArbMshNode &inode)
**
**        inode - (in)  description of the node 
**
**      Description: This function adds a new node to the mesh. 
**
**      Exceptions:
**          CArbMshCDuplNode - duplicate node id's
**
** -- */

int CArbMshSmooth2D::AddNode(const ArbMshNode &inode)
{
    ArbIntNode tmp ;
    tmp.id = inode.id ;
    tmp.coord = CArbCoord2D(inode.coord[0],inode.coord[1]) ;
    tmp.motion = ARB_FLOATING ;
    tmp.corner = false ;
    if (NodeTable->Store(inode.id,tmp) == false) {
        return(ARB_DUPLICATE_NODE_ID) ;
    }
    ZCoords->Store(inode.id,inode.coord[2]) ;
    FindBound = true ;
    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshSmooth2D::AddNodeList-void-|-int-const|-ArbMshNode-const|*const)
/* ++ ----------------------------------------------------------
**
**    AddNodeList - add an array of nodes to the mesh 
**
**      void AddNodeList(
**              const int        num_nodes,
**              const ArbMshNode *constnode_list)
**
**        num_nodes - (in)  number of nodes 
**        node_list - (in)  list of nodes to add 
**
**      Description: This function adds an array of nodes to the mesh. 
**
**      Exceptions:
**          CArbMshCDuplNode - duplicate node id's
**
** -- */

int CArbMshSmooth2D::AddNodeList(const int num_nodes,
                         const ArbMshNode *const node_list) 
{
    for (int i=0 ; i<num_nodes ; ++i) {
        int status = AddNode(node_list[i]) ;
        if (status != ARB_NORMAL_STATUS) return(status) ;
    }
    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshSmooth2D::AddElem-void-|-ArbMshElement2D-const|&)
/* ++ ----------------------------------------------------------
**
**    AddElem - add an element to the mesh 
**
**      void AddElem(const ArbMshElement2D &elem)
**
**        elem - (in)  element to add 
**
**      Description: This function adds an element to the mesh. 
**
**      Exceptions:
**          CArbMshCDuplElem - duplicate element id's
**          CArbMshCInvalidElem - invalid element
**
** -- */

int CArbMshSmooth2D::AddElem(const ArbMshElement2D &ielem)
{
    if (FirstElem) {
        if ((ielem.num_nodes == 3) || (ielem.num_nodes == 4)) {
            Order = LINEAR ;
        } else {
            Order = QUADRATIC ;
        }
        FirstElem = false ;
    } else {
        if ((ielem.num_nodes == 3) || (ielem.num_nodes == 4)) {
            if (Order != LINEAR) return(ARB_INVALID_ELEM) ;
        } else if ((ielem.num_nodes == 6) || (ielem.num_nodes == 8)) {
            if (Order != QUADRATIC) return(ARB_INVALID_ELEM) ;
        } else {
            return(ARB_INVALID_ELEM) ;
        }
    }

    // try to insert this element into the hash table

    if (ElemTable->Store(ielem.elem_id,ielem) == false) {
        return(ARB_DUPLICATE_ELEM_ID) ;
    }
    FindBound = true ;
    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshSmooth2D::AddElemList-void-|-int-const|-ArbMshElement2D-const|*const)
/* ++ ----------------------------------------------------------
**
**    AddElemList - add an array of elements to the mesh 
**
**      void AddElemList(
**              const int             num_elems,
**              const ArbMshElement2D *constelem_list)
**
**        num_elems - (in)  number of elements 
**        elem_list - (in)  elements 
**
**      Description: This function adds an array of elements to the 
**          mesh. 
**
**      Exceptions:
**          CArbMshCDuplElem - duplicate element id's
**          CArbMshCInvalidElem - invalid element
**
** -- */

int CArbMshSmooth2D::AddElemList(
               const int num_elems,
               const ArbMshElement2D *const elem_list) 
{
    for (int i=0 ; i<num_elems ; ++i) {
        int status = AddElem(elem_list[i]) ;
        if (status != ARB_NORMAL_STATUS) return(status) ;
    }
    return(ARB_NORMAL_STATUS) ;
}




// %(CArbMshSmooth2D::SmoothNodesLaplace-void-|)
/* ++ ----------------------------------------------------------
**
**    SmoothNodesLaplace - smooth nodes using a Laplace algorithm 
**
**      void SmoothNodesLaplace()
**
**      Description: This function smooths the internal nodes using a 
**          Laplace algorithm. 
**
**
** -- */

void CArbMshSmooth2D::SmoothNodesLaplace()
{
    if (FindBound) {
        FindBoundaryNodes() ;
        FindBound = false ;
    }

    CArbMshTopo2D *msh_topo = BuildMeshTopo() ;
    double min_length = FindMinCharEdgeLength(msh_topo) ;
    DoSmoothNodesLaplace(msh_topo,Max_it,min_length*It_tol) ;

    delete msh_topo ;
}




// %(CArbMshSmooth2D::SmoothNodesWinslow-void-|)
/* ++ ----------------------------------------------------------
**
**    SmoothNodesWinslow - smooth nodes using a Winslow algorithm 
**
**      void SmoothNodesWinslow()
**
**      Description: This function smooths the internal nodes using a 
**          Winslow algorithm. 
**
**
** -- */

int CArbMshSmooth2D::SmoothNodesWinslow()
{
    if (FindBound) {
        FindBoundaryNodes() ;
        FindBound = false ;
    }
    CArbMshTopo2D *msh_topo = BuildMeshTopo() ;
    double min_length = FindMinCharEdgeLength(msh_topo) ;
    if (min_length == -1)
     return 0;
    DoSmoothNodesLaplace(msh_topo,1,min_length*It_tol,true) ;
    DoSmoothNodesWinslow(msh_topo,Max_it,min_length*It_tol) ;

    delete msh_topo ;
    return 1;
}




// %(CArbMshSmooth2D::SmoothNodesConsLaplace-void-|)
/* ++ ----------------------------------------------------------
**
**    SmoothNodesConsLaplace - smooth nodes using a "constrained" 
**                             Laplace algorithm 
**
**      void SmoothNodesConsLaplace()
**
**      Description: This function smooths the internal nodes using a 
**          "constrained" Laplace algorithm. 
**
**
** -- */

void CArbMshSmooth2D::SmoothNodesConsLaplace()
{
    if (FindBound) {
        FindBoundaryNodes() ;
        FindBound = false ;
    }

    CArbMshTopo2D *msh_topo = BuildMeshTopo() ;
    DoSmoothNodesConsLaplace(msh_topo,1) ;

    delete msh_topo ;
}




// %(CArbMshSmooth2D::ConsLaplaceUpdate-CArbCoord2D-|-int-|-CArbArray-|<int>*-CArbArray-|<CArbCoord2D>*)
/* ++ ----------------------------------------------------------
**
**    ConsLaplaceUpdate - find an updated node coordinate 
**
**      CArbCoord2D ConsLaplaceUpdate(
**              int                     num_elem,
**              CArbArray<int>*         num_adj_nodes,
**              CArbArray<CArbCoord2D>* adj_coords)
**
**        num_elem      - (in)  number of elements around the node 
**        num_adj_nodes - (in)  number of nodes for each of the 
**                              adjacent elements 
**        adj_coords    - (in)  list of coordinates of nodes for each 
**                              of the adjacent elements 
**
**      Description: This function determines the updated location of a 
**          node during constrained Laplace smoothing. 
**
**      Return Value: the new coordinate location 
**
**
** -- */

CArbCoord2D CArbMshSmooth2D::ConsLaplaceUpdate(int num_elem,
                               CArbArray<int> *num_adj_nodes,
                               CArbArray<CArbCoord2D> *adj_coords)
{
    int i, j, k ;
    double tol = 0 ;

    // at this point we have a collection of all the elements
    // that surround an internal node.  We want to find the 
    // polygon that bounds the locations where we can place
    // this node and still have all valid elements.
    //
    // For a triangle this is easy, as long as the node is not
    // placed in a position that crosses the line between the other
    // two nodes the triangle is valid.
    //
    // For a quad, the node cannot cross any of the three lines
    // drawn between the other three nodes. 

    // Fist stuff all the lines that bound the region into
    // an array

    CArbArray<CArbCoord2D> lines ;
    int cur = 0 ;
    for (i=0 ; i<num_elem ; ++i) {
        if (num_adj_nodes->At(i) == 4) {
            lines.InsertAtEnd(adj_coords->At(cur)) ;
            lines.InsertAtEnd(CArbCoord2D(adj_coords->At(cur+1)-
                                          adj_coords->At(cur))) ;
            lines.InsertAtEnd(adj_coords->At(cur+3)) ;
            lines.InsertAtEnd(CArbCoord2D(adj_coords->At(cur)-
                                          adj_coords->At(cur+3))) ;
            lines.InsertAtEnd(adj_coords->At(cur+3)) ;
            lines.InsertAtEnd(CArbCoord2D(adj_coords->At(cur+1)-
                                          adj_coords->At(cur+3))) ;
            tol += (adj_coords->At(cur+1) - 
                    adj_coords->At(cur+3)).Magnitude() ;
            cur += 4 ;
        } else {
            lines.InsertAtEnd(adj_coords->At(cur)) ;
            lines.InsertAtEnd(CArbCoord2D(adj_coords->At(cur+1)-
                                          adj_coords->At(cur))) ;
            tol += (adj_coords->At(cur+1) - 
                    adj_coords->At(cur)).Magnitude() ;
            cur += 3 ;
        }
    }
    tol = 0.0000001 * tol/num_elem ;
    //huge = tol * 1e30 ;

    CArbFeasableRegion feasable(tol) ;
    for (i=0 ; i<lines.NumEntries() ; i += 2) {
        feasable.AddConstraint(lines[i],lines[i+1]) ;
    }
    CArbArray<CArbCoord2D> *poly = feasable.GetVertices() ;

    // find the average of these points

    if (poly->NumEntries() != 0) {
        CArbCoord2D center ;
        for (j=0 ; j<poly->NumEntries() ; ++j) {
            center += poly->At(j) ;
        }
        center /= poly->NumEntries() ;

        // find the centroid

        double sum_area = 0.0 ;
        double sum_prod_x = 0.0 ;
        double sum_prod_y = 0.0 ;

        for (j=0 ; j<poly->NumEntries() ; ++j) {
            k = (j+1) % poly->NumEntries() ;
            CArbCoord2D tcen = (poly->At(j)+poly->At(k)+center)/3 ;
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
        return (CArbCoord2D(center.x() + sum_prod_x/sum_area,
                            center.y() + sum_prod_y/sum_area)) ;
    } else {
        delete poly ;
        return (adj_coords->At(2)) ;
    }
}




// %(CArbMshSmooth2D::DoSmoothNodesConsLaplace-void-|-int-|)
/* ++ ----------------------------------------------------------
**
**    DoSmoothNodesConsLaplace - do constrained smoothing 
**
**      void DoSmoothNodesConsLaplace(int num_iters)
**
**        num_iters - (in)  number of smoothing iterations 
**
**      Description: This function does the actual work for the 
**          constrained Laplace smoothing. 
**
**
** -- */

int CArbMshSmooth2D::DoSmoothNodesConsLaplace(CArbMshTopo2D *msh_topo,
                                               int num_iters)
{
    int j ;

    // now the smoothing loop.

    NodeIter niter(NodeTable) ;
    for (j=0 ; j<num_iters ; ++j) {
        for (niter.First() ; niter.More() ; ++niter){

            ArbIntNode *node = niter.Entry() ;
            if ((node->motion == ARB_FIXED) || (!node->corner)) continue ;
            CArbTopoAdjVtxIterator iter(msh_topo,node->id) ;

            // fill two arrays.  The first has the number
            // of nodes in adjacent elements.  The second
            // has the acctual nodal coordinates.  Orient
            // the element so that the 3rd node is the one
            // being moved.

            CArbArray<int> num_adj_nodes ;
            CArbArray<CArbCoord2D> adj_coords ;
            int num_adj = 0 ;

            for (iter.First() ; iter.More() ; ++iter) {
                int nnodes,nodes[4] ;
                if (!msh_topo->GetElemNodes(iter.CcwElem(),
                                       node->id,&nnodes,nodes))
                 return 0;
                num_adj_nodes.InsertAtEnd(nnodes) ;
                if (nnodes == 4) {
                    for (j=0 ; j<4 ; ++j) {
                        ArbIntNode *adj =
                           NodeTable->Fetch(nodes[(j+2)%4]) ;
                        adj_coords.InsertAtEnd(adj->coord) ;
                    }
                } else {
                    for (j=0 ; j<3 ; ++j) {
                        ArbIntNode *adj =
                           NodeTable->Fetch(nodes[(j+1)%3]) ;
                        adj_coords.InsertAtEnd(adj->coord) ;
                    }
                }
                ++num_adj ;
            }

            node->coord = 
                ConsLaplaceUpdate(num_adj,&num_adj_nodes,&adj_coords) ;
        }
    }

    // if this is a quadratic order mesh, then we need to move
    // all the mid-side nodes to the proper locations

    if (Order == QUADRATIC) MoveSideNodes() ;
    return 1;
}


// position one node using the winslow update

static CArbCoord2D WinslowEval(int num_elem,CArbCoord2D node_coord,
                               int num_adj,CArbCoord2D *adj_coord,
                               double *ct,double *st)
{
    int i ;
    CArbCoord2D D_psi(0,0) ;
    CArbCoord2D D_eta(0,0) ;
    CArbCoord2D D_psi_psi(0,0) ;
    CArbCoord2D D_psi_eta(0,0) ;
    CArbCoord2D D_eta_eta(0,0) ;

    if (num_elem == 4) {
        for (i=0 ; i<4 ; ++i) {
            CArbCoord2D delta = adj_coord[i] - node_coord ;
            D_psi += delta * ct[i] ;
            D_eta += delta * st[i] ;
            D_psi_psi += delta * ct[i]*ct[i] ;
            // D_psi_eta += delta * ct[i] * st[i] ;
            D_eta_eta += delta * st[i]*st[i] ;
        }
        for (int i=4 ; i<num_adj ; ++i) {
            CArbCoord2D delta = adj_coord[i] - node_coord ;
            D_psi_eta += delta * ct[i] * st[i] ;
        }
        D_psi *= 0.5 ;
        D_eta *= 0.5 ;
        D_psi_eta *= 0.5 ;

    } else {
        for (i=0 ; i<num_adj ; ++i) {
            CArbCoord2D delta = adj_coord[i] - node_coord ;
            D_psi += delta * ct[i] ;
            D_eta += delta * st[i] ;
            D_psi_psi += delta * (4.0*ct[i]*ct[i] - 1.0) ;
            D_psi_eta += delta * ct[i] * st[i] ;
            D_eta_eta += delta * (4.0*st[i]*st[i] - 1.0) ;
        }

        D_psi *= (2.0 / num_adj) ;
        D_eta *= (2.0 / num_adj) ;
        D_psi_psi *= (2.0 / num_adj) ;
        D_psi_eta *= (8.0 / num_adj) ;
        D_eta_eta *= (2.0 / num_adj) ;
    }

    CArbCoord2D D_x = (D_eta*D_eta)*D_psi_psi -
                       2.0*(D_psi*D_eta)*D_psi_eta +
                      (D_psi*D_psi)*D_eta_eta ;
    return(D_x) ;
}



#define MAX_WITER 20

static CArbCoord2D WinslowUpdate(int num_elem,CArbCoord2D node_coord,
                                 int num_adj,CArbCoord2D *adj_coord,
                                 double *ct,double *st)
{
    // Determine a tolerance

    CArbCoord2D mean = 0.5*(adj_coord[0]+adj_coord[1]) ;
    double tol = 0.01 * (node_coord - mean).Magnitude() ;

    CArbCoord2D coord,lastf,lastc ;

    lastc = node_coord ;
//    lastf = WinslowEval(num_elem,lastf,num_adj,adj_coord,ct,st) ;
    lastf = WinslowEval(num_elem,node_coord,num_adj,adj_coord,ct,st) ;
//    if ((fabs(lastf[0]) < tol) && (fabs(lastf[1]) < tol)) return(lastc) ;

    // Compute the laplace coordinates and use them for
    // the starting guess

    coord[0] = coord[1] = 0 ;
    for (int i=0 ; i<num_elem ; ++i) coord += adj_coord[i] ;
    coord /= num_elem ;

    // deal with the case that the coord is already at the
    // centroid

    double delta = (adj_coord[0]-coord).Magnitude() * 0.0001 ;
    if (((lastc[0]-coord[0]) < delta) ||
        ((lastc[1]-coord[1]) < delta)) {
        coord[0] = coord[0] + 10000*delta ;
        coord[1] = coord[1] + 10000*delta ;
    }

    // go through a newton loop to try to update the coords

    bool done = false ;
    int iter = 0 ;
    while (!done) {
        CArbCoord2D eval = WinslowEval(num_elem,coord,num_adj,
                                       adj_coord,ct,st) ;
        if (iter >= MAX_WITER) {
            done = true ;
        } else {
            CArbCoord2D slope = (lastf - eval) / (lastc - coord) ;
            lastf = eval ;
            lastc = coord ;
            coord -= eval / slope ;
            if ((fabs(lastc[0]-coord[0]) < tol) &&
                (fabs(lastc[1]-coord[1]) < tol)) done = true ;
        }
        ++iter ;
    }

    return(coord) ;
}


// %(CArbMshSmooth2D::DoSmoothNodesWinslow-void-|)
/* ++ ----------------------------------------------------------
**
**    DoSmoothNodesWinslow - do Winslow smoothing 
**
**      void DoSmoothNodesWinslow()
**
**      Description: This function does the actual work for the Winslow 
**          smoothing. 
**
**
** -- */

void CArbMshSmooth2D::DoSmoothNodesWinslow(CArbMshTopo2D *msh_topo,
                                           int max_it,double adapt_tol)
{
    int j ;
    CArbCoord2D adj_coord[100] ;
    double ct[100], st[100] ;

    // now the smoothing loop. 

    bool first = true ;
    double max_delta = 0.0 ;
    Num_it = 0 ;

    CArbCoord2D old_max,new_max ;

#ifdef SMOOTH_DEBUG
    fprintf(stderr,"stop tolerance: %g\n",adapt_tol) ;
    fprintf(stderr,"node delta old_coords new_coords\n") ;
#endif

    NodeIter niter(NodeTable) ;
    while (first || (max_delta > adapt_tol)) {
        if (Num_it >= max_it) break ;
        Num_it++ ;
        max_delta = 0.0 ;

#ifdef SMOOTH_DEBUG
    fprintf(stderr,"%d \n",Num_it) ;
#endif

        for (niter.First() ; niter.More() ; ++niter){

            ArbIntNode *node = niter.Entry() ;
            if ((node->motion == ARB_FIXED) || (!node->corner)) continue ;
            int num = 0 ;
            CArbTopoAdjVtxIterator iter(msh_topo,node->id) ;

            for (iter.First() ; iter.More() ; ++iter) ++num ;
            if (num == 0) continue ;

            CArbCoord2D coord(0,0) ;

            if (num == 2) {
                for (iter.First() ; iter.More() ; ++iter) {
                    ArbIntNode *adj =
                       NodeTable->Fetch(iter.AdjVtx()) ;
                    coord[0] += adj->coord[0] ;
                    coord[1] += adj->coord[1] ;
                }
                if (num != 0) {
                    coord[0] /= double(num) ;
                    coord[1] /= double(num) ;
                }
            } else if (num == 3) {

                // have found that using the element far corner nodes
                // for the case where the node is adjacent to only
                // three elements can lead to trouble in some cases,
                // so we fall back to laplace.

                coord[0] = coord[1] = 0.0 ;
                for (iter.First() ; iter.More() ; ++iter) {
                    ArbIntNode *adj =
                       NodeTable->Fetch(iter.AdjVtx()) ;
                    coord += adj->coord ;
                }
                coord = coord / 3 ; 

//                 int num_adj = 0 ;
// 
//                 for (iter.First() ; iter.More() ; ++iter) {
//                     int nnodes, nodes[4] ;
//                     ArbIntNode *adj =
//                        NodeTable->Fetch(iter.AdjVtx()) ;
//                     adj_coord[num_adj] = adj->coord ;
//                     double theta = 2.0*PI*num_adj/6 ;
//                     ct[num_adj] = cos(theta) ;
//                     st[num_adj] = sin(theta) ;
//                     ++num_adj ;
// 
//                     msh_topo->GetElemNodes(iter.CcwElem(),
//                                        node->id,&nnodes,nodes) ;
//                     if (nnodes == 4) {
//                         adj = NodeTable->Fetch(nodes[2]) ;
//                         adj_coord[num_adj] = adj->coord ;
//                         theta = 2.0*PI*num_adj/6 ;
//                         ct[num_adj] = cos(theta) ;
//                         st[num_adj] = sin(theta) ;
//                         ++num_adj ;
//                     }
//                 }
// 
//                 if (num_adj == 6) {
//                     coord = WinslowUpdate(3,node->coord,num_adj,
//                                           adj_coord,ct,st) ;
//                 } else {
//                     coord[0] = coord[1] = 0.0 ;
//                     for (j=0 ; j<5 ; j+=2) coord += adj_coord[j] ;
//                     coord = coord / 3 ;
//                 }

            } else if (num == 4) {

                int num_adj = 0 ;

                for (iter.First() ; iter.More() ; ++iter) {
                    ArbIntNode *adj =
                       NodeTable->Fetch(iter.AdjVtx()) ;
                    adj_coord[num_adj] = adj->coord ;
                    double theta = 2.0*PI*num_adj/4 ;
                    ct[num_adj] = cos(theta) ;
                    st[num_adj] = sin(theta) ;
                    ++num_adj ;
                }

                for (iter.First() ; iter.More() ; ++iter) {
                    int nnodes, nodes[4] ;
                    msh_topo->GetElemNodes(iter.CcwElem(),
                                       node->id,&nnodes,nodes) ;
                    if (nnodes == 4) {
                        ArbIntNode *adj = NodeTable->Fetch(nodes[2]) ;
                        adj_coord[num_adj] = adj->coord ;
                        double theta = 2.0*PI*((num_adj-4)+0.5)/4 ;
                        ct[num_adj] = cos(theta) ;
                        st[num_adj] = sin(theta) ;
                        ++num_adj ;
                    }
                }

                if (num_adj == 8) {
                    coord = WinslowUpdate(4,node->coord,num_adj,
                                          adj_coord,ct,st) ;
                } else {
                    coord[0] = coord[1] = 0.0 ;
                    for (j=0 ; j<4 ; ++j) coord += adj_coord[j] ;
                    coord = coord / 4 ;
                }
            } else {
                int num_adj = 0 ;

                for (iter.First() ; iter.More() ; ++iter) {
                    ArbIntNode *adj =
                       NodeTable->Fetch(iter.AdjVtx()) ;
                    adj_coord[num_adj] = adj->coord ;
                    double theta = 2.0*PI*num_adj/num ;
                    ct[num_adj] = cos(theta) ;
                    st[num_adj] = sin(theta) ;
                    ++num_adj ;
                    if (num_adj == 100) break ;
                }

                coord = WinslowUpdate(num_adj,node->coord,num_adj,
                                      adj_coord,ct,st) ;
            }

            double delta = (node->coord - coord).Magnitude() ;
            if (delta > max_delta) {
                max_delta = delta ;
                old_max = node->coord ;
                new_max = coord ;

#ifdef SMOOTH_DEBUG
    fprintf(stderr,"%d %g (%g,%g) (%g,%g)\n",node->id,delta,
            node->coord.x(),node->coord.y(),coord.x(),coord.y()) ;
#endif

            }
            node->coord = coord ;
        }

        if (first) {
            first = false ;
            if (max_delta == 0.0) break ;
        }
//         fprintf(stderr,"Cnt: %d %g %d (%g,%g) (%g,%g)\n",cnt,max_delta,
//                 delt_node,old_max.x(),old_max.y(),new_max.x(),new_max.y()) ;
//         if (cnt == 1) DisplayMesh("In Winslow:") ;
//         DisplayMesh("In Winslow:") ;
    }

    // if this is a quadratic order mesh, then we need to move
    // all the mid-side nodes to the proper locations

    if (Order == QUADRATIC) MoveSideNodes() ;

#ifdef SMOOTH_DEBUG
    fprintf(stderr,"Winslow Num Iter: %d\n\n",Num_it) ;
#endif
}




// %(CArbMshSmooth2D::DoSmoothNodesLaplace-void-|)
/* ++ ----------------------------------------------------------
**
**    DoSmoothNodesLaplace - do Laplace smoothing 
**
**      void DoSmoothNodesLaplace()
**
**      Description: This function does the actual work for the Laplace 
**          smoothing. 
**
**
** -- */

int CArbMshSmooth2D::DoSmoothNodesLaplace(CArbMshTopo2D *msh_topo,
                                           int max_it,double adapt_tol,
                                           bool skip_checks)
{
    int j ;

//    DisplayMesh("Start Laplace:") ;

#ifdef SMOOTH_DEBUG
    fprintf(stderr,"stop tolerance: %g\n",adapt_tol) ;
    fprintf(stderr,"node delta old_coords new_coords\n") ;
#endif

    // now the smoothing loop.  We go through each floating
    // node and compute a new position at the centroid of it's
    // adjacent nodes.

    bool first = true ;
    double max_delta = 0.0 ;
    Num_it = 0 ;

    NodeIter niter(NodeTable) ;
    while (first || (max_delta > adapt_tol)) {
        if (Num_it >= max_it) break ;
        Num_it++ ;

#ifdef SMOOTH_DEBUG
    fprintf(stderr,"%d \n",Num_it) ;
#endif
        max_delta = 0.0 ;
        for (niter.First() ; niter.More() ; ++niter){

            ArbIntNode *node = niter.Entry() ;
            if ((node->motion == ARB_FIXED) || (!node->corner)) continue ;

            CArbTopoAdjVtxIterator iter(msh_topo,node->id) ;

            // fill two arrays.  The first has the number
            // of nodes in adjacent elements.  The second
            // has the acctual nodal coordinates.  Orient
            // the element so that the 3rd node is the one
            // being moved.

            CArbArray<int> num_adj_nodes ;
            CArbArray<CArbCoord2D> adj_coords ;
            int num_adj = 0 ;

            for (iter.First() ; iter.More() ; ++iter) {
                int nnodes,nodes[4] ;
                if (!msh_topo->GetElemNodes(iter.CcwElem(),
                                       node->id,&nnodes,nodes))
                 return 0;
                num_adj_nodes.InsertAtEnd(nnodes) ;
                if (nnodes == 4) {
                    for (j=0 ; j<4 ; ++j) {
                        ArbIntNode *adj =
                           NodeTable->Fetch(nodes[(j+2)%4]) ;
                        adj_coords.InsertAtEnd(adj->coord) ;
                    }
                } else {
                    for (j=0 ; j<3 ; ++j) {
                        ArbIntNode *adj =
                           NodeTable->Fetch(nodes[(j+1)%3]) ;
                        adj_coords.InsertAtEnd(adj->coord) ;
                    }
                }
                ++num_adj ;
            }

            CArbCoord2D coord =
                LaplaceUpdate(num_adj,&num_adj_nodes,&adj_coords,
                              skip_checks) ;

            double delta = (node->coord - coord).Magnitude() ;
            if (delta > max_delta) {

#ifdef SMOOTH_DEBUG
    fprintf(stderr,"%d %g (%g,%g) (%g,%g)\n",node->id,delta,
            node->coord.x(),node->coord.y(),coord.x(),coord.y()) ;
#endif

                max_delta = delta ;
            }
            node->coord = coord ;
        }

        if (first) {
            first = false ;
            if (max_delta == 0.0) break ;
        }
        // fprintf(stderr,"max, delta, it: %g %g %d\n",max_delta,adapt_tol,cnt) ;
        // DisplayMesh("In Laplace:") ;
    }

    // if this is a quadratic order mesh, then we need to move
    // all the mid-side nodes to the proper locations

    if (Order == QUADRATIC) MoveSideNodes() ;

//    DisplayMesh("Finish Laplace:") ;

#ifdef SMOOTH_DEBUG
    fprintf(stderr,"Laplace Num Iter: %d\n\n",Num_it) ;
#endif
    return 1;
}




// %(CArbMshSmooth2D::LaplaceUpdate-CArbCoord2D-|-int-|-CArbArray-|<int>*-CArbArray-|<CArbCoord2D>*)
/* ++ ----------------------------------------------------------
**
**    LaplaceUpdate - find an updated node coordinate 
**
**      CArbCoord2D LaplaceUpdate(
**              int                     num_elem,
**              CArbArray<int>*         num_adj_nodes,
**              CArbArray<CArbCoord2D>* adj_coords)
**
**        num_elem      - (in)  number of elements around the node 
**        num_adj_nodes - (in)  number of nodes for each of the 
**                              adjacent elements 
**        adj_coords    - (in)  list of coordinates of nodes for each 
**                              of the adjacent elements 
**
**      Description: This function determines the updated location of a 
**          node during normal Laplace smoothing 
**
**      Return Value: the new coordninate location 
**
**
** -- */

CArbCoord2D CArbMshSmooth2D::LaplaceUpdate(int num_elem,
                               CArbArray<int> *num_adj_nodes,
                               CArbArray<CArbCoord2D> *adj_coords,
                               bool skip_checks)
{
    int i, cur=0 ;
    CArbCoord2D coord(0,0) ;
    CArbArray<double> shape ;
    double sm ;

    // first find the delta to the center of the adjacent nodes

    for (i=0 ; i<num_elem ; ++i) {
        coord += adj_coords->At(cur+1) ;
        if (num_adj_nodes->At(i) == 4) {
            sm = QuadMetric(adj_coords->At(cur),
                            adj_coords->At(cur+1),
                            adj_coords->At(cur+2),
                            adj_coords->At(cur+3)) ;
            cur += 4 ;
        } else {
            sm = TriMetric(adj_coords->At(cur),
                           adj_coords->At(cur+1),
                           adj_coords->At(cur+2)) ;
            cur += 3 ;
        }
        shape.InsertAtEnd(sm) ;
    }

    CArbCoord2D delta = coord/double(num_elem) - adj_coords->At(2) ;
    if (skip_checks) return(adj_coords->At(2) + delta) ;

    // now check to see if  all elements are valid.  If not,
    // cut the motion in half and try again

    for (int step=0 ; step<10 ; ++step) {
        coord = adj_coords->At(2) + delta ;

        // loop around the element adjacent to this
        // node and see how the this delta changes
        // their shape measure

        int num_minus = 0 ;
        int num_plus = 0 ;
        int num_minus_minus = 0 ;
        int num_plus_plus = 0 ;

        for (i=0 ; i<num_elem ; ++i) {
            cur = 0 ;
            if (num_adj_nodes->At(i) == 4) {
                sm = QuadMetric(adj_coords->At(cur),
                                adj_coords->At(cur+1),
                                coord,
                                adj_coords->At(cur+3)) ;
                cur += 4 ;
            } else {
                sm = TriMetric(adj_coords->At(cur),
                               adj_coords->At(cur+1),
                               coord) ;
                cur += 3 ;
            }

            // what follows is a collection of byzantine
            // huristics used to decide if we do the update or
            // not.  These were developed by Canann, Tristano, and
            // Staten in their paper entitled "An Approach to
            // Combined Laplacian and Optimization-Based
            // Smoothing for Triangular, Quadrilateral, and
            // Quad-Dominant Meshes".  I have a copy of the
            // paper, but not a reference.  I don't know where
            // it was published

            if (sm < shape[i])
                num_minus++ ;
            else
                num_plus++ ;

            if ((shape[i] < 0.0) && (sm > shape[i])) {
                num_plus_plus++ ;
            } else if ((shape[i] < 0.05) && (sm > 0.05)) {
                num_plus_plus++ ;
            }

            if ((shape[i] > 0.0) && (sm < 0.0)) {
                num_minus_minus++ ;
            } else if ((shape[i] < 0.0) && (sm < shape[i])) {
                num_minus_minus++ ;
            } else if ((shape[i] > 0.05) && (sm < 0.05)) {
                num_minus_minus++ ;
            }
        }

        // check to see if this motion should
        // be rejected or accepted

        if ((num_minus < num_elem) &&
            (num_plus_plus >= num_minus_minus)) {
            if ((num_plus == num_elem) ||
                ((num_plus_plus > 0) && (num_minus_minus == 0)) ||
                (num_plus_plus >= num_minus_minus)) {
                return(coord) ;
            }
        }

        // reduce the size of the delta

        delta *= 0.5 ;
    }

    // if we get here, the node should not be moved

    return(adj_coords->At(2)) ;
}

double CArbMshSmooth2D::FindMinCharEdgeLength(CArbMshTopo2D *msh_topo)
{
    // find the minimum of the average lengths of all the
    // edges connected to a node

    double min_length = 0.0 ;
    bool first = true ;

    NodeIter niter(NodeTable) ;
    for (niter.First() ; niter.More() ; ++niter){

        ArbIntNode *node = niter.Entry() ;
        if ((node->motion == ARB_FIXED) || (!node->corner)) continue ;
        CArbTopoAdjVtxIterator iter(msh_topo,node->id) ;

        double sum = 0 ;
        int num = 0 ;
        for (iter.First() ; iter.More() ; ++iter) {

            ArbIntNode *anode = NodeTable->Fetch(iter.AdjVtx()) ;
            if (anode == 0)
             return -1;
            sum += (anode->coord - node->coord).Magnitude() ;
            ++num ;

        }

        // find the average length and see if this is the smallest
        // so far

        if (num > 0) {
            sum /= double(num) ;
            if (first || (sum < min_length)) min_length = sum ;
            first = false ;
        }
    }
    return(min_length) ;
}


// %(CArbMshSmooth2D::FindBoundaryNodes-void-|)
/* ++ ----------------------------------------------------------
**
**    FindBoundaryNodes - find the fixed nodes 
**
**      void FindBoundaryNodes()
**
**      Description: This function finds the nodes on the boundary and 
**          on bi-material interfaces, and flags them as being fixed. 
**
**
** -- */

void CArbMshSmooth2D::FindBoundaryNodes()
{
    CArbHashTable<ArbEdgeKey,int> edges ;

    ElemIter eiter(ElemTable) ;
    for (eiter.First() ; eiter.More() ; ++eiter) {

        ArbMshElement2D *elem = eiter.Entry() ;
        int nn ;
        if (elem->num_nodes == 6)
            nn = 3 ;
        else if (elem->num_nodes == 8)
            nn = 4 ;
        else
            nn = elem->num_nodes ;

        for (int i=0 ; i<nn ; ++i) {
            int j = (i+1) % nn ;
            ArbEdgeKey key ;
            if (elem->nodes[i] < elem->nodes[j])
                key = ArbEdgeKey(elem->nodes[i],elem->nodes[j]) ;
            else
                key = ArbEdgeKey(elem->nodes[j],elem->nodes[i]) ;

            ArbIntNode *nd = NodeTable->Fetch(elem->nodes[i]) ;
            nd->corner = true ;
            nd = NodeTable->Fetch(elem->nodes[j]) ;
            nd->corner = true ;

            // check to see if this edge is in the table

            int *value = edges.Fetch(key) ;

            // case 1: not in the table so insert it

            if (value == 0) {
                edges.Store(key,elem->mat_id) ;

            // case 2: same materials, remove table

            } else if (*value == elem->mat_id) {
                edges.Remove(key) ;
            }

            // case 3: different materials, leave in the table
        }
    }

    // All the nodes still in the table should be marked
    // as fixed nodes.
 
    CArbHashTableIterator<ArbEdgeKey,int> iter(&edges) ;

    for (iter.First() ; iter.More() ; ++iter){
        ArbIntNode *nd = NodeTable->Fetch(iter.Key().id_0) ;
        nd->motion = ARB_FIXED ;
        nd = NodeTable->Fetch(iter.Key().id_1) ;
        nd->motion = ARB_FIXED ;
    }
}




// %(CArbMshSmooth2D::BuildMeshTopo-CArbMshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    BuildMeshTopo - build a mesh topology object 
**
**      CArbMshTopo2D *BuildMeshTopo()
**
**      Description: This function builds a MshTopo2D object for the 
**          mesh described by the node and element tables. 
**
**      Return Value: a new MshTopo2D object 
**
**
** -- */

CArbMshTopo2D *CArbMshSmooth2D::BuildMeshTopo()
{
    CArbMshTopo2D *msh_topo = new CArbMshTopo2D() ;

    ElemIter eiter(ElemTable) ;
    for (eiter.First() ; eiter.More() ; ++eiter) {
        ArbMshElement2D *elem = eiter.Entry() ;
        int nn ;
        if (elem->num_nodes == 6)
            nn = 3 ;
        else if (elem->num_nodes == 8)
            nn = 4 ;
        else
            nn = elem->num_nodes ;
        msh_topo->InsertCollapsedElement(elem->elem_id,nn,elem->nodes) ;
    }
    return(msh_topo) ;
}




// %(CArbMshSmooth2D::MoveSideNodes-void-|)
/* ++ ----------------------------------------------------------
**
**    MoveSideNodes - update the coordinates of the side nodes 
**
**      void MoveSideNodes()
**
**      Description: This function updates the coordinates of all the 
**          side nodes (in the case of a quadratic order mesh) to place 
**          them at them midway between the corresponding corner nodes. 
**
**
** -- */

void CArbMshSmooth2D::MoveSideNodes()
{
    int j ;
    ElemIter eiter(ElemTable) ;
    for (eiter.First() ; eiter.More() ; ++eiter) {
        ArbMshElement2D *elem = eiter.Entry() ;
        int prev, next, side ;
        int nn, lnodes[8] ;
        for (j=0 ; j<elem->num_nodes ; ++j)
            lnodes[j] = elem->nodes[j] ;

        if (elem->num_nodes == 6)
            nn = 3 ;
        else
            nn = 4 ;

        for (j=0 ; j<nn ; ++j) {
            prev = lnodes[j] ;
            next = lnodes[(j+1) % nn] ;
            side = lnodes[j + nn] ;

            ArbIntNode *pside = NodeTable->Fetch(side) ;

            if (pside->motion == ARB_FLOATING) {
                ArbIntNode *pprev = NodeTable->Fetch(prev) ;
                ArbIntNode *pnext = NodeTable->Fetch(next) ;
                pside->coord = 0.5 * (pprev->coord + pnext->coord) ;
            }
        }
    }
}




// %(CArbMshSmooth2D::NumNodes-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumNodes - return the number of nodes in the mesh 
**
**      int NumNodes() const
**
**      Description: This function returns the current number of nodes 
**          in a mesh (typically called after smoothing). 
**
**      Return Value: number of nodes 
**
**
** -- */

int CArbMshSmooth2D::NumNodes() const
{
    return(NodeTable->NumEntries()) ;
}




// %(CArbMshSmooth2D::NumElements-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumElements - return the number of elements in the mesh 
**
**      int NumElements() const
**
**      Description: This function returns the current number of 
**          elements in a mesh (typically called after smoothing). 
**
**      Return Value: number of elements 
**
**
** -- */

int CArbMshSmooth2D::NumElements() const
{
    return(ElemTable->NumEntries()) ;
}




// %(CArbMshSmooth2D::GetNodes-ArbMshNode-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetNodes - return a list of all mesh nodes 
**
**      ArbMshNode *GetNodes() const
**
**      Description: This function returns a list of all the nodes 
**          current defined for the mesh (typically called after 
**          smoothing). 
**
**      Return Value: A list of all mesh nodes. Ownership of this 
**          memory passes to the client, which must eventually call 
**          delete []. 
**
**
** -- */

ArbMshNode *CArbMshSmooth2D::GetNodes() const
{
    int i ;
    ArbMshNode   *local_entries ;

    NodeIter iter(NodeTable) ;
    local_entries = new ArbMshNode[NodeTable->NumEntries()] ;

    for (i=0,iter.First() ; iter.More() ; ++i,++iter) {
        local_entries[i].id = iter.Key() ;
        local_entries[i].coord[0] = iter.Entry()->coord.x() ;
        local_entries[i].coord[1] = iter.Entry()->coord.y() ;
        if (ZCoords != 0) {
            local_entries[i].coord[2] =
                *(ZCoords->Fetch(iter.Entry()->id)) ;
        } else {
            local_entries[i].coord[2] = 0.0 ;
        }
    }

    return(local_entries) ;
}




// %(CArbMshSmooth2D::GetElements-ArbMshElement2D-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetElements - return a list of all mesh elements 
**
**      ArbMshElement2D *GetElements() const
**
**      Description: This function returns a list of all the elements 
**          current defined for the mesh (typically called after 
**          smoothing). 
**
**      Return Value: A list of all mesh elements. Ownership of this 
**          memory passes to the client, which must eventually call 
**          delete []. 
**
**
** -- */

ArbMshElement2D *CArbMshSmooth2D::GetElements() const
{
    CArbHashTableIterator<int,ArbMshElement2D> iter(ElemTable) ;
    ArbMshElement2D *local_elems = 
        new ArbMshElement2D[ElemTable->NumEntries()] ;

    for (int i=0 ; iter.More() ; ++i,++iter)
        local_elems[i] = *(iter.Entry()) ;

    return(local_elems) ;
}




// %(CArbMshSmooth2D::DisplayMesh-void-|)
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

void CArbMshSmooth2D::DisplayMesh(char *label)
{
//    bool do_labels = false ;
    bool do_labels = true ;

    ElemIter eiter(ElemTable) ;
    for (eiter.First() ; eiter.More() ; ++eiter) {
        ArbMshElement2D *elem = eiter.Entry() ;
        double x[8], y[8] ;
        int id[8] ;
        int nn ;
        if (elem->num_nodes == 8)
            nn = 4 ;
        else if (elem->num_nodes == 6)
            nn = 3 ;
        else
            nn = elem->num_nodes ;

        for (int j=0 ; j<elem->num_nodes ; ++j) {
            ArbIntNode *node = NodeTable->Fetch(elem->nodes[j]) ;
            x[j] = node->coord[0] ;
            y[j] = node->coord[1] ;
            id[j] = node->id ;
        }
        for (int jj=0 ; jj<nn ; ++jj) {
            int kk = (jj+1) % nn ;
            printf("l %g %g %g %g\n",x[jj],y[jj],x[kk],y[kk]) ;
        }
        if (do_labels) {
            for (int jj=0 ; jj<elem->num_nodes ; ++jj) {
                printf("t %g %g %d\n",x[jj],y[jj],id[jj]) ;
            }
        }
    }
    printf("a %s\n",label) ;
    printf("n\n") ;
    fflush(stdout) ;
}




// %(CArbMshSmooth2D::Angle-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    Angle - compute an angle 
**
**      double Angle(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j) const
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

double CArbMshSmooth2D::Angle(CArbCoord2D b,CArbCoord2D i,
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




// %(CArbMshSmooth2D::CrossProd-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    CrossProd - compute a cross product 
**
**      double CrossProd(
**              CArbCoord2D b,
**              CArbCoord2D i,
**              CArbCoord2D j) const
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

double CArbMshSmooth2D::CrossProd(CArbCoord2D b,CArbCoord2D i,
                                  CArbCoord2D j) const
{
    double cross ;

    cross = ((i[0] - b[0]) * (j[1] - b[1])) -
            ((i[1] - b[1]) * (j[0] - b[0])) ;
    return cross ;
}


/* ------------------------------------------------------------
    LenSqr - computes the square of the distance between to points
*/

inline double LenSqr(CArbCoord2D i,CArbCoord2D j)
{
    CArbCoord2D delta = i - j ;
    return(delta.x()*delta.x() + delta.y()*delta.y()) ;
}




// %(CArbMshSmooth2D::TriMetric-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
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

double CArbMshSmooth2D::TriMetric(
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
    double sum = LenSqr(i,j) + LenSqr(j,k) + LenSqr(k,i) ;
    static double factor = 3.464101615 ;  // 2 * sqrt(3)
    return(factor * cross / sum) ;
}




// %(CArbMshSmooth2D::QuadMetric-double-|^const-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|-CArbCoord2D-|)
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

double CArbMshSmooth2D::QuadMetric(
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




