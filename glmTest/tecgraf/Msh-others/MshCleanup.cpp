//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000
//   All rights reserved
//
// Revision -
//   $Revision: 1.17 $  $Date: 2004/06/01 15:46:45 $  $Author: wash $
//

#include "MshQuad2D.hpp"

namespace Msh2D {

#ifdef MEMDEBUG
#include "MemDbg.hpp"
//#define new new(__FILE__,__LINE__)
#endif

#define BIG_CTV 1000.0

/* ------------------------------------------------------------
    TopoVariance - compute a measure of topological variance
                   from an optimal connection to 4 nodes
*/

static double TopoVariance(int num_adj)
{
    switch (num_adj) {
        case 0: return(4.0) ;      break ;
        case 1: return(4.0) ;      break ;
        case 2: return(4.0) ;      break ;
        case 3: return(0.444444) ; break ;
        case 4: return(0.0) ;      break ;
        case 5: return(0.16) ;     break ;
        case 6: return(0.444444) ; break ;
    }
    double t0 = 360.0 / num_adj ;
    return(((t0*t0)/2025.0) - (t0/11.25) + 4.0) ;
}


/* ------------------------------------------------------------
    DoCollapse - do a collapse of a quad element
*/

static int DoCollapse(MshTopo2D *quad_topo,int *nodes,
                      Dict<int,IntNode> *node_table,
                      bool do_1_3)
{
    List<MshRegion2D::SeamElemCache> elem_cache ;
    int nd0,nd1,i ;

    // determine which node to remove (nd0) by checking
    // if any of the nodes are on the boundary

    if (do_1_3) {
        IntNode *ndp = node_table->Get(nodes[3]) ;
        if (ndp->motion != MSH_FIXED) {
            nd0 = nodes[3] ;
            nd1 = nodes[1] ;
        } else {
            ndp = node_table->Get(nodes[1]) ;
            if (ndp->motion != MSH_FIXED) {
                nd0 = nodes[1] ;
                nd1 = nodes[3] ;
            } else {
                return(-1) ;       
            }
        }
    } else {
        IntNode *ndp = node_table->Get(nodes[2]) ;
        if (ndp->motion != MSH_FIXED) {
            nd0 = nodes[2] ;
            nd1 = nodes[0] ;
        } else {
            ndp = node_table->Get(nodes[0]) ;
            if (ndp->motion != MSH_FIXED) {
                nd0 = nodes[0] ;
                nd1 = nodes[2] ;
            } else {
                return(-1) ;       
            }
        }
    }
    quad_topo->DeleteElement(4,nodes) ;

    TopoAdjVtxIterator iter(quad_topo,nd0) ;

    for (iter.First() ; iter.More() ; ++iter) {
        if (iter.CcwElem() != NO_ELEM) {
            MshRegion2D::SeamElemCache edata ;
            edata.elem = iter.CcwElem() ;
            quad_topo->GetElemNodes(iter.CcwElem(),
                nd0,&edata.num_nodes,edata.nodes) ;

            if (edata.num_nodes > 4) throw MshRegion2D::CleanupMeshError() ;

            elem_cache.Append(edata) ;
        }
    }
    for (i=0 ; i<elem_cache.Len() ; ++i) {
        quad_topo->DeleteElement(elem_cache[i].num_nodes,
                                 elem_cache[i].nodes) ;
    }
    for (i=0 ; i<elem_cache.Len() ; ++i) {
        elem_cache[i].nodes[0] = nd1 ;
        quad_topo->InsertElement(elem_cache[i].elem,
                                 elem_cache[i].num_nodes,
                                 elem_cache[i].nodes) ;
    }
    return(nd1) ;
}




/* ------------------------------------------------------------
    QuadTwoEdge - search for nodes that are adjacent to only
                  two edges

    This routine looks for internal nodes that are adjacent
    to only two edges.  These nodes are deleted and the two
    adjacent elements are collapsed into one.

        #             *     where:
       /|\           / \      # indicates nodes connected
      / | \         /   \       to more than 4 nodes, and
     *  o  *  ==>  *     * 
      \ | /         \   /     o indicates nodes connected
       \|/           \ /        to less than 4 nodes.
        #             *     

*/

static bool QuadTwoEdge(MshTopo2D *quad_topo,
               Dict<int,IntNode> *node_table)
{
    bool updates = false ;
    int num_nodes = quad_topo->NumNodes() ;
    int *nodes = quad_topo->GetNodeList() ;

    for (int i=0 ; i<num_nodes ; ++i) {
        int elems[3] ;
        int num = 0 ;
        TopoAdjVtxIterator iter(quad_topo,nodes[i]) ;

        for (iter.First() ; iter.More() ; ++iter) {
            if (iter.CcwElem() == NO_ELEM) {
                num = 0 ;
                break ;
            }
            elems[num] = iter.CcwElem() ;
            ++num ;
            if (num > 2) break ;
        }

        if (num == 2) {
            int num0,num1,nodes0[4],nodes1[4] ;
            quad_topo->GetElemNodes(elems[0],nodes[i],&num0,nodes0) ;
            quad_topo->GetElemNodes(elems[1],nodes[i],&num1,nodes1) ;

            if ((num0 != 4) || (num1 != 4)) continue ;

            quad_topo->DeleteElement(num0,nodes0) ;
            quad_topo->DeleteElement(num1,nodes1) ;

            IntNode *nd0 = node_table->Get(nodes0[1]) ;
            IntNode *nd1 = node_table->Get(nodes0[2]) ;
            IntNode *nd2 = node_table->Get(nodes1[2]) ;
            if ((nd0->motion == MSH_FIXED) &&
                (nd1->motion == MSH_FIXED) &&
                (nd2->motion == MSH_FIXED)) {
                quad_topo->InsertElement(elems[0],3,&nodes0[1]) ;
                quad_topo->InsertElement(elems[1],3,&nodes1[1]) ;
            } else {
                IntNode *nd0 = node_table->Get(nodes1[1]) ;
                IntNode *nd1 = node_table->Get(nodes1[2]) ;
                IntNode *nd2 = node_table->Get(nodes0[2]) ;
                if ((nd0->motion == MSH_FIXED) &&
                    (nd1->motion == MSH_FIXED) &&
                    (nd2->motion == MSH_FIXED)) {
                    quad_topo->InsertElement(elems[0],3,&nodes0[1]) ;
                    quad_topo->InsertElement(elems[1],3,&nodes1[1]) ;
                } else {
                    nodes0[0] = nodes1[2] ;
                    quad_topo->InsertElement(elems[0],4,nodes0) ;
                }
            }
        updates = true ;
        }
    }
    delete [] nodes ;
    return(updates) ;
}




/* ------------------------------------------------------------
    QuadCollapse - search for quads to collapse

    This routine looks for quads that are candidates for
    being collapsed.  If two of the nodes are connected to
    more than 4 nodes, and two are connected to less than
    4 nodes, then the element is collapsed.

           #          *     where:
          / \         |       # indicates nodes connected
         /   \        |         to more than 4 nodes, and
        o     o  ==>  * 
         \   /        |       o indicates nodes connected
          \ /         |         to less than 4 nodes.
           #          * 

*/

static bool QuadCollapse(MshTopo2D *quad_topo,
                         Dict<int,IntNode> *node_table)
{
    bool updates = false ;
    int num_elems = quad_topo->NumElements() ;
    int *elems = quad_topo->GetElemList() ;

    for (int i=0 ; i<num_elems ; ++i) {
        int num, nodes[4] ;
        quad_topo->GetElemNodes(elems[i],&num,nodes) ;

        if (num > 4) throw MshRegion2D::CleanupMeshError() ;

        if (num == 4) {
            if (quad_topo->BoundaryNode(nodes[0])) continue ;
            if (quad_topo->BoundaryNode(nodes[1])) continue ;
            if (quad_topo->BoundaryNode(nodes[2])) continue ;
            if (quad_topo->BoundaryNode(nodes[3])) continue ;
            int n0 = quad_topo->NumAdjElems(nodes[0]) ;
            if (n0 > 4) {
                int n2 = quad_topo->NumAdjElems(nodes[2]) ;
                if (n2 > 4) {
                    int n1 = quad_topo->NumAdjElems(nodes[1]) ;
                    if (n1 < 4) {
                        int n3 = quad_topo->NumAdjElems(nodes[3]) ;
                        if (n3 < 4) {
                            if (DoCollapse(quad_topo,nodes,node_table,true) != -1)
                                updates = true ;
                        }
                    }
                }
            } else if (n0 < 4) {
                int n2 = quad_topo->NumAdjElems(nodes[2]) ;
                if (n2 < 4) {
                    int n1 = quad_topo->NumAdjElems(nodes[1]) ;
                    if (n1 > 4) {
                        int n3 = quad_topo->NumAdjElems(nodes[3]) ;
                        if (n3 > 4) {
                            if (DoCollapse(quad_topo,nodes,node_table,false) != -1)
                                updates = true ;
                        }
                    }
                }
            }
        }
    }
    delete [] elems ;
    return(updates) ;
}




/* ------------------------------------------------------------
    QuadOpen - search for places to open a node and create
               a new quad

    This routine looks for nodes that are candidates for
    being opened.  If a node is connected to 5 or more
    nodes, and one or more of the adjacent nodes are
    is connected to less than 4 nodes, then the node is
    split and an new element is added.

        o          *        where:
        |         / \         # indicates nodes connected
        |        /   \          to more than 4 nodes, and
        #  ==>  *     * 
        |        \   /        o indicates nodes connected
        |         \ /           to less than 4 nodes.
        *          *     

*/

static bool QuadOpen(MshTopo2D *quad_topo,
                     MshRegion2D *region)
{
    bool updates = false ;
    int num_nodes = quad_topo->NumNodes() ;
    int *nodes = quad_topo->GetNodeList() ;

    for (int i=0 ; i<num_nodes ; ++i) {
        List<int> node_cache ;
        TopoAdjVtxIterator iter(quad_topo,nodes[i]) ;
        int num = 0 ;
        for (iter.First() ; iter.More() ; ++iter) {
            if (iter.CcwElem() == NO_ELEM) {
                num = 0 ;
                break ;
            }
            // ignore the cases with adjacent triangle elements
            // because that means we are near the edge or a messy
            // transition zone and more elements here will likely
            // make things worse.
            int nnum, lnodes[4] ;
            quad_topo->GetElemNodes(iter.CcwElem(),&nnum,lnodes) ;
            if (nnum < 4) {
                num = 0 ;
                break ;
            }
            node_cache.Append(iter.AdjVtx()) ;
            ++num ;
        }

        if (num > 4) {

            // now we have a candidate node for openning,
            // find out how many elements are adjacent to
            // to the adjacent nodes

            List<int> num_adj ;
            List<double> topo_var ;
            List<bool> bdry_flg ;
            double ctv_orig = TopoVariance(num) ;
            for (int ii=0 ; ii<num ; ++ii) {
                double tv ;
                int nadj = quad_topo->NumAdjElems(node_cache[ii]) ;
                if (quad_topo->BoundaryNode(node_cache[ii])) {
                    bdry_flg.Append(true) ;
                    nadj += 2 ;
                    // tv = TopoVariance(nadj) + 0.8 ;
                } else {
                    bdry_flg.Append(false) ;
                    // tv = TopoVariance(nadj) ;
                }
                tv = TopoVariance(nadj) ;
                num_adj.Append(nadj) ;
                topo_var.Append(tv) ;
                ctv_orig += tv ;
            }
            ctv_orig /= num+1 ;

            // here we look for a special case that can cause
            // run away refinement around a hole.  This pattern
            // is a node with 5 adjacent nodes.  4 of the nodes
            // are adjacent to 4 nodes and one is adjancent to 3.

            if (num == 5) {
                int num3 = 0 ;
                int num4 = 0 ;
                for (int k=0 ; k<5 ; ++k) {
                    if (num_adj[k] == 3) {
                        ++num3 ;
                    } else if (num_adj[k] == 4) {
                        ++num4 ;
                    } else {
                        break ;
                    }
                }
                if ((num3 == 1) && (num4 == 4)) continue ;
            }

            // now we look at all the possible split options
            // and rank them

            int top = 0, bot = 0 ;
            double ctv_min = BIG_CTV ;
            if ((num % 2) == 0) {
                for (int j=0 ; j<num ; ++j) {
                    if (bdry_flg[j] || 
                        bdry_flg[(j + num/2) % num]) continue ;
                    double ctv = TopoVariance(num_adj[j]+1) ;
                    for (int k=1 ; k<num ; ++k) {
                        if (k == num/2)
                            ctv += TopoVariance(num_adj[(j+k) % num]+1) ;
                        else
                            ctv += topo_var[(j+k) % num] ;
                    }
                    if (ctv < ctv_min) {
                        top = j ;
                        bot = (j + num/2) % num ;
                        ctv_min = ctv ;
                    }
                }
                ctv_min += 2 * TopoVariance(num-2) ;
                ctv_min /= num+2 ;
            } else {
                for (int j=0 ; j<num ; ++j) {
                    if (bdry_flg[j]) continue ;
                    double ctv0 = TopoVariance(num_adj[j]+1) ;
                    double ctv1 = ctv0 ;
                    for (int k=1 ; k<num ; ++k) {
                        if (k == num/2) {
                            ctv0 += TopoVariance(num_adj[(j+k) % num]+1) ;
                            ctv1 += topo_var[(j+k) % num] ;
                        } if (k == (num/2)+1) {
                            ctv0 += topo_var[(j+k) % num] ;
                            ctv1 += TopoVariance(num_adj[(j+k) % num]+1) ;
                        } else {
                            ctv0 += topo_var[(j+k) % num] ;
                            ctv1 += topo_var[(j+k) % num] ;
                        }
                    }
                    if ((ctv0 < ctv_min) &&
                        !bdry_flg[(j + num/2) % num]) {
                        top = j ;
                        bot = (j + num/2) % num ;
                        ctv_min = ctv0 ;
                    }
                    if ((ctv1 < ctv_min) &&
                        !bdry_flg[(j + (num/2) + 1) % num]) {
                        top = j ;
                        bot = (j + (num/2) + 1) % num ;
                        ctv_min = ctv1 ;
                    }

                }
                ctv_min += TopoVariance(num-2) + TopoVariance(num-1) ;
                ctv_min /= num+2 ;
            }

            // now we do the split

            if (ctv_min < ctv_orig) {
                int new_node = region->DuplicateNode(nodes[i]) ;
                int new_elem = region->NewElemNum() ;
                int nnum, lnodes[4] ;
                int l = top ;
                while (l != bot) {
                    int elem =
                        quad_topo->GetCCWElem(nodes[i],node_cache[l]) ;
                    quad_topo->GetElemNodes(elem,nodes[i],&nnum,lnodes) ;
                    quad_topo->DeleteElement(nnum,lnodes) ;
                    lnodes[0] = new_node ;
                    quad_topo->InsertElement(elem,nnum,lnodes) ;
                    ++l ;
                    if (l >= num) l = 0 ;
                }
                lnodes[0] = node_cache[top] ;
                lnodes[1] = new_node ;
                lnodes[2] = node_cache[bot] ;
                lnodes[3] = nodes[i] ;
                quad_topo->InsertElement(new_elem,4,lnodes) ;
                updates = true ;
            }
        }
    }
    delete [] nodes ;
    return(updates) ;
}




#if 0
/* ------------------------------------------------------------
    DiagSwap - improve the topology with a diagonal swap

    This routine looks for places where the local topological
    connectivity can be improved by a diagonal swap.  For
    example

        #-----*             *-----*   where:
       / \     \           /       \    # indicates nodes connected
      /   \     \         /         \     to more than 4 nodes, and
     o     \     o  ==>  *-----------*
      \     \   /         \         /   o indicates nodes connected
       \     \ /           \       /      to less than 4 nodes.
        *-----#             *-----*

*/

static bool DiagSwap(MshTopo2D *quad_topo)
{
    bool updates = false ;
    int num_nodes = quad_topo->NumNodes() ;
    int *nodes = quad_topo->GetNodeList() ;

    for (int i=0 ; i<num_nodes ; ++i) {

        // note that we only consider nodes with a number
        // greater than the current.  That way each edge
        // is only considered once.

        int num = quad_topo->NumAdjElems(nodes[i]) ;
        if ((num > 4) && !quad_topo->BoundaryNode(nodes[i])) {

            TopoAdjVtxIterator iter(quad_topo,nodes[i]) ;
            for (iter.First() ; iter.More() ; ++iter) {

                int nadj = quad_topo->NumAdjElems(iter.AdjVtx()) ;
                if ((nadj > 4) && !quad_topo->BoundaryNode(iter.AdjVtx())) {

                    int cw0, ccw0, cw1, ccw1 ;
                    int ncw0, nccw0, ncw1, nccw1 ;
                    cw0 = quad_topo->GetCWNode(nodes[i],iter.AdjVtx()) ;       
                    cw1 = quad_topo->GetCWNode(iter.AdjVtx(),nodes[i]) ;       
                    ccw0 = quad_topo->GetCCWNode(nodes[i],iter.AdjVtx()) ;       
                    ccw1 = quad_topo->GetCCWNode(iter.AdjVtx(),nodes[i]) ;       
                    ncw0 = quad_topo->NumAdjElems(cw0) ;
                    ncw1 = quad_topo->NumAdjElems(cw1) ;
                    nccw0 = quad_topo->NumAdjElems(ccw0) ;
                    nccw1 = quad_topo->NumAdjElems(ccw0) ;

                    double ctv_cur_cw = TopoVariance(num) +
                                        TopoVariance(nadj) +
                                        TopoVariance(ncw0) +
                                        TopoVariance(ncw1) ;

                    double ctv_cur_ccw = TopoVariance(num) +
                                         TopoVariance(nadj) +
                                         TopoVariance(nccw0) +
                                         TopoVariance(nccw1) ;

                    double ctv_cw = TopoVariance(num-1) +
                                    TopoVariance(nadj-1) +
                                    TopoVariance(ncw0+1) +
                                    TopoVariance(ncw1+1) ;
                    if (quad_topo->BoundaryNode(cw0) ||
                        quad_topo->BoundaryNode(cw1)) ctv_cw = BIG_CTV ;

                    double ctv_ccw = TopoVariance(num-1) +
                                     TopoVariance(nadj-1) +
                                     TopoVariance(nccw0+1) +
                                     TopoVariance(nccw1+1) ;
                    if (quad_topo->BoundaryNode(ccw0) ||
                        quad_topo->BoundaryNode(ccw1)) ctv_ccw = BIG_CTV ;

                    if ((ctv_cw < ctv_cur_cw) ||
                        (ctv_ccw < ctv_cur_ccw)) {

                        int elem0,elem1,num0,num1,nodes0[4],nodes1[4] ;
                        elem0 = quad_topo->GetCCWElem(nodes[i],iter.AdjVtx()) ;       
                        elem1 = quad_topo->GetCCWElem(iter.AdjVtx(),nodes[i]) ;       
                        quad_topo->GetElemNodes(elem0,nodes[i],&num0,nodes0) ;
                        quad_topo->GetElemNodes(elem1,iter.AdjVtx(),&num1,nodes1) ;

                        if ((num0 != 4) || (num1 != 4)) break ;

                        quad_topo->DeleteElement(num0,nodes0) ;
                        quad_topo->DeleteElement(num1,nodes1) ;

                        if (ctv_cw < ctv_ccw) {
                            nodes0[1] = nodes1[2] ;
                            nodes1[1] = nodes0[2] ;
                        } else {
                            nodes0[0] = nodes1[3] ;
                            nodes1[0] = nodes0[3] ;
                        }
                        quad_topo->InsertElement(elem0,num0,nodes0) ;
                        quad_topo->InsertElement(elem1,num1,nodes1) ;
                        updates = true ;
                        break ;
                    }
                }
            }
        }
    }
    delete [] nodes ;
    return(updates) ;
}

#endif


/* ------------------------------------------------------------
    ThreeEdgeInterior - improve the topology with a diagonal swap

    This routine looks for places in the topology where there
    are two adjacent three-edge interior nodes.  It then makes
    one of the modifications shown below if it will improve
    the topological variance.


        *--------*  
       / \      / \  case 0 
      /   \    /   \   
     *     *--*     *  ==>  *--*-----*--*
      \   /    \   /  
       \ /      \ /  
        *--------* 

        *--------*             *-----*  
       / \      / \  case 1   / \     \  
      /   \    /   \         /   \     \ 
     *     *--*     *  ==>  *     \     *
      \   /    \   /         \     \   / 
       \ /      \ /           \     \ / 
        *--------*             *-----*

        *--------*             *-----*  
       / \      / \  case 2   /     / \  
      /   \    /   \         /     /   \ 
     *     *--*     *  ==>  *     /     *
      \   /    \   /         \   /     / 
       \ /      \ /           \ /     / 
        *--------*             *-----*

*/

static bool ThreeEdgeInterior(MshTopo2D *quad_topo,
                              Dict<int,IntNode> *node_table)
{
    bool updates = false ;
    int num_nodes = quad_topo->NumNodes() ;
    int *nodes = quad_topo->GetNodeList() ;

    for (int i=0 ; i<num_nodes ; ++i) {

        // note that we only consider nodes with a number
        // greater than the current.  That way each edge
        // is only considered once.

        int num = quad_topo->NumAdjElems(nodes[i]) ;
        if ((num == 3) && !quad_topo->BoundaryNode(nodes[i])) {

            TopoAdjVtxIterator iter(quad_topo,nodes[i]) ;
            for (iter.First() ; iter.More() ; ++iter) {

                int nadj = quad_topo->NumAdjElems(iter.AdjVtx()) ;
                if ((nadj == 3) && !quad_topo->BoundaryNode(iter.AdjVtx())) {

                    /*  Nodes and elems are numbered as follows:
                    **
                    **         3         2
                    **        *--------* 
                    **       / \ (0)  / \ 
                    **      /(1)\    /   \ 
                    **   4 *   0 *--* 1   * 7
                    **      \   /    \(3)/ 
                    **       \ / (2)  \ / 
                    **        *--------* 
                    **         5        6
                    */

                    int lnodes[8], ncon[8] ;
                    int elems[4], ennum, enodes[4] ;

                    lnodes[0] = nodes[i] ;
                    lnodes[1] = iter.AdjVtx() ;

                    elems[0] = quad_topo->GetCCWElem(lnodes[0],lnodes[1]) ;
                    if (elems[0] == NO_ELEM) break ;
                    quad_topo->GetElemNodes(elems[0],lnodes[0],&ennum,enodes) ;
                    if (ennum != 4) break ;
                    lnodes[2] = enodes[2] ;
                    lnodes[3] = enodes[3] ;

                    elems[1] = quad_topo->GetCCWElem(lnodes[0],lnodes[3]) ;
                    if (elems[1] == NO_ELEM) break ;
                    quad_topo->GetElemNodes(elems[1],lnodes[0],&ennum,enodes) ;
                    if (ennum != 4) break ;
                    lnodes[4] = enodes[2] ;
                    lnodes[5] = enodes[3] ;

                    elems[2] = quad_topo->GetCCWElem(lnodes[0],lnodes[5]) ;
                    if (elems[2] == NO_ELEM) break ;
                    quad_topo->GetElemNodes(elems[2],lnodes[0],&ennum,enodes) ;
                    if (ennum != 4) break ;
                    lnodes[6] = enodes[2] ;

                    elems[3] = quad_topo->GetCCWElem(lnodes[1],lnodes[6]) ;
                    if (elems[3] == NO_ELEM) break ;
                    quad_topo->GetElemNodes(elems[3],lnodes[1],&ennum,enodes) ;
                    if (ennum != 4) break ;
                    lnodes[7] = enodes[2] ;

                    // find the number connected to each node and
                    // the the current ctv

                    double ctv_cur = 0.0 ;
                    for (int ii=0 ; ii<8 ; ++ii) {
                        ncon[ii] = quad_topo->NumAdjElems(lnodes[ii]) ;
                        ctv_cur += TopoVariance(ncon[ii]) ;
                    }
                    ctv_cur /= 8 ;

                    // find the ctv for the three cases

                    double ctv_0 = (TopoVariance(ncon[4]-1) +
                                    TopoVariance(ncon[7]-1) +
                                    TopoVariance(ncon[3]+ncon[5]-4) +
                                    TopoVariance(ncon[2]+ncon[6]-4)) / 4 ;

                    double ctv_1 = (TopoVariance(ncon[2]-1) +
                                    TopoVariance(ncon[3]) +
                                    TopoVariance(ncon[4]) +
                                    TopoVariance(ncon[5]-1) +
                                    TopoVariance(ncon[6]) +
                                    TopoVariance(ncon[7])) / 6 ; 

                    double ctv_2 = (TopoVariance(ncon[2]) +
                                    TopoVariance(ncon[3]-2) +
                                    TopoVariance(ncon[4]) +
                                    TopoVariance(ncon[5]) +
                                    TopoVariance(ncon[6]-2) +
                                    TopoVariance(ncon[7])) / 6 ;

                    if ((ctv_0 < ctv_cur) || (ctv_1 < ctv_cur) ||
                        (ctv_2 < ctv_cur)) {

                        // now make sure that nodes 4 and 7 are between
                        // the lines 2 -> 3 and 5 -> 6

                        IntNode *nd4 = node_table->Get(lnodes[4]) ; 
                        IntNode *nd7 = node_table->Get(lnodes[7]) ;
 
                        IntNode *nd2 = node_table->Get(lnodes[2]) ; 
                        IntNode *nd3 = node_table->Get(lnodes[3]) ;
                        Vec2D delt = (nd3->coord-nd2->coord).Normalize() ;
                        Vec2D norm = Vec2D(-delt.y(),delt.x()) ;

                        if ((norm * (nd7->coord-nd2->coord)) < 0.0) break ;
                        if ((norm * (nd4->coord-nd2->coord)) < 0.0) break ;

                        IntNode *nd5 = node_table->Get(lnodes[5]) ; 
                        IntNode *nd6 = node_table->Get(lnodes[6]) ; 
                        delt = (nd6->coord-nd5->coord).Normalize() ;
                        norm = Vec2D(-delt.y(),delt.x()) ;

                        if ((norm * (nd7->coord-nd5->coord)) < 0.0) break ;
                        if ((norm * (nd4->coord-nd5->coord)) < 0.0) break ;

                        // case 1 or 2

                        if ((ctv_0 < ctv_1) && (ctv_0 < ctv_2)) {
                            int nds[4] ;
                            nds[0] = lnodes[0] ;
                            nds[1] = lnodes[3] ;
                            nds[2] = lnodes[4] ;
                            nds[3] = lnodes[5] ;
                            lnodes[3] = DoCollapse(quad_topo,nds,node_table,true) ;
                            if (lnodes[3] == -1) break ;
                            nds[0] = lnodes[1] ;
                            nds[1] = lnodes[2] ;
                            nds[2] = lnodes[3] ;
                            nds[3] = lnodes[0] ;
                            quad_topo->DeleteElement(4,nds) ;
                            nds[0] = lnodes[1] ;
                            nds[1] = lnodes[0] ;
                            nds[2] = lnodes[3] ;
                            nds[3] = lnodes[6] ;
                            quad_topo->DeleteElement(4,nds) ;
                            nds[0] = lnodes[1] ;
                            nds[1] = lnodes[6] ;
                            nds[2] = lnodes[7] ;
                            nds[3] = lnodes[2] ;
                            quad_topo->DeleteElement(4,nds) ;
                            nds[0] = lnodes[7] ;
                            nds[1] = lnodes[2] ;
                            nds[2] = lnodes[0] ;
                            nds[3] = lnodes[6] ;
                            quad_topo->InsertElement(elems[0],4,nds) ;
                            DoCollapse(quad_topo,nds,node_table,true) ;
                        } else {

                            // case 0

                            int nds[4] ;
                            nds[0] = lnodes[0] ;
                            nds[1] = lnodes[1] ;
                            nds[2] = lnodes[2] ;
                            nds[3] = lnodes[3] ;
                            quad_topo->DeleteElement(4,nds) ;
                            nds[1] = lnodes[3] ;
                            nds[2] = lnodes[4] ;
                            nds[3] = lnodes[5] ;
                            quad_topo->DeleteElement(4,nds) ;
                            nds[1] = lnodes[5] ;
                            nds[2] = lnodes[6] ;
                            nds[3] = lnodes[1] ;
                            quad_topo->DeleteElement(4,nds) ;
                            nds[0] = lnodes[1] ;
                            nds[1] = lnodes[6] ;
                            nds[2] = lnodes[7] ;
                            nds[3] = lnodes[2] ;
                            quad_topo->DeleteElement(4,nds) ;

                            if (ctv_1 < ctv_2) {
                                nds[0] = lnodes[3] ;
                                nds[1] = lnodes[4] ;
                                nds[2] = lnodes[5] ;
                                nds[3] = lnodes[6] ;
                                quad_topo->InsertElement(elems[0],4,nds) ;
                                nds[0] = lnodes[6] ;
                                nds[1] = lnodes[7] ;
                                nds[2] = lnodes[2] ;
                                nds[3] = lnodes[3] ;
                                quad_topo->InsertElement(elems[1],4,nds) ;
                            } else {
                                nds[0] = lnodes[2] ;
                                nds[1] = lnodes[3] ;
                                nds[2] = lnodes[4] ;
                                nds[3] = lnodes[5] ;
                                quad_topo->InsertElement(elems[0],4,nds) ;
                                nds[0] = lnodes[5] ;
                                nds[1] = lnodes[6] ;
                                nds[2] = lnodes[7] ;
                                nds[3] = lnodes[2] ;
                                quad_topo->InsertElement(elems[1],4,nds) ;
                            }
                        }
                        updates = true ;
                        break ;
                    }
                }
            }
        }
    }
    delete [] nodes ;
    return(updates) ;
}

// #include "MshSmooth2D.hpp"
//  
// static void DebugCleanupSmooth(MshTopo2D *quad_topo,
//                                Dict<int,IntNode> *node_table)
// {
//     Dict<int,MshElement2D> elem_table ;
//     TopoElemIterator iter(quad_topo) ;
//     for (iter.First() ; iter.More() ; ++iter) {
//         int num,nodes[4] ;
//         quad_topo->GetElemNodes(*iter,&num,nodes) ;
//         MshElement2D elem ;
//         elem.elem_id = *iter ;
//         elem.mat_id = 1 ;
//         elem.num_nodes = num ;
//         for (int i=0 ; i<num ; ++i) elem.nodes[i] = nodes[i] ;
//         elem_table.Store(*iter,elem) ;      
//     }
//     Msh2D::MshSmooth2D smooth(node_table,&elem_table) ;
//     // smooth.SmoothNodesWinslow() ;
//     smooth.SmoothNodesLaplace() ;
// }


// %(MshRegion2D::QuadCleanup-void-|-MshTopo2D-|*)
/* ++ ----------------------------------------------------------
**
**    QuadCleanup - do topological cleanups 
**
**      void QuadCleanup(MshTopo2D *quad_topo)
**
**        quad_topo - (in)  quad mesh topology 
**
**      Description: This method performs topological cleanups on a 
**          quadrilateral mesh. 
**
**
** -- */

//extern int debug_visit ;
//extern FILE *dbfd ;

void MshQuad2D::QuadCleanup(MshTopo2D *quad_topo)
{
//if (debug_visit == 55) {
//  if (dbfd == 0) dbfd = fopen("debug.txt","w") ;
//  DebugCleanupSmooth(quad_topo,pnode_table) ;
//  DoPrint(quad_topo,dbfd,false) ;
//  fprintf(dbfd,"a Before any cleanup\n") ;
//  fprintf(dbfd,"n\n") ;
//}


int cnt = 0 ;

    bool updates = true ;
    while (updates) {
        updates = false ;
        ++cnt ;

        updates |= QuadTwoEdge(quad_topo,pnode_table) ;
            //if (dbfd != 0) {
            //    DebugCleanupSmooth(quad_topo,pnode_table) ;
            //    DoPrint(quad_topo,dbfd,false) ;
            //    fprintf(dbfd,"a After Two Edge: %d\n",cnt) ;
            //    fprintf(dbfd,"n\n") ;
            //}
            //if (debug_visit == 22) {
            //  DebugCleanupSmooth(quad_topo,pnode_table) ;
            //  DoPrint(quad_topo,dbfd,false) ;
            //  fprintf(dbfd,"a After QuadTwoEdge\n") ;
            //  fprintf(dbfd,"n\n") ;
            //}

        updates |= QuadCollapse(quad_topo,pnode_table) ;
            //if (dbfd != 0) {
            //    DebugCleanupSmooth(quad_topo,pnode_table) ;
            //    DoPrint(quad_topo,dbfd,false) ;
            //    fprintf(dbfd,"a After Collapse: %d\n",cnt) ;
            //    fprintf(dbfd,"n\n") ;
            //}
            //if (debug_visit == 22) {
            //  DebugCleanupSmooth(quad_topo,pnode_table) ;
            //  DoPrint(quad_topo,dbfd,false) ;
            //  fprintf(dbfd,"a After QuadCollapse\n") ;
            //  fprintf(dbfd,"n\n") ;
            //  DoPrint(quad_topo,dbfd,true) ;
            //  fprintf(dbfd,"a After QuadCollapse\n") ;
            //  fprintf(dbfd,"n\n") ;
            //}

        updates |= QuadOpen(quad_topo,this) ;
            //if (dbfd != 0) {
            //    DebugCleanupSmooth(quad_topo,pnode_table) ;
            //    DoPrint(quad_topo,dbfd,false) ;
            //    fprintf(dbfd,"a After Open: %d\n",cnt) ;
            //    fprintf(dbfd,"n\n") ;
            //}
            //if (debug_visit == 22) {
            //  DebugCleanupSmooth(quad_topo,pnode_table) ;
            //  DoPrint(quad_topo,dbfd,false) ;
            //  fprintf(dbfd,"a After QuadOpen\n") ;
            //  fprintf(dbfd,"n\n") ;
            //  DoPrint(quad_topo,dbfd,true) ;
            //  fprintf(dbfd,"a After QuadOpen\n") ;
            //  fprintf(dbfd,"n\n") ;
            //}

        // diagonal swaps can go on forever so we don't 
        // check the return status on these.

		// Commented out due to propagation bugs, WW 29-Oct-03
        // DiagSwap(quad_topo) ;
            // DebugCleanupSmooth(quad_topo,pnode_table) ;
            // DoPrint(quad_topo) ;
            // printf("a After DiagSwap\n") ;
            // printf("n\n") ;

        updates |= ThreeEdgeInterior(quad_topo,pnode_table) ;
            //if (dbfd != 0) {
            //    DebugCleanupSmooth(quad_topo,pnode_table) ;
            //    DoPrint(quad_topo,dbfd,false) ;
            //    fprintf(dbfd,"a After Three Edge Inter: %d\n",cnt) ;
            //    fprintf(dbfd,"n\n") ;
            //}
            //if (debug_visit == 22) {
            //  DebugCleanupSmooth(quad_topo,pnode_table) ;
            //  DoPrint(quad_topo,dbfd,false) ;
            //  fprintf(dbfd,"a After ThreeEdgeInterior\n") ;
            //  fprintf(dbfd,"n\n") ;
            //}

        //if (updates) DebugCleanupSmooth(quad_topo,pnode_table) ;
            //if (dbfd != 0) {
            //    DoPrint(quad_topo,dbfd,false) ;
            //    fprintf(dbfd,"a After Cleanup Round: %d\n",cnt) ;
            //    fprintf(dbfd,"n\n") ;
            //}
    }
}
    
/* ++ ----------------------------------------------------------
**
**    QuadAngleCheck - do quadrilateral angle checks 
**
**      void QuadAngleChecks(MshTopo2D *quad_topo)
**
**        quad_topo - (in)  quad mesh topology 
**
**      Description: Look for quadrilateral elements that have
**          angle that are out of limits and replace them with
**          two triangles. 
**
**
** -- */

#define PI 3.14159265359

void MshQuad2D::QuadAngleChecks()
{
    // loop through all elements and find the nodes

    List<int> invalid ;
    Dict<int,MshElement2D>::DictIterator elems(pelem_table) ;

    for (elems.First() ; elems.More() ; ++elems) {

        int i ;
        MshElement2D& elem = elems.Entry() ;
        int num_node = ((elem.num_nodes == 3) ||
                        (elem.num_nodes == 6)) ? 3 : 4 ;
        if (num_node != 4) continue ;

        // get the nodes coordinates

        IntNode *nd[4] ;
        for (i=0 ; i<4 ; ++i) {
            nd[i] = pnode_table->Get(elem.nodes[i]) ;
        } 

        // if we are doing node checks, loop through the nodes
        // and check the angles.

        if (QuadAngleLocation == AT_NODES) {
            for (i=0 ; i<4 ; ++i) {
                int j = i==3 ? 0 : i+1 ;
                int k = i==0 ? 3 : i-1 ;
                double angle = Angle2Pi(nd[i]->coord,nd[j]->coord,
                                        nd[k]->coord) ;
                if ((angle < MinQuadAngle) || (angle > MaxQuadAngle)) {
                    invalid.Append(elem.elem_id) ;
                    break ;
                }
            }
        } else {

            // compute the components of the jacobian at the guass
            // points (many values precomputed).

            double val0 = 0.211325 ;
            double val1 = 0.788675 ;

            double dxdr[4], dxds[4], dydr[4], dyds[4] ;
            dxdr[0] = -val1*nd[0]->coord[0] + val1*nd[1]->coord[0] +
                       val0*nd[2]->coord[0] - val0*nd[3]->coord[0] ;
            dydr[0] = -val1*nd[0]->coord[1] + val1*nd[1]->coord[1] +
                       val0*nd[2]->coord[1] - val0*nd[3]->coord[1] ;
            dxds[0] = -val1*nd[0]->coord[0] - val0*nd[1]->coord[0] +
                       val0*nd[2]->coord[0] + val1*nd[3]->coord[0] ;
            dyds[0] = -val1*nd[0]->coord[1] - val0*nd[1]->coord[1] +
                       val0*nd[2]->coord[1] + val1*nd[3]->coord[1] ;

            dxdr[1] = dxdr[0] ;
            dydr[1] = dydr[0] ;
            dxds[1] = -val0*nd[0]->coord[0] - val1*nd[1]->coord[0] +
                       val1*nd[2]->coord[0] + val0*nd[3]->coord[0] ;
            dyds[1] = -val0*nd[0]->coord[1] - val1*nd[1]->coord[1] +
                       val1*nd[2]->coord[1] + val0*nd[3]->coord[1] ;

            dxdr[2] = -val0*nd[0]->coord[0] + val0*nd[1]->coord[0] +
                       val1*nd[2]->coord[0] - val1*nd[3]->coord[0] ;
            dydr[2] = -val0*nd[0]->coord[1] + val0*nd[1]->coord[1] +
                       val1*nd[2]->coord[1] - val1*nd[3]->coord[1] ;
            dxds[2] = dxds[1] ;
            dyds[2] = dyds[1] ;

            dxdr[3] = dxdr[2] ;
            dydr[3] = dydr[2] ;
            dxds[3] = dxds[0] ;
            dyds[3] = dyds[0] ;

            for (i=0 ; i<4 ; ++i) {
                double dot = dxdr[i]*dxds[i] + dydr[i]*dyds[i] ;
                double mag0 = sqrt(dxdr[i]*dxdr[i] + dydr[i]*dydr[i]) ;
                double mag1 = sqrt(dxds[i]*dxds[i] + dyds[i]*dyds[i]) ;
                double cos_angle = dot / (mag0*mag1) ;
                if (cos_angle > 1.0) cos_angle = 1.0 ;
                if (cos_angle < -1.0) cos_angle = -1.0 ;
                double angle = acos(cos_angle) ;

                if ((angle < MinQuadAngle) || (angle > MaxQuadAngle)) {
                    invalid.Append(elem.elem_id) ;
                    break ;
                }

                angle = PI - angle ;
                if ((angle < MinQuadAngle) || (angle > MaxQuadAngle)) {
                    invalid.Append(elem.elem_id) ;
                    break ;
                }
            }
        }
    }

    // loop through the invalid elements and replace with triangles.

    int i ;
    for (i=0 ; i<invalid.Len() ; ++i) {

        MshElement2D elem = *(pelem_table->Get(invalid[i])) ;
        pelem_table->Del(invalid[i]) ;

        IntNode *nd[4] ;
        for (int j=0 ; j<4 ; ++j) {
            nd[j] = pnode_table->Get(elem.nodes[j]) ;
        } 

        // determine which of the two possible triangular
        // patterns to add

        double sm0 = TriMetric(nd[0]->coord,nd[1]->coord,
                               nd[2]->coord ) ;
        double sm1 = TriMetric(nd[0]->coord,nd[2]->coord,
                               nd[3]->coord ) ;
        double sm_case0 = sm0 < sm1 ? sm0 : sm1 ;

        double sm2 = TriMetric(nd[0]->coord,nd[1]->coord,
                               nd[3]->coord ) ;
        double sm3 = TriMetric(nd[3]->coord,nd[1]->coord,
                               nd[2]->coord ) ;
        double sm_case1 = sm2 < sm3 ? sm2 : sm3 ;

        MshElement2D new_elem ;
        new_elem.mat_id = elem.mat_id ;

        if (sm_case1 < sm_case0) {
            if (elem.num_nodes == 4) {

                new_elem.elem_id = elem.elem_id ;
                new_elem.num_nodes = 3 ;
                new_elem.nodes[0] = elem.nodes[0] ;
                new_elem.nodes[1] = elem.nodes[1] ;
                new_elem.nodes[2] = elem.nodes[2] ;
                pelem_table->Store(new_elem.elem_id,new_elem) ;
                
                new_elem.elem_id = NewElemNum() ;
                new_elem.num_nodes = 3 ;
                new_elem.nodes[0] = elem.nodes[0] ;
                new_elem.nodes[1] = elem.nodes[2] ;
                new_elem.nodes[2] = elem.nodes[3] ;
                pelem_table->Store(new_elem.elem_id,new_elem) ;

            } else {

                Vec2D pos = 0.5*(nd[0]->coord + nd[2]->coord) ;
                int new_node = NewNode(pos.x(),pos.y(),
                                       INTERIOR,MSH_FLOATING,false) ;
                new_elem.elem_id = elem.elem_id ;
                new_elem.num_nodes = 6 ;
                new_elem.nodes[0] = elem.nodes[0] ;
                new_elem.nodes[1] = elem.nodes[1] ;
                new_elem.nodes[2] = elem.nodes[2] ;
                new_elem.nodes[3] = elem.nodes[4] ;
                new_elem.nodes[4] = elem.nodes[5] ;
                new_elem.nodes[5] = new_node ;
                pelem_table->Store(new_elem.elem_id,new_elem) ;
                
                new_elem.elem_id = NewElemNum() ;
                new_elem.num_nodes = 6 ;
                new_elem.nodes[0] = elem.nodes[0] ;
                new_elem.nodes[1] = elem.nodes[2] ;
                new_elem.nodes[2] = elem.nodes[3] ;
                new_elem.nodes[3] = new_node ;
                new_elem.nodes[4] = elem.nodes[6] ;
                new_elem.nodes[5] = elem.nodes[7] ;
                pelem_table->Store(new_elem.elem_id,new_elem) ;

            }
        } else {
            if (elem.num_nodes == 4) {

                new_elem.elem_id = elem.elem_id ;
                new_elem.num_nodes = 3 ;
                new_elem.nodes[0] = elem.nodes[0] ;
                new_elem.nodes[1] = elem.nodes[1] ;
                new_elem.nodes[2] = elem.nodes[3] ;
                pelem_table->Store(new_elem.elem_id,new_elem) ;
                
                new_elem.elem_id = NewElemNum() ;
                new_elem.num_nodes = 3 ;
                new_elem.nodes[0] = elem.nodes[3] ;
                new_elem.nodes[1] = elem.nodes[1] ;
                new_elem.nodes[2] = elem.nodes[2] ;
                pelem_table->Store(new_elem.elem_id,new_elem) ;

            } else {

                Vec2D pos = 0.5*(nd[1]->coord + nd[3]->coord) ;
                int new_node = NewNode(pos.x(),pos.y(),
                                       INTERIOR,MSH_FLOATING,false) ;
                new_elem.elem_id = elem.elem_id ;
                new_elem.num_nodes = 6 ;
                new_elem.nodes[0] = elem.nodes[0] ;
                new_elem.nodes[1] = elem.nodes[1] ;
                new_elem.nodes[2] = elem.nodes[3] ;
                new_elem.nodes[3] = elem.nodes[4] ;
                new_elem.nodes[4] = new_node ;
                new_elem.nodes[5] = elem.nodes[7] ;
                pelem_table->Store(new_elem.elem_id,new_elem) ;
                
                new_elem.elem_id = NewElemNum() ;
                new_elem.num_nodes = 6 ;
                new_elem.nodes[0] = elem.nodes[3] ;
                new_elem.nodes[1] = elem.nodes[1] ;
                new_elem.nodes[2] = elem.nodes[2] ;
                new_elem.nodes[3] = new_node ;
                new_elem.nodes[4] = elem.nodes[5] ;
                new_elem.nodes[5] = elem.nodes[6] ;
                pelem_table->Store(new_elem.elem_id,new_elem) ;

            }
        }
    }
}

} // namespace



