//
// CArbMshTopo2D Class definition
//
// Description -
//   This class implements an object that maintians a dynamic
//   list of ordered lists of verticies adjacent to all
//   verticies in a mesh.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.19 $  $Date: 2002/07/26 17:56:21 $  $Author: wash $
//

#include <stdio.h>
#include <math.h>
#include "ArbMshTopo2D.hpp"
#include "assert.h"

#ifdef MEMDEBUG
#include "MemDbg.hpp"
#define new new(__FILE__,__LINE__)
#endif


// %(CArbMshTopo2D::CArbMshTopo2D-constructor-|) 
/* ++ ----------------------------------------------------------
**
**    CArbMshTopo2D - constructor 
**
**      CArbMshTopo2D()
**
**      Description: This is the ArbMshTopo2D constructor. 
**
**
** -- */

CArbMshTopo2D::CArbMshTopo2D()
{
    vtx_table = new CArbHashTable <int, ArbVtxEntry*> ( ) ;
//                       CArbMshTopo2D::ArbVtxEntry*>() ;
    elem_table = new CArbHashTable <int, int> ( ) ;

    MaxAddress = 0 ;
    MaxVtx = 0 ;
    FreeList = 0 ;
    CacheList = 0 ;
}

CArbMshTopo2D::CArbMshTopo2D(const CArbMshTopo2D &other)
{
    FreeList = 0 ;
    CacheList = 0 ;
    vtx_table = new CArbHashTable <int, ArbVtxEntry*> ( );
    CArbHashTableIterator <int, ArbVtxEntry*> viter(other.vtx_table);
    for (viter.First() ; viter.More() ; ++viter)
    {
        int vert = viter.Key ( );
        ArbVtxEntry *tmp_list = 0;
        ArbVtxEntry **next = &tmp_list;
        CArbTopoAdjVtxIterator aiter(&other,vert);
        for (aiter.First() ; aiter.More() ; ++aiter)
        {
            ArbVtxEntry *entry ;
            entry = NewVtxEntry() ;
            entry->vtx_num = aiter.AdjVtx() ;
            entry->ccw_elem = aiter.CcwElem() ;
            entry->next = 0 ;
            *next = entry ;
            next = &entry->next ;
        }
        vtx_table->Store(vert,tmp_list) ;
    }

    elem_table = new CArbHashTable<int,int>(*other.elem_table) ;

    MaxAddress = other.MaxAddress ;
    MaxVtx = other.MaxVtx ;
}

CArbMshTopo2D CArbMshTopo2D::operator = (const CArbMshTopo2D &other)
{
    DeleteEntryCache() ;
    delete vtx_table ;
    delete elem_table ;

    FreeList = 0 ;
    CacheList = 0 ;
    vtx_table = new CArbHashTable<int,ArbVtxEntry*>() ;
    CArbHashTableIterator<int,ArbVtxEntry*> viter(other.vtx_table) ;
    for (viter.First() ; viter.More() ; ++viter) {
        int vert = viter.Key() ;
        ArbVtxEntry *tmp_list = 0 ;
        ArbVtxEntry **next = &tmp_list ;
        CArbTopoAdjVtxIterator aiter(&other,vert) ;
        for (aiter.First() ; aiter.More() ; ++aiter) {
            ArbVtxEntry *entry ;
            entry = NewVtxEntry() ;
            entry->vtx_num = aiter.AdjVtx() ;
            entry->ccw_elem = aiter.CcwElem() ;
            entry->next = 0 ;
            *next = entry ;
            next = &entry->next ;
        }
        vtx_table->Store(vert,tmp_list) ;
    }

    elem_table = new CArbHashTable<int,int>(*other.elem_table) ;

    MaxAddress = other.MaxAddress ;
    MaxVtx = other.MaxVtx ;
    return(*this) ;
}

// %(CArbMshTopo2D::CArbMshTopo2D-destructor-|~)
/* ++ ----------------------------------------------------------
**
**    CArbMshTopo2D - destructor 
**
**      ~CArbMshTopo2D()
**
**      Description: THis is a destructor for ArbMshTopo2D objects. 
**
**
** -- */

CArbMshTopo2D::~CArbMshTopo2D()
{
    DeleteEntryCache() ;

//     // get a list of all the hash table entries
// 
//     CArbHashTableIterator<int,ArbVtxEntry *> iter(vtx_table) ;
// 
//     // delete all entries in all lists
// 
//     for (iter.First() ; iter.More() ; ++iter) {
//         ArbVtxEntry *ptr = *(iter.Entry()) ;
//         while (ptr != 0) {
//             ArbVtxEntry *next = ptr->next ;
//             delete ptr ;
//             ptr = next ;
//         }
//     }

    // now delete the hash table

    delete elem_table ;
    delete vtx_table ;
}


CArbMshTopo2D::ArbVtxEntry *CArbMshTopo2D::NewVtxEntry()
{
    if (FreeList == 0) {
        ArbVtxCache *ctmp = new ArbVtxCache ;
        ctmp->next = CacheList ;
        CacheList = ctmp ;
        for (int i=0 ; i<TOPO_CACHE_BLOCK_SIZE-1 ; ++i) {
            ctmp->entries[i].next = &(ctmp->entries[i+1]) ;
        }
        ctmp->entries[TOPO_CACHE_BLOCK_SIZE-1].next = 0 ;
        FreeList = &(ctmp->entries[0]) ;
    }
    ArbVtxEntry *tmp = FreeList ;
    FreeList = tmp->next ;
    return(tmp) ;
}


void CArbMshTopo2D::DeleteVtxEntry(ArbVtxEntry *entry)
{
    entry->next = FreeList ;
    FreeList = entry ;
}


void CArbMshTopo2D::DeleteEntryCache()
{
    ArbVtxCache *ptr = CacheList ;
    while (ptr != 0) {
        ArbVtxCache *tmp = ptr->next ;
        delete ptr ;
        ptr = tmp ;
    }
}

 
// %(CArbMshTopo2D::InsertElement-void-|-int-const|-int-const|-int-const|*const)
/* ++ ----------------------------------------------------------
**
**    InsertElement - insert an element 
**
**      void InsertElement(
**              const int elem_id,
**              const int num_nodes,
**              const int *constnodes)
**
**        elem_id   - (in)  element id 
**        num_nodes - (in)  number of nodes 
**        nodes     - (in)  node id's 
**
**      Description: This function inserts an element into the current 
**          mesh. 
**
**
** -- */
void CArbMshTopo2D::InsertElement(const int elem_id,
                                  const int num_nodes,
                                  const int *const nodes)
{
  int i, j, k ;
  for (i=0 ; i<num_nodes; ++i)
  {
    j = (i+1) % num_nodes;
    k = (i+2) % num_nodes;
    InsertAngle(elem_id,nodes[j],nodes[k],nodes[i]);
  }

  if (elem_table->Fetch (elem_id) == 0)
    elem_table->Store (elem_id, nodes[0]);
}




// %(CArbMshTopo2D::InsertCollapsedElement-void-|-int-const|-int-const|-int-const|*const)
/* ++ ----------------------------------------------------------
**
**    InsertCollapsedElement - insert a collapsed element 
**
**      void InsertCollapsedElement(
**              const int elem_id,
**              const int num_nodes,
**              const int *constnodes)
**
**        elem_id   - (in)  element id 
**        num_nodes - (in)  number of nodes 
**        nodes     - (in)  node id's 
**
**      Description: This function inserts an element for cases where 
**          element id's may be repeated. 
**
**
** -- */

void CArbMshTopo2D::InsertCollapsedElement(const int elem_id,
                                           const int num_nodes,
                                           const int *const nodes,
                                           const int parent_elem)
{
    int i, j, k ;
    for (i=0 ; i<num_nodes ; ++i) {
        j = (i+1) % num_nodes ;
        k = (i+2) % num_nodes ;
        if (nodes[i] != nodes[j]) {
            if (nodes[j] == nodes[k]) {
                int l = 3 ;
                k = (i+l) % num_nodes ;
                while (nodes[j] == nodes[k]) {
                    ++l ;
                    k = (i+l) % num_nodes ;
                }
            }
            InsertAngle(elem_id,nodes[j],nodes[k],nodes[i],
                        parent_elem) ;
        }
    }

    if (elem_table->Fetch(elem_id) == 0) {
        elem_table->Store(elem_id,nodes[0]) ;
    }
}




// %(CArbMshTopo2D::InsertTriangle-void-|-int-const|-int-const|-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    InsertTriangle - insert a triangular element 
**
**      void InsertTriangle(
**              const int elem_id,
**              const int nd0,
**              const int nd1,
**              const int nd2)
**
**        elem_id - (in)  element id 
**        nd0     - (in)  first node 
**        nd1     - (in)  second node 
**        nd2     - (in)  third node 
**
**      Description: This function inserts a triangular element into 
**          the current mesh. 
**
**
** -- */

int CArbMshTopo2D::InsertTriangle(const int elem_id,
                                   const int nd0,
                                   const int nd1,
                                   const int nd2)
{
    if (!InsertAngle(elem_id,nd0,nd1,nd2))
     return 0;
    if (!InsertAngle(elem_id,nd1,nd2,nd0))
     return 0;
    if (!InsertAngle(elem_id,nd2,nd0,nd1))
     return 0;
    if (elem_table->Fetch(elem_id) == 0) 
    {
        elem_table->Store(elem_id,nd0) ;
    }
    return 0;
}




// %(CArbMshTopo2D::DeleteElement-void-|-int-const|-int-const|*const)
/* ++ ----------------------------------------------------------
**
**    DeleteElement - delete an element 
**
**      void DeleteElement(
**              const int num_nodes,
**              const int *constnodes)
**
**        num_nodes - (in)  number of nodes 
**        nodes     - (in)  node id's 
**
**      Description: This function deletes an element from the mesh 
**          when the node id's are known. 
**
**
** -- */

int CArbMshTopo2D::DeleteElement(const int num_nodes,
                                  const int *const nodes)
{
    int i, j, k ;
    int elem ;
    for (i=0 ; i<num_nodes ; ++i) {
        j = (i+1) % num_nodes ;
        k = (i+2) % num_nodes ;
        elem = DeleteAngle(nodes[j],nodes[k],nodes[i]) ;
        if (elem == -1)
         return 0;
    }
    elem_table->Remove(elem) ;
    return 1;
}




// %(CArbMshTopo2D::DeleteElement-void-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    DeleteElement - delete an element 
**
**      void DeleteElement(const int elem_id)
**
**        elem_id - (in)  element id 
**
**      Description: This function deletes an element from the mesh 
**          when the element id is known. 
**
**
** -- */

#define SMALL_NUM 8

void CArbMshTopo2D::DeleteElement(const int elem)
{
    int num, small_nodes[SMALL_NUM] ;

    // get a pointer to the first node

    int *node_ptr = elem_table->Fetch(elem) ;
    assert(node_ptr != 0) ;

    int start_node ;
    small_nodes[0] = start_node = *node_ptr ;
    num = 1 ;

    // Now we traverse around the element assuming that
    // we will fit in the small nodes array

    bool too_big = false ;
    int adj_node = GetAdjVtxAlongElem(start_node,elem) ;
    while(adj_node != start_node) {
        if (num == SMALL_NUM) {
            too_big = true ;
            break ;
        }
        small_nodes[num] = adj_node ;
        ++num ;
        adj_node = GetAdjVtxAlongElem(adj_node,elem) ;
    }

    // if we get here and the too_big flag is not set,
    // delete the element

    if (!too_big) {
        DeleteElement(num,small_nodes) ;
    } else {
        // if we are too big, do the memory dynamically

        num = NumElemNodes(elem) ;
        int *nodes = new int[num] ;
        GetElemNodes(elem,&num,nodes) ;
        DeleteElement(num,nodes) ;
        delete [] nodes ;
    }
}




// %(CArbMshTopo2D::DeleteTriangle-void-|-int-const|-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    DeleteTriangle - delete a triangle 
**
**      void DeleteTriangle(
**              const int nd0,
**              const int nd1,
**              const int nd2)
**
**        nd0 - (in)  first node 
**        nd1 - (in)  second node 
**        nd2 - (in)  third node 
**
**      Description: Delete a triangular element from a mesh when the 
**          node id's are known. 
**
**
** -- */

int CArbMshTopo2D::DeleteTriangle(const int nd0,
                                   const int nd1,
                                   const int nd2)
{
    int elem ;
    elem = DeleteAngle(nd0,nd1,nd2) ;
    if (elem == -1)
     return 0;
    elem = DeleteAngle(nd1,nd2,nd0) ;
    if (elem == -1)
     return 0;
    elem = DeleteAngle(nd2,nd0,nd1) ;
    if (elem == -1)
     return 0;
    elem_table->Remove(elem) ;
    return 1;
}


void CArbMshTopo2D::SplitEdgeHelp(ArbVtxEntry *ptr0,
                                  ArbVtxEntry *ptr1,
                                  const int nd0,
                                  const int nd1,
                                  const int new_node)
{
    ArbVtxEntry **list ;

    if (nd0 != nd1) {

        // update the pointers

        ptr0->vtx_num = new_node ;
        ptr1->vtx_num = new_node ;

        // add the information for the new node

        vtx_table->Store(new_node,0) ;
        list = vtx_table->Fetch(new_node) ;

        ArbVtxEntry *entry0 = NewVtxEntry() ;
        entry0->vtx_num = nd1 ;
        entry0->next = 0 ;
        entry0->ccw_elem = ptr0->ccw_elem ;

        ArbVtxEntry *entry1 = NewVtxEntry() ;
        entry1->vtx_num = nd0 ;
        entry1->next = 0 ;
        entry1->ccw_elem = ptr1->ccw_elem ;

        if (ptr1->ccw_elem == NO_ELEM) {
            entry0->next = entry1 ;
            *list = entry0 ;
        } else {
            entry1->next = entry0 ;
            *list = entry1 ;
        }

    } else {

        // update the pointers

        ptr0->vtx_num = new_node ;

        // add the information for the new node

        vtx_table->Store(new_node,0) ;
        list = vtx_table->Fetch(new_node) ;

        ArbVtxEntry *entry0 = NewVtxEntry() ;
        entry0->vtx_num = nd0 ;
        entry0->next = 0 ;
        entry0->ccw_elem = ptr0->ccw_elem ;

        *list = entry0 ;
    }
}


void CArbMshTopo2D::SplitEdge(const int nd0,
                              const int nd1,
                              const int new_node)
{
    // new node should not exits

#ifdef _DEBUG
    ArbVtxEntry **list = vtx_table->Fetch(new_node) ;
    assert(list == 0) ;
#endif

    // make sure that the nodes exist and share and edge

    ArbVtxEntry **list0, *ptr0 ;
    ArbVtxEntry **list1, *ptr1 ;

    list0 = vtx_table->Fetch(nd0) ;
    assert(list0 != 0) ;

    for (ptr0 = *list0 ; ptr0 != 0 ; ptr0 = ptr0->next) {
        if (ptr0->vtx_num == nd1) break ;
    }
    assert(ptr0 != 0) ;

    list1 = vtx_table->Fetch(nd1) ;
    assert(list1 != 0) ;

    for (ptr1 = *list1 ; ptr1 != 0 ; ptr1 = ptr1->next) {
        if (ptr1->vtx_num == nd0) break ;
    }
    assert(ptr1 != 0) ;

    SplitEdgeHelp(ptr0,ptr1,nd0,nd1,new_node) ;
}

void CArbMshTopo2D::SplitEdgeElem(const int nd0,
                                  const int nd1,
                                  const int new_node,
                                  const int elem,
                                  const bool Ccw_flag)
{
    // this is a version of split edge that also takes an
    // element id as an argument.  This is used for the special
    // case that there may be more than one edge joining
    // nd0 and nd1

    // if Ccw_flag is set to true, then the edge to be split will
    // be the one that, haveing nd0 and the start and nd1 as the
    // adj vertex, has elem as it's ccw_elem.

    // new node should not exits

#ifdef _DEBUG
    ArbVtxEntry **list = vtx_table->Fetch(new_node) ;
    assert(list == 0) ;
#endif
    
    // look through the nodes adjacent to the first node
    // and find the 1 or 2 edges that connect to the second node

    ArbVtxEntry **list0, *ptr0, *ptr0_1 = 0, *ptr0_2 = 0 ;

    list0 = vtx_table->Fetch(nd0) ;
    assert(list0 != 0) ;

    for (ptr0 = *list0 ; ptr0 != 0 ; ptr0 = ptr0->next) {
        if (ptr0->vtx_num == nd1) {
            if (ptr0_1 == 0) {
                ptr0_1 = ptr0 ;
            } else {
                ptr0_2 = ptr0 ;
                break ;
            }
        }
    }
    assert(ptr0_1 != 0) ;

    // do the same for node 2 relative to node 1

    ArbVtxEntry **list1, *ptr1, *ptr1_1 = 0, *ptr1_2 = 0 ;

    list1 = vtx_table->Fetch(nd1) ;
    assert(list1 != 0) ;

    for (ptr1 = *list1 ; ptr1 != 0 ; ptr1 = ptr1->next) {
        if (ptr1->vtx_num == nd0) {
            if (ptr1_1 == 0) {
                ptr1_1 = ptr1 ;
            } else {
                ptr1_2 = ptr1 ;
            }
        }
    }
    assert(ptr1_1 != 0) ;

    // check to make sure that the two nodes agree on the number
    // of edges joining these nodes

    assert((ptr0_2 == 0) ? (ptr1_2 == 0) : 1) ;
    assert((ptr1_2 == 0) ? (ptr0_2 == 0) : 1) ;

    // figure out which ptrs to use to split the edge.

    if (ptr1_2 == 0) {
        SplitEdgeHelp(ptr0_1,ptr1_1,nd0,nd1,new_node) ;
    } else {
        if (Ccw_flag) {
            if (ptr0_1->ccw_elem == elem) {
                if (GetCWElemPtr(ptr1_1,list1) == elem) {
                    SplitEdgeHelp(ptr0_1,ptr1_1,nd0,nd1,new_node) ;
                } else if (GetCWElemPtr(ptr1_2,list1) == elem) {
                    SplitEdgeHelp(ptr0_1,ptr1_2,nd0,nd1,new_node) ;
                } else {
                    assert(0) ;
                }
            } else if ((ptr0_2 != 0) && (ptr0_2->ccw_elem == elem)) {
                if (GetCWElemPtr(ptr1_1,list1) == elem) {
                    SplitEdgeHelp(ptr0_2,ptr1_1,nd0,nd1,new_node) ;
                } else if (GetCWElemPtr(ptr1_2,list1) == elem) {
                    SplitEdgeHelp(ptr0_2,ptr1_2,nd0,nd1,new_node) ;
                } else {
                    assert(0) ;
                }
            } else {
                assert(0) ;
            }
        } else {
            if (GetCWElemPtr(ptr0_1,list0) == elem) {
                if (ptr1_1->ccw_elem == elem) {
                    SplitEdgeHelp(ptr0_1,ptr1_1,nd0,nd1,new_node) ;
                } else if (ptr1_2->ccw_elem == elem) {
                    SplitEdgeHelp(ptr0_1,ptr1_2,nd0,nd1,new_node) ;
                } else {
                    assert(0) ;
                }
            } else if ((ptr0_2 != 0) && GetCWElemPtr(ptr0_2,list0) == elem) {
                if (ptr1_1->ccw_elem == elem) {
                    SplitEdgeHelp(ptr0_1,ptr1_1,nd0,nd1,new_node) ;
                } else if (ptr1_2->ccw_elem == elem) {
                    SplitEdgeHelp(ptr1_1,ptr1_2,nd0,nd1,new_node) ;
                } else {
                    assert(0) ;
                }
            } else {
                assert(0) ;
            }
        }
    }
}


void CArbMshTopo2D::AddEdgeVertex(int vtx,int CW_vtx,int new_vtx)
{
    // This vtx should exist and be adjacent to the CW_vtx and
    // the CW_vtx should not be adjacent to a null element

    ArbVtxEntry **list = vtx_table->Fetch(vtx) ;
    assert(list != 0) ;

    ArbVtxEntry *ptr ;

    for (ptr = *list ; ptr != 0 ; ptr = ptr->next) {
        if (ptr->vtx_num == CW_vtx) break ;
    }
    assert(ptr != 0) ;
    assert(ptr->ccw_elem != NO_ELEM) ;

    // Get the next adjacent going in a CCW direction

    int CCW_vtx = ptr->next == 0 ? (*list)->vtx_num :
                                   ptr->next->vtx_num ;
    int elem = ptr->ccw_elem ;

    // delete the existing angle and add the two new ones

    DeleteAngle(vtx,CW_vtx,CCW_vtx) ;
    InsertAngle(elem,vtx,CW_vtx,new_vtx) ;
    InsertAngle(elem,vtx,new_vtx,CCW_vtx) ;

    InsertAngle(elem,new_vtx,vtx,vtx) ;
}

void CArbMshTopo2D::AddEdgeTwoVertex(int vtx0,int vtx1,int elem)
{
    // neither vertex should exist

    ArbVtxEntry **list0, **list1 ;

    list0 = vtx_table->Fetch(vtx0) ;
    list1 = vtx_table->Fetch(vtx1) ;
    assert(list0 == 0) ;
    assert(list1 == 0) ;

    // add the verts and angles

    InsertAngle(elem,vtx0,vtx1,vtx1) ;
    InsertAngle(elem,vtx1,vtx0,vtx0) ;

    // if necessary add the element

    if (elem_table->Fetch(elem) == 0)
        elem_table->Store(elem,vtx0) ;
}

int CArbMshTopo2D::AddFaceEdge(int vtx0,int vtx1,
                               int CW_vtx0,int CW_vtx1,
                               int new_elem)
{
    // both verticies should exist and be adjacent to their CW
    // neighbors.  Also, both CW verts should be adjacent to
    // the same element

    ArbVtxEntry **list0, *ptr0 ;
    ArbVtxEntry **list1, *ptr1 ;

    list0 = vtx_table->Fetch(vtx0) ;
    assert(list0 != 0) ;
    for (ptr0 = *list0 ; ptr0 != 0 ; ptr0 = ptr0->next) {
        if (ptr0->vtx_num == CW_vtx0) break ;
    }
    assert(ptr0 != 0) ;
//    assert(ptr0->ccw_elem != NO_ELEM) ;

    list1 = vtx_table->Fetch(vtx1) ;
    assert(list1 != 0) ;
    for (ptr1 = *list1 ; ptr1 != 0 ; ptr1 = ptr1->next) {
        if (ptr1->vtx_num == CW_vtx1) break ;
    }
    assert(ptr1 != 0) ;
//    assert(ptr1->ccw_elem != NO_ELEM) ;

    assert(ptr0->ccw_elem == ptr1->ccw_elem) ;
    int elem = ptr0->ccw_elem ;

    // update the angles

    int CCW_vtx0 = ptr0->next == 0 ? (*list0)->vtx_num :
                                     ptr0->next->vtx_num ;
    DeleteAngle(vtx0,CW_vtx0,CCW_vtx0) ;
    InsertAngle(elem,vtx0,CW_vtx0,vtx1) ;
    InsertAngle(new_elem,vtx0,vtx1,CCW_vtx0) ;

    int CCW_vtx1 = ptr1->next == 0 ? (*list1)->vtx_num :
                                     ptr1->next->vtx_num ;
    DeleteAngle(vtx1,CW_vtx1,CCW_vtx1) ;
    InsertAngle(elem,vtx1,vtx0,CCW_vtx1) ;
    InsertAngle(new_elem,vtx1,CW_vtx1,vtx0) ;

    // now loop around and update the element id's

    list0 = vtx_table->Fetch(vtx0) ;
    for (ptr0 = *list0 ; ptr0 != 0 ; ptr0 = ptr0->next) {
        if (ptr0->vtx_num == vtx1) break ;
    }
    ptr0->ccw_elem = new_elem ;

    int cur = vtx0 ;
    int nxt = ptr0->next != 0 ? ptr0->next->vtx_num :
                                (*list0)->vtx_num ;

    while (nxt != vtx0) {
        list0 = vtx_table->Fetch(nxt) ;
        for (ptr0 = *list0 ; ptr0 != 0 ; ptr0 = ptr0->next) {
            if (ptr0->vtx_num == cur) break ;
        }
        ptr0->ccw_elem = new_elem ;
        cur = nxt ;
        nxt = ptr0->next != 0 ? ptr0->next->vtx_num :
                                (*list0)->vtx_num ;
    }

    elem_table->Store(elem,vtx0) ;
    elem_table->Store(new_elem,vtx0) ;

    return(elem) ;
}

int CArbMshTopo2D::FindCrossedFace(int vtx0,int vtx1,
                                   int CW_vtx0,int CW_vtx1) const
{
    // both verticies should exist and be adjacent to their CW
    // neighbors.  Also, both CW verts should be adjacent to
    // the same element

    ArbVtxEntry **list0, *ptr0 ;
    ArbVtxEntry **list1, *ptr1 ;

    list0 = vtx_table->Fetch(vtx0) ;
    assert(list0 != 0) ;
    for (ptr0 = *list0 ; ptr0 != 0 ; ptr0 = ptr0->next) {
        if (ptr0->vtx_num == CW_vtx0) break ;
    }
    assert(ptr0 != 0) ;
//    assert(ptr0->ccw_elem != NO_ELEM) ;

    list1 = vtx_table->Fetch(vtx1) ;
    assert(list1 != 0) ;
    for (ptr1 = *list1 ; ptr1 != 0 ; ptr1 = ptr1->next) {
        if (ptr1->vtx_num == CW_vtx1) break ;
    }
    assert(ptr1 != 0) ;
//    assert(ptr1->ccw_elem != NO_ELEM) ;

    assert(ptr0->ccw_elem == ptr1->ccw_elem) ;
    return(ptr0->ccw_elem) ;
}

void CArbMshTopo2D::AddFaceTwoEdgeTwoVertex(int vtx0,int vtx1,
                               int new_elem,int in_elem)
{
    // new_elem is the (zero volume) face to be added, in_elem
    // is the face within this now loop will reside

    // neither vertex should exist

    ArbVtxEntry **list0, **list1 ;

    list0 = vtx_table->Fetch(vtx0) ;
    list1 = vtx_table->Fetch(vtx1) ;
    assert(list0 == 0) ;
    assert(list1 == 0) ;

    // add the verts and angles

    InsertAngle(in_elem,vtx0,vtx1,vtx1) ;
    InsertAngle(new_elem,vtx0,vtx1,vtx1) ;

    InsertAngle(in_elem,vtx1,vtx0,vtx0) ;
    InsertAngle(new_elem,vtx1,vtx0,vtx0) ;

    // if necessary add the element

    elem_table->Store(new_elem,vtx0) ;
}

void CArbMshTopo2D::AddFaceTwoEdgeOneVertex(int vtx0,int vtx1,
                                            int Cw_vtx0,int new_elem)
{
    // determine the face in which to do the construction

    ArbVtxEntry **list,*ptr ;

    list = vtx_table->Fetch(vtx0) ;
    assert(list != 0) ;

    for (ptr = *list ; ptr != 0 ; ptr = ptr->next) {
        if (ptr->vtx_num == Cw_vtx0) break ;
    }

    // add an edge

    int next_vtx = (ptr->next == 0) ? (*list)->vtx_num :
                                      ptr->next->vtx_num ;

    InsertAngle(ptr->ccw_elem,vtx0,Cw_vtx0,vtx1) ;
    InsertAngle(ptr->ccw_elem,vtx0,vtx1,next_vtx) ;

    // now double it

    InsertAngle(new_elem,vtx0,vtx1,vtx1) ;

    InsertAngle(ptr->ccw_elem,vtx1,vtx0,vtx0) ;
    InsertAngle(new_elem,vtx1,vtx0,vtx0) ;
}

void CArbMshTopo2D::TearEdgeAddFace(int vtx0,int vtx1,int new_elem)
{
    // find the references that defind the edge

    ArbVtxEntry **list0, *ptr0 ;
    ArbVtxEntry **list1, *ptr1 ;

    list0 = vtx_table->Fetch(vtx0) ;
    assert(list0 != 0) ;
    for (ptr0 = *list0 ; ptr0 != 0 ; ptr0 = ptr0->next) {
        if (ptr0->vtx_num == vtx1) break ;
    }
    assert(ptr0 != 0) ;

    list1 = vtx_table->Fetch(vtx1) ;
    assert(list1 != 0) ;
    for (ptr1 = *list1 ; ptr1 != 0 ; ptr1 = ptr1->next) {
        if (ptr1->vtx_num == vtx0) break ;
    }
    assert(ptr1 != 0) ;

    // add angles that will link in a second edge

    InsertAngle(new_elem,vtx0,vtx1,vtx1) ;
    InsertAngle(new_elem,vtx1,vtx0,vtx0) ;
}

// The new node will be on the side of the ccw motion
// from face0 to face1

void CArbMshTopo2D::SplitVertexMakeEdge(int vtx,int new_vtx,
                                        int face0,int face1)
{
    int i ;

    // get pointers to these faces about the vertex

    ArbVtxEntry **list, **list1 ;
    ArbVtxEntry *ptr0, *ptr1 ;

    assert(face0 != face1) ;

    list = vtx_table->Fetch(vtx) ;
    assert(list != 0) ;
    for (ptr0 = *list ; ptr0 != 0 ; ptr0 = ptr0->next) {
        if (ptr0->ccw_elem == face0) break ;
    }
    assert(ptr0 != 0) ;

    for (ptr1 = *list ; ptr1 != 0 ; ptr1 = ptr1->next) {
        if (ptr1->ccw_elem == face1) break ;
    }
    assert(ptr1 != 0) ;

    // find out how many verts we need to update

    int num = 1 ;
    ptr1 = ptr0->next == 0 ? *list : ptr0->next ;
    while (ptr1->ccw_elem != face1) {
        ptr1 = ptr1->next == 0 ? *list : ptr1->next ;
        ++num ;
    }

    // fill the adjacent list

    int *adj = new int[(num+1)*2] ;
    int cur = 0 ;

    ptr1 = ptr0->next == 0 ? *list : ptr0->next ;
    adj[cur] = ptr1->vtx_num ;  adj[cur+1] = ptr1->ccw_elem ;
    cur += 2 ;
    while (ptr1->ccw_elem != face1) {
        ptr1 = ptr1->next == 0 ? *list : ptr1->next ;
        adj[cur] = ptr1->vtx_num ;  adj[cur+1] = ptr1->ccw_elem ;
        cur += 2 ;
    }
    ptr1 = ptr1->next == 0 ? *list : ptr1->next ;
    adj[cur] = ptr1->vtx_num ;  adj[cur+1] = ptr1->ccw_elem ;

    // move around the vertex and delete the angles

    int prev_start = ptr0->vtx_num ;
    int fprev_start = ptr0->ccw_elem ;
    int prev = prev_start ;
    int fprev = fprev_start ;
    int next = adj[0] ;
    int elem ;

    DeleteAngle(vtx,prev,next) ;

    for (i=0 ; i<(num-2)+1 ; i+=2) {
        prev = adj[i] ;
        elem = adj[i+1] ;
        next = adj[i+2] ;

        if (elem != NO_ELEM) DeleteAngle(vtx,prev,next) ;
    }

    prev = adj[(num-1)*2] ;
    elem = adj[(num-1)*2+1] ;
    next = adj[num*2] ;
    if (elem != NO_ELEM) DeleteAngle(vtx,prev,next) ;

    // now insert the new angles
  
    prev = prev_start ;
    fprev = fprev_start ;
    next = adj[0] ;

    InsertAngle(fprev,vtx,prev,new_vtx) ;
    InsertAngle(fprev,new_vtx,vtx,next) ;

    // now the remote vtx

    list1 = vtx_table->Fetch(next) ;
    for (ptr1 = *list1 ; ptr1 != 0 ; ptr1 = ptr1->next) {
        if ((ptr1->vtx_num == vtx) &&
            (ptr1->ccw_elem == fprev)) {
            ptr1->vtx_num = new_vtx ;
            break ;
        }
    }

    for (i=0 ; i<(num-2)+1 ; i+=2) {
        prev = adj[i] ;
        fprev = adj[i+1] ;
        next = adj[i+2] ;

        InsertAngle(fprev,new_vtx,prev,next) ;

        list1 = vtx_table->Fetch(next) ;
        for (ptr1 = *list1 ; ptr1 != 0 ; ptr1 = ptr1->next) {
            if ((ptr1->vtx_num == vtx) &&
                (ptr1->ccw_elem == fprev)) {
                ptr1->vtx_num = new_vtx ;
                break ;
            }
        }
    }

    prev = adj[(num-1)*2] ;
    fprev = adj[(num-1)*2+1] ;
    next = adj[num*2] ;

    InsertAngle(fprev,new_vtx,prev,vtx) ;
    InsertAngle(fprev,vtx,new_vtx,next) ;

    delete [] adj ;
}

// %(CArbMshTopo2D::Splice-bool-|-int-const|-int-const|-int-const|-CArbMshTopo2D-|*-int-const|-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    Splice - splice together two element boundaries 
**
**      bool Splice(
**              const int     elem,
**              const int     start,
**              const int     finish,
**              CArbMshTopo2D *other,
**              const int     other_elem,
**              const int     other_start,
**              const int     other_finish)
**
**        elem         - (in)  this elemnent id 
**        start        - (in)  this start node 
**        finish       - (in)  this finish node 
**        other        - (in)  other mesh topology 
**        other_elem   - (in)  other element id 
**        other_start  - (in)  start node on other elem 
**        other_finish - (in)  finish node on other elem 
**
**      Description: This function splices the boundary for an element 
**          stored in another mesh to the boundary of an element stored 
**          in this mesh. This is done as follows: 
**
**               +                                +
**               |   other                        |
**               |   start                        |
**        start  +     +----+               start +----+
**               |     |    | other   ==>              |
**               |     |    | element                  |
**        finish +     +----+              other  +----+
**               |   other                 finish |
**               |   finish                       |
**               +                                +
**        element boundary
**            
**
**      Return Value: returns true if the splice was succesful 
**
**
** -- */

bool CArbMshTopo2D::Splice(const int elem,
                           const int start,
                           const int finish,
                           CArbMshTopo2D *other,
                           const int other_elem,
                           const int other_start,
                           const int other_finish)
{
    // first on the original element, find the start and
    // finish verticies, and nodes before and after these

    CArbTopoVtxOnElemCyclicIterator iter0(this,elem,start) ;
    if (!iter0.IsValid()) return(false) ;

    --iter0 ;
    int prev0 = *iter0 ;

    iter0.First() ;
    while (*iter0 != finish) ++iter0 ;

    ++iter0 ;
    int next0 = *iter0 ;
    ++iter0 ;
    int next_next0 = *iter0 ;

    // now find similar information on the other element

    CArbTopoVtxOnElemIterator iter1(other,other_elem,other_finish) ;
    if (!iter1.More()) return(false) ;

    ++iter1 ;
    while (iter1.More() &&
           (*iter1 != other_start)) ++iter1 ;
    if (!iter1.More()) return(false) ;

    // on the original element, go from the start to the
    // finish and delete all the angles

    ArbVtxEntry *ptr = FirstVtx(start) ;
    while (ptr->ccw_elem != elem) ptr = NextVtx(ptr) ;
    int nxt = ptr->vtx_num ;
    int cur = start ;
    int prv ;
    DeleteAngle(cur,nxt,prev0) ;

    if (start != finish) {
        while (nxt != finish) {
            prv = cur ;
            cur = nxt ;
            ptr = FirstVtx(cur) ;
            while (ptr->ccw_elem != elem) ptr = NextVtx(ptr) ;
            nxt = ptr->vtx_num ;
            DeleteAngle(cur,nxt,prv) ;
        }
        DeleteAngle(finish,next0,cur) ;
    }

    DeleteAngle(next0,next_next0,finish) ;

    // now splice the new angles from the other element

    InsertAngle(elem,next0,next_next0,other_finish) ;
    cur = other_finish ;
    ptr = other->FirstVtx(cur) ;
    while (ptr->ccw_elem != other_elem) ptr = other->NextVtx(ptr) ;
    cur = other_finish ;
    prv = next0 ;
    nxt = ptr->vtx_num ;
  
    while (nxt != other_start) {
        InsertAngle(elem,cur,prv,nxt) ;
        prv = cur ;
        cur = nxt ;
        ptr = other->FirstVtx(cur) ;
        while (ptr->ccw_elem != other_elem) ptr = other->NextVtx(ptr) ;
        nxt = ptr->vtx_num ;
    }

    InsertAngle(elem,cur,prv,start) ;
    InsertAngle(elem,start,cur,prev0) ;
    return(true) ;
}




// %(CArbMshTopo2D::DeleteElementRef-void-|-int-const|)
/* ++ ----------------------------------------------------------
**
**    DeleteElementRef - delete an element reference 
**
**      void DeleteElementRef(const int elem_id)
**
**        elem_id - (in)  element id 
**
**      Description: Deletes an entry in the element table 
**
**
** -- */

void CArbMshTopo2D::DeleteElementRef(const int elem)
{
    elem_table->Remove(elem) ;
}




// %(CArbMshTopo2D::HasVtx-bool-|^const-int-const|)
/* ++ ----------------------------------------------------------
**
**    HasVtx - vertex query 
**
**      bool HasVtx(const int vtx) const
**
**        vtx - (in)  vertex id 
**
**      Description: Checks to see if a vertex id is in the mesh. 
**
**      Return Value: returns true if the vetex (node) is in the mesh 
**
**
** -- */

bool CArbMshTopo2D::HasVtx(const int vtx) const
{
    ArbVtxEntry **list ;

    if (((list = vtx_table->Fetch(vtx)) != 0) && (*list !=0)) return(true) ;
    return(false) ;
}




// %(CArbMshTopo2D::FirstVtx-ArbVtxEntry-|*^const-int-const|) 
/* ++ ----------------------------------------------------------
**
**    FirstVtx - get the first adjacent vertex 
**
**      ArbVtxEntry *FirstVtx(const int vtx) const
**
**        vtx - (in)  input vertex 
**
**      Description: This function returns a pointer to the ArbVtxEntry 
**          data structure for the first vertex in the linked list of 
**          vertices adjacent to the given vertex. 
**
**      Return Value: a pointer to the first adjacent vertex 
**
**
** -- */

CArbMshTopo2D::ArbVtxEntry *CArbMshTopo2D::FirstVtx(const int vtx) const
{
    ArbVtxEntry **list ;

    // first check to see if this node is in the hash table.

    if ((list = vtx_table->Fetch(vtx)) == 0) {
        return(0) ;
    } else {
        return(*list) ;
    }
}




// %(CArbMshTopo2D::OppositeEdge-ArbEdge-|^const-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    OppositeEdge - find the edge opposite a node on an element 
**
**      ArbEdge OppositeEdge(
**              const int elem_id,
**              const int node_id) const
**
**        elem_id - (in)  element id 
**        node_id - (in)  node id 
**
**      Description: Given an element and a node on the element, this 
**          function returns the edge on the element opposite the node. 
**          This function assumes that that the element is a triangle. 
**
**      Return Value: returns the opposite edge 
**
**
** -- */

CArbMshTopo2D::ArbEdge CArbMshTopo2D::OppositeEdge(
                     const int elem_id,
                     const int node_id) const
{
    ArbEdge opp_edge ;

    // find the vertex adjacent to the given vertex which
    // is on the clockwise side of the given element

    int adj_node = GetAdjVtxAlongElem(node_id,elem_id) ;
    assert(adj_node != NO_NODE) ;

    // now we step around the element two nodes to fill
    // in the edge data

    opp_edge.nd0 = adj_node ;
    opp_edge.nd1 = GetAdjVtxAlongElem(adj_node,elem_id) ;

    // now fill in the elem data

    opp_edge.elem0 = elem_id ;
    opp_edge.elem1 = GetCWElem(opp_edge.nd0,opp_edge.nd1) ;
    return(opp_edge) ;
}




// %(CArbMshTopo2D::OppositeNode-int-|^const-int-const|-ArbEdge-const|*)
/* ++ ----------------------------------------------------------
**
**    OppositeNode - find the node opposite an edge on an element 
**
**      int OppositeNode(
**              const int     elem_id,
**              const ArbEdge *edge) const
**
**        elem_id - (in)  element id 
**        edge    - (in)  edge description 
**
**      Description: Given an element and an edge on the element, this 
**          function returns the node on the element opposite the edge. 
**          This function assumes that that the element is a triangle. 
**
**      Return Value: the opposite node 
**
**
** -- */

int CArbMshTopo2D::OppositeNode(
                     const int elem_id,
                     const ArbEdge *edge) const
{
    // check to make sure that this edge is adjacent to
    // the element

    int adj_node = GetAdjVtxAlongElem(edge->nd0,elem_id) ;
    if (adj_node == NO_NODE)
     return -1;
    //assert(adj_node == edge->nd1) ;

    // return the next node ccw about the element

    return(GetAdjVtxAlongElem(adj_node,elem_id)) ;
}




// %(CArbMshTopo2D::NextCWEdge-ArbEdge-|^const-int-const|-ArbEdge-const|*)
/* ++ ----------------------------------------------------------
**
**    NextCWEdge - get the next edge cw on an element 
**
**      ArbEdge NextCWEdge(
**              const int     elem_id,
**              const ArbEdge *edge) const
**
**        elem_id - (in)  element id 
**        edge    - (in)  current edge 
**
**      Description: Give an element and an edge, this function returns 
**          the next edge moving in a clockwise direction. 
**
**      Return Value: the next edge in a cw direction 
**
**
** -- */

CArbMshTopo2D::ArbEdge CArbMshTopo2D::NextCWEdge(
                   const int elem_id,
                   const ArbEdge *edge) const
{
    ArbEdge cw_edge ;

    // check to make sure that this edge is adjacent to
    // the element

    int adj_node = GetAdjVtxAlongElem(edge->nd0,elem_id) ;
    if (adj_node == NO_NODE)
    {
     cw_edge.nd0 = -1;
     cw_edge.nd1 = -1;
     cw_edge.elem0 = -1;
     cw_edge.elem1 = -1;
     return cw_edge;
    }
    //assert(adj_node == edge->nd1) ;

    // move around the element until we get back to the
    // starting point

    int cur_node = edge->nd1 ;
    adj_node = GetAdjVtxAlongElem(edge->nd1,elem_id) ;
    if (adj_node == NO_NODE)
    {
     cw_edge.nd0 = -1;
     cw_edge.nd1 = -1;
     cw_edge.elem0 = -1;
     cw_edge.elem1 = -1;
     return cw_edge;
    }
    while (adj_node != edge->nd0) 
    {
        cur_node = adj_node ;
        adj_node = GetAdjVtxAlongElem(adj_node,elem_id) ;
        if (adj_node == NO_NODE)
        {
         cw_edge.nd0 = -1;
         cw_edge.nd1 = -1;
         cw_edge.elem0 = -1;
         cw_edge.elem1 = -1;
         return cw_edge;
        }
    }

    // build and return the edge information

    cw_edge.nd0 = cur_node ;
    cw_edge.nd1 = adj_node ;
    cw_edge.elem0 = elem_id ;
    cw_edge.elem1 = GetCWElem(cw_edge.nd0,cw_edge.nd1) ;

    return(cw_edge) ;
}




// %(CArbMshTopo2D::NextCCWEdge-ArbEdge-|^const-int-const|-ArbEdge-const|*)
/* ++ ----------------------------------------------------------
**
**    NextCCWEdge - get the next edge cw on an element 
**
**      ArbEdge NextCCWEdge(
**              const int     elem_id,
**              const ArbEdge *edge) const
**
**        elem_id - (in)  element id 
**        edge    - (in)  current edge 
**
**      Description: Give an element and an edge, this function returns 
**          the next edge moving in a counter clockwise direction. 
**
**      Return Value: the next edge in a ccw direction 
**
**
** -- */

CArbMshTopo2D::ArbEdge CArbMshTopo2D::NextCCWEdge(
                    const int elem_id,
                    const ArbEdge *edge) const
{
    ArbEdge ccw_edge ;

    // check to make sure that this edge is adjacent to
    // the element

    int adj_node = GetAdjVtxAlongElem(edge->nd0,elem_id) ;
    if (adj_node == NO_NODE)
    {
     ccw_edge.nd0 = -1;
     ccw_edge.nd1 = -1;
     ccw_edge.elem0 = -1;
     ccw_edge.elem1 = -1;
     return(ccw_edge) ;
    }

    //assert(adj_node == edge->nd1) ;

    // build and return the edge information

    ccw_edge.nd0 = edge->nd1 ;
    ccw_edge.nd1 = GetAdjVtxAlongElem(edge->nd1,elem_id) ;
    ccw_edge.elem0 = elem_id ;
    ccw_edge.elem1 = GetCWElem(ccw_edge.nd0,ccw_edge.nd1) ;

    return(ccw_edge) ;
}




// %(CArbMshTopo2D::GetCCWBdryNode-int-|^const-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetCCWBdryNode - get the next node ccw on the boundary 
**
**      int GetCCWBdryNode(
**              const int node_0,
**              const int node_1) const
**
**        node_0 - (in)  first node 
**        node_1 - (in)  second node 
**
**      Description: This function find the node that is counter 
**          clockwise from the given edge, staying on the boundary of 
**          the mesh. 
**
**      Return Value: the ccw node id 
**
**
** -- */

int CArbMshTopo2D::GetCCWBdryNode(
                   const int node_0,
                   const int node_1) const
{
    // look at the edges adjacent to node_1 and try to
    // find node_0.

    ArbVtxEntry **list, *ptr, *next = 0 ;

    list = vtx_table->Fetch(node_1) ;
    //assert(list != 0) ;
    if (list == 0)
     return -1;

    for (ptr = *list ; ptr != 0 ; ptr = ptr->next) 
    {
        if (ptr->vtx_num == node_0) break ;
    }
    //assert(ptr != 0) ;
    if (ptr == 0)
     return -1;

    // search from this all the way around the vertex
    // keeping track of when we leave NO_ELEM regions

    bool in_void = true ;
    ptr = ptr->next ;
    if (ptr == 0) ptr = *list ;
    while (ptr->vtx_num != node_0) 
    {
        if (ptr->ccw_elem == NO_ELEM)
        {
            in_void = true ;
        } else if (in_void) {
            next = ptr ;
            in_void = false ;
        }
        ptr = ptr->next ;
        if (ptr == 0) 
         ptr = *list ;
    }
    //assert(next != 0) ;
    if (next == 0)
     return -1;
    else
    return(next->vtx_num) ;
}




// %(CArbMshTopo2D::GetCWBdryNode-int-|^const-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetCWBdryNode - get the next node cw on the boundary 
**
**      int GetCWBdryNode(
**              const int node_0,
**              const int node_1) const
**
**        node_0 - (in)  first node 
**        node_1 - (in)  second node 
**
**      Description: This function find the node that is clockwise from 
**          the given edge, staying on the boundary of the mesh. 
**
**      Return Value: the cw node id 
**
**
** -- */

int CArbMshTopo2D::GetCWBdryNode(
                   const int node_0,
                   const int node_1) const
{
    // look at the edges adjacent to node_0 and try to
    // find node_1.

    ArbVtxEntry **list, *ptr, *prev = 0 ;

    list = vtx_table->Fetch(node_0) ;
    //assert(list != 0) ;
    if (list == 0)
     return -1;

    for (ptr = *list ; ptr != 0 ; ptr = ptr->next) {
        if (ptr->vtx_num == node_1) break ;
    }
    
    //assert(ptr != 0) ;
    if (ptr == 0)
     return -1;

    // now move counter clockwise from here until we find
    // the first edge that points to NO_ELEM

    ptr = ptr->next ;
    if (ptr == 0) ptr = *list ;
    while (ptr->vtx_num != node_1) 
    {
        if (ptr->ccw_elem == NO_ELEM) {
            prev = ptr ;
            break ;
        }
        ptr = ptr->next ;
        if (ptr == 0) ptr = *list ;
    }

    //assert(prev != 0) ;
    if (prev == 0)
     return -1;
    else
    return(prev->vtx_num) ;
}




// %(CArbMshTopo2D::GetCCWNode-int-|^const-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetCCWNode - get the next node ccw on this element 
**
**      int GetCCWNode(
**              const int node_0,
**              const int node_1) const
**
**        node_0 - (in)  first node 
**        node_1 - (in)  second node 
**
**      Description: This function find the node that is counter 
**          clockwise from the given edge, staying on the same element. 
**
**      Return Value: the ccw node id 
**
**
** -- */

int CArbMshTopo2D::GetCCWNode(
                   const int node_0,
                   const int node_1) const
{
    // look at the edges adjacent to node_1 and try to
    // find node_0.

    ArbVtxEntry **list, *ptr ;

    list = vtx_table->Fetch(node_1) ;
    //assert(list != 0) ;
    if (list == 0)
     return -1;

    for (ptr = *list ; ptr != 0 ; ptr = ptr->next) 
    {
        if (ptr->vtx_num == node_0) break ;
    }

    //assert(ptr != 0) ;
    if (ptr == 0)
     return -1;

    if (ptr->next == 0) 
    {
        return((*list)->vtx_num) ;
    }
    return(ptr->next->vtx_num) ;
}




// %(CArbMshTopo2D::GetCWNode-int-|^const-int-const|-int-const|) 
/* ++ ----------------------------------------------------------
**
**    GetCWNode - get the next node cw on this element 
**
**      int GetCWNode(
**              const int node_0,
**              const int node_1) const
**
**        node_0 - (in)  first node 
**        node_1 - (in)  second node 
**
**      Description: This function find the node that is clockwise from 
**          the given edge, staying on the same element. 
**
**      Return Value: the cw node id 
**
**
** -- */

int CArbMshTopo2D::GetCWNode(
                   const int node_0,
                   const int node_1) const
{
    ArbVtxEntry **list, *ptr, *first = 0 ;
    int cw = 0 ;

    list = vtx_table->Fetch(node_0) ;
    //assert(list != 0);
    if (list == 0)
     return -1;

    first = *list ;
    for (ptr = *list ; ptr != 0 ; ptr = ptr->next) 
    {
        if ((ptr->vtx_num == node_1) &&
            (ptr != first)) return(cw) ;
        cw = ptr->vtx_num ;
    }
    return(cw) ;
}


int CArbMshTopo2D::GetDblEdgeElem(int node_0,int node_1) const
{
    // look at the edges adjacent to node_1 and try to
    // find node_0.

    ArbVtxEntry **list, *ptr ;
    int elem = -1 ;

    list = vtx_table->Fetch(node_0) ;
    //assert(list != 0) ;
    if (list == 0)
     return -1;
    for (ptr = *list ; ptr != 0 ; ptr = ptr->next) 
    {
        if (ptr->vtx_num == node_1) 
        {
            if (elem < 0) {
                elem = ptr->ccw_elem ;
            } else {
                break ;
            }
        }
    }
    return(elem) ;
}


bool CArbMshTopo2D::HasDblEdge(int node_0,int node_1) const
{
    return(GetDblEdgeElem(node_0,node_1) != NO_ELEM) ;
}


// %(CArbMshTopo2D::ElemHasNode-bool-|^const-int-const|-int-const|) 
/* ++ ----------------------------------------------------------
**
**    ElemHasNode - check for a node on an element 
**
**      bool ElemHasNode(
**              const int elem_id,
**              const int node_id) const
**
**        elem_id - (in)  element id 
**        node_id - (in)  node id 
**
**      Description: This function checks to see if a node is on an 
**          element. 
**
**      Return Value: returns true if the node is on the element 
**
**
** -- */

bool CArbMshTopo2D::ElemHasNode(const int elem_id,
                                const int node_id) const
{
    // use the element table to get the id of one of the nodes
    // on the element and if this element exists in the table,
    // check to see if the node is the one we are looking for.

    int *node_ptr = elem_table->Fetch(elem_id) ;
    if (node_ptr == 0)
     return false;
    //assert(node_ptr != 0) ;

    int start_node = *node_ptr ;
    
    if (start_node == node_id) return(true) ;

    // Now we traverse around the element looking for the node
    // and checking to see if got back to the starting node.

    int adj_node = start_node ;

    // rgd08
    int num = 0;
    do {
        adj_node = GetAdjVtxAlongElem(adj_node,elem_id) ;
        if (adj_node == NO_NODE)
         return false;
        if (adj_node == node_id) return(true) ;
        num++;
        if (num > 2000)
       break;
    } while (adj_node != start_node) ;

    return(false) ;
}




// %(CArbMshTopo2D::InsertAngle-void-|-int-const|-int-const|-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    InsertAngle - insert an angle 
**
**      void InsertAngle(
**              const int elem,
**              const int vtx,
**              const int prev,
**              const int next)
**
**        elem - (in)  element id 
**        vtx  - (in)  center vertex 
**        prev - (in)  previous vertex 
**        next - (in)  next vertex 
**
**      Description: Inserts three verticies and the included element 
**          into a mesh if they do not exist in the mesh already. The 
**          element will exist in the region moving ccw from prev to 
**          next. 
**
**
** -- */

int CArbMshTopo2D::InsertAngle(const int elem, const int vtx, const int prev,
                                const int next, const int parent_elem) 
{
    ArbVtxEntry **list, *entry ;
    ArbVtxEntry *ptr, *p_ptr, *n_ptr, *p_start = NULL, *n_end ;
    ArbVtxEntry *pp_start = NULL, *pn_ptr ;


    // first check to see if this node is already in the
    // hash table.  If not, add it

    if ((list = vtx_table->Fetch(vtx)) == 0)
    {
        vtx_table->Store(vtx,0) ;
        list = vtx_table->Fetch(vtx) ;
        if (vtx >= MaxVtx) MaxVtx = vtx+1 ;
    }


    // if the list is empty, add both nodes to the list

    if (*list == 0)
    {
        entry = NewVtxEntry() ;
        if (entry >= MaxAddress) MaxAddress = entry+1 ;
        entry->vtx_num = next ;
        entry->next = *list ;
        entry->ccw_elem = parent_elem < 0 ? NO_ELEM : parent_elem ;
        *list = entry ;
        if (prev != next)
        {
            entry = NewVtxEntry() ;
            if (entry >= MaxAddress) MaxAddress = entry+1 ;
            entry->vtx_num = prev ;
            entry->next = *list ;
            entry->ccw_elem = elem ;
            *list = entry ;
        }
        else 
        {
            entry->ccw_elem = elem ;
        }
    }
    else
    {

        // search through the list and find the
        // position of the prev and next verticies
        // also find the position of p_start, which is the
        // beginning of the contiguous block of verts that
        // preceeds the previous vert

        p_ptr = n_ptr = 0 ;
        p_start = *list ;
        for (ptr = *list ; ptr != 0 ; ptr = ptr->next)
        {
            if (ptr->vtx_num == prev) p_ptr = ptr ; 
            if (ptr->vtx_num == next) n_ptr = ptr ;
            if (p_ptr == 0)
            {
                if (p_start == 0) p_start = ptr ;
                if (ptr->ccw_elem == NO_ELEM) p_start = 0 ;
            }
        }

        // now do a case analysis depending on what's already
        // in the list.
        //
        // case 1, neither in the list, add after first
        // node that does not point to an element

        if ((p_ptr == 0) && (n_ptr == 0))
        {
            //bool found = false ;
            for (ptr = *list ; ptr != 0 ; ptr = ptr->next)
            {
                if (ptr->ccw_elem == NO_ELEM)
                {
                    entry = NewVtxEntry() ;
                    if (entry >= MaxAddress) MaxAddress = entry+1 ;
                    entry->vtx_num = next ;
                    entry->next = ptr->next ;
                    entry->ccw_elem = parent_elem < 0 ? 
                                      NO_ELEM : parent_elem ;
                    ptr->next = entry ;
                    entry = NewVtxEntry() ;
                    if (entry >= MaxAddress) MaxAddress = entry+1 ;
                    entry->vtx_num = prev ;
                    entry->next = ptr->next ;
                    entry->ccw_elem = elem ;
                    ptr->next = entry ;
                    //found = true ;
                    break ;
                } 
            }
            
            //assert(found) ;

        // case 2, prev in the list, add next after it

        } 
        else if ((p_ptr != 0) && (n_ptr == 0))
        {
            entry = NewVtxEntry() ;
            if (entry >= MaxAddress) MaxAddress = entry+1 ;
            entry->vtx_num = next ;
            entry->next = p_ptr->next ;
            entry->ccw_elem = parent_elem < 0 ? NO_ELEM : parent_elem ;
            p_ptr->next = entry ;
            p_ptr->ccw_elem = elem ;

        // case 3, next in the list, add prev before it

        } 
        else if ((p_ptr == 0) && (n_ptr != 0))
        {
            ArbVtxEntry **prev_ptr = list ; 
            for (ptr = *list ; ptr != 0 ; ptr = ptr->next)
            {
                if (ptr->vtx_num == next)
                {
                    entry = NewVtxEntry() ;
                    if (entry >= MaxAddress) MaxAddress = entry+1 ;
                    entry->vtx_num = prev ;
                    entry->next = ptr ;
                    entry->ccw_elem = elem ;
                    *prev_ptr = entry ;
                    break ;
                }
                prev_ptr = &(ptr->next) ;
            }

        // case4, both in list,  there are two sub-cases,
        //    1) the verticies are already in sequential
        //       order so we set a flag and return.
        //    2) the prev and next vert is the same.
        //    3) the verticies are not in order so we need
        //       to rearrange the list

        }
        else 
        {
            if (p_ptr == n_ptr)
            {
               if (elem != p_ptr->ccw_elem)
               {
                   entry = NewVtxEntry() ;
                   if (entry >= MaxAddress) MaxAddress = entry+1 ;
                   entry->vtx_num = next ;
                   entry->next = p_ptr->next ;
                   entry->ccw_elem = p_ptr->ccw_elem ;
                   p_ptr->ccw_elem = elem ;
                   p_ptr->next = entry ;
               }
            } 
            else if ((p_ptr->next == n_ptr) ||((p_ptr->next == 0) && (*list == n_ptr)))
            {
                p_ptr->ccw_elem = elem ;
            }
            else
            {

                // find the end of the contiguous verts starting
                // at the next ptr

                ptr = n_ptr ;
                while (ptr->ccw_elem != NO_ELEM) 
                {
                 if (ptr->next == 0)
                  return 0;
                 ptr = ptr->next ;
                }
                n_end = ptr ;

                // find the entry that points to the n_ptr

                if (n_ptr == *list)
                {
                    pn_ptr = 0 ;
                    ptr = *list ;
                    while (ptr->next != p_start) 
                    {
                     if (ptr->next == 0)
                      return 0;
                     ptr = ptr->next ;
                    }
                    pp_start = ptr ;
                } 
                else 
                {
                    ptr = *list ;
                    while (ptr->next != n_ptr) 
                    {
                     if (ptr->next == 0)
                      return 0;
                     ptr = ptr->next ;
                    }
                    pn_ptr = ptr ;
                }

                // first the case where the next block is first
                // in the list

                if (pn_ptr == 0)
                {

                    *list = p_start ;
                    pp_start->next = p_ptr->next ;
                    p_ptr->next = n_ptr ;

                    // else the more general case

                }
                else
                {

                    pn_ptr->next = n_end->next ;
                    n_end->next = p_ptr->next ;
                    p_ptr->next = n_ptr ;

                }
                p_ptr->ccw_elem = elem ;
            }
        }

    // now we want to reorder so that the entry that
    // points to NO_ELEM is at the end of the list

        n_ptr = 0 ;
        for (ptr = *list ; ptr != 0 ; ptr = ptr->next)
        {
            if (ptr->ccw_elem == NO_ELEM)
            {
                n_ptr = ptr->next ;
                ptr->next = 0 ;
                break ;
            }
        }
        ptr = n_ptr ;
        while (ptr != 0)
        {            // switch to while loop
            if (ptr->next == 0) 
            {     // because of a compiler bug
                ptr->next = *list ;
                *list = n_ptr ;
                break ;
            }
            ptr = ptr->next ;
        }
    }
    return 1;
}




// %(CArbMshTopo2D::DeleteAngle-int-|-int-const|-int-const|-int-const|) 
/* ++ ----------------------------------------------------------
**
**    DeleteAngle - delete an angle 
**
**      int DeleteAngle(
**              const int vtx,
**              const int prev,
**              const int next)
**
**        vtx  - (in)  center vertex 
**        prev - (in)  previous vertex 
**        next - (in)  next vertex 
**
**      Description: Deletes an angle from the mesh. 
**
**      Return Value: returns the formerly included element 
**
**
** -- */

int CArbMshTopo2D::DeleteAngle(const int vtx,
                                        const int prev,
                                        const int next)
{
    ArbVtxEntry **list ;
    ArbVtxEntry *ptr, *p_ptr, *n_ptr ;
    int elem = 0 ;

    // first check to see if this node is already in the
    // hash table.  If not, return

    list = vtx_table->Fetch(vtx) ;
    //assert(list != 0) ;
    if (list == 0)
     return -1;

    // if the list is empty just return

    //assert(*list != 0) ;
    if (*list == 0)
     return -1;

    // search through the list.  Find the prev vertex and
    // make sure that next is the next vertex in the list

    p_ptr = n_ptr = 0 ;
    if (prev != next) {
        for (ptr = *list ; ptr != 0 ; ptr = ptr->next) {
            if (ptr->vtx_num == prev) p_ptr = ptr ; 
            if (ptr->vtx_num == next) n_ptr = ptr ; 
        }
    } else {
        for (ptr = *list ; ptr != 0 ; ptr = ptr->next) {
            if (ptr->vtx_num == prev) {
                p_ptr = ptr ;
                break ;
            }
        }
        if (p_ptr == 0)
         return -1;
        if ((p_ptr->next != 0) && (p_ptr->next->vtx_num == next)) {
            n_ptr = p_ptr->next ;
        } else {
            n_ptr = p_ptr ;
        }
    }

    //assert((p_ptr != 0) && (n_ptr != 0)) ;
    //if ((p_ptr == 0) && (n_ptr == 0)) ;
    //return -1;

    // deal with the special case of a "dangling" edge

    if (p_ptr == n_ptr) 
    {
        if (p_ptr == 0)
         return -1;
        elem = p_ptr->ccw_elem ;
        *list = n_ptr->next ;
        int *elem_vtx = elem_table->Fetch(elem) ;
        if ((elem_vtx != 0) && (*elem_vtx == vtx)) {

            // check to make sure that the adjacent vtx
            // still exists.  If so, make this the new stored
            // node for the element.  Otherwise do a search

            elem_table->Remove(elem) ;
            ArbVtxEntry **adj = vtx_table->Fetch(p_ptr->vtx_num) ;
            if ((adj == 0) || (*adj == 0)) {
                CArbHashTableIterator<int,ArbVtxEntry *>
                                               iter(vtx_table) ;

                bool found = false ;
                for (iter.First() ; iter.More() ; ++iter) {
                    for (ptr = *(iter.Entry()) ; ptr != 0 ; ptr = ptr->next) {
                        if (ptr->ccw_elem == elem) {
                            elem_table->Store(elem,ptr->vtx_num) ;
                            found = true ;
                            break ;
                        }
                    }
                    if (found) break ;
                }

            } else {
                elem_table->Store(elem,p_ptr->vtx_num) ;
            }
        }
        DeleteVtxEntry(n_ptr) ;
        return(elem) ;
    }

    // check to see if this is the node stored for the
    // element.  If so we need to reset this

    //rgd08
    if (p_ptr==0)
     return -1;

    elem = p_ptr->ccw_elem ;
    int *elem_vtx = elem_table->Fetch(elem) ;
    if ((elem_vtx != 0) && (*elem_vtx == vtx)) {
        elem_table->Remove(elem) ;
        ArbVtxEntry **alist ;
        alist = vtx_table->Fetch(p_ptr->vtx_num) ;
        if ((alist != 0) && (*alist != 0)) {
            elem_table->Store(elem,p_ptr->vtx_num) ;
        } else {
            CArbHashTableIterator<int,ArbVtxEntry *>
                                               iter(vtx_table) ;

            bool found = false ;
            for (iter.First() ; iter.More() ; ++iter) {
                for (ptr = *(iter.Entry()) ; ptr != 0 ; ptr = ptr->next) {
                    if ((ptr->vtx_num != vtx) && (ptr->ccw_elem == elem)) {
                        alist = vtx_table->Fetch(ptr->vtx_num) ;
                        if ((alist != 0) && (*alist != 0)) {
                            elem_table->Store(elem,ptr->vtx_num) ;
                            found = true ;
                            break ;
                        }
                    }
                }
                if (found) break ;
            }
        }
    }

    // now the normal case

    if ((p_ptr->next == n_ptr) ||
        ((p_ptr->next == 0) && (*list == n_ptr))) {
        elem = p_ptr->ccw_elem ;
        p_ptr->ccw_elem = NO_ELEM ;

        // if we have two adjacent NO_ELEM regions pointed
        // to by the previous and next pointers, then remove
        // the "next" edge between them

        if (n_ptr == 0)
         return -1;

        if (n_ptr->ccw_elem == NO_ELEM) {
            if (p_ptr->next != 0) {
                p_ptr->next = n_ptr->next ;
            } else {
                *list = n_ptr->next ;
            }
            DeleteVtxEntry(n_ptr) ;
        }

        // check to see if the vertex before the prev
        // vertex points to a NO_ELEM region.  In this case
        // we remove the "prev" vertex

        if (*list != p_ptr) {
            ArbVtxEntry  *prev_ptr = 0 ;
            for (ptr = *list ; ptr != p_ptr ; ptr = ptr->next) {
                prev_ptr = ptr ;
            }
            if (prev_ptr->ccw_elem == NO_ELEM) {
                prev_ptr->next = p_ptr->next ;
                DeleteVtxEntry(p_ptr) ;
            }
        } else {
            ArbVtxEntry  *prev_ptr = 0 ;
            for (ptr = *list ; ptr != 0 ; ptr = ptr->next) {
                prev_ptr = ptr ;
            }
            if (prev_ptr->ccw_elem == NO_ELEM) {
                *list = p_ptr->next ;
                DeleteVtxEntry(p_ptr) ;
            }
        }
    }

    // now we want to reorder so that the entry that
    // points to NO_ELEM is at the end of the list

    n_ptr = 0 ;
    for (ptr = *list ; ptr != 0 ; ptr = ptr->next) {
        if (ptr->ccw_elem == NO_ELEM) {
            n_ptr = ptr->next ;
            ptr->next = 0 ;
            break ;
        }
    }
    ptr = n_ptr ;
    while(ptr != 0) {              // switch to a while loop
        if (ptr->next == 0) {      // because of a compiler bug
            ptr->next = *list ;
            *list = n_ptr ;
            break ;
        }
        ptr = ptr->next ;
    }

    return(elem) ;
}




// %(CArbMshTopo2D::GetAdjVtxAlongElem-int-|^const-int-|-int-|)
/* ++ ----------------------------------------------------------
**
**    GetAdjVtxAlongElem - get the next vertex on an element 
**
**      int GetAdjVtxAlongElem(
**              int node,
**              int elem_id) const
**
**        node    - (in)  node id 
**        elem_id - (in)  element id 
**
**      Description: Give an element and a vertex on the element, this 
**          function returns the next vertex on the element moving in a 
**          ccw direction. 
**
**      Return Value: the id of the next ccw vertex on an element 
**
**
** -- */

int CArbMshTopo2D::GetAdjVtxAlongElem(const int node,
                                               const int elem_id) const
{
    ArbVtxEntry **list, *ptr ;

    // first check to see if this node is in the hash table.

    list = vtx_table->Fetch(node) ;
    
    //assert(list != 0) ;    
    if (list == 0)
     return NO_NODE;

    for (ptr = *list ; ptr != 0 ; ptr = ptr->next) {
        if (ptr->ccw_elem == elem_id) return(ptr->vtx_num) ;
    }
    //assert(0) ;
    return(NO_NODE) ;
}




// %(CArbMshTopo2D::GetCWElem-int-|^const-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetCWElem - get the cw element 
**
**      int GetCWElem(
**              const int node_0,
**              const int node_1) const
**
**        node_0 - (in)  first node 
**        node_1 - (in)  second node 
**
**      Description: This function finds the element that is clockwise 
**          to the edge from node 0 to node 1. 
**
**      Return Value: the cw element id 
**
**
** -- */

int CArbMshTopo2D::GetCWElem(const int node_0,
                                      const int node_1) const
{
    ArbVtxEntry **list, *ptr, *first ;
    int elem = NO_ELEM ;

    // first check to see if this node is in the hash table.

    list = vtx_table->Fetch(node_0) ;
    //assert(list != 0) ;
    if (list == 0)
     return -1;

    // loop around the vertex, saving the previous ccw element

    first = *list ;
    for (ptr = *list ; ptr != 0 ; ptr = ptr->next) 
    {
        if ((ptr->vtx_num == node_1) &&
            (ptr != first)) return(elem) ;
        elem = ptr->ccw_elem ;
    }
    return(elem) ;
}

int CArbMshTopo2D::GetCWElemPtr(ArbVtxEntry *ptr,
                                ArbVtxEntry **list)
{
    ArbVtxEntry *lptr ;
    int elem ;

    elem = ptr->ccw_elem ;
    lptr = ptr->next ;
    if (lptr == 0) lptr = *list ;

    while (lptr != ptr) {
        elem = lptr->ccw_elem ;
        lptr = lptr->next ;
        if (lptr == 0) lptr = *list ;
    }
    return(elem) ;
}



// %(CArbMshTopo2D::GetCCWElem-int-|^const-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    GetCCWElem - get the ccw element 
**
**      int GetCCWElem(
**              const int node_0,
**              const int node_1) const
**
**        node_0 - (in)  first node 
**        node_1 - (in)  second node 
**
**      Description: This function finds the element that is counter 
**          clockwise to the edge from node 0 to node 1. 
**
**      Return Value: the ccw element id 
**
**
** -- */

int CArbMshTopo2D::GetCCWElem(const int node_0,
                                       const int node_1) const
{
    ArbVtxEntry **list, *ptr ;
    int elem = NO_ELEM ;

    // first check to see if this node is in the hash table.

    list = vtx_table->Fetch(node_0) ;
    //assert(list != 0) ;
    if (list == 0)
     return -1;

    // loop around the vertex, saving the previous ccw element

    for (ptr = *list ; ptr != 0 ; ptr = ptr->next) 
    {
        if (ptr->vtx_num == node_1) {
            elem = ptr->ccw_elem ;
            break ;
        }
    }
    return(elem) ;
}


void CArbMshTopo2D::GetElemsAboutEdge(const int node_0,const int node_1,
                                      int *elem_0,int *elem_1) const
{
    *elem_0 = GetCCWElem(node_0,node_1) ;
    *elem_1 = GetCCWElem(node_1,node_0) ;
}


// %(CArbMshTopo2D::NumAdjElems-int-|^const-int-const|)
/* ++ ----------------------------------------------------------
**
**    NumAdjElems - find the number of adjacent element 
**
**      int NumAdjElems(const int id) const
**
**        id - (in)  node id 
**
**      Description: Find the number of elements adjacent to a given 
**          node. 
**
**      Return Value: number of adjacent elements 
**
**
** -- */

int CArbMshTopo2D::NumAdjElems(const int id) const
{
    ArbVtxEntry **list, *ptr ;
    int num = 0 ;

    // first check to see if this node is in the hash table.

    list = vtx_table->Fetch(id) ;
    //assert(list != 0) ;
    if (list == 0)
     return -1;

    // loop around the vertex counting the elements

    for (ptr = *list ; ptr != 0 ; ptr = ptr->next)
    {
        if (ptr->ccw_elem != NO_ELEM) ++num ;
    }
    return(num) ;
}




// %(CArbMshTopo2D::NumConsecutiveAdjElems-int-|^const-int-const|)
/* ++ ----------------------------------------------------------
**
**    NumConsecutiveAdjElems - find the number of consecutive adjacent 
**                             elements 
**
**      int NumConsecutiveAdjElems(const int id) const
**
**        id - (in)  node id 
**
**      Description: Find the size of the largest sequence of elements 
**          that are adjacent to each other and the given node. 
**
**      Return Value: greatest number of adjacent elements 
**
**
** -- */

int CArbMshTopo2D::NumConsecutiveAdjElems(const int id) const
{
    ArbVtxEntry **list, *ptr ;
    int max_num = 0 ;
    int num = 0 ;

    // first check to see if this node is in the hash table.

    list = vtx_table->Fetch(id) ;
    //assert(list != 0) ;
    if (list == 0)
     return -1;

    // loop around the vertex counting the elements

    for (ptr = *list ; ptr != 0 ; ptr = ptr->next) 
    {
        if (ptr->ccw_elem == NO_ELEM) {
            if (num > max_num) max_num = num ;
            num = 0 ;
        } else {
            ++num ;
        }
    }
    return(max_num) ;
}




// %(CArbMshTopo2D::BoundaryNode-bool-|^const-int-const|)
/* ++ ----------------------------------------------------------
**
**    BoundaryNode - boundary node query 
**
**      bool BoundaryNode(const int id) const
**
**        id - (in)  node id 
**
**      Description: Checks to see if a node is on the boundary. 
**
**      Return Value: True if the node is on the boundary 
**
**
** -- */

bool CArbMshTopo2D::BoundaryNode(const int id) const
{
    ArbVtxEntry **list, *ptr ;

    // first check to see if this node is in the hash table.

    list = vtx_table->Fetch(id) ;
    //assert(list != 0) ;
    if (list == 0)
     return false;

    // loop around the vertex counting the elements

    for (ptr = *list ; ptr != 0 ; ptr = ptr->next)
    {
        if (ptr->ccw_elem == NO_ELEM) 
         return(true) ;
    }
    return(false) ;
}




// %(CArbMshTopo2D::DeleteNode-void-|-int-|)
/* ++ ----------------------------------------------------------
**
**    DeleteNode - delete a node 
**
**      void DeleteNode(int node)
**
**        node - (in)  node id 
**
**      Description: Deletes a node from the mesh. 
**
**
** -- */

void CArbMshTopo2D::DeleteNode(const int vtx) 
{
    ArbVtxEntry **list ;

    // first check to see if this node is in the hash table.

    if ((list = vtx_table->Fetch(vtx)) == 0) return ;

    // delete all entries this node's list

    ArbVtxEntry *ptr = *list ;
    while (ptr != 0) {
        ArbVtxEntry *next = ptr->next ;
        DeleteVtxEntry(ptr) ;
        ptr = next ;
    }
}




// %(CArbMshTopo2D::EdgeList-CArbQueue-|<ArbEdge>*)
/* ++ ----------------------------------------------------------
**
**    EdgeList - returns all edges in the mesh 
**
**      CArbQueue <ArbEdge>*EdgeList()
**
**      Description: This function returns all the edges in the mesh. 
**
**      Return Value: a queue object containing all the edges in the 
**          mesh 
**
**
** -- */

CArbQueue<CArbMshTopo2D::ArbEdge> *CArbMshTopo2D::EdgeList()
{
    CArbQueue<CArbMshTopo2D::ArbEdge> *elist = 
        new CArbQueue<CArbMshTopo2D::ArbEdge>() ;

    CArbHashTableIterator<int,ArbVtxEntry *> iter(vtx_table) ;

    for (iter.First() ; iter.More() ; ++iter) {
        ArbVtxEntry *ptr ;
        for (ptr = *(iter.Entry()) ; ptr != 0 ; ptr = ptr->next) {
            if (ptr->vtx_num < iter.Key()) {
                ArbEdge edge ;
                edge.nd0 = ptr->vtx_num ;
                edge.nd1 = iter.Key() ;
                elist->AddToBack(edge) ;
            }
        }
    }

    return(elist) ;
}


void CArbMshTopo2D::RebuildElemList()
{
    // because this implementation does not deal with multiply
    // connected regions properly it can be the case that
    // elements get deleted when realy some edges still point
    // to them.  This function goes through all the edges in
    // in the topology and builds a new element list.

    delete elem_table ;
    elem_table = new CArbHashTable<int,int> ;

    CArbHashTableIterator<int,ArbVtxEntry *> iter(vtx_table) ;

    for (iter.First() ; iter.More() ; ++iter) {
        ArbVtxEntry *ptr ;
        for (ptr = *(iter.Entry()) ; ptr != 0 ; ptr = ptr->next) {
            if (ptr->ccw_elem != NO_ELEM) {
                elem_table->Store(ptr->ccw_elem,iter.Key()) ;
            }
        }
    }
}


// %(CArbMshTopo2D::NumElemNodes-int-|^const-int-const|)
/* ++ ----------------------------------------------------------
**
**    NumElemNodes - number of element nodes 
**
**      int NumElemNodes(const int elem) const
**
**        elem - (in)  element id 
**
**      Description: Finds the number of nodes on a given element. 
**
**      Return Value: number of nodes on an element 
**
**
** -- */

int CArbMshTopo2D::NumElemNodes(const int elem) const
{
    int num_nodes = 0 ;

    if (elem == NO_ELEM)  return(0) ;

    // get a pointer to the first node

    int *node_ptr = elem_table->Fetch(elem) ;
    assert(node_ptr != 0) ;

    int start_node ;
    start_node = *node_ptr ;
    ++num_nodes ;

    // Now we traverse around the element

    int adj_node = GetAdjVtxAlongElem(start_node,elem) ;
    while(adj_node != start_node) {
        ++num_nodes ;
        adj_node = GetAdjVtxAlongElem(adj_node,elem) ;
    }
    return(num_nodes) ;
}




// %(CArbMshTopo2D::GetElemNodes-void-|^const-int-const|-int-|*-int-|*) 
/* ++ ----------------------------------------------------------
**
**    GetElemNodes - get an element's nodes 
**
**      void GetElemNodes(
**              const int elem,
**              int       *num_nodes,
**              int       *nodes) const
**
**        elem      - (in)  element id 
**        num_nodes - (out) number of nodes 
**        nodes     - (out) node id's 
**
**      Description: This function returns the id's of the nodes on a 
**          given element 
**
**
** -- */

int CArbMshTopo2D::GetElemNodes(const int elem,
                                 int *num_nodes,
                                 int *nodes) const
{
    *num_nodes = 0 ;

    if (elem == NO_ELEM)  
     return 0;

    // get a pointer to the first node

    int *node_ptr = elem_table->Fetch(elem) ;
    //assert(node_ptr != 0) ;
    if (node_ptr == 0)
     return 0;

    int start_node ;
    nodes[0] = start_node = *node_ptr ;
    ++(*num_nodes) ;

    // Now we traverse around the element

    int adj_node = GetAdjVtxAlongElem(start_node,elem) ;
    if (adj_node == NO_NODE)
     return 0;
    while(adj_node != start_node) 
    {
        if (*num_nodes == 4)
         return 0;
        nodes[(*num_nodes)] = adj_node ;
        ++(*num_nodes) ;
        adj_node = GetAdjVtxAlongElem(adj_node,elem) ;
        if (adj_node == NO_NODE)
         return 0;
    }
    return 1;
}




// %(CArbMshTopo2D::GetElemNodes-void-|^const-int-const|-int-const|-int-|*-int-|*)
/* ++ ----------------------------------------------------------
**
**    GetElemNodes - get an element's nodes 
**
**      void GetElemNodes(
**              const int elem,
**              const int first_node_id,
**              int       *num_nodes,
**              int       *nodes) const
**
**        elem          - (in)  element id 
**        first_node_id - (in)  first node id 
**        num_nodes     - (out) number of nodes 
**        nodes         - (out) node id's 
**
**      Description: This function returns the id's of the nodes on a 
**          given element ordered so that the given node id is first. 
**
**
** -- */

int CArbMshTopo2D::GetElemNodes(const int elem,
                                 const int node_id,
                                 int *num_nodes,
                                 int *nodes) const
{
    nodes[0] = node_id ;
    *num_nodes = 1 ;

    if (elem == NO_NODE)
     return 0;
    // Now we traverse around the element

    int adj_node = GetAdjVtxAlongElem(node_id,elem) ;
    if (adj_node == NO_NODE)
     return 0;
    while(adj_node != node_id) 
    {
        if (*num_nodes == 4)
         return 0;
        nodes[(*num_nodes)] = adj_node ;
        ++(*num_nodes) ;
        adj_node = GetAdjVtxAlongElem(adj_node,elem) ;
        if (adj_node == NO_NODE)
         return 0;
    }
    return 1;
}




// %(CArbMshTopo2D::NumElements-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumElements - number of elements int the mesh 
**
**      int NumElements() const
**
**      Description: This function returns the number of elements in 
**          the mesh. 
**
**      Return Value: number of elements 
**
**
** -- */

int CArbMshTopo2D::NumElements() const
{
    return(elem_table->NumEntries()) ;
}




// %(CArbMshTopo2D::GetElemList-int-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetElemList - get a list of element id's 
**
**      int *GetElemList() const
**
**      Description: This function returns a list of the id's for all 
**          elements in the mesh. Ownership of the memory passes to the 
**          client and must be freed with a call to delete []. 
**
**      Return Value: a list (array) of the id's of all elements in the 
**          mesh 
**
**
** -- */

int *CArbMshTopo2D::GetElemList() const
{

    return(elem_table->GetKeyList()) ;
}




// %(CArbMshTopo2D::GetElemSet-CArbElemSet-|*^const) 
/* ++ ----------------------------------------------------------
**
**    GetElemSet - get a set of the element id's 
**
**      CArbElemSet *GetElemSet() const
**
**      Description: This function returns the id's of all the elements 
**          in the mesh stored in a set object. Ownership of the set 
**          object passes to the client who must make a call to delete. 
**
**      Return Value: a set object containing the id's of the elements 
**
**
** -- */

CArbElemSet *CArbMshTopo2D::GetElemSet() const
{
    CArbElemSet *set = new CArbElemSet() ;
    CArbHashTableIterator<int,int> iter(elem_table) ;
    for (iter.First() ; iter.More() ; ++iter)
        set->Store(iter.Key(),1) ;
    return(set) ;
}




// %(CArbMshTopo2D::NumNodes-int-|^const)
/* ++ ----------------------------------------------------------
**
**    NumNodes - finds the number of nodes in the mesh 
**
**      int NumNodes() const
**
**      Description: This function returns the number of nodes in the 
**          mesh. 
**
**      Return Value: number of nodes in the mesh 
**
**
** -- */

int CArbMshTopo2D::NumNodes() const
{
    return(vtx_table->NumEntries()) ;
}




// %(CArbMshTopo2D::GetNodeList-int-|*^const) 
/* ++ ----------------------------------------------------------
**
**    GetNodeList - get a list of node id's 
**
**      int *GetNodeList() const
**
**      Description: This function returns a list of the id's for all 
**          nodes in the mesh. Ownership of the memory passes to the 
**          client and must be freed with a call to delete []. 
**
**      Return Value: a list (array) of the id's of all nodes in the 
**          mesh 
**
**
** -- */

int *CArbMshTopo2D::GetNodeList() const
{

    return(vtx_table->GetKeyList()) ;
}

/* ============================================================
   ====================== Iterator stuff ======================
   ============================================================ */

#define NULL_VTX 100000000


// %(CArbTopoAdjVtxCyclicIterator::Prev-void-|)
/* ++ ----------------------------------------------------------
**
**    Prev - point to the previous adjacent vertex 
**
**      void Prev()
**
**      Description: This function sets the iterator to point to the 
**          previous adjacent vertex. 
**
**
** -- */

void CArbTopoAdjVtxCyclicIterator::Prev()
{
    CArbMshTopo2D::ArbVtxEntry *cur = ptr ;
    while (1) {
        if (cur->next == 0) {
            CArbMshTopo2D::ArbVtxEntry *first = topo->FirstVtx(vtx) ;
            if (first == ptr) {
                ptr = cur ;
                break ;
            } else {
                cur = first ;
            }
        } else if (cur->next == ptr) {
            ptr = cur ;
            break ;
        } else {
            cur = cur->next ;
        }
    }
}



// ----------------------------------------------------------------
// vertex on element iterators


// %(CArbTopoVtxOnElemIterator::CArbTopoVtxOnElemIterator-constructor-|-CArbMshTopo2D-const|*-int-const|)
/* ++ ----------------------------------------------------------
**
**    CArbTopoVtxOnElemIterator - constructor 
**
**      CArbTopoVtxOnElemIterator(
**              const CArbMshTopo2D *aTopo,
**              const int           elem)
**
**        aTopo - (in)  an ArbMshTopo2D object 
**        elem  - (in)  input element 
**
**      Description: This is a constructor for a 
**          ArbTopoAdjVtxOnElemIterator object. 
**
**
** -- */

CArbTopoVtxOnElemIterator::
CArbTopoVtxOnElemIterator(const CArbMshTopo2D *aTopo,
                          const int ielem)
{
    elem = ielem ;
    if (elem == NO_ELEM) cur = NULL_VTX ;

    int *node_ptr = aTopo->elem_table->Fetch(elem) ;
    assert(node_ptr != 0) ;

    topo = aTopo ;
    cur = *node_ptr ;
    int tmp = aTopo->GetAdjVtxAlongElem(cur,elem) ;
    prev = aTopo->GetCCWNode(tmp,cur) ;
    first = cur ;
}




// %(CArbTopoVtxOnElemIterator::CArbTopoVtxOnElemIterator-constructor-|-CArbMshTopo2D-const|*-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    CArbTopoVtxOnElemIterator - constructor 
**
**      CArbTopoVtxOnElemIterator(
**              const CArbMshTopo2D *aTopo,
**              const int           elem,
**              const int           node)
**
**        aTopo - (in)  an ArbMshTopo2D object 
**        elem  - (in)  input element 
**        node  - (in)  first vertex 
**
**      Description: This is a constructor for a 
**          ArbTopoAdjVtxOnElemIterator object. 
**
**
** -- */

CArbTopoVtxOnElemIterator::
CArbTopoVtxOnElemIterator(const CArbMshTopo2D *aTopo,
                          const int ielem,
                          const int node)
{
    elem = ielem ;
    if (elem == NO_ELEM) cur = NULL_VTX ;

#ifdef _DEBUG
    int *node_ptr = aTopo->elem_table->Fetch(elem) ;
    assert(node_ptr != 0) ;
#endif

    topo = aTopo ;
    cur = node ;
    int tmp = aTopo->GetAdjVtxAlongElem(cur,elem) ;
    prev = topo->GetCCWNode(tmp,cur) ;
    first = cur ;
}




// %(CArbTopoVtxOnElemIterator::First-void-|) 
/* ++ ----------------------------------------------------------
**
**    First - point to the first vertex 
**
**      void First()
**
**      Description: This function sets the iterator to point to the 
**          first vertex on the element. 
**
**
** -- */

void CArbTopoVtxOnElemIterator::First()
{
    cur = first ;
    int tmp = topo->GetAdjVtxAlongElem(cur,elem) ;
    prev = topo->GetCCWNode(tmp,cur) ;
}




// %(CArbTopoVtxOnElemIterator::Next-void-|) 
/* ++ ----------------------------------------------------------
**
**    Next - point to the next vertex 
**
**      void Next()
**
**      Description: This function sets the iterator to point to the 
**          next vertex on the element. 
**
**
** -- */

void CArbTopoVtxOnElemIterator::Next()
{
    int next = topo->GetCWNode(cur,prev) ;

    if (next == first) {
        cur = NULL_VTX ;
    } else {
        prev = cur ;
        cur = next ;
    }
}




// ----------------------------------------------------------------
// vertex on element cyclic iterator


// %(CArbTopoVtxOnElemCyclicIterator::CArbTopoVtxOnElemCyclicIterator-constructor-|-CArbMshTopo2D-const|*-int-const|) 
/* ++ ----------------------------------------------------------
**
**    CArbTopoVtxOnElemCyclicIterator - constructor 
**
**      CArbTopoVtxOnElemCyclicIterator(
**              const CArbMshTopo2D *aTopo,
**              const int           elem)
**
**        aTopo - (in)  an ArbMshTopo2D object 
**        elem  - (in)  input element 
**
**      Description: This is a constructor for a 
**          ArbTopoVtxOnElemCyclicIterator object. 
**
**
** -- */

CArbTopoVtxOnElemCyclicIterator::
CArbTopoVtxOnElemCyclicIterator(const CArbMshTopo2D *aTopo,
                                const int ielem)
{
    elem = ielem ;
    if (elem == NO_ELEM) cur = NULL_VTX ;

    int *node_ptr = aTopo->elem_table->Fetch(elem) ;
    assert(node_ptr != 0) ;

    topo = aTopo ;
    cur = *node_ptr ;
    int tmp = aTopo->GetAdjVtxAlongElem(cur,elem) ;
    prev = aTopo->GetCCWNode(tmp,cur) ;
    first = cur ;
}




// %(CArbTopoVtxOnElemCyclicIterator::CArbTopoVtxOnElemCyclicIterator-constructor-|-CArbMshTopo2D-const|*-int-const|-int-const|)
/* ++ ----------------------------------------------------------
**
**    CArbTopoVtxOnElemCyclicIterator - constructor 
**
**      CArbTopoVtxOnElemCyclicIterator(
**              const CArbMshTopo2D *aTopo,
**              const int           elem,
**              const int           node)
**
**        aTopo - (in)  an ArbMshTopo2D object 
**        elem  - (in)  input element 
**        node  - (in)  first vertex 
**
**      Description: This is a constructor for a 
**          ArbTopoVtxOnElemCyclicIterator object. 
**
**
** -- */

CArbTopoVtxOnElemCyclicIterator::
CArbTopoVtxOnElemCyclicIterator(const CArbMshTopo2D *aTopo,
                                const int ielem,
                                const int node)
{
    elem = ielem ;
    if (elem == NO_ELEM) cur = NULL_VTX ;

#ifdef _DEBUG
    int *node_ptr = aTopo->elem_table->Fetch(elem) ;
    assert(node_ptr != 0) ;
#endif

    topo = aTopo ;
    cur = node ;
    int tmp = aTopo->GetAdjVtxAlongElem(cur,elem) ;
    prev = topo->GetCCWNode(tmp,cur) ;
    first = cur ;
}




// %(CArbTopoVtxOnElemCyclicIterator::First-void-|)
/* ++ ----------------------------------------------------------
**
**    First - point to the first vertex 
**
**      void First()
**
**      Description: This function sets the iterator to point to the 
**          first vertex on the element. 
**
**
** -- */

void CArbTopoVtxOnElemCyclicIterator::First()
{
    cur = first ;
    int tmp = topo->GetAdjVtxAlongElem(cur,elem) ;
    prev = topo->GetCCWNode(tmp,cur) ;
}



// %(CArbTopoVtxOnElemCyclicIterator::Next-void-|)
/* ++ ----------------------------------------------------------
**
**    Next - point to the next vertex 
**
**      void Next()
**
**      Description: This function sets the iterator to point to the 
**          next vertex on the element. 
**
**
** -- */

void CArbTopoVtxOnElemCyclicIterator::Next()
{
    int next = topo->GetCWNode(cur,prev) ;
    prev = cur ;
    cur = next ;
}




// %(CArbTopoVtxOnElemCyclicIterator::Prev-void-|)
/* ++ ----------------------------------------------------------
**
**    Prev - point to the previous vertex 
**
**      void Prev()
**
**      Description: This function sets the iterator to point to the 
**          previous vertex on the element. 
**
**
** -- */

void CArbTopoVtxOnElemCyclicIterator::Prev()
{
    int next = topo->GetCCWNode(cur,prev) ;
    cur = prev ;
    prev = next ;
}



// ----------------------------------------------------------------
// vertex on boundary iterator

// %(CArbTopoVtxOnBdryIterator::CArbTopoVtxOnBdryIterator-constructor-|-CArbMshTopo2D-const|*-int-const|)
/* ++ ----------------------------------------------------------
**
**    CArbTopoVtxOnBdryIterator - constructor 
**
**      CArbTopoVtxOnBdryIterator(
**              const CArbMshTopo2D *aTopo,
**              const int           node)
**
**        aTopo - (in)  an ArbMshTopo2D object 
**        node  - (in)  starting vertex 
**
**      Description: This is a constructor for a 
**          ArbTopoVtxOnBdryIterator object. 
**
**
** -- */

CArbTopoVtxOnBdryIterator::
CArbTopoVtxOnBdryIterator(const CArbMshTopo2D *aTopo,
                          const int inode)
{
    CArbTopoAdjVtxIterator iter(aTopo,inode) ;
    for(iter.First() ; iter.More() ; ++iter) {
        if (iter.CcwElem() == NO_ELEM) {
            owner = aTopo ;
            cur = inode ;
            first = iter.AdjVtx() ;
            prev = first ;
            node = inode ;
            return ;
        }
    }
    cur = NULL_VTX ;
}




// %(CArbTopoVtxOnBdryIterator::First-void-|) 
/* ++ ----------------------------------------------------------
**
**    First - point to the first vertex 
**
**      void First()
**
**      Description: This function sets the iterator to point to the 
**          first vertex on the boundary. 
**
**
** -- */

void CArbTopoVtxOnBdryIterator::First()
{
    cur = node ;
    prev = first ;
}




// %(CArbTopoVtxOnBdryIterator::Next-void-|)
/* ++ ----------------------------------------------------------
**
**    Next - point to the next vertex 
**
**      void Next()
**
**      Description: This function sets the iterator to point to the 
**          next vertex on the boundary. 
**
**
** -- */

void CArbTopoVtxOnBdryIterator::Next()
{
    if (cur == first) {
        cur = NULL_VTX ;
    } else {
        int next = owner->GetCCWBdryNode(prev,cur) ;
        prev = cur ;
        cur = next ;
    }
}



// ----------------------------------------------------------------
// edge on element iterator

// %(CArbTopoEdgeOnElemIterator::CArbTopoEdgeOnElemIterator-constructor-|-CArbMshTopo2D-const|*-int-const|)
/* ++ ----------------------------------------------------------
**
**    CArbTopoEdgeOnElemIterator - constructor 
**
**      CArbTopoEdgeOnElemIterator(
**              const CArbMshTopo2D *aTopo,
**              const int           elem)
**
**        aTopo - (in)  an ArbMshTopo2D object 
**        elem  - (in)  input element 
**
**      Description: This is a constructor for a 
**          ArbTopoEdgeOnElemIterator object. 
**
**
** -- */

CArbTopoEdgeOnElemIterator::
CArbTopoEdgeOnElemIterator(const CArbMshTopo2D *aTopo,
                           const int ielem)
{
    if (ielem == NO_ELEM) {
        cur = NULL_VTX ;
        return ;
    }
    elem = ielem ;
    int *node_ptr = aTopo->elem_table->Fetch(elem) ;
    if (node_ptr == 0) cur = NULL_VTX ;
    first = *node_ptr ;
    cur = aTopo->GetAdjVtxAlongElem(first,elem) ;
    prev = first ;
    adj = cur ;
    owner = aTopo ;
}

CArbTopoEdgeOnElemIterator::
CArbTopoEdgeOnElemIterator(const CArbMshTopo2D *aTopo,
                           const int ielem,
                           const CArbMshTopo2D::ElemLoop &loop)
{
    if (ielem == NO_ELEM) {
        cur = NULL_VTX ;
        return ;
    }
    elem = ielem ;
    first = loop.vtx0 ;
    cur = loop.vtx1 ;
    prev = first ;
    adj = cur ;
    owner = aTopo ;
}




// %(CArbTopoEdgeOnElemIterator::First-void-|) 
/* ++ ----------------------------------------------------------
**
**    First - point to the first edge 
**
**      void First()
**
**      Description: This function sets the iterator to point to the 
**          first edge on the element. 
**
**
** -- */

void CArbTopoEdgeOnElemIterator::First()
{
    cur = owner->GetAdjVtxAlongElem(first,elem) ;
    prev = first ;
    adj = cur ;
}




// %(CArbTopoEdgeOnElemIterator::Next-void-|)
/* ++ ----------------------------------------------------------
**
**    Next - point to the next edge 
**
**      void Next()
**
**      Description: This function sets the iterator to point to the 
**          next edge on the element. 
**
**
** -- */

void CArbTopoEdgeOnElemIterator::Next()
{
    if (cur == first) {
        cur = NULL_VTX ;
        return ;
    }
    int next = owner->GetCWNode(cur,prev) ;
    prev = cur ;
    cur = next ;
}


// ----------------------------------------------------------------
// edge on element cyclic iterator

// %(CArbTopoEdgeOnElemCyclicIterator::CArbTopoEdgeOnElemCyclicIterator-constructor-|-CArbMshTopo2D-const|*-int-const|)
/* ++ ----------------------------------------------------------
**
**    CArbTopoEdgeOnElemCyclicIterator - constructor 
**
**      CArbTopoEdgeOnElemCyclicIterator(
**              const CArbMshTopo2D *aTopo,
**              const int           elem)
**
**        aTopo - (in)  an ArbMshTopo2D object 
**        elem  - (in)  element 
**
**      Description: This is a constructor for a 
**          ArbTopoVtxOnBdryIterator object. 
**
**
** -- */

CArbTopoEdgeOnElemCyclicIterator::
CArbTopoEdgeOnElemCyclicIterator(const CArbMshTopo2D *aTopo,
                                 const int ielem)
{
    if (ielem == NO_ELEM) {
        cur = NULL_VTX ;
        return ;
    }
    elem = ielem ;
    int *node_ptr = aTopo->elem_table->Fetch(elem) ;
    if (node_ptr == 0) cur = NULL_VTX ;
    first = *node_ptr ;
    cur = aTopo->GetAdjVtxAlongElem(first,elem) ;
    prev = first ;
    adj = cur ;
    owner = aTopo ;
}

CArbTopoEdgeOnElemCyclicIterator::
CArbTopoEdgeOnElemCyclicIterator(const CArbMshTopo2D *aTopo,
                                 const int ielem,
                                 const CArbMshTopo2D::ElemLoop &loop)
{
    if (ielem == NO_ELEM) {
        cur = NULL_VTX ;
        return ;
    }
    elem = ielem ;
    first = loop.vtx0 ;
    cur = loop.vtx1 ;
    prev = first ;
    adj = cur ;
    owner = aTopo ;
}




// %(CArbTopoEdgeOnElemCyclicIterator::First-void-|)
/* ++ ----------------------------------------------------------
**
**    First - point to the first edge 
**
**      void First()
**
**      Description: This function sets the iterator to point to the 
**          first edge on the element. 
**
**
** -- */

void CArbTopoEdgeOnElemCyclicIterator::First()
{
    cur = adj ;
    prev = first ;
}




// %(CArbTopoEdgeOnElemCyclicIterator::Next-void-|) 
/* ++ ----------------------------------------------------------
**
**    Next - point to the next edge 
**
**      void Next()
**
**      Description: This function sets the iterator to point to the 
**          next edge on the element. 
**
**
** -- */

void CArbTopoEdgeOnElemCyclicIterator::Next()
{
    int next = owner->GetCWNode(cur,prev) ;
    prev = cur ;
    cur = next ;
}




// %(CArbTopoEdgeOnElemCyclicIterator::Prev-void-|) 
/* ++ ----------------------------------------------------------
**
**    Prev - point to the previous edge 
**
**      void Prev()
**
**      Description: This function sets the iterator to point to the 
**          previous edge on the element. 
**
**
** -- */

void CArbTopoEdgeOnElemCyclicIterator::Prev()
{
    int next = owner->GetCCWNode(cur,prev) ;
    cur = prev ;
    prev = next ;
}

// ----------------------------------------------------------------
// element iterator

CArbTopoElemIterator::
CArbTopoElemIterator(const CArbMshTopo2D *aTopo)
{
    Eiter = new CArbHashTableIterator<int,int>(aTopo->elem_table) ;
}

CArbTopoElemIterator::~CArbTopoElemIterator()
{
    delete Eiter ;
}

CArbTopoSortedElemIterator::
CArbTopoSortedElemIterator(const CArbMshTopo2D *aTopo)
{
    num = aTopo->elem_table->NumEntries() ;
    ids = new int[num] ;
    CArbHashTableIterator<int,int> iter(aTopo->elem_table) ;
    cur = 0 ;
    for (iter.First() ; iter.More() ; ++iter) {
        ids[cur] = iter.Key() ;
        ++cur ;
    }
    cur = 0 ;

    // dumb bubble sort for debugging

    for (int i=0 ; i<(num-1) ; ++i) {
        for (int j=i ; j<num ; ++j) {
            if (ids[j] < ids[i]) {
                int tmp = ids[i] ;
                ids[i] = ids[j] ;
                ids[j] = tmp ;
            }
        }
    }
}

CArbTopoSortedElemIterator::~CArbTopoSortedElemIterator()
{
    delete ids ;
}

// ----------------------------------------------------------------
// edge iterator

CArbTopoEdgeIterator::
CArbTopoEdgeIterator(const CArbMshTopo2D &aTopo) :
    iter(aTopo.vtx_table)
{
    ptr = *iter.Entry() ;
    while (ptr->vtx_num < iter.Key()) {
        ptr = ptr->next ;
        if (ptr == 0) {
            ++iter ;
            if (!iter.More()) return ;
            ptr = *iter.Entry() ;
        }
    }
}

void CArbTopoEdgeIterator::Next()
{
    while(1) {
        ptr = ptr->next ;
        if (ptr == 0) {
            ++iter ;
            if (!iter.More()) return ;
            ptr = *iter.Entry() ;
        }
        if (ptr->vtx_num > iter.Key()) break ;
    }
}

// ----------------------------------------------------------------
// Loop on Elem iterator

CArbTopoLoopOnElemIterator::
CArbTopoLoopOnElemIterator(const CArbMshTopo2D *aTopo,
                           const int elem)
{
    owner = aTopo ;

    CArbHashTable<EdgeKey,int> edge_hash ;

    // loop through all the verts in the topology

    CArbHashTableIterator<int,CArbMshTopo2D::ArbVtxEntry *>
        vtx_iter(owner->vtx_table) ;
    for (vtx_iter.First() ; vtx_iter.More() ; ++vtx_iter) {

        // loop through all verts adjacent to this one

        CArbMshTopo2D::ArbVtxEntry *ptr ;
        for (ptr = owner->FirstVtx(vtx_iter.Key()) ; ptr != 0 ;
             ptr = owner->NextVtx(ptr)) {

            // check to see if we are on the element.  If so
            // put this edge in the has table

            if (ptr->ccw_elem == elem) {
                edge_hash.Store(EdgeKey(vtx_iter.Key(),ptr->vtx_num),1) ;
            }
        }
    }

    // identify the loops

    EdgeKey *key ;
    int *ignore ;
    while (edge_hash.FetchAny(&key,&ignore)) {

        // create a new loop
        CArbMshTopo2D::ElemLoop loop ;
        loop.vtx0 = key->Vtx0() ;
        loop.vtx1 = key->Vtx1() ;
        loops.InsertAtEnd(loop) ;

        // set up for moving around the loop
        int first = key->Vtx0() ;
        int cur = key->Vtx1() ;
        int prev = first ;

        // delete this edge from the hash table
        edge_hash.Remove(EdgeKey(first,cur)) ;

        // move around the loop
        while (cur != first) {
            int next = owner->GetCWNode(cur,prev) ;
            prev = cur ;
            cur = next ;
            edge_hash.Remove(EdgeKey(prev,cur)) ;
        } 
    }
}
