//
// CArbMshTopo2D Class header file
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
//   $Revision: 1.15 $  $Date: 2002/07/25 13:16:14 $  $Author: wash $
//

#ifndef ArbMshTopo2D_h
#define ArbMshTopo2D_h

#include "ArbHashTable.hpp"
#include "ArbQueue.hpp"
#include "ArbArray.hpp"

#define NO_ELEM 100000000
#define NO_NODE 100000000

#define TOPO_CACHE_BLOCK_SIZE 1000

typedef CArbHashTable<int,int> CArbElemSet ;
typedef CArbHashTableIterator<int,int> CArbElemSetIterator ;

// CArbMshTopo2D
///////////////////////////////////////////////////////////////
class CArbMshTopo2D
{
    public:
    // public data structures

    struct ArbVtxEntry
    {
      int vtx_num;
      int ccw_elem;
      struct ArbVtxEntry *next;
    };

    struct ArbEdge
    {
      int nd0;
      int nd1;
      int elem0;
      int elem1;
    } ;

    struct ElemLoop
    {
      int vtx0,vtx1 ;
    } ;
           
    // constructors and destructors

    CArbMshTopo2D () ;
    CArbMshTopo2D (const CArbMshTopo2D &other) ;
    CArbMshTopo2D operator = (const CArbMshTopo2D &other) ;
    ~CArbMshTopo2D () ;

    // insert and removal routines

    void InsertElement(const int elem_id,const int num_nodes, const int *const nodes);

    int InsertTriangle(const int elem_id, const int nd0, const int nd1, const int nd2);

    void InsertCollapsedElement(const int elem_id, const int num_nodes, 
                                const int *const nodes, const int parent_elem = -1);

    int DeleteElement (const int num_nodes, const int *const nodes);

    void DeleteElement(const int elem_id);

    int DeleteTriangle (const int nd0, const int nd1, const int nd2);

    int InsertAngle (const int elem, const int vtx, const int prev,
                      const int next, const int parent_elem = -1);

    int DeleteAngle (const int vtx, const int prev, const int next);

    void DeleteElementRef(const int elem_id);

    void SplitEdge (const int nd0, const int nd1, const int new_node);

    void SplitEdgeElem (const int nd0, const int nd1, const int new_node,
                        const int elem, const bool ccw_flag);

    bool Splice (const int elem, const int start, const int finish, CArbMshTopo2D *other,
                 const int other_elem, const int other_start, const int other_finish);

    void AddEdgeVertex( int vtx,int CW_vtx,int new_vtx);

    void AddEdgeTwoVertex( int vtx0,int vtx1,int elem);

    int AddFaceEdge (int vtx0,int vtx1,int CW_vtx0,int CW_vtx1, int new_elem) ;

    int FindCrossedFace (int vtx0,int vtx1, int CW_vtx0,int CW_vtx1) const ;

    void AddFaceTwoEdgeTwoVertex (int vtx0,int vtx1, int new_elem,int in_elem) ;

    void AddFaceTwoEdgeOneVertex( int vtx0,int vtx1, int CW_vtx0,int new_elem) ;

    void TearEdgeAddFace (int vtx0,int vtx1,int new_elem) ;

    void SplitVertexMakeEdge (int vtx,int new_vtx, int face0,int face1) ;

    // insert and removal routines

    void DeleteNode (int node) ;

    // query routines

    bool HasVtx (const int vtx) const ;

    ArbEdge OppositeEdge (const int elem_id,const int node_id) const;

    int OppositeNode (const int elem_id,const ArbEdge *edge) const;

    bool ElemHasNode (const int elem_id, const int node_id) const;

    ArbEdge NextCWEdge (const int elem_id, const ArbEdge *edge) const;

    ArbEdge NextCCWEdge (const int elem_id, const ArbEdge *edge) const;

    int GetCWElem (const int node_0, const int node_1) const;

    int GetCCWElem (const int node_0, const int node_1) const;

    int GetCWBdryNode (const int node_0,const int node_1) const;

    int GetCCWBdryNode (const int node_0, const int node_1) const;

    int GetCWNode (const int node_0, const int node_1) const;

    int GetCCWNode (const int node_0, const int node_1) const;

    int GetDblEdgeElem (int node_0,int node_1) const;

    bool HasDblEdge (int node_0,int node_1) const;

    void GetElemsAboutEdge (const int node_0,const int node_1,int *elem_0,int *elem_1) const;

    int NumAdjElems (const int id) const;
    int NumConsecutiveAdjElems (const int id) const;
    bool BoundaryNode (const int id) const;

    CArbQueue<ArbEdge> *EdgeList();

    int NumElemNodes (const int elem) const;

    int GetElemNodes (const int elem, int *num_nodes, int *nodes) const;

    int GetElemNodes (const int elem, const int first_node_id, int *num_nodes, int *nodes) const;

    void RebuildElemList ( );

    int NumElements ( ) const;
    int *GetElemList ( ) const;
    CArbElemSet  *GetElemSet( ) const;

    int NumNodes ( ) const;
    int *GetNodeList ( ) const;


    private:

    struct ArbVtxCache
    {
      struct ArbVtxCache *next ;
      struct ArbVtxEntry entries[TOPO_CACHE_BLOCK_SIZE] ;
    } ;

// member variable

    CArbHashTable <int,ArbVtxEntry *> *vtx_table ;
    CArbHashTable <int,int> *elem_table ;

    ArbVtxEntry *MaxAddress ;
    int MaxVtx ;

    ArbVtxEntry *FreeList ;
    ArbVtxCache *CacheList ;

// private member functions

    ArbVtxEntry *NewVtxEntry ( ) ;
    void DeleteVtxEntry (ArbVtxEntry *entry) ;
    void DeleteEntryCache ( ) ;

    ArbVtxEntry *FirstVtx (const int vtx) const ;
    ArbVtxEntry *NextVtx (const ArbVtxEntry *ptr) const { return(ptr->next) ; } ;
    ArbVtxEntry *EolVtx() { return(0) ; } ;

    void SplitEdgeHelp (ArbVtxEntry *ptr0, ArbVtxEntry *ptr1,
                        const int nd0, const int nd1, const int new_node);

    int GetCWElemPtr(ArbVtxEntry *ptr, ArbVtxEntry **list);

    int GetAdjVtxAlongElem(int node, int elem_id) const;

    friend class CArbTopoAdjVtxIterator ;
    friend class CArbTopoAdjVtxCyclicIterator ;
    friend class CArbTopoVtxOnElemIterator ;
    friend class CArbTopoVtxOnElemCyclicIterator ;
    friend class CArbTopoEdgeOnElemIterator ;
    friend class CArbTopoEdgeOnElemCyclicIterator ;
    friend class CArbTopoElemIterator ;
    friend class CArbTopoSortedElemIterator ;
    friend class CArbTopoVtxIterator ;
    friend class CArbTopoEdgeIterator ;
    friend class CArbTopoLoopOnElemIterator ;
    friend class CArbTopoEdgeOnElemLoopIterator ;

} ;

// ----------------------------------------------------------------
// New style iterators

#define NULL_VTX 100000000

// CArbTopoAdjVtxIterator
///////////////////////////////////////////////////////////////
class CArbTopoAdjVtxIterator
{
    public:

    CArbTopoAdjVtxIterator (const CArbMshTopo2D *aTopo, const int ivtx) :
                            topo(aTopo), vtx(ivtx) 
        { ptr = topo->FirstVtx(vtx) ; } ;

    void First  ( )        { ptr = topo->FirstVtx(vtx) ; };
    void Next   ( )        { if (ptr) ptr = ptr->next ; };
    bool More   ( )        { return(ptr != 0) ; };
    void NewVtx (int ivtx) { vtx = ivtx ; ptr = topo->FirstVtx(vtx) ; };
    int AdjVtx  ( )        { if (ptr != 0)
                              return(ptr->vtx_num);
                             return -1;
                           };
    int CcwElem ( )        { if (ptr != 0)
                              return(ptr->ccw_elem) ; 
                             return -1;
                           };
    int &CcwElemRef ( )    { if (ptr != 0)
                              return (ptr->ccw_elem); 
                             static int val = -1;
                             return val;
                           };

    void operator ++ ( )     { Next() ; };
    void operator ++ (int ) { Next() ; };

    private:

    const CArbMshTopo2D *topo;
    CArbMshTopo2D::ArbVtxEntry *ptr;
    int vtx ;
};

// CArbTopoAdjVtxCyclicIterator
///////////////////////////////////////////////////////////////
class CArbTopoAdjVtxCyclicIterator
{
    public:

    CArbTopoAdjVtxCyclicIterator (const CArbMshTopo2D *aTopo, const int ivtx) :
                                  topo(aTopo), vtx(ivtx)
                            { ptr = topo->FirstVtx(vtx) ; };

    void First   ( )        { ptr = topo->FirstVtx(vtx) ; } ;
    void Next    ( )        { ptr = ptr->next ; if (!ptr) ptr = topo->FirstVtx(vtx) ; };
    void Prev    ( );
    bool IsValid ( )        { return(ptr != 0) ; };
    void NewVtx  (int ivtx) { vtx = ivtx ; ptr = topo->FirstVtx(vtx) ; };
    int AdjVtx   ( )        { if (ptr != 0)
                               return(ptr->vtx_num) ; 
                              return -1;
                            }; 
    int CcwElem  ( )        { if (ptr != 0)
                               return(ptr->ccw_elem) ; 
                              return -1;
                            };
    int &CcwElemRef ( )     { if (ptr != 0)
                               return(ptr->ccw_elem) ; 
                              static int val=-1;
                              return val;
                            };

    void operator ++ ()      { Next() ; };
    void operator ++ (int ) { Next() ; };
    void operator -- ()      { Prev() ; };
    void operator -- (int ) { Prev() ; };

    private:

    const CArbMshTopo2D *topo ;
    CArbMshTopo2D::ArbVtxEntry *ptr ;
    int vtx ;
} ;


// CArbTopoVtxOnElemIterator
///////////////////////////////////////////////////////////////
class CArbTopoVtxOnElemIterator
{
    public:

    CArbTopoVtxOnElemIterator(const CArbMshTopo2D *aTopo, const int elem);
    CArbTopoVtxOnElemIterator(const CArbMshTopo2D *aTopo, const int elem, const int node);

    void First  ( );
    void Next   ( );
    bool More   ( )            { return(cur != NULL_VTX) ; };
    int Current ( )            { return(cur) ; };

    void operator ++ ()      { Next() ; };
    void operator ++ (int ) { Next() ; };
    int operator * ()        { return(Current()) ; };

    private:

    int cur;
    int first;
    int prev;
    int elem;
    const CArbMshTopo2D *topo;
} ;

// CArbTopoVtxOnElemCyclicIterator 
///////////////////////////////////////////////////////////////
class CArbTopoVtxOnElemCyclicIterator 
{
    public:

    CArbTopoVtxOnElemCyclicIterator(const CArbMshTopo2D *aTopo, const int elem) ;
    CArbTopoVtxOnElemCyclicIterator(const CArbMshTopo2D *aTopo, const int elem, const int node);

    void First   ( );
    void Next    ( );
    void Prev    ( );
    bool IsValid ( )           { return(cur != NULL_VTX) ; };
    int Current  ( )           { return(cur) ; };

    void operator ++ ()      { Next() ; };
    void operator ++ (int ) { Next() ; };
    void operator -- ()      { Prev() ; };
    void operator -- (int ) { Prev()	;	};
    int operator * ()        { return(Current()) ; } ;

    private:

    int cur ;
    int first ;
    int prev ;
    int elem ;
    const CArbMshTopo2D *topo ;
} ;

// CArbTopoVtxOnBdryIterator 
///////////////////////////////////////////////////////////////
class CArbTopoVtxOnBdryIterator
{
    public:

    CArbTopoVtxOnBdryIterator(const CArbMshTopo2D *aTopo, const int node);

    void First   ( );
    void Next    ( );
    bool More    ( )           { return(cur != NULL_VTX) ; };
    int  Current ( )           { return(cur) ; };

    void operator ++ ()      { Next() ; };
    void operator ++ (int ) { Next() ; };
    int operator * ()        { return(Current()) ; };

    private:

    int cur ;
    int first ;
    int prev ;
    int node ;
    const CArbMshTopo2D *owner ;
} ;

// CArbTopoEdgeOnElemIterator 
///////////////////////////////////////////////////////////////
class CArbTopoEdgeOnElemIterator
{
    public:

    CArbTopoEdgeOnElemIterator (const CArbMshTopo2D *aTopo,const int elem);
    CArbTopoEdgeOnElemIterator (const CArbMshTopo2D *aTopo,const int elem,
                                const CArbMshTopo2D::ElemLoop &loop) ;

    void First    ( );
    void Next     ( );
    bool More     ( )        { return(cur != NULL_VTX ) ; };
    int  Current0 ( )        { return(prev) ; };
    int  Current1 ( )        { return(cur) ; };

    void operator ++ ()      { Next() ; };
    void operator ++ (int ) { Next() ; };
    int operator [] (int i)  { return(i==0 ? Current0() : Current1()) ; };

    private:

    int cur;
    int adj;
    int first;
    int prev;
    int elem;
    const CArbMshTopo2D *owner;
} ;

// CArbTopoEdgeOnElemCyclicIterator 
///////////////////////////////////////////////////////////////
class CArbTopoEdgeOnElemCyclicIterator
{
    public:

    CArbTopoEdgeOnElemCyclicIterator (const CArbMshTopo2D *aTopo,const int elem);
    CArbTopoEdgeOnElemCyclicIterator (const CArbMshTopo2D *aTopo,const int elem,
                                      const CArbMshTopo2D::ElemLoop &loop);

    void First    ( );
    void Next     ( );
    void Prev     ( );
    bool IsValid  ( )        { return(cur != NULL_VTX) ; };
    int  Current0 ( )        { return(prev) ; };
    int  Current1 ( )        { return(cur) ; };

    void operator ++ ()      { Next() ; } ;
    void operator ++ (int ) { Next() ; } ;
    void operator -- ()      { Prev() ; } ;
    void operator -- (int ) { Prev() ; } ;
    int operator [] (int i)  { return(i==0 ? Current0() : Current1()) ; } ;

    private:

    int cur ;
    int adj ;
    int first ;
    int prev ;
    int elem ;
    const CArbMshTopo2D *owner ;
} ;

// CArbTopoElemIterator
///////////////////////////////////////////////////////////////
class CArbTopoElemIterator 
{
    public:

    CArbTopoElemIterator (const CArbMshTopo2D *aTopo);
    ~CArbTopoElemIterator ( );

    void First   ( )         { Eiter->First() ; };
    void Next    ( )         { Eiter->Next() ; };
    bool More    ( )         { return(Eiter->More()) ; };
    int  Current ( )         { return(Eiter->Key()) ; };

    void operator ++ ()      { Next() ; };
    void operator ++ (int ) { Next() ; };
    int operator * ()        { return(Current()) ; };

    private:

    CArbHashTableIterator<int,int> *Eiter ;
} ;

// CArbTopoSortedElemIterator
///////////////////////////////////////////////////////////////
class CArbTopoSortedElemIterator 
{
    public:

    CArbTopoSortedElemIterator(const CArbMshTopo2D *aTopo);
    ~CArbTopoSortedElemIterator();

    void First   ( )         { cur = 0 ; };
    void Next    ( )         { ++cur ; };
    bool More    ( )         { return(cur < num) ; };
    int  Current ( )         { return(ids[cur]) ; };

    void operator ++ ()      { Next() ; };
    void operator ++ (int ) { Next() ; };
    int operator * ()        { return(Current()) ; };

    private:

    int num,cur,*ids ;
} ;

// CArbTopoVtxIterator
///////////////////////////////////////////////////////////////
class CArbTopoVtxIterator 
{
    public:

    CArbTopoVtxIterator (const CArbMshTopo2D &aTopo) :
                         iter(aTopo.vtx_table) {};
    void First   ( )    { iter.First() ; };
    void Next    ( )    { iter.Next() ; };
    bool More    ( )    { return(iter.More()) ; };
    int  Current ( )    { return(iter.Key()) ; };

    void operator ++ ()      { Next() ; };
    void operator ++ (int ) { Next() ; };
    int operator * ()        { return(Current()) ; };

    private:

    CArbHashTableIterator<int,CArbMshTopo2D::ArbVtxEntry *> iter ;
} ;

// CArbTopoEdgeIterator
///////////////////////////////////////////////////////////////
class CArbTopoEdgeIterator 
{
    public:

    CArbTopoEdgeIterator(const CArbMshTopo2D &aTopo) ;
    void First    ( )       { iter.First() ; };
    void Next     ( );
    bool More     ( )       { return(iter.More()) ; };
    int  Current0 ( )       { return(ptr->vtx_num) ; };
    int  Current1 ( )       { return(iter.Key()) ; };  

    void operator ++ ()      { Next() ; } ;
    void operator ++ (int ) { Next() ; } ;
    int operator [] (int i)  { return(i==0 ? Current0() : Current1()) ; } ;

    private:

    CArbHashTableIterator<int,CArbMshTopo2D::ArbVtxEntry *> iter ;
    CArbMshTopo2D::ArbVtxEntry *ptr ;
} ;

// CArbTopoLoopOnElemIterator
///////////////////////////////////////////////////////////////
class CArbTopoLoopOnElemIterator
{
    public:

    CArbTopoLoopOnElemIterator(const CArbMshTopo2D *aTopo,
                                const int elem) ;
    void First()             { cur = 0 ; } ;
    void Next()              { ++cur ; } ;
    bool More()              { return(cur < loops.NumEntries() ) ; } ;
    CArbMshTopo2D::ElemLoop Current()
                              { return(loops[cur]) ; } ;

    void operator ++ ()      { Next() ; } ;
    void operator ++ (int ) { Next() ; } ;
    CArbMshTopo2D::ElemLoop operator * ()
                              { return(Current()) ; } ;

    struct EdgeKey 
    {
      int vtx0,vtx1 ;
      EdgeKey() : vtx0(0),vtx1(0) {} ;
      EdgeKey(int v0,int v1) : vtx0(v0),vtx1(v1) {} ;
      int Vtx0() { return(vtx0) ; } ;
      int Vtx1() { return(vtx1) ; } ;
    };

    private:

    int cur ;
    CArbArray<CArbMshTopo2D::ElemLoop> loops ;
    const CArbMshTopo2D *owner ;
} ;

inline int operator == (const CArbTopoLoopOnElemIterator::EdgeKey &edge0,
                        const CArbTopoLoopOnElemIterator::EdgeKey &edge1)
{
    return((edge0.vtx0 == edge1.vtx0) && (edge0.vtx1 == edge1.vtx1)) ;
}

inline int ArbHashIndex(const CArbTopoLoopOnElemIterator::EdgeKey &edge)
{
    return(edge.vtx0 + edge.vtx1) ;
} ;
 

#undef NULL_VTX

/*
TYPEDEFS

CArbHashTable<int,int> CArbElemSet - shorthand for a hash table used as 
            an element set 

CArbHashTableIterator<int,int> CArbElemSetIterator - shorthand for an 
            iterator for an element set 


CLASS CArbMshTopo2D

  ArbMshTopo2D objects store the topology associated with a 2D mesh and 
  support queries on the mesh. 


PUBLIC INTERFACE

  Public Data Structures:

    struct ArbVtxEntry

      This structure stores information about one adjacent vertex. The 
      object maintains a linked list of these structures for all 
      adjacent verticies for all vertices in a mesh. 

      Member Variables:

        int vtx_num - adjacent vertex id 

        int ccw_elem - element that is counter clockwise from 
            the edge from this vertext to the adjacent 

        struct ArbVtxEntry *next - linked list pointer 


    struct ArbEdge

      This structure stores information about an edge. 

      Member Variables:

        int nd0 - first node 

        int nd1 - second node 

        int elem0 - adjacent element 

        int elem1 - adjacent element 


  Public Member Functions:

    CArbMshTopo2D - constructor 

      CArbMshTopo2D()

      Description: This is the ArbMshTopo2D constructor. 


    CArbMshTopo2D - destructor 

      ~CArbMshTopo2D()

      Description: THis is a destructor for ArbMshTopo2D objects. 


    InsertElement - insert an element 

      void InsertElement(
              const int elem_id,
              const int num_nodes,
              const int *constnodes)

        elem_id   - (in)  element id 
        num_nodes - (in)  number of nodes 
        nodes     - (in)  node id's 

      Description: This function inserts an element into the current 
          mesh. 


    InsertTriangle - insert a triangular element 

      void InsertTriangle(
              const int elem_id,
              const int nd0,
              const int nd1,
              const int nd2)

        elem_id - (in)  element id 
        nd0     - (in)  first node 
        nd1     - (in)  second node 
        nd2     - (in)  third node 

      Description: This function inserts a triangular element into 
          the current mesh. 


    InsertCollapsedElement - insert a collapsed element 

      void InsertCollapsedElement(
              const int elem_id,
              const int num_nodes,
              const int *constnodes)

        elem_id   - (in)  element id 
        num_nodes - (in)  number of nodes 
        nodes     - (in)  node id's 

      Description: This function inserts an element for cases where 
          element id's may be repeated. 


    DeleteElement - delete an element 

      void DeleteElement(
              const int num_nodes,
              const int *constnodes)

        num_nodes - (in)  number of nodes 
        nodes     - (in)  node id's 

      Description: This function deletes an element from the mesh 
          when the node id's are known. 


    DeleteElement - delete an element 

      void DeleteElement(const int elem_id)

        elem_id - (in)  element id 

      Description: This function deletes an element from the mesh 
          when the element id is known. 


    DeleteTriangle - delete a triangle 

      void DeleteTriangle(
              const int nd0,
              const int nd1,
              const int nd2)

        nd0 - (in)  first node 
        nd1 - (in)  second node 
        nd2 - (in)  third node 

      Description: Delete a triangular element from a mesh when the 
          node id's are known. 


    InsertAngle - insert an angle 

      void InsertAngle(
              const int elem,
              const int vtx,
              const int prev,
              const int next)

        elem - (in)  element id 
        vtx  - (in)  center vertex 
        prev - (in)  previous vertex 
        next - (in)  next vertex 

      Description: Inserts three verticies and the included element 
          into a mesh if they do not exist in the mesh already. The 
          element will exist in the region moving ccw from prev to 
          next. 


    DeleteAngle - delete an angle 

      int DeleteAngle(
              const int vtx,
              const int prev,
              const int next)

        vtx  - (in)  center vertex 
        prev - (in)  previous vertex 
        next - (in)  next vertex 

      Description: Deletes an angle from the mesh. 

      Return Value: returns the formerly included element 


    DeleteElementRef - delete an element reference 

      void DeleteElementRef(const int elem_id)

        elem_id - (in)  element id 

      Description: Deletes an entry in the element table 


    Splice - splice together two element boundaries 

      bool Splice(
              const int     elem,
              const int     start,
              const int     finish,
              CArbMshTopo2D *other,
              const int     other_elem,
              const int     other_start,
              const int     other_finish)

        elem         - (in)  this elemnent id 
        start        - (in)  this start node 
        finish       - (in)  this finish node 
        other        - (in)  other mesh topology 
        other_elem   - (in)  other element id 
        other_start  - (in)  start node on other elem 
        other_finish - (in)  finish node on other elem 

      Description: This function splices the boundary for an element 
          stored in another mesh to the boundary of an element stored 
          in this mesh. This is done as follows: 

               +                                +
               |   other                        |
               |   start                        |
        start  +     +----+               start +----+
               |     |    | other   ==>              |
               |     |    | element                  |
        finish +     +----+              other  +----+
               |   other                 finish |
               |   finish                       |
               +                                +
        element boundary
            

      Return Value: returns true if the splice was succesful 


    DeleteNode - delete a node 

      void DeleteNode(int node)

        node - (in)  node id 

      Description: Deletes a node from the mesh. 


    HasVtx - vertex query 

      bool HasVtx(const int vtx) const

        vtx - (in)  vertex id 

      Description: Checks to see if a vertex id is in the mesh. 

      Return Value: returns true if the vetex (node) is in the mesh 


    OppositeEdge - find the edge opposite a node on an element 

      ArbEdge OppositeEdge(
              const int elem_id,
              const int node_id) const

        elem_id - (in)  element id 
        node_id - (in)  node id 

      Description: Given an element and a node on the element, this 
          function returns the edge on the element opposite the node. 
          This function assumes that that the element is a triangle. 

      Return Value: returns the opposite edge 


    OppositeNode - find the node opposite an edge on an element 

      int OppositeNode(
              const int     elem_id,
              const ArbEdge *edge) const

        elem_id - (in)  element id 
        edge    - (in)  edge description 

      Description: Given an element and an edge on the element, this 
          function returns the node on the element opposite the edge. 
          This function assumes that that the element is a triangle. 

      Return Value: the opposite node 


    ElemHasNode - check for a node on an element 

      bool ElemHasNode(
              const int elem_id,
              const int node_id) const

        elem_id - (in)  element id 
        node_id - (in)  node id 

      Description: This function checks to see if a node is on an 
          element. 

      Return Value: returns true if the node is on the element 


    NextCWEdge - get the next edge cw on an element 

      ArbEdge NextCWEdge(
              const int     elem_id,
              const ArbEdge *edge) const

        elem_id - (in)  element id 
        edge    - (in)  current edge 

      Description: Give an element and an edge, this function returns 
          the next edge moving in a clockwise direction. 

      Return Value: the next edge in a cw direction 


    NextCCWEdge - get the next edge cw on an element 

      ArbEdge NextCCWEdge(
              const int     elem_id,
              const ArbEdge *edge) const

        elem_id - (in)  element id 
        edge    - (in)  current edge 

      Description: Give an element and an edge, this function returns 
          the next edge moving in a counter clockwise direction. 

      Return Value: the next edge in a ccw direction 


    GetCWElem - get the cw element 

      int GetCWElem(
              const int node_0,
              const int node_1) const

        node_0 - (in)  first node 
        node_1 - (in)  second node 

      Description: This function finds the element that is clockwise 
          to the edge from node 0 to node 1. 

      Return Value: the cw element id 


    GetCCWElem - get the ccw element 

      int GetCCWElem(
              const int node_0,
              const int node_1) const

        node_0 - (in)  first node 
        node_1 - (in)  second node 

      Description: This function finds the element that is counter 
          clockwise to the edge from node 0 to node 1. 

      Return Value: the ccw element id 


    GetCWBdryNode - get the next node cw on the boundary 

      int GetCWBdryNode(
              const int node_0,
              const int node_1) const

        node_0 - (in)  first node 
        node_1 - (in)  second node 

      Description: This function find the node that is clockwise from 
          the given edge, staying on the boundary of the mesh. 

      Return Value: the cw node id 


    GetCCWBdryNode - get the next node ccw on the boundary 

      int GetCCWBdryNode(
              const int node_0,
              const int node_1) const

        node_0 - (in)  first node 
        node_1 - (in)  second node 

      Description: This function find the node that is counter 
          clockwise from the given edge, staying on the boundary of 
          the mesh. 

      Return Value: the ccw node id 


    GetCWNode - get the next node cw on this element 

      int GetCWNode(
              const int node_0,
              const int node_1) const

        node_0 - (in)  first node 
        node_1 - (in)  second node 

      Description: This function find the node that is clockwise from 
          the given edge, staying on the same element. 

      Return Value: the cw node id 


    GetCCWNode - get the next node ccw on this element 

      int GetCCWNode(
              const int node_0,
              const int node_1) const

        node_0 - (in)  first node 
        node_1 - (in)  second node 

      Description: This function find the node that is counter 
          clockwise from the given edge, staying on the same element. 

      Return Value: the ccw node id 


    NumAdjElems - find the number of adjacent element 

      int NumAdjElems(const int id) const

        id - (in)  node id 

      Description: Find the number of elements adjacent to a given 
          node. 

      Return Value: number of adjacent elements 


    NumConsecutiveAdjElems - find the number of consecutive adjacent 
                             elements 

      int NumConsecutiveAdjElems(const int id) const

        id - (in)  node id 

      Description: Find the size of the largest sequence of elements 
          that are adjacent to each other and the given node. 

      Return Value: greatest number of adjacent elements 


    BoundaryNode - boundary node query 

      bool BoundaryNode(const int id) const

        id - (in)  node id 

      Description: Checks to see if a node is on the boundary. 

      Return Value: True if the node is on the boundary 


    EdgeList - returns all edges in the mesh 

      CArbQueue <ArbEdge>*EdgeList()

      Description: This function returns all the edges in the mesh. 

      Return Value: a queue object containing all the edges in the 
          mesh 


    NumElemNodes - number of element nodes 

      int NumElemNodes(const int elem) const

        elem - (in)  element id 

      Description: Finds the number of nodes on a given element. 

      Return Value: number of nodes on an element 


    GetElemNodes - get an element's nodes 

      void GetElemNodes(
              const int elem,
              int       *num_nodes,
              int       *nodes) const

        elem      - (in)  element id 
        num_nodes - (out) number of nodes 
        nodes     - (out) node id's 

      Description: This function returns the id's of the nodes on a 
          given element 


    GetElemNodes - get an element's nodes 

      void GetElemNodes(
              const int elem,
              const int first_node_id,
              int       *num_nodes,
              int       *nodes) const

        elem          - (in)  element id 
        first_node_id - (in)  first node id 
        num_nodes     - (out) number of nodes 
        nodes         - (out) node id's 

      Description: This function returns the id's of the nodes on a 
          given element ordered so that the given node id is first. 


    NumElements - number of elements int the mesh 

      int NumElements() const

      Description: This function returns the number of elements in 
          the mesh. 

      Return Value: number of elements 


    GetElemList - get a list of element id's 

      int *GetElemList() const

      Description: This function returns a list of the id's for all 
          elements in the mesh. Ownership of the memory passes to the 
          client and must be freed with a call to delete []. 

      Return Value: a list (array) of the id's of all elements in the 
          mesh 


    GetElemSet - get a set of the element id's 

      CArbElemSet *GetElemSet() const

      Description: This function returns the id's of all the elements 
          in the mesh stored in a set object. Ownership of the set 
          object passes to the client who must make a call to delete. 

      Return Value: a set object containing the id's of the elements 


    NumNodes - finds the number of nodes in the mesh 

      int NumNodes() const

      Description: This function returns the number of nodes in the 
          mesh. 

      Return Value: number of nodes in the mesh 


    GetNodeList - get a list of node id's 

      int *GetNodeList() const

      Description: This function returns a list of the id's for all 
          nodes in the mesh. Ownership of the memory passes to the 
          client and must be freed with a call to delete []. 

      Return Value: a list (array) of the id's of all nodes in the 
          mesh 


PRIVATE INTERFACE

  Private Member Functions:

    FirstVtx - get the first adjacent vertex 

      ArbVtxEntry *FirstVtx(const int vtx) const

        vtx - (in)  input vertex 

      Description: This function returns a pointer to the ArbVtxEntry 
          data structure for the first vertex in the linked list of 
          vertices adjacent to the given vertex. 

      Return Value: a pointer to the first adjacent vertex 


    NextVtx - get the next adjacent vertex 

      ArbVtxEntry *NextVtx(const ArbVtxEntry *ptr) const

        ptr - (in)  current vertex 

      Description: This function returns a pointer to the ArbVtxEntry 
          data structure for the next vertex in the linked list of 
          vertices adjacent to the given vertex. 

      Return Value: a pointer to the next adjacent vertex 


    EolVtx - returns an end-of-list value 

      ArbVtxEntry *EolVtx()

      Description: This function returns a value that indicates the 
          end-of-list for the linked list of adjacent verticies 
          (zero, in fact). 

      Return Value: end-of-list value 


    GetAdjVtxAlongElem - get the next vertex on an element 

      int GetAdjVtxAlongElem(
              int node,
              int elem_id) const

        node    - (in)  node id 
        elem_id - (in)  element id 

      Description: Give an element and a vertex on the element, this 
          function returns the next vertex on the element moving in a 
          ccw direction. 

      Return Value: the id of the next ccw vertex on an element 


  Private Member Variables:

    CArbHashTable<int,ArbVtxEntry*>* vtx_table - hashes node id's 
            to adjacent vertex lists 

    CArbHashTable<int,int>* elem_table - hashes element id's to one 
            node on the element 

    ArbVtxEntry *MaxAddress - always greater than the maximum 
            adjacent vertex address 

    int MaxVtx - maximum vertex id 


---------------------------------------------------------------------

CLASS CArbTopoAdjVtxIterator

  This is an iterator object that allows a client to sequence through 
  the vertices adjacent to a given vertex in ccw order. 


PUBLIC INTERFACE

  Public Member Functions:

    CArbTopoAdjVtxIterator - constructor 

      CArbTopoAdjVtxIterator(
              const CArbMshTopo2D *aTopo,
              const int           ivtx)

        aTopo - (in)  an ArbMshTopo2D object 
        ivtx  - (in)  input vertex 

      Description: This is a constructor for a ArbTopoAdjVtxIterator 
          object. 


    First - point to the first adjacent vertex 

      void First()

      Description: This function sets the iterator to point to the 
          first adjacent vertex. 


    Next - point to the next adjacent vertex 

      void Next()

      Description: This function sets the iterator to point to the 
          next adjacent vertex. 


    More - checks for more adjacent verticies 

      bool More()

      Description: This function returns true if there are more 
          adjacent vertices, false otherwise. 


    NewVtx - resets the interator for a new vertex 

      void NewVtx(int ivtx)

        ivtx - (in)  new vertex 

      Description: This function resets an iterator so that it points 
          to the first adjacent vertex for the input vertex. 


    AdjVtx - gets the adjacent vertex id 

      int AdjVtx()

      Description: This function returns the id of the adjacent 
          vertex currently being pointed to by the iterator. 

      Return Value: the adjacent vertex id 


    CcwElem - gets the ccw element's id 

      int CcwElem()

      Description: This function returns the id of the element that 
          is ccw to the edge between the vertex and the adjacent 
          vertex currently being pointed to by the iterator. 

      Return Value: the ccw element id 


    operator ++ - increment operator 

      void operator ++ ()

      Description: This is the increment operator for the iterator. 
          It is the same as iter.Next(). 


    operator ++ - increment operator 

      void operator ++ (int i)

        i - (in)  ignored parameter 

      Description: This is the postfix increment operator for the 
          iterator. It is the same as iter.Next(). 


PRIVATE INTERFACE

  Private Member Variables:

    const CArbMshTopo2D *topo - ArbMshTopo2D object being 
            referenced 

    CArbMshTopo2D::ArbVtxEntry *ptr - pointer to the current 
            adjacent vertex structure 

    int vtx - current vertex pointer 


---------------------------------------------------------------------

CLASS CArbTopoAdjVtxCyclicIterator

  This is an iterator object that allows a client to sequence through 
  the vertices adjacent to a given vertex in either direction. 
  Increments move counter clockwise, decrements move clockwise. Ccw 
  motion is more efficient 


PUBLIC INTERFACE

  Public Member Functions:

    CArbTopoAdjVtxCyclicIterator - constructor 

      CArbTopoAdjVtxCyclicIterator(
              const CArbMshTopo2D *aTopo,
              const int           ivtx)

        aTopo - (in)  an ArbMshTopo2D object 
        ivtx  - (in)  input vertex 

      Description: This is a constructor for a 
          ArbTopoAdjVtxCyclicIterator object. 


    First - point to the first adjacent vertex 

      void First()

      Description: This function sets the iterator to point to the 
          first adjacent vertex. 


    Next - point to the next adjacent vertex 

      void Next()

      Description: This function sets the iterator to point to the 
          next adjacent vertex. 


    Prev - point to the previous adjacent vertex 

      void Prev()

      Description: This function sets the iterator to point to the 
          previous adjacent vertex. 


    IsValid - checks for a valid adjacent vertex 

      bool IsValid()

      Description: This function returns true if the iterator 
          references a valid adjacent vertex. 

      Return Value: true if the iterator reference is valid 


    NewVtx - resets the interator for a new vertex 

      void NewVtx(int ivtx)

        ivtx - (in)  new vertex 

      Description: This function resets an iterator so that it points 
          to the first adjacent vertex for the input vertex. 


    AdjVtx - gets the adjacent vertex id 

      int AdjVtx()

      Description: This function returns the id of the adjacent 
          vertex currently being pointed to by the iterator. 

      Return Value: the adjacent vertex id 


    CcwElem - gets the ccw element's id 

      int CcwElem()

      Description: This function returns the id of the element that 
          is ccw to the edge between the vertex and the adjacent 
          vertex currently being pointed to by the iterator. 

      Return Value: the ccw element id 


    CcwElemRef - gets a reference to the ccw element id 

      int &CcwElemRef()

      Description: This function returns a reference to the id of the 
          element that is ccw to the edge between the vertex and the 
          adjacent vertex currently being pointed to by the iterator. 

      Return Value: reference to the ccw element id 


    operator ++ - increment operator 

      void operator ++ ()

      Description: This is the increment operator for the iterator. 
          It is the same as iter.Next(). 


    operator ++ - increment operator 

      void operator ++ (int i)

        i - (in)  ignored parameter 

      Description: This is the postfix increment operator for the 
          iterator. It is the same as iter.Next(). 


    operator -- - decrement operator 

      void operator -- ()

      Description: This is the decrement operator for the iterator. 
          It is the same as iter.Prev(). 


    operator -- - decrement operator 

      void operator -- (int i)

        i - (in)  ignored parameter 

      Description: This is the postfix decrement operator for the 
          iterator. It is the same as iter.Prev(). 


PRIVATE INTERFACE

  Private Member Variables:

        const CArbMshTopo2D *topo - ArbMshTopo2D object being 
                referenced 

        CArbMshTopo2D::ArbVtxEntry *ptr - pointer to the current 
                adjacent vertex structure 

        int vtx - current vertex pointer 


---------------------------------------------------------------------

CLASS CArbTopoVtxOnElemIterator

  This is an iterator object that allows a client to sequence through 
  the vertices adjacent to an element in ccw order. 


PUBLIC INTERFACE

  Public Member Functions:

    CArbTopoVtxOnElemIterator - constructor 

      CArbTopoVtxOnElemIterator(
              const CArbMshTopo2D *aTopo,
              const int           elem)

        aTopo - (in)  an ArbMshTopo2D object 
        elem  - (in)  input element 

      Description: This is a constructor for a 
          ArbTopoAdjVtxOnElemIterator object. 


    CArbTopoVtxOnElemIterator - constructor 

      CArbTopoVtxOnElemIterator(
              const CArbMshTopo2D *aTopo,
              const int           elem,
              const int           node)

        aTopo - (in)  an ArbMshTopo2D object 
        elem  - (in)  input element 
        node  - (in)  first vertex 

      Description: This is a constructor for a 
          ArbTopoAdjVtxOnElemIterator object. 


    First - point to the first vertex 

      void First()

      Description: This function sets the iterator to point to the 
          first vertex on the element. 


    Next - point to the next vertex 

      void Next()

      Description: This function sets the iterator to point to the 
          next vertex on the element. 


    More - checks for more verticies 

      bool More()

      Description: This function returns true if there are more 
          vertices, false otherwise. 


    Current - returns the current vertex 

      int Current()

      Description: This function returns the id of the vertex that 
          the iterator currently references. 

      Return Value: the vertex id 


    operator ++ - increment operator 

      void operator ++ ()

      Description: This is the increment operator for the iterator. 
          It is the same as iter.Next(). 


    operator ++ - increment operator 

      void operator ++ (int i)

        i - (in)  ignored parameter 

      Description: This is the postfix increment operator for the 
          iterator. It is the same as iter.Next(). 


    operator * - returns the current vertex 

      int operator * ()

      Description: This operator returns the id of the vertex that 
          the iterator currently references. 

      Return Value: the referenced id 


PRIVATE INTERFACE

  Private Member Variables:

    int cur - the current vertex id 

    int first - the first vertex id 

    int prev - the previous vertex id 

    int elem - the element currently being referenced 

    const CArbMshTopo2D *topo - the ArbMshTopo2D object being 
            referenced 


---------------------------------------------------------------------

CLASS CArbTopoVtxOnElemCyclicIterator

  This is an iterator object that allows a client to sequence through 
  the vertices on an element in either direction. Increments move 
  counter clockwise, decrements move clockwise. Ccw motion is more 
  efficient. 


PUBLIC INTERFACE

  Public Member Functions:

    CArbTopoVtxOnElemCyclicIterator - constructor 

      CArbTopoVtxOnElemCyclicIterator(
              const CArbMshTopo2D *aTopo,
              const int           elem)

        aTopo - (in)  an ArbMshTopo2D object 
        elem  - (in)  input element 

      Description: This is a constructor for a 
          ArbTopoVtxOnElemCyclicIterator object. 


    CArbTopoVtxOnElemCyclicIterator - constructor 

      CArbTopoVtxOnElemCyclicIterator(
              const CArbMshTopo2D *aTopo,
              const int           elem,
              const int           node)

        aTopo - (in)  an ArbMshTopo2D object 
        elem  - (in)  input element 
        node  - (in)  first vertex 

      Description: This is a constructor for a 
          ArbTopoVtxOnElemCyclicIterator object. 


    First - point to the first vertex 

      void First()

      Description: This function sets the iterator to point to the 
          first vertex on the element. 


    Next - point to the next vertex 

      void Next()

      Description: This function sets the iterator to point to the 
          next vertex on the element. 


    Prev - point to the previous vertex 

      void Prev()

      Description: This function sets the iterator to point to the 
          previous vertex on the element. 


    IsValid - checks for a valid vertex 

      bool IsValid()

      Description: This function returns true if the iterator 
          references a valid vertex. 

      Return Value: true if the iterator reference is valid 


    Current - returns the current vertex 

      int Current()

      Description: This function returns the id of the vertex that 
          the iterator currently references. 

      Return Value: the vertex id 


    operator ++ - increment operator 

      void operator ++ ()

      Description: This is the increment operator for the iterator. 
          It is the same as iter.Next(). 


    operator ++ - increment operator 

      void operator ++ (int i)

        i - (in)  ignored parameter 

      Description: This is the postfix increment operator for the 
          iterator. It is the same as iter.Next(). 


    operator -- - decrement operator 

      void operator -- ()

      Description: This is the decrement operator for the iterator. 
          It is the same as iter.Prev(). 


    operator -- - decrement operator 

      void operator -- (int i)

        i - (in)  ignored parameter 

      Description: This is the postfix decrement operator for the 
          iterator. It is the same as iter.Prev(). 


    operator * - returns the current vertex 

      int operator * ()

      Description: This operator returns the id of the vertex that 
          the iterator currently references. 

      Return Value: the referenced id 


PRIVATE INTERFACE

  Private Member Variables:

    int cur - the current vertex id 

    int first - the first vertex id 

    int prev - the previous vertex id 

    int elem - the element currently being referenced 

    const CArbMshTopo2D *topo - the ArbMshTopo2D object being 
            referenced 


---------------------------------------------------------------------

CLASS CArbTopoVtxOnBdryIterator

  This is an iterator object that allows a client to sequence through 
  the vertices on the boundary of the mesh. 


PUBLIC INTERFACE

  Public Member Functions:

    CArbTopoVtxOnBdryIterator - constructor 

      CArbTopoVtxOnBdryIterator(
              const CArbMshTopo2D *aTopo,
              const int           node)

        aTopo - (in)  an ArbMshTopo2D object 
        node  - (in)  starting vertex 

      Description: This is a constructor for a 
          ArbTopoVtxOnBdryIterator object. 


    First - point to the first vertex 

      void First()

      Description: This function sets the iterator to point to the 
          first vertex on the boundary. 


    Next - point to the next vertex 

      void Next()

      Description: This function sets the iterator to point to the 
          next vertex on the boundary. 


    More - checks for more verticies 

      bool More()

      Description: This function returns true if there are more 
          vertices, false otherwise. 

      Return Value: true means more verticies 


    Current - returns the current vertex 

      int Current()

      Description: This function returns the id of the vertex that 
          the iterator currently references. 

      Return Value: the vertex id 


    operator ++ - increment operator 

      void operator ++ ()

      Description: This is the increment operator for the iterator. 
          It is the same as iter.Next(). 


    operator ++ - increment operator 

      void operator ++ (int i)

        i - (in)  ignored parameter 

      Description: This is the postfix increment operator for the 
          iterator. It is the same as iter.Next(). 


    operator * - returns the current vertex 

      int operator * ()

      Description: This operator returns the id of the vertex that 
          the iterator currently references. 

      Return Value: the referenced id 


PRIVATE INTERFACE

  Private Member Variables:

    int cur - the current vertex id 

    int first - the first vertex id 

    int prev - the previous vertex id 

    int node - the element currently being referenced 

    const CArbMshTopo2D *owner - the ArbMshTopo2D object being 
            referenced 


---------------------------------------------------------------------

CLASS CArbTopoEdgeOnElemIterator

  This is an iterator object that allows a client to sequence through 
  the edges on an element in counter clockwise order. 


PUBLIC INTERFACE

  Public Member Functions:

    CArbTopoEdgeOnElemIterator - constructor 

      CArbTopoEdgeOnElemIterator(
              const CArbMshTopo2D *aTopo,
              const int           elem)

        aTopo - (in)  an ArbMshTopo2D object 
        elem  - (in)  input element 

      Description: This is a constructor for a 
          ArbTopoEdgeOnElemIterator object. 


    First - point to the first edge 

      void First()

      Description: This function sets the iterator to point to the 
          first edge on the element. 


    Next - point to the next edge 

      void Next()

      Description: This function sets the iterator to point to the 
          next edge on the element. 


    More - checks for more edges 

      bool More()

      Description: This function returns true if there are more 
          edges, false otherwise. 

      Return Value: true if there are more edges 


    Current0 - returns the first current vertex 

      int Current0()

      Description: This function returns the id of the first vertex 
          on the edge that the iterator currently references. 

      Return Value: the first vertex id 


    Current1 - returns the second current vertex 

      int Current1()

      Description: This function returns the id of the second vertex 
          on the edge that the iterator currently references. 

      Return Value: the second vertex id 


    operator ++ - increment operator 

      void operator ++ ()

      Description: This is the increment operator for the iterator. 
          It is the same as iter.Next(). 


    operator ++ - increment operator 

      void operator ++ (int i)

        i - (in)  ignored parameter 

      Description: This is the postfix increment operator for the 
          iterator. It is the same as iter.Next(). 


    operator [] - returns one of the current vertices 

      int operator [] (int i)

        i - (in)  vertex index, 0 or 1 

      Description: This operator returns the id of the first, [0], or 
          second, [1], vertex on the edge that the iterator currently 
          references. 

      Return Value: the vertex id 


PRIVATE INTERFACE

  Private Member Variables:

    int cur - the current vertex id 

    int adj - adjacent vertex id 

    int first - first vertex id 

    int prev - previous vertex id 

    int elem - element id 

    const CArbMshTopo2D *owner - the ArbMshTopo2D object being 
            referenced 


---------------------------------------------------------------------

CLASS CArbTopoEdgeOnElemCyclicIterator

  This is an iterator object that allows a client to sequence through 
  the edges on an element in either direction. Increments move counter 
  clockwise, decrements move clockwise. Ccw motion is more efficient. 


PUBLIC INTERFACE

  Public Member Functions:

    CArbTopoEdgeOnElemCyclicIterator - constructor 

      CArbTopoEdgeOnElemCyclicIterator(
              const CArbMshTopo2D *aTopo,
              const int           elem)

        aTopo - (in)  an ArbMshTopo2D object 
        elem  - (in)  element 

      Description: This is a constructor for a 
          ArbTopoVtxOnBdryIterator object. 


    First - point to the first edge 

      void First()

      Description: This function sets the iterator to point to the 
          first edge on the element. 


    Next - point to the next edge 

      void Next()

      Description: This function sets the iterator to point to the 
          next edge on the element. 


    Prev - point to the previous edge 

      void Prev()

      Description: This function sets the iterator to point to the 
          previous edge on the element. 


    IsValid - checks for a valid edge 

      bool IsValid()

      Description: This function returns true if the iterator 
          currently points to a valid edge, false otherwise. 

      Return Value: true for a valid edge 


    Current0 - returns the first current vertex 

      int Current0()

      Description: This function returns the id of the first vertex 
          on the edge that the iterator currently references. 

      Return Value: the first vertex id 


    Current1 - returns the second current vertex 

      int Current1()

      Description: This function returns the id of the second vertex 
          on the edge that the iterator currently references. 

      Return Value: the second vertex id 


    operator ++ - increment operator 

      void operator ++ ()

      Description: This is the increment operator for the iterator. 
          It is the same as iter.Next(). 


    operator ++ - increment operator 

      void operator ++ (int i)

        i - (in)  ignored parameter 

      Description: This is the postfix increment operator for the 
          iterator. It is the same as iter.Next(). 


    operator -- - decrement operator 

      void operator -- ()

      Description: This is the deccrement operator for the iterator. 
          It is the same as iter.Prev(). 


    operator -- - decrement operator 

      void operator -- (int i)

        i - (in)  ignored parameter 

      Description: This is the postfix decrement operator for the 
          iterator. It is the same as iter.Next(). 


    operator [] - returns one of the current vertices 

      int operator [] (int i)

        i - (in)  vertex index, 0 or 1 

      Description: This operator returns the id of the first, [0], or 
          second, [1], vertex on the edge that the iterator currently 
          references. 

      Return Value: the vertex id 


PRIVATE INTERFACE

  Private Member Variables:

    int cur - the current vertex id 

    int adj - adjacent vertex id 

    int first - first vertex id 

    int prev - previous vertex id 

    int elem - element id 

    const CArbMshTopo2D *owner - the ArbMshTopo2D object being 
            referenced 


*/
#endif
