//
// CArbEdgeList Class header file
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
//   $Revision: 1.12 $  $Date: 2002/09/12 18:31:10 $  $Author: wash $
//

#ifndef ArbMshEdgeList_h
#define ArbMshEdgeList_h

#include "ArbMsh.hpp"
#include "ArbHeap.hpp"
#include "ArbHashTable.hpp"
#include "ArbSet.hpp"
#include "ArbMshTopo2D.hpp"
#include "ArbCoord2D.hpp"

class CArbMshEdgeList {
    public:
    // public data structures

        struct ArbQdEdge {
            double angle_0 ;
            double angle_1 ;
            int id_0 ;
            int id_1 ;
            int level ;
            short end_code ;
        } ;

        struct ArbIntQdEdgeDesc {
            double length ;
            double angle_0 ;
            double angle_1 ;
            int id_0 ;
            int id_1 ;
            short end_code ;
            short level ;
            short adj_level_0 ;
            short adj_level_1 ;
            bool visited ;
        } ;

        struct ArbNewQuadData {
            int qcase ;
            int bl, br, tr, tl ;
            int obl, obr, otr, otl ;
            int num_front_nodes ;
            int front_nodes[8] ;
            int edge_level ;
            double ratio ;
            double base_template_ratio ;
            double top_template_ratio ;
            bool valid_transition_split ;
        } ;

        struct ArbNewSeamData {
            int id0 ;
            int id1 ;
            int id2 ;
            int cw ;
            int ccw ;
            int min_level ;
            double ratio ;
            CArbCoord2D coord ;
            bool cw_side ;
            bool do_template ;
            bool template_cw ;
            bool merge_quads ;
        } ;

        typedef CArbHeap<ArbIntQdEdgeDesc>::EntryHandle EdgeHandle ;

    // constructors and destructors

        CArbMshEdgeList(CArbHashTable<int,ArbIntNode> *inode_table) ;

        ~CArbMshEdgeList() ;

    // insert and removal routines

        void InsertEdge(const int id_0,
                        const int id_1,
                        const int cw_id,
                        const int ccw_id) ;

        int InitializeCodes() ;

        int ClassifyQuad(const int id_0,
                          const int id_1,
                          const int id_2,
                          const int id_3,
                          ArbNewQuadData *qdata) ;

        bool DoesThisCrossBoundary(const CArbCoord2D &pt0,
                                   const CArbCoord2D &pt1) const ;

        void UpdateQuad(ArbNewQuadData *qdata) ;

        void UpdateTriangle(const int id_0,
                            const int id_1,
                            const int id_2) ;

        int UpdateGeom(ArbNewQuadData *qdata) ;

        void UpdateBdryGeom(const int vtx) ;

        int UpdateSeam(ArbNewSeamData *sdata) ;

        void UpdateTransSeam(const ArbNewSeamData *sdata,
                             const int mid_id,
                             const int far_id) ;

        void UpdateTransSplit(const ArbNewQuadData *qdata,
                             const int mid_id0,
                             const int mid_id1) ;

        void UpdateTempSplit(const ArbNewQuadData *qdata,
                             const int mid_id0,
                             const int mid_id1,
                             const bool base) ;

        void UpdateSeamGeom(const ArbNewSeamData *sdata) ;

        int GetCWNode(const int id_0,
                               const int id_1) const ;
        int GetCCWNode(const int id_0,
                                const int id_1) const ;
        double GetEdgeLength(const int id_0,
                             const int id_1) const ;
        int GetEdgeLevel(const int id_0,
                                  const int id_1) const ;

        bool ContainsEdge(const int id_0,
                          const int id_1) const ;

        bool ContainsNode(const int id) const ;

        ArbQdEdge *GetNextEdge() ;

        void PushToBack(const ArbQdEdge *edge) ;

        double GetEdgeLengthRatio() const ;

        bool StagnationSweep() ;

        double GetCharNodeLength(const int id) const ;

    private:

    // member variable

        CArbHashTable<int,ArbIntNode> *node_table ;
        CArbHashTable<ArbEdgeKey,EdgeHandle> *edge_table ;
        CArbHeap<ArbIntQdEdgeDesc> *order_heap ;
        CArbSet<double> *edge_lengths ;
        CArbHashTable<int,int> *front_nodes ;
        CArbMshTopo2D *bdry_topo ;

    // member functions

        int UpdateCase1(ArbNewQuadData *qdata) ;
        int UpdateCase2(ArbNewQuadData *qdata) ;
        int UpdateCase3(ArbNewQuadData *qdata) ;
        int UpdateCase4(ArbNewQuadData *qdata) ;
        int UpdateCase5(ArbNewQuadData *qdata) ;
        int UpdateCase6(ArbNewQuadData *qdata) ;

        int UpdateGeomHelp(const int num_front_nodes,
                            const int *front_nodes,
                            const int edge_level) ;

        bool UpdateGeomOne(const int node_id,int front_nodes[]) ;

        int FindAdjacent(const int vtx,
                          const int new_vtx,
                          int *prev_vtx,
                          int *next_vtx) ;

         double Angle(CArbCoord2D b,CArbCoord2D i,
                      CArbCoord2D j) const ;
         double Angle2Pi(CArbCoord2D b,CArbCoord2D i,
                         CArbCoord2D j) const ;
         double CrossProd(CArbCoord2D b,CArbCoord2D i,
                          CArbCoord2D j) const ;

    friend class CArbEdgeListNodeIterator ;

} ;

class CArbEdgeListNodeIterator {
    public:

        CArbEdgeListNodeIterator(const CArbMshEdgeList *edge_list)
            { Niter = new CArbHashTableIterator<int,int>(
                                  edge_list->front_nodes) ; } ;
        ~CArbEdgeListNodeIterator() { delete Niter ; } ;


        void First()             { Niter->First() ; } ;
        void Next()              { Niter->Next() ; } ;
        bool More()              { return(Niter->More()) ; } ;
        int  Current()           { return(Niter->Key()) ; } ;

        void operator ++ ()      { Next() ; } ;
        void operator ++ (int ) { Next() ; } ;
        int operator * ()        { return(Current()) ; } ;

    private:

        CArbHashTableIterator<int,int> *Niter ;
} ;

//inline int operator ==(const CArbMshEdgeList::ArbIntQdEdge &op1,
//                       const CArbMshEdgeList::ArbIntQdEdge &op2)
//{
//    return((op1.id_0 == op2.id_0) && (op1.id_1 == op2.id_1)) ;
//}

inline int operator ==(const CArbMshEdgeList::ArbIntQdEdgeDesc &op1,
                       const CArbMshEdgeList::ArbIntQdEdgeDesc &op2)
{
    return(op1.id_0 == op2.id_0) ;
}

/*
CLASS CArbMshEdgeList

  An edge list is an object used during advancing front mesh
  generation. It stores a list of edges that make up the current mesh
  front.

  When using this object there is an initialization phase and then a
  dynamic update phase. In the initialization phase, the client should
  make multiple calls to InsertEdge to define the initial front. Once
  all edges are added, InitializeCodes should be call to set up the
  remaining internal data structures. After this, the Update* member
  functions can be call to perform dynamic updates.

  A number of data structures associated with this object use a
  "end_code" variable to characterize the types of angles at the ends
  of an edge. The ends can be approximately 180 degrees or bigger,
  approximately 90 degress, or a very accute angle. This information is
  encoded in the 4 low-order bits of the of the end_code.
  Unfortunately, this is done in an inconsistent manner. The bits are
  defined as follows:


              +-------- 1 means cw angle is small
              | +------ 1 means ccw angle is small
              | | +---- 1 means cw angle is ~ 90 degrees
              | | | +-- 1 means ccw angle is ~ 90 degrees
              | | | |
              v v v v
       +-----+-+-+-+-+
       | ... | | | | |
       +-----+-+-+-+-+
              3 2 1 0


        Graphically the codes look something like this

                                                      *
       code 0                      code 1             |
              *-----*=====*-----*         *-----*=====*

                *                              *     *
       code 2   |                  code 3      |     |
                *=====*-----*                  *=====*

                           *                   *   *
       code 4,5             \      code 6,7    |    \
                 *-----*=====*                 *=====*

                   *                             *   *
       code 8,10  /                code 9,11    /    |
                 *=====*-----*                 *=====*

                          * *
       code 12,13,14,15  /   \
                        *=====*



PUBLIC INTERFACE

  Public Typedefs:

        CArbHeap<ArbIntQdEdgeDesc>::EntryHandle EdgeHandle - handle
                that references a edge in the edge heap.

  Public Data Structures:

    struct ArbQdEdge

      This structure is used to return information about edges in the
      list to the client.

      Member Variables:

        double angle_0 - angle the edge makes with its cw
            neighbor

        double angle_1 - angle the edge makes with its ccw
            neighbor

        int id_0 - id of the left (cw) node

        int id_1 - id of the right (ccw) node

        int level - approximate distance (number of edges) from
            the original boundary

        short end_code - code to classify end angles


    struct ArbIntQdEdgeDesc

      This structure is used to store information about one edge in the
      list.

      Member Variables:

        double length - length of the edge

        double angle_0 - angle between this and its cw neighbor

        double angle_1 - angle between this and its ccw
            neighbor

        int id_0 - id of the left (cw) end node

        int id_1 - id of the right (ccw) end node

        short end_code - code to classify end angles

        short level - approximate distance from the original
            boundary

        short adj_level_0 - level of the cw neighbor

        short adj_level_1 - level of the ccw neighbor


    struct ArbNewQuadData

      This structure is used to store information about a new
      quadrilateral that is to be added to a mesh.

      The obl, obr, otl, and otr variabls are for client's use. They
      are not used by the object.

      Member Variables:

        int qcase - new quad case (number of edges to add), see
            the internal documentation of ClassifyQuad for
            details

        int bl - bottom left node

        int br - bottom right node

        int tr - top right node

        int tl - top left node

        int obl - original bl node

        int obr - original br node

        int otr - original tr node

        int otl - original tl node

        int num_front_nodes - number of nodes on the
            front_nodes list

        int front_nodes[8] - list of nodes on the front the go
            from one be before to one after the nodes of the
            new element that make up part of the front

        int edge_level - number of edges from the original
            front

        double ratio - ratio of a quads base to one of its
            sides

        double base_template_ratio - ratio of the length of the
            base to the smaller of its sides

        double top_template_ratio - ratio of the length of the
            top to the smaller of its sides

        bool valid_transition_split - flag to indicate that a
            transition split would be valid


    struct ArbNewSeamData

      This data structure is used so that exchange information about a
      seam update. Topologically the seaming operation looks like the
      following:

           cw *
              |                    cw *
              |                       |
          id0 *               =>      |
              |                   id0 *-----*
              |                            ccw
          id1 *-----*-----*
                   id2   ccw

      Member Variables:

        int id0 - cw node in seam

        int id1 - center node in seam

        int id2 - ccw node in seam

        int cw - node cw from id0

        int ccw - node ccw from id2

        int min_level - smaller of the levels for the two
            seamed edges

        double ratio - ratio of the lengths of the edges being
            seamed

        CArbCoord2D coord - unused

        bool cw_side - flag that indicates if the cw edge
            (true) or ccw edge (false) will be retained

        bool do_template - flag to indicate a template seam

        bool template_cw - flag to indicate if the template is
            done on the cw side (true) or ccw side (false)


  Public Member Functions:

    CArbMshEdgeList - edge list constructor

      CArbMshEdgeList(CArbHashTable<int,ArbIntNode>* inode_table)

        inode_table - (in)  node hash table

      Description: This is a constructor for an edge list. As an
          argument it takes a pointer to a hash table that maps node
          id's to ArbIntNodes.


    CArbMshEdgeList - edge list destructor

      ~CArbMshEdgeList()

      Description: This is a desctructor for an edge list.


    InsertEdge - insert an edge into the list

      void InsertEdge(
              const int id_0,
              const int id_1,
              const int cw_id,
              const int ccw_id)

        id_0   - (in)  first node id
        id_1   - (in)  second node id
        cw_id  - (in)  id of node cw from id_0
        ccw_id - (in)  if of node ccw from id_1

      Description: This function inserts a new edge into the edge
          list.


    InitializeCodes - initialize the end codes for all edges

      void InitializeCodes()

      Description: This function goes through all edges currently in
          the list and computes the angles they make with neighbors
          and store this information in the heap that will be used to
          select an edge from which an element will be constructed.


    ClassifyQuad - find information for a new quadrilateral

      void ClassifyQuad(
              const int      id_0,
              const int      id_1,
              const int      id_2,
              const int      id_3,
              ArbNewQuadData *qdata)

        id_0  - (out) first node
        id_1  - (out) second node
        id_2  - (out) third node
        id_3  - (out) forth node
        qdata - (out) quadrilateral information

      Description: Given the id's for the four corner nodes for a new
          quadrilatral, this function fills the qdata argument with
          information about the new quadrilateral.


    UpdateQuad - update the list to add a quad element

      void UpdateQuad(ArbNewQuadData *qdata)

        qdata - (in)  quad description

      Description: This function is called to update the edge list to
          add a new quadrilateral to the mesh.


    UpdateTriangle - update the list to add a quad element

      void UpdateTriangle(
              const int id_0,
              const int id_1,
              const int id_2)

        id_0 - (in)  first id
        id_1 - (in)  second id
        id_2 - (in)  third id

      Description: This function is called to update the edge list to
          add a new triangular element to the mesh. This is done only
          to close off a triangular void.


    UpdateGeom - update stored geometry parameters

      void UpdateGeom(ArbNewQuadData *qdata)

        qdata - (in)  quadrilateral description

      Description: This function updates the stored mesh front edge
          lengths and the angles associated with a new quad. This
          function is called after nodal smoothing.


    UpdateBdryGeom - update stored geometry parameters

      void UpdateBdryGeom(const int vtx)

        vtx - (in)  vertex id

      Description: This function updates the stored mesh front
          lengths and angle associated with one vertex. This function
          is called after nodal smoothing.


    UpdateSeam - update the front for a seam operation

      int UpdateSeam(ArbNewSeamData *sdata)

        sdata - (in)  seam description data

      Description: This function updates the edge list for a seam
          operation.

      Return Value: The function returns the id of the node that will
          be kept after the seam operation.


    UpdateTransSeam - update the front for a transition seam

      void UpdateTransSeam(
              const ArbNewSeamData *sdata,
              const int            mid_id,
              const int            far_id)

        sdata  - (in)  seam description data
        mid_id - (in)  id of the mid point for the transition
        far_id - (in)  id of the far point for the transition

      Description: This function updates the edge list for a
          transition seam operation.


    UpdateTransSplit - update the front for a transition split

      void UpdateTransSplit(
              const ArbNewQuadData *qdata,
              const int            mid_id0,
              const int            mid_id1)

        qdata   - (in)  quadrilateral description
        mid_id0 - (in)  id of the first mid point vertex
        mid_id1 - (in)  id of the second mid point vertex

      Description: This function updates the edge list for a
          transition split operation.


    UpdateTempSplit - update the front for a template split

      void UpdateTempSplit(
              const ArbNewQuadData *qdata,
              const int            mid_id0,
              const int            mid_id1,
              const bool           base)

        qdata   - (in)  quadrilateral data
        mid_id0 - (in)  first mid point id
        mid_id1 - (in)  second mid point id
        base    - (in)  true means split the base edge

      Description: This function updates the edge list for a template
          split operation.


    UpdateSeamGeom - update stored geometry parameters

      void UpdateSeamGeom(const ArbNewSeamData *sdata)

        sdata - (in)  seam data

      Description: This function updates the stored mesh front
          lengths and angle. This function is called after nodal
          smoothing.


    GetCWNode - find the clockwise node

      int GetCWNode(
              const int id_0,
              const int id_1) const

        id_0 - (in)  first edge node
        id_1 - (in)  second edge node

      Description: This function returns the id of the next node on
          the crack front that is clockwise from the edge defined by
          the arguments.

      Return Value: The id of the node clockwise from edge.


    GetCCWNode - find counter-clockwise node

      int GetCCWNode(
              const int id_0,
              const int id_1) const

        id_0 - (in)  first edge node
        id_1 - (in)  second edge node

      Description: This function returns the id of the next node on
          the crack front that is counter- clockwise from the edge
          defined by the arguments.

      Return Value: The id of the node counter-clockwise from edge.


    GetEdgeLength - return the length of an edge

      double GetEdgeLength(
              const int id_0,
              const int id_1) const

        id_0 - (in)  first edge node
        id_1 - (in)  second edge node

      Description: This function returns the length of the edge
          defined by the arguments.

      Return Value: The length of the edge.


    GetEdgeLevel - return the level of an edge

      int GetEdgeLevel(
              const int id_0,
              const int id_1) const

        id_0 - (in)  first edge node
        id_1 - (in)  second edge node

      Description: This function returns the level (approximate
          number of steps required to get to the original boundary)
          of the edge defined by the arguments.

      Return Value: The level of the edge.


    ContainsEdge - check to see if an edge is in the list

      bool ContainsEdge(
              const int id_0,
              const int id_1) const

        id_0 - (in)  first edge node
        id_1 - (in)  second edge node

      Description: This function checks to see if the edge define by
          the arguments is in the list.

      Return Value: True if the edge is in the list, false otherwise.


    ContainsNode - checks to see if a node is in the list

      bool ContainsNode(const int id) const

        id - (in)  node id

      Description: This function checks to see if the node specified
          by the argument is contained in in the edge list.

      Return Value: True if the node is found, false, otherwise


    GetNextEdge - get next edge to use to advance the front

      ArbQdEdge *GetNextEdge()

      Description: This function returns information about the edge
          that is ranked highest in the edge list's priority queue
          for use in advancing the mesh front.

      Return Value: An ArbQdEdge structure containing information the
          highest ranked edge.


    PushToBack - put an edge back on the edge list

      void PushToBack(const ArbQdEdge *edge)

        edge - (in)  edge to place back on the list

      Description: This function places an edge back on the edge list
          priority queue, and flags it so that it will be ranked
          lowest in the queue.


    GetEdgeLengthRatio - get the edge length ratio

      double GetEdgeLengthRatio() const

      Description: This function returns the ratio of the longest to
          shortest edges currently in the edge list.

      Return Value: The ratio of the longest to shortest edge.


PRIVATE INTERFACE

  Private Member Functions:

    UpdateCase1 - update the edge list for a case 1 quad

      int UpdateCase1(ArbNewQuadData *qdata)

        qdata - (in)  quadrilateral data

      Description: This function updates the edge list to add a case
          1 quadrilateral. A case 1 quad has one edge on the front
          and adds three new edges.

      Return Value: The edge level of the new edges.


    UpdateCase2 - update the edge list for a case 2 quad

      int UpdateCase2(ArbNewQuadData *qdata)

        qdata - (in)  quadrilateral data

      Description: This function updates the edge list to add a case
          2 quadrilateral. A case 2 quad has two edge on the front
          and adds two new edges. The two edges on the front are
          adjacent.

      Return Value: The edge level of the new edges.


    UpdateCase3 - update the edge list for a case 3 quad

      int UpdateCase3(ArbNewQuadData *qdata)

        qdata - (in)  quadrilateral data

      Description: This function updates the edge list to add a case
          3 quadrilateral. A case 3 quad has three edge on the front
          and adds one new edges.

      Return Value: The edge level of the new edges.


    UpdateCase4 - update the edge list for a case 4 quad

      int UpdateCase4(ArbNewQuadData *qdata)

        qdata - (in)  quadrilateral data

      Description: This function updates the edge list to add a case
          4 quadrilateral. A case 4 quad has all for edge on the
          front and adds no new edges.

      Return Value: The edge level of the new edges.


    UpdateCase5 - update the edge list for a case 4 quad

      int UpdateCase5(ArbNewQuadData *qdata)

        qdata - (in)  quadrilateral data

      Description: This function updates the edge list to add a case
          2 quadrilateral. A case 2 quad has two edge on the front
          and adds two new edges. The two edges on the front are not
          adjacent.

      Return Value: The edge level of the new edges.


    UpdateGeomHelp - support for UpdateGeom

      void UpdateGeomHelp(
              const int num_front_nodes,
              const int *front_nodes,
              const int edge_level)

        num_front_nodes - (in)  number of front nodes associated with
                                this quad
        front_nodes     - (in)  list of front nodes
        edge_level      - (in)  edge level to assign to updated edges

      Description: This function is a support routine for the
          UpdateGeom routine. It updates an edges end codes and edge
          lengths.


    FindAdjacent - find previous and next adjacent verticies

      void FindAdjacent(
              const int vtx,
              const int new_vtx,
              int       *prev_vtx,
              int       *next_vtx)

        vtx      - (in)  existing vertex
        new_vtx  - (in)  new vertex
        prev_vtx - (out) previous existing vertex
        next_vtx - (out) next existing vertex

      Description: This function looks around a vertex and determines
          which currently adjacent verticies would come before and
          after a new vertex.


    Angle - compute the included angle

      double Angle(
              CArbCoord2D b,
              CArbCoord2D i,
              CArbCoord2D j) const

        b - (in)  center point
        i - (in)  first point
        j - (in)  second point

      Description: This function computes the included angle between
          three points.

      Return Value: Included angle, -PI < angle < PI


    Angle2Pi - compute the included angle

      double Angle2Pi(
              CArbCoord2D b,
              CArbCoord2D i,
              CArbCoord2D j) const

        b - (in)  center point
        i - (in)  first point
        j - (in)  second point

      Description: This function computes the included angle between
          three points.

      Return Value: Included angle, 0 < angle < 2*PI


    CrossProd - cross product of two vectors

      double CrossProd(
              CArbCoord2D b,
              CArbCoord2D i,
              CArbCoord2D j) const

        b - (in)  common start point
        i - (in)  first end point
        j - (in)  second end point

      Description: This function computes the cross product of two
          vectors. The vectors are defined by one common start point
          and two separate end points.

      Return Value: the cross product


  Private Member Variables:

    CArbHashTable<int,ArbIntNode>* node_table - pointer to the node
            hash table

    CArbHashTable<ArbEdgeKey,EdgeHandle>* edge_table - pointer to
            the edge hash table

    CArbHeap<ArbIntQdEdgeDesc>* order_heap - priority queue used to
            order the edges in the list

    CArbSet<double>* edge_lengths - ordered list of the edge
            lengths currently in the edge list

    CArbHashTable<int,int>* front_nodes - hash table used to keep
            track of the id's of the nodes on the front

    CArbMshTopo2D *bdry_topo - boundary 2D topology object


NON-MEMBER OPERATORS

operator == - ArbIntQdEdgeDesc equals operator

  inline int operator == (
          const CArbMshEdgeList::ArbIntQdEdgeDesc &op1,
          const CArbMshEdgeList::ArbIntQdEdgeDesc &op2)

    op1 - (in)  first operand
    op2 - (in)  second operand

  Description: This function defines the == operator for
      ArbIntQdEdgeDesc objects.

  Return Value: 1 if the operands have the same first node id, 0
      otherwishe.

*/

#endif
