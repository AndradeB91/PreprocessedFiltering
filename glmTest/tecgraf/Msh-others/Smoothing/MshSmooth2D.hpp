//
// ArbMshSmooth2D header file
//
// Description -
//   This is the header file for the ArbMshSmooth2D objects.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2001
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.9 $  $Date: 2003/11/14 18:23:01 $  $Author: wash $
//
  
#ifndef MshSmooth2D_h
#define MshSmooth2D_h

#include <cstdio>
#include <stdexcept>

#include "Dict.hpp"
#include "List.hpp"
#include "Vec2D.hpp"

#include "Msh2DTypes.hpp"
#include "MshTopo2D.hpp"

using FTools::Dict ;
using FTools::List ;
using FTools::Vec2D ;

namespace Msh2D {

class MshSmooth2D {

    public:

        class DuplicateNodeError : public std::runtime_error {
            public:
                DuplicateNodeError() : runtime_error("duplicate node id specified") {}
        } ;

        class DuplicateElemError : public std::runtime_error {
            public:
                DuplicateElemError() : runtime_error("duplicate element id specified") {}
        } ;

        class ElementOrderError : public std::runtime_error {
            public:
                ElementOrderError() : runtime_error("missmatched element polynomial order") {}
        } ;

        class ElementTypeError : public std::runtime_error {
            public:
                ElementTypeError() : runtime_error("unknown element type") {}
        } ;

        class NodeIterator ;
        class ElemIterator ;

    // constructors and destructors

        MshSmooth2D() ;
        MshSmooth2D(Dict<int,IntNode> *node_tbl,
                    Dict<int,MshElement2D> *elem_tbl) ;
        ~MshSmooth2D() ;

    // routines to add nodes and elements to the mesh

        void AddNode(int id,const Vec2D& coord,double z_coord) ;

        void AddElem(int elem_id,int mat_id,int num_nodes,const int nodes[]) ;

    // set and query parameters effecting smoothing

        void SetStopTolerance(const double tol) { It_tol = tol ; } ;
        void SetMaxIterations(const int maxits) { Max_it = maxits ; } ;

        double GetStopTolerance() { return(It_tol) ; } ;
        int    GetMaxIterations() { return(Max_it) ; } ;

    // methods for doing the smoothing

        void SmoothNodesLaplace() ;
        void SmoothNodesWinslow() ;
        void SmoothNodesConsLaplace() ;

    // method for fetching the new node locations

        int NumNodes() const ;
        int NumElements() const ;

        NodeIterator GetNodeIterator() const { return(NodeIterator(*this)) ; }
        ElemIterator GetElemIterator() const { return(ElemIterator(*this)) ; }

    private:

        struct EdgeKey {
            int id_0 ;
            int id_1 ;

            EdgeKey() : id_0(0), id_1(0) {} ;
            EdgeKey(const int id0,const int id1) :
                id_0(id0), id_1(id1) {} ;

            int operator == (const EdgeKey &op) {
                return((id_0 == op.id_0) && (id_1 == op.id_1)) ;
            }
        } ;

        Dict<int,MshElement2D> *ElemTable ;
        Dict<int,IntNode>      *NodeTable ;
        Dict<int,double>       *ZCoords ;

        MshOrder  Order ;
        bool   FirstElem ;
        bool   DeleteTables ;
        bool   FindBound ;
        double It_tol ;
        int    Max_it ;
        int    Num_it ;


        MshTopo2D *BuildMeshTopo() ;
        void FindBoundaryNodes() ;

        Vec2D ConsLaplaceUpdate(int num_elem,
                               List<int> *num_adj_nodes,
                               List<Vec2D> *adj_coords) ;
        Vec2D LaplaceUpdate(int num_elem,
                               List<int> *num_adj_nodes,
                               List<Vec2D> *adj_coords,
                               bool skip_checks) ;

        void MoveSideNodes() ;

        void DoSmoothNodesLaplace(MshTopo2D *msh_topo,
                                  int max_iter,double adapt_tol,
                                  bool skip_checks = false) ;
        void DoSmoothNodesWinslow(MshTopo2D *msh_topo,
                                  int max_iter,double adapt_tol) ;
        void DoSmoothNodesConsLaplace(MshTopo2D *msh_topo,
                                      int num_iters) ;

        double FindMinCharEdgeLength(MshTopo2D *msh_topo) ;

        double Angle(Vec2D b,Vec2D i,
                                   Vec2D j) const ;
        double CrossProd(Vec2D b,Vec2D i,
                                       Vec2D j) const ;

        double TriMetric(Vec2D i,Vec2D j,
                         Vec2D k) const ;
        double QuadMetric(Vec2D i,Vec2D j,
                          Vec2D k,Vec2D l) const ;
 
        void DisplayMesh(char *label,FILE *fd = 0) ;

    public:

        class NodeIterator {
            public:
                NodeIterator(const MshSmooth2D &smooth) :
                    iter(smooth.NodeTable),zc(smooth.ZCoords) {} ;
                NodeIterator() {}
                NodeIterator(const NodeIterator& other) { Copy(other) ; }
                NodeIterator operator = (const NodeIterator& other) {
                    Copy(other) ;
                    return *this ;
                }
                void First() { iter.First() ; } ;
                void Next()  { iter.Next() ; } ;
                bool More()  { return(iter.More()) ; } ;
                int  Id()    { return(iter.Key()) ; } ;
                double NCoordI(int i) {
                    if (i < 2) return(iter.Entry().coord[i]) ;
                    if (zc == 0) return(0.0) ;
                    return zc->GetValue(iter.Key()) ;
                }
                Vec2D Coord2D() { return(iter.Entry().coord) ; } ;

                void operator ++ ()      { Next() ; } ;
                void operator ++ (int i) { Next() ; } ;

            private:
                Dict<int,IntNode>::DictIterator iter ;
                Dict<int,double>* zc ;
                MshNodeType ntype ;

                void Copy(const NodeIterator& other) {
                    iter = other.iter ;
                    zc = other.zc ;
                    ntype = other.ntype ;
                }
        } ;

        class ElemIterator {
            public:
                ElemIterator(const MshSmooth2D &smooth) :
                    iter(smooth.ElemTable) {}
                ElemIterator() {}
                ElemIterator(const ElemIterator& other) { Copy(other) ; }
                ElemIterator operator = (const ElemIterator& other) {
                    Copy(other) ;
                    return *this ;
                }
                void First() { iter.First() ; };
                void Next()  { iter.Next() ; };
                bool More()  { return(iter.More()) ; } ;
                int  Id()    { return(iter.Entry().elem_id) ; } ;
                int  Mat()   { return(iter.Entry().mat_id) ; } ;
                int  NumNodes() { return(iter.Entry().num_nodes) ; } ;
                int  Node(int i) { return(iter.Entry().nodes[i]) ; } ;

                void operator ++ ()      { Next() ; } ;
                void operator ++ (int i) { Next() ; } ;

            private:
                Dict<int,MshElement2D>::DictIterator iter ;

                void Copy(const ElemIterator& other) {
                    iter = other.iter ;
                }
       } ;

    friend int DictHashIndex(const EdgeKey& key) ;
} ;

inline int DictHashIndex(const MshSmooth2D::EdgeKey& key) {
    return key.id_0 ;
}

/*
CLASS MshSmooth2D

  This object provides an interface to the nodal smoothing routines. 
  There are two typical ways to use this object. In the first, the 
  client calls AddNode and AddElement to define a mesh. The Smoothing 
  routines are calles, and then the client calls GetNodes to retrieve 
  the updated nodal coordinates (a GetElements function is provided, 
  but the element information is not modified). In this mode, the 
  object determines the "floating" nodes. That is, the nodes that are 
  not part of the external boundary or on a bi-material interface. 

  The second way to use this object is to pass a node and element table 
  to the constructor. In this case it is not necessary to call the 
  Add/Get Node/Element routines, and it is the client's responsibility 
  to flag the "floating" nodes. 


PUBLIC INTERFACE

  Public Typedefs:

        Dict<int,MshElement2D> ElemHash - This table maps 
                from an element id to the element info. 

        Dict<int,IntNode> NodeHash - This table maps from a 
                node id to the node data. 

  Public Member Functions:

    MshSmooth2D - no argument constructor 

      MshSmooth2D()

      Description: This is a no argument constructor for a 
          MshSmooth2D object. 


    MshSmooth2D - constructor 

      MshSmooth2D(
              NodeHash *node_tbl,
              ElemHash *elem_tbl)

        node_tbl - (in)  node table 
        elem_tbl - (in)  element table 

      Description: This is a constructor for a MshSmooth2D object 
          with the node and element tables passed as arguments. 


    MshSmooth2D - destructor 

      ~MshSmooth2D()

      Description: This is a no destructor for a MshSmooth2D 
          object. 


    AddNode - add a node to the mesh 

      void AddNode(const ArbMshNode &inode)

        inode - (in)  description of the node 

      Description: This function adds a new node to the mesh. 

      Exceptions:
          CArbMshCDuplNode - duplicate node id's

    AddNodeList - add an array of nodes to the mesh 

      void AddNodeList(
              const int        num_nodes,
              const ArbMshNode *constnode_list)

        num_nodes - (in)  number of nodes 
        node_list - (in)  list of nodes to add 

      Description: This function adds an array of nodes to the mesh. 

      Exceptions:
          CArbMshCDuplNode - duplicate node id's

    AddElem - add an element to the mesh 

      void AddElem(const MshElement2D &elem)

        elem - (in)  element to add 

      Description: This function adds an element to the mesh. 

      Exceptions:
          CArbMshCDuplElem - duplicate element id's
          CArbMshCInvalidElem - invalid element

    AddElemList - add an array of elements to the mesh 

      void AddElemList(
              const int             num_elems,
              const MshElement2D *constelem_list)

        num_elems - (in)  number of elements 
        elem_list - (in)  elements 

      Description: This function adds an array of elements to the 
          mesh. 

      Exceptions:
          CArbMshCDuplElem - duplicate element id's
          CArbMshCInvalidElem - invalid element

    SetStopTolerance - set the smoothing stop tolerance 

      void SetStopTolerance(const double tol)

        tol - (in)  the tolerance value 

      Description: This function sets the smoothing iteration stop 
          tolerance. The smoothing loops stops when no node is moved 
          more than this factor times the maximum node movement 
          during the first iteration. 


    SetMaxIterations - set the smoothing maximum iteratons 

      void SetMaxIterations(const int maxits)

        maxits - (in)  the maximum iteratons 

      Description: This function sets the maximum number of smoothing 
          iterations. 


    GetStopTolerance - get the smoothing stop tolerance 

      double GetStopTolerance()

      Description: This function returns the smoothing stop 
          tolerance. 

      Return Value: the smoothing stop tolerance 


    GetMaxIterations - get the smoothing maximum iteratons 

      int GetMaxIterations()

      Description: This function returns the maximum number of 
          smoothing iterations. 

      Return Value: the smoothing maximum iteratons 


    SmoothNodesLaplace - smooth nodes using a Laplace algorithm 

      void SmoothNodesLaplace()

      Description: This function smooths the internal nodes using a 
          Laplace algorithm. 


    SmoothNodesWinslow - smooth nodes using a Winslow algorithm 

      void SmoothNodesWinslow()

      Description: This function smooths the internal nodes using a 
          Winslow algorithm. 


    SmoothNodesConsLaplace - smooth nodes using a "constrained" 
                             Laplace algorithm 

      void SmoothNodesConsLaplace()

      Description: This function smooths the internal nodes using a 
          "constrained" Laplace algorithm. 


    NumNodes - return the number of nodes in the mesh 

      int NumNodes() const

      Description: This function returns the current number of nodes 
          in a mesh (typically called after smoothing). 

      Return Value: number of nodes 


    NumElements - return the number of elements in the mesh 

      int NumElements() const

      Description: This function returns the current number of 
          elements in a mesh (typically called after smoothing). 

      Return Value: number of elements 


    GetNodes - return a list of all mesh nodes 

      ArbMshNode *GetNodes() const

      Description: This function returns a list of all the nodes 
          current defined for the mesh (typically called after 
          smoothing). 

      Return Value: A list of all mesh nodes. Ownership of this 
          memory passes to the client, which must eventually call 
          delete []. 


    GetElements - return a list of all mesh elements 

      MshElement2D *GetElements() const

      Description: This function returns a list of all the elements 
          current defined for the mesh (typically called after 
          smoothing). 

      Return Value: A list of all mesh elements. Ownership of this 
          memory passes to the client, which must eventually call 
          delete []. 


PRIVATE INTERFACE

  Private Typedefs:

        Dict<int,MshElement2D>::DictIterator ElemIter - Iterator 

                type for the element table. 

        Dict<int,IntNode>::DictIterator NodeIter - Iterator type 

                for the node table. 

  Private Member Functions:

    BuildMeshTopo - build a mesh topology object 

      MshTopo2D *BuildMeshTopo()

      Description: This function builds a MshTopo2D object for the 
          mesh described by the node and element tables. 

      Return Value: a new MshTopo2D object 


    FindBoundaryNodes - find the fixed nodes 

      void FindBoundaryNodes()

      Description: This function finds the nodes on the boundary and 
          on bi-material interfaces, and flags them as being fixed. 


    ConsLaplaceUpdate - find an updated node coordinate 

      Vec2D ConsLaplaceUpdate(
              int                     num_elem,
              List<int>*         num_adj_nodes,
              List<Vec2D>* adj_coords)

        num_elem      - (in)  number of elements around the node 
        num_adj_nodes - (in)  number of nodes for each of the 
                              adjacent elements 
        adj_coords    - (in)  list of coordinates of nodes for each 
                              of the adjacent elements 

      Description: This function determines the updated location of a 
          node during constrained Laplace smoothing. 

      Return Value: the new coordinate location 


    LaplaceUpdate - find an updated node coordinate 

      Vec2D LaplaceUpdate(
              int                     num_elem,
              List<int>*         num_adj_nodes,
              List<Vec2D>* adj_coords)

        num_elem      - (in)  number of elements around the node 
        num_adj_nodes - (in)  number of nodes for each of the 
                              adjacent elements 
        adj_coords    - (in)  list of coordinates of nodes for each 
                              of the adjacent elements 

      Description: This function determines the updated location of a 
          node during normal Laplace smoothing 

      Return Value: the new coordninate location 


    MoveSideNodes - update the coordinates of the side nodes 

      void MoveSideNodes()

      Description: This function updates the coordinates of all the 
          side nodes (in the case of a quadratic order mesh) to place 
          them at them midway between the corresponding corner nodes. 


    DoSmoothNodesLaplace - do Laplace smoothing 

      void DoSmoothNodesLaplace()

      Description: This function does the actual work for the Laplace 
          smoothing. 


    DoSmoothNodesWinslow - do Winslow smoothing 

      void DoSmoothNodesWinslow()

      Description: This function does the actual work for the Winslow 
          smoothing. 


    DoSmoothNodesConsLaplace - do constrained smoothing 

      void DoSmoothNodesConsLaplace(int num_iters)

        num_iters - (in)  number of smoothing iterations 

      Description: This function does the actual work for the 
          constrained Laplace smoothing. 


    Angle - compute an angle 

      double Angle(
              Vec2D b,
              Vec2D i,
              Vec2D j) const

        b - (in)  base vertex coordinate 
        i - (in)  end of first vector 
        j - (in)  end of second vector 

      Description: Compute the magnitude of the included angle 
          between two vectors that share a common starting vertex (in 
          the range -pi to pi). 

      Return Value: the magnitude of the angle (in radians) 


    CrossProd - compute a cross product 

      double CrossProd(
              Vec2D b,
              Vec2D i,
              Vec2D j) const

        b - (in)  base vertex coordinate 
        i - (in)  end of first vector 
        j - (in)  end of second vector 

      Description: This method computes the cross product of two 
          vector that share a common starting vertex. 

      Return Value: the (2D) cross product value 


    TriMetric - compute a shape metric 

      double TriMetric(
              Vec2D b,
              Vec2D i,
              Vec2D j) const

        b - (in)  first vertex 
        i - (in)  second vertex 
        j - (in)  third vertex 

      Description: This method computes a shape metric for a 
          triangular element. A zero or negative metric indicates an 
          invalid element. The maximum possible metric value is 1.0 
          for an equilateral triangle. 

      Return Value: the shape metric 


    QuadMetric - compute a shape metric 

      double QuadMetric(
              Vec2D b,
              Vec2D i,
              Vec2D j,
              Vec2D l) const

        b - (in)  first vertex 
        i - (in)  second vertex 
        j - (in)  third vertex 
        l - (in)  forth vertex 

      Description: This method computes a shape metric for a 
          quadrilateral element. A zero or negative metric indicates 
          an invalid element. The maximum possible metric value is 
          1.0 for a square. 

      Return Value: the shape metric 


    DisplayMesh - debug routine to display a mesh 

      void DisplayMesh()

      Description: This method provides debugging support. It prints 
          the commands necessary to display the mesh. 


  Private Member Variables:

        ElemHash *ElemTable - hash table of elements in the mesh 

        NodeHash *NodeTable - hash table of nodes in the mesh 

        MshOrder Order - LINEAR or QUADRATIC 

        bool FirstElem - if true there are no elements in the mesh yet 

        bool DeleteTables - if true the desctructor should delete the 
                node and element tables 

        bool FindBound - if true the floating nodes will be found 
                before smoothing 

        double It_tol - smoothing iteration tolerance 

        int Max_it - smoothing maximum number of iterations 

*/

} // namespace

#endif
