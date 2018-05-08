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
//   $Revision: 1.8 $  $Date: 2002/06/12 18:48:41 $  $Author: wash $
//
  
#ifndef ArbMshSmooth2D_h
#define ArbMshSmooth2D_h

#include "ArbMsh.hpp"
#include "ArbArray.hpp"
#include "ArbCoord2D.hpp"
#include "ArbHashTable.hpp"
#include "ArbMshTopo2D.hpp"

// status returned by this object

#ifndef ARB_NORMAL_STATUS
#define ARB_NORMAL_STATUS      0
#endif
#ifndef ARB_DUPLICATE_NODE_ID
#define ARB_DUPLICATE_NODE_ID  1
#endif
#ifndef ARB_DUPLICATE_ELEM_ID
#define ARB_DUPLICATE_ELEM_ID  2
#endif
#ifndef ARB_INVALID_ELEM
#define ARB_INVALID_ELEM       4
#endif
#ifndef ARB_INVALID_NODE
#define ARB_INVALID_NODE       5
#endif
 
// typedef CArbSet<int> CornerNodes ;

class CArbMshSmooth2D {

    public:

        typedef CArbHashTable<int,ArbMshElement2D> ElemHash ;
        typedef CArbHashTable<int,ArbIntNode>      NodeHash ;

    // constructors and destructors

        CArbMshSmooth2D() ;
        CArbMshSmooth2D(NodeHash *node_tbl,ElemHash *elem_tbl) ;
        ~CArbMshSmooth2D() ;

    // routines to add nodes and elements to the mesh

         int AddNode(const ArbMshNode &inode) ;

         int AddNodeList(const int num_nodes,
                         const ArbMshNode *const node_list) ;

         int AddElem(const ArbMshElement2D &elem) ;

         int AddElemList(const int num_elems,
                         const ArbMshElement2D *const elem_list) ;

    // set and query parameters effecting smoothing

        void SetStopTolerance(const double tol) { It_tol = tol ; } ;
        void SetMaxIterations(const int maxits) { Max_it = maxits ; } ;

        double GetStopTolerance() { return(It_tol) ; } ;
        int    GetMaxIterations() { return(Max_it) ; } ;

    // methods for doing the smoothing

        void SmoothNodesLaplace() ;
        int SmoothNodesWinslow() ;
        void SmoothNodesConsLaplace() ;

    // method for fetching the new node locations

        int NumNodes() const ;
        int NumElements() const ;

        ArbMshNode *GetNodes() const ;
        ArbMshElement2D *GetElements() const ;

    private:

        typedef CArbHashTableIterator<int,ArbMshElement2D> ElemIter ;
        typedef CArbHashTableIterator<int,ArbIntNode>      NodeIter ;

        typedef CArbHashTable<int,double>  ZCoordHash ;

        ElemHash   *ElemTable ;
        NodeHash   *NodeTable ;
        ZCoordHash *ZCoords ;

        ArbMshOrder  Order ;
        bool   FirstElem ;
        bool   DeleteTables ;
        bool   FindBound ;
        double It_tol ;
        int    Max_it ;
        int    Num_it ;


        CArbMshTopo2D *BuildMeshTopo() ;
        void FindBoundaryNodes() ;

        CArbCoord2D ConsLaplaceUpdate(int num_elem,
                               CArbArray<int> *num_adj_nodes,
                               CArbArray<CArbCoord2D> *adj_coords) ;
        CArbCoord2D LaplaceUpdate(int num_elem,
                               CArbArray<int> *num_adj_nodes,
                               CArbArray<CArbCoord2D> *adj_coords,
                               bool skip_checks) ;

        void MoveSideNodes() ;

        int DoSmoothNodesLaplace(CArbMshTopo2D *msh_topo,
                                  int max_iter,double adapt_tol,
                                  bool skip_checks = false) ;
        void DoSmoothNodesWinslow(CArbMshTopo2D *msh_topo,
                                  int max_iter,double adapt_tol) ;
        int DoSmoothNodesConsLaplace(CArbMshTopo2D *msh_topo,
                                      int num_iters) ;

        double FindMinCharEdgeLength(CArbMshTopo2D *msh_topo) ;

        double Angle(CArbCoord2D b,CArbCoord2D i,
                                   CArbCoord2D j) const ;
        double CrossProd(CArbCoord2D b,CArbCoord2D i,
                                       CArbCoord2D j) const ;

        double TriMetric(CArbCoord2D i,CArbCoord2D j,
                         CArbCoord2D k) const ;
        double QuadMetric(CArbCoord2D i,CArbCoord2D j,
                          CArbCoord2D k,CArbCoord2D l) const ;
 
        void DisplayMesh(char *label) ;
} ;

/*
CLASS CArbMshSmooth2D

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

        CArbHashTable<int,ArbMshElement2D> ElemHash - This table maps 
                from an element id to the element info. 

        CArbHashTable<int,ArbIntNode> NodeHash - This table maps from a 
                node id to the node data. 

  Public Member Functions:

    CArbMshSmooth2D - no argument constructor 

      CArbMshSmooth2D()

      Description: This is a no argument constructor for a 
          CArbMshSmooth2D object. 


    CArbMshSmooth2D - constructor 

      CArbMshSmooth2D(
              NodeHash *node_tbl,
              ElemHash *elem_tbl)

        node_tbl - (in)  node table 
        elem_tbl - (in)  element table 

      Description: This is a constructor for a CArbMshSmooth2D object 
          with the node and element tables passed as arguments. 


    CArbMshSmooth2D - destructor 

      ~CArbMshSmooth2D()

      Description: This is a no destructor for a CArbMshSmooth2D 
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

      void AddElem(const ArbMshElement2D &elem)

        elem - (in)  element to add 

      Description: This function adds an element to the mesh. 

      Exceptions:
          CArbMshCDuplElem - duplicate element id's
          CArbMshCInvalidElem - invalid element

    AddElemList - add an array of elements to the mesh 

      void AddElemList(
              const int             num_elems,
              const ArbMshElement2D *constelem_list)

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

      ArbMshElement2D *GetElements() const

      Description: This function returns a list of all the elements 
          current defined for the mesh (typically called after 
          smoothing). 

      Return Value: A list of all mesh elements. Ownership of this 
          memory passes to the client, which must eventually call 
          delete []. 


PRIVATE INTERFACE

  Private Typedefs:

        CArbHashTableIterator<int,ArbMshElement2D> ElemIter - Iterator 
                type for the element table. 

        CArbHashTableIterator<int,ArbIntNode> NodeIter - Iterator type 
                for the node table. 

  Private Member Functions:

    BuildMeshTopo - build a mesh topology object 

      CArbMshTopo2D *BuildMeshTopo()

      Description: This function builds a MshTopo2D object for the 
          mesh described by the node and element tables. 

      Return Value: a new MshTopo2D object 


    FindBoundaryNodes - find the fixed nodes 

      void FindBoundaryNodes()

      Description: This function finds the nodes on the boundary and 
          on bi-material interfaces, and flags them as being fixed. 


    ConsLaplaceUpdate - find an updated node coordinate 

      CArbCoord2D ConsLaplaceUpdate(
              int                     num_elem,
              CArbArray<int>*         num_adj_nodes,
              CArbArray<CArbCoord2D>* adj_coords)

        num_elem      - (in)  number of elements around the node 
        num_adj_nodes - (in)  number of nodes for each of the 
                              adjacent elements 
        adj_coords    - (in)  list of coordinates of nodes for each 
                              of the adjacent elements 

      Description: This function determines the updated location of a 
          node during constrained Laplace smoothing. 

      Return Value: the new coordinate location 


    LaplaceUpdate - find an updated node coordinate 

      CArbCoord2D LaplaceUpdate(
              int                     num_elem,
              CArbArray<int>*         num_adj_nodes,
              CArbArray<CArbCoord2D>* adj_coords)

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
              CArbCoord2D b,
              CArbCoord2D i,
              CArbCoord2D j) const

        b - (in)  base vertex coordinate 
        i - (in)  end of first vector 
        j - (in)  end of second vector 

      Description: Compute the magnitude of the included angle 
          between two vectors that share a common starting vertex (in 
          the range -pi to pi). 

      Return Value: the magnitude of the angle (in radians) 


    CrossProd - compute a cross product 

      double CrossProd(
              CArbCoord2D b,
              CArbCoord2D i,
              CArbCoord2D j) const

        b - (in)  base vertex coordinate 
        i - (in)  end of first vector 
        j - (in)  end of second vector 

      Description: This method computes the cross product of two 
          vector that share a common starting vertex. 

      Return Value: the (2D) cross product value 


    TriMetric - compute a shape metric 

      double TriMetric(
              CArbCoord2D b,
              CArbCoord2D i,
              CArbCoord2D j) const

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
              CArbCoord2D b,
              CArbCoord2D i,
              CArbCoord2D j,
              CArbCoord2D l) const

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

        ArbMshOrder Order - LINEAR or QUADRATIC 

        bool FirstElem - if true there are no elements in the mesh yet 

        bool DeleteTables - if true the desctructor should delete the 
                node and element tables 

        bool FindBound - if true the floating nodes will be found 
                before smoothing 

        double It_tol - smoothing iteration tolerance 

        int Max_it - smoothing maximum number of iterations 

*/

#endif
