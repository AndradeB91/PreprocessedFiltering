//
// Msh2DTypes header file
//
// Description -
//   This file defines data types common to the Msh2D routines
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2007
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//

#ifndef Msh2DTypes_h
#define Msh2DTypes_h

#include "Vec2D.hpp"

using FTools::Vec2D ;

namespace Msh2D {

// enumerated types

enum MshOrder {LINEAR, QUADRATIC} ;

enum MshElemType {TRIANGLE, QUADRILATERAL} ;

enum MshNodeType {BOUNDARY, INTERIOR, UNSPECIFIED} ;

enum MshNodeMotion {MSH_FIXED, MSH_FLOATING} ;


// structures

struct IntNode {
    int id ;                   // unique node id
    Vec2D coord ;              // nodal coordinates
    MshNodeType type ;         // boundary or interior
    MshNodeMotion motion ;     // fixed or floating
    bool corner ;              // corner node flag

    IntNode() {}
    IntNode(int iid,const Vec2D& icoord,
            MshNodeType itype = UNSPECIFIED,
            MshNodeMotion imotion = MSH_FLOATING,
            bool icorner = true) :
        id(iid),coord(icoord),type(itype),motion(imotion),corner(icorner) {}

    int operator == (const IntNode &op) {
        return(id == op.id) ;
    }
} ;

inline int DictHashIndex(const IntNode &op)
{
    return(op.id) ;
}


struct IntEdge {
    int node_id[3] ;      // node id's
    MshNodeType type ;    // boundary or interior
    double length ;       // length of the edge
    double tol ;

    IntEdge() {}
    IntEdge(int n0,int n1,int n2,MshNodeType tp) :
        type(BOUNDARY),length(-1.0),tol(-1.0) {
            node_id[0] = n0 ; node_id[1] = n1 ; node_id[2] = n2 ;
    }

    int operator == (const IntEdge &op) {
        return((node_id[0] == op.node_id[0]) &&
               (node_id[1] == op.node_id[1])) ;
    }
} ;

inline int DictHashIndex(const IntEdge &op)
{
    return(op.node_id[0]) ;
}


/* ------------------------------------------------------------
    ArbMshElement2D - a data structure used to desribe a 2D
                      element.  The nodes array contains the
                      id's of the element's nodes, specified
                      as follows:

        2         3         2       2         3    6    2
         +         +-------+         +         +---+---+
        / \        |       |        / \        |       |
       /   \       |       |     5 +   + 4   7 +       + 5
      /     \      |       |      /     \      |       |
     +-------+     +-------+     +---+---+     +---+---+
    0         1   0         1   0    3    1   0    4    1
*/

struct MshElement2D {
    int elem_id ;          // unique element id 
    int mat_id ;           // id of this element's material
    int num_nodes ;        // number of element nodes
    int nodes[8] ;         // list of node id's

    MshElement2D() {}
    MshElement2D(int ielem_id,int imat_id,int inum_nodes,const int inodes[]) :
        elem_id(ielem_id),mat_id(imat_id),num_nodes(inum_nodes) {
        for (int i=0 ; i<num_nodes ; ++i) nodes[i] = inodes[i] ;
    }

    int operator == (const MshElement2D &op) {
        return(elem_id == op.elem_id) ;
    }
} ;

inline int DictHashIndex(const MshElement2D &op)
{
    return(op.elem_id) ;
}


#if 0
/* ------------------------------------------------------------
    ArbMshNode - a data structure used to pass node descriptions
*/

struct ArbMshNode {
    int id ;               // unique node id
    double coord[3] ;      // nodal coordinates
    bool retained_flag ;   // retain this node if on a boundary
} ;


/* ------------------------------------------------------------
    ArbMshNode - a data structure used to specify the id's of
                 nodes that define an edge.  The two end nodes
                 are node_id elements 0 and 1.  Element 2 is the
                 side node if QUADRATIC elements are being
                 generated.  It is ignored otherwise.

*/

struct ArbMshEdge {
    int node_id[3] ;       // id's of nodes defining edge
} ;

/* equals operator for MshEdges */
inline int operator == (const ArbMshEdge &op1,const ArbMshEdge &op2)
{
    return((op1.node_id[0] == op2.node_id[0]) &&
           (op1.node_id[1] == op2.node_id[1])) ;
}


/* ------------------------------------------------------------
    ArbMshElement2D - a data structure used to desribe a 2D
                      element.  The nodes array contains the
                      id's of the element's nodes, specified
                      as follows:

        2         3         2       2         3    6    2
         +         +-------+         +         +---+---+
        / \        |       |        / \        |       |
       /   \       |       |     5 +   + 4   7 +       + 5
      /     \      |       |      /     \      |       |
     +-------+     +-------+     +---+---+     +---+---+
    0         1   0         1   0    3    1   0    4    1
*/

struct ArbMshElement2D {
    int elem_id ;          // unique element id 
    int mat_id ;           // id of this element's material
    int num_nodes ;        // number of element nodes
    int nodes[8] ;         // list of node id's
} ;

/* equals operator for MshElements */
inline int operator == (const ArbMshElement2D &op1,
                        const ArbMshElement2D &op2)
{
    return(op1.elem_id == op2.elem_id) ;
}


//
// Enumerated types
//

typedef enum {LINEAR, QUADRATIC} ArbMshOrder ;

typedef enum {TRIANGLE, QUADRILATERAL, MIXED} ArbMshElemType ;

typedef enum {BOUNDARY, INTERIOR, UNSPECIFIED} ArbMshNodeType ;

typedef enum {ARB_FIXED, ARB_FLOATING} ArbMshNodeMotion ;


//
// Data Structure for mesh statistics
//

struct ArbMshStats {
    int num_elements ;
    int num_nodes ;
    double mean_shape_measure ;
    double minimum_shape_measure ;
    int num_less_than[11] ;
} ;


// --------------------------------------------------
// these will move to a private header file


struct ArbIntNode {
    int id ;              // unique node id
    Vec2D coord ;              // nodal coordinates
    ArbMshNodeType type ;          // boundary or interior
    ArbMshNodeMotion motion ;      // fixed or floating
    bool corner ;                  // corner node flag
} ;

inline int operator == (const ArbIntNode &op1,const ArbIntNode &op2)
{
    return(op1.id == op2.id) ;
}

struct ArbIntEdge {
    int node_id[3] ;      // node id's
    ArbMshNodeType type ;          // boundary or interior
    double length ;                // length of the edge
} ;

inline int operator == (const ArbIntEdge &op1,const ArbIntEdge &op2)
{
    return((op1.node_id[0] == op2.node_id[0]) &&
           (op1.node_id[1] == op2.node_id[1])) ;
}

inline int DictHashIndex(const ArbIntEdge &op)
{
    return(op.node_id[0] + op.node_id[1]) ;
}


struct ArbQdEdge {
    int nd0, nd1 ;
    short end_code ;
    short level ;
    double length ;
} ;

inline int operator == (const ArbQdEdge &op1,const ArbQdEdge &op2)
{
    return((op1.nd0 == op2.nd0) && (op1.nd1 == op2.nd1)) ;
}


struct ArbEdgeKey {
    int id_0 ;
    int id_1 ;
    ArbEdgeKey() : id_0(0), id_1(0) {} ;
    ArbEdgeKey(const int id0,const int id1) :
        id_0(id0), id_1(id1) {} ;
} ;

#if 0
inline int operator == (const ArbEdgeKey &op1,const ArbEdgeKey &op2)
{
    return((op1.id_0 == op2.id_0) && (op1.id_1 == op2.id_1)) ;
}
#endif
extern int operator == (const ArbEdgeKey &op1,const ArbEdgeKey &op2) ;

inline int DictHashIndex(const ArbEdgeKey &edge)
{
    return(edge.id_0 + edge.id_1) ;
}

#endif

} // namespace

#endif

