//
// ArbMsh header file
//
// Description -
//   This is a header file for the ArbMsh objects.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.7 $  $Date: 2001/08/03 18:28:28 $  $Author: wash $
//

#ifndef ArbMsh_h
#define ArbMsh_h

#include "ArbCoord2D.hpp"

/* ------------------------------------------------------------
    ArbMshNode - a data structure used to pass node descriptions
*/

struct ArbMshNode {
    int id ;               // unique node id
    double coord[3] ;      // nodal coordinates
    bool retained_flag ;   // retain this node if on a boundary
} ;


/* ------------------------------------------------------------
    ArbMshEdge - a data structure used to specify the id's of
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
    CArbCoord2D coord ;              // nodal coordinates
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

inline int ArbHashIndex(const ArbEdgeKey &edge)
{
    return(edge.id_0 + edge.id_1) ;
} ;

#endif

