// topology_struc.hpp: structure to topology class.
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <deque>
#include <map>

#include "Vec3D.hpp"

using namespace FTools;
using namespace std;

// TopEdge_
//////////////////////////////////////////////////////////////////////
class TopEdge_
{
public:

  static TopEdge_ *newEdge ( );

  TopEdge_ ( )
  {
    id[0] = -1;
    id[1] = -1;
    type  = -1;
    size  = 0.0;
    flag  = 0;
    next  = NULL;
    prev  = NULL;
    valid = 1;
    extraInfo = 0;
  }


  int    i;                  // edge index
  int    id[2];              // node indexes
  int    type;               // 1 - bondary, 2 - Domain, 2+ - Restriction
  double size;               // size of the edge

  int    flag;               // 0 - elemento nao percorrido, 1 - ja percorrido

  int    valid;

  int    extraInfo;

  vector<TopElem_ *> AdjElms; //Adjacent elements

  TopEdge_ *next;             // ponteiro para a lista de arestas
  TopEdge_ *prev;             // ponteiro para a lista de arestas
};



#define MAX_INC  10

// TopElem_
//////////////////////////////////////////////////////////////////////
class TopElem_
{

public:

  TopElem_ ( )
  {
    id     = -1;
    type   = 0;
    flag   = 0;
    active = 0;
    next   = NULL;
    surfId = -1;
    extraInfo = 0;
  }

  int       id;                // element index
  int       type;              // element type
  int       inc[MAX_INC];      // incidence
  
  Vec3D     normal;            // normal
  int       surfId;

  vector<TopEdge_ *> AdjEdges;  // ponteiro para arestas
  vector<int>       cycles;    // 0 - ids nao trocados, 1 - ids trocados

  int       flag;              // 0 - elemento nao percorrido, 1 - ja percorrido

  int    extraInfo;

  int       active;

  int       valid;

  TopElem_   *next;      // next
};


// TopNode_
//////////////////////////////////////////////////////////////////////
class TopNode_
{

public:
  TopNode_ ( )
  {
    bound  = 0;
    flag   = 0;
    active = 0;
    next   = NULL;
    extraInfo = 0;
  }

  Vec3D      coords;  
  Vec3D      normal; 

  int        bound;           // 1-boundary node

  vector<TopEdge_ *> AdjEdges; // adjacent edges
  
  int          flag;          // in LA list
  int          active;

  int    extraInfo;


  TopNode_ *next;  // next
};


// TopBoundary
//////////////////////////////////////////////////////////////////////
class TopBoundary
{ 
public:
  vector<TopEdge_ *> Edges;  // list to edges
};

