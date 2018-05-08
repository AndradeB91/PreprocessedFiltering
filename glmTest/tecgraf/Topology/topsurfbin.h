/* topology.h:
 *
*/

#ifndef _SURF_TOPLOLOGY_BIN_H_
#define _SURF_TOPLOLOGY_BIN_H_

#include "deflib.h"

#ifdef __cplusplus
extern "C" {
#endif

MSHSURF_API void*        SurfTopInsertMesh        (int n_node, double *coords, int n_elem, int *conn,
                                       void (*func) (const char *));
MSHSURF_API void         SurfTopBoundBox          (void *surf, double min[3], double max[3]);
MSHSURF_API void         SurfTopRelease           (void *surf);

MSHSURF_API int          writeNF                  (void *surf, char *filename);

/* number of nodes, edges, bound edges, and elements */
MSHSURF_API int          SurfTopNumNodes          (void *surf);
MSHSURF_API int          SurfTopNumEdge           (void *surf);
MSHSURF_API int          SurfTopNumBoundEdge      (void *surf);
MSHSURF_API int          SurfTopNumElems          (void *surf);

/* nodes */
MSHSURF_API int          SurfTopGetCoordNode      (void *surf, int id, double coord[3]);
MSHSURF_API void         SurfTopSetCoordNode      (void *surf, int id, double *new_coord);
MSHSURF_API int          SurfTopIsBoundaryNode    (void *surf, int id);
MSHSURF_API int          SurfTopGetNormalNode     (void *surf, int id, double normal[3]);
MSHSURF_API void         SurfTopSetNormalNode     (void *surf, int id, double normal[3]);
MSHSURF_API void         SurfTopSetPtsInative     (void *surf, int id);
MSHSURF_API void         SurfTopUdateNormalNodes  (void *surf);
MSHSURF_API int          SurfTopIsActNode         (void *surf, int id);
MSHSURF_API int          SurfTopReturnNodeId      (void *surf, double x, double y, double z);

/* edges */
MSHSURF_API int          SurfTopGetEdge           (void *surf, int i, int *ei, int *ej);
MSHSURF_API int          SurfTopGetEdgeSize       (void *surf, int i, double mid[3], double *size);
MSHSURF_API int          SurfTopGetBoundEdge      (void *surf, int i, int *ei, int *ej);
MSHSURF_API int          SurfTopGetBEdgeSize      (void *surf, int i, double mid[3], double *size);
MSHSURF_API int          SurfTopGetIsEdgeValid    (void *surf, int i);
MSHSURF_API double       SurfTopGetSmalestEdge    (void *surf);

/* elements */
MSHSURF_API int          SurfTopGetElemNNodes     (void *surf, int id); 
MSHSURF_API int*         SurfTopGetElemConn       (void *surf, int id); 
MSHSURF_API int          SurfTopGetElemNorm       (void *surf, int id, double normal[3]); 
MSHSURF_API int          SurfTopGetElemIdEdge     (void *surf, int id, int pos_adj); 
MSHSURF_API int          SurfTopGetIsElemValid    (void *surf, int id);
MSHSURF_API int          SurfTopGetElemCenter     (void *surf, int id, double center[3]);
MSHSURF_API int          SurfTopIsActElem         (void *surf, int id);
MSHSURF_API void         SurfTopSetElemInfo       (void *surf, int id, int flag);
MSHSURF_API int          SurfTopGetElemInfo       (void *surf, int id);


/* loops */
MSHSURF_API int          SurfTopGetNumLoops       (void *surf);
MSHSURF_API int          SurfTopGetLoop           (void *surf, int id, int *npts, int **idpts);

#if 0
/* special function */
MSHSURF_API int          SurfTopGetNetNodes       (void *surf, int idnode, int size, int **net);
#endif

/* node - nodes */
MSHSURF_API int     SurfTopNumAdjNodeNode    (void *surf, int id);
MSHSURF_API void    SurfTopAdjNodeNode       (void *surf, int id, int *adj_nodes);
MSHSURF_API int     SurfTopOrientAdjNodeNode (void *surf, int id, int *adj_nodes);

/* node - edges */
MSHSURF_API int     SurfTopNumAdjEdgeToNode  (void *surf, int id);
MSHSURF_API int     SurfTopAdjEdgeToNode     (void *surf, int id, int pos);
MSHSURF_API int     SurfTopOppAdjEdgeToNode  (void *surf, int id, int *ne, int **edges); /* opposite adjacent edges to node */
/* node - elements */
MSHSURF_API int     SurfTopAdjElemToNode     (void *surf, int id, int *ne, int **elem);
/* edges - nodes */
MSHSURF_API int     SurfTopAdjNodeToEdge     (void *surf, int id, int *ei, int *ej);
/* edges - elements */
MSHSURF_API int     SurfTopNumAdjElemToEdge  (void *surf, int id);
MSHSURF_API int     SurfTopAdjElemToEdge     (void *surf, int id, int pos);
/* elements - edges */
MSHSURF_API int     SurfTopNumAdjEdgeToElem  (void *surf, int id);
MSHSURF_API int     SurfTopAdjEdgeToElem     (void *surf, int id, int pos);

#ifdef __cplusplus
}
#endif


#endif 
