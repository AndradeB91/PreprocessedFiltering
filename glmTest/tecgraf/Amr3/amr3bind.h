/*
* amr3bind.h
* Binding [C] for amr3bind.c [C++] implementation
*/

#ifndef rtree_h
#define rtree_h

#include "deflibRTree.h"

#ifdef __cplusplus
extern "C" {
#endif

RTREE_API void *RtreeCreate(void);
RTREE_API void RtreeDestroy(void *r);
RTREE_API void RtreeInsert(void *r,
							      void *info,
									double xmin,double xmax,
									double ymin,double ymax,
									double zmin,double zmax);
RTREE_API void RtreeDelete(void *r,
							      void *info,
									double xmin,double xmax,
									double ymin,double ymax,
									double zmin,double zmax);
RTREE_API void RtreeInitTraverse(void *r);
RTREE_API void *RtreeTraverse(void *r,
							         double *xmin,double *xmax,
										double *ymin,double *ymax,
										double *zmin,double *zmax);
RTREE_API void RtreeInitSearchBox(void *r,
							             double xmin,double xmax,
											 double ymin,double ymax,
											 double zmin,double zmax);
RTREE_API void *RtreeSearchBox(void *r,
							          double *xmin,double *xmax,
										 double *ymin,double *ymax,
										 double *zmin,double *zmax);

RTREE_API int RtreeNumber(void *r);

#ifdef __cplusplus
}
#endif

#endif
