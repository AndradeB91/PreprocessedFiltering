#ifndef msh3d_btree_h
#define msh3d_btree_h

#include "../Deflib/deflib3d.h"

//#include "..\\..\\mshaux\\include\\btree.h"				// mshaux library
//#include "btree.h"
#include "../Btree/btree.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef MSH3D_API struct btree  Msh3dBtree;

typedef MSH3D_API void (*Msh3dBtreeError) (char *);
typedef MSH3D_API int  (*Msh3dBtreeComp)  (void *,void *);

MSH3D_API struct btree
{
 Msh3dBtreeError	ferr;	/* function to handle fatal error */
 Msh3dBtreeComp	fcmp;	/* function to compare key info */
 Tbtree		*root;	/* pointer to tree's root */
};

MSH3D_API Msh3dBtree *Msh3dBTreeCreate	(Msh3dBtreeComp fc,Msh3dBtreeError fe);
MSH3D_API void   Msh3dBTreeInsert	(Msh3dBtree *bt,void *info);
MSH3D_API void  *Msh3dBTreeFind	(Msh3dBtree *bt,void *info);
MSH3D_API void   Msh3dBTreeTraverse	(Msh3dBtree *bt,void (*fn) (void *));
MSH3D_API void   Msh3dBTreeRelease	(Msh3dBtree *bt);
MSH3D_API void   Msh3dBTreeDelete	(Msh3dBtree *bt,void *info);

#ifdef __cplusplus
}
#endif

#endif	// msh3d_btree_h
