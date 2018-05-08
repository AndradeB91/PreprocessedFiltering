/*
** btree.h
**
** Implementation of B-tree of order ORDER
**
** Waldemar Celes Filho
**
** Prof. Ruy Luiz Milidiu
**
** June 16, 1992
*/


#ifndef btree_h
#define btree_h

#include "deflibMshaux.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ORDER	3

typedef struct btree_ Tbtree;
typedef struct field_ Tfield;

struct field_
{
 Tbtree *child;
 void  *info;
};
struct btree_
{
 Tfield f[ORDER];
};

MSHAUX_API void    BtreeInit (int (*fn) (void *, void*));
MSHAUX_API Tbtree *BtreeInsert (Tbtree *root, void *info);
MSHAUX_API void   *BtreeFind (Tbtree *root, void *info);
MSHAUX_API void   *BtreeDelete (Tbtree **update_root, void *info);
MSHAUX_API void    BtreeTraverse (Tbtree *root, void (*fn) (void *));


#ifdef __cplusplus
}
#endif

#endif
