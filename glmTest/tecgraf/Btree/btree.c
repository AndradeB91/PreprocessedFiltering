/*
** btree.c
**
** Implementation of B-tree of order ORDER
**
** Waldemar Celes Filho
**
** Prof. Ruy Luiz Milidiu
**
** June 16, 1992
**
**
** To simplify the implementation it is considered that all nodes (leafs and
** nonleafs) are of the same type. There are always three pointers to 
** children nodes associate with three pointers to information field.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "btree.h"


#define IS_LEAF(n)		((n)->f[0].child == NULL)

static Tbtree *BtreeNew (void);
static Tbtree *BtreeAddLeaf (Tbtree *leaf, void *info);
static Tbtree *BtreeAddNonLeaf (Tbtree *node, Tbtree *child);
static void *BtreeIntDelete (Tbtree *root, void *info);


/* Internal variables */
static Tbtree *btree_root;
static int (*fcmp) (void *, void*);  	/* function to compare key info */

/*
** Initialise tree
*/
void BtreeInit (int (*fn) (void *, void*))
{
 fcmp = fn;
}
 
/*
** Insert function
*/
Tbtree *BtreeInsert (Tbtree *root, void *info)
{
 int i;
 Tbtree *new;
 
 if( root == NULL ) root = btree_root = BtreeNew( );
  
 if (IS_LEAF (root))	/* root is a leaf */
 {
  for (i = 0; i < ORDER && root->f[i].info != NULL ; i++)
  {
   if (fcmp (info, root->f[i].info) == 0)
    return ((Tbtree *)root);			    /* already exist */
  }
  new = BtreeAddLeaf (root, info);
  if (root == btree_root && new != NULL)    /* create first non-leaf root */
  {
   btree_root = BtreeNew ();
   btree_root->f[0].child = root, btree_root->f[0].info = root->f[0].info;
   btree_root->f[1].child = new,  btree_root->f[1].info = new->f[0].info;
   return ((Tbtree *)btree_root);
  }
  return (new == NULL ? root : new);
 }
 else  			/* root is non-leaf */
 {
  for (i = 0; i < ORDER; i++)
  { 
   if (i == ORDER-1 || 
       root->f[i+1].info == NULL || 
       fcmp (info, root->f[i+1].info) < 0)
   {
    new = BtreeInsert (root->f[i].child, info);
    root->f[i].info = root->f[i].child->f[0].info;	/* update info */
    if (new != root->f[i].child)
    {
     new = BtreeAddNonLeaf (root, new);
     if (root == btree_root && new != NULL)     /* increase depth */
     {
      btree_root = BtreeNew ();
      btree_root->f[0].child = root, btree_root->f[0].info = root->f[0].info;
      btree_root->f[1].child = new,  btree_root->f[1].info = new->f[0].info;
      return ((Tbtree *)btree_root);
     }
     return (new == NULL ? root : new);
    }
    else
     return ((Tbtree *)root);
   }
  }
 }
 return NULL;
}


/*
** Find function. Return the info pointer storage in the three if it exist or
** return NULL otherwise.
*/
void *BtreeFind (Tbtree *root, void *info)
{
 int i;
 if (root == NULL)
  return NULL;
  
 if (IS_LEAF (root))
 {
  for (i = 0; i < ORDER; i++)
   if (root->f[i].info != NULL && fcmp (info, root->f[i].info) == 0)
    return root->f[i].info;
  return NULL;
 }
 else
 {
  for (i = 1; i < ORDER; i++)
   if (root->f[i].info == NULL || fcmp (info, root->f[i].info) < 0)
   return BtreeFind (root->f[i-1].child, info);
  return BtreeFind (root->f[ORDER-1].child, info);
 }
}


/*
** Traverse tree
*/
void BtreeTraverse (Tbtree *root, void (*fn) (void *))
{
 int i;
 if (IS_LEAF(root))
  for (i = 0; i < ORDER && root->f[i].info != NULL; i++)
   fn (root->f[i].info);
 else
  for (i = 0; i < ORDER && root->f[i].child != NULL; i++)
   BtreeTraverse (root->f[i].child, fn);
}


/*
** Delete function. It returns the info pointer deleted from the tree.
** It also updates the btree root if necessary.
*/

void *BtreeDelete (Tbtree **update_root, void *info)
{
 void *del;
 del = BtreeIntDelete (*update_root, info);
 *update_root = btree_root;
 return del;
}


static void *BtreeIntDelete (Tbtree *root, void *info)
{
 int   i, j, k;
 void *del;
 
 if (root == NULL)
  return NULL;	     			/* node not found */
  
 if (IS_LEAF (root))	/* root is a leaf */
 {
  for (i = 0; root->f[i].info != NULL && i < ORDER; i++) /* find info */
   if (fcmp (info, root->f[i].info) == 0)
    break;
  if (root->f[i].info == NULL || i == ORDER)
   return NULL;	     			/* node not found */
  del = root->f[i].info;
  for (j = i; j < ORDER-1; j++)		/* make shift */
   root->f[j] = root->f[j+1];
  root->f[ORDER-1].info = NULL;
  root->f[ORDER-1].child = NULL;
  if (root == btree_root && root->f[0].info == NULL)
  {
   free (root);
   btree_root = NULL;
  }
  return del;
 }
 else  			/* root is non-leaf */
 {
  for (i = 0; i < ORDER-1; i++)			/* find child */
   if (root->f[i+1].info == NULL || fcmp (info, root->f[i+1].info) < 0)
    break;
  del = BtreeIntDelete (root->f[i].child, info);
   /* check inconsistence */  
  if (root->f[i].child->f[(ORDER-1)/2].info == NULL)
  {
   Tbtree *n = root->f[i].child;
   /* verify if left brother can give a child */
   if (i != 0 && root->f[i-1].child->f[(ORDER+1)/2].info != NULL)
   {
    Tbtree *b = root->f[i-1].child;
    for (j = (ORDER+1)/2; j < ORDER-1; j++) /* find greatest info in brother */
     if (b->f[j+1].info == NULL)
      break;
    for (k = (ORDER-1)/2; k > 0; k--)
     n->f[k] = n->f[k-1];
    n->f[0] = b->f[j];
    b->f[j].info = NULL;
    b->f[j].child = NULL;
   }
   /* verify if right brother can give a child */
   else if (i != ORDER-1 && 
            root->f[i+1].child != NULL && 
            root->f[i+1].child->f[(ORDER+1)/2].info != NULL)
   {
    Tbtree *b = root->f[i+1].child;
    n->f[(ORDER-1)/2] = b->f[0];
    for (j = 0; j < ORDER-1; j++)
     b->f[j] = b->f[j+1];
    b->f[ORDER-1].info = NULL;
    b->f[ORDER-1].child = NULL;
   }
   else		/* brothers could not give a child */
   {
    if (i != 0)
    {
     Tbtree *b = root->f[i-1].child;
     for (k = 0, j = (ORDER+1)/2; k < (ORDER-1)/2; k++, j++)
      b->f[j] = n->f[k];
    }
    else
    {
     Tbtree *b = root->f[i+1].child;
     for (j = ORDER-1; j >= (ORDER-1)/2; j--)
      b->f[j] = b->f[j-(ORDER-1)/2];
     for (j = 0; j < (ORDER-1)/2; j++)
     b->f[j] = n->f[j]; 
    }
    free (n);
    for (j = i; j < ORDER-1; j++)
     root->f[j] = root->f[j+1];
    root->f[ORDER-1].child = NULL;
    root->f[ORDER-1].info = NULL;
   }
  }
  
  if (root == btree_root && root->f[1].info == NULL) /* delete root */
  {
   btree_root = root->f[0].child;
   free (root);
   return del;
  }
  for (i = 0; root->f[i].child != NULL && i < ORDER; i++)  /* update info */
   root->f[i].info = root->f[i].child->f[0].info;
  return del;
 }
}
    


/*
** Allocate a new node
*/
static Tbtree *BtreeNew (void)
{
 int i;
 Tbtree *new = (Tbtree *)calloc (1, sizeof(Tbtree));
 if (new == NULL)
 {
  puts ("Insuficient memory");
  exit (1);
 }
 for (i = 0; i < ORDER; i++)
  new->f[i].child = NULL, new->f[i].info = NULL;
 return new;
}

/*
** Add a new data info in a leaf node. If necessary, create and return a new
** leaf. Otherwise, return NULL.
*/
static Tbtree *BtreeAddLeaf (Tbtree *leaf, void *info)
{
 int i, j, k;
 for (k = 0; k < ORDER; k++)	/* find correct position to insert */
  if (leaf->f[k].info == NULL || fcmp (info, leaf->f[k].info) < 0)
   break;
 if (leaf->f[ORDER-1].info == NULL)	/* there is enough space */
 {
  for (i = ORDER-1; i > k ; i--)
   leaf->f[i].info = leaf->f[i-1].info;
  leaf->f[k].info = info;
  return NULL;
 }
 else				/* create a new leaf */
 {
  Tbtree *new = BtreeNew ();
  if (k > (ORDER-1)/2)		/* info in the new leaf */
  {
   for (j = 0, i = (ORDER+1)/2; i < k; i++, j++)
   {
    new->f[j].info = leaf->f[i].info;
    leaf->f[i].info = NULL;
   }
   new->f[j++].info = info;
   for (i = k; i < ORDER; i++, j++)
   {
    new->f[j].info = leaf->f[i].info;
    leaf->f[i].info = NULL;
   }
  }
  else				/* info in the old leaf */
  {
   for (j = 0, i = (ORDER-1)/2; i < ORDER; i++, j++)
   {
    new->f[j].info = leaf->f[i].info;
    leaf->f[i].info = NULL;
   }
   for (i = (ORDER-1)/2; i > k ; i--)
    leaf->f[i].info = leaf->f[i-1].info;
   leaf->f[k].info = info;
  }
  return new;
 }
}


/*
** Add a new node in a non-leaf node. If necessary, create and return a new
** node. Otherwise, return NULL.
*/
static Tbtree *BtreeAddNonLeaf (Tbtree *node, Tbtree *child)
{
 int i, j, k;
 for (k = 0; k < ORDER; k++)	/* find correct position to insert */
  if (node->f[k].child==NULL || fcmp (child->f[0].info,node->f[k].info)<0)
   break;
 if (node->f[ORDER-1].child == NULL)	/* there is enough space */
 {
  for (i = ORDER-1; i > k ; i--)
   node->f[i] = node->f[i-1];
  node->f[k].info = child->f[0].info;
  node->f[k].child = child;
  return NULL;
 }
 else				/* create a new node */
 {
  Tbtree *new = BtreeNew ();
  if (k > ORDER/2)		/* info in the new node */
  {
   for (j = 0, i = (ORDER+1)/2; i < k; i++, j++)
   {
    new->f[j] = node->f[i];
    node->f[i].child = NULL, node->f[i].info = NULL;
   }
   new->f[j].info = child->f[0].info;
   new->f[j++].child = child;
   for (i = k; i < ORDER; i++, j++)
   {
    new->f[j] = node->f[i];
    node->f[i].child = NULL, node->f[i].info = NULL;
   }
  }
  else				/* info in the old node */
  {
   for (j = 0, i = ORDER/2; i < ORDER; i++, j++)
   {
    new->f[j] = node->f[i];
    node->f[i].child = NULL, node->f[i].info = NULL;
   }
   for (i = ORDER/2; i > k ; i--)
    node->f[i] = node->f[i-1];
   node->f[k].info = child->f[0].info;
   node->f[k].child = child;
  }
  return new;
 }
}


