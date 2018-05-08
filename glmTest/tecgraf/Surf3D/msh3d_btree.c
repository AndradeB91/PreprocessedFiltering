
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>

/* #include "new.h" */
#include "msh3d_btree.h"

/* static Ferror fn_fatal; */        /* function to handle fatal error */

#define IS_LEAF(n)		((n)->f[0].child == NULL)

static void *Msh3dBtreeNew (void);
static Tbtree *Msh3dBtreeAddLeaf (Tbtree *leaf, void *info);
static Tbtree *Msh3dBtreeAddNonLeaf (Tbtree *node, Tbtree *child);
static void *Msh3dBtreeIntDelete (Tbtree *root, void *info);



/* Internal variables */
static Tbtree *btree_root;
static int (*fcmp) (void *, void*);  	/* function to compare key info */

/*
** Initialise tree
*/
void Msh3dBTreeInit (int (*fn) (void *, void*)/* ,Ferror f */)
{
 fcmp = fn;
/*  fn_fatal = f; */
}
 
/*
** Insert function
*/
Tbtree *Msh3dBtreeInsert (Tbtree *root, void *info)
{
 int i;
 Tbtree *new;
 
 if (root == NULL)
  root = btree_root = Msh3dBtreeNew ();
  
 if (IS_LEAF (root))	/* root is a leaf */
 {
  for (i = 0; i < ORDER && root->f[i].info != NULL; i++)
  {
   if (fcmp (info, root->f[i].info) == 0)
    return root;			    /* already exist */
  }
  new = Msh3dBtreeAddLeaf (root, info);
  if (root == btree_root && new != NULL)    /* create first non-leaf root */
  {
   btree_root = Msh3dBtreeNew ();
   btree_root->f[0].child = root, btree_root->f[0].info = root->f[0].info;
   btree_root->f[1].child = new,  btree_root->f[1].info = new->f[0].info;
   return btree_root;
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
    new = Msh3dBtreeInsert (root->f[i].child, info);
    root->f[i].info = root->f[i].child->f[0].info;	/* update info */
    if (new != root->f[i].child)
    {
     new = Msh3dBtreeAddNonLeaf (root, new);
     if (root == btree_root && new != NULL)     /* increase depth */
     {
      btree_root = Msh3dBtreeNew ();
      btree_root->f[0].child = root, btree_root->f[0].info = root->f[0].info;
      btree_root->f[1].child = new,  btree_root->f[1].info = new->f[0].info;
      return btree_root;
     }
     return (new == NULL ? root : new);
    }
    else
     return root;
   }
  }
 }
 return NULL;
}


/*
** Find function. Return the info pointer storage in the three if it exist or
** return NULL otherwise.
*/
void *Msh3dBtreeFind (Tbtree *root, void *info)
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
   return Msh3dBtreeFind (root->f[i-1].child, info);
  return Msh3dBtreeFind (root->f[ORDER-1].child, info);
 }
}


/*
** Traverse tree
*/
void Msh3dBtreeTraverse (Tbtree *root, void (*fn) (void *))
{
 int i;
 if (IS_LEAF(root))
  for (i = 0; i < ORDER && root->f[i].info != NULL; i++)
   fn (root->f[i].info);
 else
  for (i = 0; i < ORDER && root->f[i].child != NULL; i++)
   Msh3dBtreeTraverse (root->f[i].child, fn);
}


/*
** Delete function. It returns the info pointer deleted from the tree.
** It also updates the btree root if necessary.
*/

void *Msh3dBtreeDelete (Tbtree **update_root, void *info)
{
 void *del;
 del = Msh3dBtreeIntDelete (*update_root, info);
 *update_root = btree_root;
 return del;
}


static void *Msh3dBtreeIntDelete (Tbtree *root, void *info)
{
 int   i, j, k;
 void *del;
 
 if (root == NULL)
  return NULL;	     			/* node not found */
  
 if (IS_LEAF (root))	/* root is a leaf */
 {
  for (i = 0; i < ORDER && root->f[i].info != NULL; i++) /* find info */
   if (fcmp (info, root->f[i].info) == 0)
    break;
  if (i == ORDER || root->f[i].info == NULL)
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
  del = Msh3dBtreeIntDelete (root->f[i].child, info);
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
  for (i = 0; i < ORDER && root->f[i].child != NULL; i++)  /* update info */
   root->f[i].info = root->f[i].child->f[0].info;
  return del;
 }
}
    


/*
** Allocate a new node
*/
static void *Msh3dBtreeNew (void)
{
 int i;
 Tbtree *new =(Tbtree *)malloc(sizeof(Tbtree));
 if (new == NULL) return NULL; /* fn_fatal("Insuficient memory"); */
 for (i = 0; i < ORDER; i++)
  new->f[i].child = NULL, new->f[i].info = NULL;
 return new;
}

/*
** Add a new data info in a leaf node. If necessary, create and return a new
** leaf. Otherwise, return NULL.
*/
static Tbtree *Msh3dBtreeAddLeaf (Tbtree *leaf, void *info)
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
  Tbtree *new = Msh3dBtreeNew ();
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
static Tbtree *Msh3dBtreeAddNonLeaf (Tbtree *node, Tbtree *child)
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
  Tbtree *new = Msh3dBtreeNew ();
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


/*---------------------------------------------------------*/

static void BTSet(Msh3dBtree *bt)
{
 fcmp=bt->fcmp;
/*  fn_fatal=bt->ferr; */
 btree_root=bt->root;
}
 
Msh3dBtree *Msh3dBTreeCreate(Msh3dBtreeComp fc,Msh3dBtreeError fe)
{
 Msh3dBtree *bt=(Msh3dBtree *)malloc(sizeof(Msh3dBtree));
 bt->fcmp = fc;
 bt->ferr = fe;
 bt->root = NULL;
 BTSet(bt);
 return bt;
}

void Msh3dBTreeRelease(Msh3dBtree *bt)
{
 if (bt==NULL) return;
 BTSet(bt);
 while (bt->root!=NULL)
  Msh3dBtreeDelete(&(bt->root),bt->root->f[0].info);
 free(bt);
/*  fn_fatal=NULL; */
 btree_root=NULL;
 fcmp=NULL; 
}

void Msh3dBTreeDelete(Msh3dBtree *bt,void *info)
{
 if (bt==NULL||info==NULL) return;
 BTSet(bt);
 Msh3dBtreeDelete(&(bt->root),info);
}

void Msh3dBTreeInsert(Msh3dBtree *bt,void *info)
{
 BTSet(bt);
 bt->root=Msh3dBtreeInsert(bt->root,info);
}

void *Msh3dBTreeFind(Msh3dBtree *bt,void *info)
{
 BTSet(bt);
 return Msh3dBtreeFind(bt->root,info);
}

void Msh3dBTreeTraverse(Msh3dBtree *bt,void (*fn) (void *))
{
 BTSet(bt);
 if (bt->root==NULL) return;
 Msh3dBtreeTraverse(bt->root,fn);
}
