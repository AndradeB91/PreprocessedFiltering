//
// OrderedSet Template Class header file
//
// Description -
//   This class implements a set abstract data type.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.6 $  $Date: 2001/01/08 20:28:29 $  $Author: wash $
//

#ifndef OrderedSet_hh
#define OrderedSet_hh

namespace FTools {

template<class KeyType>
class OrderedSet {

    public:

        enum Dir { LEFT, RIGHT } ;

        enum Code { NO_FLIP, FLIP_GCH, FLIP_CH } ;

        class SetIterator ;

        OrderedSet(bool allow_duplicates = false) ;
        virtual ~OrderedSet() ;

        virtual int Compare(const KeyType& op0,const KeyType& op1) const = 0 ;

        bool Insert(const KeyType &key) ;
        bool Remove(const KeyType &key) ;
        bool Del(const KeyType &key) { return(Remove(key)) ; }
        bool HasKey(const KeyType &key) const ;

        KeyType *GetSmallest() const ;
        KeyType *GetLargest() const ;

        int Len() { return(NumSetEntries) ; }

        SetIterator Iterator() {
            KeyType** tmp = GetKeyList() ;
            return(SetIterator(NumSetEntries,tmp)) ;
        }

// -------DEBUG----------

#ifdef DEBUG_RBTREE
        void xCheck() ;
        void xWalkRBTree(char *name,char *mode) ;
        void Dump() ;
#endif

    private:

        struct ArbSetMember {
            KeyType      key ;
            bool         red ;
            ArbSetMember *left ;
            ArbSetMember *right ;
        } ;

        ArbSetMember *Root ;
        bool DuplicatesOK ;
        int NumSetEntries ;

        KeyType **GetKeyList() const ;

        ArbSetMember *NewMember() ;
        void DeleteRec(ArbSetMember *node) ;
        ArbSetMember *Rotate(KeyType const &key,ArbSetMember *r,
                             int flip_mode,int dir,int gdir) ;
        void Split(ArbSetMember *n,ArbSetMember **c,ArbSetMember **p,
                   ArbSetMember *g,ArbSetMember *gg) ;

        void WalkRec(ArbSetMember *n,KeyType **keys,
                     int *cur) const ;


// -------DEBUG----------

#ifdef DEBUG_RBTREE
        void xrWalk(ArbSetMember *n,int level) ;
        void xrCheck(ArbSetMember *c) ;
        void DumpRec(ArbSetMember *node) ;
#endif
// ----------------------

    public:

        class SetIterator {
            public:
                SetIterator() : cur(0),num(0),data(0) {} ;
                SetIterator(int inum,KeyType** idata) : cur(0),num(inum),data(idata) {}
                ~SetIterator() { delete [] data ; }


                void First() { cur = 0 ; }
                void Next()  { cur++ ; }
                bool More()  { return(cur < num) ; }

                KeyType &Entry()    { return(*data[cur]) ; }
                KeyType *EntryPtr() { return(data[cur]) ; }
                int Len()           { return(num) ; }

                void operator ++ ()      { Next() ; } ;
                void operator ++ (int i) { Next() ; } ;

                KeyType &  operator * () { return(Entry()) ; } ;
 
            private:
                int cur, num ;
                KeyType** data ;
        } ;



} ;

inline int ArbCmpUnsigned(const int &u0,const int &u1)
{
    if (u0 > u1) return(1) ;
    if (u0 < u1) return(-1) ;
    return(0) ;
}


// %(OrderedSet::OrderedSet-constructor-|-int-|(*^)()) 
/* ++ ----------------------------------------------------------
**
**    OrderedSet - constructor 
**
**      OrderedSet(int (*cmp_func)())
**
**        cmp_func - (in)  comparison function 
**
**      Description: This is a constructor for an ArbSet object. 
**
**
** -- */

template<class KeyType>
OrderedSet<KeyType>::OrderedSet(bool allow_duplicates)
{
    Root = NewMember() ;
    DuplicatesOK = allow_duplicates ;
    NumSetEntries = 0 ;
}




// %(OrderedSet::OrderedSet-destructor-virtual|~) 
/* ++ ----------------------------------------------------------
**
**    OrderedSet - destructor 
**
**      ~OrderedSet()
**
**      Description: This is a destructor for an ArbSet object. 
**
**
** -- */

template<class KeyType>
OrderedSet<KeyType>::~OrderedSet()
{
    DeleteRec(Root) ;
}

#if DEBUG_RBTREE
template<class KeyType>
void OrderedSet<KeyType>::Dump()
{
    DumpRec(Root) ;
}

static void dump_one(int i_val)
{
    fprintf(stderr,"%d\n",i_val) ;
}

static void dump_one(ArbIntEdge edge)
{
    fprintf(stderr,"%d %d\n",edge.node_id[0],edge.node_id[1]) ;
}

static void dump_one(ArbMshElement2D elem)
{
    fprintf(stderr,"%d\n",elem.elem_id) ;
}

static void dump_one(double d_val)
{
    fprintf(stderr,"%g\n",d_val) ;
}

template<class KeyType>
void OrderedSet<KeyType>::DumpRec(ArbSetMember *node)
{
    if (node->left != 0) DumpRec(node->left) ;
    if (node->right != 0) DumpRec(node->right) ;
    if ((node->left == 0) && (node->right == 0))
        dump_one(node->key) ;
}
#endif




// %(OrderedSet::DeleteRec-void-|-ArbSetMember-|*) 
/* ++ ----------------------------------------------------------
**
**    DeleteRec - delete all nodes 
**
**      void DeleteRec(ArbSetMember *node)
**
**        node - (in)  current node 
**
**      Description: This function recursively deletes the current node 
**          and all children. 
**
**
** -- */

template<class KeyType>
void OrderedSet<KeyType>::DeleteRec(ArbSetMember *node)
{
    if (node->left != 0) DeleteRec(node->left) ;
    if (node->right != 0) DeleteRec(node->right) ;
//    if ((node->left == 0) && (node->right == 0))
//        dump_one(node->key) ;
    delete node ;
}




// %(OrderedSet::NewMember-ArbSetMember-|*)
/* ++ ----------------------------------------------------------
**
**    NewMember - create a new node 
**
**      ArbSetMember *NewMember()
**
**      Description: This function creates a new binary tree node. 
**
**      Return Value: a new tree node 
**
**
** -- */

template<class KeyType>
typename OrderedSet<KeyType>::ArbSetMember *OrderedSet<KeyType>::NewMember()
{
    ArbSetMember *node = new ArbSetMember ;
    node->left = node->right = 0 ;
    node->red = false ;
    return(node) ;
}




// %(OrderedSet::HasKey-bool-|^const-KeyType-const|&) 
/* ++ ----------------------------------------------------------
**
**    HasKey - check for an entry 
**
**      bool HasKey(const KeyType &key) const
**
**        key - (in)  entry 
**
**      Description: This function returns true if the specified entry 
**          is in the set. 
**
**      Return Value: true if the entry is found in the set 
**
**
** -- */

template<class KeyType>
bool OrderedSet<KeyType>::HasKey(const KeyType &key) const
{
    ArbSetMember *s ;
    int dir ;

    s = Root->right ;
    while (s != NULL) {
        dir = Compare(key,s->key) ;
        if ((dir == 0) && (s->right == NULL)) return(true) ;
        if (dir < 0)
            s = s->left ;
        else
            s = s->right ;
    }
    return(false) ;
}




// %(OrderedSet::GetSmallest-KeyType-|*^const)
/* ++ ----------------------------------------------------------
**
**    GetSmallest - get the smallest entry 
**
**      KeyType *GetSmallest() const
**
**      Description: This function returns a pointer to the entry that 
**          has been ranked the smallest by the comparison function. 
**
**      Return Value: a pointer to the smallest entry 
**
**
** -- */

template<class KeyType>KeyType *OrderedSet<KeyType>::GetSmallest() const
{
    ArbSetMember *p, *s ;

    p = Root ;
    s = p->right ;
    while (s->left != NULL) {
        s = s->left ;
    }
    return(&(s->key)) ;
}




// %(OrderedSet::GetLargest-KeyType-|*^const) 
/* ++ ----------------------------------------------------------
**
**    GetLargest - get the largest entry 
**
**      KeyType *GetLargest() const
**
**      Description: This function returns a pointer to the entry that 
**          has been ranked the largest by the comparison function. 
**
**      Return Value: a pointer to the largest entry 
**
**
** -- */

template<class KeyType>KeyType *OrderedSet<KeyType>::GetLargest() const
{
    ArbSetMember *p, *s ;

    p = Root ;
    s = p->right ;
    while (s->right != NULL) {
        s = s->right ;
    }
    return(&(s->key)) ;
}




// %(OrderedSet::Insert-bool-|-KeyType-const|&)
/* ++ ----------------------------------------------------------
**
**    Insert - insert a new entry 
**
**      bool Insert(const KeyType &key)
**
**        key - (in)  entry to insert 
**
**      Description: This function inserts an entry int the set. 
**
**      Return Value: true if the entry was added, false if it was 
**          already in the set 
**
**
** -- */

template<class KeyType>
bool OrderedSet<KeyType>::Insert(const KeyType &key)
{
    int p_dir ;
    ArbSetMember *p, *s, *n ;
    ArbSetMember *g = NULL ;
    ArbSetMember *gg = NULL ;

    // make a new node containing this key

    n = NewMember() ;
    n->key = key ;

    // search until we find a leaf

    p = Root ;
    p_dir = RIGHT ;      // direction from p to s
    s = p->right ;

    if (s != 0) {
        ArbSetMember *temp ;
        int dir ;

        // look for a leaf, splitting nodes on the way down

        while (s->right != NULL) {
            if (s->left->red && s->right->red) Split(n,&s,&p,g,gg) ;
            gg = g ;
            g = p ;
            p = s ;
            if (Compare(key,s->key) < 0) {
                s = s->left ;
                p_dir = LEFT ;
            } else {
                s = s->right ;
                p_dir = RIGHT ;
            }
        }

        dir = Compare(n->key,s->key) ;
        if (!DuplicatesOK && dir == 0) {
            delete n ;
            return(false) ;
        }

        // must replace s with a new internal node that has as
        // its children s and n.  The new node gets the larger
        // of s and n as its key.  The new node gets painted red,
        // its children are black.  Coloring is done by split().

        temp = NewMember() ;
        if (dir < 0) {
            temp->key = s->key ;
            temp->red = s->red ;
            temp->left = n ;
            temp->right = s ;
        } else {
            temp->key = n->key ;
            temp->red = n->red ;
            temp->left = s ;
            temp->right = n ;
        }
        n = temp ;
    }

    // add the new node

    if (p_dir == RIGHT)
        p->right = n ;
    else
        p->left = n ;

    // color this node red and check red-black balance

    Split(n,&n,&p,g,gg) ;

    NumSetEntries++ ;

#ifdef DEBUG_RBTREE
    xCheck() ;
#endif

    return(true) ;
}




// %(OrderedSet::Remove-bool-|-KeyType-const|&)
/* ++ ----------------------------------------------------------
**
**    Remove - remove an entry 
**
**      bool Remove(const KeyType &key)
**
**        key - (in)  entry to remove 
**
**      Description: This function returns an entry from the set. 
**
**      Return Value: true if the entry was found in the set, false 
**          otherwise 
**
**
** -- */

template<class KeyType>
bool OrderedSet<KeyType>::Remove(const KeyType &key)
{
    // the goal is to arrive at a leaf with a red parent.
    // We force this by dragging a red node with us down
    // the tree, re-arranging the tree to keep the balance
    // as we go.  All the reaffangements keep the tree
    // balanced, so if we cancel the deletion or don't find
    // the specified node to delete, we can just quit.

    ArbSetMember *s, *p, *g ;
    ArbSetMember *next_node, *next_other ;
    int dir, next_dir ;

    g = NULL ;
    p = Root ;
    s = p->right ;
    dir = RIGHT ;

    // First check on the root.  It must exist, have children,
    // and either it or one of it's children must be red.  We
    // can just paint the root red, if necessary, as this will
    // affect the black height of the entire tree equally.

    if (s == NULL) return(false) ;

    // check to make sure the root isn't an only child

    if (s->left == NULL) {
        if (Compare(key,s->key) == 0) {
            p->right = NULL ;  // deleting the root
            delete s ;
            return(true) ;
        } else
            return(false) ;
    }

    // Now, either the root or one of its children must be red

    if (!s->left->red && !s->right->red)
        s->red = true ;  // just color the root red

    // March down the tree, always working to make sure the
    // current node is red.  That way, when we do arrive at a
    // leaf, its parent will be red, making the leaf very easy
    // to delete (just drop the leaf, and replace its (red)
    // parent with its (black) sibling.

    while (1) {

        // If we're at a leaf, we're done.

        if (s->left == NULL) break ;

        // decide where we are going next

        if (Compare(key,s->key) < 0) {
            next_node = s->left ;
            next_other = s->right ;
            next_dir = LEFT ;
        } else {
            next_node = s->right ;
            next_other = s->left ;
            next_dir = RIGHT ;
        }

        // If the current node or the next node is red, we
        // can advance.

        if (s->red || next_node->red) ; // do nothing

        // If the current node is black and the next node
        // is black, but the next node's sibling is red,
        // then rotate from parent towards the red child.  This
        // will lower the curent node, and give us a new
        // grandparent (the old parent) and a new parent (the
        // sibling that was red).  We then paint the current
        // node red and the new parent is painted black.

        else if (next_other->red) {
            g = p ;
            int ch_dir = (p->left == s) ? LEFT : RIGHT ;
            int gch_dir = (s->left == next_other) ? LEFT : RIGHT ;
            p = Rotate(next_node->key,p,FLIP_GCH,ch_dir,gch_dir) ;
            s->red = true ;
            p->red = false ;

        // If the current node is black and its left child is
        // black and its right child is black, then (a) the current
        // node's parent must be red (we never advance unless we
        // are leaving a red node), (b) it sibling must be black
        // (because the parent is red), and (c) we need to color
        // the current node red, the parent black, and then check
        // for tree imbalances.

        } else {
            ArbSetMember *sib ;
            //assert(p->red == true) ;
            if (dir == LEFT)
                sib = p->right ;
            else
                sib = p->left ;

            //assert(sib->red == false) ;
            //assert(sib->left != NULL) ;

            s->red = true ;
            p->red = false ;

            // First case, black sibling has two black kids.  Just
            // color the sibling red.  In effect we are reversing
            // a simple color flip.

            if (!sib->left->red && !sib->right->red)
                sib->red = true ;

            // second case, black sibling has at least on red child.
            // (it makes nod differenc if both kids are red).  We
            // need to do either a single or double rotation in
            // order to rebalance the tree.

            else {
                int redkid_dir ;

                if (sib->left->red)
                    redkid_dir = LEFT ;
                else
                    redkid_dir = RIGHT ;

                if (!dir == redkid_dir) {
                    sib->red = true ;
                    if (redkid_dir == LEFT)
                        sib->left->red = false ;
                    else
                        sib->right->red = false ;
                    g = Rotate(key,g,FLIP_GCH,-1,-1) ;
                } else {
                    Rotate(key,p,FLIP_CH+redkid_dir,-1,-1) ;
                    g = Rotate(key,g,FLIP_GCH,-1,-1) ;
                }
            }
        }

        // advance pointers

        dir = next_dir ;
        g = p ;
        p = s ;
        s = (dir == RIGHT) ? s->right : s->left ;
    }

    // make the root black

    Root->right->red = false ;

    // Delete it, if a match.  Parent is red

    if (Compare(s->key,key) == 0) {
        //assert((p->red == true) || (p == Root->right)) ;
        if ((Compare(s->key,g->key) < 0) && (g != Root)) {
            if (Compare(s->key,p->key) < 0)
                g->left = p->right ;
            else
                g->left = p->left ;
        } else {
            if (Compare(s->key,p->key) < 0)
                g->right = p->right ;
            else
                g->right = p->left ;
        }

        delete p ;
        delete s ;
        NumSetEntries-- ;

#ifdef DEBUG_RBTREE
        xCheck() ;
#endif
        return(true) ;
    } else {

#ifdef DEBUG_RBTREE
        xCheck() ;
#endif
        return(false) ;
    }
}




// %(OrderedSet::Rotate-ArbSetMember-|*-KeyType-|const&-ArbSetMember-|*-int-|-int-|-int-|)
/* ++ ----------------------------------------------------------
**
**    Rotate - do a tree rotate 
**
**      ArbSetMember *Rotate(
**              KeyType      const&key,
**              ArbSetMember *r,
**              int          flip_mode,
**              int          dir,
**              int          gdir)
**
**        key       - (in)  entry being added 
**        r         - (in)  local root 
**        flip_mode - (in)  rotation mode 
**        dir       - (in)  parent direction 
**        gdir      - (in)  grandparent direction 
**
**      Description: This function do a rotate on the binary tree. 
**
**      Return Value: new local tree root 
**
**
** -- */

template<class KeyType>
typename OrderedSet<KeyType>::ArbSetMember
*OrderedSet<KeyType>::Rotate(KeyType const &key,
                          ArbSetMember *r,int flip_mode,int dir,int gdir)
{
/*
    Rotate child and grandchild of r along the path specified
    by searching for n.  For example, if n was equal to 3,
    gc2, or 4, the following rotation occures:

               r                           r
               |                           |
               c                          gc2
            /     \         ==>          /   \
          gc1     gc2                   c     4
         /   \   /   \                 / \
        1     2 3     4              gc1  3
                                    /   \
                                   1     2

    as r may connect to c via either its left or right link,
    there are actually four symmetric variants.

    A pointe to the top of the new rotated nodes (in the case
    above, to gc2) is returned.

    This routine is complicated by the fact that the routine
    uses the value of the node n to decide which direction
    to rotate.  This may or may not be the direction the
    caller has in mind.  Rather than require the caller to specify
    the direction of rotation, it is easier to allow the
    caller to specify whether to go in the direction of n
    or away from it.  This is done by the last argument to the
    function, flip_mode.  The caller can indicate that
    either or both of the directions to child and grandchild
    should be reversed during the rotation.
*/

    ArbSetMember *ch, *gch ;
    int ch_dir, gch_dir ;

    // identify child and grandchild

    if (dir < 0) {
        if (r == Root)
            ch_dir = RIGHT ;
        else
            ch_dir = (Compare(key,r->key) < 0) ? LEFT : RIGHT ;
        if (flip_mode & FLIP_CH) ch_dir = !ch_dir ;
    } else
        ch_dir = dir ;

    if (r == Root) ch_dir = RIGHT ;  // special case
    ch = (ch_dir == RIGHT) ? r->right : r->left ;

    if (gdir < 0) {
        gch_dir = (Compare(key,ch->key) < 0) ? LEFT : RIGHT ;
        if (flip_mode) {
            if (flip_mode == FLIP_GCH)
                gch_dir = !gch_dir ;
            else
                gch_dir = flip_mode & 1 ;
        }
    } else
        gch_dir = gdir ;

    gch = (gch_dir == RIGHT) ? ch->right : ch->left ;

    // rotate: now move pointers

    if (gch_dir == RIGHT) {
        ch->right = gch->left ;
        gch->left = ch ;
    } else {
        ch->left = gch->right ;
        gch->right = ch ;
    }
    if (ch_dir == RIGHT)
        r->right = gch ;
    else
        r->left = gch ;

    return(gch) ;
}




// %(OrderedSet::Split-void-|-ArbSetMember-|*-ArbSetMember-|**-ArbSetMember-|**-ArbSetMember-|*-ArbSetMember-|*)
/* ++ ----------------------------------------------------------
**
**    Split - split a red/black tree 
**
**      void Split(
**              ArbSetMember *n,
**              ArbSetMember **c,
**              ArbSetMember **p,
**              ArbSetMember *g,
**              ArbSetMember *gg)
**
**        n  - (in)  new entry 
**        c  - (i/o) child 
**        p  - (i/o) parent 
**        g  - (in)  grandparent 
**        gg - (in)  great grandparent 
**
**      Description: This function takes care of tree colors and 
**          balances the tree. 
**
**
** -- */

template<class KeyType>
void OrderedSet<KeyType>::Split(ArbSetMember *n,ArbSetMember **c,
                             ArbSetMember **p,ArbSetMember *g,
                             ArbSetMember *gg)
{
/*
    Take care of colors and balance.  It will color the current
    location red, the current location's children black, and
    then look to see if two consecutive red nodes have been
    created.  If so, a singe or double rotation will be done
    to fix the tree.
*/

    if (Root->red) {
        fputs("dummyhead was read!!\n",stdout) ;
        Root->red = false ;
    }

    (*c)->red = true ;
    if ((*c)->left != 0)  (*c)->left->red = false ;
    if ((*c)->right != 0) (*c)->right->red = false ;

    // check to make sure we have not created two red links
    // in a row.  If so, we must rotate.

    if ((*p)->red) {
//        if (gg->red != true) g->red = true ;
        g->red = true ;

        // if the red links don't point in the same direction,
        // then we will need to do a double rotation.  The lower
        // half is around the grandparent and then the upper half
        // is around the great-grandparent.

        if (g->right == *p) {
            if ((*p)->right != *c) { 
                *p = Rotate(n->key,g,NO_FLIP,RIGHT,LEFT) ;
            }
        } else {
            if ((*p)->left != *c) { 
                *p = Rotate(n->key,g,NO_FLIP,LEFT,RIGHT) ;
            }
        }

        // same for both single and double rotations

        int ch_dir = (gg->left == g) ? LEFT : RIGHT ;
        int gch_dir = (g->left == *p) ? LEFT : RIGHT ;
        if (g != Root) *c = Rotate(n->key,gg,NO_FLIP,ch_dir,gch_dir) ;
        (*c)->red = false ;
    }
    Root->right->red = false ;
}


template<class KeyType>
void OrderedSet<KeyType>::WalkRec(ArbSetMember *n,KeyType **keys,
                               int *cur) const
{
    if (n != NULL) {
        WalkRec(n->left,keys,cur) ;
        if ((n->left == 0) && (n->right == 0)) {
            keys[*cur] = &(n->key) ;
            ++(*cur) ;
        }
        WalkRec(n->right,keys,cur) ;
    }
}




// %(OrderedSet::GetKeyList-KeyType-|**^const)
/* ++ ----------------------------------------------------------
**
**    GetKeyList - get a list of entries 
**
**      KeyType **GetKeyList() const
**
**      Description: This function returns a list of pointers to the 
**          set entries. Ownership of the list memory passes to the 
**          client, which much eventually call delete []. 
**
**      Return Value: list of pointers to entries 
**
**
** -- */

template<class KeyType>
KeyType **OrderedSet<KeyType>::GetKeyList() const
{
    if (NumSetEntries == 0) return(0) ;

    KeyType **keys = new KeyType*[NumSetEntries] ;
    int cur = 0 ;

    WalkRec(Root->right,keys,&cur) ;
    return(keys) ;
}

// --------------------------------------------

#ifdef DEBUG_RBTREE
template<class KeyType> void OrderedSet<KeyType>::xCheck()
{
    xrCheck(Root->right) ;
}

template<class KeyType> void OrderedSet<KeyType>::xrCheck(ArbSetMember *c)
{
    if (c->left == 0) {
        //assert(c->right == 0) ;
        return ;
    }
    if (c->right == 0) {
        //assert(c->left == 0) ;
        return ;
    }
    //if (c->red) {
        //assert(c->left->red == false) ;
        //assert(c->right->red == false) ;
    //}
    if (c->left != 0) xrCheck(c->left) ;
    if (c->right != 0) xrCheck(c->right) ;
}
#endif

/*
CLASS OrderedSet<class KeyType>

  This object implements a templated set data structure. A copy is made 
  of an entry when it is placed in the set. The client must provide a 
  comparison function that returns an interger less then, equal to, or 
  greater than zero depending on the relative value of it's arguments. 
  The set is implemented as a red/black binary tree. 


  Template Arguments:

    KeyType - set entry type


PUBLIC INTERFACE

  Public Member Functions:

    OrderedSet - constructor 

      OrderedSet(int (*cmp_func)())

        cmp_func - (in)  comparison function 

      Description: This is a constructor for an ArbSet object. 


    OrderedSet - destructor 

      ~OrderedSet()

      Description: This is a destructor for an ArbSet object. 


    Insert - insert a new entry 

      bool Insert(const KeyType &key)

        key - (in)  entry to insert 

      Description: This function inserts an entry int the set. 

      Return Value: true if the entry was added, false if it was 
          already in the set 


    Remove - remove an entry 

      bool Remove(const KeyType &key)

        key - (in)  entry to remove 

      Description: This function returns an entry from the set. 

      Return Value: true if the entry was found in the set, false 
          otherwise 


    HasKey - check for an entry 

      bool HasKey(const KeyType &key) const

        key - (in)  entry 

      Description: This function returns true if the specified entry 
          is in the set. 

      Return Value: true if the entry is found in the set 


    GetKeyList - get a list of entries 

      KeyType **GetKeyList() const

      Description: This function returns a list of pointers to the 
          set entries. Ownership of the list memory passes to the 
          client, which much eventually call delete []. 

      Return Value: list of pointers to entries 


    GetSmallest - get the smallest entry 

      KeyType *GetSmallest() const

      Description: This function returns a pointer to the entry that 
          has been ranked the smallest by the comparison function. 

      Return Value: a pointer to the smallest entry 


    GetLargest - get the largest entry 

      KeyType *GetLargest() const

      Description: This function returns a pointer to the entry that 
          has been ranked the largest by the comparison function. 

      Return Value: a pointer to the largest entry 


    NumEntries - number of entries 

      int NumEntries()

      Description: This function returns the number entries in the 
          set. 

      Return Value: number of set entries 


PRIVATE INTERFACE

  Private Data Structures:

    struct ArbSetMember

      This structure stores a set entry 

      Member Variables:

        KeyType key - entry data 

        bool red - flag for red/black tree 

        ArbSetMember *left - left child 

        ArbSetMember *right - right child 


  Private Member Functions:

    NewMember - create a new node 

      ArbSetMember *NewMember()

      Description: This function creates a new binary tree node. 

      Return Value: a new tree node 


    DeleteRec - delete all nodes 

      void DeleteRec(ArbSetMember *node)

        node - (in)  current node 

      Description: This function recursively deletes the current node 
          and all children. 


    Rotate - do a tree rotate 

      ArbSetMember *Rotate(
              KeyType      const&key,
              ArbSetMember *r,
              int          flip_mode,
              int          dir,
              int          gdir)

        key       - (in)  entry being added 
        r         - (in)  local root 
        flip_mode - (in)  rotation mode 
        dir       - (in)  parent direction 
        gdir      - (in)  grandparent direction 

      Description: This function do a rotate on the binary tree. 

      Return Value: new local tree root 


    Split - split a red/black tree 

      void Split(
              ArbSetMember *n,
              ArbSetMember **c,
              ArbSetMember **p,
              ArbSetMember *g,
              ArbSetMember *gg)

        n  - (in)  new entry 
        c  - (i/o) child 
        p  - (i/o) parent 
        g  - (in)  grandparent 
        gg - (in)  great grandparent 

      Description: This function takes care of tree colors and 
          balances the tree. 


  Private Member Variables:

        ArbSetMember *Root - root node 

        int (*Compare)() - comparison function 

        bool DuplicatesOK - duplicates flag 

        int NumSetEntries - number of entries 


NON-MEMBER FUNCTIONS

ArbCmpDouble - comparison function for doubles 

  inline int ArbCmpDouble(
          const double &d0,
          const double &d1)

    d0 - (in)  argument one 
    d1 - (in)  argument two 

  Description: This is set entry comparison function for doubles. 

  Return Value: -1, 0, or 1 depending on the relative values of the 
      arguments 


ArbCmpUnsigned - comparison function for integers 

  inline int ArbCmpUnsigned(
          const int &u0,
          const int &u1)

    u0 - (in)  argument one 
    u1 - (in)  argument two 

  Description: This is set entry comparison function for integers. 

  Return Value: -1, 0, or 1 depending on the relative values of the 
      arguments 

*/

} // namespace

#endif
