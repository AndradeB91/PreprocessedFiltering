//
// CArbSet Template Class implementation file
//
// Description -
//   This class implements an arbitrary size array class.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.5 $  $Date: 2000/10/31 15:20:07 $  $Author: wash $
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "ArbMsh.hpp"
#include "ArbSet.hpp"

#define LEFT       1
#define RIGHT      0

#define NO_FLIP  0
#define FLIP_GCH 1
#define FLIP_CH  2

#ifdef MEMDEBUG
#include "MemDbg.hpp"
#define new new(__FILE__,__LINE__)
#endif

template class CArbSet<int> ;
template class CArbSet<ArbIntEdge> ;
template class CArbSet<ArbMshElement2D> ;
template class CArbSet<double> ;



// %(CArbSet::CArbSet-constructor-|-int-|(*^)())
/* ++ ----------------------------------------------------------
**
**    CArbSet - constructor
**
**      CArbSet(int (*cmp_func)())
**
**        cmp_func - (in)  comparison function
**
**      Description: This is a constructor for an ArbSet object.
**
**
** -- */

template<class KeyType>
CArbSet<KeyType>::CArbSet(int ((iCompare)(const KeyType&,const KeyType&)),
                          bool allow_duplicates)
{
    Compare = iCompare ;
    DuplicatesOK = allow_duplicates ;
    Root = NewMember() ;
    NumSetEntries = 0 ;
}




// %(CArbSet::CArbSet-destructor-virtual|~)
/* ++ ----------------------------------------------------------
**
**    CArbSet - destructor
**
**      ~CArbSet()
**
**      Description: This is a destructor for an ArbSet object.
**
**
** -- */

template<class KeyType>
CArbSet<KeyType>::~CArbSet()
{
    DeleteRec(Root) ;
}




// %(CArbSet::NumEntries-int-|)
/* ++ ----------------------------------------------------------
**
**    NumEntries - number of entries
**
**      int NumEntries()
**
**      Description: This function returns the number entries in the
**          set.
**
**      Return Value: number of set entries
**
**
** -- */

template<class KeyType>
int CArbSet<KeyType>::NumEntries()
{
    return(NumSetEntries) ;
}


#if DEBUG
template<class KeyType>
void CArbSet<KeyType>::Dump()
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
void CArbSet<KeyType>::DumpRec(ArbSetMember *node)
{
    if (node->left != 0) DumpRec(node->left) ;
    if (node->right != 0) DumpRec(node->right) ;
    if ((node->left == 0) && (node->right == 0))
        dump_one(node->key) ;
}
#endif




// %(CArbSet::DeleteRec-void-|-ArbSetMember-|*)
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
void CArbSet<KeyType>::DeleteRec(ArbSetMember *node)
{
    if (node->left != 0) DeleteRec(node->left) ;
    if (node->right != 0) DeleteRec(node->right) ;
//    if ((node->left == 0) && (node->right == 0))
//        dump_one(node->key) ;
    delete node ;
}




// %(CArbSet::NewMember-ArbSetMember-|*)
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

#if 0
template<class KeyType>
CArbSet<KeyType>::ArbSetMember *CArbSet<class KeyType>::NewMember()
{
    ArbSetMember *node = new ArbSetMember ;
    node->left = node->right = 0 ;
    node->red = false ;
    return(node) ;
}
#endif



// %(CArbSet::HasKey-bool-|^const-KeyType-const|&)
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
bool CArbSet<KeyType>::HasKey(const KeyType &key) const
{
    ArbSetMember *s ;
    int dir ;

    s = Root->right ;
    while (s != NULL) {
        dir = (Compare)(key,s->key) ;
        if ((dir == 0) && (s->right == NULL)) return(true) ;
        if (dir < 0)
            s = s->left ;
        else
            s = s->right ;
    }
    return(false) ;
}




// %(CArbSet::GetSmallest-KeyType-|*^const)
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

template<class KeyType>KeyType *CArbSet<KeyType>::GetSmallest() const
{
    ArbSetMember *p, *s ;

    p = Root ;
    s = p->right ;
    while (s->left != NULL) {
        s = s->left ;
    }
    return(&(s->key)) ;
}




// %(CArbSet::GetLargest-KeyType-|*^const)
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

template<class KeyType>KeyType *CArbSet<KeyType>::GetLargest() const
{
    ArbSetMember *p, *s ;

    p = Root ;
    s = p->right ;
    while (s->right != NULL) {
        s = s->right ;
    }
    return(&(s->key)) ;
}




// %(CArbSet::Insert-bool-|-KeyType-const|&)
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
bool CArbSet<KeyType>::Insert(const KeyType &key)
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
            if ((Compare)(key,s->key) < 0) {
                s = s->left ;
                p_dir = LEFT ;
            } else {
                s = s->right ;
                p_dir = RIGHT ;
            }
        }

        dir = (Compare)(n->key,s->key) ;
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

#ifdef DEBUG
    xCheck() ;
#endif

    return(true) ;
}




// %(CArbSet::Remove-bool-|-KeyType-const|&)
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
bool CArbSet<KeyType>::Remove(const KeyType &key)
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
        if ((Compare)(key,s->key) == 0) {
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

        if ((Compare)(key,s->key) < 0) {
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
            assert(p->red == true) ;
            if (dir == LEFT)
                sib = p->right ;
            else
                sib = p->left ;

            assert(sib->red == false) ;
            assert(sib->left != NULL) ;

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

    if ((Compare)(s->key,key) == 0) {
        assert((p->red == true) || (p == Root->right)) ;
        if (((Compare)(s->key,g->key) < 0) && (g != Root)) {
            if ((Compare)(s->key,p->key) < 0)
                g->left = p->right ;
            else
                g->left = p->left ;
        } else {
            if ((Compare)(s->key,p->key) < 0)
                g->right = p->right ;
            else
                g->right = p->left ;
        }

        delete p ;
        delete s ;
        NumSetEntries-- ;

#ifdef DEBUG
        xCheck() ;
#endif
        return(true) ;
    } else {

#ifdef DEBUG
        xCheck() ;
#endif
        return(false) ;
    }
}




// %(CArbSet::Rotate-ArbSetMember-|*-KeyType-|const&-ArbSetMember-|*-int-|-int-|-int-|)
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
typename CArbSet<KeyType>::ArbSetMember
*CArbSet<KeyType>::Rotate(KeyType const &key,
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
            ch_dir = ((Compare)(key,r->key) < 0) ? LEFT : RIGHT ;
        if (flip_mode & FLIP_CH) ch_dir = !ch_dir ;
    } else
        ch_dir = dir ;

    if (r == Root) ch_dir = RIGHT ;  // special case
    ch = (ch_dir == RIGHT) ? r->right : r->left ;

    if (gdir < 0) {
        gch_dir = ((Compare)(key,ch->key) < 0) ? LEFT : RIGHT ;
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




// %(CArbSet::Split-void-|-ArbSetMember-|*-ArbSetMember-|**-ArbSetMember-|**-ArbSetMember-|*-ArbSetMember-|*)
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
void CArbSet<KeyType>::Split(ArbSetMember *n,ArbSetMember **c,
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
void CArbSet<KeyType>::WalkRec(ArbSetMember *n,KeyType **keys,
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




// %(CArbSet::GetKeyList-KeyType-|**^const)
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
KeyType **CArbSet<KeyType>::GetKeyList() const
{
    if (NumSetEntries == 0) return(0) ;

    KeyType **keys = new KeyType*[NumSetEntries] ;
    int cur = 0 ;

    WalkRec(Root->right,keys,&cur) ;
    return(keys) ;
}

// --------------------------------------------

#ifdef DEBUG
template<class KeyType> void CArbSet<KeyType>::xCheck()
{
    xrCheck(Root->right) ;
}

template<class KeyType> void CArbSet<KeyType>::xrCheck(ArbSetMember *c)
{
    if (c->left == 0) {
        assert(c->right == 0) ;
        return ;
    }
    if (c->right == 0) {
        assert(c->left == 0) ;
        return ;
    }
    if (c->red) {
        assert(c->left->red == false) ;
        assert(c->right->red == false) ;
    }
//    assert(c->key >= c->left->key) ;
//    assert(c->key <= c->right->key) ;
    if (c->left != 0) xrCheck(c->left) ;
    if (c->right != 0) xrCheck(c->right) ;
}
#endif

// --------------------------------------------


#if 0

#define TOP '+'
#define BOT '+'
#define HOR '-'
#define VRT '|'

#define RTOP '+'
#define RBOT '+'
#define RHOR '>'
#define RVRT '|'

#define DRAWBUF 100
static char draw[DRAWBUF] ;
static char work[DRAWBUF*2] ;
static int maxdepth ;
static int blackheight ;
static int maxblack ;
static FILE *outfile ;

template<class KeyType>
void CArbSet<KeyType>::xrWalk(ArbSetMember *n,int level)
{
    int i ;

    if (n != NULL) {
        if (level > maxdepth) maxdepth = level ;
        if (!n->red) blackheight++ ;

        // go right

        draw[level*2] = TOP ;
        if (n->right && n->right->red) draw[level*2] = RTOP ;
        draw[level*2+1] = ' ' ;
        xrWalk(n->right,level+1) ;

        // show the current node

        strncpy(work,draw,level*2) ;
        if (level > 0) {
            int c ;

            c = work[0] ;
            for (i=2 ; i<level*2 ; i += 2) {
                if (((c == TOP || c == RTOP) &&
                    (work[i] == TOP || work[i] == RTOP)) ||
                    ((c == BOT || c == RBOT) &&
                    (work[i] == BOT || work[i] == RBOT)))
                    work[i-1] = ' ' ;
                else
                    c = work[i] ;
            }

            work[level*2-1] = n->red ? RHOR : HOR ;

            for (i=0 ; i<level*2-2 ; i += 2) {
                if (work[i] != ' ') {
                    if (work[i] == TOP || work[i] == BOT)
                        work[i] = VRT ;
                    else
                        work[i] = RVRT ;
                }
            }
        }

        sprintf(work+level*2,"%s (%d) %c",
                n->key,level,
                n->red ? 'r' : 'b') ;
        fputs(work,outfile) ;

        if (n->left == NULL) {  // leaf
            if (maxblack < 0) maxblack = blackheight ;
            else if (maxblack != blackheight)
                fprintf(outfile,"  Leaf has black height %d!",
                        blackheight-1) ;
        }
        fputs("\n",outfile) ;

        // go left

        draw[level*2] = BOT ;
        if (n->left && n->left->red) draw[level*2] = RBOT ;
        draw[level*2+1] = ' ' ;
        xrWalk(n->left,level+1) ;
        if(!n->red) blackheight-- ;
    }
}

template<class KeyType>
void CArbSet<KeyType>::xWalkRBTree(char *name,char *mode)
{
    if (Root->right == NULL) {
        fputs("Empty tree\n",stdout) ;
    }

    maxdepth = -1 ;
    blackheight = 0 ;
    maxblack = -1 ;

    outfile = stdout ;
    if (name) {
        outfile = fopen(name,mode) ;
        if (outfile == NULL) {
            fprintf(stdout,"Can not open %s\n",name) ;
            name = NULL ;
            outfile = stdout ;
        }
    }

    xrWalk(Root->right,0) ;
    fprintf(outfile,"Max depth %d, black height %d\n",
            maxdepth,maxblack-1) ;

    if (name) fclose(outfile) ;
    else fflush(outfile) ;
}

#endif
