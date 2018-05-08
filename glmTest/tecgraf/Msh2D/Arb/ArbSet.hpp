//
// CArbSet Template Class header file
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

#ifndef ArbSet_hh
#define ArbSet_hh

//#define DEBUG

template<class KeyType>
class CArbSet {

    public:

        CArbSet(int (*cmp_func)(const KeyType&,const KeyType&),
                bool allow_duplicates = false) ;
//        CArbSet() ;
        virtual ~CArbSet() ;

        bool Insert(const KeyType &key) ;
        bool Remove(const KeyType &key) ;
        bool HasKey(const KeyType &key) const ;

        KeyType **GetKeyList() const ;

        KeyType *GetSmallest() const ;
        KeyType *GetLargest() const ;

        int NumEntries() ;

// -------DEBUG----------

#ifdef DEBUG
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
        int (*Compare)(const KeyType&,const KeyType&) ;
        bool DuplicatesOK ;
        int NumSetEntries ;

        ArbSetMember *NewMember()
        {
          ArbSetMember *node = new ArbSetMember ;
          node->left = node->right = 0 ;
          node->red = false ;
          return(node) ;
        }
        void DeleteRec(ArbSetMember *node) ;
        ArbSetMember *Rotate(KeyType const &key,ArbSetMember *r,
                             int flip_mode,int dir,int gdir) ;
        void Split(ArbSetMember *n,ArbSetMember **c,ArbSetMember **p,
                   ArbSetMember *g,ArbSetMember *gg) ;

        void WalkRec(ArbSetMember *n,KeyType **keys,
                     int *cur) const ;


// -------DEBUG----------

#ifdef DEBUG
        void xrWalk(ArbSetMember *n,int level) ;
        void xrCheck(ArbSetMember *c) ;
        void DumpRec(ArbSetMember *node) ;
#endif
} ;

inline int ArbCmpDouble(const double &d0,const double &d1)
{
    if (d0 > d1) return(1) ;
    if (d0 < d1) return(-1) ;
    return(0) ;
}

inline int ArbCmpUnsigned(const int &u0,const int &u1)
{
    if (u0 > u1) return(1) ;
    if (u0 < u1) return(-1) ;
    return(0) ;
}

/*
CLASS CArbSet<class KeyType>

  This object implements a templated set data structure. A copy is made
  of an entry when it is placed in the set. The client must provide a
  comparison function that returns an interger less then, equal to, or
  greater than zero depending on the relative value of it's arguments.
  The set is implemented as a red/black binary tree.


  Template Arguments:

    KeyType - set entry type


PUBLIC INTERFACE

  Public Member Functions:

    CArbSet - constructor

      CArbSet(int (*cmp_func)())

        cmp_func - (in)  comparison function

      Description: This is a constructor for an ArbSet object.


    CArbSet - destructor

      ~CArbSet()

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

#endif
