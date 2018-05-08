//
// ArbKDTree2D class header file
//
// Description -
//   This class implements an object that will create a 
//   2D KD tree data structure
//
// Reference - 
//   Computational Geometry and Computer Graphics in C++
//   Michael J. Laszlo, Prentice Hall, 1996
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999
//   All rights reserved
//
// Revision -
//   1.0    7-Jan-00    Wash Wawrzynek
//

#ifndef ArbKDTree2D_h
#define ArbKDTree2D_h

#include "ArbMsh.hpp"
#include "ArbArray.hpp"
#include "ArbRectangle.hpp"

class CArbKDTree2D {
    public:
        CArbKDTree2D(int n,ArbIntNode **p) ;
        ~CArbKDTree2D() ;
        CArbArray<ArbIntNode*> *RangeQuery(CArbRectangle &range) const ;
        CArbArray<ArbIntNode*> *ReturnAll() const ;

    private:

        struct KDTreeEntry {
            int         num_pnts ;
            ArbIntNode  **pnts ;
            KDTreeEntry *lchild ;
            KDTreeEntry *rchild ;
        } ;

        KDTreeEntry *Root ;
        KDTreeEntry *BuildKDTree(ArbIntNode *x[], ArbIntNode *y[],
                                 int n, int cutType) ;
        KDTreeEntry *NewEntry(int n, ArbIntNode **pt_list) ;

        void RangeQueryRec(KDTreeEntry *entry,
                           CArbRectangle &R, int cutType,
                           CArbArray<ArbIntNode*> *result) const ;
        void ReturnAllRec(KDTreeEntry *entry,
                          CArbArray<ArbIntNode*> *result) const ;
        void DeleteRec(KDTreeEntry *entry) ;

        void SplitPointSet(ArbIntNode *y[], int n, ArbIntNode *p,
                           int &nL, ArbIntNode *yL[], 
                           int &nR, ArbIntNode *yR[],
                           int &nS, ArbIntNode *yS[],
                           int (*cmp)(ArbIntNode*,ArbIntNode*)) const ;
        void MergeSort(ArbIntNode *a[], int n, 
                       int(*cmp)(ArbIntNode*,ArbIntNode*)) const ;
        void MSort(ArbIntNode *a[], int l, int r, 
                   int(*cmp)(ArbIntNode*,ArbIntNode*)) const ;
        void Merge(ArbIntNode *x[], int l, int m, int r, 
                   int(*cmp)(ArbIntNode*,ArbIntNode*)) const ;
};

/*
CLASS CArbKDTree2D

  This object implements a 2D KD tree data structure. 

  A reference for this implementation is Computational Geometry and 
  Computer Graphics in C++ Michael J. Laszlo, Prentice Hall, 1996 


PUBLIC INTERFACE

  Public Member Functions:

    CArbKDTree2D - KD tree constructor 

      CArbKDTree2D(
              int        n,
              ArbIntNode **p)

        n - (in)  number of points 
        p - (in)  points to insert into the tree 

      Description: This is the constructor for a KD tree. As 
          arguments it takes a list of points to add to the tree. 


    CArbKDTree2D - KD tree destructor 

      ~CArbKDTree2D()

      Description: This is the destructor for a KD tree. It calls a 
          private routine that recursively decends the tree and 
          deletes all entries. 


    RangeQuery - return all points in a given range 

      CArbArray <ArbIntNode*>*RangeQuery(CArbRectangle &range) const

        range - (in)  search rectangle 

      Description: This function returns an array of all the points 
          found within the specified rectangle. The caller takes 
          ownership of the array and must delete it. 

      Return Value: An array of all points found in the input 
          rectangle. 


    ReturnAll - return all points in the tree 

      CArbArray <ArbIntNode*>*ReturnAll() const

      Description: This function returns an array of all the points 
          currently stored in the tree. The caller takes ownership of 
          the array and must delete it. 

      Return Value: An array of all points in the tree. 


PRIVATE INTERFACE

  Private Data Structures:

    struct KDTreeEntry

      This is the data type used as a tree node. 

      Member Variables:

        int num_pnts - number of points stored at this node 

        ArbIntNode **pnts - array of points stored at this node 

        KDTreeEntry *lchild - left child 

        KDTreeEntry *rchild - right child 


  Private Member Functions:

    BuildKDTree - build a KD tree 

      KDTreeEntry *BuildKDTree(
              ArbIntNode *x[],
              ArbIntNode *y[],
              int        n,
              int        cutType)

        x       - (in)  input points sorted left/right 
        y       - (in)  input points sorted up/down 
        n       - (in)  number of points 
        cutType - (in)  flag for left/right or up/down node 

      Description: Given a list of points and a cut direction, this 
          function splits the points about the mid point and stores 
          the lists as left and right children of a new tree node. 

      Return Value: The newly created tree node. 


    NewEntry - create a new entry in the tree 

      KDTreeEntry *NewEntry(
              int        n,
              ArbIntNode **pt_list)

        n       - (in)  number of points 
        pt_list - (in)  points to add to this node 

      Description: This function creates a new tree node and 
          initializes it with the points passed as an argument. 

      Return Value: A pointer to a new tree node. 


    RangeQueryRec - range query routine 

      void RangeQueryRec(
              KDTreeEntry             *entry,
              CArbRectangle           &R,
              int                     cutType,
              CArbArray<ArbIntNode*>* result) const

        entry   - (in)  tree node 
        R       - (in)  query rectangle 
        cutType - (in)  flag for a left/right or up/down split 
        result  - (i/o) array to which points are added 

      Description: This function visits all tree nodes recursively 
          and adds all points found in the query rectangle to the 
          array. 


    ReturnAllRec - return all routine 

      void ReturnAllRec(
              KDTreeEntry             *entry,
              CArbArray<ArbIntNode*>* result) const

        entry  - (in)  tree node 
        result - (i/o) array to which points are added 

      Description: This function visits all tree nodes recursively 
          and adds all points in the tree to the array. 


    DeleteRec - delete all tree entries 

      void DeleteRec(KDTreeEntry *entry)

        entry - (in)  root node 

      Description: This function recursively deletes all entries in 
          the tree. 


    SplitPointSet - split the point set into three bins 

      void SplitPointSet(
              ArbIntNode *y[],
              int        n,
              ArbIntNode *p,
              int        &nL,
              ArbIntNode *yL[],
              int        &nR,
              ArbIntNode *yR[],
              int        &nS,
              ArbIntNode *yS[],
              int        (*cmp)()) const

        y   - (in)  list of input points 
        n   - (in)  number of input points 
        p   - (in)  reference point 
        nL  - (out) number in the less than bin 
        yL  - (out) less than bin 
        nR  - (out) number in the greater than bin 
        yR  - (out) greater than bin 
        nS  - (out) number in the equal bin 
        yS  - (out) equal bin 
        cmp - (in)  compare function 

      Description: This function takes a reference point and splits 
          the input point set into three bins, those less than, equal 
          to, and greater than the reference point. The definition of 
          <, =, or > is specified in the compare function. 


    MergeSort - sort a list 

      void MergeSort(
              ArbIntNode *a[],
              int        n,
              int        (*cmp)()) const

        a   - (i/o) points to sort 
        n   - (in)  number of points to sort 
        cmp - (in)  compare function 

      Description: Driver function for a merge sort routine. 


    MSort - recursive merge sort 

      void MSort(
              ArbIntNode *a[],
              int        l,
              int        r,
              int        (*cmp)()) const

        a   - (i/o) points to sort 
        l   - (in)  left index in a 
        r   - (in)  right index in a 
        cmp - (in)  compare function 

      Description: This is a recursive routine used during a merge 
          sort. 


    Merge - merge function 

      void Merge(
              ArbIntNode *x[],
              int        l,
              int        m,
              int        r,
              int        (*cmp)()) const

        x   - (in)  points to sort 
        l   - (in)  left index in x 
        m   - (in)  right index in x 
        r   - (in)  mid index in x 
        cmp - (in)  compare function 

      Description: This is a merge function that is used during the 
          merge sort. 


  Private Member Variables:

        KDTreeEntry *Root - tree data 
*/
#endif

