//
// KDTree2D class header file
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

#ifndef KDTree2D_h
#define KDTree2D_h

#include "List.hpp"

using FTools::List ;

namespace Msh2D {

template<class NodeType> class KDTree2D {

    public:

        enum { VERTICAL=0, HORIZONTAL=1 } ;

        class Rectangle {
            public:
                Rectangle(const NodeType &p0, const NodeType &p1,
                              const double &size) {
                    double xc, yc ;
                    xc = 0.5*(p0.coord[0]+p1.coord[0]);
                    yc = 0.5*(p0.coord[1]+p1.coord[1]);
                    sw.coord[0] = xc - size ;
                    sw.coord[1] = yc - size ;
                    ne.coord[0] = xc + size ;
                    ne.coord[1] = yc + size ; } ;
                Rectangle() {} ;

                bool PointInRectangle(const NodeType* p) const {
                    return ((sw.coord[0] <= p->coord[0]) &&
                            (p->coord[0] <= ne.coord[0]) &&
                            (sw.coord[1] <= p->coord[1]) &&
                            (p->coord[1] <= ne.coord[1])); } ;

                NodeType *SW() { return(&sw) ; } ;
                NodeType *NE() { return(&ne) ; } ;

            private:
                NodeType sw;
                NodeType ne;
        } ;

        KDTree2D(int n,NodeType **p) ;
        KDTree2D(const List<NodeType*>& p) ;
        KDTree2D(const List<Vec2D>& p) ;
        ~KDTree2D() ;
        List<NodeType*> *RangeQuery(Rectangle &range) const ;
        List<NodeType*> *ReturnAll() const ;

    private:

        struct KDTreeEntry {
            int         num_pnts ;
            NodeType  **pnts ;
            KDTreeEntry *lchild ;
            KDTreeEntry *rchild ;
            KDTreeEntry() : num_pnts(0),pnts(0),lchild(0),rchild(0) {}
        } ;

        KDTreeEntry *Root ;
        KDTreeEntry *BuildKDTree(NodeType *x[], NodeType *y[],
                                 int n, int cutType) {
            int nL, nR, nS;

            if (n==0)
                return 0;
            else if (n==1)
                return(NewEntry(n,x)) ;

            int m = n/2 ;
            int (*cmp)(NodeType*, NodeType*) ;
            cmp = (cutType == VERTICAL) ? LeftToRightCmp : BottomToTopCmp ;

            NodeType **yL = new NodeType*[m] ;
            NodeType **yR = new NodeType*[n-m] ;
            NodeType **yS = new NodeType*[n];
            SplitPointSet(y, n, x[m], nL, yL, nR, yR, nS, yS, cmp) ;
            KDTreeEntry *p = NewEntry(nS, yS) ;
            p->lchild = BuildKDTree(yL, x, nL, 1-cutType) ;
            p->rchild = BuildKDTree(yR, x+nL+nS, nR, 1-cutType) ;
            delete [] yL ;
            delete [] yR ;
            delete [] yS;
            return(p) ;
        }

        KDTreeEntry *NewEntry(int n, NodeType **pt_list) {
            KDTreeEntry *entry = new KDTreeEntry ;
            entry->num_pnts = n ;
            entry->pnts = new NodeType*[n] ;
            for (int i=0; i<n; i++ ) {
                entry->pnts[i] = pt_list[i];
            }
            entry->rchild = entry->lchild = 0 ;
            return(entry) ;
        }

        void RangeQueryRec(KDTreeEntry *entry,
                           Rectangle &R, int cutType,
                           List<NodeType*> *result) const ;
        void ReturnAllRec(KDTreeEntry *entry,
                          List<NodeType*> *result) const ;
        void DeleteRec(KDTreeEntry *entry) ;

        void SplitPointSet(NodeType *y[], int n, NodeType *p,
                           int &nL, NodeType *yL[], 
                           int &nR, NodeType *yR[],
                           int &nS, NodeType *yS[],
                           int (*cmp)(NodeType*,NodeType*)) const ;
        void MergeSort(NodeType *a[], int n, 
                       int(*cmp)(NodeType*,NodeType*)) const ;
        void MSort(NodeType *a[], int l, int r, 
                   int(*cmp)(NodeType*,NodeType*)) const ;
        void Merge(NodeType *x[], int l, int m, int r, 
                   int(*cmp)(NodeType*,NodeType*)) const ;

        static int LeftToRightCmp(NodeType *a, NodeType *b) ;
        static int BottomToTopCmp(NodeType *a, NodeType *b) ;
};

/* ++ ----------------------------------------------------------
**
**    KDTree2D - KD tree constructor 
**
**      KDTree2D(
**              int        n,
**              NodeType **p)
**
**        n - (in)  number of points 
**        p - (in)  points to insert into the tree 
**
**      Description: This is the constructor for a KD tree. As 
**          arguments it takes a list of points to add to the tree. 
**
**
** -- */

template<class NodeType> KDTree2D<NodeType>::KDTree2D(int n,NodeType **p)
{
    NodeType **x = new NodeType*[n];
    NodeType **y = new NodeType*[n];
    for (int i=0; i<n; i++) {
        x[i] = y[i] = p[i];
    }
// printf("%d\n",n) ;
// for (int j=0 ; j<n ; ++j)
//   printf("t %g %g %d\n",p[j]->coord[0],p[j]->coord[1],p[j]->id) ;
// fflush(stdout) ;
    MergeSort(x, n, LeftToRightCmp);
    MergeSort(y, n, BottomToTopCmp);
    Root = BuildKDTree(x, y, n, VERTICAL);
    delete [] x;
    delete [] y;
}


template<class NodeType> KDTree2D<NodeType>::KDTree2D(const List<NodeType*> &p)
{
    NodeType **x = new NodeType*[p.Len()];
    NodeType **y = new NodeType*[p.Len()];
    for (int i=0; i<p.Len(); i++) {
        x[i] = y[i] = p[i];
    }
// printf("%d\n",n) ;
// for (int j=0 ; j<n ; ++j)
//   printf("t %g %g %d\n",p[j]->coord[0],p[j]->coord[1],p[j]->id) ;
// fflush(stdout) ;
    MergeSort(x, p.Len(), LeftToRightCmp);
    MergeSort(y, p.Len(), BottomToTopCmp);
    Root = BuildKDTree(x, y, p.Len(), VERTICAL);
    delete [] x;
    delete [] y;
}


template<class NodeType> KDTree2D<NodeType>::KDTree2D(const List<Vec2D> &p)
{
    NodeType **x = new NodeType*[p.Len()];
    NodeType **y = new NodeType*[p.Len()];
    for (int i=0; i<p.Len(); i++) {
        x[i] = y[i] = p[i];
    }
// printf("%d\n",n) ;
// for (int j=0 ; j<n ; ++j)
//   printf("t %g %g %d\n",p[j]->coord[0],p[j]->coord[1],p[j]->id) ;
// fflush(stdout) ;
    MergeSort(x, p.Len(), LeftToRightCmp);
    MergeSort(y, p.Len(), BottomToTopCmp);
    Root = BuildKDTree(x, y, p.Len(), VERTICAL);
    delete [] x;
    delete [] y;
}


// %(KDTree2D::BuildKDTree-KDTreeEntry-|*-NodeType-|*^[]-NodeType-|*^[]-int-|-int-|) 
/* ++ ----------------------------------------------------------
**
**    BuildKDTree - build a KD tree 
**
**      KDTreeEntry *BuildKDTree(
**              NodeType *x[],
**              NodeType *y[],
**              int        n,
**              int        cutType)
**
**        x       - (in)  input points sorted left/right 
**        y       - (in)  input points sorted up/down 
**        n       - (in)  number of points 
**        cutType - (in)  flag for left/right or up/down node 
**
**      Description: Given a list of points and a cut direction, this 
**          function splits the points about the mid point and stores 
**          the lists as left and right children of a new tree node. 
**
**      Return Value: The newly created tree node. 
**
**
** -- */

#if 0
//template<class NodeType>
    KDTree2D::KDTreeEntry *
    KDTree2D<NodeType>::BuildKDTree(NodeType *x[],
                                        NodeType *y[],
                                        int n, int cutType)
{
    int nL, nR, nS;

    if (n==0)
        return 0;
    else if (n==1)
        return(NewEntry(n,x)) ;

    int m = n/2 ;
    int (*cmp)(NodeType*, NodeType*) ;
    cmp = (cutType == VERTICAL) ? LeftToRightCmp : BottomToTopCmp ;

    NodeType **yL = new NodeType*[m] ;
    NodeType **yR = new NodeType*[n-m] ;
    NodeType **yS = new NodeType*[n];
    SplitPointSet(y, n, x[m], nL, yL, nR, yR, nS, yS, cmp) ;
    KDTreeEntry *p = NewEntry(nS, yS) ;
    p->lchild = BuildKDTree(yL, x, nL, 1-cutType) ;
    p->rchild = BuildKDTree(yR, x+nL+nS, nR, 1-cutType) ;
    delete [] yL ;
    delete [] yR ;
    delete [] yS;
    return(p) ;
}




// %(KDTree2D::NewEntry-KDTreeEntry-|*-int-|-NodeType-|**)
/* ++ ----------------------------------------------------------
**
**    NewEntry - create a new entry in the tree 
**
**      KDTreeEntry *NewEntry(
**              int        n,
**              NodeType **pt_list)
**
**        n       - (in)  number of points 
**        pt_list - (in)  points to add to this node 
**
**      Description: This function creates a new tree node and 
**          initializes it with the points passed as an argument. 
**
**      Return Value: A pointer to a new tree node. 
**
**
** -- */

template<class NodeType>
    KDTree2D<NodeType>::KDTreeEntry*
    KDTree2D<NodeType>::NewEntry(int n, NodeType **pt_list)
{
//    fprintf(stderr,"# of points: %d\n", n);

    KDTreeEntry *entry = new KDTreeEntry ;
    entry->num_pnts = n ;
    entry->pnts = new NodeType*[n] ;
    for (int i=0; i<n; i++ ) {
        entry->pnts[i] = pt_list[i];
//        fprintf(stderr,"Point id: %d\n", entry->pnts[i]->id);
    }
    entry->rchild = entry->lchild = 0 ;
    return(entry) ;
}

#endif


// %(KDTree2D::KDTree2D-destructor-|~)
/* ++ ----------------------------------------------------------
**
**    KDTree2D - KD tree destructor 
**
**      ~KDTree2D()
**
**      Description: This is the destructor for a KD tree. It calls a 
**          private routine that recursively decends the tree and 
**          deletes all entries. 
**
**
** -- */

template<class NodeType> KDTree2D<NodeType>::~KDTree2D()
{
    DeleteRec(Root) ;
}




// %(KDTree2D::DeleteRec-void-|-KDTreeEntry-|*)
/* ++ ----------------------------------------------------------
**
**    DeleteRec - delete all tree entries 
**
**      void DeleteRec(KDTreeEntry *entry)
**
**        entry - (in)  root node 
**
**      Description: This function recursively deletes all entries in 
**          the tree. 
**
**
** -- */

template<class NodeType> void KDTree2D<NodeType>::DeleteRec(KDTreeEntry *entry)
{
    if (entry->num_pnts > 0) delete [] entry->pnts;
    if (entry->lchild) DeleteRec(entry->lchild) ;
    if (entry->rchild) DeleteRec(entry->rchild) ;
    delete entry ;
}




// %(KDTree2D::RangeQuery-List-|<NodeType*>*^const-CArbRectangle-|&)
/* ++ ----------------------------------------------------------
**
**    RangeQuery - return all points in a given range 
**
**      List <NodeType*>*RangeQuery(CArbRectangle &range) const
**
**        range - (in)  search rectangle 
**
**      Description: This function returns an array of all the points 
**          found within the specified rectangle. The caller takes 
**          ownership of the array and must delete it. 
**
**      Return Value: An array of all points found in the input 
**          rectangle. 
**
**
** -- */

template<class NodeType> List<NodeType*> *KDTree2D<NodeType>::RangeQuery(Rectangle &R) const
{
    List<NodeType*> *result = new List<NodeType*>() ;
    RangeQueryRec(Root,R,VERTICAL,result) ;
    return result;
}




// %(KDTree2D::RangeQueryRec-void-|^const-KDTreeEntry-|*-CArbRectangle-|&-int-|-List-|<NodeType*>*) 
/* ++ ----------------------------------------------------------
**
**    RangeQueryRec - range query routine 
**
**      void RangeQueryRec(
**              KDTreeEntry             *entry,
**              CArbRectangle           &R,
**              int                     cutType,
**              List<NodeType*>* result) const
**
**        entry   - (in)  tree node 
**        R       - (in)  query rectangle 
**        cutType - (in)  flag for a left/right or up/down split 
**        result  - (i/o) array to which points are added 
**
**      Description: This function visits all tree nodes recursively 
**          and adds all points found in the query rectangle to the 
**          array. 
**
**
** -- */

template<class NodeType> void KDTree2D<NodeType>::RangeQueryRec(KDTreeEntry *entry,
                                 Rectangle &R, int cutType,
                                 List<NodeType*> *result) const
{
    if (R.PointInRectangle(entry->pnts[0])) {
        for (int i=0; i<entry->num_pnts; i++)
            result->Append(entry->pnts[i]) ;
    }

    int (*cmp)(NodeType*, NodeType*) ;
    cmp = (cutType==VERTICAL) ? LeftToRightCmp : BottomToTopCmp ;

    if (entry->lchild && ((*cmp)(R.SW(), entry->pnts[0]) < 0))
        RangeQueryRec(entry->lchild,R,1-cutType,result) ;

    if (entry->rchild && ((*cmp)(R.NE(), entry->pnts[0]) > 0))
        RangeQueryRec(entry->rchild,R,1-cutType,result) ;
}




// %(KDTree2D::ReturnAll-List-|<NodeType*>*^const)
/* ++ ----------------------------------------------------------
**
**    ReturnAll - return all points in the tree 
**
**      List <NodeType*>*ReturnAll() const
**
**      Description: This function returns an array of all the points 
**          currently stored in the tree. The caller takes ownership of 
**          the array and must delete it. 
**
**      Return Value: An array of all points in the tree. 
**
**
** -- */

template<class NodeType> List<NodeType*> *KDTree2D<NodeType>::ReturnAll() const
{
    List<NodeType*> *result = new List<NodeType*>() ;
    ReturnAllRec(Root,result) ;
    return result;
}




// %(KDTree2D::ReturnAllRec-void-|^const-KDTreeEntry-|*-List-|<NodeType*>*) 
/* ++ ----------------------------------------------------------
**
**    ReturnAllRec - return all routine 
**
**      void ReturnAllRec(
**              KDTreeEntry             *entry,
**              List<NodeType*>* result) const
**
**        entry  - (in)  tree node 
**        result - (i/o) array to which points are added 
**
**      Description: This function visits all tree nodes recursively 
**          and adds all points in the tree to the array. 
**
**
** -- */

template<class NodeType> void KDTree2D<NodeType>::ReturnAllRec(KDTreeEntry *entry,
                                List<NodeType*> *result) const
{
    for (int i=0; i<entry->num_pnts; i++)
        result->Append(entry->pnts[i]) ;
    if (entry->lchild) ReturnAllRec(entry->lchild,result) ;
    if (entry->rchild) ReturnAllRec(entry->rchild,result) ;
}




// %(KDTree2D::SplitPointSet-void-|^const-NodeType-|*^[]-int-|-NodeType-|*-int-|&-NodeType-|*^[]-int-|&-NodeType-|*^[]-int-|&-NodeType-|*^[]-int-|(*^)()) 
/* ++ ----------------------------------------------------------
**
**    SplitPointSet - split the point set into three bins 
**
**      void SplitPointSet(
**              NodeType *y[],
**              int        n,
**              NodeType *p,
**              int        &nL,
**              NodeType *yL[],
**              int        &nR,
**              NodeType *yR[],
**              int        &nS,
**              NodeType *yS[],
**              int        (*cmp)()) const
**
**        y   - (in)  list of input points 
**        n   - (in)  number of input points 
**        p   - (in)  reference point 
**        nL  - (out) number in the less than bin 
**        yL  - (out) less than bin 
**        nR  - (out) number in the greater than bin 
**        yR  - (out) greater than bin 
**        nS  - (out) number in the equal bin 
**        yS  - (out) equal bin 
**        cmp - (in)  compare function 
**
**      Description: This function takes a reference point and splits 
**          the input point set into three bins, those less than, equal 
**          to, and greater than the reference point. The definition of 
**          <, =, or > is specified in the compare function. 
**
**
** -- */

template<class NodeType> void KDTree2D<NodeType>::SplitPointSet(NodeType *y[], int n, NodeType *p,
                                 int &nL, NodeType *yL[], 
                                 int &nR, NodeType *yR[],
                                 int &nS, NodeType *yS[],
                                 int (*cmp)(NodeType*,NodeType*)) const
{
  nL = 0, nR = 0, nS = 0;
  for (int i=0; i<n; i++) {
    if ((*cmp)(y[i], p) < 0)
      yL[nL++] = y[i];
    else if ((*cmp)(y[i],p) > 0)
      yR[nR++] = y[i];
    else // actual point p or equal points
      yS[nS++] = y[i];
  }
}




// %(KDTree2D::MergeSort-void-|^const-NodeType-|*^[]-int-|-int-|(*^)())
/* ++ ----------------------------------------------------------
**
**    MergeSort - sort a list 
**
**      void MergeSort(
**              NodeType *a[],
**              int        n,
**              int        (*cmp)()) const
**
**        a   - (i/o) points to sort 
**        n   - (in)  number of points to sort 
**        cmp - (in)  compare function 
**
**      Description: Driver function for a merge sort routine. 
**
**
** -- */

template<class NodeType> void KDTree2D<NodeType>::MergeSort(NodeType *a[], int n, 
                             int (*cmp)(NodeType*,NodeType*)) const
{
    MSort(a, 0, n-1, cmp) ;
}




// %(KDTree2D::MSort-void-|^const-NodeType-|*^[]-int-|-int-|-int-|(*^)())
/* ++ ----------------------------------------------------------
**
**    MSort - recursive merge sort 
**
**      void MSort(
**              NodeType *a[],
**              int        l,
**              int        r,
**              int        (*cmp)()) const
**
**        a   - (i/o) points to sort 
**        l   - (in)  left index in a 
**        r   - (in)  right index in a 
**        cmp - (in)  compare function 
**
**      Description: This is a recursive routine used during a merge 
**          sort. 
**
**
** -- */

template<class NodeType> void KDTree2D<NodeType>::MSort(NodeType *a[], int l, int r, 
                         int (*cmp)(NodeType*,NodeType*)) const
{
    if (l<r) {
        int m = (l+r)/2 ;
        MSort(a, l, m, cmp) ;
        MSort(a, m+1, r, cmp) ;
        Merge(a, l, m, r, cmp) ;
    }
}




// %(KDTree2D::Merge-void-|^const-NodeType-|*^[]-int-|-int-|-int-|-int-|(*^)()) 
/* ++ ----------------------------------------------------------
**
**    Merge - merge function 
**
**      void Merge(
**              NodeType *x[],
**              int        l,
**              int        m,
**              int        r,
**              int        (*cmp)()) const
**
**        x   - (in)  points to sort 
**        l   - (in)  left index in x 
**        m   - (in)  right index in x 
**        r   - (in)  mid index in x 
**        cmp - (in)  compare function 
**
**      Description: This is a merge function that is used during the 
**          merge sort. 
**
**
** -- */

template<class NodeType> void KDTree2D<NodeType>::Merge(NodeType *x[], int l, int m, int r, 
                         int (*cmp)(NodeType*,NodeType*)) const
{
    NodeType **a = x+l ;
    NodeType **b = x+m+1 ;
    NodeType **c = new NodeType*[r-l+1] ;
    int aindx=0, bindx=0, cindx=0 ;
    int alim=m-l+1, blim=r-m ;

    while ((aindx<alim) && (bindx<blim))
        if ((*cmp)(a[aindx], b[bindx]) < 0)
            c[cindx++] = a[aindx++] ;
        else
            c[cindx++] = b[bindx++] ;

    while (aindx<alim)              // copy rest of a
        c[cindx++] = a[aindx++] ;

    while (bindx<blim)              // copy rest of b
        c[cindx++] = b[bindx++] ;

    for (aindx=cindx=0; aindx<=r-l; a[aindx++]=c[cindx++])
        ;                          // copy back

    delete [] c ;
}


template<class NodeType> int KDTree2D<NodeType>::LeftToRightCmp(NodeType *a, NodeType *b)
{
    if ((a->coord[0] < b->coord[0]) || 
        ((a->coord[0] == b->coord[0]) && (a->coord[1] < b->coord[1]))) 
        return -1;
    if ((a->coord[0] > b->coord[0]) || 
        ((a->coord[0] == b->coord[0]) && (a->coord[1] > b->coord[1])))
        return 1;
    return 0;
}


template<class NodeType> int KDTree2D<NodeType>::BottomToTopCmp(NodeType *a, NodeType *b)
{
    if ((a->coord[1] < b->coord[1]) ||
        ((a->coord[1] == b->coord[1]) && (a->coord[0] < b->coord[0])))
        return -1;
    else if ((a->coord[1] > b->coord[1]) ||
             ((a->coord[1] == b->coord[1]) && (a->coord[0] > b->coord[0]))) 
        return 1;
    return 0;
}




/*
CLASS KDTree2D

  This object implements a 2D KD tree data structure. 

  A reference for this implementation is Computational Geometry and 
  Computer Graphics in C++ Michael J. Laszlo, Prentice Hall, 1996 


PUBLIC INTERFACE

  Public Member Functions:

    KDTree2D - KD tree constructor 

      KDTree2D(
              int        n,
              NodeType **p)

        n - (in)  number of points 
        p - (in)  points to insert into the tree 

      Description: This is the constructor for a KD tree. As 
          arguments it takes a list of points to add to the tree. 


    KDTree2D - KD tree destructor 

      ~KDTree2D()

      Description: This is the destructor for a KD tree. It calls a 
          private routine that recursively decends the tree and 
          deletes all entries. 


    RangeQuery - return all points in a given range 

      List <NodeType*>*RangeQuery(CArbRectangle &range) const

        range - (in)  search rectangle 

      Description: This function returns an array of all the points 
          found within the specified rectangle. The caller takes 
          ownership of the array and must delete it. 

      Return Value: An array of all points found in the input 
          rectangle. 


    ReturnAll - return all points in the tree 

      List <NodeType*>*ReturnAll() const

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

        NodeType **pnts - array of points stored at this node 

        KDTreeEntry *lchild - left child 

        KDTreeEntry *rchild - right child 


  Private Member Functions:

    BuildKDTree - build a KD tree 

      KDTreeEntry *BuildKDTree(
              NodeType *x[],
              NodeType *y[],
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
              NodeType **pt_list)

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
              List<NodeType*>* result) const

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
              List<NodeType*>* result) const

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
              NodeType *y[],
              int        n,
              NodeType *p,
              int        &nL,
              NodeType *yL[],
              int        &nR,
              NodeType *yR[],
              int        &nS,
              NodeType *yS[],
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
              NodeType *a[],
              int        n,
              int        (*cmp)()) const

        a   - (i/o) points to sort 
        n   - (in)  number of points to sort 
        cmp - (in)  compare function 

      Description: Driver function for a merge sort routine. 


    MSort - recursive merge sort 

      void MSort(
              NodeType *a[],
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
              NodeType *x[],
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

} // namespace

#endif

