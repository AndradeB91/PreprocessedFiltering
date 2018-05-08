//
// ArbKDTree2D class description
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
//   1.0    7-Jan-00    Bruce Carter
//

#include "ArbKDTree2D.hpp"
#include <stdio.h>

#ifdef MEMDEBUG
#include "MemDbg.hpp"
#define new new(__FILE__,__LINE__)
#endif

static int LeftToRightCmp(ArbIntNode *a, ArbIntNode *b) ;
//static int RightToLeftCmp(ArbIntNode *a, ArbIntNode *b) ;
static int BottomToTopCmp(ArbIntNode *a, ArbIntNode *b) ;
//static int TopToBottomCmp(ArbIntNode *a, ArbIntNode *b) ;

enum { VERTICAL=0, HORIZONTAL=1 } ;


// %(CArbKDTree2D::CArbKDTree2D-constructor-|-int-|-ArbIntNode-|**) 
/* ++ ----------------------------------------------------------
**
**    CArbKDTree2D - KD tree constructor 
**
**      CArbKDTree2D(
**              int        n,
**              ArbIntNode **p)
**
**        n - (in)  number of points 
**        p - (in)  points to insert into the tree 
**
**      Description: This is the constructor for a KD tree. As 
**          arguments it takes a list of points to add to the tree. 
**
**
** -- */

CArbKDTree2D::CArbKDTree2D(int n,ArbIntNode **p)
{
    ArbIntNode **x = new ArbIntNode*[n];
    ArbIntNode **y = new ArbIntNode*[n];
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




// %(CArbKDTree2D::BuildKDTree-KDTreeEntry-|*-ArbIntNode-|*^[]-ArbIntNode-|*^[]-int-|-int-|) 
/* ++ ----------------------------------------------------------
**
**    BuildKDTree - build a KD tree 
**
**      KDTreeEntry *BuildKDTree(
**              ArbIntNode *x[],
**              ArbIntNode *y[],
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

CArbKDTree2D::KDTreeEntry *CArbKDTree2D::BuildKDTree(ArbIntNode *x[],
                                                     ArbIntNode *y[],
                                                     int n, int cutType)
{
    int nL, nR, nS;

    if (n==0)
        return 0;
    else if (n==1)
        return(NewEntry(n,x)) ;

    int m = n/2 ;
    int (*cmp)(ArbIntNode*, ArbIntNode*) ;
    cmp = (cutType == VERTICAL) ? LeftToRightCmp : BottomToTopCmp ;

    ArbIntNode **yL = new ArbIntNode*[m] ;
    ArbIntNode **yR = new ArbIntNode*[n-m] ;
    ArbIntNode **yS = new ArbIntNode*[n];
    SplitPointSet(y, n, x[m], nL, yL, nR, yR, nS, yS, cmp) ;
    KDTreeEntry *p = NewEntry(nS, yS) ;
    p->lchild = BuildKDTree(yL, x, nL, 1-cutType) ;
    p->rchild = BuildKDTree(yR, x+nL+nS, nR, 1-cutType) ;
    delete [] yL ;
    delete [] yR ;
    delete [] yS;
    return(p) ;
}




// %(CArbKDTree2D::NewEntry-KDTreeEntry-|*-int-|-ArbIntNode-|**)
/* ++ ----------------------------------------------------------
**
**    NewEntry - create a new entry in the tree 
**
**      KDTreeEntry *NewEntry(
**              int        n,
**              ArbIntNode **pt_list)
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

CArbKDTree2D::KDTreeEntry *CArbKDTree2D::NewEntry(int n, ArbIntNode **pt_list)
{
//    fprintf(stderr,"# of points: %d\n", n);

    KDTreeEntry *entry = new KDTreeEntry ;
    entry->num_pnts = n ;
    entry->pnts = new ArbIntNode*[n] ;
    for (int i=0; i<n; i++ ) {
        entry->pnts[i] = pt_list[i];
//        fprintf(stderr,"Point id: %d\n", entry->pnts[i]->id);
    }
    entry->rchild = entry->lchild = 0 ;
    return(entry) ;
}




// %(CArbKDTree2D::CArbKDTree2D-destructor-|~)
/* ++ ----------------------------------------------------------
**
**    CArbKDTree2D - KD tree destructor 
**
**      ~CArbKDTree2D()
**
**      Description: This is the destructor for a KD tree. It calls a 
**          private routine that recursively decends the tree and 
**          deletes all entries. 
**
**
** -- */

CArbKDTree2D::~CArbKDTree2D()
{
    DeleteRec(Root) ;
}




// %(CArbKDTree2D::DeleteRec-void-|-KDTreeEntry-|*)
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

void CArbKDTree2D::DeleteRec(KDTreeEntry *entry)
{
    if (entry->num_pnts > 0) delete [] entry->pnts;
    if (entry->lchild) DeleteRec(entry->lchild) ;
    if (entry->rchild) DeleteRec(entry->rchild) ;
    delete entry ;
}




// %(CArbKDTree2D::RangeQuery-CArbArray-|<ArbIntNode*>*^const-CArbRectangle-|&)
/* ++ ----------------------------------------------------------
**
**    RangeQuery - return all points in a given range 
**
**      CArbArray <ArbIntNode*>*RangeQuery(CArbRectangle &range) const
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

CArbArray<ArbIntNode*> *CArbKDTree2D::RangeQuery(CArbRectangle &R) const
{
    CArbArray<ArbIntNode*> *result = new CArbArray<ArbIntNode*>() ;
    RangeQueryRec(Root,R,VERTICAL,result) ;
    return result;
}




// %(CArbKDTree2D::RangeQueryRec-void-|^const-KDTreeEntry-|*-CArbRectangle-|&-int-|-CArbArray-|<ArbIntNode*>*) 
/* ++ ----------------------------------------------------------
**
**    RangeQueryRec - range query routine 
**
**      void RangeQueryRec(
**              KDTreeEntry             *entry,
**              CArbRectangle           &R,
**              int                     cutType,
**              CArbArray<ArbIntNode*>* result) const
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

void CArbKDTree2D::RangeQueryRec(KDTreeEntry *entry,
                                 CArbRectangle &R, int cutType,
                                 CArbArray<ArbIntNode*> *result) const
{
    if (R.PointInRectangle(entry->pnts[0])) {
        for (int i=0; i<entry->num_pnts; i++)
            result->InsertAtEnd(entry->pnts[i]) ;
    }

    int (*cmp)(ArbIntNode*, ArbIntNode*) ;
    cmp = (cutType==VERTICAL) ? LeftToRightCmp : BottomToTopCmp ;

    if (entry->lchild && ((*cmp)(R.SW(), entry->pnts[0]) < 0))
        RangeQueryRec(entry->lchild,R,1-cutType,result) ;

    if (entry->rchild && ((*cmp)(R.NE(), entry->pnts[0]) > 0))
        RangeQueryRec(entry->rchild,R,1-cutType,result) ;
}




// %(CArbKDTree2D::ReturnAll-CArbArray-|<ArbIntNode*>*^const)
/* ++ ----------------------------------------------------------
**
**    ReturnAll - return all points in the tree 
**
**      CArbArray <ArbIntNode*>*ReturnAll() const
**
**      Description: This function returns an array of all the points 
**          currently stored in the tree. The caller takes ownership of 
**          the array and must delete it. 
**
**      Return Value: An array of all points in the tree. 
**
**
** -- */

CArbArray<ArbIntNode*> *CArbKDTree2D::ReturnAll() const
{
    CArbArray<ArbIntNode*> *result = new CArbArray<ArbIntNode*>() ;
    ReturnAllRec(Root,result) ;
    return result;
}




// %(CArbKDTree2D::ReturnAllRec-void-|^const-KDTreeEntry-|*-CArbArray-|<ArbIntNode*>*) 
/* ++ ----------------------------------------------------------
**
**    ReturnAllRec - return all routine 
**
**      void ReturnAllRec(
**              KDTreeEntry             *entry,
**              CArbArray<ArbIntNode*>* result) const
**
**        entry  - (in)  tree node 
**        result - (i/o) array to which points are added 
**
**      Description: This function visits all tree nodes recursively 
**          and adds all points in the tree to the array. 
**
**
** -- */

void CArbKDTree2D::ReturnAllRec(KDTreeEntry *entry,
                                CArbArray<ArbIntNode*> *result) const
{
    for (int i=0; i<entry->num_pnts; i++)
        result->InsertAtEnd(entry->pnts[i]) ;
    if (entry->lchild) ReturnAllRec(entry->lchild,result) ;
    if (entry->rchild) ReturnAllRec(entry->rchild,result) ;
}




// %(CArbKDTree2D::SplitPointSet-void-|^const-ArbIntNode-|*^[]-int-|-ArbIntNode-|*-int-|&-ArbIntNode-|*^[]-int-|&-ArbIntNode-|*^[]-int-|&-ArbIntNode-|*^[]-int-|(*^)()) 
/* ++ ----------------------------------------------------------
**
**    SplitPointSet - split the point set into three bins 
**
**      void SplitPointSet(
**              ArbIntNode *y[],
**              int        n,
**              ArbIntNode *p,
**              int        &nL,
**              ArbIntNode *yL[],
**              int        &nR,
**              ArbIntNode *yR[],
**              int        &nS,
**              ArbIntNode *yS[],
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

void CArbKDTree2D::SplitPointSet(ArbIntNode *y[], int n, ArbIntNode *p,
                                 int &nL, ArbIntNode *yL[], 
                                 int &nR, ArbIntNode *yR[],
                                 int &nS, ArbIntNode *yS[],
                                 int (*cmp)(ArbIntNode*,ArbIntNode*)) const
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




// %(CArbKDTree2D::MergeSort-void-|^const-ArbIntNode-|*^[]-int-|-int-|(*^)())
/* ++ ----------------------------------------------------------
**
**    MergeSort - sort a list 
**
**      void MergeSort(
**              ArbIntNode *a[],
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

void CArbKDTree2D::MergeSort(ArbIntNode *a[], int n, 
                             int (*cmp)(ArbIntNode*,ArbIntNode*)) const
{
    MSort(a, 0, n-1, cmp) ;
}




// %(CArbKDTree2D::MSort-void-|^const-ArbIntNode-|*^[]-int-|-int-|-int-|(*^)())
/* ++ ----------------------------------------------------------
**
**    MSort - recursive merge sort 
**
**      void MSort(
**              ArbIntNode *a[],
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

void CArbKDTree2D::MSort(ArbIntNode *a[], int l, int r, 
                         int (*cmp)(ArbIntNode*,ArbIntNode*)) const
{
    if (l<r) {
        int m = (l+r)/2 ;
        MSort(a, l, m, cmp) ;
        MSort(a, m+1, r, cmp) ;
        Merge(a, l, m, r, cmp) ;
    }
}




// %(CArbKDTree2D::Merge-void-|^const-ArbIntNode-|*^[]-int-|-int-|-int-|-int-|(*^)()) 
/* ++ ----------------------------------------------------------
**
**    Merge - merge function 
**
**      void Merge(
**              ArbIntNode *x[],
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

void CArbKDTree2D::Merge(ArbIntNode *x[], int l, int m, int r, 
                         int (*cmp)(ArbIntNode*,ArbIntNode*)) const
{
    ArbIntNode **a = x+l ;
    ArbIntNode **b = x+m+1 ;
    ArbIntNode **c = new ArbIntNode*[r-l+1] ;
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


static int LeftToRightCmp(ArbIntNode *a, ArbIntNode *b)
{
    if ((a->coord[0] < b->coord[0]) || 
        ((a->coord[0] == b->coord[0]) && (a->coord[1] < b->coord[1]))) 
        return -1;
    if ((a->coord[0] > b->coord[0]) || 
        ((a->coord[0] == b->coord[0]) && (a->coord[1] > b->coord[1])))
        return 1;
    return 0;
}


static int BottomToTopCmp(ArbIntNode *a, ArbIntNode *b)
{
    if ((a->coord[1] < b->coord[1]) ||
        ((a->coord[1] == b->coord[1]) && (a->coord[0] < b->coord[0])))
        return -1;
    else if ((a->coord[1] > b->coord[1]) ||
             ((a->coord[1] == b->coord[1]) && (a->coord[0] > b->coord[0]))) 
        return 1;
    return 0;
}


