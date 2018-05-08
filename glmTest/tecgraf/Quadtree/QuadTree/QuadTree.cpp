//
// QuadTree Class definition
//
// Description -
//   This class implements a quad tree data type.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2007
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//

#include <cstdlib>
#include <cstdio>

#include "QuadTree.hpp"

#ifdef MEMDEBUG
#include "MemDbg.hpp"
//#define new new(__FILE__,__LINE__)
#endif

namespace FTools {

/* -----------------------------------------------------------
    EW_AXIS,NS_AXIS - static member variable.  These are
                      declared as integers rather than
                      an enumerated type because they will
                      be used in bit operations and I'm not
                      sure that all compilers would be happy
                      about that. 
*/

const int QuadTree::EW_AXIS = 1 ;
const int QuadTree::NS_AXIS = 2 ;


/* -----------------------------------------------------------
    FuzzyDivideTolerance - static member variable.  This is used
                           when refining the tree for cases
                           where the refinement point is near
                           the EW or NS axis.  If it is within
                           this tolerance of the axis then both
                           sides of the axis are refined.
*/

const double QuadTree::FuzzyDivideTolerance = 0.05 ;


/* -----------------------------------------------------------
    NeighbStackQuantum - static member variable.  This is
                         the amount by which the neighbor stack
                         will grow.
*/

const int QuadTree::NeighbStackQuantum = 10 ;



// %(QuadTree::QuadTree-constructor-|-double-|-double-|-double-|)
/* ++ ----------------------------------------------------------
**
**    QuadTree - constructor 
**
**      QuadTree(
**              double _OriginX,
**              double _OriginY,
**              double _size)
**
**        _OriginX - (in)  x coordinate of the origin 
**        _OriginY - (in)  y coordinate of the origin 
**        _size    - (in)  root cell size 
**
**      Description: This is a constructor for an QuadTree object. 
**
**
** -- */

QuadTree::QuadTree(double iOriginX,double iOriginY,double iSize)
{
    OriginX = iOriginX ;
    OriginY = iOriginY ;
    RootCellSize = iSize ;
    NeighbStack = 0 ;
    NeighbStackPtr = 0 ;
    NeighbStackSize = 0 ;
    Root = new QTreeCell(0) ;
}




// %(QuadTree::QuadTree-destructor-|~) 
/* ++ ----------------------------------------------------------
**
**    QuadTree - destructor 
**
**      ~QuadTree()
**
**      Description: This is a destructor for an QuadTree object. 
**
**
** -- */

QuadTree::~QuadTree()
{
    // delete all the cells

    DeleteCellsRec(Root) ;
    delete Root ;  Root = 0 ;
    if (NeighbStack != 0) delete [] NeighbStack ;
}




// %(QuadTree::Divide-void-|-QTreeCell-|*) 
/* ++ ----------------------------------------------------------
**
**    Divide - divide a cell 
**
**      void Divide(QTreeCell *cell)
**
**        cell - (in)  cell to divide 
**
**      Description: Divide a leaf cell to add four children. 
**
**
** -- */

void QuadTree::Divide(QTreeCell *cell)
{
    if (cell->GetNW() == 0) cell->SetNW(new QTreeCell(cell)) ;
    if (cell->GetNE() == 0) cell->SetNE(new QTreeCell(cell)) ;
    if (cell->GetSW() == 0) cell->SetSW(new QTreeCell(cell)) ;
    if (cell->GetSE() == 0) cell->SetSE(new QTreeCell(cell)) ;
}




// %(QuadTree::FindCell-QTreeCell-|*^const-double-|-double-|) 
/* ++ ----------------------------------------------------------
**
**    FindCell - find a points containing cell 
**
**      QTreeCell *FindCell(
**              double x,
**              double y) const
**
**        x - (in)  point's x coordinate 
**        y - (in)  point's y coordinate 
**
**      Description: Find the quadtree cell that contains the the 
**          specified point. 
**
**      Return Value: containing cell 
**
**
** -- */

QuadTree::QTreeCell *QuadTree::FindCell(double x,double y) const
{
    double cell_size = RootCellSize ;
    double ox = OriginX ;
    double oy = OriginY ;
    return(FindCellRec(x,y,Root,&cell_size,&ox,&oy)) ;
}




// %(QuadTree::RefineToSize-void-|-double-|-double-|-double-|) 
/* ++ ----------------------------------------------------------
**
**    RefineToSize - refine the tree localy 
**
**      void RefineToSize(
**              double x,
**              double y,
**              double size)
**
**        x    - (in)  x coordinate of the refine point 
**        y    - (in)  y coordinate of the refine point 
**        size - (in)  refinement size 
**
**      Description: This function refines the tree localy about point 
**          so that the containing cell size is less than or equal to 
**          the specified size. 
**
**
** -- */

void QuadTree::RefineToSize(double x,double y,double size)
{
    RefineToSizeRec(x,y,size,Root,RootCellSize,OriginX,OriginY) ;
}


double QuadTree::ContainingCellSize(double x,double y) const
{
    double cell_size = 0.5*RootCellSize ;
    double ox = OriginX ;
    double oy = OriginY ;
    FindCellRec(x,y,Root,&cell_size,&ox,&oy) ;
    return(2*cell_size) ;
}


// %(QuadTree::RefineToSizeRec-void-|-double-|-double-|-double-|-QTreeCell-|*-double-|-double-|-double-|)
/* ++ ----------------------------------------------------------
**
**    RefineToSizeRec - refine a cell 
**
**      void RefineToSizeRec(
**              double        x,
**              double        y,
**              double        size,
**              QTreeCell *cell,
**              double        cell_size,
**              double        ox,
**              double        oy)
**
**        x         - (in)  point's x coordinate 
**        y         - (in)  point's y coordinate 
**        size      - (in)  refinement size 
**        cell      - (in)  current cell 
**        cell_size - (in)  current cell size 
**        ox        - (in)  current cell center 
**        oy        - (in)  current cell center 
**
**      Description: This function implements the recursive procedure 
**          used to do a refinement so that the cell containing a point 
**          is less than or equal to a given size. 
**
**
** -- */

void QuadTree::RefineToSizeRec(double x,double y,double size,
                                   QTreeCell *cell,double cell_size,
                                   double ox,double oy)
{
    if (cell_size <= size) return ;

    // now refine until the cell size is less than the input size.
    // Note that there is a tolerance to give us a "fuzzy" Divide
    // which allows more than one cell may be Divided.  We need
    // this when a point is near the edge of a cell.

    double tolerance = cell_size * FuzzyDivideTolerance ;
    if (tolerance > size) tolerance = size ;
    Divide(cell) ;
    cell_size *= 0.5 ;

    bool do_NW = true ;
    bool do_NE = true ;
    bool do_SW = true ;
    bool do_SE = true ;

    // check which cell(s) should be Divided

    if (x >= (ox + tolerance)) {
        do_NW = false ;   do_SW = false ;
    } else if (x <= (ox - tolerance)) {
        do_NE = false ;   do_SE = false ;
    }

    if (y >= (oy + tolerance)) {
        do_SW = false ;   do_SE = false ;
    } else if (y <= (oy - tolerance)) {
        do_NW = false ;   do_NE = false ;
    }

    // Divide all marked cells

    double half_size = 0.5*cell_size ;
    if (do_NW == true)
        RefineToSizeRec(x,y,size,cell->GetNW(),cell_size,
                        ox-half_size,oy+half_size) ;

    if (do_NE == true)
        RefineToSizeRec(x,y,size,cell->GetNE(),cell_size,
                        ox+half_size,oy+half_size) ;

    if (do_SW == true)
        RefineToSizeRec(x,y,size,cell->GetSW(),cell_size,
                        ox-half_size,oy-half_size) ;

    if (do_SE == true)
        RefineToSizeRec(x,y,size,cell->GetSE(),cell_size,
                        ox+half_size,oy-half_size) ;
}




// %(QuadTree::UniformRefine-void-|-double-|)
/* ++ ----------------------------------------------------------
**
**    UniformRefine - uniform tree refinement 
**
**      void UniformRefine(double size)
**
**        size - (in)  refinement size 
**
**      Description: This function refines the tree so that all cell 
**          sizes are less than or equal to the specified size. 
**
**
** -- */

void QuadTree::UniformRefine(double size,
                                 double minx, double maxx,
                                 double miny,double maxy)
{
    UniformRefineRec(Root,RootCellSize,size,OriginX,OriginY,
                     minx,maxx,miny,maxy) ;
}




// %(QuadTree::FindCellRec-QTreeCell-|*^const-double-|-double-|-QTreeCell-|*-double-|*-double-|*-double-|*) 
/* ++ ----------------------------------------------------------
**
**    FindCellRec - find a cell 
**
**      QTreeCell *FindCellRec(
**              double        x,
**              double        y,
**              QTreeCell *cell,
**              double        *size,
**              double        *ox,
**              double        *oy) const
**
**        x    - (in)  point's x coordinate 
**        y    - (in)  point's y coordinate 
**        cell - (in)  current cell 
**        size - (i/o) size of the current cell 
**        ox   - (i/o) origin of the current cell 
**        oy   - (i/o) origin of the current cell 
**
**      Description: This function implements the recursive procedure 
**          used to find the cell containing a point. 
**
**      Return Value: containing cell 
**
**
** -- */

QuadTree::QTreeCell *QuadTree::FindCellRec(double x,double y,
             QTreeCell *cell,double *size,double *ox,double *oy) const
{
    // if this is a leaf cell then return

    if (cell->IsLeaf() == true) return(cell) ;

    // determine which quadrant the point is in and call ourself
    // recursively.

    *size = 0.5 * *size ;

//    if (x >= *ox + *size) {
    if (x >= *ox) {
//        if (y >= *oy + *size) {
        if (y >= *oy) {
            *ox += *size ;
            *oy += *size ;
            return(FindCellRec(x,y,cell->GetNE(),size,ox,oy)) ;
        } else {
            *ox += *size ;
            *oy -= *size ;
            return(FindCellRec(x,y,cell->GetSE(),size,ox,oy)) ;
        }
    } else {
//        if (y >= *oy + *size) {
        if (y >= *oy) {
            *ox -= *size ;
            *oy += *size ;
            return(FindCellRec(x,y,cell->GetNW(),size,ox,oy)) ;
        } else {
            *ox -= *size ;
            *oy -= *size ;
            return(FindCellRec(x,y,cell->GetSW(),size,ox,oy)) ;
        }
    }
}



// %(QuadTree::FindNeighbor-QTreeCell-|*-QTreeCell-|*-int-|-QTreeDir-|-int-|*)
/* ++ ----------------------------------------------------------
**
**    FindNeighbor - find a cells neighbor 
**
**      QTreeCell *FindNeighbor(
**              QTreeCell *cell,
**              int           level,
**              QTreeDir      dir,
**              int           *n_level)
**
**        cell    - (in)  input cell 
**        level   - (in)  input cell's tree level 
**        dir     - (in)  search direction 
**        n_level - (out) neighbor's tree level 
**
**      Description: Given a cell and a search direction, this function 
**          finds the neighboring cell in the specified direction. 
**
**      Return Value: neighboring cell or zero if the input cell is on 
**          the boundary of the tree 
**
**
** -- */

QuadTree::QTreeCell *QuadTree::FindNeighbor(
                        QTreeCell *cell,
                        int level,QTreeDir dir,int *n_level)
{
    int child_flg = 0, mirror_axis = 0, levels_up ;

    // first check to see if we have reached the Root cell
    // if so, return null as there is no neighbor inside the tree

    if (cell == Root) return(0) ;

    // set a few variables based on the direction

    switch (dir) {
        case NORTH_DIR:
            child_flg = SW_DIR | SE_DIR ;
            mirror_axis = EW_AXIS ;
            break ;
        case SOUTH_DIR:
            child_flg = NW_DIR | NE_DIR ;
            mirror_axis = EW_AXIS ;
            break ;
        case EAST_DIR:
            child_flg = SW_DIR | NW_DIR ;
            mirror_axis = NS_AXIS ;
            break ;
        case WEST_DIR:
            child_flg = NE_DIR | SE_DIR ;
            mirror_axis = NS_AXIS ;
            break ;
        default:
            break ;
    }

    // traverse up the tree from the given cell until we find
    // the nearest common ancestor.  This is the first cell
    // which is reached via the child flg

    QTreeCell *cur = cell->GetParent() ;
    QTreeDir  cdir = cell->ToParentVia() ;
    NeighbStackPtr = 0 ;
    if (NeighbStackPtr == NeighbStackSize) NeighbStackGrow() ;
    NeighbStack[NeighbStackPtr] = cdir ; ++NeighbStackPtr ;

    while ((cdir & child_flg) == 0) {
        if (cur == Root) return(0) ;
        cdir = cur->ToParentVia() ;
        cur = cur->GetParent() ;
        if (NeighbStackPtr == NeighbStackSize) NeighbStackGrow() ;
        NeighbStack[NeighbStackPtr] = cdir ; ++NeighbStackPtr ;
    }
    levels_up = NeighbStackPtr ;

    // now we have a common ancestor, so traverse back down the
    // the path, but take the mirror image of the directions
    // that we took on the way up

    QTreeDir  mdir ;

    while ((NeighbStackPtr > 0) && (cur->IsLeaf() == false)) {
        --NeighbStackPtr ;
        cdir = NeighbStack[NeighbStackPtr] ;
        mdir = Mirror(mirror_axis,cdir) ;
        switch (mdir) {
            case NW_DIR: cur = cur->GetNW() ; break ;
            case NE_DIR: cur = cur->GetNE() ; break ;
            case SW_DIR: cur = cur->GetSW() ; break ;
            case SE_DIR: cur = cur->GetSE() ; break ;
            default: break ;
        }
        --levels_up ;
    }

    *n_level = level - levels_up ;
    return(cur) ;
}




// %(QuadTree::Mirror-QTreeDir-|^const-int-|-int-|)
/* ++ ----------------------------------------------------------
**
**    Mirror - find the mirror direction 
**
**      QTreeDir Mirror(
**              int mirror_axis,
**              int mirror_dir) const
**
**        mirror_axis - (in)  axis (either NorthSouth or EastWest) 
**        mirror_dir  - (in)  direction (North, South, East, or West) 
**
**      Description: given an axis and a direction finds the "mirror" 
**          direction. 
**
**      Return Value: the mirror direction 
**
**
** -- */

QuadTree::QTreeDir QuadTree::Mirror(
           int mirror_axis,int mirror_dir) const
{
    if (mirror_axis == NS_AXIS) {
        if ((mirror_dir & NW_DIR) != 0) return(NE_DIR) ;
        if ((mirror_dir & NE_DIR) != 0) return(NW_DIR) ;
        if ((mirror_dir & SW_DIR) != 0) return(SE_DIR) ;
        if ((mirror_dir & SE_DIR) != 0) return(SW_DIR) ;
    } else if (mirror_axis == EW_AXIS) {
        if ((mirror_dir & NW_DIR) != 0) return(SW_DIR) ;
        if ((mirror_dir & NE_DIR) != 0) return(SE_DIR) ;
        if ((mirror_dir & SW_DIR) != 0) return(NW_DIR) ;
        if ((mirror_dir & SE_DIR) != 0) return(NE_DIR) ;
    }
    return(ERR_DIR) ;
}




// %(QuadTree::RefineOneLevelDiff-void-|) 
/* ++ ----------------------------------------------------------
**
**    RefineOneLevelDiff - refine so that neighbors do not differ by 
**                         more than one level 
**
**      void RefineOneLevelDiff()
**
**      Description: This function refines the tree so that neighboring 
**          cells in the tree do not differ by more than one level. 
**
**
** -- */

void QuadTree::RefineOneLevelDiff()
{
    bool changed = RefineOneLevelDiffRec(Root,1) ;
    while (changed) {
        changed = RefineOneLevelDiffRec(Root,1) ;
    }
}




// %(QuadTree::UniformRefineRec-void-|-QTreeCell-|*-double-|-double-|)
/* ++ ----------------------------------------------------------
**
**    UniformRefineRec - uniform refinement 
**
**      void UniformRefineRec(
**              QTreeCell *cell,
**              double        cur,
**              double        target)
**
**        cell   - (in)  current cell 
**        cur    - (in)  current cell size 
**        target - (in)  target cell size 
**
**      Description: This function implements the recursive procedure 
**          used to do a uniform refinement of the tree. 
**
**
** -- */

void QuadTree::UniformRefineRec(QTreeCell *cell,double cur,
                                    double target,double ox,double oy,
                                    double minx,double maxx,
                                    double miny,double maxy)
{
    // Divide all cells until they are at least as small as size

    if (cur <= target) return ;

    // if we do not have children, Divide

    if (cell->GetNW() == 0) Divide(cell) ;

    // now recursively examine our children

    double half = cur * 0.5 ;
    double qtr = 0.5 * half ;

    if (((oy+half) >= miny) && (oy <= maxy)) {
        if ((ox >= minx) && ((ox-half) <= maxx)) {
            UniformRefineRec(cell->GetNW(),half,target,ox-qtr,oy+qtr,
                             minx,maxx,miny,maxy) ;
        }
        if (((ox+half) >= minx) && (ox <= maxx)) {
            UniformRefineRec(cell->GetNE(),half,target,ox+qtr,oy+qtr,
                             minx,maxx,miny,maxy) ;
        }
    }
    if ((oy >= miny) && ((oy-half) <= maxy)) {
        if ((ox >= minx) && ((ox-half) <= maxx)) {
            UniformRefineRec(cell->GetSW(),half,target,ox-qtr,oy-qtr,
                             minx,maxx,miny,maxy) ;
        }
        if (((ox+half) >= minx) && (ox <= maxx)) {
            UniformRefineRec(cell->GetSE(),half,target,ox+qtr,oy-qtr,
                             minx,maxx,miny,maxy) ;
        }
    }
}




// %(QuadTree::RefineOneLevelDiffRec-bool-|-QTreeCell-|*-int-|) 
/* ++ ----------------------------------------------------------
**
**    RefineOneLevelDiffRec - one level difference refinement 
**
**      bool RefineOneLevelDiffRec(
**              QTreeCell *cell,
**              int           level)
**
**        cell  - (in)  current cell 
**        level - (in)  current level 
**
**      Description: This function implements the recursive procedure 
**          used to do a refinement so that neighboring cells do not 
**          differ in level by more than one. 
**
**      Return Value: true if any divisions were made. 
**
**
** -- */

bool QuadTree::RefineOneLevelDiffRec(QTreeCell *cell,int level)
{
    bool changed = false ;

    // if this is not a leaf cell, call ourself recursively
    // for all children

    if (cell->IsLeaf() == false) {
        if (RefineOneLevelDiffRec(cell->GetNW(),level+1) ||
            RefineOneLevelDiffRec(cell->GetNE(),level+1) ||
            RefineOneLevelDiffRec(cell->GetSW(),level+1) ||
            RefineOneLevelDiffRec(cell->GetSE(),level+1)) {
            changed = true ;
        }

    // else we are at a leaf cell, so find our neighbors and
    // refind if necessary

    } else {
        int Nlevel ;
        QTreeCell *Nnbr ;

        Nnbr = FindNeighbor(cell,level,NORTH_DIR,&Nlevel) ;
        if (Nnbr != 0) {
            if ((Nnbr->IsLeaf() == true) &&
                ((level - Nlevel) > 1)) {
                Divide(Nnbr) ;
                changed = true ;
                (void)RefineOneLevelDiffRec(Nnbr,Nlevel+1) ;
            }
        }

        Nnbr = FindNeighbor(cell,level,SOUTH_DIR,&Nlevel) ;
        if (Nnbr != 0) {
            if ((Nnbr->IsLeaf() == true) &&
                ((level - Nlevel) > 1)) {
                Divide(Nnbr) ;
                changed = true ;
                (void)RefineOneLevelDiffRec(Nnbr,Nlevel+1) ;
            }
        }

        Nnbr = FindNeighbor(cell,level,WEST_DIR,&Nlevel) ;
        if (Nnbr != 0) {
            if ((Nnbr->IsLeaf() == true) &&
                ((level - Nlevel) > 1)) {
                Divide(Nnbr) ;
                changed = true ;
                (void)RefineOneLevelDiffRec(Nnbr,Nlevel+1) ;
            }
        }

        Nnbr = FindNeighbor(cell,level,EAST_DIR,&Nlevel) ;
        if (Nnbr != 0) {
            if ((Nnbr->IsLeaf() == true) &&
                ((level - Nlevel) > 1)) {
                Divide(Nnbr) ;
                changed = true ;
                (void)RefineOneLevelDiffRec(Nnbr,Nlevel+1) ;
            }
        }
    }
    return(changed) ;
}




// %(QuadTree::VisitLevels-void-|-void-|*-void-|(*^)()) 
/* ++ ----------------------------------------------------------
**
**    VisitLevels - visit all cells in the tree 
**
**      void VisitLevels(
**              void *ptr,
**              void (*func)())
**
**        ptr  - (in)  data to pass to the client function 
**        func - (in)  client function 
**
**      Description: This function visits all cells in the the tree and 
**          calls the passed function function at each cell. The 
**          signature of the passed function is 
**
**          void (*func)(void *ptr,double OriginX,double OriginY,
**                     double HalfSize,bool RootFlag,bool LeafFlag) -
**
**            void *ptr - client data data
**            double OriginX - x coordinate of the center of the cell
**            double OriginY - y coordinate of the center of the cell
**            double HalfSize - half size of the cell
**            bool RootFlag - true if this is the root cell
**            bool LeafFlg - true if this is a leaf cell
**           
**
**
** -- */

void QuadTree::VisitLevels(void *ptr,
         void((*func)(void *,double,double,double,bool,bool)))
{
    bool is_leaf = (Root->GetNW() == 0) ? true : false ;
    (*func)(ptr,OriginX,OriginY,RootCellSize/2,true,is_leaf) ;
    if (is_leaf) return ;
    VisitLevelsRec(Root,OriginX,OriginY,RootCellSize/2,
                     ptr,func) ;
}




// %(QuadTree::VisitRec-void-|^const-QTreeCell-|*-int-|)
/* ++ ----------------------------------------------------------
**
**    VisitRec - visit all cells 
**
**      void VisitRec(
**              QTreeCell *cell,
**              int           level) const
**
**        cell  - (in)  current cell 
**        level - (in)  current level 
**
**      Description: This function implements the recursive procedure 
**          used to visit all cells in the tree. 
**
**
** -- */

void QuadTree::VisitLevelsRec(QTreeCell *cell,
         double ox,double oy,double half,void *ptr,
         void((*func)(void *,double,double,double,bool,bool)))
{
    (*func)(ptr,ox,oy,half,false,cell->IsLeaf()) ;
    if (cell->IsLeaf()) return ;

    double hhalf = half / 2 ;
    VisitLevelsRec(cell->GetNW(),ox-hhalf,oy+hhalf,hhalf,ptr,func) ;
    VisitLevelsRec(cell->GetNE(),ox+hhalf,oy+hhalf,hhalf,ptr,func) ;
    VisitLevelsRec(cell->GetSW(),ox-hhalf,oy-hhalf,hhalf,ptr,func) ;
    VisitLevelsRec(cell->GetSE(),ox+hhalf,oy-hhalf,hhalf,ptr,func) ;
}




// %(QuadTree::NeighbStackGrow-void-|) 
/* ++ ----------------------------------------------------------
**
**    NeighbStackGrow - grow the neighbor cell 
**
**      void NeighbStackGrow()
**
**      Description: This function grows the stack used to determin 
**          neighbor cells. 
**
**
** -- */

void QuadTree::NeighbStackGrow()
{
    if (NeighbStack == 0) {
        NeighbStack = new QTreeDir[NeighbStackQuantum] ;
        NeighbStackSize = NeighbStackQuantum ;
    } else {
        QTreeDir *tmp ;
        tmp = new QTreeDir[NeighbStackSize+NeighbStackQuantum] ;
        for (int i=0 ; i<NeighbStackSize ; ++i) {
            tmp[i] = NeighbStack[i] ;
        }
        delete [] NeighbStack ;
        NeighbStack = tmp ;
        NeighbStackSize += NeighbStackQuantum ;
    }
}




// %(QuadTree::DeleteCellsRec-void-|-QTreeCell-|*)
/* ++ ----------------------------------------------------------
**
**    DeleteCellsRec - delete cells 
**
**      void DeleteCellsRec(QTreeCell *cell)
**
**        cell - (in)  current cell 
**
**      Description: Recursively delete this cell and all children. 
**
**
** -- */

void QuadTree::DeleteCellsRec(QTreeCell *cell)
{
    if (cell->GetNW() == 0) return ;

    DeleteCellsRec(cell->GetNW()) ;  delete cell->GetNW() ;
    DeleteCellsRec(cell->GetNE()) ;  delete cell->GetNE() ;
    DeleteCellsRec(cell->GetSW()) ;  delete cell->GetSW() ;
    DeleteCellsRec(cell->GetSE()) ;  delete cell->GetSE() ;
}




// %(QTreeCell::ToParentVia-QTreeDir-enum|^const)
/* ++ ----------------------------------------------------------
**
**    ToParentVia - direction from parent 
**
**      enum QTreeDir ToParentVia() const
**
**      Description: This function returns the direction along which 
**          one would need to move from the parent to this cell. 
**
**      Return Value: the direction from the parent to this cell 
**
**
** -- */

enum QuadTree::QTreeDir QuadTree::QTreeCell::ToParentVia() const
{
    if (parent == 0) return(ERR_DIR) ;
    if (this == parent->NW) return(NW_DIR) ;
    if (this == parent->NE) return(NE_DIR) ;
    if (this == parent->SW) return(SW_DIR) ;
    if (this == parent->SE) return(SE_DIR) ;
    return(ERR_DIR) ;
}




// %(QTreeCell::IsLeaf-bool-|^const) 
/* ++ ----------------------------------------------------------
**
**    IsLeaf - leaf node query 
**
**      bool IsLeaf() const
**
**      Description: This function returns true if the cell is a leaf 
**          node. 
**
**      Return Value: true if a leaf node 
**
**
** -- */

inline bool QuadTree::QTreeCell::IsLeaf() const
{
    return(NW == 0 ? true : false) ;
}

} // namespace
