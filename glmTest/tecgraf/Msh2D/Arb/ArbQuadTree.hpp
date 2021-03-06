//
// CArbQuadTree Class header file
//
// Description -
//   This class implements a quadtree object.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.5 $  $Date: 2002/06/10 17:25:49 $  $Author: wash $
//


#ifndef ArbQuadTree_hh
#define ArbQuadTree_hh

#include <stddef.h>

class CArbQuadTree {

    public:

        // constructors, destructors, and member functions

        CArbQuadTree(double _OriginX,double _OriginY,double _size) ;
        ~CArbQuadTree() ;
        void RefineToSize(double x,double y,double size) ;
        void UniformRefine(double size,double minx, double maxx,
                           double miny,double maxy) ;
        void RefineOneLevelDiff() ;
        void VisitLevels(void *ptr,
            void ((*)(void *,double,double,double,bool,bool))) ;

        // This is an enumerated type used to define directions.
        // The code relies on the fact that NE,NE,SW,SE each
        // define one bit.  N,S,E, and W are not used that way
        // even though that is how they are defined.

        enum QTreeDir { ERR_DIR = 0, NW_DIR = 1, NE_DIR = 2,
                        SW_DIR = 4, SE_DIR = 8,
                        NORTH_DIR = 16, SOUTH_DIR = 32,
                        WEST_DIR = 64, EAST_DIR = 128 } ;

    private:

        // This is a nested class that is used to implement
        // nodes in the quad tree.  Note that the new and
        // delete operators are overloaded

        class CArbQTreeCell {
            public:
                CArbQTreeCell():
                    NW(0),NE(0),SW(0),SE(0),parent(0) {} ; 
                CArbQTreeCell(CArbQTreeCell *iparent):
                    NW(0),NE(0),SW(0),SE(0),parent(iparent) {} ; 
                ~CArbQTreeCell() {} ;
                CArbQTreeCell *GetNW() { return(NW) ; } ;
                CArbQTreeCell *GetNE() { return(NE) ; } ;
                CArbQTreeCell *GetSW() { return(SW) ; } ;
                CArbQTreeCell *GetSE() { return(SE) ; } ;
                CArbQTreeCell *GetParent() { return(parent) ; } ;
                enum QTreeDir ToParentVia() const ;
                void SetNW(CArbQTreeCell *cell) { NW = cell ; } ;
                void SetNE(CArbQTreeCell *cell) { NE = cell ; } ;
                void SetSW(CArbQTreeCell *cell) { SW = cell ; } ;
                void SetSE(CArbQTreeCell *cell) { SE = cell ; } ;
                bool IsLeaf() const ;
            private:
                CArbQTreeCell *NW, *NE, *SW, *SE ;
                CArbQTreeCell *parent ;
        } ;

        // private member functions

        void Divide(CArbQTreeCell *cell) ;
        CArbQTreeCell *FindCell(double x,double y) const ;
        CArbQTreeCell *FindNeighbor(CArbQTreeCell *cell,int level,
                                    QTreeDir dir,int *n_level) ;
        QTreeDir Mirror(int mirror_axis,int mirror_dir) const ;

        // these are private member functions which call themselfs
        // recursively

        CArbQTreeCell *FindCellRec(double x,double y,CArbQTreeCell *cell,
                           double *size,double *ox,double *oy) const ;
        void UniformRefineRec(CArbQTreeCell *cell,double cur,
                              double target,double ox,double oy,
                              double minx, double maxx,
                              double miny,double maxy) ;
        bool RefineOneLevelDiffRec(CArbQTreeCell *cell,int level) ;
        void RefineToSizeRec(double x,double y,double size,
                                CArbQTreeCell *cell,double cell_size,
                                double ox, double oy) ;
        void VisitLevelsRec(CArbQTreeCell *cell,double ox,double oy,
                 double size,void *ptr,
                 void((*)(void *,double,double,double,bool,bool))) ;
        void VisitRec(CArbQTreeCell *cell,int level) const ;
        void DeleteCellsRec(CArbQTreeCell *cell) ;

        // private member variables and constants

        CArbQTreeCell *Root ;
        double OriginX, OriginY ;
        double RootCellSize ;
        static const int EW_AXIS ;
        static const int NS_AXIS ;
        static const double FuzzyDivideTolerance ;

        // these member functions and variables are used to
        // manage a dynamic stack used by find_neighbor

        QTreeDir *NeighbStack ;
        int NeighbStackPtr ;
        int NeighbStackSize ;
        static const int NeighbStackQuantum ;
        void NeighbStackGrow() ;
} ;

/*
CLASS CArbQuadTree

  This object implements a quadtree data structure used divide a 2D 
  space. 


PUBLIC INTERFACE

  Public Member Functions:

    CArbQuadTree - constructor 

      CArbQuadTree(
              double _OriginX,
              double _OriginY,
              double _size)

        _OriginX - (in)  x coordinate of the origin 
        _OriginY - (in)  y coordinate of the origin 
        _size    - (in)  root cell size 

      Description: This is a constructor for an ArbQuadTree object. 


    CArbQuadTree - destructor 

      ~CArbQuadTree()

      Description: This is a destructor for an ArbQuadTree object. 


    RefineToSize - refine the tree localy 

      void RefineToSize(
              double x,
              double y,
              double size)

        x    - (in)  x coordinate of the refine point 
        y    - (in)  y coordinate of the refine point 
        size - (in)  refinement size 

      Description: This function refines the tree localy about point 
          so that the containing cell size is less than or equal to 
          the specified size. 


    UniformRefine - uniform tree refinement 

      void UniformRefine(double size)

        size - (in)  refinement size 

      Description: This function refines the tree so that all cell 
          sizes are less than or equal to the specified size. 


    RefineOneLevelDiff - refine so that neighbors do not differ by 
                         more than one level 

      void RefineOneLevelDiff()

      Description: This function refines the tree so that neighboring 
          cells in the tree do not differ by more than one level. 


    VisitLevels - visit all cells in the tree 

      void VisitLevels(
              void *ptr,
              void (*func)())

        ptr  - (in)  data to pass to the client function 
        func - (in)  client function 

      Description: This function visits all cells in the the tree and 
          calls the passed function function at each cell. The 
          signature of the passed function is 

        void (*func)(void *ptr,double OriginX,double OriginY,
                     double HalfSize,bool RootFlag,bool LeafFlag) -

            void *ptr - client data data
            double OriginX - x coordinate of the center of the cell
            double OriginY - y coordinate of the center of the cell
            double HalfSize - half size of the cell
            bool RootFlag - true if this is the root cell
            bool LeafFlg - true if this is a leaf cell
            

PRIVATE INTERFACE

  Private Member Functions:

    Divide - divide a cell 

      void Divide(CArbQTreeCell *cell)

        cell - (in)  cell to divide 

      Description: Divide a leaf cell to add four children. 


    FindCell - find a points containing cell 

      CArbQTreeCell *FindCell(
              double x,
              double y) const

        x - (in)  point's x coordinate 
        y - (in)  point's y coordinate 

      Description: Find the quadtree cell that contains the the 
          specified point. 

      Return Value: containing cell 


    FindNeighbor - find a cells neighbor 

      CArbQTreeCell *FindNeighbor(
              CArbQTreeCell *cell,
              int           level,
              QTreeDir      dir,
              int           *n_level)

        cell    - (in)  input cell 
        level   - (in)  input cell's tree level 
        dir     - (in)  search direction 
        n_level - (out) neighbor's tree level 

      Description: Given a cell and a search direction, this function 
          finds the neighboring cell in the specified direction. 

      Return Value: neighboring cell or zero if the input cell is on 
          the boundary of the tree 


    Mirror - find the mirror direction 

      QTreeDir Mirror(
              int mirror_axis,
              int mirror_dir) const

        mirror_axis - (in)  axis (either NorthSouth or EastWest) 
        mirror_dir  - (in)  direction (North, South, East, or West) 

      Description: given an axis and a direction finds the "mirror" 
          direction. 

      Return Value: the mirror direction 


    FindCellRec - find a cell 

      CArbQTreeCell *FindCellRec(
              double        x,
              double        y,
              CArbQTreeCell *cell,
              double        *size,
              double        *ox,
              double        *oy) const

        x    - (in)  point's x coordinate 
        y    - (in)  point's y coordinate 
        cell - (in)  current cell 
        size - (i/o) size of the current cell 
        ox   - (i/o) origin of the current cell 
        oy   - (i/o) origin of the current cell 

      Description: This function implements the recursive procedure 
          used to find the cell containing a point. 

      Return Value: containing cell 


    UniformRefineRec - uniform refinement 

      void UniformRefineRec(
              CArbQTreeCell *cell,
              double        cur,
              double        target)

        cell   - (in)  current cell 
        cur    - (in)  current cell size 
        target - (in)  target cell size 

      Description: This function implements the recursive procedure 
          used to do a uniform refinement of the tree. 


    RefineOneLevelDiffRec - one level difference refinement 

      bool RefineOneLevelDiffRec(
              CArbQTreeCell *cell,
              int           level)

        cell  - (in)  current cell 
        level - (in)  current level 

      Description: This function implements the recursive procedure 
          used to do a refinement so that neighboring cells do not 
          differ in level by more than one. 

      Return Value: true if any divisions were made. 


    RefineToSizeRec - refine a cell 

      void RefineToSizeRec(
              double        x,
              double        y,
              double        size,
              CArbQTreeCell *cell,
              double        cell_size,
              double        ox,
              double        oy)

        x         - (in)  point's x coordinate 
        y         - (in)  point's y coordinate 
        size      - (in)  refinement size 
        cell      - (in)  current cell 
        cell_size - (in)  current cell size 
        ox        - (in)  current cell center 
        oy        - (in)  current cell center 

      Description: This function implements the recursive procedure 
          used to do a refinement so that the cell containing a point 
          is less than or equal to a given size. 


    VisitLevelsRec - visit all cells wrapper 

      void VisitLevelsRec(
              CArbQTreeCell *cell,
              double        ox,
              double        oy,
              double        size,
              void          *ptr,
              void          (*func)())

        cell - (in)  current cell 
        ox   - (in)  cell center 
        oy   - (in)  cell center 
        size - (in)  cell size 
        ptr  - (in)  client data 
        func - (in)  visit function 

      Description: This function is a wrapper for the function that 
          implements the recursive procedure used to visit all cells 
          in the tree. 


    VisitRec - visit all cells 

      void VisitRec(
              CArbQTreeCell *cell,
              int           level) const

        cell  - (in)  current cell 
        level - (in)  current level 

      Description: This function implements the recursive procedure 
          used to visit all cells in the tree. 


    DeleteCellsRec - delete cells 

      void DeleteCellsRec(CArbQTreeCell *cell)

        cell - (in)  current cell 

      Description: Recursively delete this cell and all children. 


    NeighbStackGrow - grow the neighbor cell 

      void NeighbStackGrow()

      Description: This function grows the stack used to determin 
          neighbor cells. 


  Private Member Variables:

    CArbQTreeCell *Root - root cell 

    double OriginX - center of tree 

    double OriginY - center of tree 

    double RootCellSize - size of the root cell 

    static const int EW_AXIS - axis constant values 

    static const int NS_AXIS - axis constant values 

    static const double FuzzyDivideTolerance - This is used when 
            refining the tree for cases where the refinement point 
            is near the EW or NS axis. If it is within this 
            tolerance of the axis then both sides of the axis are 
            refined. 

    QTreeDir *NeighbStack - stack used to find a cell's neigbor 

    int NeighbStackPtr - pointer into the neighbor stack 

    int NeighbStackSize - current size of the neighbor stack 

    static const int NeighbStackQuantum - gowth quantum for the 
            neighbor stack. 


*/

#endif
