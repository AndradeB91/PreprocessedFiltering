//
// CArbRectangle and AxisParallelEdge class header file
//
// Description -
//   This class implements an object that will create a 
//   rectangle and a rectangle with edges parallel to 
//   the Cartesian axes
//
// Reference - 
//   Computational Geometry and Computer Graphics in C++
//   Michael J. Laszlo, Prentice Hall, 1996
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000
//   All rights reserved
//
// Author -
//   Bruce Carter and Wash Wawrzynek
//
// Revision -
//   $Revision: 1.4 $  $Date: 2001/01/08 20:28:29 $  $Author: wash $
//

#ifndef ArbRectangle_h
#define ArbRectangle_h

#include "ArbMsh.hpp"

class CArbRectangle {
    public:
        CArbRectangle(const ArbIntNode &p0, const ArbIntNode &p1,
                      const double &size) {
            double xc, yc ;
            xc = 0.5*(p0.coord[0]+p1.coord[0]);
            yc = 0.5*(p0.coord[1]+p1.coord[1]);
            sw.coord[0] = xc - size ;
            sw.coord[1] = yc - size ;
            ne.coord[0] = xc + size ;
            ne.coord[1] = yc + size ; } ;
        CArbRectangle() {} ;

        bool PointInRectangle(const ArbIntNode* p) const {
            return ((sw.coord[0] <= p->coord[0]) &&
                    (p->coord[0] <= ne.coord[0]) &&
                    (sw.coord[1] <= p->coord[1]) &&
                    (p->coord[1] <= ne.coord[1])); } ;

        ArbIntNode *SW() { return(&sw) ; } ;
        ArbIntNode *NE() { return(&ne) ; } ;

    private:
        ArbIntNode sw;
        ArbIntNode ne;
};

/*
CLASS CArbRectangle

  rectangle class for a KD tree 


PUBLIC INTERFACE

  Public Member Functions:

    CArbRectangle - constructor 

      CArbRectangle(
              const ArbIntNode &p0,
              const ArbIntNode &p1,
              const double     &size)

        p0   - (in)  first corner point 
        p1   - (in)  second corner point 
        size - (in)  rectangle size 

      Description: Constructor for the rectangle class 


    CArbRectangle - constructor 

      CArbRectangle()

      Description: No argument constructor for a rectangle. 


    PointInRectangle - check for point inclusion 

      bool PointInRectangle(const ArbIntNode *p) const

        p - (in)  node point 

      Description: Check to see if a point is included in the 
          rectangle 

      Return Value: true if the rectangle contains the point 


    SW - return the southwest corner 

      ArbIntNode *SW()

      Description: Return a pointer to the node that is the southwest 
          corner of the rectangle. 

      Return Value: southwest corner 


    NE - return the northeast corner 

      ArbIntNode *NE()

      Description: Return a pointer to the node that is the northeast 
          corner of the rectangle. 

      Return Value: northeast corner 


PRIVATE INTERFACE

  Private Member Variables:

    ArbIntNode sw - southwest corner 

    ArbIntNode ne - northeast corner 

*/

#endif
