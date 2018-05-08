//
// CArbFeasableRegion Class header file
//
// Description -
//   This class implements a search for the feasable region
//   polygon due to a number of linear constraints
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2000,2001
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.3 $  $Date: 2001/01/09 14:46:34 $  $Author: wash $
//

#ifndef ArbFeasableRegion_hh
#define ArbFeasableRegion_hh

#include "ArbCoord2D.hpp"
#include "ArbArray.hpp"

class CArbFeasableRegion {
    public:

        CArbFeasableRegion(double tolerance) ;
        CArbFeasableRegion(const CArbFeasableRegion &other) ;
        CArbFeasableRegion operator = (const CArbFeasableRegion &other) ;
        ~CArbFeasableRegion() {} ;

        void AddConstraint(CArbCoord2D pt0,CArbCoord2D pt1) ;
        CArbArray<CArbCoord2D> *GetVertices() ;

        struct PolySide {
            PolySide(CArbCoord2D id,CArbCoord2D ie,double large) :
                d(id),e(ie),min(-large),max(large),prev(-1),next(-1) {} ;
            PolySide() {} ;
            CArbCoord2D d, e ;
            double min, max ;
            int prev, next ;
        } ;

    private:

        double Large, Tol ;
        CArbArray<PolySide> Sides ;
} ;

/*
CLASS CArbFeasableRegion

  This class finds the feasable region for a series of 2D linear 
  constraints. That is if you think about dividing up a 2D region with 
  a series of half spaces, this object finds the polygon that is left. 


PUBLIC INTERFACE

  Public Data Structures:

    struct PolySide

      This is a data structure used to keep track of the cutting half 
      spaces. It is public so that an array of these can be created. 


  Public Member Functions:

    PolySide - parameterized constructor 

      PolySide(
              CArbCoord2D id,
              CArbCoord2D ie,
              double      huge)

        id   - (in)  first point on the line 
        ie   - (in)  line direction vector 
        huge - (in)  relative infinity 

      Description: This is a constructor for the data structure 


    PolySide - no parameter constructor 

      PolySide()

      Description: This is a constructor with no arguments. 


  Public Member Variables:

    CArbCoord2D d - first point on the constraint line 

    CArbCoord2D e - direction vector for the line 

    double min - parametric coordinate of the greater 
        intersection point 

    double max - parametric coordinate of the lesser 
        intersection point 

    int prev - pointer to the clockwise constraint line on 
        polygon 

    int next - pointer to the counter clockwise constraint 
        on polygon 


PRIVATE INTERFACE

  Public Member Functions:

    CArbFeasableRegion - constructor 

      CArbFeasableRegion(double tolerance)

        tolerance - (in)  intersection tolerance 

      Description: This is a constructor for a FeasabilityRegion 


    CArbFeasableRegion - copy constructor 

      CArbFeasableRegion(const CArbFeasableRegion &other)

        other - (in)  object to copy 

      Description: This is the copy constructor for a 
          FeasabilityRegion 


    operator_= - assignment operator 

      CArbFeasableRegion operator = (const CArbFeasableRegion &other)

        other - (in)  object to copy 

      Description: This is the assignment operator for a 
          FeasabilityRegion 

      Return Value: the updated object 


    CArbFeasableRegion - destructor 

      ~CArbFeasableRegion()

      Description: This is the destructor for a FeasabilityRegion. 


    AddConstraint - add a linear constraint 

      void AddConstraint(
              CArbCoord2D pt0,
              CArbCoord2D pt1)

        pt0 - (in)  first point on constraint line 
        pt1 - (in)  second point on constraint line 

      Description: This function adds a linear constraint to the 
          feasable region. 


    GetVertices - get the polygon vertices 

      CArbArray <CArbCoord2D>*GetVertices()

      Description: This function returns the verticies of the 
          feasable region polygon. 

      Return Value: An array of polygon vertices. Ownership of the 
          array passes to the client, who must eventually call 
          delete. 


  Private Member Variables:

        double Huge - relative infinity 

        double Tol - intersection tolerance 

        CArbArray<PolySide> Sides - list of the linear constraints. 

*/
#endif

