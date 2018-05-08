//
// FeasableRegion Class header file
//
// Description -
//   This class implements a search for the feasable region
//   polygon due to a number of linear constraints
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2000,2001,2003
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.4 $  $Date: 2003/05/22 15:44:10 $  $Author: wash $
//

#ifndef FeasableRegion_hh
#define FeasableRegion_hh

#include "List.hpp"
#include "Vec2D.hpp"

using FTools::List ;
using FTools::Vec2D ;

namespace Msh2D {

class FeasableRegion {

    public:

        class VtxIterator ;

        FeasableRegion(double tolerance) ;
        FeasableRegion(const FeasableRegion &other) ;
        FeasableRegion operator = (const FeasableRegion &other) ;
        ~FeasableRegion() {} ;

        void AddConstraint(Vec2D pt0,Vec2D pt1) ;
        List<Vec2D> *GetVertices() ;

        VtxIterator GetVertexIterator() {
            return(VtxIterator(GetVertices())) ;
        } ;

        struct PolySide {
            PolySide(Vec2D id,Vec2D ie,double large) :
                d(id),e(ie),min(-large),max(large),prev(-1),next(-1) {} ;
            PolySide() {} ;
            Vec2D d, e ;
            double min, max ;
            int prev, next ;
        } ;

    private:

        double Large, Tol ;
        List<PolySide> Sides ;

    public:

        class VtxIterator {
            public:
                VtxIterator(List<Vec2D> *verts) :
                    vts(verts),cur(0) {} ;
                ~VtxIterator() { delete vts ; } ;
                void   First()     { cur = 0 ; } ;
                void   Next()      { ++cur ; } ;
                bool   More()      { return(cur < vts->Len()) ; } ;
                const Vec2D &Vtx() { return(vts->At(cur)) ; } ;

                void operator ++ ()      { Next() ; } ;
                void operator ++ (int i) { Next() ; } ;

                const Vec2D &operator * () { return(Vtx()) ; } ;

            private:
                List<Vec2D> *vts ;
                int cur ;
        } ;
} ;

/*
CLASS FeasableRegion

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
              Vec2D id,
              Vec2D ie,
              double      huge)

        id   - (in)  first point on the line 
        ie   - (in)  line direction vector 
        huge - (in)  relative infinity 

      Description: This is a constructor for the data structure 


    PolySide - no parameter constructor 

      PolySide()

      Description: This is a constructor with no arguments. 


  Public Member Variables:

    Vec2D d - first point on the constraint line 

    Vec2D e - direction vector for the line 

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

    FeasableRegion - constructor 

      FeasableRegion(double tolerance)

        tolerance - (in)  intersection tolerance 

      Description: This is a constructor for a FeasabilityRegion 


    FeasableRegion - copy constructor 

      FeasableRegion(const FeasableRegion &other)

        other - (in)  object to copy 

      Description: This is the copy constructor for a 
          FeasabilityRegion 


    operator_= - assignment operator 

      FeasableRegion operator = (const FeasableRegion &other)

        other - (in)  object to copy 

      Description: This is the assignment operator for a 
          FeasabilityRegion 

      Return Value: the updated object 


    FeasableRegion - destructor 

      ~FeasableRegion()

      Description: This is the destructor for a FeasabilityRegion. 


    AddConstraint - add a linear constraint 

      void AddConstraint(
              Vec2D pt0,
              Vec2D pt1)

        pt0 - (in)  first point on constraint line 
        pt1 - (in)  second point on constraint line 

      Description: This function adds a linear constraint to the 
          feasable region. 


    GetVertices - get the polygon vertices 

      List <Vec2D>*GetVertices()

      Description: This function returns the verticies of the 
          feasable region polygon. 

      Return Value: An array of polygon vertices. Ownership of the 
          array passes to the client, who must eventually call 
          delete. 


  Private Member Variables:

        double Huge - relative infinity 

        double Tol - intersection tolerance 

        List<PolySide> Sides - list of the linear constraints. 

*/

} // namespace

#endif

