//
// CArbBSpline Template Class header file
//
// Description -
//   This class implements a 2D cubic b-spline class.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2000
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.6 $  $Date: 2001/08/08 19:33:22 $  $Author: wash $
//

#ifndef ArbBSpline2D_hh
#define ArbBSpline2D_hh

#include "ArbCoord2D.hpp"
#include "ArbArray.hpp"

class CArbBSpline2D {
    public:

        CArbBSpline2D(int num,CArbCoord2D *points) ;
        CArbBSpline2D(const CArbBSpline2D &other) ;
        CArbBSpline2D operator = (const CArbBSpline2D &other) ;
        ~CArbBSpline2D() ;

        CArbCoord2D Evaluate(const double nu) const ;
        void Derivatives(const double nu,CArbCoord2D *value,
                         CArbCoord2D *der_u,CArbCoord2D *der_u2) const ;

        CArbArray<double> *IntersectLine(const CArbCoord2D &pt0,
                                         const CArbCoord2D &pt1) const ;

        bool ClosestPoint(const CArbCoord2D &pt,double *u) const ;

        CArbBSpline2D *SubSegment(double par_start,
                                  double par_stop) const ;

        double ApproxLength(const int npts = 100) const ;

        void UniformParameterization(const int npts = 10000) ;

        int GetNumSegs() { return(Num_Segs) ; } ;

        void Reverse() ;

        void Print() ;

    private:

        int Num_Segs ;
        CArbCoord2D *Coord ;
        double *U ;

        void GetControlNet(int ns,CArbCoord2D *points) ;

        void FindCross(double u0,double u1,const CArbCoord2D &norm,
                       const CArbCoord2D &tang,const CArbCoord2D &pt0,
                       double coef,
                       CArbArray<double> **results) const ;
} ;

/*
CLASS CArbBSpline2D

  This class implements a package for 2D cubic B-Splines. 


PUBLIC INTERFACE

  Public Member Functions:

    CArbBSpline2D - Build a spline interpolates the input points 

      CArbBSpline2D(
              int         num,
              CArbCoord2D *points)

        num    - (in)  number of interpolation points 
        points - (in)  list of interpolation points 

      Description: This constructor generates a spline that 
          interpolates the given input points. This is a two-step 
          process. First, a cardinal spline is fit through the points 
          (a unit parametric step between each input point). In a 
          second step, the spline is reparameterized so that the 
          mapping between spline parametric space and the arc length 
          in Cartesian space is nearly constant along the spline. 


    CArbBSpline2D - copy constructor 

      CArbBSpline2D(const CArbBSpline2D &other)

        other - (in)  object to copy 

      Description: This is the copy constructor for the class. 


    operator_= - assignment operator 

      CArbBSpline2D operator = (const CArbBSpline2D &other)

        other - (in)  object to copy for the assignment 

      Description: This function implements an assignment operator 
          for this object. 

      Return Value: returns the updated object 


    CArbBSpline2D - destructor 

      ~CArbBSpline2D()

      Description: this is the destructor function 


    Evaluate - evaluate the spline value 

      CArbCoord2D Evaluate(const double nu) const

        nu - (in)  evaluation coordinate 

      Description: This function evaluates the spline value at the 
          given parametric coordinate. 

      Return Value: The spline value at this point. 


    Derivatives - compute spline derivative 

      void Derivatives(
              const double nu,
              CArbCoord2D  *value,
              CArbCoord2D  *der_u,
              CArbCoord2D  *der_u2) const

        nu     - (in)  evaluation coordinate 
        value  - (in)  spline value 
        der_u  - (in)  first derivative 
        der_u2 - (in)  second derivative 

      Description: This functions evaluates the spline's value and 
          its first and second derivatives with respect to the 
          parametric coordinate at the specified parametric 
          coordinate. 


    IntersectLine - intersect the spline with a line 

      bool IntersectLine(
              const CArbCoord2D &pt0,
              const CArbCoord2D &pt1,
              double            *u,
              double            *t) const

        pt0 - (in)  first point on line 
        pt1 - (in)  second point on line 
        u   - (out) spline intersection coordinate 
        t   - (out) line intersection coordinate 

      Description: This function attempts to fine the point of 
          intersection between the B spline and a straight line. The 
          line is defined by two points. If succesful, the routine 
          returns the parametric coordinates of the intersection for 
          both the line and the spline. For the line, the first and 
          second points are assumed to have a parametric coordinates 
          of 0 and 1, respectively. 

      Return Value: True if an intersection point is found, false 
          otherwise. 


    ClosestPoint - find the spline point closest to the input point 

      bool ClosestPoint(
              const CArbCoord2D &pt,
              double            *u) const

        pt - (in)  search point 
        u  - (out) parametric coordinate of the close point 

      Description: Given an input point, this routine attempts to 
          find the point on the spline that is closest. 

      Return Value: A flag is returned to indicate if the search for 
          a point was successful. In theory, the search should always 
          work, but it is a nonlinear search, an in some cases the 
          search may not converge. 


    UniformParameterization - reparametrize for a uniform mapping 

      void UniformParameterization(const int npts = 10000)

        npts - (in)  number of points to use in the reevaluation, 
                     should be 'large' 

      Description: This function reparameterizes a spline so that the 
          mapping from the spline's parametric space to the arc 
          length in the Cartesian space is nearly constant along the 
          spline. 


PRIVATE INTERFACE

  Private Member Functions:

    GetControlNet - determine the control points for a spline 

      void GetControlNet(
              int         ns,
              CArbCoord2D *points)

        ns     - (in)  number of input points 
        points - (in)  input points 

      Description: This function determines a set of control points 
          (control net) for a spline. These points are determined 
          such that the resulting spline will interpolate the input 
          points. This is done by building an solving a system of 
          linear equations. 


  Private Member Variables:

        int Num_Segs - number of spline segments 

        CArbCoord2D *Coord - list of control point coordinates 

        double *U - the knot vector 

*/
#endif
