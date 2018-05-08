//
// ArbCoord2D header file
//
// Description -
//   This is the header file for the ArbCoord2D objects.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2000
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.8 $  $Date: 2002/09/06 14:56:25 $  $Author: wash $
//
 
#ifndef ArbCoord2D_h
#define ArbCoord2D_h

#include <assert.h>
#include <math.h>

class CArbCoord2D {

    public:

    // constructors and destructors

        CArbCoord2D() { coord[0] = coord[1] = 0.0 ; } ;
        CArbCoord2D(double x,double y) {
            coord[0] = x ; coord[1] = y ; } ;
        CArbCoord2D(double *x) {
            coord[0] = x[0] ; coord[1] = x[1] ; } ;

    // access
     
        double x() const { return(coord[0]) ; } ;
        double y() const { return(coord[1]) ; } ;

    // operators

        double operator [] (const int indx) const {
            assert((indx == 0) || (indx == 1)) ;
            return(coord[indx]) ;
        } ;

        double &operator[] (const int indx) {
            assert((indx == 0) || (indx == 1)) ;
            return(coord[indx]) ;
        } ;

        CArbCoord2D operator += (const CArbCoord2D &op2) {
            coord[0] += op2.coord[0] ;
            coord[1] += op2.coord[1] ;
            return(*this) ;
        } ;

        CArbCoord2D operator -= (const CArbCoord2D &op2) {
            coord[0] -= op2.coord[0] ;
            coord[1] -= op2.coord[1] ;
            return(*this) ;
        } ;

        CArbCoord2D operator *= (const double &op2) {
            coord[0] *= op2 ;
            coord[1] *= op2 ;
            return(*this) ;
        } ;

        CArbCoord2D operator /= (const double &op2) {
            coord[0] /= op2 ;
            coord[1] /= op2 ;
            return(*this) ;
        } ;

        CArbCoord2D operator /= (const CArbCoord2D &op2) {
            coord[0] /= op2.coord[0] ;
            coord[1] /= op2.coord[1] ;
            return(*this) ;
        } ;

        CArbCoord2D operator - () const {
            return(CArbCoord2D(-coord[0],-coord[1])) ;
        } ;

    // member functions

        double Magnitude() const {
            return(sqrt(coord[0]*coord[0]+coord[1]*coord[1])) ;
        } ;

        CArbCoord2D Normalize() const {
            double len = sqrt(coord[0]*coord[0]+coord[1]*coord[1]) ;
            return(CArbCoord2D(coord[0]/len,coord[1]/len)) ;
        } ;

    private:
        double coord[2] ;

        friend double operator * (const CArbCoord2D &op1,
                                  const CArbCoord2D &op2) ;
} ;
                                                 
inline int operator == (const CArbCoord2D &c1,const CArbCoord2D &c2)
{
    return((c1[0] == c2[0]) && (c1[1] == c2[1])) ;
}

inline int operator != (const CArbCoord2D &c1,const CArbCoord2D &c2)
{
    return((c1[0] != c2[0]) || (c1[1] != c2[1])) ;
}

inline CArbCoord2D operator + (const CArbCoord2D &op1,
                               const CArbCoord2D &op2)
{
    CArbCoord2D op = op1 ;
    op += op2 ;
    return(op) ;
}

inline CArbCoord2D operator - (const CArbCoord2D &op1,
                               const CArbCoord2D &op2)
{
    CArbCoord2D op = op1 ;
    op -= op2 ;
    return(op) ;
}

inline CArbCoord2D operator * (const CArbCoord2D &op1,
                               const double &op2)
{
    CArbCoord2D op = op1 ;
    op *= op2 ;
    return(op) ;
}

inline CArbCoord2D operator * (const double &op1,
                               const CArbCoord2D &op2)
{
    CArbCoord2D op = op2 ;
    op *= op1 ;
    return(op) ;
}

inline double operator * (const CArbCoord2D &op1,
                          const CArbCoord2D &op2)
{
    return(op1.coord[0]*op2.coord[0] + op1.coord[1]*op2.coord[1]) ;
}

inline CArbCoord2D operator / (const CArbCoord2D &op1,
                               const double &op2)
{
    CArbCoord2D op = op1 ;
    op /= op2 ;
    return(op) ;
}

inline CArbCoord2D operator / (const CArbCoord2D &op1,
                               const CArbCoord2D &op2)
{
    CArbCoord2D op = op1 ;
    op /= op2 ;
    return(op) ;
}

inline CArbCoord2D max(const CArbCoord2D &op1,const CArbCoord2D &op2)
{
    return(CArbCoord2D((op1.x() >= op2.x() ? op1.x() : op2.x()),
                       (op1.y() >= op2.y() ? op1.y() : op2.y()))) ;
}

inline CArbCoord2D min(const CArbCoord2D &op1,const CArbCoord2D &op2)
{
    return(CArbCoord2D((op1.x() <= op2.x() ? op1.x() : op2.x()),
                       (op1.y() <= op2.y() ? op1.y() : op2.y()))) ;
}

inline double CrossProduct(const CArbCoord2D &op1,
                           const CArbCoord2D &op2)
{
    return(op1.x()*op2.y() - op2.x()*op1.y()) ;
}

    


/*
CLASS CArbCoord2D

  This class implements a 2D coordinate object that can be used to 
  represent a point location or a vector. It also implements a number 
  of operators for perfoming basic math operations on the points. 


PUBLIC INTERFACE

  Public Member Functions:

    CArbCoord2D - no argument constructor 

      CArbCoord2D()

      Description: This is the no argument constructor for the 
          object. The coordinate value is set to 0,0. 


    CArbCoord2D - class constructor 

      CArbCoord2D(
              double x,
              double y)

        x - (in)  initial x coordinate 
        y - (in)  initial y coordinate 

      Description: This is a constructor function that takes initial 
          values for the coordinates as arguments. 


    CArbCoord2D - class constructor 

      CArbCoord2D(double *x)

        x - (in)  initial value array 

      Description: This is a constructor function that allows the 
          initial coordinate values to be specified as a two element 
          array. 


    x - x coordinate access 

      double x() const

      Description: This function returns the x coordinate value. 

      Return Value: The x coordinate value. 


    y - y coordinate access 

      double y() const

      Description: This function returns the y coordinate value. 

      Return Value: The y coordinate value. 


    operator [] - coordinate access 

      double operator [] (const int indx) const

        indx - (in)  0 for x, 1 for y 

      Description: This operator allows the x and y coordinates to be 
          accessed as if the coordinate object is a two element 
          array. 

      Return Value: The value of the x or y coordinate depending on 
          if the argument is 0 or 1, respectively. 


    operator [] - coordinate references 

      double &operator [] (const int indx)

        indx - (in)  0 for x, 1 for y 

      Description: This operator returns a reference to the x and y 
          coordinates so that they can be assigned values as if the 
          coordinate object is a two element array. 

      Return Value: A reference to the x or y coordinate depending on 
          if the argument is 0 or 1, respectively. 


    operator += - addition-assignment operator 

      CArbCoord2D operator += (const CArbCoord2D &op2)

        op2 - (in)  operand 

      Description: This operator adds the value of the operand to the 
          target and returns the updated object. 

      Return Value: (x0,y0) += (x1,y1) returns (x0+x1,y0+y1) 


    operator -= - subtraction-assignment operator 

      CArbCoord2D operator -= (const CArbCoord2D &op2)

        op2 - (in)  operand 

      Description: This operator subtracts the value of the operand 
          from the target and returns the updated object. 

      Return Value: (x0,y0) -= (x1,y1) returns (x0-x1,y0-y1) 


    operator *= - multiplication-assignment operator 

      CArbCoord2D operator *= (const double &op2)

        op2 - (in)  operand 

      Description: This operator multiplies the components of the 
          target by the value of the operand and returns the updated 
          object. 

      Return Value: (x0,y0) *= op2 returns (x0*op2,y0*op2) 


    operator /= - division-assignment operator 

      CArbCoord2D operator /= (const double &op2)

        op2 - (in)  operand 

      Description: This operator divides the components of the target 
          by the value of the operand and returns the updated object. 

      Return Value: (x0,y0) /= op2 returns (x0/op2,y0/op2) 


    operator - - unary minus operator 

      CArbCoord2D operator - ()

      Description: This is the unary minus or negation operator. 

      Return Value: -(x,y) returns (-x,-y) 


    Magnitude - vector magnitude 

      double Magnitude()

      Description: This function computes the vector magnitude (L2 
          norm) of the object. 

      Return Value: sqrt( x*x + y*y ) 


    Normalize - vector normalize 

      CArbCoord2D Normalize()

      Description: This function returns a new object that contains a 
          vector with the same direction as the taget, but normalized 
          so that it has a unit length. 

      Return Value: ( x/sqrt(x*x+y*y), y/sqrt(x*x+y*y) ) 


PRIVATE INTERFACE

  Private Member Variables:

        double coord[2] - x & y coordinates 


NON-MEMBER OPERATORS

operator == - boolian comparison 

  inline int operator == (
          const CArbCoord2D &c1,
          const CArbCoord2D &c2)

    c1 - (in)  first operand 
    c2 - (in)  second operand 

  Description: Performs an element-wise comparison of two coordinates 
      and returns true if the x and y components are equal. 

  Return Value: (x0 == x1) && (y0 == y1) 


operator + - addition operator 

  inline CArbCoord2D operator + (
          const CArbCoord2D &op1,
          const CArbCoord2D &op2)

    op1 - (in)  first operand 
    op2 - (in)  second operand 

  Description: This is an addition operator that does a 
      component-wise addition of the operands. 

  Return Value: (x0,y0) + (x1,y1) returns (x0+x1,y0+y1) 


operator - - subtraction operator 

  inline CArbCoord2D operator - (
          const CArbCoord2D &op1,
          const CArbCoord2D &op2)

    op1 - (in)  first operand 
    op2 - (in)  second operand 

  Description: This is a subtraction operator that does a 
      component-wise subtraction of the operands. 

  Return Value: (x0,y0) - (x1,y1) returns (x0-x1,y0-y1) 


operator * - scalar/vector post-multiplication operator 

  inline CArbCoord2D operator * (
          const CArbCoord2D &op1,
          const double      &op2)

    op1 - (in)  first operand 
    op2 - (in)  second operand 

  Description: This is a scalar post-multiplication operator that 
      does a component-wise multiplication of the operands. 

  Return Value: (x0,y0)*op returns (x0*op,y0*op) 


operator * - scalar/vector pre-multiplication operator 

  inline CArbCoord2D operator * (
          const double      &op1,
          const CArbCoord2D &op2)

    op1 - (in)  first operand 
    op2 - (in)  second operand 

  Description: This is a scalar pre-multiplication operator that does 
      a component-wise multiplication of the operands. 

  Return Value: op*(x0,y0) returns (x0*op,y0*op) 


operator * - scalar multiplication operator (dot product) 

  inline double operator * (
          const CArbCoord2D &op1,
          const CArbCoord2D &op2)

    op1 - (in)  first operand 
    op2 - (in)  second operand 

  Description: This is a scalar multiplication operator that returs 
      the dot- or inner- product of the operands. 

  Return Value: (x0,y0) * (x1,y1) returns x0*x1 + y0*y1 


operator / - scalar vector division operator. 

  inline CArbCoord2D operator / (
          const CArbCoord2D &op1,
          const double      &op2)

    op1 - (in)  first operand 
    op2 - (in)  second operand 

  Description: This is a scalar division operator that does a 
      component-wise division of the operands. 

  Return Value: (x0,y0)/op returns (x0/op,y0/op) 

*/

#endif
