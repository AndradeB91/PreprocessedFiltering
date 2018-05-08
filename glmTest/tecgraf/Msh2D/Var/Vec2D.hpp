//
// Vec2D object header file
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2007
//   All rights reserved
//
// ------------------------------------------------------------------
//
//  A vec2d object represents a 3 element vector, typically used
//  as a Cartesian position vector.
//
//  Public Interface:
//
//      Constructors:
//          Vec2D()
//          Vec2D(double x,double y, double z)
//          Vec2D(double *xyz)
//
//      Methods:
//          v.x()       
//          v.y()              
//          v.z()                  - returns the components of the vector
//          v.x(double x)          
//          v.y(double y)                     
//          v.z(double z)          - sets the components of the vector
//          v.Magnitude()          - returns the Cartesian length of the vector
//          v.Normalize()          - returns a normalized vector
//          v.max(const Vec2D& v2) - componentwise maximum
//          v.min(const Vec2D& v2) - comontentwise minimum
//
//      Operators:
//          v[i]                 - accesses the i'th component, 0 <= i <= 2
//          v0 += v1             - inplace vector add
//          v0 -= v1             - inplace vector subtraction
//          v0 *= float          - scalar multiplication
//          v0 /= float          - scalar division
//          -v0                  - negation
//          v0 == v1             - equal
//          v0 != v1             - not equal
//          v0 < v1              - true if all v0 components < corresponding v1 components
//          v0 > v1              - true if all v0 components > corresponding v1 components
//          v0 <= v1             - true if all v0 components <= corresponding v1 components
//          v0 >= v1             - true if all v0 components >= corresponding v1 components
//          v0 + v1              - vector addition
//          v0 - v1              - vector subtraction
//          v0 * float           - scalar multiplication
//          float * v0           - scalar multiplication
//          v0 * v1              - inner product
//          v0 / float           - scalar division
//          max(v0,v1)           - componentwise maximum
//          min(v0,v1)           - componentwise minimum
//          CrossProd(v0,v1)     - cross product
//          TripleProd(v0,v1,v2) - (v0 x v1) * v2
//          out << v0            - stream output
//          in >> v0             - stream input
//
 
#ifndef Vec2D_hpp
#define Vec2D_hpp

#include <cassert>
#include <cmath>
#include <iostream>

namespace FTools {

/// \file
///
/// \brief This file defines the Vec2D class and associated operators

/** \brief A two element vector
 *  \ingroup FundTypes
 *
 * This class represents a two element vector.  Typically used to represent
 * surface coordinates.  A number of associated operators support
 * normal vector algebra operation.
 *
 * \f[
 *  {\bf x} = \left\{ {\begin{array}{*{20}c}
 *     x  \\
 *     y  \\
 *  \end{array}} \right\}
 *  \f]
 */

class Vec2D {

    public:

        /// no argument constructor (x = y = 0)

        Vec2D() { data[0] = data[1] = 0.0 ; }

        /// construct from scalar arguments

        Vec2D(double x,double y) {
            data[0] = x ; data[1] = y ; }

        /// construct from an array argument 

        Vec2D(double *x)
        {
            if (x == 0) {
                data[0] = 0.0 ; data[1] = 0.0 ;
            } else {
                data[0] = x[0] ; data[1] = x[1] ;
            }
        }

    // access
     
        double x() const { return(data[0]) ; } ///< access to the x component
        double y() const { return(data[1]) ; } ///< access to the y component

        void x(double x) { data[0] = x ; } ///< sets the x component
        void y(double y) { data[1] = y ; } ///< sets the y component

    // operators

        /// retrieves a component by index (indx = 0 or 1)

        double operator [] (const int indx) const {
            assert((indx == 0) || (indx == 1)) ;
            return(data[indx]) ;
        }
        /// sets a component by index (indx = 0 or 1)

        double &operator[] (const int indx) {
            assert((indx == 0) || (indx == 1)) ;
            return(data[indx]) ;
        }

        /// in place vector addition, a += b => a[i] = a[i] + b[i], i = 0,1

        Vec2D operator += (const Vec2D &op2) {
            data[0] += op2.data[0] ;
            data[1] += op2.data[1] ;
            return(*this) ;
        }

        /// in place scalar addition, a += s => a[i] = a[i] + s, i = 0,1

        Vec2D operator += (const double &op2) {
            data[0] += op2 ;
            data[1] += op2 ;
            return(*this) ;
        }

        /// in place vector subtraction, a -= b => a[i] = a[i] - b[i], i = 0,1

        Vec2D operator -= (const Vec2D &op2) {
            data[0] -= op2.data[0] ;
            data[1] -= op2.data[1] ;
            return(*this) ;
        }

        /// in place scalar subtraction, a -= s => a[i] = a[i] - s, i = 0,1

        Vec2D operator -= (const double &op2) {
            data[0] -= op2 ;
            data[1] -= op2 ;
            return(*this) ;
        }
 
        /// in place multiplication by a scalar, a *= s => a[i] = a[i] * s, i = 0,1

        Vec2D operator *= (const double &op2) {
            data[0] *= op2 ;
            data[1] *= op2 ;
            return(*this) ;
        }

        /// in place division by a scalar, a /= s => a[i] = a[i] / s, i = 0,1

        Vec2D operator /= (const double &op2) {
            data[0] /= op2 ;
            data[1] /= op2 ;
            return(*this) ;
        }

        /// in place division by a vectir, a /= b => a[i] = a[i] / b[i], i = 0,1

        Vec2D operator /= (const Vec2D &op2) {
            data[0] /= op2.data[0] ;
            data[1] /= op2.data[1] ;
            return(*this) ;
        }

        /// unary negation operator (-a)

        Vec2D operator - () const {
            return(Vec2D(-data[0],-data[1])) ;
        }

    // member functions

        /// magnitude of the vector, \f$ \sqrt{x^2 + y^2} \f$

        double Magnitude() const {
            return(sqrt(data[0]*data[0]+
	                data[1]*data[1])) ;
        } ;

        /// returns a normalized version of the vector

        Vec2D Normalize() const {
            double len = sqrt(data[0]*data[0]+
	                      data[1]*data[1]) ;
            return(Vec2D(data[0]/len,data[1]/len)) ;
        } ;

        /// a.Max(b) => a[i] = max(a[i],b[i]), i = 0,1

        Vec2D Maxv(const Vec2D &op2) {
            data[0] = data[0] >= op2.x() ? data[0] : op2.x() ;
            data[1] = data[1] >= op2.y() ? data[1] : op2.y() ;
            return(*this) ;
        } ;

        /// a.Minv(b) => a[i] = min(a[i],b[i]), i = 0,1

        Vec2D Minv(const Vec2D &op2) {
            data[0] = data[0] <= op2.x() ? data[0] : op2.x() ;
            data[1] = data[1] <= op2.y() ? data[1] : op2.y() ;
            return(*this) ;
        } ;

    private:
        double data[2] ;

        friend double operator * (const Vec2D &op1,
                                  const Vec2D &op2) ;
} ;
                                                 
/// a == b => true iff a[i] == b[i], i = 0,1
/// \ingroup FundTypes

inline int operator == (const Vec2D &c1,const Vec2D &c2)
{
    return((c1[0] == c2[0]) && (c1[1] == c2[1])) ;
}

/// a != b => true if any a[i] != b[i], i = 0,1
/// \ingroup FundTypes

inline int operator != (const Vec2D &c1,const Vec2D &c2)
{
    return((c1[0] != c2[0]) || (c1[1] != c2[1])) ;
}

/// a < b => true if any a[i] <= b[i], i = 0,1
/// \ingroup FundTypes

inline int operator < (const Vec2D &c1,const Vec2D &c2)
{
    return((c1[0] < c2[0]) && (c1[1] < c2[1])) ;
}

/// a > b => true if any a[i] >= b[i], i = 0,1
/// \ingroup FundTypes

inline int operator > (const Vec2D &c1,const Vec2D &c2)
{
    return((c1[0] > c2[0]) && (c1[1] > c2[1])) ;
}

/// a <= b => true if any a[i] <= b[i], i = 0,1
/// \ingroup FundTypes

inline int operator <= (const Vec2D &c1,const Vec2D &c2)
{
    return((c1[0] <= c2[0]) && (c1[1] <= c2[1])) ;
}

/// a >= b => true if any a[i] >= b[i], i = 0,1
/// \ingroup FundTypes

inline int operator >= (const Vec2D &c1,const Vec2D &c2)
{
    return((c1[0] >= c2[0]) && (c1[1] >= c2[1])) ;
}

/// c = a + b where c[i] = a[i] + b[i], i = 0,1
/// \ingroup FundTypes

inline Vec2D operator + (const Vec2D &op1,const Vec2D &op2)
{
    Vec2D op = op1 ;
    op += op2 ;
    return(op) ;
}

/// c = a - b where c[i] = a[i] - b[i], i = 0,1
/// \ingroup FundTypes

inline Vec2D operator - (const Vec2D &op1,const Vec2D &op2)
{
    Vec2D op = op1 ;
    op -= op2 ;
    return(op) ;
}

/// c = a + s where c[i] = a[i] + s, i = 0,1
/// \ingroup FundTypes

inline Vec2D operator + (const Vec2D &op1,double op2)
{
    Vec2D op = op1 ;
    op += op2 ;
    return(op) ;
}

/// c = a - s where c[i] = a[i] - s, i = 0,1
/// \ingroup FundTypes

inline Vec2D operator - (const Vec2D &op1,double &op2)
{
    Vec2D op = op1 ;
    op -= op2 ;
    return(op) ;
}

/// c = a * s where c[i] = a[i] * s, i = 0,1
/// \ingroup FundTypes

inline Vec2D operator * (const Vec2D &op1,const double &op2)
{
    Vec2D op = op1 ;
    op *= op2 ;
    return(op) ;
}

/// c = s * a where c[i] = s * a[i], i = 0,1
/// \ingroup FundTypes

inline Vec2D operator * (const double &op1,const Vec2D &op2)
{
    Vec2D op = op2 ;
    op *= op1 ;
    return(op) ;
}

/// inner product, s = a * b where s = sum(a[i] * b[i]), i = 0,1
/// \ingroup FundTypes

inline double operator * (const Vec2D &op1,const Vec2D &op2)
{
    return(op1.x()*op2.x() + op1.y()*op2.y()) ;
}

/// c = a / s where c[i] = a[i] / s, i = 0,1
/// \ingroup FundTypes

inline Vec2D operator / (const Vec2D &op1,const double &op2)
{
    Vec2D op = op1 ;
    op /= op2 ;
    return(op) ;
}

/// c = a / b where c[i] = a[i] / b[i], i = 0,1
/// \ingroup FundTypes

inline Vec2D operator / (const Vec2D &op1,const Vec2D &op2)
{
    Vec2D op = op1 ;
    op /= op2 ;
    return(op) ;
}

/// c = Maxv(a,b) => c[i] = max(a[i],b[i]), i = 0,1
/// \ingroup FundTypes

inline Vec2D Maxv(const Vec2D &op1,const Vec2D &op2)
{
    return(Vec2D((op1.x() >= op2.x() ? op1.x() : op2.x()),
                 (op1.y() >= op2.y() ? op1.y() : op2.y()))) ;
}

/// c = Minv(a,b) => c[i] = min(a[i],b[i]), i = 0,1
/// \ingroup FundTypes

inline Vec2D Minv(const Vec2D &op1,const Vec2D &op2)
{
    return(Vec2D((op1.x() <= op2.x() ? op1.x() : op2.x()),
                 (op1.y() <= op2.y() ? op1.y() : op2.y()))) ;
}

/// vector cross product, c = a cross b
/// \ingroup FundTypes

inline double CrossProd(const Vec2D &op1,const Vec2D &op2)
{
    return((op1.x()*op2.y()) - (op1.y()*op2.x())) ; 
}

/// stream output operator
/// \ingroup FundTypes

inline std::ostream &operator << (std::ostream &out,const Vec2D &data)
{
    out <<data.x() << ' ' << data.y() ;
    return(out) ;
}

/// stream input operator
/// \ingroup FundTypes

inline std::istream &operator >> (std::istream &in, Vec2D &data)
{
    in >> data[0] >> data[1] ;
    return(in) ;
}

} // namespace

#endif
