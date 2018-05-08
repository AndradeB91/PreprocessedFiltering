//
// Vec3D object header file
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2007
//   All rights reserved
//
// ------------------------------------------------------------------
//
//  A vec3d object represents a 3 element vector, typically used
//  as a Cartesian position vector.
//
//  Public Interface:
//
//      Constructors:
//          Vec3D()
//          Vec3D(double x,double y, double z)
//          Vec3D(double *xyz)
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
//          v.MaxV(const Vec3D& v2) - componentwise maximum
//          v.MinV(const Vec3D& v2) - comontentwise minimum
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
//          MaxV(v0,v1)          - componentwise maximum
//          MinV(v0,v1)          - componentwise minimum
//          CrossProd(v0,v1)     - cross product
//          TripleProd(v0,v1,v2) - (v0 x v1) * v2
//          out << v0            - stream output
//          in >> v0             - stream input
//
 
#ifndef Vec3D_h
#define Vec3D_h

#include <cassert>
#include <cmath>
#include <iostream>

namespace FTools {

/// \file
///
/// \brief This file defines the Vec3D class and associated operators

/** \brief A three element vector
 *  \ingroup FundTypes
 *
 * This class represents a three element vector.  Typically used to represent
 * a point in Cartesian space.  A number of associated operators support
 * normal vector algebra operation.
 *
 * \f[
 *  {\bf x} = \left\{ {\begin{array}{*{20}c}
 *     x  \\
 *     y  \\
 *     z  \\
 *  \end{array}} \right\}
 *  \f]
 */

class Vec3D {

    public:

    // constructors and destructors

        /// no argument constructor (x = y = z = 0)

        Vec3D() { data[0] = data[1] = data[2] = 0.0 ; }

        /// construct from scalar arguments

        Vec3D(double x,double y,double z) {
            data[0] = x ; data[1] = y ; data[2] = z ; }

        /// construct from an array argument 

        Vec3D(double *x)
        {
            if (x == 0) {
                data[0] = 0.0 ; data[1] = 0.0 ; data[2] = 0.0 ;
            } else {
                data[0] = x[0] ; data[1] = x[1] ; data[2] = x[2] ;
            }
        }

    // access
     
        double x() const { return(data[0]) ; } ///< access to the x component
        double y() const { return(data[1]) ; } ///< access to the y component
        double z() const { return(data[2]) ; } ///< access to the z component

        void x(double x) { data[0] = x ; } ///< sets the x component
        void y(double y) { data[1] = y ; } ///< sets the y component
        void z(double z) { data[2] = z ; } ///< sets the z component

    // operators

        /// retrieves a component by index (0 <= indx <= 2)

        double operator [] (const int indx) const {
            assert((indx >= 0) && (indx <= 2)) ;
            return(data[indx]) ;
        } ;

        /// sets a component by index (0 <= indx <= 2)

        double &operator[] (const int indx) {
            assert((indx >= 0) && (indx <= 2)) ;
            return(data[indx]) ;
        } ;

        /// in place vector addition, a += b => a[i] = a[i] + b[i], i = 0..2

        Vec3D operator += (const Vec3D &op2) {
            data[0] += op2.data[0] ;
            data[1] += op2.data[1] ;
	    data[2] += op2.data[2] ;
            return(*this) ;
        } ;

        /// in place scalar addition, a += s => a[i] = a[i] + s, i = 0..2

        Vec3D operator += (const double &op2) {
            data[0] += op2 ;
            data[1] += op2 ;
            data[2] += op2 ;
            return(*this) ;
        } ; 

        /// in place vector subtraction, a -= b => a[i] = a[i] - b[i], i = 0..2

        Vec3D operator -= (const Vec3D &op2) {
            data[0] -= op2.data[0] ;
            data[1] -= op2.data[1] ;
	    data[2] -= op2.data[2] ;
            return(*this) ;
        } ;

        /// in place scalar subtraction, a -= s => a[i] = a[i] - s, i = 0..2

        Vec3D operator -= (const double &op2) {
            data[0] -= op2 ;
            data[1] -= op2 ;
            data[2] -= op2 ;
            return(*this) ;
        } ;

        /// in place multiplication by a scalar, a *= s => a[i] = a[i] * s, i = 0..2

        Vec3D operator *= (const double &op2) {
            data[0] *= op2 ;
            data[1] *= op2 ;
	    data[2] *= op2 ;
            return(*this) ;
        } ;

        /// in place division by a scalar, a /= s => a[i] = a[i] / s, i = 0..2

        Vec3D operator /= (const double &op2) {
            data[0] /= op2 ;
            data[1] /= op2 ;
	    data[2] /= op2 ;
            return(*this) ;
        } ;

        /// unary negation operator (-a)

        Vec3D operator - () const {
            return(Vec3D(-data[0],-data[1],-data[2])) ;
        } ;

    // member functions

        /// magnitude of the vector, \f$ \sqrt{x^2 + y^2 + z^2} \f$

        double Magnitude() const {
            return(sqrt(data[0]*data[0]+
	                data[1]*data[1]+
			data[2]*data[2])) ;
        } ;

        /// returns a normalized version of the vector

        Vec3D Normalize() const {
            double len = sqrt(data[0]*data[0]+
	                      data[1]*data[1]+
			      data[2]*data[2]) ;
            return(Vec3D(data[0]/len,data[1]/len,data[2]/len)) ;
        } ;

        /// a.Max(b) => a[i] = max(a[i],b[i]), i = 0..2

        Vec3D Maxv(const Vec3D &op2) {
            data[0] = data[0] >= op2.x() ? data[0] : op2.x() ;
            data[1] = data[1] >= op2.y() ? data[1] : op2.y() ;
            data[2] = data[2] >= op2.z() ? data[2] : op2.z() ;
            return(*this) ;
        } ;

        /// a.Minv(b) => a[i] = min(a[i],b[i]), i = 0..2

        Vec3D Minv(const Vec3D &op2) {
            data[0] = data[0] <= op2.x() ? data[0] : op2.x() ;
            data[1] = data[1] <= op2.y() ? data[1] : op2.y() ;
            data[2] = data[2] <= op2.z() ? data[2] : op2.z() ;
            return(*this) ;
        } ;

    private:
        double data[3] ;

        friend double operator * (const Vec3D &op1,
                                  const Vec3D &op2) ;
} ;

/// a == b => true iff a[i] == b[i], i = 0..2
/// \ingroup FundTypes

inline int operator == (const Vec3D &c1,const Vec3D &c2)
{
    return((c1[0] == c2[0]) && (c1[1] == c2[1]) && (c1[2] == c2[2])) ;
}

/// a != b => true if any a[i] != b[i], i = 0..2
/// \ingroup FundTypes

inline int operator != (const Vec3D &c1,const Vec3D &c2)
{
    return((c1[0] != c2[0]) || (c1[1] != c2[1]) || (c1[2] != c2[2])) ;
}

/// a < b => true if any a[i] <= b[i], i = 0..2
/// \ingroup FundTypes

inline int operator < (const Vec3D &c1,const Vec3D &c2)
{
    return((c1[0] <= c2[0]) && (c1[1] <= c2[1]) && (c1[2] <= c2[2])) ;
}

/// a > b => true if any a[i] >= b[i], i = 0..2
/// \ingroup FundTypes

inline int operator > (const Vec3D &c1,const Vec3D &c2)
{
    return((c1[0] >= c2[0]) && (c1[1] >= c2[1]) && (c1[2] >= c2[2])) ;
}

/// a <= b => true if any a[i] <= b[i], i = 0..2
/// \ingroup FundTypes

inline int operator <= (const Vec3D &c1,const Vec3D &c2)
{
    return((c1[0] <= c2[0]) && (c1[1] <= c2[1]) && (c1[2] <= c2[2])) ;
}

/// a >= b => true if any a[i] >= b[i], i = 0..2
/// \ingroup FundTypes

inline int operator >= (const Vec3D &c1,const Vec3D &c2)
{
    return((c1[0] >= c2[0]) && (c1[1] >= c2[1]) && (c1[2] >= c2[2])) ;
}

/// c = a + b where c[i] = a[i] + b[i], i = 0..2
/// \ingroup FundTypes

inline Vec3D operator + (const Vec3D &op1,const Vec3D &op2)
{
    Vec3D op = op1 ;
    op += op2 ;
    return(op) ;
}

/// c = a + s where c[i] = a[i] + s, i = 0..2
/// \ingroup FundTypes

inline Vec3D operator + (const Vec3D &op1,double op2)
{
    Vec3D op = op1 ;
    op += op2 ;
    return(op) ;
}

/// c = a - b where c[i] = a[i] - b[i], i = 0..2
/// \ingroup FundTypes

inline Vec3D operator - (const Vec3D &op1,const Vec3D &op2)
{
    Vec3D op = op1 ;
    op -= op2 ;
    return(op) ;
}

/// c = a - s where c[i] = a[i] - s, i = 0..2
/// \ingroup FundTypes

inline Vec3D operator - (const Vec3D &op1,double op2)
{
    Vec3D op = op1 ;
    op -= op2 ;
    return(op) ;
}

/// c = a * s where c[i] = a[i] * s, i = 0..2
/// \ingroup FundTypes

inline Vec3D operator * (const Vec3D &op1,const double &op2)
{
    Vec3D op = op1 ;
    op *= op2 ;
    return(op) ;
}

/// c = s * a where c[i] = s * a[i], i = 0..2
/// \ingroup FundTypes

inline Vec3D operator * (const double &op1,const Vec3D &op2)
{
    Vec3D op = op2 ;
    op *= op1 ;
    return(op) ;
}

/// inner product, s = a * b where s = sum(a[i] * b[i]), i = 0..2
/// \ingroup FundTypes

inline double operator * (const Vec3D &op1,const Vec3D &op2)
{
    return(op1.x()*op2.x() + op1.y()*op2.y() + op1.z()*op2.z()) ;
}

/// c = a / s where c[i] = a[i] / s, i = 0..2
/// \ingroup FundTypes

inline Vec3D operator / (const Vec3D &op1,const double &op2)
{
    Vec3D op = op1 ;
    op /= op2 ;
    return(op) ;
}

/// c = Maxv(a,b) => c[i] = max(a[i],b[i]), i = 0..2
/// \ingroup FundTypes

inline Vec3D Maxv(const Vec3D &op1,const Vec3D &op2)
{
    return(Vec3D((op1.x() >= op2.x() ? op1.x() : op2.x()),
                 (op1.y() >= op2.y() ? op1.y() : op2.y()),
                 (op1.z() >= op2.z() ? op1.z() : op2.z()))) ;
}

/// c = Minv(a,b) => c[i] = min(a[i],b[i]), i = 0..2
/// \ingroup FundTypes

inline Vec3D Minv(const Vec3D &op1,const Vec3D &op2)
{
    return(Vec3D((op1.x() <= op2.x() ? op1.x() : op2.x()),
                 (op1.y() <= op2.y() ? op1.y() : op2.y()),
                 (op1.z() <= op2.z() ? op1.z() : op2.z()))) ;
}

/// vector cross product, c = a cross b
/// \ingroup FundTypes

inline Vec3D CrossProd(const Vec3D &op1,const Vec3D &op2)
{
    return(Vec3D((op1.y()*op2.z()) - (op1.z()*op2.y()), 
                 (op1.z()*op2.x()) - (op1.x()*op2.z()), 
                 (op1.x()*op2.y()) - (op1.y()*op2.x()))) ; 
}

/// vector triple product, s = (a cross b) * c
/// \ingroup FundTypes

inline double TripleProd(const Vec3D &op1,const Vec3D &op2,
                         const Vec3D &op3)
{
    Vec3D tmp = CrossProd(op1,op2) ;
    return(tmp * op3) ; 
}

/// stream output operator
/// \ingroup FundTypes

inline std::ostream &operator << (std::ostream &out,const Vec3D &data)
{
    out <<data.x() << ' ' << data.y() << ' ' << data.z() ;
    return(out) ;
}

/// stream input operator
/// \ingroup FundTypes

inline std::istream &operator >> (std::istream &in, Vec3D &data)
{
    in >> data[0] >> data[1] >> data[2] ;
    return(in) ;
}

} // namespace

#endif
