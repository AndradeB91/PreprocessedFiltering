//
// Geom2DMixIn MixIn Class header file
//
// Description -
//   This is a MixIn class that supports common geometrical
//   computations
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2007
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//

#ifndef Geom2DMixIn_hh
#define Geom2DMixIn_hh

#include "Vec2D.hpp"

namespace FTools {

/// \file
///
/// \brief This file defines the Geom2DMixIn class

/** \brief A MixIn class for 2D geometrical operations
 *
 *  This class is a mixin class that defines many common 2D
 *  geometrical computations. 
 */


class Geom2DMixIn {
    public:

        Geom2DMixIn() {} ;

        /// return the angle ibj in the range -pi to pi

        double Angle(const Vec2D &b,const Vec2D &i,
                     const Vec2D &j) const ;

        /// return the angle ibj in the range 0 to 2pi

        double Angle2Pi(const Vec2D &b,const Vec2D &i,
                        const Vec2D &j) const ;

        /// return the angle between prev and next in the range -pi to pi

        double AngleVect(const Vec2D &prev,
                         const Vec2D &next) const ;

        /// return the angle between prev and next in the range 0 to 2pi

        double Angle2PiVect(const Vec2D &prev,
                            const Vec2D &next) const ;

        /// return the normalized vector that bisects angle ibj

        Vec2D BisectNorm(const Vec2D &b,
                               const Vec2D &i,
                               const Vec2D &j) const ;

        /// return the normalized vector that bisects angle between prev and next

        Vec2D BisectNormVect(const Vec2D &prev,
                                   const Vec2D &next) const ;

        /// return the perpendicular distance between b and line segment ij

        double PerpendicularDist(const Vec2D &b,
                                 const Vec2D &i,
                                 const Vec2D &j,
                                 double *blen) const ;

        /// return a "score" for the angle between ref and vec in the range -1 to 1

        double ScoreRelativeAngle(const Vec2D &ref,
                                  const Vec2D &vect) const ;

        /// true if a line from b to minus infinity crosses line segment ij

        bool ScanCross(const Vec2D &b,
                       const Vec2D &i,
                       const Vec2D &j) const ;

        /// true if pt falls on line segment line0-line1 within and internal tolearance

        bool IsPointOnLine(const Vec2D &pt,
                           const Vec2D &line0,
                           const Vec2D &line1,
                           double *line_param) const ;

        /// true if pt falls within the polygon described by pts

        bool IsPointInPoly(const Vec2D &pt,
                           int num_verts,
                           const Vec2D *poly_verts,
                           double tol) const ;
} ;

} //namespace

#endif
