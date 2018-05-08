//
// CArbGeom2DMixIn MixIn Class header file
//
// Description -
//   This is a MixIn class that supports common geometrical
//   computations
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2001
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.4 $  $Date: 2002/07/25 19:23:49 $  $Author: wash $
//

#ifndef ArbGeom2DMixIn_hh
#define ArbGeom2DMixIn_hh

#include "ArbCoord2D.hpp"

class CArbGeom2DMixIn {
    public:

        CArbGeom2DMixIn() {} ;

        double Angle(const CArbCoord2D &b,const CArbCoord2D &i,
                     const CArbCoord2D &j) const ;

        double Angle2Pi(const CArbCoord2D &b,const CArbCoord2D &i,
                        const CArbCoord2D &j) const ;

        double AngleVect(const CArbCoord2D &prev,
                         const CArbCoord2D &next) const ;

        double Angle2PiVect(const CArbCoord2D &prev,
                            const CArbCoord2D &next) const ;

        CArbCoord2D BisectNorm(const CArbCoord2D &b,
                               const CArbCoord2D &i,
                               const CArbCoord2D &j) const ;

        CArbCoord2D BisectNormVect(const CArbCoord2D &prev,
                                   const CArbCoord2D &next) const ;

        double PerpendicularDist(const CArbCoord2D &b,
                                 const CArbCoord2D &i,
                                 const CArbCoord2D &j,
                                 double *blen) const ;

        double ScoreRelativeAngle(const CArbCoord2D &ref,
                                  const CArbCoord2D &vect) const ;

        bool ScanCross(const CArbCoord2D &b,
                       const CArbCoord2D &i,
                       const CArbCoord2D &j) const ;

        bool IsPointOnLine(const CArbCoord2D &pt,
                           const CArbCoord2D &line0,
                           const CArbCoord2D &line1,
                           double *line_param) const ;

        bool IsPointInPoly(const CArbCoord2D &pt,
                           int num_verts,
                           const CArbCoord2D *poly_verts,
                           double tol) const ;

    private:
} ;

#endif
