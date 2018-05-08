//
// Geom2DMixIn Class definition
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

#include "Geom2DMixIn.hpp"

#ifdef MEMDEBUG
#include "MemDbg.hpp"
//#define new new(__FILE__,__LINE__)
#endif

#define PI      3.14159265359
#define HALF_PI 1.57079632680
#define TWO_PI  6.28318530718


namespace FTools {

/// \file
///
/// \brief This file implements the Geom2DMixIn class


double Geom2DMixIn::Angle(const Vec2D &b,
                              const Vec2D &i,
                              const Vec2D &j) const
{
    return(AngleVect(i-b,j-b)) ;
}


double Geom2DMixIn::AngleVect(const Vec2D &prev,
                                  const Vec2D &next) const
{
    Vec2D bi = prev ;
    Vec2D bj = next ;

    double tmp = ((bi.x() * bj.x()) + (bi.y() * bj.y())) /
                 (sqrt(bi.x()*bi.x() + bi.y()*bi.y()) *
                  sqrt(bj.x()*bj.x() + bj.y()*bj.y())) ;
    if (tmp >  1.0) tmp =  1.0 ;
    if (tmp < -1.0) tmp = -1.0 ;
    double cross = bi.x()*bj.y() - bi.y()*bj.x() ;

    if (cross == 0.0) {
        return((tmp > 0) ? 0.0 : acos(-1.0)) ;
    }

    return(acos(tmp) * (cross/fabs(cross))) ;
}


double Geom2DMixIn::Angle2Pi(const Vec2D &b,
                                 const Vec2D &i,
                                 const Vec2D &j) const
{
    double angle = Angle(b,i,j) ;
    if (angle < 0.0) return(TWO_PI + angle) ;
    return(angle) ;
}


double Geom2DMixIn::Angle2PiVect(const Vec2D &prev,
                                     const Vec2D &next) const
{
    double angle = AngleVect(prev,next) ;
    if (angle < 0.0) return(TWO_PI + angle) ;
    return(angle) ;
}


/** This is the assumed node configuration:
 *
 *       * i
 *        \
 *         \
 *          \
 *           \
 *            *----------*
 *           b          j
 *
 *  We want to find the normalized vector that is the bisector of
 *  angle ibj.  The procedure is as follows:
 *
 *  1. normalize the vectors bi and bj.
 *  2. find the "right hand" normal to bj (that is, bj cross bj_normal
 *     is positive).
 *  3. find the "left hand" normal to bi
 *  4. Average the normals and renormalize.
 *
 *  There is only one degenerate case, were the include angle is 0.
 *  for this case return the normal at 180.
 */ 

Vec2D Geom2DMixIn::BisectNorm(
                         const Vec2D &b,
                         const Vec2D &i,
                         const Vec2D &j) const
{
    Vec2D n, d, bi, bj, nbi, nbj, nmean ;
    double len ;

    d = i - b ;
    len = d.Magnitude() ;
    bi = d / len ;
    nbi[0] = bi[1] ;  nbi[1] = -bi[0] ;

    d = j - b ;
    len = d.Magnitude() ;
    bj = d / len ;
    nbj[0] = -bj[1] ;  nbj[1] = bj[0] ;

    nmean = 0.5 * (nbi + nbj) ;
    len = nmean.x()*nmean.x() + nmean.y()*nmean.y() ;

    if (len < 0.0001) {   // degenerate case
        n = (bi + bj) / 2.0 ;
        n = n.Normalize() ;
        double angle = Angle(b,j,i) ;
        if (angle < 0.0) n = -n ;
    } else {
        len = sqrt(len) ;
        n = nmean / len ;
    }
    return(n) ;
}

/** This is the assumed node configuration:
 *
 *       ^
 *        \
 *         \next
 *          \
 *           \
 *            *---------->
 *                prev
 *
 *  We want to find the normalized vector that is the bisector of
 *  angle between prev and next.  We do this as follows
 *
 *  1. normalize the vectors.
 *  2. find the "right hand" normal to prev
 *        (that is, prev cross prev_normal)
 *  3. find the "left hand" normal to next
 *  4. Average the normals and renormalize.
 *
 *  There is only one degenerate case, were the include angle is 0.
 *  for this case return the normal at 180.
 */ 

Vec2D Geom2DMixIn::BisectNormVect(
                         const Vec2D &prev,
                         const Vec2D &next) const
{
    Vec2D tmp ;

    tmp = next.Normalize() ;
    Vec2D next_norm(tmp.y(),-tmp.x()) ;

    tmp = prev.Normalize() ;
    Vec2D prev_norm(-tmp.y(),tmp.x()) ;

    Vec2D nmean = 0.5 * (next_norm + prev_norm) ;

    if (nmean.Magnitude() < 0.0001) {   // degenerate case
        nmean = 0.5 * (prev + next) ;
        return(-nmean.Normalize()) ;
    }
    
    return(nmean.Normalize()) ;
}

/** This routine finds the perpendicular minimum distance
 *  between a point (b) and a line segment (i,j).  It also
 *  returns the length of the triangle base (i,j) in *blen
 *
 *  This is the picture:
 *
 *                   b
 *                 *        
 *                /|\       
 *               / | \  
 *              /  |  \     
 *             *---+---*    
 *            i    k    j   
 * 
 * the variables are defined as follows:
 *
 * blength = |ij|
 * area    = area of triangle ijb
 * pdist   = |bk|
 * adist1  = |bi|
 * adist2  = |bj|
 */

double Geom2DMixIn::PerpendicularDist(const Vec2D &b,
                     const Vec2D &i,const Vec2D &j,
                     double *blen) const
{
    double twice_area, pdist ;

    Vec2D base = j - i ;
    *blen = base.Magnitude() ;

    // check for b between i and j

    Vec2D ib = b - i ;
    Vec2D jb = b - j ;

    if ((*blen != 0.0) && (ib*base * jb*base <= 0.0)) {

        twice_area = fabs((i[0]-b[0])*(j[1]-b[1]) - 
                          (j[0]-b[0])*(i[1]-b[1])) ;
        pdist = twice_area / *blen ;
        return(pdist) ;

    }

    // otherwise find return the shorter radial distance

    double ibmag = ib.Magnitude() ;
    double jbmag = jb.Magnitude() ;
    return(ibmag <= jbmag ? ibmag : jbmag) ;
}

/** given two (assumed to be) normalized vectors, this function
 *  gives a "score" to the angle between them.  If the angle from
 *  the reference to the vector is differentially small in the
 *  counter clockwise direction the rank will be 1.0.  If the
 *  angle is pi, the rank is 0.0.  If the angle is 2*pi
 *  differentially small in the clockwise direction) the rank is -1.0
 *
 * The score function is:
 *
 * score = cos(theta/2)*sin(theta)/|sin(theta)|
 *
 * cos(theta/2) is equivalent to sqrt((1+cos(theta))/2)
 *
 * and because things are normalized this then becomes:
 *
 * score = sqrt((1+(ref*vect))/2)*CrossProd(ref,vect)/
 *              fabs(CrossProd(ref,vect)))
 */

double Geom2DMixIn::ScoreRelativeAngle(const Vec2D &ref,
                                           const Vec2D &vect) const
{
    if (CrossProd(ref,vect) < 0.0) {
        double dot = ref*vect ;
        if (dot <= -1.0) {  // can happen due to roundoff
            return(0.0) ;
        } else {
            return(-(sqrt(0.5 * (1.0 + dot)))) ;
        }
    } else {
        double dot = ref*vect ;
        if (dot <= -1.0) {  // can happen due to roundoff
            return(0.0) ;
        } else {
            return(sqrt(0.5 * (1.0 + ref*vect))) ;
        }
    }
}


/**  Given the coordinates of the vertices of a straight 
 *   edge and a point location, this function returns a 
 *   status to indicate if a scan line drawn from the
 *   point to minus infinity crosses the line segement.
 */

bool Geom2DMixIn::ScanCross(const Vec2D &b,
                                const Vec2D &i,
                                const Vec2D &j) const
{
    double y_min,y_max,x_min ;

    // first reject if the line segment is to the right of
    // the point

    if ((i[0] > b[0]) && (j[0] > b[0])) return(false) ;

    // find the minimum and maximum values of the edge

    if (i[1] == j[1]) {                // horizontal line segment
        return(false) ;
    } else if (i[1] <= j[1]) {
        y_max = j[1] ;   y_min = i[1] ;   x_min = i[0] ;
    } else {
        y_max = i[1] ;   y_min = j[1] ;   x_min = j[0] ;
    }

    if (i[0] == j[0]) {                // vertical line segment
        if ((b[1] >= y_min) && (b[1] < y_max))
            return(true) ;
        else
            return(false) ;
    }

    if ((b[1] == y_min) && (b[0] > x_min)) return(true) ;

    // now reject if the line segment is above or below
    // the point

    if ((b[1] > y_max) || (b[1] < y_min)) return(false) ;

    // Find the x coordinates of the intersection and see
    // if this is to the left of the point

    double m = (i[1] - j[1]) / (i[0] - j[0]) ;
    double c = i[1] - m * i[0] ;
    double xint = (b[1] - c)/m ;

    if (xint < b[0]) return(true) ; 

    return(false) ;
}


bool Geom2DMixIn::IsPointOnLine(const Vec2D &pt,
                                    const Vec2D &line0,
                                    const Vec2D &line1,
                                    double *line_param) const
{
    // find a tolerance based on a fraction of the length
    // of the line

    Vec2D dline = line1 - line0 ;
    double tol = 0.0001 * (dline).Magnitude() ;

    // find the normal to the line

    Vec2D dlinen = dline.Normalize() ;
    Vec2D norm(-dlinen.y(),dlinen.x()) ;

    // now determine if we are on the line

    if (fabs(norm*line0 - norm*pt) > tol) return(false) ;

    // if the point is on the line, determine the parametric
    // coordinate

    if (fabs(dline.x()) > fabs(dline.y())) {
        *line_param = (pt.x() - line0.x()) / dline.x() ;
    } else {
        *line_param = (pt.y() - line0.y()) / dline.y() ;
    }
    return(true) ;
}


bool Geom2DMixIn::IsPointInPoly(const Vec2D &pt,
                                    int num_verts,
                                    const Vec2D *poly_verts,
                                    double tol) const
{
    int num_cross = 0 ;

    // check to see if the given point is inside the polygon
    // described by the input verticies

    for (int i=0 ; i<num_verts ; ++i) {

        int j = (i+1) % num_verts ;

        // if the scan line goes through this vertex then
        // check to see if we cross or are changing direction
        // here.

        if ((pt.y() < poly_verts[i].y()+tol) &&
            (pt.y() > poly_verts[i].y()-tol)) {
            if (pt.x() < poly_verts[i].x()) continue ;
            if ((pt.y() > poly_verts[j].y()+tol) ||
                (pt.y() < poly_verts[j].y()-tol)) {
                int k = i==0 ? num_verts-1 : i-1 ;
                while ((pt.y() < poly_verts[k].y()+tol) &&
                       (pt.y() > poly_verts[k].y()-tol)) {
                    k = k==0 ? num_verts-1 : k-1 ;
                }
                double delt0 = poly_verts[j].y() - pt.y() ;
                double delt1 = poly_verts[k].y() - pt.y() ;
                if (delt0*delt1 < 0.0) ++num_cross ;
            } else {
                continue ;
            }

        // if the scan line goes through the next point then skip

        } else if ((pt.y() < poly_verts[j].y()+tol) &&
                   (pt.y() > poly_verts[j].y()-tol)) {

            continue ;            

        // else check for a scan line cross

        } else {

            if (ScanCross(pt,poly_verts[i],poly_verts[j])) ++num_cross ;

        }
    }

    // if we have an odd number of crossings we are inside

    return ((num_cross%2) == 1) ;
}

} // namespace


