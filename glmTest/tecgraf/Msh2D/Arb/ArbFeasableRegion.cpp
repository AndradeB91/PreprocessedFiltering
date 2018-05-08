//
// CArbFeasableRegion Class implementation file
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
//   $Revision: 1.4 $  $Date: 2001/02/15 14:19:40 $  $Author: wash $
//

#include "ArbFeasableRegion.hpp"

// %(CArbFeasableRegion::CArbFeasableRegion-constructor-|-double-|)
/* ++ ----------------------------------------------------------
**
**    CArbFeasableRegion - constructor 
**
**      CArbFeasableRegion(double tolerance)
**
**        tolerance - (in)  intersection tolerance 
**
**      Description: This is a constructor for a FeasabilityRegion 
**
**
** -- */

CArbFeasableRegion::CArbFeasableRegion(double tolerance)
{
    Large = tolerance * 1e30 ;
    Tol = tolerance ;
}




// %(CArbFeasableRegion::CArbFeasableRegion-constructor-|-CArbFeasableRegion-const|&)
/* ++ ----------------------------------------------------------
**
**    CArbFeasableRegion - copy constructor 
**
**      CArbFeasableRegion(const CArbFeasableRegion &other)
**
**        other - (in)  object to copy 
**
**      Description: This is the copy constructor for a 
**          FeasabilityRegion 
**
**
** -- */

CArbFeasableRegion::CArbFeasableRegion(const CArbFeasableRegion &other)
{
    Large = other.Large ;
    Tol = other.Tol ;
    for (int i=0 ; i<other.Sides.NumEntries() ; ++i)
        Sides[i] = other.Sides[i] ;
}




// %(CArbFeasableRegion::operator_=-CArbFeasableRegion-|-CArbFeasableRegion-const|&)
/* ++ ----------------------------------------------------------
**
**    operator_= - assignment operator 
**
**      CArbFeasableRegion operator = (const CArbFeasableRegion &other)
**
**        other - (in)  object to copy 
**
**      Description: This is the assignment operator for a 
**          FeasabilityRegion 
**
**      Return Value: the updated object 
**
**
** -- */

CArbFeasableRegion CArbFeasableRegion::operator =
     (const CArbFeasableRegion &other)
{
    Large = other.Large ;
    Tol = other.Tol ;
    Sides.Clear() ;
    for (int i=0 ; i<other.Sides.NumEntries() ; ++i)
        Sides[i] = other.Sides[i] ;
    return(*this) ;
}




// %(CArbFeasableRegion::AddConstraint-void-|-CArbCoord2D-|-CArbCoord2D-|)
/* ++ ----------------------------------------------------------
**
**    AddConstraint - add a linear constraint 
**
**      void AddConstraint(
**              CArbCoord2D pt0,
**              CArbCoord2D pt1)
**
**        pt0 - (in)  first point on constraint line 
**        pt1 - (in)  second point on constraint line 
**
**      Description: This function adds a linear constraint to the 
**          feasable region. 
**
**
** -- */

#define PARALLEL_TOL 1.0e-12

void CArbFeasableRegion::AddConstraint(CArbCoord2D pt0,CArbCoord2D pt1)
{
    int j, k ;

    // create a new side

    PolySide nside(pt0,pt1,Large) ;
    Sides.InsertAtEnd(nside) ;
    k = Sides.NumEntries() - 1 ;

    // check for intersections

    for (j=0 ; j<Sides.NumEntries()-1 ; ++j) {

        // find the intersection 

        double denom = -Sides[k].e.y()*Sides[j].e.x() +
                        Sides[k].e.x()*Sides[j].e.y() ;

        if (fabs(denom) > PARALLEL_TOL) {
            double t0 = -(Sides[j].d.y()*Sides[j].e.x() -
                          Sides[k].d.y()*Sides[j].e.x() +
                          Sides[j].e.y()*Sides[k].d.x() -
                          Sides[j].e.y()*Sides[j].d.x()) / denom ;
            double t1 = -(Sides[k].e.y()*Sides[k].d.x() -
                          Sides[k].e.y()*Sides[j].d.x() -
                          Sides[k].e.x()*Sides[k].d.y() +
                          Sides[k].e.x()*Sides[j].d.y()) / denom ;
            double cross = Sides[k].e.x()*Sides[j].e.y() -
                           Sides[k].e.y()*Sides[j].e.x() ;

            // update the polygon side information checking all
            // the special cases

            int upd = 0 ;

            if ((cross > 0.0) && (t0 <= Sides[k].max-Tol)) {
                if ((t1 > Sides[j].min+Tol) &&
                    (t1 < Sides[j].max-Tol)) {
                    upd = 1 ;
                } else if ((Sides[j].prev >= 0) &&
                           (t1 > Sides[j].min-Tol) &&
                           (t1 < Sides[j].min+Tol)) {
                    double cross2 = 
                        Sides[k].e.x()*Sides[Sides[j].prev].e.y() -
                        Sides[k].e.y()*Sides[Sides[j].prev].e.x() ;
                        if (cross2 > 0.0) upd = 1 ;
                } else if ((Sides[j].next >= 0) &&
                           (t1 > Sides[j].max-Tol) &&
                           (t1 < Sides[j].max+Tol)) {
                    double cross2 = 
                        Sides[k].e.x()*Sides[Sides[j].next].e.y() -
                        Sides[k].e.y()*Sides[Sides[j].next].e.x() ;
                        if (cross2 > 0.0) upd = 3 ;
                }
            } else if ((cross < 0.0) && (t0 >= Sides[k].min+Tol)) {
                if ((t1 > Sides[j].min+Tol) &&
                    (t1 < Sides[j].max-Tol)) {
                    upd = 2 ;
                } else if ((Sides[j].prev >= 0) &&
                           (t1 > Sides[j].min-Tol) &&
                           (t1 < Sides[j].min+Tol)) {
                    double cross2 = 
                        Sides[k].e.x()*Sides[Sides[j].prev].e.y() -
                        Sides[k].e.y()*Sides[Sides[j].prev].e.x() ;
                        if (cross2 < 0.0) upd = 4 ;
                } else if ((Sides[j].next >= 0) &&
                           (t1 > Sides[j].max-Tol) &&
                           (t1 < Sides[j].max+Tol)) {
                    double cross2 = 
                        Sides[k].e.x()*Sides[Sides[j].next].e.y() -
                        Sides[k].e.y()*Sides[Sides[j].next].e.x() ;
                        if (cross2 < 0.0) upd = 2 ;
                }
            }

            if (upd == 1) {

                Sides[k].max = t0 ;
                Sides[j].min = t1 ;
                if ((Sides[k].next != -1) &&
                    (Sides[Sides[k].next].prev == k))
                    Sides[Sides[k].next].prev = -1 ;
                Sides[k].next = j ;
                if ((Sides[j].prev != -1) &&
                    (Sides[Sides[j].prev].next == j))
                    Sides[Sides[j].prev].next = -1 ;
                Sides[j].prev = k ;

            } else if (upd == 2) {

                Sides[k].min = t0 ;
                Sides[j].max = t1 ;
                if ((Sides[k].prev != -1) &&
                    (Sides[Sides[k].prev].next == k))
                    Sides[Sides[k].prev].next = -1 ;
                Sides[k].prev = j ;
                if ((Sides[j].next != -1) &&
                    (Sides[Sides[j].next].prev == j))
                    Sides[Sides[j].next].prev = -1 ;
                Sides[j].next = k ;

            } else if (upd == 3) {

                Sides[k].max = t0 ;
                if ((Sides[k].next != -1) &&
                    (Sides[Sides[k].next].prev == k))
                    Sides[Sides[k].next].prev = -1 ;
                Sides[k].next = Sides[j].next ;
                Sides[Sides[j].next].prev = k ;
                Sides[j].prev = -1 ;
                Sides[j].next = -1 ;

            } else if (upd == 4) {

                Sides[k].min = t0 ;
                if ((Sides[k].prev != -1) &&
                    (Sides[Sides[k].prev].next == k))
                    Sides[Sides[k].prev].next = -1 ;
                Sides[k].prev = Sides[j].prev ;
                Sides[Sides[j].prev].next = k ;
                Sides[j].prev = -1 ;
                Sides[j].next = -1 ;

            }
        }
    }
}




// %(CArbFeasableRegion::GetVertices-CArbArray-|<CArbCoord2D>*)
/* ++ ----------------------------------------------------------
**
**    GetVertices - get the polygon vertices 
**
**      CArbArray <CArbCoord2D>*GetVertices()
**
**      Description: This function returns the verticies of the 
**          feasable region polygon. 
**
**      Return Value: An array of polygon vertices. Ownership of the 
**          array passes to the client, who must eventually call 
**          delete. 
**
**
** -- */

CArbArray<CArbCoord2D> *CArbFeasableRegion::GetVertices()
{
    // extract the points that make up the polygon

    CArbArray<CArbCoord2D> *poly = new CArbArray<CArbCoord2D> ;

    for (int j=0 ; j<Sides.NumEntries() ; ++j) {
        if ((Sides[j].next != -1) && (Sides[j].prev != -1)) {
            bool invalid = false ;
            int first = j ;
            int cur = j ;
            poly->InsertAtEnd(CArbCoord2D(Sides[cur].d +
                                 Sides[cur].max*Sides[cur].e)) ;
            while (Sides[cur].next != first) {
                cur = Sides[cur].next ;
                if ((Sides[cur].next == -1) || (Sides[cur].prev == -1)) {
                    invalid = true ;
                    break ;
                } else {
                    poly->InsertAtEnd(
                        CArbCoord2D(Sides[cur].d + 
                            Sides[cur].max*Sides[cur].e)) ;
                }
            }
            if (!invalid) break ;
        }
        poly->Clear() ;
    }
    return(poly) ;
}
