//
// CArbBdrySeg Template Class header file
//
// Description -
//   This class implements a boundary segments spline class.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2000
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.3 $  $Date: 2001/08/08 19:33:22 $  $Author: wash $
//

#ifndef ArbBdrySeg_hh
#define ArbBdrySeg_hh

#include "ArbBSpline2D.hpp"

class CArbBdrySeg {
    public:

        // constructor for a spline segment

        CArbBdrySeg(int num,CArbCoord2D *points,
                    int inode0,double ichar_len0,
                    int inode1,double ichar_len1,
                    int iparent_id,
                    bool bi_mat_interface) :
                node0(inode0),node1(inode1),
                char_len0(ichar_len0),
                char_len1(ichar_len1),
                parent_id(iparent_id),
                interface(bi_mat_interface),
                mate(0)
                { spline = new CArbBSpline2D(num,points) ; } ;

        CArbBdrySeg(const CArbCoord2D &icoord,
                    const CArbCoord2D &itang,
                    int inode0,double ichar_len0,
                    int inode1,double ichar_len1,
                    int iparent_id,
                    bool bi_mat_interface) :
                node0(inode0),node1(inode1),
                char_len0(ichar_len0),
                char_len1(ichar_len1),
                parent_id(iparent_id),
                interface(bi_mat_interface),
                coord(icoord),tang(itang),
                spline(0),mate(0) {} ;

        CArbBdrySeg(const CArbBdrySeg &other) :
                node0(other.node0),node1(other.node1),
                char_len0(other.char_len0),
                char_len1(other.char_len1),
                parent_id(other.parent_id),
                interface(other.interface),
                coord(other.coord),tang(other.tang),mate(0)
                { spline = new CArbBSpline2D(*other.spline) ; } ;

        ~CArbBdrySeg() { if (spline != 0) delete spline ; } ;

        int GetStartNode() { return(node0) ; } ;
        int GetStopNode() { return(node1) ; } ;

        double GetStartCharLen() { return(char_len0) ; } ;
        double GetStopCharLen() { return(char_len1) ; } ;

        bool IsBiMatInterface() { return(interface) ; } ;

        int GetParentId() { return(parent_id) ; } ;

        void SetStartNode(int nd) { node0 = nd ; } ;
        void SetStopNode(int nd) { node1 = nd ; } ;

        void SetStartCharLen(double cl) { char_len0 = cl ; } ;
        void SetStopCharLen(double cl) { char_len1 = cl ; } ;

        void SetMate(CArbBdrySeg *imate) { mate = imate ; } ;

        CArbBSpline2D *GetSpline() { return(spline) ; } ;
        CArbBdrySeg *GetMate() { return(mate) ; } ;

        CArbCoord2D Evaluate(double u) const
            { return((spline == 0) ? coord : spline->Evaluate(u)) ; } ;

        CArbCoord2D Tangent(double u) const
            { if (spline != 0) {
                 CArbCoord2D v,du,d2u ;
                 spline->Derivatives(u,&v,&du,&d2u) ;
                 return((u != 1.0) ? du.Normalize() : -du.Normalize()) ;
              } else {
                 return((u != 1.0) ? tang : -tang) ;
              }
            } ;

        double TangentMag(double u) const
            { if (spline != 0) {
                 CArbCoord2D v,du,d2u ;
                 spline->Derivatives(u,&v,&du,&d2u) ;
                 return(du.Magnitude()) ;
              } else {
                 return(0.0) ;
              }
            } ;

        bool ClosestPoint(const CArbCoord2D &pt,double *u) const
            { if (spline != 0) return(spline->ClosestPoint(pt,u)) ;
              *u = 0.0 ; return(true) ; } ;

        double ApproxLength(const int npts = 100) const
            { return((spline != 0) ? spline->ApproxLength(npts) : 0.0) ; } ;

        CArbArray<double> *IntersectLine(const CArbCoord2D &pt0,
                                         const CArbCoord2D &pt1) const
            { if (spline != 0) return(spline->IntersectLine(pt0,pt1)) ;
              // currently only intersect lines
              return(0) ; }
         

        CArbBdrySeg *SubSegment(
                               double par_start,double par_stop,
                               int inode0,double ichar_len0,
                               int inode1, double ichar_len1,
                               bool bi_mat_interface) const
            {
                int Num_Segs = spline->GetNumSegs() ;
                CArbCoord2D *points = new CArbCoord2D[Num_Segs+1] ;
                for (int i=0 ; i<=Num_Segs ; ++i) {
                    double f = double(i)/double(Num_Segs) ;
                    double u = par_start + f*(par_stop - par_start) ;
                    points[i] = Evaluate(u) ;
                }
                CArbBdrySeg *nspl = new CArbBdrySeg(
                                             Num_Segs+1,points,
                                             inode0,ichar_len0,
                                             inode1,ichar_len1,
                                             parent_id,
                                             bi_mat_interface) ;
                delete [] points ;
                return(nspl) ;
            } ;


        void Reverse()
            { int tmp0 = node0 ; node0 = node1 ; node1 = tmp0 ;
              double tmp1 = char_len0 ; char_len0 = char_len1 ;
              char_len1 = tmp1 ;
              if (spline != 0) spline->Reverse() ;
              else tang = -tang ;
            } ;

        void Print()
            { if (spline != 0) spline->Print() ; } ;

    private:

        int node0, node1 ;
        double char_len0, char_len1 ;
        int parent_id ;
        bool interface ;
        CArbCoord2D coord ;
        CArbCoord2D tang ;
        CArbBSpline2D *spline ;
        CArbBdrySeg *mate ;
} ;

#endif

