
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "ArbBSpline2D.hpp"


// %(CArbBSpline2D::CArbBSpline2D-constructor-|-int-|-CArbCoord2D-|*)
/* ++ ----------------------------------------------------------
**
**    CArbBSpline2D - Build a spline interpolates the input points
**
**      CArbBSpline2D(
**              int         num,
**              CArbCoord2D *points)
**
**        num    - (in)  number of interpolation points 
**        points - (in)  list of interpolation points 
**
**      Description: This constructor generates a spline that 
**          interpolates the given input points. This is a two-step 
**          process. First, a cardinal spline is fit through the points 
**          (a unit parametric step between each input point). In a 
**          second step, the spline is reparameterized so that the 
**          mapping between spline parametric space and the arc length in 
**          Cartesian space is nearly constant along the spline. 
**
**
** -- */


CArbBSpline2D::CArbBSpline2D(int num,CArbCoord2D *points)
{
    Coord = new CArbCoord2D[num+2] ;
    U = new double[num+4] ;
    Num_Segs = num - 1 ;
    GetControlNet(Num_Segs,points) ;
    UniformParameterization() ;
}




// %(CArbBSpline2D::CArbBSpline2D-constructor-|-CArbBSpline2D-const|&)
/* ++ ----------------------------------------------------------
**
**    CArbBSpline2D - copy constructor
**
**      CArbBSpline2D(const CArbBSpline2D &other)
**
**        other - (in)  object to copy 
**
**      Description: This is the copy constructor for the class. 
**
**
** -- */


CArbBSpline2D::CArbBSpline2D(const CArbBSpline2D &other)
{
    Num_Segs = other.Num_Segs ;

    Coord = new CArbCoord2D[Num_Segs+3] ;
    for (int i=0 ; i<Num_Segs+3 ; ++i)
        Coord[i] = other.Coord[i] ;

    U = new double[Num_Segs+5] ;
    for (int j=0 ; j<Num_Segs+5 ; ++j) U[j] = other.U[j] ;
}




// %(CArbBSpline2D::operator_=-CArbBSpline2D-|-CArbBSpline2D-const|&)
/* ++ ----------------------------------------------------------
**
**    operator_= - assignment operator
**
**      CArbBSpline2D operator = (const CArbBSpline2D &other)
**
**        other - (in)  object to copy for the assignment 
**
**      Description: This function implements an assignment operator for 
**          this object. 
**
**      Return Value: returns the updated object 
**
**
** -- */


CArbBSpline2D CArbBSpline2D::operator = (const CArbBSpline2D &other)
{
    if (Num_Segs > 0) {
        delete [] Coord ;
        delete [] U ;
    }

    Num_Segs = other.Num_Segs ;

    Coord = new CArbCoord2D[Num_Segs+3] ;
    for (int i=0 ; i<Num_Segs+3 ; ++i)
        Coord[i] = other.Coord[i] ;

    U = new double[Num_Segs+5] ;
    for (int j=0 ; j<Num_Segs+5 ; ++j) U[j] = other.U[j] ;
    return(*this) ;
}




// %(CArbBSpline2D::CArbBSpline2D-destructor-|~) 
/* ++ ----------------------------------------------------------
**
**    CArbBSpline2D - destructor
**
**      ~CArbBSpline2D()
**
**      Description: this is the destructor function 
**
**
** -- */


CArbBSpline2D::~CArbBSpline2D()
{
    delete [] Coord ;
    delete [] U ;
}




// %(CArbBSpline2D::Evaluate-CArbCoord2D-|^const-double-const|)
/* ++ ----------------------------------------------------------
**
**    Evaluate - evaluate the spline value
**
**      CArbCoord2D Evaluate(const double nu) const
**
**        nu - (in)  evaluation coordinate 
**
**      Description: This function evaluates the spline value at the 
**          given parametric coordinate. 
**
**      Return Value: The spline value at this point. 
**
**
** -- */


CArbCoord2D CArbBSpline2D::Evaluate(const double nu) const
{
    //  Reference:
    //  Bohm, W., Farin G., and Kahmann J., "A survay of curve and 
    //  surface methods in CAGD," Computer Aided Geometric Design,
    //  Vol. 1, 1984, pp. 1-60.

    int i, j ;

    // Check make sure that the evaluation point is within the
    // supported section of the spline

    double uevl = nu * Num_Segs ;
    assert((uevl >= U[2]) && (uevl <= U[Num_Segs+2])) ;

    // Find the interval that contains ueval

    int u_int = 2 ;
    while (uevl >= U[u_int+1]) u_int++ ;
    if (u_int >= Num_Segs+2) u_int-- ;
    int evl_s = u_int-2 ;

    // Evaluate the spline using the deBoor recurance formulat
    // Compute the three first level coordinates

    double alpha ;
    CArbCoord2D x1[3] ;

    for (i=0; i<3; i++) {
        if ((U[evl_s+i+3]-U[evl_s+i]) != 0.0) 
            alpha = (uevl - U[evl_s+i]) / (U[evl_s+i+3]-U[evl_s+i]) ;
        else
            alpha = 0.0;
        x1[i] = alpha * Coord[evl_s+i+1]  +  (1.0 - alpha) * Coord[evl_s+i]  ;
    }

    // Compute the two second level coordinates

    CArbCoord2D x2[2] ;

    for (i=1,j=0; i<3; i++,j++) {
        if ((U[evl_s+i+2]-U[evl_s+i]) != 0.0)	
            alpha = (uevl - U[evl_s+i]) / (U[evl_s+i+2]-U[evl_s+i]) ;
        else
            alpha = 0.0 ;
        x2[j] = alpha * x1[j+1]  +  (1.0 - alpha) * x1[j] ;
    }

    // Compute the spline value

    if ((U[evl_s+3]-U[evl_s+2]) != 0.0)
        alpha = (uevl - U[evl_s+2]) / (U[evl_s+3]-U[evl_s+2]) ;
    else
        alpha = 0.0;

    return(alpha * x2[1]  +  (1.0 - alpha) * x2[0]) ;
}




// %(CArbBSpline2D::Derivatives-void-|^const-double-const|-CArbCoord2D-|*-CArbCoord2D-|*-CArbCoord2D-|*)
/* ++ ----------------------------------------------------------
**
**    Derivatives - compute spline derivative
**
**      void Derivatives(
**              const double nu,
**              CArbCoord2D  *value,
**              CArbCoord2D  *der_u,
**              CArbCoord2D  *der_u2) const
**
**        nu     - (in)  evaluation coordinate 
**        value  - (in)  spline value 
**        der_u  - (in)  first derivative 
**        der_u2 - (in)  second derivative 
**
**      Description: This functions evaluates the spline's value and its 
**          first and second derivatives with respect to the parametric 
**          coordinate at the specified parametric coordinate. 
**
**
** -- */


void CArbBSpline2D::Derivatives(const double nu,CArbCoord2D *value,
                        CArbCoord2D *der_u,CArbCoord2D *der_u2) const
{
    int i, j ;

    // Check make sure that the evaluation point is within the
    // supported section of the spline

    double uevl = nu * Num_Segs ;
    assert((uevl >= U[2]) && (uevl <= U[Num_Segs+2])) ;

    // Find the interval that contains ueval

    int u_int = 2 ;
    while (uevl > U[u_int+1]) u_int++ ;
    if (u_int >= Num_Segs+2) u_int-- ;
    int evl_s = u_int-2 ;

    //  First reduction of the coordinates in u direction

    double alpha, beta ;
    CArbCoord2D x1[3],xdu[3] ;
    for (i=0; i<3; i++) {
        if ((U[evl_s+i+3]-U[evl_s+i]) != 0.0) {
            alpha = (uevl - U[evl_s+i]) / (U[evl_s+i+3]-U[evl_s+i]);
            beta  = 3.0 / (U[evl_s+i+3]-U[evl_s+i]);
        } else {
            alpha = 0.0;
            beta = 0.0;
        }
        x1[i] = alpha * Coord[evl_s+i+1]  +  (1.0 - alpha) * Coord[evl_s+i];
        xdu[i] = beta * (Coord[evl_s+i+1] - Coord[evl_s+i]);
    }

    // Second reduction of the coordinates

    CArbCoord2D x2[2],xdu1[2],xddu[2] ;
    for (i=1,j=0; i<3; i++,j++) {
        if ((U[evl_s+i+2]-U[evl_s+i]) != 0.0) {
            alpha = (uevl - U[evl_s+i]) / (U[evl_s+i+2]-U[evl_s+i]) ;
            beta = 2.0 / (U[evl_s+i+2]-U[evl_s+i]) ;
        } else {
            alpha = 0.0 ;
            beta = 0.0 ;
        }
        x2[j] = alpha * x1[j+1]  +  (1.0 - alpha) * x1[j] ;
        xdu1[j] = alpha * xdu[j+1] + (1.0 - alpha) * xdu[j] ;
        xddu[j] = beta * (xdu[j+1] - xdu[j]) ;
    }


    // Third reduction (spline values and derivatives)

    if ((U[evl_s+3]-U[evl_s+2]) != 0.0)
        alpha = (uevl - U[evl_s+2]) / (U[evl_s+3]-U[evl_s+2]) ;
    else
        alpha = 0.0 ;

    *value = alpha * x2[1]  +  (1.0 - alpha) * x2[0] ;
    *der_u  = (alpha * xdu1[1] + (1.0 - alpha) * xdu1[0]) * Num_Segs ;
    *der_u2 = (alpha * xddu[1] + (1.0 - alpha) * xddu[0]) * Num_Segs ;
}




// %(CArbBSpline2D::GetControlNet-void-|-int-|-CArbCoord2D-|*)
/* ++ ----------------------------------------------------------
**
**    GetControlNet - determine the control points for a spline
**
**      void GetControlNet(
**              int         ns,
**              CArbCoord2D *points)
**
**        ns     - (in)  number of input points 
**        points - (in)  input points 
**
**      Description: This function determines a set of control points 
**          (control net) for a spline. These points are determined such 
**          that the resulting spline will interpolate the input points. 
**          This is done by building an solving a system of linear 
**          equations. 
**
**
** -- */


void CArbBSpline2D::GetControlNet(int ns,CArbCoord2D *points)
{
    /* Barsky, B. A., "A Method Describing Curved Surfaces by 
       Transforming between Interpolatory and B-spline Representations,"
       MSc Thesis, Cornell University, 1979.

       Barsky, B. A. and Greenberg, D. P., " Determining a Set of B-spline
       Control Vertices to Generate an Interpolating Surface," Computer 
       Graphics and Image Processing, No. 14, 1979, pp. 203-226.

       this routine finds an appropriate set of cubic 
       B-spline curve control points by interpolating a set of given points 
       on the curve.  

       Since an (n+1) array of control points of a cubic B-spline curve
       defines (n-2) segments of curve, this routine extends the
       control polygon at the end of the curve.  In this way the user does 
       not need to extend the control polygon or knot vectors.
       The control polygon is extended with pseudo control points ("phantom
       vertices") at the end of the curve so that the curve tangent directions
       agree with the directions of the end control segments.  The pseudo
       control points are defined in terms of the original set of control
       points in such a manner so as to satisfy the boundary condition of
       having the appropriate second parametric derivative at the end points
       equal to zero.
       An uniform (periodic) knot vector is generated and the knot values 
       are -2,-1,0,1,...,np,np+1,np+2.
    */

    int i ;
    int n = ns + 2 ;  // # of control segments

    // generate the knot vector

    for (i=0 ; i < (ns+5) ; i++) U[i] = double(i-2) ;

    // set the end points of the curve

    Coord[1]   = points[0] ;
    Coord[n-1] = points[ns] ;

    // if we have only one segment, we need only to extend it

    if (ns == 1) {
        Coord[0] = 2.0*Coord[1]   - Coord[2] ;
        Coord[n] = 2.0*Coord[n-1] - Coord[n-2] ;
        return ;
    }

    // get space for and fill an auxiliar array (beta[0] is never used)

    double *beta = new double[n] ;
    beta[1] = 0.25 ;
    for (i = 2 ; i<int(n) ; i++)
	beta[i] = 1.0 / (4.0 - beta[i-1]) ;

    // get the rhs terms

    double coeff = 6.0 ;

    if (ns > 2) {
	Coord[2] = coeff * points[1] - Coord[1] ;
	for (i=2 ; i<(ns-1) ; i++) {
	    Coord[i+1] = coeff * points[i] ;
	}
	Coord[n-2] = coeff * points[ns-1] - Coord[n-1] ;
    } else {
	Coord[2] = coeff * points[1] - Coord[1] - Coord[3] ;
    }

    // transform the rhs

    Coord[2] = beta[1] * Coord[2] ;
    for (i=2 ; i<ns ; i++) {
	Coord[i+1] = beta[i]*(Coord[i+1] - Coord[i]) ;
    }

    // solve (back-substitution)

    for (i=(ns-2) ; i>0 ; i--) {
	Coord[i+1] = Coord[i+1] - beta[i] * Coord[i+2] ;
    }
    delete [] beta;

    // now extend the control values once in outward direction.

    Coord[0] = 2.0*Coord[1]   - Coord[2] ;
    Coord[n] = 2.0*Coord[n-1] - Coord[n-2] ;
}




// %(CArbBSpline2D::IntersectLine-bool-|^const-CArbCoord2D-const|&-CArbCoord2D-const|&-double-|*-double-|*) 
/* ++ ----------------------------------------------------------
**
**    IntersectLine - intersect the spline with a line
**
**      bool IntersectLine(
**              const CArbCoord2D &pt0,
**              const CArbCoord2D &pt1,
**              double            *u,
**              double            *t) const
**
**        pt0 - (in)  first point on line 
**        pt1 - (in)  second point on line 
**        u   - (out) spline intersection coordinate 
**        t   - (out) line intersection coordinate 
**
**      Description: This function attempts to fine the point of 
**          intersection between the B spline and a straight line. The 
**          line is defined by two points. If succesful, the routine 
**          returns the parametric coordinates of the intersection for 
**          both the line and the spline. For the line, the first and 
**          second points are assumed to have a parametric coordinates of 
**          0 and 1, respectively. 
**
**      Return Value: True if an intersection point is found, false 
**          otherwise. 
**
**
** -- */


#define MAX_ITER 100

CArbArray<double> *CArbBSpline2D::IntersectLine(const CArbCoord2D &pt0,
                                       const CArbCoord2D &pt1) const
{
    /* 
       First we divide the spline into two times the number of
       control points segments.  We step through these segments.
       For each segment we check to see two end points are on
       opposite sides of the line.  If so we have bracketed a
       crossing and we use a newton algorithm to find the crossing.

       We also check the tangent to the spline on each end of
       the segment to see if it changes direction relative to
       the line.  If so there may be two crossings in the
       segment.  Here we use bisection algorithm until we
       find the point where the segment tangent is parallel to
       the line or until we have bracketed a root.
     */

    CArbArray<double> *results = 0 ;
    
    // determine the tangent and normal to the line

    CArbCoord2D tang = pt1 - pt0 ;
    CArbCoord2D tmp(-tang.y(),tang.x()) ;
    CArbCoord2D norm = tmp.Normalize() ;

    // find the third coefficient for the line equation

    double coef = -norm.x()*pt0.x() - norm.y()*pt0.y() ;
    double on_tol = fabs(0.001 * coef) ;

    // now loop through the spline looking for a segment
    // that crosses the line

    double u1, u0 = 0.0 ;
    CArbCoord2D sp1, sp0 = Evaluate(0.0) ;
    double side1, side0 = norm * sp0 + coef ;

    for (int i=1 ; i<=Num_Segs*2 ; ++i) {
        u1 = double(i) / double(Num_Segs*2) ;
        sp1 = Evaluate(u1) ;
        side1 = norm * sp1 + coef ;

        if (side0*side1 <= 0.0) {

            // we have bracketed a crossing
            // so find it using a Newton algorithm

            if ((fabs(side0) > on_tol) || (fabs(side1) > on_tol))
                FindCross(u0,u1,norm,tang,pt0,coef,&results) ;
        } else {

            // check the direction of the spline tangent relative
            // to the line

            CArbCoord2D ntan = tang.Normalize() ;
            CArbCoord2D s0, ds0, ds20 ;
            CArbCoord2D s1, ds1, ds21 ;
            Derivatives(u0,&s0,&ds0,&ds20) ;
            Derivatives(u1,&s1,&ds1,&ds21) ;
            ds0 = ds0.Normalize() ;
            ds1 = ds1.Normalize() ;

            if (CrossProduct(ds0,ntan)*CrossProduct(ds1,ntan) < 0.0) {
                double cross0 = CrossProduct(ds0,ntan) ;
                double tol = 0.001 ;
                double u_low, u_high, u_mid ;
                if (cross0 < 0.0) {
                    u_low = u0 ;  u_high = u1 ;
                } else {
                    u_low = u1 ;  u_high = u0 ;
                }
                u_mid = 0.5*(u_low + u_high) ;
                CArbCoord2D sm, dsm, ds2m ;
                Derivatives(u_mid,&sm,&dsm,&ds2m) ;
                dsm = dsm.Normalize() ;
                double cross = CrossProduct(dsm,tang) ;

                int iter = 0 ;
                while ((fabs(cross) > tol) && (iter < MAX_ITER)) {
                    side1 = norm * Evaluate(u_mid) + coef ;

                    // check to see if we have a crossing

                    if (side0*side1 <= 0) {
                        FindCross(u0,u_mid,norm,tang,pt0,coef,&results) ;
                        FindCross(u_mid,u1,norm,tang,pt0,coef,&results) ;
                        break ;
                    }

                    if (cross < 0.0) {
                        u_low = u_mid ;
                    } else {
                        u_high = u_mid ;
                    }
                    u_mid = 0.5*(u_low + u_high) ;
                    Derivatives(u_mid,&sm,&dsm,&ds2m) ;
                    cross = CrossProduct(dsm,ntan) ;
                    ++iter ;
                }
            }
        }

        u0 = u1 ;
        sp0 = sp1 ;
        side0 = side1 ;
    } 

    return(results) ;
}


void CArbBSpline2D::FindCross(double u0,double u1,
                              const CArbCoord2D &norm,
                              const CArbCoord2D &tang,
                              const CArbCoord2D &pt0,
                              double coef,
                              CArbArray<double> **results) const
{
    double u_new = u0 ;
    double u_tol = 0.001 ;
    double u = .0f;
    CArbCoord2D s ;

    for (int iter=0 ; iter < MAX_ITER ; ++iter) {
        u = u_new ;

        // evaluate the spline at this point

        CArbCoord2D ds, ds2 ;
        Derivatives(u,&s,&ds,&ds2) ;

        double delta_u = (norm * s + coef)/(norm * ds) ;
        u_new = u - delta_u ;

        if (u_new < u0) u_new = u0 ;
        if (u_new > u1) u_new = u1 ;

        // check for convergance

        if (fabs(delta_u) < u_tol) break ;
    }

    if (*results == 0) {
        *results = new CArbArray<double> ;
    } else {
        if (((*results)->At((*results)->NumEntries()-2) < u+u_tol ) &&
            ((*results)->At((*results)->NumEntries()-2) > u-u_tol))
           return ;
    }

    (*results)->InsertAtEnd(u) ;
    double par ;
    if (fabs(tang.x()) > fabs(tang.y()))
        par = (s.x() - pt0.x()) / tang.x() ;
    else
        par = (s.y() - pt0.y()) / tang.y() ;
    (*results)->InsertAtEnd(par) ;
}



// %(CArbBSpline2D::ClosestPoint-bool-|^const-CArbCoord2D-const|&-double-|*) 
/* ++ ----------------------------------------------------------
**
**    ClosestPoint - find the spline point closest to the input point
**
**      bool ClosestPoint(
**              const CArbCoord2D &pt,
**              double            *u) const
**
**        pt - (in)  search point 
**        u  - (out) parametric coordinate of the close point 
**
**      Description: Given an input point, this routine attempts to find 
**          the point on the spline that is closest. 
**
**      Return Value: A flag is returned to indicate if the search for a 
**          point was successful. In theory, the search should always 
**          work, but it is a nonlinear search, an in some cases the 
**          search may not converge. 
**
**
** -- */


bool CArbBSpline2D::ClosestPoint(const CArbCoord2D &pt,double *uf) const
{
    // loop through the spline looking for a point that
    // is closest to the given point

    CArbCoord2D sp = Evaluate(0.0) ;
    CArbCoord2D delt = sp - pt ;
    double dist = delt.Magnitude() ;
    double u_min = 0.0 ;
    double d_min = dist ;

    for (int i=1 ; i<=10 ; ++i) {
        double u = i / 10.0 ;
        sp = Evaluate(u) ;
        delt = sp - pt ;
        dist = delt.Magnitude() ;

        if (dist < d_min) {
            u_min = u ;
            d_min = dist ;
        }
    } 

    // using the min as the starting point, iterate until
    // we converge

    double u = .0f, u_new = u_min ;
    double u_tol = 0.001 ;

    for (int iter=0 ; iter < MAX_ITER ; ++iter) {

        u = u_new ;

        // evaluate the spline at this point

        CArbCoord2D s, ds, ds2 ;
        Derivatives(u,&s,&ds,&ds2) ;
        delt = s - pt ;

        double delta_u = (delt * ds) / ((ds * ds)+(delt * ds2)) ;
        u_new = u - delta_u ;

        if (u_new < 0.0) u_new = 0.0 ;
        if (u_new > 1.0) u_new = 1.0 ;

        // check for convergance

        if (fabs(delta_u) < u_tol) {
            *uf = u_new ;
//            return(fabs(delt*ds) < u_tol ? true : false) ;
            return(true) ;
        }
    }
    *uf = u ;
    return(false) ;
}


CArbBSpline2D *CArbBSpline2D::SubSegment(double par_start,
                                         double par_stop) const
{
    CArbCoord2D *points = new CArbCoord2D[Num_Segs+1] ;
    for (int i=0 ; i<=Num_Segs ; ++i) {
        double u = par_start + (double(i)/double(Num_Segs)) *
                       (par_stop - par_start) ;
        points[i] = Evaluate(u) ;
    }
    CArbBSpline2D *nspl = new CArbBSpline2D(Num_Segs+1,points) ;
    delete [] points ;
    return(nspl) ;
}


double CArbBSpline2D::ApproxLength(const int npts) const
{
    double len = 0.0 ;
    CArbCoord2D pt0 = Evaluate(0.0) ;
    for (int i=0 ; i<npts-1 ; ++i) {
        double u = double(i+1)/double(npts-1) ;
        CArbCoord2D pt1 = Evaluate(u) ;
        len += (pt1-pt0).Magnitude() ;
        pt0 = pt1 ;
    }
    return(len) ;
}


// %(CArbBSpline2D::UniformParameterization-void-|-int-const|) 
/* ++ ----------------------------------------------------------
**
**    UniformParameterization - reparametrize for a uniform mapping
**
**      void UniformParameterization(const int npts = 10000)
**
**        npts - (in)  number of points to use in the reevaluation, 
**                     should be 'large' 
**
**      Description: This function reparameterizes a spline so that the 
**          mapping from the spline's parametric space to the arc length 
**          in the Cartesian space is nearly constant along the spline. 
**
**
** -- */


void CArbBSpline2D::UniformParameterization(const int npts)
{
    int i, cur ;
    if (Num_Segs == 1) return ;

    // compute the approximate total arc length

    double arc_len = 0.0 ;
    CArbCoord2D pt1, pt0 = Evaluate(0.0) ;

    for (i=1 ; i<npts ; ++i) {
        double u = double(i) / double(npts - 1) ;
        pt1 = Evaluate(u) ;
        CArbCoord2D delta = pt1 - pt0 ;
        arc_len += delta.Magnitude() ;
        pt0 = pt1 ;
    }

    double thresh_inc = arc_len / Num_Segs ;
    double cur_threshold = thresh_inc ;

    CArbCoord2D *new_pts = new CArbCoord2D[Num_Segs + 1] ;

    // loop through and compute the new curve points

    arc_len = 0.0 ;
    pt0 = new_pts[0] = Evaluate(0.0) ;
    for (i=1,cur=1 ; i<npts ; ++i) {
        double u = double(i) / double(npts - 1) ;
        pt1 = Evaluate(u) ;
        CArbCoord2D delta = pt1 - pt0 ;
        arc_len += delta.Magnitude() ;
        pt0 = pt1 ;
        if (arc_len >= cur_threshold) {
            new_pts[cur] = pt0 ;
            ++cur ;
            if (cur >= Num_Segs) break ;
            cur_threshold += thresh_inc ;
        }
    }
    new_pts[Num_Segs] = Evaluate(1.0) ;

    GetControlNet(Num_Segs,new_pts) ;
    delete [] new_pts ;
}


void CArbBSpline2D::Reverse()
{
    // reverse the control points

    int num = Num_Segs+3 ;
    for (int i=0 ; i<num/2 ; ++i) {
        int j = num - i - 1 ;
        CArbCoord2D tmp = Coord[i] ;
        Coord[i] = Coord[j] ;
        Coord[j] = tmp ;
    }

    // reverse the knot vector

    num = Num_Segs + 5 ;
    double range = U[num-1] - U[0] ;
    double *new_U = new double[num] ;

    for (int j=0 ; j<num ; ++j) {
        new_U[j] = U[0] + range*(1.0 - (U[num-1]-U[j])/range) ;
    }
    delete U ;
    U = new_U ;
}

void CArbBSpline2D::Print()
{
    fprintf(stderr,"\nNum: %d\n",Num_Segs) ;
    for (int i=0 ; i<Num_Segs+3 ; ++i) {
        fprintf(stderr,"    %21.15g %21.15g\n",Coord[i].x(),Coord[i].y()) ;
    }
    for (int j=0 ; j<Num_Segs+5 ; ++j) {
        fprintf(stderr,"    %21.15g\n",U[j]) ;
    }
}

