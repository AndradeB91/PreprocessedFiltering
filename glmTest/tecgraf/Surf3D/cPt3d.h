// cPt3d.h
// Cpoint3d and cBbox3d definition and methods 

#ifndef cPt3d_h
#define cPt3d_h

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>

class cPoint3d
{
 public:
  double x,y,z;
  cPoint3d() {x=y=z=0;}
  cPoint3d(const cPoint3d& p){x=p.x;y=p.y;z=p.z;}
  cPoint3d(const cPoint3d *p){x=p->x;y=p->y;z=p->z;}
  cPoint3d(const double x0, const double y0, const double z0) {x=x0;y=y0;z=z0;}
  ~cPoint3d(){}

  double Distance(const cPoint3d& p){return sqrt((p.x-x)*(p.x-x)+
                                                (p.y-y)*(p.y-y)+
                                                (p.z-z)*(p.z-z));}
  double Dot(const cPoint3d& p0) const { return(x * p0.x + y * p0.y + z * p0.z); }
  cPoint3d Cross(const cPoint3d& p0) const
   { return cPoint3d((y * p0.z) - (z * p0.y),
                    (z * p0.x) - (x * p0.z),
                    (x * p0.y) - (y * p0.x)); 
  }
  double Mix(const cPoint3d& p0,const cPoint3d& p1) const
  { 
   return(Dot(p0.Cross(p1))); 
  }
  cPoint3d& operator +=(const cPoint3d& p) {x+=p.x;y+=p.y;z+=p.z;return *this;}
  cPoint3d& operator -=(const cPoint3d& p) {x-=p.x;y-=p.y;z-=p.z;return *this;}
  cPoint3d& operator *=(double v) {x*=v;y*=v;z*=v;return *this;}
  double   Len() const                   {return sqrt(Dot(*this));}
  int      Normalize() { double l=Len(); if(l>0) { x/=l; y/=l; z/=l; return 1; } return 0; }
  int      Normalize(double tol) 
  { 
   double l=Len();
   if (l == 0) return 0;
   if (l>tol)
   { 
    x/=l; 
    y/=l; 
    z/=l; 
    return 1;
   } 
   return 0; 
  }
  friend cPoint3d operator -(const cPoint3d& p1, const cPoint3d& p2)
   { return cPoint3d(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z); } 
  friend cPoint3d operator +(const cPoint3d& p1, const cPoint3d& p2)
   { return cPoint3d(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z); } 

  bool Colinear(cPoint3d &b, cPoint3d &c,double tol)
  {
   cPoint3d ab=b-this;
   cPoint3d bc=b-c;
   cPoint3d _cross=ab.Cross(bc);
   return (fabs(_cross.x)<tol&&fabs(_cross.y)<tol&&fabs(_cross.z)<tol);
  }
};

class cBbox3d
{
 public:
  double xmin,xmax,ymin,ymax,zmin,zmax;
  cBbox3d(){xmin=xmax=ymin=ymax=zmin=zmax=0;}
  ~cBbox3d(){}
  void Init(const cPoint3d& p)
           {xmin=xmax=p.x;ymin=ymax=p.y;zmin=zmax=p.z;}
  void Init(const double xminI, const double yminI, const double zminI)
           {xmin=xmax=xminI; ymin=ymax=yminI; zmin=zmax=zminI;}
  void Init(const cBbox3d& b)
           {*this=b;}
  void Update(const cPoint3d& p)
           {if (p.x<xmin) xmin=p.x;if (p.x>xmax) xmax=p.x;
            if (p.y<ymin) ymin=p.y;if (p.y>ymax) ymax=p.y;
            if (p.z<zmin) zmin=p.z;if (p.z>zmax) zmax=p.z;}
  void Update(const cBbox3d& b)
           {if (b.xmin<xmin) xmin=b.xmin;if (b.xmax>xmax) xmax=b.xmax;
            if (b.ymin<ymin) ymin=b.ymin;if (b.ymax>ymax) ymax=b.ymax;
            if (b.zmin<zmin) zmin=b.zmin;if (b.zmax>zmax) zmax=b.zmax;}
  void Update(const double x,const double y,const double z)
           {if (x<xmin) xmin=x;if (x>xmax) xmax=x;
            if (y<ymin) ymin=y;if (y>ymax) ymax=y;
            if (z<zmin) zmin=z;if (z>zmax) zmax=z;}
  void reset()  
           {double maxValue = (double) DBL_MAX;
            xmin=ymin=zmin =   maxValue;
            xmax=ymax=zmax =  -maxValue;}      
  cBbox3d& operator =(const cBbox3d& b)
     {xmin=b.xmin;xmax=b.xmax;
      ymin=b.ymin;ymax=b.ymax;
      zmin=b.zmin;zmax=b.zmax;
      return *this;}
  cBbox3d& operator =(const cPoint3d& p)
     {xmin=xmax=p.x;
      ymin=ymax=p.y;
      zmin=zmax=p.z;
      return *this;}
  
};

#endif
