/*
* amr3bind.c
* Binding [C] for amr3tree.c [C++] implementation
*/

#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "amr3tree.h"
#include "amr3bind.h"

void *RtreeCreate(void)
{
 amR3tree *r3;
 r3=new amR3tree(4,10);
 r3->Create();
 return (void*)r3;
}

void RtreeDestroy(void *r)
{
 amR3tree *r3=(amR3tree*)r;
 r3->Delete();
 delete r3;
}

void RtreeInsert(void *r,
                 void *info,
                 double xmin,double xmax,
                 double ymin,double ymax,
                 double zmin,double zmax)
{
 amR3tree *r3=(amR3tree*)r;
 amBoxId bid;
 bid.P=info;
 bid.Xmin=xmin, bid.Xmax=xmax,
 bid.Ymin=ymin, bid.Ymax=ymax,
 bid.Zmin=zmin, bid.Zmax=zmax;
 r3->Insert(&bid);
}

void RtreeDelete(void *r,
                 void *info,
                 double xmin,double xmax,
                 double ymin,double ymax,
                 double zmin,double zmax)
{
 amR3tree *r3=(amR3tree*)r;
 amBoxId bid;
 bid.P=info;
 bid.Xmin=xmin, bid.Xmax=xmax,
 bid.Ymin=ymin, bid.Ymax=ymax,
 bid.Zmin=zmin, bid.Zmax=zmax;
 r3->Remove(&bid);
}

void RtreeInitTraverse(void *r)
{
 amR3tree *r3=(amR3tree*)r;
 r3->SearchAll();
}

void *RtreeTraverse(void *r,
                    double *xmin,double *xmax,
                    double *ymin,double *ymax,
                    double *zmin,double *zmax)
{
 amR3tree *r3=(amR3tree*)r;
 amBoxId bid;
 if (!r3->ResultAll(&bid))
 {
  *xmin=bid.Xmin,
  *xmax=bid.Xmax,
  *ymin=bid.Ymin,
  *ymax=bid.Ymax,
  *zmin=bid.Zmin,
  *zmax=bid.Zmax;
  return bid.P;
 }
 return NULL;
}

void RtreeInitSearchBox(void *r,
                        double xmin,double xmax,
                        double ymin,double ymax,
                        double zmin,double zmax)
{
 amR3tree *r3=(amR3tree*)r;
 amBox b;
 b.Xmin=xmin, b.Ymin=ymin, b.Zmin=zmin,
 b.Xmax=xmax, b.Ymax=ymax, b.Zmax=zmax;
 r3->Search(&b);
}

void *RtreeSearchBox(void *r,
                     double *xmin,double *xmax,
                     double *ymin,double *ymax,
                     double *zmin,double *zmax)
{
 amR3tree *r3=(amR3tree*)r;
 amBoxId bid;
 if(!r3->Result(&bid))
 {
  *xmin=bid.Xmin,
  *xmax=bid.Xmax,
  *ymin=bid.Ymin,
  *ymax=bid.Ymax,
  *zmin=bid.Zmin,
  *zmax=bid.Zmax;
  return bid.P;
 }
 return NULL;
}

int RtreeNumber( void *r )
{
 amR3tree *r3=(amR3tree*)r;

 return r3->Number();
}
