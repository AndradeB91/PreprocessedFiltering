#include <stdio.h>

#include "seamLib.h"
#include "../Surf3D/cPt3d.h"
#include "msh2d.h"

typedef struct point	point_3d;
struct point
{
 double x,y,z;
};

float p_vetorial(float Ax,float Ay,float Bx,float By)
{
  return(Ax*By - Ay*Bx);
}

// when a face with collinear edge is found
// the check connectivity routine skip this face to the next one
int ColinearEdge(int pos,int *con, double *pts)
{
 if (con != NULL)
 {
  if (con[pos]!=4) return -1;
  int id0=con[pos+1];
  int id1=con[pos+2];
  int id2=con[pos+3];
  int id3=con[pos+4];

  if (id0!=id1 && id1!=id2 && id2!=id3 && id3!=id0)
  {
   cPoint3d v01(pts[id1*2+0]-pts[id0*2+0],pts[id1*2+1]-pts[id0*2+1],0);
   cPoint3d v12(pts[id2*2+0]-pts[id1*2+0],pts[id2*2+1]-pts[id1*2+1],0);
   v01.Normalize();
   v12.Normalize();
   cPoint3d vt=v12.Cross(v01);
   if (!vt.Normalize(0.1)) return 1;

   cPoint3d v23(pts[id3*2+0]-pts[id2*2+0],pts[id3*2+1]-pts[id2*2+1],0);
            v23.Normalize();
            vt=v23.Cross(v12);
   if (!vt.Normalize(0.1)) return 2;

   cPoint3d v30(pts[id0*2+0]-pts[id3*2+0],pts[id0*2+1]-pts[id3*2+1],0);
            v30.Normalize();
            vt=v30.Cross(v23);
   if (!vt.Normalize(0.1)) return 3;

            vt=v01.Cross(v30);
   if (!vt.Normalize(0.1)) return 0;
  }
  return -1;
 }
 return -1;
}

int CheckConnectivity(int ne, int *con, double *pts)
{
 if (con != NULL)
 {
  int pos=0;
  if (ne == 1) return 1;
  for (int i=0;i<ne;++i,pos+=con[pos]+1)
  {
   if (ColinearEdge(pos,con,pts)> 0)
    continue;
   for (int j=0;j<con[pos];++j)
   {
    int a=con[pos+j+1];
    int b= con[pos+((j+1)%con[pos])+1];
    // reta é um ponto
    if (a == b)
     continue;

    for (int k=0; k <con[pos]; k++)
    {
     if (k == j)
      continue;
     int c= con[pos+k+1];
     int d= con[pos+((k+1)%con[pos])+1];
     // reta é um ponto
     if (c == d)
      continue;
     double z1,z2;
     point_3d p1,p2,p3,p4;
     p1.x = pts[a*2+0];
     p1.y = pts[a*2+1];
     p1.z = 0;
     p2.x = pts[b*2+0];
     p2.y = pts[b*2+1];
     p2.z = 0;
     p3.x = pts[c*2+0];
     p3.y = pts[c*2+1];
     p3.z = 0;
     p4.x = pts[d*2+0];
     p4.y = pts[d*2+1];
     p4.z = 0;
     z1 = p_vetorial((p2.x-p1.x),(p2.y-p1.y),(p3.x-p1.x),(p3.y-p1.y));
     z2 = p_vetorial((p2.x-p1.x),(p2.y-p1.y),(p4.x-p1.x),(p4.y-p1.y));
     //if ((z1 < 0 && z1 < 1) || (z1 > 0 && z2 < 1))
     if (z1 > 0 && z2 > 0)
      return 0;
    }
   }//end for
  }// end for all ne
  return 1;
 }
 return 0;
}

int MSH2DQuadSeam (int n_loops, int *loop_segs, double *bdry_pts, int type_mesh,
                   int *n_node, double **coords, int *n_elem, int **conn)
{
 if (type_mesh == 4)
 {
  // rtree containing all edges of the mesh
  int ne = *n_elem;
  int np = *n_node;
  int *con = *conn;
  double *pts = *coords;

  double i;

  for (i=1; i >=0.1; i=i-0.1)
  {
   Msh2DSetRefFactor (i);
   printf("Mesh factor %f\n",i);
   if (Msh2DQuadSeam(n_loops,loop_segs, bdry_pts,type_mesh,&np,&pts,&ne,&con))
   {
    if (!CheckConnectivity(ne,con,pts))
    {
     ne=0;
     np=0;
     con = NULL;
     pts = NULL;
     //printf("con errado!!!\n");
     continue;
    }
    else
    {
     //printf("Mesh factor %f ok\n",i);
     *n_elem = ne;
     *n_node = np;
     *conn = con;
     *coords = pts;
     return i;
   }
   }
   //printf("mudando fator da malha!!!\n");
   continue;
  }

  for (i=1.1; i <=15; i=i+0.1)
  {
   Msh2DSetRefFactor (i);
   printf("Mesh factor %f\n",i);
   if (Msh2DQuadSeam(n_loops,loop_segs, bdry_pts,type_mesh,&np,&pts,&ne,&con))
   {
    if (!CheckConnectivity(ne,con,pts))
    {
     ne=0;
     np=0;
     con = NULL;
     pts = NULL;
     //printf("con errado!!!\n");
     continue;
    }
    else
    {
     //printf("Mesh factor %f ok\n",i);
     *n_elem = ne;
     *n_node = np;
     *conn = con;
     *coords = pts;
     return i;
   }
   }
   //printf("mudando fator da malha!!!\n");
   continue;
  }
 }
 return -1;
}

int MSH2DBilinear (double *bry, int m, int n, int /*elem_type*/, int /*diagtype*/,
                   int *nno, int *nel, double **pt, int **conn)
 {
  // rtree containing all edges of the mesh
  int ne = *nel;
  int np = *nno;
  int *con = *conn;
  double *pts = *pt;
  double i;

  for (i=1.0; i >=0.1; i=i-0.1)
  {
   Msh2DSetRefFactor (i);
   printf("Mesh factor %f\n",i);
   if (Msh2DBilinear(bry,m,n,4,4,&np, &ne, &pts, &con))
   {
    int ne = *nel;
    int np = *nno;
    int *con = *conn;
    double *pts = *pt;
    if (!CheckConnectivity(ne,con,pts))
    {
     ne=0;
     np=0;
     con = NULL;
     pts = NULL;
     printf("con errado!!!\n");
     continue;
    }
    else
    {
     //printf("Mesh factor %f ok\n",i);
     *nel = ne;
     *nno = np;
     *conn = con;
     *pt = pts;
     return 1;
    }
   }
   //printf("mudando fator da malha!!!\n");
   continue;
  }

  for (i=1.1; i <=15; i=i+0.1)
  {
   Msh2DSetRefFactor (i);
   printf("Mesh factor %f\n",i);
   if (Msh2DBilinear(bry,m,n,4,4,&np, &ne, &pts, &con))
   {
    if (!CheckConnectivity(ne,con,pts))
    {
     ne=0;
     np=0;
     con = NULL;
     pts = NULL;
     //printf("con errado!!!\n");
     continue;
    }
    else
    {
     //printf("Mesh factor %f ok\n",i);
     *nel = ne;
     *nno = np;
     *conn = con;
     *pt = pts;
     return 1;
    }
   }
   //printf("mudando fator da malha!!!\n");
   continue;
  }
 return 0;
}
