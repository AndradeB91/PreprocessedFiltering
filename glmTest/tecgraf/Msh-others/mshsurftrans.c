
/*
** MshSurfTransf.c
**
**
** This module provides functions to build finite element mesh based
** on transfinite maping. The goal here is to allow the user to specify
** linear or quadratic mesh type with any boundary edge (linear or quadratic).
**
** TecGraf - Grupo de Tecnologia em Computacao Grafica - PUC-Rio
**
** Antonio Miranda - Out 98
**
** Last review: 
*/

#include <string.h>
#include <stdlib.h>			/* memcpy */
#include <stdio.h>

#include "mshsurfmapp.h"
#include "mshsurfmap.h"
// #include "mshmurf.h" 


#ifndef NULL
#define NULL 0L
#endif

/*
** Sort the nodes and connectivity of elements
** Update number of nodes
*/
static void SortNodesAndConn(
int     *nno,     /* number of nodes                 */
int      nel,     /* number of elements              */
double  *points,  /* elements coordinates            */
int     *conn     /* elements connectivity           */
)
{
  int *list1, *list2, index, i;

  list1=(int *)calloc(*nno,sizeof(int));
  list2=(int *)calloc(*nno,sizeof(int));
  for(i=0;i<*nno;i++)
   list1[i]=-1;

  index=0;
  for(i=0;i<nel;i++)
  {
   int j;
   for(j=0;j<conn[index];j++)
   {
    int id=conn[index+j+1];
    list1[id]=1;
   }
   index=index+conn[index]+1;
  }
  
  index=0;
  for(i=0;i<*nno;i++)
  {
   if(list1[i]==1)
   {
    list2[index]=i;
    list1[i]=index;
    index++;
   }
  }
  *nno=index;

  /* Sort nodes */
  for(i=0;i<*nno;i++)
  {
    points[i*3+0]=points[list2[i]*3+0];
    points[i*3+1]=points[list2[i]*3+1];
    points[i*3+2]=points[list2[i]*3+2];
  }
  /* Sort conn */
  index=0;
  for(i=0;i<nel;i++)
  {
   int j;
   for(j=0;j<conn[index];j++)
   {
    int id=conn[index+j+1];
    conn[index+j+1]=list1[id];
   }
   index=index+conn[index]+1;
  }
}
  



/*
** Bilinear adaptative maping
**
** The function computes number of nodes and elments, and alocates and 
** fill the mesh node coordinates and mesh connectivity.
** The function returns 0 if a error occur and 1 on success.
**
** Pararmeters:
** In:
**	bry    -> boundary coordinate vector given in clockwise order
**	m      -> number of nodes in first direction
**	n      -> number of nodes in the other direction
**	type   -> mesh type:
			  3 - T3
			  6 - T6
			  4 - Q4
			  8 - Q8
			  
**	diagtype -> 1-rigth, 2-left, 3-union jack, 4-optimal 
** Out:
**	nno    -> number of nodes
**	nel    -> number of elements
**	pt     -> node coordinates
**	conn   -> elements connectivity
*/
int MshSurfBilinear (double *bry, int m, int n, int elem_type, int diagtype, 
                   int *nno, int *nel, double **pt, int **conn)
{
 /* double *points;
  double *bdry;

  bdry=(double *)bry; */

 
 /* compute number of nodes, elements and incidence */
 *nno=m*n; /* number of nodes */
 if(elem_type == 3 )
 {
    *nel = 2*(m-1)*(n-1);
 }
 else if( elem_type == 4)
 {
    *nel = (m-1)*(n-1);
 }
 else if( elem_type == 6)
 {
    *nel = (m-1)*(n-1)/2;
 }
 else if( elem_type == 8)
 {
    *nel = (m-1)*(n-1)/4;
 }


 /* alocate memory */
 *pt = (double *)calloc ((*nno)*3, sizeof(double));
 if (*pt == NULL)
  return 0;
 *conn = (int *)calloc ((*nel)*(elem_type+1), sizeof(int));
 if (*conn == NULL)
  return 0;
 
 MshSurfTransfin(bry,m,n,elem_type,diagtype,*pt,*conn);

 SortNodesAndConn(nno,*nel,*pt,*conn);

 return 1;
} 


/*
 *  MshSurfTryBilinear
 *************************************************************************
 */
int MshSurfTryBilinear (double *bry, int np, int elem_type, int diagtype, 
                        int *nno, int *nel, double **pt, int **conn)
{
  int m, n;
  
  if (!MshSurfAutomaticBilinearCornes(np, bry, &m, &n))
    return 0;

  if (!MshSurfBilinear (bry, m, n, elem_type, diagtype, nno, nel, pt, conn))
    return 0;

  return 1;
}


/*
 *  MshSurfTryTrilinear
 *************************************************************************
 */
int MshSurfTryTrilinear (double *bry, int np, int elem_type,  
                         int *nno, int *nel, double **pt, int **conn)
{
  int n;

  if (!MshSurfAutomaticTrilinearCornes(np, bry, &n))
    return 0;

  if (!MshSurfTrilinear(bry, n, elem_type, nno, nel, pt, conn))
    return 0;
  
  return 1;
}

/*
** Bilinear mapping meshing with straigth boundary (Lofting)
**
** The function computes number of nodes and elments, and alocates and
** fill the mesh node coordinates and mesh connectivity.
** The function returns 0 if a error occur and 1 on success.
**
** Pararmeters:
** In:
**      bry    -> boundary coordinate vector given in clockwise order
**      m      -> number of nodes in first direction
**      n      -> number of nodes in the other direction
**      dir    -> Direction of lofting (0 => first direction, 1 => other)
**      weight -> Lofting weight
**      type   -> mesh type:
                          3 - T3
                          6 - T6
                          4 - Q4
                          8 - Q8

**      diagtype -> 1-rigth, 2-left, 3-union jack, 4-optimal
** Out:
**      nno    -> number of nodes
**      nel    -> number of elements
**      pt     -> node coordinates
**      conn   -> elements connectivity
*/
int MshSurfLoft (double *bry, int m, int n, int dir, double weight,
               int elem_type, int diagtype,
               int *nno, int *nel, double **pt, int **conn)
{
 /* double *points;
  double *bdry;

  bdry=(double *)bry; */


 /* compute number of nodes, elements and incidence */
 *nno=m*n; /* number of nodes */
 if(elem_type == 3 )
 {
    *nel = 2*(m-1)*(n-1);
 }
 else if( elem_type == 4)
 {
    *nel = (m-1)*(n-1);
 }
 else if( elem_type == 6)
 {
    *nel = (m-1)*(n-1)/2;
 }
 else if( elem_type == 8)
 {
    *nel = (m-1)*(n-1)/4;
 }


 /* alocate memory */
 *pt = (double *)calloc ((*nno)*3, sizeof(double));
 if (*pt == NULL)
  return 0;
 *conn = (int *)calloc ((*nel)*(elem_type+1), sizeof(int));
 if (*conn == NULL)
  return 0;

 MshSurfLofting(bry,m,n,dir,weight,elem_type,diagtype,*pt,*conn);

 SortNodesAndConn(nno,*nel,*pt,*conn);

 return 1;
}



/*
** Trilinear adaptative maping
**
** The function computes number of nodes and elments, and alocates and 
** fill the mesh node coordinates and mesh connectivity.
** The function returns 0 if a error occur and 1 on success.
**
** Pararmeters:
** In:
**	bry    -> boundary coordinate vector given in clockwise order
**	m      -> number of nodes in one direction
**	type   -> mesh type: 0 - linear or 1 - quadratic
**	cod[3] -> boundaries type: 0 - linear or 1 - quadratic
** Out:
**	nno    -> number of nodes
**	nel    -> number of elements
**	pt     -> node coordinates
**	conn   -> elements connectivity
*/
int MshSurfTrilinear (double *bry, int m, int elem_type,  
                  int *nno, int *nel, double **pt, int **conn)
{
  int     i;

  
  /* compute number of nodes, elements and incidence */
 *nno=m*m; /* number of nodes */
 if(elem_type == 3 )
 {
    *nel = (m-1)*(m-1);
 }
 else if( elem_type == 4)
 {
	 *nel=0;
	 for(i=(m-1); i>=0; i--)
		 *nel=*nel+i;
 }
 else if( elem_type == 6)
 {
    *nel = (m-1)*(m-1)/4;
 }
 else if( elem_type == 8)
 {
	 *nel=0;
	 for(i=((m-1)/2); i>=0; i--)
		 *nel=*nel+i;
 }

 /* alocate memory */
 *pt = (double *)calloc ((*nno)*3, sizeof(double));
 if (pt == NULL)
  return 0;
 *conn = (int *)calloc ((*nel)*(elem_type+1), sizeof(int));
 if (*conn == NULL)
  return 0;
 

 MshSurfTrimap(bry,m,elem_type,*pt,*conn);

 SortNodesAndConn(nno,*nel,*pt,*conn);

 return 1;
} 


/*
** Bilinear adaptative maping with collapsed node
**
** The function computes number of nodes and elments, and alocates and 
** fill the mesh node coordinates and mesh connectivity.
** The function returns 0 if a error occur and 1 on success.
**
** Pararmeters:
** In:
**	bry    -> boundary coordinate vector given in clockwise order
**	m      -> number of nodes in the side opposite to the collapsed one
**	n      -> number of nodes in the other direction
**	type   -> mesh type:
			  3 - T3
			  6 - T6
			  4 - Q4
			  8 - Q8
			  
**	diagtype -> 1-rigth, 2-left, 3-union jack, 4-optimal 
** Out:
**	nno    -> number of nodes
**	nel    -> number of elements
**	pt     -> node coordinates
**	conn   -> elements connectivity
*/

int MshSurfCollBilinear (double *bry, int m, int n, int elem_type, int diagtype,  
                        int *nno, int *nel, double **pt, int **conn)
{
  int     nsegs = 0; /* number of segments int the opposite to the collapsed node */ 
  int     blank = 0; /* number of repeated nodes in the collapsed node */

  int     i,index;
  

  /* compute number of nodes, elements and incidence */
 *nno=m*n; /* number of nodes */
 if(elem_type == 3 )
 {
    *nel = 2*(m-1)*(n-2)+(m-1);
	nsegs = m-1;
	blank=nsegs;
 }
 else if( elem_type == 4)
 {
	 *nel = (m-1)*(n-1);
	 nsegs = m-1;
	 blank=nsegs;
 }
 else if( elem_type == 6)
 {
    *nel = (m-1)*(n-3)/2+(m-1)/2;
	nsegs = (m-1)/2;
	blank=nsegs*2;
 }
 else if( elem_type == 8)
 {
	 *nel = (m-1)*(n-1)/4;
	 nsegs = (m-1)/2;
	 blank=nsegs*2;
 }

 
 /* alocate memory */
 *pt = (double *)calloc ((*nno)*3, sizeof(double));
 if (*pt == NULL)
  return 0;
 *conn = (int *)calloc ((*nel)*(elem_type+1), sizeof(int));
 if (*conn == NULL)
  return 0;

   MshSurfTrsfncoll(bry,m,n,elem_type,diagtype,*pt,*conn); 
 
  /* Corrections on the Mesh */
  /* repeated nodes in quadratic edges*/
  index=0;
  if((elem_type==6)||(elem_type==8))
  {
   int i,j,cont=1;
   for(i=0;i<*nel;i++)
   {
    for(j=0;j<(*conn)[index];j++)
    {
     int id=(*conn)[index+j+1];
     if((nsegs*2+1+(cont*2))==id)
     {
      (*conn)[index+j+1]=nsegs*2+cont*2;
      cont++;
     }
    }
    index=index+(*conn)[index]+1;
    if(cont==nsegs) break;
   }
  }
  /* blank nodes in the collapsed node */
  *nno=*nno-blank;
  for(i=1;i<*nno;i++)
  {
	  (*pt)[i*3+0]=(*pt)[(i+blank)*3+0];
	  (*pt)[i*3+1]=(*pt)[(i+blank)*3+1];
	  (*pt)[i*3+2]=(*pt)[(i+blank)*3+2];  }
  /* Sort conn */
  index=0;
  for(i=0;i<*nel;i++)
  {
   int j;
   for(j=0;j<(*conn)[index];j++)
   {
	   int id=(*conn)[index+j+1];
	   if(id<=blank)
		   (*conn)[index+j+1]=0;
	   else
		   (*conn)[index+j+1]=id-blank;
   }
   index=index+(*conn)[index]+1;
  }

  SortNodesAndConn(nno,*nel,*pt,*conn);
 
 return 1;
} 



/*
** Bilinear collapsed mapping meshing with straigth
** boundary shape (Lofting)
**
** The function computes number of nodes and elments, and alocates and
** fill the mesh node coordinates and mesh connectivity.
** The function returns 0 if a error occur and 1 on success.
**
** Pararmeters:
** In:
**      bry    -> boundary coordinate vector given in clockwise order
**      m      -> number of nodes in the side opposite to the collapsed one
**      n      -> number of nodes in the other direction
**      weight -> Lofting weight
**      type   -> mesh type:
                          3 - T3
                          6 - T6
                          4 - Q4
                          8 - Q8

**      diagtype -> 1-rigth, 2-left, 3-union jack, 4-optimal
** Out:
**      nno    -> number of nodes
**      nel    -> number of elements
**      pt     -> node coordinates
**      conn   -> elements connectivity
*/

int MshSurfCollLoft (double *bry, int m, int n, double weight,
                   int elem_type, int diagtype,
                   int *nno, int *nel, double **pt, int **conn)
{
  int     nsegs = 0; /* number of segments int the opposite to the collapsed node */
  int     blank = 0; /* number of repeated nodes in the collapsed node */

  int     i,index;

  /* compute number of nodes, elements and incidence */
 *nno=m*n; /* number of nodes */
 if(elem_type == 3 )
 {
    *nel = 2*(m-1)*(n-2)+(m-1);
        nsegs = m-1;
        blank=nsegs;
 }
 else if( elem_type == 4)
 {
         *nel = (m-1)*(n-1);
         nsegs = m-1;
         blank=nsegs;
 }
 else if( elem_type == 6)
 {
    *nel = (m-1)*(n-3)/2+(m-1)/2;
        nsegs = (m-1)/2;
        blank=nsegs*2;
 }
 else if( elem_type == 8)
 {
         *nel = (m-1)*(n-1)/4;
         nsegs = (m-1)/2;
         blank=nsegs*2;
 }
 /* alocate memory */
 *pt = (double *)calloc ((*nno)*3, sizeof(double));
 if (*pt == NULL)
  return 0;
 *conn = (int *)calloc ((*nel)*(elem_type+1), sizeof(int));
 if (*conn == NULL)
  return 0;

  MshSurfLoftcoll(bry,m,n,weight,elem_type,diagtype,*pt,*conn); 

  /* Corrections on the Mesh */
  /* repeated nodes in quadratic edges*/
  index=0;
  if((elem_type==6)||(elem_type==8))
  {
   int i,j,cont=1;
   for(i=0;i<*nel;i++)
   {
    for(j=0;j<(*conn)[index];j++)
    {
     int id=(*conn)[index+j+1];
     if((nsegs*2+1+(cont*2))==id)
     {
      (*conn)[index+j+1]=nsegs*2+cont*2;
      cont++;
     }
    }
    index=index+(*conn)[index]+1;
    if(cont==nsegs) break;
   }
  }
  /* blank nodes in the collapsed node */
  *nno=*nno-blank;
  for(i=1;i<*nno;i++)
  {
          (*pt)[i*3+0]=(*pt)[(i+blank)*3+0];
          (*pt)[i*3+1]=(*pt)[(i+blank)*3+1];
          (*pt)[i*3+2]=(*pt)[(i+blank)*3+2];
  }
  /* Sort conn */
  index=0;
  for(i=0;i<*nel;i++)
  {
   int j;
   for(j=0;j<(*conn)[index];j++)
   {
           int id=(*conn)[index+j+1];
           if(id<=blank)
                   (*conn)[index+j+1]=0;
           else
                   (*conn)[index+j+1]=id-blank;
   }
   index=index+(*conn)[index]+1;
  }

  SortNodesAndConn(nno,*nel,*pt,*conn);

 return 1;
}

