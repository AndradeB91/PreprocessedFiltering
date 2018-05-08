
/*
** msh2d.c
*/
#include <stdlib.h>
#include "msh2d.h" 

#ifndef NULL
#define NULL 0L
#endif

/* ref. factor */
static double ref_factor = 1.0;

/*
**
*/
void   Msh2DSetRefFactor (double factor)
{
  if (factor < 0.1)
    ref_factor = 0.1;
  else
    ref_factor = factor;
}

/*
**
*/
double Msh2DGetRefFactor (void)
{
   return (ref_factor); 
}


/*
** Free memory of allocate points into Msh2D
*/
void Msh2DFreeNodes(double *points)
{
	free(points);
	points=NULL;
}

/*
** Free memory of allocate connectivity into Msh2D
*/
void Msh2DFreeConn(int *Conn)
{
	free(Conn);
	Conn=NULL;
}
