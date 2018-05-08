/*
** ----------------------------------------------------------------------
**
** msh_bdrsurf.h - Header file for advancing front routine. 
**
** ----------------------------------------------------------------------
**
** Created:      Fev-2008      Antonio C.O. Miranda
**
** ----------------------------------------------------------------------
**
*/

#ifndef _MSH_BOUND_SURF_H_
#define _MSH_BOUND_SURF_H_

#include "mshsurf3d.h"

/*
** -----------------------------------------------------------------
** Public Functions:
*/

int SurfBdryContraction ( int num_org_nodes, 
                          int  num_org_edges,
                          double original_nodes[][2], 
                          int original_edges[][2],
                          int *num_int_nodes, 
                          double **internal_nodes, 
                          int *num_gen_elements, 
                          int **generated_elements,
                          void (*f_surf)(double, double, double *));

/* set a function to consult the size of elements */
void MshSurfBdryRegFunc (MshSurfSizeElement *mshsurf_size);

#endif
