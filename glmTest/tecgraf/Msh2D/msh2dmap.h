#ifndef _MSHMAP2D_H
#define _MSHMAP2D_H

void Msh2DTransfin( double *bdynodes, int nu, int nv, int elemtype, 
                    int diagtype, double *gennodes, int *genelems );
void Msh2DLofting( double *bdynodes, int nu, int nv, int dir, double weight,
                  int elemtype, int diagtype, double *gennodes, int *genelems );
void Msh2DTrsfncoll( double *bdynodes, int nu, int nv, int elemtype, 
                    int diagtype, double *gennodes, int *genelems );
void Msh2DLoftcoll( double *bdynodes, int nu, int nv, double weight,
                  int elemtype, int diagtype, double *gennodes, int *genelems );
void Msh2DTrimap( double *bdynodes, int n, int elemtype, 
                  double *gennodes, int *genelems );

#endif
