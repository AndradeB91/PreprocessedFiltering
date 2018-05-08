#ifndef _MSHMAPSURF_H
#define _MSHMAPSURF_H

void MshSurfTransfin( double *bdynodes, int nu, int nv, int elemtype, 
                    int diagtype, double *gennodes, int *genelems );
void MshSurfLofting( double *bdynodes, int nu, int nv, int dir, double weight,
                  int elemtype, int diagtype, double *gennodes, int *genelems );
void MshSurfTrsfncoll( double *bdynodes, int nu, int nv, int elemtype, 
                    int diagtype, double *gennodes, int *genelems );
void MshSurfLoftcoll( double *bdynodes, int nu, int nv, double weight,
                  int elemtype, int diagtype, double *gennodes, int *genelems );
void MshSurfTrimap( double *bdynodes, int n, int elemtype, 
                  double *gennodes, int *genelems );

#endif
