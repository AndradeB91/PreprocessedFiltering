/*
** ---------------------------------------------------------------------------
**
** sll.h:   Definitions for singly linked lists
**
** TeCGraf - Grupo de Tecnologia em Computacao Grafica, PUC-Rio
**
** ---------------------------------------------------------------------------
*/

#ifndef _SLL_H
#define _SLL_H 

#include "deflib2d.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef MSH2D_API struct _sllrec {
  struct _sllrec *next;
  int            id;
} *Sll, Sllrec;

/*
** ---------------------------------------------------------------------------
** Public Functions:
*/

MSH2D_API int  SllAddTop( Sll *, unsigned int );
MSH2D_API int  SllAddEnd( Sll *, unsigned int, Sll * );
MSH2D_API int  SllAddBef( Sll *, unsigned int, Sll, Sll * );
MSH2D_API int  SllAddOrd( Sll *, unsigned int, int, Sll * );
MSH2D_API int  SllDelete( Sll *, Sll );
MSH2D_API void SllDelAll( Sll * );
MSH2D_API void SllRenum ( Sll * );
MSH2D_API int  SllNumElm( Sll * );
MSH2D_API int  SllGetElm( Sll *, int, Sll * );

#ifdef __cplusplus
}
#endif

#endif

