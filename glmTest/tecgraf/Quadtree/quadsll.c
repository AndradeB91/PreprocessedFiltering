/*
** ---------------------------------------------------------------------------
**
** sll.c - Package to manage singly linked lists.
**
**	remarks:
** ---------------------------------------------------------------------------
**
**	(0) Functionality:
**	This is general package for generic handling of singly linked
**	lists.  It is generic in the sense that it performs the basic
**	operations on singly linked lists (add at the top, add at the
**	end, etc.).  To do this it is required that the first field
**	of a list element be a pointer to the next element in the list.
**	As usual, this field is null for the last element of the list.
**	The allocation of memory for a list element is done in this
**	package.  For this, the insertion operators require the size
**	of the list element (including the next field) to be created
**	as a parameter.
**
**	(1) Id mark:
**	The functionality of this package was increased so that it also
**	numbers the list elements.  Therefore, one field was added to the
**	generic structure of a singly list element: an id field.  This
**	generic structure is shown below:
**
**		typedef struct _sllrec {
**			struct _sllrec	*next;
**			int		id;
**		   } *Sll, Sllrec;
**
**
** ---------------------------------------------------------------------------
**
**	entry pts:
**
** ---------------------------------------------------------------------------
**
** 	int    SllAddTop( Sll *head, unsigned int size )
**
**		    head	- head of list 			(in/out)
**		    size	- size of element of list	(in)
**
**		This function inserts an element at the top of a list.
**		The new element created is returned as the new head.
**
** ---------------------------------------------------------------------------
**
**	int	SllAddEnd( Sll *head, unsigned int size, Sll *last )
**
**		    head	- head of list 			(in/out)
**		    size	- size of element of list	(in)
**		    last	- last element created		(out)
**
**		This function inserts an element at the end of a list.
**		If the list is empty (null head) the new element created
**		is also returned as the new head.
**
** ---------------------------------------------------------------------------
**
** 	int    SllAddBef( Sll *head, unsigned int size, Sll after, Sll *elem )
**
**		    head	- head of list 			(in/out)
**		    size	- size of element of list	(in)
**		    after	- elem after the new one	(in)
**		    elem	- new element created		(out)
**
**		This function inserts an element before a given element
**		of a list.
**		If the list is empty (null head) or if the given after
**		element was not found in the list, the new element is
**		put at the end of the list and the after parameter is
**		disregarded.
**		If the list is empty (null head), the returned head is
**		the same as the new element
**
** ---------------------------------------------------------------------------
**
** 	int    SllAddOrd( Sll *head, unsigned int size, int id, Sll *elem )
**
**		    head	- head of list 			(in/out)
**		    size	- size of element of list	(in)
**		    id		- id of new element created	(in)
**		    elem	- new element created		(out)
**
**		This function inserts an element in order of id's in
**		a list.
**		If already exists an element with the given id, it inserts
**		the new element just before the element which had this id.
**		If not, the element is inserted at the end of the list.
**		If the list is empty (null head), the returned head is
**		the same as the new element
**
** ---------------------------------------------------------------------------
**
**	int	SllDelete( Sll *head, Sll elem )
**
**		    head	- head of list 			(in/out)
**		    elem	- element to delete		(in)
**
**		This function deletes the element of a list.
**		Delete here means release the memory used by the element.
**		If the list is empty (null head) it returns a false (0)
**		status.  If the given element does not belong to the list,
**		also returns a false (0) status.  In any other case,
**		it returns a true (1) status.
**
** ---------------------------------------------------------------------------
**
**	void	SllDelAll( Sll *head )
**
**		    head	- head of list 			(in/out)
**
**		This function deletes all the elements of a list.
**		The return value of the head of the list is null.
**
** ---------------------------------------------------------------------------
**
**	void	SllRenum( Sll *head )
**
**		    head	- head of list 			(in)
**
**		This function renumbers consecutively, starting from
**		one (1), the elements of a list.
**
** ---------------------------------------------------------------------------
**
**	int	SllNumElm( Sll *head )
**
**		    head	- head of list 			(in)
**
**		This function returns as the function value the
**		total number of elements in a list.
**
** ---------------------------------------------------------------------------
**
**	int	SllGetElm( Sll *head, int id, Sll *elem )
**
**		    head	- head of list 			(in)
**		    id		- id of target element		(in)
**		    elem	- retrieved element		(out)
**
**		This function returns an element of a list which has
**		a given id.  If the list is empty (null head) or
**		no element has the given id, the function returns a
**		false (0) flag.  Otherwise, it returns a true (1) flag.
**
** ---------------------------------------------------------------------------
**
**	Version:	0-001
**
**	History:	Created    28-Jan-92
**					Eduardo Corseuil
**					Adriane Cavalieri
**					Vinicius Samu
**
**      TeCGraf - Grupo de Tecnologia em Computacao Grafica, PUC-Rio
**
** ---------------------------------------------------------------------------
*/

/*
** ---------------------------------------------------------------------------
**	Global variables and symbols
*/
#include <stdio.h>
#include <stdlib.h>
#include "quadsll.h"

/*
** ---------------------------------------------------------------------------
**	Local variables and symbols: (none)
*/

/*
** ---------------------------------------------------------------------------
**	Privite functions: (none)
*/

/*
** ---------------------------------------------------------------------------
**	Entry points begin here:
*/

/* ========================  SllAddTop  =================================== */

int SllAddTop( Sll *head, unsigned int size )
{
 Sll nw;

 nw = (Sll)calloc( 1, size );

 if( nw == NULL ) return( 0 );

 nw->next = *head;
 nw->id   = 0;

 *head = nw;

 return( 1 );
}

/* ========================  SllAddEnd  =================================== */

int SllAddEnd( Sll *head, unsigned int size, Sll *last )
{
 Sll nw;
 Sll curr;

 nw = (Sll)calloc( 1, size );

 if( nw == NULL ) return( 0 );

 nw->next = NULL;
 nw->id   = 0;

 if( *head == NULL )
 {
  *head = nw;
 }
 else
 {
  curr = *head;
  while( curr->next != NULL )
  {
   curr = curr->next;
  }
  curr->next = nw;
 }

 *last = nw;

 return( 1 );
}

/* ========================  SllAddBef  =================================== */

int SllAddBef( Sll *head, unsigned int size, Sll after, Sll *elem )
{
 Sll nw;
 Sll curr;

 if( *head == NULL) return( SllAddEnd( head, size, elem ) );

 nw = (Sll)calloc( 1, size );

 if( nw == NULL ) return( 0 );

 curr = *head;
 while( (curr->next != NULL) && (curr->next != after) )
 {
  curr = curr->next;
 }

 nw->next = curr->next;
 nw->id   = 0;

 curr->next = nw;

 *elem = nw;

 return( 1 );
}

/* ========================  SllAddOrd  =================================== */

int SllAddOrd( Sll *head, unsigned int size, int id, Sll *elem )
{
 Sll nw;
 Sll curr;

 if( *head == NULL) return( SllAddEnd( head, size, elem ) );

 nw = (Sll)calloc( 1, size );

 if( nw == NULL ) return( 0 );

 curr = *head;
 while( (curr->next != NULL) && (curr->id < (id-1)) )
 {
  curr = curr->next;
 }

 nw->next = curr->next;
 nw->id   = id;

 curr->next = nw;

 *elem = nw;

 return( 1 );
}

/* ========================  SllDelete  =================================== */

int SllDelete( Sll *head, Sll elem )
{
 Sll prev;

 if( *head == NULL ) return( 0 );

 if( *head == elem )
 {
  *head = elem->next;
 }
 else
 {
  prev = *head;

  while( ( prev->next != elem ) && ( prev->next != NULL ) )
  {
   prev = prev->next;
  }

  if( prev->next == NULL ) return( 0 );

  prev->next = elem->next;
 }

 free ( elem );

 return( 1 );
}

/* ========================  SllDelAll  =================================== */

void SllDelAll( Sll *head )
{
 while( *head != NULL )
 {
  SllDelete( head, *head );
 }
 return;
}

/* ========================  SllRenum  ==================================== */

void SllRenum( Sll *head )
{
 Sll elem;
 int cur = 0;			/* current id value		*/

 for( elem = *head; elem != NULL; elem = elem->next )
 {
  elem->id = ++cur;
 }
 return;
}

/* ========================  SllNumElm  =================================== */

int SllNumElm( Sll *head )
{
 Sll elem;
 int total = 0;

 for( elem = *head; elem != NULL; elem = elem->next )
 {
  total++;
 }
 return( total );
}

/* ========================  SllGetElm  =================================== */

int SllGetElm( Sll *head, int id, Sll *elem )
{
 Sll cur_elem;

 *elem = NULL;

 for( cur_elem = *head; cur_elem != NULL; cur_elem = cur_elem->next )
 {
  if( cur_elem->id == id )
  {
   *elem = cur_elem;
   break;
  }
 }

 if( *elem == NULL ) return( 0 );
 else                return( 1 );
}

