//
// CArbQueue Template Class header file
//
// Description -
//   This class implements a FIFO stack.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.5 $  $Date: 2001/06/15 13:30:03 $  $Author: wash $
//

#ifndef ArbQueue_hh
#define ArbQueue_hh

#define SMALL_HASH_TABLE 4

template<class EntryType>
class CArbQueue {
    public:

        CArbQueue() ;
        ~CArbQueue() ;

        void AddToBack(const EntryType &value) ;
        bool RemoveFromFront(EntryType *value) ;

        int NumEntries() { return(NumQueueEntries) ; } ;

    private:

        struct ArbQueueEntry {
            EntryType    entry ;
            ArbQueueEntry *prev ;
            ArbQueueEntry *next ;
        } ;

        ArbQueueEntry *head ;
        ArbQueueEntry *tail ;
        int NumQueueEntries ;
} ;

/*
TEMPLATE CLASS CArbQueue<class EntryType>

  This object implements a first in first out queue template class. A 
  copy is made of items placed in the queue 

  Template Arguments:

    EntryType - type for the queueu entries


PUBLIC INTERFACE

  Public Member Functions:

    CArbQueue - constructor 

      CArbQueue()

      Description: This is the constructor for an ArbQueue object. 


    CArbQueue - destructor 

      ~CArbQueue()

      Description: This is the destructor for an ArbQueue object. 


    AddToBack - add an entry 

      void AddToBack(const EntryType &value)

        value - (in)  entry data 

      Description: This function adds an entry to the queue. 


    RemoveFromFront - remove and item from the queue 

      bool RemoveFromFront(EntryType *value)

        value - (out) place to copy the entry data 

      Description: This function removes an element from the "front" 
          of the queue. The queue entry is copied to the memory 
          pointed to by the argument 

      Return Value: false if the queue is empty, true otherwise 


    NumEntries - number of queue entries 

      int NumEntries()

      Description: This routine returns the number of entries 
          currently in the queue. 

      Return Value: number of queue entries 


PRIVATE INTERFACE

  Private Data Structures:

    struct ArbQueueEntry

      data structure for a queue entry 

      Member Variables:

        EntryType entry - entry data 

        ArbQueueEntry *prev - previous doubly linked list 
            pointer 

        ArbQueueEntry *next - next doubly linked list pointer 


  Private Member Variables:

    ArbQueueEntry *head - pointer to the head entry 

    ArbQueueEntry *tail - pointer to the tail entry 

    int NumQueueEntries - current number of entries in the queue 

*/

// %(CArbQueue::CArbQueue-constructor-|) 
/* ++ ----------------------------------------------------------
**
**    CArbQueue - constructor 
**
**      CArbQueue()
**
**      Description: This is the constructor for an ArbQueue object. 
**
**
** -- */

template<class EntryType> CArbQueue<EntryType>::CArbQueue() :
    head(0), tail(0), NumQueueEntries(0)
{ } ;




// %(CArbQueue::CArbQueue-destructor-|~)
/* ++ ----------------------------------------------------------
**
**    CArbQueue - destructor 
**
**      ~CArbQueue()
**
**      Description: This is the destructor for an ArbQueue object. 
**
**
** -- */

template<class EntryType> CArbQueue<EntryType>::~CArbQueue()
{
    if (NumQueueEntries > 0) {
        ArbQueueEntry *nextPtr, *hPtr = head ;
        while (hPtr != 0) {
            nextPtr = hPtr->next ;
	    delete hPtr ;
	    hPtr = nextPtr ;
	}
    }
}




// %(CArbQueue::AddToBack-void-|-EntryType-const|&) 
/* ++ ----------------------------------------------------------
**
**    AddToBack - add an entry 
**
**      void AddToBack(const EntryType &value)
**
**        value - (in)  entry data 
**
**      Description: This function adds an entry to the queue. 
**
**
** -- */

template<class EntryType>
void CArbQueue<EntryType>::AddToBack(const EntryType &entry)
{
    ArbQueueEntry *hPtr = new ArbQueueEntry ;
    hPtr->entry = entry ;

    if (NumQueueEntries == 0) {
        hPtr->prev = hPtr->next = 0 ;
        head = tail = hPtr ;
    } else {
        hPtr->prev = tail ;
        hPtr->next = 0 ;
        tail->next = hPtr ;
        tail = hPtr ;
    }
    NumQueueEntries++ ;
}




// %(CArbQueue::RemoveFromFront-bool-|-EntryType-|*) 
/* ++ ----------------------------------------------------------
**
**    RemoveFromFront - remove and item from the queue 
**
**      bool RemoveFromFront(EntryType *value)
**
**        value - (out) place to copy the entry data 
**
**      Description: This function removes an element from the "front" 
**          of the queue. The queue entry is copied to the memory 
**          pointed to by the argument 
**
**      Return Value: false if the queue is empty, true otherwise 
**
**
** -- */

template<class EntryType>
bool CArbQueue<EntryType>::RemoveFromFront(EntryType *entry)
{
    if (NumQueueEntries == 0) return(false) ;

    *entry = head->entry ;
    NumQueueEntries-- ;

    if (NumQueueEntries == 0) {
        delete head ;
        head = tail = 0 ;
    } else if (NumQueueEntries == 1) {
        delete head ;
        head = tail ;
        head->prev = 0 ;
    } else {
        ArbQueueEntry *nextPtr = head->next ;
        nextPtr->prev = 0 ;
        delete head ;
        head = nextPtr ;
    }
    return(true) ;
}


#endif
