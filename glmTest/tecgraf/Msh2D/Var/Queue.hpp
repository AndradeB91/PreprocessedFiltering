//
// Queue Template Class header file
//
// Description -
//   This class implements a FIFO stack.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2007
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//

#ifndef Queue_hh
#define Queue_hh

namespace FTools {

#define SMALL_HASH_TABLE 4

/// \file
///
/// \brief This file defines the Queue class

/** \brief A templated variable FIFO queue
 *
 *  This class implements a templated a FIFO queue
 */

template<class EntryType>
class Queue {
    public:

        Queue() : head(0), tail(0), NumQueueEntries(0) {}
        ~Queue() ;

        /// append an item to the back of the queue

        void AddToBack(const EntryType &value) ;

        /// remove an item from the front of the queue

        bool RemoveFromFront(EntryType *value) ;

        /// return the number of entries in the queue

        int Len() { return(NumQueueEntries) ; } ;

    private:

        struct QueueEntry {
            EntryType    entry ;
            QueueEntry *prev ;
            QueueEntry *next ;
        } ;

        QueueEntry *head ;
        QueueEntry *tail ;
        int NumQueueEntries ;
} ;


template<class EntryType> Queue<EntryType>::~Queue()
{
    if (NumQueueEntries > 0) {
        QueueEntry *nextPtr, *hPtr = head ;
        while (hPtr != 0) {
            nextPtr = hPtr->next ;
	    delete hPtr ;
	    hPtr = nextPtr ;
	}
    }
}

template<class EntryType>
void Queue<EntryType>::AddToBack(const EntryType &entry)
{
    QueueEntry *hPtr = new QueueEntry ;
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

template<class EntryType>
bool Queue<EntryType>::RemoveFromFront(EntryType *entry)
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
        QueueEntry *nextPtr = head->next ;
        nextPtr->prev = 0 ;
        delete head ;
        head = nextPtr ;
    }
    return(true) ;
}

} // namespace

#endif
