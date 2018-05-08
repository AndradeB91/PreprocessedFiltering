//
// Heap Template Class header file
//
// Description -
//   This class implements a "heap" based priority queue.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2007
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//

#ifndef Heap_hh
#define Heap_hh

//#define DEBUG_HEAP
#ifdef DEBUG_HEAP
#include <iostream>
#endif

#define HEAP_CACHE_BLOCK_SIZE 100

namespace FTools {

/// \file
///
/// \brief This file defines the Heap class

/** \brief An abstract templated Heap data structure
 *
 *  This class implements a templated Heap data structure that can be
 *  used to implement priority queues or a heap sort.  To use this
 *  class a client must create a derived class and overide the Compare method.
 */

template<class EntryType>
class Heap {
    public:

        struct HandleData {
            int handle_data ;
            HandleData *next ;
        } ;

        typedef HandleData *EntryHandle ;
 
        Heap(int initial_size_hint = 100) ;
        Heap(const Heap<EntryType> &other) ;
        virtual ~Heap() ;

        /// comparison method should return -1, 0, or 1

        virtual int Compare(const EntryType&,const EntryType&) = 0 ;

        /// insert a new entry into the heap

        void Insert(const EntryType &entry) ;

        /// removes the minimum item fron the heap and returns a pointer to it

        EntryType *GetMin() ;

        /// returns a pointer to the minimum item in the heap without removing it

        EntryType *ViewMin() ;

        /// true if the item is in the heap

        bool Contains(EntryType &entry) ;

        /// remove an item from the heap

        bool Remove(EntryType &entry) ;

        /// remove all items from the heap

        void Clear() ;

        /// insert an entry and return a handle to it

        EntryHandle InsertWithHandle(EntryType &entry) ;

        /// remove an entry by giving a handle

        bool RemoveWithHandle(EntryHandle handle) ;

        /// view the item with the given handle

        EntryType *ViewWithHandle(EntryHandle handle) ;

        /// returns a pointer to the list of all heap entries

        EntryType **GetEntryList() ;

        /// return the number of entries in the heap

        int Len() { return(NumHeapEntries) ; } ;

#ifdef DEBUG_HEAP

        void StartDebugTrace(EntryType &entry) {
            doing_trace = true ;
            debug_entry = entry ;
        }

#endif

    private:

        struct HeapEntry {
            EntryType    entry ;
            EntryHandle handle ;
        } ;

        struct HandleCache {
            struct HandleCache *next ;
            HandleData entries[HEAP_CACHE_BLOCK_SIZE] ;
        } ;

        HeapEntry *Data ;
        int NumHeapEntries ;
        int NumAllocated ;

        HandleData  *FreeList ;
        HandleCache *CacheList ;
 
        void Exchange(HeapEntry *data,
                      const int i,
                      const int j) ;
        void FixUp(int k) ;
        void FixDown(int k) ;

        EntryHandle NewHandle() ;
        void DeleteHandle(EntryHandle handle) ;

#ifdef DEBUG_HEAP

        bool doing_trace ;
        EntryType debug_entry ;
        void DebugTrace() ;

#endif
} ;


#ifdef MEMDEBUG
#include "MemDbg.hpp"
//#define new new(__FILE__,__LINE__)
#endif

template<class EntryType>
Heap<EntryType>::
Heap(int initial_size_hint)
{
    NumHeapEntries = 0 ;
    NumAllocated = unsigned(1.1 * (initial_size_hint+1)) ;
    Data = new HeapEntry[NumAllocated] ;
    FreeList = 0 ;
    CacheList = 0 ;
#ifdef DEBUG_HEAP
    doing_trace = false ;
#endif
}

template<class EntryType>
Heap<EntryType>::Heap(const Heap<EntryType> &other)
{
    NumHeapEntries = other.NumHeapEntries ;
    NumAllocated = int(1.1 * NumHeapEntries) ;
    if (NumAllocated < NumHeapEntries+2)
        NumAllocated = NumHeapEntries+2 ;
    Data = new HeapEntry[NumAllocated] ;
    for (int i=1 ; i<=NumHeapEntries ; ++i) {
        Data[i] = other.Data[i] ;
        if (Data[i].handle != 0) {
            EntryHandle handle = NewHandle() ;
            handle->handle_data = other.Data[i].handle->handle_data ;
            Data[i].handle = handle ;
        }
    }
}

template<class EntryType>
Heap<EntryType>::~Heap()
{
    HandleCache *ptr = CacheList ;
    while (ptr != 0) {
        HandleCache *tmp = ptr->next ;
        delete ptr ;
        ptr = tmp ;
    }
    delete [] Data ;
}

template<class EntryType>
void Heap<EntryType>::Insert(const EntryType &entry)
{
    // Check to see if there is room.  If not the double the
    // size of the heap.

    if (NumHeapEntries+1 == NumAllocated) {
        HeapEntry *tmp = new HeapEntry[2*NumAllocated] ;
        for (int i=1 ; i<=NumHeapEntries ; ++i) {
            tmp[i] = Data[i] ;
        }
        delete [] Data ;
        Data = tmp ;
        NumAllocated *= 2 ;
    }

    // Add this entry at the end of the heap

    ++NumHeapEntries ;
    Data[NumHeapEntries].entry = entry ;
    Data[NumHeapEntries].handle = 0 ;

    // bubble this up the heap to the proper position.

    FixUp(NumHeapEntries) ;

#ifdef DEBUG_HEAP
    if (doing_trace) DebugTrace() ;
#endif
}

template<class EntryType>
EntryType *Heap<EntryType>::GetMin()
{
    if (NumHeapEntries == 0) return(0) ;

    // temporarily store the entry at the end of the heap. Note
    // that entry 0 of the Data array is not used.

    Exchange(Data,1,NumHeapEntries) ;

    // reduce the size of the heap and fix it up

    --NumHeapEntries ;
    FixDown(1) ;

    if (Data[NumHeapEntries+1].handle != 0) {
        DeleteHandle(Data[NumHeapEntries+1].handle) ;
        Data[NumHeapEntries+1].handle = 0 ;
    }

    // return pointer to what was the smallest entry

    return(&Data[NumHeapEntries+1].entry) ;
}

template<class EntryType>
EntryType *Heap<EntryType>::ViewMin()
{
    if (NumHeapEntries == 0) return(0) ;
    return(&Data[1].entry) ;
}

template<class EntryType>
bool Heap<EntryType>::Contains(EntryType &entry)
{
    // Search all elements for this entry.

    for (int i=1 ; i <= NumHeapEntries ; ++i) {
        if (Data[i].entry == entry) return(true) ;
    }
    return(false) ;
}

template<class EntryType>
bool Heap<EntryType>::Remove(EntryType &entry)
{
    // Search all elements for this entry.

    for (int i=1 ; i <= NumHeapEntries ; ++i) {
        if (Data[i].entry == entry) {

            if (Data[i].handle != 0) {
                DeleteHandle(Data[i].handle) ;
                Data[i].handle = 0 ;
            }

            // found a match, so copy the last item to this location

            Exchange(Data,i,NumHeapEntries) ;
            --NumHeapEntries ;

            // either fix up or fix down

            if ((i > 1) &&
                (Compare(Data[i].entry,Data[i/2].entry) < 0)) {
                FixUp(i) ;
            } else {
                FixDown(i) ;
            }

#ifdef DEBUG_HEAP
            if (doing_trace) DebugTrace() ;
#endif
            return(true) ;
        }
    }
    return(false) ;
}

template<class EntryType>
void Heap<EntryType>::Clear()
{
    NumHeapEntries = 0 ;
}

template<class EntryType>
EntryType **Heap<EntryType>::GetEntryList()
{
    EntryType **entries = new EntryType*[NumHeapEntries] ;

    for (int i=1 ; i <= NumHeapEntries ; ++i) {
        entries[i-1] = &Data[i].entry ;
    }
    return(entries) ;
}

template<class EntryType>
void Heap<EntryType>::Exchange(HeapEntry *data,
                                   const int i,
                                   const int j)
{
    HeapEntry tmp ;
    if (data[i].handle != 0) data[i].handle->handle_data = j ;
    if (data[j].handle != 0) data[j].handle->handle_data = i ;
    tmp = data[i] ; data[i] = data[j] ; data[j] = tmp ;

}

template<class EntryType>
void Heap<EntryType>::FixUp(int k)
{
    while ((k > 1) && (Compare(Data[k/2].entry,Data[k].entry) > 0)) {
//    while ((k > 1) && (Compare(Data[k/2].entry,Data[k].entry) < 0)) {
        Exchange(Data,k,k/2) ;
        k = k/2 ;
    }
}

template<class EntryType>
void Heap<EntryType>::FixDown(int k)
{
    while (2*k <= NumHeapEntries) {
        int j = 2*k ;
        if ((j < NumHeapEntries) &&
            (Compare(Data[j].entry,Data[j+1].entry) > 0)) ++j ;
        if (Compare(Data[k].entry,Data[j].entry) <= 0) break ;
        Exchange(Data,k,j) ;
        k = j ;
    }
}

template<class EntryType>
typename Heap<EntryType>::EntryHandle
Heap<EntryType>::InsertWithHandle(EntryType &entry)
{
#ifdef DEBUG_INSERT
    std::cout << "Inserting: " << entry << std::endl ;
#endif

    // Check to see if there is room.  If not the double the
    // size of the heap.

    if (NumHeapEntries+1 == NumAllocated) {
        HeapEntry *tmp = new HeapEntry[2*NumAllocated] ;
        for (int i=1 ; i<=NumHeapEntries ; ++i) {
            tmp[i] = Data[i] ;
        }
        delete [] Data ;
        Data = tmp ;
        NumAllocated *= 2 ;
    }

    // Add this entry at the end of the heap

    EntryHandle handle = NewHandle() ;

    ++NumHeapEntries ;
    Data[NumHeapEntries].entry = entry ;
    Data[NumHeapEntries].handle = handle ;
    handle->handle_data = NumHeapEntries ;

    // bubble this up the heap to the proper position.

    FixUp(NumHeapEntries) ;

#ifdef DEBUG_HEAP
    if (doing_trace) DebugTrace() ;
#endif

    return(handle) ;
}

template<class EntryType>
bool Heap<EntryType>::RemoveWithHandle(EntryHandle handle)
{
    int i = handle->handle_data ;

    if (i <= NumHeapEntries) {
        if (Data[i].handle != 0) {
            DeleteHandle(Data[i].handle) ;
            Data[i].handle = 0 ;
        }

        // copy the last item to this location

        if (i == NumHeapEntries) {
            --NumHeapEntries ;
            return(true) ;
        }

        Exchange(Data,i,NumHeapEntries) ;
        --NumHeapEntries ;

        // either fix up or fix down

        if ((i > 1) &&
            (Compare(Data[i].entry,Data[i/2].entry) < 0)) {
            FixUp(i) ;
        } else {
            FixDown(i) ;
        }
#ifdef DEBUG_HEAP
        if (doing_trace) DebugTrace() ;
#endif
        return(true) ;
    }
    return(false) ;
}

template<class EntryType>
EntryType *Heap<EntryType>::ViewWithHandle(EntryHandle handle)
{
    int i = handle->handle_data ;

    if (i <= NumHeapEntries) return(&Data[i].entry) ;
    return(0) ;
}

template<class EntryType>
typename Heap<EntryType>::EntryHandle 
Heap<EntryType>::NewHandle()
{
    if (FreeList == 0) {
        HandleCache *ctmp = new HandleCache ;
        ctmp->next = CacheList ;
        CacheList = ctmp ;
        for (int i=0 ; i<HEAP_CACHE_BLOCK_SIZE-1 ; ++i) {
            ctmp->entries[i].next = &(ctmp->entries[i+1]) ;
        }
        ctmp->entries[HEAP_CACHE_BLOCK_SIZE-1].next = 0 ;
        FreeList = &(ctmp->entries[0]) ;
    }
    HandleData *tmp = FreeList ;
    FreeList = tmp->next ;
    return(tmp) ;
}

template<class EntryType>
void Heap<EntryType>::DeleteHandle(EntryHandle handle)
{
    handle->next = FreeList ;
    FreeList = handle ;
}

#ifdef DEBUG_HEAP
template<class EntryType>
void Heap<EntryType>::DebugTrace()
{
    for (int i=1 ; i <= NumHeapEntries ; ++i) {
        if (Data[i].entry == debug_entry) {
            std::cerr << "Item found at: " << i << std::endl ;
            return ;
        }
    }
    return ;
}
#endif

/*
TEMPLATE CLASS Heap<class EntryType>

  This template class implements heap data structure. This data 
  structure can be used as a priority queue. Clients define a function 
  that ranks heap entries. Entries can be added to the heap in any 
  order. The GetMin function returns the entry in the heap currently 
  ranked the lowest. 

  Entries can inserted into the heap with or without a handle. A handle 
  is a token that is returned to the client. The client can use the 
  handle to delete entries in the heap that are not currently ranked 
  lowest. 


  Template Arguments:

    EntryType - Type of data to be stored in the heap.


PUBLIC INTERFACE

  Public Typedefs:

        int *EntryHandle - Defines the data type that is returned as an 
                entry handle. 

  Public Member Functions:

    Heap - heap constructor 

      Heap(
              int (*compare_func)(),
              int initial_size_hint = 100)

        compare_func      - (in)  a comparison function 
        initial_size_hint - (in)  a hint for the initial allocation 
                                  size 

      Description: This is a constructor for a heap object. As an 
          parameter it takes the address of a function that ranks 
          heap entries. The function has the following interface: 

          int (*func)(const EntryType&,const EntryType&) 

          This function returns -1, 0, or 1 depending on if it's 
          first argument is less than, equal to, or greater than its 
          second. 

          A copy is made of all entries and the copy is stored in the 
          heap. Therefore if the values are large and stored 
          elseware, a reasonable strategy is to build a heap of 
          pointers to these values. 


    Heap - heap destructor 

      ~Heap()

      Description: This is a destructor for the heap object. 


    Insert - insert an entry into the heap 

      void Insert(EntryType &entry)

        entry - (in)  entry to be added. 

      Description: This function inserts a new entry into the heap. 


    GetMin - get the lowest ranked entry in the heap 

      EntryType *GetMin()

      Description: This function returns a pointer to the the entry 
          in the table that is currently ranked lowest. This entry is 
          removed from the heap. 

      Return Value: A pointer to the entry in the table that is 
          currently ranked lowest. 


    ViewMin - return a pointer to the lowest ranked entry 

      EntryType *ViewMin()

      Description: This function returns a pointer to the entry in 
          the table that is currently ranked lowest. This entry is 
          not removed from the heap. 

      Return Value: A pointer to the entry in the table that is 
          currently ranked lowest. 


    Contains - check to see if an entry is in the heap 

      bool Contains(EntryType &entry)

        entry - (in)  entry to search for 

      Description: This function checks to see if the entry specified 
          as an argument is currently stored in the heap. 

      Return Value: True if the entry is found in the heap, false 
          otherwise. 


    Remove - remove an entry from the heap 

      bool Remove(EntryType &entry)

        entry - (in)  entry to delete 

      Description: This function checks to see if the entry specified 
          as an argument is currently stored in the heap. If found 
          the entry is removed from the heap. If one wants to remove 
          an item from heap but does not have the full entry 
          information, the RemoveWithHandle interface may be useful. 

      Return Value: True if an item was found and removed, false 
          otherwise. 


    InsertWithHandle - insert an entry and return a handle 

      EntryHandle InsertWithHandle(EntryType &entry)

        entry - (in)  entry to be added 

      Description: This function inserts an item into the heap and 
          returns a handle, or token, that can be used to reference 
          the entry. 

      Return Value: A handle that can be used to referenc the item. 


    RemoveWithHandle - remove an entry referenced by a handle 

      bool RemoveWithHandle(EntryHandle handle)

        handle - (in)  handle of the item to remove 

      Description: This function removes an entry from the heap where 
          the entry is referenced by its handle. 

      Return Value: True if an item is removed, false otherwise. 


    ViewWithHandle - view an entry referenced by a handle 

      EntryType *ViewWithHandle(EntryHandle handle)

        handle - (in)  handle of the item to be viewed 

      Description: This function returns a pointer to a heap entry 
          that is referenced by a handle. 

      Return Value: A pointer to the reference entry. 


    GetEntryList - returns a list of all entries in the table. 

      EntryType **GetEntryList()

      Description: This function returns a list (array) of all the 
          entries in the table. The client takes ownership of the 
          list and must delete it by a call to delete []. 

      Return Value: The address of an array containing all the heap 
          entries. 


    NumEntries - number of entries 

      int NumEntries()

      Description: This function returns the number of entries 
          currently in the heap. 

      Return Value: Number of entries currently in the heap. 


PRIVATE INTERFACE

  Private Data Structures:

    struct HeapEntry

      this structure used to store entries in the heap 

      Member Variables:

        EntryType entry - the entry's value 

        EntryHandle handle - the entry's handle 


  Private Member Functions:

    Exchange - exchange two entries 

      void Exchange(
              HeapEntry *data,
              const int    i,
              const int    j)

        data - (i/o) the heap's data 
        i    - (in)  index of the first item 
        j    - (in)  index of the second item 

      Description: This function is call to exhange the location of 
          two entries in the heap. 


    FixUp - bubble a new entry up the tree 

      void FixUp(int k)

        k - (in)  index of the item to bubble up 

      Description: This function "bubbles" an entry up the heap's 
          binary tree until it is ranked smaller than all its 
          decendants. 


    FixDown - bubble an entry down the tree 

      void FixDown(int k)

        k - (in)  index of the item to bubble down 

      Description: This function "bubbles" an entry down the heap's 
          binary tree until it is ranked larger than its parent. 


  Private Member Variables:

    HeapEntry *pData - the heap's data 

    int NumHeapEntries - current number of entries 

    int NumAllocated - current number of entries allocated 

    int (*Compare)() - pointer to the compare function 

*/

} // namespace

#endif
