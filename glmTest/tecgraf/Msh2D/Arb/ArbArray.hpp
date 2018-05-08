//
// CArbArray Template Class header file
//
// Description -
//   This class implements an arbitrary size array class.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.7 $  $Date: 2001/06/15 13:30:02 $  $Author: wash $
//

#ifndef ArbArray_hh
#define ArbArray_hh

template<class EntryType>
class CArbArray {
    public:

        CArbArray(int initial_size_hint = 10) ;
        CArbArray(const CArbArray &other) ;
        CArbArray operator = (const CArbArray &other) ;
        virtual ~CArbArray() ;

        void InsertAtEnd(const EntryType &entry) ;
        void Clear() ;
        int NumEntries() const { return(NumArrayEntries) ; } ;

        EntryType *AsVector() const ;
        EntryType &At(const int i) ;
        EntryType At(const int i) const ;

        EntryType &operator[] (const int i) ;
        EntryType operator[] (const int i) const ;

    private:

        EntryType *data ;
        int NumArrayEntries ;
        int NumAllocated ;
} ;


#ifdef MEMDEBUG
#include "MemDbg.hpp"
#define new new(__FILE__,__LINE__)
#endif

// force templated code to be generated

// %(CArbArray::CArbArray-constructor-|-int-|)
/* ++ ----------------------------------------------------------
**
**    CArbArray - constructor for array objects
**
**      CArbArray(int initial_size_hint = 10)
**
**        initial_size_hint - (in)  this parameters gives a hint about 
**                                  the initial number of array slots for 
**                                  which memory should be allocated 
**
**      Description: This is the constructor method for an array object. 
**
**
** -- */


template<class EntryType>
CArbArray<EntryType>::
CArbArray(int initial_size_hint)
{
    NumArrayEntries = 0 ;
    NumAllocated = initial_size_hint ;
    data = new EntryType[NumAllocated] ;
}




// %(CArbArray::CArbArray-constructor-|-CArbArray-const|&)
/* ++ ----------------------------------------------------------
**
**    CArbArray - copy constructor
**
**      CArbArray(const CArbArray &other)
**
**        other - (in)  object to be copied 
**
**      Description: This is the copy constructor for array objects. 
**
**
** -- */


template<class EntryType>
CArbArray<EntryType>::
CArbArray(const CArbArray &other)
{
    NumArrayEntries = other.NumArrayEntries ;
    NumAllocated    = other.NumAllocated ;
    data = new EntryType[NumAllocated] ;
    for (int i=0 ; i<NumArrayEntries ; ++i)
        data[i] = other.data[i] ;
}




// %(CArbArray::operator_=-CArbArray-|-CArbArray-const|&)
/* ++ ----------------------------------------------------------
**
**    operator_= - assignment operator
**
**      CArbArray operator = (const CArbArray &other)
**
**        other - (in)  arry to be copied 
**
**      Description: This implements the assignment operator for array 
**          objects. 
**
**      Return Value: This function returns a copy of the array. 
**
**
** -- */


template<class EntryType> CArbArray<EntryType>
CArbArray<EntryType>::operator =(const CArbArray &other)
{
    delete [] data ;
    NumArrayEntries = other.NumArrayEntries ;
    NumAllocated    = other.NumAllocated ;
    data = new EntryType[NumAllocated] ;
    for (int i=0 ; i<NumArrayEntries ; ++i)
        data[i] = other.data[i] ;
    return(*this) ;
}




// %(CArbArray::CArbArray-destructor-virtual|~)
/* ++ ----------------------------------------------------------
**
**    CArbArray - destructor
**
**      ~CArbArray()
**
**      Description: This is the destructor function for array objects. 
**
**
** -- */


template<class EntryType> CArbArray<EntryType>::~CArbArray()
{
    delete [] data ;
}




// %(CArbArray::InsertAtEnd-void-|-EntryType-const|&)
/* ++ ----------------------------------------------------------
**
**    InsertAtEnd - insert an entry at the end of the array
**
**      void InsertAtEnd(const EntryType &entry)
**
**        entry - (in)  Entity to be added at the end of the array. 
**
**      Description: This function adds an entry to the end of the Array. 
**          The array object makes a copy of the argument data. 
**
**
** -- */


template<class EntryType>
void CArbArray<EntryType>::InsertAtEnd(const EntryType &entry)
{
    // Check to see if there is room.  If not the double the
    // size of the Array.

    if (NumArrayEntries == NumAllocated) {
        EntryType *tmp = new EntryType[2*NumAllocated] ;
        for (int i=0 ; i<NumArrayEntries ; ++i) {
            tmp[i] = data[i] ;
        }
        delete [] data ;
        data = tmp ;
        NumAllocated *= 2 ;
    }

    // Add this entry at the end of the Array

    data[NumArrayEntries] = entry ;
    ++NumArrayEntries ;
}




// %(CArbArray::Clear-void-|) 
/* ++ ----------------------------------------------------------
**
**    Clear - clear all entries in the array
**
**      void Clear()
**
**      Description: This function clears (deletes) all entries in the 
**          array. 
**
**
** -- */


template<class EntryType> void CArbArray<EntryType>::Clear()
{
    NumArrayEntries = 0 ;
}




// %(CArbArray::AsVector-EntryType-|*^const) 
/* ++ ----------------------------------------------------------
**
**    AsVector - return the array data in a vector
**
**      EntryType *AsVector() const
**
**      Description: This function return the data in the array as a 
**          vector (a standard C++ array). The calling routine takes 
**          ownership of the pointer. 
**
**      Return Value: The data in the array stored in a standard C++ 
**          array. 
**
**
** -- */


template<class EntryType>
EntryType *CArbArray<EntryType>::AsVector() const
{
    EntryType *vector = new EntryType[NumArrayEntries] ;
    for (int i=0 ; i<NumArrayEntries ; ++i)
        vector[i] = data[i] ;
    return(vector) ;
}




// %(CArbArray::At-EntryType-|&-int-const|)
/* ++ ----------------------------------------------------------
**
**    At - return an array entry
**
**      EntryType &At(const int i)
**
**        i - (in)  index of the entry to return 
**
**      Description: This function returns a reference to one of the 
**          entries of the array. This reference can be used as an 
**          l-value. 
**
**      Return Value: A reference to one of the array entries. 
**
**
** -- */


template<class EntryType>
EntryType &CArbArray<EntryType>::At(const int i)
{
    int indx = (i>=0) ? i : NumArrayEntries+i ;

    if (indx > NumAllocated-1) {
        int new_size = int(indx * 1.5) ;
        EntryType *tmp = new EntryType[new_size] ;
        int j ;
        for (j=0 ; j<NumArrayEntries ; ++j) {
            tmp[j] = data[j] ;
        }
        delete [] data ;
        data = tmp ;
        NumAllocated = new_size ;
    }
    if (indx+1 > NumArrayEntries) NumArrayEntries = indx+1 ; 
    return(data[indx]) ;
}




// %(CArbArray::At-EntryType-|^const-int-const|) 
/* ++ ----------------------------------------------------------
**
**    At - return an array entry
**
**      EntryType At(const int i) const
**
**        i - (in)  index of the entry to return 
**
**      Description: This function returns one of the entries of the 
**          array. 
**
**      Return Value: One of the array entries. 
**
**
** -- */


template<class EntryType>
EntryType CArbArray<EntryType>::At(const int i) const
{
    // assert((i >= -NumArrayEntries) && (i < NumArrayEntries)) ; 
    return((i>=0) ? data[i] : data[NumArrayEntries+i]) ;
}




// %(CArbArray::operator_[]-EntryType-|&-int-const|) 
/* ++ ----------------------------------------------------------
**
**    operator_[] - access operator for the array
**
**      EntryType &operator [] (const int i)
**
**        i - (in)  index of the array entry 
**
**      Description: This operator returns a reference to one of the 
**          entries in the array. This reference can be used as an 
**          l-value. 
**
**      Return Value: A reference to an array entry. 
**
**
** -- */


template<class EntryType>
EntryType &CArbArray<EntryType>::operator[] (const int i)
{
    int indx = (i>=0) ? i : NumArrayEntries+i ;

    if (indx > NumAllocated-1) {
        int new_size = int(indx * 1.5) ;
        EntryType *tmp = new EntryType[new_size] ;
        int j ;
        for (j=0 ; j<NumArrayEntries ; ++j) {
            tmp[j] = data[j] ;
        }
        delete [] data ;
        data = tmp ;
        NumAllocated = new_size ;
    }
    if (indx+1 > NumArrayEntries) NumArrayEntries = indx+1 ; 
    return(data[indx]) ;
}




// %(CArbArray::operator_[]-EntryType-|^const-int-const|) 
/* ++ ----------------------------------------------------------
**
**    operator_[] - access operator for the array
**
**      EntryType operator [] (const int i) const
**
**        i - (in)  index of the array entry 
**
**      Description: This operator returns one of the entries in the 
**          array. 
**
**      Return Value: an array entry 
**
**
** -- */


template<class EntryType>
EntryType CArbArray<EntryType>::operator[] (const int i) const
{
    // assert((i >= -NumArrayEntries) && (i < NumArrayEntries)) ; 
    return((i>=0) ? data[i] : data[NumArrayEntries+i]) ;
}

#endif
