//
// CArbSmallSet Template Class header file
//
// Description -
//   This class implements a set class for a sets with only
//   a handful of members
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2001,2002
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.3 $  $Date: 2002/06/11 14:20:00 $  $Author: wash $
//

#ifndef ArbSmallSet_hh
#define ArbSmallSet_hh

template<class EntryType,const int def_size> class CArbConstSmallSetIterator ;

template<class EntryType,const int def_size>
class CArbSmallSet {
    public:

        CArbSmallSet()
        {
            NumAllocated = def_size ;
            NumElems = 0 ;
            data = &StaticData[0] ;
        } ;

        CArbSmallSet(const CArbSmallSet &other)
        {
            NumAllocated = other.NumAllocated ;
            NumElems = other.NumElems ;
            data = &StaticData[0] ;
            if (other.data != &other.StaticData[0]) {
                data = new EntryType[NumAllocated] ;
            }
            for (int i=0 ; i<NumElems ; ++i) {
                data[i] = other.data[i] ;
            }
        } ;

        CArbSmallSet operator = (const CArbSmallSet &other)
        {
            if (data != &StaticData[0]) delete [] data ;
            NumAllocated = other.NumAllocated ;
            NumElems = other.NumElems ;
            data = &StaticData[0] ;
            if (other.data != &other.StaticData[0]) {
                data = new EntryType[NumAllocated] ;
            }
            for (int i=0 ; i<NumElems ; ++i) {
                data[i] = other.data[i] ;
            }
            return(*this) ;
        } ;

        ~CArbSmallSet() 
        {
            if (data != &StaticData[0]) delete [] data ;
        } ;

        void Insert(const EntryType &entry)
        {
            if (NumElems >= NumAllocated) {
                EntryType *tmp = new EntryType[2*NumAllocated] ;
                for (int i=0 ; i<NumElems ; ++i)
                    tmp[i] = data[i] ;
                if (data != &StaticData[0]) delete [] data ;
                data = tmp ;
                NumAllocated *= 2 ;
            }
            data[NumElems] = entry ;
            ++NumElems ;
        } ;

        void Remove(const EntryType &entry)
        {
            for (int i=0 ; i<NumElems ; ++i) {
                if (data[i] == entry) {
                    for (int j=i ; j<NumElems-1 ; ++j) {
                        data[j] = data[j+1] ;
                    }
                    return ;
                }
            }
        } ;

        int NumElements() const { return(NumElems) ; } ;

        bool HasElement(const EntryType &entry)
        {
            for (int i=0 ; i<NumElems ; ++i) {
                if (data[i] == entry) return(true) ;
            }
            return(false) ;
        } ;

        EntryType &Element(const int i) const { return(data[i]) ; } ;

    private:

        EntryType StaticData[def_size], *data ;
        int NumElems ;
        int NumAllocated ;

    friend class CArbConstSmallSetIterator<EntryType,def_size> ;

#ifdef USE_STL_IO
#if defined(WIN32) || defined(WIN64)
    friend ostream &operator << (ostream &out,
                const CArbSmallSet<EntryType,def_size> &set) ;
#else
    friend ostream &operator << <EntryType,def_size> (ostream &out,
                const CArbSmallSet<EntryType,def_size> &set) ;
#endif
#endif
} ;


template<class EntryType,const int def_size>
class CArbConstSmallSetIterator {
    public:

        CArbConstSmallSetIterator(
             const CArbSmallSet<EntryType,def_size> *const aSet) :
             set(aSet),cur(0) {} ;
        void First()        { cur = 0 ; } ;
        void Next()         { ++cur ; } ;
        bool More()         { return(cur < set->NumElems) ; } ;
        EntryType Current() { return(cur < set->NumElems ?
                                     set->data[cur] : 0) ; } ;
        EntryType operator * ()  { return(Current()) ; } ;
        void operator ++ ()      { Next() ; } ;
        void operator ++ (int i) { Next() ; } ;

    private:
        const CArbSmallSet<EntryType,def_size> *set ;
        int cur ;
} ;

#ifdef USE_STL_IO
#include <iostream.h>

template<class EntryType,const int def_size>
inline ostream &operator << (ostream &out,
                    const CArbSmallSet<EntryType,def_size> &set)
{
    out << "[" ;
    for (int i=0 ; i<set.NumElems ; ++i) {
        if (i > 0) out << ' ' ;
        out << set.data[i] ;
    }
    out << ']' ;
    return(out) ;
}
#endif

/* ----------------------------------------------------------------
class CArbSmallSet<class EntryType,class int>

  This object implements a set data structure. The elements of the set 
  are maintained in an array with linear searching, which is 
  inefficient fall all but small sets. 


  Template Arguments:

    EntryType - type of set elements

    int - Default set allocation size.  If the set is
           <= to this number no additional allocations are
           performed.

  Public Member Functions:

    CArbSmallSet - constructor 

      CArbSmallSet()

      Description: Constructor method for ArbSmallSet 


    CArbSmallSet - copy constructor 

      CArbSmallSet(const CArbSmallSet &other)

        other - (in)  object to copy 

      Description: Copy constructor for an ArbSmallSet object 

      Return Value: the new object 


    operator_= - assignment operator 

      CArbSmallSet operator = (const CArbSmallSet &other)

        other - (in)  object to copy 

      Description: Assignment operator for an ArbSmallSet 

      Return Value: the updated object 


    CArbSmallSet - destructor 

      ~CArbSmallSet()

      Description: Destructor method for ArbSmallSet 


    Insert - insert an element 

      void Insert(const EntryType &entry)

        entry - (in)  element to add 

      Description: Inserts a new element into the set. 


    Remove - remove an element 

      void Remove(const EntryType &entry)

        entry - (in)  element to remove 

      Description: Removes an element from the set. 


    NumElements - number of elements 

      int NumElements() const

      Description: Returns the number of elements in the set. 

      Return Value: number of elements 


    HasElement - check for element 

      bool HasElement(const EntryType &entry)

        entry - (in)  element to check 

      Description: Returns true if the given element is in the set. 

      Return Value: true => element is in the set 


    Element - return an element 

      EntryType &Element(const int i) const

        i - (in)  element index 

      Description: Returns the element at the given index 

      Return Value: one element of the set 


  Private Member Variables:

        EntryType StaticData[def_size] - data storage if the number of 
                elements is small enough 

        EntryType *data - pointer to the set data (may be the address 
                of StaticData) 

        int NumElems - number of elements in the set 

        int NumAllocated - number of slots currently allocated 



class CArbConstSmallSetIterator<>

  Iterator object for SmallSets 


  Template Arguments:

  Public Member Functions:

    CArbConstSmallSetIterator - constructor 

      CArbConstSmallSetIterator(const CArbSmallSet<EntryType,def_size>*const aSet)

        aSet - (in)  pointer to the set to interate through 

      Description: Constructor method. 


    First - point to first element 

      void First()

      Description: Sets (resets) the iterator to point to the first 
          element. 


    Next - point to next element 

      void Next()

      Description: Moves the iterator to point to the next element. 


    More - check for available data 

      bool More()

      Description: Checks to see if the iterator points to a valid 
          element. 

      Return Value: true => valid data, false => end of data 


    Current - return the current element 

      EntryType Current()

      Description: Returns the currently referenced element. 

      Return Value: the current element 


    operator_* - dereference operator 

      EntryType operator * ()

      Description: Dereference operator, returns the current element. 

      Return Value: the current element 


    operator_++ - increment operator 

      void operator ++ ()

      Description: Calls the Next method. 


    operator_++ - increment operator 

      void operator ++ (int i)

        i - (in)  ignored 

      Description: Calls the Next method. 

  Private Member Variables:

        const CArbSmallSet<EntryType,def_size>* set - pointer to the 
                set that is being iterated over 

        int cur - current index in the referenced set 


*/

#endif

