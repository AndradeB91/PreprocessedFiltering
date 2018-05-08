//
// List.hpp
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2007
//   All rights reserved
//
// ------------------------------------------------------------------
//
//  A List object is a templated dynamically allocated list of objects.
//     A derived class that supports the observer pattern along with
//     an observer mixin class.
//
//  Objects are copied when they are placed in the list.
//
//  Negative indices count back from the end of the list.  That is,
//  l[-1] is the last item in the list, l[-2] the second last, etc.
//
//  Public Interface:
//
//      Constructors:
//          List<type>()
//          List<type>(int initial_size)
//          List<type>(const type* in,int len)
//          List<type>(const List<type>& other)
//          ~List()
//
//      Methods:
//          l.Append(const type& in)       - append to the end of the list
//          l.Insert(int i,const type& in) - insert at position i moving existing 
//          l.RemoveFromEnd()              - removes the last item from the list
//          l.Remove(int i)                - removes the i'th item from the list
//          l.Clear()                      - removes all items from the list
//          l.Len()                        - returns the number items in the list
//          l.At(int i)                    - accesses the i'th item in the list
//          l.Reverse()                    - reverses the order of items in the list
//          l.Extend(const List<type>& l2) - appends the items in list l2
//          l.Iterator()
//          l.ConstIterator()              - returns and iterator (see below)
//
//      Operators:
//          l0 = l1        - assignment
//          l[i]           - accesses the i'th item in the list
//
//
//  Lists can be accessed with iterators.  Typical access with an iterator
//  looks like
//
//  List<type>::ListIterator iter = l.Iterator()
//  for (iter.First() ; iter.More() ; ++iter) {
//      ...
//      do something
//      ...
//  }
//
//  Public Iterator Interface:
//
//      Constructor:
//          ListIterator(List<type>& list)
//
//      Methods:
//          i.First()      - positions the iterator to the first item
//          i.Next()       - moves the iterator to the next item
//          i.More()       - true if the iterator points to a valid item
//          i.Index()      - returns the index of the current item
//          i.Entry()      - returns a reference to the current item
//
//      Operators:
//          ++i, i++       - same as i.Next()
//          &i             - same as i.Index()
//          *i             - same as i.Entry()
//
//
//
//  Observed derived class:
//
//      Constructors:
//          ObsList<type>()
//          ObsList<type>(int initial_size)
//          ObsList<type>(const type* in,int len)
//          ObsList<type>(const List<type>& other) ;
//          ~ObsList<type>()
//
//      Methods:
//          ol.Attach(ListObsrvr<type>* observer) - add an observer
//          ol.Detach(ListObsrvr<type>* observer) - remove an observer
//
//
//
//  Public Interface for an ListObserver:
//
//      Constructors:
//          no constructures, should be a mixin class
//
//      Methods:
//          lo.NotifyUpdate(ObsList<type>* ol,
//                           int index,type value)
//          lo.NotifyAppend(ObsList<type>* ol,type value)
//          lo.NotifyStartBatch(ObsList<type>* ol)
//          lo.NotifyEndBatch(ObsList<type>* ol)
//
//


#ifndef List_hpp
#define List_hpp

#include <cassert>
#include <cstdlib>

namespace FTools {

#define List_STATIC_DATA_SIZE 6

/// \file
///
/// \brief This file defines the List class

/** \brief A templated variable length list
 *
 *  This class implements a templated variable length list
 */


template<class EntryType>
class List {
    public:

        class ListIterator ;
        class ConstListIterator ;

        List(int initial_size_hint = List_STATIC_DATA_SIZE) ;
        List(const EntryType* in,int num,int initial_size_hint = 10) ;
        List(EntryType* in,int num,int initial_size_hint = 10) ;
        List(const List &other) ;
        //List(List &other) ;
        List operator = (const List &other) ;
        //List operator = (List &other) ;
        virtual ~List() ;

        /// appends a value to the end of the list

        void Append(const EntryType &entry) ;
        void Append(EntryType &entry) ;

        /// appends a value to the end of the list only if not already present

        void AddUnique(const EntryType &entry) ;
        void AddUnique(EntryType &entry) ;

        /// insert a value at the index, existing values are shuffled down

        void Insert(int i,const EntryType &entry) ;
        void Insert(int i,EntryType &entry) ;

        /// removes the last value in the list

        void RemoveFromEnd() ;

        /// removes a value at the index, other values are shuffled up

        void Remove(int i) ;

        /// removes a value at the index, other values are shuffled up

        void Del(int i) ;

        /// remove all entries

        void Clear() ;

        /// returns the number of entries in the list

        int Len() const { return(NumArrayEntries) ; } ;

        /// returns a reference to the entry at the index

        EntryType& At(const int i) ;

        /// returns the value of the entry

        const EntryType& At(const int i) const ;

        /// returns a reference to the entry at the index

        EntryType& operator[] (const int i) ;

        /// returns the value of the entry

        const EntryType& operator[] (const int i) const ;

        /// returns a pointer to the list data as an array

        const EntryType* Data() const { return(data) ; } ;

        /// reverses the order of the entries in the list

        void Reverse() ;

        /// append the values from another list

        void Extend(const List<EntryType> &other) ;
        void Extend(List<EntryType> &other) ;

        // This version of extend increases the size of the list by
        // one but does not copy anything into the new location.  This
        // optimization should only be used when the new location will
        // imediately initialized.  Use this method with caution

        void ExtendOne() ;

        /// returns a list iterator

        ListIterator Iterator() { return(ListIterator(*this)) ; } ;

        /// returns a constant list iterator

        ConstListIterator ConstIterator() { return(ConstListIterator(*this)) ; } ;

        /// sort the list entries

        void Sort(int(*compare)(const void*,const void*)) ;

    protected:

        EntryType *data ;
        EntryType StaticData[List_STATIC_DATA_SIZE] ;
        int NumArrayEntries ;
        int NumAllocated ;

    public:

        class ListIterator {
            public:
                ListIterator(
                    List<EntryType>& alist) :
                        list(alist) {
                    index = 0 ;
                } ;

                void First()       { index = 0 ; } ;
                void Next()        { ++index ; } ;
                bool More()        { return(index < list.NumArrayEntries) ; } ;
                int  Index()       { return(index) ; } ;
                EntryType &Entry() { return(list.data[index]) ; } ; 

                void operator ++ ()      { Next() ; } ;
                void operator ++ (int i) { Next() ; } ;

                int         operator & () { return(Index()) ; } ;
                EntryType  &operator * () { return(Entry()) ; } ;

            private:
                List<EntryType>& list ;
                int index ;
        } ;

        class ConstListIterator {
            public:
                ConstListIterator(
                    const List<EntryType>& alist) :
                        list(alist) {
                    index = 0 ;
                } ;

                void First()       { index = 0 ; } ;
                void Next()        { ++index ; } ;
                bool More()        { return(index < list.NumArrayEntries) ; } ;
                int  Index()             { return(index) ; } ;
                const EntryType &Entry() { return(list.data[index]) ; } ; 

                void operator ++ ()      { Next() ; } ;
                void operator ++ (int i) { Next() ; } ;

                int               operator & () { return(Index()) ; } ;
                const EntryType  &operator * () { return(Entry()) ; } ;

            private:
                const List<EntryType>& list ;
                int index ;
        } ;
} ;


template<class EntryType>
    List<EntryType>::List(
        int initial_size_hint)
{
    NumArrayEntries = 0 ;
    NumAllocated = initial_size_hint ;
    if (initial_size_hint > List_STATIC_DATA_SIZE) {
        data = new EntryType[NumAllocated] ;
    } else {
        data = &StaticData[0] ;
        NumAllocated = List_STATIC_DATA_SIZE ;
    }
}


template<class EntryType>
    List<EntryType>::List(
        const EntryType *in,
        int num,
        int initial_size_hint)
{
    NumArrayEntries = num ;
    NumAllocated = (initial_size_hint > num) ? initial_size_hint : num ;
    if (NumAllocated > List_STATIC_DATA_SIZE) {
        data = new EntryType[NumAllocated] ;
    } else {
        data = &StaticData[0] ;
        NumAllocated = List_STATIC_DATA_SIZE ;
    }
    NumArrayEntries = num ;
    for (int i=0 ; i<num ; ++i) data[i] = in[i] ;
}

template<class EntryType>
    List<EntryType>::List(
        EntryType *in,
        int num,
        int initial_size_hint)
{
    NumArrayEntries = num ;
    NumAllocated = (initial_size_hint > num) ? initial_size_hint : num ;
    if (NumAllocated > List_STATIC_DATA_SIZE) {
        data = new EntryType[NumAllocated] ;
    } else {
        data = &StaticData[0] ;
        NumAllocated = List_STATIC_DATA_SIZE ;
    }
    NumArrayEntries = num ;
    for (int i=0 ; i<num ; ++i) data[i] = in[i] ;
}


template<class EntryType>
    List<EntryType>::List(
        const List &other)
{
    NumArrayEntries = other.NumArrayEntries ;
    if (NumArrayEntries > List_STATIC_DATA_SIZE) {
        NumAllocated    = other.NumAllocated ;
        data = new EntryType[NumAllocated] ;
    } else {
        data = &StaticData[0] ;
        NumAllocated = List_STATIC_DATA_SIZE ;
    }
    for (int i=0 ; i<NumArrayEntries ; ++i)
        data[i] = other.data[i] ;
}

// template<class EntryType>
//     List<EntryType>::List(
//     List &other)
// {
//     NumArrayEntries = other.NumArrayEntries ;
//     if (NumArrayEntries > List_STATIC_DATA_SIZE) {
//         NumAllocated    = other.NumAllocated ;
//         data = new EntryType[NumAllocated] ;
//     } else {
//         data = &StaticData[0] ;
//         NumAllocated = List_STATIC_DATA_SIZE ;
//     }
//     for (int i=0 ; i<NumArrayEntries ; ++i)
//         data[i] = other.data[i] ;
// }

template<class EntryType>
    List<EntryType> List<EntryType>::operator =(
        const List &other)
{
    if (data != 0 && data != &StaticData[0]) 
        delete [] data ;
    NumArrayEntries = other.NumArrayEntries ;
    if (NumArrayEntries > List_STATIC_DATA_SIZE) {
        NumAllocated    = other.NumAllocated ;
        data = new EntryType[NumAllocated] ;
    } else {
        data = &StaticData[0] ;
        NumAllocated = List_STATIC_DATA_SIZE ;
    }
    for (int i=0 ; i<NumArrayEntries ; ++i)
        data[i] = other.data[i] ;
    return(*this) ;
}

// template<class EntryType>
//     List<EntryType> List<EntryType>::operator =(
//     List &other)
// {
//     if (data != 0 && data != &StaticData[0]) 
//         delete [] data ;
//     NumArrayEntries = other.NumArrayEntries ;
//     if (NumArrayEntries > List_STATIC_DATA_SIZE) {
//         NumAllocated    = other.NumAllocated ;
//         data = new EntryType[NumAllocated] ;
//     } else {
//         data = &StaticData[0] ;
//         NumAllocated = List_STATIC_DATA_SIZE ;
//     }
//     for (int i=0 ; i<NumArrayEntries ; ++i)
//         data[i] = other.data[i] ;
//     return(*this) ;
// }

template<class EntryType>
    List<EntryType>::~List()
{
    if (data != 0 && data != &StaticData[0]) 
        delete [] data ;
}


template<class EntryType>
    void List<EntryType>::Append(
        const EntryType &entry)
{
    // Check to see if there is room.  If not the double the
    // size of the Array.

    EntryType *realloc_data = 0;

    if ((data == &StaticData[0]) &&
        (NumArrayEntries == List_STATIC_DATA_SIZE)) {
        NumAllocated = 2 * List_STATIC_DATA_SIZE ;
        EntryType *tmp = new EntryType[NumAllocated] ;
        for (int i=0 ; i<List_STATIC_DATA_SIZE ; ++i) {
            tmp[i] = data[i] ;
        }
        data = tmp ;
    //} else if (NumArrayEntries > List_STATIC_DATA_SIZE) {
    } else {
        if (NumArrayEntries == NumAllocated) {
            EntryType *tmp = new EntryType[2*NumAllocated] ;
            for (int i=0 ; i<NumArrayEntries ; ++i) {
                tmp[i] = data[i] ;
            }
            realloc_data = data ;
            //delete [] data ;
            data = tmp ;
            NumAllocated *= 2 ;
        }
    }

    // Add this entry at the end of the Array

    data[NumArrayEntries] = entry ;
    ++NumArrayEntries ;

    if (realloc_data != 0) delete [] realloc_data ;
}

template<class EntryType>
    void List<EntryType>::Append(
        EntryType &entry)
{
    // Check to see if there is room.  If not the double the
    // size of the Array.

    EntryType *realloc_data = 0;

    if ((data == &StaticData[0]) &&
        (NumArrayEntries == List_STATIC_DATA_SIZE)) {
        NumAllocated = 2 * List_STATIC_DATA_SIZE ;
        EntryType *tmp = new EntryType[NumAllocated] ;
        for (int i=0 ; i<List_STATIC_DATA_SIZE ; ++i) {
            tmp[i] = data[i] ;
        }
        data = tmp ;
    //} else if (NumArrayEntries > List_STATIC_DATA_SIZE) {
    } else {
        if (NumArrayEntries == NumAllocated) {
            EntryType *tmp = new EntryType[2*NumAllocated] ;
            for (int i=0 ; i<NumArrayEntries ; ++i) {
                tmp[i] = data[i] ;
            }
            realloc_data = data ;
            //delete [] data ;
            data = tmp ;
            NumAllocated *= 2 ;
        }
    }

    // Add this entry at the end of the Array

    data[NumArrayEntries] = entry ;
    ++NumArrayEntries ;

    if (realloc_data != 0) delete [] realloc_data ;
}

template<class EntryType>
    void List<EntryType>::AddUnique(
        const EntryType &entry)
{
    for (int j=0 ; j<NumArrayEntries ; ++j) {
        if (data[j] == entry) return ;
    }
    Append(entry) ;
}

template<class EntryType>
    void List<EntryType>::AddUnique(
        EntryType &entry)
{
    for (int j=0 ; j<NumArrayEntries ; ++j) {
        if (data[j] == entry) return ;
    }
    Append(entry) ;
}

template<class EntryType>
    void List<EntryType>::Insert(
        int i,
        const EntryType &entry)
{
    // get a local copy of the entry for the special case where it is
    // a reference to a current element and we realloc

    EntryType tmp = entry ;

    Append(data[NumArrayEntries-1]) ;

    // move the current entries

    for (int j=NumArrayEntries-1 ; j>i ; --j) {
        data[j] = data[j-1] ;
    }

    // Add this entry

    data[i] = tmp ;
}

template<class EntryType>
    void List<EntryType>::Insert(
        int i,
        EntryType &entry)
{
    // get a local copy of the entry for the special case where it is
    // a reference to a current element and we realloc

    EntryType tmp = entry ;

    Append(data[NumArrayEntries-1]) ;

    // move the current entries

    for (int j=NumArrayEntries-1 ; j>i ; --j) {
        data[j] = data[j-1] ;
    }

    // Add this entry

    data[i] = tmp ;
}

template<class EntryType>
    void List<EntryType>::RemoveFromEnd()
{
    // Remove the last item in the array

    if (NumArrayEntries > 0) --NumArrayEntries ;
}

template<class EntryType>
    void List<EntryType>::Remove(int i)
{
    // Remove one item in the array

    assert(i < NumArrayEntries) ;
    for (int j=i ; j<NumArrayEntries-1 ; ++j) data[j] = data[j+1] ;
    if (NumArrayEntries > 0) --NumArrayEntries ;
}

template<class EntryType>
    void List<EntryType>::Del(int i)
{
    // Remove one item in the array

    assert(i < NumArrayEntries) ;
    for (int j=i ; j<NumArrayEntries-1 ; ++j) data[j] = data[j+1] ;
    if (NumArrayEntries > 0) --NumArrayEntries ;
}

template<class EntryType>
    void List<EntryType>::Clear()
{
    NumArrayEntries = 0 ;
}


// template<class EntryType>
//     EntryType &List<EntryType>::At(
//         const int i)
// {
//     int indx = (i>=0) ? i : NumArrayEntries+i ;
// 
//     if (indx > NumAllocated-1) {
//         int new_size = int(indx * 1.5) ;
//         EntryType *tmp = new EntryType[new_size] ;
//         int j ;
//         for (j=0 ; j<NumArrayEntries ; ++j) {
//             tmp[j] = data[j] ;
//         }
//         delete [] data ;
//         data = tmp ;
//         NumAllocated = new_size ;
//     }
//     if (indx+1 > NumArrayEntries) NumArrayEntries = indx+1 ; 
//     return(data[indx]) ;
// }


template<class EntryType>
    EntryType &List<EntryType>::At(
        const int i)
{
    int indx = (i>=0) ? i : NumArrayEntries+i ;
    assert(indx < NumArrayEntries) ;
    return(data[indx]) ;
}

template<class EntryType>
    const EntryType& List<EntryType>::At(
        const int i) const
{
    assert((i >= -NumArrayEntries) && (i < NumArrayEntries)) ; 
    return((i>=0) ? data[i] : data[NumArrayEntries+i]) ;
}


// template<class EntryType>
//     EntryType &List<EntryType>::operator[] (
//         const int i)
// {
//     int indx = (i>=0) ? i : NumArrayEntries+i ;
// 
//     if (indx > NumAllocated-1) {
//         int new_size = int(indx * 1.5) ;
//         EntryType *tmp = new EntryType[new_size] ;
//         int j ;
//         for (j=0 ; j<NumArrayEntries ; ++j) {
//             tmp[j] = data[j] ;
//         }
//         delete [] data ;
//         data = tmp ;
//         NumAllocated = new_size ;
//     }
//     if (indx+1 > NumArrayEntries) NumArrayEntries = indx+1 ; 
//     return(data[indx]) ;
// }

template<class EntryType>
    EntryType& List<EntryType>::operator[] (
        const int i)
{
    assert((i >= -NumArrayEntries) && (i < NumArrayEntries)) ; 
    return((i>=0) ? data[i] : data[NumArrayEntries+i]) ;
}


template<class EntryType>
    const EntryType& List<EntryType>::operator[] (
        const int i) const
{
    assert((i >= -NumArrayEntries) && (i < NumArrayEntries)) ; 
    return((i>=0) ? data[i] : data[NumArrayEntries+i]) ;
}


template<class EntryType>
    void List<EntryType>::Reverse()
{
    // Reverse the order of the entries in the array

    for (int i=0 ; i<(NumArrayEntries/2) ; ++i) {
        int j = NumArrayEntries - i - 1 ;
        EntryType tmp = data[i] ;
        data[i] = data[j] ;
        data[j] = tmp ;
    }
}

template<class EntryType>
    void List<EntryType>::Extend(
        const List<EntryType> &other)
{
    for (int i=0 ; i<other.NumArrayEntries ; ++i)
        Append(other.data[i]) ;
}

template<class EntryType>
    void List<EntryType>::Extend(
        List<EntryType> &other)
{
    for (int i=0 ; i<other.NumArrayEntries ; ++i)
        Append(other.data[i]) ;
}

template<class EntryType>
    void List<EntryType>::ExtendOne()
{
    // Check to see if there is room.  If not the double the
    // size of the Array.

    EntryType *realloc_data = 0;

    if ((data == &StaticData[0]) &&
        (NumArrayEntries == List_STATIC_DATA_SIZE)) {
        NumAllocated = 2 * List_STATIC_DATA_SIZE ;
        EntryType *tmp = new EntryType[NumAllocated] ;
        for (int i=0 ; i<List_STATIC_DATA_SIZE ; ++i) {
            tmp[i] = data[i] ;
        }
        data = tmp ;
    //} else if (NumArrayEntries > List_STATIC_DATA_SIZE) {
    } else {
        if (NumArrayEntries == NumAllocated) {
            EntryType *tmp = new EntryType[2*NumAllocated] ;
            for (int i=0 ; i<NumArrayEntries ; ++i) {
                tmp[i] = data[i] ;
            }
            realloc_data = data ;
            //delete [] data ;
            data = tmp ;
            NumAllocated *= 2 ;
        }
    }

    // Add this entry at the end of the Array

    ++NumArrayEntries ;

    if (realloc_data != 0) delete [] realloc_data ;
}


template<class EntryType>
    void List<EntryType>::Sort(
        int(*compare)(const void*,const void*))
{
    qsort(data,NumArrayEntries,sizeof(EntryType),compare) ;
}

} // namespace

#endif
