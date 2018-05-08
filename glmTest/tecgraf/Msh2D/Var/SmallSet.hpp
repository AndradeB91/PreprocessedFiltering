//
// SmallSet
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2004
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision$  $Date$  $Author$
//

#ifndef SmallSet_hh
#define SmallSet_hh

#ifdef MEMDEBUG
#include "MemDbg.hpp"
#endif

namespace FTools {

/// \file
///
/// \brief This file implements a SmallSet class

/// \brief A templated small set class
///
/// The second template argument is an integer that gives the default
/// entry slots.  The set size will be resized automatically if more
/// entries are added, but at the cost of some efficiency

template<class EntryType,const int def_size>
class SmallSet {
    public:

        class SmallSetIterator ;
        class ConstSmallSetIterator ;

        SmallSet()
        {
            NumAllocated = def_size ;
            NumElems = 0 ;
            data = &StaticData[0] ;
        } ;

        SmallSet(const SmallSet &other)
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

        SmallSet operator = (const SmallSet &other)
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

        ~SmallSet() 
        {
            if (data != &StaticData[0]) delete [] data ;
        } ;

        /// insert an item into the set

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

        /// insert an item into the set if it is not already there

        bool InsertUnique(const EntryType &entry)
        {
            for (int i=0 ; i<NumElems ; ++i) {
                if (data[i] == entry) return(false) ;
            }
            Insert(entry) ;
            return(true) ;
        } ;

        /// delete an item from the set

        void Del(const EntryType &entry)
        {
            for (int i=0 ; i<NumElems ; ++i) {
                if (data[i] == entry) {
                    for (int j=i ; j<NumElems-1 ; ++j) {
                        data[j] = data[j+1] ;
                    }
                    --NumElems ;
                    return ;
                }
            }
        } ;

        /// number of items in the set

        int Len() const { return(NumElems) ; } ;

        /// true if the entry is in the set

        bool HasElement(const EntryType &entry) const
        {
            for (int i=0 ; i<NumElems ; ++i) {
                if (data[i] == entry) return(true) ;
            }
            return(false) ;
        } ;

        /// return a reference to the item at the given index

        EntryType &Element(const int i) const { return(data[i]) ; } ;

        /// return a SmallSetIterator

        SmallSetIterator Iterator() { return(SmallSetIterator(this)) ; } ;

        /// return a ConstSmallSetIterator

        ConstSmallSetIterator ConstIterator() { 
            return(ConstSmallSetIterator(this)) ;
        } ;

        /// return the intersection of two sets

        SmallSet Intersect(const SmallSet &other) {
            SmallSet tmp ;
            for (int i=0 ; i<NumElems ; ++i) {
                if (other.HasElement(data[i])) tmp.Insert(data[i]) ;
            }
            return tmp ;
        }

        /// return the union of two sets

        SmallSet Union(const SmallSet &other) {
            SmallSet tmp ;
            for (int i=0 ; i<NumElems ; ++i) {
                tmp.Insert(data[i]) ;
            }
            for (int i=0 ; i<other.NumElems ; ++i) {
                tmp.InsertUnique(other.data[i]) ;
            }
            return tmp ;
        }

    private:

        EntryType StaticData[def_size], *data ;
        int NumElems ;
        int NumAllocated ;

    public:

        class SmallSetIterator {
            public:

                 SmallSetIterator(
                     SmallSet<EntryType,def_size> * aSet) :
                         set(aSet),cur(0) {} ;

                 /// set the cursor to point to the first entry

                 void First()        { cur = 0 ; } ;

                 /// set the cursor to point to the next entry

                 void Next()         { ++cur ; } ;
 
                 /// true if the cursor points to a valid entry

                 bool More()         { return(cur < set->NumElems) ; } ;

                 /// return the current set entry

                 EntryType Current() { return(cur < set->NumElems ?
                                       set->data[cur] : 0) ; } ;

                 /// return the current set entry

                 EntryType operator * ()  { return(Current()) ; } ;

                 /// set the cursor to point to the next entry

                 void operator ++ ()      { Next() ; } ;

                 /// set the cursor to point to the next entry

                 void operator ++ (int i) { Next() ; } ;

            private:
                SmallSet<EntryType,def_size> *set ;
                int cur ;
        } ;

        class ConstSmallSetIterator {
            public:

                 ConstSmallSetIterator(
                     const SmallSet<EntryType,def_size> * aSet) :
                         set(aSet),cur(0) {} ;

                 /// set the cursor to point to the first entry

                 void First()        { cur = 0 ; } ;

                 /// set the cursor to point to the next entry

                 void Next()         { ++cur ; } ;
 
                 /// true if the cursor points to a valid entry

                 bool More()         { return(cur < set->NumElems) ; } ;

                 /// return the current set entry

                 EntryType Current() { return(cur < set->NumElems ?
                                       set->data[cur] : 0) ; } ;

                 /// return the current set entry

                 EntryType operator * ()  { return(Current()) ; } ;

                 /// set the cursor to point to the next entry

                 void operator ++ ()      { Next() ; } ;

                 /// set the cursor to point to the next entry

                 void operator ++ (int i) { Next() ; } ;

            private:
                const SmallSet<EntryType,def_size> *set ;
                int cur ;
        } ;


// #ifdef USE_STL_IO
// #if defined(WIN32) || defined(WIN64)
//     friend std::ostream &operator << (std::ostream &out,
//                 const CArbSmallSet<EntryType,def_size> &set) ;
// #else
//     friend std::ostream &operator << <EntryType,def_size> (std::ostream &out,
//                 const CArbSmallSet<EntryType,def_size> &set) ;
// #endif
// #endif
} ;


// #ifdef USE_STL_IO
// #include <iostream>
// 
// template<class EntryType,const int def_size>
// inline std::ostream &operator << (std::ostream &out,
//                     const CArbSmallSet<EntryType,def_size> &set)
// {
//     out << "[" ;
//     for (int i=0 ; i<set.NumElems ; ++i) {
//         if (i > 0) out << ' ' ;
//         out << set.data[i] ;
//     }
//     out << ']' ;
//     return(out) ;
// }
// #endif

} // namespace

#endif

