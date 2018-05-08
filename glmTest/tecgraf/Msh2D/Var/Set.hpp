//
// Set
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
// ------------------------------------------------------------------
//
//  A Set object is a templated set object based on a hash table.
//  The Hash table gives efficient lookups, but has a relatively
//  high overhead.  If you need lots of small sets then look at the
//  SmallSet object.  If you need the information in the set to
//  be ordered (sorted) look at the OrderedSet object.
//
//  Objects are copied when they are placed in the Set.
//
//  In order to use used defined types in the Set the types
//  must have define a default (no argument) constructor.
//  There must also be defined an equivalence operator (==) and
//  a DictHashIndex(const mytype& key), which generates an integer
//  to be used as a hash key.
//
//  Public Interface:
//
//      Constructors:
//          Set<type>(bool allow_duplicates = false)
//          Set<type>(const Set<type>& other)
//          ~Set()
//
//      Methods:
//          s.Store(const type& value)      - insert an item into the set
//          s.HasElement(const type& value) - returns true if value is in the set
//          s.Len()                         - returns the number of items in the set
//          s.Del(const type& value)        - removes a value from the set
//          s.Clear()                       - removes all items from the set
//          s.Iterator()
//          s.ConstIterator()               - returns and iterator (see below)
//
//      Operators:
//          s0 = s1     - assignment
//
//
//  Dictionaries can be accessed with iterators.  Typical access with an iterator
//  looks like
//
//  Set<type>::SetIterator iter = s.Iterator()
//  for (iter.First() ; iter.More() ; ++iter) {
//      ...
//      do something
//      ...
//  }
//
//  Public Iterator Interface:
//
//      Constructors:
//          SetIterator(Set<type>& set)
//          SetIterator(Set<type>* set)
//
//      Methods:
//          i.First()   - positions the iterator to the first item
//          i.Next()    - moves the iterator to the next item
//          i.More()    - true if the iterator points to a valid item
//          i.Current() - returns the current value
//
//      Operators:
//          ++i, i++    - same as i.Next()
//          *i          - same as i.Current()
//

#ifndef Set_hh
#define Set_hh

#ifdef MEMDEBUG
#include "MemDbg.hpp"
#endif

#include "Dict.hpp"

namespace FTools {

/// \file
///
/// \brief This file implements a Set class using a Dictionary

/// \brief A Set class
///
/// This class implements a Set class using a Dictionary.
/// Objects are copied when they are placed in the Set.
/// In order to use used defined types in the Set the types
/// must have define a default (no argument) constructor and
/// there must also be defined an equivalence operator (==)
/// and a DictHashIndex(const mytype& key)
/// function, which generates an integer to be used as a hash key.

template<class EntryType> class Set {
    public:

        class SetIterator ;
        class ConstSetIterator ;

        Set(bool allow_duplicates = false) :
            Duplicates(allow_duplicates),
            Hash(allow_duplicates) {}

        Set(const Set &other)
        {
            Duplicates = other.Duplicates ;
            Hash = other.Hash ;
        }

        Set operator = (const Set &other)
        {
            Duplicates = other.Duplicates ;
            Hash = other.Hash ;
            return(*this) ;
        }

        ~Set() {}

        /// insert an item into the set

        void Store(const EntryType &entry)
        {
            if (!Duplicates) {
                if (Hash.HasKey(entry)) return ;
            }
            Hash.Store(entry,1) ;
        }

        /// delete an item from the set

        void Del(const EntryType &entry) { Hash.Del(entry) ; }

        /// number of items in the set

        int Len() const { return(Hash.Len()) ; }

        /// delete all items from the set

        void Clear() { Hash.Clear() ; }

        /// true if the entry is in the set

        bool HasElement(const EntryType &entry) const
        {
            return(Hash.HasKey(entry)) ;
        }

        /// return a SetIterator

        SetIterator Iterator() {
            return(SetIterator(Hash.Iterator())) ;
        }

        /// return a ConstSetIterator

        ConstSetIterator ConstIterator() const { 
            return(ConstSetIterator(Hash.ConstIterator())) ;
        }

    private:

        bool Duplicates ;
        Dict<EntryType,int> Hash ;

    public:

        /// \brief A Set Iterator

        class SetIterator {
            public:

                 SetIterator() {} ;
                 SetIterator(
                     typename Dict<EntryType,int>::DictIterator iter) : Iter(iter) {} ; 

                 /// set the cursor to point to the first entry

                 void First()             { Iter.First() ; } ;

                 /// set the cursor to point to the next entry

                 void Next()              { ++Iter ; } ;
 
                 /// true if the cursor points to a valid entry

                 bool More()              { return(Iter.More()) ; } ;

                 /// return the current set entry

                 EntryType Current()      { return(Iter.Key()) ; } ;

                 /// return the current set entry

                 EntryType operator * ()  { return(Current()) ; } ;

                 /// set the cursor to point to the next entry

                 void operator ++ ()      { Next() ; } ;

                 /// set the cursor to point to the next entry

                 void operator ++ (int i) { Next() ; } ;

            private:
                 typename Dict<EntryType,int>::DictIterator Iter ;
        } ;

        /// \brief A constant Set Iterator

        class ConstSetIterator {
            public:

                 ConstSetIterator() {} ;
                 ConstSetIterator(
                     typename Dict<EntryType,int>::ConstDictIterator iter) : Iter(iter) {} ; 


                 /// set the cursor to point to the first entry

                 void First()                   { Iter.First() ; } ;
 
                 /// set the cursor to point to the next entry

                 void Next()                    { ++Iter ; } ;
 
                 /// true if the cursor points to a valid entry

                 bool More()                    { return(Iter.More()) ; } ;

                 /// return the current set entry

                 const EntryType Current()      { return(Iter.Key()) ; } ;

                 /// return the current set entry

                 const EntryType operator * ()  { return(Current()) ; } ;

                 /// set the cursor to point to the next entry

                 void operator ++ ()            { Next() ; } ;

                 /// set the cursor to point to the next entry

                 void operator ++ (int i)       { Next() ; } ;

            private:
                 typename Dict<EntryType,int>::ConstDictIterator Iter ;
        } ;
} ;

} // namespace

#endif

