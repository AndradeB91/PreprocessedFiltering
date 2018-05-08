//
// Dict.hpp
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2007
//   All rights reserved
//
// ------------------------------------------------------------------
//
//  A Dict object is a templated hash table.
//
//  Objects are copied when they are placed in the table.
//
//  In order to use used defined types in the table the types
//  must have define a default (no argument) constructor.  If
//  they are to be used as keys, there must also be defined an
//  equivalence operator (==) and a DictHashIndex(const mytype& key),
//  which generates an integer to be used as a hash key.
//
//  Public Interface:
//
//      Constructors:
//          Dict<type0,type1>()
//          Dict<type0,type1>(const Dict<type0,type1>& other)
//          ~Dict()
//
//      Methods:
//          d.Store(const type0& key,
//                  const type1& value) - insert an item into the dictionary
//          d.Get(const type0& key)     - returns a pointer to the value
//                                        associated with the key or zero
//                                        if key not in the dictionary
//          d.GetValue(const type0& key) - returns a reference to the value
//                                        associated with the key
//          d.HasKey(const type0& key)  - returns true if key is in the dictionary
//          d.Len()                     - returns the number of items in the dictionary
//          d.Del(const type0& key)     - removes an entry from the dictionary
//          d.Clear()                   - removes all items from the dictionary
//          d.Iterator()
//          d.ConstIterator()           - returns and iterator (see below)
//
//      Operators:
//          d0 = d1     - assignment
//          d[key]      - accesses an entry in the dictionary
//
//
//  Dictionaries can be accessed with iterators.  Typical access with an iterator
//  looks like
//
//  Dict<type0,type1>::DictIterator iter = d.Iterator()
//  for (iter.First() ; iter.More() ; ++iter) {
//      ...
//      do something
//      ...
//  }
//
//  Public Iterator Interface:
//
//      Constructor:
//          DictIterator(Dict<type0,type1>& dict)
//
//      Methods:
//          i.First()   - positions the iterator to the first item
//          i.Next()    - moves the iterator to the next item
//          i.More()    - true if the iterator points to a valid item
//          i.Key()     - returns the current key
//          i.Entry()   - returns a reference to the current value
//
//      Operators:
//          ++i, i++    - same as i.Next()
//          &i          - same as i.Key()
//          *i          - same as i.Entry()
//


#ifndef Dict_hpp
#define Dict_hpp

#ifdef MEMDEBUG
#include "MemDbg.hpp"
#endif

#include <iostream>
#include <cassert>

namespace FTools {

#define RebuildMultiplier 3

#define SMALL_HASH_TABLE 4
#define CACHE_BLOCK_SIZE 20

/// \file
///
/// \brief This file implements a dictionary (associated array) class

/// \brief A dictionary (associated array) class
///
/// This class implements a dictionary (associated array).
/// Objects are copied when they are placed in the dictionary.
/// In order to use used defined types in the dictionary the types
/// must have define a default (no argument) constructor.  If
/// they are to be used as keys, there must also be defined an
/// equivalence operator (==) and a DictHashIndex(const mytype& key)
/// function, which generates an integer to be used as a hash key.

template<class KeyType,class EntryType>
class Dict {
    public:

        class DictIterator ;
        class ConstDictIterator ;

        Dict(bool allow_duplicates = false) ;
        Dict(int num,KeyType keys[],EntryType entries[],
             bool allow_duplicates=false) ;
        Dict(const Dict &other) { Copy(other) ; }
        Dict<KeyType,EntryType>
             operator = (const Dict &other) ;
         ~Dict() ;

        /// insert a key/value pair

        bool Store(const KeyType &key,const EntryType &value) ;

        /// return a reference to a dictionary entry

        EntryType &operator[] (const KeyType &key) const ;

        /// return a pointer to an entry (zero if not present)

        EntryType *Get(const KeyType &key) const ;

        /// return a reference to a ditionary entry

        EntryType &GetValue(const KeyType &key) const ;

        /// true if this key is in the dictionary

        bool HasKey(const KeyType &key) const ;

        /// number of entries in the dictionary

        int Len() const { return(NumHashEntries) ; } ;

        /// delete an entry

        void Del(const KeyType &key) ;

        /// unused

        void DelEntry(const KeyType &key,const EntryType &value) ;

        /// delete all entries

        void Clear() ;

        /// return an iterator

        DictIterator Iterator() { return(DictIterator(*this)) ; } ;

        /// return a constant iterator

        ConstDictIterator ConstIterator() const { return(ConstDictIterator(*this)) ; } ;

        /// (Debug) print a list of keys

        void PrintKeys(std::ostream &out) const ;

        /// (Debug) print a list of entries

        void PrintEntrys(std::ostream &out) const ;

#ifdef DEBUG_PRINT
        void DebugInfo(const KeyType &key) const ;
#endif

    private:

        // static const int RebuildMultiplier ;

        struct DictEntry {
            KeyType      key ;
            EntryType    entry ;
            DictEntry *next ;
        } ;

        struct DictCache {
             struct DictCache *next ;
             struct DictEntry entries[CACHE_BLOCK_SIZE] ;
        } ;

        bool AllowDuplicates ;
        DictEntry **Buckets ;
        DictEntry *StaticBuckets[SMALL_HASH_TABLE] ;
        int NumBuckets ;
        int NumHashEntries ;
        int RebuildSize ;
        int DownShift ;
        int Mask ;

        int BucketIndex(const int i) const
            { return (int)(((((long)i)*1103515245) >>
                      DownShift) & Mask) ; } ;

        void Copy(const Dict &other) ;
        void RebuildTable() ;

        DictEntry *FreeList ;
        DictCache *CacheList ;

        DictEntry *NewHashEntry() ;
        void DeleteHashEntry(DictEntry *entry) ;

    public:

        /// \brief A dictionary Iterator

        class DictIterator {
            public:
                DictIterator() {} ;

                DictIterator(
                    Dict<KeyType,EntryType>& adict) : ptr(0)
                {
                    dict = &adict ;
                    bucket = 0 ;
                    while ((bucket < dict->NumBuckets) &&
                           (dict->Buckets[bucket] == 0)) ++bucket ;
                    if (bucket < dict->NumBuckets) ptr = dict->Buckets[bucket] ;
                } ;

                DictIterator(
                    Dict<KeyType,EntryType>* adict) : ptr(0)
                {
                    dict = adict ;
                    bucket = 0 ;
                    while ((bucket < dict->NumBuckets) &&
                           (dict->Buckets[bucket] == 0)) ++bucket ;
                    if (bucket < dict->NumBuckets) ptr = dict->Buckets[bucket] ;
                } ;

                /// set the cursor to point to the first entry

                void First() {
                    bucket = 0 ;
                    while ((bucket < dict->NumBuckets) &&
                           (dict->Buckets[bucket] == 0)) ++bucket ;
                    if (bucket < dict->NumBuckets) ptr = dict->Buckets[bucket] ;
                } ;

                /// set the cursor to point to the next entry

                void Next() {
                    ptr = ptr->next ;
                    while (ptr == 0) {
                        ++bucket ;
                        if (bucket >= dict->NumBuckets) return ;
                        ptr = dict->Buckets[bucket] ;
                    }
                } ;

                /// true if the cursor points to a valid entry

                bool More()        { return(bucket < dict->NumBuckets) ; } ;

                /// return the current key

                KeyType Key()      { return(ptr->key) ; } ;

                /// return a reference to the current entry

                EntryType &Entry() { return(ptr->entry) ; } ;

                /// return a pointer to the current entry

                EntryType *EntryPtr() { return(&(ptr->entry)) ; } ;

                /// return the number of entries in the dictionary

                int Len()          { return(dict->Len()) ; } ;

                /// set the cursor to point to the next entry

                void operator ++ ()      { Next() ; } ;

                /// set the cursor to point to the next entry

                void operator ++ (int /*i*/) { Next() ; } ;

                /// return the current key

                KeyType    operator & () { return(Key()) ; } ;

                /// return a reference to the current entry

                EntryType &operator * () { return(Entry()) ; } ;

            private:
                Dict<KeyType,EntryType>* dict ;
                int bucket ;
                typename Dict<KeyType,EntryType>::DictEntry *ptr ;
        } ;

        /// \brief A constant dictionary Iterator

        class ConstDictIterator {
            public:
                ConstDictIterator() {} ;

                ConstDictIterator(
                    const Dict<KeyType,EntryType>& adict) : ptr(0)
                {
                    dict = &adict ;
                    bucket = 0 ;
                    while ((bucket < dict->NumBuckets) &&
                           (dict->Buckets[bucket] == 0)) ++bucket ;
                    if (bucket < dict->NumBuckets) ptr = dict->Buckets[bucket] ;
                } ;

                ConstDictIterator(
                    const Dict<KeyType,EntryType>* adict) : ptr(0)
                {
                    dict = adict ;
                    bucket = 0 ;
                    while ((bucket < dict->NumBuckets) &&
                           (dict->Buckets[bucket] == 0)) ++bucket ;
                    if (bucket < dict->NumBuckets) ptr = dict->Buckets[bucket] ;
                } ;

                /// set the cursor to point to the first entry

                void First() {
                    bucket = 0 ;
                    while ((bucket < dict->NumBuckets) &&
                           (dict->Buckets[bucket] == 0)) ++bucket ;
                    if (bucket < dict->NumBuckets) ptr = dict->Buckets[bucket] ;
                } ;

                /// set the cursor to point to the next entry

                void Next() {
                    ptr = ptr->next ;
                    while (ptr == 0) {
                        ++bucket ;
                        if (bucket >= dict->NumBuckets) return ;
                        ptr = dict->Buckets[bucket] ;
                    }
                } ;

                /// true if the cursor points to a valid entry

                bool More()        { return(bucket < dict->NumBuckets) ; } ;

                /// return the current key

                KeyType Key()      { return(ptr->key) ; } ;

                /// return a constant reference to the current entry

                const EntryType &Entry() { return(ptr->entry) ; } ;

                /// return a pointer to the current entry

                const EntryType *EntryPtr() { return(&(ptr->entry)) ; } ;

                /// return the number of entries in the dictionary

                int Len()          { return(dict->Len()) ; } ;

                /// set the cursor to point to the next entry

                void operator ++ ()      { Next() ; } ;

                /// set the cursor to point to the next entry

                void operator ++ (int /*i*/) { Next() ; } ;

                /// return the current key

                KeyType   operator & () { return(Key()) ; } ;

                /// return a reference to the current entry

                const EntryType &operator * () { return(Entry()) ; } ;

            private:
                const Dict<KeyType,EntryType>* dict ;
                int bucket ;
                typename Dict<KeyType,EntryType>::DictEntry *ptr ;
        } ;

// #if defined(WIN32) || defined(WIN64)
//     friend std::ostream &operator << (std::ostream &out,
//                 const Dict<KeyType,EntryType> &tbl) ;
// #else
//     friend std::ostream &operator << <KeyType,EntryType> (std::ostream &out,
//                 const Dict<KeyType,EntryType> &tbl) ;
// #endif
} ;

// template<class KeyType,class EntryType>
// std::ostream &operator << (std::ostream &out,
//             const Dict<KeyType,EntryType> &tbl) ;

// here define a hash index routine for int types
// to keep the compiler happy (even though it will not get called)

inline int DictHashIndex(int i) { return(static_cast<int>(i)) ; }


template<class KeyType,class EntryType>
    Dict<KeyType,EntryType>::Dict(bool allow_duplicates)
{
    AllowDuplicates = allow_duplicates ;
    Buckets = StaticBuckets ;
    StaticBuckets[0] = StaticBuckets[1] = 0 ;
    StaticBuckets[2] = StaticBuckets[3] = 0 ;
    NumBuckets = SMALL_HASH_TABLE ;
    NumHashEntries = 0 ;
    RebuildSize = SMALL_HASH_TABLE*RebuildMultiplier ;
    DownShift = 28 ;
    Mask = 3 ;
    FreeList = 0 ;
    CacheList = 0 ;
}

template<class KeyType,class EntryType>
    Dict<KeyType,EntryType>::Dict(int num,KeyType keys[],
                                  EntryType entries[],bool allow_duplicates)
{
    AllowDuplicates = allow_duplicates ;
    Buckets = StaticBuckets ;
    StaticBuckets[0] = StaticBuckets[1] = 0 ;
    StaticBuckets[2] = StaticBuckets[3] = 0 ;
    NumBuckets = SMALL_HASH_TABLE ;
    NumHashEntries = 0 ;
    RebuildSize = SMALL_HASH_TABLE*RebuildMultiplier ;
    DownShift = 28 ;
    Mask = 3 ;
    FreeList = 0 ;
    CacheList = 0 ;

    for (int i=0 ; i<num ; ++i) {
        Store(keys[i],entries[i]) ;
    }
}

// template<class KeyType,class EntryType>
//     Dict<KeyType,EntryType>::Dict(const Dict &other)
// {
//     AllowDuplicates = other.AllowDuplicates ;
//     Buckets = StaticBuckets ;
//     StaticBuckets[0] = StaticBuckets[1] = 0 ;
//     StaticBuckets[2] = StaticBuckets[3] = 0 ;
//     NumBuckets = SMALL_HASH_TABLE ;
//     NumHashEntries = 0 ;
//     RebuildSize = SMALL_HASH_TABLE*RebuildMultiplier ;
//     DownShift = 28 ;
//     Mask = 3 ;
//     FreeList = 0 ;
//     CacheList = 0 ;
//
//     DictEntry *hPtr ;
//
//     for (int i=0 ; i<other.NumBuckets ; ++i) {
//         hPtr = other.Buckets[i] ;
//         while (hPtr != 0) {
//             Store(hPtr->key,hPtr->entry) ;
//             hPtr = hPtr->next ;
//         }
//     }
// }

template<class KeyType,class EntryType>
    void Dict<KeyType,EntryType>::Copy(const Dict &other)
{
    AllowDuplicates = other.AllowDuplicates ;
    if (other.Buckets == other.StaticBuckets) {
        Buckets = StaticBuckets ;
        for (int i=0 ; i < SMALL_HASH_TABLE ; i++) {
            Buckets[i] = 0 ;
        }
    } else {
        Buckets = new DictEntry*[other.NumBuckets] ;
        for (int i=0 ; i < other.NumBuckets ; i++) {
            Buckets[i] = 0 ;
        }
    }

    NumBuckets = other.NumBuckets ;
    NumHashEntries = other.NumHashEntries ;
    RebuildSize = other.RebuildSize ;
    DownShift = other.DownShift ;
    Mask = other.Mask ;

    FreeList = 0 ;
    CacheList = 0 ;

    DictEntry *hPtr,**insert ;
    for (int i=0 ; i<other.NumBuckets ; ++i) {
        hPtr = other.Buckets[i] ;
        insert = &(Buckets[i]) ;
        while (hPtr != 0) {
            // new entry
            DictEntry* nPtr = NewHashEntry() ;
            nPtr->key = hPtr->key ;
            nPtr->entry = hPtr->entry ;
            nPtr->next = 0 ;
            // insert
            *insert = nPtr ;
            insert = &(nPtr->next) ;
            hPtr = hPtr->next ;
        }
    }
}

template<class KeyType,class EntryType>
    Dict<KeyType,EntryType> Dict<KeyType,EntryType>::operator = (
        const Dict &other)
{
    DictCache *ptr = CacheList ;
    while (ptr != 0) {
        DictCache *tmp = ptr->next ;
        delete ptr ;
        ptr = tmp ;
    }

    // Free up the bucket array, if it was dynamically allocated.

    if (Buckets != 0 && Buckets != StaticBuckets)
        delete [] Buckets ;

    Copy(other) ;

//     AllowDuplicates = other.AllowDuplicates ;
//     Buckets = StaticBuckets ;
//     StaticBuckets[0] = StaticBuckets[1] = 0 ;
//     StaticBuckets[2] = StaticBuckets[3] = 0 ;
//     NumBuckets = SMALL_HASH_TABLE ;
//     NumHashEntries = 0 ;
//     RebuildSize = SMALL_HASH_TABLE*RebuildMultiplier ;
//     DownShift = 28 ;
//     Mask = 3 ;
//     FreeList = 0 ;
//     CacheList = 0 ;
//
//     DictEntry *hPtr ;
//
//     for (int i=0 ; i<other.NumBuckets ; ++i) {
//         hPtr = other.Buckets[i] ;
//         while (hPtr != 0) {
//             Store(hPtr->key,hPtr->entry) ;
//             hPtr = hPtr->next ;
//         }
//     }

    return(*this) ;
}


template<class KeyType,class EntryType>
    Dict<KeyType,EntryType>::~Dict()
{
    // delete all entries in the vtx cache

    DictCache *ptr = CacheList ;
    while (ptr != 0) {
        DictCache *tmp = ptr->next ;
        delete ptr ;
        ptr = tmp ;
    }
    CacheList = 0 ;

    // Free up the bucket array, if it was dynamically allocated.

    if (Buckets != 0 && Buckets != StaticBuckets)
        delete [] Buckets ;
}

template<class KeyType,class EntryType>
    typename Dict<KeyType,EntryType>::DictEntry
         *Dict<KeyType,EntryType>::NewHashEntry()
{
    if (FreeList == 0) {
        DictCache *ctmp = new DictCache ;
        ctmp->next = CacheList ;
        CacheList = ctmp ;
        for (int i=0 ; i<CACHE_BLOCK_SIZE-1 ; ++i) {
            ctmp->entries[i].next = &(ctmp->entries[i+1]) ;
        }
        ctmp->entries[CACHE_BLOCK_SIZE-1].next = 0 ;
        FreeList = &(ctmp->entries[0]) ;
    }
    DictEntry *tmp = FreeList ;
    FreeList = tmp->next ;
    return(tmp) ;
}

template<class KeyType,class EntryType> void
    Dict<KeyType,EntryType>::DeleteHashEntry(DictEntry *entry)
{
    entry->next = FreeList ;
    FreeList = entry ;
}

template<class KeyType,class EntryType>
    bool Dict<KeyType,EntryType>::Store(
        const KeyType &key,
        const EntryType &entry)
{
    DictEntry *hPtr ;
    int  index = BucketIndex(DictHashIndex(key)) ;

    // Search all of the entries in this bucket to see if we
    // already have this key.

    if (!AllowDuplicates) {
        for (hPtr = Buckets[index] ; hPtr != 0 ; hPtr = hPtr->next) {
            if (hPtr->key == key) {
                hPtr->entry = entry ;
                return(false) ;
            }
        }
    }

    // Entry not found.  Create a new table element and store
    // the key and the entry value.  Now link this into the
    // singly linked list rooted at this bucket.

    hPtr = NewHashEntry() ;
    hPtr->key = key ;
    hPtr->entry = entry ;
    hPtr->next = Buckets[index] ;
    Buckets[index] = hPtr ;
    NumHashEntries++ ;

    // If the table has exceeded a decent size, rebuild it
    // with many more buckets.

    if (NumHashEntries >= RebuildSize) {
        RebuildTable();
    }
    return(true) ;
}


template<class KeyType,class EntryType>
    EntryType &Dict<KeyType,EntryType>::operator[] (const KeyType &key) const
{
    int  index = BucketIndex(DictHashIndex(key)) ;

    // Search all of the entries in this bucket.  If we find the
    // key return a pointer to the entry.

    DictEntry *hPtr ;
    for (hPtr = Buckets[index] ;
         hPtr != 0 ; hPtr = hPtr->next) {
        if (hPtr->key == key) {
            return(hPtr->entry) ;
        }
    }
    assert(0) ;
    hPtr = Buckets[0] ;
    return(hPtr->entry) ;
}

template<class KeyType,class EntryType>
    EntryType &Dict<KeyType,EntryType>::GetValue(const KeyType &key) const
{
    int  index = BucketIndex(DictHashIndex(key)) ;

    // Search all of the entries in this bucket.  If we find the
    // key return a pointer to the entry.

    DictEntry *hPtr ;
    for (hPtr = Buckets[index] ;
         hPtr != 0 ; hPtr = hPtr->next) {
        if (hPtr->key == key) {
            return(hPtr->entry) ;
        }
    }
    assert(0) ;
    hPtr = Buckets[0] ;
    return(hPtr->entry) ;
}


template<class KeyType,class EntryType>
    EntryType *Dict<KeyType,EntryType>::Get(const KeyType &key) const
{
    int  index = BucketIndex(DictHashIndex(key)) ;

    // Search all of the entries in this bucket.  If we find the
    // key return a pointer to the entry.

    for (DictEntry *hPtr = Buckets[index] ;
         hPtr != 0 ; hPtr = hPtr->next) {
        if (hPtr->key == key) {
            return(&hPtr->entry) ;
        }
    }
    return(0) ;
}

template<class KeyType,class EntryType>
    bool Dict<KeyType,EntryType>::HasKey(const KeyType &key) const
{
    int  index = BucketIndex(DictHashIndex(key)) ;

    // Search all of the entries in this bucket.  If we find the
    // key return a pointer to the entry.

    for (DictEntry *hPtr = Buckets[index] ;
         hPtr != 0 ; hPtr = hPtr->next) {
        if (hPtr->key == key) {
            return(true) ;
        }
    }
    return(false) ;
}

template<class KeyType,class EntryType>
    void Dict<KeyType,EntryType>::Del(const KeyType &key)
{
    int  index = BucketIndex(DictHashIndex(key)) ;
    DictEntry *hPtr ;
    DictEntry **pPtr ;
    bool found = false ;

    // Search all of the entries in this bucket.  If we find
    // the key then delete the entry and fix up the linked list

    for (hPtr = Buckets[index], pPtr = &Buckets[index] ;
         hPtr != 0 ;
         pPtr = &hPtr->next, hPtr = hPtr->next) {
        if (hPtr->key == key) {
            *pPtr = hPtr->next ;
            DeleteHashEntry(hPtr) ;
            NumHashEntries-- ;
            found = true ;
//            if (!AllowDuplicates) break ;
            break ;
        }
    }

//    if (found) return(true) ;
    assert(found) ;
}

template<class KeyType,class EntryType>
    void Dict<KeyType,EntryType>::DelEntry(const KeyType &key,
                                           const EntryType &value)
{
    int  index = BucketIndex(DictHashIndex(key)) ;
    DictEntry *hPtr ;
    DictEntry **pPtr ;
    bool found = false ;

    // Search all of the entries in this bucket.  If we find
    // the key then delete the entry and fix up the linked list

    for (hPtr = Buckets[index], pPtr = &Buckets[index] ;
         hPtr != 0 ;
         pPtr = &hPtr->next, hPtr = hPtr->next) {
        if ((hPtr->key == key) && (hPtr->entry == value)) {
            *pPtr = hPtr->next ;
            DeleteHashEntry(hPtr) ;
            NumHashEntries-- ;
            found = true ;
//            if (!AllowDuplicates) break ;
            break ;
        }
    }

//    if (found) return(true) ;
    assert(found) ;
}

template<class KeyType,class EntryType>
    void Dict<KeyType,EntryType>::Clear()
{
    DictEntry *hPtr ;
    DictEntry *nPtr ;

    for (int i=0 ; i<NumBuckets ; ++i) {
        hPtr = Buckets[i] ;
        while (hPtr != 0) {
            nPtr = hPtr->next ;
            DeleteHashEntry(hPtr) ;
            hPtr = nPtr ;
        }
        Buckets[i] = 0 ;
    }
    NumHashEntries = 0 ;
}

template<class KeyType,class EntryType>
    void Dict<KeyType,EntryType>::RebuildTable()
{
    int oldSize = NumBuckets ;
    DictEntry **oldBuckets = Buckets ;
    DictEntry **newChainPtr ;
    int count ;

    // Allocate and initialize the new bucket array, and set up
    // hashing constants for new array size.

    NumBuckets *= 4;
    Buckets = new DictEntry*[NumBuckets] ;
    for (count = NumBuckets, newChainPtr = Buckets ;
         count > 0 ; count--, newChainPtr++) {
        *newChainPtr = 0 ;
    }
    RebuildSize *= 4 ;
    DownShift -= 2 ;
    Mask = (Mask << 2) + 3 ;

    // Rehash all of the existing entries into the new bucket array.

    for (DictEntry **oldChainPtr = oldBuckets ;
         oldSize > 0 ; oldSize--, oldChainPtr++) {
        for (DictEntry *hPtr = *oldChainPtr ;
             hPtr != 0 ; hPtr = *oldChainPtr) {
            *oldChainPtr = hPtr->next ;
            int index = BucketIndex(DictHashIndex(hPtr->key)) ;
            hPtr->next = Buckets[index] ;
            Buckets[index] = hPtr ;
        }
    }

    // Free up the old bucket array, if it was dynamically allocated.

    if (oldBuckets != StaticBuckets) {
        delete [] oldBuckets ;
    }
}

#ifdef DEBUG_PRINT
template<class KeyType,class EntryType>
    void Dict<KeyType,EntryType>::DebugInfo(const KeyType &key) const
{
    int  index = BucketIndex(DictHashIndex(key)) ;
    std::cerr << "Dict Debug: Hash Index: " << index << std::endl ;

    // Search all of the entries in this bucket.  If we find the
    // key return a pointer to the entry.

    int num = 1 ;
    for (DictEntry *hPtr = Buckets[index] ;
         hPtr != 0 ; hPtr = hPtr->next) {
        if (hPtr->key == key) {
            std::cerr << "Dict Debug: position in list " << num << std::endl ;
            return ;
        }
        num++ ;
    }
    std::cerr << "Dict Debug: Could Not Find Entry" << std::endl ;
    return ;
}
#endif

// #ifdef USE_STL_IO
//
// template<class KeyType,class EntryType>
//     void Dict<KeyType,EntryType>::PrintKeys(std::ostream &out) const
// {
//     // loop through the full hash table.
//
//     out << '{' << endl ;
//     for (int bkt=0 ; bkt < NumBuckets ; bkt++) {
//         DictEntry *chainPtr = Buckets[bkt] ;
//         while (chainPtr != 0) {
//             out << "  " << chainPtr->key << endl ;
// 	    chainPtr = chainPtr->next ;
// 	}
//     }
//     out << '}' << endl ;
// }
//
//
// template<class KeyType,class EntryType>
//     void Dict<KeyType,EntryType>::PrintEntrys(std::ostream &out) const
// {
//     // loop through the full hash table.
//
//     out << '{' << endl ;
//     for (int bkt=0 ; bkt < NumBuckets ; bkt++) {
//         DictEntry *chainPtr = Buckets[bkt] ;
//         while (chainPtr != 0) {
//             out << "  " << chainPtr->entry << endl ;
// 	    chainPtr = chainPtr->next ;
// 	}
//     }
//     out << '}' << endl ;
// }
//
// template<class KeyType,class EntryType>
//     std::ostream &operator << (std::ostream &out,
//             const Dict<KeyType,EntryType> &tbl)
// {
//     // loop through the full hash table.
//
//     out << '{' << std::endl ;
//     for (int bkt=0 ; bkt < tbl.NumBuckets ; bkt++) {
//         typename Dict<KeyType,EntryType>::DictEntry
//             *chainPtr = tbl.Buckets[bkt] ;
//         while (chainPtr != 0) {
//             out << "  " << chainPtr->key << " : " ;
//             out << chainPtr->entry << std::endl ;
// 	    chainPtr = chainPtr->next ;
// 	}
//     }
//     out << '}' << std::endl ;
//     return(out) ;
// }
// #endif

} // namespace

#endif
