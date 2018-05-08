//
// CArbHashTable Template Class header file
//
// Description -
//   This class implements a templated hash table object.
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 1999,2000,2001,2002
//   All rights reserved
//
// Author -
//   Wash Wawrzynek
//
// Revision -
//   $Revision: 1.9 $  $Date: 2002/06/11 14:36:08 $  $Author: wash $
//

#ifndef ArbHashTable_hh
#define ArbHashTable_hh

#ifdef USE_STL_IO
#include <iostream.h>
#endif

#define RebuildMultiplier 3

template <class KeyType, class EntryType> class CArbHashTableIterator;
template <class KeyType, class EntryType> class CArbConstHashTableIterator;

#define SMALL_HASH_TABLE 4
#define CACHE_BLOCK_SIZE 20
 
template <class KeyType, class EntryType>
class CArbHashTable
{
    public:

        CArbHashTable  (bool allow_duplicates = false);
        CArbHashTable  (const CArbHashTable &other);
        CArbHashTable  <KeyType, EntryType> operator = (const CArbHashTable &other);
        ~CArbHashTable ( ) ;

        bool Store           (const KeyType &key,const EntryType &value);
        EntryType *Fetch     (const KeyType &key) const;
        EntryType FetchValue (const KeyType &key) const;
        bool Remove          (const KeyType &key);
        bool RemoveEntry     (const KeyType &key,const EntryType &value);
        bool FetchAny        (KeyType **key,EntryType **value) const;
        void Clear           ( ) ;

        KeyType   *GetKeyList    ( ) const;
        EntryType **GetEntryList ( ) const;

        int NumEntries ( ) const { return(NumHashEntries) ; };

#ifdef USE_STL_IO
        void PrintKeys(ostream &out) const ;
        void PrintEntrys(ostream &out) const ;
#endif

    private:

        // static const int RebuildMultiplier ;

        struct ArbHashEntry
        {
            KeyType      key;
            EntryType    entry;
            ArbHashEntry *next;
        };

        struct ArbHashCache
        {
             struct ArbHashCache *next ;
             struct ArbHashEntry entries[CACHE_BLOCK_SIZE] ;
        };
 
        bool          AllowDuplicates ;
        ArbHashEntry  **ppBuckets ;
        ArbHashEntry  *pStaticBuckets[SMALL_HASH_TABLE] ;
        int           NumBuckets ;
        int           NumHashEntries ;
        int           RebuildSize ;
        int           DownShift ;
        int           Mask ;

        int BucketIndex(const int i) const
            { return (((((long)i)*1103515245) >> DownShift) & Mask) ; };

        void RebuildTable ( );

        ArbHashEntry *FreeList ;
        ArbHashCache *CacheList ;

        ArbHashEntry  *NewHashEntry    ( );
        void           DeleteHashEntry (ArbHashEntry *entry);

    friend class CArbHashTableIterator      <KeyType, EntryType>;
    friend class CArbConstHashTableIterator <KeyType, EntryType>;

#ifdef USE_STL_IO
#if defined(WIN32) || defined(WIN64)
    friend ostream &operator << 
      (ostream &out, const CArbHashTable <KeyType, EntryType> &tbl);
#else
    friend ostream &operator << <KeyType, EntryType> 
      (ostream &out, const CArbHashTable <KeyType, EntryType> &tbl);
#endif
#endif
};  // class CArbHashTable


#ifdef USE_STL_IO
template <class KeyType, class EntryType>
ostream &operator << (ostream &out, const CArbHashTable <KeyType, EntryType> &tbl);
#endif


template <class KeyType, class EntryType>
class CArbHashTableIterator
{
    public:

        CArbHashTableIterator (const CArbHashTable <KeyType, EntryType> *const aTable);
        void       First ( );
        void       Next  ( );
        bool       More  ( ) { return(bucket < table->NumBuckets) ; };
        KeyType    Key   ( ) { return(ptr->key) ; };
        EntryType *Entry ( ) { return(&ptr->entry) ; }; 

        void operator ++ ( )     { Next( ) ; };
        void operator ++ (int ) { Next( ) ; };

    private:

        const CArbHashTable <KeyType, EntryType> *table;
        int   bucket ;
        typename CArbHashTable <KeyType, EntryType>::ArbHashEntry *ptr;
};



template <class KeyType, class EntryType>
class CArbConstHashTableIterator
{
    public:

        CArbConstHashTableIterator (const CArbHashTable <KeyType, EntryType> *const aTable);
        void       First ( );
        void       Next  ( );
        bool       More  ( )  { return(bucket < table->NumBuckets) ; };
        KeyType    Key   ( )  { return(ptr->key) ; };
        EntryType *Entry ( )  { return(&ptr->entry) ; }; 

        void operator ++ ()      { Next() ; } ;
        void operator ++ (int i) { Next() ; } ;

    private:

        const CArbHashTable <KeyType, EntryType> *table;
        int   bucket;
        typename CArbHashTable <KeyType, EntryType>::ArbHashEntry *ptr;
} ;



#ifdef MEMDEBUG
#include "MemDbg.hpp"
#define new new(__FILE__,__LINE__)
#endif

//#include "ArbTopo3D.hpp"

// force templated code to be generated

// here define a hash index routine for int types
// to keep the compiler happy (even though it will not get called)

inline int ArbHashIndex(int i) { return(static_cast<int>(i)) ; }


// %(CArbHashTable::CArbHashTable-constructor-|-bool-|) 
/* ++ ----------------------------------------------------------
**
**    CArbHashTable - hash table constructor
**
**      CArbHashTable(bool allow_duplicates = false)
**
**        allow_duplicates - (in)  If this flag is set, entries with 
**                                 duplicate keys are allowed. 
**
**      Description: This is the constructor for a hash table. 
**
**
** -- */


template<class KeyType,class EntryType>
CArbHashTable<KeyType,EntryType>::CArbHashTable(bool allow_duplicates)
{
    AllowDuplicates = allow_duplicates ;
    ppBuckets = pStaticBuckets ;
    pStaticBuckets[0] = pStaticBuckets[1] = 0 ;
    pStaticBuckets[2] = pStaticBuckets[3] = 0 ;
    NumBuckets = SMALL_HASH_TABLE ;
    NumHashEntries = 0 ;
    RebuildSize = SMALL_HASH_TABLE*RebuildMultiplier ;
    DownShift = 28 ;
    Mask = 3 ;
    FreeList = 0 ;
    CacheList = 0 ;
}

template<class KeyType,class EntryType>
CArbHashTable<KeyType,EntryType>::
CArbHashTable(const CArbHashTable &other)
{
    AllowDuplicates = other.AllowDuplicates ;
    ppBuckets = pStaticBuckets ;
    pStaticBuckets[0] = pStaticBuckets[1] = 0 ;
    pStaticBuckets[2] = pStaticBuckets[3] = 0 ;
    NumBuckets = SMALL_HASH_TABLE ;
    NumHashEntries = 0 ;
    RebuildSize = SMALL_HASH_TABLE*RebuildMultiplier ;
    DownShift = 28 ;
    Mask = 3 ;
    FreeList = 0 ;
    CacheList = 0 ;

    ArbHashEntry *hPtr ;

    for (int i=0 ; i<other.NumBuckets ; ++i) {
        hPtr = other.ppBuckets[i] ;
        while (hPtr != 0) {
            Store(hPtr->key,hPtr->entry) ;
            hPtr = hPtr->next ;
        }
    }
}


template<class KeyType,class EntryType>
CArbHashTable<KeyType,EntryType> 
CArbHashTable<KeyType,EntryType>::
operator = (const CArbHashTable &other)
{
    ArbHashCache *ptr = CacheList ;
    while (ptr != 0) {
        ArbHashCache *tmp = ptr->next ;
        delete ptr ;
        ptr = tmp ;
    }

    // Free up the bucket array, if it was dynamically allocated.

    if (ppBuckets != pStaticBuckets) delete [] ppBuckets ;

    AllowDuplicates = other.AllowDuplicates ;
    ppBuckets = pStaticBuckets ;
    pStaticBuckets[0] = pStaticBuckets[1] = 0 ;
    pStaticBuckets[2] = pStaticBuckets[3] = 0 ;
    NumBuckets = SMALL_HASH_TABLE ;
    NumHashEntries = 0 ;
    RebuildSize = SMALL_HASH_TABLE*RebuildMultiplier ;
    DownShift = 28 ;
    Mask = 3 ;
    FreeList = 0 ;
    CacheList = 0 ;

    ArbHashEntry *hPtr ;

    for (int i=0 ; i<other.NumBuckets ; ++i) {
        hPtr = other.ppBuckets[i] ;
        while (hPtr != 0) {
            Store(hPtr->key,hPtr->entry) ;
            hPtr = hPtr->next ;
        }
    }

    return(*this) ;
}



// %(CArbHashTable::CArbHashTable-destructor-|~) 
/* ++ ----------------------------------------------------------
**
**    CArbHashTable - hash table destructor
**
**      ~CArbHashTable()
**
**      Description: This is a destructor for a hash table. 
**
**
** -- */


template<class KeyType,class EntryType>
CArbHashTable<KeyType,EntryType>::~CArbHashTable()
{
//    ArbHashEntry *hPtr, *nextPtr;

    // Delete all entries in the table

//     for (int i = 0; i < NumBuckets; i++) {
// 	hPtr = ppBuckets[i] ;
// 	while (hPtr != 0) {
// 	    nextPtr = hPtr->next ;
// 	    delete hPtr ;
// 	    hPtr = nextPtr ;
// 	}
//     }

    // delete all entries in the vtx cache

    ArbHashCache *ptr = CacheList ;
    while (ptr != 0) {
        ArbHashCache *tmp = ptr->next ;
        delete ptr ;
        ptr = tmp ;
    }

    // Free up the bucket array, if it was dynamically allocated.

    if (ppBuckets != pStaticBuckets) delete [] ppBuckets ;
}



template<class KeyType,class EntryType>
typename CArbHashTable<KeyType,EntryType>::ArbHashEntry
*CArbHashTable<KeyType,EntryType>::NewHashEntry()
{
    if (FreeList == 0) {
        ArbHashCache *ctmp = new ArbHashCache ;
        ctmp->next = CacheList ;
        CacheList = ctmp ;
        for (int i=0 ; i<CACHE_BLOCK_SIZE-1 ; ++i) {
            ctmp->entries[i].next = &(ctmp->entries[i+1]) ;
        }
        ctmp->entries[CACHE_BLOCK_SIZE-1].next = 0 ;
        FreeList = &(ctmp->entries[0]) ;
    }
    ArbHashEntry *tmp = FreeList ;
    FreeList = tmp->next ;
    return(tmp) ;
}

template<class KeyType,class EntryType> void
CArbHashTable<KeyType,EntryType>::DeleteHashEntry(ArbHashEntry *entry)
{
    entry->next = FreeList ;
    FreeList = entry ;
}

// /* -----------------------------------------------------------
//     RebuildMultiplier - static member variable.  When we grow
//                         the hash table we will be this much bigger.
// */
// 
// template<class KeyType,class EntryType>
// const int CArbHashTable<KeyType,EntryType>::RebuildMultiplier = 3 ;


// %(CArbHashTable::Store-bool-|-KeyType-const|&-EntryType-const|&)
/* ++ ----------------------------------------------------------
**
**    Store - store an entry in the hash table
**
**      bool Store(
**              const KeyType   &key,
**              const EntryType &value)
**
**        key   - (in)  Key associated with this entry in the hash table. 
**                      Keys must be unique for all entries. 
**        value - (in)  Entry to be added the hash table. A copy is make 
**                      of the entry when it is placed in the table. 
**
**      Description: This function stores a new entry in the table. 
**
**      Return Value: This function returns true if the entry is stored 
**          in the hash table succesfully. A false status indicates that 
**          this key is already in the table and the duplicates flag was 
**          not set at construction. 
**
**
** -- */


template<class KeyType,class EntryType>
bool CArbHashTable<KeyType,EntryType>::
Store(const KeyType &key,const EntryType &entry)
{
    ArbHashEntry *hPtr ;
    int  index = BucketIndex(ArbHashIndex(key)) ;

    // Search all of the entries in this bucket to see if we
    // already have this key.
    if (!AllowDuplicates)
    {
      for (hPtr = ppBuckets[index] ; hPtr != 0 ; hPtr = hPtr->next)
      {
        if (hPtr->key == key)
        {
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
    hPtr->next = ppBuckets[index] ;
    ppBuckets[index] = hPtr ;
    NumHashEntries++ ;

    // If the table has exceeded a decent size, rebuild it
    // with many more buckets.

    if (NumHashEntries >= RebuildSize)
      RebuildTable();

    return(true) ;
}




// %(CArbHashTable::Fetch-EntryType-|*^const-KeyType-const|&) 
/* ++ ----------------------------------------------------------
**
**    Fetch - retrieve an entry from the table
**
**      EntryType *Fetch(const KeyType &key) const
**
**        key - (in)  search key 
**
**      Description: Given a key, fetch the associated entry from the 
**          hash table. 
**
**      Return Value: This function returns a pointer to the table entry 
**          associated with the givn key. A null pointer is returned if 
**          this key is not found in the table. 
**
**
** -- */


template<class KeyType,class EntryType>
EntryType *CArbHashTable<KeyType,EntryType>::
Fetch(const KeyType &key) const
{
    int  index = BucketIndex(ArbHashIndex(key)) ;

    // Search all of the entries in this bucket.  If we find the
    // key return a pointer to the entry.

    for (ArbHashEntry *hPtr = ppBuckets[index] ;
         hPtr != 0 ; hPtr = hPtr->next) {
	if (hPtr->key == key) {
            return(&hPtr->entry) ;
        }
    }
    return(0) ;
}


template<class KeyType,class EntryType>
EntryType CArbHashTable<KeyType,EntryType>::
FetchValue(const KeyType &key) const
{
    int  index = BucketIndex(ArbHashIndex(key)) ;

    // Search all of the entries in this bucket.  If we find the
    // key return a pointer to the entry.

    for (ArbHashEntry *hPtr = ppBuckets[index] ;
         hPtr != 0 ; hPtr = hPtr->next) {
	if (hPtr->key == key) {
            return(hPtr->entry) ;
        }
    }
    return(0) ;
}


// %(CArbHashTable::Remove-bool-|-KeyType-const|&)
/* ++ ----------------------------------------------------------
**
**    Remove - remove an entry from the table
**
**      bool Remove(const KeyType &key)
**
**        key - (in)  search key 
**
**      Description: Given a key, this function removes the key and 
**          associated entry from the table. If the duplicates flag is 
**          set and there are multiple entries associated with the search 
**          key, only one entry (selected abitrarily) is deleted. 
**
**      Return Value: This function returns true if the entry is found 
**          and removed from the table. A false status indicates that 
**          this key was not found in the table. 
**
**
** -- */


template<class KeyType,class EntryType>
bool CArbHashTable<KeyType,EntryType>::
Remove(const KeyType &key)
{
    int  index = BucketIndex(ArbHashIndex(key)) ;
    ArbHashEntry *hPtr ;
    ArbHashEntry **pPtr ;
    bool found = false ;

    // Search all of the entries in this bucket.  If we find
    // the key then delete the entry and fix up the linked list

    for (hPtr = ppBuckets[index], pPtr = &ppBuckets[index] ;
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

    if (found)
        return(true) ;
    else
        return(false) ;
}




// %(CArbHashTable::RemoveEntry-bool-|-KeyType-const|&-EntryType-const|&) 
/* ++ ----------------------------------------------------------
**
**    RemoveEntry - remove an entry from the table
**
**      bool RemoveEntry(
**              const KeyType   &key,
**              const EntryType &value)
**
**        key   - (in)  search key 
**        value - (in)  search value 
**
**      Description: Given both a key and an entry value, this function 
**          removes the corresponding entry from the table. In practice, 
**          this function is useful for hash tables for which duplicates 
**          are allowed. 
**
**      Return Value: This function returns true if the entry is found 
**          and removed from the table. A false status indicates that 
**          this entry was not found in the table. 
**
**
** -- */


template<class KeyType,class EntryType>
bool CArbHashTable<KeyType,EntryType>::
RemoveEntry(const KeyType &key,const EntryType &entry)
{
    int  index = BucketIndex(ArbHashIndex(key)) ;
    ArbHashEntry *hPtr ;
    ArbHashEntry **pPtr ;

    // Search all of the entries in this bucket.  If we find
    // the key then delete the entry and fix up the linked list

    for (hPtr = ppBuckets[index], pPtr = &ppBuckets[index] ;
         hPtr != 0 ;
         pPtr = &hPtr->next, hPtr = hPtr->next) {
	if ((hPtr->key == key) && (hPtr->entry == entry)) {
            *pPtr = hPtr->next ;
            DeleteHashEntry(hPtr) ;
            NumHashEntries-- ;
            return(true) ;
        }
    }
    return(false) ;
}




// %(CArbHashTable::FetchAny-bool-|^const-KeyType-|**-EntryType-|**)
/* ++ ----------------------------------------------------------
**
**    FetchAny - return an entry from the table
**                      selected arbitrarily
**
**      bool FetchAny(
**              KeyType   **key,
**              EntryType **value) const
**
**        key   - (out) Address of a pointer to a key type. This is set 
**                      to point to the key of the returned entry. 
**        value - (out) Address of a pointer to an entry type. This is 
**                      set to point to the returned entry. 
**
**      Description: This function returns one entry fron the table. The 
**          entry is selected abitrarily. 
**
**      Return Value: This function returns true if an entry was found. A 
**          false status indicates that the hash table is empty. 
**
**
** -- */


template<class KeyType,class EntryType>
bool CArbHashTable<KeyType,EntryType>::
FetchAny(KeyType **key,EntryType **entry) const
{
    // return first the entry in the table.

    for (int i=0 ; i<NumBuckets ; ++i) {
        if (ppBuckets[i] != 0) {
            *key = &(ppBuckets[i]->key) ;
            *entry = &(ppBuckets[i]->entry) ;
            return(true) ;
        }
    }
    return(false) ;
}


template<class KeyType,class EntryType>
void CArbHashTable<KeyType,EntryType>::Clear()
{
    ArbHashEntry *hPtr ;
    ArbHashEntry *nPtr ;

    for (int i=0 ; i<NumBuckets ; ++i) {
        hPtr = ppBuckets[i] ;
        while (hPtr != 0) {
            nPtr = hPtr->next ;
            DeleteHashEntry(hPtr) ;
            hPtr = nPtr ;
        }
        ppBuckets[i] = 0 ;
    }
    NumHashEntries = 0 ;
}


// %(CArbHashTable::GetKeyList-KeyType-|*^const) 
/* ++ ----------------------------------------------------------
**
**    GetKeyList - gets a list of all the keys in the table
**
**      KeyType *GetKeyList() const
**
**      Description: This function returns a pointer to a list (array) of 
**          all the keys in the table. The caller takes ownership of the 
**          pointer and must call delete []. 
**
**      Return Value: A pointer to a list of the keys in the table. 
**
**
** -- */


template<class KeyType,class EntryType>
KeyType *CArbHashTable<KeyType,EntryType>::GetKeyList() const
{
    int cur ;
    KeyType *keys = new KeyType[NumHashEntries] ;

    // loop through the full hash table.

    cur = 0 ;
    for (int bkt=0 ; bkt < NumBuckets ; bkt++) {
        ArbHashEntry *chainPtr = ppBuckets[bkt] ; 
        while (chainPtr != 0) {
            keys[cur] = chainPtr->key ;
            ++cur ;
	    chainPtr = chainPtr->next ;
	}
    }
    return(keys) ;
}




// %(CArbHashTable::GetEntryList-EntryType-|**^const) 
/* ++ ----------------------------------------------------------
**
**    GetEntryList - gets a list of all values in the table
**
**      EntryType **GetEntryList() const
**
**      Description: This function returns a pointer to a list (array) of 
**          pointer to all the values in the table. The caller takes 
**          ownership of the pointer and must call delete []. 
**
**      Return Value: A pointer to a list of the values in the table. 
**
**
** -- */


template<class KeyType,class EntryType>
EntryType **CArbHashTable<KeyType,EntryType>::GetEntryList() const
{
    int cur ;
    EntryType **entries = new EntryType*[NumHashEntries] ;

    // loop through the full hash table.

    cur = 0 ;
    for (int bkt=0 ; bkt < NumBuckets ; bkt++) {
        ArbHashEntry *chainPtr = ppBuckets[bkt] ; 
        while (chainPtr != 0) {
            entries[cur] = &chainPtr->entry ;
            ++cur ;
	    chainPtr = chainPtr->next ;
	}
    }
    return(entries) ;
}




// %(CArbHashTable::RebuildTable-void-|) 
/* ++ ----------------------------------------------------------
**
**    RebuildTable - rebuild the hash table to increase its size
**
**      void RebuildTable()
**
**      Description: This function increases the number of buckets by the 
**          factor RebuildMultiplier and rehashes all table entries. 
**
**
** -- */


template<class KeyType,class EntryType>
void CArbHashTable<KeyType,EntryType>::RebuildTable()
{
    int oldSize = NumBuckets ;
    ArbHashEntry **oldBuckets = ppBuckets ;
    ArbHashEntry **newChainPtr ;
    int count ;

    // Allocate and initialize the new bucket array, and set up
    // hashing constants for new array size.

    NumBuckets *= 4;
    ppBuckets = new ArbHashEntry*[NumBuckets] ;
    for (count = NumBuckets, newChainPtr = ppBuckets ;
	 count > 0 ; count--, newChainPtr++) {
	*newChainPtr = 0 ;
    }
    RebuildSize *= 4 ;
    DownShift -= 2 ;
    Mask = (Mask << 2) + 3 ;

    // Rehash all of the existing entries into the new bucket array.

    for (ArbHashEntry **oldChainPtr = oldBuckets ;
         oldSize > 0 ; oldSize--, oldChainPtr++) {
        for (ArbHashEntry *hPtr = *oldChainPtr ;
             hPtr != 0 ; hPtr = *oldChainPtr) {
	    *oldChainPtr = hPtr->next ;
            int index = BucketIndex(ArbHashIndex(hPtr->key)) ;
	    hPtr->next = ppBuckets[index] ;
	    ppBuckets[index] = hPtr ;
	}
    }

    // Free up the old bucket array, if it was dynamically allocated.

    if (oldBuckets != pStaticBuckets) {
	delete [] oldBuckets ;
    }
}

#ifdef USE_STL_IO

template<class KeyType,class EntryType>
void CArbHashTable<KeyType,EntryType>::PrintKeys(ostream &out) const
{
    // loop through the full hash table.

    out << '{' << endl ;
    for (int bkt=0 ; bkt < NumBuckets ; bkt++) {
        ArbHashEntry *chainPtr = ppBuckets[bkt] ; 
        while (chainPtr != 0) {
            out << "  " << chainPtr->key << endl ;
	    chainPtr = chainPtr->next ;
	}
    }
    out << '}' << endl ;
}


template<class KeyType,class EntryType>
void CArbHashTable<KeyType,EntryType>::PrintEntrys(ostream &out) const
{
    // loop through the full hash table.

    out << '{' << endl ;
    for (int bkt=0 ; bkt < NumBuckets ; bkt++) {
        ArbHashEntry *chainPtr = ppBuckets[bkt] ; 
        while (chainPtr != 0) {
            out << "  " << chainPtr->entry << endl ;
	    chainPtr = chainPtr->next ;
	}
    }
    out << '}' << endl ;
}

template<class KeyType,class EntryType>
ostream &operator << (ostream &out,
            const CArbHashTable<KeyType,EntryType> &tbl)
{
    // loop through the full hash table.

    out << '{' << endl ;
    for (int bkt=0 ; bkt < tbl.NumBuckets ; bkt++) {
        CArbHashTable<KeyType,EntryType>::ArbHashEntry
            *chainPtr = tbl.ppBuckets[bkt] ; 
        while (chainPtr != 0) {
            out << "  " << chainPtr->key << " : " ;
            out << chainPtr->entry << endl ;
	    chainPtr = chainPtr->next ;
	}
    }
    out << '}' << endl ;
    return(out) ;
}
#endif


// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------

// %(CArbHashTableIterator::CArbHashTableIterator-constructor-|-CArbHashTable-const|<KeyType,EntryType>*) 
/* ++ ----------------------------------------------------------
**
**    CArbHashTableIterator - constructor for hash table iterators
**
**      CArbHashTableIterator(const CArbHashTable<KeyType,EntryType>* aTable)
**
**        aTable - (in)  address of the hash table to be visited 
**
**      Description: This is a constructor for HashTableIterator objects. 
**
**
** -- */


template<class KeyType,class EntryType>
CArbHashTableIterator<KeyType,EntryType>::
CArbHashTableIterator(const CArbHashTable<KeyType,EntryType> *aTable)
{
    table = aTable ;
    bucket = 0 ;
    while ((bucket < table->NumBuckets) && 
           (table->ppBuckets[bucket] == 0)) ++bucket ;
    if (bucket < table->NumBuckets) ptr = table->ppBuckets[bucket] ;
}




// %(CArbHashTableIterator::First-void-|) 
/* ++ ----------------------------------------------------------
**
**    First - point to the first element
**
**      void First()
**
**      Description: This function sets the iterator to point to the 
**          first entry in the table. 
**
**
** -- */


template<class KeyType,class EntryType>
void CArbHashTableIterator<KeyType,EntryType>::First()
{
    bucket = 0 ;
    while ((bucket < table->NumBuckets) && 
           (table->ppBuckets[bucket] == 0)) ++bucket ;
    if (bucket < table->NumBuckets) ptr = table->ppBuckets[bucket] ;
}



// %(CArbHashTableIterator::Next-void-|) 
/* ++ ----------------------------------------------------------
**
**    Next - move to the next item in the table
**
**      void Next()
**
**      Description: This function sets the iterator to point to the next 
**          entry in the table. 
**
**
** -- */


template<class KeyType,class EntryType>
void CArbHashTableIterator<KeyType,EntryType>::Next()
{
    ptr = ptr->next ;
    while (ptr == 0) {
        ++bucket ;
        if (bucket >= table->NumBuckets) return ;
        ptr = table->ppBuckets[bucket] ;
    }
}



template<class KeyType,class EntryType>
CArbConstHashTableIterator<KeyType,EntryType>::
CArbConstHashTableIterator(const CArbHashTable<KeyType,EntryType> *aTable)
{
    table = aTable ;
    bucket = 0 ;
    while ((bucket < table->NumBuckets) && 
           (table->ppBuckets[bucket] == 0)) ++bucket ;
    if (bucket < table->NumBuckets) ptr = table->ppBuckets[bucket] ;
}


template<class KeyType,class EntryType>
void CArbConstHashTableIterator<KeyType,EntryType>::First()
{
    bucket = 0 ;
    while ((bucket < table->NumBuckets) && 
           (table->ppBuckets[bucket] == 0)) ++bucket ;
    if (bucket < table->NumBuckets) ptr = table->ppBuckets[bucket] ;
}


template<class KeyType,class EntryType>
void CArbConstHashTableIterator<KeyType,EntryType>::Next()
{
    ptr = ptr->next ;
    while (ptr == 0) {
        ++bucket ;
        if (bucket >= table->NumBuckets) return ;
        ptr = table->ppBuckets[bucket] ;
    }
}


/*
TEMPLATE CLASS CArbHashTable<class KeyType,class EntryType>

  This is a general hash table template object. It maps a key to a 
  value stored in the table. Copies are made of both the key and value 
  and are stored in the table. Therefore, if the values are large and 
  stored elseware, a reasonable strategy is to build a hash table of 
  pointers to these values. 

  The caller must define a function of the form: 

  int ArbHashIndex(KeyType i) 

  that transforms a key of abitrary type to an integer. A hash index 
  function for integers is already built in. 


  Template Arguments:

    KeyType - key type

    EntryType - value type


PUBLIC INTERFACE

  Public Member Functions:

    CArbHashTable - hash table constructor 

      CArbHashTable(bool allow_duplicates = false)

        allow_duplicates - (in)  If this flag is set, entries with 
                                 duplicate keys are allowed. 

      Description: This is the constructor for a hash table. 


    CArbHashTable - hash table destructor 

      ~CArbHashTable()

      Description: This is a destructor for a hash table. 


    Store - store an entry in the hash table 

      bool Store(
              const KeyType   &key,
              const EntryType &value)

        key   - (in)  Key associated with this entry in the hash 
                      table. Keys must be unique for all entries. 
        value - (in)  Entry to be added the hash table. A copy is 
                      make of the entry when it is placed in the 
                      table. 

      Description: This function stores a new entry in the table. 

      Return Value: This function returns true if the entry is stored 
          in the hash table succesfully. A false status indicates 
          that this key is already in the table and the duplicates 
          flag was not set at construction. 


    Fetch - retrieve an entry from the table 

      EntryType *Fetch(const KeyType &key) const

        key - (in)  search key 

      Description: Given a key, fetch the associated entry from the 
          hash table. 

      Return Value: This function returns a pointer to the table 
          entry associated with the givn key. A null pointer is 
          returned if this key is not found in the table. 


    Remove - remove an entry from the table 

      bool Remove(const KeyType &key)

        key - (in)  search key 

      Description: Given a key, this function removes the key and 
          associated entry from the table. If the duplicates flag is 
          set and there are multiple entries associated with the 
          search key, only one entry (selected abitrarily) is 
          deleted. 

      Return Value: This function returns true if the entry is found 
          and removed from the table. A false status indicates that 
          this key was not found in the table. 


    RemoveEntry - remove an entry from the table 

      bool RemoveEntry(
              const KeyType   &key,
              const EntryType &value)

        key   - (in)  search key 
        value - (in)  search value 

      Description: Given both a key and an entry value, this function 
          removes the corresponding entry from the table. In 
          practice, this function is useful for hash tables for which 
          duplicates are allowed. 

      Return Value: This function returns true if the entry is found 
          and removed from the table. A false status indicates that 
          this entry was not found in the table. 


    FetchAny - return an entry from the table selected arbitrarily 

      bool FetchAny(
              KeyType   **key,
              EntryType **value) const

        key   - (out) Address of a pointer to a key type. This is set 
                      to point to the key of the returned entry. 
        value - (out) Address of a pointer to an entry type. This is 
                      set to point to the returned entry. 

      Description: This function returns one entry fron the table. 
          The entry is selected abitrarily. 

      Return Value: This function returns true if an entry was found. 
          A false status indicates that the hash table is empty. 


    GetKeyList - gets a list of all the keys in the table 

      KeyType *GetKeyList() const

      Description: This function returns a pointer to a list (array) 
          of all the keys in the table. The caller takes ownership of 
          the pointer and must call delete []. 

      Return Value: A pointer to a list of the keys in the table. 


    GetEntryList - gets a list of all values in the table 

      EntryType **GetEntryList() const

      Description: This function returns a pointer to a list (array) 
          of pointer to all the values in the table. The caller takes 
          ownership of the pointer and must call delete []. 

      Return Value: A pointer to a list of the values in the table. 


    NumEntries - number of entries in the table 

      int NumEntries()

      Description: This function returns the number of entries in the 
          table. 

      Return Value: number of entries in the table 


PRIVATE INTERFACE

  Private Data Structures:

    struct ArbHashEntry

      this structure used to store entries in the table 

      Member Variables:

        KeyType key - entry's key 

        EntryType entry - entry's value 

        ArbHashEntry *next - pointer to the next entry in a 
            linked list 

  Private Member Functions:

    BucketIndex - map a key to a bucket 

      int BucketIndex(const int i) const

        i - (in)  key 

      Description: This function takes an integer key and maps it to 
          a hash bucket, i.e., this is the hash function. 

      Return Value: the bucket associated with the input key 


    RebuildTable - rebuild the hash table to increase its size 

      void RebuildTable()

      Description: This function increases the number of buckets by 
          the factor RebuildMultiplier and rehashes all table 
          entries. 


  Private Member Variables:

        static const int RebuildMultiplier - factor by which a table 
                will grow when it is rebuilt. 

        bool AllowDuplicates - flag that indicates if duplicate keys 
                are allowed 

        ArbHashEntry **ppBuckets - list of hash table buckets 

        ArbHashEntry *pStaticBuckets[SMALL_HASH_TABLE] - when the table 
                is small, this statically allocated array of buckets 
                are used 

        int NumBuckets - current number of buckets 

        int NumHashEntries - total number of entries in the table 

        int RebuildSize - number of entries at which the table will be 
                rebuilt 

        int DownShift - variable used by the hash function that depends 
                on the current table size 

        int Mask - variable use by the hash function that depends on 
                the current table size 

-----------------------------------------------------------------------

TEMPLATE CLASS CArbHashTableIterator<class KeyType,class EntryType>

  This is an object that allows one to sequence through all the entries 
  in the table. The entries are visited in arbitrary order. 

  A typical use of the object could look something like: 

  CArbHashTableIterator<KeyType,EntryType> iterator ; 
  for (iterator.First() ; iterator.More() ; ++iterator()) { 
  ... 
  } 


  Template Arguments:

    KeyType - key type

    EntryType - value type


PUBLIC INTERFACE

  Public Member Functions:

    CArbHashTableIterator - constructor for hash table iterators 

      CArbHashTableIterator(const CArbHashTable<KeyType,EntryType>* aTable)

        aTable - (in)  address of the hash table to be visited 

      Description: This is a constructor for HashTableIterator 
          objects. 


    First - point to the first element 

      void First()

      Description: This function sets the iterator to point to the 
          first entry in the table. 


    Next - move to the next item in the table 

      void Next()

      Description: This function sets the iterator to point to the 
          next entry in the table. 


    More - query to see if there are more entries in the table 

      bool More()

      Description: This function tells if there are more entries in 
          the table or if the end of the table has been reached. 

      Return Value: True indicates that the iterator points to a 
          valid table entry. False indicates that the end of the 
          table has been reached and the iterator does not point to a 
          valid table. 


    Key - return the entry key 

      KeyType Key()

      Description: Returns the key of the hash table entry that the 
          iterator is currently pointing to. 

      Return Value: the key of the current entry 


    Entry - return the entry value 

      EntryType *Entry()

      Description: Returns a pointer to the value of the hash table 
          entry that the iterator is currently pointing to. 

      Return Value: the value of the current entry 


    operator ++ - increment operator 

      void operator ++ ()

      Description: This operator increments the iterator so that it 
          points to the next entry in the table. 


    operator ++ - increment operator 

      void operator ++ (int i)

        i - (in)  dummy argument 

      Description: This is the postfix form of the increment 
          operator. 


PRIVATE INTERFACE

  Private Member Variables:

    const CArbHashTable<KeyType,EntryType>* table - pointer to the 
            hash table that we are iterating over. 

    int bucket - bucket we are currently visiting in the hash table 

    CArbHashTable<KeyType,EntryType>::ArbHashEntry* ptr - entry we 
            are currently visiting in the hash table 


*/

#endif
