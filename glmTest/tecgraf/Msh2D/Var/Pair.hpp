//
// Pair.hpp
//
// Copyright -
//   (c) Fracture Analysis Consultants, Inc. 2004
//   All rights reserved
//
// Revision -
//   $Revision$  $Date$  $Author$
//
// ------------------------------------------------------------------
//

#ifndef Pair_hpp
#define Pair_hpp

namespace FTools {

/// \file
///
/// \brief This file defines the Pair and UPair class

/// \brief A templated class to store a pair of values

template<class Type0,class Type1> class Pair {
    public:

        Pair(Type0 item0,Type1 item1) : Item0(item0),Item1(item1) {} ;
        Pair() { Item0 = Type0() ; Item1 = Type1() ; } ;

        Pair(const Pair& other) { Item0 = other.Item0 ; Item1 = other.Item1 ; }

        Pair operator = (const Pair& other) {
            Item0 = other.Item0 ; Item1 = other.Item1 ; return *this ;
        }

        Type0 Item0 ;  ///< the first item
        Type1 Item1 ;  ///< the second item
} ;


/// \brief A templated class to store a pair of values of a uniform type
 
template<class Type> class UPair {
    public:

        UPair(Type item0,Type item1) : Item0(item0),Item1(item1) {} ;
        UPair() { Item0 = Type() ; Item1 = Type() ; } ;

        UPair(const UPair& other) { Item0 = other.Item0 ; Item1 = other.Item1 ; }

        UPair operator = (const UPair& other) {
            Item0 = other.Item0 ; Item1 = other.Item1 ; return *this ;
        }

        /// return the indexed item (0 or 1)

        Type operator[] (int indx) const 
        {
            assert(indx<2) ;

            return((indx==0) ? Item0 : Item1) ;
        }

        /// return r reference to the indexed item (0 or 1)

        Type& operator[] (int indx)
        {
            assert(indx<2) ;

            return((indx==0) ? Item0 : Item1) ;
        }

        Type Item0 ;  ///< the first item
        Type Item1 ;  ///< the second item
} ;

} // namespace

#endif
