
#ifndef MSHSURF_IDSMAP_HPP
#define MSHSURF_IDSMAP_HPP

#include <algorithm>
#include <iostream>
#include <map>
using namespace std;

class mshsurfIndexMap
{
private:

  typedef struct key
  {
    int id[3];
  } Key;


  struct ltstr
  {
    bool operator()(const Key& s1, const Key& s2) const
    {
      if      (s1.id[0] < s2.id[0]) return 1;
      else if (s1.id[0] > s2.id[0]) return 0;
      else if (s1.id[1] < s2.id[1]) return 1;
      else if (s1.id[1] > s2.id[1]) return 0;
      else if (s1.id[2] < s2.id[2]) return 1;
      else if (s1.id[2] > s2.id[2]) return 0;
      else return 0;
    }
  };

  typedef map<const Key, void *, ltstr> TypeLocal;

  TypeLocal LocalMap;
  TypeLocal::iterator LocalInterator;

public:
   mshsurfIndexMap ( ) { };
  ~mshsurfIndexMap ( ) {clear ( );};
 
  // find method
  void *find (int id1, int id2, int id3)
  {
    Key copy;


    if (id3 == -1)
    {
      if (id1 < id2) {
        copy.id[1] = id1;
        copy.id[2] = id2;
      }
      else {
        copy.id[1] = id2;
        copy.id[2] = id1;
      }
      copy.id[0] = id3;
    }
    else
    {
      if (id1 < id2) {
        copy.id[0] = id1;
        copy.id[1] = id2;
      }
      else {
        copy.id[0] = id2;
        copy.id[1] = id1;
      }
      copy.id[2] = id3;
    }

    
    LocalInterator = LocalMap.find (copy);

    if (LocalInterator == LocalMap.end ())
    {
      //cout << "not found " << copy.id[0] << "-" 
      //  << copy.id[1] << "-" << copy.id[2] << endl; 
      return NULL;
    }
    else
    {
      //cout << "found " << copy.id[0] << "-" 
      //  << copy.id[1] << "-" << copy.id[2] <<  LocalInterator->second << endl; 
      return (LocalInterator->second);
    }
  }

  // insert method
  void insert (int id1, int id2, int id3, void *info)
  {
    Key copy;

    if (id3 == -1)
    {
      if (id1 < id2) {
        copy.id[1] = id1;
        copy.id[2] = id2;
      }
      else {
        copy.id[1] = id2;
        copy.id[2] = id1;
      }
      copy.id[0] = id3;
    }
    else
    {
      if (id1 < id2) {
        copy.id[0] = id1;
        copy.id[1] = id2;
      }
      else {
        copy.id[0] = id2;
        copy.id[1] = id1;
      }
      copy.id[2] = id3;
    }

    LocalMap.insert (make_pair(copy, info));
  }

  // get the number of entities
  int size (void) { return ((int) LocalMap.size ()); }

  // clear all
  void clear (void) { LocalMap.clear (); }

  // first entity
  void *first ( )
  {
    LocalInterator = LocalMap.begin( );
    return (LocalInterator->second);
  }

  // next entity
  void *next ( )
  {
    ++LocalInterator;
    if (LocalInterator == LocalMap.end ())
      return NULL;
    else
      return (LocalInterator->second);
  }

};

#endif
