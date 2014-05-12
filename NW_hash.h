// -*- Mode: C++;  -*-
// File: hash.h
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:
//*
//* CLASSES:
//* 
//* RELATED PACKAGES:
//*
//* HISTORY:
//*  Apr 19 10:59 1994 (sjr): hash --> hash()
//*                           tested okay.
//*  Apr  7 17:23 1994 (sjr): Tested Okay.
//* Created: Thu Apr  7 12:25:06 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//  Ref:  code heavily borrowed from the Sather 0.5 libraries
//  -- Author: Stephen M. Omohundro (om@ICSI.Berkeley.EDU)
//-- Copyright (C) International Computer Science Institute, 1991, 1992, 1993
//--
//-- COPYRIGHT NOTICE: This code is provided "AS IS" WITHOUT ANY WARRANTY
//-- and is subject to the terms of the SATHER LIBRARY GENERAL PUBLIC
//-- LICENSE contained in the file: "sather/doc/license.txt" of the Sather
//-- distribution. The license is also available from ICSI, 1947 Center
//-- St., Suite 600, Berkeley CA 94704, USA.
//--
//-- Changes: Heinz W. Schmidt (hws@csis.dit.csiro.au)
//-- (c) Commonwealth Scientific and Industrial Research Organisation (CSIRO),
//-- Australia, 1992, 1993.

// TYPE must have routine is_equal and UINT function hash() defined

#ifndef NW_HASH_H
#define NW_HASH_H

#include <assert.h>
#include <memory.h>
#include <stdio.h>
#include <string.h>
#include "NW_type.h"

template <class TYPE>
class HASH_TABLE {
  TYPE** table;
  UINT    mask;
  UINT    alloc;

  VOID   double_size() {
    TYPE *entry;
    UINT new_alloc = (alloc - 1) * 2 + 1;
    mask = (mask << 1) + 1;
    TYPE** new_table = new TYPE*[new_alloc];
    memset(new_table, 0, new_alloc*sizeof(TYPE*));
    size = 0;
    for(UINT i = 0; i < alloc; i++) {
      entry = table[i];
      if(entry != NULL) {
	UINT h = entry->hash() & mask;
	while(new_table[h] != NULL){
	  if(h == new_alloc - 2)
	    h = 0;
	  else
	    h++;
	}
	new_table[h] = entry;
	size++;
      }
    }
    delete table;
    table = new_table;
    alloc = new_alloc;
  }

public:
  HASH_TABLE() { 
    size = 0; alloc = 9; mask = 7; 
    table = new TYPE*[alloc]; 
    memset(table, 0, alloc*sizeof(TYPE*));
  }

  ~HASH_TABLE() { 
    for(UINT i = 0; i < alloc; i++)
      delete table[i];
    delete table; 
  }
  
  UINT size;
  UINT get_size() { return(size); }
  UINT get_alloc() { return(alloc); }

  TYPE *get(TYPE*entry) {
    UINT h = entry->hash() & mask;
    TYPE *res;
    while(1) {
      while((res = table[h]) != NULL) {
	if(res->is_equal(entry)) {
	  return(res); // found it
	}
	h++;
      }
      if(h != alloc-1) break;	// failed
      h = 0;			// hit the end so back to the start
    }
    return(NULL);			// if we get here we've failed...
  }

  TYPE* get(UINT i) { return(table[i]); }


  BOOL test(TYPE* e) {
    return(get(e) != NULL);
  }

  VOID insert(TYPE* e) {
    if (((size + 1) << 2) > alloc)  double_size();
    UINT h = e->hash() & mask;
    while(table[h] != NULL) {
      if(table[h]->is_equal(e)){
	table[h] = e; 
	return;
      }
      if (h == alloc-2) 
	h = 0;
      else 
	h++;
    }
    table[h] = e;
    size++;
  }

  VOID remove(TYPE* e) {
    // Delete element e if present
    UINT h = e->hash() & mask;
    UINT hole;
    while(!table[h]->is_equal(e)) {
      if(table[h] == NULL) 
	return;			// Not in the table
      else if (h == alloc-2)
	h = 0;
      else
	h++;
    }
    hole = h; table[hole] = NULL; size--; // Delete it
  
    // Now move up and fill in the hole if necessary
    UINT index = hole;
    while(1){
      if(index == alloc - 2) index = 0;
      else                  index++;
      if(table[index] == NULL) break;
      h = table[index]->hash() & mask;
      if(h <= index) {		// No wraparound
	if(hole < index && hole >= h){ // Fill in the hole
	  table[hole] = table[index];
	  table[index] = NULL;
	  hole = index;
	}
      }
      else {			// Wraparound
	if(hole >= h || hole < index) { // Fill in the hole
	  table[hole] = table[index];
	  table[index] = NULL;
	  hole = index;
	}
      }
    }
  }

  VOID clear() {
    delete table;
    size = 0; alloc = 9; mask = 7; 
    table = new TYPE*[alloc]; 
    memset(table, 0, alloc*sizeof(TYPE*));
  }

  BOOL is_empty() { return (size==0); }

};

template <class T>
class STR_HASH_ENTRY {
private:
  UINT hash_cache;
public:
  STR_HASH_ENTRY(STRING s, T v) 
  { key = s;  val = v; set_hash(); }

  STRING key;
  T      val;
  UINT    hash() { return(hash_cache); }

  //  Gives a 22 bit hash value (rightmost bits)
  //  okay for up to 1 million elements
//  VOID set_hash() {
//    UINT i = 0; 
//    hash_cache = 0;
//    while(key[i] != '\0'){
//      hash_cache =  hash_cache ^ ((UINT)key[i] << (i & 0xf));
//      i++; 
//    }
//    hash_cache &= 0x3fffff;
//    if((hash_cache = hash_cache - (((hash_cache << (i % 10)) | 1) & 0x3fffff)) < 0)
//      hash_cache += 0x3fffff;
//  }

  // 
  // Sstrhash: returns a 31 bit hash of the bottom 7 bits in the given string 
  // uses the "Additive Congruential Method" as described by Robert
  //   Sedgewick in "Algorithms", (2nd edition), page 515.  This is used
  //   over any "modulo" method as it is independent of the sizeof the hash
  //   table.  The feedback tap has been chosen to be as small as possible
  //   and so help randomise very short strings 
  // (Thanks Tony.)
  VOID set_hash() {
    CCHAR* cptr = key;

    hash_cache = 0;
    while(*cptr != '\0')
      hash_cache = (hash_cache << 7 & 0x7ffffff8) |
	((hash_cache >> 6 ^ hash_cache >> 24 ^ *cptr++) & 0x7f);
  }
  

#if 0
  //  Gives a 31 bit hash value (rightmost bits)
  VOID set_hash() {
    UINT i = 0; 
    hash_cache = 0;
    while(key[i] != '\0'){
      hash_cache =  hash_cache ^ ((UINT)key[i] << (i % 25));
      i++; 
    }
    hash_cache &= 0x7fffffff;
    if((hash_cache = hash_cache - (((hash_cache << (i % 10)) | 1) & 0x7fffffff)) < 0)
      hash_cache += 0x7fffffff;
  }  
#endif

  BOOL is_equal(STR_HASH_ENTRY<T> *e) 
  { return( hash() == e->hash() && !strcmp(key, e->key)); }
};

template <class T>
class STR_HASH_TABLE {
  HASH_TABLE < STR_HASH_ENTRY <T> > *tbl;
  STR_HASH_ENTRY<T> *entry;
  T null_value;

public:
  STR_HASH_TABLE(T nv = NULL) { 
    null_value = nv;
    tbl = new HASH_TABLE<STR_HASH_ENTRY <T> >(); 
    entry = new STR_HASH_ENTRY <T> ("", null_value);
  }

  ~STR_HASH_TABLE() { delete tbl; delete entry; }

  T get(STRING s) {		
    entry->key = s;  
    entry->set_hash();
    STR_HASH_ENTRY<T> *res = tbl->get(entry);
    if(res != NULL)
      return res->val;
    else
      return null_value;
  }

  UINT size() { return(tbl->size); }

  BOOL test(STRING s) {
    entry->key = s;  
    entry->set_hash();
    STR_HASH_ENTRY<T> *res = tbl->get(entry);
    return(res != NULL);
  }


  BOOL is_empty() { return (size()==0); }

  VOID insert(STRING s, T v) { 
    tbl->insert(new STR_HASH_ENTRY<T>(s, v)); 
  }

  VOID remove(STRING s) { 
    entry->key = s; entry->set_hash();
    tbl->remove(entry);
  }

  VOID clear() { tbl->clear(); }

  VOID dump() { 
    printf("size = %d, alloc = %d\n", tbl->get_size(), tbl->get_alloc()); 
  } 

  T get_null() { return null_value; }

};


#endif

