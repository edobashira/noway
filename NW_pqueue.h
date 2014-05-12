// -*- Mode: C++;  -*-
// File: pqueue.h
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:
//*
//* CLASSES:  template <class TYPE> class PQUEUE 
//* 
//* RELATED PACKAGES:
//*
//* HISTORY:
//*  May 24 15:43 1994 (sjr): Fixed delete bug in replace
//*  May 21 10:49 1994 (sjr): Added replace(INT indx, TYPE t)
//*  Apr 22 16:32 1994 (sjr): Changed the array from TYPE** to TYPE*, so
//*                           indirection is specified by caller as
//*                           PQUEUE<T*>, similar to LIST, QUEUE, etc.
//*  Apr  6 14:25 1994 (sjr): separate (non-recursive) heapify routines for 
//*                           going up (insert) and down (pop) the heap.  
//*                           Tested okay.
//* Created: Tues Apr  5 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Priority queue implemented using a heap data structure
// A heap is an array which conforms to the heap property:
//    qitem[i]->is_less_than(qitem[parent(i)])
// So the largest item is always qitem[1] (data[0] is not used).
// class TYPE must have routine is_less_than(TYPE*) defined

// Note that this is not efficient for base types (int, float, etc)
// as the heap is stored as an array of pointers

// ref:  "Introduction to Algorithms", Cormen, Leiserson + Rivest
//       MIT Press 1990, Chapter 7


// First edited 5/4/94
// 6/4/94:  template <class TYPE> class PQUEUE with separate
//           (non-recursive) heapify routines for going up (insert)
//           and down (pop) the heap.  Tested okay.

#ifndef NW_PQUEUE_H
#define NW_PQUEUE_H

#include <assert.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include "NW_type.h"

#ifndef MIN
#define MIN(X,Y) ((X)>(Y) ? (Y) : (X))
#endif
#ifndef MAX
#define MAX(X,Y) ((X)>(Y) ? (X) : (Y))
#endif

template <class TYPE>
class PQUEUE {
protected:
  // heapified array of TYPEs;  data[0] is unused
  TYPE* data;
  INT alloc;
  INT default_alloc;

  VOID resize(INT n) {
    assert(n > size);
    alloc = n;
    TYPE* old_data = data;
    data = new TYPE[n];
    memcpy(data+1, old_data+1, size*sizeof(TYPE));
    delete[] old_data;
  }

  // Use this routine when we add an item to the end of the array (insert)
  virtual VOID heapify_up(INT l, INT u) {
    //    assert(u >= l && u >= 1 && l >= 1 && u <= size);
    INT i = u; 
    while (i > l) {
      INT j = parent(i);
      // If item is less than a parent then it's a heap
      if (data[i]->is_less_than(data[j])) 
	break;
      // Otherwise swap and compare the parent with its parent
      else {
	TYPE tmp = data[j];
	data[j] = data[i];
	data[i] = tmp;
	i = j;
      }
    }
  }

  // Use this when adding an item to the head of the array (pop, when 
  //  the final item is moved to the start)
  virtual VOID heapify_down(INT l, INT u) {
    //    assert(u >= l && u >= 1 && l >= 1 && u <= size);
    INT c, le, ri;
    INT i = l; 
    while(1) {
      le = left_child(i);
      ri = right_child(i);
      if(le > u) break;

      // set c to largest child
      if (ri <= u && data[le]->is_less_than(data[ri]))
	c = ri;
      else
	c = le;

      if(data[i]->is_less_than(data[c])) {
	// Not heapified as largest child is bigger, so swap
	TYPE tmp = data[c];
	data[c] = data[i];
	data[i] = tmp;
	i = c;
      }  else {
	// Heapified as largest child isn't bigger
	break;
      }
    }
  }

  // Tree functions
  INT parent(int i) { return(i>>1); }
  INT left_child(int i) { return(i<<1); }
  INT right_child(int i) { return(i<<1 | 1); }

public:

  PQUEUE(INT n = 7) {
    size = 0;
    n = MAX(n+1, 8);
    alloc = n;
    default_alloc = n;
    data = new TYPE[n];
  }

  // This doesn't delete the objects pointed to by the elements of data
  // works around compiler errors in 2.6.3 and 2.7.0
  virtual ~PQUEUE() { 
    if(data != NULL) {
      delete[] data;
      data = NULL;
    }  
  }


  INT size;
  BOOL is_empty() { return(size == 0); }

  VOID clear() { size = 0; resize(default_alloc); }
  VOID fast_clear() { size = 0; }

  TYPE pop() { 
    assert(size > 0);

    // return the top of the heap
    TYPE res = data[1]; 

    // Move the bottom to the top and reestablish the heap property
    data[1] = data[size--];
    if(size > 1) heapify_down(1,size);

    return(res);
  }

  TYPE top() { return(data[1]); }

  // assumes indexing from 0
  TYPE get(INT i) { return(data[i+1]); }

  VOID insert(TYPE h) {
    if(size >= alloc-2)
      resize(alloc*2);
    data[++size] = h;
    //    assert(size < alloc);
    heapify_up(1,size);
  }

  VOID push(TYPE h) { insert(h); }

  TYPE remove(INT i) {
    TYPE res;
    i = i+1;
    assert(i <= size && i > 0);
    if (i < size) {
      res = data[i];
      data[i] = data[size--];
      if(size > 1) heapify_down(i, size);
    } else { 
      // i == size
      res = data[size--];
    }
    return(res);
  }

  PQUEUE *copy() {
    PQUEUE *res = new PQUEUE(alloc-1);
    res->size = size;
    memcpy(&res->data[1], &data[1], size*sizeof(TYPE));
    return(res);    
  }

  // Arg assumes indexing from 0
  // Delete item i, replace with h, re-establish the heap property
  virtual VOID replace(INT i, TYPE h) {
    TYPE oldh = data[++i];
    data[i] = h;
    if(h->is_less_than(oldh)) {
      // heap property is okay above i, check below
      heapify_down(i, size);
    } else {
      // heap property is okay below i, check above
      heapify_up(1, i);
    }
    delete oldh;
  }

  // Arg assumes indexing from 0
  // Return item i, replace with h, re-establish the heap property
  virtual TYPE swap(INT i, TYPE h) {
    TYPE oldh = data[++i];
    data[i] = h;
    if(h->is_less_than(oldh)) {
      // heap property is okay above i, check below
      heapify_down(i, size);
    } else {
      // heap property is okay below i, check above
      heapify_up(1, i);
    }
    return oldh;
  }

};

// reverse PQUEUE (i.e. smallest items come out first)
// requires TYPE to have method is_greater_than
// not inherited for efficiency reasons (don't want to
//  use virtual functions).  Done due to gprof results - 
//  better to be safe than sorry...
//
template <class TYPE>
class MINPQUEUE : public PQUEUE<TYPE> {
protected:

  // Use this routine when we add an item to the end of the array (insert)
  VOID heapify_up(INT l, INT u) {
    //    assert(u >= l && u >= 1 && l >= 1 && u <= size);
    INT i = u; 
    while (i > l) {
      INT j = this->parent(i);
      if (this->data[i]->is_less_than(this->data[j])){
	// If item is less than a parent then it's not a min-heap,
	//  so swap and compare the parent with its parent
	TYPE tmp = this->data[j];
	this->data[j] = this->data[i];
	this->data[i] = tmp;
	i = j;
      }
      else {
	break;
      }
    }
  }

  // Use this when adding an item to the head of the array (pop, when 
  //  the final item is moved to the start)
  VOID heapify_down(INT l, INT u) {
    //    assert(u >= l && u >= 1 && l >= 1 && u <= size);
    INT c, le, ri;
    INT i = l; 
    while(1) {
      le = this->left_child(i);
      ri = this->right_child(i);
      if(le > u) break;

      // set c to smallest child
      if (ri <= u && this->data[le]->is_greater_than(this->data[ri]))
	c = ri;
      else
	c = le;

      if(this->data[i]->is_greater_than(this->data[c])) {
	// Not min-heapified as smallest child is smaller
	TYPE tmp = this->data[c];
	this->data[c] = this->data[i];
	this->data[i] = tmp;
	i = c;
      } else {
	// Min-heapified as smaller than smallest child
	break;
      }
    }
  }

public:

  MINPQUEUE(INT n = 7) : PQUEUE<TYPE>(n) { }

  // works around compiler errors in 2.6.3 and 2.7.0
  ~MINPQUEUE() {     
    if(this->data != NULL) {
      delete[] this->data;
      this->data = NULL;
    }  
  }


  TYPE last() {
    assert(this->size > 0);
    return(this->data[this->size]);
  }
  // Arg assumes indexing from 0
  // Delete item i, replace with h, re-establish the heap property
  VOID replace(INT i, TYPE h) {
    TYPE oldh = this->data[++i];
    this->data[i] = h;
    if(h->is_greater_than(oldh)) {
      // heap property is okay above i, check below
      heapify_down(i, this->size);
    } else {
      // heap property is okay below i, check above
      heapify_up(1, i);
    }
    delete oldh;
  }

  // Arg assumes indexing from 0
  // Return item i, replace with h, re-establish the heap property
  TYPE swap(INT i, TYPE h) {
    TYPE res = this->data[++i];
    this->data[i] = h;
    if(h->is_greater_than(res)) {
      // heap property is okay above i, check below
      heapify_down(i, this->size);
    } else {
      // heap property is okay below i, check above
      heapify_up(1, i);
    }
    return res;
  }

};

#endif
