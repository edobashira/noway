// -*- Mode: C++;  -*-
// File: list.h
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
//*  Jan 27 12:59 1995 (sjr): removed fast_clear() from append; added
//*                            reverse_append()
//*  Apr 19 10:11 1994 (sjr): added append() second() put_and_grow().  
//*                           Tested okay.
//*  Apr 12 15:35 1994 (sjr): added remove() and index()
//*  Apr 12 10:11 1994 (sjr): get(), top() and pop() of empty list now
//                             caught by assert, don't return NULL.
//*  Apr  7 11:44 1994 (sjr):  Tested okay
//* Created: Wed Apr  6 15:37:07 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef NW_LIST_H
#define NW_LIST_H

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
class LIST {
 protected:
  TYPE *data;
  INT alloc;
  INT default_alloc;

  VOID resize(INT n) {
    if(n <= alloc)
      return;
    TYPE* old_data = data;
    data = new TYPE[n];
    memcpy(data, old_data, size*sizeof(TYPE));
    delete[] old_data;    
    alloc = n;
  }

  INT cursor;

 public:
  INT size;
  BOOL is_empty() { return(size == 0); }

  LIST(INT n = 5){ 
    size = 0; 
    cursor = -1;
    if(n > 0) {
      alloc = n; 
      default_alloc = n; 
      data = new TYPE[n]; 
    }
  }

  ~LIST() { delete[] data; }

  VOID fast_clear() { size = 0; }
  VOID clear() { size = 0;  resize(default_alloc); }

  LIST<TYPE> *copy() { 
    LIST<TYPE> *res = new LIST<TYPE>(alloc);
    res->size = size;
    memcpy(res->data, data, size*sizeof(TYPE));
    return(res);
  }

  VOID copy(LIST *l) {
    INT sz = l->size;
    if(sz > alloc) {
      INT a = alloc;
      while(sz > a)
	a *= 2;
      resize(a);
    }
    size = sz;
    memcpy(data, l->data, size*sizeof(TYPE));    
  }

  VOID init_cursor() { 
    cursor = size;
  }

  BOOL next() {
    return(--cursor >= 0);
  }

  TYPE read() {
    return(data[cursor]);
  }

  INT get_cursor() {
    return cursor;
  }

  VOID write(TYPE val) {
    data[cursor] = val;
  }

  BOOL end_cursor() {
    return(cursor < 0);
  }


  VOID push(TYPE val) { 
    if (size == alloc)
      resize(alloc*2);
    data[size++] = val;
  }

  VOID append(LIST *l) {
    if(l == NULL || l->size == 0) return;
    INT a = alloc;
    while(size+l->size > a)
      a *= 2;
    resize(a);
    memcpy(&data[size], l->data, l->size*sizeof(TYPE));
    size += l->size;
  }

  VOID append(TYPE *newdata, INT n) {
    if(newdata == NULL || n == 0) return;
    INT a = alloc;
    while(size+n > a)
      a *= 2;
    resize(a);
    memcpy(&data[size], newdata, n*sizeof(TYPE));
    size += n;
  }

  VOID reverse_append(LIST *l) {
    if(l == NULL || l->size == 0) return;
    INT a = alloc;
    while(size+l->size > a)
      a *= 2;
    resize(a);
    for(INT i = l->size-1; i >= 0; i++)
      data[size++] = l->get(i);
  }

  TYPE pop() { 
    assert(size > 0);
    return(data[--size]); 
  }

  TYPE top() { 
    assert(size>0);
    return(data[size-1]); 
  }

  TYPE nth(INT i) {
    assert(size >= i);
    return(data[size-i]);
  }

  TYPE get(INT i) { 
    assert(i >= 0 && i < size);
    return(data[i]); 
  }

  inline TYPE& operator[](INT i) {
    assert(i >= 0 && i < size);
    return(data[i]); 
  }

  VOID put(INT i, TYPE val) { 
    assert(i >= 0 && i < size);
    data[i] = val; 
  }

  VOID put_and_grow(INT i, TYPE val, TYPE nullval) { 
    assert(i >= 0);
    if(i < size){
      data[i] = val; 
    } else if (i == size) {
      push(val);
    } else { 
      if(i >= alloc) {
	INT a = alloc;
	while(i >= a)
	  a *= 2;
	resize(a);
      }
      for(INT j = size; j < i; j++)
	data[j] = nullval;
      data[i] = val;
      size = i+1;
    }
  }

  VOID remove(INT i) {
    assert(i >= 0 && i < size);
    size--;
    for(; i < size; i++)
      data[i] = data[i+1];
  }

  INT index(TYPE val) {
    INT i = 0;
    while(i < size){
      if(data[i] == val) return(i);
      i++;
    }
    return(-1);
  }

  BOOL index_remove(TYPE val) {
    INT i;
    TYPE* p;
    for(i = 0, p = data; i < size; i++, p++){
      if(*p == val)
	break;
    }
    if(i == size)
      return(FALSE);
    size--;
    for(; i < size; i++)
      data[i] = data[i+1];
    return(TRUE);
  }  

  VOID set(TYPE val, INT sz = -1, INT start = 0) {
    if(sz == -1) 
      sz = size;
    assert(start >= 0 && start <= size);
    INT end  = sz + start;
    sz = MAX(size, end);
    INT a = alloc;
    while(sz > a)
      a *= 2;
    resize(a);
    for(INT i = start; i < end; i++)
      data[i] = val;
    size = sz;
  }

  VOID squeeze() { resize(size); }

  VOID sort(INT (*compar)(const VOID *, const VOID *)) {
    qsort(data, size, sizeof(TYPE), compar);
  }

  VOID info(STRING title = "LIST", FILE *fp = stderr) {
    fprintf(fp, "List %s (0x%p):\n", title, this);
    fprintf(fp, "  size = %d (%d allocated)  data = 0x%p\n", size, alloc, data);
    //    for(INT i=0; i<size; i++)
    //      fprintf(fp, "  [%d] = 0x%p\n", i, (void*)data[i]);
  }

  TYPE* get_data_ptr() { return(data); }
};


#endif
