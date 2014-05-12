// -*- Mode: C++;  -*-
// File: queue.h
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:  FIFO queue 
//*
//* CLASSES:
//* 
//* RELATED PACKAGES:
//*
//* HISTORY:
//*  Apr 19 10:47 1994 (sjr):  added cursor routines, append() copy()
//*                            tested okay.
// *  Apr 12 10:10 1994 (sjr): fixed bug in resize (> --> >=), works in
//                             breadth-first search of lexicon
// *  Apr 11 11:01 1994 (sjr): tested okay
//* Created: Mon Apr 11 10:06:21 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef NW_QUEUE_H
#define NW_QUEUE_H

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

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif


template <class TYPE>
class QUEUE {
protected:
  TYPE *data;
  INT alloc;
  INT default_alloc;
  INT head;
  INT tail;
  INT cursor;
  INT count;

  VOID resize(INT n) {
    assert(n >= size);
    TYPE* new_data = new TYPE[n];
    copy_data(new_data);
    head = 0;
    tail = size;
    alloc = n;
    delete[] data; 
    data = new_data;
  }

public:
  INT size;
  BOOL is_empty() { return(size == 0); }

  QUEUE(INT n = 5){ 
    size = 0; head = 0; tail = 0; 
    alloc = n; data = new TYPE[n]; 
    default_alloc = n;
  }

  ~QUEUE() { delete[] data; }

  /*virtual*/ VOID clear() { size = 0; head = 0; tail = 0; resize(default_alloc); }
  /*virtual*/ VOID fast_clear() { size = 0; head = 0; tail = 0; }

  VOID enqueue(TYPE val) { 
    if (size == alloc)
      resize(alloc*2);
    data[tail] = val;
    if(++tail == alloc)
      tail = 0;
    size++;
  }

  VOID push(TYPE val) { enqueue(val); }

  /*virtual*/ TYPE dequeue() { 
    assert(size > 0);
    TYPE res = data[head];
    if(++head == alloc)
      head = 0;
    size--;
    return(res);
  }

  TYPE pop() { return(dequeue()); }

  TYPE top() { 
    assert(size > 0);
    return(data[head]); 
  }

  TYPE second() { 
    assert(size > 1);
    INT i = (head+1)%alloc;
    return(data[i]); 
  }

  TYPE last() {
    assert(size > 0);
    if(tail == 0)
      return(data[alloc-1]);
    else
      return(data[tail-1]);
  }

  /*virtual*/ TYPE get(INT i) {
    assert(i >= 0 && i < size);
    return(data[(head+i)%alloc]);
  }

  /*virtual*/ VOID put(INT i, TYPE val) {
    assert(i >= 0 && i < size);
    data[(head+i)%alloc] = val;
  }

  /*virtual*/ INT index(TYPE val) {
    INT i = 0;
    while(i < size){
      if(data[(head+i)%alloc] == val) return(i);
      i++;
    }
    return(-1);
  }

  VOID init_cursor() { 
    cursor = (head)%alloc - 1; 
    count = 0; 
  }

  VOID init_reverse_cursor() { 
    assert(size > 0); 
    count = size-1; 
    if(tail == 0)
      cursor = alloc-1;
    else
      cursor = tail-1;
  }


  BOOL next() {
    if(++count > size) {
      return(FALSE);
    } else {
      cursor = (cursor+1)%alloc;
      return(TRUE);
    }
  }

  TYPE read() {
    return(data[cursor]);
  }
  
  VOID write(TYPE val) {
    data[cursor] = val;
  }

  INT previous() {
    if(--count < 0){
      return(-1);
    }  else {
      if(--cursor < 0)
	cursor += alloc;
      return(count);
    }
  }

  /*virtual*/ QUEUE* copy() {
    QUEUE *res = new QUEUE(alloc);
    res->size = size;
    res->head = head;
    res->tail = tail;
    memcpy(res->data, data, size*sizeof(TYPE));
    return(res);
  }

  VOID copy_data(TYPE* mem) {
    if(size > 0){
      if(head >= tail) {
	INT j = alloc - head;
	memcpy(mem, &data[head], j*sizeof(TYPE));
	memcpy(&mem[j], data, (size-j)*sizeof(TYPE));
      } else {
	memcpy(mem, &data[head], size*sizeof(TYPE));
      }
    }
  }

  QUEUE* copy_insert(INT i, TYPE val) {
    assert(i < size);
    if(size == 1) { // so i == 0
      QUEUE *res = new QUEUE();
      res->enqueue(val);
      res->enqueue(data[head]);
      return(res);
    }
    INT a = alloc;
    INT sz = size+1;
    if(a <= sz)
      a *= 2;
    QUEUE *res = new QUEUE(a);
    res->size = sz;
    res->head = 0;
    res->tail = sz;

    // copy over the data
    INT nleft;
    INT indx = head;
    INT n;

    if(i > 0) {
      // copy chunk 1
      n = MIN(i, alloc-head);
      memcpy(res->data, &data[indx], n*sizeof(TYPE));
      indx += n;
      if(indx==alloc) 
	indx = 0;
    
      // loopback if necessary
      if((nleft = (i-n)) > 0) {
	memcpy(&res->data[n], data, nleft*sizeof(TYPE));
	indx += nleft;
      }
    }

    // insert i
    res->data[i] = val;

    // copy chunk 2
    nleft = size-i;
    n = MIN(nleft, alloc-indx);
    memcpy(&res->data[i+1], &data[indx], n*sizeof(TYPE));

    // loopback if necessary
    if(n < nleft) 
      memcpy(&res->data[i+1+n], data, (nleft-n)*sizeof(TYPE));

    return(res);
    
  }

  VOID append(QUEUE *q) {
    INT a = alloc;
    while(a <= size + q->size)
      a *= 2;
    resize(a);			// this has the side effect of making head = 0
    q->copy_data(&data[tail]);
    tail += q->size;
    size += q->size;
    q->fast_clear();
  }

  VOID squeeze() { resize(size); }

  /*virtual*/ VOID info(STRING title = "QUEUE") {
    int i, j;
    printf("Queue %s %p:\n", title, this);
    printf("  size = %d (%d allocated)  data = %p\n", 
	   size, alloc, data);
    printf("  head = %d  tail = %d\n", head, tail);
    for(j = 0, i = head; j < size; j++){
      printf("  [%d] = %p\n", i, (void*)data[i]);
      i = (i+1) % alloc;
    }
  }

  TYPE* get_data_ptr() { return(data); }

};

template <class TYPE>
class QUEUE_OFFSET : public QUEUE<TYPE> {

protected:
  INT offset;

public:
  QUEUE_OFFSET(INT n = 32) : QUEUE<TYPE>(n) {
    offset = 0;
  }

  VOID clear() { this->size = 0; this->head = 0; this->tail = 0; offset = 0; resize(this->default_alloc); }
  VOID fast_clear() { this->size = 0; this->head = 0; this->tail = 0; offset =0; }

  
  TYPE dequeue() { 
    assert(this->size > 0);
    TYPE res = this->data[this->head];
    if(++(this->head) == this->alloc)
      this->head = 0;
    this->size--;
    offset++;
    return(res);
  }

  TYPE get(INT i) {
    i -= offset;
    //    assert(i >= 0 && i < size);
    return(this->data[(this->head+i)%this->alloc]);
  }

  VOID put(INT i, TYPE val) {
    i -= offset;
    assert(i >= 0 && i < this->size);
    this->data[(this->head+i)%this->alloc] = val;
  }

  INT index(TYPE val) {
    INT i = 0;
    while(i < this->size){
      if(this->data[(this->head+i)%this->alloc] == val) return(i+offset);
      i++;
    }
    return(-1);
  }


  QUEUE_OFFSET* copy() {
    QUEUE_OFFSET* res = new QUEUE_OFFSET(this->alloc);
    res->size = this->size;
    res->head = this->head;
    res->tail = this->tail;
    memcpy(res->data, this->data, this->size*sizeof(TYPE));
    res->offset = offset;
    return(res);
  }

  VOID info(STRING title = "QUEUE_OFFSET") {
    int i, j;
    printf("Queue %s %p:\n", title, this);
    printf("  size = %d offset = %d (%d allocated)  data = %p\n", 
	   this->size, offset, this->alloc, this->data);
    printf("  head = %d  tail = %d\n", this->head, this->tail);
    for(j = 0, i = this->head; j < this->size; j++){
      printf("  [%d] = %p\n", i, (void*)this->data[i]);
      i = (i+1) % this->alloc;
    }
  }
  
};


#endif
