// -*- Mode: C++;  -*-
// File: type.h
// Author: Steve Renals (sjr@eng.cam.ac.uk)
// Copyright (C) Cambridge University Engineering Department, 1994
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//* FUNCTION:  various typedefs
//*
//* CLASSES:
//* 
//* RELATED PACKAGES:
//*
//* HISTORY:
//* Created: Tues Apr 5 1994 (sjr)
//*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#ifndef NW_TYPE_H
#define NW_TYPE_H

typedef int BOOL;
typedef unsigned int UINT;
typedef int INT;
typedef unsigned long ULONG;
typedef long LONG;
typedef unsigned short USHORT;
typedef short SHORT;
typedef unsigned char UCHAR;
typedef char CHAR;
typedef const char CCHAR;
typedef unsigned const char CUCHAR;
typedef const char *STRING;
typedef float FLOAT;
typedef double DOUBLE;
typedef void VOID;

// constants
#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#endif
