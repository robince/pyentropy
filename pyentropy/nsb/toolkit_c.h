/*
 *  Copyright 2009, Weill Medical College of Cornell University
 *  All rights reserved.
 *
 *  This software is distributed WITHOUT ANY WARRANTY
 *  under license "license.txt" included with distribution and
 *  at http://neurodatabase.org/src/license.
 */

/* Make header C++ compatible */
#ifdef __cplusplus
extern "C" {
#endif

/* Some useful constants */
#define MAXCHARS 256
#define MAXPATH 260
#define BITS_IN_A_BYTE 8

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
//#include "../input/input_c.h"
#include "entropy_c.h"
//#include "hist_c.h"
#include "gen_c.h"
#include "sort_c.h"

#ifdef TOOLKIT
#include "mex.h"
#endif

/* Some useful macros */
#define NAT2BIT(x) x/log(2.0)
#define LOGZ(x) (x <= 0 ? 0 : log(x))
#define LOG2Z(x) NAT2BIT(LOGZ(x))
#define XLOGX(x) (-x*LOGZ(x))
#define XLOG2X(x) (-x*LOG2Z(x))
#define MAX(a,b) (a > b ? a : b) 
#define MIN(a,b) (a < b ? a : b) 
#define MIN3(a,b,c) MIN(MIN(a,b),c)

/* Define round and lround for compilers lacking these math.h functions (e.g., Microsoft Visual C++).
 * See http://www.velocityreviews.com/forums/t532986-rounding-functions-in-microsoft-visual-cc.html.
 */
#ifndef round
#define round(x) floor(x + 0.5)
#define lround(x) (long)floor(x + 0.5)
#endif

/* Define INFINITY/NAN for compilers lacking symbol (e.g., Microsoft Visual C++).
 * See http://www.gamedev.net/community/forums/topic.asp?topic_id=465682.
 */
#ifndef INFINITY
union MSVC_EVIL_FLOAT_HACK
{
	unsigned char Bytes[4];
	float Value;
};
static union MSVC_EVIL_FLOAT_HACK INFINITY_HACK = {{0x00, 0x00, 0x80, 0x7F}};
static union MSVC_EVIL_FLOAT_HACK NAN_HACK = {{0x00, 0x00, 0xC0, 0x7F}};
#define INFINITY (INFINITY_HACK.Value)
#define NAN (NAN_HACK.Value)
#endif

/* Make header C++ compatible */
#ifdef __cplusplus
}
#endif

