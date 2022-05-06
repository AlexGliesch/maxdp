#ifndef MYSYS
#define MYSYS 1
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

namespace btdsatur {
using namespace std;

/* time information */
extern long seconds, microsecs;
extern struct rusage tmp;

extern int getrusage(int who, struct rusage* rusage);

#include <signal.h>
extern int setrlimit(), getrlimit();

/* useful to help link from output files to res files */
// extern int getpid();

// #ifndef linux
// extern void ungetc();
// 
// extern int printf();
// extern int fprintf();
// extern int scanf();
// extern int fscanf();
// 
// extern int fgetc();
// 
// extern int fread();
// extern void fflush();
// extern void fclose();
// extern void rewind();
// #endif

/* the following required to make qsort calling args
silent */
typedef int (*compfunc)(const void*, const void*);
} // namespace btdsatur
#endif
