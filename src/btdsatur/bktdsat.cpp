#include "../util.h"
#include "mysys.h"

#include <float.h>
#include <fstream>
#include <iomanip>
#include <limits.h>

#include "bktdsat.h"
#include "colorrtns.h"
#include "graph.h"
#include "maxclique.h"

namespace btdsatur {

extern unsigned long long numConfChecks;
extern unsigned long long maxChecks;
// extern ofstream timeStream, checksStream;
extern int verbose;
extern clock_t startTime;

/*choose minsat */

/*#define DEBUG */
/*#define IMPACT */

#define MAXCLR 1280
#define LONG_BITS BITS(long)
#define MOD_MASK (LONG_BITS - 1)

#ifdef ANZAC
#define SHIFT_VALUE 6
#else
#define SHIFT_VALUE 5
#endif

/* prototypes */
void BlockColor(vertextype v, colortype c, colortype maxclr, int branch,
                popmembertype* m, Timer t);
void FindPair(colortype maxclr, vertextype* v, colortype* c, int* impval);
void move(vertextype v, int newsatur);
void fix(void);
void ApplyColor(vertextype v, colortype c, colortype maxclr, int branch,
                popmembertype* m, Timer t);
void Color(colortype maxclr, int branch, popmembertype* m, Timer t);
#ifndef IMPACT
int impact(vertextype v, colortype c);
#endif

/* global variables */

vertextype nextv[MAXCLR + MAXVERTEX];
vertextype Prev[MAXCLR + MAXVERTEX];
vertextype lclindex[MAXCLR + MAXVERTEX];

/* how many colors are conflicting */
vertextype satur[MAXVERTEX];
/* pointer to position in lists */
vertextype current[MAXCLR];

/* total of each adjacent color to vertex */
short clrset[MAXVERTEX][MAXCLR];

#ifdef IMPACT
/* color impacts */
short impact[MAXVERTEX][MAXCLR];
#endif

colortype bestcolor, maxsat, minsat;
vertextype numcolored;
popmembertype bestmember;
int Fixed;
colortype target;
int maxbranch;
int minlimit, maxlimit;
int MinMax;

/* choose min or max sat? */

#ifdef IMPACT
/* code to print out impact array */
#define USECLR 10
void pimpact(void) {
  int i, j;

  printf("## ");
  for (i = 0; i < USECLR; i++)
    printf("%4d ", i);
  printf("\n");

  for (i = 0; i < order; i++) {
    printf("%2d ", i);
    for (j = 0; j < USECLR; j++)
      printf("%4d ", impact[i][j]);
    printf("\n");
  }
}
#endif

/* interrupt control */
int stopflag = 0;
void cpulimbk() {
  printf("CPU TIME EXCEEDED -- let me clean up\n");
  stopflag = 1;
}
void intfcbk() {
  printf("INTERRUPT DETECTED -- cleaning up\n");
  stopflag = 1;
}

/* m is an ineffecient way to organize data  NOTE */
void bktdsat(popmembertype* m, int branch, colortype targetclr, int min,
             int max, Timer t) {
  /* vertices of equal saturation are kept in a circular doubly linked
  list. There are MAXCLR such lists. Vertex i is represented 0 <= i <
  MAXVERTEX by the ith position. nextv and Prev indicate the nextv and
  Previous vertices in the list.
  lclindex indicates the lclindex number of the vertex
  in the permutation handed to brelaz. Each list is kept in its
  lclindex order. The positions MAXVERTEX <= i < MAXCLR represent the
  reference positions in the circular lists; lclindex for these
  positions is ENDLIST to make it easy to detect end of list, and
  by keeping ENDLIST large makes it easy to detect insertion
  conditions at end of list.
  */
  vertextype v, w;
  vertextype i, j;
  int cliquesize;

  stopflag = 0;

  // Here are some other parameters that have been hard coded from the original
  // version of Culberson This means at each iteration a vertex of maximum
  // saturation is chosen (conforming to the Dsatur algorithm)
  MinMax = 1;

  // Find a large clique
  cliquesize = maxclique(m);
  if (cliquesize > targetclr) {
    printf("WARNING: a clique > target color was found\n");
    target = cliquesize;
  } else {
    target = targetclr;
  }
  bestcolor = order + 1;

  /* initialize */
  maxsat = 0;
  minsat = 0;
  numcolored = 0;
  bestmember = *m;
  Fixed = 0;
  maxbranch = order;
  minlimit = min;
  maxlimit = max;

#ifdef IMPACT
  /* initialize impact values to degree */
  for (i = 0; i < order; i++) {
    for (j = 0; j < i; j++) {
      numConfChecks++;
      if (edge(i, j)) {
        impact[i][0]++;
        impact[j][0]++;
      }
    }
  }

  for (i = 0; i < order; i++) {
    for (j = 1; j < MAXCLR; j++) {
      impact[i][j] = impact[i][0];
    }
  }

  printf("initially\n");
#ifdef DEBUG
  pimpact();
#endif
#endif

  /* initially no color saturation */
  for (i = 0; i < order; i++) {
    for (j = 0; j < MAXCLR; j++)
      clrset[i][j] = 0;
    satur[i] = 0;
  }

  /* all vertices are on zero conflict list */
  for (i = MAXVERTEX; i < MAXCLR + MAXVERTEX; i++) {
    nextv[i] = Prev[i] = i;
    lclindex[i] = ENDLIST;
  }

  /* the 0th conflict list is anchored in MAXVERTEX position of array */
  w = MAXVERTEX;
  for (i = 0; i < order; i++) {
    lclindex[v = m->vc[i].vertex] = i;
    /* vertices are coming in from smallest to largest,
    so insert at end of list  thus never changing w.
    */
    nextv[v] = w;
    Prev[v] = Prev[w];
    nextv[Prev[w]] = v;
    Prev[w] = v;
  }

  Color(0, branch, m, t);
  *m = bestmember;
  m->clrdata.numcolors = bestcolor;
}

#ifndef IMPACT
int impact(vertextype v, colortype c) {
  adjacencytype* x;
  vertextype w;
  int impval = 0;

  numConfChecks++;
  initnbr(x, v);
  for (w = 0; w < order; w++) {
    numConfChecks++;
    if (isnbr(x, w) && lclindex[w] != ENDLIST && clrset[w][c] == 0) impval++;
  }
  return (impval);
}
#endif

void fix(void) {
  int j;

  /* initialize pointers to reference positions */
  for (j = 0; j < MAXCLR; j++)
    current[j] = j + MAXVERTEX;

  /* scan for maximum saturation list */
  while (nextv[current[maxsat]] == current[maxsat] && maxsat > 0)
    maxsat--;
  /* scan for min saturation list */
  while (nextv[current[minsat]] == current[minsat] && minsat < maxsat)
    minsat++;

  Fixed = 1;
}

void Color(colortype maxclr, int branch, popmembertype* m, Timer t) {
  /* maxcnf is maxsat */
  vertextype v;
  colortype c;
  int impval;

#ifdef DEBUG
  printf("maxsat =%d\n", maxsat);
#endif

  if (stopflag) {
    //     cout << "stopflag\n";
    return;
  }

  if (numcolored >= order) {
    if (maxclr < bestcolor) {
      // At this point we have reduced the number of colors by one So we record
      // the time and number of checks
      //       clock_t colDuration =
      //           (int)((double)(clock() - startTime) / CLOCKS_PER_SEC * 1000);
      //       if (verbose >= 1)
      //         cout << setw(5) << maxclr << setw(11) << colDuration << "ms\t"
      //              << numConfChecks << endl;
      // timeStream<<maxclr<<"\t"<<colDuration<<"\n";
      // checksStream<<maxclr<<"\t"<<numConfChecks<<"\n";
      if (maxclr <= target) {
        if (verbose >= 1)
          cout << "\nSolution with <= maxclr (" << maxclr
               << ") colours has been found. Ending..." << endl;
        stopflag = 1;
      }
      bestcolor = maxclr;
      bestmember = *m;
    } else {
      //       cout << "End of branch, no better coloring\n";
    }
  } else if (bestcolor <= target) {
    if (verbose >= 1)
      cout << "\nSolution with <= bestcolor (" << bestcolor
           << ") colours has been found. Ending..." << endl;
    //     clock_t colDuration =
    //         (int)((double)(clock() - startTime) / CLOCKS_PER_SEC * 1000);
    // timeStream<<bestcolor<<"\t"<<colDuration<<"\n";
    // checksStream<<bestcolor<<"\t"<<numConfChecks<<"\n";
    stopflag = 1;
  } else if (maxclr >= bestcolor) {
    if (verbose >= 1) cout << "Worse or equal coloring, returning\n";
  } else if (/*numConfChecks > maxChecks*/ t.timedOut()) {
    if (verbose >= 1) cout << "TimeoutException\n";
    throw TimeoutException(t);
    //     clock_t colDuration =
    //         (int)((double)(clock() - startTime) / CLOCKS_PER_SEC * 1000);
    //     if (verbose >= 1)
    //       cout << "\nRun limit exceeded. No solution using " << bestcolor - 1
    //            << " colours was achieved (Checks = " << numConfChecks << ", "
    //            << colDuration << "ms)" << endl;
    // timeStream<<bestcolor-1<<"\tX\t"<<colDuration<<"\n";
    // checksStream<<bestcolor-1<<"\tX\t"<<numConfChecks<<"\n";
    //     stopflag = 1;
  } else {
    fix();
    if (maxsat == maxclr) {
      /* some vertex is completely saturated */
      v = nextv[current[maxsat]];
      if (maxclr + 1 < bestcolor) {
        ApplyColor(v, maxclr + 1, maxclr, branch, m, t);
      }
#ifdef DEBUG
      else
        cout << bestcolor << ":" << numcolored << " Worse or equal coloring\n";
#endif
    } else {
      FindPair(maxclr, &v, &c, &impval);
      ApplyColor(v, c, maxclr, branch, m, t);
      if (maxclr >= bestcolor) return;
      /* if impact==0 then this color was not at fault */
      if (branch > 0 && impval > 0 &&
          (numcolored < minlimit || numcolored > maxlimit)) {
        if (numcolored <= maxbranch) {
          maxbranch = numcolored;
          fflush(stdout);
        }
        BlockColor(v, c, maxclr, branch - 1, m, t);
      }
    }
  }
}

/* #### use impact size array to save list of nbrs? */
void ApplyColor(vertextype v, colortype c, colortype maxclr, int branch,
                popmembertype* m, Timer t) {
  vertextype oldlclindex, w;
  int oldmaxsat, oldminsat;
  int j;
  adjacencytype* x;

  oldmaxsat = maxsat;
  oldminsat = minsat;

  /* pull v off its list */
  nextv[Prev[v]] = nextv[v];
  Prev[nextv[v]] = Prev[v];

  if (c > maxclr) maxclr = c;

  numcolored++;
  m->vc[lclindex[v]].color = c;
  oldlclindex = lclindex[v];
  lclindex[v] = ENDLIST; /* no longer on any list */

  /*#### use impact size array to save list of nbrs? */

  /* update saturation and impact lists */
  numConfChecks++;
  initnbr(x, v);
  for (j = 0; j < order; j++) {
    w = m->vc[j].vertex;
    numConfChecks++;
    if (isnbr(x, w) && lclindex[w] != ENDLIST) {
#ifdef IMPACT
      /* do impact */
      for (i = 1; i < MAXCLR; i++)
        if (clrset[v][i] == 0) impact[w][i]--;
#endif

      /* mark color in colorset and check if */
      /* color not Previously adjacent to w */
      if (0 == (clrset[w][c]++)) {
#ifdef IMPACT
        /* do impact */
        numConfChecks++;
        initnbr(x2, w);
        for (i = 0; i < order; i++) {
          w2 = m->vc[i].vertex;
          numConfChecks++;
          if (isnbr(x2, w2) && lclindex[w2] != ENDLIST) impact[w2][c]--;
        }
#endif

        /* move vertex to nextv list */
        move(w, satur[w] + 1);
        satur[w]++;
        if (maxsat < satur[w]) maxsat = satur[w];
#ifdef DEBUG
        printf("newmaxsat=%d\n", maxsat);
#endif
      }
    }
  }
  Fixed = 0;

#ifdef IMPACT
#ifdef DEBUG
  pimpact();
#endif
#endif

  Color(maxclr, branch, m, t);

  if (!Fixed) {
    fix();
#ifdef DEBUG
    printf("ERROR: current[] was not Fixed\n");
#endif
  }

  /* restore saturation and impact lists */
  for (j = 0; j < order; j++) {
    w = m->vc[j].vertex;
    numConfChecks++;
    if (isnbr(x, w) && lclindex[w] != ENDLIST) {
#ifdef IMPACT
      /* do impact */
      for (i = 0; i < MAXCLR; i++)
        if (clrset[v][i] == 0) impact[w][i]++;
#endif

      /* unmark color in colorset and check if */
      /* color now not adjacent to w */
      if (0 == (--clrset[w][c])) {
#ifdef IMPACT
        /* do impact */
        numConfChecks++;
        initnbr(x2, w);
        for (i = 0; i < order; i++) {
          w2 = m->vc[i].vertex;
          numConfChecks++;
          if (isnbr(x2, w2) && lclindex[w2] != ENDLIST) impact[w2][c]++;
        }
#endif

        /* assume satur[w]>0 */
        /* return vertex to Prev list */
        move(w, satur[w] - 1);
        satur[w]--;
        if (minsat > satur[w]) minsat = satur[w];
      }
    }
  }
  Fixed = 0;

#ifdef IMPACT
#ifdef DEBUG
  pimpact();
#endif
#endif

  /* put v back on its list */
  if (Prev[nextv[v]] != Prev[v])
    printf("ERROR: Prev nextv %d != Prev %d\n", v, v);
  if (nextv[Prev[v]] != nextv[v])
    printf("ERROR: nextv Prev %d != nextv %d\n", v, v);

  Prev[nextv[v]] = v;
  nextv[Prev[v]] = v;

  numcolored--;
  lclindex[v] = oldlclindex;
  m->vc[lclindex[v]].color = 0;

  maxsat = oldmaxsat;
  minsat = oldminsat;
}

void move(vertextype v, int newsatur) {
  vertextype z, zp;

  nextv[Prev[v]] = nextv[v];
  Prev[nextv[v]] = Prev[v];
  if (current[satur[v]] == v) current[satur[v]] = Prev[v];

  /* insert v into new list */
  /* note use of maximal ENDLIST lclindex */
  z = current[newsatur];
  while (lclindex[zp = nextv[z]] < lclindex[v])
    z = zp;

  nextv[v] = zp;
  Prev[v] = z;
  nextv[z] = v;
  Prev[zp] = v;
  current[newsatur] = v;

#ifdef DEBUG
  if (current[newsatur] == nextv[current[newsatur]] &&
      current[newsatur] < MAXVERTEX) {
    printf("ERROR: move(): current==nextv[current]\n");
  }
#endif
}

void FindPair(colortype maxclr, vertextype* v, colortype* c, int* impval) {
  int w, i, t;

  *impval = order;
  *c = 1;
  *v = 0;

  if (MinMax == 0)
    w = nextv[current[minsat]];
  else
    w = nextv[current[maxsat]];

  for (i = 1; i <= maxclr; i++) {
    if (clrset[w][i] == 0) {
#ifdef IMPACT
      t = impact[w][i];
#else
      t = impact(w, i);
#endif
      if (t < *impval) {
        *impval = t;
        *c = i;
        *v = w;
      }
    }
    if (*impval == 0) break;
  }
}

void BlockColor(vertextype v, colortype c, colortype maxclr, int branch,
                popmembertype* m, Timer t) {
  clrset[v][c] = order;

  fix();

  move(v, satur[v] + 1);
  satur[v]++;
  if (maxsat < satur[v]) maxsat = satur[v];
  Fixed = 0;

  Color(maxclr, branch, m, t);

  if (!Fixed) {
    fix();
#ifdef DEBUG
    fprintf(stderr, "ERROR: BlockColor: current[] was not Fixed\n");
#endif
  }

  move(v, satur[v] - 1);
  satur[v]--;
  if (minsat > satur[v]) minsat = satur[v];
  Fixed = 0;

  clrset[v][c] = 0;
}
} // namespace btdsatur
