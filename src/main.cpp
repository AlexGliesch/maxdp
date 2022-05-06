/*
 * A Hybrid Heuristic for the Maximum Dispersion Problem
 * Copyright (c) 2020 Alex Gliesch, Marcus Ritt
 *
 * Permission is hereby granted, free of charge, to any person (the "Person")
 * obtaining a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software, including the rights to use, copy, modify,
 * merge, publish, distribute the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * 1. The above copyright notice and this permission notice shall be included in
 *    all copies or substantial portions of the Software.
 * 2. Under no circumstances shall the Person be permitted, allowed or authorized
 *    to commercially exploit the Software.
 * 3. Changes made to the original Software shall be labeled, demarcated or
 *    otherwise identified and attributed to the Person.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
#include "main.h"
#include "bal.h"
#include "btdsatur/bktdsat.h"
#include "btdsatur/colorrtns.h"
#include "btdsatur/graph.h"
#include "btdsatur/mysys.h"
#include "cmdline.h"
#include "color.h"
#include "coudert.h"
#include "mcmf.h"
#include "oscillate.h"
#include "readinstance.h"
#include "solution.h"
#include "ub.h"
#include "umdp.h"
#include "umdpfernandez.h"
bool exited = false;
void exitFun(int) {
  if (not exited) {
    exited = true;
    if (lastExitCode == EXIT_SUCCESS) {
      stats::printStats();
      if (lbSol.isComplete()) {
        ofstream of(outputFilename);
        if (of.good()) {
          for (int i = 0; i < n; ++i) {
            if (i > 0) of << ' ';
            of << lbSol.a[i];
          }
        }
      }
    }
  }
  exit(lastExitCode);
}
void exitFun2() {
  if (not exited) exitFun(0);
}
namespace btdsatur {
extern int colorsearch(int targetnumcolors, Timer t);
extern unsigned long long numConfChecks;
extern unsigned long long maxChecks;
extern int verbose;
} // namespace btdsatur
int main(int argc, char** argv) {
  signal(SIGINT, exitFun);
  atexit(exitFun2);
  try {
    parseCommandLine(argc, argv);
    readInstance(inputFilename);
    globalTimer.reset(timeLimit);
    ran = true;
    if (tabuColTest) {
      assert(tabuColDValue > 0);
      double minD = NLD::max();
      int minI = -1;
      for (uint i = 0; i < Rd.size(); ++i) {
        double d = abs(tabuColDValue - Rd[i]);
        if (d < ff::EPS and d < minD) minD = d, minI = i;
      }
      if (minI == -1) {
        print("d value ({}) for tabuCol is not valid for this instance.\n",
              tabuColDValue);
        exit(EXIT_SUCCESS);
      }
      [[maybe_unused]] VI coloring;
      int colors = colorTabuCol(minI, coloring, globalTimer);
      print("{}", colors);
      exit(EXIT_SUCCESS);
    }
    pr("|R|: {}\n", Rd.size());
    if (testType == "ub") {
      computeUpperBound(Timer(ubTimeLimit, globalTimer), true);
    } else if (testType == "umdp") {
      computeUpperBound(Timer(ubTimeLimit, globalTimer), true);
      umdpSol = umdp(Timer(umdpTimeLimit, globalTimer), true);
      assert(umdpSol.isComplete());
      pr("UMDP ended!\n");
    } else if (testType == "full") {
      oscillate(lbSol, globalTimer, true);
    }
  } catch (std::exception& e) {
    print("Exception caught (instance {}): {}\n", inputFilename, e.what());
    exit(EXIT_FAILURE);
  }
}
