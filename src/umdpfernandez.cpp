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
#include "umdpfernandez.h"
#include "color.h"
#include "main.h"
#include "solution.h"
#include "tabucol/PartialColAndTabuCol.h"
#include "tabucol/util/inputGraph.h"
int colorTabuCol(int d, VI &bestColoring, Timer t) {
  bestColoring.resize(n);
  gCol::Graph g(n);
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      if (di[i][j] < d) {
        ++g.nbEdges;
        g[i][j] = 1;
        g[j][i] = 1;
      }
  unsigned long long maxChecks = 100000000LL;
  int algorithm = 1, tenure = 0, targetCols = m, constructiveAlg = 1,
      verbose = 0;
  [[maybe_unused]] int timeLimitSeconds = t.secsLeft();
  confTimeLogger log;
  Timer colTimer;
  for (int iter = 1; iter <= umdpFernandezTabuColRepl; ++iter) {
    int seed = randInt(1, n);
    unsigned bestK = gCol::pctc::pctc(
        g, bestColoring, log, algorithm, tenure, maxChecks, targetCols, seed,
        verbose, constructiveAlg, umdpFernandezTabuColMaxIter);
    stats::ferTotalColTime += colTimer.elapsedSecs();
    stats::ferMaxColTime = max(stats::ferMaxColTime, colTimer.elapsedSecs());
    return bestK;
  }
  return n;
}
bool colorableTabuCol(int d, VI &bestColoring, Timer t) {
  int colors = colorTabuCol(d, bestColoring, t);
  return colors <= m;
}
Solution getSolutionFromColoring(VI coloring, [[maybe_unused]] int expectedD) {
  VI num(m, 0);
  for (int i : coloring)
    ++num[i];
  for (int g = 0; g < m; ++g)
    if (num[g] == 0) {
      int h = max_element(begin(num), end(num)) - begin(num);
      assert(h != g);
      for (int i = 0; i < n and num[g] < num[h]; ++i)
        if (coloring[i] == h)
          coloring[i] = g, --num[h], ++num[g];
    }
#ifndef NDEBUG
  assert((int)size(coloring) == n);
  assert(all_of(begin(coloring), end(coloring),
                [&](int i) { return inrange(i, 0, m - 1); }));
  for (int i = 0; i < m; ++i)
    assert(linearIn(coloring, i));
#endif
  Solution s;
  s.dd.setConsiderEqualDisp(nDupl >= n);
  s.populate(coloring);
  if (s.dispInt() < expectedD) {
    pr("dispInt: {}, expectedD: {}\n", s.dispInt(), expectedD);
  }
  assert(s.isComplete() and s.dispInt() >= expectedD);
  return s;
}
Solution solveUmdpFernandez(const Solution &lb, Timer t, bool verbose) {
  if (umdpFernandezDir == "down") {
    int step = umdpFernandezSteps == "linear" ? 1 : ub - lb.dispInt();
    int lastDNotOk = ub + 1;
    for (uint d = ub; d >= lb.dispInt() + 1 and not t.timedOut(); d -= step) {
      pr("#{}; d: {}/{}, step: {}\n", stats::ferIter + 1, d, ub, step);
      ++stats::ferIter;
      VI coloring;
      bool colorable = colorableTabuCol(d, coloring, t);
      if (colorable) {
        if (verbose)
          pr("d {}(r{})/{}: Ok, colorable with <= m.\n", d, Rd[d], Rd.size());
        if (step != 1) {
          step = floor(step / 2.0);
          d = lastDNotOk;
        } else {
          return getSolutionFromColoring(coloring, d);
        }
      } else if (verbose) {
        pr("d {}(r{})/{}: NOT colorable with <= m!\n", d, Rd[d], Rd.size());
        lastDNotOk = d;
      }
    }
    return lb;
  } else {
    assert(umdpFernandezDir == "up");
    int step = umdpFernandezSteps == "linear" ? 1 : ub - lb.dispInt() - 1;
    uint lastDOk = lb.dispInt();
    VI lastColoringOk, coloring;
    for (uint d = lb.dispInt() + 1; d <= ub and not t.timedOut(); d += step) {
      ++stats::ferIter;
      pr("#{}; d: {}/{}, step: {}\n", stats::ferIter + 1, d, ub, step);
      bool colorable = colorableTabuCol(d, coloring, t);
      if (colorable) {
        if (verbose)
          pr("d {}(r{})/{}: Ok, colorable with <= m.\n", d, Rd[d], Rd.size());
        lastDOk = d;
        lastColoringOk = coloring;
      } else {
        if (verbose)
          pr("d {}(r{})/{}: NOT colorable with <= m!\n", d, Rd[d], Rd.size());
        if (step != 1) {
          step = floor(step / 2.0);
          d = lastDOk;
        } else
          break;
      }
    }
    return (lastDOk == lb.dispInt())
               ? lb
               : getSolutionFromColoring(lastColoringOk, lastDOk);
  }
}
