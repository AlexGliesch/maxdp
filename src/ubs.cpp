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
#include "color.h"
#include "ub.h"
VI computeUbs(Timer t, bool verbose) {
  if (verbose) pr("Computing ub^S...\n");
  if (ubAlg == "all") {
    if (ubi == NLI::max()) computeUbi(t, verbose);
    if (ubrb == NLI::max()) computeUbrb(t, verbose);
    ubs = min(ubi, ubrb);
  } else if (ubAlg == "ubsblind") {
    ubs = Rd.size();
  } else {
    assert(ubAlg == "ubs");
    if (ubi == NLI::max()) computeUbi(t, verbose);
    ubs = ubi;
  }
  Timer ubsTimer;
  if (ubsSigma <= 0) throw ubsSigma;
  ubsSubsetSize =
      sigmaFixed ? m + ubsSigma : min(n, max(m + 3, (int)floor(ubsSigma * m)));
  VI ss(ubsSubsetSize), deg, order;
  [[maybe_unused]] VI ssBest;
  ubsim::computeDeg(ubs, deg, order);
  VB tried(n, false);
  int j = 0;
  while (stats::ubsIter < (ubsFewer ? ubsSubsetSize * 2 : n) and
         not t.timedOut()) {
    for (; j < n; ++j)
      if (not tried[order[j]]) break;
    int i = order[j];
    tried[i] = true;
    ++stats::ubsIter;
    ubsim::initialSubsetGreedy(ss, i, ubs, deg);
    [[maybe_unused]] int ssEdges = ubsim::optimizeSubset(ss, ubs);
    if (verbose)
      pr("#{}: Trying to color subset of node {}, with {} nodes and {} edges "
         "at ub^S = {}... ",
         stats::ubsIter, i, ubsSubsetSize, ssEdges, ubs);
    assert(ssEdges == ubsim::numEdges(ss, ubs));
    constexpr bool useTLOnColoring = false;
    constexpr double coloringTL = 2.0;
    Timer colorTimer(
        min(t.secsLeft(), useTLOnColoring ? coloringTL : NLD::max()));
    bool isLb = isLbColor(ubs, ss, colorTimer);
    if (verbose) {
      if (colorTimer.timedOut())
        pr("timed out at {}s!\n", coloringTL);
      else
        pr("{}ok!\n", isLb ? "" : "not ");
    }
    if (not isLb and not colorTimer.timedOut()) {
      if (verbose) pr("Not colorable, binary searching for new ub...\n");
      int lo = 0, hi = ubs - 1, mid = ubs;
      [[maybe_unused]] int old = ubs;
      while (lo <= hi) {
        ++stats::ubsBsIter;
        if (t.timedOut()) break;
        mid = (lo + hi) / 2;
        if (isLbColor(mid, ss, t)) {
          if (verbose)
            pr("Subset from {} of size {} IS {}-colorable at ub^S = {}(r{})\n",
               i, ubsSubsetSize, m, mid, Rd[mid]);
          lo = mid + 1;
        } else {
          if (verbose)
            pr("Subset from {} of size {} IS NOT {}-colorable at ub^S = "
               "{}(r{})\n",
               i, ubsSubsetSize, m, mid, Rd[mid]);
          hi = mid - 1;
          ubs = hi;
        }
      }
      ssBest = ss;
      if (ubs == old) --ubs;
      assert(ubs <= old);
      ubsim::computeDeg(ubs, deg, order), j = 0;
      stats::ubsIterToBest = stats::ubsIter;
      ++stats::ubsNumImproves;
      stats::ubsTtb = globalTimer.elapsedSecs();
    }
  }
  stats::ubsTime = ubsTimer.elapsedSecs();
  if (t.timedOut()) {
    if (verbose) pr("Time limit ({}s) reached. Exiting.\n", t.tmLim);
  }
  if (verbose)
    pr("ub^S: {} (r{}), subsets tried: {}, calls to coloring: {}, iters. to "
       "best: {}, num. improves: {}, time: {}\n\n",
       ubs, Rd[ubs], stats::ubsIter, stats::ubsNumColorCalls,
       stats::ubsIterToBest, stats::ubsNumImproves, stats::ubsTime);
  return (ubs < min(ubi, ubc)) ? move(ss) : VI();
}
