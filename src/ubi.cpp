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
#include "main.h"
#include "ub.h"
namespace {
struct UBISearchSolution {
  UBISearchSolution(int i) {
    S.resize(m + 1);
    copy(begin(duu[i]), begin(duu[i]) + m + 1, begin(S));
    compute();
  }
  void compute() {
    compMaxDist(maxDist, mi, mj, -1);
    compMaxDist(maxDist_i, mi_i, mj_i, mi);
    compMaxDist(maxDist_j, mi_j, mj_j, mj);
  }
  int valSwap(int i, int u) {
    int val = (i == mi ? maxDist_i : (i == mj ? maxDist_j : maxDist));
    for (int j = 0; j < m + 1; ++j) {
      if (j == i) continue;
      if (S[j] != u) {
        val = max(val, di[S[j]][u]);
      } else {
        val = NLI::max();
        break;
      }
    }
    return val;
  }
  void swap(int i, int u) {
    assert(i >= 0 and i < m + 1);
    S[i] = u;
    compute();
    assert((int)size(S) == m + 1);
  }
  void compMaxDist(int& maxDist, int& mi, int& mj, int ignore) {
    mi = mj = -1;
    maxDist = NLI::min();
    for (int i = 0; i < m + 1; ++i)
      if (i != ignore)
        for (int j = i + 1; j < m + 1; ++j)
          if (j != ignore) {
            if (S[i] == S[j]) {
              maxDist = NLI::max(), mi = i, mj = j;
              return;
            }
            if (di[S[i]][S[j]] > maxDist)
              mi = i, mj = j, maxDist = di[S[i]][S[j]];
          }
  }
  int maxDist, mi, mj;
  int mi_j, mj_j, maxDist_j;
  int mi_i, mj_i, maxDist_i;
  VI S;
};
}
void computeUbi(Timer t, bool verbose) {
  if (verbose) pr("Computing ub^I...\n");
  Timer ubiTimer;
  ubi = NLI::max();
  double totUb = 0.0, totUbCons = 0.0;
  for (int k = 0; k < n and not t.timedOut(); ++k) {
    UBISearchSolution S(k);
    [[maybe_unused]] int maxDistB = S.maxDist;
    totUbCons += maxDistB;
    stats::ubiMinCons = min(stats::ubiMinCons, maxDistB);
    bool improved = true;
    while (improved and not t.timedOut()) {
      improved = false;
      for (int i : {S.mi, S.mj}) {
        if (improved) break;
        for (int u = 0; u < n; ++u)
          if (S.valSwap(i, u) < S.maxDist) {
            improved = true;
            S.swap(i, u);
            break;
          }
      }
    }
    totUb += S.maxDist;
    assert(S.maxDist <= maxDistB);
    if (S.maxDist < ubi) {
      ubi = S.maxDist;
      S.S.erase(S.S.begin() + S.mi);
      ubiSubset = S.S;
      stats::ubiIterToBest = k + 1;
      stats::ubiTtb = globalTimer.elapsedSecs();
    }
  }
  stats::ubiAvgCons = totUbCons / double(n);
  stats::ubiAvg = totUb / double(n);
  stats::ubiTime = ubiTimer.elapsedSecs();
  if (t.timedOut()) {
    if (verbose) pr("Time limit ({}s) reached. Exiting.\n", t.tmLim);
  }
  if (verbose)
    pr("ub^I: {} (r{}), time: {}\n", ubi,
       inrange(ubi, 0, (int)size(Rd) - 1) ? Rd[ubi] : -1.0, stats::ubiTime);
}
