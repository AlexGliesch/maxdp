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
int ukValue(const VI& S, VII& e) {
  throw logic_error(
      "Deprecated; this function considers only edges where all previous "
      "edges are disjoint. Remove later? Code might be useful somewhere.");
  TIME_BLOCK("ukValue");
  using namespace ubkim;
  int k = S.size() - m;
  assert(k >= 1 and k <= m);
  sortSubsetEdgeList(S, e);
  unordered_set<int> vs;
  for (int i = 0; i < k; ++i)
    vs.insert(e[i].first), vs.insert(e[i].second);
  if ((int)vs.size() == k * 2) {
    return di[S[e[k - 1].first]][S[e[k - 1].second]];
  } else if ((int)vs.size() < k * 2 - 2) {
    return -1;
  } else {
    int r1 = -1, r2 = -1;
    if ((int)vs.size() == k * 2 - 2) {
      int num0 = 0, num1 = 0, num2 = 0;
      for (int i = 0; i < k; ++i) {
        int num = 0;
        for (int j = 0; j < k; ++j)
          num += not edgesAreDisjoint(i, j, e);
        if (num == 2) r1 = i;
        num0 += (num == 0), num1 += (num == 1), num2 += (num == 2);
      }
      if (not(num0 == k - 3 and num1 == 2 and num2 == 1)) {
        return -1;
      }
    } else {
      for (int i = 0; i < k and r1 == -1; ++i)
        for (int j = i + 1; j < k and r1 == -1; ++j)
          if (not edgesAreDisjoint(i, j, e)) r1 = i, r2 = j;
    }
    bool kp1Dis = true;
    bool uRem = false;
    for (int i = 0; i < k and kp1Dis; ++i) {
      bool dis = edgesAreDisjoint(i, k, e);
      if (not dis and not uRem and (i == r1 or i == r2))
        uRem = true;
      else
        kp1Dis = dis;
    }
    if (kp1Dis) {
      return di[S[e[k].first]][S[e[k].second]];
    } else {
      return -1;
    }
  }
}
int ukValueNew(const VI& S, VII& e) {
  using namespace ubkim;
  sortSubsetEdgeList(S, e);
  int n = size(S);
  int r = n % m;
  int b1 = floor(double(n) / double(m)), b2 = ceil(double(n) / double(m));
  b1 = (b1 == 0 ? 0 : (b1 * (b1 - 1)) / 2);
  assert(b2 != 0);
  b2 = (b2 * (b2 - 1)) / 2;
  int i = (m - r) * b1 + r * b2;
  assert(inrange(i, 1, (int)e.size()));
  return di[S[e[i - 1].first]][S[e[i - 1].second]];
}
void computeUbk(int k, Timer t, bool verbose) {
  using namespace ubkim;
  VI deg, order;
  VI S(m + k);
  [[maybe_unused]] int itb = -1;
  int numImp = 0;
  ubsim::computeDeg(ubk, deg, order);
  auto e = subsetEdgeList(k);
  VB tried(n, false);
  int j = 0;
  int iter = 0;
  const bool fewer = false;
  while (iter < (fewer ? (m + k) * 2 : n) and not t.timedOut()) {
    for (; j < n; ++j)
      if (not tried[order[j]]) break;
    int i = order[j];
    tried[i] = true;
    ++iter;
    {
      ubsim::initialSubsetGreedy(S, i, ubk, deg);
      ubsim::optimizeSubset(S, ubk);
    }
    int uk = ukValueNew(S, e);
    if (uk != -1 and uk < ubk) {
      ubk = uk;
      itb = iter;
      ++numImp;
      ubsim::computeDeg(ubs, deg, order), j = 0;
    }
  }
  if (verbose) pr("UB^k, k: {}: {} (r{})\n", k, ubk, Rd[ubk]);
}
