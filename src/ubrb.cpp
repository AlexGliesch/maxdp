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
III b2b3(const VI& S) {
  using namespace ubkim;
  static VII e = subsetEdgeList(2);
  sortSubsetEdgeList(S, e);
  auto ret = [&](int i) {
    return III(di[S[e[i].first]][S[e[i].second]], e[i].first, e[i].second);
  };
  if (edgesAreDisjoint(0, 1, e)) return ret(1);
  if (edgesAreTriangle(0, 1, 2, e) or edgesAreDisjoint(0, 2, e) or
      edgesAreDisjoint(1, 2, e))
    return ret(2);
  int u = commonVertex(0, 1, e);
  assert(u == commonVertex(0, 2, e) and u == commonVertex(1, 2, e));
  for (uint i = 3; i < e.size(); ++i)
    if (commonVertex(i, 0, e) != u) return ret(i);
  return ret(e.size() - 1);
}
III b2(const VI& S) {
  using namespace ubkim;
  static VII e = subsetEdgeList(2);
  sortSubsetEdgeList(S, e);
  auto getNextDisjoint = [&](int e1, int e2) {
    for (int i = e1 + 1; i < e2; ++i)
      if (edgesAreDisjoint(e1, i, e)) return i;
    return e2;
  };
  int e2 = getNextDisjoint(0, e.size());
  assert(e2 != (int)e.size());
  for (int i = 1; i < e2; ++i)
    e2 = min(e2, getNextDisjoint(i, e2));
  int ans = di[S[e[e2].first]][S[e[e2].second]];
  return mt(ans, e[e2].first, e[e2].second);
}
III b3(const VI& S) {
  III ans(NLI::min(), -1, -1);
  for (uint i = 0; i < S.size(); ++i)
    for (uint j = i + 1; j < S.size(); ++j)
      for (uint k = j + 1; k < S.size(); ++k)
        ans =
            max(ans,
                min(III(di[S[i]][S[j]], i, j),
                    min(III(di[S[i]][S[k]], i, k), III(di[S[j]][S[k]], j, k))));
  return ans;
}
III urbs(const VI& S) {
  assert((int)S.size() == m + 2);
  return b2b3(S);
}
int LS([[maybe_unused]] VI& S, [[maybe_unused]] Timer t) {
  int i = 0, noImp = 0;
  auto [Svalue, m1, m2] = urbs(S);
  for (; noImp <= n; i = (i + 1) % n, ++noImp) {
    if (linearIn(S, i)) continue;
    for (int j : {m1, m2}) {
      if (t.timedOut()) return Svalue;
      int oldSj = S[j];
      S[j] = i;
      auto [nVal, nm1, nm2] = urbs(S);
      if (nVal < Svalue) {
        Svalue = nVal, m1 = nm1, m2 = nm2;
        noImp = 0;
        break;
      }
      S[j] = oldSj;
    }
  }
  return Svalue;
}
void computeUbrb(Timer t, bool verbose) {
  if (verbose) pr("Computing ub^RB...\n");
  ubrb = NLI::max();
  Timer ubrbTimer;
  double totUb = 0.0, totUbCons = 0.0;
  const int k = m;
  [[maybe_unused]] VI ssBest;
  [[maybe_unused]] vector<pair<int, VI>> kBest;
  for (int i = 0; i < n and not t.timedOut(); ++i) {
    VI S;
    S.push_back(i);
    while ((int)S.size() != m + 2) {
      int minDist = NLI::max(), bk = -1;
      for (int k = 0; k < n; ++k)
        if (not linearIn(S, k)) {
          int d = NLI::min();
          for (int l : S)
            d = max(d, di[l][k]);
          if (d < minDist) minDist = d, bk = k;
        }
      assert(bk != -1 and not linearIn(S, bk));
      S.push_back(bk);
    }
    int bound = get<0>(urbs(S));
    if (ubrbDoLS) {
      if ((int)kBest.size() < k or mp(bound, S) < kBest.back()) {
        kBest.emplace_back(bound, S);
        if ((int)kBest.size() > k) {
          sort(begin(kBest), end(kBest));
          kBest.erase(kBest.begin() + k, kBest.end());
        }
      }
    }
    totUbCons += bound;
    if (bound < ubrb) {
      ubrb = bound;
      stats::ubrbIterToBest = i + 1 + n;
      stats::ubrbTtb = globalTimer.elapsedSecs();
      ssBest = S;
    }
  }
  stats::ubrbMinCons = ubrb;
  if (ubrbDoLS) {
    for (int i = 0; i < k and not t.timedOut(); ++i) {
      int lsBound = LS(kBest[i].second, t);
      if (verbose)
        pr("ub^RB LS({}): {}, time: {}\n", i, lsBound, ubrbTimer.elapsedSecs());
      totUb += lsBound;
      if (lsBound < ubrb) {
        ubrb = lsBound;
        stats::ubrbIterToBest = i + 1 + n;
        stats::ubrbTtb = globalTimer.elapsedSecs();
        ssBest = kBest[i].second;
      }
    }
  }
  stats::ubrbAvg = totUb / double(k);
  stats::ubrbAvgCons = totUbCons / double(n);
  stats::ubrbTime = ubrbTimer.elapsedSecs();
  if (verbose)
    pr("ub^RB: {} (r{}), time: {}\n", ubrb, Rd[ubrb], stats::ubrbTime);
}
