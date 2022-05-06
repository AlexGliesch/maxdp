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
#include "ec.h"
#include "balvns.h"
#include "dynamicdispersion.h"
#include "solution.h"
#include "stats.h"
#include "vnsmove.h"
namespace ec {
VVI nodeGroupFailed;
VI failHist;
struct SimpleSolution {
  SimpleSolution() = default;
  SimpleSolution(Solution& s) { a = s.a, disp = s.dd.disp(); }
  void ejectBatch(const VI& P) {
    for (int i : P)
      a[i] = -1;
    if (a[disp.n1] == -1 or a[disp.n2] == -1) recomputeDispBF();
  }
  void eject(int i) {
    assert(a[i] != -1);
    a[i] = -1;
    if (i == disp.n1 or i == disp.n2) recomputeDispBF();
  }
  void insert(int i, int g) {
    assert(a[i] == -1);
    for (int u : duu[i]) {
      if (di[u][i] >= disp.val) break;
      if (a[u] == g) {
        disp = Disp(i, u);
        break;
      }
    }
    a[i] = g;
  }
  void recomputeDispBF([[maybe_unused]] bool computeAmt = false) {
    [[maybe_unused]] int oldDisp = disp.val;
    disp = Disp();
    for (auto [i, j] : R)
      if (a[i] == a[j] and a[i] != -1) {
        disp = Disp(i, j);
        break;
      }
    assert(disp.val >= oldDisp);
  }
  Disp disp;
  VI a;
  int numInserts = 0;
};
struct Conflicts {
  int u = -1;
  VVI L;
  VB feasible;
};
int ls(Solution& s, Timer t, bool verbose) {
  Timer lsTimer;
  int lsSteps = 0;
  while (not t.timedOut()) {
    ++stats::ecNumLsSteps;
    int n1 = s.dispN1(), n2 = s.dispN2();
    assert(di[n1][n2] == s.disp().val);
    assert(s.groupSize(s.a[n1]) >= 2);
    VNSMove mvSh(VNSMove::shift, NLI::min());
    for (int i : VI{n1, n2}) {
      static int g = 0;
      for (int gc = 0; gc < m and not t.timedOut(); ++gc, g = (g + 1) % m) {
        if (g == s.a[i]) continue;
        mvSh.considerDisp(s.shiftCostDisp(i, g), i, g);
        if (mvSh.valDisp > s.disp()) goto doShift;
      }
    }
    if (mvSh.valDisp > s.disp()) {
    doShift:
      ++stats::ecNumShifts;
      if (verbose)
        pr("shift {}({})->{}\n", mvSh.u, s.a[mvSh.u], mvSh.g, mvSh.valDisp);
      s.shift(mvSh.u, mvSh.g, true, true);
      s.checkCorrect(t);
      assert(s.disp() == mvSh.valDisp);
      ++lsSteps;
      continue;
    }
    VNSMove mvSwp(VNSMove::swap, NLI::min());
    for (int i : VI{n1, n2}) {
      static int j = 0;
      for (int jc = 0; jc < n and not t.timedOut(); ++jc, j = (j + 1) % n) {
        if (s.a[i] == s.a[j]) continue;
        mvSwp.considerDisp(s.swpCostDisp(i, j), i, j);
        if (mvSwp.valDisp > s.disp()) goto doSwp;
      }
    }
    if (mvSwp.valDisp > s.disp()) {
    doSwp:
      ++stats::ecNumSwaps;
      int g1 = s.a[mvSwp.u1], g2 = s.a[mvSwp.u2];
      if (verbose) pr("swap {}({})<>{}({})\n", mvSwp.u1, g1, mvSwp.u2, g2);
      s.swp(mvSwp.u1, mvSwp.u2, true, true);
      s.checkCorrect(t);
      assert(s.disp() == mvSwp.valDisp);
      ++lsSteps;
      continue;
    }
    break;
  }
  s.checkCorrect(t);
  stats::ecLsTime += lsTimer.elapsedSecs();
  return lsSteps;
}
void findConflicts(int u, int notG, const VB& fixed, const SimpleSolution& s,
                   Conflicts& c, int maxConf = NLI::max()) {
  assert(s.a[u] == -1);
  assert(maxConf >= 0);
  c.u = u;
  c.L.assign(m, VI());
  c.feasible.assign(m, true);
  c.feasible[notG] = false;
  for (int j : duu[u]) {
    if (di[j][u] > s.disp.val) break;
    int k = s.a[j];
    if (k == -1) continue;
    if (fixed[j]) c.feasible[k] = false;
    if ((int)size(c.L[k]) >= maxConf) c.feasible[k] = false;
    if (not c.feasible[k]) continue;
    c.L[k].push_back(j);
  }
}
void countConflicts(int u, int notG, const VB& fixed, const SimpleSolution& s,
                    VI& numConfsPerGroup) {
  numConfsPerGroup.assign(m, 0);
  numConfsPerGroup[notG] = NLI::max();
  for (int j : duu[u]) {
    if (di[j][u] > s.disp.val) break;
    int k = s.a[j];
    if (k == -1) continue;
    if (fixed[j]) numConfsPerGroup[k] = NLI::max();
    numConfsPerGroup[k] += (numConfsPerGroup[k] != NLI::max());
  }
}
void getGroupOrder(const Conflicts& c, VI& ordGr) {
  ordGr.resize(m);
  iota(begin(ordGr), end(ordGr), 0);
  if (groupOrder == ordRandom) {
    throw logic_error("Deprecated.");
    shuffle(begin(ordGr), end(ordGr), rng);
  } else if (groupOrder == ordConflicts) {
    sort(begin(ordGr), end(ordGr), [&](int i, int j) {
      return mt(size(c.L[i]), nodeGroupFailed[c.u][i], i) <
             mt(size(c.L[j]), nodeGroupFailed[c.u][j], j);
    });
  } else if (groupOrder == ordHistory) {
    sort(begin(ordGr), end(ordGr), [&](int i, int j) {
      return mt(nodeGroupFailed[c.u][i], size(c.L[i]), i) <
             mt(nodeGroupFailed[c.u][j], size(c.L[j]), j);
    });
  }
}
void sortNodes(VI& L, int g, const VB& fixed, const SimpleSolution& s) {
  if (nodeOrder == ExpansionOrder::ordRandom) {
    throw logic_error("Deprecated.");
    shuffle(begin(L), end(L), rng);
  } else if (nodeOrder == ExpansionOrder::ordConflicts) {
    static VVI nodeConfList;
    if (nodeConfList.empty()) nodeConfList.assign(n, VI(m));
    for (int u : L) {
      countConflicts(u, g, fixed, s, nodeConfList[u]);
      sort(begin(nodeConfList[u]), end(nodeConfList[u]));
    }
    sort(begin(L), end(L), [&](int u, int w) {
      return tie(nodeConfList[u], failHist[u]) >
             tie(nodeConfList[w], failHist[w]);
    });
  } else if (nodeOrder == ExpansionOrder::ordHistory) {
    sort(begin(L), end(L),
         [&](int u, int w) { return failHist[u] > failHist[w]; });
  }
}
int maxConflictsAllowed(int depth) {
  return depth > p3 ? 0 : (depth > p1 ? 1 : n - 1);
}
bool insert(SimpleSolution& s, int u, int notG, VB& fixed, int depth, Timer t) {
  ++s.numInserts;
  ++stats::ecNodeExp;
  assert(s.a[u] == -1);
  if (t.timedOut()) return false;
  int maxC = maxConflictsAllowed(depth);
  Conflicts c;
  findConflicts(u, notG, fixed, s, c, maxC);
  VI grOrd;
  getGroupOrder(c, grOrd);
  assert((int)size(grOrd) == m);
  int grTried = 0;
  for (int g : grOrd) {
    if (t.timedOut()) break;
    if (depth > p2 and grTried > 0)
      break;
    if (not c.feasible[g]) continue;
    auto& L = c.L[g];
    assert((int)size(L) <= maxC);
    ++grTried;
    if (L.empty()) {
      s.insert(u, g);
      fixed[u] = true;
      return true;
    }
    SimpleSolution sTmp = s;
    VB fxTmp = fixed;
    sTmp.ejectBatch(L);
    sTmp.insert(u, g);
    fxTmp[u] = true;
    assert(sTmp.disp.val >= s.disp.val);
    bool ok = true;
    sortNodes(L, g, fxTmp, sTmp);
    for (int w : L) {
      if (t.timedOut()) {
        ok = false;
        break;
      }
      ok = insert(sTmp, w, g, fxTmp, depth + 1, t);
      if (not ok) {
        ++failHist[w];
        break;
      }
    }
    if (ok) {
      assert(sTmp.disp.val >= s.disp.val);
      fixed = fxTmp;
      s = sTmp;
      return true;
    } else {
      ++nodeGroupFailed[u][g];
    }
  }
  return false;
}
bool improve(Solution& s, Timer t) {
  s.checkCorrect(t);
  VB fixed;
  SimpleSolution st;
  VI us;
  if (s.dd.considerEqualDisp and s.dd.nearSet.size() > 1) {
    VI numInc(n, 0);
    for (const auto& d : s.dd.nearSet) {
      ++numInc[d.n1], ++numInc[d.n2];
      us.push_back(d.n1);
      us.push_back(d.n2);
    }
    removeDuplicates(us);
    sort(begin(us), end(us), [&](int i, int j) {
      return mt(numInc[i], -failHist[i]) > mt(numInc[j], -failHist[i]);
    });
  } else {
    us = VI{s.dispN1(), s.dispN2()};
    sortNodes(us, s.a[us[0]], fixed, st);
  }
  for (int u : us) {
    if (t.timedOut()) break;
    fixed.assign(n, false);
    st = SimpleSolution(s);
    st.eject(u);
    Timer insertTimer;
    bool inserted = insert(st, u, s.a[u], fixed, 0, t);
    stats::ecInsertTime += insertTimer.elapsedSecs();
    if (inserted and st.disp > s.dd.disp()) {
      s.populate(st.a);
      assert(s.disp().val == st.disp.val);
      return true;
    }
  }
  return false;
}
void solve(Solution& s, int maxAlt, Timer t, bool verbose) {
  Timer ecTimer, lastIterTimer;
  static bool init = false;
  if (not init) {
    failHist.assign(n, 0);
    nodeGroupFailed.assign(n, VI(m, 0));
    init = true;
  }
  if (verbose) pr("Solving UMDP by EC...:\n");
  int numImproves = 0;
  Disp dBeforeLS, dBeforeEC;
  for (int i = 0; not t.timedOut(); ++i) {
    lastIterTimer.reset();
    ++stats::ecNumAlt;
    dBeforeLS = s.disp();
    Timer lsTimer;
    int lsSteps = ls(s, t, false);
    assert(s.disp() >= dBeforeLS);
    bool lsImproved = s.disp() > dBeforeLS;
    stats::ecLsAltImp += lsImproved;
    if (verbose)
      pr("#{} LS: r{} {} -> r{} {}, steps: {}, time: {}\n", i,
         Rd[dBeforeLS.val], dBeforeLS, s.disp().real(), s.disp(), lsSteps,
         lsTimer.elapsedSecs());
    if (s.disp().val > dBeforeLS.val) {
      stats::ecTtb = globalTimer.elapsedSecs();
      stats::ecIterTb = i;
    }
    if (i > 0 and s.disp().val > dBeforeEC.val) {
      ++numImproves;
      if (numImproves >= maxAlt) break;
    }
    if (s.dispOptimal()) break;
    Timer ecImprTimer;
    dBeforeEC = s.disp();
    bool ecImproved = improve(s, t);
    if (verbose)
      pr("#{} EC: r{} {}, visited: {}, time: {}\n", i, s.disp().real(),
         s.disp(), stats::ecNodeExp, ecImprTimer.elapsedSecs());
    if (s.disp().val > dBeforeEC.val) {
      stats::ecTtb = globalTimer.elapsedSecs();
      stats::ecIterTb = i;
    }
    assert(s.disp() >= dBeforeEC);
    if (s.dispOptimal() or not(lsImproved or ecImproved)) break;
  }
  stats::ecTimeLstEc = lastIterTimer.elapsedSecs();
  stats::ecTime += ecTimer.elapsedSecs();
  s.b.init(s);
  s.checkCorrect(t);
}
}
