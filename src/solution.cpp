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
#include "solution.h"
#include "main.h"
#include "ub.h"
void Solution::init() {
  a.assign(n, -1);
  ga.resize(m);
  numAssigned = 0;
  iiga.assign(n, -1);
  dd.init(*this);
  b.init(*this);
}
void Solution::populate(const VI& assigned) {
  a = assigned;
  ga.assign(m, VI());
  iiga.assign(n, -1);
  numAssigned = 0;
  for (int i = 0; i < n; ++i)
    if (a[i] != -1) {
      ++numAssigned;
      iiga[i] = ga[a[i]].size();
      ga[a[i]].push_back(i);
    }
  dd.init(*this);
  b.init(*this);
}
Solution Solution::trivial() {
  Solution s;
  s.dd.setConsiderEqualDisp(nDupl >= n);
  s.init();
  for (int i = 0; i < n; ++i)
    s.shift(i, i % m, false, false);
  assert(s.isComplete());
  dynamicDispBf(s, true);
  balanceBF(s, true);
  return s;
}
Solution Solution::random() {
  pr("---- We should not be using random initial solutions. Consider using a "
     "greedy procedure.\n");
  Solution s;
  s.dd.setConsiderEqualDisp(nDupl >= n);
  s.init();
  VI indn = iotaed(n);
  shuffle(begin(indn), end(indn), rng);
  for (int i = 0; i < m; ++i)
    s.shift(indn[i], i, false, false);
  for (int i = m; i < n; ++i)
    s.shift(indn[i], randInt(0, m - 1), false, false);
  assert(s.isComplete());
  dynamicDispBf(s, true);
  balanceBF(s, true);
  return s;
}
void Solution::shift(int u, int g, bool updBal, bool updDsp) {
  int p = a[u];
  if (p == g) return;
  if (updDsp) {
    dd.shift(*this, u, g, p);
  }
  if (updBal) {
    b.shift(*this, u, g, p);
  }
  a[u] = g;
  if (p != -1) {
    iiga[ga[p].back()] = iiga[u];
    swap(ga[p][iiga[u]], ga[p].back());
    ga[p].pop_back();
  } else {
    ++numAssigned;
  }
  if (g != -1) {
    ga[g].push_back(u);
    iiga[u] = ga[g].size() - 1;
  } else {
    --numAssigned;
    iiga[u] = -1;
  }
}
void Solution::swp(int u1, int u2, bool updBal, bool updDsp) {
  if (updDsp) dd.swp(*this, u1, u2);
  if (updBal) b.swp(*this, u1, u2);
  int g1 = a[u1], g2 = a[u2];
  int &i1 = iiga[u1], &i2 = iiga[u2];
  assert(i1 != -1 and i2 != -1);
  swap(a[u1], a[u2]);
  swap(ga[g1][i1], ga[g2][i2]);
  swap(i1, i2);
}
int Solution::randomWalk(int numMoves, int minDisp) {
  static VI indn = iotaed(n), indm = iotaed(m);
  shuffle(begin(indn), end(indn), rng), shuffle(begin(indm), end(indm), rng);
  [[maybe_unused]] int movesTested = 0;
  int movesDone = 0, round = 0;
  for (int i = 0; round < 2 and movesDone < numMoves; i = (i + 1) % m) {
    int g = indm[i];
    for (int j = 0; j < n and movesDone < numMoves; ++j) {
      int u = indn[j];
      ++movesTested;
      if (a[u] != g and shiftCostDisp(u, g).val >= minDisp)
        shift(u, g, true, true), ++movesDone;
    }
    if (i == m - 1) ++round;
  }
  return movesDone;
}
void Solution::checkCorrect([[maybe_unused]] Timer t) const {
  if constexpr (dbgSol or dbgBal or dbgDisp)
    if (empty(a)) return;
  if constexpr (dbgSol) {
    [[maybe_unused]] int numAssigned = (n - count(begin(a), end(a), -1));
    assert(numAssigned == this->numAssigned);
    for (int g = 0; g < m; ++g)
      assert(groupSize(g) >= 1);
    vector<vector<int>> ga(m), ga2 = this->ga;
    for (int u = 0; u < n; ++u)
      if (a[u] != -1) ga[a[u]].push_back(u);
    for (auto& g : ga)
      sort(begin(g), end(g));
    for (auto& g : ga2)
      sort(begin(g), end(g));
    assert(ga == ga2);
  }
  if constexpr (dbgBal) {
    Solution o = *this;
    balanceBF(o, true);
    assert(ff(bal()) == ff(o.bal()));
    for (int g = 0; g < m; ++g)
      assert(ffuple(b.gb[g], b.gw[g]) == ffuple(o.b.gb[g], o.b.gw[g]));
  }
  if constexpr (dbgDisp) {
    VDisp gdBF;
    Disp globalDisp;
    bool cEqD = dd.considerEqualDisp;
    for (int g = 0; g < m; ++g) {
      Disp d;
      assert(groupSize(g) == (int)size(ga[g]));
      for (uint i = 0; i < size(ga[g]); ++i)
        for (uint j = i + 1; j < size(ga[g]); ++j)
          d.consider(Disp(ga[g][i], ga[g][j]), cEqD);
      globalDisp.consider(d, cEqD);
      gdBF.push_back(d);
      assert(gdBF[g].val == dd.gd[g].val);
    }
    for (auto& d : dd.nearSet)
      if (inrange(d.n1, 0, n - 1) and inrange(d.n2, 0, n - 1))
        assert(a[d.n1] == a[d.n2]);
    if (not cEqD)
      assert(dd.nearSet.size() == 1 and dd.gd[dd.gds[0]] == dd.nearSet[0]);
    assert(dd.disp().val == globalDisp.val);
    assert(not cEqD or linearIn(dd.nearSet, globalDisp));
    if (t.timedOut()) return;
    vector<VDisp> gid;
    vector<vector<i16>> numRd;
    for (int g = 0; g < m and not t.timedOut(); ++g) {
      VDisp gidG;
      vector<i16> numRdG(Rd.size(), 0);
      assert(not cEqD or size(dd.numRd[g]) == size(numRdG));
      for (uint i = 0; i < size(ga[g]); ++i)
        for (uint j = i + 1; j < size(ga[g]); ++j) {
          gidG.emplace_back(ga[g][i], ga[g][j]);
          ++numRdG[gidG.back().val];
        }
      sort(begin(gidG), end(gidG));
      if (cEqD) {
        numRdG[0] = size(ga[g]);
        for (Disp& d : gidG)
          d.amt = numRdG[d.val];
        assert(numRdG == dd.numRd[g]);
      }
      assert(equal(begin(gidG), begin(gidG) + size(dd.gid[g]), begin(dd.gid[g]),
                   end(dd.gid[g])));
      gid.emplace_back(move(gidG));
      numRd.emplace_back(move(numRdG));
    }
    VI amt;
    amt.assign(Rd.size() + 1, 0);
    for (int u = 0; u < n and not t.timedOut(); ++u)
      for (int g = 0; g < m; ++g) {
        amt.assign(Rd.size() + 1, 0);
        VDisp utgd;
        for (int w : ga[g])
          if (w != u) {
            Disp duw(u, w);
            if (cEqD) ++amt[duw.val];
            utgd.push_back(duw);
          }
        sort(begin(utgd), end(utgd));
        while (utgd.size() > 2 and utgd.back().val != utgd[1].val)
          utgd.pop_back();
        if (cEqD)
          for (auto& i : utgd)
            i.amt = amt[i.val];
        VDisp ddUtgd = dd.utgd[u][g];
        sort(begin(ddUtgd), end(ddUtgd));
        assert(mp(utgd[0], utgd[1]) == mp(ddUtgd[0], ddUtgd[1]));
      }
    if (cEqD and not t.timedOut()) {
      VDisp nearSet;
      for (int g = 0; g < m; ++g)
        for (const Disp& d : gid[g]) {
          if (d.val == globalDisp.val)
            nearSet.push_back(d);
          else
            break;
        }
      for (auto& d : nearSet)
        d.amt = size(nearSet);
      sort(begin(nearSet), end(nearSet));
      assert(is_sorted(begin(dd.nearSet), end(dd.nearSet)));
      assert(nearSet == dd.nearSet);
    }
  }
}
void Solution::writeToFile(const string& filename) {
  ofstream f(filename);
  if (f.fail()) throw logic_error("invalid output filename.");
  f << format("{} {} {}\n", instName, n, m);
  for (int g = 0; g < m; ++g) {
    f << groupSize(g) << ' ';
    for (int u : ga[g])
      f << u << ' ';
    f << '\n';
  }
}
Solution Solution::readFromFile(const string& filename) {
  ifstream f(filename);
  if (f.fail()) throw logic_error("invalid input filename.");
  string iNameI;
  int nI, mI;
  f >> iNameI >> nI >> mI;
  if (mt(nI, mI) != mt(n, m))
    throw logic_error(
        format("inconsistent solution size in input file; {} != {}", mt(nI, mI),
               mt(n, m)));
  VI a(n, -1);
  for (int g = 0; g < m; ++g) {
    int sz;
    f >> sz;
    for (int i = 0; i < sz; ++i) {
      int u;
      f >> u;
      a[u] = g;
    }
  }
  if (count(begin(a), end(a), -1) != 0)
    throw logic_error("input solution is incomplete");
  Solution s;
  s.populate(a);
  pr("Read solution from file: {}\n", s);
  return s;
}
