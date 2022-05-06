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
#include "dynamicdispersion.h"
#include "solution.h"
namespace {
void gidTrim(vector<Disp>& gid) {
  if (size(gid) < 2) return;
  uint i = 1;
  for (; i < size(gid); ++i) {
    bool isfree = false;
    for (uint j = 0; j < i and not isfree; ++j)
      isfree = not gid[j].intersects(gid[i]);
    if (isfree ) {
      int k = i;
      while (gid[i].val == gid[k].val and i < size(gid))
        ++i;
      --i;
      break;
    }
  }
  assert(i <= size(gid));
  if (i < size(gid)) gid.erase(begin(gid) + i + 1, end(gid));
}
void gidBF(VDisp& gid, const vector<i16>& numRd, const VI& ga, int ignore,
           bool considerEqualDisp) {
  gid.clear();
  Disp v1, v2;
  for (uint i = 0; i < size(ga); ++i)
    if (ga[i] != ignore)
      for (uint j = i + 1; j < size(ga); ++j)
        if (ga[j] != ignore) {
          Disp dij(ga[i], ga[j]);
          if (dij < v1) {
            if (not dij.intersects(v1))
              v2 = v1;
            else if (dij.intersects(v2)) {
              v2 = Disp();
              for (const Disp& v : gid)
                if (not dij.intersects(v) and v < v2) v2 = v;
            }
            v1 = dij;
          } else if (dij < v2 and not dij.intersects(v1)) {
            v2 = dij;
          }
          if (dij.val <= v2.val) gid.push_back(dij);
          assert((v1.val == NLI::max() and v2.val == NLI::max()) or
                 (not v1.intersects(v2) and v1.val <= v2.val));
        }
  sort(begin(gid), end(gid));
  gidTrim(gid);
  if (considerEqualDisp and size(numRd) == size(Rd))
    for (Disp& d : gid)
      if (inrange(d.val, 0, (int)size(Rd) - 1)) d.amt = numRd[d.val];
}
Disp distToGr(const Solution& s, int u, int g, int ignore,
              bool considerEqualDisp) {
  auto& ud = s.dd.utgd[u][g];
  if (considerEqualDisp and ud[0].val == ud[1].val) {
    bool red = ud[0].contains(ignore);
    Disp r = ud[red];
    for (uint i = 1; i < ud.size() and not red; ++i)
      red = ud[i].contains(ignore);
    r.amt -= red;
    return r;
  } else {
    return ud[0].contains(ignore) ? ud[1] : ud[0];
  }
}
Disp grInDisp(const Solution& s, int g, int ignore, bool considerEqualDisp) {
  if (considerEqualDisp) {
    int igCount = 0;
    for (uint i = 0; i < size(s.dd.gid[g]); ++i) {
      const Disp& w = s.dd.gid[g][i];
      if (i > 0 and w.val != s.dd.gid[g][i - 1].val) igCount = 0;
      if (w.contains(ignore)) {
        ++igCount;
      } else {
        Disp r = w;
        for (uint j = i + 1;
             j < size(s.dd.gid[g]) and s.dd.gid[g][j].val == w.val; ++j)
          igCount += s.dd.gid[g][j].contains(ignore);
        r.amt -= igCount;
        return r;
      }
    }
  } else {
    for (const Disp& w : s.dd.gid[g])
      if (not w.contains(ignore)) return w;
  }
  return EmptyDisp;
}
}
void DynamicDispersion::init(Solution& s) {
  gd.resize(m);
  gds.resize(m);
  gid.resize(m);
  utgd.assign(n, vector<VDisp>(m));
  numRd.resize(m);
  if (considerEqualDisp) {
    numRd.assign(m, vector<i16>(Rd.size(), 0));
  }
  dynamicDispBf(s, true);
}
Disp DynamicDispersion::shiftCost(const Solution& s, int u, int g) const {
  int p = s.a[u];
  const bool e = considerEqualDisp;
  if (p == g) return s.disp();
  if (e) {
    Disp d = g == -1 ? Disp() : gd[g];
    d.consider(distToGr(s, u, g, -1, e), e);
    if (p != -1) d.consider(grInDisp(s, p, u, e), e);
    for (int i : gds)
      if (i != p and i != g) {
        d.consider(gd[i], e);
      }
    return d;
  } else {
    Disp dg = g == -1 ? Disp() : min(gd[g], distToGr(s, u, g, -1, e));
    Disp dp = p == -1 ? Disp() : grInDisp(s, p, u, e);
    Disp d = min(dg, dp);
    for (int i : gds)
      if (i != p and i != g) return min(d, gd[i]);
  }
  assert(false);
  return {};
}
Disp DynamicDispersion::swpCost(const Solution& s, int u1, int u2) const {
  int g1 = s.a[u1], g2 = s.a[u2];
  bool e = considerEqualDisp;
  if (e) {
    Disp d = grInDisp(s, g1, u1, e);
    d.consider(distToGr(s, u2, g1, u1, e), e);
    d.consider(grInDisp(s, g2, u2, e), e);
    d.consider(distToGr(s, u1, g2, u2, e), e);
    for (int i : gds)
      if (i != g1 and i != g2) d.consider(gd[i], e);
    return d;
  } else {
    Disp d1 = min(grInDisp(s, g1, u1, e), distToGr(s, u2, g1, u1, e));
    Disp d2 = min(grInDisp(s, g2, u2, e), distToGr(s, u1, g2, u2, e));
    Disp d = min(d1, d2);
    for (int i : gds)
      if (i != g1 and i != g2) return min(d, gd[i]);
  }
  return Disp{};
}
void DynamicDispersion::shift(const Solution& s, int u, int g, int p) {
  if (considerEqualDisp) {
    assert(inrange(g, -1, m - 1) and inrange(p, -1, m - 1));
    if (g != -1) {
      for (int w : s.ga[g])
        ++numRd[g][di[u][w]];
      ++numRd[g][0];
    }
    if (p != -1)
      for (int w : s.ga[p])
        --numRd[p][di[u][w]];
  }
  if (p != -1) gidLose(p, u, s.ga[p]);
  if (g != -1) gidGain(g, u, s.ga[g], -1);
  if (g != -1) gd[g] = gid[g].empty() ? Disp() : gid[g][0];
  if (p != -1) gd[p] = gid[p].empty() ? Disp() : gid[p][0];
  sortGds();
  for (int w = 0; w < n; ++w)
    if (w != u) {
      if (g != -1)
        utgdGain(g, u, w);
      if (p != -1)
        utgdLose(p, u, w, s.ga[p]);
    }
  if (considerEqualDisp)
    rcmpNearSet();
  else
    nearSet = VDisp(1, gd[gds[0]]);
}
void DynamicDispersion::swp(const Solution& s, int u1, int u2) {
  int g1 = s.a[u1], g2 = s.a[u2];
  if (considerEqualDisp) {
    for (int w : s.ga[g1]) {
      if (w != u1) ++numRd[g1][di[u2][w]];
      --numRd[g1][di[u1][w]];
    }
    for (int w : s.ga[g2]) {
      if (w != u2) ++numRd[g2][di[u1][w]];
      --numRd[g2][di[u2][w]];
    }
    ++numRd[g1][0];
    ++numRd[g2][0];
  }
  gidGainLose(g1, u2, u1, s.ga[g1]);
  gidGainLose(g2, u1, u2, s.ga[g2]);
  gd[g1] = gid[g1].empty() ? Disp() : gid[g1][0];
  gd[g2] = gid[g2].empty() ? Disp() : gid[g2][0];
  sortGds();
  for (int w = 0; w < n; ++w) {
    utgdGainLose(g1, u2, u1, w, s.ga[g1]);
    utgdGainLose(g2, u1, u2, w, s.ga[g2]);
  }
  if (considerEqualDisp)
    rcmpNearSet();
  else
    nearSet = VDisp(1, gd[gds[0]]);
}
void DynamicDispersion::rcmpNearSet() {
  int curDisp = min_element(begin(gd), end(gd))->val;
  nearSet.clear();
  for (int g = 0; g < m; ++g) {
    for (const Disp& i : gid[g]) {
      if (i.val == curDisp)
        nearSet.push_back(i);
      else
        break;
    }
  }
  for (auto& d : nearSet)
    d.amt = size(nearSet);
  sort(begin(nearSet), end(nearSet));
}
void DynamicDispersion::gidGain(int g, int u, const VI& ga, int ignore) {
  auto last = empty(gid[g]) ? Disp() : gid[g].back();
  for (int w : ga) {
    if (w != ignore and w != u and Disp(u, w).val <= last.val)
      gid[g].emplace_back(u, w);
  }
  if (size(gid[g]) and gid[g].back() != last) {
    if (considerEqualDisp)
      for (auto& d : gid[g])
        d.amt = numRd[g][d.val];
    sort(begin(gid[g]), end(gid[g]));
    gidTrim(gid[g]);
  }
}
void DynamicDispersion::gidLose(int g, int u, const VI& ga) {
  if (any_of(begin(gid[g]), end(gid[g]),
             [&](auto& d) { return d.contains(u); })) {
    assert(not empty(gid[g]));
    assert((int)size(numRd) == m);
    gidBF(gid[g], numRd[g], ga, u, considerEqualDisp);
  }
}
void DynamicDispersion::gidGainLose(int g, int uGain, int uLose, const VI& ga) {
  gidLose(g, uLose, ga);
  gidGain(g, uGain, ga, uLose);
}
void DynamicDispersion::utgdGain(int g, int u, int w) {
  udConsider(utgd[w][g], Disp(u, w), considerEqualDisp);
}
void DynamicDispersion::utgdLose(int g, int u, int w, const VI& ga) {
  auto& ud = utgd[w][g];
  if (ud.size() < 2 or ud[0].contains(u) or ud[1].contains(u)) {
    ud.clear();
    for (int v : ga)
      if (v != u and v != w) udConsider(ud, Disp(w, v), considerEqualDisp);
  } else if (considerEqualDisp) {
    auto it = find_if(ud.begin(), ud.end(),
                      [&](const Disp& d) { return d.contains(u); });
    if (it != ud.end()) {
      int v = it->val;
      ud.erase(it);
      for (auto& d : ud)
        if (d.val == v) --d.amt;
    }
  }
}
void DynamicDispersion::utgdGainLose(int g, int uGain, int uLose, int w,
                                     const VI& ga) {
  if (w != uLose) utgdLose(g, uLose, w, ga);
  if (w != uGain) utgdGain(g, uGain, w);
}
Disp dynamicDispBf(Solution& s, bool updateS) {
  bool cEqD = s.dd.considerEqualDisp;
  if (updateS) {
    assert((int)size(s.dd.gid) == m);
    assert((int)size(s.dd.gd) == m);
    auto& ga = s.ga;
    auto& gid = s.dd.gid;
    [[maybe_unused]] auto& numRd = s.dd.numRd;
    [[maybe_unused]] auto& nearSet = s.dd.nearSet;
    nearSet.clear();
    for (int g = 0; g < m; ++g) {
      if (empty(ga[g])) continue;
      if (cEqD) {
        assert(size(numRd[g]) == size(Rd));
        fill(begin(numRd[g]), end(numRd[g]), 0);
        for (int i = 0; i < (int)size(ga[g]); ++i)
          for (int j = i + 1; j < (int)size(ga[g]); ++j) {
            assert(ga[g][i] != ga[g][j]);
            Disp dij(ga[g][i], ga[g][j]);
            ++numRd[g][dij.val];
            if (nearSet.empty() or dij.val <= nearSet[0].val) {
              if (size(nearSet) and dij.val < nearSet[0].val) nearSet.clear();
              nearSet.push_back(dij);
            }
          }
        numRd[g][0] = size(s.ga[g]);
        for (Disp& d : nearSet)
          d.amt = size(nearSet);
      }
      gidBF(gid[g], numRd[g], ga[g], -1, cEqD);
      assert((int)size(s.dd.gd) > g);
      s.dd.gd[g] = gid[g].empty() ? Disp() : gid[g][0];
    }
    sort(begin(nearSet), end(nearSet));
    auto& utgd = s.dd.utgd;
    assert((int)size(utgd) == n and
           all_of(begin(utgd), end(utgd),
                  [&](auto& x) { return (int)size(x) == m; }));
    for (int u = 0; u < n; ++u)
      for (int g = 0; g < m; ++g) {
        auto& ud = utgd[u][g];
        ud.clear();
        for (int w : s.ga[g])
          if (w != u) udConsider(ud, Disp(w, u), cEqD);
      }
    assert((int)size(s.dd.gds) == m);
    iota(begin(s.dd.gds), end(s.dd.gds), 0);
    s.dd.sortGds();
    if (not cEqD) nearSet = VDisp(1, s.dd.gd[s.dd.gds[0]]);
    return s.dd.disp();
  } else {
    Disp d;
    for (int g = 0; g < m; ++g) {
      Disp v;
      for (uint i = 0; i < size(s.ga[g]); ++i)
        for (uint j = i + 1; j < size(s.ga[g]); ++j)
          v.consider(Disp(s.ga[g][i], s.ga[g][j]), cEqD);
      d.consider(v, cEqD);
    }
    if (d.n1 == -1 or d.n2 == -1) d.val = NLI::min(), d.amt = 0;
    return d;
  }
}
