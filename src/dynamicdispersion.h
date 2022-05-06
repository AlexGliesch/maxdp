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
#pragma once
       
#include "main.h"
#include "util.h"
struct Disp {
  Disp() : n1(-1), n2(-1), amt(1), val(NLI::max()) {}
  Disp(int n1, int n2)
      : n1(min(n1, n2)), n2(max(n1, n2)), amt(1), val(di[n1][n2]) {}
  bool operator<(const Disp& o) const {
    return mt(val, -amt, n1, n2) < mt(o.val, -o.amt, o.n1, o.n2);
  }
  bool operator>(const Disp& o) const { return o < *this; }
  bool operator>=(const Disp& o) const { return not(*this < o); }
  bool operator<=(const Disp& o) const { return not(*this > o); }
  bool operator==(const Disp& o) const {
    return mt(val, amt, n1, n2) == mt(o.val, o.amt, o.n1, o.n2);
  }
  bool operator!=(const Disp& o) const { return not(*this == o); }
  bool contains(int u) const { return n1 == u or n2 == u; }
  bool intersects(const Disp& o) const {
    return contains(o.n1) or contains(o.n2);
  }
  double real() const {
    return inrange(val, 0, (int)Rd.size() - 1) ? Rd[val] : -1.0;
  }
  void consider(const Disp& o, bool considerEqualDisp) {
    if (considerEqualDisp) {
      if (o.val == val)
        amt += o.amt;
      else if (o.val < val)
        amt = o.amt;
    }
    if (mt(val, n1, n2) > mt(o.val, o.n1, o.n2))
      tie(val, n1, n2) = mt(o.val, o.n1, o.n2);
  }
  int n1, n2;
  int amt;
  int val;
};
static const Disp EmptyDisp;
using VDisp = vector<Disp>;
using DispDisp = pair<Disp, Disp>;
inline ostream& operator<<(ostream& o, const Disp& di) {
  o << "(" << (inrange(di.val, 0, (int)Rd.size() - 1) ? Rd[di.val] : -1)
    << "r, " << di.amt << ", " << di.n1 << ", " << di.n2 << ")";
  return o;
}
struct Solution;
Disp dynamicDispBf(Solution& s, bool updateS);
struct DynamicDispersion {
  friend Disp dynamicDispBf(Solution&, bool);
  void init(Solution& s);
  Disp shiftCost(const Solution& s, int u, int g) const;
  Disp swpCost(const Solution& s, int u1, int u2) const;
  void shift(const Solution& s, int u, int g, int p);
  void swp(const Solution& s, int u1, int u2);
  const Disp& disp() const {
    return considerEqualDisp ? (nearSet.empty() ? EmptyDisp : nearSet.front())
                             : gd[gds.front()];
  }
  bool isComplete() const {
    return gds.size() and inrange(gds.front(), 0, (int)gd.size() - 1) and
           gid.size() and utgd.size() and
           inrange(disp().val, 0, (int)Rd.size() - 1);
  }
  VDisp gd;
  vector<int> gds;
  vector<VDisp> gid;
  vector<vector<VDisp>> utgd;
  vector<vector<i16>> numRd;
  VDisp nearSet;
  void setConsiderEqualDisp(bool c) { considerEqualDisp = c; }
  bool considerEqualDisp = false;
  void sortGds() {
    assert((int)size(gd) == m);
    sort(begin(gds), end(gds), [&](int i, int j) { return gd[i] < gd[j]; });
  }
  void rcmpNearSet();
  void gidGain(int g, int u, const VI& ga, int ignore);
  void gidLose(int g, int u, const VI& ga);
  void gidGainLose(int g, int uGain, int uLose, const VI& ga);
  void utgdGain(int g, int u, int w);
  void utgdLose(int g, int u, int w, const VI& ga);
  void utgdGainLose(int g, int uGain, int uLose, int w, const VI& ga);
};
inline void udConsider(VDisp& p, Disp t, bool considerEqualDisp) {
  if (p.size() < 2) {
    p.push_back(t);
    if (p.size() == 2) {
      if (p[0] > p[1]) swap(p[0], p[1]);
      if (considerEqualDisp and p[0].val == p[1].val) p[0].amt = p[1].amt = 2;
    }
  } else {
    if (considerEqualDisp) {
      if (t.val == p[0].val) t.amt = ++p[0].amt;
      if (t.val == p[1].val) t.amt = ++p[1].amt;
    }
    assert(p[0] <= p[1]);
    if (t < p[0])
      p.insert(begin(p), t);
    else if (t < p[1])
      p.insert(begin(p) + 1, t);
    else if (t.val == p[1].val)
      p.push_back(t);
    while (p.size() > 2 and p.back().val != p[1].val)
      p.pop_back();
  }
}
