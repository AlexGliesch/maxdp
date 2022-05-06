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
#include "dynamicbalance.h"
#include "solution.h"
void DynamicBalance::init(Solution& s) {
  gw.assign(m, 0);
  gb.assign(m, 1.0 - alpha);
  balanceBF(s, true);
}
double DynamicBalance::shiftCost(const Solution& s, int u, int g) const {
  int p = s.a[u];
  double b = vio;
  if (g != -1) b += deltaGain(g, u) - gb[g];
  if (p != -1) b += deltaLose(p, u) - gb[p];
  ff::fixZero(b);
  if constexpr (dbgBal) {
    auto& r = noConst(s);
    r.a[u] = g;
    assert(ff(b) == ff(balanceBF(r, false)));
    r.a[u] = p;
  }
  return b;
}
double DynamicBalance::swpCost(const Solution& s, int u1, int u2) const {
  int g1 = s.a[u1], g2 = s.a[u2];
  assert(g1 != g2 and g1 != -1 and g2 != -1);
  double b = vio;
  b -= (gb[g1] + gb[g2]);
  b += deltaSwap(g1, u2, u1) + deltaSwap(g2, u1, u2);
  ff::fixZero(b);
  if constexpr (dbgBal) {
    auto& r = noConst(s);
    swap(r.a[u1], r.a[u2]);
    assert(ff(b) == ff(balanceBF(r, false)));
    swap(r.a[u1], r.a[u2]);
  }
  return b;
}
void DynamicBalance::shift(const Solution&, int u, int g, int p) {
  if (g != -1) {
    vio -= gb[g];
    gb[g] = deltaGain(g, u);
    ff::fixZero(gb[g]);
    vio += gb[g];
    gw[g] += obw[u];
    ff::fixZero(gw[g]);
  }
  if (p != -1) {
    vio -= gb[p];
    gb[p] = deltaLose(p, u);
    ff::fixZero(gb[p]);
    vio += gb[p];
    gw[p] -= obw[u];
    ff::fixZero(gw[p]);
  }
  ff::fixZero(vio);
}
void DynamicBalance::swp(const Solution& s, int u1, int u2) {
  int g1 = s.a[u1], g2 = s.a[u2];
  vio -= (gb[g1] + gb[g2]);
  gb[g1] = deltaSwap(g1, u2, u1);
  gb[g2] = deltaSwap(g2, u1, u2);
  vio += gb[g1] + gb[g2];
  gw[g1] += obw[u2] - obw[u1];
  gw[g2] += obw[u1] - obw[u2];
  ff::fixZero(gb[g1]), ff::fixZero(gb[g2]);
  ff::fixZero(gw[g1]), ff::fixZero(gw[g2]);
  ff::fixZero(vio);
}
double balanceBF(Solution& s, bool updateS) {
  vector<double> gw(m, 0.0), gb(m);
  double bal = 0.0;
  for (int u = 0; u < n; ++u)
    if (s.a[u] != -1) gw[s.a[u]] += obw[u];
  for (int g = 0; g < m; ++g) {
    gb[g] = max(0.0, abs(gw[g] - tw[g]) / tw[g] - alpha);
    bal += gb[g];
  }
  if (updateS) {
    swap(gw, s.b.gw);
    swap(gb, s.b.gb);
    s.b.vio = bal;
  }
  ff::fixZero(bal);
  return bal;
}
