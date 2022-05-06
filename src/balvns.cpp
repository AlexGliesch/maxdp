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
#include "balvns.h"
#include "2ge.h"
#include "bal.h"
#include "main.h"
#include "solution.h"
#include "vnsmove.h"
namespace {
VB ci;
VNSMove findShift(const Solution& s, int dmin) {
  if (ci.empty()) ci.assign(m, true);
  static int g1 = 0, g2 = 0;
  double sBal = s.bal();
  VNSMove mv(VNSMove::shift, NLD::max());
  for (int gc = 0; gc < m; ++gc, g1 = (g1 + 1) % m) {
    for (int g2c = 0; g2c < m; ++g2c, g2 = (g2 + 1) % m)
      if (g2 != g1 and (s.b.gw[g1] < tw[g1] and s.b.gw[g2] > tw[g2]) and
          (ci[g2] == true or ci[g1] == true) and s.groupSize(g2) > 1)
        for (int u : s.ga[g2]) {
          if (s.shiftCostDisp(u, g1).val < dmin) continue;
          double valBal = s.shiftCostBal(u, g1);
          if (valBal < mv.valBal) {
            mv.valBal = valBal;
            mv.u = u, mv.g = g1;
            if (mv.valBal < sBal) {
              ci[g1] = ci[g2] = true;
              g1 = (g1 + 1) % m;
              g2 = (g2 + 1) % m;
              return mv;
            }
          }
        }
  }
  ci.assign(m, false);
  return mv;
}
VNSMove findSwap(const Solution& s, int dmin) {
  static int u1 = 0, u2 = 0;
  double sBal = s.bal();
  VNSMove mv(VNSMove::swap, NLD::max());
  for (int u1c = 0; u1c < n; ++u1c, u1 = (u1 + 1) % n) {
    for (int u2c = 0; u2c < n; ++u2c, u2 = (u2 + 1) % n)
      if (s.a[u1] != s.a[u2]) {
        if (s.swpCostDisp(u1, u2).val < dmin) continue;
        double valBal = s.swpCostBal(u1, u2);
        if (valBal < mv.valBal) {
          mv.valBal = valBal;
          mv.u1 = u1, mv.u2 = u2;
          if (mv.valBal < sBal) {
            ci[s.a[u1]] = ci[s.a[u2]] = true;
            u1 = (u1 + 1) % n;
            u2 = (u2 + 1) % n;
            return mv;
          }
        }
      }
  }
  return mv;
}
}
void balVNS(Solution& s, int dmin, Timer t, [[maybe_unused]] bool verbose) {
  TIME_BLOCK("balVNS");
  Solution bst = s;
  int numMvs = 0;
  int numShakes = 0;
  while (not t.timedOut() and not bst.isBalanced()) {
    VNSMove mv = findShift(s, dmin);
    if (mv.valBal >= s.bal()) {
      VNSMove mvSw = findSwap(s, dmin);
      if (mvSw.valBal < mv.valBal) mv = mvSw;
    }
    bool mvImproves = ff(mv.valBal) < ff(s.bal());
    if (mv.valid() and mvImproves) {
      [[maybe_unused]] int p = s.a[mv.u];
      if (mv.type == VNSMove::shift) {
        s.shift(mv.u, mv.g, true, true);
        ++stats::balVnsShifts;
      } else if (mv.type == VNSMove::swap) {
        s.swp(mv.u1, mv.u2, true, true);
        ++stats::balVnsSwaps;
      }
      assert(s.dispInt() >= dmin);
      s.checkCorrect(t);
      if (ff(s.bal()) < ff(bst.bal())) {
        bst = s;
        numShakes = 0;
        stats::balTtb = globalTimer.elapsedSecs();
      }
      ++numMvs;
      if (verbose and numMvs % 1000 == 0) {
        pr("#{}: {} ({}/{}--{}); shakes: {}, s: {}, bst: {}, time: {}\n",
           numMvs, mv.type == VNSMove::shift ? "shift" : "swap", mv.u1, p,
           mv.u2, numShakes, s, bst, t.elapsedSecs());
      }
      if (bst.isBalanced()) break;
    } else {
      if (numShakes >= balMaxShakes) break;
      [[maybe_unused]] double balBefore = s.bal();
      s = bst;
      [[maybe_unused]] int shakeSteps = s.randomWalk(balShakeAmt, dmin);
      ++numShakes;
      ++stats::balTotalShakes;
      if (verbose)
        pr("Shaking... moves: {}, shakes: {}, shake steps: {}, s before: {}, "
           "bst: {}\n",
           numMvs, numShakes, shakeSteps, balBefore, bst);
    }
  }
  s = move(bst);
}
