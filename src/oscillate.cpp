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
#include "oscillate.h"
#include "bal.h"
#include "constructive.h"
#include "ec.h"
#include "lowerbound.h"
#include "solution.h"
#include "umdp.h"
void oscillate(Solution& lb, Timer t, bool verbose) {
  if (umdpAlg != "ec") {
    pr("Option --umdpalg was set to {}, but the oscillation "
       "procedure uses EC. Using EC.",
       umdpAlg);
    umdpAlg = "ec";
  }
  computeUpperBound(Timer(ubTimeLimit, t), true);
  auto checkStoppingCriteria = [&](const Solution& s) {
    if (s.dispOptimal() and s.isBalanced()) {
      pr("Optimal solution found; breaking.\n");
      stats::oscTtb = globalTimer.elapsedSecs();
      stats::oscIterTb = stats::oscIter;
      lb = s;
      return true;
    } else if (t.timedOut()) {
      pr("Timed out; breaking.\n");
      return true;
    }
    return false;
  };
  Timer oscTm;
  pr("Computing initial lb...\n");
  lb = lowerBound(t, true);
  pr("lbSol: {}\n", lb);
  if (lbSol.isBalanced()) {
    stats::oscTtb = globalTimer.elapsedSecs();
    stats::oscIterTb = 0;
  } else {
    pr("Lower bound is not balanced, ending.\n");
  }
  Solution s = lb;
  set<Solution> visited;
  visited.insert(s);
  while (s.isBalanced()) {
    if (checkStoppingCriteria(s)) break;
    ++stats::oscIter;
    Timer tm;
    Disp dispBefore = s.disp();
    ec::solve(s, oscMaxEcIter, t, true);
    if (checkStoppingCriteria(s)) break;
    if (visited.count(s)) {
      pr("Duplicate found, breaking...\n");
      break;
    }
    visited.insert(s);
    bool ecImproved = s.disp() >= dispBefore;
    int dispValDiff = s.dispInt() - dispBefore.val;
    stats::oscTotDispStep += dispValDiff;
    stats::oscMaxDispStep = max(stats::oscMaxDispStep, dispValDiff);
    stats::oscTotBalStep += s.bal();
    stats::oscMaxBalStep = max(stats::oscMaxBalStep, s.bal());
    stats::oscEcTime += tm.elapsedSecs();
    stats::oscMaxEcTime = max(stats::oscMaxEcTime, tm.elapsedSecs());
    if (stats::oscIter == 1) {
      stats::oscDispFirstEC = s.dispInt();
      stats::oscDispFirstECReal = s.dispReal();
      stats::oscBalFirstEC = s.bal();
    }
    if (verbose)
      pr("#{} ec.: {}, time: {}\n", stats::oscIter, s, oscTm.elapsedSecs());
    tm.reset();
    balanceSolution(s, s.dispInt(), t, true);
    pr("After 2GE: {}\n", s);
    stats::oscMaxBalTime = max(stats::oscMaxBalTime, tm.elapsedSecs());
    if (verbose)
      pr("#{} bal.: {}, time: {}\n", stats::oscIter, s, oscTm.elapsedSecs());
    if (checkStoppingCriteria(s)) break;
    if (visited.count(s)) {
      pr("Duplicate found, breaking...\n");
      break;
    }
    visited.insert(s);
    if (s.isBalanced()) {
      if (s.dispInt() > lb.dispInt()) {
        if (verbose) pr("New lower bound {}\n", s);
        stats::oscTtb = globalTimer.elapsedSecs();
        stats::oscIterTb = stats::oscIter;
      }
      lb = s;
      if (not ecImproved) {
        pr("EC did not improve; breaking.");
        break;
      }
    } else {
      if (verbose)
        pr("Not balanced, binary searching between dr={} and dr={}.\n",
           lb.dispReal(), s.dispReal());
      stats::oscDispFirstFail = s.dispInt();
      stats::oscBalFirstFail = s.bal();
      int lo = lb.dispInt(), hi = s.dispInt();
      while (lo <= hi and not t.timedOut()) {
        ++stats::oscBsSteps;
        int mid = (lo + hi) / 2;
        Solution sp = lb.dispInt() > mid ? lb : s;
        assert(Rd[mid] <= sp.dispReal());
        balanceSolution(sp, mid, t, true);
        pr("mid: {} (r{}), sp: {}, lb: {}\n", mid, Rd[mid], sp, lb);
        if (sp.isBalanced()) {
          if (lb.disp() < sp.disp()) lb = sp;
          lo = sp.dispInt() + 1;
          stats::oscTtb = globalTimer.elapsedSecs();
          stats::oscIterTb = stats::oscIter;
          if (lb.dispOptimal()) break;
        } else {
          hi = mid - 1;
        }
      }
      break;
    }
  }
  if (verbose) pr("Final lb: {}\n", lb);
}
