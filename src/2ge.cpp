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
#include "2ge.h"
#include "bal.h"
#include "main.h"
#include "mcmf.h"
#include "solution.h"
#include "stats.h"
namespace twoge {
constexpr int geIterPrintStep = 100;
bool improveRR(Solution& s, int dmin, int maxNodesBB, Timer t, bool verbose);
bool improveFlow(Solution& s, int dmin, int maxNodesBB, Timer t, bool verbose);
struct BB {
  BB(Solution& s, int g1, int g2, int dmin, int maxNodes);
  int solve(Timer t, bool verbose);
private:
  struct Node {
    Node() = default;
    Node(const BB& b) : f(b.n1 + b.n2, false), g(b.n1 + b.n2, true) {
      numFixed = v[0] = v[1] = rem[0] = rem[1] = 0;
      fill(begin(g), begin(g) + b.n1, false);
      computeValue(b);
      computeLb(b);
    }
    bool operator<(const Node& o) const {
      if (nodeSelectionStrategy == dfs) {
        return id < o.id;
      } else if (nodeSelectionStrategy == valueLb) {
        return ffuple(value, lb) > ffuple(o.value, o.lb);
      } else {
        assert(nodeSelectionStrategy == lbValue);
        return ffuple(lb, value) > ffuple(o.lb, o.value);
      }
    }
    void fix(int u, int val, const BB& b) {
      ++depth;
      fixRec(u, val, b);
      value = gBal(v[0] + rem[0], tw[b.g1]) + gBal(v[1] + rem[1], tw[b.g2]);
      computeLb(b);
    }
    void apply(BB& b) const;
    bool isComplete() const { return numFixed == (int)f.size(); }
    double value = -1.0, lb = -1.0;
    u64 id = 0;
    int numFixed = 0;
    int depth = 0;
    VB f;
    VB g;
    double v[2];
    double rem[2];
  private:
    void computeValue(const BB& b);
    void computeLb(const BB& b);
    void fixRec(int u, int val, const BB& b);
  };
  int selectNextBranchingVar(const BB::Node& n) const;
  Solution& s;
  int maxNodes;
  int g1, g2;
  int n1, n2;
  VI o;
  VVI conf;
  double BStar;
  boost::heap::fibonacci_heap<Node> open;
};
void BB::Node::computeValue(const BB& b) {
  numFixed = 0;
  v[0] = v[1] = rem[0] = rem[1] = 0.0;
  for (int u = 0; u < b.n1 + b.n2; ++u)
    if (f[u])
      v[g[u]] += obw[b.o[u]], ++numFixed;
    else
      rem[g[u]] += obw[b.o[u]];
  value = gBal(v[0] + rem[0], tw[b.g1]) + gBal(v[1] + rem[1], tw[b.g2]);
}
void BB::Node::computeLb(const BB& b) {
  double totRem = rem[0] + rem[1];
  auto SCost = [&](double s) {
    return inrange(ff(s), ff(0.0), ff(totRem))
               ? gBal(v[0] + s, tw[b.g1]) + gBal(v[1] + totRem - s, tw[b.g2])
               : NLD::max();
  };
  double T1 = tw[b.g1], T2 = tw[b.g2], c = v[0] - T1, d = v[1] + totRem - T2;
  lb = vmin(SCost(0.0), SCost(totRem), SCost(-c - T1 * alpha),
            SCost(-c + T1 * alpha), SCost(d - T2 * alpha),
            SCost(d + T2 * alpha));
}
void BB::Node::fixRec(int u, int val, const BB& b) {
  assert(inrange(u, 0, (int)max(f.size(), g.size()) - 1));
  assert(f[u] == false);
  rem[g[u]] -= obw[b.o[u]], v[val] += obw[b.o[u]];
  f[u] = true, g[u] = val;
  ++numFixed;
  for (int i : b.conf[u]) {
    assert((int)f.size() > i and (int) g.size() > i);
    if (f[i])
      assert(g[i] != g[u]);
    else
      fixRec(i, 1 - val, b);
  }
}
void BB::Node::apply(BB& b) const {
  int shiftsDone = 0;
  for (int i = 0; i < (int)f.size(); ++i) {
    int ogi = i >= b.n1;
    if (f[i] and g[i] != ogi) {
      ++shiftsDone;
      b.s.shift(b.o[i], ogi ? b.g1 : b.g2, true, true);
    }
  }
  stats::geShiftsDone += shiftsDone;
}
BB::BB(Solution& s, int g1, int g2, int dmin, int maxNodes)
    : s(s), maxNodes(maxNodes), g1(g1), g2(g2), n1(s.groupSize(g1)),
      n2(s.groupSize(g2)), conf(n1 + n2), BStar(s.b.gb[g1] + s.b.gb[g2]) {
  assert(g1 != g2);
  assert(inrange(g1, 0, m - 1) and inrange(g2, 0, m - 1));
  copy(begin(s.ga[g1]), end(s.ga[g1]), back_inserter(o));
  copy(begin(s.ga[g2]), end(s.ga[g2]), back_inserter(o));
  for (int i = 0; i < n1; ++i)
    for (int j = n1; j < n1 + n2; ++j)
      if (di[o[i]][o[j]] < dmin) conf[i].push_back(j), conf[j].push_back(i);
}
int BB::solve(Timer t, bool verbose) {
  TIME_BLOCK("twoge::BB::solve");
  double tmBef = t.elapsedSecs();
  u64 nodesExpanded = 0;
  u64 nodeId = 0;
  int maxDepth = 0;
  Node bestSoFar(*this);
  bestSoFar.id = nodeId++;
  open.emplace(bestSoFar);
  bestSoFar.value = bestSoFar.lb = NLD::max();
  auto hopeless = [&](const Node& c) {
    return ff(c.lb) >= ff(bestSoFar.value);
  };
  while (not open.empty() and not t.timedOut() and
         nodesExpanded < (u64)maxNodes) {
    Node n = open.top();
    open.pop();
    ++nodesExpanded;
    maxDepth = max(maxDepth, n.depth);
    if (verbose and nodesExpanded % 1000 == 0)
      pr("{}: open: {}, fixed: {}/{}, val: {}, heu: {}, bSF: {}, B*: {}\n",
         nodesExpanded, open.size(), n.numFixed, n1 + n2, n.value, n.lb,
         bestSoFar.value, BStar);
    if (hopeless(n)) continue;
    int j = selectNextBranchingVar(n);
    assert(inrange(j, 0, n1 + n2 - 1) and n.f[j] == false);
    auto n0 = n;
    n0.fix(j, 0, *this);
    n.fix(j, 1, *this);
    if (nodeSelectionStrategy == dfs) {
      if (n0.lb > n.lb) {
        n0.id = nodeId++;
        n.id = nodeId++;
      } else {
        n.id = nodeId++;
        n0.id = nodeId++;
      }
    }
    auto considerCandidate = [&](const Node& c) {
      if (hopeless(c)) return;
      if (ff(c.value) < ff(bestSoFar.value)) bestSoFar = c;
      if (not c.isComplete()) open.emplace(move(c));
    };
    considerCandidate(n0);
    considerCandidate(n);
  }
  if (ff(bestSoFar.value) < ff(BStar)) bestSoFar.apply(*this);
  stats::geBBTime += t.elapsedSecs() - tmBef;
  stats::geTotDepth += maxDepth;
  return nodesExpanded;
}
int BB::selectNextBranchingVar(const BB::Node& n) const {
  int j = -1;
  for (int i = 0; i < n1 + n2; ++i)
    if (n.f[i] == false and (j == -1 or obw[o[i]] > obw[o[j]])) j = i;
  return j;
}
II selectNextPairBipartiteFlow(const Solution& s, const TabuList& tabu,
                               const VVB& canImprove) {
  const int source = 2 * m, sink = source + 1;
  const double LARGE = pow(10.0, 8.0);
  constexpr bool tbEdges = true;
  constexpr bool notCiEdges = false;
  constexpr bool bothBalEdges = false;
  const double sinkCap = 1.0;
  const double sameCap =
      1.0 - alpha;
  const double sameCost = 1;
  const double difCost = 2;
  using Flow = MinCostMaxFlow<double, double>;
  Flow f(m * 2 + 2);
  [[maybe_unused]] int numTb = 0;
  for (int g1 = 0; g1 < m; ++g1) {
    for (int g2 = 0; g2 < m; ++g2) {
      bool tb = tabu.isTabu(g1 * m + g2), ci = canImprove[g1][g2],
           oneImb = s.b.isGrImb(g1) or s.b.isGrImb(g2);
      numTb += g1 != g2 and tb;
      if (g1 == g2)
        f.addEdge(g1, g1 + m, sameCost, sameCap * tw[g1]);
      else if (oneImb and ci and not tb)
        f.addEdge(g1, g2 + m, difCost, LARGE);
      else if ((tbEdges and tb) or (notCiEdges and not ci) or
               (bothBalEdges and not oneImb)) {
        f.addEdge(g1, g2 + m, LARGE, LARGE);
      }
    }
    f.addEdge(source, g1, 0, s.b.gw[g1]);
    f.addEdge(g1 + m, sink, 0, sinkCap * tw[g1]);
  }
  double totalFlow, totalCost;
  tie(totalFlow, totalCost) = f.solve(source, sink);
  Flow::Edge be;
  bool choseTabu = true;
  for (auto& e : f.edge) {
    if (e.from == source or e.to == sink) continue;
    if (inrange(e.from, 0, m - 1) and inrange(e.to, m, 2 * m - 1) and
        e.from != e.to - m and canImprove[e.from][e.to - m] and
        (s.b.isGrImb(e.from) or s.b.isGrImb(e.to - m))) {
      bool tb = tabu.isTabu(e.from * m + e.to - m);
      if (be.from == -1 or (be.flow < e.flow and tb <= choseTabu))
        be = e, choseTabu = tb;
    }
  }
  if (be.to != -1) be.to -= m;
  if (choseTabu and be.from != -1) {
    assert(tabu.isTabu(be.from * m + be.to));
  }
  if (tbEdges == false and be.from == -1) {
    for (int g1 = 0; g1 < m; ++g1)
      for (int g2 = g1 + 1; g2 < m; ++g2)
        if (tabu.isTabu(g1 * m + g2) and
            (s.b.isGrImb(g1) or s.b.isGrImb(g2)) and canImprove[g1][g2]) {
          be.from = g1, be.to = g2;
          break;
        }
  }
#ifndef NDEBUG
  if (be.from == -1)
    for (int g1 = 0; g1 < m; ++g1)
      for (int g2 = g1 + 1; g2 < m; ++g2)
        assert((not s.b.isGrImb(g1) and not s.b.isGrImb(g2)) or
               not canImprove[g1][g2]);
#endif
  return II(be.from, be.to);
}
bool improveFlow(Solution& s, int dmin, int maxNodesBB, Timer t, bool verbose) {
  const uint tabuTenure = tau * m;
  TabuList tabu(m * m, tabuTenure);
  bool improved = false;
  VVB ci(m, VB(m, true));
  int numCi = m * m - m, g1, g2;
  int expandedMaxNodes = 0;
  while (not s.isBalanced() and not t.timedOut()) {
    ++stats::geIter;
    double balBef = s.bal();
    tie(g1, g2) = selectNextPairBipartiteFlow(s, tabu, ci);
    if (g1 == -1 and g2 == -1) {
      if ((u64)maxNodesBB < th2 and expandedMaxNodes) {
        ++stats::geNumMaxNodeIncr;
        return improveFlow(s, dmin, maxNodesBB * 2, t, verbose);
      } else
        break;
    }
    assert(g1 != g2 and inrange(g1, 0, m - 1) and inrange(g2, 0, m - 1));
    assert(ci[g1][g2] and ci[g1][g2] == ci[g2][g1] and
           (s.b.isGrImb(g1) or s.b.isGrImb(g2)));
    int expanded = BB(s, g1, g2, dmin, maxNodesBB).solve(t, false);
    bool expMax = (expanded >= maxNodesBB);
    expandedMaxNodes += expMax;
    stats::geDidReachTh2 += expMax;
    stats::geNodesExp += expanded;
    improved = ff(s.bal()) < ff(balBef);
    tabu.add(g1 * m + g2), tabu.add(g2 * m + g1);
    if (verbose and stats::geIter % geIterPrintStep == 0 ) {
      pr("#{}: g1: {}, g2: {}, B*: {}, B: {}, dr: {}, exp: {}, c.i.: {}, "
         "stats::balfex2gNodesExp: {}, maxNodes: {}, time: {}\n",
         stats::geIter, g1, g2, balBef, s.bal(), s.dispReal(), expanded, numCi,
         stats::geNodesExp, maxNodesBB, t.elapsedSecs());
    }
    if (improved) {
      stats::balTtb = globalTimer.elapsedSecs();
      tabu.advanceIter();
      for (int g = 0; g < m; ++g) {
        assert(ci[g1][g] == ci[g][g1]);
        numCi += (not ci[g1][g]) + (not ci[g2][g]) + (not ci[g][g1]) +
                 (not ci[g][g2]);
        ci[g1][g] = ci[g][g1] = ci[g2][g] = ci[g][g2] = true;
      }
    } else {
      ci[g1][g2] = ci[g2][g1] = false;
      numCi -= 2;
    }
  }
  if (improved and verbose)
    pr("Improved! s: {}, tot. exp.: {}\n", s, stats::geNodesExp);
  return improved;
}
bool improveRR(Solution& s, int dmin, int maxNodesBB, Timer t, bool verbose) {
  bool improved = true;
  while (improved and not s.isBalanced() and not t.timedOut()) {
    double balBef = s.bal();
    improved = false;
    static int g1 = 0, g2 = 1;
    for (int g1c = 0; g1c < m and not improved and not t.timedOut();
         ++g1c, g1 = (g1 + 1) % m)
      for (int g2c = 0; g2c < m and not improved and not t.timedOut();
           ++g2c, g2 = (g2 + 1) % m)
        if (g1 != g2) {
          ++stats::geIter;
          int expanded =
              BB(s, g1, g2, dmin, maxNodesBB).solve(t, false);
          stats::geNodesExp += expanded;
          if (verbose and stats::geIter % geIterPrintStep == 0)
            pr("#{}: g1: {}, g2: {}, B*: {}, B: {}, dr: {}, exp: {}, c.i.: {}, "
               "stats::balfex2gNodesExp: {}, maxNodes: {}, time: {}\n",
               stats::geIter, g1, g2, balBef, s.bal(), s.dispReal(), expanded,
                         0, stats::geNodesExp, maxNodesBB, t.elapsedSecs());
          improved = ff(s.bal()) < ff(balBef);
        }
    if (improved and verbose)
      pr("Improved! s: {}, tot. exp.: {}\n", s, stats::geNodesExp);
  }
  return improved;
}
int deterministicShuffle(Solution& s, int dmin, Timer t,
                         [[maybe_unused]] bool verbose) {
  TIME_BLOCK("deterministicShuffle");
  int shiftsDone = 0;
  if (true) {
    bool doneShift = true;
    static int g1 = 0, g2 = 0;
    for (int iter = 0; iter < balShakeAmt and doneShift and not t.timedOut();
         ++iter) {
      doneShift = false;
      for (int g1c = 0; g1c < m and not doneShift; ++g1c, g1 = (g1 + 1) % m)
        for (int i = 0; i < (int)s.ga[g1].size() and not doneShift; ++i)
          for (int g2c = 0; g2c < m and not doneShift; ++g2c, g2 = (g2 + 1) % m)
            if (g2 != g1 and s.shiftCostDisp(s.ga[g1][i], g2).val >= dmin) {
              doneShift = true;
              ++shiftsDone;
              s.shift(s.ga[g1][i], g2, true, true);
            }
    }
  }
  return shiftsDone;
}
bool improveBalance(Solution& sIn, int dmin, Timer t,
                    bool verbose ) {
  TIME_BLOCK("twoge::improveBalance");
  int numShakes = 0;
  Solution s = sIn;
  for (; numShakes <= balMaxShakes and not t.timedOut(); ++numShakes) {
    [[maybe_unused]] bool ok =
        (strategy == "rr" ? improveRR(s, dmin, th1, t, verbose)
                          : improveFlow(s, dmin, th1, t, verbose));
    if (ff(s.bal()) < ff(sIn.bal())) sIn = s;
    if (sIn.isBalanced()) break;
    if (s.bal() > balShakeThreshold) break;
    if (numShakes != balMaxShakes) {
      if (verbose) pr("Shuffling solution #{}: {}...\n", numShakes, s);
      int shiftsDone = deterministicShuffle(s, dmin, t, verbose);
      if (verbose) pr("Done {} shifts, {}...\n", shiftsDone, s);
    }
  }
  stats::balTotalShakes += numShakes;
  return sIn.isBalanced();
}
}
