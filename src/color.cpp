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
#include "btdsatur/bktdsat.h"
#include "btdsatur/colorrtns.h"
#include "btdsatur/graph.h"
#include "btdsatur/mysys.h"
#include "coudert.h"
#include "ec.h"
#include "main.h"
#include "ub.h"
void writeDimacsSAT(int ub, const VI& ss, const string& filename) {
  string s;
  int numClauses = size(ss);
  for (uint i = 1; i <= size(ss); ++i) {
    for (int j = 1; j <= m; ++j)
      s += format("{} ", (i - 1) * m + j);
    s += "0\n";
  }
  for (int g1 = 1; g1 <= m; ++g1)
    for (int g2 = g1 + 1; g2 <= m; ++g2) {
      numClauses += size(ss);
      for (uint i = 1; i <= size(ss); ++i)
        s += format("-{} -{} 0\n", (i - 1) * m + g1, (i - 1) * m + g2);
    }
  for (uint i = 1; i <= size(ss); ++i)
    for (uint j = i + 1; j <= size(ss); ++j)
      if (di[ss[i - 1]][ss[j - 1]] < ub) {
        numClauses += m;
        for (int g = 1; g <= m; ++g)
          s += format("-{} -{} 0\n", (i - 1) * m + g, (j - 1) * m + g);
      }
  ofstream f(filename);
  f << "p cnf " << size(ss) * m << " " << numClauses << endl;
  f << s;
}
void writeDimacsGraph(int ub, const VI& ss, const string& filename) {
  string s;
  int numEdges = 0;
  for (uint i = 0; i < size(ss); ++i)
    for (uint j = i + 1; j < size(ss); ++j)
      if (di[ss[i]][ss[j]] < ub)
        ++numEdges, s += format("e {} {}\n", i + 1, j + 1);
  ofstream f(filename);
  assert(f.good());
  f << format(
      "c ub(={})-induced subgraph of MaxDP instance {}. Is it {}-colorable?\n",
      Rd[ub], instName, m);
  f << format("p edge {} {}\n", size(ss), numEdges);
  f << s;
}
void writeReducedInstance(const VI& ss, const string& filename) {
  ofstream f(filename);
  string s;
  s += format("{} {}\n{} {} {}\n", ss.size(), m, instPrm.type, instPrm.rndSeed,
              instPrm.beta);
  s += format("{}\n", tw);
  for (int i : ss)
    s += format("{} ", obw[i]);
  s += "\n";
  for (int i : ss) {
    if (instPrm.type == "weee")
      s += format("{} {}\n", obx[i], oby[i]);
    else
      s += format("{}\n", oblik[i]);
  }
  f << s;
}
bool isLbColorSATGlucose(int ub, const VI& ss, Timer t) {
  throw logic_error("isLbColorSATGlucose() is deprecated.");
  const string glucosePath = "../coloring/glucose/simp/glucose_static";
  if (not fileExists(glucosePath)) {
    print("\ERROR: glucose does not exist at {}. "
          "Please build the program before running MaxDP.\n\n",
          glucosePath);
    exit(EXIT_SUCCESS);
  }
  string filename = format("sat-{}-{}", instName, uniqueRandomSeed());
  writeDimacsSAT(ub, ss, filename);
  int secsLeft = int(t.secsLeft() + 1);
  assert(secsLeft > 0);
  string output = exec(
      format("timeout {} ./{} -verb=0 {}", secsLeft, glucosePath, filename));
  exec(format("rm {}", filename));
  if (t.timedOut()) throw TimeoutException(t);
  return boost::contains(output, "s SATISFIABLE");
}
bool isLbColorCoudert(int ub, const VI& ss, int& chi, VI& C, Timer t) {
  coudert::Graph g;
  for (uint i = 0; i < size(ss); ++i)
    add_vertex(coudert::VertexInformation(i), g);
  for (uint i = 0; i < size(ss); ++i) {
    assert(inrange(ss[i], 0, n - 1));
    for (uint j = i + 1; j < size(ss); ++j) {
      assert(inrange(ss[j], 0, n - 1));
      if (di[ss[i]][ss[j]] < ub) add_edge(i, j, g);
    }
  }
  brandes_betweenness_centrality(
      g, get(&coudert::VertexInformation::centrality, g));
  C = coudert::maxClique(g);
  for (auto c : C)
    g[c].clique = true;
  if ((int)size(C) > m) {
    chi = size(C);
    return false;
  }
  chi = coudert::seqColor(g, C, m, t).first;
  return chi != -1 and chi <= m;
}
bool isLbColorHEA(int ub, const VI& ss, Timer t) {
  throw logic_error("isLbColorHEA() is deprecated.");
  const string HEAPath = "../coloring/hybrid-ea/HybridEA";
  if (not fileExists(HEAPath)) {
    print("\ERROR: HEA program does not exist at {}. "
          "Please build the program before running MaxDP.\n\n",
          HEAPath);
    exit(EXIT_SUCCESS);
  }
  string filename = format("hea-{}-{}", instName, uniqueRandomSeed());
  writeDimacsGraph(ub, ss, filename);
  if (t.secsLeft() <= 0) return true;
  string output = exec(format("./{} {} -s 100000000000000 -T {} -t {}", HEAPath,
                              filename, m, t.secsLeft()));
  exec(format("rm {}", filename));
  return not output.empty() and stoi(output) <= m;
}
bool isLbColorBTDSatur(int ub, const VI& ss, Timer t) {
  throw logic_error("isLbColorBTDSatur() is deprecated.");
  const string DSaturPath =
      "../coloring/backtracking-d-satur/BacktrackingDSatur";
  if (not fileExists(DSaturPath)) {
    print("\nError: backtracking-DSatur program does not exist at {}. "
          "Please build the program before running MaxDP.\n\n",
          DSaturPath);
    exit(EXIT_SUCCESS);
  }
  string filename = format("dsatur-{}-{}", instName, uniqueRandomSeed());
  writeDimacsGraph(ub, ss, filename);
  int secsLeft = int(t.secsLeft() + 1);
  assert(secsLeft > 0);
  string output = exec(format("timeout {} ./{} "
                              "{} -s 100000000000000 -T {}",
                              secsLeft, DSaturPath, filename, m));
  exec(format("rm {}", filename));
  if (t.timedOut()) {
    return true;
  }
  return not output.empty() and stoi(output) <= m;
}
namespace btdsatur {
extern int colorsearch(int targetnumcolors, Timer t);
extern unsigned long long numConfChecks;
extern unsigned long long maxChecks;
extern int verbose;
}
bool isLbColorBTDSaturIntegrated(int ub, const VI& ss, Timer t) {
  srand(n + m );
  btdsatur::partitionflag = 0;
  btdsatur::cheatflag = 0;
  btdsatur::order = size(ss);
  btdsatur::maxChecks = 100000000000000LL;
  btdsatur::numConfChecks = 0;
  btdsatur::verbose = 0;
  memset(btdsatur::graph, 0, GRAPHSIZE);
  for (uint i = 0; i < size(ss); ++i)
    for (uint j = i + 1; j < size(ss); ++j)
      if (di[ss[i]][ss[j]] < ub) {
        btdsatur::setedge(i, j);
        btdsatur::setedge(j, i);
      }
  if (t.secsLeft() <= 0) return true;
  try {
    int colors = btdsatur::colorsearch(m, t);
    return colors <= m;
  } catch (TimeoutException& e) {
    return true;
  }
}
bool isLbColor(int ub, const VI& ss, Timer t) {
  bool ans;
  Timer colorTimer;
  if (lbColorAlg == lbColorGlucose) {
    ans = isLbColorSATGlucose(ub, ss, t);
  } else if (lbColorAlg == lbColorCoudert) {
    int chi;
    VI clique;
    ans = isLbColorCoudert(ub, ss, chi, clique, t);
  } else if (lbColorAlg == lbColorHEA) {
    ans = isLbColorHEA(ub, ss, t);
  } else {
    assert(lbColorAlg == lbColorBTDSatur);
    ans = isLbColorBTDSaturIntegrated(ub, ss, t);
  }
  ++stats::ubsNumColorCalls;
  stats::ubsTotColorTime += colorTimer.elapsedSecs();
  stats::ubsAvgColorTime =
      stats::ubsTotColorTime / double(stats::ubsNumColorCalls);
  return ans;
}
