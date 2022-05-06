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
namespace ubkim {
VII subsetEdgeList(int k) {
  VII e;
  for (int i = 0; i < m + k; ++i)
    for (int j = i + 1; j < m + k; ++j)
      e.emplace_back(i, j);
  assert((int)e.size() == ((m + k) * (m + k - 1)) / 2);
  return e;
}
}
namespace ubsim {
int numEdges(const VI& ss, int ub) {
  int E = 0;
  for (uint i = 0; i < size(ss); ++i)
    for (uint j = i + 1; j < size(ss); ++j)
      E += isEdge(ss[i], ss[j], ub);
  return E;
}
[[maybe_unused]] void initialSubsetNN(VI& ss, int u, int ) {
  assert(not ss.empty());
  copy(begin(duu[u]), begin(duu[u]) + size(ss), begin(ss));
}
void computeDeg(int ub, VI& deg, VI& order) {
  deg.assign(n, 0);
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      if (isEdge(i, j, ub)) ++deg[i], ++deg[j];
  order.resize(n);
  iota(begin(order), end(order), 0);
  stable_sort(begin(order), end(order), [&](int i, int j) {
    return deg[i] > deg[j];
  });
}
void initialSubsetGreedy(VI& ss, int u, int ub, const VI& deg) {
  uint k = 0;
  ss[k++] = u;
  VI valC(n);
  for (int i = 0; i < n; ++i)
    valC[i] = i == u ? -1 : isEdge(i, u, ub);
  while (k < size(ss)) {
    int me = -1;
    for (int i = 0; i < n; ++i)
      if (me == -1 or mp(valC[i], deg[i]) > mp(valC[me], deg[me])) me = i;
    assert(inrange(me, 0, n - 1));
    ss[k++] = me;
    for (int i = 0; i < n; ++i)
      if (valC[i] != -1) {
        if (i == me)
          valC[i] = -1;
        else
          valC[i] += isEdge(i, me, ub);
      }
  }
}
int optimizeSubset(VI& ss, int ub) {
  VB inSS(n, false);
  for (int i : ss)
    inSS[i] = true;
  int val = numEdges(ss, ub);
  VI valC(n, 0);
  for (int i = 0; i < n; ++i)
    for (int j : ss)
      valC[i] += isEdge(i, j, ub);
  bool improved = true;
  const int sz = size(ss);
  int i = 0, k = 0;
  while (improved) {
    improved = false;
    for (int ict = 0; ict < sz and not improved; ++ict, i = (i + 1) % sz)
      for (int kct = 0; kct < n; ++kct, k = (k + 1) % n)
        if (not inSS[k]) {
          int ol = ss[i], nw = k;
          int val2 = val - valC[ol] + valC[nw] - isEdge(ol, nw, ub);
          if (val2 > val) {
            val = val2;
            inSS[ol] = false;
            ss[i] = nw;
            inSS[nw] = true;
            improved = true;
            for (int j = 0; j < n; ++j)
              valC[j] = valC[j] - isEdge(ol, j, ub) + isEdge(nw, j, ub);
            break;
          }
        }
  }
  assert(isUnique(ss));
  return val;
}
}
void computeUbiUbRbOnBestSubset(const VI& ss) {
  if ((int)ss.size() < m) {
    stats::ubRbOnBestSubset = stats::ubiOnBestSubset = -1.0;
    return;
  }
  string s = exec("mktemp");
  s.erase(remove_if(begin(s), end(s), [&](char c) {
    return c == ' ' or c == '\n' or c == '\r';
  }));
  writeReducedInstance(ss, s);
  string cmdUbi = format("./maxdp-r --in {} --test ub --ub ubi --irace", s);
  stats::ubiOnBestSubset = stod(exec(cmdUbi));
  if (ubAlg == "ubrb") {
    stats::ubRbOnBestSubset = Rd[ubrb];
  } else {
    string cmdUbRb = format("./maxdp-r --in {} --test ub --ub ubrb --irace", s);
    stats::ubRbOnBestSubset = stod(exec(cmdUbRb));
  }
  exec(format("rm {}", s));
}
