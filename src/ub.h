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
inline int ubi = NLI::max();
inline VI ubiSubset;
void computeUbi(Timer t, bool verbose);
inline int ubc = NLI::max();
void computeUbc(Timer t, bool verbose);
inline int ubrb = NLI::max();
inline bool ubrbDoLS = false;
void computeUbrb(Timer t, bool verbose);
inline int ubs = NLI::max();
inline int ubsSubsetSize = -1;
inline bool ubsFewer = false;
inline double ubsSigma = 2.0;
inline bool sigmaFixed = false;
VI computeUbs(Timer t, bool verbose);
inline int ubk = NLI::max();
void computeUbk(int k, Timer t,
                bool verbose);
inline string ubAlg;
inline double ubTimeLimit;
void computeUpperBound(Timer t, bool verbose);
inline int ub;
namespace ubkim {
VII subsetEdgeList(int k);
inline void sortSubsetEdgeList(const VI& S, VII& e) {
  sort(begin(e), end(e), [&](II a, II b) {
    return di[S[a.first]][S[a.second]] > di[S[b.first]][S[b.second]];
  });
}
inline int commonVertex(int e1, int e2, const VII& e) {
  if (e[e1].first == e[e2].first or e[e1].second == e[e2].first)
    return e[e2].first;
  if (e[e1].first == e[e2].second or e[e1].second == e[e2].second)
    return e[e2].second;
  return -1;
}
inline bool edgesAreDisjoint(int e1, int e2, const VII& e) {
  return commonVertex(e1, e2, e) == -1;
}
inline bool edgesAreTriangle(int e1, int e2, int e3, const VII& e) {
  auto c12 = commonVertex(e1, e2, e), c23 = commonVertex(e2, e3, e),
       c13 = commonVertex(e1, e3, e);
  return c12 != -1 and c23 != -1 and c13 != -1 and c12 != c23 and c12 != c13 and
         c13 != c23;
};
}
namespace ubsim {
inline bool isEdge(int i, int j, int ub) { return i != j and di[i][j] < ub; }
int numEdges(const VI& ss, int ub);
void initialSubsetNN(VI& ss, int u, int );
void computeDeg(int ub, VI& deg, VI& order);
void initialSubsetGreedy(VI& ss, int u, int ub, const VI& deg);
int optimizeSubset(VI& ss, int ub);
}
void computeUbiUbRbOnBestSubset(const VI& ss);
