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
       
#include "stats.h"
#include "util.h"
inline constexpr bool dbgDisp = false;
inline constexpr bool dbgBal = false;
inline constexpr bool dbgSol = false;
inline string inputFilename;
inline string outputFilename;
inline string instName;
inline string instFmt;
inline size_t rndSeed;
inline double timeLimit;
inline Timer globalTimer;
inline bool iraceTest = false;
inline bool tabuColTest = false;
inline double tabuColDValue;
inline string testType;
inline int n, m;
inline double alpha;
inline VD tw;
inline VD obw;
inline VD obx, oby;
inline vector<array<int, 25>>
    oblik;
inline vector<VD> d;
inline VVI di;
inline VVI duu;
struct InstanceParameters {
  string type;
  double beta;
  size_t rndSeed;
};
inline InstanceParameters instPrm;
inline VII R;
inline VD Rd;
inline int nDupl = 0;
inline int dispIndex(double disp) {
  if ((ff)disp == (ff)0.0) return 0;
  int i = lower_bound(Rd.begin(), Rd.end(), disp) - Rd.begin();
  for (int j : {i - 1, i, i + 1}) {
    if (j >= 0 and j < (int)Rd.size() and Rd[j] == disp) return j;
  }
  return -1;
}
template <typename... A> void pr([[maybe_unused]] A... a) {
#ifndef NDEBUG
  fmt::print(a...);
#endif
}
inline bool ran = false;
