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
struct Solution;
inline double gBal(double w, double tw) {
  return max(0.0, abs(w - tw) / tw - alpha);
}
struct DynamicBalance {
  void init(Solution& s);
  double shiftCost(const Solution& s, int u, int g) const;
  double swpCost(const Solution& s, int u1, int u2) const;
  void shift(const Solution& s, int u, int g, int p);
  void swp(const Solution& s, int u1, int u2);
  double bal() const { return vio; }
  bool isGrImb(int g) const { return ff(gb[g]) > ff(0.0); }
  double deltaGain(int g, int u) const { return gBal(gw[g] + obw[u], tw[g]); }
  double deltaLose(int g, int u) const { return gBal(gw[g] - obw[u], tw[g]); }
  inline double deltaSwap(int g, int ug, int ul) const {
    return gBal(gw[g] + obw[ug] - obw[ul], tw[g]);
  }
  VD gw;
  VD gb;
  double vio;
};
double balanceBF(Solution& s, bool updateS);
