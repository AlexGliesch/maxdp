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
       
#include "dynamicbalance.h"
#include "dynamicdispersion.h"
#include "main.h"
#include "ub.h"
#include "util.h"
struct Solution {
  bool operator<(const Solution &s) const {
    return tie(disp(), b.vio, a) < tie(s.disp(), s.b.vio, s.a);
  }
  void init();
  void populate(const VI &assigned);
  static Solution random();
  static Solution trivial();
  bool isComplete() const { return numAssigned == n; }
  bool isDispComplete() const { return dd.isComplete(); }
  Disp shiftCostDisp(int u, int g) const { return dd.shiftCost(*this, u, g); }
  Disp swpCostDisp(int u1, int u2) const { return dd.swpCost(*this, u1, u2); }
  double shiftCostBal(int u, int g) const { return b.shiftCost(*this, u, g); }
  double swpCostBal(int u1, int u2) const { return b.swpCost(*this, u1, u2); }
  void shift(int u, int g, bool updBal, bool updDsp);
  void swp(int u1, int u2, bool updBal, bool updDsp);
  int randomWalk(int numMoves, int minDisp);
  void checkCorrect(Timer t) const;
  void setConsiderEqualDisp(bool c) {
    dd.setConsiderEqualDisp(c);
    dynamicDispBf(*this, true);
  }
  void writeToFile(const string &filename);
  static Solution readFromFile(const string &filename);
  int groupSize(int g) const { return size(ga[g]); }
  const Disp &disp() const { return dd.disp(); }
  int dispInt() const { return disp().val; }
  double dispReal() const { return disp().real(); }
  int dispN1() const { return dd.disp().n1; }
  int dispN2() const { return dd.disp().n2; }
  double bal() const { return b.bal(); }
  bool isBalanced() const { return isComplete() and ff(bal()) == ff(0.0); }
  bool dispOptimal() const { return dispInt() == ub; }
  int numAssigned;
  VI a;
  VVI ga;
  VI iiga;
  DynamicDispersion dd;
  DynamicBalance b;
};
inline ostream &operator<<(ostream &o, const Solution &s) {
  o << "(b " << s.bal() << ", dr "
    << (inrange(s.dispInt(), 0, (int)size(Rd) - 1) ? s.dispReal() : -1);
  return o << (s.dispOptimal()
                   ? (s.isBalanced() ? ", optimal)" : ", disp. optimal)")
                   : ")");
}
inline Solution
    umdpSol,
    lbSol,
    consSol;
