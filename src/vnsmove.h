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
       
#include "util.h"
struct VNSMove {
  enum Type { shift, swap } type;
  ReservoirSampling rs;
  union {
    double valBal;
    Disp valDisp;
  };
  union {
    int u = -1;
    int u1;
  };
  union {
    int g = -1;
    int u2;
  };
  VNSMove(Type tp, double v) : type(tp), valBal(v) {}
  VNSMove(Type tp, int v) : type(tp), valDisp() { valDisp.val = v; }
  void considerBal(double v, int u1OrU, int u2OrG) {
    rs.considerV(ff(v), ff(valBal),
                 [&]() { tie(u1, u2, valBal) = tie(u1OrU, u2OrG, v); });
  }
  void considerDisp(Disp v, int u1OrU, int u2OrG) {
    rs.considerV(valDisp, v,
                 [&]() { tie(u1, u2, valDisp) = tie(u1OrU, u2OrG, v); });
  }
  bool valid() const { return u1 != -1 and u2 != -1; }
};
