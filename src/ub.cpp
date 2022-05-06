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
#include "ub.h"
#include "color.h"
#include "main.h"
void computeUpperBound(Timer t, bool verbose) {
  if (ubAlg == "ubi") {
    computeUbi(t, verbose);
    ub = ubi;
  } else if (ubAlg == "ubc") {
    throw logic_error("Deprecated.");
    computeUbc(t, verbose);
    ub = ubc;
  } else if (ubAlg == "ubrb") {
    computeUbrb(t, verbose);
    ub = ubrb;
  } else if (ubAlg == "ubs" or ubAlg == "all" or
             ubAlg == "ubsblind") {
    [[maybe_unused]] VI ss = computeUbs(t, verbose);
    assert(ubAlg != "all" or (ub <= ubrb and ub <= ubi));
    ub = ubs;
  } else if (ubAlg == "ubk") {
    print("summary_line ");
    print("instance={} ", instName);
    print("instanceType={} ", instPrm.type);
    print("n={} ", n);
    print("m={} ", m);
    print("beta={} ", instPrm.beta);
    print("ubAlg={} ", ubAlg);
    VI ks;
    for (int i = 1; i <= 9; ++i)
      ks.push_back(i);
    for (double i = 10; i <= 50; i += 10)
      ks.push_back(i);
    for (int k : ks) {
      Timer ubkTimer;
      ubk = (int)Rd.size() - 1;
      computeUbk(k, t, false);
      print("k{}ub={} ", k, ubk);
    }
    print("\n");
    exit(EXIT_SUCCESS);
  } else {
    assert(ubAlg == "none");
    ub = size(Rd) - 1;
  }
  assert(inrange(ub, 0, (int)size(Rd) - 1));
}
