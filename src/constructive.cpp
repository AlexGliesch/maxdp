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
#include "constructive.h"
#include "ub.h"
Solution constructive(Timer t, [[maybe_unused]] bool verbose) {
  if (consAlg == "random")
    return Solution::random();
  else if (consAlg == "trivial")
    return Solution::trivial();
  else {
    assert(consAlg == "greedy");
    if (ubiSubset.empty())
      computeUbi(t, false);
    assert((int)size(ubiSubset) == m);
    Solution s;
    s.dd.setConsiderEqualDisp(nDupl >= n);
    s.init();
    for (int i = 0; i < m; ++i) {
      s.shift(ubiSubset[i], i, true, true);
    }
    while (not s.isComplete()) {
      if (t.timedOut())
        return Solution::trivial();
      int g = 0;
      for (int i = 1; i < m; ++i)
        if (s.b.gw[i] / tw[i] < s.b.gw[g] / tw[g])
          g = i;
      int bu = -1, buD;
      double buB;
      for (int u = 0; u < n; ++u)
        if (s.a[u] == -1) {
          int uD = s.shiftCostDisp(u, g).val;
          double uB = s.shiftCostBal(u, g);
          if (bu == -1 or mp(buD, -buB) < mp(uD, -uB))
            bu = u, buD = uD, buB = ub;
        }
      assert(inrange(bu, 0, n - 1));
      s.shift(bu, g, true, true);
    }
    assert(s.isComplete());
    dynamicDispBf(s, true);
    balanceBF(s, true);
    return s;
  }
}
