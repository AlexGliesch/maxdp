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
#include "umdp.h"
#include "constructive.h"
#include "ec.h"
#include "main.h"
#include "solution.h"
#include "umdpfernandez.h"
Solution umdp(Timer t, bool verbose) {
  Solution best = constructive(t, verbose);
  consSol = best;
  pr("Initial: {}\n", best);
  best.checkCorrect(t);
  if (consAlg != "random")
    umdpRepl = 1;
  Timer umdpTimer;
  for (int i = 0; i < umdpRepl and not best.dispOptimal() and not t.timedOut();
       ++i) {
    ++stats::umdpReplDone;
    Solution s = umdpRepl == 1 ? consSol : constructive(t, verbose);
    if (s.disp() > consSol.disp())
      consSol = s;
    if (umdpAlg == "vns") {
      throw logic_error("UMDP VNS is deprecated.");
    } else if (umdpAlg == "fer") {
      s = solveUmdpFernandez(best, t, true);
      assert(s.dispInt() >= best.dispInt());
    } else {
      assert(umdpAlg == "ec");
      ec::solve(s, NLI::max(), t, verbose);
    }
    s.checkCorrect(t);
    pr("UMDP replication #{}: r{} ({}){}\n", i + 1, s.dispReal(), s.dispInt(),
       s.dispOptimal() ? ", disp. optimal" : "");
    if (s.disp() > best.disp())
      best = s;
  }
  stats::umdpTime = umdpTimer.elapsedSecs();
  return best;
}
