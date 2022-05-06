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
#include "lowerbound.h"
#include "bal.h"
#include "constructive.h"
#include "solution.h"
#if !defined(__CYGWIN__)
#define USE_FEASIBILITY_MODEL 
#include <ilcplex/ilocplex.h>
#endif
#ifdef USE_FEASIBILITY_MODEL
Solution solveFeasibilityModel(Timer t, bool) {
  IloEnv env;
  IloModel model(env);
  IloNumVarArray xij(env, n * m, 0, 1, ILOBOOL);
  model.add(xij);
  IloAdd(model, IloMinimize(env));
  for (int j = 0; j != n; j++) {
    IloExpr e(env);
    for (int i = 0; i != m; i++) {
      e += xij[i * n + j];
    }
    model.add(e == 1);
  }
  for (int i = 0; i != m; i++) {
    IloExpr e(env);
    for (int j = 0; j != n; j++)
      e += obw[j] * xij[i * n + j];
    model.add((1 - alpha) * tw[i] <= e <= (1 + alpha) * tw[i]);
  }
  IloCplex solver(model);
  solver.setParam(IloCplex::Param::Threads, 1);
  solver.setParam(IloCplex::TiLim, t.secsLeft());
  solver.setParam(IloCplex::Param::WorkMem, 2048);
  solver.setOut(env.getNullStream());
  solver.setWarning(env.getNullStream());
  solver.solve();
  pr("Solver status: {}\n", solver.getStatus());
  Solution s;
  if (not solver.isPrimalFeasible()) {
    s.init();
    return s;
  } else {
    IloNumArray val(env);
    solver.getValues(val, xij);
    VI a(n, 0);
    for (int j = 0; j != n; j++) {
      for (int i = 0; i != m; ++i)
        if (val[i * n + j]) a[j] = i;
      assert(inrange(a[j], 0, m - 1));
    }
    s.populate(a);
    return s;
  }
}
#endif
Solution lowerBound(Timer t, bool verbose) {
  Solution s;
  s = constructive(t, verbose);
  balanceSolution(s, 0, t, verbose);
  assert(s.isComplete());
#ifdef USE_FEASIBILITY_MODEL
  if (not s.isBalanced()) {
    pr("No balanced solution found, going to run feasibility model.\n");
    s = solveFeasibilityModel(t, verbose);
    pr("Solution value: {}\n", s.isComplete() ? s.dispReal() : -1.0);
  }
#endif
  stats::feasible = s.isBalanced();
  return s;
}
