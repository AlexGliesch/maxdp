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
#include "readinstance.h"
#include "bal.h"
#include "cmdline.h"
#include "dynamicdispersion.h"
#include "ec.h"
#include "main.h"
#include "ub.h"
template <typename T>
void readInstanceData(ifstream& f, T* data, const char* desc) {
  if (not(f >> *data))
    throw logic_error(
        format("failed while parsing \"{}\". The instance could be in "
               "an invalid format.\n",
               desc));
}
void readInstanceAlex(ifstream& f) {
  readInstanceData(f, &n, "number of vertices");
  readInstanceData(f, &m, "number of groups");
  readInstanceData(f, &instPrm.type, "instance type");
  readInstanceData(f, &instPrm.rndSeed, "seed");
  readInstanceData(f, &instPrm.beta, "beta");
  boost::to_lower(instPrm.type);
  if (instPrm.type != "study" and instPrm.type != "weee")
    throw logic_error(format("invalid instance type {}.", instPrm.type));
  tw.resize(m);
  obw.resize(n);
  d.assign(n, vector<double>(n));
  for (auto& w : tw)
    readInstanceData(f, &w, "target weights");
  for (auto& w : obw)
    readInstanceData(f, &w, "object weights");
  if (instPrm.type == "weee") {
    obx.resize(n), oby.resize(n);
    for (int i = 0; i < n; ++i) {
      readInstanceData(f, &obx[i], "object coordinates");
      readInstanceData(f, &oby[i], "object coordinates");
    }
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        d[i][j] = sqrt(pow(obx[i] - obx[j], 2) + pow(oby[i] - oby[j], 2));
  } else if (instPrm.type == "study") {
    oblik.resize(n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < 25; ++j)
        readInstanceData(f, &oblik[i][j], "object Likert points");
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        d[i][j] = 0;
        if (i != j)
          for (int k = 0; k < 25; ++k)
            d[i][j] += abs(oblik[i][k] - oblik[j][k]);
      }
  }
}
void readInstanceMarcus(ifstream& f) {
  readInstanceData(f, &n, "n");
  readInstanceData(f, &m, "m");
  readInstanceData(f, &instPrm.beta, "beta");
  alpha = 0.05;
  instPrm.type = "weee";
  instPrm.rndSeed = 0;
  tw.resize(m);
  obw.resize(n);
  d.assign(n, vector<double>(n));
  obx.resize(n), oby.resize(n);
  for (int i = 0; i < n; ++i) {
    readInstanceData(f, &obw[i], "object weights");
    readInstanceData(f, &obx[i], "object coordinates");
    readInstanceData(f, &oby[i], "object coordinates");
  }
  for (auto& w : tw)
    readInstanceData(f, &w, "target weights");
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      if (i == j)
        d[i][j] = 0.0;
      else
        d[i][j] = sqrt(pow(obx[i] - obx[j], 2) + pow(oby[i] - oby[j], 2));
    }
}
void readInstance(const string& inputFilename) {
  ifstream inputFile(inputFilename);
  if (not inputFile) {
    fmt::print("ERROR: could not open input file.\n");
    exit(EXIT_FAILURE);
  }
  instName = boost::filesystem::path(inputFilename).filename().string();
  instFmt == "alex" ? readInstanceAlex(inputFile)
                    : readInstanceMarcus(inputFile);
  pr("Alpha: {}\n", alpha);
  assert(empty(R));
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      R.emplace_back(i, j);
  sort(begin(R), end(R), [&](II& a, II& b) {
    return d[a.first][a.second] < d[b.first][b.second];
  });
  assert(empty(Rd));
  Rd.push_back(0.0);
  transform(begin(R), end(R), back_inserter(Rd),
            [&](const II& p) { return d[p.first][p.second]; });
  Rd.erase(unique(begin(Rd), end(Rd)), end(Rd));
  nDupl = size(R) + 1 - size(Rd);
  if (nDupl > 0) {
    pr("There are {} duplicate values.\n", nDupl);
    pr("Size(Rd): {}, size(R): {}\n", size(Rd), size(R));
    assert(size(R) + 1 >= size(Rd));
  }
  assert(isUnique(Rd));
  assert(Rd.size() >= 1);
  assert(Rd.front() <= Rd.back());
  assert((ff)Rd.front() == (ff)0.0);
  di.assign(n, VI(n));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      di[i][j] = dispIndex(d[i][j]);
  duu.assign(n, VI(n));
  for (int i = 0; i < n; ++i) {
    iota(begin(duu[i]), end(duu[i]), 0);
    sort(begin(duu[i]), end(duu[i]),
         [&](int a, int b) { return mp(di[i][a], a) < mp(di[i][b], b); });
    for (int j = 0; j < n - 1; ++j)
      assert(di[i][duu[i][j]] <= di[i][duu[i][j + 1]]);
    assert(duu[i][0] == i);
  }
  pr("\n");
  if (ubAlg == "auto") ubAlg = instPrm.type == "weee" ? "ubs" : "none";
  if (balAlg == "vns") {
    if (balMaxShakes < 0) balMaxShakes = 50;
    if (balShakeAmtDbl < 0) balShakeAmtDbl = 0.03;
  } else {
    if (balMaxShakes < 0) balMaxShakes = 35;
    if (balShakeAmtDbl < 0) balShakeAmtDbl = 0.01;
  }
  balShakeAmt = int(balShakeAmtDbl * n);
}
