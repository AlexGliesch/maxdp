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
#include "stats.h"
#include "2ge.h"
#include "bal.h"
#include "ec.h"
#include "main.h"
#include "oscillate.h"
#include "solution.h"
#include "ub.h"
#include "umdp.h"
namespace stats {
void printUbStats() {
  if (iraceTest) {
    print("{} ", Rd[ub]);
    return;
  }
  if (ubAlg == "ubi" or ubAlg == "all") {
    print("ubiTime={:.5f} ", ubiTime);
    print("ubiIterToBest={} ", ubiIterToBest);
    print("ubiTtb={} ", ubiTtb);
    print("ubiAvgCons={:.5f} ", ubiAvgCons);
    print("ubiMinCons={} ", ubiMinCons);
    print("ubiAvg={:.5f} ", ubiAvg);
    print("ubi={} ", ubi);
    print("ubiReal={} ", not inrange(ubi, 0, (int)Rd.size() - 1)
                             ? "NA"
                             : format("{:.5f}", Rd[ubi]));
  } else {
    print("ubiTime=NA ubiIterToBest=NA ubiTtb=NA ubiAvgCons=NA ubiMinCons=NA "
          "ubiAvg=NA ubi=NA ubiReal=NA ");
  }
  if (ubAlg == "ubrb" or ubAlg == "all") {
    print("ubRbTime={:.5f} ", ubrbTime);
    print("ubRbIterToBest={} ", ubrbIterToBest);
    print("ubrbTtb={} ", ubrbTtb);
    print("ubRbAvgCons={:.5f} ", ubrbAvgCons);
    print("ubRbMinCons={} ", ubrbMinCons);
    print("ubRbAvg={:.5f} ", ubrbAvg);
    print("ubRb={} ", ubrb);
    print("ubRbReal={} ",
          not inrange(ubrb, 0, (int)Rd.size() - 1)
              ? "NA"
              : format("{:.5f}", Rd[ubrb]));
  } else {
    print("ubRbTime=NA ubRbIterToBest=NA ubrbTtb=NA ubRbAvgCons=NA "
          "ubRbMinCons=NA ubRbAvg=NA ubRb=NA ubRbReal=NA ");
  }
  if (boost::starts_with(ubAlg, "ubs") or ubAlg == "all") {
    print("ubsTime={:.5f} ", ubsTime);
    print("ubsSigma={} ", ubsSigma);
    print("ubsSubsetSize={} ", ubsSubsetSize);
    print("ubsIter={} ", ubsIter);
    print("ubsTtb={} ", ubsTtb);
    print("ubsBsIter={} ", ubsBsIter);
    print("ubsIterToBest={} ", ubsIterToBest);
    print("ubsNumImproves={} ", ubsNumImproves);
    print("ubsNumColorCalls={} ", ubsNumColorCalls);
    print("ubsAvgColorTime={:.5f} ", ubsAvgColorTime);
    print("ubs={} ", ubs);
    print("ubsReal={} ", ubs == NLI::max()
                             ? "NA"
                             : format("{:.5f}", Rd[ubs]));
  } else {
    print(
        "ubsTime=NA ubsSigma=NA ubsSubsetSize=NA ubsIter=NA ubsTtb=NA "
        "ubsBsIter=NA ubsIterToBest=NA ubsNumColorCalls=NA ubsAvgColorTime=NA "
        "ubs=NA ubsReal=NA ");
  }
  print("ubiOnBestSubset={} ",
        valOrNA(ubiOnBestSubset != -1.0, ubiOnBestSubset));
  print("ubRbOnBestSubset={} ",
        valOrNA(ubRbOnBestSubset != -1.0, ubRbOnBestSubset));
}
void printUmdpStats() {
  if (iraceTest) {
    print("{} ", -umdpSol.dispReal());
    return;
  }
  if (not outputFilename.empty()) {
    umdpSol.writeToFile(outputFilename);
  }
  print("umdpAlg={} ", umdpAlg);
  print("ubTime={} ", ubsTime + ubiTime + ubrbTime);
  print("ub={} ", ub);
  print("lb={} ", consSol.dispInt());
  print("p1={} ", ec::p1);
  print("p2={} ", ec::p2);
  print("p3={} ", ec::p3);
  print("disp={} ", umdpSol.disp().val);
  print("dispReal={} ", Rd[umdpSol.disp().val]);
  print("umdpTime={} ", umdpTime);
  print("timeLastEC={} ", ecTimeLstEc);
  print("ecTtb={} ", ecTtb);
  print("ecIterTb={} ", ecIterTb);
  print("ecNumAlt={} ", ecNumAlt);
  print("ecTime={} ", ecTime);
  print("ecInsertTime={} ", ecInsertTime);
  print("ecNodes={} ", ecNodeExp);
  print("lsTime={} ", ecLsTime);
  print("lsAltImpr={} ", ecLsAltImp);
  print("lsSteps={} ", ecNumLsSteps);
  print("lsShifts={} ", ecNumShifts);
  print("lsSwaps={} ", ecNumSwaps);
  print("ferIter={} ", ferIter);
  print("ferAvgColTime={} ",
        valOrNA(ferIter > 0, ferTotalColTime / double(ferIter)));
  print("ferMaxColTime={} ", ferMaxColTime);
}
void printFullStats() {
  if (iraceTest) {
    print("{} ", -lbSol.dispReal());
    return;
  }
  print("ttb={} ", oscTtb);
  print("iter={} ", oscIter);
  print("iterTb={} ", oscIterTb);
  print("bsSteps={} ", oscBsSteps);
  print("disp={} ", lbSol.dispInt());
  print("dispReal={} ", lbSol.disp().real());
  print("bal={} ", lbSol.bal());
  print("dispFirstFail={} ", valOrNA(oscDispFirstFail != -1, oscDispFirstFail));
  print("balFirstFail={} ", valOrNA(oscBalFirstFail != -1, oscBalFirstFail));
  print("dispFirstEC={} ", oscDispFirstEC);
  print("dispFirstECReal={} ", oscDispFirstECReal);
  print("balFirstEC={} ", oscBalFirstEC);
  print("ubTime={} ", ubsTime + ubiTime + ubrbTime);
  print("ub={} ", ub);
  print("ubReal={} ", Rd[ub]);
  print("ecTime={} ", oscEcTime);
  print("ecTimeMax={} ", oscMaxEcTime);
  print("ecNodes={} ", ecNodeExp);
  print("maxEcIter={} ", oscMaxEcIter);
  print("avgDispStep={} ",
        valOrNA(oscIter > 0, oscTotDispStep / double(oscIter)));
  print("maxDispStep={} ", oscMaxDispStep);
  print("feasible={} ", int(feasible));
  print("balAlg={} ", balAlg);
  print("th1={} ", twoge::th1);
  print("th2={} ", twoge::th2);
  print("balTime={} ", balTime);
  print("balTtb={} ", balTtb);
  print("bbTime={} ", geBBTime);
  print("balTimeMax={} ", oscMaxBalTime);
  print("balNodes={} ", geNodesExp);
  print("balCalls={} ", callsToBalancing);
  print("avgBalStep={} ",
        valOrNA(oscIter > 0, oscTotBalStep / double(oscIter)));
  print("maxBalStep={} ", oscMaxBalStep);
  print("geNodes={} ", geNodesExp);
  print("geNumBB={} ", geIter);
  print("avgBbNodes={} ", valOrNA(geIter > 0, geNodesExp / double(geIter)));
  print("avgBbMoves={} ", valOrNA(geIter > 0, geShiftsDone / double(geIter)));
  print("avgBbDepth={} ", valOrNA(geIter > 0, geTotDepth / double(geIter)));
  print("numMaxNodeIncr={} ", geNumMaxNodeIncr);
  print("geDidReachTh2={} ", geDidReachTh2);
  print("balShakes={} ", balTotalShakes);
  print("balVnsShifts={} ", valOrNA(balAlg == "vns", balVnsShifts));
  print("balVnsSwaps={} ", valOrNA(balAlg == "vns", balVnsSwaps));
}
inline void printTimedBlocks() {
#ifdef USE_TIMED_BLOCKS
  vector<pair<double, string>> tbS;
  for (const auto& p : timedBlocks)
    tbS.emplace_back(-p.second, p.first);
  sort(begin(tbS), end(tbS));
  if (!timedBlocks.empty()) {
    VS s;
    size_t maxSize = 0;
    for (const auto& p : tbS) {
      s.push_back(fmt::format("{}: {}", p.second, -p.first));
      maxSize = std::max(maxSize, s.back().size());
    }
    fmt::print("{}\n| Timed blocks:{}|\n", string(maxSize + 4, '='),
               string(maxSize - 12, ' '));
    for (const auto& i : s)
      fmt::print("| {}{} |\n", i, string(maxSize - i.size(), ' '));
    fmt::print("{}\n", string(maxSize + 4, '='));
  }
#endif
}
void printBaseStats() {
  if (not iraceTest) {
    print("instance={} ", instName);
    print("instanceType={} ", instPrm.type);
    print("n={} ", n);
    print("m={} ", m);
    print("R={} ", Rd.size());
    print("beta={} ", instPrm.beta);
    print("duplDist={} ", nDupl);
    print("seed={} ", rndSeed);
    print("alpha={} ", alpha);
    print("ubAlg={} ", ubAlg);
    print("testType={} ", testType);
    print("tau={} ", twoge::tau);
    print("time={} ", globalTimer.elapsedSecs());
  }
}
void printStats() {
  if (ubAlg == "ubk") {
    return;
  }
  if (not ran) return;
  pr("\n"), printTimedBlocks(), pr("\n");
  if (tabuColTest) return;
  if (not iraceTest) {
    print("summary_line ");
    printBaseStats();
  }
  if (testType == "ub")
    printUbStats();
  else if (testType == "umdp")
    printUmdpStats();
  else if (testType == "full")
    printFullStats();
  if (not iraceTest) print("\n");
}
}
