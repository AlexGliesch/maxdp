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
namespace stats {
inline double ubiTime = 0.0;
inline int ubiIterToBest = -1;
inline double ubiTtb = 0.0;
inline double ubiAvgCons = -1;
inline int ubiMinCons = NLI::max();
inline double ubiAvg = -1;
inline double ubrbTime = 0.0;
inline int ubrbIterToBest = -1;
inline double ubrbTtb = 0.0;
inline double ubrbAvgCons = -1;
inline int ubrbMinCons = NLI::max();
inline double ubrbAvg = -1;
inline double ubsTime = 0.0;
inline int ubsIter = 0;
inline double ubsTtb = 0.0;
inline int ubsBsIter = 0;
inline int ubsIterToBest = -1;
inline int ubsNumColorCalls = 0;
inline double ubsTotColorTime = 0.0;
inline double ubsAvgColorTime = 0.0;
inline int ubsNumImproves = 0;
inline double ubiOnBestSubset = -1.0;
inline double ubRbOnBestSubset = -1.0;
}
namespace stats {
inline u64 ecNodeExp = 0;
inline u64 ecNumAlt = 0;
inline u64 ecNumLsSteps = 0;
inline u64 ecNumShifts = 0;
inline u64 ecNumSwaps = 0;
inline double ecTimeLstEc = 0;
inline double umdpTime = 0;
inline u64 umdpReplDone = 0;
inline double ecInsertTime = 0;
inline double ecLsTime = 0;
inline u64 ecLsAltImp = 0;
inline double ecTtb = 0;
inline int ecIterTb = 0;
inline double ecTime = 0;
inline u64 ferIter = 0;
inline double ferTotalColTime = 0;
inline double ferMaxColTime = 0;
}
namespace stats {
inline double balTime = 0.0;
inline u64 callsToBalancing = 0;
inline u64 geNodesExp = 0;
inline u64 geIter = 0;
inline u64 geShiftsDone = 0;
inline u64 geNumMaxNodeIncr = 0;
inline double geBBTime = 0.0;
inline u64 geTotDepth = 0;
inline u64 geDidReachTh2 = 0;
inline u64 balVnsShifts = 0;
inline u64 balVnsSwaps = 0;
inline u64 balTotalShakes = 0;
inline double balTtb = 0;
inline bool feasible = false;
}
namespace stats {
inline double oscMaxBalTime = 0.0;
inline double oscEcTime = 0.0;
inline double oscMaxEcTime = 0.0;
inline double oscTtb = 0.0;
inline int oscIterTb = 0;
inline int oscBsSteps = 0;
inline int oscIter = 0;
inline int oscDispFirstFail = -1;
inline double oscBalFirstFail = -1;
inline int oscDispFirstEC = -1;
inline double oscDispFirstECReal = -1.0;
inline double oscBalFirstEC = -1;
inline int oscTotDispStep = 0;
inline int oscMaxDispStep = 0;
inline double oscTotBalStep = 0;
inline double oscMaxBalStep = 0;
}
namespace stats {
inline int numBsSteps = 0;
void printStats();
}
