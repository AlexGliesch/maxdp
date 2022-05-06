#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
/**********************************************************************
 * \file PartialColAndTabuCol.h
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 *   \version $Id: emacs 7968 2017-05-17 23:03:17Z ritt $
 *   \date Time-stamp: <2018-06-21 20:55:00 ritt>
 **********************************************************************/
#include "util/Graph.h"
#include "util/util.h"
#include <limits>

namespace gCol {
using namespace std;
namespace pctc {
int pctc(const Graph& g, vector<int>& best, confTimeLogger& log, int algorithm,
         int tenure, unsigned long long maxChecks, int targetCols,
         int randomSeed, int verbose, int constructiveAlg,
         int maxIterations = std::numeric_limits<int>::max());
int main(int argc, char** argv);
} // namespace pctc
} // namespace gCol
