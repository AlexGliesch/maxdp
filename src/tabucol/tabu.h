#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma once

#include "util/Graph.h"
#include "util/manipulateArrays.h"
#include "util/util.h"

namespace gCol {
namespace pctc {
long tabu(const Graph& g, int* c, int k, unsigned long long maxChecks,
          unsigned long long maxIterations, int tenure, int verbose,
          int frequency, int increment, int** neighbors);
}
} // namespace gCol
