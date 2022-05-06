#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma once

#include <vector>

#include "util/Graph.h"
#include "util/util.h"

namespace gCol {
  namespace pctc {
    int generateInitialK(const Graph& g, int alg, int *bestColouring);
    void initializeColoring(const Graph& g, int *c, int k);
    void initializeColoringForTabu(const Graph& g, int *c, int k);
  }
}
