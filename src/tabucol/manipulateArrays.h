#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma once

#include "util/util.h"
#include "util/Graph.h"

namespace gCol {
  namespace pctc {
    void moveNodeToColor(int bestNode, int bestColor, const Graph& g, int *c, int **nodesByColor, int **conflicts, int *nbcPosition, int **neighbors,
			 int **tabuStatus, long totalIterations, int tabuTenure);
  }
}
