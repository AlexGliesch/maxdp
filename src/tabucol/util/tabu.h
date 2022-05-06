#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma once

#include <vector>

#include "util.h"
#include "Graph.h"

namespace gCol {  
  using namespace std;

  int tabu(const Graph & g, vector < int >&c, int k, int maxIterations, int verbose, int **neighbors);
}
