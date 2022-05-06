#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma once
#include "../util/Graph.h"

#include <istream>

namespace gCol {
  using namespace std;
  void inputDimacsGraph(Graph& g, istream&);
  void inputDimacsGraph(Graph& g, const char *filename);
  void readInputFile(istream& inStream, int& numNodes, int& numEdges, vector<vector<bool>>& adjacent, vector<int>& degree, vector<vector<int>>& adjList);
}
