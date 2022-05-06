#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma once

#include <vector>

namespace gCol {
  using namespace std;
  
  class Graph {

  public:

    Graph();
    Graph(int n);
    ~Graph();

    void resize(int n);

    int *matrix;
    int n;	  // number of nodes
    int nbEdges;  // number of edges

    int *operator[] (int index) const;
  };

  void makeAdjList(int** neighbors, const Graph& g);
  void makeAdjList(vector<vector<int>>& neighbors, const Graph& g);
}
