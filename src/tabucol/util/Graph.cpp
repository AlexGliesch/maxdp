#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "Graph.h"
#include <fstream>
#include <string.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>

namespace gCol {  
  using namespace std;

  Graph::Graph() {
    matrix = nullptr;
    n = 0;
    nbEdges = 0;
  }

  Graph::Graph(int n) {
    matrix = nullptr;
    resize(n);
  }

  int *Graph::operator[] (int index) const {
    if (index < 0 || index >= this->n) {
      cerr << "First node index out of range: " << index << "\n";
      matrix[-1] = 0;		//Make it crash.
    }
    return matrix + this->n * index;
  }

  void Graph::resize(int _n) {
    if (matrix != nullptr)
      delete [] matrix;

    if (_n > 0) {
      n = _n;
      nbEdges = 0;
      matrix = new int[n * n];
      for (int i = 0; i < n*n; i++)
	matrix[i] = 0;
    }
  }

  Graph::~Graph() {
    resize(0);
  }

  // Makes the adjacency list corresponding to G
  void makeAdjList(int **neighbors, const Graph& g) {
    for (int i = 0; i < g.n; i++) {
      neighbors[i] = new int[g.n + 1];
      neighbors[i][0] = 0;
    }
    for (int i = 0; i < g.n; i++)
      for (int j = 0; j < g.n; j++)
	if (g[i][j] && i != j)
	  neighbors[i][++neighbors[i][0]] = j;
  }

  // Makes the adjacency list corresponding to G
  void makeAdjList(vector<vector<int>>& neighbors, const Graph& g) {
    for (int i = 0; i < g.n; i++)
      for (int j = 0; j < g.n; j++)
	if (g[i][j] && i != j)
	  neighbors[i].push_back(j);
  }
}
