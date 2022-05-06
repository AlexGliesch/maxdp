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
       
#include "main.h"
#include "pre.h"
namespace coudert {
using namespace std;
using namespace boost;
struct VertexInformation {
  VertexInformation() : id(0), color(0) {}
  typedef int NodeID;
  VertexInformation(NodeID _id) : id(_id), color(0) {}
  NodeID id;
  unsigned color;
  double centrality;
  bool clique;
  unsigned rank;
};
typedef adjacency_list<vecS, vecS, undirectedS, VertexInformation> Graph;
typedef graph_traits<Graph>::vertex_descriptor Node;
typedef unsigned Color;
inline Graph readDIMACS(ifstream& f) {
  char type;
  string s;
  Graph G;
  while (!f.eof()) {
    f >> type;
    if (f.eof()) break;
    switch (type) {
    case 'c':
      getline(f, s);
      break;
    case 'p':
      unsigned n, e;
      f >> s >> n >> e;
      for (unsigned i = 0; i < n; i++)
        add_vertex(VertexInformation(i), G);
      break;
    case 'n':
      getline(f, s);
      break;
    case 'e':
      unsigned src, dst;
      f >> src >> dst;
      if (!is_adjacent(G, src - 1, dst - 1)) add_edge(src - 1, dst - 1, G);
      break;
    default:
      throw("File format error");
    }
  }
  return G;
}
class vertex_writer {
  const Graph& G;
public:
  enum GraphType { none, queen };
  GraphType t;
  vertex_writer(const Graph& _G, GraphType _t = none) : G(_G), t(_t) {}
  template <typename VertexOrEdge>
  void operator()(ostream& out, const VertexOrEdge& v) const {
    out << "[";
    if (G[v].clique)
      out << "shape=box";
    else
      out << "shape=circle";
    out << ",label=\"" << G[v].rank
        << "\",fontcolor=\"/white\",fillcolor=" << G[v].color;
    if (t == queen) {
      const unsigned scale = 60;
      unsigned n = sqrt(num_vertices(G));
      out << ",pos=\"" << scale * (G[v].id % n) << "," << scale * (G[v].id / n)
          << "\"";
    }
    out << "]";
  }
};
struct graph_writer {
  const Graph& G;
public:
  graph_writer(const Graph& _G) : G(_G) {}
  void operator()(std::ostream& out) const {
    out << "graph [outputorder=edgesfirst]" << std::endl;
    out << "node  [colorscheme=spectral11,style=filled]" << std::endl;
  }
};
vector<VertexInformation::NodeID> maxClique(const Graph& G);
pair<unsigned, unsigned> seqColor(Graph& G,
                                  const vector<VertexInformation::NodeID>& C,
                                  unsigned m, Timer timer);
}
