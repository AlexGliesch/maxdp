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
#include "coudert.h"
namespace coudert {
vector<VertexInformation::NodeID>
maxCliqueRec(const Graph& G, vector<VertexInformation::NodeID>& C,
             vector<VertexInformation::NodeID>& B, unsigned ub, unsigned level);
bool remove_first_non_adjacent(Graph& G, Node& v) {
  assert(v < num_vertices(G));
  graph_traits<Graph>::vertex_iterator vp, ve;
  for (tie(vp, ve) = vertices(G); vp != ve; vp++)
    if (*vp != v && !is_adjacent(G, v, *vp)) {
      clear_vertex(*vp, G);
      remove_vertex(*vp, G);
      return true;
    }
  return false;
}
void keep_neighborhood(Graph& G, Node v) {
  VertexInformation::NodeID id = G[v].id;
  while (remove_first_non_adjacent(G, v)) {
    if (G[v].id != id || !(v < num_vertices(G))) v--;
    assert(G[v].id == id && v < num_vertices(G));
  }
  clear_vertex(v, G);
  remove_vertex(v, G);
}
Node maxDegree(const Graph& G) {
  assert(num_vertices(G) > 0);
  graph_traits<Graph>::vertex_iterator vp, ve;
  Node r = 0;
  unsigned max_degree = 0;
  double max_central = 0;
  for (tie(vp, ve) = vertices(G); vp != ve; vp++) {
    if (degree(*vp, G) > max_degree ||
        (degree(*vp, G) == max_degree && G[*vp].centrality > max_central)) {
      max_degree = degree(*vp, G);
      max_central = G[*vp].centrality;
      r = *vp;
    }
  }
  return r;
}
pair<unsigned, unsigned> saturation(const Graph& G, const vector<Color>& color,
                                    Node v) {
  set<Color> c;
  unsigned uncolored = 0;
  graph_traits<Graph>::adjacency_iterator up, ue;
  for (tie(up, ue) = adjacent_vertices(v, G); up != ue; up++)
    if (color[*up] != 0)
      c.insert(color[*up]);
    else
      uncolored++;
  return make_pair(c.size(), uncolored);
}
Node maxSatur(const Graph& G, vector<Color>& color) {
  assert(num_vertices(G) > 0);
  assert(num_vertices(G) <= color.size());
  graph_traits<Graph>::vertex_iterator vp, ve;
  Node r = 0;
  int max_satur = -1;
  unsigned max_uncolored = 0;
  double max_central = 0;
  for (tie(vp, ve) = vertices(G); vp != ve; vp++) {
    if (color[*vp] == 0) {
      int satur;
      unsigned uncolored;
      tie(satur, uncolored) = saturation(G, color, *vp);
      if (satur > max_satur ||
          (satur == max_satur && uncolored > max_uncolored) ||
          (satur == max_satur && uncolored == max_uncolored &&
           G[*vp].centrality > max_central)) {
        max_satur = satur;
        max_uncolored = uncolored;
        max_central = G[*vp].centrality;
        r = *vp;
      }
    }
  }
  return r;
}
unsigned color_H3(const Graph& G, vector<Color>& color) {
  assert(color.size() >= num_vertices(G));
  assert(all_of(color.begin(), color.begin() + num_vertices(G),
                [](unsigned c) { return c == 0; }));
  unsigned colored = 0;
  unsigned k = 0;
  while (colored < num_vertices(G)) {
    Node v = maxSatur(G, color);
    for (unsigned c = 1; c <= k + 1; c++) {
      graph_traits<Graph>::adjacency_iterator up, ue;
      for (tie(up, ue) = adjacent_vertices(v, G); up != ue; up++)
        if (color[*up] == c) break;
      if (up != ue) continue;
      color[v] = c;
      colored++;
      k = max(c, k);
      break;
    }
  }
  assert( k <= num_vertices(G));
  assert(all_of(color.begin(), color.begin() + num_vertices(G),
                [](unsigned c) { return c != 0; }));
  return k;
}
unsigned color_H3(const Graph& G) {
  vector<Color> color(num_vertices(G));
  return color_H3(G, color);
}
bool rule_B(Graph& G, int limit) {
  graph_traits<Graph>::vertex_iterator vp, ve;
  for (tie(vp, ve) = vertices(G); vp != ve; vp++)
    if (int(degree(*vp, G)) < limit) {
      clear_vertex(*vp, G);
      remove_vertex(*vp, G);
      return true;
    }
  return false;
}
bool rule_C(Graph& G, vector<VertexInformation::NodeID>& C) {
  graph_traits<Graph>::vertex_iterator vp, ve;
  for (tie(vp, ve) = vertices(G); vp != ve; vp++)
    if (degree(*vp, G) >= num_vertices(G) - 2) {
      C.push_back(G[*vp].id);
      keep_neighborhood(G, *vp);
      return true;
    }
  return false;
}
bool rule_Q(Graph& G, const vector<Color>& color, const unsigned k,
            const int limit) {
  graph_traits<Graph>::vertex_iterator vp, ve;
  for (tie(vp, ve) = vertices(G); vp != ve; vp++)
    if (int(k) - int(saturation(G, color, *vp).first) > limit) {
      clear_vertex(*vp, G);
      remove_vertex(*vp, G);
      return true;
    }
  return false;
}
unsigned maxClique_backtracks;
vector<VertexInformation::NodeID> maxClique(const Graph& G) {
  TIME_BLOCK("coudert::maxClique");
  maxClique_backtracks = 0;
  vector<VertexInformation::NodeID> C, B;
  return maxCliqueRec(G, C, B, numeric_limits<unsigned>::max(), 0);
}
vector<VertexInformation::NodeID>
maxCliqueRec(const Graph& G, vector<VertexInformation::NodeID>& C_,
             vector<VertexInformation::NodeID>& B, unsigned ub,
             unsigned level) {
  Graph G_(G);
  vector<VertexInformation::NodeID> C(C_);
  unsigned k;
  vector<Color> color(num_vertices(G));
  do {
    if (num_vertices(G_) == 0) return C;
    if (C.size() + num_vertices(G_) <= B.size()) return B;
    fill(color.begin(), color.end(), 0);
    k = color_H3(G_, color);
    ub = min(ub, unsigned(C.size() + k));
    if (ub <= B.size()) return B;
  } while (rule_B(G_, int(B.size()) - int(C.size())) || rule_C(G_, C) ||
           rule_Q(G_, color, k, int(C.size()) - int(B.size()) + int(k)));
  unsigned Bold = B.size();
  Node v = maxDegree(G_);
  {
    C.push_back(G_[v].id);
    Graph G1(G_);
    assert(num_vertices(G1) == num_vertices(G_));
    assert(G_[v].id == G1[v].id);
    keep_neighborhood(G1, v);
    assert(num_vertices(G1) == degree(v, G_));
    B = maxCliqueRec(G1, C, B, ub, level + 1);
    C.pop_back();
  }
  if (ub == B.size()) return B;
  Graph G0(G_);
  clear_vertex(v, G0);
  remove_vertex(v, G0);
  B = maxCliqueRec(G0, C, B, ub, level + 1);
  if (B.size() == Bold) maxClique_backtracks++;
  return B;
}
unsigned seqColorRec(Graph& G, unsigned k, unsigned B, const unsigned lb,
                     vector<Color>& color, unsigned colors, unsigned& back,
                     unsigned m, Timer timer);
pair<unsigned, unsigned> seqColor(Graph& G,
                                  const vector<VertexInformation::NodeID>& C,
                                  unsigned m, Timer timer) {
  TIME_BLOCK("coudert::seqColor");
  unsigned chi = 0;
  unsigned back = 0;
  vector<Color> color(num_vertices(G));
  for (unsigned i = 0; i < C.size(); i++)
    color[C[i]] = i + 1;
  chi = seqColorRec(G, C.size(), num_vertices(G) + 1, C.size(), color, C.size(),
                    back, m, timer);
  return make_pair(chi, back);
}
unsigned seqColorRec(Graph& G, unsigned k, unsigned B, const unsigned lb,
                     vector<Color>& color, unsigned colored, unsigned& back,
                     unsigned m, Timer timer) {
  assert(num_vertices(G) == color.size());
  if (timer.timedOut()) throw TimeoutException(timer);
  if (B <= m) return B;
  if (num_vertices(G) == colored) {
    if (k < B) {
      graph_traits<Graph>::vertex_iterator vp, ve;
      for (tie(vp, ve) = vertices(G); vp != ve; vp++)
        G[*vp].color = color[*vp];
    }
    return k;
  }
  if (k >= B) return B;
  Node v = maxSatur(G, color);
  G[v].rank = colored + 1;
  for (Color c = 1; c <= min(k + 1, B - 1); c++) {
    graph_traits<Graph>::adjacency_iterator up, ue;
    for (tie(up, ue) = adjacent_vertices(v, G); up != ue; up++)
      if (color[*up] == c) break;
    if (up != ue) continue;
    color[v] = c;
    B = seqColorRec(G, max(c, k), B, lb, color, colored + 1, back, m, timer);
    color[v] = 0;
    if (timer.timedOut()) throw TimeoutException(timer);
    if (lb == B || k >= B || B <= m) return B;
  }
  back++;
  return B;
}
}
