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
       
#include "util.h"
template <typename FlowT = double, typename CostT = double>
struct MinCostMaxFlow {
  MinCostMaxFlow(int n) : g(n + 1) {}
  pair<FlowT, CostT> solve(int s, int t);
  void addEdge(int u, int v, CostT wt, FlowT cap) {
    g[u].push_back(edge.size());
    edge.push_back({u, v, cap, 0, wt});
    g[v].push_back(edge.size());
    edge.push_back({v, u, 0, 0, -wt});
  }
  struct Edge {
    int from = -1, to = -1;
    FlowT cap = 0, flow = 0;
    CostT wt = 0;
  };
  vector<Edge> edge;
private:
  bool dijkstra(int s, int t);
  VI pred;
  vector<CostT> dist;
  VVI g;
};
template <typename FlowT, typename CostT>
pair<FlowT, CostT> MinCostMaxFlow<FlowT, CostT>::solve(int s, int t) {
  FlowT totFlow = 0;
  CostT totCost = 0;
  while (dijkstra(s, t)) {
    FlowT f = numeric_limits<FlowT>::max();
    for (int e = pred[t]; e != -1; e = pred[edge[e].from])
      f = min(f, edge[e].cap - edge[e].flow);
    if (f == 0) continue;
    for (int e = pred[t]; e != -1; e = pred[edge[e].from]) {
      edge[e].flow += f;
      edge[e ^ 1].flow -= f;
    }
    totFlow += f;
    totCost += (CostT)(f * dist[t]);
  }
  return make_pair(totFlow, totCost);
}
template <typename FlowT, typename CostT>
bool MinCostMaxFlow<FlowT, CostT>::dijkstra(int s, int t) {
  dist.assign(g.size(), -1);
  pred.assign(g.size(), -1);
  priority_queue<pair<CostT, int>> pq;
  dist[s] = 0.0;
  pq.emplace(0.0, s);
  while (pq.size()) {
    int v = pq.top().second;
    CostT w = -pq.top().first;
    pq.pop();
    if (dist[v] != w) continue;
    if (v == t) break;
    for (int e : g[v]) {
      int u = edge[e].to;
      CostT d = edge[e].wt + w;
      if (ff(edge[e].cap) <= ff(edge[e].flow)) continue;
      if (dist[u] == -1 or dist[u] > d) {
        pred[u] = e;
        dist[u] = d;
        pq.emplace(-d, u);
      }
    }
  }
  return pred[t] != -1;
}
