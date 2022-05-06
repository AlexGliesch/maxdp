#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
/**********************************************************************
 * \file solution.h
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 *   \version $Id: emacs 7968 2017-05-17 23:03:17Z ritt $
 *   \date Time-stamp: <2018-06-29 17:50:01 ritt>
 **********************************************************************/
#pragma once

#include <cassert>
#include <iostream>
#include <vector>

namespace gCol {
using namespace std;
void prettyPrintSolution(vector<vector<int>>& candSol);
void writeSolution(vector<int>&, string file = "solution.txt");
void writeGrouping(int n, vector<vector<int>>& candSol,
                   string file = "solution.txt");
void group2sol(vector<vector<int>>& candSol, vector<int>& coloring);

template <typename AdjList>
void checkSolution(vector<vector<int>>& candSol, AdjList& adjacent,
                   int numNodes) {
  assert(numNodes >= 0);
  unsigned i, count = 0, group;
  bool valid = true;

  // first check that the permutation is the right length
  for (group = 0; group < candSol.size(); group++)
    count = count + candSol[group].size();

  if (count != unsigned(numNodes)) {
    cout << "Error: Permutations length is not equal to the problem size\n";
    valid = false;
  }

  // Now check that all the nodes are in the permutation (EXPENSIVE)
  vector<int> a(numNodes, 0);
  for (group = 0; group < candSol.size(); group++)
    for (i = 0; i < candSol[group].size(); i++)
      a[candSol[group][i]]++;

  for (i = 0; i < unsigned(numNodes); i++) {
    if (a[i] != 1) {
      cout << "Error: Vertex " << i
           << " should be in the solution exactly once.\n";
      valid = false;
    }
  }

  // Finally, check for illegal colourings: I.e. check that each colour class
  // contains non conflicting nodes
  for (group = 0; group < candSol.size(); group++) {
    for (i = 0; i < candSol[group].size() - 1; i++) {
      for (unsigned j = i + 1; j < candSol[group].size(); j++) {
        if (adjacent[candSol[group][i]][candSol[group][j]]) {
          cout << "Error: Nodes " << candSol[group][i] << " and "
               << candSol[group][j] << " are in the same group, but they clash"
               << endl;
          valid = false;
        }
      }
    }
  }

  if (!valid) cout << "This solution is not valid" << endl;
}
} // namespace gCol
