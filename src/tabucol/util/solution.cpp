#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
/**
 * \file solution.cpp
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br>
 *   \version $Id: emacs 7968 2017-05-17 23:03:17Z ritt $
 *   \date Time-stamp: <2018-06-29 17:49:43 ritt>
 */
#include "solution.h"

#include <fstream>
#include <iostream>

namespace gCol {
void prettyPrintSolution(vector<vector<int>>& candSol) {
  unsigned i, count = 0, group;
  cout << endl << endl;
  for (group = 0; group < candSol.size(); group++) {
    cout << "C- " << group << " = {";
    if (candSol[group].size() == 0)
      cout << "empty}\n";
    else {
      for (i = 0; i < candSol[group].size() - 1; i++)
        cout << candSol[group][i] << ", ";
      cout << candSol[group][candSol[group].size() - 1]
           << "} |G| = " << candSol[group].size() << endl;
      count = count + candSol[group].size();
    }
  }
  cout << "Total Number of Nodes = " << count << endl;
}

// output the solution to a text file
void writeSolution(vector<int>& S, string file) {
  ofstream solStrm(file);
  solStrm << S.size() << endl;
  for (int color : S)
    solStrm << color << endl;
  solStrm.close();
}

void writeGrouping(int n, vector<vector<int>>& candSol, string file) {
  vector<int> coloring(n);
  group2sol(candSol, coloring);
  writeSolution(coloring, file);
}

void group2sol(vector<vector<int>>& candSol, vector<int>& coloring) {
  for (unsigned i = 0; i < candSol.size(); i++)
    for (unsigned j = 0; j < candSol[i].size(); j++)
      coloring[candSol[i][j]] = i;
}
} // namespace gCol
