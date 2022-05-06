#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
/**********************************************************************
 * \file color.h
 *   \author Marcus Ritt <marcus.ritt@inf.ufrgs.br> 
 *   \version $Id: emacs 7968 2017-05-17 23:03:17Z ritt $
 *   \date Time-stamp: <2018-06-22 00:03:01 ritt>
 **********************************************************************/
#pragma once

#include <vector>

namespace gCol {
  using namespace std;
  
  // check if vertex `v` can be colored `c` in solution `sol`
  template <typename AdjMatrix>
    inline bool colourIsFeasible(int v, const vector<vector<int>>& sol, int c, const vector<int>& colNode, const vector<vector<int>>& adjList, const AdjMatrix& adj) {
    unsigned i;
    numConfChecks++;
    if (sol[c].size() > adjList[v].size()) {
      //check if any neighbours of v are currently in colour c
      for (i = 0; i < adjList[v].size(); i++) {
	numConfChecks++;
	if (colNode[adjList[v][i]] == c)
	  return false;
      }
      return true;
    } else {
      //check if any vertices in colour c are adjacent to v
      for (i = 0; i < sol[c].size(); i++) {
	numConfChecks++;
	if (adj[v][sol[c][i]])
	  return false;
      }
      return true;
    }
  }

  template <class AdjList>
  inline
    void assignAColourDSatur(bool& foundColour, vector<vector<int>>& candSol, const AdjList& adjacent,
			     const vector<int>& permutation, int nodePos, vector<int>& satDeg, const vector<vector<int>>& adjList, vector<int>& colNode) {
    int i, j, c = 0, v = permutation[nodePos];
    bool alreadyAdj;

    while (c < candSol.size() && !foundColour) {
      //check if colour c is feasible for vertex v
      if (colourIsFeasible(v, candSol, c, colNode, adjList, adjacent)) {
	//v can be added to this colour
	foundColour = true;
	candSol[c].push_back(v);
	colNode[v] = c;
	//We now need to update satDeg. To do this we identify the uncoloured nodes i that are adjacent to
	//this newly coloured node v. If i is already adjacent to a node in colour c we do nothing, 
	//otherwise its saturation degree is increased...
	for (i = 0; i < satDeg.size(); i++) {
	  numConfChecks++;
	  if (adjacent[v][permutation[i]]) {
	    alreadyAdj = false;
	    j = 0;
	    while (j < candSol[c].size() - 1 && !alreadyAdj) {
	      numConfChecks++;
	      if (adjacent[candSol[c][j]][permutation[i]])
		alreadyAdj = true;
	      j++;
	    }
	    if (!alreadyAdj)
	      satDeg[i]++;
	  }
	}
      }
      c++;
    }
  }
}
