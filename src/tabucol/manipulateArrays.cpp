#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "manipulateArrays.h"
#include <iostream>

#include "util/manipulateArrays.h"

namespace gCol {
  namespace pctc {
    void moveNodeToColor(int bestNode, int bestColor, const Graph& g, int *c, int **nodesByColor, int **conflicts, int *nbcPosition, int **neighbors,
			 int **tabuStatus, long totalIterations, int tabuTenure) {

      // move bestNodes to bestColor
      c[bestNode] = bestColor;
      // Replace bestNode by the last node in the nodesByColor[0] array and shorten it
      nodesByColor[0][nbcPosition[bestNode]] = nodesByColor[0][nodesByColor[0][0]--];
      // Update the nbcPosition array the node that has taken the place of best node
      nbcPosition[nodesByColor[0][nbcPosition[bestNode]]] = nbcPosition[bestNode];
      // Insert bestNode into the nodesByColor[bestColor] array, increase its lenght
      // and update the nbcPosition for bestNode
      nodesByColor[bestColor][(nbcPosition[bestNode] = ++nodesByColor[bestColor][0])] = bestNode;

      // Update the conflicts array and remove conflicting nodes
      numConfChecks++;
      for (int j = 1; j <= neighbors[bestNode][0]; j++) {
	int i = neighbors[bestNode][j];
	numConfChecks++;

	// Do not move neighbors to bestColor for a couple of iterations in order to
	// avoid bestNode from dropping back out too soon
	tabuStatus[i][bestColor] = totalIterations + tabuTenure;

	// Increase the conflicts for bestColor
	conflicts[bestColor][i]++;
	numConfChecks++;

	// Check for conflict created by moving bestNode to bestColor
	if (c[i] == bestColor) {
	  // If so, remove the node and put it back to the uncolored nodes
	  nodesByColor[bestColor][nbcPosition[i]] = nodesByColor[bestColor][nodesByColor[bestColor][0]--];
	  nbcPosition[nodesByColor[bestColor][nbcPosition[i]]] = nbcPosition[i];
	  nodesByColor[0][(nbcPosition[i] = ++nodesByColor[0][0])] = i;
	  c[i] = 0;
	  // Reduce the conflicts of all neighbors.
	  numConfChecks++;
	  for (int k = 1; k <= neighbors[i][0]; k++) {
	    conflicts[bestColor][neighbors[i][k]]--;
	    numConfChecks += 2;
	  }
	}
      }
    }
  }
}
