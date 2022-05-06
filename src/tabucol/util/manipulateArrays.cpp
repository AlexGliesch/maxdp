#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "manipulateArrays.h"

#include <iostream>

namespace gCol {
  void freeArrays(int **&nodesByColor, int **&conflicts, int **&tabuStatus, int *&nbcPosition, int k, int n) {
    for (int i = 0; i <= k; i++) {
      delete [] nodesByColor[i];
      delete [] conflicts[i];
    }
    for (int i = 0; i < n; i++)
      delete [] tabuStatus[i];

    delete [] nodesByColor;
    delete [] conflicts;
    delete [] tabuStatus;
    delete [] nbcPosition;
  }
}
