#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
/******************************************************************************/
//  This code implements the PartialCol, and Tabucol graph coloring
//  heuristics described in "A Reactive Tabu Search Using Partial Solutions for
//  the Graph Coloring Problem" by Ivo Bloechliger and Nicolas Zuffery.
//
//  The code was originally written by Ivo Bloechliger
//  http://rose.epfl.ch/~bloechli/coloring/
//  but has been modified for running with multiple values of k and for
//  algorithm perfromance analysis
//
//      See: Lewis, R. (2015) A Guide to Graph Colouring: Algorithms and
//      Applications. Berlin, Springer.
//       ISBN: 978-3-319-25728-0. http://www.springer.com/us/book/9783319257280
//
//      for further details
/******************************************************************************/

#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <string.h>

#include <sys/ioctl.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "PartialColAndTabuCol.h"
#include "initializeColoring.h"
#include "manipulateArrays.h"
#include "reactcol.h"
#include "tabu.h"
#include "util/inputGraph.h"
#include "util/solution.h"

namespace gCol {
using namespace std;
namespace pctc {
int pctc(const Graph& g, vector<int>& best, confTimeLogger& log, int algorithm,
         int tenure, unsigned long long maxChecks, int targetCols,
         int randomSeed, int verbose, int constructiveAlg, int maxIterations) {

  int k, frequency = 0, increment = 0, cost, duration;

  // Produce some output
  vstream(1) << " COLS     CPU-TIME\tCHECKS" << endl;

  // Make the adjacency list structure
  int** neighbors = new int*[g.n];
  makeAdjList(neighbors, g);

  // The solution is held in the following array
  int* coloring = new int[g.n];
  int* bestColouring = new int[g.n];

  // Now start the timer
  timer tm;

  // Generate the initial value for k using greedy or dsatur algorithm
  k = generateInitialK(g, constructiveAlg, bestColouring);
  //..and write the results to the output file
  duration = tm.elapsed_ms();
  vstream(1) << setw(5) << k << setw(11) << duration << "ms\t" << numConfChecks
             << " (via constructive)" << endl;
  log.found(k, numConfChecks, duration);

  // MAIN ALGORITHM
  int bestK = k;
  k--;
  while (numConfChecks < maxChecks && k + 1 > targetCols) {

    // Initialise the solution array
    fill(&coloring[0], &coloring[g.n], 0);

    // Do the algorithm for this value of k, either until a slution is found, or
    // maxChecks is exceeded
    if (algorithm == 1)
      cost = reactcol(g, coloring, k, maxChecks, tenure, verbose, frequency,
                      increment, neighbors);
    else
      cost = tabu(g, coloring, k, maxChecks, maxIterations, tenure, verbose,
                  frequency, increment, neighbors);

    // Algorithm has finished at this k
    duration = tm.elapsed_ms();
    if (cost == 0) {
      vstream(1) << setw(5) << k << setw(11) << duration << "ms\t"
                 << numConfChecks << endl;
      log.found(k, numConfChecks, duration);
      bestK = min(bestK, k);
      // Copy the current solution as the best solution
      for (int i = 0; i < g.n; i++)
        bestColouring[i] = coloring[i] - 1;
      // Check if the target has been met
      if (k <= targetCols) {
        vstream(1) << "\nSolution with <=" << k
                   << " colours has been found. Ending..." << endl;
        log.target(1);
        break;
      }
    } else {
      vstream(1) << "\nRun limit exceeded. No solution using " << k
                 << " colours was achieved (Checks = " << numConfChecks << ", "
                 << duration << "ms)" << endl;
      log.fail(k, numConfChecks, duration);
    }

    // Decrement k (if the run time hasn't been reached, we'll carry on with
    // this new value)
    k--;
  }
  copy(&bestColouring[0], &bestColouring[g.n], best.begin());
  delete[] coloring;
  delete[] bestColouring;
  // Delete the arrays and end
  for (int i = 0; i < g.n; i++)
    delete[] neighbors[i];
  delete[] neighbors;

  return bestK;
}

int main(int argc, char** argv) {
  // (0) process commandline options

  unsigned long long maxChecks;
  int algorithm = 1, tenure = 0, randomSeed, targetCols, constructiveAlg,
      timeLimitSeconds;
  bool printColors;
  string solution;

  struct winsize wsize;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &wsize);
  options_counter verbosec;

  po::options_description desc("Options", wsize.ws_col);
  desc.add_options()("help", "Show help")(
      "instance", po::value<string>(), "Coloring instance in DIMACS format.")(
      "tabucol,t", po::bool_switch(),
      "If present, TabuCol is used. Else PartialCol is used.")(
      "dynamic,d", po::bool_switch(),
      "If present, a "
      "dynamic tabu "
      "tenure is used "
      "(i.e. "
      "tabuTenure = "
      "(int)(0.6*nc) "
      "+ rand(0,9)). "
      "Otherwise a "
      "reactive "
      "tenure is "
      "used.") // TBD:
               // differs!
      ("stop,s",
       po::value<unsigned long long>(&maxChecks)->default_value(100000000),
       "Maximum number of constraint checks.")(
          "random,r", po::value<int>(&randomSeed)->default_value(1),
          "Random seed.")(
          "target,T", po::value<int>(&targetCols)->default_value(2),
          "Target number of colours. Algorithm halts if this is reached.")(
          "verbose,v", po::value(&verbosec)->zero_tokens(),
          "Verbosity. If present, output is sent to screen. If -v is repeated, "
          "more output is given.")(
          "algorithm,a", po::value<int>(&constructiveAlg)->default_value(1),
          "Choice of construction algorithm to determine initial value for k. "
          "DSsatur = 1, Greedy = 2.")(
          "timelimit", po::value<int>(&timeLimitSeconds)->default_value(3600),
          "Time limit in seconds.") // fix: only code where -t is not for the
                                    // timelimit
      ("printcolors", po::bool_switch(&printColors)->default_value(false),
       "Print the best number of colors found in standard output.")(
          "solution",
          po::value<string>(&solution)->default_value("solution.txt"),
          "Name of solution file.");

  po::positional_options_description pod;
  pod.add("instance", 1);

  po::variables_map vm;
  po::store(
      po::command_line_parser(argc, argv).options(desc).positional(pod).run(),
      vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("instance")) {
    cout << desc << endl;
    return 1;
  }

  int verbose = verbosec.count;

  if (vm["tabucol"].as<bool>()) algorithm = 2;
  if (vm["dynamic"].as<bool>()) tenure++;
  // Read in program parameters
  Graph g;
  vstream(1) << "PartialCol/TabuCol Algorithm* using <"
             << vm["instance"].as<string>() << ">" << endl
             << endl;
  inputDimacsGraph(g, vm["instance"].as<string>().c_str());

  if (targetCols < 2 || targetCols > g.n) targetCols = 2;

  // This variable keeps count of the number of times information about the
  // instance is looked up
  numConfChecks = 0;

  // Seed
  srand(randomSeed);

  // Now set up some output files
  confTimeLogger log;

  // Do a check to see if we have the empty graph. If so, end immediately.
  if (g.nbEdges <= 0) {
    log.found(1, 0, 0);
    log.fail(0, 0, 0);
    vstream(1) << "Graph has no edges. Optimal solution is obviously using one "
                  "colour. Exiting."
               << endl;
    exit(1);
  }

  vector<int> bestColouring(g.n);
  unsigned bestK = pctc(g, bestColouring, log, algorithm, tenure, maxChecks,
                        targetCols, randomSeed, verbose, constructiveAlg);

  writeSolution(bestColouring, solution);

  if (printColors) cout << bestK << endl;

  return 0;
}
} // namespace pctc
} // namespace gCol
