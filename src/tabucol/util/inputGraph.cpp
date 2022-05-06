#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "inputGraph.h"

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>

namespace gCol {  
  using namespace std;

  void inputDimacsGraph(Graph& g, const char *file) {
    ifstream in(file, ios::in);
    inputDimacsGraph(g,in);
    in.close();
  }
    
  void inputDimacsGraph(Graph& g, istream& in) {
    char c;
    char str[1000];
    int line = 0;
    g.nbEdges = 0;
    int edges = -1;
    int blem = 1;
    int multiple = 0;
    
    while (!in.eof()) {
      line++;
      in.get(c);
      if (in.eof())
	break;
      switch (c) {
      case 'p':
	in.get(c);
	in.getline(str, 999, ' ');
	if (strcmp(str, "edge") && strcmp(str, "edges")) {
	  cerr << " line " << line << ":\n";
	  cerr << "Error reading 'p' line: no 'edge' keyword found.\n";
	  cerr << "'" << str << "' found instead\n";
	  exit(-1);
	}
	in >> g.n;
	in >> edges;
	blem = 0;
	g.resize(g.n);
	break;
      case 'n':
	if (blem) {
	  cerr << " line " << line << ":\n";
	  cerr << "Found 'n' line before a 'p' line.\n";
	  exit(-1);
	}
	int node;
	in >> node;
	if (node < 1 || node > g.n) {
	  cerr << " line " << line << ":\n";
	  cerr << "Node number " << node << " is out of range!\n";
	  exit(-1);
	}
	node--;
	cout << "Tags (n Lines) not implemented in g object\n";
	break;
      case 'e':
	int node1, node2;
	in >> node1 >> node2;
	if (node1 < 1 || node1 > g.n || node2 < 1 || node2 > g.n) {
	  cerr << " line " << line << ":\n";
	  cerr << "Node " << node1 << " or " << node2 << " is out of range!\n";
	  exit(-1);
	}
	node1--;
	node2--;
	// Undirected graph
	// cout << "Read edge from " << node1 << " to " << node2 << endl;
	if (g[node1][node2] == 0) {
	  g.nbEdges++;
	} else {
	  multiple++;
	  if (multiple < 5) {
	    cerr << "Warning: " << " at line " << line << ": edge is defined more than once.\n";
	    if (multiple == 4) {
	      cerr << "  No more multiple edge warnings will be issued\n";
	    }
	  }
	}
	g[node1][node2] = 1;
	g[node2][node1] = 1;
	break;
      case 'd':
      case 'v':
      case 'x':
	cerr << " line " << line << ":\n";
	cerr << "'" << c << "' lines are not implemented yet...\n";
	in.getline(str, 999, '\n');
	break;
      case 'c':
	in.putback('c');
	in.get(str, 999, '\n');
	break;
      default:
	cerr << " line " << line << ":\n";
	cerr << "'" << c << "' is an unknown line code\n";
	exit(-1);
      }
      in.get();			// Kill the newline;
    }
    if (multiple)
      cerr << multiple << " multiple edges encountered\n";
  }

  // simple wrapper around existing function, TBD: improve
  void readInputFile(istream& inStream, int& numNodes, int& numEdges, vector<vector<bool>>& adjacent, vector<int>& degree, vector<vector<int>>& adjList) {
    Graph g;
    inputDimacsGraph(g, inStream);
    numNodes = g.n;
    numEdges = g.nbEdges;

    // copy over
    adjacent.clear();
    adjacent.resize(numNodes, vector<bool>(numNodes));
    for (int i = 0; i < numNodes; i++)
      for (int j = 0; j < numNodes; j++)
	adjacent[i][j]=g[i][j];
    for (int i = 0; i < numNodes; i++)
      adjacent[i][i]=true;
      
    //Now use the adjacency matrix to construct the degree array and adjacency list
    adjList.clear();
    adjList.resize(numNodes, vector<int>());
    makeAdjList(adjList,g);
      
    degree.resize(numNodes, 0);
    for (int i = 0; i < numNodes; i++) {
      for (int j = 0; j < numNodes; j++) {
	if (adjacent[i][j] && i != j) {
	  degree[i]++;
	}
      }
    }
  }
}

