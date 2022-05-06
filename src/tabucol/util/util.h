#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma once

#include <chrono>
#include <limits>
#include <utility>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

#include <boost/any.hpp>

#define vstream(l) if (verbose>=(l)) cout

extern unsigned long long numConfChecks;

//// Type aliases
using nld = std::numeric_limits<double>;

// helper class for output
struct confTimeLogger {
  std::ofstream conf, time;

  confTimeLogger() {
    time.open("teffort.txt");
    conf.open("ceffort.txt");
  }
  ~confTimeLogger() {
    time.close();
    conf.close();
  }
  void found(int k, unsigned long long numConfChecks, int duration) {
    conf << k << "\t" << numConfChecks << std::endl;
    time << k << "\t" << duration << std::endl;
  }
  void target(int k) {
    conf << k << "\tX" << std::endl;
    time << k << "\tX" << std::endl;
  }
  void fail(int k, unsigned long long numConfChecks, int duration) {
    conf << k << "\tX\t" << numConfChecks << std::endl;
    time << k << "\tX\t" << duration << std::endl;
  }
};

// Time-tracking structure
class timer {
  using clock = std::chrono::steady_clock;
  using timepoint = clock::time_point;
  timepoint tpstart;
  double tmlim;

 public:
  // Initialize a timer with a given time limit
  timer(double time_lim_secs = nld::max()) { reset(time_lim_secs); }
  // Initialize a timer with a given time limit or with a "parent timer": the
  // time limit will be min(tmlim, parent.secs_left())
  timer(double time_lim_secs, const timer& parent) {
    reset(std::min(time_lim_secs, parent.secs_left()));
  }
  // Resets a timer with a given time limit
  void reset(double time_lim_secs = nld::max()) {
    tmlim = time_lim_secs, tpstart = clock::now();
  }
  // Elapsed seconds since the start
  double elapsed_secs() const {
    return std::chrono::duration_cast<std::chrono::duration<double>>(clock::now()-tpstart).count();
  }
  unsigned elapsed_ms() const {
    return std::chrono::duration_cast<std::chrono::milliseconds>(clock::now()-tpstart).count();
  }
  // Number of seconds left until timeout
  double secs_left() const { return tmlim - elapsed_secs(); }
  // Whether it timed out
  bool timed_out() const { return elapsed_secs() >= tmlim; }
};

struct options_counter { int count = 0; };
inline void validate(boost::any& v, std::vector<std::string> const& xs, options_counter*, long) {
  if (v.empty()) v = options_counter{1};
  else ++boost::any_cast<options_counter&>(v).count;
}

