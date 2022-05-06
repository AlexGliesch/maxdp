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
#include <ctime>
#include <random>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif
using namespace std;
inline mt19937 rng;
inline int randInt(int from, int to) {
  static uniform_int_distribution<int> d;
  return d(rng, decltype(d)::param_type{from, to});
}
inline double randDouble(double from, double to) {
  static std::uniform_real_distribution<double> d;
  return d(rng, decltype(d)::param_type{from, to});
}
inline size_t uniqueRandomSeed() {
  size_t a = (size_t)clock(), b = (size_t)time(nullptr);
#ifdef _WIN32
  size_t c = (size_t)GetCurrentProcessId();
#else
  size_t c = (size_t)getpid();
#endif
  a = (a - b - c) ^ (c >> 13);
  b = (b - c - a) ^ (a << 8);
  c = (c - a - b) ^ (b >> 13);
  a = (a - b - c) ^ (c >> 12);
  b = (b - c - a) ^ (a << 16);
  c = (c - a - b) ^ (b >> 5);
  a = (a - b - c) ^ (c >> 3);
  b = (b - c - a) ^ (a << 10);
  c = (c - a - b) ^ (b >> 15);
  return c;
}
template <typename T>
inline void randomChoice(const vector<T>& v, vector<T>& result, int k) {
  int n = v.size();
  result.resize(k);
  for (int i = 0; i < k; ++i)
    result[i] = v[i];
  for (int i = k + 1; i < n; ++i) {
    int j = randInt(0, i);
    if (j < k) result[j] = v[i];
  }
}
struct ReservoirSampling {
  bool consider() { return (1.0 / ++num) >= randDouble(0.0, 1.0); }
  template <typename T, typename F>
  void considerV(const T& v1, const T& v2, const F& f) {
    bool eq = v1 == v2, le = v1 < v2;
    if (le or (eq and consider())) {
      f();
      if (le) reset();
    }
  }
  void reset(double num = 0.0) { this->num = num; }
  double num = 0.0;
};
