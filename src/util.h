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
       
#include "pre.h"
#include "random.h"
using u64 = uint64_t;
using i64 = int64_t;
using u32 = uint32_t;
using uint = u32;
using i32 = int32_t;
using u16 = uint16_t;
using i16 = int16_t;
using u8 = uint8_t;
using i8 = int8_t;
using ld = long double;
using II = pair<int, int>;
using DD = pair<double, double>;
using SS = pair<string, string>;
using VS = vector<string>;
using VD = vector<double>;
using VI = vector<int>;
using VB = vector<bool>;
using NLD = numeric_limits<double>;
using NLI = numeric_limits<int>;
using ILI = initializer_list<int>;
using VII = vector<II>;
using VVI = vector<VI>;
using VVB = vector<VB>;
using VVVI = vector<VVI>;
using III = tuple<int, int, int>;
using VIII = vector<III>;
using VVII = vector<VII>;
using VVVII = vector<VVII>;
using USI = unordered_set<int>;
#define mt make_tuple
#define mp make_pair
#define COMBINE1(X,Y) X ##Y
#define COMBINE(X,Y) COMBINE1(X, Y)
inline int lastExitCode = EXIT_SUCCESS;
#define exit(x) \
  (exit)(lastExitCode = (x))
struct ff {
  static constexpr double EPS = 1e-7;
  static bool fEq(double a, double b) { return abs(a - b) < EPS; }
  ff() = default;
  ff(double v) : v(v) {}
  bool operator==(const ff &o) const { return fEq(v, o.v); }
  bool operator!=(const ff &o) const { return not fEq(v, o.v); }
  bool operator<=(const ff &o) const { return v < o.v or fEq(v, o.v); }
  bool operator>=(const ff &o) const { return v > o.v or fEq(v, o.v); }
  bool operator<(const ff &o) const { return v + EPS < o.v; }
  bool operator>(const ff &o) const { return o < *this; }
  operator double() const { return v; }
  static void fixZero(double &v) {
    if (ff(v) == ff(0.0))
      v = 0.0;
  }
  double v;
};
inline auto ffuple(double f) { return tuple(ff(f)); }
template <typename... Ts> auto ffuple(double f, Ts... v) {
  static_assert((is_same<Ts, double>::value && ...));
  return tuple_cat(tuple(ff(f)), ffuple(v...));
}
inline string exec(const char *cmd) {
  array<char, 128> buffer;
  string result;
  shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
  if (!pipe)
    throw std::runtime_error("popen() failed!");
  while (!feof(pipe.get())) {
    if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
      result += buffer.data();
  }
  return result;
}
inline string exec(const string &cmd) { return exec(cmd.c_str()); }
using fmt::format;
using fmt::print;
template <typename T> class IsStreamable {
  template <typename U> static auto t(const U *u) -> decltype(std::cout << *u);
  static auto t(...) -> std::false_type;
public:
  enum {
    value = !std::is_same<decltype(t((T *)nullptr)), std::false_type>::value
  };
};
template <typename T>
typename enable_if<!is_void<decltype(begin(T()))>::value &&
                       !IsStreamable<T>::value,
                   ostream &>::type
operator<<(ostream &o, const T &v) {
  for (auto it = begin(v); it != end(v); ++it) {
    if (next(it) != end(v))
      o << *it << " ";
    else
      o << *it;
  }
  return o;
}
template <typename T1, typename T2>
ostream &operator<<(ostream &o, const std::pair<T1, T2> &p) {
  o << "(" << p.first << ", " << p.second << ")";
  return o;
}
namespace aux {
template <std::size_t...> struct seq {};
template <std::size_t N, std::size_t... Is>
struct gen_seq : gen_seq<N - 1, N - 1, Is...> {};
template <std::size_t... Is> struct gen_seq<0, Is...> : seq<Is...> {};
template <class Ch, class Tr, class Tuple, std::size_t... Is>
void print_tuple(std::basic_ostream<Ch, Tr> &os, Tuple const &t, seq<Is...>) {
  using swallow = int[];
  (void)swallow{0,
                (void(os << (Is == 0 ? "" : ", ") << std::get<Is>(t)), 0)...};
}
}
template <class Ch, class Tr, class... Args>
auto operator<<(std::basic_ostream<Ch, Tr> &os, std::tuple<Args...> const &t)
    -> std::basic_ostream<Ch, Tr> & {
  os << "(";
  aux::print_tuple(os, t, aux::gen_seq<sizeof...(Args)>());
  return os << ")";
}
struct Timer {
  using Clock = chrono::steady_clock;
  using TimePoint = Clock::time_point;
  TimePoint tpStart;
  double tmLim;
  Timer(double timeLimitSecs, const Timer &parent) {
    reset(min(timeLimitSecs, parent.secsLeft()));
  }
  Timer(double timeLimitSecs = NLD::max()) { reset(timeLimitSecs); }
  void reset(double timeLimitSecs = NLD::max()) {
    tmLim = timeLimitSecs, tpStart = Clock::now();
  }
  double elapsedSecs() const {
    return chrono::duration_cast<chrono::duration<double>>(Clock::now() -
                                                           tpStart)
        .count();
  }
  double secsLeft() const { return tmLim - elapsedSecs(); }
  bool timedOut() const { return elapsedSecs() >= tmLim; }
};
#define USE_TIMED_BLOCKS 
#ifdef NDEBUG
#undef USE_TIMED_BLOCKS
#endif
#ifdef USE_TIMED_BLOCKS
inline unordered_map<string, double> timedBlocks;
#endif
struct TimedBlock {
#ifdef USE_TIMED_BLOCKS
  TimedBlock(const string &name) : name(name) {}
  ~TimedBlock() { timedBlocks[name] += timer.elapsedSecs(); }
  string name;
  Timer timer;
#else
  TimedBlock(const string &) {}
#endif
};
#define TIME_BLOCK(s) TimedBlock COMBINE(tbDummy, __LINE__)(s)
template <typename T> T vmin(T &&t) { return std::forward<T>(t); }
template <typename T0, typename T1, typename... Ts>
typename std::common_type<T0, T1, Ts...>::type vmin(T0 &&val1, T1 &&val2,
                                                    Ts &&... vs) {
  return val2 < val1 ? vmin(val2, std::forward<Ts>(vs)...)
                     : vmin(val1, std::forward<Ts>(vs)...);
}
double ternarySearchReal(double lo, double hi, const auto &f, auto cmp) {
  while (ff(lo) < ff(hi)) {
    double leftThird = lo + (hi - lo) / 3.0;
    double rightThird = hi - (hi - lo) / 3.0;
    if (!cmp(f(leftThird), f(rightThird)))
      lo = leftThird;
    else
      hi = rightThird;
  }
  return (lo + hi) / 2.0;
}
double binarySearchReal(double lo, double hi, double target, const auto &f,
                        auto cmp) {
  while (ff(lo) < ff(hi)) {
    double mid = (lo + hi) / 2.0, fmid = f(mid);
    if (ff(fmid) == ff(target))
      return mid;
    else if (cmp(ff(fmid), ff(target)))
      lo = mid + ff::EPS;
    else
      hi = mid - ff::EPS;
  }
  return lo;
}
inline string valOrNA(bool yes, double val) {
  return yes ? format("{}", val) : "NA";
}
template <typename ContainerType>
void insertToSortedVector(ContainerType &v,
                          const typename ContainerType::value_type &t,
                          bool allowDuplicates) {
  auto ub = upper_bound(begin(v), end(v), t);
  if (allowDuplicates or ub == begin(v) or *(ub - 1) != t)
    v.insert(ub, t);
}
template <typename ContainerType>
void eraseFromSortedVector(ContainerType &v,
                           const typename ContainerType::value_type &t) {
  auto it = lower_bound(begin(v), end(v), t);
  if (it != end(v) and *it == t)
    v.erase(it);
}
template <typename T> auto &noConst(const T &s) {
  return const_cast<remove_const_t<remove_reference_t<decltype(s)>> &>(s);
}
inline void removeDuplicates(auto &v) {
  sort(begin(v), end(v));
  v.erase(unique(v.begin(), v.end()), v.end());
}
inline bool isUnique(auto v) {
  sort(begin(v), end(v));
  return adjacent_find(begin(v), end(v)) == end(v);
}
inline auto sorted(auto r) {
  sort(begin(r), end(r));
  return r;
}
inline VI iotaed(size_t n, int start = 0) {
  VI ind(n);
  iota(begin(ind), end(ind), start);
  return ind;
}
template <typename ContainerType, typename ValueType>
size_t index(const ContainerType &c, const ValueType &v) {
  return find(begin(c), end(c), v) - begin(c);
}
template <typename ContainerType, typename ValueType>
bool linearIn(const ContainerType &c, const ValueType &v) {
  return index(c, v) != size(c);
}
template <typename ContainerType, typename ValueType>
bool binaryIn(const ContainerType &c, const ValueType &v) {
  return binary_search(begin(c), end(c), v);
}
template <typename T> void vecErase(vector<T> &v, const T &elem) {
  auto it = find(begin(v), end(v), elem);
  if (it != end(v))
    v.erase(it);
}
inline double interp(double val, double min1, double max1, double min2,
                     double max2) {
  return min2 + ((max2 - min2) * (val - min1)) / double(max1 - min1);
}
inline bool fileExists(const std::string &name) {
  struct stat buffer;
  return (stat(name.c_str(), &buffer) == 0);
}
template <typename V> inline bool inrange(V v, V lo, V hi) {
  return v >= lo and v <= hi;
}
struct TabuList {
  TabuList(int size, int tenure)
      : tabu(size, -tenure), iter(1), ten(tenure), sz(size) {}
  bool isTabu(int u) const { return not tabu.empty() and tabu[u] + ten > iter; }
  void add(int u) { tabu[u] = iter; }
  void remove(int u) { tabu[u] = -ten; }
  void advanceIter() { ++iter; }
  int since(int u) const { return isTabu(u) ? tabu[u] : NLI::max(); }
  void reset() { tabu.assign(sz, -ten), iter = 1; }
  VI tabu;
  int iter, ten, sz;
};
namespace std {
template <typename T1, typename T2> struct hash<pair<T1, T2>> {
  size_t operator()(const pair<T1, T2> &p) const {
    size_t seed(0);
    boost::hash_combine(seed, p.first);
    boost::hash_combine(seed, p.second);
    return seed;
  }
};
}
struct TimeoutException {
  Timer t;
  TimeoutException(Timer timer) : t(timer) {}
  const char *what() const {
    return format("Timeout exception (elapsed: {}, time limit: {})",
                  t.elapsedSecs(), t.tmLim)
        .c_str();
  }
};
