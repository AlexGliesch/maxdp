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
#include "util.h"
using VD = vector<double>;
int n, m;
double bbeta;
string type;
size_t rndSeed;
string outputFilename, suffix;
string outputFmt;
VD gw;
VD ow;
VD obx, oby;
vector<array<int, 25>>
    oblik;
void parseCommandLine(int argc, char** argv);
void generateDistances(const string& type);
void generateObjectWeights(VD& ow);
void generateGroupTargetWeights(VD& gw, double A);
void writeInstanceAlex(const VD& ow, const VD gw, ostream& o);
void writeInstanceMarcus(const VD& ow, const VD gw, ostream& o);
void parseCommandLine(int argc, char** argv) {
  namespace po = boost::program_options;
  po::options_description desc("MaxDP instance generator.");
  desc.add_options()("help", "");
  desc.add_options()("n", po::value<int>(&n)->required(), "Number of objects");
  desc.add_options()("m", po::value<int>(&m)->default_value(-1),
                     "Number of groups; default: floor(5+3n/200)");
  desc.add_options()("type", po::value<string>(&type)->default_value("weee"),
                     "Instance type, in [weee,study].");
  desc.add_options()("beta", po::value<double>(&bbeta)->default_value(0.25),
                     "(see Fernandez et al. (2013) for details on this "
                     "parameter).");
  desc.add_options()("out",
                     po::value<string>(&outputFilename)->default_value(""),
                     "Output filename. If unset, a default name will be used");
  desc.add_options()("suffix", po::value<string>(&suffix)->default_value(""),
                     "Suffix to filename. If \"seed\", the suffix will be "
                     "the seed number.");
  desc.add_options()("format",
                     po::value<string>(&outputFmt)->default_value("marcus"),
                     "Output format, in {alex,marcus}.");
  desc.add_options()("seed", po::value<size_t>(&rndSeed)->default_value(0),
                     "Random seed. If =0, a random value will be used.");
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
      cout << desc << endl;
      exit(EXIT_SUCCESS);
    } else {
      po::notify(vm);
    }
    boost::to_lower(type);
    boost::to_lower(outputFmt);
    size_t seed = rndSeed;
    if (seed == 0) seed = uniqueRandomSeed();
    seed ^= n;
    seed ^= m;
    seed ^= hash<string>()(type);
    rng.seed(seed);
    if (m == -1) m = 5 + (3 * n) / 200.0;
    if (type != "weee" and type != "study") {
      fmt::print("ERROR: invalid type \"{}\". See --help.\n", type);
      exit(EXIT_FAILURE);
    }
    if (outputFmt != "alex" and outputFmt != "marcus") {
      fmt::print("ERROR: invalid instance format \"{}\". Possible values are: "
                 "[alex,marcus].\n",
                 outputFmt);
      exit(EXIT_FAILURE);
    }
    if (type == "study" and outputFmt == "marcus") {
      fmt::print(
          "Warning: Marcus' instance format is not suitable for \'study\' "
          "type instances. Using Alex's format instead.\n");
      outputFmt = "alex";
    }
  } catch (po::error& e) {
    fmt::print("ERROR parsing command line: {}.\n", e.what());
    cout << desc << endl;
    exit(EXIT_FAILURE);
  }
}
void generateDistances(const string& type) {
  if (type == "weee") {
    assert(size(obx) + size(oby) == 0);
    for (int i = 0; i < n; ++i) {
      obx.push_back(randDouble(0, 10));
      oby.push_back(randDouble(0, 10));
    }
  } else if (type == "study") {
    oblik.resize(n);
    for (auto& i : oblik)
      for (auto& j : i)
        j = randInt(0, 4);
  }
}
void generateObjectWeights(VD& ow) {
  ow.resize(n);
  for (auto& w : ow)
    w = randDouble(1000.0, 4000.0);
}
void generateGroupTargetWeights(VD& gw, double A) {
  gw.resize(m);
  for (auto& w : gw)
    w = (A / (double)m) * randDouble(1.0 - bbeta, 1.0 + bbeta);
  double A2 = accumulate(begin(gw), end(gw), 0.0);
  for (auto& w : gw)
    w = (w * A) / A2;
  assert(ff(accumulate(begin(gw), end(gw), 0.0)) == ff(A));
}
void writeInstanceAlex(const VD& ow, const VD gw, ostream& o) {
  o << fmt::format("{} {}\n{} {} {}\n", n, m, type, rndSeed, bbeta);
  o << gw << "\n";
  o << ow << "\n";
  if (type == "weee") {
    for (int i = 0; i < n; ++i)
      o << obx[i] << " " << oby[i] << endl;
  } else if (type == "study") {
    for (auto& i : oblik)
      o << i << endl;
  }
}
void writeInstanceMarcus(const VD& ow, const VD gw, ostream& o) {
  if (type != "weee") throw false;
  o << n << " " << m << " " << bbeta << "\n";
  for (int i = 0; i < n; ++i)
    o << ow[i] << " " << obx[i] << " " << oby[i] << "\n";
  for (int i = 0; i < m; ++i)
    o << gw[i] << "\n";
}
int main(int argc, char** argv) {
  parseCommandLine(argc, argv);
  generateDistances(type);
  generateObjectWeights(ow);
  generateGroupTargetWeights(gw, accumulate(begin(ow), end(ow), 0));
  if (empty(outputFilename)) {
    outputFilename = fmt::format("{}-{}-{}-{:03}", type, n, m, bbeta * 100);
  }
  if (size(suffix))
    outputFilename += "-" + (suffix == "seed" ? to_string(rndSeed) : suffix);
  ofstream f(outputFilename);
  if (not f) {
    fmt::print("ERROR: could not open output file.\n");
    exit(EXIT_FAILURE);
  }
  if (outputFmt == "alex") {
    writeInstanceAlex(ow, gw, f);
  } else if (outputFmt == "marcus") {
    writeInstanceMarcus(ow, gw, f);
  }
}
