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
#include "cmdline.h"
#include "2ge.h"
#include "bal.h"
#include "balvns.h"
#include "constructive.h"
#include "ec.h"
#include "oscillate.h"
#include "readinstance.h"
#include "ub.h"
#include "umdp.h"
#include "umdpfernandez.h"
#include "util.h"
namespace {
string ecnodeorderS, ecgrouporderS, balNodeSelStratS;
void buildOptionsDesc(po::options_description& desc) {
  desc.add_options()("help,h", "Show help menu");
  desc.add_options()("in,i", po::value<string>(&inputFilename)->required(),
                     "Input filename");
  desc.add_options()("out,o", po::value<string>(&outputFilename),
                     "Output filename, to which solution will be printed. If "
                     "unset, no output file will be produced.");
  desc.add_options()(
      "time,t",
      po::value<double>(&timeLimit)->default_value(double(NLI::max())),
      "Time limit, in seconds");
  desc.add_options()("seed,s", po::value<size_t>(&rndSeed)->default_value(0),
                     "Random seed. If 0, a random value will be used");
  desc.add_options()("test",
                     po::value<string>(&testType)->default_value("full"),
                     "Test type, in [full,ub,umdp,bal].");
  desc.add_options()(
      "format", po::value<string>(&instFmt)->default_value( "alex"),
      "Instance format, in [alex,marcus].");
  desc.add_options()("irace", "Output a single value, for irace tests.");
  desc.add_options()("tabucol",
                     po::value<double>(&tabuColDValue)->default_value(-1.0),
                     "Output a single value, which is the solution of tabuCol "
                     "with the given d-value.");
  desc.add_options()("cons",
                     po::value<string>(&consAlg)->default_value("greedy"),
                     "Constructive heuristic to generate intial solution, in "
                     "[random,greedy,trivial]");
  desc.add_options()(
      "ub", po::value<string>(&ubAlg)->default_value("auto"),
      "(Upper bound to use, in "
      "[ubi,ubs,ubsblind,ubrb,ubk,all,none,auto]. Default (auto): "
      "ubs for weee instances, none for study instances. If "
      "option \'all\' is selected, will use min(ubi,ubs,ubrb). "
      "If \"none\" is selected, d^R will be used as upper "
      "bound. If \"ubsblind\" is selected, will compute ub^S "
      "starting from d^R instead of ub^I.");
  desc.add_options()(
      "ubtime",
      po::value<double>(&ubTimeLimit)->default_value(double(NLI::max())),
      "Time limit for computing the upper bound.");
  desc.add_options()(
      "ubssigma", po::value<double>(&ubsSigma)->default_value(1.5),
      "Size of subsets (relative to m) to be considered by ub^S.");
  desc.add_options()(
      "sigmafixed",
      "If ub=ubs uses ubsssigma as a fixed value, rather than relative to m.");
  desc.add_options()("ubrbls", "If ub=ubrb, if set will perform a local search "
                               "on the subsets selected by U^RB_h. ");
  desc.add_options()("ubsall", "If ub=ubs, will consider all n subsets, "
                               "instead of only the first 2n\\sigma subsets.");
  desc.add_options()("umdpalg",
                     po::value<string>(&umdpAlg)->default_value("ec"),
                     "Algorithm to use for UMDP, in [ec,fer]");
  desc.add_options()(
      "umdptime",
      po::value<double>(&umdpTimeLimit)->default_value(double(NLI::max())),
      "Time limit for UMDP.");
  desc.add_options()("umdprepl", po::value<int>(&umdpRepl)->default_value(1),
                     "Maximum number replications of UMDP. Only used if option "
                     "--cons is random.");
  desc.add_options()(
      "umdpferdir", po::value<string>(&umdpFernandezDir)->default_value("down"),
      "Optimization direction for Fernandez's UMDP, in [up,down].");
  desc.add_options()(
      "umdpfermaxitertabucol",
      po::value<int>(&umdpFernandezTabuColMaxIter)->default_value(20000),
      "Maximum TabuCol iterations for Fernandez's UMDP.");
  desc.add_options()(
      "umdpfersteps",
      po::value<string>(&umdpFernandezSteps)->default_value("linear"),
      "Search strategy for Fernandez's UMDP, in [linear,bs].");
  desc.add_options()("ecp1", po::value<int>(&ec::p1)->default_value(4),
                     "Lower depth (p1) for EC.");
  desc.add_options()("ecp2", po::value<int>(&ec::p2)->default_value(8),
                     "Middle depth (p2) for EC.");
  desc.add_options()("ecp3", po::value<int>(&ec::p3)->default_value(8),
                     "Upper depth (p3) for EC.");
  desc.add_options()("ecnodeorder",
                     po::value<string>(&ecnodeorderS)->default_value("hist"),
                     "EC group ordering, in [rand,conf,hist]");
  desc.add_options()("ecgrouporder",
                     po::value<string>(&ecgrouporderS)->default_value("conf"),
                     "EC group ordering, in [rand,conf,hist]");
  desc.add_options()("alpha", po::value<double>(&alpha)->default_value(0.001),
                     "Balance tolerance parameter \alpha.");
  desc.add_options()("balalg", po::value<string>(&balAlg)->default_value("2ge"),
                     "Algorithm to use for balancing solutions, in [2ge,vns]");
  desc.add_options()(
      "baltime", po::value<double>(&balTimeLimit)->default_value(NLI::max()),
      "Time limit for balancing algorithm.");
  desc.add_options()(
      "balnodesel",
      po::value<string>(&balNodeSelStratS)->default_value("value"),
      "Node selection strategy for balncing B&B, in [dfs,value,lb].");
  desc.add_options()("balmaxshakes",
                     po::value<int>(&balMaxShakes)->default_value(-1.0),
                     "Maximum number of shakes per call to balancing. Default: "
                     "35, or 50 for VNS.");
  desc.add_options()(
      "balshakeamt", po::value<double>(&balShakeAmtDbl)->default_value(-1.0),
      "Number of random moves to perform (relative to n) "
      "when shaking a balancing solution. Default: 0.01, or 0.03 for VNS.");
  desc.add_options()(
      "balshakethr", po::value<double>(&balShakeThreshold)->default_value(0.5),
      "Threshold for balancing solutions on BB: if a solution's imbalance is "
      "larger than this value, shakes will not be performed.");
  desc.add_options()("tau", po::value<double>(&twoge::tau)->default_value(1.0),
                     "Balancing tabu tenure, relative to m (i.e., the actual "
                     "tenure is \tau*m).");
  desc.add_options()(
      "th1", po::value<u64>(&twoge::th1)->default_value(pow(2, 11)),
      "Minimum (maximum) number of nodes for the balancing branch & bound.");
  desc.add_options()(
      "th2", po::value<u64>(&twoge::th2)->default_value(pow(2, 18)),
      "Maximum (maximum) number of nodes for the balancing branch & bound.");
  desc.add_options()("osceciter",
                     po::value<int>(&oscMaxEcIter)->default_value(1),
                     "Maximum EC iterations per oscillation steps.");
}
void readStringOpt(string optName, string& optValue, VS values) {
  boost::to_lower(optValue);
  if (not linearIn(values, optValue))
    throw logic_error(
        format("{} ({}) must be in [{}]\n", optName, optValue, values));
}
}
void parseCommandLine(int argc, char** argv) {
  po::options_description desc(
      "Heuristic solutions for the maximum dispersion problem");
  buildOptionsDesc(desc);
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
      cout << desc << endl;
      exit(EXIT_SUCCESS);
    } else {
      po::notify(vm);
    }
    iraceTest = vm.count("irace");
    if (rndSeed == 0) rndSeed = uniqueRandomSeed();
    rng.seed(rndSeed);
    tabuColTest = tabuColDValue > 0.0;
    ubrbDoLS = vm.count("ubrbls");
    ubsFewer = not vm.count("ubsall");
    sigmaFixed = vm.count("sigmafixed");
    readStringOpt(
        "ub", ubAlg,
        VS{"ubi", "ubs", "ubsblind", "ubrb", "ubk", "all", "none", "auto"});
    readStringOpt("format", instFmt, VS{"marcus", "alex"});
    readStringOpt("cons", consAlg, VS{"greedy", "random", "trivial"});
    readStringOpt("balalg", balAlg, VS{"2ge", "vns"});
    readStringOpt("umdpalg", umdpAlg, VS{"fer", "ec"});
    readStringOpt("test", testType, VS{"full", "ub", "umdp"});
    readStringOpt("umdpferdir", umdpFernandezDir, VS{"up", "down"});
    readStringOpt("umdpfersteps", umdpFernandezSteps, VS{"linear", "bs"});
    readStringOpt("balnodesel", balNodeSelStratS, VS{"dfs", "lb", "value"});
    readStringOpt("ecnodeorder", ecnodeorderS, VS{"rand", "conf", "hist"});
    readStringOpt("ecgrouporder", ecgrouporderS, VS{"rand", "conf", "hist"});
    if (floor(ubsSigma * m) > n)
      throw logic_error(
          "invalid value for ubssigma: must be within 1.0 and n/m.");
    map<string, ec::ExpansionOrder> ecOrderNames = {{"rand", ec::ordRandom},
                                                    {"conf", ec::ordConflicts},
                                                    {"hist", ec::ordHistory}};
    ec::nodeOrder = ecOrderNames[ecnodeorderS];
    ec::groupOrder = ecOrderNames[ecgrouporderS];
    map<string, twoge::NodeSelectionStrategy> nodeSelNames = {
        {"dfs", twoge::dfs}, {"lb", twoge::lbValue}, {"value", twoge::valueLb}};
    twoge::nodeSelectionStrategy = nodeSelNames[balNodeSelStratS];
    if (ec::p1 > ec::p2 or ec::p1 > ec::p3 or ec::p2 > ec::p3)
      throw logic_error(
          format("ecp1 ({}), ecp2 ({}), ecp3 ({}) should be such that "
                 "p1<=p2<=p3.\n",
                 ec::p1, ec::p2, ec::p3));
    if (testType != "umdp" and umdpAlg == "fer") {
      throw logic_error("\"fer\" umdp algorithm is only valid for umdp tests.");
    }
  } catch (po::error& e) {
    print("ERROR parsing command line: {}.\n", e.what());
    cout << desc << endl;
    exit(EXIT_SUCCESS);
  }
  pr("Input file: {}\n", inputFilename);
  pr("Instance format: {}\n", instFmt);
  pr("Seed: {}\n", rndSeed);
  pr("Time limit: {}\n", timeLimit);
}
