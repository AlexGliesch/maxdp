#!/usr/bin/env julia
################################################################################
# mdp: An implementation of the method of Fernandez et al. (2013) for the 
# maximum dispersion problem.
#
# The overall idea is the following:
#
# 1. First solve the unrestricted maxDP (i.e. ignore balance constraints), to
# obtain an upper bound u^*.
# 2. Next solve a sequence of coloring problems to obtain the chromatic 
# lower bound LC^C.
# 3. Then, repeatedly solve the maxDP (in the covering formulation) using the 
# bounds above to restrict the variables.
#
#
# Euclidean instances are represented as the number of objects and groups
#   n m
# followed by a list of objects
#   weight x y
# and a list of groups
#  demands
#
################################################################################
using ArgParse
using JuMP
using CPLEX
using DataStructures

##### configuration
# for CPLEX
global timelimit = 600
global instanceName
global threads = 1

include("instance.jl")

## build the umaxDPcov model
function buildUmaxDPcov(I::Instance, lb::Float64, ub::Float64, verbose::Bool=false)
    m = Model(solver=CplexSolver(CPX_PARAM_TILIM=timelimit-(time()-start),
                                 CPX_PARAM_THREADS=threads,
                                 CPX_PARAM_EPGAP=1e-6, 
                                 CPX_PARAM_SCRIND=(verbose ? 1 : 0)))
    R = size(I.R,1)
    @variable(m, w[1:R], Bin) 
    @variable(m, x[1:I.n,1:I.m], Bin)
    @objective(m, Max, I.R[R][3] + sum((I.R[r][3] - I.R[r+1][3]) * w[r] for r=1:R-1))
    @constraint(m, [i=1:I.n], sum(x[i,k] for k=1:I.m)==1) ## assignment constraint (2)
    @constraint(m, [k=1:I.m,r=1:R], x[I.R[r][1],k]+x[I.R[r][2],k] <= 1+w[r]) ## exclusion constraint (7)
    @constraint(m, [r=2:R], w[r-1]<=w[r]) ## weight constraint (8)
    @constraint(m, [r=1:R; I.R[r][3]<lb], w[r]==0) ## (simple) lower and upper bound
    @constraint(m, [r=1:R; I.R[r][3]>ub], w[r]==1) # may be >=

    m,w,x
end

## solve the covering formulation of the maxDP
function buildmaxDPcov(I::Instance, α::Float64, lb::Float64, ub::Float64)
    m,w,x = buildUmaxDPcov(I,lb,ub)

    ## balance lower and upper bounds (3) and (4)
    @constraint(m, [k=1:I.m], (1-α)*I.M[k]<=sum(I.a[i]*x[i,k] for i=1:I.n))
    @constraint(m, [k=1:I.m], sum(I.a[i]*x[i,k] for i=1:I.n)<=(1+α)*I.M[k])
    m,w,x
end

## solve the unrestricted maxDP
function umaxDPcov(I::Instance, lb::Float64, ub::Float64)
    m = buildUmaxDPcov(I,lb,ub)[1]
    status = solve(m)
    status == :Optimal && return getobjectivevalue(m)
    return -1
    error("Unrestricted problem umaxDP could not be solved.")
end

## helper for collecting statistics on maxDP
mutable struct dpInfo
  nsolved :: Int # number of solved integer programs
  lb :: Float64 # lower bound
  ub :: Float64 # upper bound
  lbColor :: Float64 # lower bound obtained by repeated coloring
  lbUMaxDPCov :: Float64 # lower bound obtained by umaxDPcov, if lbColor != ubi (otherwise = lbColor)
  ubI :: Float64 # UB^I, in the paper
  memlim :: Int
  timelim :: Int
end

## solve the maxDP, return a coloring
function maxDP(I::Instance, ub::Int, lb::Int, α::Float64, verbose::Bool, justUmaxDP::Bool)
  stat = dpInfo(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0) 
  ## (1) compute UB^I and LB^C
  ubReal = I.rUnique[ub]
  println("Upper bound $ub (r$ubReal).")
  lbReal = I.rUnique[lb]
  println("Lower bound $lb (r$lbReal).")  
  stat.ub = stat.ubI = ubReal
  stat.lb = stat.lbColor = stat.lbUMaxDPCov = lbReal
  r = lb
  ## (2) find u^*
  if time()-start > timelimit
    return zeros(I.n),-1,stat
  end
  if lb != ub
    verbose && println("lb != ub, solving umaxDPcov...")
    # ubReal = 0
    try 
      ubReal = umaxDPcov(I,lbReal,ubReal)
    catch ex
      if isa(ex, OutOfMemoryError)
        println("Caught OutOfMemoryError!!!")  
        stat.memlim=1
        rethrow(ex)
      end
      println("Exception: $ex")
      bt = catch_backtrace()
      msg = sprint(showerror, ex, bt)
      println("Message: $msg")
      if contains(msg, "memory") || contains(msg, "OutOfMemory")
        stat.memlim=1
      end
      return zeros(I.n),-1,stat
    end 
    ubV = findall(x->abs(x-ubReal) < 10e-5, I.rUnique)
    if length(ubV)==0
      verbose && println("length(ubV)==0; ubReal: $(ubReal), rUnique: $(I.rUnique)")
      return zeros(I.n),-1,stat
    end
    ub = ubV[1]
    stat.lbUMaxDPCov = ubReal
    verbose && println("Upper bound UMaxDPCov is $ub (r$ubReal).")
  end
  justUmaxDP && return zeros(I.n),ub,stat

  ## (3) repeated solve maxDP
  println("(3) repeatedly solve MaxDP...")
  m,w,x = buildmaxDPcov(I,α,lbReal,ubReal)
  foundopt = false
  while r >= 0 && time()-start < timelimit
    verbose && println("Solving maxDP in [$lb,$ub]...")
    status = 0
    try
      status = solve(m)
    catch ex
      if isa(ex, OutOfMemoryError)
        println("Caught OutOfMemoryError!")  
        stat.memlim=1
        rethrow(ex)        
      end
      println("Exception: $ex")
      bt = catch_backtrace()
      msg = sprint(showerror, ex, bt)
      println("Message: $msg")
      if contains(msg, "memory") || contains(msg, "OutOfMemory")
        stat.memlim=1
      end
      return zeros(I.n),-1,stat
    end 
    stat.nsolved += 1
    if status == :Optimal
      foundopt=true
      break
    end
    ## not optimal: next round
    r -= 1
    if r == 0
      println("Infeasible.")
      break
    end
    lb = ub = r 
    lbReal = ubReal = I.rUnique[r]
    stat.lb = stat.ub = lbReal
    m,w,x = buildmaxDPcov(I,α,lbReal,ubReal)
  end
  ## (4) process the solution
  C=zeros(I.n)
  if foundopt
    verbose && println("Optimal solution of value $(getobjectivevalue(m)) found.")
    for i=1:I.n, k=1:I.m
      if getvalue(x[i,k])>0.5
          C[i]=k
      end
    end
  end
  C, foundopt ? getobjectivevalue(m) : -1, stat
end

# (0) check the input
s = ArgParseSettings()
@add_arg_table s begin
  "instance"
    help = "MDP instance"
    required = true
  "--alpha" 
    arg_type = Float64
    help = "Maximum relative deviation from target weight"
    default = 0.05
  "--timelimit"
    help = "time limit for the solver"
    arg_type = Int
    default = 600
  "--threads"
    help = "number of threads for solving mathematical programs"
    arg_type = Int
    default = 1
  "--ub"
    help = "upper bound"
    arg_type = Float64
    required = true
  "--lb"
    help = "lower bound"
    arg_type = Float64
    required = true
  "--wdot"
    help="write dot graph output file"
  "--verbose"
    help="verbose"
end

args = parse_args(s)

# verbose=args["verbose"]
verbose=true
timelimit=args["timelimit"]
threads=args["threads"]
justUmaxDP=false
instanceName = args["instance"]

# (1) read instance
inst=open(args["instance"])
I=readInstance(inst)
close(inst)

ubIntV = findall(x -> abs(x-args["ub"]) < 10.0^(-5), I.rUnique)
if length(ubIntV) == 0 || (args["ub"] == 0)
  println("Invalid upper bound, should be within $(I.rUnique[2]) and $(I.rUnique[length(I.rUnique)]).")
  exit()
end 
ubInt = ubIntV[1]
ubReal = I.rUnique[ubInt]

lbIntV = findall(x -> abs(x-args["lb"]) < 10.0^(-5), I.rUnique)
if length(lbIntV) == 0 || (args["lb"] == 0)
  println("Invalid lower bound, should be within $(I.rUnique[2]) and $(I.rUnique[length(I.rUnique)]).")
  exit()
end 
lbInt = lbIntV[1]
lbReal = I.rUnique[lbInt]

verbose && println("alpha: $(args["alpha"])")

# (2) solve the problem
start=time()
C,v,stat=zeros(I.n),-1,dpInfo(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0) 
try
  C,v,stat = maxDP(I, ubInt, lbInt, args["alpha"], verbose, justUmaxDP)
catch ex
  if isa(ex, OutOfMemoryError)
    println("OOM error reached main function")
  end
  println("summary_line instance=$(basename(args["instance"])) n=$(I.n) m=$(I.m) beta=$(I.β) time=$(time()-start) alpha=$(args["alpha"]) nsolved=0 ubi=$(ubInt) lbColor=$(lbInt) ub=$(ubInt) lb=$(lbInt) vs=NA timelim=0 memlim=1")
  exit()
end
vS = (v == -1 ? "NA" : "$v")
# instance n m beta time alpha nsolved ubI lb ub vS
stat.timelim = Int(time()-start > timelimit)
println("summary_line instance=$(basename(args["instance"])) n=$(I.n) m=$(I.m) beta=$(I.β) time=$(time()-start) alpha=$(args["alpha"]) nsolved=$(stat.nsolved) ubi=$(stat.ubI) lbColor=$(stat.lbColor) ub=$(stat.ub) lb=$(stat.lb) vs=$vS timelim=$(stat.timelim) memlim=$(stat.memlim)")

# (3) optionally, provide some output
if args["wdot"] != nothing
  dots=open(args["wdot"],"w")
  println(dots,"graph {")
  println(dots,"node [colorscheme=\"paired12\",style=filled,shape=circle,margin=\"0,0\"]")
  for i=1:I.n
      println(dots,"v$i [pos=\"$(250*I.x[i]),$(250*I.y[i])\",fillcolor=$(C[i]),label=\"$i\"]")
  end
  println(dots,"}")
  close(dots)
end