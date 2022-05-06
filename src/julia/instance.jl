######################################################################
# Representation of a problem instance, and basic functions
######################################################################
using Distributions

mutable struct Instance
  n :: Int               # number of objects
  m :: Int               # number of groups
  a :: Array{Float64}    # weights of the objects
  M :: Array{Float64}    # target weights for the groups
  # x :: Array{Float64}    # coordinates
  # y :: Array{Float64}    #
  d :: Array{Float64,2}  # alternative: distances
  R :: Array{Tuple{Int,Int,Float64}}      # edge list, in increasing distance
  rUnique :: Array{Float64}             # edge list, unique values
  β :: Float64
  instanceType :: String         # type, in [study,weee]
end
    
function distance(I::Instance, u::Int, v::Int)
  I.d[u,v]
end

## Compute the maximum distance points in `S` (given by their indices
## in instance `I`). Return maximum distance, a pair of maximum
## distance, and the number of pairs of maximum distance.
function maxDistance(I::Instance, S::Array{Int,1})
  D = 0.0
  nD = 0
  Dp = (0,0)
  for i=1:length(S), j=i+1:length(S)
    d = distance(I,S[i],S[j])
    if d>D
      D = d
      nD = 1
      Dp = (i,j)
    elseif d==D
      nD += 1
    end
  end
  D,Dp,nD
end

## Generate a random WEEE instance.
function genWEEE(n,m,β)
  a=[rand(Uniform(1000,4000)) for i=1:n]
  x=[rand(Uniform(0,10)) for i=1:n]
  y=[rand(Uniform(0,10)) for i=1:n]
  A=sum(a)
  M=[rand(Uniform((1-β)*A/m,(1+β)*A/m)) for j=1:m]
  M=A*M/sum(M)
  R=Tuple{Int,Int,Float64}[]
  I=Instance(n,m,a,M,x,y,R,β)
  computeR(I)
  I
end

## Compute all pairwise distances, order by non-decreasing distance
function computeR(I::Instance)
  for i=1:I.n, j=i+1:I.n
    push!(I.R, (i,j,distance(I,i,j)))
    push!(I.rUnique, distance(I,i,j))
  end
  push!(I.rUnique, 0) # to be compatible with C++ src
  I.R=sort(I.R,lt=(a,b)->a[3]<b[3])
  I.rUnique=sort(unique(I.rUnique))
  # println("rUnique: $(I.rUnique)")
end

function readInstance(io::IO)
  # # (1) read size
  (n,m) = map(x->parse(Int,x), split(readline(io)))
  (instanceType,seeds,betas) = split(readline(io))
  β = parse(Float64, betas)
  # (2) read objects
  M = map(x->parse(Float64,x), split(readline(io)))
  a = map(x->parse(Float64,x), split(readline(io)))
  d = zeros(Float64, n, n)
  if instanceType == "weee"
    x = zeros(n)
    y = zeros(n)
    for i=1:n
      (x[i], y[i]) = map(x->parse(Float64, x), split(readline(io)))
    end    
    for u=1:n,v=1:n
      dx = x[u]-x[v]
      dy = y[u]-y[v]
      d[u,v] = sqrt(dx*dx+dy*dy)
    end
  elseif instanceType == "study"
    lik = zeros(Int, n, 25)
    for i=1:n
      vals = map(x->parse(Int,x),split(readline(io)))
      for j in 1:25
        lik[i,j] = vals[j]
      end      
    end
    for u=1:n, v=1:n, k=1:25
      d[u,v]+=abs(lik[u,k]-lik[v,k])
    end    
  end 
  # (4) create distances
  R=Tuple{Int,Int,Float64}[]
  rUnique = Float64[]
  I=Instance(n,m,a,M,d,R,rUnique,β,instanceType)
  computeR(I)
  I
end

function writeInstance(io::IO,I)
  println(io, "$(I.n) $(I.m) $(I.β)")
  for i = 1:I.n
    println(io,"$(I.a[i]) $(I.x[i]) $(I.y[i])")
  end
  #println("# A=$(sum(I.a)) A/m=$(sum(I.a)/I.m)")
  for j = 1:I.m
    println(io,"$(I.M[j])")
  end
end
