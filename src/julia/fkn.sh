#!/usr/bin/env bash
# $1: instance
# $2: time limit, in seconds (default: 3600 seconds)
# $3: alpha (default: 0.001)

TIMELIMIT=${2:-3600}
ALPHA=${3:-0.001}
UB="$(./../maxdp --test ub --ub ubi --irace --in $1)"
UMDP="$(./../maxdp --test umdp --ub ubi --irace --in $1 --time $TIMELIMIT --umdpalg fer)"
UMDP="${UMDP//-/}" # remove negative sign 
echo "UB = $UB"
echo "UMDP = $UMDP"
julia mdp.jl $1 --timelimit $TIMELIMIT --alpha $ALPHA --ub $UB --lb $UMDP