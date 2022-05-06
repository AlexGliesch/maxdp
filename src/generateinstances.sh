#!/usr/bin/env bash

# This script generates the test instances of our paper
# $1: num. processors (default: 1)

P=${1:-1} # Num. processors
mkdir inst

# Calibration set, seeds 11-20
parallel -j $P ./generateinstance --n {1} --type study --beta 0.1 --suffix {2} --seed {2} --format alex ::: `seq 400 400 4000` ::: `seq 11 1 20`
mkdir inst/calib
rm inst/calib/*
mv study-* inst/calib

parallel -j $P ./generateinstance --n {1} --type weee --beta 1.0 --suffix {2} --seed {2} --format alex ::: `seq 400 400 4000` ::: `seq 11 1 20`
mv weee-* inst/calib

# Test set, seeds 1-10
parallel -j $P ./generateinstance --n {1} --type weee --beta {2} --suffix {3} --seed {3}  --format alex ::: `seq 400 400 4000` ::: 0.25 0.5 0.75 1.0 ::: `seq 10` 
parallel -j $P ./generateinstance --n {1} --type study --beta {2} --suffix {3} --seed {3}  --format alex ::: `seq 400 400 4000` ::: 0.0 0.1 ::: `seq 10` 
mkdir inst/test
rm inst/test/*
mv weee-* study-* inst/test

# Instances with same sizes as Fernandez, seeds 1-10
parallel -j $P ./generateinstance --n 200 --m {1} --type weee --beta {2} --suffix {3} --seed {3}  --format alex ::: 4 6 8 ::: 0.25 0.5 0.75 1.0 ::: `seq 10` 
parallel -j $P ./generateinstance --n 300 --m {1} --type weee --beta {2} --suffix {3} --seed {3}  --format alex ::: 5 7 9 ::: 0.25 0.5 0.75 1.0 ::: `seq 10` 
parallel -j $P ./generateinstance --n 400 --m {1} --type weee --beta {2} --suffix {3} --seed {3}  --format alex ::: 7 9 11 ::: 0.25 0.5 0.75 1.0 ::: `seq 10` 
parallel -j $P ./generateinstance --n 500 --m {1} --type weee --beta {2} --suffix {3} --seed {3}  --format alex ::: 8 10 12 ::: 0.25 0.5 0.75 1.0 ::: `seq 10` 
parallel -j $P ./generateinstance --n 600 --m {1} --type weee --beta {2} --suffix {3} --seed {3}  --format alex ::: 10 12 14 ::: 0.25 0.5 0.75 1.0 ::: `seq 10` 
parallel -j $P ./generateinstance --n 700 --m {1} --type weee --beta {2} --suffix {3} --seed {3}  --format alex ::: 12 14 16 ::: 0.25 0.5 0.75 1.0 ::: `seq 10` 
parallel -j $P ./generateinstance --n 100 --m {1} --type study --beta {2} --suffix {3} --seed {3}  --format alex ::: 4 5 ::: 0.0 0.1 ::: `seq 10` 
parallel -j $P ./generateinstance --n 200 --m {1} --type study --beta {2} --suffix {3} --seed {3}  --format alex ::: 4 6 8 ::: 0.0 0.1 ::: `seq 10` 
parallel -j $P ./generateinstance --n 300 --m {1} --type study --beta {2} --suffix {3} --seed {3}  --format alex ::: 5 7 ::: 0.0 0.1 ::: `seq 10` 
mkdir inst/fer
rm inst/fer/*
mv weee-* study-* inst/fer

# Quick test instances 
parallel -j $P ./generateinstance --n {1} --type weee --beta 1.0 --suffix 1 --seed 99 --format alex ::: `seq 400 800 4000`
parallel -j $P ./generateinstance --n {1} --type study --beta 0.1 --suffix 1 --seed 99 --format alex ::: `seq 400 800 4000`
mkdir inst/quick
rm inst/quick/*
mv weee-* study-* inst/quick