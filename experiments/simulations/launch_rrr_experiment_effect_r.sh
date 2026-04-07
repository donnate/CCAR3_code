#!/bin/bash

# Define the values for the variables
theta_strengths="high medium low"
n_values="500"
#p_values="80 100 300 500 800"
p_values="1000"
r_values="2 3 5 7 10"

for theta in $theta_strengths; do
  for n in $n_values; do
    for p in $p_values; do
      for r in $r_values; do
          sbatch rrr_experiment_effect_r.sh "$n" "$theta" "$p" "$r" "30"
    #echo "$theta"
      done
    done
  done
done
#1: n
#2: strength theta
#3: p
#4: r
#5: q    
