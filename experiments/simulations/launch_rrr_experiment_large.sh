#!/bin/bash

# Define the values for the variables
theta_strengths="high medium low"
#theta_strengths="medium"
n_values="500"
#p_values="100 300 500 800"
p_values="1000 2000 3000"
r_values="3"
#q_values="80"
q_values="5 10 30 50 80"
for theta in $theta_strengths; do
  for n in $n_values; do
    for p in $p_values; do
      for r in $r_values; do
	 for q in $q_values; do
      sbatch rrr_experiment_large.sh "$n" "$theta" "$p" "$r" "$q"
    #echo "$theta"
    done
      done
    done
  done
done     
