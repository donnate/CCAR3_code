#!/bin/sh

set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
batch_script="$script_dir/rrr_experiment.sh"

mkdir -p "$script_dir/logs"

# Define the values for the variables
theta_strengths="high medium low"
n_values="500"
#p_values="100 300 500 800"
p_values="50 100 300 500 750 1000 2500 5000 10000"
r_values="3"
q_values="10 20"

for theta in $theta_strengths; do
  for n in $n_values; do
    for p in $p_values; do
      for r in $r_values; do
        for q in $q_values; do
          sbatch --chdir="$script_dir" "$batch_script" "$n" "$theta" "$p" "$r" "$q"
        done
      done
    done
  done
done
