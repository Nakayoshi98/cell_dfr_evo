#!/bin/bash
#for文 betaは任意の変数一般的に「beta」「state_sigma」「cost」「n」なとが使われることが多い
#inより後ろを順に代入してゆくことになる。
export CPATH=:~/include/


for delta in 0.4 0.5 #0.3
do
for genes in 4 5 6 
do
for num_noise_type in 1000
do
for max_fitness in 10 15 20
do
for seed in `seq 2 10`
do
mkdir -p ../data_random_noise${delta}/genes_${genes}_MC_step_fitness=pro_ignoring_lc/num_noise_type_${num_noise_type}/max=${max_fitness}/samples/seed_${seed}
bsub -e ../Error/error_MC_${delta}_${num_noise_type}_${max_fitness}_${seed}.out -q "normalCPU" "
./sampling_MC_step_fitness=pro $genes $num_noise_type $max_fitness $delta $seed &&
"
done
done
done
done
done