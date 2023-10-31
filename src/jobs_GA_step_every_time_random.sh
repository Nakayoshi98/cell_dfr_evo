#!/bin/bash
#for文 betaは任意の変数一般的に「beta」「state_sigma」「cost」「n」なとが使われることが多い
#inより後ろを順に代入してゆくことになる。
export CPATH=:~/include/

for genes in 5
do
for delta in 0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
for seed in `seq 6 10`
do
for num_noise_type in 1000
do
mkdir -p ../data_random_noise${delta}/genes_${genes}_step/num_noise_type_${num_noise_type}_every_time_random/network_info_step/seed_${seed}/GA_ignoring_lc_10^6gen
./GA_step_ignoring_lc2_including_e_every_time_random $genes $num_noise_type $seed $delta &
done
done
done
done