#!/bin/bash
#for文 betaは任意の変数一般的に「beta」「state_sigma」「cost」「n」なとが使われることが多い
#inより後ろを順に代入してゆくことになる。
export CPATH=:~/include/

for genes in 5 #3 4 5 6 
do
for delta in   0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
for seed in `seq 58 62` #`seq 53 57` #`seq 51 52` 
do
for L in 1000
do
for epsilon in 0.01
do
for num_noise_type in 1000
do
# mkdir -p ../data_random_noise${delta}/genes_${genes}_step/epsilon_${epsilon}/num_noise_type_${num_noise_type}_without_e/network_info_step/seed_${seed}/GA_ignoring_lc
# mkdir -p ../data_random_noise${delta}/genes_${genes}_step/epsilon_${epsilon}/num_noise_type_${num_noise_type}/network_info_step/seed_${seed}/GA_ignoring_lc_10^6gen
mkdir -p ../data_random_noise${delta}/genes_${genes}_step/epsilon_${epsilon}/num_noise_type_${num_noise_type}/network_info_step/seed_${seed}/GA_ignoring_lc_10^7gen
bsub -e ../Error/error_GAstep10^7gen_${genes}_${num_noise_type}_${seed}_${delta}_${epsilon}_${L}.out  -q "jobLimit60" "
# ./GA_step $genes $num_noise_type $seed $delta $epsilon $L&&
python L_noise_write.py $genes -1 $seed $delta -1 $L &&
# ./GA_step_ignoring_lc2_without_e $genes $num_noise_type $seed $delta $epsilon $L&&
# ./GA_step_ignoring_lc2_including_e $genes $num_noise_type $seed $delta $epsilon $L&&
./GA_step_ignoring_lc2_including_e_10^7 $genes $num_noise_type $seed $delta $epsilon $L&&
"
done
done
done
done
done
done