#!/bin/bash
#for文 betaは任意の変数一般的に「beta」「state_sigma」「cost」「n」なとが使われることが多い
#inより後ろを順に代入してゆくことになる。
export CPATH=:~/include/

for delta in 0.3 0.4 0.5
do
for genes in 4 5 6
do
# python L_noise_write.py $genes -1 1 $delta -1 1000
for num_noise_type in 1000
do
for max_fitness in 10 15 20 
do
mkdir -p ../data_random_noise${delta}/genes_${genes}_MC_step_fitness=pre/num_noise_type_${num_noise_type}/max=${max_fitness}/
# mv ../data_random_noise${delta}/genes_${genes}_MC_step_fitness=pro_ignoring_lc/ ../data_random_noise${delta}/NG_genes_${genes}_MC_step_fitness=pro_ignoring_lc/
bsub -e ../Error/error_MC_pro_${delta}_${num_noise_type}_${max_fitness}.out -q "normalCPU" "
./MC_step_fitness=pro_debugged $genes $num_noise_type $max_fitness $delta &&
"
done
done
done
done