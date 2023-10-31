#!/bin/bash
#for文 betaは任意の変数一般的に「beta」「state_sigma」「cost」「n」なとが使われることが多い
#inより後ろを順に代入してゆくことになる。
export CPATH=:~/include/

# g++ -O3 sampling_multican.cpp -o sampling_multican
# ~/intel/oneapi/compiler/latest/linux/bin/intel64/icpc -fast MC_step.cpp -o MC_step -I ~/include -I ~/include/boost/

for delta in 0.05 0.1 0.2 0.3 0.4 #0.9 #0.8 #0.7 #0.6 #0.01 0.5 #0.6 0.7 0.8 0.9 1.0 0.01 #0.05 0.1 0.2 0.3 0.4
do
for genes in 4 5 6 #2 #4 #3 #4 5 6
do
# python L_noise_write.py $genes -1 1 $delta -1 1000
for num_noise_type in 1000
do
for max_fitness in 10 15 20 25 #4 5 #10 15 20 25 #24 25 26 #5 7 10 15 #20 #15 #3 4 5 6 7 #10 15 20 # 6 #5 6 8 9 12 #7 10 15 20
do
mkdir -p ../data_random_noise${delta}/genes_${genes}_MC_step/num_noise_type_${num_noise_type}/max=${max_fitness}/
bsub -e ../Error/error_MC_debugged_${delta}_${num_noise_type}_${max_fitness}.out  -m "cell6" -q "jobLimit60" "
# ./MC_step $genes $num_noise_type $max_fitness $delta &&
./MC_step_debugged $genes $num_noise_type $max_fitness $delta &&
"
# mv ../data_random_noise${delta}/genes_${genes}_MC_step_debugged/ ../data_random_noise${delta}/NG_genes_${genes}_MC_step_debugged/
# ./MC_step_debugged $genes $num_noise_type $max_fitness $delta &
done
done
done
done