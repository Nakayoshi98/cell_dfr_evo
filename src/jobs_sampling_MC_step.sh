#!/bin/bash
#for文 betaは任意の変数一般的に「beta」「state_sigma」「cost」「n」なとが使われることが多い
#inより後ろを順に代入してゆくことになる。
export CPATH=:~/include/

# g++ -O3 sampling_multican.cpp -o sampling_multican
# ~/intel/oneapi/compiler/latest/linux/bin/intel64/icpc -fast MC_step.cpp -o MC_step -I ~/include -I ~/include/boost/

for delta in 0.05 0.4 #0.1 0.2 0.3 0.4 0.45 0.46 0.47 0.48 #0.5 #0.1 0.2 #0.3 0.4 0.5
do
for genes in 5 #2 #4 #5
do
for num_noise_type in 1000
do
for max_fitness in 20  #6 7 #15 #3 #6 #3 4 5 6 7 #7 10 15 20 
do
for seed in `seq 51 150` #`seq 21 50` #`seq 1 3`
do
mkdir -p ../data_random_noise${delta}/genes_${genes}_MC_step/num_noise_type_${num_noise_type}/max=${max_fitness}/samples/seed_${seed}
bsub -e ../Error/error_MC_${delta}_${num_noise_type}_${max_fitness}_${seed}.out -m "hpc" -q "normalCPU" "
./sampling_MC_step $genes $num_noise_type $max_fitness $delta $seed &&
"
done
done
done
done
done