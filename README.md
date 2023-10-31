Evolution of hierarchy and irreversibility in theoretical cell differentiation model
====

Overview

## Description

## Requirement
C++: Eigen 3.3.8 

python: itertools, numpy
## Usage
cd src

mkdir -p ../perturbation_lib

python noise_write.py $genes $num_noise_type $seed $delta  

### Genetic algorithm
mkdir -p ../data_random_noise\${delta}/genes_\${genes}_step/num_noise_type_\${num_noise_type}/seed_\${seed}/

./GA $genes $num_noise_type $seed $delta

### Multicanonical Monte Carlo (Wang-Landau method)
mkdir -p ../data_random_noise\${delta}/genes_\${genes}_MC_step/num_noise_type_\${num_noise_type}/max=\${max_fitness}/

./MC $genes $num_noise_type $max_fitness $delta

### Multicanonical Monte Carlo (sampling)
mkdir -p ../data_random_noise\${delta}/genes_\${genes}_MC_step/num_noise_type_\${num_noise_type}/max=\${max_fitness}/samples/seed_\${seed}

./sampling_MC_step $genes $num_noise_type $max_fitness $delta $seed 

## Licence

[MIT](https://github.com/Nakayoshi98/cell_dfr_evo/blob/main/LICENSE)

## Reference
[1], 'Evolution of hierarchy and irreversibility in theoretical cell differentiation model' 

