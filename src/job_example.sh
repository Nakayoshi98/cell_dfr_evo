mkdir -p ../perturbation_lib
mkdir -p ../data_random_noise${delta}/genes_${genes}_step/epsilon_${epsilon}/num_noise_type_${num_noise_type}_without_e/network_info_step/seed_${seed}/GA_ignoring_lc
python noise_write.py $genes -1 $seed $delta -1 $L 
./GA $genes $num_noise_type $seed $delta $epsilon $L
