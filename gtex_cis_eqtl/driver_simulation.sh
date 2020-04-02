#######################
# Used Directories
#######################
# Simulated results directory
simulated_results_dir="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/gtex_cis_eqtl/simulation/simulated_eqtl_results/"
# visualize simulated results directory
visualize_simulated_results_dir="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/gtex_cis_eqtl/simulation/visualize_simulated_eqtl_results/"


##############
# Parameters
##############
num_individuals="100"
num_samples_per_individual="10"
num_tests="1000"
num_latent_factors="7"

module load python/3.7.4-anaconda


###############
# Run model
################
svi="True"
parrallel="False"
seeds=("0")
for seed in "${seeds[@]}"; do
	echo "Seed: "$seed
	output_stem=$simulated_results_dir$num_individuals"_individuals_"$num_samples_per_individual"_samples_per_individual_"$num_tests"_tests_"$num_latent_factors"_latent_factors_"$svi"_svi_boolean_"$seed"_seed_parrallelized_"
	python simulate_and_run_factorization.py $output_stem $svi $num_individuals $num_samples_per_individual $num_tests $num_latent_factors $seed $parrallel
done





if false; then
svi="True"
seeds=("0")
for seed in "${seeds[@]}"; do
	echo "Seed: "$seed
	output_stem=$simulated_results_dir$num_individuals"_individuals_"$num_samples_per_individual"_samples_per_individual_"$num_tests"_tests_"$num_latent_factors"_latent_factors_"$svi"_svi_boolean_"$seed"_seed_parrallelized_"
	python simulate_and_run_factorization.py $output_stem $svi $num_individuals $num_samples_per_individual $num_tests $num_latent_factors $seed
done
fi


###############
# Visualize Results
################
if false; then
Rscript visualize_simulated_eqtl_factorization_results.R $simulated_results_dir $visualize_simulated_results_dir
fi
