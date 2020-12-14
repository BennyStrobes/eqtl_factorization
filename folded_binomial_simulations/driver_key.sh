# For pystan
module load python/2.7-anaconda
module load gcc/6.4.0

output_directory="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/folded_binomial_simulations/"

python folded_binomial_simulations.py $output_directory