args = commandArgs(trailingOnly=TRUE)
library(flexmix)



generate_eqtl_mixture_input <- function(G, Y) {
	num_tests = dim(G)[1]
	df = data.frame(G1=as.numeric(G[1,]), Y1=as.numeric(Y[1,]))
	for (test_num in 2:num_tests) {
		df[paste0("G", test_num)] = as.numeric(G[test_num,])
		df[paste0("Y", test_num)] = as.numeric(Y[test_num,])
	}
	return(df)
}

run_flexmix <- function(eqtl_mixture_input_df, num_tests, num_latent_factors) {
	model_list = list()

	for (test_num in 1:num_tests) {

		model_list[[test_num]] = FLXMRglm(as.formula(paste0("Y", test_num, " ~ G", test_num)))
	}

	# Run model

	res <- flexmix(formula = ~1, data=eqtl_mixture_input_df, k=num_latent_factors, model=model_list)
	return(res)
}


genotype_training_file = args[1]
expression_training_file = args[2]
sample_overlap_file = args[3]
num_latent_factors = as.numeric(args[4])
output_root = args[5]

#print(num_latent_factors)


if (FALSE) {
G = read.table(genotype_training_file, header=FALSE)
Y = read.table(expression_training_file, header=FALSE)

num_tests <- dim(G)[1]

eqtl_mixture_input_df <- generate_eqtl_mixture_input(G, Y)
print("START RUN")
flex_result <- run_flexmix(eqtl_mixture_input_df, num_tests, num_latent_factors)

output_file <- paste0(output_root, "flexmix_model.RDS")
saveRDS(flex_result, file=output_file)
}
output_file <- paste0(output_root, "flexmix_model.RDS")

flex_result <- readRDS(output_file)


posterior_probs <- posterior(flex_result)
output_file <- paste0(output_root, "flexmix_model_posterior_prob.txt")

write.table(posterior_probs, file=output_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

print(summary(posterior(flex_result)))

