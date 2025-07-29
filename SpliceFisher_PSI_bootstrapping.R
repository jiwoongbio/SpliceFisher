args <- commandArgs(TRUE)

input_file <- args[1]
output_file <- args[2]
comma_separated_controls <- args[3]
comma_separated_tests <- args[4]
boot_replications <- 1000

input_df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "")
id_colname <- colnames(input_df)[1]
rownames(input_df) <- input_df[, 1]
input_df <- input_df[, -1]
colnames <- colnames(input_df)
colnames(input_df) <- 1:length(colnames)

controls <- colnames(input_df)[colnames %in% unlist(strsplit(comma_separated_controls, ","))]
tests <- colnames(input_df)[colnames %in% unlist(strsplit(comma_separated_tests, ","))]

output_df <- as.data.frame(t(apply(input_df, 1, function(x) {
	pvalue <- wilcox.test(x[controls], x[tests])$p.value
	mean_diff <- mean(x[tests]) - mean(x[controls])
	boot_diffs <- replicate(boot_replications, {
		control_boot <- sample(x[controls], size = length(x[controls]), replace = TRUE)
		test_boot <- sample(x[tests], size = length(x[tests]), replace = TRUE)
		mean(test_boot) - mean(control_boot)
	})
	pvalue <- 1.0
	if(mean_diff > 0) {
		pvalue <- 2 * (sum(boot_diffs <= 0) / length(boot_diffs))
	}
	if(mean_diff < 0) {
		pvalue <- 2 * (sum(boot_diffs >= 0) / length(boot_diffs))
	}
	return(c(
		mean_diff = mean_diff,
		pvalue = pvalue
	))
})))
output_df[[id_colname]] <- rownames(output_df)
output_df$fdr <- p.adjust(output_df$pvalue, method = "fdr")

write.table(output_df[, c(id_colname, "mean_diff", "pvalue", "fdr")], file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
