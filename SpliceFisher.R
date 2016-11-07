args <- commandArgs(TRUE)

table <- read.table(args[1], header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)

parts <- c("Head", "Tail", "Body")
for(part in parts) {
	count.table <- table[, grep(paste("^count", part, sep = ""), colnames(table))]
	count.table[, ] <- unlist(lapply(as.character(unlist(count.table)), function(x) {sum(as.numeric(unlist(strsplit(x, ','))))}))
	tests <- apply(count.table, 1, function(x) {fisher.test(matrix(x, nrow = 2))})
	table[, paste("pvalue",    part, sep = "")] <- unlist(lapply(tests, function(x) {x$p.value}))
	table[, paste("oddsratio", part, sep = "")] <- unlist(lapply(tests, function(x) {x$estimate}))
	table[apply(count.table, 1, function(x) {any(c(rowSums(matrix(x, nrow = 2)), colSums(matrix(x, nrow = 2))) == 0)}), paste(c("pvalue", "oddsratio"), part, sep = "")] <- NA
}

table[, paste("padjust", c("Head", "Tail"), sep = "")] <- p.adjust(data.matrix(table[, paste("pvalue", c("Head", "Tail"), sep = "")]), method = "fdr")
table[, paste("padjust", "Body", sep = "")] <- p.adjust(table[, paste("pvalue", "Body", sep = "")], method = "fdr")

table[, "change"] <- 0
table[apply(table[, paste("oddsratio", parts, sep = "")], 1, function(x) {all(!is.na(x) & x > 1)}), "change"] <- -1
table[apply(table[, paste("oddsratio", parts, sep = "")], 1, function(x) {all(!is.na(x) & x < 1)}), "change"] <- 1

colnames <- c("chromosome", "start", "end", "strand", "gene", "change")
colnames <- c(colnames, apply(expand.grid(c("pvalue", "padjust", "oddsratio"), parts), 1, function(x) {paste(x, collapse = "")}))
colnames <- c(colnames, colnames(table)[grep("^count", colnames(table))])
table <- table[, colnames]

write.table(table, file = args[2], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
