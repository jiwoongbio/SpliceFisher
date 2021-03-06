args <- commandArgs(TRUE)

table <- read.table(args[1], header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)

parts <- c("Head", "Tail", "Body")
for(part in parts) {
	count.table <- table[, grep(paste("^count", part, sep = ""), colnames(table))]
	table[, paste("kruskal_pvalue", part, sep = "")] <- apply(count.table, 1, function(x) {
		ratiosList <- lapply(as.list(data.frame(matrix(as.character(unlist(x)), nrow = 2), stringsAsFactors = FALSE)), function(y) {
			counts <- as.numeric(unlist(strsplit(y[2], ',')))
			ratios <- counts / (counts + as.numeric(unlist(strsplit(y[1], ','))))
			ratios[!is.na(ratios)]
		})
		lengths <- unlist(lapply(ratiosList, length))
		if(all(lengths > 0)) {
			kruskal.test(unlist(ratiosList), rep(1:length(lengths), lengths))$p.value
		} else {
			NA
		}
	})
	count.table[, ] <- unlist(lapply(as.character(unlist(count.table)), function(x) {sum(as.numeric(unlist(strsplit(x, ','))))}))
	tests <- apply(count.table, 1, function(x) {fisher.test(matrix(x, nrow = 2))})
	table[, paste("pvalue",    part, sep = "")] <- unlist(lapply(tests, function(x) {x$p.value}))
	table[, paste("oddsratio", part, sep = "")] <- unlist(lapply(tests, function(x) {x$estimate}))
	table[apply(count.table, 1, function(x) {any(c(rowSums(matrix(x, nrow = 2)), colSums(matrix(x, nrow = 2))) == 0)}), paste(c("pvalue", "oddsratio"), part, sep = "")] <- NA
}

table[, paste("padjust", c("Head", "Tail"), sep = "")] <- p.adjust(data.matrix(table[, paste("pvalue", c("Head", "Tail"), sep = "")]), method = "fdr")
table[, paste("padjust", "Body", sep = "")] <- p.adjust(table[, paste("pvalue", "Body", sep = "")], method = "fdr")

table[, paste("kruskal_padjust", c("Head", "Tail"), sep = "")] <- p.adjust(data.matrix(table[, paste("kruskal_pvalue", c("Head", "Tail"), sep = "")]), method = "fdr")
table[, paste("kruskal_padjust", "Body", sep = "")] <- p.adjust(table[, paste("kruskal_pvalue", "Body", sep = "")], method = "fdr")

table[, "change"] <- 0
table[apply(table[, paste("oddsratio", parts, sep = "")], 1, function(x) {all(!is.na(x) & x > 1)}), "change"] <- -1
table[apply(table[, paste("oddsratio", parts, sep = "")], 1, function(x) {all(!is.na(x) & x < 1)}), "change"] <- 1

colnames <- c("change")
colnames <- c(colnames, apply(expand.grid(c("pvalue", "padjust", "oddsratio", "kruskal_pvalue", "kruskal_padjust"), parts), 1, function(x) {paste(x, collapse = "")}))
colnames <- c(colnames, colnames(table)[grep("^count", colnames(table))])
if("pairType" %in% colnames(table)) {
	colnames <- c("chromosome", "start1", "end1", "start2", "end2", "strand", "pairType", "gene", colnames)
} else {
	colnames <- c("chromosome", "start", "end", "strand", "gene", colnames)
}
table <- table[, colnames]

write.table(table, file = args[2], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
