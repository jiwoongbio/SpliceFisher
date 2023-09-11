args <- commandArgs(TRUE)

table <- read.table(args[1], header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)

parts <- c("Head", "Tail", "Body")
sampleNumbers <- unique(as.numeric(sub("^.*_", "", grep("^count", colnames(table), value = TRUE))))
for(rowNumber in 1:nrow(table)) {
	test.list <- list()
	for(part in parts) {
		countPairs.list <- list()
		for(sampleNumber in sampleNumbers) {
			countPairs.list[[sampleNumber]] <- apply(sapply(table[rowNumber, paste("count", part, 1:2, "_", sampleNumber, sep = "")], function(x) {as.numeric(unlist(strsplit(x, ',')))}), 1, function(x) {paste(x, collapse = ",")})
		}
		test.list <- c(test.list, lapply(apply(expand.grid(countPairs.list), 1, function(x) {paste(x, collapse = ",")}), function(x) {fisher.test(matrix(as.numeric(unlist(strsplit(x, ","))), nrow = 2))}))
	}
	oddsratios <- unlist(lapply(test.list, function(x) {x$estimate["odds ratio"]}))
	change <- 0
	if(all(oddsratios > 1)) {
		change <- 1
	}
	if(all(oddsratios < 1)) {
		change <- -1
	}
	table[rowNumber, "change"] <- change
	table[rowNumber, "pvalue"] <- max(unlist(lapply(test.list, function(x) {x$p.value})))
}

colnames <- c("change", "pvalue")
colnames <- c(colnames, grep("^count", colnames(table), value = TRUE))
if("pairType" %in% colnames(table)) {
	colnames <- c("chromosome", "start1", "end1", "start2", "end2", "strand", "pairType", "gene", colnames)
} else {
	colnames <- c("chromosome", "start", "end", "strand", "gene", colnames)
}
table <- table[, colnames]

write.table(table, file = args[2], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
