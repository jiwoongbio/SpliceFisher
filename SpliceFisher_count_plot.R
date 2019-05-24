args <- commandArgs(TRUE)

table <- read.table(args[1], header = FALSE, sep = "\t")
rownames(table) <- table[, 1]
table <- table[, -1]

png(filename = args[2], width = as.numeric(args[3]), height = as.numeric(args[4]), type = "cairo")
par(mar = c(par()$mar[1] + as.numeric(args[5]), par()$mar[2], par()$mar[3], par()$mar[4]))
barplot(t(table), beside = TRUE, las = 2)
dev.off()
