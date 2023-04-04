#!/usr/bin/env Rscript

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if (!requireNamespace("CMplot", quietly = TRUE)) install.packages("CMplot")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")


library("ggplot2")
library("RColorBrewer")
library("CMplot")
library("dplyr")

prepare_snp <- function(FILE) {
  table <- read.table(FILE, header = T, sep = "\t")
  
  table$chr[table$chr == "X"] <- 23
  table$chr[table$chr == "Y"] <- 24
  
  snp <- table$rsid
  chr <- as.numeric(table$chr)
  pos <- as.numeric(table$pos)
  pval <- as.numeric(table$pval)
  final_table <- data.frame(snp, chr, pos, pval)
  return(final_table)
  
}

draw_mhplot_qq <- function(table, name, format, color, SNPs, genes, max_pval = 8) {
  CMplot(table,
         plot.type = "m", col = c("grey40", "grey70"), # chr colors
         highlight = SNPs,
         highlight.text = genes,
         highlight.col = c("#fb8072"),
         highlight.cex = 1, highlight.pch = c(16),
         LOG10 = TRUE, ylim = c(0, max_pval), # limits of log pval
         threshold = c(0.05 / nrow(table), 1e-4), # cut-offs of pval
         threshold.lty = c(1, 2), threshold.lwd = c(1, 1),
         threshold.col = c("black", "grey"), # threshold colors
         amplify = TRUE, chr.den.col = NULL,
         signal.col = c("#fb8072", "#b3de69"), # colors of significant
         signal.cex = c(1.5, 1.5), signal.pch = c(19, 19),
         file = format, # file format
         memo = name, # file postfix
         dpi = 500, file.output = TRUE, verbose = TRUE,
         width = 14, height = 5
  )
  CMplot(table,
         plot.type = "q", col = color, box = FALSE, file = format, memo = name, dpi = 500,
         conf.int = TRUE, conf.int.col = NULL, threshold.col = "red", threshold.lty = 2,
         file.output = TRUE, verbose = TRUE,
         width = 5, height = 3.5
  )
}


args = commandArgs(trailingOnly=TRUE)
PHENO_NAME <- args[1]
FILE <- args[2]
cat(PHENO_NAME)
cat("\n")
cat(FILE)
cat("\n")
# ./DRAW.R metal /media/array/pregnancy/data/extended_attempt1.TBL
# PHENO_NAME <- 'metal'
# FILE <- "/media/array/pregnancy/data/extended_attempt1.TBL"

colors <- c("#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")
COLOR <- sample(colors, 1)

pheno_name <- PHENO_NAME
color <- COLOR

table <- prepare_snp(FILE)

max_pval <- max(-log10(table$pval))
max_pval <- ceiling(max_pval)
cutoff <-  round(max_pval*0.8)

SNPs <- table[-log10(table$pval) > cutoff, 1]

draw_mhplot_qq(table,
           pheno_name,
           "pdf",
           color,
           SNPs,
           SNPs,
           max_pval)

quit(save='no')
