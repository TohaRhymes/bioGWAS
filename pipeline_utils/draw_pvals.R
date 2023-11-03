#!/usr/bin/env Rscript


if (!requireNamespace("svMisc", quietly = TRUE)) install.packages("svMisc")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if (!requireNamespace("CMplot", quietly = TRUE)) install.packages("CMplot")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library("svMisc")
library("ggplot2")
library("RColorBrewer")
library("CMplot")
library("dplyr")

# points settings
CEX=0.7

is_true_or_false <- function(check_var) {
  check_var_ <- as.character(check_var)
  
  check_var_ <- toupper(check_var_)
  
  if (check_var_ %in% c("1", "TRUE", "T")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

prepare_snp <- function(FILE, delete_sex) {
  print('Reading table...')
  table <- read.table(FILE, header = T, sep = "\t")
  print('Table imported!')
  if(delete_sex){
      table <- table[table$CHR != 23,]
      table <- table[table$CHR != 24,]
      table <- table[table$CHR != "X",]
      table <- table[table$CHR != "Y",]
  }else{
      table$CHR[table$CHR == "X"] <- 23
      table$CHR[table$CHR == "Y"] <- 24
  }
  
  snp <- table$SNP
  chr <- as.numeric(table$CHR)
  pos <- as.numeric(table$BP)
  pval <- as.numeric(table$P)
  final_table <- na.omit(data.frame(snp, chr, pos, pval))
  return(final_table)
  
}

try_read_file <- function(file_pval, CLUMP_FILE){
    SNPs <- tryCatch({
            print('tryin to clump')
            CLUMP_FILE = paste0(gsub(".qassoc", "", file_pval), '_clumps.txt')
            SNPs <- read.table(CLUMP_FILE, header = T, sep = "\t")
            SNPs <- SNPs$SNP
            print('Markers clumped')
            return(SNPs)
        }, error=function(e) {
            SNPs <- NA
            return(SNPs)
        }
    )
    return(SNPs)
}


str2bool <- function(input_str)
{
  if(input_str == "0" || input_str == "FALSE" || input_str == "False" || input_str == "false" || input_str == "F"  || input_str == "f"){
    input_str = FALSE
  }else{
    input_str = TRUE
  }
  return(input_str)
}

get_markers <- function(mutations, CUT_OFF, WINDOW) {
  if (sum(mutations$pval<=CUT_OFF, na.rm=TRUE)==0){
    return(c())
  }
  t_mutations <- mutations[mutations$pval<=CUT_OFF & !is.na(mutations$pval),]
  t_mutations$index <- rownames(t_mutations)
  rownames(t_mutations) <- 1:nrow(t_mutations)
  
  difference <- t_mutations$pos[-1] - t_mutations$pos[-nrow(t_mutations)]
  pair_index <- difference > 0 & difference < WINDOW
  
  index_to_take <- c()
  prev_diff <- FALSE
  prev_pval <- 1
  prev_index <- -1
  n_pairs <- nrow(t_mutations)
  for (i in seq_along(pair_index)){
    progress(i, n_pairs)
    if(pair_index[i]){
      prev_diff <- TRUE
      if (t_mutations[i,'pval'] < prev_pval){
        prev_pval <- t_mutations[i,'pval']
        prev_index <- t_mutations[i,'index']
      }
    }else{
      if(prev_diff){
        if (t_mutations[i,'pval'] < prev_pval){
          prev_pval <- t_mutations[i,'pval']
          prev_index <- t_mutations[i,'index']
        }
      }else{
        prev_index <- t_mutations[i,'index']
      }
      index_to_take <- c(index_to_take, c(prev_index))
      prev_diff <- FALSE
      prev_pval <- 1
      prev_index <- -1
    }
  }
  i <- nrow(t_mutations)
  progress(i, n_pairs)
  if(prev_diff){
    if (t_mutations[i,'pval'] < prev_pval){
      prev_pval <- t_mutations[i,'pval']
      prev_index <- t_mutations[i,'index']
    }
  }else{
    prev_index <- t_mutations[i,'index']
  }
  index_to_take <- c(index_to_take, c(prev_index))
  
  flag <- rownames(mutations) %in% index_to_take
  print('')
  return(mutations[flag,]$snp)
}


draw_qq <- function(table, name, format, color, SNPs, genes, max_pval = 8) {
  print("Q-Q printing...")
  CMplot(table,
         plot.type = "q", 
         col = color, 
         box = FALSE, 
         file = format, 
         memo = name, 
         dpi = 40,
         conf.int = TRUE, 
         conf.int.col = NULL, 
         threshold.col = "red", 
         threshold.lty = 2,
         cex=CEX,
         file.output = TRUE, 
         verbose = TRUE,
         width = 7, 
         height = 5
  )
  print("Q-Q printed")
}

draw_mh <- function(table, name, format, color, SNPs, genes, max_pval = 8) {
  print("MH printing...")
  CMplot(table,
         plot.type = "m", col = c("grey40", "grey70"), # chr colors
         highlight = SNPs,
         highlight.text = genes,
         highlight.col = c("#fb8072"),
         highlight.cex = CEX, 
         highlight.pch = c(16),
         cex=CEX,
         LOG10 = TRUE, 
         ylim = c(0, max_pval), # limits of log pval
         threshold = c(0.05 / nrow(table), 1e-4), # cut-offs of pval
         threshold.lty = c(1, 2), threshold.lwd = c(1, 1),
         threshold.col = c("black", "grey"), # threshold colors
         amplify = TRUE, chr.den.col = NULL,
         signal.col = c("#fb8072", "#b3de69"), # colors of significant
         signal.cex = c(CEX+1, CEX+1), 
         signal.pch = c(19, 19),
         file = format, # file format
         memo = name, # file postfix
         dpi = 40, 
         file.output = TRUE, 
         verbose = TRUE,
         width = 14, 
         height = 5
  )
  print("MH printed")
}


colors <- c("#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")
color <- sample(colors, 1)


args = commandArgs(trailingOnly=FALSE)

file_path = ""
for (arg in args){
    if (grepl('--file', arg)){
        file_path=dirname(gsub("--file=", "", arg))
        break
    }
}


args = commandArgs(trailingOnly=TRUE)
pheno_name <- args[1]
file_pval <- args[2]
clump <- str2bool(args[3])
file_binary <- args[4]
PLINK_PATH <- args[5]
delete_sex <- args[6]

# file_path <- "../../1000genomes/pipeline_utils"
# pheno_name <- "M06"
# file_pval <- "data/M06.gwas.imputed_v3.both_sexes.qassoc"
# clump <- TRUE
# file_binary <- "../1000genomes/data3/PATTERN_filt_sim"
# PLINK_PATH <- "/home/achangalidi/tools/plink"
# delete_sex <- TRUE


if(is.na(delete_sex)){
  delete_sex <- FALSE
}else{
  delete_sex <- is_true_or_false(delete_sex)
}



        
print(pheno_name)
print(file_pval)

table <- prepare_snp(file_pval, delete_sex)

head(table)

# There are some times pvals = 0 => -log10(pval) = inf. 
# We make these pvals equal to min pval that are not 0 (~e-200 -- e-300) 
zero_flag <- table$pval==0 
if (sum(zero_flag, na.rm=TRUE)>0){
  non_zero_min <- min(table[!zero_flag,]$pval, na.rm=TRUE)
  table[zero_flag,]$pval <- non_zero_min
}

max_pval <- max(-log10(table$pval), na.rm=TRUE)
max_pval <- ceiling(max_pval)
# p-value cut off for significant SNPs
cutoff <-  0.05/nrow(table)



print('getting markers...')
##  getting array of boolean: markers to sign: ##
if(clump){
    #R2 significance
    R2 <- 0.1
    # loci interval (only one SNP with max pval will be signed)
    WINDOW <- 5000
        
    get_clumped <- paste0(file_path, "/./clumped.py ", 
                          file_binary, " ", 
                          file_pval, " ", 
                          cutoff, " ",
                          R2, " ",
                          WINDOW, " ",
                          PLINK_PATH,  " ")
    
    print(get_clumped)
    system(get_clumped)
    
    SNPs <- try_read_file(file_pval, CLUMP_FILE)

    if(is.null(SNPs)){
            clump <- FALSE
    }
            
}

if (! clump){
    # loci interval (only one SNP with max pval will be signed)
    WINDOW <- 500000
    # maximum SNPs to sign (per chr)
    N_MARKERS <- 7
    # p-value cut off for significant SNPs
    CUT_OFF <- cutoff

    SNPs <- c()
    for (chr in unique(table$chr)){
      cur_window <- WINDOW
      flag <- table$chr == chr
      markers <- get_markers(table[flag,], CUT_OFF, cur_window)
      while(length(markers) > N_MARKERS){
        cur_window <- 2 * cur_window
        markers <- get_markers(table[flag,], CUT_OFF, cur_window)
      }
      SNPs <- c(SNPs, markers)
    }
    print('Markers got')
    
}
    
               
draw_qq(table[,c('snp', 'chr', 'pos', 'pval')],
           pheno_name,
           "pdf",
           color,
           NULL,
           NULL,
           max_pval)       

draw_mh(table[,c('snp', 'chr', 'pos', 'pval')],
           pheno_name,
           "pdf",
           color,
           SNPs,
           SNPs,
           max_pval)

quit(save='no')
