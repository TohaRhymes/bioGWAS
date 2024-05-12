#!/usr/bin/env Rscript


if (!requireNamespace("CMplot", quietly = TRUE)) install.packages("CMplot")

library("CMplot")

# points settings
CEX=0.7
# colors for q-q
colors <- c("#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")



# function checks check_var for true/false value
is_true_or_false <- function(check_var) {
  if(is.na(check_var)){
    return(FALSE)
  }
    
  check_var_ <- as.character(check_var)
  check_var_ <- toupper(check_var_)
  if (check_var_ %in% c("1", "TRUE", "T")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# open summstats file and prepare for CMplot
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

# try read clumps file, if there is no - return NA
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

# from the table of SNPs get significant ones (one per window)
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

# draw Q-Q and Manhattan using CMPlot
draw_qq <- function(table, name, format, color, SNPs, genes, max_pval = 8, width=7, height=5, dpi=40) {
  print("Q-Q printing...")
  CMplot(table,
         plot.type = "q", 
         col = color, 
         box = FALSE, 
         file = format, 
         memo = name, 
         conf.int = TRUE, 
         conf.int.col = NULL, 
         threshold.col = "red", 
         threshold.lty = 2,
         cex=CEX,
         file.output = TRUE, 
         verbose = TRUE,
         width = width, 
         height = height,
         dpi = dpi
  )
  print("Q-Q printed")
}

draw_mh <- function(table, name, format, color, SNPs, genes, max_pval = 8, width=14, height=5, dpi=40) {
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
         file.output = TRUE, 
         verbose = TRUE,
         width = width, 
         height = height,
         dpi = dpi
  )
  print("MH printed")
}


# === script started ===

#pick color for mh
color <- sample(colors, 1)


#parse args
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
clump <- is_true_or_false(args[3])
file_binary <- args[4]
PLINK_PATH <- args[5]
delete_sex <- is_true_or_false(args[6])

qq_width <- args[7]
qq_height <- args[8]
qq_dpi <- args[9]
mh_width <- args[10]
mh_height <- args[11]
mh_dpi <- args[12]


#read and parse table with snps

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



##  getting array of boolean: markers to sign: ##
print('getting markers...')
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

#now we have significant SNPs, we can draw them
    
               
draw_qq(table[,c('snp', 'chr', 'pos', 'pval')],
           pheno_name,
           "pdf",
           color,
           NULL,
           NULL,
           max_pval,
           qq_width, 
           qq_height, 
           qq_dpi)       

draw_mh(table[,c('snp', 'chr', 'pos', 'pval')],
           pheno_name,
           "pdf",
           color,
           SNPs,
           SNPs,
           max_pval,
           mh_width, 
           mh_height, 
           mh_dpi)

quit(save='no')
