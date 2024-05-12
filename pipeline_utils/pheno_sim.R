#!/usr/bin/env Rscript

count_rows_in_file <- function(filename) {
  full_filename <- paste0(filename, ".fam")
  if (!file.exists(full_filename)) {
    return("File not found.")
  }
  count <- as.numeric(length(readLines(full_filename)))
  return(count)
}

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


binarize_vector <- function(Y, case_fraction) {
  # Length of the vector Y
  n <- length(Y)
  
  # Calculate the number of elements for the case portion
  n_case <- ceiling(case_fraction * n)
  
  # Sort Y to find the threshold
  sorted_Y <- sort(Y)
  
  # Find the cutoff point for the largest case_fraction part
  if (n_case == 0) {
    # All are 0 if case_fraction is 0
    cutoff_case <- Inf  # So no values are >= Inf
  } else {
    cutoff_case <- sorted_Y[n - n_case + 1]
  }
  
  # Transform Y to 0 or 1 based on the cutoff
  transformed_Y <- ifelse(Y < cutoff_case, 0, 1)
  
  return(transformed_Y)
}


args = commandArgs(trailingOnly=TRUE)

library(PhenotypeSimulator)

# parse args (more info -- in PhenotypeSimulator documentation)
WD <- args[1]
geno_file <- args[2]
pheno_file <- args[3]
N_genes_file <-args[4]


mBeta <- as.numeric(args[5])
sdBeta <- as.numeric(args[6])

genVar <- as.numeric(args[7]) # noiseVar = 1-genVar 

# Genetic
h2s <- as.numeric(args[8]) # a fraction of the genvar that the effects of my snips are.
pIndependentGenetic <- as.numeric(args[9]) #proportion of genetic variant effects to have a trait-independent fixed effect
theta  <- as.numeric(args[10]) #  proportion of variance of shared genetic variant effects

# Noise
phi <- as.numeric(args[11])
alpha <- as.numeric(args[12])


seed <- as.numeric(args[13])
if (is.na(seed)){
    seed <- 566
}

binary <- is_true_or_false(args[14])
case_fraction <- as.numeric(args[15])

if(genVar==0){
  h2s <- NULL
  pIndependentGenetic <- NULL
  theta <- NULL
}
if(genVar==1){
  phi <- NULL
  alpha <- NULL
}

setwd(WD)

N_samples <- count_rows_in_file(geno_file)   # args[4]
N_genes <- nrow(read.table(N_genes_file, header = F, sep = "\t"))
print(N_samples)
print(N_genes)

# Use Phenotype Simulator functions to simulate phenotypes
genotypes <- readStandardGenotypes(N=N_samples, 
					filename=geno_file,
					format="plink")
phenotype <- runSimulation(N = N_samples, 
                           P = 1, 
                           genotypefile = geno_file,
                           cNrSNP=N_genes, 
                           format="plink",
                           genVar = genVar, 
                           h2s = h2s, 
                           mBetaGenetic = mBeta, 
                           sdBetaGenetic = sdBeta,
                           pIndependentGenetic = pIndependentGenetic,
                           theta=theta,
                           phi = phi,
                           alpha=alpha,
                           verbose = TRUE,
                          seed=seed)

# Save to file
print("Saving...")
Y <- phenotype$phenoComponentsFinal$Y



row.names(Y) <- genotypes[["id_samples"]]
Y <- as.data.frame(Y)

if (binary) {
  Y[["t1"]] <- binarize_vector(Y[["Trait_1"]], case_fraction)
} else {
  Y[["t1"]] <- Y[["Trait_1"]]
}


Y[["#FID"]] <- row.names(Y) 
Y[["IID"]] <- sub("id1", "id2", Y[["#FID"]])

write.table(Y[c("#FID", "IID", "t1")], file=pheno_file, quote=FALSE, sep='\t', col.names = TRUE, row.names=FALSE)

print(paste0('Pheno wrote to: ', pheno_file))

# todo delete .RData
image_name <- paste0(pheno_file, ".RData")

save.image(file = image_name)

