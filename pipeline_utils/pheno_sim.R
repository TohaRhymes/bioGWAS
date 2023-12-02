#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(PhenotypeSimulator)

# parse args (more info -- in PhenotypeSimulator documentation)
WD <- args[1]
geno_file <- args[2]
pheno_file <- args[3]
N_samples <- as.numeric(args[4])
N_genes_file <-args[5]

N_genes <- nrow(read.table(N_genes_file, header = F, sep = "\t"))

print(N_samples)
print(N_genes)

mBeta <- as.numeric(args[6])
sdBeta <- as.numeric(args[7])

genVar <- as.numeric(args[8]) # noiseVar = 1-genVar 

# Genetic
h2s <- as.numeric(args[9]) # a fraction of the genvar that the effects of my snips are.
pIndependentGenetic <- as.numeric(args[10]) #proportion of genetic variant effects to have a trait-independent fixed effect
theta  <- as.numeric(args[11]) #  proportion of variance of shared genetic variant effects

# Noise
phi <- as.numeric(args[12])
alpha <- as.numeric(args[13])


seed <- as.numeric(args[14])
if (is.na(seed)){
    seed <- 566
}

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
Y[["t1"]] <- Y[["Trait_1"]]
Y[["#FID"]] <- row.names(Y) 
Y[["IID"]] <- sub("id1", "id2", Y[["#FID"]])

write.table(Y[c("#FID", "IID", "t1")], file=pheno_file, quote=FALSE, sep='\t', col.names = TRUE, row.names=FALSE)

print(paste0('Pheno wrote to: ', pheno_file))

image_name <- paste0(pheno_file, ".RData")

save.image(file = image_name)

