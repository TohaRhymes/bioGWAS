#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (!require("PhenotypeSimulator")) install.packages("PhenotypeSimulator")  #0.3.4

library(PhenotypeSimulator)

# args<-c("./", 
# 	"data2/chr_ph_ss_K10_filt_sim_snps",
# 	"data2/chr_ph_ss_K10_m0.5_sd0.3_gv0.5_h2s0.5_phenos.tsv",
# 	"5000", 
# 	"10", 
# 	"0.5", 
# 	"0.5", 
# 	"1.0",
#        "0.5",
#        "0.3")
#	/home/achangalidi/ukb_finngen/1000genomes chr_EUR_sim_snps phenos.tsv 5000 20 0.5 0.5 0.5


WD <- args[1]
geno_file <- args[2]
pheno_file <- args[3]
N_samples <- as.numeric(args[4])
N_genes <-as.numeric(args[5])

mBeta <- as.numeric(args[6])
sdBeta <- as.numeric(args[7])

genVar <- as.numeric(args[8]) # noiseVar = 1-genVar 

# Genetic
h2s <- as.numeric(args[9]) # доля от генвар, которую составля.т эффекты моих снипов
pIndependentGenetic <- as.numeric(args[10]) #proportion of genetic variant effects to have a trait-independent fixed effect
theta  <- as.numeric(args[11]) #  proportion of variance of shared genetic variant effects

# Noise
phi <- as.numeric(args[12])
alpha <- as.numeric(args[13])


seed <- as.numeric(args[14])
if (is.na(seed)){
    seed <- 219453
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

