#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(PhenotypeSimulator)

# plink2 --bfile chr_EUR_sim_snps --recode vcf --out chr_EUR_sim_snps
# plink2 --vcf chr_EUR_sim_snps.vcf --make-bed --out _chr_EUR_sim_snps

#args[1] <- wd
#args[2] <- readgenotypes
#args[3] <- save to
#args[4] <- # samples
#args[5] <- # of snps
#args[6] <- genvar
#args[7] <- h2s


#args<-c("./", 
#	"data/chr_filt_sim_snps",
#	"data/phenos_ph1.tsv",
#	"5000", 
#	"10", 
#	"0.5", 
#	"0.5", 
#	"0.5")
#	/home/achangalidi/ukb_finngen/1000genomes chr_EUR_sim_snps phenos.tsv 5000 20 0.5 0.5 0.5


WD <- args[1]
geno_file <- args[2]
pheno_file <- args[3]
N_samples <- as.numeric(args[4])
N_genes <-as.numeric(args[5])
genVar <- as.numeric(args[6]) #!!! 
h2s <- as.numeric(args[7]) # доля от генвар, которую составля.т эффекты моих снипов
shared <- as.numeric(args[8])
mBeta <- as.numeric(args[9])
sdBeta <- as.numeric(args[10])
totalSNPeffect <- genVar*h2s
noiseVar <- 1 - genVar
independent <- 1 - shared


setwd(WD)


total <- independent * noiseVar + shared * noiseVar +
  shared * h2s *genVar + independent * h2s * genVar +
  shared * (1-h2s) * genVar + independent * (1-h2s) * genVar


print("KEK")
print(total)
print(abs(total-1) < 0.00001)
print("LOL")
stopifnot(abs(total-1) < 0.00001)



genotypes <- readStandardGenotypes(N=N_samples, 
					filename=geno_file,
					format="plink")
causalSNPs <-getCausalSNPs(N=N_samples, 
				genotypes = genotypes$genotypes,
				NrCausalSNPs = N_genes, 
				verbose = TRUE)
genFixed <-geneticFixedEffects(N = N_samples, 
				P = 1, 
				X_causal = causalSNPs, 
				pIndependentGenetic = 1,
                pTraitIndependentGenetic = 1, 
                mBeta = mBeta, 
                sdBeta = sdBeta)


genotypes_sd <-standardiseGenotypes(genotypes$genotypes)
kinship <- getKinship(N=N_samples, 
                      X=genotypes_sd, 
                      verbose = FALSE)

genBg <-geneticBgEffects(N=N_samples, 
                         kinship = kinship, 
                         P = 1)

noiseBg <- noiseBgEffects(N=N_samples, 
                          P=1)


genFixed_shared_scaled <- rescaleVariance(genFixed$shared, 
                                          shared * h2s * genVar)
genFixed_independent_scaled <- rescaleVariance(genFixed$independent, 
                                               independent * h2s * genVar)

genBg_shared_scaled <- rescaleVariance(genBg$shared, 
                                       shared * (1-h2s) * genVar)
genBg_independent_scaled <- rescaleVariance(genBg$independent, 
                                            independent * (1-h2s) * genVar)

noiseBg_shared_scaled <- rescaleVariance(noiseBg$shared, 
                                         shared * noiseVar)
noiseBg_independent_scaled <- rescaleVariance(noiseBg$independent,
                                              independent *noiseVar)



Y <- scale(genBg_independent_scaled$component +
             noiseBg_independent_scaled$component + 
             genFixed_independent_scaled$component)

row.names(Y) <- genotypes[["id_samples"]]
Y <- as.data.frame(Y)
Y[["t1"]] <- Y[["Trait_1"]]
Y[["#FID"]] <- row.names(Y) 
Y[["IID"]] <- sub("id1", "id2", Y[["#FID"]])
write.table(Y[c("#FID", "IID", "t1")], file=pheno_file, quote=FALSE, sep='\t', col.names = TRUE, row.names=FALSE)

#write.table(Y[c("#FID", "Trait_1")], file=pheno_file, quote=FALSE, sep='\t', col.names = TRUE, row.names=FALSE) 


