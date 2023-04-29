## GWAS data simulation

Human GWAS data simulator using 1000 human genomes project data.

### 1. Getting simulated data
* [one-time] Download files from 1khg web-site: `./download_script.sh`
    * ungzip, beautify headers
* [one-time] Preprocess files: `./1_script_make_eur_haplo.sh ./data chr 0.05 ./EUR_SAMPLES_ID.txt  > 1_script_make_eur_haplo.log`:
    * Filter: only 2-alleles SNPs, EUR samples, MAF > 0.05
    * --export ped + --export hapslegend
* Run hapgen to make samplings: `./2_script_hapgen.sh ./data chr 5000`
* Correct files: `./3_script_rename_hap.sh ./data chr > 3_script_rename_hap.log`
### 2. Merge this data
Merged all bfiles to one and made vcf & ucsc’s bed.
These 2 steps in one script: `./4_script_make_sim_bfile.sh ./data chr 5010 > 4_script_make_sim_bfile.log`

### 3. Annotation 
* [one-time] `./download_gtf.sh ./data` - download gtf file with genes
* `./5_script_make_gtf.sh ./data gencode.v37.annotation` - prepare gtf file with genes
* `./6_bedtools_closest.sh ./data chr gencode.v37.annotation` -- annotate bfile with bedtools closest.

### 4. Getting set of snps

Using annotated SNPs and gmt-file with list of pathways we can select list of causal SNPs. Script:  `./get_snps_set.py data chr data/chr_filt_sim_anno.bed data/h.all.v2023.1.Hs.symbols.gmt data/pathways.csv 20 10`

### 5. Extracting specific SNPs from bfile for simulation

Just use fiule from previous step and plink: `plink2 --bfile chr_EUR_sim --extract range chr_chosen_snps.tsv --make-bed --out chr_EUR_sim_snps`

### 6. Simulation

Use PhenotypeSimulator and R-script `./pheno_sim.R`:
`./pheno_sim.R /home/achangalidi/ukb_finngen/1000genomes chr_EUR_sim_snps phenos.tsv 5000 20 0.5 0.5 0.5`

### 7. GWAS and PCA 

Script with GWAS and PCA (plink):
`./gwas_analysis.sh`

Draw:
* `./draw_pvals.R chr_EUR_sim_p1 chr_EUR_sim_p1.tsv` -- Q-Q plot and Manhattan plot;
* `./draw_pca.py chr_EUR_sim` -- PCA plots for full set of SNPS;
* `./draw_pca.py chr_EUR_sim_indep` -- PCA plots for independent set of SNPS.

## Results:
Pictures:
* Q-Q: `QQplot.pval_chr_EUR_sim_p1.pdf`
* Manhattan: `Rectangular-Manhattan.pval_chr_EUR_sim_p1.pdf`
* PCAs:
    * Full set: `chr_EUR_sim_indep_pca`
    * Independent set: `chr_EUR_sim_pca.pdf`

