## GWAS data simulation

Human GWAS data simulator using 1000 human genomes project data.

### 1. Getting simulated data
* [one-time] Downloaded files from 1khg web-site: `./download_script.sh`
    * Before filtering by MAF, there were a lot of snips (~90% were gone), and hapgen worked super long, so it was a good idea to filter right before generation (if you need another cutoff, it's better to filter right away) -- next step.
* [one-time] Preprocess files: `./script_make_eur_haplo.sh > script_make_eur_haplo.log`:
    * ungzip, beautify headers
    * Filter: only 2-alleles SNPs, EUR samples, MAF > 0.05
    * --export ped + --export hapslegend
* Run hapgen: `./script_hapgen.sh`
* Correct files: `./script_rename_hap.sh > script_rename_hap.log`
### 2. Merge this data
Merged all bfiles to one and made vcf & ucsc’s bed.
These 2 steps in one script: `./script_make_sim_bfile.sh > script_make_sim_bfile.log`

### 3. Annotation 
* [one-time] `./script_download_gtf.sh` - download gtf file with genes
* `./script_make_gtf.sh` - prepare gtf file with genes
* `bedtools closest -d -a filtered_chr_EUR_sim_ucsc.bed -b gen_gencode.v37.annotation.gene.sorted.gtf  > chr_EUR_sim.sorted.annotated.bed` -- annotate bfile with bedtools closest.

### 4. Getting set of snps

Using annotated SNPs and gmt-file with list of pathways we can select list of causal SNPs. Script:  `./get_snps_set.py`

### 5. Extracting specific SNPs from bfile for simulation

Just use fiule from previous step and plink: `plink2 --bfile chr_EUR_sim --extract range snps_chr_EUR.tsv --make-bed --out chr_EUR_sim_snps`

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

