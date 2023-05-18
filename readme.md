## GWAS data simulation

Human GWAS data simulator using 1000 human genomes project data.

### Setup with Snakemake

#### 1. Configuration:

Configuration parameters are in `config.yaml`


#### 2. Launch

To launch pipeline, write: `snakemake --cores <X> -s Snakefile`

To pring DAG of rules: `snakemake --rulegraph -s Snakefile | dot -Tpdf > images/rg.pdf`


#### 3. Results:
Images:
* Q-Q plot: `images/<PATTERN>_<PHENO>_gwas_qq.pdf`
* Manhattan plot: `images/<PATTERN>_<PHENO>_gwas_mh.pdf`
* Genomes' PCAs:
    * Full set: `images/<PATTERN>_filt_pca.pdf`
    * Independent set: `images/<PATTERN>_filt_indep_pca`
    
    
Data:
* `data/...`