## GWAS data simulation

Human GWAS data simulator using 1000 human genomes project data.

### Setup with Snakemake

#### 1. Configuration:

Configuration parameters are in `config.yaml`


#### 2. Launch

To launch pipeline, write: `snakemake --cores <X> -s Snakefile`

To pring DAG of rules: `snakemake --rulegraph -s Snakefile | dot -Tpdf > images/rg.pdf`


#### 3. Results:

Results are stored in:
* Images -- `images/`:
   * Q-Q plot: `images/<PATTERN>_<PHENOS_ID>_<SIM_ID>_gwas_qq.pdf`
   * Manhattan plot: `images/<PATTERN>_<PHENOS_ID>_{SIM_ID}_gwas_mh.pdf`
   * Genomes' PCAs:
       * Full set: `images/<PATTERN>_filt_pca.pdf`
       * Independent set: `images/<PATTERN>_filt_indep_pca`
* Data -- `data/` (data is not uploaded, since files are too large):
   * `data/<PATTERN>_filt_sim.vcf` -- simulated genomes file;
   * `data/<PATTERN>_<PHENOS_ID>_<SIM_ID>_phenos.vcf` -- simulated phenotypes file;
   * `data/<PATTERN>_<PHENOS_ID>_<SIM_ID>_gwas.tsv` -- association data (GWAS summary statistics).

## Extra-scripts for article

### Iterate over parameters

* `iterate_params.py` -- script for iteration over parameters and save all the results to `DATA_DIR`
* `iterate_results.py` -- script for iteration over results of parameters iteratin (previous script). It generates: 
   * `data/<PATTERN>_<PHENOS_ID>_compare_results.tsv` -- table with summary of comparison (with number of found causal SNPS, ...)
* `results_1_check_best_parameters.ipynb` -- visualization of the results of parameters enumeration.

### LSEA launching

* `iterate_lsea.py` -- script for launching LSEA a lot of times to validate results.
* `results_1_check_best_parameters.ipynb` -- visualization of the results of validating LSEA.
