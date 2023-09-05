## GWAS data simulation

Human GWAS data simulator using 1000 human genomes project data.


### Launch

Example of launching simulation:

```
./biogwas.py \
--vcf_dir data2  \
--data_dir data3  \
--img_dir images  \
--vcfs_list  data2/chr.list \
--ids_file  data2/EUR_SAMPLES_ID.txt \
--anno_file data2/gencode.v37.annotation.gtf \
--gmt_file data2/h.all.v2023.1.Hs.symbols.gmt  \
--causal_pathways data2/pathways.csv \
--pattern PATTERN  \
--causal_id CAUSAL_ID  \
--sim_id SIM_ID 
```

Or shortly:

```
./biogwas.py \
-vd data2  \
-dd data3  \
-id images  \
-vl  data2/chr.list \
-if  data2/EUR_SAMPLES_ID.txt \
-af data2/gencode.v37.annotation.gtf \
-gf data2/h.all.v2023.1.Hs.symbols.gmt  \
-cp data2/pathways.csv \
-p PATTERN  \
-cid CAUSAL_ID  \
-sid SIM_ID 
```


Full list of parameters with its description and default values can be found here:

```
./biogwas.py --help
```


### Launch with Snakemake

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

### Extra-scripts for article

#### Iterate over parameters

* `iterate_params.py` -- script for iteration over parameters and save all the results to `DATA_DIR`
* `iterate_results.py` -- script for iteration over results of parameters iteratin (previous script). It generates: 
   * `data/<PATTERN>_<PHENOS_ID>_compare_results.tsv` -- table with summary of comparison (with number of found causal SNPS, ...)
* `results_1_check_best_parameters.ipynb` -- visualization of the results of parameters enumeration.

#### LSEA launching

* `iterate_lsea.py` -- script for launching LSEA a lot of times to validate results.
* `results_2_check_LSEA.ipynb` -- visualization of the results of validating LSEA.
