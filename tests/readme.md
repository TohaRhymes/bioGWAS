# Parameters validation

(Initial datasets as well, as all experiments are described in the Changalidis et al.,2023)

## 1. Default parameter fitting

Directory: `tests/1_validate_params`.

#### Files

* First round of validation and default parameter set fitting
    * `utils_1.py` - file to be importes in other files (with param sets and K sets, names of files);
    * `1.0_iterate_params.py` - script to iterate over params and run simulation;
    * `1.1_iterate_results.py` - script to iterate over results of simulation, clump SNPs and summarize results. 
    * `1.2_check_parameters.ipynb` -- vizualization of these results.
* Second round of the revalidation of selected sets, same 4 files:
    * `utils_2.py` - file to be importes in other files (with param sets and K sets, names of files);
    * `2.0_iterate_best_params.py` - script to iterate over params and run simulation;
    * `2.1_iterate_best_results.py` - script to iterate over results of simulation, clump SNPs and summarize results. 
    * `2.2_validate_parameters.ipynb` -- vizualization of these results.
* `utils.py` - common utilities.

#### Images

All images are located in `tests/1_validate_params/images`:

* `SF3.1.pdf` - `SF3.6.pdf` - parameters distribution for the first round of validation (code in `1.2_check_parameters.ipynb` )
* `extra_SF3.1.pdf` - `extra_SF3.6.pdf` - parameters distribution for the second round of validation (code in `2.2_validate_parameters.ipynb`)
* `SF4.1.pdf` and `SF4.2.pdf` - distributions of precisions and recalls for all experiments (code in `1.2_check_parameters.ipynb` and `2.2_validate_parameters.ipynb` respectively)


#### Data

All datasets are located in `tests/1_validate_params/data`:

* `test10000_hallmark__compare_results.tsv` and `test10000_best_hallmark__compare_results.tsv` - results after simulation and clumping (produced by `1.0_iterate_params.py`+`1.1_iterate_results.py` and `2.0_iterate_best_params.py`+`2.1_iterate_best_results.py` respectively)
* `aggregated_sim1.csv` and `aggregated_sim_best.csv` - aggregated results (aggregated in corresponding notebooks).



## 2. GWAS replication 

Directory: `tests/2_simulate_gwas`.

#### Files

* `fg_to_plink.py` and `ukb_to_plink.py` - scripts to map GWAS summary statistics from FinnGen and UKB respectively to PLINK format of summary statistics.
* `draw_mh_qq_for_initial.sh` - scripts to draw Manhattan and Q-Q plots for original datasets
* Simulation was conducted in standard manner (see other examples) with drawing option.


#### Images

All images are located in `tests/2_simulate_gwas/images`:

* `QQplot.pval_<X>.pdf`, `Rectangular-Manhattan.pval_<X>.pdf` - Q-Q and Manhatta plots for original UKB GWAS (M06), original FinnGen GWAS (I9), simulated UKB GWAS (SIM_M06), simulated FinnGen GWAS (SIM_I9)  

#### Data

All datasets are located in `tests/2_simulate_gwas/data`:

* Original significant SNPs in corresponding GWAS results from UKB ind Finngen: `I9_SNPs.txt` and `M06_SNPs.txt`
* Significant clumps in our simulation: `SIM_I9_S3_gwas_clumps.clumped` and `SIM_M06_S2_gwas_clumps.clumped`



## 3. Enrichment analysis simulation 

Directory: `tests/3_pathways`.

#### Files

* `0_path_pick.ipynb` - pick pathways for the analysis;
* `1_iterate_pathways.py` - iterate through pathways and simulate GWAS datasets;
* Launch MAGMA:
    * `2.0_prepare_for_magma.sh` - prepare files for MAGMA;
    * `2.1_iterate_magma.py` - make bash-scripts (in `magma_scripts` directory) to launch three models on all simulated data;
    * `2.2_launch_ea.sh` - launch magma using scripts, prepared in previous file.
* Launch Pascal (Pathway scoring algorithm):
    * `3.0_mapping_rsid.py` - map simulated GWAS datasets' positions to rsids.
    * `3.1_batch_pascal.sh` - launch Pascal.
* Aggregation of all enrichment data: `4.0_pascal_process.ipynb` and `4.1_path_results.ipynb`
    


#### Images

All images are located in `tests/3_pathways/images`:

* `X`

#### Data

All datasets are located in `tests/3_pathways/data`:

* `X`


