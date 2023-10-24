# GWAS data simulation

Human GWAS data simulator.


# Installation

To make this tool working, you need:
* `Python 3.9` (+ packages are listed in `requirements.txt`)
*  `R 4.1` (+ packages are listed in `requirements-R.txt`)
* `PLINK`, `PLINK2`, `bedtools` and `HAPGEN2` be installed. Pathes to these tools have to be presented in `dependencies.yaml` (use `-d` parameter to specify). The example is in `dependencies.yaml` file.

# Launch

Example of launching simulation:

```
./biogwas.py \
--input_dir data2  \
--data_dir data3  \
--img_dir images  \
--vcf_in_flag \
--input_list  data2/chr.list \
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
-id data2  \
-dd data3  \
-imd images  \
-vcf \
-il  data2/chr.list \
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


# Examples

Examples of usage as well as validation, ment in the paper is located in `tests` directory.




