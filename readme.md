# GWAS data simulation

Human GWAS data simulator.


# Installation

To make this tool working, you need:
* `Python 3.9` (+ packages are listed in `requirements.txt`, installation: `pip install -r requirements.txt`)
*  `R 4.1` (+ packages are listed in `requirements-R.txt`, quick installation: `install_r_reqs.R`)
* `PLINK`, `PLINK2`, `bedtools` and `HAPGEN2` be installed. Pathes to these tools have to be presented in `dependencies.yaml` (use `-d` parameter to specify). The example is in `dependencies.yaml` file.

# Launch

Example of launching simulation:

```bash
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

```bash
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

```bash
./biogwas.py --help
```

# Launch with Docker

In order to build image use Docker file:
```bash
docker build -t biogwas .
```


After that you can easily run our tool using `docker-compose`. Example of configuration is in `docker-compose.yaml` file. You have to change:

```
volumes:
  - "/media/MIRROR/ukb_finngen:/data" 
```

To point to the working dir.

Also you have to work with the command part:
```
command: >
/bioGWAS/biogwas.py
-d /dependencies.yml
-id "/data/1000genomes/data2"
-dd "/data/gwassim_check/attempt_docker/data"
-imd "/data/gwassim_check/attempt_docker/images"
--vcf_in_flag
-il "/data/1000genomes/data2/chr.list"
-if "/data/1000genomes/data2/EUR_SAMPLES_ID.txt"
-af "/data/1000genomes/data2/gencode.v37.annotation.gtf"
--gmt_file "/data/1000genomes/data2/h.all.v2023.1.Hs.symbols.gmt"
--causal_pathways "/data/1000genomes/data2/pathways.csv"
-p "dpat"
-cid "dcid"
-sid "dsid"
```

You have to leave first three strings as they are:
```
command: >
/bioGWAS/biogwas.py
-d /dependencies.yml
```

All other parameters should be changed to your settings. 

Note: not all settings and flags are shown here, to display the full list, change all parameters to: 
```
command: >
/bioGWAS/biogwas.py
--help
```

After editing `docker-compose.yaml`, run simulation using command: `docker-compose up`



# Examples

Examples of usage as well as validation, ment in the paper is located in `tests` directory.




