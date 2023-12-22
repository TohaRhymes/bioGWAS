# GWAS data simulation

Human GWAS data simulator from Changalidis et al., 2023.

Examples of usage as well as experiments described in the article are presented in the `tests/` directory (with descriptions in [tests/readme.md](tests/readme.md)).


# Installation

To make this tool working, you need:
* `Python 3.9` (+ packages are listed in `requirements.txt`, installation: `pip install -r requirements.txt`)
* `R 4.1` (+ packages are listed in `requirements-R.txt`, quick installation: `install_r_reqs.R`)
* `PLINK`, `PLINK2`, `bedtools` and `HAPGEN2` be installed. Pathes to these tools have to be presented in `dependencies.yaml` (use `-d` parameter to specify). The example is in `dependencies.yaml` file.

# Launch

The only __required__ parameters are:
* `-id/--input_dir` - path to the directory which contains input genotypes.
* `-il/--input_list` - csv-file that describes one input file per line with structure: `<file>,<chromosome_number>`.
* `-if/--ids_file` - file with samples ids (from input files) to use for sampling (one id per line).
* `-af/--anno_file` -- file with chromosome annotations in gtf-format.


Full list of parameters with its description and default values can be found here:

```bash
./biogwas.py --help
```


Example of launching simulation:

```bash
./biogwas.py \
--input_dir <input_data_dir>  \
--data_dir <output_data_dir>  \
--img_dir <output_images_dir>  \
--vcf_in_flag \
--input_list  <path_to_genofile_list>.list \
--ids_file  <path_to_samples_list>.txt \
--anno_file <path_to_gtf>.gtf \
--gmt_file <path_to_gmt>.gmt  \
--causal_pathways <path_to_pathways>.csv \
--pattern PATTERN  \
--causal_id CAUSAL_ID  \
--sim_id SIM_ID 
```

Or shortly:

```bash
./biogwas.py \
-id <input_data_dir>  \
-dd <output_data_dir>  \
-imd <output_images_dir>  \
-vcf \
-il  <path_to_genofile_list>.list \
-if  <path_to_samples_list>.txt \
-af <path_to_gtf>.gtf \
-gf <path_to_gmt>.gmt  \
-cp <path_to_pathways>.csv \
-p PATTERN  \
-cid CAUSAL_ID  \
-sid SIM_ID 
```


# Docker

In order to build image use Docker file:
```bash
docker build -t biogwas .
```

### Launch from CLI

After that you can easily run docker by:
```bash
docker run \
-v "<work_dir>:<work_dir_inside_container>" \
biogwas \
/bioGWAS/biogwas.py \
-d /dependencies.yaml \
-id <input_data_dir_inside_container>  \
-dd <output_data_dir_inside_container>  \
-imd <output_images_dir_inside_container>  \
-vcf \
-il  <path_to_genofile_list_inside_container>.list \
-if  <path_to_samples_list_inside_container>.txt \
-af <path_to_gtf_inside_container>.gtf \
-gf <path_to_gmt_inside_container>.gmt  \
-cp <path_to_pathways_inside_container>.csv \
-p "dpat" \
-cid "dcid" \
-sid "dsid"
```

Your parameters goes after 5th string, i.e. you have to leave first 5 strings as they are:
```bash
docker run \
-v "<work_dir>:<work_dir_inside_container>" \
biogwas \
/bioGWAS/biogwas.py \
-d /dependencies.yaml
```

All other flags can be changed. To read manual, simply run with flag `--help`:

```bash
docker run \
-v "<work_dir>:<work_dir_inside_container>" \
biogwas \
/bioGWAS/biogwas.py \
--help
```


### Launch with docker-compose.yaml

Alternatively you can launch bioGWAS using `docker-compose`. Example of configuration is in `docker-compose.yaml` file. You have to change:

```
volumes:
  - "<path_to_workdir>:<path_inside_container>" 
```

To point to the working dir.

You have to work with the command part:
```
command: >
/bioGWAS/biogwas.py
-d /dependencies.yaml
-id <input_data_dir_inside_container>  \
-dd <output_data_dir_inside_container>  \
-imd <output_images_dir_inside_container>  \
-vcf \
-il  <path_to_genofile_list_inside_container>.list \
-if  <path_to_samples_list_inside_container>.txt \
-af <path_to_gtf_inside_container>.gtf \
-gf <path_to_gmt_inside_container>.gmt  \
-cp <path_to_pathways_inside_container>.csv \
-p "dpat"
-cid "dcid"
-sid "dsid"
```

You have to leave first three strings as they are:
```
command: >
/bioGWAS/biogwas.py
-d /dependencies.yaml
```

All other parameters should be changed to your settings. 

Note: not all settings and flags are shown here, to display the full list, change all parameters to: 
```
command: >
/bioGWAS/biogwas.py
--help
```

After editing `docker-compose.yaml`, run simulation using command: 
```bash
docker-compose up
```

You can also go inside container (not recommended):
```bash
docker run --rm -it --entrypoint /bin/bash biogwas
```


# Examples

Examples of usage as well as validation steps, described in the paper is located in `tests` directory.

## Step-by-step tutorial


### 1. Install Docker
The simple way to launch our tool is using [Docker](https://docs.docker.com/get-docker/).
If you use docker, there is no need to install all other packages, everything will be installed for you!


### 2. Download source code
Download this repo source code:

```bash
git clone https://github.com/TohaRhymes/bioGWAS.git
```

Go to the directory with the source code:

```bash
cd bioGWAS
```

Build an image:

```bash
docker build -t biogwas .
```


### 3. Make working directory

Create a working directory and browse to it, in linux:

```bash
mkdir biogwas_test
cd biogwas_test
```


### 4. Prepare essential files

For this example, let all necessary files be in `./data/` directory (otherwise you have to mount all directories with data when launching docker).

Essential files are:

* genotypes file(s) (in `.vcf` format, we will add support of initial `.bed` files in the nearest future).
	* you need to store them in one directory (let it be `./data/genotypes`)
	* you need to create a list of all necessary file(s) to be included (it can be just names, or full path), and the chromosome corresponding to this file, e.g.[the example](data/chr.list):
        ```txt
        ./data/genotypes/file1.vcf,chr1
        ./data/genotypes/file2.vcf,chr2
        ...
        ```
        Let this file be `./data/genotypes.list` in our example. In the nearest future we are planning to add support for direct passing file names to the command line. 

* txt-file with samples ids to be included in the analysis (one per line), [the example](data/EUR_SAMPLES_ID.txt):
    ```txt
    HG00096
    HG00097
    HG00099
    HG00100
    ...
    ```
    Let this file be `./data/samles.txt` in our example.

* Annotation file in gtf format. In our test we used comprehensive gene annotation downloaded from [gencode site](https://www.gencodegenes.org/human/release_37.html). Let this file be `./data/gencodes.gtf` in our example.
* You also need to provide:
	* If you want to use specific causal SNPs, a list of these SNPs, one per line, example:
        ```txt
        1:172643220
        3:128435895
        4:76045432
        4:87976387
        5:7891402
        ```
        Let this file be `./data/snps.txt` in our example.
	* If you want to use specific pathways, you need to:
		* use GMT files with your pathways (you can download hallmark, and other gmt files from [gsea-msigdb website](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)). Let this file be `./data/pathways.gmt` in our example.
		* A list of causal pathway(s), e.g.:
            ```txt
            KEGG_PPAR_SIGNALING_PATHWAY
            KEGG_LONG_TERM_DEPRESSION
            
            ```
            Let this file be `./data/pathways_list.txt` in our example.


So, in total the necessary files are:
* `./data/genotypes/` - directory with genotypes files
* `./data/genotypes.list` - list with genotypes files
* `./data/samles.txt` - list with samples
* `./data/gencodes.gtf` - annotation file
* One of two:
	* `./data/snps.txt` - list of causal SNPs
	* `./data/pathways.gmt` and `./data/pathways_list.txt` - gmt-file with pathways, and your causal pathways. 


### 5. Run simulation

(Let working directory inside container be `/wd`, and our data will be in `/wd/data`, however it can be anything.)


If you want to run using causal pathways:

```bash
docker run \
-v "./:/wd" \
biogwas \
/bioGWAS/biogwas.py \
--input_dir /path/to/genotypes  \
--data_dir /wd \
--img_dir /wd  \
--vcf_in_flag \
--input_list  /wd/data/genotypes.list \
--ids_file  /wd/data/samles.txt \
--anno_file /wd/data/gencodes.gtf \
--gmt_file /wd/data/pathways.gmt  \
--causal_pathways /wd/data/pathways_list.txt
```


If you want to run using causal SNPs:

```bash
docker run \
-v "./:/wd" \
biogwas \
/bioGWAS/biogwas.py \
--input_dir /path/to/genotypes  \
--data_dir /wd \
--img_dir /wd  \
--vcf_in_flag \
--input_list  /wd/data/genotypes.list \
--ids_file  /wd/data/samles.txt \
--anno_file /wd/data/gencodes.gtf \
--gmt_file /wd/data/pathways.gmt  \
--use_causal_snps \
--causal_snps /wd/data/snps.txt
```


### 6. Wait!

After finishing, all output data will be in `./data` dir.
