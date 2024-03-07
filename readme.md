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
* Genotypes, one of these:
    * `-il/--input_list` - csv-file that describes one input file per line with structure: `<abs_path_to_file>,<chromosome_number>`.
    * `-gel/--geno_list` and  `-chl/--chr_list` and specify genotype file(s) and chromosome(s) respectively:
```text
    --geno_list /wd/data/genotypes/file1.vcf,/wd/data/genotypes/file2.vcf
    --chr_list chr1,chr2
```
* `-af/--anno_file` -- file with chromosome annotations in gtf-format.
* Dependencies file with tools pathes (set by `-d/--dependencies`).
* Output directories for data and images (default: `./data` and `./images`) -- they should be created.
* For choosing causal SNPs/pathways:
    * `-ucs/--use_causal_snps` - has to be set or not set.
    * then `-cp/--causal_pathways` OR `-cs/--causal_snps` for specifying causal pathways/SNPs respectively.
    * if we use causal pathways (i.e. `use_causal_snps` is not set), we also have to specify `--gmt_file`.


Full list of parameters with its description and default values can be found here:

```bash
./biogwas.py --help
```


Example of launching simulation:

```bash
./biogwas.py \
--data_dir <output_data_dir>  \
--img_dir <output_images_dir>  \
--input_list  <path_to_genofile_list>.list \
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
-dd <output_data_dir>  \
-imd <output_images_dir>  \
-il  <path_to_genofile_list>.list \
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

After that you can easily run docker by something like:
```bash
docker run \
-v "<work_dir>:<work_dir_inside_container>" \
biogwas \
/bioGWAS/biogwas.py \
-d /dependencies.yaml \
-dd <output_data_dir_inside_container>  \
-imd <output_images_dir_inside_container>  \
-il  <path_to_genofile_list_inside_container>.list \
-af <path_to_gtf_inside_container>.gtf \
-gf <path_to_gmt_inside_container>.gmt  \
-cp <path_to_pathways_inside_container>.csv \
-p "dpat" \
-cid "dcid" \
-sid "dsid"
```

Your parameters for biogwas goes after 5th string, i.e. you have to leave first 5 strings as they are (except `-v` flag which shows how to mount directories inside container, you **have** to specify all directories you're using: the easiest way  is to make just one directory with subdirectories for I/O, and specify this one directory):


All other flags can be changed. To read manual, simply run with flag `--help`:

```bash
docker run \
-v "<work_dir>:<work_dir_inside_container>" \
biogwas \
/bioGWAS/biogwas.py \
--help
```


### Launch with docker-compose.yaml

Alternatively you can launch bioGWAS using `docker-compose`. Example of configuration is in `docker-compose.yaml` file. You have to change (again, it shows how to mount directories inside container, you **have** to specify all directories you're using: the easiest way  is to make just one directory with subdirectories for I/O, and specify this one directory):

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
-dd "/data/gwassim_check/attempt_docker/data"
-imd "/data/gwassim_check/attempt_docker/images"
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


#### Genotypes files

There are 2 ways to set genotypes files:

__Method 1: using file__

* genotypes file(s) (in `.vcf` or PLINK's bfile format (the former is available only with skipping genotypes simulation)).
	* you need to create a list of all necessary file(s) to be included, and the chromosome corresponding to this file, e.g.[the example](data/chr.list):
        ```txt
        /wd/data/genotypes/file1.vcf,chr1
        /wd/data/genotypes/file2.vcf,chr2
        ...
        ```
        __Note__: path should be written as files will appear in docker container. All other files path will be specified later, when we launch docker.  
        Let this file be `./data/genotypes.list` in our example.

__Method 2: using flags__        

Instead of creating a file with genotypes and specifying it using `--input_list`,
we can run bioGWAS using `--geno_list` and  `--chr_list` and specify genotype file(s) and chromosome(s) respectively:
```bash
    ...
    --geno_list /wd/data/genotypes/file1.vcf,/wd/data/genotypes/file2.vcf
    --chr_list chr1,chr2
    ...
```

#### Samples to be included (optionally)

It should be txt-file with samples ids to be included in the analysis (one per line), [the example](data/EUR_SAMPLES_ID.txt):
    ```txt
    HG00096
    HG00097
    HG00099
    HG00100
    ...
    ```
    Let this file be `./data/samples.txt` in our example.

  In case you don't have this file, you can create it using `bcftools`: ```bcftools query -l your/vcf/file.vcf >  ./data/samles.txt```
  
#### Annotations

It should be annotation file in gtf format. In our test we used comprehensive gene annotation downloaded from [gencode site](https://www.gencodegenes.org/human/release_37.html). Let this file be `./data/gencodes.gtf` in our example.


#### Specifying causal SNPs/pathways

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

#### Summary

So, in total the necessary files are:
* `./data/genotypes/` - directory with genotypes files
* `./data/genotypes.list` - list with genotypes files
* `./data/samples.txt` - list with samples
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
--data_dir /wd \
--img_dir /wd  \
--input_list  /wd/data/genotypes.list \
--ids_file  /wd/data/samples.txt \
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
--data_dir /wd \
--img_dir /wd  \
--input_list  /wd/data/genotypes.list \
--ids_file  /wd/data/samles.txt \
--anno_file /wd/data/gencodes.gtf \
--gmt_file /wd/data/pathways.gmt  \
--use_causal_snps \
--causal_snps /wd/data/snps.txt
```


### 6. Wait!

After finishing, all output data will be in `./data` dir.
