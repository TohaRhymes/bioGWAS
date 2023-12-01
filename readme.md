# GWAS data simulation

Human GWAS data simulator from Changalidis et al., 2023.

Examples of usage as well as experiments described in the article are presented in the `tests/` directory (with descriptions in [tests/readme.md](tests/readme.md)).


# Installation

To make this tool working, you need:
* `Python 3.9` (+ packages are listed in `requirements.txt`, installation: `pip install -r requirements.txt`)
* `R 4.1` (+ packages are listed in `requirements-R.txt`, quick installation: `install_r_reqs.R`)
* `PLINK`, `PLINK2`, `bedtools` and `HAPGEN2` be installed. Pathes to these tools have to be presented in `dependencies.yaml` (use `-d` parameter to specify). The example is in `dependencies.yaml` file.

# Launch

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


Full list of parameters with its description and default values can be found here:

```bash
./biogwas.py --help
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




