#!/bin/bash

# download gzip-file
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
gzip -d gencode.v37.annotation.gtf.gz 

mv gencode.v37.annotation.gtf data

