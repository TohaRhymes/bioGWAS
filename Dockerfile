# Use an official Python runtime as a parent image
FROM python:3.9.12-slim


# Install system dependencies for R
RUN apt-get update && apt-get install -y --no-install-recommends \
    gnupg \
    dirmngr \
    software-properties-common \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    r-base \
    r-base-dev \
    rename

# Download and install bedtools
RUN apt-get install -y --no-install-recommends wget \
    && wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools-2.31.0.tar.gz \
    && tar -zxvf bedtools-2.31.0.tar.gz \
    && cd bedtools2 \
    && make \
    && mv bin/* /usr/local/bin/

# Download and install HAPGEN
RUN wget https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/download/builds/x86_64/v2.2.0/hapgen2_x86_64.tar.gz \
    && tar -zxvf hapgen2_x86_64.tar.gz \
    && mv hapgen2 /usr/local/bin/

# Download and install PLINK 1.90b6.26
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231018.zip \
    && unzip plink_linux_x86_64_20231018.zip \
    && mv plink /usr/local/bin/

# Download and install PLINK 2.00a4.9LM
RUN wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20231030.zip \
    && unzip plink2_linux_x86_64_20231030.zip \
    && mv plink2 /usr/local/bin/

# Write paths to dependencies.yaml
RUN echo "plink: /usr/local/bin/plink" > dependencies.yaml \
    && echo "plink2: /usr/local/bin/plink2" >> dependencies.yaml \
    && echo "bedtools: /usr/local/bin/bedtools" >> dependencies.yaml \
    && echo "hapgen2: /usr/local/bin/hapgen2" >> dependencies.yaml

# Clone the desired GitHub repository
RUN apt-get install -y --no-install-recommends git \
    && git clone https://github.com/TohaRhymes/bioGWAS.git


# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r bioGWAS/requirements.txt

# Install R packages from requirements-R.txt
RUN Rscript bioGWAS/install_r_reqs.R

# Clean up APT when done
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Define the data directory as a volume
VOLUME /data

