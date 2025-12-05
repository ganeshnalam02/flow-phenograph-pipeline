############################################################
# OPTION 2: Build directly from Bioconductor base image
# This removes the need for flowcore-base
############################################################
FROM bioconductor/bioconductor_docker:RELEASE_3_17 AS builder

# Install system libs needed for R packages
RUN apt-get update && apt-get install -y \
    libssl-dev libcurl4-openssl-dev libxml2-dev libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

############################################################
# Install required R packages
############################################################
RUN R -e "install.packages(c( \
    'aws.s3','RcppAnnoy','igraph','uwot','Hmisc','Rtsne', \
    'RColorBrewer','ggrepel','FNN','lars','pheatmap', \
    'magrittr','plyr','inflection','tibble','dplyr','tidyr' \
), repos='https://cloud.r-project.org', dependencies=TRUE)"

# Install devtools + Rphenograph
RUN R -e "install.packages('devtools', repos='https://cloud.r-project.org')"
RUN R -e "devtools::install_github('JinmiaoChenLab/Rphenograph', upgrade='never')"

############################################################
# Add pipeline script
############################################################
WORKDIR /pipeline
COPY phenograph_pipeline.R /pipeline/phenograph_pipeline.R
RUN chmod +x /pipeline/phenograph_pipeline.R

############################################################
# Run test inside image build to confirm required libs load
############################################################
RUN R -e "library(flowCore); library(Rphenograph); library(uwot); cat('✓ Image build successful.\n')"

ENTRYPOINT ["Rscript", "/pipeline/phenograph_pipeline.R"]
