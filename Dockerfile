############################################################
## Use the WORKING flowCore base image you built earlier  ##
############################################################
FROM flowcore-base AS builder

# Install remaining CRAN packages
RUN R -e "install.packages(c( \
    'aws.s3','RcppAnnoy','igraph','uwot','Hmisc','Rtsne', \
    'RColorBrewer','ggrepel','FNN','lars','pheatmap', \
    'magrittr','plyr','inflection','tibble','dplyr','tidyr' \
), repos='https://cloud.r-project.org/')"

# Install Rphenograph (GitHub)
RUN R -e "install.packages('devtools', repos='https://cloud.r-project.org/')"
RUN R -e "devtools::install_github('JinmiaoChenLab/Rphenograph', upgrade='never')"

# Copy pipeline
WORKDIR /pipeline
COPY phenograph_pipeline.R /pipeline/phenograph_pipeline.R
RUN chmod +x /pipeline/phenograph_pipeline.R

# Final verification
RUN R -e "library(flowCore); library(Rphenograph); library(uwot); \
          cat('✓ flowCore + Rphenograph + uwot OK\n')"

ENTRYPOINT ["Rscript", "/pipeline/phenograph_pipeline.R"]
