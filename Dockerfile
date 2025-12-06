FROM flowcore-base:latest AS builder

RUN R -e "install.packages(c(
    'aws.s3','RcppAnnoy','igraph','uwot','Hmisc','Rtsne',
    'RColorBrewer','ggrepel','FNN','lars','pheatmap',
    'magrittr','plyr','inflection','tibble','dplyr','tidyr'
), repos='https://cloud.r-project.org/')"

RUN R -e "devtools::install_github('JinmiaoChenLab/Rphenograph', upgrade='never')"

WORKDIR /pipeline
COPY phenograph_pipeline.R /pipeline/phenograph_pipeline.R
RUN chmod +x /pipeline/phenograph_pipeline.R

ENTRYPOINT ["Rscript", "/pipeline/phenograph_pipeline.R"]
