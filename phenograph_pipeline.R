#!/usr/bin/env Rscript

###############################################
## Phenograph + UMAP Automated AWS Pipeline  ##
## Fully automated version for AWS Batch     ##
###############################################

suppressPackageStartupMessages({
  #library(aws.s3)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(Rtsne)
  library(igraph)
  library(FNN)
  library(lars)
  library(pheatmap)
  library(Rphenograph)
  library(inflection)
  library(uwot)
  library(flowCore)
  library(Biobase)
  library(Hmisc)
  library(magrittr)
  library(plyr)
})

options(scipen = 999)
set.seed(1)

###############################################
### 1. Read environment variables
###############################################

input_s3  <- Sys.getenv("input_s3")   # REQUIRED
output_s3 <- Sys.getenv("output_s3")  # REQUIRED
k_value   <- as.numeric(Sys.getenv("k", "60"))
events    <- as.numeric(Sys.getenv("events", "500"))

if (input_s3 == "" | output_s3 == "") {
  stop("ERROR: input_s3 and output_s3 must be provided as environment variables.")
}

message("Input S3: ", input_s3)
message("Output S3: ", output_s3)
message("k for Phenograph: ", k_value)
message("Events per file: ", events)

###############################################
### 2. Create local directories
###############################################

dir.create("InputFiles", showWarnings = FALSE)
dir.create("OutputFiles", showWarnings = FALSE)
dir.create("OutputFiles/Linear_MFIs", showWarnings = FALSE)
dir.create("OutputFiles/Linear_MFIs/linear_UMAP", showWarnings = FALSE)
dir.create("OutputFiles/BiExp_MFIs", showWarnings = FALSE)
dir.create("OutputFiles/BiExp_MFIs/hi_biex_transform_UMAP", showWarnings = FALSE)
dir.create("OutputFiles/BiExp_MFIs/mid_biex_transform_UMAP", showWarnings = FALSE)
dir.create("OutputFiles/BiExp_MFIs/lo_biex_transform_UMAP", showWarnings = FALSE)

###############################################
### 3. Download all CSVs from S3
###############################################

message("Downloading CSV files from S3…")

objs <- get_bucket_df(bucket = input_s3)

csv_files <- objs %>% filter(grepl("\\.csv$", Key)) %>% pull(Key)

if (length(csv_files) == 0)
  stop("ERROR: No CSV files found in the S3 input folder.")

for (f in csv_files) {
  save_object(
    object = f,
    bucket = input_s3,
    file = file.path("InputFiles", basename(f))
  )
  message("Downloaded: ", f)
}

###############################################
### 4. Load CSV files
###############################################

all_files <- list.files("InputFiles", pattern = "*.csv$", full.names = TRUE)

# Manifest file
manifest_file <- all_files[grepl("Manifest", all_files, ignore.case = TRUE)]

if (length(manifest_file) != 1) {
  stop("ERROR: Manifest.csv must exist exactly once in the input folder!")
}

manifest <- read.csv(manifest_file)

# Flow cytometry CSV files
fc_files <- all_files[!all_files %in% manifest_file]

if (length(fc_files) == 0)
  stop("ERROR: No flow cytometry CSV files found.")

###############################################
### 5. Read and downsample flow cytometry files
###############################################

flowLS <- lapply(fc_files, function(FILE) {
  
  fName <- gsub(".*_.*_(.*?_\\d+hrs)_.*", "\\1", basename(FILE))
  df <- read.csv(FILE) %>% as.data.frame()
  
  df_sample <- df[c(sample.int(min(events, nrow(df)))), ] %>%
    mutate(DonorID = fName)
  
  message("Loaded: ", basename(FILE), " (Donor: ", fName, ")")
  return(df_sample)
})

flowDF0 <- do.call(rbind, flowLS)

###############################################
### 6. Remove markers not used for clustering
###############################################

Columns_not_for_clustering <- c(
  "FSC.A", "SSC.A", "FSC.H", "SSC.H", "FSC.W", "SSC.W",
  "LD", "CD3"
)

flowDF <- do.call(rbind, flowLS) %>%
  select(-Columns_not_for_clustering, -Time, -DonorID)

###############################################
### 7. Scale markers & run UMAP
###############################################

flowDF_scaled <- scale(flowDF)
um <- umap(flowDF_scaled)

x1 <- as.data.frame(um)

pdf("OutputFiles/UMAP_All_events.pdf")
print(ggplot(x1, aes(V1, V2)) + geom_point(size = .1, alpha = .5))
dev.off()

###############################################
### 8. Phenograph clustering
###############################################

pg <- Rphenograph(flowDF_scaled, k = k_value)
clusters <- membership(pg[[2]])

x0 <- cbind(x1, flowDF0)
x0$louvain <- as.factor(clusters)

centroids <- x0 %>% group_by(louvain) %>% summarise(V1 = mean(V1), V2 = mean(V2))

pdf("OutputFiles/UMAP_All_events_clusters.pdf")
print(
  ggplot(x0, aes(V1, V2, color = louvain)) +
    geom_point(size = .1, alpha = .5) +
    geom_label_repel(data = centroids, aes(label = louvain)) +
    guides(color = FALSE)
)
dev.off()

###############################################
### 9. MFI per cluster + heatmaps
###############################################

x3 <- cbind(x0[,3:ncol(x0)], x0[,1:2])

mfi <- aggregate(x3[,1:(ncol(x3)-5)], list(x3$louvain), median)
rownames(mfi) <- mfi[,1]
mfi <- mfi[,-1]

write.csv(mfi, "OutputFiles/MFI_per_cluster.csv")

###############################################
### 10. Full folder of BiExp and Linear plots
###############################################

fcsData <- as.matrix(flowDF)

metaData <- data.frame(name = colnames(fcsData),
                       desc = paste("column", colnames(fcsData)))

ff <- new("flowFrame", exprs = fcsData, parameters = AnnotatedDataFrame(metaData))

# Three biexponential transforms (hi, mid, low)
biex_values <- list(
  high = c(a=5000,c=5000),
  mid  = c(a=50,c=50),
  low  = c(a=0.05,c=0.05)
)

for (level in names(biex_values)) {
  
  bpar <- biex_values[[level]]
  biex <- biexponentialTransform("BE", a=bpar["a"], b=1, c=bpar["c"],
                                 d=1, f=0, w=0)
  
  BET <- transform(ff, transformList(colnames(exprs(ff)), biex))
  df <- as.data.frame(exprs(BET))
  df2 <- cbind(df, x0[,c("V1","V2","louvain")])
  
  out_dir <- paste0("OutputFiles/BiExp_MFIs/", level, "_biex_transform_UMAP")
  
  for (i in 1:(ncol(df2)-3)) {
    
    pdf(sprintf("%s/MFI_UMAP_%s_%s.pdf",
                out_dir, level, colnames(df2)[i]),
        width = 6, height = 5)
    
    print(
      ggplot(df2, aes(V1, V2, col = df2[,i])) +
        geom_point(size = .1, alpha = .5) +
        scale_colour_gradient2(low="green", mid="yellow", high="red")
    )
    
    dev.off()
  }
}

###############################################
### 11. Cluster frequency tables
###############################################

x6 <- as.matrix(table(x3$DonorID, x3$louvain))
write.csv(x6, "OutputFiles/lovain_cluster_table_absolute_counts.csv")

props <- prop.table(x6, margin = 1)
write.csv(props, "OutputFiles/lovain_cluster_table_proportions.csv")

###############################################
### 12. Statistics + correlations
###############################################

show <- merge(manifest, as.data.frame(props),
              by.x="SCOPE_ID", by.y="Var1")

pdf("OutputFiles/ClusterFreq_group.pdf")
print(
  ggplot(show, aes(Var2, Freq, color = Group)) +
    geom_boxplot(outlier.shape = NA)
)
dev.off()

###############################################
### 13. Upload entire OutputFiles folder to S3
###############################################

message("Uploading results to S3…")

out_files <- list.files("OutputFiles", recursive = TRUE, full.names = TRUE)

for (f in out_files) {
  rel <- gsub("OutputFiles/", "", f)
  put_object(
    file = f,
    bucket = output_s3,
    object = paste0("OutputFiles/", rel)
  )
  message("Uploaded: ", rel)
}

message("===== PIPELINE COMPLETE =====")
