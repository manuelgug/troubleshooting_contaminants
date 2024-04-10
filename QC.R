
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(stringr)


#INPUT: raw unfiltered allele_data.txt from mad4hatter


problem_run <- basename("runs/ASINT_NextSeq01_RESULTS/") #USER INPUT (folder name)
problem_run <- sub("_RESULTS.*", "", problem_run)


#### 1) Contextualize the run in terms of previous ones

# import allele_data.txt files
directories <- list.dirs("runs", full.names = TRUE, recursive = FALSE)

alelle_data_list <- list()

for (directory in directories) {

  folder_name <- basename(directory)
  
  allele_files <- list.files(directory, pattern = "allele_data.txt", full.names = TRUE)
  
  for (file in allele_files) {
    df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    alelle_data_list[[folder_name]] <- df
  }
}

# Clean the names of the data_frames_list
names_list <- names(alelle_data_list)
names_list_clean <- sub("_RESULTS.*", "", names_list)
names(alelle_data_list) <- names_list_clean

# Add a new column named "run" filled with the name of each run
for (i in seq_along(alelle_data_list)) {
  
  df_name <- names(alelle_data_list)[i]
  alelle_data_list[[i]]$run <- df_name
}

# Combine all data frames from allele_data_list into a single data frame
allele_data_df <- do.call(rbind, alelle_data_list)

# Create locus_allele column
allele_data_df$locus_allele <- paste0(allele_data_df$locus, "___", allele_data_df$pseudo_cigar)

# TEHERE IS A BUG WITH THE MASK THAT GENERATES MORE ALLELES THAN THERE REALLY ARE. THIS SNIPPET OF CODE COLLAPSES REPETITIONSAND SUMS THE READS AND FREQS. CRITICAL!!!
allele_data_df <- allele_data_df %>%
  group_by(run, sampleID, locus_allele) %>%
  summarize(reads = sum(reads))

# filter by a specific pool?
pool <- c("2") # USER INPUT: pool

if (pool == "2"){
  
    allele_data_df <- allele_data_df[grep("-2__", allele_data_df$locus_allele), ]

  }else if (pool == "1A"){
  
    allele_data_df <- allele_data_df[grep("-1A__", allele_data_df$locus_allele), ]

  }else if (pool == "1B"){
  
    allele_data_df <- allele_data_df[grep("-1B__", allele_data_df$locus_allele), ]
  
  }else if (pool == "1AB"){
    
    allele_data_df <- allele_data_df[grep("-1AB__", allele_data_df$locus_allele), ]
    
  }else if (pool == "1B2"){
    
    allele_data_df <- allele_data_df[grep("-1B2__", allele_data_df$locus_allele), ]
  }

# wide df format
allele_data_df <- pivot_wider(allele_data_df, names_from = locus_allele, values_from = reads)

# Replace NA with 0 in allele_data_df
allele_data_df[is.na(allele_data_df)] <- 0

allele_data_df_pres_abs <- allele_data_df
allele_data_df_pres_abs[allele_data_df_pres_abs > 1] <- 1



#### 2) Find outlier samples within a run

#subset problem run
problem_run_df <-  allele_data_df[allele_data_df$run == problem_run,]

#split counts and metadata, calculate allele proportions
problem_run_df_meta <- problem_run_df[,1:2]
problem_run_df_counts <- problem_run_df[,-1:-2]

#calculate allele proportions
problem_run_df_proprotions <- problem_run_df_counts/rowSums(problem_run_df_counts)
problem_run_df_proprotions[is.na(problem_run_df_proprotions)] <- 0 # change NaN for 0

#create presence/absence df
problem_run_df_pres_abs <- problem_run_df_counts
problem_run_df_pres_abs[problem_run_df_pres_abs > 1] <- 1


### pca ###
perform_pca <- function(problem_run_df_proportions, problem_run_df_meta, alpha = 0.01) {
  # Remove columns with constant or zero variance
  constant_cols <- sapply(problem_run_df_proportions, function(x) length(unique(x)) == 1)
  zero_variance_cols <- apply(problem_run_df_proportions, 2, function(x) var(x) == 0)
  cols_to_remove <- constant_cols | zero_variance_cols
  problem_run_df_proportions <- problem_run_df_proportions[, !cols_to_remove]
  
  # Perform PCA
  pca_result <- prcomp(problem_run_df_proportions, scale. = FALSE)
  
  # Extract PC scores
  pc_scores <- pca_result$x
  
  # Convert PC scores to dataframe and add SampleID column
  pc_scores_df <- as.data.frame(pc_scores)
  pc_scores_df$sampleID <- problem_run_df_meta$sampleID
  
  # Calculate Mahalanobis distance
  mah_dist <- mahalanobis(pc_scores_df[, c("PC1", "PC2")], colMeans(pc_scores_df[, c("PC1", "PC2")]), cov(pc_scores_df[, c("PC1", "PC2")]))
  mah_threshold <- qchisq(1 - alpha, df = 2)  # Chi-square threshold
  
  # Identify outliers
  pc_scores_df <- pc_scores_df %>%
    mutate(is_outlier = ifelse(mah_dist > mah_threshold, TRUE, FALSE))
  
  # Calculate percentage variance explained by each PC
  variance_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)
  
  # Plot PCA results with outliers colored in red and percentage variance explained
  pca_plot <- ggplot(pc_scores_df, aes(x = PC1, y = PC2, color = is_outlier, label = sampleID)) +
    geom_point(alpha = 0.5) +
    geom_text_repel(data = subset(pc_scores_df, is_outlier), aes(label = sampleID), size = 3, color = "red", segment.color = "red", segment.size = 0.5, max.overlaps = 1000) + 
    scale_color_manual(values = c("black", "red")) + 
    labs(title = "",
         x = paste0("PC1 (", variance_explained[1], "%)"),
         y = paste0("PC2 (", variance_explained[2], "%)")) +
    theme_minimal()
  
  return(pca_plot)
}

#perform pca analyses on allele proportions and presence/absence 
pca_plot_proportions <- perform_pca(problem_run_df_proprotions, problem_run_df_meta)
pca_plot_proportions <- pca_plot_proportions + labs(title = "Allele proportions")
pca_plot_pres_abs <- perform_pca(problem_run_df_pres_abs, problem_run_df_meta)
pca_plot_pres_abs <- pca_plot_pres_abs + labs(title = "Allele presence/absence")

pca_merged <- grid.arrange(pca_plot_proportions, pca_plot_pres_abs, ncol = 2)

ggsave(paste0(problem_run, "_samples_PCA.png"), pca_merged, width = 12, height = 6)


#### 3) Check presence of REF alleles from 3d7 controls

#subset problem run again because fuck efficiency
problem_run_df_again <- alelle_data_list[[problem_run]]

#extract 3D7 controls
controls_3D7 <- unique(problem_run_df_again[(grepl("(?i)3D7", problem_run_df_again$sampleID) & !grepl("(?i)(Dd2|PM|HB3|RDT)", problem_run_df_again$sampleID)),])

#turn masked pseudocigars into REF allele; those with N and without A, T, C or G into "."
controls_3D7$pseudo_cigar <- ifelse(grepl("N", controls_3D7$pseudo_cigar) & !grepl("[ATCG]", controls_3D7$pseudo_cigar), ".", controls_3D7$pseudo_cigar)

#coun ref alleles found
controls_3D7_ref_counts <- controls_3D7 %>%
  group_by(sampleID) %>%
  summarise(ref_alleles_found = sum(pseudo_cigar == "."))

#total ref alleles that SHOULD be found
total_amplicons <- length(unique(controls_3D7$locus))

#percentage ref alleles found
controls_3D7_ref_counts$perc_ref_alleles_found <- controls_3D7_ref_counts$ref_alleles_found / total_amplicons

write.csv(controls_3D7_ref_counts, paste0(problem_run, "_controls_3D7.csv"), row.names = F)



#### 4) Compare alleles from neg controls against all samples from previous and current run to detect potential sources of contaminants