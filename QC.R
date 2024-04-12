
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(stringr)
library(vegan)
library(ape)


#INPUT: raw unfiltered allele_data.txt from mad4hatter


problem_run <- basename("runs/ASINT_NextSeq01_RESULTS/") #USER INPUT (folder name)
problem_run <- sub("_RESULTS.*", "", problem_run)


#### 1) Contextualize the run with previous ones

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
allele_data_df_wide <- pivot_wider(allele_data_df, names_from = locus_allele, values_from = reads)

# Replace NA with 0 in allele_data_df_wide
# allele_data_df_wide[is.na(allele_data_df_wide)] <- 0
allele_data_df_wide <- replace(allele_data_df_wide, is.na(allele_data_df_wide), 0)




## PcoA of allele pres/abs of all samples from all runs
allele_data_df_pres_abs <- allele_data_df_wide
allele_data_df_metadata <- allele_data_df_pres_abs[,1:2] #separate metadata from values
allele_data_df_pres_abs <- allele_data_df_pres_abs[,-1:-2]
# allele_data_df_pres_abs[allele_data_df_pres_abs > 1] <- 1 #convert data to presence/absence
allele_data_df_pres_abs <- replace(allele_data_df_pres_abs, allele_data_df_pres_abs > 1, 1) #convert data to presence/absence

# Compute Bray-Curtis dissimilarity matrix
#bray_curtis_dist <- vegdist(allele_data_df_pres_abs, method = "bray") #SLOW AF
bray_curtis_dist <- designdist(allele_data_df_pres_abs, method = "(A+B-2*J)/(A+B)", terms = "binary") # braycurtis distance formula (faster than vegdist(allele_data_df_pres_abs, method = "bray") #SLOW AF)

# Perform PCoA
pcoa_result <- pcoa(bray_curtis_dist)

# Extract PCoA scores
pc_scores <- as.data.frame(pcoa_result$vectors)
pc_scores <- cbind(pc_scores, allele_data_df_metadata)

# # Plot PCoA
variance_explained <- round(pcoa_result$values / sum(pcoa_result$values) * 100, 2)
variance_explained_axis1 <- variance_explained$Eigenvalues[1]
variance_explained_axis2 <- variance_explained$Eigenvalues[2]

set.seed(123)

# Select the unique runs
unique_runs <- unique(pc_scores$run)
n_runs <- length(unique_runs)

# Number of colors to generate
n_colors <- 50

# Generate random colors in hexadecimal format
random_colors <- sample(rgb(runif(n_colors), runif(n_colors), runif(n_colors)), n_colors)

# Repeat the sampled colors to match the number of runs
color_palette <- rep(random_colors, length.out = n_runs)

# Plot with the custom color palette
pc_scores$problem_run <- ifelse(pc_scores$run == problem_run, "problem_run", "other")

# Plot with different shapes for problem_run and others
pcoa_runs<- ggplot(pc_scores, aes(x = Axis.1, y = Axis.2, color = run, label = sampleID, shape = problem_run)) +
  geom_point(alpha = 0.5, size = 4) +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = c(problem_run = 15, other = 20)) +
  labs(title = "PCoA Plot of PRESENCE/ABSENCE of Alleles",
       x = paste0("PCo 1: ", variance_explained_axis1, "%\n"),
       y = paste0("PCo 2: ", variance_explained_axis2, "%")) +
  theme_minimal() +
  guides(color = guide_legend(ncol = 1))

ggsave(paste0(problem_run, "_runs_PCoA.png"), pcoa_runs, width = 18, height = 10, bg ="white")


#### 2) Find outlier samples within a run

#subset problem run
problem_run_df <-  allele_data_df_wide[allele_data_df_wide$run == problem_run,]

#split counts and metadata, calculate allele proportions
problem_run_df_meta <- problem_run_df[,1:2]
problem_run_df_counts <- problem_run_df[,-1:-2]

#calculate allele proportions
problem_run_df_proprotions <- problem_run_df_counts/rowSums(problem_run_df_counts)
problem_run_df_proprotions[is.na(problem_run_df_proprotions)] <- 0 # change NaN for 0

#create presence/absence df
problem_run_df_pres_abs <- problem_run_df_counts
problem_run_df_pres_abs[problem_run_df_pres_abs > 1] <- 1


### pca to spot outlier samples ###
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
#find neg control names in problem_run
problem_run_df_again_again <- allele_data_df[allele_data_df$run == problem_run, ]
neg_controls <- unique(problem_run_df_again_again$sampleID[grepl("(?i)Neg",problem_run_df_again_again$sampleID)])

#extract data from neg controls
problem_run_df_again_again_neg_controls <- problem_run_df_again_again[problem_run_df_again_again$sampleID %in% neg_controls,]

#count unique alleles on each neg control
n_alleles_in_neg_controls <- problem_run_df_again_again_neg_controls %>% 
  group_by(sampleID) %>%
  summarize(length(unique(locus_allele)))


# Create an empty list to store sorted_combined_df for each neg_control
sorted_combined_dfs <- list()

# Loop through each element in neg_controls_asntmal
for (neg_control in neg_controls) {
  # Subset problem_run_df_again_again_neg_controls for the current neg_control
  problem_run_df_again_again_neg_controls_subset <- problem_run_df_again_again_neg_controls[grepl(neg_control, problem_run_df_again_again_neg_controls$sampleID), ]
  
  l <- length(unique(problem_run_df_again_again_neg_controls_subset$locus_allele))
  
  # Initialize an empty list to store results for each run
  run_results <- list()
  
  # Loop through each unique run
  for (run in unique(allele_data_df$run)) {
    # Subset the data for the current run
    run_data <- subset(allele_data_df, run == run)
    
    # Perform the required operations for the current run
    result <- run_data %>% 
      group_by(sampleID) %>% 
      summarize(
        n_alleles_shared = sum(unique(locus_allele) %in% unique(problem_run_df_again_again_neg_controls_subset$locus_allele)),
        perc_alleles_shared = sum(unique(locus_allele) %in% unique(problem_run_df_again_again_neg_controls_subset$locus_allele)) / l,
        run = run
      ) %>%
      arrange(desc(perc_alleles_shared))
    
    # Store the result for the current run
    run_results[[run]] <- result
  }
  
  # Combine the results for all runs into a single dataframe
  combined_df <- do.call(rbind, run_results)
  
  # Rename columns and adjust percentages
  combined_df$perc_alleles_shared <- ifelse(combined_df$perc_alleles_shared > 1, 1, combined_df$perc_alleles_shared)
  colnames(combined_df)[2:3] <- paste0(colnames(combined_df)[2:3], "_", neg_control)
  
  # Add the result to the list with the name of the neg_control
  sorted_combined_dfs[[neg_control]] <- combined_df
}



#nrow(sorted_combined_dfs[[1]]) #holy shit
sorted_combined_dfs <- lapply(sorted_combined_dfs, distinct) #lazy cope, shouldn't generate that amount of repeated rows!!!



# Merge all data frames in sorted_combined_dfs by sampleID and run
merged_df <- Reduce(function(x, y) merge(x, y, by = c("sampleID", "run"), all = TRUE), sorted_combined_dfs)

# Extract columns percentage
selected_cols <- c("sampleID", "run", grep("perc", colnames(merged_df), value = TRUE))

# Subset the dataframe
selected_df_perc <- merged_df[, selected_cols]

selected_df_perc$mean_perc_shared_alleles<- rowMeans(selected_df_perc[,3:ncol(selected_df_perc)])

selected_df_perc <- selected_df_perc %>%
  arrange(desc(mean_perc_shared_alleles))

selected_df_perc[,c(1,2,ncol(selected_df_perc))]

# Remove duplicates from sampleID column
unique_sampleID <- unique(selected_df_perc$sampleID)

# Reorder levels based on mean_perc_shared_alleles
sorted_levels <- unique_sampleID[order(-selected_df_perc$mean_perc_shared_alleles[match(unique_sampleID, selected_df_perc$sampleID)])]

# Set sampleID as factor with reordered levels
selected_df_perc$sampleID <- factor(selected_df_perc$sampleID, levels = sorted_levels)

#remove undetermined since they don't provide any info on source of contamination
selected_df_perc <- selected_df_perc[selected_df_perc$sampleID !="Undetermined_S0_L001",]

write.csv(selected_df_perc, paste(problem_run, "_percentages_contaminants.csv"), row.names = F)


## TAILPLOT

#keep 100 more likely samples for easier visualization
selected_df_perc_subset <- selected_df_perc[1:100,]

#calcualte q99 for mean percentages:
q99 <- quantile(selected_df_perc$mean_perc_shared_alleles, 0.99)

set.seed(123)

# Select the unique runs
unique_runs <- unique(selected_df_perc_subset$run)
n_runs <- length(unique_runs)

# Generate random colors in hexadecimal format
color_palette <- rainbow(n_runs)

shared_neg <- ggplot(selected_df_perc_subset, aes(x = 1:length(sampleID), y = mean_perc_shared_alleles, color = run)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = q99, linetype = "dashed", color = "red") +  # Add horizontal line
  scale_color_manual(values = color_palette) +
  labs(title = "",
       x = "Sample ID",
       y = "Mean Percentage of Shared Alleles",
       color = "Run") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  guides(color = guide_legend(ncol = 1))

ggsave(paste0(problem_run, "_shared_alleles_w_neg_controls.png"), shared_neg, width = 14, height = 8, bg = "white")


## BARPLOT

# Extract runs with mean_perc_shared_alleles > q99
selected_df_perc_q99 <- selected_df_perc[selected_df_perc$mean_perc_shared_alleles > q99,]

# Create a bar plot of the frequency of each unique run
run_freq <- selected_df_perc_q99 %>%
  group_by(run) %>%
  summarise(count = n())

barplot <- ggplot(run_freq, aes(x = reorder(run, -count), y = count,)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "# of samples with mean_perc_shared_alleles > q99",
       x = "Run",
       y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(paste0(problem_run, "_contam_runs.png"), barplot, width = 8, height = 6, bg = "white")
