
library(dplyr)

#import allele data from ASINTMAL and ICAE runs
ASINTMAL <- read.table("ASINT_NextSeq01_RESULTS/allele_data.txt", header = T) 
ICAE <- read.table("ICAE_NextSeq01_demux__RESULTS/allele_data.txt", header = T)
BOH <- read.table("BOH22_Nextseq01_RESULTS_v0.1.8/allele_data.txt", header = T) # good run, just for comparison
HFS1 <- read.table("HFS_NextSeq01_RESULTS_v0.1.8//allele_data.txt", header = T) # good old run, just for comparison

# extract loci from pool 2
A_p2<- ASINTMAL[grep("-2$", ASINTMAL$locus), ]
unique(A_p2$locus)

I_p2<- ICAE[grep("-2$", ICAE$locus), ]
unique(I_p2$locus)

B_p2 <- BOH[grep("-2$", BOH$locus), ]
unique(B_p2$locus)

H_p2 <- HFS1[grep("-2$", HFS1$locus), ]
unique(H_p2$locus)


#remove 3d7 controls from ASINTMAL run because they DO have pool 2
A_p2_no3D7 <- A_p2[!grepl("3D7", A_p2$sampleID, ignore.case = TRUE), ]
A_p2_no3D7 <- A_p2[!grepl("dd2", A_p2$sampleID, ignore.case = TRUE), ]
I_p2 <- I_p2[!grepl("3D7", I_p2$sampleID, ignore.case = TRUE), ]
I_p2 <- I_p2[!grepl("dd2", I_p2$sampleID, ignore.case = TRUE), ]
B_p2 <- B_p2[!grepl("3D7", B_p2$sampleID, ignore.case = TRUE), ]
B_p2 <- B_p2[!grepl("dd2", B_p2$sampleID, ignore.case = TRUE), ]
H_p2 <- H_p2[!grepl("3D7", H_p2$sampleID, ignore.case = TRUE), ]
H_p2 <- H_p2[!grepl("dd2", H_p2$sampleID, ignore.case = TRUE), ]


#how many alleles of pool 2 from ASINTMAL are in ICAE?
alleles_A <- length(unique(A_p2_no3D7$pseudo_cigar))
alleles_I <- length(unique(I_p2$pseudo_cigar))
alleles_B <- length(unique(B_p2$pseudo_cigar))
alleles_H <- length(unique(H_p2$pseudo_cigar))

alleles_A_in_I <- length(unique(A_p2_no3D7$pseudo_cigar[unique(I_p2$pseudo_cigar) %in% unique(I_p2$pseudo_cigar)]))
alleles_A_in_B <- length(unique(A_p2_no3D7$pseudo_cigar[unique(B_p2$pseudo_cigar) %in% unique(B_p2$pseudo_cigar)]))
alleles_A_in_H <- length(unique(A_p2_no3D7$pseudo_cigar[unique(H_p2$pseudo_cigar) %in% unique(H_p2$pseudo_cigar)]))

print(paste0(alleles_A/alleles_A_in_I *100, "% of alleles in ASINTMAL run (excluding 3D7 controls) from pool 2 are also in ICAE run"))
print(paste0(alleles_A/alleles_A_in_B *100, "% of alleles in ASINTMAL run (excluding 3D7 controls) from pool 2 are also in BOH run"))
print(paste0(alleles_A/alleles_A_in_H *100, "% of alleles in ASINTMAL run (excluding 3D7 controls) from pool 2 are also in HFS1 run"))
print(paste0(round(alleles_A_in_I/alleles_I *100, 2), "% of alleles in ICAE from pool 2 are also in ASINTMAL run (excluding 3D7 controls)"))
print(paste0(round(alleles_A_in_B/alleles_B *100, 2), "% of alleles in BOH from pool 2 are also in ASINTMAL run (excluding 3D7 controls)"))
print(paste0(round(alleles_A_in_H/alleles_H *100, 2), "% of alleles in HFS1 from pool 2 are also in ASINTMAL run (excluding 3D7 controls)"))

#how many/ which amplicons of pool 2 are shared between runs?
unique_amps_ICAB <- setdiff(unique(I_p2$locus), unique(A_p2_no3D7$locus))
print(paste0(unique_amps_ICAB, " is present in ICAB but not in ASINTMAL run"))




library(tidyr)
library(ggplot2)
library(ggrepel)

#### PCA ASINTMAL pool2
reformatted_df_A <- A_p2_no3D7 %>%
  group_by(sampleID, locus) %>%
  summarize(reads = sum(reads)) %>%
  pivot_wider(names_from = locus, values_from = reads, values_fill = 0)

# Calculate proportions
row_sums <- rowSums(reformatted_df_A[,-1], na.rm = TRUE)
proportions_df_A <- reformatted_df_A[,-1] / row_sums

#### PCA ICAE pool2
reformatted_df_I <- I_p2 %>%
  group_by(sampleID, locus) %>%
  summarize(reads = sum(reads)) %>%
  pivot_wider(names_from = locus, values_from = reads, values_fill = 0)

#remove amplicon not present in ASINTMAL for a more fair comparison
reformatted_df_I <- reformatted_df_I[, !names(reformatted_df_I) %in% unique_amps_ICAB]

# Calculate proportions
row_sums <- rowSums(reformatted_df_I[,-1], na.rm = TRUE)
proportions_df_I <- reformatted_df_I[,-1] / row_sums


#### PCA BOH pool2
reformatted_df_B <- B_p2 %>%
  group_by(sampleID, locus) %>%
  summarize(reads = sum(reads)) %>%
  pivot_wider(names_from = locus, values_from = reads, values_fill = 0)

#remove amplicon not present in ASINTMAL for a more fair comparison
reformatted_df_B <- reformatted_df_B[, !names(reformatted_df_B) %in% unique_amps_ICAB]

# Calculate proportions
row_sums <- rowSums(reformatted_df_B[,-1], na.rm = TRUE)
proportions_df_B <- reformatted_df_B[,-1] / row_sums


#### PCA HFS1 pool2
reformatted_df_H <- H_p2 %>%
  group_by(sampleID, locus) %>%
  summarize(reads = sum(reads)) %>%
  pivot_wider(names_from = locus, values_from = reads, values_fill = 0)

#remove amplicon not present in ASINTMAL for a more fair comparison
reformatted_df_H <- reformatted_df_H[, !names(reformatted_df_H) %in% unique_amps_ICAB]

# Calculate proportions
row_sums <- rowSums(reformatted_df_H[,-1], na.rm = TRUE)
proportions_df_H <- reformatted_df_H[,-1] / row_sums





##### JOINED PCA pool2

# Merge data frames by column names
proportions_df_A <- as.data.frame(cbind(sampleID = reformatted_df_A$sampleID, proportions_df_A))
proportions_df_A$run <- "ASINTMAL"
proportions_df_I <- as.data.frame(cbind(sampleID = reformatted_df_I$sampleID, proportions_df_I))
proportions_df_I$run <- "ICAB"
proportions_df_B <- as.data.frame(cbind(sampleID = reformatted_df_B$sampleID, proportions_df_B))
proportions_df_B$run <- "BOH"
proportions_df_H <- as.data.frame(cbind(sampleID = reformatted_df_H$sampleID, proportions_df_H))
proportions_df_H$run <- "HFS1"

dfs <- list(proportions_df_A, proportions_df_I, proportions_df_B, proportions_df_H)
proportions_joined <- Reduce(function(x, y) merge(x, y, by = intersect(names(x), names(y)), all = TRUE), dfs)

runs <- proportions_joined[,c("run")]
sample_names <- proportions_joined[,c("sampleID")]

proportions_joined <- proportions_joined[,c(-1, -ncol(proportions_joined))]

# Perform PCA
pca_result <- prcomp(proportions_joined)
pc_scores <- as.data.frame(pca_result$x)

# Calculate outlier samples
mah_dist <- mahalanobis(pc_scores[, c("PC1", "PC2")], colMeans(pc_scores[, c("PC1", "PC2")]), cov(pc_scores[, c("PC1", "PC2")]))^2
alpha <- 0.05
n_components <- min(dim(pc_scores)[2], dim(pc_scores)[1])
chi_sq_crit <- qchisq(1 - alpha, df = n_components)
outliers <- mah_dist > chi_sq_crit

# Add outliers to PC scores dataframe
pc_scores$outlier <- outliers
pc_scores <- as.data.frame(cbind(pc_scores, sampleID = sample_names))
pc_scores$runs <- runs

# Plot PCA results with outliers labeled
ggplot(pc_scores, aes(x = PC1, y = PC2, color = runs, label = sampleID)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(data = subset(pc_scores, outlier), aes(label = sampleID), size = 3, color = "grey60", segment.color = "grey60", segment.size = 0.5) + 
  scale_color_manual(values = c("red", "purple", "limegreen", "orange")) +
  labs(title = "PCA Plot of Amplicon Proportions (Pool 2 only)",
       x = paste0("Principal Component 1 (", round(100 * pca_result$sdev[1]^2 / sum(pca_result$sdev^2), 2), "%)"), 
       y = paste0("Principal Component 2 (", round(100 * pca_result$sdev[2]^2 / sum(pca_result$sdev^2), 2), "%)")) +
  theme_minimal()



library(vegan)
library(ape)

# Compute Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(proportions_joined, method = "bray")

# Perform PCoA
pcoa_result <- pcoa(bray_curtis_dist)

# Extract PCoA scores
pc_scores <- as.data.frame(pcoa_result$vectors)

# Calculate outlier samples
mah_dist <- mahalanobis(pc_scores, colMeans(pc_scores), cov(pc_scores))^2
alpha <- 0.05
n_components <- min(dim(pc_scores)[2], dim(pc_scores)[1])
chi_sq_crit <- qchisq(1 - alpha, df = n_components)
outliers <- mah_dist > chi_sq_crit

# Add outliers to PC scores dataframe
pc_scores$outlier <- outliers
pc_scores <- cbind(pc_scores, sampleID = sample_names)

# # Plot PCoA
variance_explained <- round(pcoa_result$values / sum(pcoa_result$values) * 100, 2)
variance_explained_axis1 <- variance_explained$Eigenvalues[1]
variance_explained_axis2 <- variance_explained$Eigenvalues[2]

# Plot PCoA results with outliers labeled
ggplot(pc_scores, aes(x = Axis.1, y = Axis.2, color = runs, label = sampleID)) +
  geom_point(alpha = 0.5) +
  #geom_text_repel(data = subset(pc_scores, outlier), aes(label = sampleID), size = 3, color = "grey60", segment.color = "grey60", segment.size = 0.5) + 
  scale_color_manual(values = c("red", "purple", "limegreen", "orange")) +
  labs(title = "PCoA Plot of Amplicon Proportions (Pool 2 only)",
       x = paste0("PCo 1: ", variance_explained_axis1, "%\n"),
       y = paste0("PCo 2: ", variance_explained_axis2, "%")) +
  theme_minimal()

