
library(dplyr)

#import allele data from ASINTMAL and ICAE runs
ASINTMAL <- read.table("ASINT_NextSeq01_RESULTS/allele_data.txt", header = T) 
ICAE <- read.table("ICAE_NextSeq01_demux__RESULTS/allele_data.txt", header = T)
BOH <- read.table("BOH22_Nextseq01_RESULTS_v0.1.8/allele_data.txt", header = T) # good run, just for comparison
HFS1 <- read.table("HFS_NextSeq01_RESULTS_v0.1.8/allele_data.txt", header = T) # good old run, just for comparison
SMC2_1 <- read.table("SMC2_NextSeq01_211022_RESULTS_v0.1.8/allele_data.txt", header = T) # good old run, just for comparison
TES <-read.table("TES22_NextSeq01_RESULTS_v0.1.8/allele_data.txt", header = T) # good old run, just for comparison

# RENAME ALLELES
ASINTMAL$locus_allele <- paste0(ASINTMAL$locus, "___", ASINTMAL$pseudo_cigar)
ICAE$locus_allele <- paste0(ICAE$locus, "___", ICAE$pseudo_cigar)
BOH$locus_allele <- paste0(BOH$locus, "___", BOH$pseudo_cigar)
HFS1$locus_allele <- paste0(HFS1$locus, "___", HFS1$pseudo_cigar)
SMC2_1$locus_allele <- paste0(SMC2_1$locus, "___", SMC2_1$pseudo_cigar)
TES$locus_allele <- paste0(TES$locus, "___", TES$pseudo_cigar)

# extract loci from pool 2
A_p2<- ASINTMAL[grep("-2$", ASINTMAL$locus), ]
I_p2<- ICAE[grep("-2$", ICAE$locus), ]
B_p2 <- BOH[grep("-2$", BOH$locus), ]
H_p2 <- HFS1[grep("-2$", HFS1$locus), ]
S_p2 <- SMC2_1[grep("-2$", SMC2_1$locus), ]
T_p2 <- TES[grep("-2$", TES$locus), ]

# #no pool filtering when enabled
# A_p2 <- ASINTMAL
# I_p2 <- ICAE
# B_p2 <- BOH
# H_p2 <- HFS1

#remove 3d7 controls from ASINTMAL run because they DO have pool 2
A_p2 <- A_p2[!grepl("3D7", A_p2$sampleID, ignore.case = TRUE), ]
A_p2 <- A_p2[!grepl("dd2", A_p2$sampleID, ignore.case = TRUE), ]
# I_p2 <- I_p2[!grepl("3D7", I_p2$sampleID, ignore.case = TRUE), ]
# I_p2 <- I_p2[!grepl("dd2", I_p2$sampleID, ignore.case = TRUE), ]
# B_p2 <- B_p2[!grepl("3D7", B_p2$sampleID, ignore.case = TRUE), ]
# B_p2 <- B_p2[!grepl("dd2", B_p2$sampleID, ignore.case = TRUE), ]
# H_p2 <- H_p2[!grepl("3D7", H_p2$sampleID, ignore.case = TRUE), ]
# H_p2 <- H_p2[!grepl("dd2", H_p2$sampleID, ignore.case = TRUE), ]

#remove reference alleles (".") from runs
A_p2 <- A_p2[A_p2$pseudo_cigar != ".",]
I_p2 <- I_p2[I_p2$pseudo_cigar != ".",]
B_p2 <- B_p2[B_p2$pseudo_cigar != ".",]
H_p2 <- H_p2[A_p2$pseudo_cigar != ".",]
S_p2 <- S_p2[S_p2$pseudo_cigar != ".",]
T_p2 <- T_p2[T_p2$pseudo_cigar != ".",]

#how many alleles of pool 2 from ASINTMAL are in ICAE?
alleles_A <- length(unique(A_p2$locus_allele))
alleles_I <- length(unique(I_p2$locus_allele))
alleles_B <- length(unique(B_p2$locus_allele))
alleles_H <- length(unique(H_p2$locus_allele))
alleles_S <- length(unique(S_p2$locus_allele))
alleles_T <- length(unique(T_p2$locus_allele))

A_I_shared_alleles <- unique(A_p2$locus_allele[unique(A_p2$locus_allele) %in% unique(I_p2$locus_allele)])
A_B_shared_alleles <- unique(A_p2$locus_allele[unique(A_p2$locus_allele) %in% unique(B_p2$locus_allele)])
A_H_shared_alleles <- unique(A_p2$locus_allele[unique(A_p2$locus_allele) %in% unique(H_p2$locus_allele)])
A_S_shared_alleles <- unique(A_p2$locus_allele[unique(A_p2$locus_allele) %in% unique(S_p2$locus_allele)])
A_T_shared_alleles <- unique(A_p2$locus_allele[unique(A_p2$locus_allele) %in% unique(T_p2$locus_allele)])

n_alleles_A_in_I <- length(A_I_shared_alleles)
n_alleles_A_in_B <- length(A_B_shared_alleles)
n_alleles_A_in_H <- length(A_H_shared_alleles)
n_alleles_A_in_S <- length(A_S_shared_alleles)
n_alleles_A_in_T <- length(A_T_shared_alleles)

print(paste0(round(n_alleles_A_in_I/alleles_A *100, 2), "% of alleles in ASINTMAL run (n = ", n_alleles_A_in_I, ") are also in ICAE run"))
print(paste0(round(n_alleles_A_in_B/alleles_A *100, 2), "% of alleles in ASINTMAL run (n = ", n_alleles_A_in_B, ") are also in BOH run"))
print(paste0(round(n_alleles_A_in_H/alleles_A *100, 2), "% of alleles in ASINTMAL run (n = ", n_alleles_A_in_H, ") are also in HFS1 run"))
print(paste0(round(n_alleles_A_in_S/alleles_A *100, 2), "% of alleles in ASINTMAL run (n = ", n_alleles_A_in_S, ") are also in SMC2_1 run"))
print(paste0(round(n_alleles_A_in_T/alleles_A *100, 2), "% of alleles in ASINTMAL run (n = ", n_alleles_A_in_T, ") are also in TES run"))

#samples that share alleles between runs
unique(A_p2[A_p2$locus_allele %in% A_I_shared_alleles, ]$sampleID)
unique(A_p2[A_p2$locus_allele %in% A_B_shared_alleles, ]$sampleID)
unique(A_p2[A_p2$locus_allele %in% A_H_shared_alleles, ]$sampleID)
unique(A_p2[A_p2$locus_allele %in% A_S_shared_alleles, ]$sampleID)
unique(A_p2[A_p2$locus_allele %in% A_T_shared_alleles, ]$sampleID)


#how many/ which amplicons of pool 2 are shared between runs?
unique_amps_ICAE <- setdiff(unique(I_p2$locus), unique(A_p2$locus))
unique_amps_BOH <- setdiff(unique(B_p2$locus), unique(A_p2$locus))
unique_amps_H1 <- setdiff(unique(H_p2$locus), unique(A_p2$locus))
print(paste0(unique_amps_ICAE, " is present in ICAE but not in ASINTMAL run"))
print(paste0(unique_amps_BOH, " is present in BOH but not in ASINTMAL run"))
print(paste0(unique_amps_H1, " is present in HFS1 but not in ASINTMAL run"))




library(tidyr)
library(ggplot2)
library(ggrepel)

#### PCA ASINTMAL pool2
reformatted_df_A <- A_p2 %>%
  group_by(sampleID, locus_allele) %>%
  summarize(reads = sum(reads)) %>%
  pivot_wider(names_from = locus_allele, values_from = reads, values_fill = 0)

# Calculate proportions
row_sums <- rowSums(reformatted_df_A[,-1], na.rm = TRUE)
proportions_df_A <- reformatted_df_A[,-1] / row_sums

#### PCA ICAE pool2
reformatted_df_I <- I_p2 %>%
  group_by(sampleID, locus_allele) %>%
  summarize(reads = sum(reads)) %>%
  pivot_wider(names_from = locus_allele, values_from = reads, values_fill = 0)

#remove amplicon not present in ASINTMAL for a more fair comparison
reformatted_df_I <- reformatted_df_I[, !names(reformatted_df_I) %in% unique_amps_ICAE]

# Calculate proportions
row_sums <- rowSums(reformatted_df_I[,-1], na.rm = TRUE)
proportions_df_I <- reformatted_df_I[,-1] / row_sums


#### PCA BOH pool2
reformatted_df_B <- B_p2 %>%
  group_by(sampleID, locus_allele) %>%
  summarize(reads = sum(reads)) %>%
  pivot_wider(names_from = locus_allele, values_from = reads, values_fill = 0)

#remove amplicon not present in ASINTMAL for a more fair comparison
reformatted_df_B <- reformatted_df_B[, !names(reformatted_df_B) %in% unique_amps_ICAE]

# Calculate proportions
row_sums <- rowSums(reformatted_df_B[,-1], na.rm = TRUE)
proportions_df_B <- reformatted_df_B[,-1] / row_sums


#### PCA HFS1 pool2
reformatted_df_H <- H_p2 %>%
  group_by(sampleID, locus_allele) %>%
  summarize(reads = sum(reads)) %>%
  pivot_wider(names_from = locus_allele, values_from = reads, values_fill = 0)

#remove amplicon not present in ASINTMAL for a more fair comparison
reformatted_df_H <- reformatted_df_H[, !names(reformatted_df_H) %in% unique_amps_ICAE]

# Calculate proportions
row_sums <- rowSums(reformatted_df_H[,-1], na.rm = TRUE)
proportions_df_H <- reformatted_df_H[,-1] / row_sums


#### PCA SMC pool2
reformatted_df_S <- S_p2 %>%
  group_by(sampleID, locus_allele) %>%
  summarize(reads = sum(reads)) %>%
  pivot_wider(names_from = locus_allele, values_from = reads, values_fill = 0)

#remove amplicon not present in ASINTMAL for a more fair comparison
reformatted_df_S <- reformatted_df_S[, !names(reformatted_df_S) %in% unique_amps_ICAE]

# Calculate proportions
row_sums <- rowSums(reformatted_df_S[,-1], na.rm = TRUE)
proportions_df_S <- reformatted_df_S[,-1] / row_sums


#### PCA TES pool2
reformatted_df_T <- T_p2 %>%
  group_by(sampleID, locus_allele) %>%
  summarize(reads = sum(reads)) %>%
  pivot_wider(names_from = locus_allele, values_from = reads, values_fill = 0)

#remove amplicon not present in ASINTMAL for a more fair comparison
reformatted_df_T <- reformatted_df_T[, !names(reformatted_df_T) %in% unique_amps_ICAE]

# Calculate proportions
row_sums <- rowSums(reformatted_df_T[,-1], na.rm = TRUE)
proportions_df_T <- reformatted_df_T[,-1] / row_sums





##### JOINED PCA pool2

# Merge data frames by column names
proportions_df_A <- as.data.frame(cbind(sampleID = reformatted_df_A$sampleID, proportions_df_A))
proportions_df_A$run <- "ASINTMAL"
proportions_df_I <- as.data.frame(cbind(sampleID = reformatted_df_I$sampleID, proportions_df_I))
proportions_df_I$run <- "ICAE"
proportions_df_B <- as.data.frame(cbind(sampleID = reformatted_df_B$sampleID, proportions_df_B))
proportions_df_B$run <- "BOH"
proportions_df_H <- as.data.frame(cbind(sampleID = reformatted_df_H$sampleID, proportions_df_H))
proportions_df_H$run <- "HFS1"
proportions_df_S <- as.data.frame(cbind(sampleID = reformatted_df_S$sampleID, proportions_df_S))
proportions_df_S$run <- "SMC2_1"
proportions_df_T <- as.data.frame(cbind(sampleID = reformatted_df_T$sampleID, proportions_df_T))
proportions_df_T$run <- "TES"

dfs <- list(proportions_df_A, proportions_df_I, proportions_df_B, proportions_df_H, proportions_df_S, proportions_df_T)
proportions_joined <- Reduce(function(x, y) merge(x, y, by = intersect(names(x), names(y)), all = TRUE), dfs)

runs <- proportions_joined[,c("run")]
sample_names <- proportions_joined[,c("sampleID")]

# Remove sampleID and runs columns
proportions_joined <- proportions_joined[, !(names(proportions_joined) %in% c("sampleID", "run"))]

#replace NA with 0
proportions_joined[is.na(proportions_joined)] <- 0




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
  scale_color_manual(values = c("red", "purple", "limegreen", "orange", "blue", "pink")) +
  labs(title = "PCA Plot of non-ref Allele Proportions (Pool 2 only)",
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
  geom_point(alpha = 0.5, size =2) +
  #geom_text_repel(data = subset(pc_scores, outlier), aes(label = sampleID), size = 3, color = "grey60", segment.color = "grey60", segment.size = 0.5) +
  scale_color_manual(values = c("red", "purple", "limegreen", "orange", "blue", "pink")) +
  labs(title = "PCoA Plot of non-ref Allele Proportions (Pool 2 only)",
       x = paste0("PCo 1: ", variance_explained_axis1, "%\n"),
       y = paste0("PCo 2: ", variance_explained_axis2, "%")) +
  theme_minimal()




# ##### PCA VISUALIZATION WITHOUT OUTLIERS
# # Perform PCA
# pca_result <- prcomp(proportions_joined)
# pc_scores <- as.data.frame(pca_result$x)
# 
# # Calculate outlier samples
# mah_dist <- mahalanobis(pc_scores[, c("PC1", "PC2")], colMeans(pc_scores[, c("PC1", "PC2")]), cov(pc_scores[, c("PC1", "PC2")]))^2
# alpha <- 0.05
# n_components <- min(dim(pc_scores)[2], dim(pc_scores)[1])
# chi_sq_crit <- qchisq(1 - alpha, df = n_components)
# outliers <- mah_dist > chi_sq_crit
# 
# # Filter out outliers
# pc_scores <- pc_scores[!outliers, ]
# sample_names <- sample_names[!outliers]
# runs <- runs[!outliers]
# 
# # Add sample names and runs
# pc_scores$sampleID <- sample_names
# pc_scores$runs <- runs
# 
# # Plot PCA results without outliers labeled
# ggplot(pc_scores, aes(x = PC1, y = PC2, color = runs, label = sampleID)) +
#   geom_point(alpha = 0.5) +
#   scale_color_manual(values = c("red", "purple", "limegreen", "orange")) +
#   labs(title = "PCA Plot of non-ref Allele Proportions (Pool 2 only)",
#        x = paste0("Principal Component 1 (", round(100 * pca_result$sdev[1]^2 / sum(pca_result$sdev^2), 2), "%)"), 
#        y = paste0("Principal Component 2 (", round(100 * pca_result$sdev[2]^2 / sum(pca_result$sdev^2), 2), "%)")) +
#   theme_minimal()





sample_names_A <- proportions_df_A["sampleID"]
proportions_df_A_pca <- proportions_df_A[, !(names(proportions_df_A) %in% c("sampleID", "run"))]

# Compute Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(proportions_df_A_pca, method = "bray")

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
pc_scores <- cbind(pc_scores, sampleID = sample_names_A)

# # Plot PCoA
variance_explained <- round(pcoa_result$values / sum(pcoa_result$values) * 100, 2)
variance_explained_axis1 <- variance_explained$Eigenvalues[1]
variance_explained_axis2 <- variance_explained$Eigenvalues[2]

# Plot PCoA results with outliers labeled
ggplot(pc_scores, aes(x = Axis.1, y = Axis.2, color = outlier, label = sampleID)) +
  geom_point(alpha = 0.5, size =2, color = "black") +
  geom_text_repel(data = subset(pc_scores, outlier), aes(label = sampleID), size = 3, color = "grey60", segment.color = "grey60", segment.size = 0.5,  max.overlaps = 1000) +
  #scale_color_manual(values = c("red", "purple", "limegreen", "orange")) +
  labs(title = "PCoA Plot of non-ref Allele Proportions (Pool 2 only)",
       x = paste0("PCo 1: ", variance_explained_axis1, "%\n"),
       y = paste0("PCo 2: ", variance_explained_axis2, "%")) +
  theme_minimal()

