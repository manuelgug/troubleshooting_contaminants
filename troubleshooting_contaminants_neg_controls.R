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
# S_p2 <- SMC2_1
# T_p2 <- TES

#remove reference alleles (".") from runs
A_p2 <- A_p2[A_p2$pseudo_cigar != ".",]
I_p2 <- I_p2[I_p2$pseudo_cigar != ".",]
B_p2 <- B_p2[B_p2$pseudo_cigar != ".",]
H_p2 <- H_p2[A_p2$pseudo_cigar != ".",]
S_p2 <- S_p2[S_p2$pseudo_cigar != ".",]
T_p2 <- T_p2[T_p2$pseudo_cigar != ".",]

#keep only neg controls from asintmal
neg_controls_asntmal <- unique(A_p2$sampleID[grepl("NN",A_p2$sampleID)])

A_p2 <- A_p2[A_p2$sampleID %in% neg_controls_asntmal,]

A_p2 %>% 
  group_by(sampleID) %>%
  summarize(length(unique(locus_allele)))

#how many alleles of pool 2 from ASINTMAL are in ICAE (and other runs)?
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

#samples from ASINTMAL that share alleles between runs
unique(A_p2[A_p2$locus_allele %in% A_I_shared_alleles, ]$sampleID)
unique(A_p2[A_p2$locus_allele %in% A_B_shared_alleles, ]$sampleID)
unique(A_p2[A_p2$locus_allele %in% A_H_shared_alleles, ]$sampleID)
unique(A_p2[A_p2$locus_allele %in% A_S_shared_alleles, ]$sampleID)
unique(A_p2[A_p2$locus_allele %in% A_T_shared_alleles, ]$sampleID)


## eprcentage of reads from neg controls shared with other samples
l <- length(unique(A_p2$locus_allele))
l

Y <- ASINTMAL %>% 
  group_by(sampleID) %>% 
  summarize(n_alleles_shared = sum(locus_allele %in% A_p2$locus_allele),
            perc_alleles_shared = sum(locus_allele %in% A_p2$locus_allele)/l,
            run = "ASINTMAL") %>%
  arrange(desc(perc_alleles_shared))

X <- ICAE %>% 
  group_by(sampleID) %>% 
  summarize(n_alleles_shared = sum(locus_allele %in% A_p2$locus_allele),
            perc_alleles_shared = sum(locus_allele %in% A_p2$locus_allele)/l,
            run = "ICAE") %>%
  arrange(desc(perc_alleles_shared))

W <- BOH %>% 
  group_by(sampleID) %>% 
  summarize(n_alleles_shared = sum(locus_allele %in% A_p2$locus_allele),
            perc_alleles_shared = sum(locus_allele %in% A_p2$locus_allele)/l,
            run = "BOH") %>%
  arrange(desc(perc_alleles_shared))

Z <- HFS1 %>% 
  group_by(sampleID) %>% 
  summarize(n_alleles_shared = sum(locus_allele %in% A_p2$locus_allele),
            perc_alleles_shared = sum(locus_allele %in% A_p2$locus_allele)/l,
            run = "HFS1") %>%
  arrange(desc(perc_alleles_shared))

Z1 <- SMC2_1 %>% 
  group_by(sampleID) %>% 
  summarize(n_alleles_shared = sum(locus_allele %in% A_p2$locus_allele),
            perc_alleles_shared = sum(locus_allele %in% A_p2$locus_allele)/l,
            run = "SMC2_1") %>%
  arrange(desc(perc_alleles_shared))

Z2 <- TES %>% 
  group_by(sampleID) %>% 
  summarize(n_alleles_shared = sum(locus_allele %in% A_p2$locus_allele),
            perc_alleles_shared = sum(locus_allele %in% A_p2$locus_allele)/l,
            run = "TES") %>%
  arrange(desc(perc_alleles_shared))

combined_df <- rbind(Y, X, W, Z, Z1, Z2)

sorted_combined_df <- combined_df %>%
  arrange(desc(perc_alleles_shared))

sorted_combined_df


#check is there is a single sample in ICAB and ASINTMAL that has all neg control alleles.
# Create an empty list to store sorted_combined_df for each neg_control
sorted_combined_dfs <- list()

# Loop through each element in neg_controls_asntmal
for (neg_control in neg_controls_asntmal) {
  # Subset A_p2 for the current neg_control
  A_p2_subset <- A_p2[grepl(neg_control, A_p2$sampleID), ]
  
  l <- length(unique(A_p2_subset$locus_allele))
  l
  
  Y <- ASINTMAL %>% 
    group_by(sampleID) %>% 
    summarize(n_alleles_shared = sum(locus_allele %in% A_p2_subset$locus_allele),
              perc_alleles_shared = sum(locus_allele %in% A_p2_subset$locus_allele)/l,
              run = "ASINTMAL") %>%
    arrange(desc(perc_alleles_shared))
  
  X <- ICAE %>% 
    group_by(sampleID) %>% 
    summarize(n_alleles_shared = sum(locus_allele %in% A_p2_subset$locus_allele),
              perc_alleles_shared = sum(locus_allele %in% A_p2_subset$locus_allele)/l,
              run = "ICAE") %>%
    arrange(desc(perc_alleles_shared))
  
  W <- BOH %>% 
    group_by(sampleID) %>% 
    summarize(n_alleles_shared = sum(locus_allele %in% A_p2_subset$locus_allele),
              perc_alleles_shared = sum(locus_allele %in% A_p2_subset$locus_allele)/l,
              run = "BOH") %>%
    arrange(desc(perc_alleles_shared))
  
  Z <- HFS1 %>% 
    group_by(sampleID) %>% 
    summarize(n_alleles_shared = sum(locus_allele %in% A_p2_subset$locus_allele),
              perc_alleles_shared = sum(locus_allele %in% A_p2_subset$locus_allele)/l,
              run = "HFS1") %>%
    arrange(desc(perc_alleles_shared))
  
  Z1 <- SMC2_1 %>% 
    group_by(sampleID) %>% 
    summarize(n_alleles_shared = sum(locus_allele %in% A_p2_subset$locus_allele),
              perc_alleles_shared = sum(locus_allele %in% A_p2_subset$locus_allele)/l,
              run = "SMC2_1") %>%
    arrange(desc(perc_alleles_shared))
  
  Z2 <- TES %>% 
    group_by(sampleID) %>% 
    summarize(n_alleles_shared = sum(locus_allele %in% A_p2_subset$locus_allele),
              perc_alleles_shared = sum(locus_allele %in% A_p2_subset$locus_allele)/l,
              run = "TES") %>%
    arrange(desc(perc_alleles_shared))
  
  combined_df <- rbind(Y, X, W, Z, Z1, Z2)
  
  sorted_combined_df <- combined_df %>%
    arrange(desc(perc_alleles_shared))
  
  sorted_combined_df$perc_alleles_shared <- ifelse(sorted_combined_df$perc_alleles_shared > 1, 1, sorted_combined_df$perc_alleles_shared) #don't now why i have more than 100% in some samples, just changing to 100
  
  colnames(sorted_combined_df)[2:3] <- paste0(colnames(sorted_combined_df)[2:3], "_", neg_control)
  
  # Add Y to the list with the name of the neg_control converted to character
  sorted_combined_dfs[[neg_control]] <- sorted_combined_df
}

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

ggplot(selected_df_perc, aes(x = 1:length(sampleID), y = mean_perc_shared_alleles, color = run)) +
  geom_point(size =3, alpha = 0.7) +
  labs(title = "Mean Percentage of Shared Alleles by Sample ID",
       x = "Sample ID",
       y = "Mean Percentage of Shared Alleles",
       color = "Run") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  # Rotate x-axis labels for better readability


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





##### JOINED PCoA pool2

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



library(vegan)
library(ape)

### PcoA presence/absence

proportions_joined_PA <- ifelse(proportions_joined == 0 , 0, 1)

# Compute Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(proportions_joined_PA, method = "bray")

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
  labs(title = "PCoA Plot of non-ref PRESENCE/ABSENCE of Alleles (Pool 2 only)",
       x = paste0("PCo 1: ", variance_explained_axis1, "%\n"),
       y = paste0("PCo 2: ", variance_explained_axis2, "%")) +
  theme_minimal()



