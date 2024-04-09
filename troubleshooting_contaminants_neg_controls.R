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

#keep only neg controls from asintmal
neg_controls_asntmal <- unique(A_p2$sampleID[grepl("NN",A_p2$sampleID)])

A_p2 <- A_p2[A_p2$sampleID %in% neg_controls_asntmal,]

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



