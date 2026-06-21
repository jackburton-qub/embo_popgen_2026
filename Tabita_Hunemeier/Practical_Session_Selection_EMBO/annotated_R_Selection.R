install.packages("rehh")
library("rehh")

EDAR_LWK <- data2haplohh(hap_file = "Chr2_EDAR_LWK_500K.recode.vcf", 
                         polarize_vcf = FALSE)
EDAR_CHS <- data2haplohh(hap_file = "Chr2_EDAR_CHS_500K.recode.vcf", 
                         polarize_vcf = FALSE)

# estimate EHH decay for rs3827760 in both populations

EDAR_LWK_ehh <- calc_ehh(EDAR_LWK, mrk = "rs3827760")
EDAR_CHS_ehh <- calc_ehh(EDAR_CHS, mrk = "rs3827760")

# Plot the EHH decay and furcation trees for both AFR and EAS

library(ggplot2)

# Plot EHH decay for LWK

plot(EDAR_LWK_ehh)

# Plot EHH decay for CHS

plot(EDAR_CHS_ehh)

# scan for genome wide homozygosity in both populations using scan_hh function

scan_LWK <- scan_hh(EDAR_LWK, polarized = FALSE)
scan_CHS <- scan_hh(EDAR_CHS, polarized = FALSE)

# get ihs scores, check score at rs3827760, and generate a single-site ihs plot

ihs_LWK <- ihh2ihs(scan_LWK)
ihs_CHS <- ihh2ihs(scan_CHS)

# check score at rs3827760
ihs_LWK[which(ihs_LWK$CHR = "rs3827760"), ]

# plot the at one position

plot(ihs_LWK$ihs$POSITION, ihs_LWK$ihs$IHS,
     col = ifelse(ihs_LWK$ihs$POSITION == "109513601", "red", "black"),
     pch = 19)
abline(h = 0, col = "blue", lty = 2)

plot(ihs_CHS$ihs$POSITION, ihs_CHS$ihs$IHS,
     col = ifelse(ihs_CHS$ihs$POSITION == "109513601", "red", "black"),
     pch = 19)
abline(h = 0, col = "blue", lty = 2)


# create a function to estimate the average absolute 
# iHS in sliding windows (50 SNPs/40 step) and plot the results

slide <- function(ihs_data, window_size = 50, step_size = 40) {
  positions <- ihs_data$ihs$POSITION
  ihs_values <- ihs_data$ihs$IHS
  
  # Create a data frame to store the results
  results <- data.frame(Window_Start = numeric(),
                        Window_End = numeric(),
                        Avg_Abs_IHS = numeric())
  
  # Loop through the positions in sliding windows
  for (start in seq(1, length(positions) - window_size + 1, by = step_size)) {
    end <- start + window_size - 1
    if (end > length(positions)) break
    
    window_positions <- positions[start:end]
    window_ihs_values <- ihs_values[start:end]
    
    avg_abs_ihs <- mean(abs(window_ihs_values), na.rm = TRUE)
    
    results <- rbind(results, data.frame(Window_Start = min(window_positions),
                                         Window_End = max(window_positions),
                                         Avg_Abs_IHS = avg_abs_ihs))
  }
  
  return(results)
}

# Apply the sliding window function to both populations

slide_LWK <- slide(ihs_LWK)
slide_CHS <- slide(ihs_CHS)

# Plot the sliding window results for both populations in separate plots

ggplot(slide_LWK, aes(x = Window_Start, y = Avg_Abs_IHS)) +
  geom_point(color = "blue") +
  labs(title = "Sliding Window Average Absolute iHS - LWK",
       x = "Genomic Position",
       y = "Average Absolute iHS") +
  theme_minimal()

ggplot(slide_CHS, aes(x = Window_Start, y = Avg_Abs_IHS)) +
  geom_point(color = "red") +
  labs(title = "Sliding Window Average Absolute iHS - CHS",
       x = "Genomic Position",
       y = "Average Absolute iHS") +
  theme_minimal()

# Estimate cross-population XP-EHH between LWK and CHS using ies2xpehh(), 
# calculate window-based averages, and plot them.

xp_ehh <- ies2xpehh(scan_LWK, scan_CHS)

window_xp_ehh <- slide(xp_ehh)

# Plot the sliding window results for XP-EHH

ggplot(window_xp_ehh, aes(x = Window_Start, y = Avg_Abs_IHS)) +
  geom_point(color = "purple") +
  labs(title = "Sliding Window Average Absolute XP-EHH",
       x = "Genomic Position",
       y = "Average Absolute XP-EHH") +
  theme_minimal()

#####################################################################

# Write the R code necessary to: 1. Read the pairwise $F_{ST}$ files involving Native Americans (input/Part_1_HumanDiversity/Chr2_NAM_EAS.weir.fst and input/Part_1_HumanDiversity/Chr2_NAM_EUR.weir.fst) 
# and Europeans-East Asians (input/Part_1_HumanDiversity/Chr2_EUR_EAS.weir.fst). 
# 2. Filter duplicates, exclude NA values, and align positions. 
# 3. Convert negative $F_{ST}$ values to zero. 
# 4. Estimate $PBS_{NAM}$ using NAM, EAS, and EUR. 
# 5. Check PBS value at rs3827760, check quantiles, and plot the PBS scan.

# first read FST files

NAM_EAS_fst <- read.table("Chr2_NAM_EAS.weir.fst", header = TRUE)
EUR_EAS_fst <- read.table("Chr2_EUR_EAS.weir.fst", header = TRUE)
NAM_EUR_fst <- read.table("Chr2_NAM_EUR.weir.fst", header = TRUE)

# filter duplicates and exclude NA FST values

NAM_EAS_filt <- NAM_EAS_fst[!duplicated(NAM_EAS_fst$POS) & !is.na(NAM_EAS_fst$WEIR_AND_COCKERHAM_FST), ]
EUR_EAS_filt <- EUR_EAS_fst[!duplicated(EUR_EAS_fst$POS) & !is.na(EUR_EAS_fst$WEIR_AND_COCKERHAM_FST), ]
NAM_EUR_filt <- NAM_EUR_fst[!duplicated(NAM_EUR_fst$POS) & !is.na(NAM_EUR_fst$WEIR_AND_COCKERHAM_FST), ]

NAM_EAS_filt$WEIR_AND_COCKERHAM_FST <- pmax(NAM_EAS_filt$WEIR_AND_COCKERHAM_FST, 0)
EUR_EAS_filt$WEIR_AND_COCKERHAM_FST <- pmax(EUR_EAS_filt$WEIR_AND_COCKERHAM_FST, 0)
NAM_EUR_filt$WEIR_AND_COCKERHAM_FST <- pmax(NAM_EUR_filt$WEIR_AND_COCKERHAM_FST, 0)

# align positions across the three datasets

aligned_positions <- Reduce(intersect, list(NAM_EAS_filt$POS, EUR_EAS_filt$POS, NAM_EUR_filt$POS))

# filter the datasets to keep only aligned positions

aligned_filtered <- function(df, aligned_positions) {
  df[df$POS %in% aligned_positions, ]
}

NAM_EAS_filt <- aligned_filtered(NAM_EAS_filt, aligned_positions)
EUR_EAS_filt <- aligned_filtered(EUR_EAS_filt, aligned_positions)
NAM_EUR_filt <- aligned_filtered(NAM_EUR_filt, aligned_positions)

# estimate PBS

T_nam_eas <- -log(1 - NAM_EAS_filt$WEIR_AND_COCKERHAM_FST)   # T(AFR, EAS)
T_nam_eur <- -log(1 - NAM_EUR_filt$WEIR_AND_COCKERHAM_FST)   # T(AFR, EUR)
T_eur_eas <- -log(1 - EUR_EAS_filt$WEIR_AND_COCKERHAM_FST)   # T(EAS, EUR)

PBS_EAS_raw <- ((T_nam_eas + T_eur_eas) - T_nam_eur) / 2

# ── 4. Convert negative branch lengths to zero ───────────────
PBS_EAS <- pmax(PBS_EAS_raw, 0)

# Assemble results data frame
pbs_df <- data.frame(
  CHROM   = NAM_EAS_filt$CHROM,
  POS     = NAM_EAS_filt$POS,
  PBS_EAS = PBS_EAS
)

cat("\nPBS_EAS summary:\n")
print(summary(pbs_df$PBS_EAS))

#  plot the PBS scan.

plot(pbs_df$POS, pbs_df$PBS_EAS, type = "l", col = "blue",
     xlab = "Genomic Position", ylab = "PBS (EAS)",
     main = "PBS Scan for EAS Population")

# restrict x-axis to the region around rs3827760 (109,513,601) using points instead of a line

library(ggplot2)

ggplot(pbs_df, aes(x = POS, y = PBS_EAS)) +
  geom_point(color = "blue") +
  geom_vline(xintercept = 109513601, color = "red", linetype = "dashed") +
  labs(title = "PBS Scan for EAS Population",
       x = "Genomic Position",
       y = "PBS (EAS)") +
  theme_minimal()


#####################################################################

# Write the R code necessary to: 
# 1. Load the polarized Small Dogs and Large Dogs VCFs into rehh using data2haplohh(). 
# 2. Compute the haplotype homozygosity scan for both populations using scan_hh(). 
# 3. Compute the cross-population Rsb statistic using ines2rsb(). 
# 4. Plot the Rsb Manhattan plot using ggplot2, highlighting the IGF1 region between 41 Mb and 45.5 Mb in red.

library("rehh")

# Load the polarized Small Dogs and Large Dogs VCFs into rehh

install.packages("R.utils")
library("R.utils")

small_dogs <- data2haplohh(hap_file = "small_dogs_polarized.vcf.gz", 
                           polarize_vcf = FALSE)
large_dogs <- data2haplohh(hap_file = "large_dogs_polarized.vcf.gz",
                           polarize_vcf = FALSE)

scan_small <- scan_hh(small_dogs, polarized = FALSE)
scan_large <- scan_hh(large_dogs, polarized = FALSE)

# calculate Rsb statistics


