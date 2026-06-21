# EMBO Practical Course: Genomic Diversity & Natural Selection Scan (Worksheet)

This worksheet guides you through the analysis of genomic data in humans and canines to identify genomic regions under natural selection. The tutorial is divided into two main parts:

1.  **Part 1: Human Genomic Diversity and Natural Selection**: Investigating selection signatures in the candidate gene ***EDAR*** (associated with ectodermal traits in East Asians and Native Americans) using population differentiation ($F_{ST}$, PBS) and haplotype-based metrics (EHH, iHS, XP-EHH).
2.  **Part 2: Genomic Selection Scan in Canines**: Identifying the selective sweep at the ***IGF1*** body-size locus by comparing small vs. large dog breeds using PCA, PCAdapt (outlier scan), and haplotype homozygosity methods (XP-nSL and Rsb).

------------------------------------------------------------------------

# Part 1: Human Genomic Diversity and Natural Selection

## 1. Background and Dataset

### Goal

Our goal is to explore approaches and methods which seek to identify regions of the genome with signatures of natural selection. We will use real genomic data and two classes of tests: one based on population differentiation ($F_{ST}$ / PBS) and another based on extended haplotype homozygosity (EHH / iHS / XP-EHH).

### Dataset

Whole-genome sequencing data from the 1000 Genomes Project Phase III. The full database can be accessed via: <ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>

### Data Pre-processing

We will analyze a pre-processed dataset for chromosome 2 corresponding to individuals sampled from the African (AFR: 504 individuals), European (EUR: 503 individuals), and East Asian (EAS: 504 individuals) populations. In this dataset, INDELs, singletons, and SNPs with MAF \< 0.05 have been removed. The pairwise $F_{ST}$ was then estimated using `vcftools`.

All data files are located in the `input/` directory: - `input/Part_1_HumanDiversity/AFR_EAS.weir.fst` (Fst between Africans and East Asians) - `input/Part_1_HumanDiversity/AFR_EUR.weir.fst` (Fst between Africans and Europeans) - `input/Part_1_HumanDiversity/EAS_EUR.weir.fst` (Fst between East Asians and Europeans) - `input/Part_1_HumanDiversity/Chr2_EDAR_LWK_500K.recode.vcf` (Phased African haplotypes around *EDAR*) - `input/Part_1_HumanDiversity/Chr2_EDAR_CHS_500K.recode.vcf` (Phased East Asian haplotypes around *EDAR*) - `input/Part_1_HumanDiversity/Chr2_NAM_EAS.weir.fst` (Fst between Native Americans and East Asians in candidate region) - `input/Part_1_HumanDiversity/Chr2_NAM_EUR.weir.fst` (Fst between Native Americans and Europeans in candidate region) - `input/Part_1_HumanDiversity/Chr2_EUR_EAS.weir.fst` (Fst between Europeans and East Asians in candidate region)

------------------------------------------------------------------------

## 2. Genetic Differentiation (FST and PBS)

### Investigating the Candidate Gene *EDAR*

The human Ectodysplasin A receptor gene, or ***EDAR***, is part of the EDA signaling pathway which specifies prenatally the location, size, and shape of ectodermal appendages (such as hair follicles, teeth, and glands). *EDAR* is a textbook example of positive selection in East Asians. A specific non-synonymous variant, **rs3827760** (chr2:109,513,601 A\>G), results in a Val370Ala substitution and is strongly associated with thicker hair shafts and shovel-shaped incisors. Another hypothesis states that *EDAR* acted along with *FADS* and *VDR* in the Beringia Standstill, allowing Native American ancestors to survive in extreme arctic environments.

### Questions for Students

1.  **The estimate of** $F_{ST}$ by the Weir and Cockerham metric can sometimes generate negative values and "NA". What does that mean? How can this interfere with the results?
    -   *Answer*: Missing data values
2.  **The** $F_{ST}$ values observed between pairs of populations for the SNP rs3827760 (position 109,513,601) fall within which distribution quantiles of $F_{ST}$ values for the studied chromosome? Can they be considered outliers?
    -   *Answer*: yes
3.  **From the observed** $F_{ST}$ values between population pairs and the significance estimates, what can we say about the rs3827760 SNP differentiation between populations?
    -   *Answer*:
4.  **Discuss how these results justify performing another type of analysis based on PBS (Population Branch Statistics).**
    -   *Answer*:
5.  **What does the PBS analysis reveal? What is the difference between PBS and** $F_{ST}$ analysis?
    -   *Answer*:

------------------------------------------------------------------------

### R Code Exercise: Pairwise $F_{ST}$ Calculation

Write the R code necessary to perform the following: 1. Read the pairwise $F_{ST}$ files from `input/`. 2. Filter duplicate SNP positions and exclude NA values. 3. Align the datasets by overlapping positions. 4. Set negative $F_{ST}$ values to zero. 5. Check Fst values at position `109513601`. 6. Calculate distribution quantiles to determine if rs3827760 is an outlier. 7. Plot pairwise Fst around `109513601` in a 10kb window, highlighting the candidate SNP.

**Write your R code here:**

``` r
#   write code to read Fst files from input/. folder

infile_1 <- read.table(file = "EAS_EUR.weir.fst", header=T)
infile_2 <- read.table(file = "AFR_EUR.weir.fst", header=T)
infile_3 <- read.table(file = "EAS_EUR.weir.fst", header=T)



# ============================================================
# Pairwise FST Analysis – rs3827760 (chr2:109513601)
# Comparisons: AFR_EUR | AFR_EAS | EAS_EUR
# ============================================================

# ── 0. Libraries ────────────────────────────────────────────
library(ggplot2)
library(dplyr)

# ── 1. Read input files ──────────────────────────────────────
infile_1 <- read.table(file = "AFR_EAS.weir.fst", header = TRUE)
infile_2 <- read.table(file = "AFR_EUR.weir.fst", header = TRUE)
infile_3 <- read.table(file = "EAS_EUR.weir.fst", header = TRUE)

cat("Rows read — AFR_EAS:", nrow(infile_1),
    "| AFR_EUR:", nrow(infile_2),
    "| EAS_EUR:", nrow(infile_3), "\n")

# ── 2. Filter duplicate positions & remove NAs ───────────────
clean_fst <- function(df, label) {
  df <- df[!duplicated(df$POS), ]          # drop duplicate SNP positions
  df <- df[!is.na(df$WEIR_AND_COCKERHAM_FST), ]  # exclude NA FST values
  cat("After filtering [", label, "] :", nrow(df), "SNPs\n")
  df
}

fst1 <- clean_fst(infile_1, "AFR_EAS")
fst2 <- clean_fst(infile_2, "AFR_EUR")
fst3 <- clean_fst(infile_3, "EAS_EUR")

# ── 3. Align datasets by overlapping positions ───────────────
shared_pos <- Reduce(intersect, list(fst1$POS, fst2$POS, fst3$POS))
cat("Overlapping positions across all three comparisons:", length(shared_pos), "\n")

fst1 <- fst1[fst1$POS %in% shared_pos, ]
fst2 <- fst2[fst2$POS %in% shared_pos, ]
fst3 <- fst3[fst3$POS %in% shared_pos, ]

# Ensure consistent ordering by position
fst1 <- fst1[order(fst1$POS), ]
fst2 <- fst2[order(fst2$POS), ]
fst3 <- fst3[order(fst3$POS), ]

# ── 4. Set negative FST values to zero ──────────────────────
fst1$WEIR_AND_COCKERHAM_FST <- pmax(fst1$WEIR_AND_COCKERHAM_FST, 0)
fst2$WEIR_AND_COCKERHAM_FST <- pmax(fst2$WEIR_AND_COCKERHAM_FST, 0)
fst3$WEIR_AND_COCKERHAM_FST <- pmax(fst3$WEIR_AND_COCKERHAM_FST, 0)

# ── 5. Check FST values at rs3827760 (position 109513601) ────
target_pos <- 109513601

check_pos <- function(df, label) {
  row <- df[df$POS == target_pos, ]
  if (nrow(row) == 0) {
    cat(label, ": position", target_pos, "NOT FOUND\n")
  } else {
    cat(label, "FST at", target_pos, "=",
        round(row$WEIR_AND_COCKERHAM_FST, 4), "\n")
  }
}

cat("\n── FST at rs3827760 (pos", target_pos, ") ──\n")
check_pos(fst1, "AFR_EAS")
check_pos(fst2, "AFR_EUR")
check_pos(fst3, "EAS_EUR")

# ── 6. Quantile analysis – is rs3827760 an outlier? ─────────
cat("\n── Quantile distributions ──\n")

quantile_summary <- function(df, label) {
  fst_val <- df$WEIR_AND_COCKERHAM_FST
  q        <- quantile(fst_val, probs = c(0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.999))
  target   <- df$WEIR_AND_COCKERHAM_FST[df$POS == target_pos]
  pctile   <- if (length(target) > 0)
                  round(mean(fst_val <= target) * 100, 2)
              else NA

  cat("\n[", label, "]\n")
  print(round(q, 4))
  if (!is.na(pctile))
    cat("  rs3827760 FST =", round(target, 4),
        "→ percentile:", pctile, "%",
        if (pctile >= 99) "*** OUTLIER ***" else "", "\n")
}

quantile_summary(fst1, "AFR_EAS")
quantile_summary(fst2, "AFR_EUR")
quantile_summary(fst3, "EAS_EUR")

# ── 7. Plot pairwise FST in ±10 kb window around 109513601 ──
window_size <- 10000
win_min     <- target_pos - window_size
win_max     <- target_pos + window_size

# Combine all three comparisons into one tidy data frame
make_window <- function(df, label) {
  sub <- df[df$POS >= win_min & df$POS <= win_max, ]
  sub$comparison <- label
  sub
}

plot_df <- bind_rows(
  make_window(fst1, "AFR vs EAS"),
  make_window(fst2, "AFR vs EUR"),
  make_window(fst3, "EAS vs EUR")
)

# Candidate SNP annotation layer
candidate <- data.frame(
  POS        = target_pos,
  comparison = c("AFR vs EAS", "AFR vs EUR", "EAS vs EUR")
)
candidate <- left_join(
  candidate,
  bind_rows(
    data.frame(POS = target_pos, comparison = "AFR vs EAS",
               WEIR_AND_COCKERHAM_FST = fst1$WEIR_AND_COCKERHAM_FST[fst1$POS == target_pos]),
    data.frame(POS = target_pos, comparison = "AFR vs EUR",
               WEIR_AND_COCKERHAM_FST = fst2$WEIR_AND_COCKERHAM_FST[fst2$POS == target_pos]),
    data.frame(POS = target_pos, comparison = "EAS vs EUR",
               WEIR_AND_COCKERHAM_FST = fst3$WEIR_AND_COCKERHAM_FST[fst3$POS == target_pos])
  ),
  by = c("POS", "comparison")
)

p <- ggplot(plot_df, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, colour = comparison)) +

  # Background shading for candidate region
  annotate("rect",
           xmin = target_pos - 500, xmax = target_pos + 500,
           ymin = -Inf, ymax = Inf,
           fill = "gold", alpha = 0.25) +

  # Vertical line at candidate SNP
  geom_vline(xintercept = target_pos,
             linetype = "dashed", colour = "grey40", linewidth = 0.6) +

  # FST values for all SNPs in window
  geom_point(size = 1.8, alpha = 0.7) +
  geom_line(alpha = 0.4, linewidth = 0.5) +

  # Highlight candidate SNP
  geom_point(data = candidate,
             aes(x = POS, y = WEIR_AND_COCKERHAM_FST),
             colour = "red", size = 4, shape = 18) +

  # Label
  geom_label(data = candidate,
             aes(x = POS, y = WEIR_AND_COCKERHAM_FST,
                 label = "rs3827760"),
             nudge_y = 0.03, size = 3, colour = "red",
             fill = "white", fontface = "bold") +

  facet_wrap(~comparison, ncol = 1) +
  scale_colour_manual(values = c("AFR vs EAS" = "#E64B35",
                                  "AFR vs EUR" = "#4DBBD5",
                                  "EAS vs EUR" = "#00A087")) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title    = expression(paste("Pairwise ", F[ST], " around rs3827760 (chr2:109,513,601)")),
    subtitle = paste0("±", window_size / 1000, " kb window | Candidate SNP highlighted in red"),
    x        = "Chromosomal position (bp)",
    y        = expression(F[ST]),
    colour   = "Comparison"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey40", size = 10),
    strip.text    = element_text(face = "bold"),
    legend.position = "none",         
    panel.grid.minor = element_blank()
  )

print(p)
```

------------------------------------------------------------------------

### R Code Exercise: Population Branch Statistics (PBS)

Write the R code necessary to: 1. Estimate the Population Branch Statistic for East Asians ($PBS_{EAS}$) using the AFR, EAS, and EUR populations. 2. Convert negative branch lengths to zero. 3. Check the PBS value for the candidate SNP rs3827760. 4. Calculate distribution quantiles to determine if it is an outlier. 5. Plot PBS values around the candidate SNP in a 10kb window.

**Write your R code here:**

``` r
# # ============================================================
# Population Branch Statistic – PBS_EAS
# Populations: AFR | EAS | EUR
# Candidate SNP: rs3827760 (chr2:109,513,601)
# ============================================================
# PBS is computed from three pairwise FST values:
#   T(x,y) = -log(1 - FST(x,y))   [branch-length transformation]
#   PBS_EAS = ( T(EAS,AFR) + T(EAS,EUR) - T(AFR,EUR) ) / 2
# ============================================================

# ── 0. Libraries ─────────────────────────────────────────────
library(ggplot2)
library(dplyr)

# ── 1. Read & clean pairwise FST files ───────────────────────
# (mirrors pairwise_fst.R; repeated here so the script is self-contained)

read_fst <- function(path, label) {
  df <- read.table(file = path, header = TRUE)
  df <- df[!duplicated(df$POS), ]
  df <- df[!is.na(df$WEIR_AND_COCKERHAM_FST), ]
  # Floor FST at 0 before the log transform (avoids complex numbers)
  df$WEIR_AND_COCKERHAM_FST <- pmax(df$WEIR_AND_COCKERHAM_FST, 0)
  # Cap FST just below 1 to avoid log(0) = -Inf
  df$WEIR_AND_COCKERHAM_FST <- pmin(df$WEIR_AND_COCKERHAM_FST, 0.9999)
  cat("Loaded [", label, "] :", nrow(df), "SNPs\n")
  df
}

afr_eas <- read_fst("AFR_EAS.weir.fst", "AFR_EAS")
afr_eur <- read_fst("AFR_EUR.weir.fst", "AFR_EUR")
eas_eur <- read_fst("EAS_EUR.weir.fst", "EAS_EUR")

# ── 2. Align to overlapping positions ────────────────────────
shared_pos <- Reduce(intersect, list(afr_eas$POS, afr_eur$POS, eas_eur$POS))
cat("Overlapping positions:", length(shared_pos), "\n")

afr_eas <- afr_eas[afr_eas$POS %in% shared_pos, ]
afr_eur <- afr_eur[afr_eur$POS %in% shared_pos, ]
eas_eur <- eas_eur[eas_eur$POS %in% shared_pos, ]

afr_eas <- afr_eas[order(afr_eas$POS), ]
afr_eur <- afr_eur[order(afr_eur$POS), ]
eas_eur <- eas_eur[order(eas_eur$POS), ]

# ── 3. Estimate PBS_EAS ───────────────────────────────────────
# Branch-length transformation: T = -log(1 - FST)
T_afr_eas <- -log(1 - afr_eas$WEIR_AND_COCKERHAM_FST)   # T(AFR, EAS)
T_afr_eur <- -log(1 - afr_eur$WEIR_AND_COCKERHAM_FST)   # T(AFR, EUR)
T_eas_eur <- -log(1 - eas_eur$WEIR_AND_COCKERHAM_FST)   # T(EAS, EUR)

PBS_EAS_raw <- (T_afr_eas + T_eas_eur - T_afr_eur) / 2

# ── 4. Convert negative branch lengths to zero ───────────────
PBS_EAS <- pmax(PBS_EAS_raw, 0)

# Assemble results data frame
pbs_df <- data.frame(
  CHROM   = afr_eas$CHROM,
  POS     = afr_eas$POS,
  PBS_EAS = PBS_EAS
)

cat("\nPBS_EAS summary:\n")
print(summary(pbs_df$PBS_EAS))

# ── 5. Check PBS value at rs3827760 (pos 109513601) ──────────
target_pos <- 109513601

cat("\n── PBS_EAS at rs3827760 (pos", target_pos, ") ──\n")
candidate_row <- pbs_df[pbs_df$POS == target_pos, ]

if (nrow(candidate_row) == 0) {
  cat("Position", target_pos, "NOT FOUND in aligned dataset.\n")
} else {
  cat("PBS_EAS =", round(candidate_row$PBS_EAS, 6), "\n")
}

# ── 6. Quantile analysis – outlier detection ─────────────────
cat("\n── PBS_EAS quantile distribution ──\n")

q_probs <- c(0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.999)
q_vals  <- quantile(pbs_df$PBS_EAS, probs = q_probs)
print(round(q_vals, 6))

if (nrow(candidate_row) > 0) {
  pctile <- round(mean(pbs_df$PBS_EAS <= candidate_row$PBS_EAS) * 100, 3)
  cat("\nrs3827760 PBS_EAS =", round(candidate_row$PBS_EAS, 6),
      "→ percentile:", pctile, "%")
  if (pctile >= 99.9) {
    cat("  *** TOP 0.1% OUTLIER ***\n")
  } else if (pctile >= 99) {
    cat("  *** TOP 1% OUTLIER ***\n")
  } else if (pctile >= 95) {
    cat("  * Top 5%\n")
  } else {
    cat("  (not a strong outlier)\n")
  }
}

# ── 7. Plot PBS_EAS in ±10 kb window around rs3827760 ────────
window_size <- 10000
win_min     <- target_pos - window_size
win_max     <- target_pos + window_size

win_df <- pbs_df[pbs_df$POS >= win_min & pbs_df$POS <= win_max, ]
cat("\nSNPs in ±10 kb window:", nrow(win_df), "\n")

# 99th-percentile threshold line for reference
p99 <- quantile(pbs_df$PBS_EAS, 0.99)

# Candidate SNP row for annotation
cand_df <- win_df[win_df$POS == target_pos, ]

p <- ggplot(win_df, aes(x = POS, y = PBS_EAS)) +

  # Shaded band around candidate SNP (±500 bp)
  annotate("rect",
           xmin = target_pos - 500, xmax = target_pos + 500,
           ymin = -Inf, ymax = Inf,
           fill = "gold", alpha = 0.30) +

  # 99th-percentile reference line
  geom_hline(yintercept = p99,
             linetype = "dotted", colour = "grey50", linewidth = 0.7) +
  annotate("text",
           x = win_min + (win_max - win_min) * 0.02,
           y = p99 * 1.05,
           label = "99th pctile",
           colour = "grey50", size = 3, hjust = 0) +

  # Vertical line at candidate position
  geom_vline(xintercept = target_pos,
             linetype = "dashed", colour = "grey30", linewidth = 0.6) +

  # All SNPs in window
  geom_point(colour = "#00A087", size = 2.2, alpha = 0.75) +
  geom_line(colour  = "#00A087", alpha = 0.35, linewidth = 0.5) +

  # Candidate SNP highlighted
  geom_point(data   = cand_df,
             aes(x = POS, y = PBS_EAS),
             colour = "red", size = 5, shape = 18) +

  # Label for candidate SNP
  geom_label(data  = cand_df,
             aes(x = POS, y = PBS_EAS, label = "rs3827760"),
             nudge_y   = max(win_df$PBS_EAS, na.rm = TRUE) * 0.07,
             size       = 3.2,
             colour     = "red",
             fill       = "white",
             fontface   = "bold") +

  scale_x_continuous(labels = scales::comma) +
  labs(
    title    = expression(paste("Population Branch Statistic (", PBS[EAS], ") around rs3827760")),
    subtitle = paste0("chr2:109,513,601  |  ±", window_size / 1000,
                      " kb window  |  candidate SNP in red"),
    x        = "Chromosomal position (bp)",
    y        = expression(PBS[EAS])
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(colour = "grey40", size = 10),
    panel.grid.minor = element_blank()
  )

print(p)
```

------------------------------------------------------------------------

## 3. Extended Haplotype Homozygosity (EHH)

### Extended Haplotype Homozygosity (EHH) and Haplotype Sweeps

Different approaches can detect genomic signatures of selection at different timescales. More recent selection signals can be detected from haplotype-based tests. Positive selection causes a rapid rise in the frequency of the selected allele, such that recombination does not have enough time to break down the haplotype on which the mutation arose. This creates a signature of **Extended Haplotype Homozygosity (EHH)** extending over a long physical distance.

### Questions for Students

1.  **How is the haplotype profile of genetic variants under recent positive selection?**
    -   *Answer*:
2.  **What is the profile of ancestral and derived haplotypes of the rs3827760 SNP in AFR and EAS?**
    -   *Answer*:
3.  **The iHS score observed for the SNP rs3827760 falls within which distribution quantiles of iHS values for the studied chromosome? Can it be considered an outlier? How can we make this analysis more robust?**
    -   *Answer*:
4.  **What information does the XP-EHH analysis add about natural selection in the candidate SNP?**
    -   *Answer*:

------------------------------------------------------------------------

### R Code Exercise: EHH & Furcation Trees

Write the R code necessary to: 1. Convert the VCF databases to `rehh` format using `data2haplohh()`. 2. Estimate the EHH decay for rs3827760 in both populations. 3. Plot the EHH decay and furcation trees for both AFR and EAS.

**Write your R code here:**

``` r
# 
```

------------------------------------------------------------------------

### R CodeExercise: iHS & XP-EHH (Window-based)

Write the R code necessary to: 1. Perform a genome-wide scan of homozygosity using `scan_hh()` for AFR and EAS. 2. Calculate iHS scores for both populations using `ihh2ihs()`. 3. Check the iHS score at rs3827760 and generate a single-site iHS plot in EAS. 4. Create a function to estimate the average absolute iHS in sliding windows (50 SNPs/40 step) and plot the results. 5. Estimate cross-population XP-EHH between EAS and AFR using `ies2xpehh()`, calculate window-based averages, and plot them.

**Write your R code here:**

``` r
# 
```

------------------------------------------------------------------------

## 4. Native American Selection Analysis

### Background

Hlusko et al. (2018), using morphological data, found a strong selection signal in the *EDAR* gene in Native Americans. Using the additional database from the 1000 Genomes Project (Peruvian samples with over 95% Native American Ancestry, represented as **NAM**), we evaluate genomic signatures of selection at the functional variant rs3827760.

### Questions for Students

1.  **Is the functional allele in East Asian at high frequency in other human populations (e.g. Native Americans)?**
    -   *Answer*:
2.  **Can we identify signatures of natural selection on EDAR in Native Americans using PBS?**
    -   *Answer*:
3.  **Is selection targeting the same functional variant?**
    -   *Answer*:
4.  **What is your conclusion based on the data generated?**
    -   *Answer*:

------------------------------------------------------------------------

### R Code Exercise: PBS in Native Americans (NAM)

Write the R code necessary to: 1. Read the pairwise $F_{ST}$ files involving Native Americans (`input/Part_1_HumanDiversity/Chr2_NAM_EAS.weir.fst` and `input/Part_1_HumanDiversity/Chr2_NAM_EUR.weir.fst`) and Europeans-East Asians (`input/Part_1_HumanDiversity/Chr2_EUR_EAS.weir.fst`). 2. Filter duplicates, exclude NA values, and align positions. 3. Convert negative $F_{ST}$ values to zero. 4. Estimate $PBS_{NAM}$ using NAM, EAS, and EUR. 5. Check PBS value at rs3827760, check quantiles, and plot the PBS scan.

**Write your R code here:**

``` r
# 
library(ggplot2)
library(dplyr)

# ── 1. Read input files ──────────────────────────────────────
infile_1 <- read.table(file = "Chr2_NAM_EAS.weir.fst", header = TRUE)
infile_2 <- read.table(file = "Chr2_NAM_EUR.weir.fst", header = TRUE)
infile_3 <- read.table(file = "Chr2_EUR_EAS.weir.fst", header = TRUE)

cat("Rows read — AFR_EAS:", nrow(infile_1),
    "| AFR_EUR:", nrow(infile_2),
    "| EAS_EUR:", nrow(infile_3), "\n")

# ── 2. Filter duplicate positions & remove NAs ───────────────
clean_fst <- function(df, label) {
  df <- df[!duplicated(df$POS), ]          # drop duplicate SNP positions
  df <- df[!is.na(df$WEIR_AND_COCKERHAM_FST), ]  # exclude NA FST values
  cat("After filtering [", label, "] :", nrow(df), "SNPs\n")
  df
}

fst1 <- clean_fst(infile_1, "AFR_EAS")
fst2 <- clean_fst(infile_2, "AFR_EUR")
fst3 <- clean_fst(infile_3, "EAS_EUR")

# ── 3. Align datasets by overlapping positions ───────────────
shared_pos <- Reduce(intersect, list(fst1$POS, fst2$POS, fst3$POS))
cat("Overlapping positions across all three comparisons:", length(shared_pos), "\n")

fst1 <- fst1[fst1$POS %in% shared_pos, ]
fst2 <- fst2[fst2$POS %in% shared_pos, ]
fst3 <- fst3[fst3$POS %in% shared_pos, ]

# Ensure consistent ordering by position
fst1 <- fst1[order(fst1$POS), ]
fst2 <- fst2[order(fst2$POS), ]
fst3 <- fst3[order(fst3$POS), ]

# ── 4. Set negative FST values to zero ──────────────────────
fst1$WEIR_AND_COCKERHAM_FST <- pmax(fst1$WEIR_AND_COCKERHAM_FST, 0)
fst2$WEIR_AND_COCKERHAM_FST <- pmax(fst2$WEIR_AND_COCKERHAM_FST, 0)
fst3$WEIR_AND_COCKERHAM_FST <- pmax(fst3$WEIR_AND_COCKERHAM_FST, 0)
```

---
# Part 2: Genomic Selection Sweep Scan in Canines
The dataset is sourced from the **Dog10K** consortium ([Download Link](https://dog10k.kiz.ac.cn/Home/Download)). The original genomic dataset is a high-coverage phased BCF file containing 1,929 individuals and over 29 million SNPs:
- Original file: `AutoAndXPAR.Dog10K.phased.bcf`
- Metadata table: `dog10K-alignment-sample-table.2022-02-23-v7.txt`
For this analysis, we will use a subset of **130 individuals** representing body size extremes:
| :--- | :--- | :---: |
editor_options: 
  markdown: 
    wrap: sentence
---

## 2. Preprocessing and Filtering

We filter the massive BCF file to include only our 130 samples and chromosome 15, while removing low-frequency SNPs (MAF \< 0.05) that are not informative for breed/size differentiation. The coordinates correspond to the CanFam3 reference genome, where the ***IGF1*** gene is located approximately in the **41.2 Mb to 44.5 Mb** region.

### Bash Code Exercise 1: Extraction & Format Conversion

``` bash
# 1. Extract chr15 and filter for samples and MAF >= 0.05
bcftools view \
  -S input/Part_2_CanidDiversity/subset_dogs.txt \
  -r chr15 \
  -q 0.05:minor \
  -O b \
  -o input/Part_2_CanidDiversity/subset_chr15.bcf \
  input/Part_2_CanidDiversity/AutoAndXPAR.Dog10K.phased.bcf

# 2. Index the subset BCF
bcftools index input/Part_2_CanidDiversity/subset_chr15.bcf

# 3. Convert BCF to PLINK binary format (.bed/.bim/.fam)
plink1.9 \
  --bcf input/Part_2_CanidDiversity/subset_chr15.bcf \
  --dog \
  --keep-allele-order \
  --make-bed \
  --out output/subset_chr15
```

> **Premise**: The filtered BCF file contains **177,953 SNPs** on chromosome 15 across the 130 samples. The coordinates correspond to the CanFam3 reference genome, where the ***IGF1*** gene is located approximately in the **41.2 Mb to 44.5 Mb** region.

------------------------------------------------------------------------

## 3. Population Structure Analysis (PCA)

Before checking for selection outliers, we must examine the genetic structure of our subset. We will run a PCA using **PLINK 1.9** and visualize it in R.

### R Code Exercise 2: PCA Visualization in R

Write the R code necessary to: 1. Load the eigenvectors and eigenvalues generated by PLINK. 2. Merge them with the sample metadata (`input/Part_2_CanidDiversity/sample_info.txt`). 3. Calculate the percentage of variance explained by PC1 and PC2. 4. Generate a scatter plot of PC1 vs. PC2, coloring the points by breed and shaping them by size group (small vs. large).

**Write your R code here:**

``` r
# 
```

### Questions for Students

1.  **What pattern do you observe along the first principal component (PC1)?**
    -   *Answer*:
2.  **Why does PC1 capture body size differences in this particular dataset?**
    -   *Answer*:

------------------------------------------------------------------------

## 4. Genomic Outlier Detection using PCAdapt

**PCAdapt** is a method designed to find SNPs that are exceptionally related to population structure (PCs) rather than neutral drift.

### The Role of Linkage Disequilibrium (LD) and Clumping

> [!IMPORTANT] **Key Concept**: PCA is highly sensitive to Linkage Disequilibrium (LD). If a region contains many highly correlated markers (due to a selective sweep or low recombination), that single region will dominate the principal components, biasing the PCA and masking other genomic signals (such as the *IGF1* gene sweep).
>
> To resolve this, we must enable **LD Clumping** in PCAdapt. Thinning out redundant SNPs in strong LD allows the global genomic structure to be correctly computed and helps locate narrow selection sweeps.

### R Code Exercise 3: PCAdapt with LD Clumping

Write the R script necessary to: 1. Load the genotypes into PCAdapt. 2. Execute `pcadapt()` using $K = 2$ components and enable **LD clumping** (for example, with a window size of 500 SNPs and an $r^2$ threshold of 0.1). 3. Merge the resulting p-values with the physical genomic positions from the `.bim` file. 4. Generate a Manhattan Plot of the results, plotting physical position (in Mb) on the X-axis and $-\log_{10}(\text{p-value})$ on the Y-axis.

**Write your R code here:**

``` r
# 
```

### Questions for Students

1.  **If you ran PCAdapt without enabling LD clumping, a single massive peak at \~61 Mb would dominate the plot, hiding other regions. What is the effect of LD clumping on outlier detection and why is it necessary here?**
    -   *Answer*:
2.  **Did you detect the outlier peak at the *IGF1* locus? What is the approximate coordinate of the peak and what is its biological significance?**
    -   *Answer*:

------------------------------------------------------------------------

## 5. Cross-Population Selection Scan using XP-nSL

To confirm that the outlier peak on chromosome 15 is indeed driven by a selective sweep, we will perform a haplotype-based selection scan. Specifically, we will run **XP-nSL (Cross-Population Number of Segregating Sites by Length)** to compare the haplotype homozygosity decay between small dogs and large dogs.

### Key Concepts: nSL and XP-nSL

-   **nSL**: A within-population selection scan metric similar to iHS. However, instead of measuring haplotype decay in terms of genetic distance (which requires a genetic map), nSL measures distance by counting the number of segregating sites (segregating site count by length). This makes it highly robust to recombination rate variation and suitable for genomes without well-defined genetic maps.
-   **XP-nSL**: A cross-population statistic that compares nSL profiles between a target population and a reference population. A high positive score indicates a selective sweep specific to the target population (longer haplotypes around the derived allele).
-   **Phased Mode**: Since our input Dog10K BCF file is already phased (containing haplotype data formatted as `0|0`, `1|0`, etc.), we will perform a phased XP-nSL scan. This utilizes the precise haplotype sequences, which provides a significantly stronger selection signal compared to unphased analyses.

### Outgroup Allele Polarization

Haplotype selection scans require knowing which allele is **ancestral** (original) and which is **derived** (new mutant). `selscan` expects a VCF file where `0` is the ancestral allele and `1` is the derived allele. To polarize our dataset, we use the **gray wolves** in the Dog10K metadata as an outgroup: - Gray wolves are the evolutionary ancestor of domestic dogs. - For each SNP, the most common (major) allele in the gray wolf population is designated as the ancestral allele. - If the ALT allele in the original VCF is the major allele in wolves, we must physically swap the REF/ALT alleles and swap genotypes (`0` becomes `1`, and `1` becomes `0`) for all individuals.

------------------------------------------------------------------------

### Bash Code Exercise 4: Extracting and Polarizing Alleles

``` bash
# 1. Extract wolf samples and combine with dog samples
bcftools query -l input/Part_2_CanidDiversity/AutoAndXPAR.Dog10K.phased.bcf | grep -E '^CLUP' > input/Part_2_CanidDiversity/wolves.txt
cat input/Part_2_CanidDiversity/subset_dogs.txt input/Part_2_CanidDiversity/wolves.txt > input/Part_2_CanidDiversity/dogs_and_wolves.txt

# 2. Extract polymorphic sites for both dogs and wolves
bcftools query -f '%CHROM\t%POS\n' input/Part_2_CanidDiversity/subset_chr15.bcf > input/Part_2_CanidDiversity/subset_chr15_positions.txt
bcftools view \
  -S input/Part_2_CanidDiversity/dogs_and_wolves.txt \
  -T input/Part_2_CanidDiversity/subset_chr15_positions.txt \
  -O z \
  -o input/Part_2_CanidDiversity/subset_chr15_with_wolves.vcf.gz \
  input/Part_2_CanidDiversity/AutoAndXPAR.Dog10K.phased.bcf

# 3. Polarize alleles using Python (major allele in wolves = 0)
python3 scripts/polarize_by_wolves.py

# 4. Re-compress to block gzip format and index polarized VCF
bcftools view input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz -O z -o output/subset_chr15_polarized_bgzf.vcf.gz
mv output/subset_chr15_polarized_bgzf.vcf.gz input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz
bcftools index input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz

# 5. Extract polarized VCFs for small and large dogs separately
bcftools view -S input/Part_2_CanidDiversity/small_dogs.txt -O z -o input/Part_2_CanidDiversity/small_dogs_polarized.vcf.gz input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz
bcftools index input/Part_2_CanidDiversity/small_dogs_polarized.vcf.gz

bcftools view -S input/Part_2_CanidDiversity/large_dogs.txt -O z -o input/Part_2_CanidDiversity/large_dogs_polarized.vcf.gz input/Part_2_CanidDiversity/subset_chr15_polarized.vcf.gz
bcftools index input/Part_2_CanidDiversity/large_dogs_polarized.vcf.gz
```

------------------------------------------------------------------------

### Selection Scan Execution and Normalization

We will execute the selection scan using the compiled `selscan` binary, comparing small dogs (target) against large dogs (reference) in phased mode, and normalize the raw scores.

### Bash Code Exercise 5: Running XP-nSL and Normalizing

Write the bash commands necessary to: 1. Execute the `selscan` program in phased XP-nSL mode using the polarized target (small dogs) and reference (large dogs) VCFs. Use 4 threads and output results to `output/xpnsl_phased`. 2. Normalize the raw XP-nSL scores using the `selscan norm` utility. 3. Calculate window-based statistics (fraction of outliers) in 100 Kb non-overlapping windows.

**Write your Bash code here:**

``` bash
# 
```

------------------------------------------------------------------------

### R Code Exercise 6: XP-nSL Haplotype Manhattan Plots

Write the R code necessary to: 1. Load the normalized SNP-level XP-nSL scores from `input/Part_2_CanidDiversity/xpnsl_phased.xpnsl.out.norm`. 2. Load the 100 Kb window-level data from `input/Part_2_CanidDiversity/xpnsl_phased.xpnsl.out.norm.100kb.windows`. 3. Highlight the *IGF1* region (between 41 Mb and 45.5 Mb). 4. Generate two Manhattan plots: one for the raw SNP-level scores (calculating $-\log_{10}(p\text{-value})$ from the normalized Z-scores), and another for the window-based fraction of extreme positive SNPs (`frac_top`).

**Write your R code here:**

``` r
# 
```

### Questions for Students

1.  **Why do we need to polarize alleles using an outgroup like the gray wolf? What does the `0` vs `1` coding represent in `selscan`?**
    -   *Answer*:
2.  **Why does the raw SNP-level XP-nSL scan look like a noisy cloud at individual sites? What is the effect of calculating window-based scores (e.g. 100 Kb)?**
    -   *Answer*:

------------------------------------------------------------------------

## 6. Alternative Haplotype Selection Scan: Rsb using `rehh`

As a complementary approach to XP-nSL, we will run **Rsb**, another widely used cross-population EHH-based statistic.

### Key Concepts & Comparison

-   **Rsb**: Compares the integrated Extended Haplotype Homozygosity (iHH) between two populations. It is calculated as $\ln(iES_{pop1} / iES_{pop2})$, where $iES$ is the integrated EHH over physical distance (bp). A high positive score indicates selection in population 1 (small dogs), while a negative score indicates selection in population 2 (large dogs).
-   **Difference from XP-nSL**:
    -   **XP-nSL** integrates the nSL metric over the number of segregating sites (SNP count). This makes it highly robust to local recombination rate variation.
    -   **Rsb** integrates EHH over physical distance (bp). In species with very strong selective sweeps and long-range Linkage Disequilibrium (like domestic dogs), Rsb can produce exceptionally high, clear peaks at sweep loci like *IGF1*.

### Phase Preservation during Outgroup Polarization

A common question in haplotype selection scans is: **Does swapping the REF and ALT alleles (and flipping 0 to 1 and 1 to 0) to align with the outgroup corrupt or destroy the phasing information?** The answer is **No**. In a phased VCF, genotypes are represented as `0|1` or `1|0` to denote which allele lies on which homologous chromosome (haplotype). Swapping REF and ALT alleles and swapping `0` and `1` (converting `0|1` to `1|0` and vice-versa) is a mathematically symmetric operation. It maintains the exact same physical haplotype alignment across all sites on each chromosome, merely updating the label of which allele is ancestral and which is derived. Thus, outgroup polarization is fully compatible with phased haplotype scans and no data is lost.

------------------------------------------------------------------------

### R Code Exercise 7: Running Rsb in R using `rehh`

Write the R code necessary to: 1. Load the polarized Small Dogs and Large Dogs VCFs into `rehh` using `data2haplohh()`. 2. Compute the haplotype homozygosity scan for both populations using `scan_hh()`. 3. Compute the cross-population Rsb statistic using `ines2rsb()`. 4. Plot the Rsb Manhattan plot using `ggplot2`, highlighting the *IGF1* region between 41 Mb and 45.5 Mb in red.

**Write your R code here:**

``` r
# 
```

------------------------------------------------------------------------

### Questions for Students

1.  **Explain the physical and mathematical difference between Rsb and XP-nSL. Why does Rsb show a much higher, less noisy peak at the *IGF1* locus in dogs compared to XP-nSL?**
    -   *Answer*:
2.  **Does outgroup polarization of a phased VCF file corrupt or destroy the phasing information? Why or why not?**
    -   *Answer*:
