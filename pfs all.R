################################################################################
# IPD META-ANALYSIS: survival workflow
#
# - Checks proportional hazards (global + Schoenfeld plots + per-study checks)
# - If PH plausible -> stratified Cox (recommended) for pooled HR
# - If PH violated -> alternatives: time-varying Cox, piecewise (landmark),
#   flexible parametric (Royston-Parmar), and RMST
# - Plots: KM with HR annotation, RMST curves, RMST difference shaded plot
# - Clear didactic comments explain "why", "when", and "how to interpret"
#
# NOTE: Change the import path to your IPD file if different.
################################################################################

# ---- Packages ---------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  rio, survival, survminer, tidyverse, survRM2, broom,
  patchwork, glue
)

#import custom functions
source("scripts/functions/custom_km_plot.R")
source("scripts/functions/tvhr_plot.R")
source("scripts/functions/compute_followup_stats.R")
source("scripts/functions/rmst_plot.R")

#make output dir
dir.create("output/km plots", recursive = T)
dir.create("output/rmst plots", recursive = T)
dir.create("output/tvhr plots", recursive = T)
dir.create("output/schoenfeld residuals plots", recursive = T)


# ---- Import data ------------------------------------------------------------
# Replace path with your file. The dataset must contain:
#   - time: follow-up time (months)
#   - status: event indicator (1=event, 0=censored)
#   - treat: treatment indicator (0 = control, 1 = treatment) or factor with those labels
#   - study: study identifier (for stratification/cluster)
# Import digitized KM data (pfs_all endpoint)

#dual ici + ct
poseidon_pfs_all <- import("ipd/dual ici plus chemo/pfs/poseidon pfs all.csv") |>
  mutate(study = "poseidon")

cm9la_pfs_all <- import("ipd/dual ici plus chemo/pfs/cm 9LA 6y pfs all.csv") |> 
  mutate(study = "cm9la")

#dual ici alone
mystic_pfs_all <- import("ipd/dual ici/pfs/mystic pfs all.csv") |> 
  mutate(study = "mystic")

cm227_pfs_pdl1pos <- import("ipd/dual ici/pfs/cm227 pfs pdl1pos.csv") |> 
  mutate(study = "cm227")

cm227_pfs_pdl1neg <- import("ipd/dual ici/pfs/cm227 pfs pdl1neg.csv") |> 
  mutate(study = "cm227")

# Combine individual patient data (IPD)
combined_pfs_all <- rbind(poseidon_pfs_all,
                         cm9la_pfs_all,
                         mystic_pfs_all,
                         cm227_pfs_pdl1pos,
                         cm227_pfs_pdl1neg) |> 
  mutate(treat_num = case_when(treat == "Dual ICI" ~ as.numeric(0),
                               treat == "Dual ICI + CT" ~ as.numeric(1)))

# ---- Basic sanity checks ----------------------------------------------------
glimpse(combined_pfs_all)
stopifnot(all(c("time", "status", "treat", "study") %in% names(combined_pfs_all)))

# Make sure treat and study are factors and create numeric treat for some models
combined_pfs_all <- combined_pfs_all |> 
  mutate(
    treat = factor(treat, levels = c("Dual ICI", "Dual ICI + CT")),
    study = factor(study)
  )

# Quick per-study descriptives (why: to check sample sizes and event balance)
combined_pfs_all |> 
  group_by(study, treat) |> 
  summarise(
    n = n(),
    events = sum(status),
    median_follow_up = median(time[status == 0], na.rm = TRUE),
    .groups = "drop"
  ) |>  print(n = Inf)

# ---- 1) Cox models for diagnostics ------------------------------------------
# 1.a Pooled Cox
coxph_pooled <- coxph(Surv(time, status) ~ treat, data = combined_pfs_all)
summary(coxph_pooled)

# ---- 2) PH diagnostics (Grambsch-Therneau + plots) ---------------------------
# Why: Cox model requires proportional hazards (PH) for the HR to be constant.
# If PH is violated, pooled HR is not fully adequate — consider time-varying or RMST.

# 2.a cox.zph on pooled model (approximate)
zph_pooled <- cox.zph(coxph_pooled)
print(zph_pooled)

png("output/schoenfeld residuals plots/pfs_all.png", width = 6000, height = 4000, res = 300)  
ggcoxzph(
  zph_pooled
)
dev.off() 

# Short guidance:
# - If zph > 0.05 and plots show no trend -> PH plausible.
# - If p < 0.05 or clear trend -> PH likely violated; proceed to alternatives.

# ---- 3) If PH plausible -> report stratified Cox HR and KM -------------------
# Kaplan-Meier for visualization (KM does not account for study stratification)

tiff("output/km plots/pfs_all.tiff", width = 6000, height = 4000, res = 300)  
custom_km_plot(combined_pfs_all)
dev.off()  


# Median follow-up time
compute_followup_stats(combined_pfs_all)


# ---- 4) If PH violated -> alternatives (we run them as sensitivity analyses)-
# A) Use time-varying coefficient in Cox (treatment * g(time))
# B) Use RMST (robust, interpretable difference in mean survival up to tau)
################################################################################

# ----- A Time-varying Cox (tt) ---------------------------------------------
# Why: estimate HR(t) when effect changes over time. Use tt() with a function g(t).

png("output/tvhr plots/pfs_all.png", width = 6000, height = 4000, res = 300) 
tvhr_plot(combined_pfs_all)
dev.off()

# ---- C) RMST analysis (robust, interpretable) -------------------------------
#  USER PARAMETERS

arm0_label <- "Dual ICI"           # control
arm1_label <- "Dual ICI + CT"      # treated
tau <- 60

png("output/rmst plots/pfs_all.png", width = 6000, height = 4000, res = 300) 
rmst_plot(data = combined_pfs_all, arm0_label =arm0_label, arm1_label = arm1_label, tau = tau,
          ypos = 0.9, xposfactor = 0.5) |> 
  print()
dev.off()

# Reporting recommendations (what to include in paper / report) ----------
# - Report which PH tests were used: cox.zph results (global and per-covariate), and
#   show Schoenfeld plots. Mention caveats: clustered data -> interpret cautiously.
# - If PH holds: report Cox estimate (HR), 95% CI, p-value.
# - If PH violated: present one main alternative (justify choice) and show sensitivity analyses:
#     * time-varying HR curve
#     * RMST difference at pre-specified tau(s)
# - Always show Kaplan-Meier curves for visualization (with number at risk).
# - If using random-effects Cox (coxme), report estimated between-study variance (sigma) and interpret heterogeneity accordingly.
#

################################################################################
# END OF SCRIPT
#
# Notes:
# - Tweak tau and landmark_time according to clinical context.
# - When reporting choose 1 primary approach (e.g.,Cox if PH holds,
#   RMST if PH violated) and present others as sensitivity analyses.
################################################################################