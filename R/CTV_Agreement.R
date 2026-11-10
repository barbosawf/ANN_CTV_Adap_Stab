###########################################################################
# Agreement between the ANN and the classical CTV method (real data)
#
# 1) Runs the classical CTV analysis (Cruz, Torres and Vencovsky 1989) on
#    the real data and classifies each genotype with t tests (beta1 = 1,
#    beta1 + beta2 = 1) and the F test (sigma2_d = 0) at alpha = 0.05.
# 2) Classifies the same genotypes with the trained network used in the
#    article (Trained_Nets/2025_11_11_trained_mlp_ctv.rds).
# 3) Reports the agreement between both methods for beta1 (unfavorable
#    environments), beta1 + beta2 (favorable environments) and stability.
#
# The function ctv() returns only p-values (in %); the classification at
# alpha is done here. alpha = 0.05 refers to the classical CTV tests used in
# the comparison; the simulation used alpha = 0.01 (see R/Adap_Stab_CTV.R).
#
# Requires Data/new_data.xlsx (not distributed). Run from the project root:
#   Rscript R/CTV_Agreement.R
###########################################################################

alpha <- 0.05
net_file <- "Trained_Nets/2025_11_11_trained_mlp_ctv.rds"

suppressMessages({library(RSNNS); library(tidyverse)})

# Functions
for (f in list.files("Functions", pattern = ".rds")) {
  assign(gsub(".rds", "", f), readRDS(file.path("Functions", f)))
}

# Real data (missing plots replaced by the genotype x environment mean,
# as in R/Adap_Stab_CTV.R)
data <- readxl::read_xlsx("Data/new_data.xlsx")
data <- data |>
  group_by(LINHAGEM, LOCAL) |>
  mutate(REND = ifelse(is.na(REND), mean(REND, na.rm = TRUE), REND)) |>
  ungroup()

data_ctv <- ctv(data, genotype = "LINHAGEM", env = "LOCAL",
                block = "REP", y = "REND")
cp <- data_ctv$curve_parameters

# Classical CTV classification at alpha (p-values from ctv() are in %)
classify_beta <- function(p_percent, estimate) {
  ifelse(p_percent >= 100 * alpha, "=1", ifelse(estimate < 1, "<1", ">1"))
}
ctv_class <- tibble(
  Genotype  = as.character(cp[[1]]),
  CTV_b1    = classify_beta(cp$prob_b1, cp$B1),
  CTV_b1_b2 = classify_beta(cp$prob_b2_b2, cp$B1...B2),
  CTV_stab  = ifelse(cp$prob_dev >= 100 * alpha, "High", "Low")
)

# ANN classification (net.adaptability.analysis reads data_for_classification
# from the global environment)
data_for_classification <- data_ctv$data_for_classification
net <- readRDS(net_file)
ann <- net.adaptability.analysis(data_for_classification, net, norm_type = "norm")
ann_class <- tibble(
  Genotype  = as.character(ann[[1]]),
  ANN_b1    = substr(ann$Adaptability, 2, 3),
  ANN_b1_b2 = substr(ann$Adaptability, 7, 8),
  ANN_stab  = ann$Stability
)

both <- inner_join(ctv_class, ann_class, by = "Genotype")
stopifnot(nrow(both) == nrow(ctv_class))

# Cohen's kappa: agreement corrected for the agreement expected by chance
cohen_kappa <- function(a, b) {
  lev <- union(a, b)
  tab <- table(factor(a, levels = lev), factor(b, levels = lev))
  n <- sum(tab)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / n^2
  (po - pe) / (1 - pe)
}

agreement <- tibble(
  Parameter = c("beta1 (unfavorable environments)",
                "beta1 + beta2 (favorable environments)",
                "sigma2_d (stability)"),
  Agreements = c(sum(both$CTV_b1 == both$ANN_b1),
                 sum(both$CTV_b1_b2 == both$ANN_b1_b2),
                 sum(both$CTV_stab == both$ANN_stab)),
  Genotypes = nrow(both),
  Kappa = c(cohen_kappa(both$CTV_b1, both$ANN_b1),
            cohen_kappa(both$CTV_b1_b2, both$ANN_b1_b2),
            cohen_kappa(both$CTV_stab, both$ANN_stab))
) |>
  mutate(Percent = round(100 * Agreements / Genotypes, 1),
         Kappa = round(Kappa, 2)) |>
  relocate(Percent, .before = Kappa)

cat(sprintf("\nAgreement between the ANN (%s) and the classical CTV method (alpha = %.2f)\n\n",
            basename(net_file), alpha))
print(as.data.frame(agreement), row.names = FALSE)

cat("\nCross-tabulations (rows: classical CTV; columns: ANN):\n")
for (p in list(c("CTV_b1", "ANN_b1", "beta1"), c("CTV_b1_b2", "ANN_b1_b2", "beta1 + beta2"),
               c("CTV_stab", "ANN_stab", "stability"))) {
  cat("\n", p[3], "\n"); print(table(CTV = both[[p[1]]], ANN = both[[p[2]]]))
}

cat("\nGenotypes classified differently by the two methods:\n")
both |>
  filter(CTV_b1 != ANN_b1 | CTV_b1_b2 != ANN_b1_b2 | CTV_stab != ANN_stab) |>
  as.data.frame() |>
  print(row.names = FALSE)

# Uncomment to save the genotype-level comparison
# writexl::write_xlsx(both, "Results/ctv_vs_ann_agreement.xlsx")
