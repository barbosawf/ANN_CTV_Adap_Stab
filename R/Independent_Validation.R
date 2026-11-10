###########################################################################
# Independent validation of the selected network (reviewer comment 1)
#
# Evaluates the trained network (Trained_Nets/2025_11_11_trained_mlp_ctv.rds,
# architecture 14-9-17-18) WITHOUT retraining on:
#   1) the regenerated original simulated set (seed 4321, train + test);
#   2) three new independent simulated sets (seeds 2026, 7, 99; 9,000
#      genotypes each), never used for training or model selection;
#   3) the real data, comparing the paper's normalization (own z-score)
#      with the training-set normalization parameters.
# Normalization parameters are taken from the regenerated original set.
# Run from the project root: Rscript R/Independent_Validation.R
###########################################################################

suppressMessages({library(RSNNS); library(tidyverse)})
ff <-list.files("Functions", pattern = ".rds")
for (f in ff) assign(gsub(".rds", "", f), readRDS(file.path("Functions", f)))

data <- readxl::read_xlsx("Data/new_data.xlsx")
data |> group_by(LINHAGEM, LOCAL) |>
  mutate_at(vars(REND), ~ case_match(.x, NA ~ mean(.x, na.rm = T), .default = .x)) |>
  ungroup() -> data
data_ctv <- ctv(data, genotype = "LINHAGEM", env = "LOCAL", block = "REP", y = "REND")
parameters <- data_ctv$parameters_for_simulation

net <- readRDS("Trained_Nets/2025_11_11_trained_mlp_ctv.rds")
cat("Net size:", net$archParams$size, "\n")

acc_fun <- function(net, X, y, normParams) {
  Xn <- normalizeData(X, normParams)
  p <- apply(predict(net, Xn), 1, which.max)
  cm <- table(factor(y, 1:18), factor(p, 1:18))
  k <- caret::confusionMatrix(cm)
  adap <- mean(((y - 1) %% 9) == ((p - 1) %% 9))
  stab <- mean((y > 9) == (p > 9))
  c(acc = unname(k$overall["Accuracy"]), kappa = unname(k$overall["Kappa"]),
    adap_acc = adap, stab_acc = stab)
}

# 1) Regenerate the original simulated set (same seed as the paper's script)
set.seed(4321)
t0 <- Sys.time()
orig <- train.set_ctv(parameters, n_training = 500, alpha = 0.01)$data2
cat("Sim time (s):", as.numeric(Sys.time() - t0, units = "secs"), "\n")
Xo <- as.matrix(orig[, -1]); yo <- as.numeric(as.character(orig[[1]]))
Xo_n <- normalizeData(Xo, type = "norm")
np <- getNormParameters(Xo_n)
cat("\nOriginal simulated set (train+test, 9000):\n"); print(round(acc_fun(net, Xo, yo, np), 4))

# 2) Independent new simulated sets, never used for training/selection
res <- NULL
for (s in c(2026, 7, 99)) {
  set.seed(s)
  ind <- train.set_ctv(parameters, n_training = 500, alpha = 0.01)$data2
  Xi <- as.matrix(ind[, -1]); yi <- as.numeric(as.character(ind[[1]]))
  res <- rbind(res, c(seed = s, acc_fun(net, Xi, yi, np)))
}
cat("\nIndependent simulated sets (9000 each):\n"); print(round(res, 4))

# 3) Real data: as in the paper (own z-score) vs training normalization params
dc <- data_ctv$data_for_classification
Xr <- as.matrix(dc[, -1])
p_own <- apply(predict(net, normalizeData(Xr, type = "norm")), 1, which.max)
p_trn <- apply(predict(net, normalizeData(Xr, np)), 1, which.max)
cat("\nReal data: agreement own-norm vs train-norm class:", mean(p_own == p_trn), "\n")
rk <- readxl::read_xlsx("Results/ranking_2025_11_11.xlsx")
print(head(rk, 5))
cat("Input column names:", colnames(dc)[-1], "\n")
