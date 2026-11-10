

###########################################################################
################################ RUNNING ##################################
###########################################################################


# Reading the functions in the Functions folder ---------------------------

functions_files <- list.files("Functions", pattern = ".rds")

functions_names <-
  gsub(pattern = ".rds", replacement = "", functions_files)
  
names(functions_names) <- functions_names

functions_path <- paste0("Functions/", functions_files)

for (i in 1:length(functions_path)) {
  
  assign(functions_names[i], readRDS(functions_path[i]))
  
}


# Installing (if necessary) and requiring the packages --------------------

packages = c(
  "RSNNS",
  "tidyverse",
  "gt",
  "gtExtras",
  "NeuralNetTools",
  "caret",
  "ggformula",
  "flextable"
)

package.check(packages)


# Reading the real data ---------------------------------------------------

data <- readxl::read_xlsx("Data/new_data.xlsx")

# Replace NAs with the average
data |>
  group_by(LINHAGEM, LOCAL) |>
  mutate_at(vars(REND), ~
              case_match(.x,
                         NA ~ mean(.x, na.rm = T),
                         .default = .x)) |>
  ungroup() -> data

data


# Running the analysis of Cruz, Torres & Vencovsky ------------------------

data_ctv <-
  ctv(
    data,
    genotype = "LINHAGEM",
    env = "LOCAL",
    block = "REP",
    y = "REND"
  )

data_ctv$curve_parameters |>
  view()


# Graphical visualization of the real data --------------------------------

cr <- length(data_ctv$data_for_classification)
real_data <- data_ctv$data_for_classification[,-cr]

show_graphic(
  real_data,
  jump = 1,
  y_lab = "Production (Kg/ha)",
  env_index = c(),
  x_lab = "Environmental index",
  obs_name = "Lineage",
  group_name = "Lineage:",
  obs_column = 1,
  y_lim = c(2000, 7000),
  y_break = 1000,
  sec = 3
)


# Collecting the parameters for simulation --------------------------------

parameters <- data_ctv$parameters_for_simulation


# Running the data simulation ---------------------------------------------

set.seed(4321)

simulated_data <-
  train.set_ctv(parameters, n_training = 500, alpha = 0.01)

 #simulated_data <- readRDS("Simulated_Data/simulated_data.rds")

sim_data <- simulated_data$data3

colMeans(sim_data[, -1]) - mean(as.matrix(sim_data[, -1]))


# Graphical visualization of the simulated data ---------------------------

show_graphic(
  sim_data,
  jump = 500,
  y_lab = "Production (Kg/ha)",
  env_index = c(),
  x_lab = "Environmental index",
  obs_name = "Genotype",
  group_name = "Group",
  obs_column = 1,
  y_lim = c(2000, 7000),
  y_break = 1000,
  sec = 3
)


# Training the neural network ---------------------------------------------

trained.mlp_ctv <- train.mlp.net(
  data = simulated_data$data2,
  decay = 0.0005,
  maxit = 8000,
  error.max = 0.01,
  neurons = seq(5, 20),
  layers = 2, 
  column_group = 1,
  p_train = 0.8,
  norm_type = "norm"
)



# Saving the trained neural network ---------------------------------------

trained_mlp_ctv <- trained.mlp_ctv$mlp_with_smallest_error

# save if necessary read it later (Do not forget to change date before saving!)

saveRDS(trained_mlp_ctv, file = "Trained_Nets/2025_11_11_trained_mlp_ctv.rds") 

list.files("Trained_Nets")
#trained_mlp_ctv <- readRDS(file = "Trained_Nets/2023_17_10_trained_mlp_ctv.rds")


# Classification of real genotypes ----------------------------------------

data_for_classification <- data_ctv$data_for_classification

ranking <- net.adaptability.analysis(data_for_classification,
                                     trained_mlp_ctv,
                                     norm_type = "norm")

ranking

ranking |>
  writexl::write_xlsx("Results/ranking_2025_11_11.xlsx")

# Visualization of the classification of the real genotypes ---------------

gt_table(ranking)
ft_table(ranking) # Maybe it works! Just another representation of ranking table.
