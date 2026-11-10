# Adap_Stab_CTV — ANN based on bisegmented regression (CTV) for adaptability and stability

R code and the trained artificial neural network (ANN) used in:

> *Artificial neural network based on bisegmented regression for adaptability and stability analysis.*
> Crop Breeding and Applied Biotechnology (manuscript CBAB-2026-5660, under review).

The ANN classifies genotypes into 18 adaptability × stability classes defined by the
Cruz, Torres and Vencovsky (1989) bisegmented regression method (CTV). It is trained on
genotypes simulated from the CTV model and then applied to real multi-environment trial data.

## Repository contents

| Path | Description |
|---|---|
| `R/Functions.R` | Source of all functions (CTV analysis, data simulation, ANN training, classification, tables). Running it regenerates the files in `Functions/`. |
| `R/Adap_Stab_CTV.R` | Main pipeline: CTV analysis of the real data → simulation → ANN training → classification of real genotypes. |
| `R/Independent_Validation.R` | Evaluates the trained network, without retraining, on new independent simulated datasets. |
| `R/CTV_Agreement.R` | Classifies the real genotypes with the classical CTV method (t and F tests, α = 0.05) and computes the agreement with the ANN (percentage and Cohen's kappa for β1, β1 + β2 and stability). |
| `Functions/*.rds` | Serialized versions of the functions in `R/Functions.R` (read by the main script). |
| `Trained_Nets/2025_11_11_trained_mlp_ctv.rds` | **The trained network used in the article** (RSNNS `mlp` object). |
| `Data/` | Empty. The field data are not distributed (see below). |
| `Results/` | Output folder for the classification tables. |

## Network used in the article

| Item | Setting |
|---|---|
| Type | Multilayer perceptron (feedforward), RSNNS `mlp()` |
| Architecture | 14 inputs – 9 – 17 – 18 outputs |
| Inputs | 13 genotype mean yields (kg ha⁻¹), ordered by increasing environmental index, + estimated variance of the regression deviations (σ²δ) |
| Outputs | 18 classes (9 adaptability × 2 stability), 1-of-18 coding; class = highest output |
| Activation | Logistic (hidden and output layers) |
| Learning | Standard backpropagation, online, shuffled patterns; learning rate 0.0005; 8,000 epochs |
| Weight initialization | Uniform in [−0.3, 0.3] (RSNNS default) |
| Normalization | z-score (`type = "norm"`) |
| Training data | 9,000 simulated genotypes (500 per class; seed 4321); 80% training / 20% testing, stratified |

Load and use it:

```r
library(RSNNS)
net <- readRDS("Trained_Nets/2025_11_11_trained_mlp_ctv.rds")
# x: matrix with 14 columns (13 environment means ordered by environmental index + S2d)
pred <- predict(net, normalizeData(x, type = "norm"))
class <- apply(pred, 1, which.max)
```

## Data availability

The soybean yield data (33 genotypes, 13 environments, 2020/2021, Embrapa Soja) are not
publicly distributed. They can be requested from the corresponding author, subject to
Embrapa Soja's authorization. To run the full pipeline, place the data file at
`Data/new_data.xlsx` with columns `LINHAGEM` (genotype), `LOCAL` (environment), `REP` (block)
and `REND` (yield, kg ha⁻¹).

## How to run

Open `Adap_Stab_CTV.Rproj` (the scripts use paths relative to the project root), then:

```r
source("R/Adap_Stab_CTV.R")          # full pipeline (requires Data/new_data.xlsx)
source("R/Independent_Validation.R")  # independent validation of the trained network
source("R/CTV_Agreement.R")           # agreement between the ANN and the classical CTV method (requires the data)
```

Simulation uses seed 4321, 500 genotypes per class and α = 0.01 for the t and F tests that
retain simulated genotypes. The comparison with the classical CTV method uses α = 0.05.
Note that `ctv()` returns p-values (in %) only; the classification at a given α is done in
`R/CTV_Agreement.R`.

Note: `R/Adap_Stab_CTV.R` retrains a network and saves it under `Trained_Nets/`. Rename the
output file to avoid overwriting the network used in the article.

Packages: `RSNNS`, `tidyverse`, `readxl`, `writexl`, `caret`, `abind`, `gt`, `gtExtras`,
`flextable`, `NeuralNetTools`, `ggformula`.

## License

Code released under the MIT License (see `LICENSE`). If you use this code or the trained
network, please cite the article above.

## References

- Cruz CD, Torres RA and Vencovsky R (1989) An alternative approach to the stability analysis proposed by Silva and Barreto. Revista Brasileira de Genética 12: 567-580.
- Bergmeir C and Benítez JM (2012) Neural networks in R using the Stuttgart Neural Network Simulator: RSNNS. Journal of Statistical Software 46: 1-26.
