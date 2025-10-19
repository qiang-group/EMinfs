# EMinfs: Gamma frailty proportional hazards model for Regression analysis of current status data subject to informative censoring

This repository provides the R implementation for the manuscript  
**"Gamma frailty proportional hazards model for Regression analysis of current status data subject to informative censoring"**  
by *Qiang Zheng, Lianming Wang, Xiaoyan Lin, and Dengdeng Yu*.

---

## 🧭 Overview

**EMinfs** implements the algorithm of the **Gamma frailty proportional hazards (PH)** model for analyzing **current status data** subject to **mixed censoring** (both informative and non-informative).  

The method introduces a **three-step data-augmentation EM algorithm** using multinomial and Poisson latent variables.  
It provides:
- Closed-form variance estimation  
- Robust convergence and numerical stability  
- Fast computation and scalability  
- Flexible baseline hazard estimation via **monotone I-splines**

This approach extends traditional PH frailty models by introducing a latent Gamma variable that captures the dependence between failure and censoring times without requiring additional parameters.

---

## ⚙️ Installation

You can install the latest version directly from GitHub:

```r
devtools::install_github("https://github.com/qiang-group/EMinfs")
library(EMinfs)

# Load the dataset included in the package
data("chloroprene_data")

# Prepare data
d  <- chloroprene_data$Lesion.Found                  # 1 = tumor found, 0 = none
s  <- chloroprene_data$Informative.C.s               # 1 = informative censoring
Ci <- chloroprene_data$Days.On.Study                 # observed study time
Xp <- as.matrix(chloroprene_data[, c(
  "Group.12.8", "Group.30", "Group.80", "Gender.1.male.0.female"
)])

# Fit the Gamma-Frailty PH model using the developed EM algorithm
model_results <- estimate_em_model(d = d, s = s, Ci = Ci, Xp = Xp)

# Summarize results
print(model_results$coefficients_t)   # Effects on failure time (tumor onset)
print(model_results$coefficients_c)   # Effects on censoring time
