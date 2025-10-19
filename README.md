# EMinfs: Gamma Frailty Proportional Hazards Model for Current Status Data

This repository provides the R implementation for the manuscript  
**"Regression Analysis of Current Status Data Under Informative Censoring"**  
by *Qiang Zheng, Lianming Wang, Xiaoyan Lin, and Dengdeng Yu*.

---

## 🧭 Overview

**EMinfs** implements a **Gamma frailty proportional hazards (PH)** model for analyzing **current status data** subject to **mixed censoring** (both informative and non-informative).  

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
