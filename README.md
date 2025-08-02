# Robust Estimation of Change Points in Linear Spline Models with Missing Data: Simulation and Real Data Application

This repository contains R code for change point estimation in linear spline models with potentially missing outcomes. The code implements and compares three approaches:

- **Complete Case (CC)**
- **Inverse Probability Weighting (IPW)**
- **Doubly Robust-Augmented Inverse Probability Weighting (DR-AIPW)**


---

## Contents

- `Change_point.R`:  
  Contains code for simulation studies with user-defined sample sizes, change point locations, and other parameters. It outputs results including bias, variance estimate, Monte-Carlo variance and confidence interval coverage of change point estimate for each method.

- `ACTG_application.R`:  
  Applies the above methods to the **ACTG175** clinical trial dataset, estimating the change point in CD4 counts at 96 weeks versus CD4 counts at baseline. The estimate, variance estimate and confidence interval of change point for CC, IPW and DR-AIPW methods are recorded. The resulting figure shows estimated piecewise linear curves and confidence intervals for CC and DR-AIPW methods.

- `LICENSE`:  
  MIT License.

---

## How to Run

Make sure you have installed the following required R packages:

```r
install.packages(c("ggplot2", "foreach", "parallel","dplyr", "speff2trial"))
```

Then source and run either of the main scripts:

```r
source("Change_point.R")       # For simulation study
source("ACTG_application.R")   # For real data application
```

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
