# Modeling Co-Infections of Tuberculosis with HIV

## Overview

This repository contains the code and data for the project **Modeling Co-Infections of Tuberculosis with HIV**. The project aims to develop a mathematical model that captures the interactions between Tuberculosis (TB) and HIV across different countries. We use a compartmental **SLIRS model** with separate compartments for HIV-positive and HIV-negative individuals to analyze the transmission dynamics of TB and its co-infections with HIV. The modeling approach is based on **Ordinary Differential Equations (ODE)**, and we perform **Bayesian parameter optimization** to estimate key parameters governing the disease dynamics. Refer to the complete **[Project Report](https://drive.google.com/file/d/1T0NR-2HJ2N8Zw2Wt9Q28wn-mKfUp_M8r/view?usp=sharing)** for more details.



## Key Objectives

- Develop a mathematical model that describes the transmission of TB and HIV.
- Compare the model results across different countries (Canada, Norway, Australia, UK, and Japan).
- Perform **sensitivity analysis** to identify key parameters that influence TB and HIV spread.
- Analyze **intervention strategies** to determine effective policies for disease control.

## Methodology

The **modeling approach** classifies the population into **Susceptible (S), Latent Infected (L), Infected (I)**, and **Recovered (R)** compartments. Each compartment is further divided into **HIV-positive** and **HIV-negative** subgroups, leading to a total of 8 compartments. Transmission, activation, recovery, and immunity loss rates are considered for both populations.

The **data sources** used in this pr**oject include datasets from the World Health Organization (WHO) on TB and HIV and statistics from Our World in Data. These sources provide reliable global health information necessary for our modeling and analysis.

For **parameter estimation**, we employ **Bayesian inference using Markov Chain Monte Carlo (MCMC) methods**. The model parameters are estimated separately for each country, allowing for an inter-country comparison of TB and HIV transmission dynamics.

The **sensitivity analysis** is performed using **Global Sensitivity Analysis with Latin Hypercube Sampling (LHS)**. This approach helps us evaluate the impact of key model parameters on TB spread and identify which parameters are most influential.

Our **intervention strategies** focus on two primary aspects: reducing the HIV transmission rate $\delta$  through awareness programs and reducing the TB transmission rate $\beta$ through improved healthcare policies. We simulate these interventions to assess their effectiveness in controlling TB and HIV spread.

## Figures

 ![Model Diagram](Figures/Model%20Diagram.png)
 *Figure 1. Model Diagram*

   

     

![EPV](Figures/Estimated%20Parameter%20Values%20(Log%20Scale).png)
*Figure 2: Estimated Parameter Values (Log Scale)*


![SA](Figures/Sesitivity%20Analysis%20(Log%20Scale,%20Australia).png)
*Figure 3: Sensitivity Analysis (Log Scale, Australia)*

![SIB](Figures/Intervention%20in%20Beta.png)
*Figure 4: Simulations of Intervention on Beta (Australia)*

![SID](Figures/Intervention%20in%20Delta.png)
*Figure 5: Simulations of Intervention on Delta (Australia)*

## Results

- The **TB transmission rate $\beta$** is the most critical factor in determining the spread of TB.
- **HIV-positive individuals have a higher TB transmission rate**, but their TB activation rate is lower due to medical care.
- Norway has **lower TB activation rates** for HIV-positive individuals but **higher TB transmission rates**.
- Japan has the **lowest TB transmission rates** despite its high population density.
- **UK has the highest proportion of latent TB cases**, suggesting a need for increased screening.
- **Reducing TB transmission $\beta$ is the most effective strategy**, whereas reducing HIV transmission $\delta$ has limited impact on TB dynamics.

## Limitations and Future Work

- **Assumptions:** The model assumes minimal interaction between HIV-positive and HIV-negative individuals, which may not always hold.
- **Parameter Estimation Uncertainty:** Some parameters show high variability across countries.
- **Future Enhancements:**
  - Extend the model to include **vaccination and treatment interventions**.
  - Perform **basic reproduction number $R_0$ analysis**.
  - Implement **Disease-Free Equilibrium (DFE) analysis** to assess eradication thresholds.

## Getting Started

### Requirements

- MATLAB

### Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/tb-hiv-model.git
cd tb-hiv-model
```

### Running the Model

The project contains three primary folders with MATLAB source code for performing the relevant operations:

```
├── InterventionAnalysis
├── ParameterEstimation
└── SensitivityAnalysis
```

Each folder contains MATLAB scripts that need to be opened and run in the MATLAB environment.

#### Running Parameter Estimation

```matlab
run('ParameterEstimation/parameterestimation.m')
```

#### Running Intervention Analysis

```matlab
run('InterventionAnalysis/InterventionAnalysis.m')
```

#### Running Sensitivity Analysis

```matlab
run('SensitivityAnalysis/SensitivityAnalysis.m')
```

Ensure that all dependencies are correctly set up before running the scripts.

## Data References

- [WHO Global Tuberculosis Report](https://www.who.int/teams/global-tuberculosis-programme/tb-reports/global-tuberculosis-report-2022)
- [Our World in Data - HIV Statistics](https://ourworldindata.org/hiv-aids)

## Acknowledgments

This project was developed in collaboration with **[Ali Parsaee](https://www.linkedin.com/in/ali-parsaee/)** and **[Rajkumar Patel](https://www.linkedin.com/in/rajkumarpatel96/)**, under the guidance of **Prof. Marie Varughese** at the University of Alberta.



