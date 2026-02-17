# Metadiet

# Supplementary Data and Code

## Fecal trnL Metabarcoding Reliably Reconstructs Plant Intake and Reveals Multi-Omics Links to Gut Microbiota and Host Metabolism

**Status:** Submitted to *Food & Function* **Authors:** 
[Nubia Andrea Villota-Salazar1#*, Oscar J. Lara-Guzmán1,2#, Vanessa Corrales-Agudelo1,2, Diego A. Rivera1, O. Lucía Ortega1,3, Jelver A. Sierra1,3, Katalina Muñoz-Durango1,3, Juan S. Escobar1,3*

(#)These authors contributed equally to this work. 
1 Vidarium–Nutrition, Health, and Wellness Research Center, Grupo Empresarial Nutresa. Carrera 52 #2-38, 050023 Medellin, Colombia.
2 Current address: School of Nutrition and Dietetics, Universidad de Antioquia, Calle 70 No. 52-21, 050010 Medellín, Colombia.
3 Current address: Faculty of Medicine, Universidad de Antioquia, Carrera 51D #62-29, 050010 Medellin, Colombia.]  

**Contact:** [* Corresponding author:
Nubia Andrea Villota-Salazar, andreavillota172024@outlook.com]

---

### 📌 Overview

This repository contains the source code, statistical scripts, and processed datasets required to reproduce the analysis presented in the manuscript **"Fecal trnL Metabarcoding Reliably Reconstructs Plant Intake and Reveals Multi-Omics Links to Gut Microbiota and Host Metabolism"**.

The core of this analysis focuses on the multi-omics integration of **plant intake (trnL metabarcoding)**, **gut microbiota (16S rRNA)**, and **plasma metabolome** using the **DIABLO** framework from the `mixOmics` R package.

### 📂 Repository Structure

```text
├── data/
│   ├── raw/                 # Raw count tables (trnL, 16S, metabolites)
│   ├── processed/           # CLR and Z-score transformed data
│   └── metadata.csv         # Participant metadata (T1/T2, Groups)
├── scripts/
│   ├── 1_preprocessing.R         # Normalization (CLR) & Alpha/Beta diversity
│   ├─────  1.1_correlations.R    # Procrustes & sPLS covariation analysis
│   ├── 2_diablo_main.R           # Normalization (Z-score) & CORE SCRIPT: Multi-omics DIABLO integration
│   └─────  2.1_pairwise.R        # univariate correlations & Biomarker validation
├── results/                      # Output figures and tables
└── README.md
