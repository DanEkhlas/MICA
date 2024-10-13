# MICA

MICA is a workflow for minimal inhibitory concentration (MIC) analysis. The goal of this workflow is to allow users the easy assessment of antimicrobial resistance profiles of their lab strains. More information on MICs can be found in the STAR protocol by [Barnes et al., 2023](https://www.sciencedirect.com/science/article/pii/S2666166723004793).

Currently the following assay formats are supported with this script:

- **Plate size:**                               96-well plates
- **Antimicrobial concentrations (columns):**   2-fold and customizable
- **Assay layout:**                             multiple strains & one antimicrobial / one strain & multiple antimicrobials

## Import of Absorbance Data from Spectrophotometer

In the first step, it will be required to load the actual absorbance data, obtained via a spectrophotometer after microbial growth. In this workflow, we assume that the spectrophotometer read-out will be from a 96-well plate.

For the purpose of this workflow example, I will provide the following example data, which can be loaded using the following command:

```
absorbance_data <- as.data(c(
  0.156, 0.200, 1.000, 0.950, 0.880, 0.930, 1.050, 1.020, 0.890, 0.980, 0.850, 0.060,  # Row A
  0.172, 0.123, 0.102, 0.088, 0.112, 0.200, 1.000, 1.030, 0.950, 0.900, 0.880, 0.087,  # Row B
  0.094, 0.132, 0.150, 0.970, 1.050, 0.980, 1.020, 0.910, 0.890, 0.880, 1.060, 0.043,  # Row C
  0.095, 0.103, 0.100, 0.980, 1.020, 0.980, 0.890, 0.870, 0.910, 0.950, 1.030, 0.101,  # Row D
  0.092, 0.095, 0.097, 0.094, 0.300, 1.030, 0.870, 0.910, 0.900, 1.050, 0.890, 0.098,  # Row E
  0.098, 0.250, 1.020, 1.000, 0.950, 1.010, 1.030, 0.940, 0.880, 0.870, 1.060, 0.093,  # Row F
  0.095, 1.020, 1.000, 0.940, 1.050, 0.900, 0.880, 0.920, 0.890, 0.980, 0.860, 0.081,  # Row G
  0.099, 0.093, 0.180, 0.910, 0.970, 0.880, 0.860, 0.930, 1.050, 0.890, 1.020, 0.078   # Row H
), nrow = 8, ncol = 12, byrow = TRUE)

df_plate <- as.data.frame(absorbance_data)
```