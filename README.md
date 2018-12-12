Replicating the results for the data application found in "Causal inference with interfering units for cluster and population level treatment allocation programs" by Papadogeorgou G, Mealli F, Zilger C

```Data_analysis.R```: The main code to replicate the results of the paper. In this script, we load the data set, form 50 clusters based on Ward's agglomerative method, and estimate population average potential outcomes, direct and indirect effects.

```Sensitivity_analysis.R```: Estimating the population average potential outcomes, direct and indirect effects using alternative clustering approaches.

```Data_analysis_Fa_indirect.R```: Estimating the population level intervention indirect effect, where the distribution of treatment cluster coverage shifts from within the 20-80th percentiles to within the 50-80th percentiles of the observed cluster coverage distribution.
