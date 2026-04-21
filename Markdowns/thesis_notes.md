---
title: 'Thesis results'
author: "Manuel Rivera"
date: "2026-04-19"
---

# Description
This markdown is for summarizing all principal results about research. Plots are up to date and data, script and plot paths will be included.

# Data Analysis
First step for generating data was analysis the parameter space available to simulate. 

## Interaction matrix
To analyze how to generate the interaction matrix required for the gLV model, 
we simulated several communities at different sizes, sparsity levels, and 
proportions of positive and negative interactions. Results suggested that the 
model successfully simulated communities in which all off-diagonal interactions 
were set to values <= 0.

The diagonal of the interaction matrix was always set to −0.5, while 
off-diagonal elements were assigned negative or null interactions with 
proportions *p_neg* and *p_null*, respectively.

![Figure 1: Parameter analysis](./PEA/params_space_heatmap.png)

A sensitivity analysis of the solver tolerances showed no effect on model 
performance among simulations that converged successfully. Cases that failed 
to run under a given tolerance setting were excluded from this analysis.

![Figure 2: Solver tolerance analysis](./PEA/tol_analysis/tolerance_heatmap.png)

# Control generation
Once we identified which cases were succesful to run we tried several approaches to generate in silico communities with a specie of the community being the keystone of it. Aproaches tested were:

## Keystone by cascading effect
In silico communities were generated in which one species was selected at random to initiate the cascading effect. 
For example, in a community of 5 species where Species 1 affects Species 2, Species 2 affects Species 3, and 
so on, the resulting interaction chain would be:

S1 -> S2 -> S3 -> S4 -> S5

The interaction matrix was set as follows: diagonal elements were fixed at −0.5, while off-diagonal elements maintained a defined proportion of null and negative interactions, with negative values drawn from −U(0,1). 
Off-diagonal cascade elements did not contribute to this proportion and were instead drawn from −U(0,1) × 10, reflecting the stronger effect of cascading interactions.

Data path: `/mnt/data/sur/users/mrivera/Data/Controls/Cascade_keystone`

![Figure 3: Metrics distributions]( "./Controls/Cascade_keystone/all_metrics_distribution.png")

![Figure 4: Relative abundance and metrics]( "./Controls/Cascade_keystone/keystone_relative_metrics.png")

![Figure 5: Metric distributions for survival species]( "./Controls/Cascade_keystone/alive_metrics_distribution.png")

## Different magnitude keystone column (neg_mag)
This experiment was designed to generate control communities that explicitly include a keystone (k) species. The interaction matrix was set as follows: diagonal elements were fixed at −0.5, while off-diagonal non-zero elements excluding column k were drawn from −U(1,2). Column k was drawn from -U(0,1), representing less strong interactions from the keystone species on the rest of the community.

The rationale for this definition was that removal of the keystone species k would destabilise the community by releasing previously buffered strong interactions, allowing them to dominate community dynamics.

Data path: `/mnt/data/sur/users/mrivera/Data/clean_controls/Negative_magnitude_d1e776c23c74`

![Figure 6: Metrics distributions]("./clean_controls/Negative_magnitude_d1e776c23c74/all_metrics_distribution.png")

![Figure 7: Keystone species relative abundance and metrics]("./clean_controls/Negative_magnitude_d1e776c23c74/keystone_relative_metrics.png")

![Figure 8: Metric distributions for survival species]("./clean_controls/Negative_magnitude_d1e776c23c74/alive_metrics_distribution.png")

## Boosting keystone column (Kboost)
This experiment was designed to generate control communities that explicitly include a keystone (k) species. 
The interaction matrix was set as follows: diagonal elements were fixed at −0.5, off-diagonal had a proportion of being zero or otherwise negative drawn from a -U(0,1). Column k was boosted by a constant of 10, excepting by interaction with itself.

This was operationalised as follows: species k's effect on other species was either zero or negative, with non-zero values drawn from −U(0,1) and subsequently multiplied by 10, making the keystone species the strongest interactor in the system. Consequently, its removal would destabilise the community by allowing those strong interactions to go unchecked.

Data path: `/mnt/data/sur/users/mrivera/Data/Controls/Boosted_keystone`

![Figure 9: Metrics distributions]("./Controls/Boosted_keystone/all_metrics_distribution.png")

![Figure 10: Keystone species relative abundance and metrics]("./Controls/Boosted_keystone/keystone_relative_metrics.png")

![Figure 11: Metric distributions for survival species]("./Controls/Boosted_keystone/alive_metrics_distribution.png")

## Impact quantification
Community-level impact was assessed for both the full community and a sub-community restricted to surviving species. Relative abundance, dissimilarity, and keystoneness were computed independently for each case.

Data path: `/mnt/data/sur/users/mrivera/Data/clean_controls/ImpactAnalysis_9ee93e`

![Figure 12: Metric distributions for survival species with full-community extinctions impact]("./clean_controls/ImpactAnalysis_9ee93e/full_metrics_distribution.png")

![Figure 13: Metric distributions for survival species with sub-community extinctions impact]("./clean_controls/ImpactAnalysis_9ee93e/sub_metrics_distribution.png")