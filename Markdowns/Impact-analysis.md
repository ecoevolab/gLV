---
title: "Experiment *91074c4e25b4* notes"
author: "Manuel Rivera"
date: "2026-03-22"
---

# Description
Previous GNN models trained on extinctions of survival nodes showed a negative correlation, 
suggesting that computing impact at the full community level is critical. To test this hypothesis, 
we simulated the extinction of survival species only and compared two impact calculations: 
one considering all nodes in the community, and one restricted to survival nodes only.

## Impact of extinctions at full community
Distribution of impact metrics when a species is removed from the full community, which includes all species regardless of whether they survived.

![Figure 1: Metrics distribution at full impact](/mnt/data/sur/users/mrivera/Plots/91074c4e25b4/full_impact_distr.png)

## Impact of extinctions at sub-community
Distribution of impact metrics when a species is removed from the sub-community, which includes only the species that survived.

![Figure 2: Metrics distribution at sub-community impact](/mnt/data/sur/users/mrivera/Plots/91074c4e25b4/sub_impact_distr.png)

