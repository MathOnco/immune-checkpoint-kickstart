# immune-checkpoint-kickstart
The code (MATLAB) for the paper "[The Immune Checkpoint Kick Start: Optimization of Neoadjuvant Combination Therapy Using Game Theory](https://ascopubs.org/doi/full/10.1200/CCI.18.00078)" published in JCO Clinical Cancer Informatics.

## Figure functions:
Included in the repository are the MATLAB code used to make the following figures:
- figure 3A
- figure 3B
- figure 4
- figure 5A
- figure 5B

## Helper functions:
Used inside the above functions are several helper functions, listed below.
- calcMetScore.m (calculates time-weighted average of met score, in fig 5)
- gradientOfMetPotential.m (calculates the map of metastatic potential, used as the background color of fig 3B)
- nice_plot.m (resizes fonts for line plots)
- payoff.m (returns the payoff matrix for given parameter set)
- rep_ode (the set of replicator ordinary differential equations; typically solved using ode45.m)
