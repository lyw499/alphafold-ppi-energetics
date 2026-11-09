# Exploratory same-species BSA-adjusted comparison



The audited reference contains 23 M. genitalium–M. genitalium chain pairs from 7 PDB entries and 8 assemblies. Mixed-species and human-only pairs were excluded from this primary scope, but retained in the audit.

Using positive BSA and common support, a parsimonious NB2 model with dataset-specific log1p(BSA) slopes estimated a Pedro/experimental standardized mean-count ratio of 0.748 (approximate PDB-cluster-t 95% CI 0.545–1.025; exploratory p=0.0652). The comparison included 37069 Pedro observations and 21 reference interfaces from 7 PDB entries.

The original FoldBench and SAFES estimates were 0.548 and 0.528. The Mgen point estimate is closer to one than both original estimates, descriptively; this is not a formal between-comparison contrast.

Leave-one-PDB-out point estimates ranged from 0.653 to 0.869.

With very few independent reference PDBs, clustered intervals and p-values remain approximate and potentially overconfident. Sensitivities, low-resolution entries, shared protein dependence and changing standardization targets should be considered. A nonsignificant result would not demonstrate equivalence, and these unpaired comparisons cannot establish AF3 prediction error or a causal species effect.

Additional diagnostics: the QC-adjusted fit failed its Hessian check and is not reported as a valid estimate. The reference Pearson variance diagnostic was approximately 3.19, indicating residual heterogeneity under the shared-dispersion model. Excluding 8PBZ/8PC1/8PC0 gave a point estimate of 0.622 but left only 13 interfaces from four PDBs; PDB-equal weighting gave 0.767. These sensitivities should not replace the primary result.
