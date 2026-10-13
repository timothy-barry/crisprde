# crisprde

Analysis tools for Caliper-GW and amplicon sequencing.

- Caliper-GW tuning, NB/NB-mixture models, occupancy models, and homology weighting.
- GUIDE-seq clustering, annotation, replicate count tables, and UMI deduplication.
- Weighted/unweighted QQ plots, pipeline QC plots, and local UMI plots.
- Frequentist and Bayesian amplicon inference and sample-size planning.
- Multivariate amplicon processing and modeling in `R/multivariate_amplicon_seq.R`.
- Base-editing count-table import, simulations, and synthetic data generation.
- CFD lookup-table generation in `inst/scripts/load_pkl_files.R`.

Multivariate Bayesian fitting uses CmdStanR and the bundled
`inst/stan/reg_model.stan` source; the compiled executable is generated locally.

## Cleanup compatibility

Legacy DM GUIDE-seq methods, old alignment utilities, all amplicon plotting
functions, and superseded GUIDE-seq plots have been removed. Historical analyses
using those functions can use pre-cleanup commit
`56cef6800bb1ba76fc6d11eecc904c5fcd03d02e`.

Retained plots are organized in `plotting_guideseq_results.R`, `weighted_qq_plot.R`,
`qq_plot_helpers.R`, `plotting_guideseq_qc.R`, and `plotting_local_umi.R`.
