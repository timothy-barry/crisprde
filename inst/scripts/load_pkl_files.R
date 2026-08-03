crispritz_score_dir <- "/Users/timbarry/research_code/CRISPRitz/sourceCode/Python_Scripts/Scores"
mismatch_score <- reticulate::py_load_object(file.path(crispritz_score_dir, "mismatch_score.pkl"))
mismatch_score_v <- as.numeric(mismatch_score) |> setNames(names(mismatch_score))
pam_scores <- reticulate::py_load_object(file.path(crispritz_score_dir, "PAM_scores.pkl"))
pam_scores_v <- as.numeric(pam_scores) |> setNames(names(pam_scores))

extdata_dir <- file.path("inst", "extdata")
dir.create(extdata_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(mismatch_score_v, file.path(extdata_dir, "cfd_mismatch_scores.rds"))
saveRDS(pam_scores_v, file.path(extdata_dir, "cfd_pam_scores.rds"))
