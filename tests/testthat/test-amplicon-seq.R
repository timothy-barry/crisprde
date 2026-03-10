testthat::test_that("low-information amplicons remain in the output but skip Wald inference", {
  data_list <- list(
    amplicon_ids = c("boundary", "lowinfo"),
    n_mat_trt = matrix(c(100, 100, 100, 100, 100, 100), nrow = 3),
    n_mat_cntrl = matrix(c(100, 100, 100, 100, 100, 100), nrow = 3),
    k_mat_trt = matrix(c(10, 12, 11, 1, 0, 0), nrow = 3),
    k_mat_cntrl = matrix(c(0, 0, 0, 0, 0, 0), nrow = 3)
  )

  res <- run_amplicon_seq_analysis(data_list, min_mutation_count = 20L)
  result_df <- res$result_df

  testthat::expect_equal(result_df$amplicon_id, c("boundary", "lowinfo"))
  testthat::expect_equal(result_df$total_mutated_reads, c(33, 1))
  testthat::expect_equal(result_df$passes_mutation_count_qc, c(TRUE, FALSE))

  lowinfo_row <- result_df[result_df$amplicon_id == "lowinfo", ]
  testthat::expect_false(is.na(lowinfo_row$theta_hat))
  testthat::expect_true(is.na(lowinfo_row$theta_hat_se))
  testthat::expect_true(is.na(lowinfo_row$theta_hat_lower_ci))
  testthat::expect_true(is.na(lowinfo_row$theta_hat_upper_ci))
  testthat::expect_true(is.na(lowinfo_row$p_value))
  testthat::expect_false(lowinfo_row$significant)
  testthat::expect_true(is.na(lowinfo_row$rho_hat))
  testthat::expect_true(is.na(lowinfo_row$pilot_rho_hat))
})


testthat::test_that("Jeffreys regularization yields finite Wald inference at the boundary", {
  data_list <- list(
    amplicon_ids = "boundary",
    n_mat_trt = matrix(c(100, 100, 100), nrow = 3),
    n_mat_cntrl = matrix(c(100, 100, 100), nrow = 3),
    k_mat_trt = matrix(c(10, 12, 11), nrow = 3),
    k_mat_cntrl = matrix(c(0, 0, 0), nrow = 3)
  )

  res <- run_amplicon_seq_analysis(data_list, min_mutation_count = 20L)
  result_df <- res$result_df

  testthat::expect_true(result_df$passes_mutation_count_qc)
  testthat::expect_gt(result_df$theta_hat_se, 0)
  testthat::expect_false(is.na(result_df$theta_hat_lower_ci))
  testthat::expect_false(is.na(result_df$theta_hat_upper_ci))
  testthat::expect_false(is.na(result_df$p_value))
  testthat::expect_equal(res$dispersion_diagnostics$shared_rho_hat, 0)
})


testthat::test_that("amplicon IDs fall back to matrix column names when needed", {
  data_list <- list(
    n_mat_trt = matrix(c(100, 100, 100, 100, 100, 100), nrow = 3,
                       dimnames = list(NULL, c("amp_a", "amp_b"))),
    n_mat_cntrl = matrix(c(100, 100, 100, 100, 100, 100), nrow = 3,
                         dimnames = list(NULL, c("amp_a", "amp_b"))),
    k_mat_trt = matrix(c(5, 6, 4, 0, 0, 0), nrow = 3,
                       dimnames = list(NULL, c("amp_a", "amp_b"))),
    k_mat_cntrl = matrix(c(0, 0, 0, 0, 0, 0), nrow = 3,
                         dimnames = list(NULL, c("amp_a", "amp_b")))
  )

  fisher_res <- run_fisher_exact_test(data_list)
  analysis_res <- run_amplicon_seq_analysis(data_list, min_mutation_count = 20L)

  testthat::expect_equal(fisher_res$amplicon_id, c("amp_a", "amp_b"))
  testthat::expect_equal(analysis_res$result_df$amplicon_id, c("amp_a", "amp_b"))
})
