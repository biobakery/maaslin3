context("collapse migration regressions")

# Regressions introduced by the dplyr -> collapse migration, plus a
# separately-discovered duplicate-feature-name bug in maaslin_read_data.

test_that("maaslin_process_metadata handles all-categorical metadata", {
    # collapse::num_vars(metadata) returns an empty data frame when there
    # are no numeric columns, and collapse::fscale on that empty frame
    # previously triggered:
    #   Error in `get_vars_ind<-`(x, .Call(C_vtypes, x, 1L), value) :
    #     NROW(value) must match nrow(x)
    metadata <- data.frame(
        treatment = c("A", "B", "A", "B", "A", "B"),
        group = c("X", "Y", "X", "Y", "X", "Y"),
        stringsAsFactors = FALSE,
        row.names = paste0("S", 1:6)
    )

    result <- maaslin3:::maaslin_process_metadata(
        metadata = metadata,
        fixed_effects = c("treatment", "group"),
        reference = NULL,
        feature_specific_covariate_name = NULL,
        standardize = TRUE
    )

    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 6)
    expect_true(all(c("treatment", "group") %in% colnames(result)))
})

test_that("maaslin_process_metadata standardizes numeric and preserves categorical", {
    metadata <- data.frame(
        age = c(25, 30, 35, 40, 45, 50),
        treatment = c("A", "B", "A", "B", "A", "B"),
        stringsAsFactors = FALSE,
        row.names = paste0("S", 1:6)
    )

    result <- maaslin3:::maaslin_process_metadata(
        metadata = metadata,
        fixed_effects = c("age", "treatment"),
        reference = NULL,
        feature_specific_covariate_name = NULL,
        standardize = TRUE
    )

    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 6)
    # numeric column should be z-scored (mean ~0, sd ~1)
    expect_equal(mean(result$age), 0, tolerance = 1e-8)
    expect_equal(sd(result$age), 1, tolerance = 1e-8)
    # categorical column unchanged
    expect_equal(as.character(result$treatment), metadata$treatment)
})

test_that("collapse::fscale on empty data frame is guarded", {
    # Document the underlying collapse behavior that motivated the guard.
    empty_df <- data.frame(row.names = paste0("S", 1:6))
    expect_equal(NCOL(empty_df), 0)
    # Without the guard, assigning fscale(empty_df) back into num_vars<-
    # errors. The maaslin_process_metadata code must therefore skip the
    # assignment when there are no numeric columns.
})

test_that("maaslin_read_data makes duplicate feature names unique", {
    tmpdir <- tempfile("maaslin_dup_")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

    data_path <- file.path(tmpdir, "data.tsv")
    meta_path <- file.path(tmpdir, "meta.tsv")

    set.seed(1)
    n_samples <- 10
    data_mat <- matrix(
        abs(rnorm(n_samples * 4, mean = 100, sd = 30)),
        nrow = n_samples,
        ncol = 4
    )
    data_df <- as.data.frame(data_mat)
    # Two duplicate column names - exactly the condition that blew up
    # downstream subsetting with "undefined columns selected".
    colnames(data_df) <- c("Bacteroides", "Prevotella", "Bacteroides", "Firmicutes")
    rownames(data_df) <- paste0("S", 1:n_samples)

    # write.table would mangle dup names; use a manual writer to preserve them.
    header <- paste(c("", colnames(data_df)), collapse = "\t")
    body <- apply(data_df, 1, function(row) paste(row, collapse = "\t"))
    writeLines(c(header, paste(rownames(data_df), body, sep = "\t")), data_path)

    meta_df <- data.frame(
        treatment = rep(c("A", "B"), length.out = n_samples),
        row.names = paste0("S", 1:n_samples),
        stringsAsFactors = FALSE
    )
    write.table(meta_df, meta_path, sep = "\t", quote = FALSE,
                col.names = NA, row.names = TRUE)

    result <- maaslin3:::maaslin_read_data(
        input_data = data_path,
        input_metadata = meta_path,
        feature_specific_covariate = NULL,
        unscaled_abundance = NULL
    )

    expect_true(!any(duplicated(colnames(result$data))))
    expect_true("Bacteroides" %in% colnames(result$data))
    expect_true("Bacteroides.1" %in% colnames(result$data))
})

test_that("add_qvals handles prevalence/abundance fits with all-NA pvals", {
    # data.table::fifelse is strict about yes/no types. When a fit object's
    # $results$pval is all-NA (logical), fifelse(..., NA_real_, pval)
    # previously errored with
    #   'no' is of type logical but 'yes' is double
    # The fix coerces pvals to numeric first.
    fake_fit <- list(
        results = data.frame(
            pval = NA,
            error = NA,
            stringsAsFactors = FALSE
        )
    )

    result <- maaslin3:::add_qvals(
        fit_data_abundance = fake_fit,
        fit_data_prevalence = fake_fit,
        correction = "BH"
    )

    expect_true(is.list(result))
    expect_true(all(is.na(result$fit_data_abundance$results$qval_individual)))
    expect_true(all(is.na(result$fit_data_prevalence$results$qval_individual)))
})

test_that("full maaslin3 run works with only categorical metadata", {
    tmpdir <- tempfile("maaslin_cat_only_")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

    set.seed(42)
    n_samples <- 20
    n_features <- 8
    data_mat <- matrix(
        abs(rnorm(n_samples * n_features, mean = 100, sd = 30)),
        nrow = n_samples,
        ncol = n_features
    )
    data_df <- as.data.frame(data_mat)
    colnames(data_df) <- paste0("Feature", 1:n_features)
    rownames(data_df) <- paste0("S", 1:n_samples)

    meta_df <- data.frame(
        treatment = rep(c("A", "B"), each = n_samples / 2),
        row.names = paste0("S", 1:n_samples),
        stringsAsFactors = FALSE
    )

    out_dir <- file.path(tmpdir, "out")

    expect_error(
        suppressWarnings(maaslin3::maaslin3(
            input_data = data_df,
            input_metadata = meta_df,
            output = out_dir,
            fixed_effects = c("treatment"),
            normalization = "TSS",
            transform = "LOG",
            standardize = TRUE,
            plot_summary_plot = FALSE,
            plot_associations = FALSE,
            save_models = FALSE,
            verbosity = "ERROR"
        )),
        NA
    )
})
