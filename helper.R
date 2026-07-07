library(tidyverse)

##################################

exclude_collinear_vars <- function(
  pred_matrix,
  corr_mat,
  lower_threshold = 0.1,
  upper_threshold = 0.99
) {
  # That's the function of mice for removing collinear vars
  # https://github.com/amices/mice/blob/7df3487a56cd49dffc9192a8ebe9697bfa7258ca/R/internal.R#L84-L91
  # Inclusion of low correlated varvs is not useful
  # See: Auxiliary_variables_in_multiple_imputation_in_regr
  # See also a practical guide
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC8017730/
  # It suggests keeping auxiliary vars with cors >=0.5

  # Remove diagonal pairs (self-correlated variables)
  diag(corr_mat) <- 0

  # Just use the lower triangle of the matrix
  remove_idx <- lower.tri(corr_mat) &
    !is.na(corr_mat) &
    (abs(corr_mat) > upper_threshold |
      abs(corr_mat) < lower_threshold)

  idx <- which(remove_idx, arr.ind = TRUE)

  removed_pairs <- data.frame(
    row_var = rownames(corr_mat)[idx[, 1]],
    col_var = colnames(corr_mat)[idx[, 2]],
    corr = corr_mat[remove_idx]
  )
  message(
    glue::glue(
      "Removing; `{NROW(removed_pairs)}` pair(s) "
    )
  )

  # Remove from lower triangle
  pred_matrix[remove_idx] <- 0

  # Mirror to upper triangle
  pred_matrix[upper.tri(pred_matrix)] <-
    t(pred_matrix)[upper.tri(pred_matrix)]

  # Attach the removed pairs as an attribute!
  attr(pred_matrix, "removed_pairs") <- removed_pairs

  return(pred_matrix)
}


test_df <- data.frame(
  age_child_12_1 = c(0, 1, 1, 1, 1, 1),
  age_parent_14 = c(1, 0, 1, 1, 1, 1),
  age_teach_14_1 = c(1, 1, 0, 1, 1, 1),
  age_child_14_1 = c(1, 1, 1, 0, 1, 1),
  # age_child_16_1 & age_parent_16 are highly correlated
  age_child_16_1 = c(1, 1, 1, 1, 0, 0),
  age_parent_16 = c(1, 1, 1, 1, 0, 0)
)

predMatrix <- mice::make.predictorMatrix(data = test_df)

corr_mat <- cor(test_df, method = "spearman", use = "pairwise.complete.obs")
# ΝΑ are produced because one var does not vary when their pair does, so
# their SD ==0 and the correlation is NA

corr_mat <- as.data.frame(corr_mat)
s <- exclude_collinear_vars(
  pred_matrix = predMatrix,
  corr_mat = cor(test_df, use = "pairwise.complete.obs")
)


test_s <- data.frame(
  age_child_12_1 = c(0, 1, 1, 1, 1, 1),
  age_parent_14 = c(1, 0, 1, 1, 1, 1),
  age_teach_14_1 = c(1, 1, 0, 1, 1, 1),
  age_child_14_1 = c(1, 1, 1, 0, 1, 1),
  # age_child_16_1 & age_parent_16 are highly correlated
  age_child_16_1 = c(1, 1, 1, 1, 0, 0),
  age_parent_16 = c(1, 1, 1, 1, 0, 0)
)
names(test_s) <- rownames(predMatrix)
rownames(test_s) <- rownames(predMatrix)
test_s <- as.matrix(test_s)
attr(test_s, "removed_pairs") <- data.frame(
  row_var = "age_parent_16",
  col_var = "age_child_16_1",
  corr = 1
)
stopifnot(
  all.equal(
    s,
    test_s
  )
)
rm(list = c("test_s", "s", "corr_mat", "predMatrix", "test_df"))

####################################

modify_pred_matrix_scales <- function(
  pred_matrix,
  item_pattern,
  total_pattern = "total"
) {
  vars <- colnames(pred_matrix)

  item_vars <- vars[grepl(item_pattern, vars)]
  message(glue::glue("Found `{length(item_vars)}`"))
  message(paste(item_vars, collapse = " "))
  total_vars <- vars[grepl(total_pattern, vars)]

  message(glue::glue("Found `{length(total_vars)}`"))
  message(paste(total_vars, collapse = " "))
  non_scale_vars <- setdiff(vars, c(item_vars, total_vars))

  non_scale_vars_items <- non_scale_vars[grepl(
    pattern = "item",
    non_scale_vars
  )]
  non_scale_vars_totals <- non_scale_vars[grepl(
    pattern = "total",
    non_scale_vars
  )]
  non_scale_vars_other <- setdiff(
    non_scale_vars,
    c(non_scale_vars_items, non_scale_vars_totals)
  )

  # See Van Buuren p.181 #

  ####################################################################
  # 1.
  # Impute variables that are NOT items or totals from any scale given
  # totals and other variables (not items!)
  ####################################################################
  pred_matrix[non_scale_vars_other, ] <- 0
  pred_matrix[non_scale_vars_other, total_vars] <- 1
  pred_matrix[non_scale_vars_other, non_scale_vars_totals] <- 1
  pred_matrix[non_scale_vars_other, non_scale_vars_other] <- 1

  #######################################################
  # 2
  # Items from a scale can predict ONLY their co-items.
  # And they can not be predicted by other items
  # and their own total!
  # Thus, items:
  #    can be predicted by:
  #      - non-scale vars
  #      - other items
  #      - other totals!
  #######################################################

  pred_matrix[, item_vars] <- 0

  # Items impute given the scale items
  pred_matrix[item_vars, item_vars] <- 1

  # Items can not be predicted by other scales' items
  pred_matrix[item_vars, non_scale_vars_items] <- 0

  # Do not use their own total
  pred_matrix[item_vars, total_vars] <- 0

  #######################################################
  # 3. Totals are passive:
  #    not predicted by anything
  #######################################################

  pred_matrix[total_vars, ] <- 0
  # The rest of `non_scale_vars_totals` will be handled seperately
  # It is not necessary to handle them in this iteration
  # pred_matrix[non_scale_vars_totals, ] <- 0

  # IMPORTANT #
  # if there are any subtotals, pls modify the predMatrix accordingly
  # Subtotals should not predict or be predicted by anything, because
  # they are linear transformation of other variables and collinearity issues
  # will arise if you leave both the subtotals and the scale total.

  # ------------------------------------------------------
  # 4. Non-scale vars:
  #    predicted only by non-scale vars + totals
  # ------------------------------------------------------

  pred_matrix[non_scale_vars_other, item_vars] <- 0
  # The rest of `non_scale_vars_items` will be handled seperately
  # It is not necessary to handle them in this iteration
  # pred_matrix[non_scale_vars_other, non_scale_vars_items] <- 0

  # ------------------------------------------------------
  # No self-prediction
  # ------------------------------------------------------

  diag(pred_matrix) <- 0

  return(pred_matrix)
}

test_df <- data.frame(
  PANSS1p_item1 = 1,
  PANSS1p_item2 = 1,
  PANSS1n_item1 = 1,
  PANSS1g_item1 = 1,
  PANSS1p.total = 1,
  PANSS1n.total = 1,
  PANSS1g.total = 1,
  PANSS1.total = 1,
  PANSS2p_item1 = 1,
  PANSS2p_item2 = 1,
  PANSS2n_item1 = 1,
  PANSS2g_item1 = 1,
  PANSS2p.total = 1,
  PANSS2n.total = 1,
  PANSS2g.total = 1,
  PANSS2.total = 1,
  age = 50
)
pred <- mice::make.predictorMatrix(test_df)
pred <- modify_pred_matrix_scales(
  pred,
  item_pattern = "^PANSS1[png]_item",
  total_pattern = "^PANSS1[png]?.*total"
) %>%
  modify_pred_matrix_scales(
    item_pattern = "^PANSS2[png]_item",
    total_pattern = "^PANSS2[png]?.*total"
  ) %>%
  modify_pred_matrix_scales(
    item_pattern = "mpvs_item",
    total_pattern = "mpvs_total"
  )

testthat::test_that("items can predict within wave but not self", {
  items_1 <- grep("^PANSS1", colnames(pred), value = TRUE)
  items_1 <- items_1[grepl("[png]_item\\d+$", items_1)]

  sub <- as.matrix(pred[items_1, items_1])

  testthat::expect_true(all(sub[lower.tri(sub) | upper.tri(sub)] == 1))
  testthat::expect_true(all(diag(sub) == 0))
})

testthat::test_that("same-wave totals do not predict items", {
  # Wave 1
  testthat::expect_true(
    all(
      pred[
        grepl("^PANSS1p|^PANSS1n|^PANSS1g", rownames(pred)),
        grepl("^PANSS1.*total", colnames(pred))
      ] ==
        0
    )
  )
  # Wave 2
  testthat::expect_true(
    all(
      pred[
        grepl("^PANSS2p|^PANSS2n|^PANSS2g", rownames(pred)),
        grepl("^PANSS2.*total", colnames(pred))
      ] ==
        0
    )
  )
})

testthat::test_that("totals are not predicted by any variable", {
  total_rows <- pred[grepl("total", rownames(pred)), ]

  testthat::expect_true(all(total_rows == 0))
})


testthat::test_that("age is predicted by all totals and items only by other-wave totals", {
  total_cols <- grep("total", colnames(pred), value = TRUE)

  wave1_totals <- grep("^PANSS1.*total", total_cols, value = TRUE)
  wave2_totals <- grep("^PANSS2.*total", total_cols, value = TRUE)

  wave1_items <- grep("^PANSS1[png]\\d+", rownames(pred), value = TRUE)
  wave2_items <- grep("^PANSS2[png]\\d+", rownames(pred), value = TRUE)

  # Age predicted by all totals
  testthat::expect_true(
    all(pred["age", total_cols] == 1)
  )

  # Wave 1 items predicted by Wave 2 totals
  testthat::expect_true(
    all(pred[wave1_items, wave2_totals] == 1)
  )

  # Wave 1 items NOT predicted by Wave 1 totals
  testthat::expect_true(
    all(pred[wave1_items, wave1_totals] == 0)
  )

  # Wave 2 items predicted by Wave 1 totals
  testthat::expect_true(
    all(pred[wave2_items, wave1_totals] == 1)
  )

  # Wave 2 items NOT predicted by Wave 2 totals
  testthat::expect_true(
    all(pred[wave2_items, wave2_totals] == 0)
  )
})

testthat::test_that("items can be predicted by other wave totals", {
  testthat::expect_true(
    all(
      pred[
        grepl("^PANSS1", rownames(pred)) &
          grepl("[png]\\d+", rownames(pred)),
        grepl("^PANSS2.*total", colnames(pred))
      ] ==
        1
    )
  )

  testthat::expect_true(
    all(
      pred[
        grepl("^PANSS2", rownames(pred)) &
          grepl("[png]\\d+", rownames(pred)),
        grepl("^PANSS1.*total", colnames(pred))
      ] ==
        1
    )
  )
})

testthat::test_that("predictor matrix is valid binary matrix", {
  testthat::expect_true(all(pred %in% c(0, 1)))
})

rm(test_df)
rm(pred)
#######################################


color_corr_residuals <- function(resid_df, limit, bg_color) {
  ft <- flextable::flextable(data = resid_df) |>
    flextable::bg(
      # https://stackoverflow.com/a/11672966
      j = base::setdiff(
        x = colnames(resid_df),
        y = "vars"
      ),
      bg = function(x) {
        dplyr::case_when(
          abs(x) > limit ~ bg_color,
          # x < -0.1 ~ "gray",
          TRUE ~ "white"
        )
      }
    )
  return(ft)
}


add_labels_cols <- function(df) {
  df <- df |>
    mutate(
      to = case_when(
        lhs == "PANSS1.total" ~ "PANSS Total (T1)",
        lhs == "PANSS2.total" ~ "PANSS Total (T2)",
        lhs == "PANSS3.total" ~ "PANSS Total (T3)",
        lhs == "gender" ~ "Gender",
        lhs == "age" ~ "Age",
        lhs == "mpvs_total_score" ~ "MPVS",
        lhs == "employment_ever" ~ "Ever employed"
      )
    ) |>
    mutate(
      from = case_when(
        rhs == "PANSS1.total" ~ "PANSS Total (T1)",
        rhs == "PANSS2.total" ~ "PANSS Total (T2)",
        rhs == "PANSS3.total" ~ "PANSS Total (T3)",
        rhs == "gender" ~ "Gender",
        rhs == "age" ~ "Age",
        rhs == "mpvs_total_score" ~ "MPVS",
        rhs == "employment_ever" ~ "Ever employed"
      )
    )
}
##################
replace_target <- function(data, columns, mapping) {
  data[columns] <- lapply(data[columns], function(x) {
    idx <- x %in% names(mapping)
    x[idx] <- unname(mapping[x[idx]])
    return(x)
  })

  return(data)
}
library(testthat)

test_that("mapped values are replaced in a single column", {
  df <- data.frame(
    outcome = c("mpvs_total_score", "height"),
    stringsAsFactors = FALSE
  )

  mapping <- c(
    mpvs_total_score = "MPVS"
  )

  result <- replace_target(df, "outcome", mapping)

  expect_equal(result$outcome, c("MPVS", "height"))
})

test_that("mapped values are replaced in multiple columns", {
  df <- data.frame(
    outcome = c("mpvs_total_score", "height"),
    predictor = c("age_years", "weight"),
    stringsAsFactors = FALSE
  )

  mapping <- c(
    mpvs_total_score = "MPVS",
    age_years = "Age"
  )

  result <- replace_target(df, c("outcome", "predictor"), mapping)

  expect_equal(result$outcome, c("MPVS", "height"))
  expect_equal(result$predictor, c("Age", "weight"))
})

test_that("unmapped values remain unchanged", {
  df <- data.frame(
    outcome = c("height", "weight"),
    predictor = c("sex", "site"),
    stringsAsFactors = FALSE
  )

  mapping <- c(
    mpvs_total_score = "MPVS"
  )

  result <- replace_target(df, c("outcome", "predictor"), mapping)

  expect_identical(result, df)
})

test_that("columns not listed are left unchanged", {
  df <- data.frame(
    outcome = c("mpvs_total_score", "height"),
    predictor = c("age_years", "weight"),
    stringsAsFactors = FALSE
  )

  mapping <- c(
    mpvs_total_score = "MPVS",
    age_years = "Age"
  )

  result <- replace_target(df, "outcome", mapping)

  expect_equal(result$outcome, c("MPVS", "height"))
  expect_equal(result$predictor, c("age_years", "weight"))
})

test_that("original data frame is not modified", {
  df <- data.frame(
    outcome = c("mpvs_total_score"),
    stringsAsFactors = FALSE
  )

  mapping <- c(
    mpvs_total_score = "MPVS"
  )

  original <- df

  result <- replace_target(df, "outcome", mapping)

  expect_identical(df, original)
  expect_equal(result$outcome, "MPVS")
})

test_that("works with an empty vector of columns", {
  df <- data.frame(
    outcome = c("mpvs_total_score"),
    stringsAsFactors = FALSE
  )

  mapping <- c(
    mpvs_total_score = "MPVS"
  )

  result <- replace_target(df, character(0), mapping)

  expect_identical(result, df)
})
###################
