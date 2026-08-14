library(gtsummary)
library(tidyverse)


needed <- df_imp |> select(!all_of("PANSS3.total"))
labels_list <- list(
  "PANSS1.total" ~ "PANSS Total (T1)",
  "PANSS2.total" ~ "PANSS Total (T2)",
  #"gender" ~ "Gender",
  "age" ~ "Age",
  "mpvs_total_score" ~ "MPVS",
  "employment_ever" ~ "Ever employed"
)

labels_list_deparsed <- lapply(
  labels_list,
  function(f) {
    lhs_str <- deparse(f[[2]])
    rhs_str <- deparse(f[[3]])

    # Change to known env
    d <- as.formula(paste(lhs_str, '~', rhs_str))
    environment(d) <- .GlobalEnv
    return(d)
  }
)
# test
# is_list(test)
# names(test)
# str(test)
################
# Descriptives #
################

summary_all2 <- needed |>
  gtsummary::tbl_summary(
    by = gender,
    statistic = list(
      all_categorical() ~ "{n}    ({p}%)",
      all_continuous() ~ "{mean} ({sd})"
    ),
    type = list(
      # gender ~ "categorical",
      employment_ever ~ "categorical"
    ),
    missing_text = "Missing",
    missing_stat = "{N_miss} ({p_miss}%)",
    label = labels_list_deparsed
    #   lapply(
    #   list_to_individual_formulas(VAR_LABELS),
    #   function(f) deparse(f)
    # )
  ) |>
  bold_labels() |>
  modify_header(
    # See https://www.danieldsjoberg.com/gtsummary/reference/modify.html#ref-examples
    all_stat_cols() ~ "**{level}** (n = {n})"
  ) |>
  # Include an "overall" column
  add_overall(
    last = T,
    # The ** make it bold
    col_label = "**All participants**<br> \n (N = {N})"
  ) |>
  # modify_caption("Table 1. Participant characteristics") |>
  add_p(
    test = list(all_continuous() ~ "t.test", all_categorical() ~ "chisq.test"),
    pvalue_fun = label_style_pvalue(digits = 3)
  ) |>
  # See https://stackoverflow.com/q/74705598
  modify_footnote(everything() ~ NA) |>
  # This functions converts gtsummary to gt
  bstfun::bold_italicize_group_labels(bold = T) |>
  gt::tab_header(
    title = "Table 1. Participants characteristics",
    subtitle = ""
  ) |>
  gt::tab_options(heading.subtitle.font.size = "20px") |>
  gt::tab_footnote(
    footnote = "",
    locations = NULL
  )

summary_all2


######################
# Correlation Matrix #
######################
needed_num <- needed |> mutate(across(everything(), as.numeric))
corr_mat <- cor(
  needed_num,
  use = "pairwise.complete.obs",
  method = "spearman"
)
vars <- c(
  "PANSS1.total",
  "PANSS2.total",
  "mpvs_total_score",
  "age",
  "gender",
  "employment_ever"
)

# P-values
p_mat <- matrix(
  NA,
  nrow = ncol(needed_num),
  ncol = ncol(needed_num),
  dimnames = list(vars, vars)
)


for (i in seq_along(vars)) {
  for (j in seq_along(vars)) {
    p_mat[i, j] <- cor.test(
      needed_num[[i]],
      needed_num[[j]],
      method = "spearman",
      exact = FALSE
    )$p.value
  }
}

# extract raw p-values
pvals <- p_mat[lower.tri(p_mat)]

pvals_adj <- p.adjust(
  pvals,
  method = "BH"
)

# Recreate the matrix with the adjusted p-values
p_mat_adj <- p_mat
p_mat_adj[lower.tri(p_mat)] <- pvals_adj
p_mat_adj[upper.tri(p_mat_adj)] <- NA
diag(p_mat_adj) <- NA

# The correlation matrix to be displayed
# Converted to characters
corr_display <- matrix(
  sprintf("%.2f", corr_mat),
  nrow = nrow(corr_mat),
  dimnames = dimnames(corr_mat)
)
# Add significance
# corr_display[pvals_adj < .05] <- paste0(corr_display[pvals_adj < .05], "*")
# corr_display[pvals_adj < .01] <- paste0(sprintf("%.2f", corr_mat[pvals_adj < .01]), "**")
# corr_display[pvals_adj < .001] <- paste0(sprintf("%.2f", corr_mat[pvals_adj < .001]), "***")

corr_display[upper.tri(corr_display)] <- ""
diag(corr_display) <- ""
list(
  "PANSS1.total" = "PANSS Total (T1)",
  "PANSS2.total" = "PANSS Total (T2)",
  "gender" = "Gender",
  "age" = "Age",
  "mpvs_total_score" = "MPVS",
  "employment_ever" = "Ever employed"
)

table_corr_mat_mpvs_dcq <- corr_display |>
  as.data.frame() |>
  gt::gt(
    # https://stackoverflow.com/questions/75260770/how-to-add-the-row-names-to-my-gt-table
    rownames_to_stub = T,
    caption = "" # Table 2. Correlation matrix of DCQ and MPVS across waves
  ) |>
  gt::cols_label(
    .list = list_to_individual_formulas(
      VAR_LABELS
    )
  ) |>
  gt::text_case_match(
    "PANSS1.total" ~ "PANSS Total (T1)",
    "PANSS2.total" ~ "PANSS Total (T2)",
    "gender" ~ "Gender",
    "age" ~ "Age",
    "mpvs_total_score" ~ "MPVS",
    "employment_ever" ~ "Ever employed",
    .locations = gt::cells_stub()
  ) |>
  # gt::fmt_number(
  #   decimals = 2
  # ) |>
  gt::tab_options(heading.subtitle.font.size = "20pt") |>
  gt::tab_header(
    title = "Table 2. Spearman correlation matrix",
    subtitle = ""
  ) |>
  gt::tab_footnote(
    footnote = paste(
      "PANSS: Positive and Negative Syndrome Scale;",
      "MPVS: Multidimensional peer victimization scale;",
      "y: year(s).",
      "Bold coefficients indicate p-value <0.05",
      # "* p < .05, ** p < .01, *** p < .001",
      "(Adjusted for multiple testing using Benjamini & Hochberg procedure)."
    ),
    locations = NULL
  )


# Bold cells for p-value <0.05
sig_idx <- which(
  p_mat_adj < .05 &
    lower.tri(p_mat_adj),
  arr.ind = TRUE
)

for (k in seq_len(nrow(sig_idx))) {
  row_name <- rownames(p_mat_adj)[sig_idx[k, "row"]]
  col_name <- colnames(p_mat_adj)[sig_idx[k, "col"]]

  table_corr_mat_mpvs_dcq <- table_corr_mat_mpvs_dcq |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(
        rows = row_name,
        columns = col_name
      )
    )
}
table_corr_mat_mpvs_dcq
