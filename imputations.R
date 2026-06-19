library(mice)
library(miceadds)
library(micemd)
library(beepr)

SEED <- 123
set.seed(seed = SEED)

df_imp <- df_named %>%
  dplyr::select(
    all_of(
      c(
        "gender",
        "migration",
        "work_now",
        "age",
        "employment_ever",
        "income_montly_now",
        # "income_monthly_1y_ago",
        # "income_monthly_5y_ago",
        grep("^PANSS1p_item\\d+$", names(df_named), value = TRUE),
        grep("^PANSS1n_item\\d+$", names(df_named), value = TRUE),
        grep("^PANSS1g_item\\d+$", names(df_named), value = TRUE),
        "PANSS1g.total",
        "PANSS1p.total",
        # 1 NA
        "PANSS1n.total",
        "PANSS1.total",
        # 2 NA
        grep("^PANSS2p_item\\d+$", names(df_named), value = TRUE),
        grep("^PANSS2n_item\\d+$", names(df_named), value = TRUE),
        grep("^PANSS2g_item\\d+$", names(df_named), value = TRUE),
        "PANSS2g.total",
        "PANSS2p.total",
        "PANSS2n.total",
        "PANSS2.total",
        # 4
        # CTQ1-25
        grep("^ctq_", names(df_named), value = TRUE),
        # >5
        colnames(df_named)[grepl(
          pattern = "mpvs_item", x = colnames(df_named)
        )],
        "mpvs_total_score"
      )
    )
  )

#####################
# Predictors Matrix #
#####################
predMatrix <- make.predictorMatrix(df_imp)


predMatrix <- modify_pred_matrix_scales(
  pred_matrix = predMatrix,
  item_pattern = "mpvs_item",
  total_pattern = "mpvs_total"
) %>%
  modify_pred_matrix_scales(
    item_pattern = "^PANSS1[png]_item",
    total_pattern = "^PANSS1[png]?.*total"
  ) %>%
  modify_pred_matrix_scales(
    item_pattern = "^PANSS2[png]_item",
    total_pattern = "^PANSS2[png]?.*total"
  )

# If we want to include the PANSS subscales, then these subscales should not
# predict anything because collinearity issues arised in mice!
predMatrix[, grep("^PANSS[12][png].total", names(df_named), value = TRUE)] <- 0



predMatrix <- suppressWarnings( #
  exclude_collinear_vars(
    pred_matrix = predMatrix,
    corr_mat = as.data.frame(
      cor(
        df_imp,
        use = "pairwise.complete.obs"
      )
    ),
    lower_threshold = 0.1,
    upper_threshold = 0.99
  )
)

######################
# Imputation Methods #
######################
impMethod <- make.method(data = df_imp)

for (var in colnames(df_imp)) {
  impMethod[var] <- "pmm"
}

##############################################
# The sum score needs to be properly handled #
##############################################

# `mice: Multivariate Imputation by Chained Equations in R`
# https://www.jstatsoft.org/article/view/v045i03
# See `Sum scores` section

# https://stefvanbuuren.name/fimd/sec-knowledge.html
# See `Sum scores` section
# Sum scores will be added as mentioned in url above
impMethod["mpvs_total_score"] <- "~I(mpvs_item1 + mpvs_item2 +
mpvs_item3 + mpvs_item4 +
mpvs_item5 + mpvs_item6 +
mpvs_item7 + mpvs_item8 +
mpvs_item9 + mpvs_item10 +
mpvs_item11 + mpvs_item12 +
mpvs_item13 + mpvs_item14 +
mpvs_item15 + mpvs_item16)"

impMethod["PANSS1p.total"] <- "~I(PANSS1p_item1 + PANSS1p_item2 +
PANSS1p_item3 + PANSS1p_item4 +
PANSS1p_item5 + PANSS1p_item6 + PANSS1p_item7)"

impMethod["PANSS1n.total"] <- "~I(PANSS1n_item1 + PANSS1n_item2 +
PANSS1n_item3 + PANSS1n_item4 +
PANSS1n_item5 + PANSS1n_item6 + PANSS1n_item7)"

impMethod["PANSS1g.total"] <- "~I(PANSS1g_item1 + PANSS1g_item2 +
PANSS1g_item3 + PANSS1g_item4 +
PANSS1g_item5 + PANSS1g_item6 + PANSS1g_item7 + PANSS1g_item8 +
PANSS1g_item9 + PANSS1g_item10 + PANSS1g_item11 + PANSS1g_item12 +
PANSS1g_item13 + PANSS1g_item14 + PANSS1g_item15 + PANSS1g_item16)"

impMethod["PANSS1.total"] <- "~I(PANSS1p_item1 + PANSS1p_item2 +
 PANSS1p_item3 + PANSS1p_item4 +
 PANSS1p_item5 + PANSS1p_item6 + PANSS1p_item7+PANSS1n_item1 + PANSS1n_item2 +
 PANSS1n_item3 + PANSS1n_item4 +
 PANSS1n_item5 + PANSS1n_item6 + PANSS1n_item7+PANSS1g_item1 + PANSS1g_item2 +
 PANSS1g_item3 + PANSS1g_item4 +
 PANSS1g_item5 + PANSS1g_item6 + PANSS1g_item7 + PANSS1g_item8 +
 PANSS1g_item9 + PANSS1g_item10 + PANSS1g_item11 + PANSS1g_item12 +
 PANSS1g_item13 + PANSS1g_item14 + PANSS1g_item15 + PANSS1g_item16)"


impMethod["PANSS2p.total"] <- "~I(PANSS2p_item1 + PANSS2p_item2 + PANSS2p_item3 + PANSS2p_item4 +
PANSS2p_item5 + PANSS2p_item6 + PANSS2p_item7)"

impMethod["PANSS2n.total"] <- "~I(PANSS2n_item1 + PANSS2n_item2 + PANSS2n_item3 + PANSS2n_item4 +
PANSS2n_item5 + PANSS2n_item6 + PANSS2n_item7)"

impMethod["PANSS2g.total"] <- "~I(PANSS2g_item1 + PANSS2g_item2 + PANSS2g_item3 + PANSS2g_item4 +
PANSS2g_item5 + PANSS2g_item6 + PANSS2g_item7 + PANSS2g_item8 +
PANSS2g_item9 + PANSS2g_item10 + PANSS2g_item11 + PANSS2g_item12 +
PANSS2g_item13 + PANSS2g_item14 + PANSS2g_item15 + PANSS2g_item16)"

impMethod["PANSS2.total"] <- "~I(PANSS2p_item1 + PANSS2p_item2 + PANSS2p_item3 + PANSS2p_item4 +
 PANSS2p_item5 + PANSS2p_item6 + PANSS2p_item7+PANSS2n_item1 + PANSS2n_item2 + PANSS2n_item3 + PANSS2n_item4 +
 PANSS2n_item5 + PANSS2n_item6 + PANSS2n_item7+PANSS2g_item1 + PANSS2g_item2 + PANSS2g_item3 + PANSS2g_item4 +
 PANSS2g_item5 + PANSS2g_item6 + PANSS2g_item7 + PANSS2g_item8 +
 PANSS2g_item9 + PANSS2g_item10 + PANSS2g_item11 + PANSS2g_item12 +
 PANSS2g_item13 + PANSS2g_item14 + PANSS2g_item15 + PANSS2g_item16)"


###########
# Impute! #
###########
start_time <- Sys.time()
imp <- mice(
  df_imp,
  method = impMethod,
  predictorMatrix = predMatrix,
  maxit = 5,
  m = 5,
  donors = 5
)
View(imp$loggedEvents)
print(imp$visitSequence)
print(Sys.time() - start_time)
