library(lavaan)
library(Lmisc)
source("helper.R")
source("constants.R")
df_selected <- df_named |>
  select(
    all_of(
      c(
        "gender_fct",
        "age",
        "employment_ever",
        "employment_ever_fct",
        "PANSS1.total",
        "PANSS2.total",
        "PANSS3.total",
        "mpvs_total_score"
      )
    )
  )


#####################
# Linear regression #
#####################

model_lm1 <- lm(
  mpvs_total_score ~ PANSS1.total + age + gender_fct + employment_ever_fct,
  data = df_selected
)

summary(model_lm1)

model_lm2 <- lm(
  PANSS2.total ~
    mpvs_total_score + PANSS1.total + age + gender_fct + employment_ever_fct,
  data = df_selected
)

summary(model_lm2)

#######
# SEM #
#######
model <- "
PANSS1.total ~ gender_fct
PANSS1.total ~ age
PANSS1.total ~ employment_ever_fct

PANSS2.total ~ gender_fct
PANSS2.total ~ age
PANSS2.total ~ employment_ever_fct
PANSS2.total ~ PANSS1.total
PANSS2.total ~ mpvs_total_score


PANSS3.total ~ gender_fct
PANSS3.total ~ age
PANSS3.total ~ employment_ever_fct
PANSS3.total ~ PANSS1.total
PANSS3.total ~ PANSS2.total
PANSS3.total ~ mpvs_total_score

mpvs_total_score ~ gender_fct
mpvs_total_score ~ age
mpvs_total_score ~ PANSS1.total
"


model_ml <- sem(
  model = model,
  data = df_selected
)
summary(model_ml)


model_ml_residuals <- round(
  change_df_labels(
    df = extract_cov_residuals(
      residual_obj = resid(
        model_ml,
        type = "cor.bollen"
      )
    ),
    labels = MAPPING
  ),
  digits = 3
)


color_corr_residuals(
  resid_df = model_ml_residuals,
  limit = 0.1,
  bg_color = "gray"
)

# To adjust the model
# modindices(model_ml)

parameters_model_ml_unstd <- Lmisc::modify_parameter_estimates(
  df = parameterestimates(
    model_ml,
    standardized = F
  ),
  round_digits = 2
)


########
# Fiml #
########
model_fiml <- sem(
  model = model,
  data = df_selected,
  missing = "fiml"
)
summary(model_fiml)

modificationIndices(model_fiml, sort. = T)

model_fiml_residuals <- round(
  change_df_labels(
    df = extract_cov_residuals(
      residual_obj = resid(
        model_fiml,
        type = "cor.bollen"
      )
    ),
    labels = MAPPING
  ),
  digits = 3
)


color_corr_residuals(
  resid_df = model_fiml_residuals,
  limit = 0.1,
  bg_color = "gray"
)


parameters_model_fiml_unstd <- Lmisc::modify_parameter_estimates(
  df = parameterestimates(
    model_fiml,
    standardized = F
  ),
  round_digits = 2
)


# #
# model_modified <- "
# PANSS1.total ~ gender_fct
# PANSS1.total ~ age
# PANSS1.total ~ employment_ever_fct
#
#
# PANSS2.total ~ gender_fct
# PANSS2.total ~ age
# PANSS2.total ~ employment_ever_fct
# PANSS2.total ~ PANSS1.total
#
# PANSS3.total ~ gender_fct
# PANSS3.total ~ age
# PANSS3.total ~ employment_ever_fct
# PANSS3.total ~ PANSS2.total
# PANSS3.total ~ PANSS1.total
#
#
# mpvs_total_score ~ gender_fct
# mpvs_total_score ~ age
# mpvs_total_score ~ PANSS1.total
# mpvs_total_score ~ PANSS2.total
#
# mpvs_total_score ~~ 0*PANSS3.total
#
# "
#
#
# model_fiml_modified <- sem(
#   model = model_modified,
#   data = df_selected,
#   missing = "fiml"
# )
# summary(model_fiml_modified)
#
# modificationIndices(model_fiml_modified, sort. = T)
#
# color_corr_residuals(
#   resid_df = round(
#     change_df_labels(
#       df = extract_cov_residuals(
#         residual_obj = resid(
#           model_fiml_modified,
#           type = "cor.bollen"
#         )
#       ),
#       labels = MAPPING
#     ),
#     digits = 3
#   ),
#   limit = 0.1,
#   bg_color = "gray"
# )
#
#
# parameters_model_fiml_modified <- modify_parameter_estimates(
#   df = parameterestimates(
#     model_fiml_modified,
#     standardized = F
#   ),
#   round_digits = 2
# )

####################
# Use imputed data #
####################

# Do not use semTools, it's deprecated
library(lavaan.mi)
print(sys.nframe())
# if (sys.nframe() %in% c(0, 4)) {
if (exists("imp_data") == F) {
  if (file.exists("imp.Rdata")) {
    load("imp.Rdata")
  } else {
    source("imputations.R")
  }
}
# }

model <- "
PANSS1.total ~ gender
PANSS1.total ~ age
PANSS1.total ~ employment_ever

PANSS2.total ~ gender
PANSS2.total ~ age
PANSS2.total ~ employment_ever
PANSS2.total ~ PANSS1.total
PANSS2.total ~ mpvs_total_score

PANSS3.total ~ gender
PANSS3.total ~ age
PANSS3.total ~ employment_ever
PANSS3.total ~ PANSS1.total
PANSS3.total ~ PANSS2.total
PANSS3.total ~ mpvs_total_score

mpvs_total_score ~ gender
mpvs_total_score ~ age
mpvs_total_score ~ PANSS1.total
"


fit_mi <- sem.mi(
  model = model,
  data = imp_data
)
# modindices.mi(fit_mi)
summary(fit_mi)
parameterEstimates.mi(fit_mi)


parameters_model_mi <- modify_parameter_estimates(
  df = parameterEstimates.mi(
    fit_mi,
    standardized = F
  ),
  round_digits = 2
)


# Raw
# resid_fit_mi <- as.data.frame(resid(fit_mi))

# Corralated residuals
model_mi_residuals <- round(
  change_df_labels(
    df = extract_cov_residuals(
      residual_obj = resid(
        fit_mi,
        type = "cor.bollen"
      )
    ),
    labels = MAPPING
  ),
  digits = 3
)


color_corr_residuals(
  resid_df = model_mi_residuals,
  limit = 0.1,
  bg_color = "gray"
)
