library(ggplot2)
library(ggdag)
library(tidyverse)
library(Lmisc)

source("helper.R")
source("constants.R")
dag <- ggdag::dagify(
  mpvs_total_score ~ PANSS1.total,
  mpvs_total_score ~ age,
  mpvs_total_score ~ gender,

  PANSS1.total ~ age,
  PANSS1.total ~ gender,
  PANSS1.total ~ employment_ever,

  PANSS2.total ~ age,
  PANSS2.total ~ gender,
  PANSS2.total ~ employment_ever,
  PANSS2.total ~ PANSS1.total,
  PANSS2.total ~ mpvs_total_score,

  PANSS3.total ~ age,
  PANSS3.total ~ gender,
  PANSS3.total ~ employment_ever,
  PANSS3.total ~ PANSS1.total,
  PANSS3.total ~ PANSS2.total,
  PANSS3.total ~ mpvs_total_score,

  coords = list(
    x = c(
      PANSS1.total = 0,
      PANSS2.total = 6,
      PANSS3.total = 12,
      mpvs_total_score = 6,
      gender = 0,
      age = 3,
      employment_ever = 3
    ),
    y = c(
      PANSS1.total = 4,
      PANSS2.total = 4,
      PANSS3.total = 4,
      mpvs_total_score = 0,
      gender = 0,
      age = 7,
      employment_ever = 0.5
    )
  ),
  labels = c(
    PANSS1.total = "PANSS Total (T1)",
    PANSS2.total = "PANSS Total (T2)",
    PANSS3.total = "PANSS Total (T3)",
    mpvs_total_score = "MPVS",
    gender = "Gender",
    age = "Age",
    employment_ever = "Ever employed"
  )
)


tidy_modified <- tidy_dagitty(dag)

nodes <- tidy_modified$data %>%
  filter(!duplicated(name)) %>%
  dplyr::transmute(
    node_id = label,
    x = x,
    y = y,
    label = label
  )

edges <- tidy_modified$data

# Replace names with their corresponding labels
edges <- replace_target(
  data = edges,
  column = c("name", "to"),
  mapping = MAPPING
)

temp <- data.frame(
  # Extract only the regression paths
  parameters_model_ml_unstd[
    parameters_model_ml_unstd$op == "~",
  ]
) %>%
  replace_target(
    columns = c("lhs", "rhs"),
    mapping = MAPPING
  )

edges <- left_join(
  x = edges,
  y = temp,
  by = join_by(
    name == rhs,
    to == lhs
  )
)

edges$curvature <- 0

edges <- edges %>%
  filter(!is.na(to)) %>%
  transmute(
    from = name,
    to = to,
    curvature = curvature,
    pvalue = pvalue,
    est = est,
    ci.lower = ci.lower,
    ci.upper = ci.upper
  )
edges$curvature_amount <- 0
edges$gap <- -1


edges <- edges %>%
  mutate(
    gap = case_when(
      from == "PANSS Total (T1)" & to == "PANSS Total (T3)" ~ 3,
      .default = gap
    )
  ) |>
  mutate(
    curvature = case_when(
      from == "PANSS Total (T1)" & to == "PANSS Total (T3)" ~ 1,
      .default = curvature
    )
  ) |>
  mutate(
    curvature_amount = case_when(
      from == "PANSS Total (T1)" & to == "PANSS Total (T3)" ~ -0.1,
      .default = curvature_amount
    )
  ) |>
  mutate(
    hjust = case_when(
      curvature == 1 & pvalue <= 0.05 ~ 0.5,
      curvature == 1 & pvalue > 0.05 ~ 0.5,
      pvalue <= 0.05 ~ 0.5,
      TRUE ~ 0.5,
    )
  ) %>%
  mutate(
    hjust = case_when(
      # from == "MPVS12" & to == "MPVS16" ~ 0.8,
      from == "PANSS Total (T1)" & to == "PANSS Total (T3)" ~ 0.6,
      from == "Gender" & to == "PANSS Total (T2)" ~ 0.75,
      from == "Gender" & to == "PANSS Total (T3)" ~ 0.4,
      from == "Ever employed" & to == "PANSS Total (T1)" ~ 0.6,
      from == "PANSS Total (T1)" & to == "MPVS" ~ 0.3,
      from == "Ever employed" & to == "PANSS Total (T2)" ~ 0.5,
      from == "Age" & to == "MPVS" ~ 0.15,
      .default = 0.5
    )
  ) %>%
  mutate(
    vjust = case_when(
      # from == "MPVS12" & to == "MPVS16" ~ 0.8,
      # from == "MPVS14" & to == "MPVS21" ~ 0.8,
      # from == "Gender" & to == "PANSS Total (T2)" ~ 1,
      from == "MPVS" & to == "PANSS Total (T2)" ~ -1.5,
      from == "Gender" & to == "PANSS Total (T1)" ~ 0.8,
      from == "Ever employed" & to == "PANSS Total (T1)" ~ 0.2,
      .default = 0.5
    )
  ) %>%
  mutate(
    vertical_label_position = case_when(
      from == "Age" & to == "PANSS Total (T1)" ~ 0,
      .default = 0.5
    )
  )


p <- plot_dag(
  nodes = nodes,
  edges = edges,
  label_size = 10,
  label_size_unit = "pt",
  text_size = 4,
  xlim = c(-0.5, 12.2),
  ylim = c(-0.5, 7)
  # footnote = foonote_acronymns_for_dag_plots,
  # footnote_size = 14
)
p

save_dag("img\\model_ml_unstd.tiff", p)


########
# Fiml #
########

edges <- tidy_modified$data
edges <- wrap_edges(
  edges = edges,
  parameters = parameters_model_fiml_unstd
)

p <- plot_dag(
  nodes = nodes,
  edges = edges,
  label_size = 10,
  label_size_unit = "pt",
  text_size = 4,
  xlim = c(-0.5, 12.2),
  ylim = c(-0.5, 7)
  # footnote = foonote_acronymns_for_dag_plots,
  # footnote_size = 14
)
p

save_dag("img\\model_fiml_unstd.tiff", p)

#######################
# Multiple Imputation #
#######################

edges <- tidy_modified$data
edges <- wrap_edges(edges = edges, parameters = parameters_model_mi)

p <- plot_dag(
  nodes = nodes,
  edges = edges,
  label_size = 10,
  label_size_unit = "pt",
  text_size = 4,
  xlim = c(-0.5, 12.2),
  ylim = c(-0.5, 7)
  # footnote = foonote_acronymns_for_dag_plots,
  # footnote_size = 14
)
p

save_dag("img\\model_mi_unstd.tiff", p)
