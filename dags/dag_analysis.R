library(ggplot2)
library(ggdag)
library(tidyverse)
library(Lmisc)

source("helper.R")

dag <- ggdag::dagify(
  mpvs_total_score ~ PANSS1.total,
  mpvs_total_score ~ PANSS2.total,
  mpvs_total_score ~ age,
  mpvs_total_score ~ gender,
  PANSS1.total ~ age,
  PANSS1.total ~ gender,
  PANSS1.total ~ employment_ever,
  PANSS2.total ~ age,
  PANSS2.total ~ gender,
  PANSS2.total ~ employment_ever,
  PANSS2.total ~ PANSS1.total,
  outcome = "mpvs_total_score",
  coords = list(
    x = c(
      PANSS1.total = 0,
      PANSS2.total = 6,
      mpvs_total_score = 7,
      gender = 0,
      age = 3,
      employment_ever = 3
    ),
    y = c(
      PANSS1.total = 4,
      PANSS2.total = 4,
      mpvs_total_score = 0,
      gender = 0,
      age = 7,
      employment_ever = 0.5
    )
  ),
  labels = c(
    PANSS1.total = "PANSS Total (T1)",
    PANSS2.total = "PANSS Total (T2)",
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
    node_id = name,
    x = x,
    y = y,
    label = label
  ) %>%
  mutate(node_id = label)


edges <- tidy_modified$data

# Replace names with their corresponding labels
edges <- replace_target(
  data = edges,
  column = c("name", "to"),
  mapping = c(
    "PANSS1.total" = "PANSS Total (T1)",
    "PANSS2.total" = "PANSS Total (T2)",
    "gender" = "Gender",
    "age" = "Age",
    "mpvs_total_score" = "MPVS",
    "employment_ever" = "Ever employed"
  )
)

temp <- data.frame(
  # Extract only the regression paths
  parameters_model_ml_unstd[
    parameters_model_ml_unstd$op == "~",
  ]
) %>%
  replace_target(
    columns = c("lhs", "rhs"),
    mapping = c(
      "PANSS1.total" = "PANSS Total (T1)",
      "PANSS2.total" = "PANSS Total (T2)",
      "gender" = "Gender",
      "age" = "Age",
      "mpvs_total_score" = "MPVS",
      "employment_ever" = "Ever employed"
    )
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
    hjust = case_when(
      curvature == 1 & pvalue <= 0.05 ~ 0.6,
      curvature == 1 & pvalue > 0.05 ~ 0.5,
      pvalue <= 0.05 ~ 0.5,
      TRUE ~ 0.5,
    )
  ) %>%
  mutate(
    hjust = case_when(
      # from == "MPVS12" & to == "MPVS16" ~ 0.8,
      from == "Gender" & to == "PANSS Total (T2)" ~ 0.8,
      from == "Ever employed" & to == "PANSS Total (T1)" ~ 0.6,
      from == "PANSS Total (T1)" & to == "MPVS" ~ 0.3,
      from == "Ever employed" & to == "PANSS Total (T2)" ~ 0.55,
      from == "Age" & to == "MPVS" ~ 0.3,
      .default = 0.5
    )
  ) %>%
  mutate(
    vjust = case_when(
      # from == "MPVS12" & to == "MPVS16" ~ 0.8,
      # from == "MPVS14" & to == "MPVS21" ~ 0.8,
      # from == "Gender" & to == "PANSS Total (T2)" ~ 1,
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
  xlim = c(-0.5, 7.2),
  ylim = c(-0.5, 7)
  # footnote = foonote_acronymns_for_dag_plots,
  # footnote_size = 14
)
p

save_dag("model_ml_unstd.tiff", p)


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
  xlim = c(-0.5, 7.2),
  ylim = c(-0.5, 7)
  # footnote = foonote_acronymns_for_dag_plots,
  # footnote_size = 14
)
p

save_dag("model_fiml_unstd.tiff", p)

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
  xlim = c(-0.5, 7.2),
  ylim = c(-0.5, 7)
  # footnote = foonote_acronymns_for_dag_plots,
  # footnote_size = 14
)
p

save_dag("model_unstd_mi.tiff", p)
