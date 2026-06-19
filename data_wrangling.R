library(tidyverse)


df_raw <- readxl::read_excel(path = "Athens FEP_lampros.xlsx")

name_map <- c(
  AGE = "age",
  GENDER = "gender",
  DATE_BIRTH = "date_of_birth",
  AGE_OF_ONSET = "age_of_onset",
  EMPLOYMENT = "employment_ever",
  M_INCOME_NOW = "income_montly_now",
  M_INCOME_1YR = "income_monthly_1y_ago",
  M_INCOME_5YR = "income_monthly_5y_ago",
  TOTAL_SCORE = "mpvs_total_score",
  GAF1st = "gaf_baseline",
  GAF2nd = "gaf_1month",
  LANGUAGE = "language",
  MIGRATION = "migration",
  WORK_NOW = "work_now",
  PHYSICAL_ABUSE = "ctq_physical_abuse",
  SEXUAL_ABUSE = "ctq_sexual_abuse",
  EMOTIONAL_ABUSE = "ctq_emotional_abuse",
  PHYSICAL_NEGLECT = "ctq_physical_neglect",
  EMOTIONAL_NEGLECT = "ctq_emotional_neglect"
)

df_named <- df_raw %>%
  rename(!!!setNames(names(name_map), name_map))

# PANSS1p1 -> PANSS1p_item1
names(df_named) <- gsub(
  "^(PANSS\\d+[pngt])(\\d+)$",
  "\\1_item\\2",
  names(df_named)
)

# peervictimization1 <- peervictimization_item1
names(df_named) <- gsub(
  "^(peervictimization)(\\d+)$",
  "\\1_item\\2",
  names(df_named)
)

# Rename peervictimization to mpvs
names(df_named) <- gsub(
  "^peervictimization",
  "mpvs",
  names(df_named)
)

df_named[df_named == 777] <- NA # This value appears in AGE_OF_MIGRATION
df_named[df_named == 999] <- NA # This indicates NA
df_named[df_named == 9999] <- NA # This probably is also NA

# View(explore::describe(df_named))

all.equal(
  df_named$PANSS1.total,
  df_named$PANSS1p.total + df_named$PANSS1n.total + df_named$PANSS1g.total
)
all.equal(
  df_named$PANSS2.total,
  df_named$PANSS2p.total + df_named$PANSS2n.total + df_named$PANSS2g.total
)
