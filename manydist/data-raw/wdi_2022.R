# Rebuild the fixed 2022 World Development Indicators snapshot distributed
# with manydist. The packaged snapshot was retrieved on 2026-05-05. Because
# historical WDI values and country metadata can be revised, inspect any
# differences before replacing data/wdi_2022.rda.

library(dplyr)
library(tidyr)
library(WDI)

wdi_raw <- WDI(
  country = "all",
  indicator = c(
    gdp_pc = "NY.GDP.PCAP.CD",
    life_exp = "SP.DYN.LE00.IN",
    unemployment = "SL.UEM.TOTL.ZS",
    urban_pop = "SP.URB.TOTL.IN.ZS",
    pop_growth = "SP.POP.GROW"
  ),
  start = 2022,
  end = 2022,
  extra = TRUE
)

wdi_2022 <- wdi_raw |>
  filter(.data$region != "Aggregates") |>
  transmute(
    country = .data$country,
    region = factor(.data$region),
    income = factor(.data$income),
    lending = factor(.data$lending),
    gdp_pc_kusd = round(.data$gdp_pc / 1000, digits = 1),
    life_exp = round(.data$life_exp, digits = 1),
    unemployment = round(.data$unemployment, digits = 1),
    urban_pop = round(.data$urban_pop, digits = 1),
    pop_growth = round(.data$pop_growth, digits = 2)
  ) |>
  drop_na() |>
  arrange(.data$country)

# The row-count check prevents an upstream revision from silently replacing
# the documented snapshot with a materially different dataset.
stopifnot(
  nrow(wdi_2022) == 181L,
  ncol(wdi_2022) == 9L
)

save(
  wdi_2022,
  file = "data/wdi_2022.rda",
  version = 3,
  compress = "xz"
)
