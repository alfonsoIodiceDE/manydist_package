install.packages("WDI")

library(WDI)
library(tidyverse)

wdi_raw <- WDI(
  country = "all",
  indicator = c(
    gdp_pc       = "NY.GDP.PCAP.CD",
    life_exp     = "SP.DYN.LE00.IN",
    unemployment = "SL.UEM.TOTL.ZS",
    urban_pop    = "SP.URB.TOTL.IN.ZS",
    pop_growth   = "SP.POP.GROW"
  ),
  start = 2022,
  end = 2022,
  extra = TRUE
)

wdi_data <- wdi_raw |>
  filter(region != "Aggregates") |>
  transmute(
    country,
    region  = factor(region),
    income  = factor(income),
    lending = factor(lending),
    gdp_pc_kusd  = round(gdp_pc / 1000, 1),
    life_exp     = round(life_exp, 1),
    unemployment = round(unemployment, 1),
    urban_pop    = round(urban_pop, 1),
    pop_growth   = round(pop_growth, 2)
  ) |>
  drop_na()


wdi_labels <- c(
  country      = "Country",
  region       = "Region",
  income       = "Income group",
  lending      = "World Bank lending category",
  gdp_pc_kusd  = "GDP per capita (k USD)",
  life_exp     = "Life expectancy (years)",
  unemployment = "Unemployment (%)",
  urban_pop    = "Urban population (% total)",
  pop_growth   = "Population growth (%)"
)

save(wdi_data, wdi_labels, file = "../zakopane_talk/data/wdi_socio_2022.RData")
