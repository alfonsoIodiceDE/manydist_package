#' Selected World Development Indicators for 2022
#'
#' A fixed snapshot of selected economic, demographic, and labour-market
#' indicators for 181 countries, together with World Bank regional, income,
#' and lending classifications. It provides a compact mixed-type dataset for
#' illustrating distance construction and distance-based learning.
#'
#' @format A data frame with 181 rows and 9 variables:
#' \describe{
#'   \item{country}{Country or economy name.}
#'   \item{region}{World Bank region, as a factor with 7 levels.}
#'   \item{income}{World Bank income group, as a factor with 5 levels.}
#'   \item{lending}{World Bank lending category, as a factor with 4 levels.}
#'   \item{gdp_pc_kusd}{GDP per capita in thousands of current US dollars.}
#'   \item{life_exp}{Life expectancy at birth, in years.}
#'   \item{unemployment}{Unemployment as a percentage of the total labour
#'     force, based on the modelled ILO estimate.}
#'   \item{urban_pop}{Urban population as a percentage of the total
#'     population.}
#'   \item{pop_growth}{Annual population growth, in percent.}
#' }
#'
#' @details
#' The reference year is 2022 and the snapshot was retrieved from the World
#' Bank Indicators API on 2026-05-05. Aggregate regions and observations with
#' missing values in any selected field were removed. GDP per capita was
#' divided by 1,000 and rounded to one decimal place; life expectancy,
#' unemployment, and urban population were rounded to one decimal place; and
#' population growth was rounded to two decimal places.
#'
#' The numerical indicators and their World Bank codes are:
#' \itemize{
#'   \item GDP per capita: \code{NY.GDP.PCAP.CD};
#'   \item life expectancy: \code{SP.DYN.LE00.IN};
#'   \item unemployment: \code{SL.UEM.TOTL.ZS};
#'   \item urban population: \code{SP.URB.TOTL.IN.ZS};
#'   \item population growth: \code{SP.POP.GROW}.
#' }
#'
#' @section License:
#' The original World Development Indicators data are distributed by the
#' World Bank under the Creative Commons Attribution 4.0 International
#' license, subject to the World Bank Dataset Terms. Attribution: The World
#' Bank: World Development Indicators.
#'
#' @source
#' The World Bank, \emph{World Development Indicators},
#' \url{https://datacatalog.worldbank.org/search/dataset/0037712/world-development-indicators}.
#'
#' @references
#' World Bank Dataset Terms,
#' \url{https://www.worldbank.org/ext/en/legal/terms-conditions/datasets}.
#'
#' @examples
#' data("wdi_2022", package = "manydist")
#' dplyr::glimpse(wdi_2022)
#'
#' @docType data
#' @keywords datasets
#' @name wdi_2022
NULL
