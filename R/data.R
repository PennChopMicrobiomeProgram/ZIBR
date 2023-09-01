#' Longitudinal human microbiome data
#'
#' A dataset containing the bacterial abundance and clinical information from
#' a longitudinal human microbiome study
#'
#' @docType data
#' @keywords datasets
#'
#' @format A data frame with 236 rows and 5 variables:
#' \describe{
#'   \item{Sample}{Sample IDs}
#'   \item{Subject}{Subject IDs}
#'   \item{Time}{Time points}
#'   \item{Treatment}{Treatment, 0 for antiTNF, 1 for EEN}
#'   \item{Abundance}{Abundance for Eubacterium}
#'   ...
#' }
#' @references Lewis and Chen et al. (2016) Cell Host & Microbe 18 (4), 489-500
"ibd"
