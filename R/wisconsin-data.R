#' Breast Cancer Wisconsin Diagnostic
#'
#' Features of this dataset are computed from a digitized image of a
#' fine needle aspirate (FNA), describing characteristics of a breast
#' mass cell nuclei present in the image.
#'
#' @docType data
#'
#' @usage data(wisconsin)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
#' @references W.H. Wolberg, W.N. Street, O.L. Mangasarian, Breast Cancer Wisconsin
#' (Diagnostic) Data Set, UCI Mach. Learn. Repos. 1992.
#'
#' @source \href{https://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Diagnostic)}{UCI ML Repository}
#'
#' @examples
#' data(wisconsin)
#' area_mean <- wisconsin$area_mean
"wisconsin"
