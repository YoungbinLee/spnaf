#' Calculate spatial autocorrelation with OD data and corresponding flows.
#'
#' @param df A data.frame that contains your Origin-Destination data. The df must consist of "oid" (origin id), "did" (destination id), and "n" (flow weight).
#' @param shape A shapefile (in a polygon type) that matches your OD dataframe. The shape must have an "id" column to match your ids in df.
#' @param method A string value among "queen" (spatial contiguity), "KNN" (k-nearest neighbors), and "fixed_distance" (fixed distance).
#' @param snap A parameter used to calculate \code{spdep}'s spatial contiguity (please refer to the \link[spdep]{poly2nb} documentation for more information).
#' @param k An integer value to define the number of nearest neighbors for the K-nearest neighbors method.
#' @param d An integer value to define the distance for the fixed distance spatial weight matrix.
#' @param idw A logical value indicating whether to use inverse distance weighting.
#' @param row_standardized A logical value indicating whether to row-standardize the spatial weights.
#' @param OD A string value among "o" (origin-based), "d" (destination-based), and "t" (both ways), which determines how to generate spatial weights. The default value is "t".
#' @param R An integer value to define how many times you want to execute bootstrapping.
#' @return The result is a list containing a dataframe and an \code{sf} object. Both contain Gij statistics and p-value columns merged with your input df. The geometry type of the latter is linestring.
#' @examples
#' # Data manipulation
#' CA <- spnaf::CA
#' OD <- cbind(CA$FIPS.County.Code.of.Geography.B, CA$FIPS.County.Code.of.Geography.A)
#' OD <- cbind(OD, CA$Flow.from.Geography.B.to.Geography.A)
#' OD <- data.frame(OD)
#' names(OD) <- c("oid", "did", "n")
#' OD$n <- as.numeric(OD$n)
#' OD <- OD[order(OD[,1], OD[,2]),]
#' head(OD) # check the input df's format
#'
#' # Load sf polygon
#' CA_polygon <- spnaf::CA_polygon
#' head(CA_polygon) # it has a geometry column
#'
#' # Execution of Gij.flow with data above and given parameters
#' \dontrun{
#' result <- Gij.flow(df = OD, shape = CA_polygon, method = 'queen', snap = 1, OD = 't', R = 1000)
#' }
#'
#' # check the results
#' \dontrun{
#' head(result[[1]])
#' head(result[[2]])
#' }
#' @references
#' Berglund, S., & Karlström, A. (1999). Identifying local spatial association in flow data. Journal of Geographical Systems, 1(3), 219-236. https://doi.org/10.1007/s101090050013
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export Gij.flow


Gij.flow <- function(df, shape, method = "queen", snap = 1, k = NULL, d = NULL,
                     idw = FALSE, row_standardized = FALSE, OD = 't', R = 1000){
  oid <- did <- Gij <- NULL

  sw <- SpatialWeight(df, shape, snap = snap, method = method, k = k, d = d)

  # result_frame: OD data + G statistic
  result_frame <- Gstat(SpatialWeights = sw, method = OD) %>%
    dplyr::select(oid, did, .data$n, Gij)

  # n > 0
  result_frame <- result_frame %>%
    dplyr::filter(.data$n > 0)

  # result_frame: OD data + G statistic + pval
  result_frame <- Boot(rf = result_frame, R = R)

  # result_lines: OD data + G statistic + pval + WKT(lines)
  result_lines <- Resultlines(shape, result_frame)

  return(list(result_frame, result_lines))
}
