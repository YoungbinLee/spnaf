#' Calculate spatial autocorrelation with OD data and corresponding polygons.
#'
#' @param df A data.frame that contains your Origin-Destination data. The df must consist of "oid" (origin id), "did" (destination id), "n" (flow weight).
#' @param shape A shapefile (in a polygon type) that matches to your OD dataframe. The shape must have an "id" column to match your ids in df.
#' @param queen A T/F input that is used to calculate \code{spdep}'s spatial contingency (Please view documents of \link[spdep]{poly2nb} for more information).
#' @param snap A parameter that is also used to calculate \code{spdep}'s spatial contingency (Please view documents of \link[spdep]{poly2nb} for more information).
#' @param method A string value among "o" (origin based), "d" (destination based), and "t" (both way) which determines the way to generate Spatial Weights. The default value is "t".
#' @param n An integer value to define how many times you want to execute bootstrapping.
#' @return The result is in the form of a list which includes a dataframe and a \code{sf} object.
#' Both contain Gij statistics and p-value columns merged to your input df. The geometry type of the latter is linestring.
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
#' # Execution of Gij.polygon with data above and given parameters
#' result <- Gij.polygon(df = OD, shape = CA_polygon, queen = TRUE, snap = 1,
#' method = 't', n = 1000)
#'
#' # check the results
#' head(result[[1]])
#' head(result[[2]])
#'
#' @references
#' Berglund, S., & Karlström, A. (1999). Identifying local spatial association in
#' flow data, Journal of Geographical Systems, 1(3), 219-236. https://doi.org/10.1007/s101090050013

#' @importFrom magrittr %>%
#' @export Gij.polygon

# Functions ---------------------------------------------------------------


Gij.polygon <- function(df, shape,
                  queen = TRUE, snap = 1,
                  method = 't', n = 1000){
    oid <- did <- Gij <- NULL

    sw <- SpatialWeight(df, shape, snap, queen)
    # result_frame: OD data + G statistic
    result_frame <- Gstat(SpecialWeights = sw, method = method) %>%
        dplyr::select(oid, did, n, Gij)
    # result_frame: OD data + G statistic + pval
    result_frame <- Boot(rf = result_frame, n = n)
    # result_lines: OD data + G statistic + pval + WKT(lines)
    result_lines <- Resultlines(shape, result_frame)

    return(list(result_frame, result_lines))

}
