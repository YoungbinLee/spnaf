#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom sf st_geometry_type
#' @importFrom sf st_centroid
#' @importFrom sf st_as_sf
#' @importFrom sf st_coordinates
#' @importFrom sf st_polygon
#' @importFrom sf st_sfc
#' @importFrom sf st_sf
#' @importFrom spdep nbdists

SpatialWeight <- function(df, shape, snap = NULL, method = "queen", k = NULL,
                          d = NULL, idw=FALSE, row_standardized= FALSE) {
    cat("(1) Creating Spatial Weights... ")

    w <- NULL
    original_id <- as.character(shape$id)

    if(method == "queen") {
        if (!inherits(shape, "SpatialPolygons") && !any(st_geometry_type(shape) == "POLYGON")) {
            centroids <- st_centroid(st_as_sf(shape))
            shape <- create_thiessen_polygons(centroids)
        }
        nb <- unclass(spdep::poly2nb(shape, snap = snap, queen = TRUE)) # create nb data from input shape
    } else if (method == "KNN") {
        coo <- st_centroid(shape)
        nb <- unclass(spdep::knn2nb(spdep::knearneigh(coo, k = k))) # create nb data from input shape
    } else if (method == "fixed_distance") {
        if (!inherits(shape, "SpatialPolygons") && !any(st_geometry_type(shape) == "POLYGON")) {
            centroids <- st_centroid(st_as_sf(shape))
            shape <- create_thiessen_polygons(centroids)
        }
        nb <- spdep::poly2nb(shape, snap = snap, queen = TRUE)
        nb2 <- spdep::nblag(neighbours = nb, maxlag = 2)
        nb <- unclass(nb2[[2]])
    }

    nb <- lapply(nb, function(data) {
        return(original_id[data]) # substitute adjacent index for original id
    })
    names(nb) <- original_id # reference index for original id

    if (idw) {
        if (!inherits(shape, "SpatialPolygons") && !any(st_geometry_type(shape) == "POLYGON")) {
            centroids <- st_centroid(st_as_sf(shape))
            shape <- create_thiessen_polygons(centroids)
        }
        coordinates_matrix <- as.matrix(st_coordinates(st_centroid(st_as_sf(shape))))
        dist_matrix <- nbdists(nb, coordinates_matrix, longlat = TRUE)
        inv_dist_list <- lapply(dist_matrix, function(x) 1 / x)
        inv_dist_list <- lapply(inv_dist_list, function(x) ifelse(is.infinite(x), 0, x))
        if (row_standardized) {
            inv_dist_list <- lapply(inv_dist_list, function(x) x / sum(x))
        }
    } else {
        weight_list <- lapply(nb, function(x) rep(1, length(x)))
        if (row_standardized) {
            weight_list <- lapply(weight_list, function(x) x / sum(x))
        }
    }

    if (idw) {
        nb_frame <- data.frame(
            oid = rep(names(nb), sapply(nb, length)),
            did = unlist(nb),
            w = unlist(inv_dist_list)
        )
    } else {
        nb_frame <- data.frame(
            oid = rep(names(nb), sapply(nb, length)),
            did = unlist(nb),
            w = unlist(weight_list)
        )
    }

    # convert list to data.frame
    did <- utils::stack(nb)
    did <- stats::setNames(as.character(did$values), as.character(did$ind))
    oid <- names(did)
    nb_frame <- data.frame(oid = as.character(oid))
    nb_frame$did <- as.character(did)
    nb_frame$w <- 1 # add w as 1 for contiguous polygons

    # join spatial weights to input dataframe
    df <- df %>%
        dplyr::mutate(oid = as.character(oid), did = as.character(did))
    result <- df %>%
        dplyr::left_join(nb_frame, by = c("oid", "did")) %>%
        dplyr::filter(!is.na(.data$n)) # filter observations with valid n only

    # Union is created by input polygons^2
    U <- expand.grid(oid = original_id, did = original_id)
    U <- U %>% dplyr::filter(oid != did) # remove flow from i to i

    # Join OD data with spatial weights to Union
    result <- U %>%
        dplyr::left_join(result, by = c("oid", "did")) %>%
        dplyr::mutate(n = ifelse(is.na(.data$n), 0, .data$n),
                      w = ifelse(is.na(.data$w), 0, 1))

    cat("Done!\n")
    # Alarm the size of the union set
    cat(paste0("note: Total ", nrow(U), " network combinations are ready to be analyzed\n"))
    cat(paste0("note: which is calculated by the input polygon data\n"))

    return(result)
}
