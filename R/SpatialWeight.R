#' @importFrom magrittr %>%
#' @importFrom rlang .data

SpatialWeight <- function(df, shape, snap, queen){

    cat("(1) Creating Spatial Weights... ")

    w <- NULL
    nb <- unclass(spdep::poly2nb(shape, snap = snap, queen =queen)) # create nb data from input shape
    original_id <- shape$id
    nb <- lapply(nb, function(data){
        return(original_id[data]) # substitute adjacent index for original id
    })
    names(nb) <- original_id # reference index for original id

    # convert list to data.frame
    did <- utils::stack(nb)
    did <- stats::setNames(did$values, did$ind)
    oid <- names(did)
    nb_frame <- data.frame(oid = oid)
    nb_frame$did <- did
    nb_frame$w <- 1 # add w as 1 for contiguous polygons

    # join spatial weights to input dataframe
    result <- df %>%
        dplyr::left_join(nb_frame, by = c("oid", "did")) %>%
        dplyr::filter(!is.na(.data$n)) # filter observations with valid n only

    ## Union is created by input polygons^2
    U <- merge(shape$id, shape$id)
    names(U) <- c("oid", "did")
    U <- U %>% dplyr::filter(oid != did) # remove flow from i to i

    ### Join OD data with spatial weights to Union
    result <- U %>%
        dplyr::left_join(result, by = c("oid" = "oid", "did" = "did")) %>%
        dplyr::mutate(n = ifelse(is.na(.data$n), 0, .data$n),
                      w = ifelse(is.na(w), 0, 1))

    cat("Done !\n")
    ## Alarm the size of the union set
    cat(paste0("note: Total ", nrow(U), " network combinations are ready to be analyzed\n"))
    cat(paste0("note: which is calculated by the input polygon data\n"))

    return(result)
}
