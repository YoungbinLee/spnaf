#' Calculate spatial weights for networks based on input polygons.
#'
#' @param shape A shapefile (in a polygon type) that matches to your OD dataframe. The shape must have an "id" column to match your ids in df.
#' @param queen A TRUE/FALSE input that is used to calculate \code{spdep}'s spatial contingency (Please view documents of \link[spdep]{poly2nb} for more information).
#' @param snap A parameter that is also used to calculate \code{spdep}'s spatial contingency (Please view documents of \link[spdep]{poly2nb} for more information).
#' @param method A string value among "o" (origin based), "d" (destination based), and "t" (both way) which determines the way to generate Spatial Weights. The default value is "t".
#' @return The result is in the form of a list which includes combinations of origin ids and destination ids.
#' @examples
#' # Data manipulation
#' # Load sf polygon
#' \donttest{
#' CA_polygon <- spnaf::CA_polygon
#' }
#' # Execution of Networknb with data above and given parameters
#' \donttest{
#' nnb <- Networknb(shape = CA_polygon, queen = TRUE, snap = 1, method = 'o')
#' }
#'
#' # check the results
#' \donttest{
#' head(nnb)
#' }


#' @importFrom magrittr %>%
#' @export Networknb

Networknb <- function(shape, snap = 1, queen = TRUE, method = "t"){

    if(!method %in% c("t", "o", "d")){
        stop("method must be one of c('t', 'o', 'd') \n")
    }

    cat("(1) Creating base spatial weights... ")
    w <- NULL
    nb <- unclass(spdep::poly2nb(shape, snap = snap, queen = queen)) # create nb data from input shape
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

    ## Union is created by input polygons^2
    U <- merge(shape$id, shape$id)
    names(U) <- c("oid", "did")
    result <- U %>% dplyr::filter(oid != did) %>% # remove flow from i to i
        dplyr::left_join(nb_frame, by = c("oid", "did")) %>%
        dplyr::mutate(w = ifelse(is.na(w), 0, 1)) %>%
        dplyr::arrange(oid, did)

    rm(nb_frame)
    rm(did)
    rm(oid)
    rm(w)
    rm(original_id)
    rm(nb)
    cat("Done !\n")
    ## Alarm the size of the union set
    cat(paste0("note: Total ", nrow(U), " network combinations are ready to be analyzed\n"))
    cat(paste0("note: which is calculated by the input polygon data\n"))
    rm(U)

    cat("(2) Creating network spatial weights... ")

    SpatialWeightsl <- result %>%
        dplyr::mutate(seq = paste(oid, did, sep="-"))

    SWL <- split(SpatialWeightsl,
                 f = SpatialWeightsl$seq)

    oid <- did <- w <- NULL

    subframe <- function(l, ref = result, m = method){
        o <- l$oid
        d <- l$did

        if(m == 't'){
            origins <- ref %>%
                dplyr::filter(oid == o, w == 1) %>%
                dplyr::select(did) %>% unlist()
            origins <- unique(c(o, origins)) # include o
            destinations <- ref %>%
                dplyr::filter(oid == d, w == 1) %>%
                dplyr::select(did) %>% unlist()
            destinations <- unique(c(d, destinations)) # include d
        }else if(m == 'o'){
            origins <- o
            destinations <- ref %>%
                dplyr::filter(oid == d, w == 1) %>%
                dplyr::select(did) %>% unlist()
            destinations <- unique(c(d, destinations)) # include d
        }else if(m == 'd'){
            origins <- ref %>%
                dplyr::filter(oid == o, w == 1) %>%
                dplyr::select(did) %>% unlist()
            origins <- unique(c(o, origins)) # include o
            destinations <- d
        }

        ## Merge valid networks
        set1 <- ref %>%
            dplyr::filter(oid %in% origins, did == d)
        set2 <- ref %>%
            dplyr::filter(oid == o, did %in% destinations)
        set <- rbind(set1, set2) %>%
            dplyr::distinct() %>%
            dplyr::mutate(ntwkid = paste(oid, did, sep = '-'))
        return(c(set$ntwkid))

    }

    ntwknb <- lapply(SWL, subframe)

    cat("Done !\n")
    ## Alarm the size of the union set
    return(ntwknb)
}


