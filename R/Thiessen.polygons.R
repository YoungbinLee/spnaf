
#' @importFrom deldir deldir
#' @importFrom deldir tile.list
# Functions ---------------------------------------------------------------

#
# library(deldir)
# library(sp)
# library(dplyr)

create_thiessen_polygons <- function(points) {
  voronoi <- deldir(st_coordinates(points)[,1], st_coordinates(points)[,2])
  voronoi_polygons <- tile.list(voronoi)

  close_polygon <- function(coords) {
    if (!all(coords[1,] == coords[nrow(coords),])) {
      coords <- rbind(coords, coords[1,])
    }
    coords
  }

  polys <- lapply(1:length(voronoi_polygons), function(i) {
    coords <- cbind(voronoi_polygons[[i]]$x, voronoi_polygons[[i]]$y)
    coords <- close_polygon(coords)
    st_polygon(list(coords))
  })
  voronoi_sf <- st_sfc(polys, crs = 4326)
  voronoi_sf <- st_sf(id = points$id, geometry = voronoi_sf)
  return(voronoi_sf)

}
