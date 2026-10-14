make_par_image <- function(binary_raster) {
  
  nr <- nrow(binary_raster)
  nc <- ncol(binary_raster)
  
  list(
    as.matrix(binary_raster, wide = TRUE),
    nc/2,
    nr/2,
    min(nr, nc)/2
  )
}