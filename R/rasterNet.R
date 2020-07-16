#' @import sp raster
rasterNet <- function (x, resolution = NULL, xbin = NULL, ybin = NULL, mask = FALSE,
                       degree = 111325, xOffset = NULL, yOffset = NULL, checkerboard = FALSE,
                       maxpixels = 250000)
{
  ext <- raster::extent(x)
  extRef <- raster::extent(x)
  if (is.na(sp::proj4string(x))) {
    mapext <- raster::extent(x)[1:4]
    if (mapext >= -180 && mapext <= 180) {
      resolution <- resolution/degree
      # warning("The input layer has no CRS defined. Based on the extent of the input map it is assumed to have an un-projected reference system")
    }
    else {
      resolution <- resolution
      # warning("The input layer has no CRS defined. Based on the extent of the input map it is assumed to have a projected reference system")
    }
  }
  else {
    if (sp::is.projected(sp::SpatialPoints((matrix(1:10,
                                                   5, byrow = FALSE)), proj4string = raster::crs(x)))) {
      resolution <- resolution
    }
    else {
      resolution <- resolution/degree
    }
  }
  if (!is.null(xbin) && is.null(ybin)) {
    rasterNet <- raster::raster(ext, nrow = 1, ncol = xbin,
                                crs = raster::crs(x))
  }
  else if (is.null(xbin) && !is.null(ybin)) {
    rasterNet <- raster::raster(ext, nrow = ybin, ncol = 1,
                                crs = raster::crs(x))
  }
  else if (!is.null(xbin) && !is.null(ybin)) {
    rasterNet <- raster::raster(ext, nrow = ybin, ncol = xbin,
                                crs = raster::crs(x))
  }
  else if (is.null(xbin) && is.null(ybin) && !is.null(resolution)) {
    xrange <- raster::xmax(x) - raster::xmin(x)
    yrange <- raster::ymax(x) - raster::ymin(x)
    xPix <- ceiling(xrange/resolution)
    yPix <- ceiling(yrange/resolution)
    xdif <- ((xPix * resolution) - xrange)/2
    ydif <- ((yPix * resolution) - yrange)/2
    ext@xmin <- raster::xmin(x) - xdif
    ext@xmax <- raster::xmax(x) + xdif
    ext@ymin <- raster::ymin(x) - ydif
    ext@ymax <- raster::ymax(x) + ydif
    if (!is.null(xOffset)) {
      if (xOffset > 1 || xOffset < 0) {
        stop("xOffset should be between 0 and 1")
      }
      ext@xmin <- ext@xmin + (resolution * xOffset)
      ext@xmax <- ext@xmax + (resolution * xOffset)
    }
    if (!is.null(yOffset)) {
      if (yOffset > 1 || yOffset < 0) {
        stop("yOffset should be between 0 and 1")
      }
      ext@ymin <- ext@ymin + (resolution * yOffset)
      ext@ymax <- ext@ymax + (resolution * yOffset)
    }
    if (ext@xmin > extRef@xmin) {
      ext@xmin <- ext@xmin - resolution
      xPix <- xPix + 1
    }
    if (ext@ymin > extRef@ymin) {
      ext@ymin <- ext@ymin - resolution
      yPix <- yPix + 1
    }
    rasterNet <- raster::raster(ext, nrow = yPix, ncol = xPix,
                                crs = raster::crs(x))
  }
  else stop("A value should be specified for the block size")
  if (checkerboard == TRUE) {
    values(rasterNet) <- 1:ncell(rasterNet)
    m <- as.matrix(rasterNet)
    for (i in 1:ncol(rasterNet)) {
      if (i%%2 == 0) {
        m[, i] <- rep(1:2, nrow(m))[1:nrow(m)]
      }
      else {
        m[, i] <- rep(2:1, nrow(m))[1:nrow(m)]
      }
    }
    rasterNet[] <- m
  }
  else {
    values(rasterNet) <- 1:ncell(rasterNet)
  }
  rasterNet <- raster::rasterToPolygons(rasterNet)
  if (mask == TRUE) {
    if (methods::is(x, "Raster")) {
      points <- raster::rasterToPoints(x[[1]], spatial = TRUE)
      if (nrow(points) > 750000) {
        points2 <- points[sample(1:nrow(points), maxpixels,
                                 replace = FALSE), ]
        rasterNet <- raster::intersect(rasterNet, points2)
      }
      else rasterNet <- raster::intersect(rasterNet, points)
    }
    else {
      rasterNet <- raster::intersect(rasterNet, x)
    }
  }
  return(rasterNet)
}
