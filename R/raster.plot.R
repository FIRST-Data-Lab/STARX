#' STARX: spatiotemporal autoregressive model with exogenous variables
#'
#' This function generates raster plot using ggplot2.
#' @import ggplot2
#' @importFrom colorRamps matlab.like
#' @param data 	The data to be displayed with columns:
#' \code{x}, \code{y}, and \code{value}.
#' @param value.range The range of value displayed on the plot.
#' @param boundary Boundary of the plot. The value out of boundary is \code{NA}.
#' @param coord.ratio Aspect ratio, expressed as y / x.
#' @export
raster.plot <- function (data, value.range, boundary, coord.ratio = 1){
  # input check
  if( sum(c("x", "y", "value") %in% names(data)) != 3) {
    stop("Please check the column names of data.")
  }
  if(sum(c("V1", "V2") %in% names(boundary)) != 2) {
    stop("Please check the column names of boundary.")
  }

  # generate plots
  p <- ggplot(data = data, aes(x = x, y = y)) +
    geom_raster(aes(fill = value)) +
    scale_fill_gradientn(colours = colorRamps::matlab.like(104),
                         limits = value.range, na.value = 'transparent') +
    coord_fixed(ratio = coord.ratio)+
    theme_bw()+theme(axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), legend.title=element_blank(),
                     text=element_text(size = 20),legend.key.size = unit(1.5, "cm"),
                     legend.key.width = unit(0.6,"cm"),
                     legend.box.spacing = unit(0.4, "cm"),
                     legend.spacing.x = unit(0.2, 'cm')) +
    geom_polygon(data = boundary, aes(V1, V2), fill=NA, colour="black", size = 1)

  p
}
