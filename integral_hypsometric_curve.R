HypsoIntCurve <- function(basins, dem, labelfield, nrow, manexcl){
  # Generate Hypsometric Curve and Hypsometric Integral
  # of each stream of a network
  # Args:
  #   basins:       One SpatialPolygons* object
  #                 If a string is provided, it will be used
  #                 as the name for a vector map in the GRASS location
  #                 containing the basins as areas.
  #   labelfield:   One string with the column field to label the plots.
  #   dem:          One Raster* object.
  #                 If a string is provided, it will be used
  #                 as the name for a raster map in the GRASS location  
  #                 containing the DEM.
  #   nrow:         Number of rows of the facet_wrap of stream profiles.
  #   manexcl:      Vector of elements to manually exclude from the analysis.
  
  # Returns:
  #   GRASS GIS maps of the longest flow path and its tributaries
  # Note:
  #   A GRASS session must be initiated using rgrass7 package
  
  # Error handling
  if (!is.character(labelfield)) {
    stop("Argument prefix must be a character string.")
  }
  if (!length(nrow) == 1) {
    stop("Argument nrow must be a vector with one single positive.")
  }

  # Packages
  require(rgrass7)
  require(sp)
  require(raster)
  require(ggplot2)
  require(DescTools)
  require(directlabels)
  require(scales)
  require(gtools)
  
  # Read the sources
  if(class(basins)=="SpatialPolygonsDataFrame") {
    basins <- basins
  } else {
    if(is.character(basins)) {
      basins <- readVECT(basins)
    }
  }
  if(class(dem)=='RasterLayer') {
    dem <- dem
  } else {
    if(is.character(dem)) {
      dem <- raster(readRAST(dem))
    }
  }
  
  # Exclude fake basins and artifacts
  excl <- which(sapply(area(basins), function(x)  x > prod(res(dem))))
  basins <- basins[excl, ]
  
  # Generate DEMs and data.frames of dimensionless A/Ao and H/Ho, and Hypsometric Integral (AUC)
  index <- gtools::mixedsort(as.character(basins@data[, labelfield]))
  if (exists('manexcl')) {
    index <- index[!index %in% manexcl]
  }
  hypsodfl <- sapply(
    index,
    function(x){
      foobasin <- basins[basins@data[,labelfield] == x,]
      foodem <- mask(crop(dem, extent(foobasin)), foobasin)
      z <- sort(na.omit(foodem[]), decreasing = T)
      df <- data.frame(
        cumarea = rescale(
          cumsum(as.numeric(1:length(z)))*prod(res(dem))
          ),
        height = rescale(z)
      )
      return(df = df)
    },
    simplify = F
  )
  hypsodf <- ldply(hypsodfl, data.frame, .id = labelfield)
  HypsoInt <- ldply(
    sapply(
      index,
      function(x){
        data.frame(hypsoint = AUC(
          hypsodf[hypsodf[,labelfield] == x, 'cumarea'],
          hypsodf[hypsodf[,labelfield] == x, 'height'])
        )
      },
      simplify = F
    ),
    data.frame,
    .id = labelfield
  )
  # Generate the Hypsometric Curve
  p <- ggplot(hypsodf, aes(x = cumarea, y = height)) +
    geom_line(col = 'red', lwd = 1) +
    coord_equal() +
    theme(
      legend.position = "none",
      text = element_text(size = 18),
      panel.background = element_rect(fill = 'white', colour = 'black'),
      panel.grid.major.y = element_line(colour = "grey", linetype = "dashed", size = 0.25),
      strip.background = element_rect(colour = "black", fill = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      strip.text.x = element_text(colour = "white", face = "bold")
    ) +
    scale_x_continuous(breaks = c(0,0.5,1), labels = c(0,0.5,1)) +
    scale_y_continuous(breaks = c(0,0.5,1), labels = c(0,0.5,1)) +
    annotate(
      "text",
      x = 0.9, y = 0.9,
      label = paste0('HI==', round(HypsoInt[,'hypsoint'],2)),
      size = 5,
      hjust = 1,
      parse = T
    ) +
    facet_wrap(paste0('~', labelfield), nrow = nrow) +
    labs(x = 'A/Ao', y = 'H/Ho')
  
  # Returns
  return(
    list(
      DataFrame = hypsodf,
      HypsoInt = HypsoInt,
      HypsoCurve = p
    )
  )
}
