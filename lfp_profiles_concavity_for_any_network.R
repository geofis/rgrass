AnyNetworkProfilesConcavity <- function(xycoords, network, prefix, dem, direction,
                                 crs = '+init=epsg:32619', smns = 0.5, nrow = 4){
  # Generate the profile curves and concavity indices
  # of each stream of a network
  # Args:
  #   xycoords:     One vector with the coordinates of the basin outlet.
  #   network:      The stream network for which profiles may be generated.
  #   prefix:       One string for the prefix of the profile labels.
  #   dem:          DEM from which to extract z for profiles.
  #   direction:    Flow direction raster map. May be generated with r.stream*.
  #   crs:          CRS of the xycoords object.
  #   smns:         Smoothing parameter for stats::spline function.
  #   nrow:         Number of rows of the facet_wrap of stream profiles.
  # Returns:
  #   GRASS GIS maps of the longest flow path and its tributaries
  # Note:
  #   A GRASS session must be initiated using rgrass7 package
  
  # Error handling
  if (!length(xycoords)==2) {
    stop("Argument xycoords must be a vector with two values.")
  }
  if (TRUE %in% is.na(xycoords)) {
    stop("Argument xycoords must not have missing values.")
  }
  if (!is.character(prefix)) {
    stop("Argument prefix must be a character string.")
  }
  if (!is.character(network) || !is.character(dem) || !is.character(direction)) {
    stop("Arguments network, dem and direction must be character strings.")
  }
  if (!length(smns) == 1) {
    stop("Argument smns must be a vector with one single value.")
  }
  if (!(smns >= 0 & smns <= 1)) {
    stop("Argument smns must be a vector with a positive value between 0 and 1.")
  }
  if (!length(nrow) == 1) {
    stop("Argument nrow must be a vector with one single positive.")
  }
  
  # Packages
  require(rgrass7)
  require(sp)
  require(raster)
  require(plyr)
  require(ggplot2)
  require(RColorBrewer)
  require(DescTools)
  require(directlabels)
  require(scales)
  
  # Generate outlet raster from coordinates
  sp <- data.frame(x=xycoords[1], y=xycoords[2], prefix = prefix)
  coordinates(sp) <- ~x+y
  proj4string(sp) <- CRS(crs)
  writeVECT(
    sp,
    paste0('AnyNetwork_outlet_', prefix),
    v.in.ogr_flags = 'overwrite'
  )
  execGRASS(
    "v.to.rast",
    flags='overwrite',
    parameters = list(
      input = paste0('AnyNetwork_outlet_', prefix),
      output = paste0('AnyNetwork-outlet-', prefix),
      use = 'cat'
    )
  )
  
  # Generate the stream distance raster map
  execGRASS(
    "r.stream.distance",
    flags='overwrite',
    parameters = list(
      stream_rast = paste0('AnyNetwork-outlet-', prefix),
      direction = direction,
      distance = paste0('AnyNetwork-flds-', prefix),
      method = 'downstream'
    )
  )
  
  #Bring to R both rasters the DEM and the distance map, and also the LFP vector map
  lfpall <- readVECT(network)
  length <- raster(readRAST(paste0('AnyNetwork-flds-', prefix)))
  names(length) <- 'length'
  z <- raster(readRAST(dem))
  names(z) <- 'z'
  stk <- stack(length, z)
  
  # Extract values from raster stack
  listofdfs <- raster::extract(stk, lfpall)
  names(listofdfs) <- paste0(prefix, '-', as.character(lfpall$cat))
  dfs <- plyr::ldply(listofdfs, data.frame, .id = 'stream')
  dfs <- na.omit(dfs)
  dfs <- droplevels(dfs[dfs$stream %in% names(which(table(dfs$stream)>3)),])
  
  # Spline smoothing
  listofdfs2 <- sapply(
    as.character(unique(dfs$stream)),
    function(x)
      as.data.frame(
        smooth.spline(dfs[dfs$stream==x,c('length','z')], spar = smns)[c('x','y')]),
    simplify = F,
    USE.NAMES = T
  )
  # names(listofdfs2) <- names(listofdfs)
  # return(listofdfs2)
  dfs2 <- plyr::ldply(listofdfs2, data.frame, .id = 'stream')
  names(dfs2)[2:3] <- c('length', 'z')

  # Delete points which surpassed the limits of their own subbasin
  # Not an elegant solution, but for now it works
  dfs2[,'rowname'] <- as.integer(rownames(dfs2))
  lengthdif <- plyr::ldply(
    sapply(
      levels(dfs2[,'stream']),
      function(x)
        t(
          sapply(
            1:nrow(dfs2[dfs2[,'stream']==x,]),
            function(y){
              dif <- ifelse(
                y==1,
                0,
                dfs2[dfs2[,'stream']==x,][y,'length'] -
                  dfs2[dfs2[,'stream']==x,][y-1,'length']
              )
              rowname <- dfs2[dfs2[,'stream']==x,][y,'rowname']
              return(c(dif = dif, rowname = rowname))
            }
          )
        )
    ),
    data.frame,
    .id = 'stream'
  )
  deldfs2 <- as.vector(
    unlist(
      sapply(
        levels(lengthdif[,'stream']),
        function(x){
          lengthdifx <- lengthdif[lengthdif[,'stream']==x, ]
          rownamex <- lengthdifx[,'rowname']
          rownamegreater <- lengthdifx[, 'dif'] > sqrt(gmeta()$nsres^2 + gmeta()$ewres^2)
          rownamegreater2 <- lengthdifx[rownamegreater,'rowname']
          rownamex[rownamex >= rownamegreater2]
        }
      )
    )
  )
  if(length(deldfs2)>0) dfs3 <- dfs2[-deldfs2,] else dfs3 <- dfs2

  # Profiles in one single plot
  longeststream <- dfs3[which.max(dfs3[,'length']),'stream']
  ncolors <- length(levels(dfs3[,'stream']))
  p <- ggplot(dfs3, aes(x=length, y=z, group=stream)) +
    geom_line(lwd = 1) +
    geom_line(
      data = dfs3[dfs3[,'stream']==longeststream,],
      aes(x=length, y=z), col = 'black', lwd = 2) +
    scale_color_manual(values = rep(brewer.pal(8,"Dark2"), length.out = ncolors)) +
    # geom_dl(aes(label=stream), method="last.points") +
    theme(
      legend.position="none",
      text = element_text(size = 18),
      panel.background = element_rect(fill = 'white', colour = 'black'),
      panel.grid.major.y = element_line(colour = "grey", linetype = "dashed", size = 0.25),
      strip.background = element_rect(colour = "black", fill = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) +
    labs(x = 'Length (m)', y = 'Altitude (m)')

  # Dimensionless tables, profile graphs
  dmnls <- sapply(
    levels(dfs3[,'stream']),
    function(x){
      dmnls <- dfs3[dfs3[,'stream']==x,]
      dmnls[,c('length.dmnls','z.dmnls')] <- sapply(dmnls[,c('length','z')], rescale)
      dmnls
    },
    simplify = F
  )
  dmnls <- plyr::ldply(dmnls, data.frame, .id = 'stream')
  ci <- plyr::ldply(sapply(
    levels(dmnls[,'stream']),
    function(x) {
      data.frame(ci = (0.5 -
                         AUC(
                           dmnls[dmnls['stream']==x, 'length.dmnls'],
                           dmnls[dmnls['stream']==x, 'z.dmnls'])) / 0.5
      )
    },
    simplify = F
  ), data.frame, .id = 'stream')
  pdmnls <- ggplot(dmnls, aes(x = length.dmnls, y = z.dmnls)) +
    geom_line(col = 'red', lwd = 1) +
    coord_equal() +
    facet_wrap(~stream, nrow = nrow) +
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
      x = 0.1, y = 0.9,
      label = paste0('C[a]==', round(ci[,'ci'],2)),
      size = 5,
      hjust = 0,
      parse = T
    ) +
    labs(x = 'L/Lo', y = 'H/Ho')
  # Returns
  return(
    list(
      lengthzdata = dfs3,
      profiles = p,
      lengthzdatadmnls = dmnls,
      concavityindex = ci,
      dimensionlessprofiles = pdmnls
    )
  )
}
