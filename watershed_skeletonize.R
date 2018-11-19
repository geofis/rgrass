WatershedSkeletonize <- function(suffix, 
                                 direction = NULL, xycoords = NULL,
                                 watershed = NULL, thin = NULL){
  # Skeletonize watershed. Watershed may be provided by the user
  # or may be generated with a direction map and a pair of coordinates
  # Args:
  #   direction:    Name of the flow direction GRASS raster map. May be generated with r.stream*.
  #                 If direction map is not provided, the function attempts
  #                 to skeletonize a watershed supplied
  #   xycoords:     One vector with the coordinates of the basin outlet.
  #                 X and Y coordinates must be placed in positions 
  #                 1 and 2 of the vector, respectively
  #                 If no coordinates are provided, the function attempts
  #                 to skeletonize a watershed supplied
  #   watershed:    Name of the watershed GRASS vector polygon map name.
  #                 If no watershed is provided, the function attemps
  #                 to create a watershed from direction map and xycoords
  #   suffix:       One string for the names of the GRASS maps and the shapefile.
  #   thin:         Maximum dangle length of skeletons. Argument passed to GRASS
  # Returns:
  #   GRASS GIS maps of the skeleton of the basin
  # Note:
  #   A GRASS session must be initiated using rgrass7 package
  
  # Error handling
  # if (!is.character(direction) || !is.character(watershed) || !is.character(suffix)) {
  #   stop("Arguments direction, watershed and suffix must be character strings.")
  # }
  if (!length(xycoords)==2) {
    stop("Argument xycoords must be a vector with two values.")
  }
  if (TRUE %in% is.na(xycoords)) {
    stop("Argument xycoords must not have missing values.")
  }
  if (is.null(thin) || length(thin)>1 || thin<=0) {
    stop("Argument thin must be one positive real number.")
  }
  if (!is.null(watershed)) {
    watershedprovided <- TRUE
  } else {
    if  (!is.null(direction) & !is.null(xycoords)) {
      watershedprovided <- FALSE
    } else {
      stop('A watershed polygon map name or a direction map and a pair of coordinates must be provided')
    }
  }
  
  #Packages
  require(rgrass7)
  require(sp)
  
  if (!watershedprovided) {
    # Generate the watershed basin from outlet coordinates and the drainage direction
    execGRASS(
      "r.water.outlet",
      flags='overwrite',
      parameters = list(
        input = direction,
        coordinates = xycoords,
        output = "skeleton_tmp_watershed"
      )
    )
    
    # Convert the watershed basin to vector format
    execGRASS(
      "r.to.vect",
      flags = 'overwrite',
      parameters = list(
        input = 'skeleton_tmp_watershed',
        output = 'skeleton_tmp_watershed',
        type = 'area'
      )
    )
  } else {
    if (watershedprovided) {
      execGRASS(
        "g.copy",
        flags = 'overwrite',
        parameters = list(
          vector = paste0(watershed,',skeleton_tmp_watershed')
        )
      )
    }
  }
  
  # Convert the watershed basin to line vector format
  execGRASS(
    "v.to.lines",
    flags = 'overwrite',
    parameters = list(
      input = 'skeleton_tmp_watershed',
      output = 'skeleton_tmp_watershed_lines'
    )
  )
  
  # Skeletonize the watershed basin
  execGRASS(
    "v.voronoi",
    flags = c('overwrite','s'),
    parameters = list(
      input = 'skeleton_tmp_watershed',
      output = 'skeleton_tmp_skeleton',
      thin = thin
    )
  )
  
  # Generalize the longest flow path line vector
  execGRASS(
    "v.generalize",
    flags='overwrite',
    parameters = list(
      input = 'skeleton_tmp_skeleton',
      output = 'skeleton_tmp_skeleton_generalized',
      method = 'reumann',
      threshold = gmeta()$nsres
    )
  )
  
  #Rename GRASS maps to be saved
  ##The skeleton
  execGRASS(
    "g.rename",
    flags = 'overwrite',
    parameters = list(
      vector = paste0(
        'skeleton_tmp_skeleton_generalized,',
        'skltnfun_the_skeleton_',
        suffix
      )
    )
  )
  ##The watershed polygon
  if (!watershedprovided) {
    execGRASS(
      "g.rename",
      flags = 'overwrite',
      parameters = list(
        vector = paste0(
          'skeleton_tmp_watershed,',
          'skltnfun_watershed_',
          suffix
        )
      )
    )
  }
  
  ##The watershed boundary
  execGRASS(
    "g.rename",
    flags = 'overwrite',
    parameters = list(
      vector = paste0(
        'skeleton_tmp_watershed_lines,',
        'skltnfun_watershed_boundary_',
        suffix
      )
    )
  )
  
  #Remove tmp maps
  execGRASS(
    "g.remove",
    flags = 'f',
    parameters = list(
      type = c('raster','vector'),
      pattern = 'skeleton_tmp_*'
    )
  )
}
