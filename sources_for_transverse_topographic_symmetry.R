SourcesForTransTopoSym <- function(direction, xycoords, prefix,
                                   smthpathf = 3, smthwatshdf = 3, thin = -1,
                                   tspacing = 100, fdistft = 1){
  # Generates the sources for calculating the Transverse Topographic Basin Symmetry (T)
  # vector of a stream path given the coordinates of a basin
  # Args:
  #   direction:    Flow direction GRASS raster map. May be generated with r.stream*.
  #   xycoords:     One vector with the coordinates of the basin outlet.
  #                 X and Y coordinates must be placed in positions 
  #                 1 and 2 of the vector, respectively
  #   prefix:       One string for the names of the GRASS maps and the shapefile.
  #   smthpathf:    Numeric. Factor to which multiply region resolution.
  #                 resolution*smthpathf is passed a to v.generalize as threshold
  #                 parameter. The threshold parameter, given in map units is the degree of
  #                 simplification increases with the increasing value of threshold.
  #                 This important parameter controls whether the transects (perpendicular to
  #                 the path) will be meaningful or not.
  #   smthwatshdf:   Numeric. Factor to which multiply region resolution.
  #                 resolution*smthwatshd is passed a to v.generalize as threshold
  #                 parameter. The threshold parameter, given in map units is the degree of
  #                 simplification increases with the increasing value of threshold.
  #   ptrim:        Proportion of vertices with respect of the total to trim on each endpoint
  #   thin:         Maximum dangle length of skeletons. Argument passed to GRASS > v.voronoi
  #                 The default value (-1) extracts the centerline
  #   tspacing:     Separation between transects (in map units)
  #   fdistft:      Factor to multiply the length of the longest flow path. The product is used
  #                 as the distance transect extends to the left and the right of the flow path.
  # Returns:
  #   Values and GRASS GIS maps of the magnitude and bearing of T vector of a basin,
  #   as well as the centerline of the basin and a smoothed version of the longest flow path
  # Note:
  #   A GRASS session must be initiated using rgrass7 package
  
  # Error handling
  if (!is.character(direction) || !is.character(prefix)) {
    stop("Arguments direction and prefix must be character strings.")
  }
  if (!length(xycoords)==2) {
    stop("Argument xycoords must be a vector with two values.")
  }
  if (TRUE %in% is.na(xycoords)) {
    stop("Argument xycoords must not have missing values.")
  }

  #Packages
  require(rgrass7)
  require(rgeos)
  require(sp)

  # Generates the longest flow path, which is considered the active meanderbelt
  execGRASS(
    "r.lfp",
    flags='overwrite',
    parameters = list(
      input = direction,
      coordinates = xycoords,
      output = 'tts_tmp_lfp'
    )
  )
  
  #If r.lfp generates more than one path, randomly select one of them
  if(vInfo('tts_tmp_lfp')['lines']>1){
    execGRASS(
      'v.category',
      flags = 'overwrite',
      parameters = list(
        input = 'tts_tmp_lfp',
        output = 'tts_tmp_lfp_catdel',
        option = 'del',
        cat = -1
      )
    )
    execGRASS(
      'v.category',
      flags = 'overwrite',
      parameters = list(
        input = 'tts_tmp_lfp_catdel',
        output = 'tts_tmp_lfp_catadd',
        option = 'add',
        cat = 1
      )
    )
    execGRASS(
      'v.extract',
      flags = 'overwrite',
      parameters = list(
        input = 'tts_tmp_lfp_catadd',
        output = 'tts_tmp_lfp',
        random = 1
      )
    )
  }
  
  # Generalize the longest flow path line vector
  execGRASS(
    "v.generalize",
    flags='overwrite',
    parameters = list(
      input = 'tts_tmp_lfp',
      output = 'tts_tmp_lfp_smoothed',
      method = 'reumann',
      threshold = gmeta()$nsres*smthpathf
    )
  )

  # Generate the watershed basin from outlet coordinates and the drainage direction
  execGRASS(
    "r.water.outlet",
    flags='overwrite',
    parameters = list(
      input = direction,
      coordinates = xycoords,
      output = "tts_tmp_watershed"
    )
  )

  # Convert the watershed basin to vector format
  execGRASS(
    "r.to.vect",
    flags = 'overwrite',
    parameters = list(
      input = 'tts_tmp_watershed',
      output = 'tts_tmp_watershed',
      type = 'area'
    )
  )

  # Convert the watershed basin to line vector format
  execGRASS(
    "v.to.lines",
    flags = 'overwrite',
    parameters = list(
      input = 'tts_tmp_watershed',
      output = 'tts_tmp_watershed_lines'
    )
  )
  
  # Convert the watershed basin to line vector format
  execGRASS(
    "v.to.lines",
    flags = 'overwrite',
    parameters = list(
      input = 'tts_tmp_watershed',
      output = 'tts_tmp_watershed_lines'
    )
  )
  
  # Generalize the longest flow path line vector
  execGRASS(
    "v.generalize",
    flags='overwrite',
    parameters = list(
      input = 'tts_tmp_watershed_lines',
      output = 'tts_tmp_watershed_lines_smoothed',
      method = 'chaiken',
      threshold = gmeta()$nsres*smthwatshdf
    )
  )
  # Convert to polygon the smoothed watershed polyline
  execGRASS(
    "v.type",
    flags='overwrite',
    parameters = list(
      input = 'tts_tmp_watershed_lines_smoothed',
      output = 'tts_tmp_watershed_lines_smoothed_boundary',
      from_type = 'line',
      to_type = 'boundary'
    )
  )
  # Add the centroid to the smoothed watershed polygon
  execGRASS(
    "v.centroids",
    flags='overwrite',
    parameters = list(
      input = 'tts_tmp_watershed_lines_smoothed_boundary',
      output = 'tts_tmp_watershed_lines_smoothed_polygon_with_centroid'
    )
  )
  
  # Skeletonize the watershed basin
  execGRASS(
    "v.voronoi",
    flags = c('overwrite','s'),
    parameters = list(
      input = 'tts_tmp_watershed_lines_smoothed_polygon_with_centroid',
      output = 'tts_tmp_skeleton',
      thin = thin
    )
  )
  
  lfp <- readVECT('tts_tmp_lfp_smoothed')
  distfortransects <- gLength(lfp)
  
  # Generate the transects
  execGRASS(
    "v.transects",
    flags = 'overwrite',
    parameters = list(
      input = 'tts_tmp_lfp_smoothed',
      output = 'tts_tmp_transects',
      transect_spacing = tspacing,
      dleft = distfortransects*fdistft,
      dright = distfortransects*fdistft
    )
  )
  
  
  
  #Rename GRASS maps to be saved
  ##The transects
  execGRASS(
    "g.rename",
    flags = 'overwrite',
    parameters = list(
      vector = paste0(
        'tts_tmp_transects,',
        'tts_', prefix, '_transects'
      )
    )
  )
  
  ##The watershed divide
  execGRASS(
    "g.rename",
    flags = 'overwrite',
    parameters = list(
      vector = paste0(
        'tts_tmp_watershed_lines_smoothed,',
        'tts_', prefix, '_wshed_divide'
      )
    )
  )
  
  ##The centerline
  execGRASS(
    "g.rename",
    flags = 'overwrite',
    parameters = list(
      vector = paste0(
        'tts_tmp_skeleton,',
        'tts_', prefix, '_centerline'
      )
    )
  )
  
  ##The LFP or stream
  execGRASS(
    "g.rename",
    flags = 'overwrite',
    parameters = list(
      vector = paste0(
        'tts_tmp_lfp_smoothed,',
        'tts_', prefix, '_stream'
      )
    )
  )
  
  #Remove tmp maps
  execGRASS(
    "g.remove",
    flags = 'f',
    parameters = list(
      type = c('raster','vector'),
      pattern = 'tts_tmp*'
    )
  )
  
}
