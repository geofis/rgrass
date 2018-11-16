TrasnTopoSym <- function(direction, xycoords, suffix){
  # Calculate the Transverse Topographic Basin Symmetry (T) vector
  # of each stream of a network
  # Args:
  #   direction:    Flow direction GRASS raster map. May be generated with r.stream*.
  #   xycoords:     One vector with the coordinates of the basin outlet.
  #                 X and Y coordinates must be placed in positions 
  #                 1 and 2 of the vector, respectively
  #   suffix:       One string for the names of the GRASS maps and the shapefile.
  # Returns:
  #   Values and GRASS GIS maps of the magnitude and bearing of T vector of a basin,
  #   as well as the centerline of the basin and a smoothed version of the longest flow path
  # Note:
  #   A GRASS session must be initiated using rgrass7 package
  
  # Error handling
  if (!is.character(direction) || !is.character(suffix)) {
    stop("Arguments direction and suffix must be character strings.")
  }
  if (!length(xycoords)==2) {
    stop("Argument xycoords must be a vector with two values.")
  }
  if (TRUE %in% is.na(xycoords)) {
    stop("Argument xycoords must not have missing values.")
  }

  #Packages
  require(rgrass7)
  require(sp)
  require(rgdal)
  require(devtools)
  if (!'pT' %in% ls(envir = .GlobalEnv)) {
    source_url('https://raw.githubusercontent.com/geofis/rgrass/master/pT.R')
  }
  
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
  
  # Generate the watershed basin from outlet coordinates and the drainage direction
  execGRASS(
    "r.water.outlet",
    flags='overwrite',
    parameters = list(
      input = direction,
      coordinates = as.vector(
        ParraSubasins[ParraSubasins$name=='Naranjal',c('x','y')],
        mode = 'numeric'),
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
  
  # Generate the centerline of the watershed basin
  execGRASS(
    "v.voronoi",
    flags = c('overwrite','s'),
    parameters = list(
      input = 'tts_tmp_watershed',
      output = 'tts_tmp_watershed_centerline'
    )
  )
  
  # Generalize the longest flow path line vector
  execGRASS(
    "v.generalize",
    flags='overwrite',
    parameters = list(
      input = 'tts_tmp_lfp',
      output = 'tts_tmp_lfp_smoothed',
      method = 'reumann',
      threshold = gmeta()$nsres
    )
  )
  
  # Convert the longest flow path line to a point vector
  # based on the vertices of the generalized flow path
  execGRASS(
    "v.to.points",
    flags = c('overwrite','i','t'),
    parameters = list(
      input = 'tts_tmp_lfp_smoothed',
      output = 'tts_tmp_lfp_pnt',
      use = 'vertex',
      dmax = gmeta()$nsres
    )
  )
  
  # Add an attribute table to the LFP point vector with the fields
  # that will be later populated with the distances to
  # basin centerline and basin margin
  execGRASS(
    "v.db.addtable",
    parameters = list(
      map = 'tts_tmp_lfp_pnt',
      layer = 2,
      columns = paste(
        'to_centerline_x double precision',
        'to_centerline_y double precision',
        'to_basinmargin_x double precision',
        'to_basinmargin_y double precision',
        sep = ','
      )
    )
  )
  
  # Calculate shortest distances from vertices of LFP point layer 
  # to the basin centerline
  execGRASS(
    "v.distance",
    flags = 'overwrite',
    parameters = list(
      from = 'tts_tmp_lfp_pnt',
      from_layer = '2',
      from_type = 'point',
      to = 'tts_tmp_watershed_centerline',
      to_type = 'line',
      upload = paste('to_x', 'to_y', sep = ','),
      column = paste('to_centerline_x', 'to_centerline_y', sep = ',')
    )
  )
  
  # Calculate shortest distances from vertices of LFP point layer 
  # to the basin margin
  execGRASS(
    "v.distance",
    flags = 'overwrite',
    parameters = list(
      from = 'tts_tmp_lfp_pnt',
      from_layer = '2',
      from_type = 'point',
      to = 'tts_tmp_watershed_lines',
      to_type = 'line',
      upload = paste('to_x', 'to_y', sep = ','),
      column = paste('to_basinmargin_x', 'to_basinmargin_y', sep = ',')
    )
  )
  
  # Importing data into R
  distxynear <- readVECT('tts_tmp_lfp_pnt', layer = 2)
  distxynear@data[,c('from_x','from_y')] <- distxynear@coords
  
  #Magnitude
  ##Distance from active meander to the basin midline (Da)
  distxynear@data$Da <- with(
    distxynear@data,
    sqrt((to_centerline_x - from_x) ^ 2 + (to_centerline_y - from_y) ^ 2)
  )
  
  ##Distance from the active meander to the basin divide
  distxynear@data$distbasinmargin <- with(
    distxynear@data,
    sqrt((to_basinmargin_x - from_x) ^ 2 + (to_basinmargin_y - from_y) ^ 2)
  )
  
  ##Distance from the basin divide to the basin midline (Dd)
  distxynear@data$Dd <- with(
    distxynear@data,
    Da + distbasinmargin
  )
  
  ##Transverse Topography Symmetry (T or TTS)
  distxynear@data$TTS <- with(
    distxynear@data,
    Da/Dd
  )
  
  #Bearing of the Direction
  #Radians
  distxynear@data$bearingrad <- with(
    distxynear@data,
    atan2((to_centerline_y - from_y),(to_centerline_x - from_x))
  )
  #Degrees
  distxynear@data$bearingdeg <- with(
    distxynear@data,
    atan2((to_centerline_y - from_y),(to_centerline_x - from_x)) * 180 / pi
  )
  
  #Assign an iterative sequence 1 to 5 based on cat numbers
  #as an index for subsetting points based on a somewhat constant distance
  distxynear@data$seq <- (distxynear@data$cat - 1) %% 5 + 1
  
  #Mean value and probabality of change by chance of the magnitude of
  #Transverse Topographic Basin Symmetry (T) for each index number
  mandpTbygrp <- sapply(
    1:5,
    function(x){
      tts <- distxynear@data[distxynear@data$seq == x, 'TTS']
      p <- pT(tts)$p
      m <- mean(tts)
      paste0(
        'Group ', x, ': ',
        'm=', round(m*100, 2), '%. ',
        'p=', round(p, 2), '\n'
      )
    }
  )
  
  #Overal mean value and probabality of change by chance of the magnitude of
  #Transverse Topographic Basin Symmetry (T)
  ttsov <- distxynear@data$TTS
  pov <- pT(ttsov)$p
  mov <- mean(ttsov)
  mandpTov <- paste0(
    'm=', round(mov*100, 2), '%. ',
    'p=', round(pov, 2), '\n'
  )
  
  #Distance between consecutive points
  disconpnt <- sapply(
    unique(distxynear@data$seq),
    function(x){
      subsetted <- distxynear@coords[distxynear@data$seq == x,]
      mins <- sapply(
        2:nrow(subsetted),
        function(i)
          dist(subsetted[c(i-1,i),])
      )
      paste0(
        'Group ', x, ': ',
        'min=', round(range(mins)[1], 2), ' m. ',
        'max=', round(range(mins)[2], 2), ' m\n'
      )
    }
  )
  
  #Exporting results
  ##Entire sample
  writeVECT(
    distxynear,
    paste0(
      'tts_TTS_magnitude_bearing_entire_sample_',
      suffix
    ),
    v.in.ogr_flags = 'overwrite'
  )
  
  ##By groups
  sapply(
    unique(distxynear@data$seq),
    function(x){
      writeVECT(
        distxynear[distxynear@data$seq == x,],
        paste0(
          'tts_TTS_magnitude_bearing_group_',
          x,
          '_',
          suffix
        ),
        v.in.ogr_flags = 'overwrite'
      )
    }
  )
  
  #Rename GRASS maps to be saved
  execGRASS(
    "g.rename",
    flags = 'overwrite',
    parameters = list(
      vector = paste0(
        'tts_tmp_lfp_smoothed,',
        'tts_lfp_smoothed_',
        suffix
      )
    )
  )
  
  execGRASS(
    "g.rename",
    flags = 'overwrite',
    parameters = list(
      vector = paste0(
        'tts_tmp_lfp_pnt,',
        'tts_lfp_pnt_',
        suffix
      )
    )
  )
  
  execGRASS(
    "g.rename",
    flags = 'overwrite',
    parameters = list(
      vector = paste0(
        'tts_tmp_watershed_centerline,',
        'tts_lfp_watershed_centerline_',
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
      pattern = 'tts_tmp*'
    )
  )
  
  #Returns:
  return(
    list(
      mandpTov = capture.output(
        cat(
          'Overall mean value and probability of change by chance of the magnitude of T',
          '\n',
          mandpTov,
          sep = ''
        )
      ),
      mandpTbygrp = capture.output(
        cat(
          'Mean value and probability of change by chance of the magnitude of T by group',
          '\n',
          mandpTbygrp,
          sep = ''
          )
        ),
      disconpnt = capture.output(
        cat(
          'Minimum and maximum distances between consecutive points by group',
          '\n',
          disconpnt,
          sep = ''
          )
        ),
      SpatialObjectWithAttributes = distxynear
    )
  )
}
