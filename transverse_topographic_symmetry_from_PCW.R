TransTopoSymfromPCW <- function(path, centerline, watershed, suffix,
                                grps = 5, smthpath = F, smthcnt = F, ptrim = 0.04){
  # Calculate the Transverse Topographic Basin Symmetry (T) vector
  # of a stream path providing path, centerline and watershed polygon
  # Args:
  #   path:         GRASS vector map name of the stream path(s).
  #                 May be generated with r.lfp or r.stream*
  #   centerline:   GRASS vector map name of the centerline of the watershed basin divide.
  #                 May be generated with v.voronoi or with the
  #                 custom function WatershedSkeletonize
  #   watershed:    GRASS vector map name of the watershed basin divide as a line
  #   suffix:       One string for the names of the GRASS maps and the shapefile.
  #   smthpath:     Logical. Smoothes the path.
  #   smthcnt:      Logical. Smoothes the centerline.
  #   grps:         Integer. Number of groups in which to split the sample.
  #                 Must be greater than 1
  #   ptrim:        Proportion of vertices with respect of the total to trim on each endpoint
  # Returns:
  #   Values and GRASS GIS maps of the magnitude and bearing
  #   of T vector of a watershed basin
  # Note:
  #   A GRASS session must be initiated using rgrass7 package
  
  # Error handling
  if (!is.character(path) || !is.character(centerline) ||
      !is.character(watershed) || !is.character(suffix)) {
    stop("Arguments path, centerline, watershed and suffix must be character strings.")
  }
  if (grps<=1) {
    stop("Argument grps must be greater than 1.")
  }
  
  
  #Packages
  require(rgrass7)
  require(sp)
  require(rgdal)
  require(devtools)
  if (!'pT' %in% ls(envir = .GlobalEnv)) {
    source_url('https://raw.githubusercontent.com/geofis/rgrass/master/pT.R')
  }
  
  # Groups as integers
  grps <- as.integer(grps)
  
  # Smooth path
  if (smthpath) {
    execGRASS(
      "v.generalize",
      flags='overwrite',
      parameters = list(
        input = path,
        output = 'ttsfpcw_tmp_path_smoothed',
        method = 'reumann',
        threshold = gmeta()$nsres
      )
    )
    path <- 'ttsfpcw_tmp_path_smoothed'
  }
  
  # Smooth centerline
  if (smthcnt) {
    execGRASS(
      "v.generalize",
      flags='overwrite',
      parameters = list(
        input = centerline,
        output = 'ttsfpcw_tmp_centerline_smoothed',
        method = 'reumann',
        threshold = gmeta()$nsres
      )
    )
    centerline <- 'ttsfpcw_tmp_centerline_smoothed'
  }
  
  # Convert the path line(s) to a point vector
  # based on the vertices of the path
  execGRASS(
    "v.to.points",
    flags = c('overwrite','i','t'),
    parameters = list(
      input = path,
      output = 'ttsfpcw_tmp_path_pnt',
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
      map = 'ttsfpcw_tmp_path_pnt',
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
  # to the watershed basin centerline
  execGRASS(
    "v.distance",
    flags = 'overwrite',
    parameters = list(
      from = 'ttsfpcw_tmp_path_pnt',
      from_layer = '2',
      from_type = 'point',
      to = centerline,
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
      from = 'ttsfpcw_tmp_path_pnt',
      from_layer = '2',
      from_type = 'point',
      to = watershed,
      to_type = 'line',
      upload = paste('to_x', 'to_y', sep = ','),
      column = paste('to_basinmargin_x', 'to_basinmargin_y', sep = ',')
    )
  )
  
  # Importing data into R
  distxynear <- readVECT('ttsfpcw_tmp_path_pnt', layer = 2)

  # Trim vertices
  distxynear <- distxynear[
    seq(
      ceiling(quantile(1:nrow(distxynear), ptrim)),
      floor(quantile(1:nrow(distxynear), 1-ptrim)),
      by = 1
    ),
    ]
  
  # XY coords
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
  distxynear@data$seq <- (distxynear@data$cat - 1) %% grps + 1
  
  #Mean value and probabality of change by chance of the magnitude of
  #Transverse Topographic Basin Symmetry (T) for each index number
  mandpTbygrp <- sapply(
    unique(distxynear@data$seq),
    function(x){
      tts <- distxynear@data[distxynear@data$seq == x, 'TTS']
      p <- pT(tts)$p
      m <- pT(tts)$m
      n <- pT(tts)$n
      paste0(
        'Group ', x, ': ',
        'm=', round(m*100, 2), '%, ',
        'p=', round(p, 4), ', ',
        'n=', n, '.\n'
      )
    }
  )
  
  #Overal mean value and probabality of change by chance of the magnitude of
  #Transverse Topographic Basin Symmetry (T)
  ttsov <- distxynear@data$TTS
  pov <- pT(ttsov)$p
  mov <- pT(ttsov)$m
  nov <- pT(ttsov)$n
  mandpTov <- paste0(
    'm=', round(mov*100, 2), '%, ',
    'p=', round(pov, 4), ', ',
    'n=', nov, '.\n'
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
      'ttsfpcw_points_with_magnitude_bearing_entire_sample_',
      suffix
    ),
    v.in.ogr_flags = 'overwrite'
  )
  
  # ##By groups
  # sapply(
  #   unique(distxynear@data$seq),
  #   function(x){
  #     writeVECT(
  #       distxynear[distxynear@data$seq == x,],
  #       paste0(
  #         'ttsfpcw_points_with_magnitude_bearing_group_',
  #         x,
  #         '_',
  #         suffix
  #       ),
  #       v.in.ogr_flags = 'overwrite'
  #     )
  #   }
  # )
  
  #Rename GRASS maps to be saved
  ##The LFP as a point map
  execGRASS(
    "g.rename",
    flags = 'overwrite',
    parameters = list(
      vector = paste0(
        'ttsfpcw_tmp_path_pnt,',
        'ttsfpcw_path_pnt_',
        suffix
      )
    )
  )
  
  ##The smoothed maps
  if (smthpath) {
    execGRASS(
      "g.rename",
      flags = 'overwrite',
      parameters = list(
        vector = paste0(
          path,
          ',',
          'ttsfpcw_path_smoothed',
          suffix
        )
      )
    )
  }
  
  if (smthcnt) {
    execGRASS(
      "g.rename",
      flags = 'overwrite',
      parameters = list(
        vector = paste0(
          centerline,
          ',',
          'ttsfpcw_centerline_smoothed',
          suffix
        )
      )
    )
  }
  
  #Remove tmp maps
  execGRASS(
    "g.remove",
    flags = 'f',
    parameters = list(
      type = c('raster','vector'),
      pattern = 'ttsfpcw_tmp*'
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
