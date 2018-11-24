TransTopoSym <- function(gtransects, gdivide, gcenterline,
                                    gstream, prefix, ptrim = 0.02){
  # Calculate the Transverse Topographic Basin Symmetry (T)
  # vector given the sources generated from the SourcesForTransTopoSym function
  # Args:
  #   gtransects:     String. Name of the GRASS vector map for the transects.
  #   gdivide:        String. Name of the GRASS vector map for the waterhsed divide.
  #   gcenterline:    String. Name of the GRASS vector map for the basin centerline.
  #   gstream:        String. Name of the GRASS vector map for the longest flow path.
  #   prefix:         One string for the names of the GRASS map.
  #   ptrim:          Numeric. Proportion of transects with respect of the total to
  #                   trim on each endpoint
  # Returns:
  #   Values and GRASS GIS maps of the magnitude and bearing of T vector of a basin,
  #   as well as the centerline of the basin and a smoothed version of the longest flow path
  # Note:
  #   A GRASS session must be initiated using rgrass7 package
  
  # Error handling
  if (!is.character(gtransects) || !is.character(gdivide) || 
      !is.character(gcenterline) || !is.character(gstream)) {
    stop("Arguments gtransects, gdivide, gcenterline and gstream must be character strings.")
  }

  #Packages
  require(rgrass7)
  require(sp)
  require(rgeos)
  require(rgdal)
  require(devtools)
  if (!'pT' %in% ls(envir = .GlobalEnv)) {
    source_url('https://raw.githubusercontent.com/geofis/rgrass/master/pT.R')
  }
  
  #Read the transects vector
  transects <- readVECT(gtransects)
  
  #Read the divide vector
  divide <- readVECT(gdivide)
  
  #Read the centerline vector
  centerline <- readVECT(gcenterline)
  
  #Read the stream or longest flow path
  stream <- readVECT(gstream)
  
  #Rename the slots of...
  ##...the transects SpatialLinesDataFrame
  transects@lines <- sapply(
    transects@lines,
    function(x) {
      x@ID <- paste0('TRANSECT', x@ID); return(x)
    }
  )
  ##...the divide SpatialLinesDataFrame
  divide@lines <- sapply(
    divide@lines,
    function(x) {
      x@ID <- paste0('DIVIDE', x@ID); return(x)
    }
  )
  ##...the centerline SpatialLinesDataFrame
  centerline@lines <- sapply(
    centerline@lines,
    function(x) {
      x@ID <- paste0('CENTERLINE', x@ID); return(x)
    }
  )
  ##...the stream SpatialLinesDataFrame
  stream@lines <- sapply(
    stream@lines,
    function(x) {
      x@ID <- paste0('STREAM', x@ID); return(x)
    }
  )
  
  #Pass the ID to a column. Useless by now, so commented.
  ##This step may be useful if some value column from transect is needed
  ##furtherly, in which case a merge through TRANSECT field is suitable
  # transects$TRANSECT <- sapply(transects@lines, function(x) x@ID)
  
  
  #Intersection of transects and stream
  ##Create the SpatialPointsDataFrame (SPDF)
  transectsstream <- gIntersection(transects, stream, byid = T)
  ##Save the rownames, which contains the attributes of the intersected lines
  IDts <- rownames(transectsstream@coords)
  ##Remove rownames, which prevents a warning message in further steps
  rownames(transectsstream@coords) <- c()
  ##Create a SPDF form the previously created SPPF
  transectsstream2 <- SpatialPointsDataFrame(
    transectsstream,
    data.frame(
      UNIQUE = 1:length(transectsstream),
      ID = IDts,
      x = as.vector(transectsstream@coords[,'x']),
      y = as.vector(transectsstream@coords[,'y'])
    )
  )
  ##Dummy object containing the points of the intersection
  foo3 <- data.frame(
    UNIQUE = transectsstream2$UNIQUE,
    stringr::str_split_fixed(transectsstream2$ID, ' ', 2)
  )
  colnames(foo3)[2:3] <- c('TRANSECT','STREAM')
  transectsstream2 <- merge(transectsstream2, foo3, by = 'UNIQUE', sort = F)
  transectsstream2@data[,c('x', 'y', 'TRANSECT')]
  
  
  #Intersection of transects and divide
  ##Create the SpatialPointsDataFrame (SPDF)
  transectsdivide <- gIntersection(transects, divide, byid = T)
  ##Save the rownames, which contains the attributes of the intersected lines
  IDtd <- rownames(transectsdivide@coords)
  ##Remove rownames, which prevents a warning message in further steps
  rownames(transectsdivide@coords) <- c()
  ##Create a SPDF form the previously created SPPF
  transectsdivide2 <- SpatialPointsDataFrame(
    transectsdivide,
    data.frame(
      UNIQUE = 1:length(transectsdivide),
      ID = IDtd,
      x = as.vector(transectsdivide@coords[,'x']),
      y = as.vector(transectsdivide@coords[,'y'])
    )
  )
  ##Dummy object containing the points of the intersection
  foo1 <- data.frame(
    UNIQUE = transectsdivide2$UNIQUE,
    stringr::str_split_fixed(transectsdivide2$ID, ' ', 2)
  )
  ##Update the colnames of the dummy object
  colnames(foo1)[2:3] <- c('TRANSECT','DIVIDE')
  ##Merge the dummy object with the data frame of the slot of SPDF object
  transectsdivide2 <- merge(transectsdivide2, foo1, by = 'UNIQUE', sort = F)
  ##Calculate distance between stream points and divide points
  
  #Calculate distance between divide and stream points
  distdivide <- sapply(
    as.character(transectsstream2$TRANSECT),
    function(x) {
      xytstream <- transectsstream2@data[transectsstream2$TRANSECT==x,]
      xytdivide <- transectsdivide2@data[transectsdivide2$TRANSECT==x,]
      dist <- min(
        sqrt((xytstream$x-xytdivide$x)^2 + (xytstream$y-xytdivide$y)^2),
        na.rm = T)
      return(dist)
    }
  )
  distdivide2 <- data.frame(
    DIVIDESTREAMDIST = as.vector(distdivide, mode = 'numeric'),
    TRANSECT = names(distdivide))
  
  #Merge results with SPDF
  transectsstream3 <- merge(transectsstream2, distdivide2, by = 'TRANSECT', sort = F)
  

  
  #Intersection of transects and centerline
  ##Create the SpatialPointsDataFrame (SPDF)
  transectscenterline <- gIntersection(transects, centerline, byid = T)
  ##Save the rownames, which contains the attributes of the intersected lines
  IDtc <- rownames(transectscenterline@coords)
  ##Remove rownames, which prevents a warning message in further steps
  rownames(transectscenterline@coords) <- c()
  ##Create a SPDF form the previously created SPPF
  transectscenterline2 <- SpatialPointsDataFrame(
    transectscenterline,
    data.frame(
      UNIQUE = 1:length(transectscenterline),
      ID = IDtc,
      x = as.vector(transectscenterline@coords[,'x']),
      y = as.vector(transectscenterline@coords[,'y'])
    )
  )
  ##Dummy object containing the points of the intersection
  foo1 <- data.frame(
    UNIQUE = transectscenterline2$UNIQUE,
    stringr::str_split_fixed(transectscenterline2$ID, ' ', 2)
  )
  ##Update the colnames of the dummy object
  colnames(foo1)[2:3] <- c('TRANSECT','CENTERLINE')
  ##Merge the dummy object with the data frame of the slot of SPDF object
  transectscenterline2 <- merge(transectscenterline2, foo1, by = 'UNIQUE', sort = F)
  
  #Calculate distance between centerline and stream points
  distcenterline <- sapply(
    as.character(transectsstream2$TRANSECT),
    function(x) {
      xytstream <- transectsstream2@data[transectsstream2$TRANSECT==x,]
      xytcenterline <- transectscenterline2@data[transectscenterline2$TRANSECT==x,]
      dist <- sqrt((xytstream$x-xytcenterline$x)^2 + (xytstream$y-xytcenterline$y)^2)
      return(dist)
    }
  )
  distcenterline2 <- data.frame(
    Da_CENTERLINESTREAMDIST = as.vector(distcenterline, mode = 'numeric'),
    TRANSECT = names(distcenterline))

  #Merge results with SPDF
  transectsstream4 <- merge(transectsstream3, distcenterline2, by = 'TRANSECT', sort = F)
  
  #Calculate Dd
  transectsstream4$Dd <- transectsstream4$Da_CENTERLINESTREAMDIST / 
    (transectsstream4$Da_CENTERLINESTREAMDIST + transectsstream4$DIVIDESTREAMDIST)
  
  #Bearing between stream points and centerline points
  bearing <- sapply(
    as.character(transectsstream2$TRANSECT),
    function(x) {
      xytstream <- transectsstream2@data[transectsstream2$TRANSECT==x,]
      xytcenterline <- transectscenterline2@data[transectscenterline2$TRANSECT==x,]
      bearing <- atan2(xytcenterline$y-xytstream$y, xytcenterline$x-xytstream$x)
      return(bearing)
    }
  )
  bearing2 <- data.frame(
    BEARINGRAD = as.vector(bearing, mode = 'numeric'),
    TRANSECT = names(bearing))
  
  #Merge results with SPDF
  transectsstream5 <- merge(transectsstream4, bearing2, by = 'TRANSECT', sort = F)
  
  transectsstream5$BEARINGDEG <- with(
    transectsstream5@data,
    BEARINGRAD * 180 / pi
  )
  
  #Trim the sample
  transectsstream6 <- transectsstream5[seq(
    ceiling(quantile(1:nrow(transectsstream5), ptrim)),
    floor(quantile(1:nrow(transectsstream5), 1 - ptrim)),
    by = 1
  ),]
  
  #Mean value and probabality of change by chance of the magnitude of
  #Transverse Topographic Basin Symmetry (T) 
  mandpT <- pT(transectsstream6$Dd)
  
  #Write to GRASS database
  # writeVECT(
  #   transectsstream6,
  #   paste0('tts_', prefix, '_points'),
  #   v.in.ogr_flags = 'overwrite'
  #   )
  
  #Write to SHP
  writeOGR(
    transectsstream6,
    dsn = paste0(wd, '/TTSresults/'),
    layer = paste0('tts_', prefix, '_points'),
    driver = 'ESRI Shapefile',
    overwrite_layer = T
  )
  
  #Returns:
  return(
    list(
      mandp = mandpT,
      SpatialObjectWithAttributes = transectsstream6
    )
  )
}
