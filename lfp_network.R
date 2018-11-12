LfpNetwork <- function(xycoords, suffix, stream_vect, direction){
  # Generate the longest flow path of a basin and its tributaries
  # Args:
  #   xycoords:     One vector with the coordinates of the basin outlet
  #   suffix:       One string for the suffix of the GRASS GIS maps to be generated
  #   stream_vect:  One string of the existing stream network in GRASS GIS
  #   direction:    Flow direction raster map. May be generated with r.stream*.
  # Returns:
  #   GRASS GIS maps of the longest flow path and its tributaries
  #   of the basin with outlet at xycoords
  # Note:
  #   A GRASS session must be initiated using rgrass7 package
  # Error handling
  if (!length(xycoords)==2) {
    stop("Argument xycoords must be a vector with two values.")
  }
  if (TRUE %in% is.na(xycoords)) {
    stop("Argument xycoords must not have missing values.")
  }
  if (!is.character(suffix)) {
    stop("Argument suffix must be a character string.")
  }
  if (!is.character(direction)) {
    stop("Argument direction must be a character string.")
  }
  if (!is.character(stream_vect)) {
    stop("Argument stream_vect must be a character string.")
  }
  # Packages
  require(rgrass7)
  require(sp)
  # Generate the main LFP
  execGRASS(
    "r.lfp",
    flags = 'overwrite',
    parameters = list(
      input = direction,
      coordinates = xycoords,
      output = paste0('LfpNetwork_lfp_', suffix)
    )
  )
  # Select tributaries that touches the main LFP
  execGRASS(
    "v.select",
    flags = c('c','overwrite'),
    parameters = list(
      ainput = stream_vect,
      binput = paste0('LfpNetwork_lfp_', suffix),
      output = paste0('LfpNetwork_tributaries_', suffix),
      operator = 'touches'
    )
  )
  # Generate the prelink point of each tributary and write them to GRASS
  trib <- readVECT(paste0('LfpNetwork_tributaries_', suffix))
  tribpts <- as(trib, 'SpatialPointsDataFrame')
  tribpreconf <- tribpts[sapply(
    unique(tribpts@data$Lines.ID),
    function(x){
      n <- nrow(tribpts@data[tribpts@data$Lines.ID==x,])-1
      rownames(tribpts@data[tribpts@data$Lines.ID==x,][n,])
    }),]
  writeVECT(
    tribpreconf,
    paste0('LfpNetwork_tributaries_preconf_', suffix),
    v.in.ogr_flags='overwrite'
  )
  # Generate the LFP of each tributary.
  # In the case of more than one LFP, select one randomly
  sapply(
    tribpreconf$stream,
    function(x){
      execGRASS(
        "r.lfp",
        flags = 'overwrite',
        parameters = list(
          input = direction,
          coordinates = unlist(
            unlist(tribpreconf[tribpreconf@data$stream == x, ]@coords)
          ),
          output = paste0('LfpNetwork_lfp_tmp_', suffix, '_stream_', x)
        )
      )
      if(vInfo(paste0('LfpNetwork_lfp_tmp_', suffix, '_stream_', x))['lines']>1){
        execGRASS(
          'v.category',
          flags = 'overwrite',
          parameters = list(
            input = paste0('LfpNetwork_lfp_tmp_', suffix, '_stream_', x),
            output = paste0(
              'LfpNetwork_lfp_tmp_',
              suffix,
              '_stream_',
              x,
              '_onesinglepath'),
            option = 'del',
            cat = -1
          )
        )
        execGRASS(
          'v.category',
          flags = 'overwrite',
          parameters = list(
            input = paste0(
              'LfpNetwork_lfp_tmp_',
              suffix, '_stream_',
              x,
              '_onesinglepath'),
            output = paste0('LfpNetwork_lfp_tmp_', suffix, '_stream_', x),
            option = 'add',
            cat = 1
          )
        )
        execGRASS(
          'v.extract',
          flags = 'overwrite',
          parameters = list(
            input = paste0('LfpNetwork_lfp_tmp_', suffix, '_stream_', x),
            output = paste0(
              'LfpNetwork_lfp_tmp_',
              suffix,
              '_stream_',
              x,
              '_onesinglepath'),
            random = 1
          )
        )
        execGRASS(
          "g.remove",
          flags = 'f',
          parameters = list(
            type = 'vector',
            name = paste0('LfpNetwork_lfp_tmp_', suffix, '_stream_', x)
          )
        )
      }
    }
  )
  # Merge all the tributaries in a single vector map
  lfptriblist <- execGRASS(
    'g.list',
    parameters = list(
      type = 'vector',
      pattern = 'LfpNetwork_lfp_tmp_*'
    )
  )
  execGRASS(
    'v.patch',
    flags = 'overwrite',
    parameters = list(
      input = attr(lfptriblist, 'resOut'),
      output = paste0('LfpNetwork_lfp_tmp_all_tributaries_', suffix)
    )
  )
  execGRASS(
    'v.patch',
    flags = 'overwrite',
    parameters = list(
      input = paste0(
        paste0('LfpNetwork_lfp_', suffix),
        ',',
        paste0('LfpNetwork_lfp_tmp_all_tributaries_', suffix)
      ),
      output = paste0('LfpNetwork_lfp_all_final_', suffix)
    )
  )
  # Assign unique categories to objects in the cat field
  execGRASS(
    'v.category',
    flags = 'overwrite',
    parameters = list(
      input = paste0('LfpNetwork_lfp_all_final_', suffix),
      output = paste0('LfpNetwork_lfp_tmp_', suffix, 'all_final_'),
      option = 'del',
      cat = -1
    )
  )
  execGRASS(
    'v.category',
    flags = 'overwrite',
    parameters = list(
      input = paste0('LfpNetwork_lfp_tmp_', suffix, 'all_final_'),
      output = paste0('LfpNetwork_lfp_all_final_', suffix),
      option = 'add',
      cat = 1
    )
  )
  # Remove temporary vector maps
  execGRASS(
    "g.remove",
    flags = 'f',
    parameters = list(
      type = 'vector',
      pattern = 'LfpNetwork_lfp_tmp_*'
    )
  )
}
