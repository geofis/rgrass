LfpNetworkMerge <- function(main, net.to.merge, suffix){
  # Bring stream networks together in one vector map
  # and assign a unique category to each stream in the cat field
  # Args:
  #   main:             One string with one vector map name
  #   net.to.merge:     One string vector with two or more vector maps names
  #                     to merge to the main network
  # Returns:
  #   GRASS GIS vector map of the merged network
  # Notes:
  #   1) A GRASS session must be initiated using rgrass7 package
  # Error handling
  if (!length(main)==1) {
    stop("Argument main must be a vector with one single value.")
  }
  if (TRUE %in% is.na(main) || TRUE %in% is.na(net.to.merge)) {
    stop("Arguments main and net.to.merge must not have missing values.")
  }
  if (!is.character(main) || !is.character(net.to.merge)) {
    stop("Arguments main and net.to.merge must be character strings.")
  }
  if (!is.character(suffix)) {
    stop("Argument suffix must be a character string.")
  }
  # Packages
  require(rgrass7)
  # Merge main and net.to.merge
  execGRASS(
    'v.patch',
    flags = 'overwrite',
    parameters = list(
      input = c(main, net.to.merge),
      output = paste0('LfpNetwork_lfp_tmp_merged_', suffix)
    )
  )
  execGRASS(
    'v.category',
    flags = 'overwrite',
    parameters = list(
      input = paste0('LfpNetwork_lfp_tmp_merged_', suffix),
      output = paste0('LfpNetwork_lfp_tmp_merged2_', suffix),
      option = 'del',
      cat = -1
    )
  )
  execGRASS(
    'v.clean',
    flags = 'overwrite',
    parameters = list(
      input = paste0('LfpNetwork_lfp_tmp_merged2_', suffix),
      output = paste0('LfpNetwork_lfp_tmp_merged2_clean_', suffix),
      type = 'line',
      tool = 'rmdupl'
    )
  )
  execGRASS(
    'v.category',
    flags = 'overwrite',
    parameters = list(
      input = paste0('LfpNetwork_lfp_tmp_merged2_clean_', suffix),
      output = paste0('LfpNetwork_lfp_all_final_', suffix),
      option = 'add',
      cat = 1
    )
  )
  execGRASS(
    "g.remove",
    flags = 'f',
    parameters = list(
      type = 'vector',
      pattern = 'LfpNetwork_lfp_tmp_*'
    )
  )
}
