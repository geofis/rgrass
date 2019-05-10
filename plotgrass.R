plotgrass <- function(
  gl = NULL, #GRASS GIS layer
  scaledist = 4, #Distance for scale bar (in km)
  qual = F,
  legpos = 'right', #'none' to remove
  cols = scale_fill_viridis()  #Palette for ggplot
  #Other options: scale_fill_gradientn(colours = terrain_hcl(3))
  #Other options: scale_fill_viridis(values = rescale(c(1,5,30,100,1000,40775)), direction = -1)
){
  library(rgrass7)
  library(ggplot2)
  library(ggsn)
  library(raster)
  rlayer <- as.data.frame(
    as(
      raster(
        readRAST(gl)
      ),
      "SpatialPixelsDataFrame")
  )
  datacolname <- colnames(rlayer)[which(!colnames(rlayer) %in% c('x','y'))]
  if(qual) {rlayer[,datacolname] <- as.factor(rlayer[,datacolname])}
  p <- ggplot() +  
    geom_raster(
      data=rlayer,
      aes_string(x = 'x', y = 'y', fill = datacolname),
      interpolate = ifelse(qual, F, T),
      alpha=1) +
    cols +
    coord_equal() +
    theme_bw() +
    theme(legend.position = legpos) +
    labs(title = gl) +
    ggsn::scalebar(
      dist = scaledist,
      st.size=3,
      height=0.02,
      x.min=min(rlayer$x),
      x.max=max(rlayer$x),
      y.min = min(rlayer$y),
      y.max=max(rlayer$y),
      model = 'WGS84') +
    ggsn::north(
      scale = .08,
      x.min=min(rlayer$x),
      x.max=max(rlayer$x),
      y.min = min(rlayer$y),
      y.max=max(rlayer$y))
  return(p)
}
