PolarPlotforTTS <- function(bearings, tts, maintitle){
  # Generate a polar plot of the Transverse Topographic Basin Symmetry (T) vector
  # Args:
  #   bearings:    One numeric vector with bearing observations. Angles are in radians
  #               between -Pi and Pi. May be generated with the TransTopoSym function.
  #               Positive values measure angles in CCW direction
  #               from the +x axis (horizontal axis pointing East).
  #               Negative values measure angles in CW direction.
  #               Must be the same length as tts.
  #   tts:        One numeric vector with the TTS values
  #   maintitle   Character string. Main title for the plot,
  #               which will apppear in the topmost margin.
  # Returns:
  #   Polar plot of TTS vectors as point geometries in ggplot2 style
  
  # Error handling:
  if(!length(bearings)==length(tts)){
    stop('Bearings and tts must be the same length.')
  }
  
  #Packages:
  require(ggplot2)
  
  #Transform bearings to ggplot2 style
  bearings <- ifelse(bearings > 0, bearings, (2*pi) + bearings)
  mbearings <- mean(bearings, na.rm = T)
  mtts <- mean(tts, na.rm = T)
  
  #Generate the plot
  p <- ggplot() +
    aes(
      x=bearings,
      y=tts) +
    geom_point(alpha = 0.4, size = 3) +
    coord_polar(start = -pi/2, direction = -1) +
    scale_x_continuous(
      labels = c('E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'),
      breaks = pi*(seq(0,7/4,by=1/4)),
      limits = c(0, 2*pi)) +
    scale_y_continuous(labels = seq(0,1, by=0.2),
                       breaks = seq(0,1, by=0.2),
                       limits = c(0,1)) +
    theme(
      legend.position="none",
      text = element_text(size = 22),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'black'), 
      panel.grid = element_line(
        colour = "black",
        linetype = "dashed",
        size = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) +
    annotate(
      'label',
      x = pi/2, y = seq(0,1, by=0.2),
      label = seq(0,1, by=0.2),
      size = 3,
      hjust = 0.5,
      alpha = 0.5
    ) +
    annotate(
      'point',
      x = mbearings, y = mtts, alpha = 0.7,
      shape = 24, size = 5, fill = 'red') +
    annotate(
      'label',
      x = mbearings, y = mtts, alpha = 0.7,
      label = paste0(round(mtts,2)*100, '%, ', round(mbearings*180/pi,0),'°'),
      size = 4, colour = 'red', hjust = -0.1) +
    labs(title= maintitle)
  return(p)
}
