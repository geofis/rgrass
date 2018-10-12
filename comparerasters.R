comparerasters <- function(r1, r2, n=100, nsig = 0.05, seed = 131){
  set.seed(seed)
  p1 <- spsample(
    as(extent(r1), 'SpatialPolygons'),
    n,
    type = 'random')
  proj4string(p1) <- CRS(proj4string(r1))
  p2 <- spTransform(p1, CRS = proj4string(r2))
  sr1 <- raster::extract(r1, p1)
  sr2 <- raster::extract(r2, p2)
  ttest <- t.test(sr1, sr2)
  pvalue <- ttest$p.value
  result <- ifelse(
    pvalue>=nsig,
    'there were NO significant differences',
    'there were significant differences'
    )
  RMSE <- format(round(sqrt(mean(sr1-sr2, na.rm = T)^2),2), nsmall = 2)
  return(paste0('RMSE = ', RMSE, ', and t.test suggests that ', result))
}
