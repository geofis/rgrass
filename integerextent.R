intext <- function(e, r, type = c('inner', 'outer')){#Integer numbers for extent object
  #e    Extent object
  #r    Resolution
  #type Inner or outer extent object
  xyfloor <- function(coord, res=r){
    m <- as.integer(floor(coord))
    n <- as.integer(round(res,0))
    i <- 0
    for (i in 0:n){
      if((m-i)%%n==0){
        break
      }
    }
    mu <- m-i
    return(mu)
  }
  xyceiling <- function(coord, res=r){
    m <- as.integer(ceiling(coord))
    n <- as.integer(round(res,0))
    i <- 0
    for (i in 0:n){
      if((m+i)%%n==0){
        break
      }
    }
    mu <- m+i
    return(mu)
  }
  if(type=='inner'){
    outextent <- extent(
      xyceiling(slot(e, 'xmin')),
      xyfloor(slot(e, 'xmax')),
      xyceiling(slot(e, 'ymin')),
      xyfloor(slot(e, 'ymax'))
    )
  } else {
    outextent <- extent(
      xyfloor(slot(e, 'xmin')),
      xyceiling(slot(e, 'xmax')),
      xyfloor(slot(e, 'ymin')),
      xyceiling(slot(e, 'ymax'))
    )
  }
  return(outextent)
}
