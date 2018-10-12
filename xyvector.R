xyvector <- function(eorb = NULL){
  #Generates a vector of type "xmin xmax ymin ymax"
  #eorb Extent object or bbox object. If bbox object is supplied, columns contain min and max values (in that order), and rows contain x and y (in that order)
  if(class(eorb)=='Extent'){
    z <- c(eorb@xmin, eorb@xmax, eorb@ymin, eorb@ymax) 
  } else if(class(eorb)=='matrix') {
    z <- c(eorb[1,1], eorb[1,2], eorb[2,1], eorb[2,2]) 
  } else {
    print('An object of class Extent or a matrix must be supplied')
  }
  return(z)
}
