pT <- function(tts){
  # Probabality of magnitude of Transverse Topographic Basin Symmetry (T)
  # Description
  # Calculates the probability (p) of Ho: "obtaining a greater mean vector magnitude
  # by pure chance combination of random vectors", proposed
  # by Cox (1994) based on Curray (1956)
  # Args:
  #   tts:          Vector of Transverse Topographic Basin Symmetry (T) values
  
  # Returns:
  #   The probability value
  # Error handling
  if (length(tts) < 5) {
    stop("Argument tts must be a vector with at least 5 values")
  }
  if (!is.numeric(tts)) {
    stop("Argument tts must be a numeric vector")
  }
  p <- exp(
    -((mean(tts)*100)^2)*
      length(tts)*
      (1e-04)
  )
  return(p)
}
